"""Utility helpers shared by Cellpose workflow notebooks.

This module keeps notebook cells simple and beginner-friendly by moving
implementation details into reusable, documented Python functions.

The helpers are intentionally focused on the two notebooks:
- notebooks/part_1_segmentation_workflows/NEW_infer_masks_cellpose.ipynb
- notebooks/part_1_segmentation_workflows/NEW_cellpose-prep.ipynb
"""

from __future__ import annotations

import json
import threading
import time
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
import torch
from cellpose import models
from skimage.transform import resize

from infer_subc.organelles.membrane import membrane_composite


REQUIRED_CONFIG_KEYS: tuple[str, ...] = (
    "IM_TYPE",
    "SEGMENT_CELLS",
    "SEGMENT_NUCLEI",
    "CELL_SOURCE",
    "PM_CH",
    "NUC_CH",
    "USE_NUCLEI_FOR_CELLS",
    "XY_DOWNSAMPLE",
    "CELL_DIAMETER",
    "NUC_DIAMETER",
    "FLOW_THRESHOLD",
    "CELLPROB_THRESHOLD",
    "MIN_SIZE_CELL",
    "MIN_SIZE_NUC",
    "BATCH_SIZE",
    "USE_GPU",
)


def ensure_cellpose_v3_models_available() -> str:
    """Validate that the current Cellpose version still includes `cyto3`/`nuclei`.

    Returns:
        Installed Cellpose version string.

    Raises:
        ImportError: If Cellpose major version is 4 or above.
    """
    cp_version = version("cellpose")
    if int(cp_version.split(".")[0]) >= 4:
        raise ImportError(
            'Cellpose >= 4 no longer includes the cyto3/nuclei models. '
            'Run: pip install "cellpose<4" and restart the kernel.'
        )
    return cp_version


def prepare_volume_zyx(volume: np.ndarray, xy_downsample: int = 1) -> np.ndarray:
    """Scale a 3D volume to [0, 1] and optionally shrink XY.

    Args:
        volume: Input image volume in ZYX order.
        xy_downsample: XY shrink factor. `1` keeps original size.

    Returns:
        Float32 volume in ZYX order.
    """
    vol = np.asarray(volume, dtype=np.float32)
    p1, p99 = np.percentile(vol, (1, 99.8))
    vol = np.clip((vol - p1) / (p99 - p1 + 1e-8), 0, 1)

    if int(xy_downsample) != 1:
        vol = resize(
            vol,
            (
                vol.shape[0],
                round(vol.shape[1] / int(xy_downsample)),
                round(vol.shape[2] / int(xy_downsample)),
            ),
            order=1,
            anti_aliasing=True,
        ).astype(np.float32)

    return vol


def restore_label_volume(labels: np.ndarray, output_shape_zyx: Sequence[int]) -> np.ndarray:
    """Resize a label volume back to target shape with nearest-neighbor interpolation.

    Args:
        labels: Input labels in ZYX order.
        output_shape_zyx: Target ZYX shape.

    Returns:
        Restored labels as `uint16`.
    """
    restored = labels
    output_shape = tuple(int(v) for v in output_shape_zyx)
    if tuple(labels.shape) != output_shape:
        restored = resize(
            labels,
            output_shape,
            order=0,
            anti_aliasing=False,
            preserve_range=True,
        )
    return restored.astype(np.uint16)


def get_process_rss_gib() -> float:
    """Read current process RSS memory in GiB (Linux `/proc` path).

    Returns:
        RSS in GiB, or `nan` if unavailable.
    """
    try:
        with open("/proc/self/status", "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    kb = float(line.split()[1])
                    return kb / (1024**2)
    except Exception:
        pass
    return float("nan")


def start_resource_monitor(
    *,
    enabled: bool,
    sample_sec: float,
    gpu_enabled: bool,
) -> tuple[threading.Event | None, threading.Thread | None, dict[str, Any]]:
    """Start background RAM/GPU sampling for one inference run.

    Args:
        enabled: If `False`, sampling is skipped and defaults are returned.
        sample_sec: Sampling period in seconds.
        gpu_enabled: Whether to sample CUDA memory values.

    Returns:
        Tuple `(stop_event, thread, stats)`.

        - `stop_event`/`thread` are `None` when monitoring is disabled.
        - `stats` contains `ram_gb_peak`, `gpu_alloc_gb_peak`,
          `gpu_reserved_gb_peak`, and `gpu_device`.
    """
    stats: dict[str, Any] = {
        "ram_gb_peak": get_process_rss_gib(),
        "gpu_alloc_gb_peak": 0.0,
        "gpu_reserved_gb_peak": 0.0,
        "gpu_device": None,
    }

    if not enabled:
        return None, None, stats

    if gpu_enabled:
        torch.cuda.reset_peak_memory_stats()
        stats["gpu_device"] = torch.cuda.get_device_name(torch.cuda.current_device())

    stop_event = threading.Event()

    def _worker() -> None:
        while not stop_event.is_set():
            rss = get_process_rss_gib()
            if not np.isnan(rss):
                stats["ram_gb_peak"] = max(stats["ram_gb_peak"], rss)

            if gpu_enabled:
                alloc = torch.cuda.memory_allocated() / (1024**3)
                reserv = torch.cuda.memory_reserved() / (1024**3)
                stats["gpu_alloc_gb_peak"] = max(stats["gpu_alloc_gb_peak"], alloc)
                stats["gpu_reserved_gb_peak"] = max(stats["gpu_reserved_gb_peak"], reserv)

            time.sleep(max(0.05, float(sample_sec)))

    thread = threading.Thread(target=_worker, daemon=True)
    thread.start()
    return stop_event, thread, stats


def stop_resource_monitor(
    *,
    stop_event: threading.Event | None,
    thread: threading.Thread | None,
    stats: dict[str, Any],
    gpu_enabled: bool,
) -> None:
    """Stop background resource monitoring and finalize peak GPU stats.

    Args:
        stop_event: Event returned from `start_resource_monitor`.
        thread: Thread returned from `start_resource_monitor`.
        stats: Mutable stats dictionary returned from `start_resource_monitor`.
        gpu_enabled: Whether CUDA memory peaks should be finalized.
    """
    if stop_event is not None:
        stop_event.set()
    if thread is not None:
        thread.join(timeout=2)

    if gpu_enabled:
        stats["gpu_alloc_gb_peak"] = max(
            stats["gpu_alloc_gb_peak"], torch.cuda.max_memory_allocated() / (1024**3)
        )
        stats["gpu_reserved_gb_peak"] = max(
            stats["gpu_reserved_gb_peak"], torch.cuda.max_memory_reserved() / (1024**3)
        )


def to_jsonable(value: Any) -> Any:
    """Recursively convert common scientific-Python values to JSON-safe values.

    Args:
        value: Any scalar, array, sequence, or mapping.

    Returns:
        JSON-serializable representation where practical.
    """
    if value is None:
        return None

    if isinstance(value, (np.integer, np.floating, np.bool_)):
        value = value.item()

    if isinstance(value, np.ndarray):
        return [to_jsonable(x) for x in value.tolist()]

    if isinstance(value, (list, tuple)):
        return [to_jsonable(x) for x in value]

    if isinstance(value, dict):
        return {str(k): to_jsonable(v) for k, v in value.items()}

    try:
        if pd.isna(value):
            return None
    except Exception:
        pass

    return value


def append_profile_record_upsert(
    log_dir: Path,
    record: Mapping[str, Any],
    upsert_keys: Sequence[str] | None = None,
) -> tuple[Path, Path]:
    """Upsert one benchmark record into CSV and JSONL logs.

    Existing rows with the same values for every key in `upsert_keys` are removed,
    then the new row is appended. This keeps logs compact while preserving the
    latest run per dataset (or other chosen key set).

    Args:
        log_dir: Output directory for log files.
        record: Record to append.
        upsert_keys: Key columns used to identify replaceable old rows.

    Returns:
        Tuple `(csv_path, jsonl_path)`.
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    csv_path = log_dir / "cellpose_resource_log.csv"
    jsonl_path = log_dir / "cellpose_resource_log.jsonl"

    keys = list(upsert_keys) if upsert_keys else ["dataset_id"]
    new_row = pd.DataFrame([dict(record)])

    if csv_path.exists():
        existing = pd.read_csv(csv_path)
        for key in keys:
            if key not in existing.columns:
                existing[key] = None
            if key not in new_row.columns:
                new_row[key] = None

        keep_mask = np.ones(len(existing), dtype=bool)
        for key in keys:
            keep_mask &= existing[key].astype(str).values != str(new_row.iloc[0][key])

        updated = pd.concat([existing.loc[keep_mask], new_row], ignore_index=True)
    else:
        updated = new_row

    updated.to_csv(csv_path, index=False)

    with open(jsonl_path, "w", encoding="utf-8") as fh:
        for _, row in updated.iterrows():
            payload = {k: to_jsonable(v) for k, v in row.to_dict().items()}
            fh.write(json.dumps(payload) + "\n")

    return csv_path, jsonl_path


def assign_nuclei_to_cells(nuclei_labels: np.ndarray, cell_labels: np.ndarray) -> np.ndarray:
    """Map each nucleus label to the most-overlapping cell label.

    This keeps only nuclei that overlap with a cell and relabels them so nucleus
    IDs match their paired cell IDs.

    Args:
        nuclei_labels: Nuclei instance labels in ZYX order.
        cell_labels: Cell instance labels in ZYX order.

    Returns:
        Nucleus label volume where each kept nucleus has its paired cell label.
    """
    mapped = np.zeros_like(nuclei_labels, dtype=np.uint16)
    inside = (nuclei_labels > 0) & (cell_labels > 0)
    if not inside.any():
        return mapped

    pairs, counts = np.unique(
        np.stack([nuclei_labels[inside], cell_labels[inside]]),
        axis=1,
        return_counts=True,
    )
    best_cell: dict[int, int] = {}
    best_count: dict[int, int] = {}
    for (nuc_id, cell_id), overlap_count in zip(pairs.T, counts):
        if int(overlap_count) > best_count.get(int(nuc_id), 0):
            best_cell[int(nuc_id)] = int(cell_id)
            best_count[int(nuc_id)] = int(overlap_count)

    lut = np.zeros(int(nuclei_labels.max()) + 1, dtype=np.uint16)
    for nuc_id, cell_id in best_cell.items():
        lut[nuc_id] = np.uint16(cell_id)

    return lut[nuclei_labels]


def run_cellpose_single_image(
    *,
    img_data: np.ndarray,
    scale_zyx: Sequence[float],
    use_composite: bool,
    composite_zyx: np.ndarray | None,
    pm_channel: int,
    nuc_channel: int | None,
    segment_nuclei: bool,
    use_nuclei_for_cells: bool,
    xy_downsample: int,
    cell_diameter: float,
    nuc_diameter: float,
    flow_threshold: float,
    cellprob_threshold: float,
    min_size_cell: int,
    min_size_nuc: int,
    batch_size: int,
    cell_model_path: str | None,
    nuc_model_path: str | None,
    profile_resources: bool,
    profile_sample_sec: float,
) -> dict[str, Any]:
    """Run one Cellpose 3D segmentation pass for cells and optionally nuclei.

    Args:
        img_data: Multi-channel image in CZYX order.
        scale_zyx: Physical voxel spacing `(z, y, x)`.
        use_composite: Whether to use `composite_zyx` instead of one PM channel.
        composite_zyx: Pre-built composite in ZYX order when `use_composite=True`.
        pm_channel: PM channel index in `img_data`.
        nuc_channel: Nuclei channel index in `img_data` when needed.
        segment_nuclei: If `True`, run the nuclei model.
        use_nuclei_for_cells: If `True`, pass nuclei as helper channel to `cyto3`.
        xy_downsample: XY downsample factor before inference.
        cell_diameter: Cell diameter for Cellpose.
        nuc_diameter: Nucleus diameter for Cellpose.
        flow_threshold: Cellpose flow threshold.
        cellprob_threshold: Cellpose probability threshold.
        min_size_cell: Minimum accepted cell size.
        min_size_nuc: Minimum accepted nucleus size.
        batch_size: Cellpose batch size.
        cell_model_path: Optional custom cell model path.
        nuc_model_path: Optional custom nucleus model path.
        profile_resources: If `True`, track peak RAM/GPU during model calls.
        profile_sample_sec: Sampling interval for resource monitor.

    Returns:
        Dictionary with segmentation outputs and runtime stats:
            - `cellmask_obj`
            - `nuc_all`
            - `nuc_obj`
            - `elapsed_sec`
            - `resource_stats`
            - `anisotropy`
            - `pm_vol_shape`
            - `use_gpu`
            - `need_nuc_input`
    """
    ensure_cellpose_v3_models_available()

    use_gpu = torch.cuda.is_available()
    use_nuc = bool(use_nuclei_for_cells) and bool(segment_nuclei)
    need_nuc_input = bool(segment_nuclei) or use_nuc

    if use_composite:
        if composite_zyx is None:
            raise ValueError("composite_zyx must be provided when use_composite=True")
        pm_vol = prepare_volume_zyx(composite_zyx, xy_downsample)
    else:
        pm_vol = prepare_volume_zyx(img_data[int(pm_channel)], xy_downsample)

    if need_nuc_input:
        if nuc_channel is None:
            raise ValueError("nuc_channel must be set when nuclei are required")
        nuc_vol = prepare_volume_zyx(img_data[int(nuc_channel)], xy_downsample)
    else:
        nuc_vol = None

    output_shape = img_data[int(pm_channel)].shape
    anisotropy = float(scale_zyx[0]) / (float(scale_zyx[1]) * int(xy_downsample))

    common = dict(
        do_3D=True,
        z_axis=0,
        anisotropy=anisotropy,
        flow_threshold=float(flow_threshold),
        cellprob_threshold=float(cellprob_threshold),
        batch_size=int(batch_size),
    )

    t0 = time.perf_counter()
    stop_event, monitor_thread, resource_stats = start_resource_monitor(
        enabled=bool(profile_resources),
        sample_sec=float(profile_sample_sec),
        gpu_enabled=use_gpu,
    )

    nuc_raw = None
    cell_raw = None
    try:
        if segment_nuclei:
            if nuc_model_path:
                nuc_model = models.CellposeModel(gpu=use_gpu, pretrained_model=str(nuc_model_path))
                nuc_diam_arg = None
            else:
                nuc_model = models.CellposeModel(gpu=use_gpu, model_type="nuclei")
                nuc_diam_arg = float(nuc_diameter)

            nuc_raw = nuc_model.eval(
                nuc_vol,
                diameter=nuc_diam_arg,
                channels=[0, 0],
                min_size=int(min_size_nuc),
                **common,
            )[0]

        if cell_model_path:
            cell_model = models.CellposeModel(gpu=use_gpu, pretrained_model=str(cell_model_path))
            cell_diam_arg = None
        else:
            cell_model = models.CellposeModel(gpu=use_gpu, model_type="cyto3")
            cell_diam_arg = float(cell_diameter)

        if use_nuc:
            if nuc_vol is None:
                raise ValueError("nuc_vol is required when use_nuclei_for_cells=True")
            cell_input = np.stack([pm_vol, nuc_vol], axis=-1)
            cell_raw = cell_model.eval(
                cell_input,
                diameter=cell_diam_arg,
                channels=[1, 2],
                channel_axis=3,
                min_size=int(min_size_cell),
                **common,
            )[0]
        else:
            cell_raw = cell_model.eval(
                pm_vol,
                diameter=cell_diam_arg,
                channels=[0, 0],
                min_size=int(min_size_cell),
                **common,
            )[0]

    finally:
        stop_resource_monitor(
            stop_event=stop_event,
            thread=monitor_thread,
            stats=resource_stats,
            gpu_enabled=use_gpu,
        )

    elapsed_sec = time.perf_counter() - t0

    cellmask_obj = restore_label_volume(cell_raw, output_shape)

    if segment_nuclei and nuc_raw is not None:
        nuc_all = restore_label_volume(nuc_raw, output_shape)
        nuc_obj = assign_nuclei_to_cells(nuc_all, cellmask_obj)
    else:
        nuc_all = None
        nuc_obj = None

    return {
        "cellmask_obj": cellmask_obj,
        "nuc_all": nuc_all,
        "nuc_obj": nuc_obj,
        "elapsed_sec": float(elapsed_sec),
        "resource_stats": resource_stats,
        "anisotropy": float(anisotropy),
        "pm_vol_shape": list(pm_vol.shape),
        "use_gpu": bool(use_gpu),
        "need_nuc_input": bool(need_nuc_input),
        "effective_use_nuclei_for_cells": bool(use_nuc),
    }


def canonical_settings_keys(raw: Mapping[str, Any]) -> dict[str, Any]:
    """Allow config payloads to be either flat or nested under `settings`.

    Args:
        raw: Loaded JSON object from config file.

    Returns:
        Flat dictionary of settings and selected metadata fields.

    Raises:
        ValueError: If input is not a mapping.
    """
    if not isinstance(raw, Mapping):
        raise ValueError("Config content must be a JSON object.")

    if "settings" in raw and isinstance(raw["settings"], Mapping):
        merged = dict(raw["settings"])
        for key in ("dataset_id", "selected_on", "selected_by", "notes", "config_version"):
            if key in raw and key not in merged:
                merged[key] = raw[key]
        return merged

    return dict(raw)


def validate_and_normalize_settings(
    settings: Mapping[str, Any],
    *,
    default_settings: Mapping[str, Any],
    required_keys: Sequence[str] = REQUIRED_CONFIG_KEYS,
    source_label: str = "settings",
) -> dict[str, Any]:
    """Merge defaults with user settings, type-cast values, and validate constraints.

    Args:
        settings: Incoming settings from config or notebook inputs.
        default_settings: Fallback values for missing fields.
        required_keys: Keys that must exist after merge.
        source_label: Label included in validation error messages.

    Returns:
        Normalized settings dictionary.

    Raises:
        ValueError: If settings are invalid.
    """
    merged = dict(default_settings)
    merged.update(dict(settings))

    if merged.get("CELL_SOURCE") not in ("composite", "single"):
        raise ValueError(f"{source_label}: CELL_SOURCE must be 'composite' or 'single'.")

    merged["SEGMENT_CELLS"] = bool(merged["SEGMENT_CELLS"])
    merged["SEGMENT_NUCLEI"] = bool(merged["SEGMENT_NUCLEI"])
    merged["USE_NUCLEI_FOR_CELLS"] = bool(merged["USE_NUCLEI_FOR_CELLS"])
    merged["USE_GPU"] = bool(merged.get("USE_GPU", True))

    merged["PM_CH"] = int(merged["PM_CH"])
    merged["NUC_CH"] = int(merged["NUC_CH"])
    merged["XY_DOWNSAMPLE"] = int(merged.get("XY_DOWNSAMPLE", 1))
    merged["CELL_DIAMETER"] = int(merged.get("CELL_DIAMETER", 0))
    merged["NUC_DIAMETER"] = int(merged.get("NUC_DIAMETER", 0))
    merged["MIN_SIZE_CELL"] = int(merged.get("MIN_SIZE_CELL", 0))
    merged["MIN_SIZE_NUC"] = int(merged.get("MIN_SIZE_NUC", 0))
    merged["BATCH_SIZE"] = int(merged.get("BATCH_SIZE", 4))

    merged["FLOW_THRESHOLD"] = float(merged.get("FLOW_THRESHOLD", 0.2))
    merged["CELLPROB_THRESHOLD"] = float(merged.get("CELLPROB_THRESHOLD", -3))

    upsert_keys = merged.get("PROFILE_UPSERT_KEYS", ["dataset_id"])
    if isinstance(upsert_keys, str):
        upsert_keys = [upsert_keys]
    if not isinstance(upsert_keys, (list, tuple)) or len(upsert_keys) == 0:
        upsert_keys = ["dataset_id"]
    merged["PROFILE_UPSERT_KEYS"] = [str(v) for v in upsert_keys]

    if merged["XY_DOWNSAMPLE"] < 1:
        raise ValueError(f"{source_label}: XY_DOWNSAMPLE must be >= 1.")

    if merged["CELL_SOURCE"] == "composite":
        weights = merged.get("COMPOSITE_WEIGHTS")
        if not isinstance(weights, Mapping) or not any(float(w) > 0 for w in weights.values()):
            raise ValueError(
                f"{source_label}: COMPOSITE_WEIGHTS must include at least one weight > 0."
            )

    missing = [key for key in required_keys if key not in merged]
    if missing:
        raise ValueError(f"{source_label}: missing required keys: {missing}")

    return merged


def write_dataset_config(
    dataset_folder: Path,
    *,
    config_filename: str,
    settings: Mapping[str, Any],
    config_folder: Path | None = None,
    dataset_id: str | None = None,
    selected_by: str,
    notes: str,
    generated_by: str,
    generation_reason: str,
    write_history: bool,
    history_dirname: str,
) -> tuple[Path, Path | None]:
    """Write active dataset config JSON and optional timestamped history snapshot.

    Args:
        dataset_folder: Dataset directory path.
        config_filename: Active config filename.
        settings: Normalized settings payload.
        config_folder: Optional folder where config/history files are written.
            If omitted, `dataset_folder` is used.
        dataset_id: Optional explicit dataset ID for metadata.
            If omitted, `dataset_folder.name` is used.
        selected_by: User or process label.
        notes: Optional free-text notes.
        generated_by: Source notebook/script name.
        generation_reason: Short reason label.
        write_history: If `True`, also write a timestamped snapshot.
        history_dirname: Name of history subfolder.

    Returns:
        Tuple `(active_path, history_path)` where `history_path` may be `None`.
    """
    dataset_folder = Path(dataset_folder)
    target_folder = Path(config_folder) if config_folder is not None else dataset_folder
    target_folder.mkdir(parents=True, exist_ok=True)

    payload = {
        "dataset_id": dataset_id or dataset_folder.name,
        "config_version": "1.0",
        "selected_on": datetime.now().isoformat(timespec="seconds"),
        "selected_by": selected_by,
        "notes": notes,
        "generated_by": generated_by,
        "generation_reason": generation_reason,
        "settings": dict(settings),
    }

    text = json.dumps(payload, indent=2)

    active_path = target_folder / config_filename
    active_path.write_text(text)

    history_path: Path | None = None
    if write_history:
        hist_dir = target_folder / history_dirname
        hist_dir.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        history_path = hist_dir / f"cellpose_config_{stamp}.json"
        history_path.write_text(text)

    return active_path, history_path


def update_config_im_type(config_path: Path, image_suffix: str) -> None:
    """Update the `IM_TYPE` value in an existing dataset config file.

    This is used after extension auto-detection when the configured extension
    no longer matches files on disk.

    Args:
        config_path: Path to active `cellpose_config.json`.
        image_suffix: File suffix to store (example: `.tiff`).
    """
    raw = json.loads(Path(config_path).read_text())
    if "settings" in raw and isinstance(raw["settings"], dict):
        raw["settings"]["IM_TYPE"] = image_suffix
    else:
        raw["IM_TYPE"] = image_suffix
    Path(config_path).write_text(json.dumps(raw, indent=2))


def resolve_dataset_settings(
    dataset_folder: Path,
    *,
    config_folder: Path | None = None,
    settings_source: str,
    config_filename: str,
    strict_config: bool,
    fallback_to_manual_if_missing_config: bool,
    default_settings: Mapping[str, Any],
    manual_settings: Mapping[str, Any],
    selected_by: str,
    notes: str,
    write_config_history: bool,
    history_dirname: str,
    generated_by: str,
) -> tuple[dict[str, Any], Path, str]:
    """Resolve settings for one dataset from config or manual defaults.

    Args:
        dataset_folder: Dataset folder.
        config_folder: Optional folder where config/history files are read/written.
            If omitted, `dataset_folder` is used.
        settings_source: Either `config` or `manual`.
        config_filename: Config file name in each dataset folder.
        strict_config: If `True`, missing config in config mode is an error.
        fallback_to_manual_if_missing_config: In config mode, create config from
            manual settings when file is missing.
        default_settings: Default values for missing keys.
        manual_settings: Notebook manual inputs.
        selected_by: Value for config metadata.
        notes: Value for config metadata.
        write_config_history: Whether to store timestamped snapshots.
        history_dirname: Snapshot folder name.
        generated_by: Notebook/script name for metadata.

    Returns:
        `(settings, config_path, source_label)`.

    Raises:
        ValueError: For invalid mode or invalid settings.
        FileNotFoundError: If strict config mode requires a missing config.
    """
    dataset_folder = Path(dataset_folder)
    config_base = Path(config_folder) if config_folder is not None else dataset_folder
    config_base.mkdir(parents=True, exist_ok=True)

    mode = str(settings_source).lower().strip()
    config_path = config_base / config_filename

    if mode == "manual":
        settings = validate_and_normalize_settings(
            manual_settings,
            default_settings=default_settings,
            source_label="manual inputs",
        )
        config_path, _ = write_dataset_config(
            dataset_folder,
            config_filename=config_filename,
            settings=settings,
            config_folder=config_base,
            dataset_id=dataset_folder.name,
            selected_by=selected_by,
            notes=notes,
            generated_by=generated_by,
            generation_reason="manual-settings",
            write_history=write_config_history,
            history_dirname=history_dirname,
        )
        return settings, config_path, "manual"

    if mode != "config":
        raise ValueError("SETTINGS_SOURCE must be 'config' or 'manual'.")

    if config_path.exists():
        loaded = canonical_settings_keys(json.loads(config_path.read_text()))
        settings = validate_and_normalize_settings(
            loaded,
            default_settings=default_settings,
            source_label=str(config_path),
        )
        return settings, config_path, "config-file"

    if fallback_to_manual_if_missing_config:
        settings = validate_and_normalize_settings(
            manual_settings,
            default_settings=default_settings,
            source_label=f"manual fallback for {dataset_folder}",
        )
        config_path, _ = write_dataset_config(
            dataset_folder,
            config_filename=config_filename,
            settings=settings,
            config_folder=config_base,
            dataset_id=dataset_folder.name,
            selected_by=selected_by,
            notes=notes,
            generated_by=generated_by,
            generation_reason="missing-config-manual-fallback",
            write_history=write_config_history,
            history_dirname=history_dirname,
        )
        return settings, config_path, "manual-fallback"

    if strict_config:
        raise FileNotFoundError(f"Missing {config_filename} in {dataset_folder}")

    settings = validate_and_normalize_settings(
        default_settings,
        default_settings=default_settings,
        source_label=f"default fallback for {dataset_folder}",
    )
    config_path, _ = write_dataset_config(
        dataset_folder,
        config_filename=config_filename,
        settings=settings,
        config_folder=config_base,
        dataset_id=dataset_folder.name,
        selected_by=selected_by,
        notes=notes,
        generated_by=generated_by,
        generation_reason="missing-config-default-fallback",
        write_history=write_config_history,
        history_dirname=history_dirname,
    )
    return settings, config_path, "default-fallback"


def normalize_volume_01(volume: np.ndarray) -> np.ndarray:
    """Normalize any numeric volume to [0, 1] using robust percentiles.

    Args:
        volume: Input image volume.

    Returns:
        Float32 normalized volume.
    """
    v = np.asarray(volume, dtype=np.float32)
    p1, p99 = np.percentile(v, (1, 99.8))
    return np.clip((v - p1) / (p99 - p1 + 1e-8), 0, 1).astype(np.float32)


def downsample_zyx(volume: np.ndarray, xy_downsample: int) -> np.ndarray:
    """Downsample a ZYX volume in XY while keeping Z unchanged.

    Args:
        volume: Input volume in ZYX order.
        xy_downsample: Downsample factor (>=1).

    Returns:
        Downsampled volume.
    """
    if int(xy_downsample) == 1:
        return volume
    shape = (
        volume.shape[0],
        round(volume.shape[1] / int(xy_downsample)),
        round(volume.shape[2] / int(xy_downsample)),
    )
    return resize(volume, shape, order=1, anti_aliasing=True).astype(np.float32)


def build_cell_signal(image_czyx: np.ndarray, config: Mapping[str, Any]) -> np.ndarray:
    """Build the Cellpose cell-signal input from one channel or a weighted composite.

    Args:
        image_czyx: Source image in CZYX order.
        config: Dataset settings dictionary.

    Returns:
        Float32 cell-signal volume in ZYX order.
    """
    if str(config["CELL_SOURCE"]) == "single":
        return normalize_volume_01(image_czyx[int(config["PM_CH"])])

    weights = {
        f"weight_ch{int(ch)}": float(w)
        for ch, w in dict(config["COMPOSITE_WEIGHTS"]).items()
        if float(w) > 0
    }
    comp = membrane_composite(
        image_czyx,
        Invert_PM=bool(config.get("INVERT_PM", False)),
        PM_Channel=int(config["PM_CH"]),
        **weights,
    )
    return normalize_volume_01(comp)


def build_cellpose_cli(
    *,
    indir_or_file: Path,
    savedir: Path,
    model: str,
    chan: int,
    chan2: int,
    diameter: int,
    min_size: int,
    anisotropy: float,
    flow_threshold: float,
    cellprob_threshold: float,
    batch_size: int,
    use_gpu: bool,
    channel_axis: int | None = None,
    is_file: bool = False,
) -> str:
    """Create a ready-to-run Cellpose CLI command string.

    Args:
        indir_or_file: Input directory or single TIFF path.
        savedir: Output directory for masks.
        model: Model name (`cyto3`/`nuclei`) or model path.
        chan: Main input channel index for Cellpose CLI.
        chan2: Secondary input channel index for Cellpose CLI.
        diameter: Object diameter parameter.
        min_size: Minimum object size filter.
        anisotropy: Z/XY voxel-size ratio.
        flow_threshold: Cellpose flow threshold.
        cellprob_threshold: Cellpose cellprob threshold.
        batch_size: Cellpose batch size.
        use_gpu: Whether to add `--use_gpu`.
        channel_axis: Optional channel axis for multi-channel TIFF input.
        is_file: If `True`, use `--image_path`; else use `--dir`.

    Returns:
        Command string suitable for a `.bat` script line.
    """
    model_str = str(model)
    base = [
        "python -m cellpose",
        f'--image_path "{indir_or_file}"' if is_file else f'--dir "{indir_or_file}"',
        f'--pretrained_model "{model_str}"' if model_str.startswith("/") else f"--pretrained_model {model_str}",
        "--do_3D",
        "--z_axis 0",
    ]

    if channel_axis is not None:
        base.append(f"--channel_axis {int(channel_axis)}")

    base.extend(
        [
            f"--chan {int(chan)}",
            f"--chan2 {int(chan2)}",
            f"--diameter {int(diameter)}",
            f"--anisotropy {float(anisotropy):.4f}",
            f"--flow_threshold {float(flow_threshold)}",
            f"--cellprob_threshold {float(cellprob_threshold)}",
            f"--min_size {int(min_size)}",
            f"--batch_size {int(batch_size)}",
        ]
    )

    if use_gpu:
        base.append("--use_gpu")

    base.extend(["--save_tif", "--no_npy", f'--savedir "{savedir}"', "--verbose"])
    return " ".join(base)


def build_nonzero_weight_dict(max_channels: int, weight_values: Sequence[Any]) -> dict[int, float]:
    """Convert indexed weight values to a compact dictionary of positive weights.

    Args:
        max_channels: Maximum number of expected channels.
        weight_values: Sequence where index = channel number.

    Returns:
        Dictionary `{channel_index: weight}` including only weights > 0.
    """
    out: dict[int, float] = {}
    for idx in range(int(max_channels)):
        if idx >= len(weight_values):
            break
        try:
            w = float(weight_values[idx])
        except Exception:
            w = 0.0
        if w > 0:
            out[idx] = w
    return out


def optional_path_str(value: Any) -> str | None:
    """Convert a model-path-like value to `str | None`.

    Args:
        value: Any incoming path-like object.

    Returns:
        Trimmed string path or `None` for empty values.
    """
    if value is None:
        return None
    text = str(value).strip()
    return text if text else None


def build_dataset_config_payload(
    *,
    dataset_folder: Path,
    selected_by: str,
    notes: str,
    generation_reason: str,
    settings: Mapping[str, Any],
) -> dict[str, Any]:
    """Build a full dataset config payload with metadata and nested settings.

    Args:
        dataset_folder: Dataset path used to derive `dataset_id`.
        selected_by: User label.
        notes: Free-text notes.
        generation_reason: Short reason value.
        settings: Normalized settings dictionary.

    Returns:
        JSON-ready payload dictionary.
    """
    folder = Path(dataset_folder)
    return {
        "dataset_id": folder.name,
        "config_version": "1.0",
        "selected_on": datetime.now().isoformat(timespec="seconds"),
        "selected_by": selected_by,
        "notes": notes,
        "generated_by": "NEW_infer_masks_cellpose.ipynb",
        "generation_reason": generation_reason,
        "settings": dict(settings),
    }
