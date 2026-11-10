"""Worker for NEW_waveorder_cytoland.ipynb.

Runs in a dedicated Python >= 3.12 environment (waveorder + cytoland + torch):
  1. waveorder 3D quantitative phase imaging (QPI) reconstruction of a through-focus
     transmitted-light stack
  2. Cytoland VSCyto3D virtual staining of nuclei and plasma membrane from the phase

Usage: python waveorder_cytoland_worker.py <params.json>
All inputs/outputs are exchanged as ZYX TIFF files described in the params file.
"""
import json
import sys
from pathlib import Path

import numpy as np
import tifffile
import torch
import xarray as xr


def reconstruct_phase(zyx, p):
    from waveorder.api import phase

    z_um, y_um, x_um = p["scale_zyx_um"]
    if not np.isclose(y_um, x_um):
        raise ValueError(
            "This WaveOrder version accepts only scalar yx_pixel_size; "
            f"received anisotropic Y/X spacing ({y_um}, {x_um}) um."
        )
    settings = phase.Settings(
        transfer_function=phase.TransferFunctionSettings(
            yx_pixel_size=float(y_um),
            z_pixel_size=float(z_um),
            wavelength_illumination=p["wavelength_um"],
            numerical_aperture_detection=p["na_detection"],
            numerical_aperture_illumination=p["na_illumination"],
            index_of_refraction_media=p["ri_media"],
            z_padding=p["z_padding"],
            invert_phase_contrast=p["invert_phase_contrast"],
        ),
        apply_inverse=phase.ApplyInverseSettings(
            regularization_strength=p["regularization_strength"],
        ),
    )
    czyx = xr.DataArray(
        zyx[None].astype(np.float32),
        dims=("c", "z", "y", "x"),
        coords={
            "c": ["Brightfield"],
            "z": np.arange(zyx.shape[0]) * z_um,
            "y": np.arange(zyx.shape[1]) * y_um,
            "x": np.arange(zyx.shape[2]) * x_um,
        },
    )
    out = phase.reconstruct(czyx, recon_dim=3, settings=settings, device=p["device"])
    return np.asarray(out.values[0], dtype=np.float32)


def load_vscyto3d(ckpt, device):
    from cytoland.engine import AugmentedPredictionVSUNet, VSUNet

    model_config = {
        "in_channels": 1,
        "out_channels": 2,
        "encoder_blocks": [3, 3, 9, 3],
        "dims": [96, 192, 384, 768],
        "decoder_conv_blocks": 2,
        "stem_kernel_size": [5, 4, 4],
        "in_stack_depth": 15,
        "pretraining": False,
    }
    model = VSUNet.load_from_checkpoint(
        ckpt,
        architecture="fcmae",
        model_config=model_config,
        map_location="cpu",
    )
    model.eval().to(device)
    predictor = AugmentedPredictionVSUNet(model=model.model).to(device).eval()
    return predictor, model_config["in_stack_depth"]


@torch.inference_mode()
def virtual_stain(phase_zyx, p):
    device = torch.device(p["device"] if p["device"] != "auto" else ("cuda" if torch.cuda.is_available() else "cpu"))
    predictor, depth = load_vscyto3d(p["checkpoint"], device)
    step = int(p.get("z_stride", 1))
    if not 1 <= step <= depth:
        raise ValueError(f"z_stride must be between 1 and the model depth ({depth}), got {step}")
    if p.get("tile_yx", 0):
        raise ValueError("Cytoland's package predictor requires full XY planes; set tile_yx=0.")

    # Match the authors' field-of-view median/IQR normalization.
    median = np.median(phase_zyx)
    q25, q75 = np.percentile(phase_zyx, [25, 75])
    iqr = q75 - q25
    if iqr <= 1e-8:
        raise ValueError("Phase volume has near-zero IQR; cannot normalize for VSCyto3D inference.")
    vol = (phase_zyx - median) / iqr
    original_z, _, _ = vol.shape
    padded_z = max(depth, depth + ((max(original_z - depth, 0) + step - 1) // step) * step)
    if padded_z > original_z:
        vol = np.pad(vol, ((0, padded_z - original_z), (0, 0), (0, 0)), mode="reflect")
    volume = torch.from_numpy(vol[None, None].astype(np.float32)).to(device)
    prediction = predictor.predict_sliding_windows(volume, out_channel=2, step=step)
    return prediction[0, :, :original_z].cpu().numpy().astype(np.float32)


def resample(vol, src_zyx, dst_zyx, shape=None, order=1):
    """Resample a ZYX volume from voxel size src to dst (micrometers) with linear interpolation."""
    from scipy.ndimage import gaussian_filter, zoom

    if shape is None:
        shape = [max(int(round(n * s / d)), 1) for n, s, d in zip(vol.shape, src_zyx, dst_zyx)]
    factors = [t / n for t, n in zip(shape, vol.shape)]
    antialias_sigma = [max((1 / factor - 1) / 2, 0) for factor in factors]
    if any(sigma > 0 for sigma in antialias_sigma):
        vol = gaussian_filter(vol, sigma=antialias_sigma, mode="nearest")
    return zoom(vol, factors, order=order, mode="nearest", grid_mode=True).astype(np.float32)


def main(params_path):
    p = json.loads(Path(params_path).read_text())
    out_dir = Path(p["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    zyx = np.asarray(tifffile.imread(p["input_tiff"]), dtype=np.float32)
    print("input", zyx.shape, flush=True)

    print("waveorder phase reconstruction (native voxel size)...", flush=True)
    phase_zyx = reconstruct_phase(zyx, p)
    tifffile.imwrite(out_dir / "phase.tiff", phase_zyx)

    # VSCyto3D was trained at ~0.2 x 0.1 x 0.1 um (ZYX); resample the phase to match
    src = p["scale_zyx_um"]
    dst = p.get("model_voxel_zyx_um")
    model_in = phase_zyx
    if dst and not np.allclose(src, dst, rtol=0.05):
        model_in = resample(phase_zyx, src, dst)
        print(f"resampled phase {phase_zyx.shape} -> {model_in.shape} for the model ({src} -> {dst} um)", flush=True)

    print("Cytoland VSCyto3D virtual staining...", flush=True)
    try:
        vs = virtual_stain(model_in, p)
    except torch.cuda.OutOfMemoryError:
        if p["device"] == "cpu" or not torch.cuda.is_available():
            raise
        torch.cuda.empty_cache()
        print("GPU out of memory; retrying VSCyto3D inference on CPU...", flush=True)
        vs = virtual_stain(model_in, dict(p, device="cpu"))
    if vs.shape[1:] != phase_zyx.shape:
        vs = np.stack([resample(c, None, None, shape=phase_zyx.shape) for c in vs])
    tifffile.imwrite(out_dir / "vs_nucleus.tiff", vs[0].astype(np.float32))
    tifffile.imwrite(out_dir / "vs_membrane.tiff", vs[1].astype(np.float32))
    print("done", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
