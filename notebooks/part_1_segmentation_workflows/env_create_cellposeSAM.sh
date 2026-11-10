#!/usr/bin/env bash
# Recreates the "infer-subc-cellposeSAM" environment used by
# notebooks/part_1_segmentation_workflows/NEW_infer_masks_cellposeSAM.ipynb
# (Cellpose 4 / Cellpose-SAM + infer-subc + napari).
#
# Run from the infer-subc repository root (the folder containing pyproject.toml), e.g. in Git Bash
# or Anaconda Prompt with bash:   bash env_create_cellposeSAM.sh
#
# GPU: installs PyTorch built for CUDA 12.6 (needs an NVIDIA driver supporting CUDA 12.6 or newer).
# For a different CUDA version change TORCH_INDEX (e.g. .../whl/cu128), or use .../whl/cpu for CPU only.

set -e

ENV_NAME="infer-subc-cpdino"
TORCH_INDEX="https://download.pytorch.org/whl/cu126"

# Python 3.11+ is required by dinov3 (the cpdino / cpdino-vitb models)
conda create -y -n "$ENV_NAME" python=3.11

# `conda run` executes inside the new environment without needing `conda activate`.
RUN="conda run --no-capture-output -n $ENV_NAME"

$RUN python -m pip install --upgrade pip

# PyTorch first, from the CUDA wheel index (installs the GPU build instead of the CPU build from PyPI)
$RUN python -m pip install torch torchvision --index-url "$TORCH_INDEX"

# Cellpose 4 (Cellpose-SAM) and its segment_anything dependency
$RUN python -m pip install "cellpose>=4" segment_anything

# DINOv3 backbone for the cpdino / cpdino-vitb models (needs git)
$RUN python -m pip install git+https://github.com/facebookresearch/dinov3

# infer-subc from this repository (editable) -- also installs napari 0.5.5, numpy 1.26.4, aicsimageio, etc.
$RUN python -m pip install -e .

# Napari GUI backend and Jupyter kernel support
$RUN python -m pip install pyqt5 ipykernel
$RUN python -m ipykernel install --user --name "$ENV_NAME" --display-name "$ENV_NAME"

# Check the install
$RUN python -c "import torch, cellpose, napari, numpy; from importlib.metadata import version; print('cellpose', version('cellpose'), '| torch', torch.__version__, '| CUDA available:', torch.cuda.is_available(), '| napari', napari.__version__, '| numpy', numpy.__version__)"

echo
echo "Done. In VS Code / Jupyter choose the kernel \"$ENV_NAME\"."
echo "The first Cellpose-SAM run downloads ~1.2 GB of model weights (internet required) to ~/.cellpose/models."
