#!/usr/bin/env bash
# Create the environment for NEW_infer_masks_cellpose.ipynb,
# NEW_cellpose-prep.ipynb, NEW_cellpose-training.ipynb, and the Cellpose CLI.
# Run: bash create_cellpose_env.sh
# Optional: CUDA=cpu bash create_cellpose_env.sh
# For another supported PyTorch wheel, for example: CUDA=cu126 bash create_cellpose_env.sh
set -euo pipefail

ENV_NAME="cellpose"
CUDA="${CUDA:-cu128}"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "$SCRIPT_DIR/../.." && pwd)"

# conda create -y -n "$ENV_NAME" python=3.12 pip
PY=(conda run --no-capture-output -n "$ENV_NAME" python)

# Install a GPU-enabled PyTorch build before installing Cellpose.
"${PY[@]}" -m pip install "torch==2.11.0" --index-url "https://download.pytorch.org/whl/$CUDA"
"${PY[@]}" -m pip install -e "$PROJECT_DIR" "cellpose==3.1.1.3" ipykernel jupyterlab pyqt5
"${PY[@]}" -m pip check
"${PY[@]}" -m ipykernel install --user --name "$ENV_NAME" --display-name "$ENV_NAME"

"${PY[@]}" -c "import torch, napari, tifffile; from cellpose import models, train, metrics; from infer_subc.core.file_io import read_czi_image; from infer_subc.organelles.membrane import membrane_composite; from importlib.metadata import version; print('Cellpose:', version('cellpose'), '| PyTorch:', torch.__version__, '| CUDA available:', torch.cuda.is_available())"
echo
echo "Select the '$ENV_NAME' kernel in Jupyter or VS Code."
echo "For the CLI, run: conda activate $ENV_NAME"
echo "The first use of pretrained models may download weights; Napari requires a graphical display."
