#!/usr/bin/env bash
# Recreates the "waveorder-cytoland" conda environment used by NEW_waveorder_cytoland.ipynb
# (waveorder QPI reconstruction + Cytoland VSCyto3D virtual staining).
# Run from Git Bash / WSL / Linux / macOS with conda on PATH:  bash create_waveorder_cytoland_env.sh
# Optional: CUDA=cu124 bash create_waveorder_cytoland_env.sh   (or CUDA=cpu)
set -euo pipefail

ENV_NAME="waveorder-cytoland"
CUDA="${CUDA:-cu128}"   # PyTorch wheel flavor: cu128, cu124, cpu
CKPT_URL="https://public.czbiohub.org/comp.micro/viscy/VS_models/VSCyto3D/epoch=83-step=14532-loss=0.492.ckpt"
CKPT_PATH="$HOME/.cache/cytoland/VSCyto3D.ckpt"

# waveorder 3.x and cytoland require Python >= 3.12
conda create -y -n "$ENV_NAME" python=3.12
PY=(conda run --no-capture-output -n "$ENV_NAME" python)

"${PY[@]}" -m pip install --upgrade pip
"${PY[@]}" -m pip install torch torchvision --index-url "https://download.pytorch.org/whl/$CUDA"
"${PY[@]}" -m pip install waveorder cytoland tifffile scikit-image scipy seaborn napari[all] infer-subc

# Installing cytoland can swap in a CPU-only torch; force the chosen build back (without touching deps)
"${PY[@]}" -m pip install --force-reinstall --no-deps torch torchvision --index-url "https://download.pytorch.org/whl/$CUDA"

# VSCyto3D checkpoint (~420 MB), used by the notebook
if [ ! -f "$CKPT_PATH" ]; then
    mkdir -p "$(dirname "$CKPT_PATH")"
    curl -L -o "$CKPT_PATH" "$CKPT_URL"
fi

"${PY[@]}" -c "import torch, waveorder, cytoland; print('torch', torch.__version__, '| CUDA available:', torch.cuda.is_available())"
echo "Done. Set worker_python in the notebook to this environment's python executable."
