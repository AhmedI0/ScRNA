#!/usr/bin/env bash
# Reproduce the full analysis on Linux/macOS: environment -> data -> notebooks -> figures.
#   bash run_pipeline.sh              # everything
#   SKIP_SETUP=1 bash run_pipeline.sh # env already built
#   FROM=05 bash run_pipeline.sh      # resume at notebook 05
# R must be on PATH (Rscript); set r.rscript_linux in config/config.yaml otherwise.
set -euo pipefail
cd "$(dirname "$0")"
ENV_NAME="${ENV_NAME:-scrna-ipf}"
CUDA="${CUDA:-cu126}"   # cu126 | cu124 | cpu

if [[ "${SKIP_SETUP:-0}" != "1" ]]; then
    conda env create -n "$ENV_NAME" -f envs/environment.yml || conda env update -n "$ENV_NAME" -f envs/environment.yml
    conda run -n "$ENV_NAME" python -m pip install torch --index-url "https://download.pytorch.org/whl/$CUDA"
    conda run -n "$ENV_NAME" python -m pip install -r envs/requirements.txt
    conda run -n "$ENV_NAME" python -m pip install -e . --no-deps
    conda run -n "$ENV_NAME" python -m ipykernel install --sys-prefix --name scrna-ipf --display-name "Python (scrna-ipf)"
    Rscript envs/r_packages.R
fi
run() { conda run --no-capture-output -n "$ENV_NAME" "$@"; }
run python -m pytest -q tests
run python workflow/scripts/download_data.py
if [[ -n "${FROM:-}" ]]; then run python workflow/scripts/run_notebooks.py --from "$FROM"
else run python workflow/scripts/run_notebooks.py; fi
echo "Done. Figures: results/figures  Tables: results/tables"
