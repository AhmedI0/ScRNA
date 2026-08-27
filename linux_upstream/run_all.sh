#!/usr/bin/env bash
# ============================================================================
#  IPF scRNA-seq - Linux upstream pipeline (STARsolo, CellBender, pySCENIC)
#
#  usage:  bash run_all.sh                 # all steps 00-10
#          bash run_all.sh 04 05 06        # selected steps
#          TEST_MODE=1 LIBRARIES="F01214" bash run_all.sh 00 02 03 04 05 06 09 10
#
#  Long runs: start inside tmux/screen or with nohup, e.g.
#          nohup bash run_all.sh > run_all.log 2>&1 &
# ============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/scripts/common.sh"

declare -A STEP=(
  [00]=00_check_system.sh   [01]=01_install_envs.sh   [02]=02_get_reference.sh
  [03]=03_build_star_index.sh [04]=04_download_fastq.sh [05]=05_read_qc.sh
  [06]=06_starsolo.sh       [07]=07_cellbender.sh     [08]=08_pyscenic.sh
  [09]=09_export_h5ad.sh    [10]=10_collect_results.sh
)
steps=("$@")
(( ${#steps[@]} )) || steps=(00 01 02 03 04 05 06 07 08 09 10)

info "project dir: $PROJECT_DIR | threads: $THREADS | RAM: ${MEM_GB} GB | test mode: $TEST_MODE"
for s in "${steps[@]}"; do
    script="${STEP[$s]:-}"
    [[ -n "$script" ]] || die "unknown step $s"
    if [[ "$s" == "08" && ! -s "$SCENIC_INPUT" ]]; then
        warn "step 08 skipped: $SCENIC_INPUT not found"; continue
    fi
    info "================ step $s: $script ================"
    t0=$(date +%s)
    bash "$HERE/scripts/$script" 2>&1 | tee "$LOG_DIR/step_${s}.log"
    info "step $s finished in $(( ($(date +%s) - t0) / 60 )) min"
done
info "all requested steps finished"
