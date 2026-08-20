#!/usr/bin/env bash
# Step 07 - ambient RNA removal with CellBender remove-background (Fleming et al. 2023).
#
# Every droplet contains some "soup": cell-free mRNA released by lysed cells.
# In lung tissue the soup is dominated by highly expressed epithelial and
# macrophage transcripts (SFTPC, SFTPB, SCGB1A1, APOC1...), which then appear at
# low levels in unrelated cells and can masquerade as "intermediate" states.
# CellBender fits a deep generative model to ALL droplets (including empty
# ones, so it needs the RAW matrix) and returns background-corrected counts.
#
#   --expected-cells           STARsolo EmptyDrops cell count
#   --total-droplets-included  cells + a large set of surely-empty droplets
#   --fpr 0.01                 tolerated false-positive rate of removed signal
source "$(dirname "$0")/common.sh"

if ! has_gpu; then
    if [[ "$USE_GPU" == "no" ]]; then
        warn "USE_GPU=no: skipping CellBender; step 10 returns the raw matrices for workflow/scripts/run_cellbender.py"
        exit 0
    fi
    warn "no GPU: CellBender will run on CPU and may take several hours per library"
    cuda=()
else
    cuda=(--cuda)
fi

for lib in $(library_ids); do
    raw="$SOLO_DIR/$lib/Solo.out/Gene/raw"
    out="$CELLBENDER_DIR/$lib"
    require_file "$raw/matrix.mtx.gz" "run 06_starsolo.sh"
    if [[ -s "$out/cellbender_filtered.h5" ]]; then info "[skip] $lib"; continue; fi
    mkdir -p "$out"
    ncell=$(grep -E '^Estimated Number of Cells' "$SOLO_DIR/$lib/Solo.out/Gene/Summary.csv" | cut -d, -f2)
    ndrop=$(( ncell * 3 > 25000 ? ncell * 3 : 25000 ))
    info "CellBender: $lib (expected cells $ncell, droplets included $ndrop)"
    ( cd "$out" && in_env "$ENV_CELLBENDER" cellbender remove-background \
        --input "$raw" \
        --output "$out/cellbender.h5" \
        --expected-cells "$ncell" \
        --total-droplets-included "$ndrop" \
        --fpr 0.01 \
        --epochs 150 \
        "${cuda[@]}" ) 2>&1 | tee "$LOG_DIR/cellbender_$lib.log"
done
info "CellBender finished"
