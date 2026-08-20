#!/usr/bin/env bash
# Step 08 - gene regulatory network inference with pySCENIC (Aibar et al. 2017;
# Van de Sande et al. 2020).
#
# Input: $SCENIC_INPUT, a loom file exported by notebook 09 on Windows
# (raw counts of epithelial + mesenchymal cells, genes detected in >= 1 % cells).
#
#   1. GRN    : GRNBoost2 regresses every gene on all TFs -> TF-target co-expression modules
#   2. cisTarget (ctx): keep only targets whose promoters are enriched for the TF's motif
#                       (removes indirect co-expression) -> regulons
#   3. AUCell : per-cell regulon activity = area under the recovery curve of the
#               regulon's genes in the cell's expression ranking
#
# GRNBoost2 is stochastic, so it is run $SCENIC_RUNS times with different seeds;
# regulons recovered in >= 2/3 of runs are reported as robust (aggregate_regulons.py).
source "$(dirname "$0")/common.sh"
require_file "$SCENIC_INPUT" "copy scenic_input.loom from the Windows project (notebook 09) into $PKG_DIR/input/"
db="$SCENIC_DIR/resources"
mkdir -p "$db" "$SCENIC_DIR/runs"

for f in "$SCENIC_DB1" "$SCENIC_DB2" "$SCENIC_MOTIFS" "$SCENIC_TFS"; do
    dest="$db/$(basename "$f")"
    if [[ ! -s "$dest" ]]; then
        info "downloading $(basename "$f")"
        curl -fL --retry 5 -C - -o "$dest" "$SCENIC_BASE/$f"
    fi
done
( cd "$db" && sha256sum ./* > "$SCENIC_DIR/resources.sha256" )

loom="$SCENIC_INPUT"
for run in $(seq 1 "$SCENIC_RUNS"); do
    rd="$SCENIC_DIR/runs/run$run"; mkdir -p "$rd"
    seed=$(( SCENIC_SEED + run ))
    if [[ ! -s "$rd/adjacencies.tsv" ]]; then
        info "[run $run/$SCENIC_RUNS] GRNBoost2 (seed $seed)"
        in_env "$ENV_SCENIC" arboreto_with_multiprocessing.py "$loom" "$db/$(basename "$SCENIC_TFS")" \
            --method grnboost2 --output "$rd/adjacencies.tsv" --num_workers "$THREADS" --seed "$seed" \
            2>&1 | tee "$LOG_DIR/scenic_grn_run$run.log"
    fi
    if [[ ! -s "$rd/regulons.csv" ]]; then
        info "[run $run/$SCENIC_RUNS] cisTarget motif pruning"
        in_env "$ENV_SCENIC" pyscenic ctx "$rd/adjacencies.tsv" \
            "$db/$(basename "$SCENIC_DB1")" "$db/$(basename "$SCENIC_DB2")" \
            --annotations_fname "$db/$(basename "$SCENIC_MOTIFS")" \
            --expression_mtx_fname "$loom" \
            --output "$rd/regulons.csv" --mask_dropouts --num_workers "$THREADS" \
            2>&1 | tee "$LOG_DIR/scenic_ctx_run$run.log"
    fi
    if [[ ! -s "$rd/auc.loom" ]]; then
        info "[run $run/$SCENIC_RUNS] AUCell"
        in_env "$ENV_SCENIC" pyscenic aucell "$loom" "$rd/regulons.csv" \
            --output "$rd/auc.loom" --num_workers "$THREADS" --seed "$seed" \
            2>&1 | tee "$LOG_DIR/scenic_aucell_run$run.log"
    fi
done

info "aggregating regulons across runs"
in_env "$ENV_SCENIC" python "$PKG_DIR/python/aggregate_regulons.py" "$SCENIC_DIR/runs" "$SCENIC_DIR/consensus" \
    --min-fraction 0.66
info "pySCENIC finished: $SCENIC_DIR/consensus"
