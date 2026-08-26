#!/usr/bin/env bash
# Step 10 - collect everything the Windows analysis needs into ONE archive.
#   results_to_return/
#     starsolo/   per-library Summary.csv, Log.final.out, Features.stats, filtered matrices
#                 (+ the raw Gene matrix when CellBender was skipped here)
#     h5ad/       per-library AnnData with spliced/unspliced layers
#     cellbender/ background-removed matrices, reports and metrics
#     scenic/     consensus regulons, AUCell matrix, per-run regulon tables
#     qc/         MultiQC report (FastQC + STAR)
#     logs/       all logs and exported environment locks
#     MANIFEST.sha256
source "$(dirname "$0")/common.sh"
R="$RESULTS_DIR"
rm -rf "$R"; mkdir -p "$R"/{starsolo,h5ad,cellbender,scenic,qc,logs}

for lib in $(library_ids); do
    s="$SOLO_DIR/$lib"
    [[ -d "$s/Solo.out" ]] || continue
    mkdir -p "$R/starsolo/$lib"
    cp "$s/Log.final.out" "$R/starsolo/$lib/" 2>/dev/null || true
    cp "$s/Solo.out/Barcodes.stats" "$R/starsolo/$lib/" 2>/dev/null || true
    for feat in Gene GeneFull; do
        mkdir -p "$R/starsolo/$lib/$feat"
        cp "$s/Solo.out/$feat/Summary.csv" "$s/Solo.out/$feat/Features.stats" "$R/starsolo/$lib/$feat/" 2>/dev/null || true
        cp -r "$s/Solo.out/$feat/filtered" "$R/starsolo/$lib/$feat/" 2>/dev/null || true
    done
    # no NVIDIA GPU here: CellBender runs on the analysis PC (workflow/scripts/run_cellbender.py),
    # which needs the raw droplet matrix (all barcodes, empty droplets included)
    if [[ ! -s "$CELLBENDER_DIR/$lib/cellbender_filtered.h5" ]]; then
        cp -r "$s/Solo.out/Gene/raw" "$R/starsolo/$lib/Gene/"
    fi
    [[ -s "$SOLO_DIR/h5ad/$lib.h5ad" ]] && cp "$SOLO_DIR/h5ad/$lib.h5ad" "$R/h5ad/"
    if [[ -d "$CELLBENDER_DIR/$lib" ]]; then
        mkdir -p "$R/cellbender/$lib"
        cp "$CELLBENDER_DIR/$lib"/cellbender_filtered.h5 "$CELLBENDER_DIR/$lib"/*.html \
           "$CELLBENDER_DIR/$lib"/*metrics.csv "$CELLBENDER_DIR/$lib"/*cell_barcodes.csv \
           "$CELLBENDER_DIR/$lib"/*.pdf "$R/cellbender/$lib/" 2>/dev/null || true
    fi
done
[[ -s "$SOLO_DIR/starsolo_summary.tsv" ]] && cp "$SOLO_DIR/starsolo_summary.tsv" "$R/starsolo/"

if [[ -d "$SCENIC_DIR/consensus" ]]; then
    cp -r "$SCENIC_DIR/consensus" "$R/scenic/"
    for run in "$SCENIC_DIR"/runs/run*; do
        mkdir -p "$R/scenic/$(basename "$run")"
        cp "$run/regulons.csv" "$R/scenic/$(basename "$run")/" 2>/dev/null || true
    done
    cp "$SCENIC_DIR/resources.sha256" "$R/scenic/" 2>/dev/null || true
fi

if conda_cmd >/dev/null; then
    in_env "$ENV_UPSTREAM" multiqc -f -q -o "$R/qc" "$QC_DIR/fastqc" "$SOLO_DIR" || warn "MultiQC failed"
fi
cp -r "$LOG_DIR"/. "$R/logs/"
cp "$SAMPLES_TSV" "$R/"
{ echo "host: $(hostname)"; echo "date: $(date -Iseconds)"; echo "threads: $THREADS"; echo "mem_gb: $MEM_GB"; \
  echo "test_mode: $TEST_MODE"; } > "$R/RUN_INFO.txt"

( cd "$R" && find . -type f ! -name MANIFEST.sha256 -print0 | sort -z | xargs -0 sha256sum > MANIFEST.sha256 )
archive="$(dirname "$R")/ipf_upstream_results_$(date +%Y%m%d).tar"   # next to results_to_return/
tar -cf "$archive" -C "$(dirname "$R")" "$(basename "$R")"
sha256sum "$archive" > "$archive.sha256"
info "archive ready: $archive ($(du -h "$archive" | cut -f1))"
info "copy it to the Windows project folder data/linux_results/ and extract it there (see README, section 6)"
