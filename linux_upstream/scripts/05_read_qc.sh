#!/usr/bin/env bash
# Step 05 - read-level QC and read-layout detection.
#   * detect_read_layout.py : which FASTQ of each pair holds the cell barcode/UMI
#   * FastQC on a 1M-read subsample of every file, summarised with MultiQC
source "$(dirname "$0")/common.sh"
mkdir -p "$QC_DIR/fastqc"

for lib in $(library_ids); do
    info "read layout: $lib"
    in_env "$ENV_UPSTREAM" python "$PKG_DIR/python/detect_read_layout.py" "$FASTQ_DIR/$lib" "$WHITELIST" \
        --cb-len "$CB_LEN" | tee "$LOG_DIR/read_layout_$lib.log"
done

info "FastQC on 1M-read subsamples"
tmp="$QC_DIR/subsample"; mkdir -p "$tmp"
for f in "$FASTQ_DIR"/*/*.fastq.gz; do
    out="$tmp/$(basename "$f")"
    [[ -s "$out" ]] || in_env "$ENV_UPSTREAM" seqkit head -n 1000000 "$f" -o "$out"
done
in_env "$ENV_UPSTREAM" fastqc -q -t "$THREADS" --outdir "$QC_DIR/fastqc" "$tmp"/*.fastq.gz
rm -rf "$tmp"
info "read QC finished (MultiQC report is built in step 10)"
