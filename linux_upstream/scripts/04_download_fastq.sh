#!/usr/bin/env bash
# Step 04 - download raw FASTQ files from ENA and verify their MD5 checksums.
#
# ENA mirrors SRA and serves gzipped FASTQ over HTTPS together with MD5 sums,
# so no SRA toolkit conversion (fasterq-dump) is needed. Downloads resume after
# interruption and run $DOWNLOAD_JOBS in parallel.
# TEST_MODE=1 streams only the first $TEST_READS read pairs of one run per library.
source "$(dirname "$0")/common.sh"
mkdir -p "$FASTQ_DIR"

fetch() {  # fetch <url> <dest> <md5>
    local url="$1" dest="$2" md5="$3"
    if [[ -s "$dest" ]] && [[ "$(md5sum "$dest" | cut -d' ' -f1)" == "$md5" ]]; then
        echo "[skip] $(basename "$dest")"; return 0
    fi
    for attempt in 1 2 3 4 5; do
        if curl -fsSL --retry 5 --retry-delay 20 -C - -o "$dest" "$url"; then
            if [[ "$(md5sum "$dest" | cut -d' ' -f1)" == "$md5" ]]; then echo "[ok]   $(basename "$dest")"; return 0; fi
            echo "[md5 mismatch] $(basename "$dest") attempt $attempt"; rm -f "$dest"
        fi
        sleep 30
    done
    echo "[FAILED] $url" >&2; return 1
}
export -f fetch

fetch_head() {  # first N reads only (test mode); no md5 possible on a partial file
    local url="$1" dest="$2" n="$3"
    [[ -s "$dest" ]] && { echo "[skip] $(basename "$dest")"; return 0; }
    (curl -fsSL "$url" | gzip -dc | head -n $(( n * 4 )) | gzip -1 > "$dest") || true
    echo "[test] $(basename "$dest"): $(( $(gzip -dc "$dest" | wc -l) / 4 )) reads"
}
export -f fetch_head

jobs_file="$FASTQ_DIR/download_jobs.tsv"
: > "$jobs_file"
for lib in $(library_ids); do
    mkdir -p "$FASTQ_DIR/$lib"
    paste <(sheet_col "$lib" run_accession) <(sheet_col "$lib" fastq_1) <(sheet_col "$lib" md5_1) \
          <(sheet_col "$lib" fastq_2) <(sheet_col "$lib" md5_2) | \
    { if [[ "$TEST_MODE" == "1" ]]; then head -1; else cat; fi; } | \
    while IFS=$'\t' read -r run u1 m1 u2 m2; do
        printf "%s\t%s\t%s\n" "$u1" "$FASTQ_DIR/$lib/${run}_1.fastq.gz" "$m1" >> "$jobs_file"
        printf "%s\t%s\t%s\n" "$u2" "$FASTQ_DIR/$lib/${run}_2.fastq.gz" "$m2" >> "$jobs_file"
    done
done
info "$(wc -l < "$jobs_file") FASTQ files to fetch for: $(library_ids | paste -sd ' ')"

if [[ "$TEST_MODE" == "1" ]]; then
    warn "TEST_MODE: only the first $TEST_READS reads of one run per library"
    while IFS=$'\t' read -r u d _md5; do fetch_head "$u" "$d" "$TEST_READS"; done < "$jobs_file"
else
    tr '\t' '\n' < "$jobs_file" | xargs -P "$DOWNLOAD_JOBS" -n 3 -d $'\n' bash -c 'fetch "$0" "$1" "$2"'
fi
info "download finished: $(du -sh "$FASTQ_DIR" | cut -f1) in $FASTQ_DIR"
