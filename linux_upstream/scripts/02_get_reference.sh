#!/usr/bin/env bash
# Step 02 - download the reference genome/annotation and the 10x barcode whitelist.
#
# We use exactly the reference of the original study (Cell Ranger 3.0.2,
# refdata-cellranger-GRCh38-3.0.0), so differences between our STARsolo counts
# and the published Cell Ranger counts reflect the aligner, not the annotation.
# The 737K-august-2016 list contains every valid barcode of 10x 3' v2 and 5' v1 gel beads.
source "$(dirname "$0")/common.sh"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

tarball="$REF_DIR/$REF_NAME.tar.gz"
if [[ ! -d "$REF_DIR/$REF_NAME/fasta" ]]; then
    info "downloading $REF_URL (11.4 GB, resumable)"
    curl -fL --retry 10 --retry-delay 30 -C - -o "$tarball" "$REF_URL"
    size=$(stat -c %s "$tarball")
    [[ "$size" == "$REF_BYTES" ]] || die "reference tarball has $size bytes, expected $REF_BYTES"
    info "extracting reference"
    tar -xzf "$tarball" -C "$REF_DIR"
    sha256sum "$tarball" > "$REF_DIR/$REF_NAME.tar.gz.sha256"
    rm -f "$tarball"
fi
require_file "$REF_DIR/$REF_NAME/fasta/genome.fa" "reference FASTA"
require_file "$REF_DIR/$REF_NAME/genes/genes.gtf" "reference GTF"
info "reference: $(grep -c '^>' "$REF_DIR/$REF_NAME/fasta/genome.fa") sequences, \
$(awk '$3=="gene"' "$REF_DIR/$REF_NAME/genes/genes.gtf" | wc -l) genes"

if [[ ! -s "$WHITELIST" ]]; then
    info "downloading 10x barcode whitelist"
    curl -fsSL "$WHITELIST_URL" | gzip -dc > "$WHITELIST"
fi
md5=$(md5sum "$WHITELIST" | cut -d' ' -f1)
[[ "$md5" == "$WHITELIST_MD5" ]] || die "whitelist md5 $md5 != $WHITELIST_MD5"
info "whitelist OK: $(wc -l < "$WHITELIST") barcodes"
