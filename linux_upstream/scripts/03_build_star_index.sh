#!/usr/bin/env bash
# Step 03 - build the STAR genome index (one-off, ~1 h with 16 threads).
#
# --sjdbOverhang 100 (read length - 1 for typical 101-151 bp cDNA reads) is the
# STAR default and works for all read lengths in this data set.
# With < 40 GB RAM a sparse suffix array (--genomeSAsparseD 2) halves the
# memory needed at the cost of somewhat slower mapping.
source "$(dirname "$0")/common.sh"

if [[ -s "$STAR_INDEX/SAindex" ]]; then
    info "STAR index already present: $STAR_INDEX"; exit 0
fi
mkdir -p "$STAR_INDEX"
extra=()
if (( MEM_GB < 40 )); then
    warn "only ${MEM_GB} GB RAM: building a sparse index"
    extra=(--genomeSAsparseD 2 --limitGenomeGenerateRAM $(( (MEM_GB - 4) * 1024 * 1024 * 1024 )))
fi
info "building STAR index with $THREADS threads"
in_env "$ENV_UPSTREAM" STAR \
    --runMode genomeGenerate \
    --runThreadN "$THREADS" \
    --genomeDir "$STAR_INDEX" \
    --genomeFastaFiles "$REF_DIR/$REF_NAME/fasta/genome.fa" \
    --sjdbGTFfile "$REF_DIR/$REF_NAME/genes/genes.gtf" \
    --sjdbOverhang 100 \
    --outFileNamePrefix "$STAR_INDEX/" \
    "${extra[@]}"
info "STAR index ready: $(du -sh "$STAR_INDEX" | cut -f1)"
