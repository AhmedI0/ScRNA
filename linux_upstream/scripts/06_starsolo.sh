#!/usr/bin/env bash
# Step 06 - align and count with STARsolo (Kaminow, Yunusov & Dobin 2021).
#
# STARsolo re-implements the Cell Ranger algorithm inside the STAR aligner and,
# with the "_CR" options below, reproduces Cell Ranger counts while being
# open-source and several times faster.
#
# Key settings
#   CB_UMI_Simple, CB 16 bp + UMI 10 bp   10x 5' v1 bead structure
#   --soloBarcodeReadLength 0             barcode read may be longer than 26 bp
#   --soloStrand Reverse                  5' kits: the cDNA read is antisense
#   --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts, --soloUMIfiltering MultiGeneUMI_CR,
#   --soloUMIdedup 1MM_CR                 Cell Ranger >= 4 barcode/UMI rules
#   --soloCellFilter EmptyDrops_CR        Cell Ranger 3+ cell calling (EmptyDrops)
#   --soloFeatures Gene GeneFull Velocyto
#       Gene      exonic reads (= Cell Ranger default)
#       GeneFull  exonic + intronic reads (pre-mRNA)
#       Velocyto  spliced / unspliced / ambiguous counts for RNA velocity
#   --outSAMtype None                     no BAM (saves ~30 GB per library)
source "$(dirname "$0")/common.sh"
require_file "$STAR_INDEX/SAindex" "run 03_build_star_index.sh"
require_file "$WHITELIST" "run 02_get_reference.sh"

for lib in $(library_ids); do
    out="$SOLO_DIR/$lib"
    if [[ -s "$out/Solo.out/Gene/Summary.csv" ]]; then info "[skip] $lib already counted"; continue; fi
    layout="$FASTQ_DIR/$lib/read_layout.tsv"
    require_file "$layout" "run 05_read_qc.sh"
    cdna=$(awk -F'\t' -v d="$FASTQ_DIR/$lib" 'NR>1 {printf "%s%s/%s", (n++ ? "," : ""), d, $3}' "$layout")
    bc=$(awk -F'\t' -v d="$FASTQ_DIR/$lib" 'NR>1 {printf "%s%s/%s", (n++ ? "," : ""), d, $2}' "$layout")
    mkdir -p "$out"
    info "STARsolo: $lib ($(awk 'NR>1' "$layout" | wc -l) runs)"
    in_env "$ENV_UPSTREAM" STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --genomeLoad NoSharedMemory \
        --readFilesIn "$cdna" "$bc" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$out/" \
        --outSAMtype None \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "$WHITELIST" \
        --soloCBstart 1 --soloCBlen "$CB_LEN" \
        --soloUMIstart $(( CB_LEN + 1 )) --soloUMIlen "$UMI_LEN" \
        --soloBarcodeReadLength 0 \
        --soloStrand "$SOLO_STRAND" \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR \
        --soloCellFilter EmptyDrops_CR \
        --soloFeatures Gene GeneFull Velocyto \
        --soloOutFileNames Solo.out/ features.tsv barcodes.tsv matrix.mtx \
        --outFilterScoreMin 30 \
        2>&1 | tee "$LOG_DIR/starsolo_$lib.log"

    # gzip matrices so that the folders look like Cell Ranger v3 output (read by CellBender/scanpy)
    for feat in Gene GeneFull; do
        for kind in raw filtered; do
            d="$out/Solo.out/$feat/$kind"
            [[ -d "$d" ]] && in_env "$ENV_UPSTREAM" pigz -f -p "$THREADS" "$d"/*.tsv "$d"/*.mtx
        done
    done
    in_env "$ENV_UPSTREAM" pigz -f -p "$THREADS" "$out"/Solo.out/Velocyto/raw/*.tsv "$out"/Solo.out/Velocyto/raw/*.mtx || true
    rm -rf "$out/_STARtmp"
    info "$lib: $(grep -E '^Estimated Number of Cells' "$out/Solo.out/Gene/Summary.csv" | cut -d, -f2) cells, \
$(grep -E '^Median UMI per Cell' "$out/Solo.out/Gene/Summary.csv" | cut -d, -f2) median UMI/cell"
done

# one table with the key mapping/counting statistics of every library
summary="$SOLO_DIR/starsolo_summary.tsv"
first=1
for lib in $(library_ids); do
    s="$SOLO_DIR/$lib/Solo.out/Gene/Summary.csv"; [[ -s "$s" ]] || continue
    if (( first )); then { printf "library_id\t"; cut -d, -f1 "$s" | paste -sd '\t'; } > "$summary"; first=0; fi
    { printf "%s\t" "$lib"; cut -d, -f2 "$s" | paste -sd '\t'; } >> "$summary"
done
info "summary written: $summary"
