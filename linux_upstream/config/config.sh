#!/usr/bin/env bash
# ============================================================================
#  Central configuration for the Linux upstream package.
#  Every script sources this file. Override any variable from the shell, e.g.
#      PROJECT_DIR=/data/ipf THREADS=32 bash run_all.sh
# ============================================================================

# Location of this package (resolved automatically)
PKG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PKG_DIR

# ---- working directory: needs ~200 GB free (FASTQ ~96 GB, reference ~45 GB) ----
export PROJECT_DIR="${PROJECT_DIR:-$HOME/ipf_scrna_upstream}"

# ---- compute ----
export THREADS="${THREADS:-$(nproc)}"
export MEM_GB="${MEM_GB:-$(awk '/MemTotal/ {printf "%d", $2/1024/1024}' /proc/meminfo)}"
export DOWNLOAD_JOBS="${DOWNLOAD_JOBS:-4}"     # parallel FASTQ downloads
export USE_GPU="${USE_GPU:-auto}"              # auto | yes | no (CellBender)

# ---- conda ----
# Uses mamba/micromamba/conda if present; otherwise 01_install_envs.sh installs
# Miniforge into $PROJECT_DIR/miniforge3.
export CONDA_ROOT="${CONDA_ROOT:-$PROJECT_DIR/miniforge3}"
# Optional: create the environments by path in this folder (e.g. on a data drive)
# instead of conda's default envs location, and keep conda's package cache in PKGS_DIR.
export ENVS_DIR="${ENVS_DIR:-}"
export PKGS_DIR="${PKGS_DIR:-}"
if [[ -n "$PKGS_DIR" ]]; then export CONDA_PKGS_DIRS="$PKGS_DIR"; fi
export ENV_UPSTREAM="ipf-upstream"
export ENV_CELLBENDER="ipf-cellbender"
export ENV_SCENIC="ipf-scenic"

# ---- inputs ----
export SAMPLES_TSV="${SAMPLES_TSV:-$PKG_DIR/config/samples.tsv}"
# Restrict to some libraries (space-separated library_id), empty = all in SAMPLES_TSV
export LIBRARIES="${LIBRARIES:-}"
# TEST_MODE=1: stream only the first TEST_READS reads of ONE run per library
export TEST_MODE="${TEST_MODE:-0}"
export TEST_READS="${TEST_READS:-2000000}"

# ---- reference: identical to the one used by the original study ----
# Cell Ranger 3.0.2 reference GRCh38-3.0.0 (Ensembl 93 annotation, filtered GTF)
export REF_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz"
export REF_BYTES=11387238361
export WHITELIST_URL="https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/737K-august-2016.txt.gz"
export WHITELIST_MD5="f62a276e262fdd85262a889d0f48556b"   # md5 of the *uncompressed* 737K list

# ---- 10x Chromium Single Cell 5' v1 read structure ----
export CB_LEN=16
export UMI_LEN=10
export SOLO_STRAND="Reverse"     # 5' kits: Read 2 is antisense to the mRNA

# ---- pySCENIC resources (Aerts lab, hg38, motif collection v10 clustered) ----
export SCENIC_BASE="https://resources.aertslab.org/cistarget"
export SCENIC_DB1="databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
export SCENIC_DB2="databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
export SCENIC_MOTIFS="motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
export SCENIC_TFS="tf_lists/allTFs_hg38.txt"
export SCENIC_INPUT="${SCENIC_INPUT:-$PKG_DIR/input/scenic_input.loom}"
export SCENIC_RUNS="${SCENIC_RUNS:-3}"       # independent GRNBoost2 runs (stochastic)
export SCENIC_SEED=777

# ---- derived paths ----
export REF_DIR="$PROJECT_DIR/reference"
export REF_NAME="refdata-cellranger-GRCh38-3.0.0"
export STAR_INDEX="$REF_DIR/star_index_GRCh38-3.0.0"
export WHITELIST="$REF_DIR/737K-august-2016.txt"
# TEST_MODE reads carry the real file names, and every step skips finished work, so test
# reads, counts and logs go to their own folder; only the reference is shared
RUN_DIR="$PROJECT_DIR"
if [[ "$TEST_MODE" == "1" ]]; then RUN_DIR="$PROJECT_DIR/test_mode"; fi
export FASTQ_DIR="$RUN_DIR/fastq"
export SOLO_DIR="$RUN_DIR/starsolo"
export CELLBENDER_DIR="$RUN_DIR/cellbender"
export SCENIC_DIR="$RUN_DIR/scenic"
export QC_DIR="$RUN_DIR/qc"
export RESULTS_DIR="$RUN_DIR/results_to_return"
export LOG_DIR="$RUN_DIR/logs"
