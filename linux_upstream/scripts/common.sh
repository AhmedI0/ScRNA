#!/usr/bin/env bash
# Shared helpers sourced by every step script.
set -euo pipefail

# shellcheck source=../config/config.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"
mkdir -p "$LOG_DIR"

ts()   { date "+%Y-%m-%d %H:%M:%S"; }
info() { echo -e "[$(ts)] \033[1;34mINFO\033[0m  $*"; }
warn() { echo -e "[$(ts)] \033[1;33mWARN\033[0m  $*" >&2; }
die()  { echo -e "[$(ts)] \033[1;31mERROR\033[0m $*" >&2; exit 1; }

# Locate a conda front-end (mamba > micromamba > conda), preferring our own Miniforge.
conda_cmd() {
    if [[ -x "$CONDA_ROOT/bin/mamba" ]]; then echo "$CONDA_ROOT/bin/mamba"
    elif [[ -x "$CONDA_ROOT/bin/conda" ]]; then echo "$CONDA_ROOT/bin/conda"
    elif command -v mamba >/dev/null 2>&1; then command -v mamba
    elif command -v micromamba >/dev/null 2>&1; then command -v micromamba
    elif command -v conda >/dev/null 2>&1; then command -v conda
    else return 1; fi
}

# Address an environment by path when ENVS_DIR is set, otherwise by name
env_target() {
    if [[ -n "${ENVS_DIR:-}" ]]; then echo "-p $ENVS_DIR/$1"; else echo "-n $1"; fi
}

# Run a command inside a named conda environment:  in_env ipf-upstream STAR --version
# conda needs --no-capture-output to stream output; mamba >= 2 streams by default and rejects it
in_env() {
    local env="$1"; shift
    local c; c="$(conda_cmd)" || die "no conda/mamba found - run scripts/01_install_envs.sh first"
    if [[ -z "${RUN_FLAGS+x}" ]]; then
        RUN_FLAGS=""
        "$c" run --help 2>&1 | grep -q -- "--no-capture-output" && RUN_FLAGS="--no-capture-output"
    fi
    "$c" run $RUN_FLAGS $(env_target "$env") "$@"
}

# Library ids from the sample sheet (optionally restricted by $LIBRARIES)
library_ids() {
    local all
    all=$(awk -F'\t' 'NR>1 {print $1}' "$SAMPLES_TSV" | sort -u)
    if [[ -n "$LIBRARIES" ]]; then
        for l in $LIBRARIES; do grep -qx "$l" <<<"$all" && echo "$l" || warn "library $l not in $SAMPLES_TSV"; done
    else
        echo "$all"
    fi
}

# Column value from the sample sheet by header name, for all rows of a library
sheet_col() {  # sheet_col <library_id> <column_name>
    awk -F'\t' -v lib="$1" -v col="$2" '{sub(/\r$/, "")} NR==1 {for (i=1;i<=NF;i++) if ($i==col) c=i; next} $1==lib {print $c}' "$SAMPLES_TSV"
}

require_file() { [[ -s "$1" ]] || die "required file missing: $1 ($2)"; }

has_gpu() {
    case "$USE_GPU" in
        yes) return 0 ;;
        no)  return 1 ;;
        *)   command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1 ;;
    esac
}
