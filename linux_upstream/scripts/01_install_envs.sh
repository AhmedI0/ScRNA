#!/usr/bin/env bash
# Step 01 - create the three conda environments.
#   ipf-upstream   : STAR 2.7.11b, samtools, seqkit, FastQC, MultiQC, anndata
#   ipf-cellbender : CellBender 0.3.2 (+ PyTorch)
#   ipf-scenic     : pySCENIC 0.12.1
# If no conda is available, Miniforge is installed into $CONDA_ROOT (no root needed).
source "$(dirname "$0")/common.sh"

if ! conda_cmd >/dev/null; then
    info "installing Miniforge into $CONDA_ROOT"
    mkdir -p "$PROJECT_DIR"
    inst="$PROJECT_DIR/Miniforge3-Linux-x86_64.sh"
    curl -fsSL -o "$inst" "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh"
    bash "$inst" -b -p "$CONDA_ROOT"
    rm -f "$inst"
fi
C="$(conda_cmd)"
info "using $C"

env_exists() {
    if [[ -n "$ENVS_DIR" ]]; then [[ -d "$ENVS_DIR/$1/conda-meta" ]]
    else "$C" env list | awk '{print $1}' | grep -qx "$1"; fi
}

# --yes where the front-end has it: mamba >= 2 asks for confirmation, which would stall
# an unattended run; older conda neither asks nor knows the flag
yes_flag() { "$C" env "$1" --help 2>&1 | grep -q -- "--yes" && echo "--yes"; }

create_env() {  # create_env <name> <yml>
    if env_exists "$1"; then
        info "environment $1 exists - updating"
        "$C" env update $(env_target "$1") -f "$2" --prune $(yes_flag update)
    else
        info "creating environment $1${ENVS_DIR:+ in $ENVS_DIR}"
        if [[ -n "$ENVS_DIR" ]]; then mkdir -p "$ENVS_DIR"; fi
        "$C" env create $(env_target "$1") -f "$2" $(yes_flag create)
    fi
}

create_env "$ENV_UPSTREAM"   "$PKG_DIR/envs/upstream.yml"
create_env "$ENV_SCENIC"     "$PKG_DIR/envs/scenic.yml"
if [[ "${SKIP_CELLBENDER:-0}" != "1" ]]; then
    create_env "$ENV_CELLBENDER" "$PKG_DIR/envs/cellbender.yml"
fi

info "tool versions:"
in_env "$ENV_UPSTREAM" STAR --version
in_env "$ENV_UPSTREAM" samtools --version | head -1
in_env "$ENV_SCENIC" python -c "import pyscenic; print('pySCENIC', pyscenic.__version__)"
[[ "${SKIP_CELLBENDER:-0}" == "1" ]] || in_env "$ENV_CELLBENDER" cellbender --version || true

mkdir -p "$LOG_DIR/envs"
for e in "$ENV_UPSTREAM" "$ENV_SCENIC" "$ENV_CELLBENDER"; do
    env_exists "$e" || continue
    "$C" env export $(env_target "$e") > "$LOG_DIR/envs/$e.lock.yml" 2>/dev/null || true
done
info "environments ready (exported locks in $LOG_DIR/envs)"
