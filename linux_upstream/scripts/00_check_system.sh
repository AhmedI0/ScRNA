#!/usr/bin/env bash
# Step 00 - check that the workstation can run the pipeline.
# Prints CPU / RAM / disk / GPU / conda / network status and the resource
# estimate for the selected samples. Nothing is installed or changed.
source "$(dirname "$0")/common.sh"

info "Host: $(hostname)  |  $(uname -srm)"
[[ -r /etc/os-release ]] && info "OS: $(. /etc/os-release; echo "$PRETTY_NAME")"
info "CPU threads: $(nproc)  -> THREADS=$THREADS"
info "RAM: ${MEM_GB} GB  (STAR genome index needs ~32 GB; with <40 GB a sparse index is built)"

mkdir -p "$PROJECT_DIR"
avail=$(df -BG --output=avail "$PROJECT_DIR" | tail -1 | tr -dc '0-9')
need=$(awk -F'\t' '{sub(/\r$/, "")} NR==1 {for (i=1;i<=NF;i++) if ($i=="gb") c=i; next} {s+=$c} END {printf "%d", s + 110}' "$SAMPLES_TSV")
info "Disk free at $PROJECT_DIR: ${avail} GB  |  estimated need: ~${need} GB (FASTQ + reference + index + outputs)"
(( avail > need )) || warn "not enough free disk space for all samples - set PROJECT_DIR or LIBRARIES"

if command -v nvidia-smi >/dev/null 2>&1; then
    info "GPU: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | paste -sd ';')"
else
    warn "no NVIDIA GPU detected: CellBender will run on CPU (slow, hours per library) or set USE_GPU=no to skip"
fi

if c=$(conda_cmd); then info "conda front-end: $c"; else warn "no conda found: 01_install_envs.sh will install Miniforge in $CONDA_ROOT"; fi

for url in "https://www.ebi.ac.uk/ena/portal/api/" "$REF_URL" "$SCENIC_BASE/$SCENIC_TFS"; do
    if curl -sfIL --max-time 20 "$url" >/dev/null; then info "reachable: $url"; else warn "NOT reachable: $url"; fi
done

info "Samples in $SAMPLES_TSV:"
awk -F'\t' '{sub(/\r$/, "")} NR>1 {gb[$1]+=$NF; d[$1]=$2"\t"$3} END {for (l in gb) printf "    %-10s %-10s %6.1f GB\n", l, d[l], gb[l]}' "$SAMPLES_TSV" | sort -k3
info "system check finished"
