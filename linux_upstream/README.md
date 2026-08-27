# Linux upstream package — FASTQ → counts, ambient-RNA removal and gene-regulatory networks

This folder contains the parts of the IPF single-cell project that need Linux:

| step | tool | why Linux |
|---|---|---|
| raw reads → count matrices | **STARsolo** 2.7.11b | STAR is Linux/macOS only and needs ~32 GB RAM |
| ambient RNA removal | **CellBender** 0.3.2 | needs the *raw* droplet matrix (only available after re-alignment) and a CUDA GPU |
| gene regulatory networks | **pySCENIC** 0.12.1 | multiprocessing-heavy; officially supported on Linux only |

The main analysis (QC → integration → annotation → statistics → figures) runs on
Windows from the notebooks. This package produces **additional evidence** that
the notebooks pick up automatically once the results are copied back:

* **Notebook 14** compares our STARsolo counts with the published Cell Ranger
  matrix, quantifies ambient RNA per cell type (CellBender) and runs RNA velocity
  on the reprocessed libraries.
* **Notebook 15** analyses SCENIC regulons along the alveolar epithelial lineage (AT2 → transitional → AT1 / aberrant basaloid).

---

## 0. Getting this package onto the Linux machine

The package is sent from the Windows PC with Tailscale file sharing (Taildrop),
which needs no SSH server. On `ahmedb`:

```bash
mkdir -p ~/incoming && sudo tailscale file get ~/incoming      # receives the waiting files
cd ~/incoming && sha256sum -c ipf_linux_upstream.tar.gz.sha256  # integrity check
tar -xzf ipf_linux_upstream.tar.gz -C ~                         # creates ~/ipf_linux_upstream
cd ~/ipf_linux_upstream && chmod +x run_all.sh scripts/*.sh
```

(`sudo tailscale set --operator=$USER` once removes the need for `sudo`.)

`input/scenic_input.loom` (for step 08) is sent separately once Notebook 09
has run on Windows; move it into `~/ipf_linux_upstream/input/`.

## 1. What is processed

`config/samples.tsv` lists **8 donors with one 10x library each** from GSE135893
(Habermann *et al.* 2020): 4 controls and 4 IPF, ~96 GB of FASTQ from ENA.

| library | donor | diagnosis | instrument | FASTQ |
|---|---|---|---|---|
| F01394 | VUHD71 | Control | NovaSeq 6000 | 5.9 GB |
| F00409 | VUHD69 | Control | HiSeq 4000 | 6.6 GB |
| F01174 | THD0002 | Control | NovaSeq 6000 | 9.9 GB |
| F01157 | THD0001 | Control | NovaSeq 6000 | 17.3 GB |
| F01214 | TILD010 | IPF | NovaSeq 6000 | 4.1 GB |
| F01390 | VUILD64 | IPF | NovaSeq 6000 | 6.0 GB |
| F00431 | TILD001 | IPF | HiSeq 4000 | 7.3 GB |
| ILD63 | VUILD63 | IPF | NovaSeq 6000 | 38.6 GB |

Donors with two libraries were left out on purpose: their runs are filed under
one GEO sample, and merging two 10x channels in one STARsolo run would mix
barcodes from different captures.

## 2. Requirements

* 64-bit Linux, bash ≥ 4, `curl`, internet access to `ebi.ac.uk`, `10xgenomics.com`, `aertslab.org`
* ≥ 16 CPU threads recommended, **≥ 40 GB RAM** (32 GB works with a sparse STAR index)
* **~210 GB free disk** in `PROJECT_DIR` (FASTQ 96 GB + reference 11 GB + index 30 GB + outputs)
* NVIDIA GPU with ≥ 8 GB for CellBender (optional: without one, set `USE_GPU=no` and run
  CellBender on the analysis PC instead, see section 6; CPU mode works but takes hours per library)
* conda/mamba (optional — Miniforge is installed automatically if missing, no root needed).
  An existing installation is used when `CONDA_ROOT` points to it; conda ≥ 23 and mamba ≥ 2 both work.

Run `bash scripts/00_check_system.sh` first; it reports all of the above.

## 3. Quick start

```bash
# unpack the package (see section 0) and enter it
cd ~/ipf_linux_upstream

# where the big files go (default: ~/ipf_scrna_upstream)
export PROJECT_DIR=/data/ipf_scrna_upstream
export THREADS=16

# optional: an existing conda/mamba, with environments and package cache on a data
# drive instead of the home directory
export CONDA_ROOT=$HOME/anaconda3
export ENVS_DIR=/data/conda/envs PKGS_DIR=/data/conda/pkgs

# a 10-minute smoke test on 2 M reads of one library first; its reads, counts and
# logs go to $PROJECT_DIR/test_mode/, only the reference is shared with the full run
TEST_MODE=1 LIBRARIES="F01214" bash run_all.sh 00 01 02 03 04 05 06 09 10

# then the full run, detached from the terminal (tmux or nohup)
tmux new-session -d -s ipf 'bash run_all.sh 2>&1 | tee -a run_all.log'
tail -f run_all.log
```

Without an NVIDIA GPU, add `USE_GPU=no SKIP_CELLBENDER=1` to both commands: step 07
is then skipped, step 01 does not build the CellBender environment, and step 10
returns the raw droplet matrices for CellBender on the analysis PC (section 6).

Steps can be run one by one (`bash run_all.sh 06`) or directly
(`bash scripts/06_starsolo.sh`). Every step skips work that is already done,
so a crashed run is resumed by starting it again.

## 4. Steps in detail

| # | script | what it does | approx. time (16 threads) |
|---|---|---|---|
| 00 | `00_check_system.sh` | CPU, RAM, disk, GPU, conda and network check | seconds |
| 01 | `01_install_envs.sh` | creates `ipf-upstream`, `ipf-scenic` and (unless `SKIP_CELLBENDER=1`) `ipf-cellbender` conda envs, by path in `ENVS_DIR` if set | 10–20 min |
| 02 | `02_get_reference.sh` | 10x GRCh38-3.0.0 reference (same as the original study) + 737K barcode whitelist, size/MD5 verified | 15 min |
| 03 | `03_build_star_index.sh` | STAR genome index (`--sjdbOverhang 100`) | ~1 h |
| 04 | `04_download_fastq.sh` | FASTQ from ENA over HTTPS, MD5-verified, resumable, 4 in parallel | 1–3 h |
| 05 | `05_read_qc.sh` | detects which FASTQ holds the barcode (whitelist hit rate), FastQC on subsamples | 15 min |
| 06 | `06_starsolo.sh` | STARsolo: Gene, GeneFull and Velocyto counts, EmptyDrops cell calling | 20–60 min / library |
| 07 | `07_cellbender.sh` | CellBender remove-background on the raw droplet matrix (skipped with `USE_GPU=no`) | 20–40 min / library (GPU) |
| 08 | `08_pyscenic.sh` | GRNBoost2 → cisTarget → AUCell, ×3 seeds, consensus regulons | 3–8 h |
| 09 | `09_export_h5ad.sh` | one `.h5ad` per library with spliced/unspliced layers | minutes |
| 10 | `10_collect_results.sh` | MultiQC report, checksums, single `.tar` archive to return | minutes |

### Key parameter choices

* **`--soloStrand Reverse`**: in 10x 5′ libraries, Read 2 (cDNA) is antisense to
  the transcript. Using `Forward` here would assign almost no reads to genes.
* **Barcode-read detection** (`python/detect_read_layout.py`): ENA does not say
  which FASTQ is Read 1. The script tests the first 16 bases of 20,000 reads
  against the whitelist; the barcode read has ~90 % hits, the cDNA read ~0 %.
* **`_CR` options + `EmptyDrops_CR`**: reproduce Cell Ranger ≥ 4 barcode
  correction, UMI deduplication and cell calling, so counts are comparable with GEO.
* **Velocyto counts from 5′ data**: 5′ libraries capture fewer intronic reads
  than 3′ libraries. RNA velocity on these data is treated as **exploratory**
  (Notebook 14 reports the unspliced fraction per cell type before interpreting anything).
* **pySCENIC ×3**: GRNBoost2 is stochastic; only regulons recovered in ≥ 2 of 3
  runs are kept (`python/aggregate_regulons.py`), following Van de Sande *et al.* 2020.

## 5. pySCENIC input

Step 08 needs `input/scenic_input.loom`, which is exported by
**Notebook 09** on Windows (raw counts of the subsampled alveolar epithelial lineage, genes detected in ≥ 1 % of cells).
If the file is missing, `run_all.sh` skips step 08 with a warning. To run it
later, copy the loom into `input/` and run `bash run_all.sh 08 10`.

## 6. Bringing the results back to Windows

Step 10 writes `$PROJECT_DIR/ipf_upstream_results_<date>.tar` (+ `.sha256`).
Choose one of these to copy it:

```bash
# a) Tailscale file transfer (Taildrop) to the Windows PC "ahmed"
tailscale file cp $PROJECT_DIR/ipf_upstream_results_*.tar ahmed:

# b) scp from the Windows side (PowerShell), if SSH is enabled on this machine
scp <user>@ahmedb:/data/ipf_scrna_upstream/ipf_upstream_results_*.tar J:\Portfolio_projects\Projects\scRNA\data\linux_results\
```

On Windows, extract the archive into `data/linux_results/` (so that
`data/linux_results/results_to_return/` exists), then run notebooks 14 and 15:

```powershell
tar -xf ipf_upstream_results_<date>.tar -C J:\Portfolio_projects\Projects\scRNA\data\linux_results
python workflow\scripts\run_notebooks.py 14 15
```

If CellBender was skipped on Linux (`USE_GPU=no`), the archive holds the raw droplet
matrices (`starsolo/<library>/Gene/raw/`). Run CellBender on the analysis PC's GPU
first, with the parameters of step 07, in its own environment (`envs/cellbender.yml`;
CellBender pins Jupyter < 7, which would clash with the analysis environment):

```powershell
conda env create -p J:\Programsv1\Anaconda_envs\cellbender -f envs\cellbender.yml   # once
python workflow\scripts\run_cellbender.py        # writes results_to_return\cellbender\<library>\
python workflow\scripts\run_notebooks.py 14 15
```

The environment location is set in `config/config.yaml` (`cellbender: env_windows`).

## 7. Folder layout

```
linux_upstream/
├── run_all.sh                 # orchestrator (step selection, logging, timing)
├── config/
│   ├── config.sh              # ALL paths, URLs, threads, chemistry settings
│   └── samples.tsv            # libraries, ENA run accessions, FASTQ URLs + MD5
├── envs/                      # conda environment definitions
│   ├── upstream.yml  cellbender.yml  scenic.yml
├── scripts/                   # one script per step (00-10) + common.sh helpers
├── python/                    # helper programs used by the steps
│   ├── detect_read_layout.py  starsolo_to_h5ad.py  aggregate_regulons.py
└── input/                     # place scenic_input.loom here
```

## 8. Troubleshooting

| symptom | fix |
|---|---|
| `EXITING because of FATAL ERROR: not enough memory for BAM sorting` | not used here (`--outSAMtype None`); if you enable BAM, add `--limitBAMsortRAM 30000000000` |
| STAR killed during index build | RAM too small → `MEM_GB=30 bash scripts/03_build_star_index.sh` builds a sparse index |
| `md5 mismatch` repeatedly | ENA mirror hiccup; rerun step 04 (it resumes) |
| very low "Reads Mapped to Gene" in `Summary.csv` | wrong strand — check `SOLO_STRAND` (Reverse for 5′, Forward for 3′) |
| CellBender `CUDA out of memory` | add `--projected-ambient-count-threshold 1` or lower `--total-droplets-included` in step 07 |
| pySCENIC `AttributeError: module 'numpy' has no attribute 'object'` | env built with numpy ≥ 1.24 → recreate `ipf-scenic` from `envs/scenic.yml` |
