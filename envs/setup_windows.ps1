<#
.SYNOPSIS
    Build the `scrna-ipf` analysis environment on Windows.

.DESCRIPTION
    1. creates the conda environment from envs/environment.yml
    2. installs PyTorch with CUDA support (used by scvi-tools on an NVIDIA GPU)
    3. installs the scverse / downstream stack from envs/requirements.txt
    4. installs the project package (src/ipfatlas) in editable mode
    5. registers a Jupyter kernel and writes a frozen lock file

    All caches are redirected to the project / J: drive so nothing large is
    written to C:.

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File envs\setup_windows.ps1
    powershell -ExecutionPolicy Bypass -File envs\setup_windows.ps1 -EnvPrefix D:\envs\scrna-ipf -Cuda cpu
#>
param(
    [string]$Conda     = "J:\Programsv1\Anaconda\Scripts\conda.exe",
    [string]$EnvPrefix = "J:\Programsv1\Anaconda_envs\scrna-ipf",
    [ValidateSet("cu126", "cu124", "cpu")]
    [string]$Cuda      = "cu126"
)
$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $PSScriptRoot

# keep caches off C:
$env:PIP_CACHE_DIR = Join-Path $Root ".cache\pip"
$env:TMP  = Join-Path $Root ".cache\tmp"
$env:TEMP = $env:TMP
New-Item -ItemType Directory -Force $env:PIP_CACHE_DIR, $env:TMP | Out-Null

Write-Host "==> [1/5] conda environment at $EnvPrefix"
if (-not (Test-Path "$EnvPrefix\python.exe")) {
    & $Conda env create --prefix $EnvPrefix --file "$Root\envs\environment.yml"
    if ($LASTEXITCODE -ne 0) { throw "conda env create failed" }
}
$Py = "$EnvPrefix\python.exe"

Write-Host "==> [2/5] PyTorch ($Cuda)"
& $Py -m pip install --upgrade pip
& $Py -m pip install torch --index-url "https://download.pytorch.org/whl/$Cuda"
if ($LASTEXITCODE -ne 0) { throw "torch install failed" }

Write-Host "==> [3/5] scverse stack"
& $Py -m pip install -r "$Root\envs\requirements.txt"
if ($LASTEXITCODE -ne 0) { throw "requirements install failed" }

Write-Host "==> [4/5] project package (editable)"
& $Py -m pip install -e $Root --no-deps

Write-Host "==> [5/5] Jupyter kernel + lock file"
& $Py -m ipykernel install --sys-prefix --name scrna-ipf --display-name "Python (scrna-ipf)"  # kernel spec lives inside the env (no C: writes)
New-Item -ItemType Directory -Force "$Root\envs\lock" | Out-Null
& $Py -m pip freeze | Out-File -Encoding utf8 "$Root\envs\lock\pip_freeze_windows.txt"
& $Conda list --prefix $EnvPrefix --explicit | Out-File -Encoding utf8 "$Root\envs\lock\conda_explicit_windows.txt"

& $Py -c "import torch, scanpy, scvi; print('torch', torch.__version__, 'CUDA', torch.cuda.is_available()); print('scanpy', scanpy.__version__, '| scvi-tools', scvi.__version__)"
Write-Host "Environment ready: $EnvPrefix"
