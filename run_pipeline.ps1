<#
.SYNOPSIS
    Reproduce the full Windows analysis: environment -> data -> notebooks -> figures.

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File run_pipeline.ps1                  # everything
    powershell -ExecutionPolicy Bypass -File run_pipeline.ps1 -SkipSetup       # env already built
    powershell -ExecutionPolicy Bypass -File run_pipeline.ps1 -From 05         # resume at notebook 05
#>
param(
    [string]$EnvPrefix = "J:\Programsv1\Anaconda_envs\scrna-ipf",
    [string]$Rscript   = "K:\R_setup\R\R-4.6.1\bin\x64\Rscript.exe",
    [switch]$SkipSetup,
    [string]$From = ""
)
$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

if (-not $SkipSetup) {
    Write-Host "==> Python environment"; & "$PSScriptRoot\envs\setup_windows.ps1" -EnvPrefix $EnvPrefix
    Write-Host "==> R packages";         & $Rscript "envs\r_packages.R"
}
. "$PSScriptRoot\envs\activate_env.ps1" -EnvPrefix $EnvPrefix

Write-Host "==> Unit tests";             python -m pytest -q -p no:cacheprovider tests
Write-Host "==> Download public data";   python workflow\scripts\download_data.py
Write-Host "==> Notebooks"
if ($From) { python workflow\scripts\run_notebooks.py --from $From } else { python workflow\scripts\run_notebooks.py }
Write-Host "Done. Figures: results\figures  Tables: results\tables  Notebooks: notebooks\*.ipynb"
