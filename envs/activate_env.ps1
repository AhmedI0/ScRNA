<#
.SYNOPSIS
    Activate the scrna-ipf environment in the current PowerShell session.
.DESCRIPTION
    Dot-source this file:  . .\envs\activate_env.ps1
    It prepends the environment folders to PATH (equivalent to `conda activate`),
    which is required on Windows so that MKL loads its DLLs from the environment,
    and redirects temporary files to the project's .cache folder.
#>
param([string]$EnvPrefix = "J:\Programsv1\Anaconda_envs\scrna-ipf")
$Root = Split-Path -Parent $PSScriptRoot
$dirs = @($EnvPrefix, "$EnvPrefix\Library\mingw-w64\bin", "$EnvPrefix\Library\usr\bin",
          "$EnvPrefix\Library\bin", "$EnvPrefix\Scripts", "$EnvPrefix\bin")
$env:PATH = ($dirs -join ";") + ";" + $env:PATH
$env:CONDA_PREFIX = $EnvPrefix
$env:PYTHONUTF8 = "1"
$env:TMP = Join-Path $Root ".cache\tmp"; $env:TEMP = $env:TMP
New-Item -ItemType Directory -Force $env:TMP | Out-Null
Write-Host "scrna-ipf active: $(python --version)  ($EnvPrefix)"
