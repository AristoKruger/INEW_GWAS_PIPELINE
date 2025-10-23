<#
Runs the PLINK sequence described in the README (LD pruning, binary conversion,
pruned subset, and PED/MAP export). Execute from the project root.

Usage examples:
  # Uses PLINK on PATH
  powershell -ExecutionPolicy Bypass -File scripts/cli_examples/run_plink_stage.ps1

  # Or define an explicit PLINK path for this session
  $env:PLINK_PATH = "C:\Tools\plink\plink.exe"
  powershell -ExecutionPolicy Bypass -File scripts/cli_examples/run_plink_stage.ps1
#>

$plink = $env:PLINK_PATH
if ([string]::IsNullOrWhiteSpace($plink)) {
    $plink = "plink"
}

if (-not (Get-Command $plink -ErrorAction SilentlyContinue)) {
    throw "PLINK executable not found. Either add it to PATH or set the PLINK_PATH environment variable."
}

$pythonCmd = Get-Command python -ErrorAction SilentlyContinue
if (-not $pythonCmd) {
    throw "Python executable not found on PATH. Required to read config/pipeline.yaml."
}
$pythonExe = $pythonCmd.Source
$yamlCmd = @"
import json
from pathlib import Path
import yaml

cfg = yaml.safe_load(Path('config/pipeline.yaml').read_text())
params = cfg.get('tools', {}).get('plink_ld_prune', {})
step11 = cfg.get('steps', {}).get('step11_plink_pruning', {})
print(json.dumps({
    'window_kb': params.get('window_kb', 50),
    'step': params.get('step', 5),
    'r2': params.get('r2', 0.2),
    'panel_prefix': step11.get('panel_prefix', 'outputs/step11_plink_pruning/pruned_panel'),
}))
"@
$rawParams = & $pythonExe "-c" $yamlCmd
if ($LASTEXITCODE -ne 0 -or -not $rawParams) {
    throw "Failed to read PLINK LD prune parameters from config/pipeline.yaml."
}
$params = ConvertFrom-Json $rawParams
$windowKb = [int]$params.window_kb
$stepSize = [int]$params.step
$r2Threshold = [double]$params.r2
$panelPrefix = $params.panel_prefix
if ([string]::IsNullOrWhiteSpace($panelPrefix)) {
    $panelPrefix = "outputs\step11_plink_pruning\pruned_panel"
}
$panelPrefix = $panelPrefix -replace '/', '\'
$panelDir = Split-Path -Parent $panelPrefix
if (-not (Test-Path $panelDir)) {
    New-Item -ItemType Directory -Path $panelDir -Force | Out-Null
}
$tasselLogDir = Join-Path $panelDir "logs"
if (-not (Test-Path $tasselLogDir)) {
    New-Item -ItemType Directory -Path $tasselLogDir -Force | Out-Null
}

$base = "outputs\step09_csv_to_tped\filtered_genotypes_strict"
$mapPath = "outputs\step08_create_plink_map\filtered_genotypes_strict.map"
$required = @("$base.tped", "$base.tfam", $mapPath)
$missing = $required | Where-Object { -not (Test-Path $_) }
if ($missing) {
    $list = ($missing | ForEach-Object { "  - $_" }) -join "`n"
    throw "Missing required input file(s):`n$list`nRun `python -m gwas_pipeline --steps step08_create_plink_map step09_csv_to_tped` before this helper."
}

$ldPrefix = $panelPrefix
$qcPrefix = Join-Path $panelDir "filtered_qc"

function Invoke-PlinkStep {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Message,

        [Parameter(Mandatory = $true)]
        [string[]]$Arguments,

        [string]$ErrorNote
    )

    Write-Host $Message -ForegroundColor Cyan
    & $plink @Arguments
    if ($LASTEXITCODE -ne 0) {
        if ($ErrorNote) {
            throw $ErrorNote
        }
        throw "$Message failed."
    }
}

function Move-PlinkLogs {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Source,

        [Parameter(Mandatory = $true)]
        [string]$Destination,

        [string]$Pattern = "*.log"
    )

    Get-ChildItem -Path $Source -Filter $Pattern -ErrorAction SilentlyContinue |
        ForEach-Object {
            Move-Item -Path $_.FullName -Destination (Join-Path $Destination $_.Name) -Force
        }
}

Invoke-PlinkStep -Message "Running LD pruning..." -Arguments @(
    "--tfile", $base,
    "--allow-extra-chr",
    "--indep-pairwise", $windowKb, $stepSize, $r2Threshold,
    "--out", $ldPrefix
) -ErrorNote "LD pruning failed."

Invoke-PlinkStep -Message "Creating binary files from QC dataset..." -Arguments @(
    "--tfile", $base,
    "--allow-extra-chr",
    "--make-bed",
    "--out", $qcPrefix
) -ErrorNote "QC binary conversion failed."

Invoke-PlinkStep -Message "Restricting to pruned SNPs..." -Arguments @(
    "--bfile", $qcPrefix,
    "--extract", ($panelPrefix + ".prune.in"),
    "--make-bed",
    "--allow-extra-chr",
    "--out", $panelPrefix
) -ErrorNote "Pruned binary conversion failed."

Invoke-PlinkStep -Message "Generating PED/MAP for TASSEL..." -Arguments @(
    "--bfile", $panelPrefix,
    "--recode",
    "--allow-extra-chr",
    "--tab",
    "--out", $panelPrefix
) -ErrorNote "TASSEL PED/MAP export failed."

Move-PlinkLogs -Source $panelDir -Destination $tasselLogDir -Pattern "*.log"
Move-PlinkLogs -Source $panelDir -Destination $tasselLogDir -Pattern "*.nosex"

Write-Host "Done. Files saved under $panelDir" -ForegroundColor Green
