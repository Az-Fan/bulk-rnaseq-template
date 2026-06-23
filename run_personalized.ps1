param(
    [Parameter(Mandatory = $true)][ValidatePattern('^[A-Za-z0-9_.-]+$')]
    [string]$Name,
    [string]$Config = "config/project.psd1",
    [string]$Rscript = "C:\Program Files\R\R-4.4.2\bin\Rscript.exe"
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location -LiteralPath $root
$settings = Import-PowerShellDataFile -LiteralPath (Join-Path $root $Config)
$script = Join-Path $root "personalized/code/$Name.R"
if (-not (Test-Path -LiteralPath $script)) { throw "Personalized script not found: $script" }

$profile = [string]$settings.AnalysisProfile
$standardWork = Join-Path (Join-Path (Join-Path $root $settings.WorkRoot) $settings.ProjectName) $profile
$standardResults = Join-Path (Join-Path (Join-Path $root $settings.ResultsRoot) $settings.ProjectName) $profile
$personalizedWork = Join-Path $root "personalized/work/$Name"
$personalizedResults = Join-Path $root "personalized/results/$Name"
New-Item -ItemType Directory -Force -Path $personalizedWork, $personalizedResults | Out-Null

$env:R_LIBS_USER = Join-Path $root ".r-library"
$env:PERSONALIZED_WORK_DIR = $personalizedWork
$env:PERSONALIZED_OUTPUT_DIR = $personalizedResults
$env:STANDARD_WORK_DIR = $standardWork
$env:STANDARD_RESULTS_DIR = $standardResults

& $Rscript $script
if ($LASTEXITCODE -ne 0) { throw "Personalized analysis failed: $Name" }
Write-Host "Personalized results: $personalizedResults"
