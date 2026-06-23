param(
    [Parameter(Mandatory = $true)][string]$ProjectRoot,
    [Parameter(Mandatory = $true)][string]$Destination
)

$ErrorActionPreference = "Stop"
$ProjectRoot = [IO.Path]::GetFullPath($ProjectRoot)
$Destination = [IO.Path]::GetFullPath($Destination)

if (Test-Path -LiteralPath $Destination) {
    Remove-Item -LiteralPath $Destination -Recurse -Force
}

New-Item -ItemType Directory -Force -Path `
    $Destination,
    (Join-Path $Destination "INPUT"),
    (Join-Path $Destination "config") | Out-Null

foreach ($File in @(
    "run.ps1", "run_personalized.ps1",
    "README.md", "README.txt", ".gitignore"
)) {
    Copy-Item -LiteralPath (Join-Path $ProjectRoot $File) -Destination $Destination -Force
}
foreach ($Directory in @("scripts", "resources", "personalized")) {
    $Path = Join-Path $ProjectRoot $Directory
    if (Test-Path -LiteralPath $Path) {
        Copy-Item -LiteralPath $Path -Destination $Destination -Recurse -Force
    }
}
# Keep reusable downloaded databases and network caches, but omit transient logs.
$TemplateCacheLogs = Join-Path $Destination "resources/cache/logs"
if (Test-Path -LiteralPath $TemplateCacheLogs) {
    Remove-Item -LiteralPath $TemplateCacheLogs -Recurse -Force
}
$ProjectStringCache = Join-Path $Destination "resources/cache/stringdb/string_edges_project.csv"
if (Test-Path -LiteralPath $ProjectStringCache) {
    Remove-Item -LiteralPath $ProjectStringCache -Force
}
foreach ($GeneratedDir in @("personalized/work", "personalized/results")) {
    $Path = Join-Path $Destination $GeneratedDir
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Force -Path $Path | Out-Null
}
Copy-Item -LiteralPath (Join-Path $ProjectRoot "config/theme.yml") `
    -Destination (Join-Path $Destination "config/theme.yml") -Force
Copy-Item -LiteralPath (Join-Path $ProjectRoot "config/plots.yml") `
    -Destination (Join-Path $Destination "config/plots.yml") -Force
Copy-Item -LiteralPath (Join-Path $ProjectRoot "config/stage_checkpoint.example.json") `
    -Destination (Join-Path $Destination "config/stage_checkpoint.example.json") -Force

$ReadmeMd = Join-Path $Destination "README.md"
(Get-Content -LiteralPath $ReadmeMd -Raw).Replace(
    "# FKBP1A-PAI bulk RNA-seq workflow",
    "# Reusable bulk RNA-seq workflow template"
) | Set-Content -LiteralPath $ReadmeMd -Encoding UTF8
$ReadmeTxt = Join-Path $Destination "README.txt"
(Get-Content -LiteralPath $ReadmeTxt -Raw).Replace(
    "FKBP1A-PAI BULK RNA-SEQ WORKFLOW",
    "REUSABLE BULK RNA-SEQ WORKFLOW TEMPLATE"
) | Set-Content -LiteralPath $ReadmeTxt -Encoding UTF8

$TemplateConfig = @'
@{
    ProjectName    = "CHANGE_ME"
    AnalysisProfile = "Primary_padj0.05_LFC1.0"
    CountMatrix    = "INPUT/matrix_gene.count.xls"
    Metadata       = "INPUT/sample_metadata.csv"
    WorkRoot       = "work"
    ResultsRoot    = "results"
    SharedCache    = "resources/cache"
    DefaultPadjThreshold = 0.05
    DefaultLfcThreshold  = 1.0

    ControlGroup   = "Control"
    TreatGroup     = "Treatment"
    ExcludeSamples = ""
    BatchColumn    = ""
    MinCount       = 10
    MinSamples     = 3
    HeatmapTopPerDirection = 20

    PrepareCache   = $false
    EnablePathview = $false
}
'@
$TemplateConfig | Set-Content -LiteralPath (Join-Path $Destination "config/project.psd1") -Encoding UTF8
$TemplateConfig.Replace(
    'AnalysisProfile = "Primary_padj0.05_LFC1.0"',
    'AnalysisProfile = "Exploratory_padj0.05_LFC0.5"'
).Replace(
    'DefaultLfcThreshold  = 1.0',
    'DefaultLfcThreshold  = 0.5'
) | Set-Content -LiteralPath (Join-Path $Destination "config/project.exploratory.psd1") -Encoding UTF8

$Metadata = @'
sample,group,batch
Control_1,Control,
Control_2,Control,
Control_3,Control,
Treatment_1,Treatment,
Treatment_2,Treatment,
Treatment_3,Treatment,
'@
$Metadata | Set-Content -LiteralPath (Join-Path $Destination "INPUT/sample_metadata.csv") -Encoding UTF8

New-Item -ItemType File -Force -Path (Join-Path $Destination "INPUT/PUT_COUNT_MATRIX_HERE.txt") | Out-Null

Write-Host "Reusable template created: $Destination"

