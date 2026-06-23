param(
    [ValidateSet(
        "check", "01_qc", "02_differential", "03_enrichment",
        "04_regulation", "05_ppi", "06_custom", "organize"
    )]
    [string]$Stage = "01_qc",
    [string]$Config = "config/project.psd1",
    [string]$Rscript = "",
    [switch]$Overwrite
)

$ErrorActionPreference = "Stop"
$ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location -LiteralPath $ProjectRoot

function Resolve-ProjectPath([string]$Path) {
    if ([IO.Path]::IsPathRooted($Path)) { return [IO.Path]::GetFullPath($Path) }
    return [IO.Path]::GetFullPath((Join-Path $ProjectRoot $Path))
}

function Remove-SafeDirectory([string]$Path, [string]$Parent) {
    $ResolvedPath = [IO.Path]::GetFullPath($Path)
    $ResolvedParent = [IO.Path]::GetFullPath($Parent)
    if (-not $ResolvedPath.StartsWith($ResolvedParent + [IO.Path]::DirectorySeparatorChar)) {
        throw "Refusing to remove directory outside its configured root: $ResolvedPath"
    }
    if (Test-Path -LiteralPath $ResolvedPath) {
        Remove-Item -LiteralPath $ResolvedPath -Recurse -Force
    }
}

function Optimize-PipelineArtifacts([string]$PipelineRoot, [switch]$ForCache) {
    if (-not (Test-Path -LiteralPath $PipelineRoot)) { return }

    # Remove non-deliverable and duplicate representations before copying.
    Get-ChildItem -LiteralPath $PipelineRoot -Recurse -File -Filter *.md |
        Remove-Item -Force
    Get-ChildItem -LiteralPath $PipelineRoot -Recurse -File -Filter *.pdf |
        Where-Object {
            Test-Path -LiteralPath ([IO.Path]::ChangeExtension($_.FullName, ".png"))
        } | Remove-Item -Force
    Get-ChildItem -LiteralPath $PipelineRoot -Recurse -File |
        Where-Object {
            ($_.Name -eq "Excluded_Samples.txt" -and $_.Length -eq 0) -or
            $_.Name -in @(
                "recipe.json", "run_status.tsv", "session_info_command.log",
                "figure_registry.csv", "Analysis_Summary_Report.md"
            )
        } | Remove-Item -Force
    Get-ChildItem -LiteralPath $PipelineRoot -Recurse -Directory |
        Sort-Object FullName -Descending |
        Where-Object { -not (Get-ChildItem -LiteralPath $_.FullName -Force) } |
        Remove-Item -Force
    if (-not $ForCache) { return }

    # Final visual deliverables already live in results. Work keeps only
    # reusable numeric/intermediate objects needed to redraw without refitting.
    Get-ChildItem -LiteralPath $PipelineRoot -Recurse -File |
        Where-Object {
            $_.Extension.ToLowerInvariant() -in @(
                ".png", ".pdf", ".html", ".htm", ".svg", ".jpg", ".jpeg",
                ".tif", ".tiff"
            )
        } | Remove-Item -Force
    Get-ChildItem -LiteralPath $PipelineRoot -Recurse -Directory |
        Where-Object { $_.Name -like "*_files" } |
        Sort-Object FullName -Descending |
        Remove-Item -Recurse -Force
    Get-ChildItem -LiteralPath (Join-Path $PipelineRoot "functional_enrichment") `
        -Recurse -File -Filter *.xlsx -ErrorAction SilentlyContinue |
        Remove-Item -Force

    # Keep only upstream objects actually read by later stages or needed to
    # redraw from stored numerical results.
    $dataKeep = @(
        "Analysis_Config_Used.csv",
        "metadata.Rdata",
        "metadata.csv",
        "exprSet_vst.Rdata",
        "exprSet_vst_filtered.Rdata",
        "exprSet_vst_unfiltered.Rdata",
        "exprSet_for_heatmap.Rdata"
    )
    Get-ChildItem -LiteralPath (Join-Path $PipelineRoot "data_processed") `
        -File -ErrorAction SilentlyContinue |
        Where-Object {
            $_.Name -notin $dataKeep -and
            $_.Name -notlike "exprSet_vst_batch_corrected_*.Rdata"
        } | Remove-Item -Force

    $diffKeep = @(
        "DEseq2_Diff_Annotated.Rdata",
        "DEG_results_annotated.csv"
    )
    Get-ChildItem -LiteralPath (Join-Path $PipelineRoot "diff_analysis") `
        -File -ErrorAction SilentlyContinue |
        Where-Object { $_.Name -notin $diffKeep } |
        Remove-Item -Force

    # Detailed logs are delivered to results/99_Run. Work retains compact
    # stage status, recipe, completion markers and session info only.
    Get-ChildItem -LiteralPath (Join-Path $PipelineRoot "logs") `
        -File -Filter *.log -ErrorAction SilentlyContinue |
        Remove-Item -Force
    Remove-Item -LiteralPath (Join-Path $PipelineRoot "summary") `
        -Recurse -Force -ErrorAction SilentlyContinue

    # Raw extracted counts can always be reconstructed from INPUT and are not
    # required by downstream plotting.
    Remove-Item -LiteralPath (Join-Path $PipelineRoot "data_processed/counts_sub_confwe.csv") `
        -Force -ErrorAction SilentlyContinue
    Remove-Item -LiteralPath (Join-Path $PipelineRoot "data_processed/plot_data_long_for_selected_boxplots.Rdata") `
        -Force -ErrorAction SilentlyContinue

    Get-ChildItem -LiteralPath $PipelineRoot -Recurse -Directory |
        Sort-Object FullName -Descending |
        Where-Object { -not (Get-ChildItem -LiteralPath $_.FullName -Force) } |
        Remove-Item -Force
}

function Restore-PipelineWorkspace([string]$ProfileDir, [string]$RuntimeDir) {
    $LegacyDir = Join-Path $ProfileDir "_pipeline"
    if ((Test-Path -LiteralPath $LegacyDir) -and -not (Test-Path -LiteralPath $RuntimeDir)) {
        Move-Item -LiteralPath $LegacyDir -Destination $RuntimeDir
    }
    New-Item -ItemType Directory -Force -Path $RuntimeDir | Out-Null

    foreach ($StageDir in @(
        "01_QC", "02_Differential", "03_Enrichment", "04_Regulation",
        "05_Network", "06_Custom", "99_Run"
    )) {
        $Container = if ($StageDir -eq "99_Run") { "Records" } else { "Intermediate" }
        $CacheRoot = Join-Path (Join-Path $ProfileDir $StageDir) $Container
        if (-not (Test-Path -LiteralPath $CacheRoot)) { continue }
        Get-ChildItem -LiteralPath $CacheRoot -Recurse -File | ForEach-Object {
            $Relative = $_.FullName.Substring($CacheRoot.Length).TrimStart('\', '/')
            $Target = Join-Path $RuntimeDir $Relative
            New-Item -ItemType Directory -Force -Path (Split-Path -Parent $Target) | Out-Null
            Move-Item -LiteralPath $_.FullName -Destination $Target -Force
        }
    }
    Get-ChildItem -LiteralPath $ProfileDir -Directory |
        Where-Object { $_.Name -ne "_runtime" } |
        Sort-Object FullName -Descending |
        ForEach-Object {
            Get-ChildItem -LiteralPath $_.FullName -Recurse -Directory |
                Sort-Object FullName -Descending |
                Where-Object { -not (Get-ChildItem -LiteralPath $_.FullName -Force) } |
                Remove-Item -Force
            if (-not (Get-ChildItem -LiteralPath $_.FullName -Force)) {
                Remove-Item -LiteralPath $_.FullName -Force
            }
        }
}

function Archive-PipelineWorkspace([string]$ProfileDir, [string]$RuntimeDir) {
    Optimize-PipelineArtifacts $RuntimeDir -ForCache
    if (-not (Test-Path -LiteralPath $RuntimeDir)) { return }

    Get-ChildItem -LiteralPath $RuntimeDir -Recurse -File | ForEach-Object {
        $Relative = $_.FullName.Substring($RuntimeDir.Length).TrimStart('\', '/')
        $Normalized = $Relative.Replace('\', '/')
        $StageDir = switch -Regex ($Normalized) {
            '^data_processed/exprSet_for_heatmap\.Rdata$' { "02_Differential"; break }
            '^data_processed/' { "01_QC"; break }
            '^diff_analysis/' { "02_Differential"; break }
            '^functional_enrichment/custom_gene_sets/' { "06_Custom"; break }
            '^functional_enrichment/ppi_analysis/' { "05_Network"; break }
            '^functional_enrichment/(tf_analysis|gsva_analysis)/' { "04_Regulation"; break }
            '^functional_enrichment/gsea_analysis/tables/GSEA_Full_Table_GTRD_TF\.csv$' {
                "04_Regulation"; break
            }
            '^functional_enrichment/' { "03_Enrichment"; break }
            '^logs/' { "99_Run"; break }
            default { "99_Run" }
        }
        $Container = if ($StageDir -eq "99_Run") { "Records" } else { "Intermediate" }
        $Target = Join-Path (Join-Path (Join-Path $ProfileDir $StageDir) $Container) $Relative
        New-Item -ItemType Directory -Force -Path (Split-Path -Parent $Target) | Out-Null
        Move-Item -LiteralPath $_.FullName -Destination $Target -Force
    }

    Get-ChildItem -LiteralPath $RuntimeDir -Recurse -Directory |
        Sort-Object FullName -Descending |
        Where-Object { -not (Get-ChildItem -LiteralPath $_.FullName -Force) } |
        Remove-Item -Force
    if (-not (Get-ChildItem -LiteralPath $RuntimeDir -Force)) {
        Remove-Item -LiteralPath $RuntimeDir -Force
    }
}

$ConfigPath = Resolve-ProjectPath $Config
if (-not (Test-Path -LiteralPath $ConfigPath)) { throw "Config not found: $ConfigPath" }
$Settings = Import-PowerShellDataFile -LiteralPath $ConfigPath
if (-not $Settings.ProjectName -or $Settings.ProjectName -match '[\\/:*?"<>|]') {
    throw "ProjectName must be a simple Windows folder name."
}

if (-not $Rscript) {
    $Rscript = @(
        "C:\Program Files\R\R-4.4.2\bin\Rscript.exe",
        "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"
    ) | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
}
if (-not $Rscript) { throw "Rscript.exe not found." }

$CountMatrix = Resolve-ProjectPath $Settings.CountMatrix
$Metadata = Resolve-ProjectPath $Settings.Metadata
$WorkRoot = Resolve-ProjectPath $Settings.WorkRoot
$ResultsRoot = Resolve-ProjectPath $Settings.ResultsRoot
$SharedCache = Resolve-ProjectPath $Settings.SharedCache
$AnalysisProfile = if ($Settings.AnalysisProfile) {
    [string]$Settings.AnalysisProfile
} else {
    "Primary_padj0.05_LFC1.0"
}
if ($AnalysisProfile -match '[\\/:*?"<>|]') {
    throw "AnalysisProfile must be a simple Windows folder name."
}
$AnalysisTag = if ([string]::IsNullOrWhiteSpace([string]$Settings.BatchColumn)) {
    "nobatch"
} else {
    "batch_" + (([string]$Settings.BatchColumn) -replace '[^A-Za-z0-9_.-]', '_')
}
$WorkProfileDir = Join-Path (Join-Path $WorkRoot $Settings.ProjectName) $AnalysisProfile
$WorkDir = Join-Path $WorkProfileDir "_runtime"
$ResultDir = Join-Path (Join-Path $ResultsRoot $Settings.ProjectName) $AnalysisProfile
$LogDir = Join-Path $WorkDir "logs"
$LocalRLibrary = Join-Path $ProjectRoot ".r-library"
New-Item -ItemType Directory -Force -Path $LocalRLibrary | Out-Null
$env:R_LIBS_USER = $LocalRLibrary
$env:THEME_CONFIG = Join-Path $ProjectRoot "config/theme.yml"
$env:PLOT_REGISTRY = Join-Path $ProjectRoot "config/plots.yml"

$CheckpointPath = Join-Path $ProjectRoot "config/stage_checkpoint.json"
if ($Stage -match '^\d{2}_') {
    if (-not (Test-Path -LiteralPath $CheckpointPath)) {
        throw @"
Stage checkpoint required before running '$Stage'.
Ask the user whether this stage should run and confirm its stage-specific options.
Then create config/stage_checkpoint.json from config/stage_checkpoint.example.json.
"@
    }
    $Checkpoint = Get-Content -LiteralPath $CheckpointPath -Raw | ConvertFrom-Json
    if (-not $Checkpoint.approved) { throw "Agent checkpoint is not approved." }
    if ([string]$Checkpoint.stage -ne $Stage) {
        throw "Checkpoint is for '$($Checkpoint.stage)', not '$Stage'."
    }
    if ($Stage -eq "02_differential" -and -not $Checkpoint.qc_reviewed) {
        throw "Stage 02 requires qc_reviewed=true after reviewing PCA, batch and outliers."
    }
    Remove-Item -LiteralPath $CheckpointPath -Force

    $env:DEG_PADJ_THRESHOLD = [string]$(if ($null -ne $Checkpoint.padj_threshold) {
        $Checkpoint.padj_threshold
    } else { $Settings.DefaultPadjThreshold })
    $env:DEG_LFC_THRESHOLD = [string]$(if ($null -ne $Checkpoint.lfc_threshold) {
        $Checkpoint.lfc_threshold
    } else { $Settings.DefaultLfcThreshold })
    $env:VOLCANO_LABEL_COUNT = [string]$(if ($null -ne $Checkpoint.volcano_label_count) {
        $Checkpoint.volcano_label_count
    } else { 20 })
    $env:SELECTED_GENES = if ($Checkpoint.selected_genes) {
        ($Checkpoint.selected_genes -join ",")
    } else { "" }
    $env:GSEA_HEATMAP_DATABASE = [string]$(if ($Checkpoint.gsea_database) {
        $Checkpoint.gsea_database
    } else { "Hallmark" })
    $env:GSEA_HEATMAP_PATHWAY = [string]$Checkpoint.gsea_pathway
    $env:ORA_PLOT_TOP_N = [string]$(if ($null -ne $Checkpoint.ora_top_terms) {
        $Checkpoint.ora_top_terms
    } else { 10 })
    $env:PPI_TOP_N = [string]$(if ($null -ne $Checkpoint.ppi_top_genes) {
        $Checkpoint.ppi_top_genes
    } else { 200 })
}

Write-Host ""
Write-Host "Windows bulk RNA-seq pipeline" -ForegroundColor Cyan
Write-Host "Project : $($Settings.ProjectName)"
Write-Host "Profile : $AnalysisProfile"
Write-Host "Stage   : $Stage"
Write-Host "Work    : $WorkProfileDir"
Write-Host "Results : $ResultDir"
Write-Host ""

& $Rscript (Join-Path $ProjectRoot "scripts/support/workflow_utils.R") check-dependencies
if ($LASTEXITCODE -ne 0) { throw "Required R packages are missing." }
if ($Stage -eq "check") { Write-Host "Environment check passed." -ForegroundColor Green; exit 0 }

if (-not (Test-Path -LiteralPath $CountMatrix)) { throw "Count matrix not found: $CountMatrix" }
if (-not (Test-Path -LiteralPath $Metadata)) { throw "Metadata not found: $Metadata" }

if ($Overwrite -and $Stage -eq "01_qc") {
    Remove-SafeDirectory $WorkProfileDir $WorkRoot
    Remove-SafeDirectory $ResultDir $ResultsRoot
}
New-Item -ItemType Directory -Force -Path $WorkProfileDir, $SharedCache | Out-Null
Restore-PipelineWorkspace $WorkProfileDir $WorkDir
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null

if ($Stage -eq "organize") {
    Optimize-PipelineArtifacts $WorkDir
    $padj = if ($Settings.DefaultPadjThreshold) { [double]$Settings.DefaultPadjThreshold } else { 0.05 }
    $lfc = if ($Settings.DefaultLfcThreshold) { [double]$Settings.DefaultLfcThreshold } else { 1.0 }
    & (Join-Path $ProjectRoot "scripts/finalize_results.ps1") `
        -Source $WorkDir -Destination $ResultDir `
        -ProjectName $Settings.ProjectName -AnalysisLabel $AnalysisProfile `
        -ActiveStage $Stage -PadjThreshold $padj -LfcThreshold $lfc -Overwrite
    Archive-PipelineWorkspace $WorkProfileDir $WorkDir
    exit $LASTEXITCODE
}

function Save-Recipe {
    param([hashtable]$Recipe, [string]$Path)
    $Recipe | ConvertTo-Json -Depth 12 | Set-Content -LiteralPath $Path -Encoding UTF8
}

$GitCommit = ""
try {
    $GitCommit = (& git rev-parse HEAD 2>$null).Trim()
} catch {}

$RecipePath = Join-Path $LogDir ("recipe_" + $Stage + ".json")
$Recipe = [ordered]@{
    schema_version = "1.0"
    project = [string]$Settings.ProjectName
    analysis_profile = $AnalysisProfile
    stage = $Stage
    started_at = (Get-Date).ToString("o")
    git_commit = $GitCommit
    inputs = [ordered]@{
        count_matrix = [ordered]@{
            path = $CountMatrix
            sha256 = (Get-FileHash -LiteralPath $CountMatrix -Algorithm SHA256).Hash
        }
        metadata = [ordered]@{
            path = $Metadata
            sha256 = (Get-FileHash -LiteralPath $Metadata -Algorithm SHA256).Hash
        }
    }
    design = [ordered]@{
        control = [string]$Settings.ControlGroup
        treatment = [string]$Settings.TreatGroup
        batch_column = [string]$Settings.BatchColumn
        excluded_samples = [string]$Settings.ExcludeSamples
    }
    checkpoint = if ($null -ne $Checkpoint) { $Checkpoint } else { $null }
    theme_config = $env:THEME_CONFIG
    plot_registry = $env:PLOT_REGISTRY
    steps = @()
}
Save-Recipe -Recipe $Recipe -Path $RecipePath

$env:OUTPUT_DIR = $WorkDir
$env:PROJECT_NAME = [string]$Settings.ProjectName
$env:ANALYSIS_MODE_TAG = $AnalysisTag
$env:BATCH_COLUMN = [string]$Settings.BatchColumn
$env:COUNT_MATRIX_PATH = $CountMatrix
$env:SAMPLE_METADATA_PATH = $Metadata
$env:SHARED_CACHE_DIR = $SharedCache
$env:CONTROL_GROUP = [string]$Settings.ControlGroup
$env:TREAT_GROUP = [string]$Settings.TreatGroup
$env:EXCLUDE_SAMPLES = [string]$Settings.ExcludeSamples
$env:MIN_COUNT = [string]$(if ($null -ne $Settings.MinCount) { $Settings.MinCount } else { 10 })
$env:MIN_SAMPLES = [string]$(if ($null -ne $Settings.MinSamples) { $Settings.MinSamples } else { 3 })
$env:HEATMAP_TOP_PER_DIRECTION = [string]$(if ($null -ne $Settings.HeatmapTopPerDirection) {
    $Settings.HeatmapTopPerDirection
} else {
    20
})

if ($Settings.PrepareCache) {
    & $Rscript (Join-Path $ProjectRoot "scripts/support/prepare_resource_cache.R") *>&1 |
        Set-Content -LiteralPath (Join-Path $LogDir "00_prepare_resource_cache.log")
}

$QcMarker = Join-Path $LogDir "stage_01_qc.complete.json"
$DegMarker = Join-Path $LogDir "stage_02_differential.complete.json"
if ($Stage -eq "02_differential" -and -not (Test-Path -LiteralPath $QcMarker)) {
    throw "Run and review Stage 01 QC first."
}
if ($Stage -in @("03_enrichment", "04_regulation", "05_ppi", "06_custom") -and
    -not (Test-Path -LiteralPath $DegMarker)) {
    throw "Run Stage 02 differential expression first."
}

switch ($Stage) {
    "01_qc" {
        $env:PIPELINE_PHASE = "qc"
        $Scripts = @("scripts/pipeline/01_deseq2_qc_diff.R")
    }
    "02_differential" {
        $env:PIPELINE_PHASE = "differential"
        $Scripts = @(
            "scripts/pipeline/01_deseq2_qc_diff.R",
            "scripts/pipeline/02_prepare_heatmap_matrix.R",
            "scripts/pipeline/03_plot_volcano_ma_heatmap.R",
            "scripts/support/export_deg_workbook.R",
            "scripts/pipeline/12_gene_boxplots.R"
        )
    }
    "03_enrichment" {
        $Scripts = @(
            "scripts/pipeline/04_ora_enrichment.R",
            "scripts/pipeline/04b_ora_sankey_bubble.R",
            "scripts/pipeline/05_gsea_core.R",
            "scripts/pipeline/07_gsea_pathway_heatmap.R"
        )
        if ($Settings.EnablePathview) { $Scripts += "scripts/pipeline/08_kegg_pathview_maps.R" }
        $Scripts += "scripts/support/postprocess_enrichment.R"
    }
    "04_regulation" {
        $Scripts = @(
            "scripts/pipeline/06_gsea_tf_databases.R",
            "scripts/pipeline/10_tf_activity_decoupler.R",
            "scripts/pipeline/11_gsva_tf_correlation.R",
            "scripts/pipeline/13_progeny_decoupler.R"
        )
    }
    "05_ppi" {
        $Scripts = @(
            "scripts/pipeline/09_ppi_string_network.R",
            "scripts/support/postprocess_ppi.R"
        )
    }
    "06_custom" {
        $Scripts = @("scripts/pipeline/14_custom_gene_sets.R")
    }
}

$StatusFile = Join-Path $LogDir ("run_status_" + $Stage + ".tsv")
"step`tscript`tstatus`tseconds`tlog" | Set-Content -LiteralPath $StatusFile

for ($Index = 0; $Index -lt $Scripts.Count; $Index++) {
    $Script = $Scripts[$Index]
    $Name = [IO.Path]::GetFileNameWithoutExtension($Script)
    $Log = Join-Path $LogDir "$Stage`__$Name.log"
    $Stdout = Join-Path $LogDir "$Name.stdout.tmp"
    $Stderr = Join-Path $LogDir "$Name.stderr.tmp"
    $Start = Get-Date

    Write-Host "[$($Index + 1)/$($Scripts.Count)] $Name"
    $Process = Start-Process -FilePath $Rscript `
        -ArgumentList (Join-Path $ProjectRoot $Script) `
        -WorkingDirectory $ProjectRoot `
        -RedirectStandardOutput $Stdout `
        -RedirectStandardError $Stderr `
        -WindowStyle Hidden -Wait -PassThru

    $Lines = @()
    if (Test-Path $Stdout) { $Lines += Get-Content $Stdout }
    if (Test-Path $Stderr) { $Lines += Get-Content $Stderr }
    $Lines | Set-Content -LiteralPath $Log
    Remove-Item $Stdout, $Stderr -Force -ErrorAction SilentlyContinue
    $Elapsed = [math]::Round(((Get-Date) - $Start).TotalSeconds, 1)
    $Status = if ($Process.ExitCode -eq 0) { "SUCCESS" } else { "FAILED" }
    "$($Index + 1)`t$Script`t$Status`t$Elapsed`t$Log" | Add-Content $StatusFile

    $Outputs = Get-ChildItem -LiteralPath $WorkDir -Recurse -File |
        Where-Object { $_.LastWriteTime -ge $Start } |
        ForEach-Object { $_.FullName.Substring($WorkDir.Length).TrimStart('\', '/') }
    $WarningCount = @($Lines | Where-Object { $_ -match '(?i)warning|warn|skipped|璺宠繃' }).Count
    $Recipe.steps += [ordered]@{
        id = $Name
        script = $Script
        status = $Status
        seconds = $Elapsed
        warning_count = $WarningCount
        log = $Log.Substring($WorkDir.Length).TrimStart('\', '/')
        outputs = @($Outputs)
    }
    Save-Recipe -Recipe $Recipe -Path $RecipePath

    if ($Process.ExitCode -ne 0) {
        Write-Host ($Lines | Select-Object -Last 30)
        throw "Step failed: $Name. See $Log"
    }
    Remove-Item -LiteralPath (Join-Path $ProjectRoot "Rplots.pdf") -Force -ErrorAction SilentlyContinue
    Remove-Item -LiteralPath (Join-Path $ProjectRoot "omnipathr-log") -Recurse -Force -ErrorAction SilentlyContinue
}

& $Rscript (Join-Path $ProjectRoot "scripts/support/workflow_utils.R") write-session-info | Out-Null

if ($Stage -in @("02_differential", "03_enrichment", "04_regulation", "05_ppi", "06_custom")) {
    & $Rscript (Join-Path $ProjectRoot "scripts/support/generate_summary_report.R") *>&1 |
        Set-Content -LiteralPath (Join-Path $LogDir "generate_summary_report.log")
}

$Recipe.completed_at = (Get-Date).ToString("o")
$Recipe.result_directory = $ResultDir
Save-Recipe -Recipe $Recipe -Path $RecipePath

$StageMarker = Join-Path $LogDir ("stage_" + ($Stage -replace '_.*$', '') + "_" +
    ($Stage -replace '^\d{2}_', '') + ".complete.json")
[ordered]@{
    stage = $Stage
    completed_at = (Get-Date).ToString("o")
    checkpoint = $Checkpoint
} | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $StageMarker -Encoding UTF8

Optimize-PipelineArtifacts $WorkDir

$FinalPadj = [double]$env:DEG_PADJ_THRESHOLD
$FinalLfc = [double]$env:DEG_LFC_THRESHOLD
& (Join-Path $ProjectRoot "scripts/finalize_results.ps1") `
    -Source $WorkDir -Destination $ResultDir `
    -ProjectName $Settings.ProjectName -AnalysisLabel $AnalysisProfile `
    -ActiveStage $Stage -PadjThreshold $FinalPadj -LfcThreshold $FinalLfc -Overwrite
Archive-PipelineWorkspace $WorkProfileDir $WorkDir

Write-Host ""
Write-Host "Pipeline completed." -ForegroundColor Green
Write-Host "Browse results: $ResultDir"

