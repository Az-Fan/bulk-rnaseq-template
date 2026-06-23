param(
    [Parameter(Mandatory = $true)][string]$Source,
    [Parameter(Mandatory = $true)][string]$Destination,
    [string]$ProjectName = "Project",
    [string]$AnalysisLabel = "Analysis",
    [string]$ActiveStage = "",
    [double]$PadjThreshold = 0.05,
    [double]$LfcThreshold = 1.0,
    [switch]$Overwrite
)

$ErrorActionPreference = "Stop"
$Source = [IO.Path]::GetFullPath($Source)
$Destination = [IO.Path]::GetFullPath($Destination)
if (-not (Test-Path -LiteralPath $Source)) { throw "Pipeline output not found: $Source" }

if ((Test-Path -LiteralPath $Destination) -and $Overwrite) {
    $stageCategory = @{
        "01_qc" = "01_QC"
        "02_differential" = "02_Differential"
        "03_enrichment" = "03_Enrichment"
        "04_regulation" = "04_Regulation"
        "05_ppi" = "05_Network"
        "06_custom" = "06_Custom"
    }[$ActiveStage]
    if ($stageCategory) {
        $target = Join-Path $Destination $stageCategory
        if (Test-Path -LiteralPath $target) {
            Remove-Item -LiteralPath $target -Recurse -Force
        }
    }
}
New-Item -ItemType Directory -Force -Path $Destination | Out-Null

function Test-ActiveStage([string]$StageName) {
    return [string]::IsNullOrWhiteSpace($ActiveStage) -or
        $ActiveStage -eq "organize" -or
        $ActiveStage -eq $StageName
}

function Get-Purpose([IO.FileInfo]$File) {
    switch ($File.Extension.ToLowerInvariant()) {
        ".png"  { "Figures" }
        ".jpg"  { "Figures" }
        ".jpeg" { "Figures" }
        ".svg"  { "Figures" }
        ".pdf"  { "Figures" }
        ".csv"  { "Tables" }
        ".tsv"  { "Tables" }
        ".xlsx" { "Tables" }
        ".xls"  { "Tables" }
        ".html" { "Interactive" }
        ".md"   { "Methods" }
        ".txt"  { "Tables" }
        ".rdata" { "Intermediate" }
        ".rds"   { "Intermediate" }
        default { "Records" }
    }
}

function Copy-ArtifactGroup(
    [string]$RelativeRoot,
    [string]$Target,
    [string[]]$Patterns = @("*")
) {
    $root = Join-Path $Source $RelativeRoot
    if (-not (Test-Path -LiteralPath $root)) { return }
    Get-ChildItem -LiteralPath $root -Recurse -File | Where-Object {
        $name = $_.Name
        $matches = ($Patterns | Where-Object { $name -like $_ }).Count -gt 0
        $hasPng = $_.Extension -ieq ".pdf" -and
            (Test-Path -LiteralPath ([IO.Path]::ChangeExtension($_.FullName, ".png")))
        $matches -and -not $hasPng
    } | ForEach-Object {
        $purpose = Get-Purpose $_
        if ($purpose -eq "Methods") { return }
        if ($purpose -eq "Intermediate") { return }
        if ($purpose -eq "Records" -and -not $Target.StartsWith("99_Run")) { return }
        $relative = $_.FullName.Substring($root.Length).TrimStart('\', '/')
        $relativeDir = Split-Path -Parent $relative
        $targetBase = Join-Path $Destination $Target
        $targetLeaf = (($Target -split '[\\/]') | Select-Object -Last 1)
        if ($targetLeaf -eq $purpose) {
            $targetDir = $targetBase
        } else {
            $targetDir = Join-Path $targetBase $purpose
        }
        if ($relativeDir) { $targetDir = Join-Path $targetDir $relativeDir }
        New-Item -ItemType Directory -Force -Path $targetDir | Out-Null
        Copy-Item -LiteralPath $_.FullName -Destination (Join-Path $targetDir $_.Name) -Force
    }
}

if (Test-ActiveStage "01_qc") {
    Copy-ArtifactGroup "plots/qc_plots" "01_QC"
    Copy-ArtifactGroup "data_processed" "01_QC" @(
        "Sample_*", "metadata.*", "Gene_Filtering_Summary.csv",
        "Analysis_Config_Used.csv", "Analysis_Limitations.txt", "Excluded_Samples.txt",
        "QC_REVIEW_REQUIRED.txt"
    )
}

# 02 Differential expression
if (Test-ActiveStage "02_differential") {
    Copy-ArtifactGroup "diff_analysis" "02_Differential"
    Copy-ArtifactGroup "plots/volcano_plots" "02_Differential"
    Copy-ArtifactGroup "plots/heatmaps" "02_Differential"
    Copy-ArtifactGroup "plots/individual_gene_boxplots" "02_Differential"
}

# 03 Enrichment, grouped by database/analysis type
if (Test-ActiveStage "03_enrichment") {
    Copy-ArtifactGroup "functional_enrichment/go_analysis" "03_Enrichment/ORA/GO"
    Copy-ArtifactGroup "functional_enrichment/kegg_analysis" "03_Enrichment/ORA/KEGG"
    Copy-ArtifactGroup "functional_enrichment/hallmark_analysis" "03_Enrichment/ORA/Hallmark"
    Copy-ArtifactGroup "functional_enrichment/advanced_plots" "03_Enrichment/ORA/Advanced"
    Copy-ArtifactGroup "functional_enrichment/apear_network" "03_Enrichment/ORA/Networks"
    Copy-ArtifactGroup "functional_enrichment/sankey_bubble" "03_Enrichment/ORA/Sankey"
    Copy-ArtifactGroup "functional_enrichment/gsea_analysis/GO_BP_plots" "03_Enrichment/GSEA/GO_BP"
    Copy-ArtifactGroup "functional_enrichment/gsea_analysis/KEGG_plots" "03_Enrichment/GSEA/KEGG"
    Copy-ArtifactGroup "functional_enrichment/gsea_analysis/Hallmark_plots" "03_Enrichment/GSEA/Hallmark"
    Copy-ArtifactGroup "functional_enrichment/gsea_analysis/Reactome_plots" "03_Enrichment/GSEA/Reactome"
    Copy-ArtifactGroup "functional_enrichment/gsea_analysis/pathway_heatmaps" "03_Enrichment/Selected_Pathways"
    Copy-ArtifactGroup "functional_enrichment/gsea_analysis/tables" "03_Enrichment/GSEA/Tables" @(
        "GSEA_Full_Table_GO_BP.csv", "GSEA_Full_Table_KEGG.csv",
        "GSEA_Full_Table_Hallmark.csv", "GSEA_Full_Table_Reactome.csv",
        "GSEA_ranked_gene_list.csv", "GSEA_Gene_Set_Provenance.csv"
    )
    Copy-ArtifactGroup "functional_enrichment" "03_Enrichment/Overview" @(
        "ORA_universe_entrez.csv", "ORA_Gene_Set_Provenance.csv"
    )
    Copy-ArtifactGroup "functional_enrichment/gsea_analysis" "03_Enrichment/Overview" @("GSEA_Summary_Dotplots.*")
}

# 04 Regulation
if (Test-ActiveStage "04_regulation") {
    Copy-ArtifactGroup "functional_enrichment/tf_analysis/boxplots" "04_Regulation/boxplots"
    Copy-ArtifactGroup "functional_enrichment/tf_analysis/decoupler_inference" "04_Regulation/TF_Activity"
    Copy-ArtifactGroup "functional_enrichment/tf_analysis/ENCODE_2015" "04_Regulation/TF_Target_GSEA/ENCODE"
    Copy-ArtifactGroup "functional_enrichment/tf_analysis/ChEA_Consensus" "04_Regulation/TF_Target_GSEA/ChEA"
    Copy-ArtifactGroup "functional_enrichment/gsea_analysis/GTRD_TF_plots" "04_Regulation/TF_Target_GSEA/GTRD"
    Copy-ArtifactGroup "functional_enrichment/gsea_analysis/tables" "04_Regulation/TF_Target_GSEA/GTRD" @(
        "GSEA_Full_Table_GTRD_TF.csv"
    )
    Copy-ArtifactGroup "functional_enrichment/gsva_analysis" "04_Regulation/GSVA"
    Copy-ArtifactGroup "functional_enrichment/tf_analysis/decouple_progeny" "04_Regulation/PROGENy"
}

# 05 Network
if (Test-ActiveStage "05_ppi") {
    Copy-ArtifactGroup "functional_enrichment/ppi_analysis" "05_Network/PPI"
}

# 06 Custom gene sets
if (Test-ActiveStage "06_custom") {
    Copy-ArtifactGroup "functional_enrichment/custom_gene_sets" "06_Custom"
}

# 99 Run records
Copy-ArtifactGroup "logs" "99_Run/Logs" @(
    "run_status_*.tsv", "recipe_*.json", "stage_*.complete.json",
    "session_info.txt", "*__*.log"
)
Copy-ArtifactGroup "summary" "99_Run/Reports" @("*.html")

function Write-Readme([string]$RelativeDir, [string[]]$Lines) {
    $dir = Join-Path $Destination $RelativeDir
    New-Item -ItemType Directory -Force -Path $dir | Out-Null
    $Lines | Set-Content -LiteralPath (Join-Path $dir "README.txt") -Encoding UTF8
}
function Write-ReadmeIfExists([string]$RelativeDir, [string[]]$Lines) {
    if (Test-Path -LiteralPath (Join-Path $Destination $RelativeDir)) {
        Write-Readme $RelativeDir $Lines
    }
}

$rule = "padj < $PadjThreshold and |ashr-shrunken log2FC| > $LfcThreshold"
Write-Readme "" @(
    "$ProjectName - $AnalysisLabel",
    "",
    "差异表达阈值：$rule",
    "基因预过滤：至少 3 个样本中 count >= 10。",
    "GSEA 使用全部过滤后基因，并按 DESeq2 Wald statistic 排序。",
    "ORA 使用进入检验的基因作为背景，Up/Down 分开分析。",
    "",
    "查看结果前请先阅读每个编号目录中的 README.txt。",
    "Figures 为图片；Tables 为数据表；Interactive 为 HTML；",
    "Intermediate 仅存在于 work。每个阶段只创建实际有内容的目录。"
)
Write-ReadmeIfExists "01_QC" @(
    "01_QC", "", "样本层面 count QC、PCA、相关性、距离、聚类和文库大小。",
    "建议先查看 PCA 和 Sample_QC_Metrics.csv。CHECK 表示需要人工复核，不代表自动排除。"
)
Write-ReadmeIfExists "02_Differential" @(
    "02_Differential", "", "DESeq2 差异表达模型，使用 ashr 收缩 log2FC。",
    "当前 DEG 阈值：$rule",
    "差异表达阶段完成时立即生成 Excel 工作簿。",
    "建议先查看 Tables/00_${ProjectName}_Differential_Expression.xlsx。"
)
Write-ReadmeIfExists "03_Enrichment" @(
    "03_Enrichment", "", "GO、KEGG、Hallmark 和 Reactome 的 ORA 与 preranked GSEA。",
    "ORA 使用的 DEG 阈值：$rule",
    "GSEA 不使用 DEG 子集，而是按 Wald statistic 排序全部检验基因。",
    "Pathview 及部分数据库在没有缓存时可能需要网络。"
)
Write-ReadmeIfExists "04_Regulation" @(
    "04_Regulation", "", "TF-target GSEA、CollecTRI TF 活性、GSVA 相关性和 PROGENy。",
    "boxplots 中的 Top20 TF 表达图：优先展示 GSEA 显著且 TF 活性排名靠前的 TF；不足 20 个时用 TF 活性 Top 排名补齐。",
    "基因名颜色：绿色=TF 自身表达方向与预测活性一致；黄色=表达变化弱；红色=方向相反。",
    "这些属于探索性调控推断，不能单独建立机制结论。"
)
Write-ReadmeIfExists "05_Network" @(
    "05_Network", "", "STRING PPI、hub 排名、模块和 Cytoscape 输入表。",
    "网络中心性属于探索性结果；本地缓存缺失时 STRING 需要网络。"
)
Write-ReadmeIfExists "06_Custom" @(
    "06_Custom", "", "基于 resources/gene_sets/custom_gene_sets.xlsx 的热图和 GSEA。"
)
Write-Readme "99_Run" @(
    "99_Run", "", "日志、recipe、session information、报告和可复现记录。",
    "使用 run_status_<stage>.tsv 和 recipe_<stage>.json 审计各阶段。"
)

$manifest = Get-ChildItem -LiteralPath $Destination -Recurse -File | ForEach-Object {
    [pscustomobject]@{
        category = $_.DirectoryName.Substring($Destination.Length).TrimStart('\', '/')
        file = $_.Name
        bytes = $_.Length
    }
}
$manifest | Export-Csv -LiteralPath (Join-Path $Destination "manifest.csv") -NoTypeInformation -Encoding UTF8
Write-Host "Results view updated: $Destination"
