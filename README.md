# Reusable bulk RNA-seq workflow template

Windows 原生 bulk RNA-seq 分析项目。输入为基因层面的原始整数 count matrix。

## 目录

```text
INPUT/                  用户输入：count matrix、sample metadata
resources/              自定义基因集、TF 数据库、MSigDB 文件和网络缓存
config/                 项目配置、图形配置和一次性运行参数
scripts/
  pipeline/             六阶段调用的标准分析脚本
  support/              共享函数、报告、Excel 和后处理
  finalize_results.ps1  按当前阶段整理 results，不复制 work
work/<project>/<profile>/
results/<project>/<profile>/
```

`work/<project>/<profile>/` 与 results 使用相同的阶段目录：`01_QC`、
`02_Differential`、`03_Enrichment`、`04_Regulation`、`05_Network`、
`06_Custom` 和 `99_Run`。各阶段只在 `Intermediate` 中保存可复用缓存；
不会复制最终图片。`results` 保存 Figures、Tables、Interactive 和 README。

## 分析配置

- 正式分析：`Primary_padj0.05_LFC1.0`
- 探索性分析：`Exploratory_padj0.05_LFC0.5`
- 基因预过滤：至少 3 个样本中 count ≥ 10
- DESeq2 设计：`~ group`（有批次列时为 `~ batch + group`）
- fold change：ashr shrinkage
- 正式 DEG：padj < 0.05 且 |ashr log2FC| > 1
- 探索性 DEG：padj < 0.05 且 |ashr log2FC| > 0.5
- GSEA：全部通过预过滤的基因，按 DESeq2 Wald statistic 排序
- ORA：进入检验的基因作为背景，Up/Down 分开分析

## 已执行分析

样本 QC、PCA、相关性、距离和聚类；DESeq2 差异表达；火山图、MA 图、热图和
单基因箱线图；GO/KEGG/Hallmark ORA；GO/KEGG/Hallmark/Reactome GSEA；
TF-target GSEA；CollecTRI TF activity；GSVA 与 TF–pathway correlation；
PROGENy；STRING PPI、hub 和模块；自定义基因集热图及 GSEA；KEGG Pathview。

## 六阶段运行

安装依赖：

```powershell
& "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" scripts\setup_packages.R
```

正式分析使用 `config/project.psd1`，探索性分析使用
`config/project.exploratory.psd1`。禁止默认 full 连跑。每个阶段开始前询问用户，
并从示例创建只针对当前阶段的一次性 `config/stage_checkpoint.json`。

```powershell
.\run.ps1 -Stage 01_qc -Config config/project.psd1 -Overwrite
.\run.ps1 -Stage 02_differential -Config config/project.psd1
.\run.ps1 -Stage 03_enrichment -Config config/project.psd1
.\run.ps1 -Stage 04_regulation -Config config/project.psd1
.\run.ps1 -Stage 05_ppi -Config config/project.psd1
.\run.ps1 -Stage 06_custom -Config config/project.psd1
```

01 完成后必须停止并检查 PCA、批次和离群样本，再修改 `BatchColumn`、
`ExcludeSamples` 并批准 02。每个后续阶段也可选择跳过。

差异表达阶段结束时会立即生成：

```text
02_Differential/Tables/00_<Project>_Differential_Expression.xlsx
```

Excel 只保留在 `02_Differential/Tables/`，避免重复副本。

## 个性化分析

特殊绘图和项目专属分析放入 `personalized/`，使用 `run_personalized.ps1`。
其代码、work 和 results 与标准 Step 01–06 完全隔离，不进入标准 manifest。

## 网络依赖

- STRING PPI：本地缓存缺失时访问 STRING。
- CollecTRI/PROGENy：缓存缺失时访问 OmniPath/decoupleR 数据源。
- KEGG Pathview：未缓存通路需要访问 KEGG。
- R 包安装、数据库首次安装或升级需要网络。
- DESeq2、QC、已有缓存的调控分析、本地 TF 数据库和自定义基因集可离线运行。

## 解释边界

TF activity、PROGENy、GSVA correlation、PPI hub/module 和数据库推断的通路方向
属于探索性结果，需要独立实验验证。本项目从 count matrix 开始，不能替代
FASTQ 层面的 FastQC/MultiQC、接头、比对率、链特异性和污染检查。
