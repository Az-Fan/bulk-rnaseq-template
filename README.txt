REUSABLE BULK RNA-SEQ WORKFLOW TEMPLATE
=======================================

输入
----
INPUT 只保存用户输入：
1. 原始整数 count matrix
2. sample_metadata.csv

资源
----
resources/gene_sets       自定义基因集
resources/tf_databases    ENCODE、ChEA 等本地 TF 数据库
resources/msigdb          本地 GMT 文件
resources/cache           CollecTRI、PROGENy、STRING、KEGG 缓存

分析配置
--------
正式分析：Primary_padj0.05_LFC1.0
探索分析：Exploratory_padj0.05_LFC0.5

基因预过滤：至少 3 个样本中 count >= 10
正式 DEG：padj < 0.05 且 |ashr-shrunken log2FC| > 1
探索 DEG：padj < 0.05 且 |ashr-shrunken log2FC| > 0.5
GSEA：全部过滤后基因按 DESeq2 Wald statistic 排序，不使用 DEG 子集
ORA：进入检验的基因作为背景，Up 和 Down 分开分析

完成的分析
----------
样本 QC、PCA、相关性、距离、聚类、测序深度
DESeq2 差异表达、ashr 收缩、火山图、MA 图、热图、单基因箱线图
GO、KEGG、Hallmark ORA
GO、KEGG、Hallmark、Reactome GSEA
ENCODE、ChEA、GTRD TF-target GSEA
CollecTRI TF 活性
GSVA 和 TF-pathway correlation
PROGENy 通路活性
STRING PPI、hub、Louvain 模块和 Cytoscape 表
自定义基因集热图和 GSEA
KEGG Pathview

六阶段执行
----------
01_QC：仅做 QC，完成后停止；人工判断批次和离群样本。
02_Differential：按最终 batch/outlier 决定运行 DESeq2 并生成 Excel。
03_Enrichment：确认 DEG 阈值后运行 ORA/GSEA。
04_Regulation：询问是否进行 TF、GSVA 和 PROGENy。
05_PPI：询问是否进行 STRING PPI。
06_Custom：询问是否进行自定义基因集分析。

每个阶段开始前都必须创建只针对该阶段的 config/stage_checkpoint.json。
禁止默认 full 连跑，允许跳过 03–06 中不需要的阶段。

结果结构
--------
results/<项目>/<分析配置>/ 使用编号目录：
01_QC
02_Differential
03_Enrichment
04_Regulation
05_Network
06_Custom
99_Run

Figures     图片
Tables      CSV、Excel 等表格
Interactive HTML 交互图
work/<项目>/<分析配置>/ 与 results 使用相同的 01_QC、02_Differential、
03_Enrichment、04_Regulation、05_Network、06_Custom、99_Run 目录。
各阶段仅在 Intermediate 中保存下次绘图或续跑需要的缓存，不复制最终图片。

每个编号目录均有 README.txt，建议先阅读后查看结果。

差异表达 Excel
--------------
差异表达阶段结束时立即生成：
02_Differential/Tables/00_<Project>_Differential_Expression.xlsx
Excel 只保留在 02_Differential/Tables，避免重复副本。

网络风险
--------
高风险：首次安装 R 包、KEGG Pathview、无缓存的 STRING、CollecTRI、PROGENy。
本地缓存存在时，上述分析可复用缓存。DESeq2、QC、本地 TF 数据库和自定义基因集
不依赖网络。

个性化分析
----------
特殊绘图或项目专属分析放在 personalized/code、personalized/work 和
personalized/results。它们不进入标准结果，也不进入标准 manifest。

解释边界
--------
TF、PROGENy、GSVA、PPI 和数据库通路推断属于探索性结果，不能单独作为机制结论。
本项目从 count matrix 开始，不能替代 FASTQ 层面的质量控制。
