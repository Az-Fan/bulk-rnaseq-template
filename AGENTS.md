# Bulk RNA-seq v4：Agent 强制执行合同

所有 Agent 在本仓库行动前必须完整阅读本文件。这里是“一仓库、多项目、固定结构、
单 Pixi、结果易查”的通用分析平台，不是一次性项目目录。

## 1. 不可违反的规则

1. 不修改、移动、覆盖或删除用户原始数据、v3、历史项目或已发布结果。
2. 不引入 renv、BiocManager::install、install.packages、remotes::install_*、系统
   Rscript 或第二套下游包管理流程；用户只需 `pixi install --locked --all`。
3. 正式运行不联网、不更新数据库、不安装依赖。`resources-sync` 是唯一资源联网阶段，
   且必须得到明确批准。
4. 不静默跳过模块，不静默猜测设计、物种、comparison、reference level 或 counts
   来源；未确认必须停止。
5. 每一次分析前都必须重新向用户展示模块清单并获得本次确认。旧配置中的确认不能
   代替本次询问。
6. 不因精简目录、报告或图件数量而删除分析能力、完整数值结果或 provenance。
7. 不用 TPM/FPKM/CPM、批次校正值或其他非 raw counts 做 DESeq2 检验。
8. motif、TF activity、PROGENy、PPI hub/module 只可解释为关联或探索证据，不写成
   直接结合或因果机制。
9. 原始输入只读；导入时复制并保存来源、字节数和 SHA256。
10. 不手工编辑正式 `results/` 或 `work/history/published/`。
11. 新图必须实际渲染检查；不能只确认代码成功。
12. 统计方法变更必须补单元测试、通用集成测试，并在用户提供历史基线时做文件级和
    数值回归。
13. 删除缓存、运行或资源前必须展示解析后的绝对目标、预计空间和确认令牌；默认 dry-run。
14. 个性化代码只能读取标准产物，不得复制或改写核心统计链。
15. 未经用户明确批准，不得改变本文件定义的仓库或项目结构。

## 2. 固定目录树与职责

```text
bulk-rnaseq-v4/
├── README.md                         # 生信使用者入口
├── AGENTS.md                         # 本合同
├── pixi.toml                         # 唯一公开命令和软件声明
├── pixi.lock                         # 唯一锁文件
├── workflow/                         # 用户可审查的科学方法层
│   ├── run.R                         # 只编排，不隐藏统计决策
│   ├── functions.R                   # 共享科学、I/O、主题和绘图函数
│   ├── 01_qc.R                       # counts/sample QC、PCA、相关性、聚类
│   ├── 02_differential.R             # DESeq2、contrast、LFC 收缩、DEG
│   ├── 03_enrichment.R               # ORA、GSEA、Pathview
│   ├── 04_activity.R                 # GSVA、TF、PROGENy、相关性
│   ├── 05_network.R                  # PPI、module、hub、module enrichment
│   ├── 06_motif.R                    # promoter motif / peak-aware motif
│   └── 07_exploratory.R              # 指定基因、自定义基因集、适用探索
├── projects/
│   ├── _template/                    # 唯一空白模板；故意未确认
│   └── <project_id>/
│       ├── project.yml               # 来源、设计、模块、格式和资源
│       ├── samples.tsv               # 样本、协变量、排除及理由
│       ├── contrasts.tsv             # numerator 相对 denominator 的方向
│       ├── qc_approval.yml            # 用户 QC/排除确认
│       ├── input/                    # counts 校验副本与上游配置
│       ├── results/                  # 唯一当前正式结果
│       │   └── index.html            # 唯一浏览总入口
│       └── work/                     # staging、cache、logs、history
├── resources/
│   ├── registry.yml                  # release/source/hash/license 元数据
│   ├── sources/                      # 唯一允许联网的同步声明
│   └── data/                         # 本地冻结资源，不提交 Git
├── tests/                            # 通用结构、禁令、数值和渲染测试
└── internal/                         # 不放新的科学决策
    ├── cli/                          # init/doctor/check/analyze/clean
    ├── lib/                          # schema、哈希、manifest、控制逻辑
    ├── schemas/                      # 项目和资源 schema
    ├── reporting/                    # HTML、发布和视觉 QA
    ├── audit/                        # 可选历史结果盘点/比较工具
    ├── motif/                        # MEME/STREME 执行包装
    ├── upstream/
    │   ├── adapter.py                # 显式上游服务器适配器
    │   └── omics-pipelines/          # 固定提交 Git submodule
    ├── vendor/r-apear/               # 固定源码、Pixi 构建依赖
    ├── utilities/                    # 可审计维护工具
    └── windows/                      # Windows→WSL 导入入口
```

仓库根目录不得新增 `runs/`、`analysis/`、`scripts/`、`R/` 等平级实现。所有项目
内容只能进入 `projects/<project_id>/`；共享资源只能进入 `resources/`。

## 3. 收到 counts 后的强制流程

### 3.1 先问，不运行

先确认项目 ID；同名项目不得覆盖。然后逐项询问并记录：

- counts 来自 featureCounts、STAR、Salmon/tximport、其他还是未知；
- gene-level 还是 transcript-level；
- 物种、assembly、annotation release 和 gene ID 类型；
- 是否做过 TPM/FPKM/CPM/RPKM、批次校正或其他归一化；
- 技术重复、上游样本/基因过滤、strandedness；
- 样本到组/批次/配对/协变量的映射；
- 设计公式、参考水平、每个 comparison 的 numerator/denominator；
- 是否有 FASTQ 等上游数据，是否本次需要上游；
- 下面列出的每一个分析模块是否启用、非适用或由用户跳过。

`unknown` 可以记录，但必须在报告中暴露限制。明确为标准化表达值时禁止 DESeq2；
transcript counts 必须先确认 tximport/聚合策略。

### 3.2 只读检查和导入

检查格式、分隔符、行列数、样本列、ID、重复 ID、NA/Inf、负值、小数、library
size 和零值比例。整数不等于已证明的 raw counts。确认后用 `import-counts` 复制到
项目 `input/` 并校验哈希；不得修改源文件。

### 3.3 设计与 QC 双重确认

生成 `project.yml`、`samples.tsv`、`contrasts.tsv` 草案。任何模块或上游决定为
`confirmed: false` 时必须停止。先做预分析 QC；样本排除只能由用户确认，并写入
`samples.tsv` 与 `qc_approval.yml`，不得因 PCA 离群而自动删除。

### 3.4 每次运行重新确认模块

必须把 11 个细分决定逐项展示给用户：

- QC；差异表达；ORA/GSEA 富集；GSVA/TF/PROGENy 调控；PPI 网络；
- promoter motif；peak-aware motif；自定义基因集；Pathview；个性化分析；WGCNA。

得到本次确认后运行 `pixi run modules -- --project ...`，并把其完整
`confirmation_token` 原样交给 `analyze --confirm-modules`。交互/非交互、Python/R
入口都不得绕过这一条件。

### 3.5 staging、QA、发布

正式数值分析只计算一次：一个 DESeq2 共享模型供所有 contrasts 使用；GSEA 使用
完整 Wald 排序；报告详细度和导出格式不得触发统计重算。结果先进入 staging。
启用模块失败必须 `failed_explicit`，非适用或跳过必须有理由。

发布前必须通过：配置/资源/输入哈希、单元与集成测试、PDF 渲染、HTML 链接、模块
状态和 manifest。若用户提供旧结果，还要生成 `legacy_analysis_matrix.tsv`、
`legacy_output_inventory.tsv`；`missing` 必须阻止发布。旧正式结果移入 history，
不得覆盖或删除。

## 4. 模块、报告与导出是独立维度

```yaml
analysis:
  modules: {}       # 是否计算
report:
  detail: comprehensive
export:
  figures: all      # 发布图件范围
  formats: [pdf]    # 文件格式
```

报告精简不代表分析被跳过。默认只生成 PDF；只有用户明确说需要 PNG 时，才允许加入
`png`。不得默认生成 SVG。火山图固定为用户选定的经典单面板样式；风格调整不能改变
点集合、阈值、方向、标签选择规则或统计值。

## 5. 可选上游的硬边界

每个新项目必须先问“是否有上游 FASTQ，以及本次是否需要上游分析”。默认模板的
`upstream.status: unconfirmed` 会阻止分析，不能静默当作 disabled。

明确不需要时记录 `status: disabled, confirmed: true`。明确需要时才可改为 enabled，
准备项目内的 `input/upstream/` 配置，并再次运行 `upstream-plan` 展示 provider、固定
提交、执行主机、目标和确认令牌。

- 上游来源固定为 `xuzhougeng/omics-pipelines` 提交
  `ce4e2ec88da6663a32b7099c5850e1a51ad66952`。
- 仅调用 FASTQ QC、fastp、STAR、BAM/coverage 依赖和 raw gene-count merge。
- 目标必须是 `results/04-quant/counts.tsv`；不得调用该上游的 DESeq2 或 enrichment。
- 默认 `execution_host: server`、`allow_wsl: false`。本机 WSL 不建 index、不跑 STAR。
- STAR index 必须由服务器管理员预先提供；适配器不创建、不下载、不猜测版本。
- 上游输出仍须通过标准 counts 导入、来源确认、QC 和下游模块确认。
- submodule 是第三方边界，不得直接改其源码。更新提交必须用户批准、重新审计并更新测试。

## 6. 软件、资源与计算

软件由一个根 Pixi workspace 和一个锁文件管理。锁内可以有兼容的 default/motif
环境，但用户没有第二套安装流程。大型数据库、FASTA/GTF、MSigDB、JASPAR、STRING、
OmniPath/CollecTRI/PROGENy、KEGG/Pathview、ChEA/ENCODE/GTRD 和自定义基因集属于
资源层，每个资源必须记录 ID、物种、assembly、release、来源、时间、SHA256 和许可。

普通任务最多 3 worker；DESeq2 共享模型单任务；BLAS/OMP/data.table 每 worker 1
线程。运行前显示 CPU、内存与磁盘估计。大表优先压缩 TSV/Parquet。

## 7. 结果与科学解释

结果按 7 个清晰模块组织，共享 QC 和模型只出现一次，多 comparison 结果进入
`Comparisons/<contrast_id>/`。表放 Tables，图放 Figures；没有显著结果仍保存完整表并
明确说明，不生成空图冒充完成。

ORA 依赖阈值及明确 universe；GSEA 使用完整排序，两者不可互相替代。WGCNA、SVA、
细胞比例等仅在适用条件满足且用户确认后运行。图形不得隐藏异常值、调换上下调方向、
改变坐标语义或重算既有布局。

## 8. Agent 交付时必须报告

- 唯一正式结果/或 staging 路径；
- 本次用户确认的模块清单和上游决定；
- 启用、非适用、失败、用户跳过状态；
- tests、PDF 实际渲染、资源与哈希 QA；
- 主要 Figures/Tables/HTML 入口和已知限制；
- 发生的任何清理及其可恢复性。
