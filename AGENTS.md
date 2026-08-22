# Bulk RNA-seq v4：Agent 强制执行合同

所有 Agent 在本仓库中行动前，必须完整阅读本文件。本项目是“一仓库、多项目、固定结构、单 Pixi、正式运行离线”的通用 Bulk RNA-seq 平台，不是一次性分析目录。

## 1. 不可违反的规则

1. 不修改、移动、覆盖或删除原始数据、历史项目、v3、已发布结果及用户提供的旧结果。
2. 不引入 renv、BiocManager、`install.packages()`、`remotes::install_*`、系统 R 或第二套包管理流程。用户只执行 `pixi install --locked --all`。不得删除或忽略仓库跟踪的 `.pixi/config.toml`；它是 Bioconda annotation data 包可重建所必需的安装策略。
3. 正式分析不联网、不安装依赖、不更新数据库。`resources-sync` 是唯一资源联网阶段，且必须单独获批。
4. 不猜测 counts 来源、物种、assembly、设计、参考水平、contrast、模块或分析参数；未确认即停止。
5. 每次正式运行都必须重新展示完整分析计划并取得本次确认。项目中过去的 `confirmed: true` 不能代替当前对话中的确认。
6. 不静默跳过启用模块；不完整或未实现的模块必须阻止运行，不能标记为完成。
7. 不因精简目录、报告或图件数量而删除分析能力、完整数值表或 provenance。
8. 不使用 TPM、FPKM、CPM、批次校正值或其他非 raw counts 做 DESeq2 检验。
9. 技术重复不得当作生物学重复；必须先明确其合并策略。
10. 原始输入只读；导入时复制并保存来源、字节数和 SHA256。
11. 不手工编辑正式 `results/` 或 `work/history/published/`。
12. motif、TF/PROGENy activity、PPI hub/module 只能解释为关联或探索证据，不能写成直接结合或因果机制。
13. 新图必须实际渲染检查；不能只看代码是否成功。
14. 统计方法变化必须补单元测试、通用集成测试；存在历史基线时还必须做文件级和数值回归。
15. 删除缓存、运行或资源前必须先 dry-run，展示解析后的绝对路径、文件数、预计空间和确认令牌。
16. 个性化代码只能读取标准产物，不能复制或改写核心统计链。
17. 未经用户明确批准，不改变本文件定义的仓库和项目结构。

## 2. 固定目录树

```text
bulk-rnaseq-v4/
├── README.md                         # 生信使用者唯一入口
├── AGENTS.md                         # Agent 执行合同
├── pixi.toml / pixi.lock             # 唯一软件环境和锁文件
├── .pixi/config.toml                  # 锁定环境的 post-link 策略；其余 .pixi 内容不提交
├── workflow/                         # 可审查的科学方法层
│   ├── Snakefile                     # 正式 DAG：依赖、资源和断点续跑
│   ├── functions.R                   # 共享验证、I/O、图形与导出函数
│   ├── 01_qc.R                       # counts/sample QC、PCA、相关性、聚类
│   ├── 02_differential.R             # DESeq2、contrast、LFC 收缩、DEG
│   ├── 03_enrichment.R               # ORA、GSEA
│   ├── 04_activity.R                 # GSVA、TF、PROGENy
│   ├── 05_network.R                  # PPI、module、hub、module enrichment
│   ├── 06_motif.R                    # promoter motif；peak motif 仍受硬门禁
│   ├── 07_exploratory.R              # boxplot、自定义基因集、Pathview 等
│   └── snakemake/                    # 单模块执行与最终封存
├── projects/
│   ├── _template/                    # 唯一空白模板，故意不可直接运行
│   └── <project_id>/
│       ├── project.yml               # 来源、设计、参数、模块、资源、导出
│       ├── samples.tsv               # 样本、组、批次、配对、排除及理由
│       ├── contrasts.tsv             # numerator 相对 denominator 的方向
│       ├── qc_approval.yml            # 用户 QC/排除审批
│       ├── input/                    # counts 校验副本及来源清单
│       ├── results/                  # 唯一当前正式结果；index.html 为入口
│       └── work/                     # staging、cache、logs、history
├── resources/
│   ├── registry.yml                  # release/source/hash/license 元数据
│   ├── sources/                      # 获批联网同步声明
│   └── data/                         # 本地冻结资源，不提交 Git
├── tests/                            # 结构、禁令、数值、回归、渲染测试
└── internal/                         # CLI、schema、发布、审计、上游适配基础设施
```

根目录不得新增 `runs/`、`analysis/`、`scripts/`、`R/` 等平级实现。项目内容只能进入 `projects/<project_id>/`；共享资源只能进入 `resources/`。

## 3. 收到 counts 后必须先问的问题

在导入或计算前逐项确认并记录：

- 项目 ID；是否会与现有项目同名。
- counts 来自 featureCounts、STAR、Salmon/tximport、其他还是未知。
- gene-level 还是 transcript-level；gene ID 类型。
- 物种、assembly、annotation release。
- 是否做过 TPM/FPKM/CPM/RPKM、批次校正或其他归一化。
- 是否包含技术重复；若包含，是否已在上游按相同生物样本求和。
- 上游样本/基因过滤、strandedness 和已知限制。
- 样本到组、批次、配对和协变量的映射；任何排除及理由。
- 设计公式、每个因子的参考水平。
- 每个 contrast 的 factor、numerator、denominator 和方向解释。
- 是否有 FASTQ；本次是否需要服务器上游分析。
- 11 个细分模块的 enabled / not_applicable / skipped_by_user 状态及理由。

`unknown` 可以作为来源限制保留，但正式 DESeq2 前 `normalization` 必须明确为 `raw_counts`。transcript counts 必须先确认 gene-level 汇总方法。整数值只能证明格式可能兼容，不能证明其一定是 raw counts。

## 4. 每次运行必须重新确认的决策

Agent 必须执行：

```bash
pixi run plan -- --project projects/<project_id>
```

并把输出中的每组决策完整展示给用户。至少包括：

| 决策组 | 必须询问/展示的内容 |
|---|---|
| 输入与设计 | 来源、物种/assembly、样本、排除、设计、reference、所有 contrasts、上游选择 |
| 过滤与 DE | `min_count`、`min_samples`、screening 或 LFC-threshold test、每个 padj/abs-LFC profile、主 profile |
| DE 图形 | 火山图标签数、热图每方向基因数；boxplot 基因在 exploratory 中单独列出 |
| ORA/GSEA | ORA 使用哪个 DEG profile、ORA/GSEA 数据库、ORA FDR、GSEA 报告 FDR、set size、冗余阈值、曲线数 |
| 调控 | GSVA 集合、TF 网络、目标 TF、关注通路、最小 targets、TF-GSEA rank/FDR/曲线、报告数量 |
| PPI | DEG profile、最大输入基因数、STRING confidence、module FDR、标签数量、固定布局 seed |
| motif | DEG profile、promoter 窗口、最少 foreground、匹配背景、STREME 宽度/P 阈值、展示 motif 数；是否有 peaks |
| 探索 | boxplot 目标基因、自定义基因集及 rank、Pathview ID、个性化任务、WGCNA 适用性 |
| 报告与导出 | 实际计算模块、报告深度、发布图件范围、格式；默认仅 PDF |

确认后同时使用 `plan` 输出的两个令牌：

```bash
pixi run analyze -- --project projects/<project_id> \
  --confirm-modules '<module_confirmation_token>' \
  --confirm-plan '<plan_confirmation_token>'
```

`plan_confirmation_token` 绑定 `project.yml`、`samples.tsv`、`contrasts.tsv`、`qc_approval.yml` 和 `input/source_manifest.yml`。任一输入或参数变化后旧令牌自动失效。Agent 不得替用户生成“同意”；令牌只能在当前计划已被用户确认后提交。

## 5. 哪些方法是科学锁定项，不作为任意选择

- DESeq2 只使用 raw integer counts；批次/协变量进入设计公式，不在校正表达矩阵上检验。
- QC PCA 使用 blind VST；模型相关表达图使用 `blind = FALSE` 的 VST。
- 一个项目只拟合一个共享 DESeq2 模型，所有已确认 contrasts 复用该模型。
- 展示效应量使用 ashr 收缩 LFC；Wald P 值和 BH padj 不因收缩而替换。
- ORA 背景固定为“通过过滤且成功映射到相应 ID 空间的全部 tested genes”。Up/Down 分开。
- GSEA 使用完整有限排序，不用 DEG 子集；常规通用项目默认推荐 Wald statistic。
- 完整结果表永远保留；top-N、冗余约简和图件数量只是报告层。
- 任何样本排除只能由用户审批，不能依据 PCA 自动删除。
- 随机布局和 permutation 使用配置中记录的固定 seed。

若用户要改变这些方法，Agent 必须先说明影响，建立新的方法版本和回归测试，不能把变更混进一次普通运行。

## 6. 当前硬门禁与支持边界

- `motif_peaks`：当前没有通过通用验证的 peak-aware executor；启用会失败，不能伪装完成。
- `wgcna`：当前没有通过通用验证的执行器；启用会失败。样本不足时应标为 not_applicable。
- `personalized`：必须先提供单独审核的、只读标准产物的执行器，否则启用会失败。
- 通用 `Pathview`：当前通用计算器未完成；历史迁移项目只能复制并校验冻结旧图。新项目启用会失败。
- enrichment、regulation、network 当前仅对已冻结的人类资源完成验证；非人项目启用这些模块会失败，而不是静默套用人类数据库。

这份边界必须在计划阶段向用户说明。只有真正实现并通过测试后，才能删除相应门禁。

## 7. 历史迁移项目

当 `migration.requires_legacy_parity: true` 时：

1. 旧结果是不可退化的回归基线；新增分析不能替换旧能力。
2. 先写 staging，不覆盖 `results/`。
3. 发布前生成 `legacy_analysis_matrix.tsv`、逐文件 `legacy_output_inventory.tsv` 和数值回归报告。
4. 任一旧输出为 `missing`，或 DE、ORA/GSEA、TF/GSVA/PROGENy、PPI、自定义分析回归失败，均阻止发布。
5. 不同 gene-set release、ID、universe、rank、min/max size、seed 或算法都属于统计变更，不是可视化优化。
6. motif 是新增探索模块，不得反向改写旧结果的生物学含义。

## 8. QC、staging 与发布

1. counts 导入后保存 `source_manifest.yml` 与 SHA256。
2. 只读预检格式、行列、样本列、重复 ID、NA/Inf、负值、小数、library size 和零值比例。
3. 输入/设计/过滤和 QC 模块确认后，使用当前 plan token 执行 `qc-preview`。该阶段只生成 QC，状态必须为 `awaiting_user_qc_approval`，不能发布。
4. 用户查看全部 QC/PCA 后，审批排除并更新 `samples.tsv` 与 `qc_approval.yml`；随后必须重新生成 formal plan token。
5. 正式结果先进入 `work/staging/<run_id>/`。
6. 启用模块失败必须为 `failed_explicit`；无显著结果也保留完整表和明确状态，不能制造空图。
7. 发布前验证 `runtime-check`、配置/资源/输入哈希、禁令、单元与集成测试、PDF 签名和实际渲染、HTML 链接、模块状态及 manifest。
8. 已发布结果移入 history 后才可发布新结果；不得覆盖或删除历史。

## 9. 图件与报告

- 默认只输出可编辑字体的矢量 PDF；只有用户明确要求时才增加 PNG。
- 火山图保持用户选定的经典单面板风格，完整显示所有点、阈值线和固定规则标签。
- 不隐藏极端值、异常值或不符合预期的结果，不调换上下调方向。
- 热图必须记录输入矩阵、缩放方向、距离、聚类方法和最终顺序。
- ORA 与 GSEA 不得互相替代；GeneRatio 不是 effect size，NES 的正负必须跟 contrast 一致。
- 每张图应有对应绘图数据或可追溯的标准输入。

## 10. Agent 交付时必须报告

- 唯一 staging 或正式结果路径。
- 本次用户确认的输入、设计、参数、模块和上游选择。
- enabled、not_applicable、skipped_by_user、failed_explicit 状态。
- tests、forbidden check、schema、资源/哈希、PDF/HTML QA 结果。
- 主要 Figures、Tables、HTML 入口和已知限制。
- 任何科学方法变化及对应回归结果。
- 任何清理操作、目标范围和可恢复性。
