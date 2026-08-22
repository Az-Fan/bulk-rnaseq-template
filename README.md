# Bulk RNA-seq analysis platform v4

面向 Linux x86-64 的“一仓库、多项目”Bulk RNA-seq 平台。WSL2 适合下游分析；可选 FASTQ→STAR→counts 上游默认只允许在 Linux 服务器上显式运行。项目根目录由代码位置自动解析，不绑定用户名或个人绝对路径。

GitHub 保存的是可迁移的分析平台源码：七个科学模块、控制层、空白项目模板、固定 Pixi 环境和资源同步声明。真实 counts、项目结果、运行缓存及受许可约束的冻结资源不属于源码仓库，也不需要随模板迁移。

## 使用者只需关注

```text
bulk-rnaseq-v4/
├── README.md                 # 本入口
├── AGENTS.md                 # Agent 的硬约束、决策清单和完整目录树
├── pixi.toml / pixi.lock     # 唯一软件环境
├── .pixi/config.toml         # 允许锁定 Bioconda 数据包完成校验式 post-link 安装
├── workflow/                 # Snakemake + 7 个可审查科学模块
├── projects/
│   ├── _template/            # 唯一空白模板
│   └── <project_id>/         # 每个真实项目的数据、配置、结果和工作区
├── resources/                # 多项目共享的冻结资源
├── tests/                    # 结构、统计、回归、渲染和禁令测试
└── internal/                 # CLI、schema、发布、审计和上游适配
```

核心科学代码只有七个编号模块：QC、差异表达、富集、调控、网络、motif、探索分析。详细职责和完整目录树见 [AGENTS.md](AGENTS.md)。
当前科学支持边界和逐模块审计见 [workflow/SCIENTIFIC_AUDIT.md](workflow/SCIENTIFIC_AUDIT.md)。

## 新项目

```bash
pixi install --locked --all
pixi run init-project -- --project-id MY_PROJECT
```

## 在 Linux 服务器重建平台

```bash
git clone --recurse-submodules https://github.com/Az-Fan/bulk-rnaseq-template.git
cd bulk-rnaseq-template
pixi install --locked --all
pixi run runtime-check
pixi run forbidden-check
pixi run test
pixi run doctor
```

仓库跟踪的 `.pixi/config.toml` 允许锁文件中的 Bioconda annotation data 包执行 post-link。该步骤会在首次环境安装时联网下载固定版本的 `org.Hs.eg.db`、`GO.db` 和 `reactome.db`，并使用 Bioconda 随包元数据中的 MD5 校验；正式分析仍不得联网。`runtime-check` 会验证这些包和本地 Pixi 构建的 aPEAR 是否真实存在且版本正确，防止出现“Pixi 显示安装成功但数据库包实际缺失”。

只有启用依赖外部数据库的模块时才同步对应资源。例如 promoter motif：

```bash
pixi run resources-sync -- \
  --manifest resources/sources/hg38_motif.yml \
  --destination resources/data/shared/hg38-motif
```

资源文件与软件环境严格分开；同步声明固定 release、URL 和 SHA256。项目特有或受许可限制的资源由用户导入并登记，不进入 GitHub。

新项目故意不可直接运行。Agent 必须先确认 counts 来源、物种/annotation、样本设计、contrast、QC、模块，以及会影响结果或展示的参数。

导入 counts、填写样本/设计并确认输入过滤参数后，先生成计划并只运行 QC 预览：

```bash
pixi run plan -- --project projects/MY_PROJECT
pixi run qc-preview -- --project projects/MY_PROJECT \
  --confirm-plan '<plan_confirmation_token>'
```

QC 预览不会运行 DE/ORA，也不能发布。查看全部 QC/PCA 后，由用户明确更新 `samples.tsv` 排除项和 `qc_approval.yml`。随后完成其他参数选择，再执行 `check-project` 和正式计划。

## 每次运行都要重新确认

```bash
pixi run plan -- --project projects/MY_PROJECT
```

该命令列出本次必须由用户确认的完整计划：输入与设计、过滤/DEG 阈值、ORA/GSEA、TF/GSVA/PROGENy、PPI、motif、boxplot/自定义分析、报告和导出。

用户在当前对话确认后，使用输出中的两个令牌：

```bash
pixi run analyze -- --project projects/MY_PROJECT \
  --confirm-modules '<module_confirmation_token>' \
  --confirm-plan '<plan_confirmation_token>'
```

计划令牌绑定项目配置、样本表、contrasts、QC 审批和 counts 来源清单；任一内容变化都会使旧令牌失效。结果先写入 `projects/<id>/work/staging/<run_id>/`，通过 QA 后才能显式发布到唯一的 `projects/<id>/results/`。

## 科学默认与边界

- DESeq2 只接受明确确认的 raw integer counts。
- 批次和协变量进入设计公式；不在校正表达值上做 DESeq2。
- QC PCA 使用 blind VST；模型相关表达图使用 `blind = FALSE` VST。
- 一个共享 DESeq2 模型供所有 contrasts 使用；效应量展示使用 ashr 收缩。
- ORA 使用 tested mapped genes 作为 universe，并分开 Up/Down。
- GSEA 使用完整排序；通用项目默认推荐 Wald statistic。
- top-N、冗余约简和曲线数只改变展示，不删完整结果表。
- 默认仅导出矢量 PDF；只有用户明确要求才增加 PNG。
- motif、TF/PROGENy activity 和 PPI 只能解释为关联或探索证据。

当前尚未通过通用验证的 peak-aware motif、WGCNA、personalized executor 和新项目 Pathview 会被硬门禁阻止，不能静默标记完成。人类以外物种若启用尚无冻结物种资源的模块，也会明确失败。

## 可选服务器上游

`internal/upstream/omics-pipelines` 固定为 `xuzhougeng/omics-pipelines` 的指定提交。适配器只调用 FASTQ QC、fastp、STAR 和 raw gene-count merge，不调用其下游 DESeq2 或 enrichment，不在 WSL 建立或下载 STAR index。

只有用户明确选择上游后才执行：

```bash
pixi run upstream-plan -- --project projects/MY_PROJECT
```

服务器克隆应使用：

```bash
git clone --recurse-submodules <repo-url>
pixi install --locked --all
```

FASTQ、STAR index、冻结数据库、项目输入、缓存和结果默认均不提交 GitHub。

## 验证

```bash
pixi run doctor
pixi run runtime-check
pixi run forbidden-check
pixi run test
```

`doctor` 检查路径可迁移性和本机资源；禁令测试确保没有第二套 R 安装流程；测试覆盖 schema、确认门禁、通用多 contrast、历史输出契约和图形发布。

## GitHub CLI 授权

`gh` 已纳入同一个 Pixi 锁文件。不要把 Personal Access Token 发到聊天或写进 remote URL。在 WSL 仓库根目录执行：

```bash
pixi run gh auth login --hostname github.com --git-protocol https --web
pixi run gh auth refresh --hostname github.com --scopes repo,workflow
pixi run gh auth status
```

浏览器中只授权你希望本仓库使用的 GitHub 账户。当前 remote 为 `https://github.com/Az-Fan/bulk-rnaseq-template.git`；授权完成后仍需用户明确批准，Agent 才能 commit/push。
