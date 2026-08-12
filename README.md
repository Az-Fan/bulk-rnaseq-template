# Bulk RNA-seq analysis platform v4

这是一个面向 Linux x86-64 的“一仓库、多项目”Bulk RNA-seq 分析模板。WSL2
适合下游分析；上游 FASTQ→STAR→counts 适配器默认只允许在 Linux 服务器上显式
运行。项目根目录由代码位置自动解析，不绑定用户名或绝对路径。

## 使用者只需关注

```text
bulk-rnaseq-v4/
├── README.md                 # 本入口
├── AGENTS.md                 # Agent 必须遵守的合同和完整目录树
├── pixi.toml / pixi.lock     # 唯一软件环境
├── workflow/                 # 7 个可审查的科学模块
├── projects/
│   └── _template/            # 空白项目模板
├── resources/                # 多项目共享的冻结资源
├── tests/                    # 通用结构、统计与运行测试
└── internal/                 # CLI、schema、发布、QA、上游适配等基础设施
```

新项目不要复制仓库，只复制项目模板：

```bash
pixi install --locked --all
pixi run init-project -- --project-id MY_PROJECT
```

新项目故意处于不可运行状态。Agent 必须先向用户确认 counts 来源、物种/注释、
样本设计、contrast、QC 决定、是否有上游 FASTQ，以及每个分析模块。任何未确认项
都会使校验失败。

## 每次分析都重新确认模块

即使 `project.yml` 已经填写，每次正式分析前 Agent 仍必须把模块清单展示给用户并
获得本次确认。随后执行：

```bash
pixi run modules -- --project projects/MY_PROJECT
pixi run analyze -- --project projects/MY_PROJECT \
  --confirm-modules '<上一条命令打印的完整确认串>'
```

直接调用 R 入口也需要同一确认串，不能绕过。分析先写入
`projects/<id>/work/staging/<run_id>/`；只有完整 QA 通过后才能显式发布到唯一的
`projects/<id>/results/`。

## 图件规则

默认只输出矢量 PDF。只有用户明确要求 PNG 时，才在该项目中加入
`export.formats: [pdf, png]`。火山图使用经典单面板风格，完整显示全部数据点、
阈值线与固定规则标签；不为美观隐藏极端值或改变统计分类。

## 可选上游

`internal/upstream/omics-pipelines` 固定为上游项目的特定 Git 提交。适配器只调用
FASTQ QC、fastp、STAR 和合并 raw gene counts，不调用上游差异分析或富集，也不
创建或下载 STAR index。

默认模板为：

```yaml
upstream:
  status: unconfirmed
  execution_host: server
  allow_wsl: false
```

用户明确选择上游后，Agent 才能准备服务器配置，并仍需提交 `upstream-plan` 打印
的精确确认令牌。服务器应使用：

```bash
git clone --recurse-submodules <repo-url>
```

FASTQ、STAR index、冻结数据库、项目输入、运行缓存和结果均不提交 GitHub。

完整初始化问题、目录职责、科学边界和删除规则见 [AGENTS.md](AGENTS.md)。
