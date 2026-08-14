# Scientific workflow

本目录是用户和 Agent 审查科学方法的核心入口：

```text
QC → Differential → Enrichment → Activity → Network → Motif → Exploratory
```

`Snakefile` 只负责编排、资源、缓存和断点续跑。统计逻辑保留在七个编号 R 文件；共享验证、I/O 和绘图函数位于 `functions.R`。`run.R` 仅为兼容串行入口。

## 运行前门禁

1. `project.yml` 通过 schema 和科学条件校验。
2. counts、samples、contrasts、QC 审批及来源清单存在且哈希一致。
3. 所有模块和八组细参数均已在项目级明确确认。
4. 当前运行必须提交 `pixi run plan` 生成的 module token 与 plan token。
5. plan token 绑定所有科学输入；任何编辑后必须重新展示计划并询问用户。

新项目先使用 `qc-preview` 生成不可发布的 QC-only staging。用户检查 PCA/相关性/距离/聚类后审批样本排除；正式分析必须使用审批后重新生成的 plan token。

## 方法摘要

- QC：raw counts 格式验证、明确基因过滤、blind VST PCA、library size、检测率、相关性、距离和聚类；不自动排除样本。
- Differential：检查设计变量、full rank、残差自由度和每个 contrast 的生物学重复；拟合一个 DESeq2 模型；保存 raw 与 ashr-shrunken LFC、Wald P、BH padj、Cook's、size factor、dispersion 和 P-value 诊断。
- Enrichment：ORA 使用 tested mapped universe、Up/Down 分开；GSEA 使用完整 Wald 排序。数据库、阈值、set size、报告 FDR、冗余和曲线数来自配置。
- Activity：冻结 CollecTRI/PROGENy/MSigDB；目标 TF 和关注通路必须存在，不能静默换成其他对象。
- Network：冻结 STRING；DEG profile、输入数量、confidence、module FDR 和布局 seed 显式配置。
- Motif：promoter foreground 使用声明的 DEG profile；背景匹配 GC、长度和 baseMean；窗口和 STREME 参数显式配置。peak-aware executor 尚未验证，启用会失败。
- Exploratory：目标基因 boxplot、自定义 gene-set GSEA/heatmap 和历史 Pathview；rank、FDR、曲线和聚类均显式配置。

默认只导出 Cairo PDF。完整数值表不受报告精简、top-N 或图片格式改变影响。启用模块缺少实现、输入或资源时必须失败，不能静默跳过。
