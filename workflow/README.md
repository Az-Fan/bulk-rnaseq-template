# Scientific workflow

这里是用户和 Agent 审查科学方法的唯一核心入口。运行顺序固定为：

QC → Differential → Enrichment → Activity → Network → Motif → Exploratory

`run.R` 只负责编排，7 个编号文件保存相应统计逻辑，`functions.R` 保存共享的科学、
I/O 与绘图函数。设计公式、factor levels、contrasts、模块状态、目标基因和冻结资源均
来自项目配置；启用模块缺少实现或输入时必须明确失败。

每一次分析都要求与当前配置完全一致的模块确认串。图件默认只导出 Cairo PDF；只有
项目记录了用户的 PNG 请求时才同时生成 PNG。火山图保持经典单面板样式。改变报告或
导出选择不得重拟合 DESeq2、重跑 GSEA 或 motif。
