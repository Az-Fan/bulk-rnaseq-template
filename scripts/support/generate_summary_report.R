#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results/final")
padj_threshold <- as.numeric(Sys.getenv("DEG_PADJ_THRESHOLD", unset = "0.05"))
lfc_threshold <- as.numeric(Sys.getenv("DEG_LFC_THRESHOLD", unset = "1"))
summary_dir <- file.path(output_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
report_html_path <- file.path(summary_dir, "Analysis_Summary_Report.html")

read_csv_safe <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read.csv(path, check.names = FALSE, ...), error = function(e) NULL)
}

read_tsv_safe <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read.delim(path, check.names = FALSE, ...), error = function(e) NULL)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", format(round(x, digits), nsmall = digits, trim = TRUE))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "e", digits = 3))
}

to_md_table <- function(df) {
  if (is.null(df) || nrow(df) == 0) return("（无数据）")
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      out <- ifelse(abs(col) >= 1000, format(round(col, 2), big.mark = ","), as.character(signif(col, 4)))
      out[is.na(col)] <- "NA"
      return(out)
    }
    as.character(col)
  })
  header <- paste0("| ", paste(colnames(df), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

to_html_table <- function(df) {
  if (is.null(df) || nrow(df) == 0) return("<p>（无数据）</p>")
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      out <- ifelse(abs(col) >= 1000, format(round(col, 2), big.mark = ","), as.character(signif(col, 4)))
      out[is.na(col)] <- "NA"
      return(out)
    }
    as.character(col)
  })
  th <- paste0("<th>", html_escape(colnames(df)), "</th>", collapse = "")
  tr <- apply(df, 1, function(r) {
    paste0("<tr>", paste0("<td>", html_escape(r), "</td>", collapse = ""), "</tr>")
  })
  paste0(
    "<table class='summary-table'><thead><tr>", th, "</tr></thead><tbody>",
    paste0(tr, collapse = ""),
    "</tbody></table>"
  )
}

status_files <- list.files(
  file.path(output_dir, "logs"),
  pattern = "^run_status_.*\\.tsv$",
  full.names = TRUE
)
status_file <- if (length(status_files)) {
  status_files[which.max(file.info(status_files)$mtime)]
} else {
  NA_character_
}
status_df <- if (!is.na(status_file)) read_tsv_safe(status_file) else NULL
qc_df <- read_csv_safe(file.path(output_dir, "data_processed", "Sample_QC_Metrics.csv"))
meta_df <- read_csv_safe(file.path(output_dir, "data_processed", "metadata.csv"))
deg_df <- read_csv_safe(file.path(output_dir, "diff_analysis", "DEG_results_annotated.csv"))
model_df <- read_csv_safe(file.path(output_dir, "diff_analysis", "DEG_model_summary.csv"))

go_bp_ora <- read_csv_safe(file.path(output_dir, "functional_enrichment", "go_analysis", "Enrichment_Full_GO_BP.csv"))
hallmark_gsea <- read_csv_safe(file.path(output_dir, "functional_enrichment", "gsea_analysis", "tables", "GSEA_Full_Table_Hallmark.csv"))
reactome_gsea <- read_csv_safe(file.path(output_dir, "functional_enrichment", "gsea_analysis", "tables", "GSEA_Full_Table_Reactome.csv"))
tf_activity <- read_csv_safe(file.path(output_dir, "functional_enrichment", "tf_analysis", "decoupler_inference", "TF_Activity_Contrast_Full.csv"))
custom_gsea <- read_csv_safe(file.path(output_dir, "functional_enrichment", "custom_gene_sets", "gsea", "GSEA_CustomGeneSets_FullTable.csv"))

sample_note <- "（缺少样本质控文件）"
sample_qc_tbl <- NULL
if (!is.null(qc_df)) {
  sample_qc_tbl <- qc_df[, intersect(
    c("sample", "group", "library_size", "detected_genes", "zero_fraction", "mean_pearson_corr", "qc_flag"),
    colnames(qc_df)
  ), drop = FALSE]
  sample_note <- paste0("共 ", nrow(qc_df), " 个样本，", sum(qc_df$qc_flag == "CHECK", na.rm = TRUE), " 个样本被标记为 CHECK。")
}

within_corr <- NA_real_
between_corr <- NA_real_
if (!is.null(qc_df) && !is.null(meta_df)) {
  cor_path <- file.path(output_dir, "data_processed", "Sample_Correlation_Pearson.csv")
  cor_mat <- tryCatch(as.matrix(read.csv(cor_path, row.names = 1, check.names = FALSE)), error = function(e) NULL)
  if (!is.null(cor_mat)) {
    grps <- meta_df$group
    names(grps) <- meta_df$sample
    common <- intersect(colnames(cor_mat), names(grps))
    if (length(common) >= 3) {
      cor_mat <- cor_mat[common, common, drop = FALSE]
      grps <- grps[common]
      within <- c()
      between <- c()
      for (i in seq_along(common)) {
        for (j in seq_along(common)) {
          if (i < j) {
            if (grps[i] == grps[j]) {
              within <- c(within, cor_mat[i, j])
            } else {
              between <- c(between, cor_mat[i, j])
            }
          }
        }
      }
      if (length(within) > 0) within_corr <- mean(within, na.rm = TRUE)
      if (length(between) > 0) between_corr <- mean(between, na.rm = TRUE)
    }
  }
}

deg_summary <- data.frame(
  Metric = c(
    "Total tested genes (non-NA)",
    paste0("Up genes (adj.P<", padj_threshold, " & logFC>", lfc_threshold, ")"),
    paste0("Down genes (adj.P<", padj_threshold, " & logFC<-", lfc_threshold, ")")
  ),
  Value = c(NA, NA, NA),
  stringsAsFactors = FALSE
)
top_up <- NULL
top_down <- NULL
if (!is.null(deg_df)) {
  d2 <- deg_df[!is.na(deg_df$adj.P.Val) & !is.na(deg_df$logFC), , drop = FALSE]
  up <- d2[d2$adj.P.Val < padj_threshold & d2$logFC > lfc_threshold, , drop = FALSE]
  down <- d2[d2$adj.P.Val < padj_threshold & d2$logFC < -lfc_threshold, , drop = FALSE]
  deg_summary$Value <- c(nrow(d2), nrow(up), nrow(down))
  top_up <- head(up[order(up$adj.P.Val, -up$logFC), c("gene", "logFC", "adj.P.Val")], 5)
  top_down <- head(down[order(down$adj.P.Val, down$logFC), c("gene", "logFC", "adj.P.Val")], 5)
}

top_ora <- NULL
if (!is.null(go_bp_ora)) {
  cols <- intersect(c("Description", "p.adjust", "Count", "GeneRatio"), colnames(go_bp_ora))
  top_ora <- head(go_bp_ora[order(go_bp_ora$p.adjust), cols, drop = FALSE], 5)
}

top_hallmark <- NULL
if (!is.null(hallmark_gsea)) {
  cols <- intersect(c("Description", "NES", "p.adjust"), colnames(hallmark_gsea))
  top_hallmark <- head(hallmark_gsea[order(hallmark_gsea$p.adjust), cols, drop = FALSE], 5)
}

top_reactome <- NULL
if (!is.null(reactome_gsea)) {
  cols <- intersect(c("Description", "NES", "p.adjust"), colnames(reactome_gsea))
  top_reactome <- head(reactome_gsea[order(reactome_gsea$p.adjust), cols, drop = FALSE], 5)
}

tf_pos <- NULL
tf_neg <- NULL
if (!is.null(tf_activity)) {
  tf_pos <- head(tf_activity[order(-tf_activity$score), c("source", "score", "p_value")], 5)
  tf_neg <- head(tf_activity[order(tf_activity$score), c("source", "score", "p_value")], 5)
}

custom_top <- NULL
if (!is.null(custom_gsea) && nrow(custom_gsea) > 0) {
  custom_top <- head(custom_gsea[order(custom_gsea$p.adjust), c("ID", "NES", "p.adjust")], 5)
}

kegg_map_dir <- file.path(output_dir, "functional_enrichment", "kegg_analysis", "pathview_maps")
kegg_map_count <- if (dir.exists(kegg_map_dir)) length(list.files(kegg_map_dir, recursive = TRUE)) else 0

lines <- c(
  "# RNA-seq 结果摘要报告",
  "",
  paste0("- 生成时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("- 输出目录: `", output_dir, "`"),
  "",
  "## 1) 运行状态",
  "",
  if (!is.null(status_df)) to_md_table(status_df[, c("step", "script", "status", "seconds"), drop = FALSE]) else "（缺少阶段运行状态表）",
  "",
  "## 2) 样本质控",
  "",
  paste0("- ", sample_note),
  paste0("- 组内平均相关性: ", fmt_num(within_corr, 4)),
  paste0("- 组间平均相关性: ", fmt_num(between_corr, 4)),
  "",
  if (!is.null(sample_qc_tbl)) to_md_table(sample_qc_tbl) else "（无样本 QC 表）",
  "",
  "## 3) 差异表达",
  "",
  if (!is.null(model_df)) to_md_table(model_df) else "（缺少 DEG_model_summary.csv）",
  "",
  to_md_table(deg_summary),
  "",
  "### Top 5 上调基因",
  "",
  if (!is.null(top_up)) to_md_table(top_up) else "（无数据）",
  "",
  "### Top 5 下调基因",
  "",
  if (!is.null(top_down)) to_md_table(top_down) else "（无数据）",
  "",
  "## 4) 富集与调控摘要",
  "",
  "### GO BP ORA Top 5",
  "",
  if (!is.null(top_ora)) to_md_table(top_ora) else "（无数据）",
  "",
  "### Hallmark GSEA Top 5",
  "",
  if (!is.null(top_hallmark)) to_md_table(top_hallmark) else "（无数据）",
  "",
  "### Reactome GSEA Top 5",
  "",
  if (!is.null(top_reactome)) to_md_table(top_reactome) else "（无数据）",
  "",
  "### TF 活性 Top 5 (Activated)",
  "",
  if (!is.null(tf_pos)) to_md_table(tf_pos) else "（无数据）",
  "",
  "### TF 活性 Top 5 (Inhibited)",
  "",
  if (!is.null(tf_neg)) to_md_table(tf_neg) else "（无数据）",
  "",
  "### 自定义基因集 GSEA Top 5",
  "",
  if (!is.null(custom_top)) to_md_table(custom_top) else "（无数据）",
  "",
  "## 5) 关键文件",
  "",
  "- `plots/qc_plots/`：PCA、相关性热图、样本聚类图",
  "- `diff_analysis/DEG_results_annotated.csv`：差异基因全表",
  "- `functional_enrichment/gsea_analysis/tables/`：GSEA 主表",
  "- `functional_enrichment/tf_analysis/decoupler_inference/`：TF 活性结果",
  "- `functional_enrichment/custom_gene_sets/gsea/`：自定义基因集 GSEA 结果",
  "",
  "## 6) 备注",
  "",
  paste0("- KEGG pathview 图文件数量: ", kegg_map_count,
         if (kegg_map_count > 0) "（Step08 已完成）" else "（Step08 未生成图）"),
  "- 报告为自动摘要，结论建议结合原始图件人工复核。"
)

html_lines <- c(
  "<!doctype html>",
  "<html lang='zh-CN'>",
  "<head>",
  "  <meta charset='utf-8'/>",
  "  <meta name='viewport' content='width=device-width, initial-scale=1'/>",
  "  <title>RNA-seq 结果摘要报告</title>",
  "  <style>",
  "    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'PingFang SC', 'Hiragino Sans GB', 'Microsoft YaHei', sans-serif; margin: 28px auto; max-width: 1100px; padding: 0 18px; color: #1f2937; line-height: 1.55; }",
  "    h1 { margin-bottom: 8px; }",
  "    h2 { margin-top: 28px; border-bottom: 1px solid #e5e7eb; padding-bottom: 6px; }",
  "    h3 { margin-top: 20px; }",
  "    ul { margin-top: 4px; }",
  "    code { background: #f3f4f6; padding: 2px 6px; border-radius: 4px; }",
  "    .meta { color: #4b5563; margin: 2px 0; }",
  "    .summary-table { border-collapse: collapse; width: 100%; margin: 8px 0 16px 0; font-size: 14px; }",
  "    .summary-table th, .summary-table td { border: 1px solid #d1d5db; padding: 6px 8px; text-align: left; vertical-align: top; }",
  "    .summary-table thead th { background: #f9fafb; }",
  "  </style>",
  "</head>",
  "<body>",
  "  <h1>RNA-seq 结果摘要报告</h1>",
  paste0("  <p class='meta'>生成时间: ", html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "</p>"),
  paste0("  <p class='meta'>输出目录: <code>", html_escape(output_dir), "</code></p>"),
  "  <h2>1) 运行状态</h2>",
  if (!is.null(status_df)) to_html_table(status_df[, c("step", "script", "status", "seconds"), drop = FALSE]) else "<p>（缺少阶段运行状态表）</p>",
  "  <h2>2) 样本质控</h2>",
  paste0("  <p>", html_escape(sample_note), "</p>"),
  paste0("  <p>组内平均相关性: <b>", html_escape(fmt_num(within_corr, 4)), "</b><br/>组间平均相关性: <b>", html_escape(fmt_num(between_corr, 4)), "</b></p>"),
  if (!is.null(sample_qc_tbl)) to_html_table(sample_qc_tbl) else "<p>（无样本 QC 表）</p>",
  "  <h2>3) 差异表达</h2>",
  if (!is.null(model_df)) to_html_table(model_df) else "<p>（缺少 DEG_model_summary.csv）</p>",
  to_html_table(deg_summary),
  "  <h3>Top 5 上调基因</h3>",
  if (!is.null(top_up)) to_html_table(top_up) else "<p>（无数据）</p>",
  "  <h3>Top 5 下调基因</h3>",
  if (!is.null(top_down)) to_html_table(top_down) else "<p>（无数据）</p>",
  "  <h2>4) 富集与调控摘要</h2>",
  "  <h3>GO BP ORA Top 5</h3>",
  if (!is.null(top_ora)) to_html_table(top_ora) else "<p>（无数据）</p>",
  "  <h3>Hallmark GSEA Top 5</h3>",
  if (!is.null(top_hallmark)) to_html_table(top_hallmark) else "<p>（无数据）</p>",
  "  <h3>Reactome GSEA Top 5</h3>",
  if (!is.null(top_reactome)) to_html_table(top_reactome) else "<p>（无数据）</p>",
  "  <h3>TF 活性 Top 5 (Activated)</h3>",
  if (!is.null(tf_pos)) to_html_table(tf_pos) else "<p>（无数据）</p>",
  "  <h3>TF 活性 Top 5 (Inhibited)</h3>",
  if (!is.null(tf_neg)) to_html_table(tf_neg) else "<p>（无数据）</p>",
  "  <h3>自定义基因集 GSEA Top 5</h3>",
  if (!is.null(custom_top)) to_html_table(custom_top) else "<p>（无数据）</p>",
  "  <h2>5) 关键文件</h2>",
  "  <ul>",
  "    <li><code>plots/qc_plots/</code>：PCA、相关性热图、样本聚类图</li>",
  "    <li><code>diff_analysis/DEG_results_annotated.csv</code>：差异基因全表</li>",
  "    <li><code>functional_enrichment/gsea_analysis/tables/</code>：GSEA 主表</li>",
  "    <li><code>functional_enrichment/tf_analysis/decoupler_inference/</code>：TF 活性结果</li>",
  "    <li><code>functional_enrichment/custom_gene_sets/gsea/</code>：自定义基因集 GSEA 结果</li>",
  "  </ul>",
  "  <h2>6) 备注</h2>",
  "  <ul>",
  paste0("    <li>KEGG pathview 图文件数量: ", kegg_map_count,
         if (kegg_map_count > 0) "（Step08 已完成）</li>" else "（Step08 未生成图）</li>"),
  "    <li>报告为自动摘要，结论建议结合原始图件人工复核。</li>",
  "  </ul>",
  "</body>",
  "</html>"
)

writeLines(html_lines, con = report_html_path, useBytes = TRUE)
cat("Summary report generated:\n", report_html_path, "\n", sep = "")
