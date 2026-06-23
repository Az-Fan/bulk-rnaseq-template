options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
})

output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results")
project_name <- Sys.getenv("PROJECT_NAME", unset = "Project")
padj_cutoff <- as.numeric(Sys.getenv("DEG_PADJ_THRESHOLD", unset = "0.05"))
lfc_cutoff <- as.numeric(Sys.getenv("DEG_LFC_THRESHOLD", unset = "1"))
diff_dir <- file.path(output_dir, "diff_analysis")
input_path <- file.path(diff_dir, "DEG_results_annotated.csv")
if (!file.exists(input_path)) stop("Missing DEG table: ", input_path)

x <- read.csv(input_path, check.names = FALSE)
x$Direction <- ifelse(
  !is.na(x$adj.P.Val) & x$adj.P.Val < padj_cutoff & x$logFC > lfc_cutoff,
  "Up",
  ifelse(
    !is.na(x$adj.P.Val) & x$adj.P.Val < padj_cutoff & x$logFC < -lfc_cutoff,
    "Down",
    "Not significant"
  )
)
x$DEG_Status <- ifelse(x$Direction %in% c("Up", "Down"), "DEG", "Not DEG")
x$P_value_display <- ifelse(
  is.na(x$P.Value),
  NA_character_,
  ifelse(x$P.Value == 0, "<1e-300", format(x$P.Value, scientific = TRUE, digits = 5))
)
x$Adjusted_P_display <- ifelse(
  is.na(x$adj.P.Val),
  NA_character_,
  ifelse(x$adj.P.Val == 0, "<1e-300", format(x$adj.P.Val, scientific = TRUE, digits = 5))
)
x <- x %>%
  arrange(adj.P.Val, desc(abs(stat))) %>%
  mutate(Rank = row_number()) %>%
  select(
    Rank, gene_id, gene, ENTREZID, AveExpr, logFC_raw, logFC, lfcSE, stat,
    P.Value, adj.P.Val, P_value_display, Adjusted_P_display,
    Direction, DEG_Status
  )
sig <- x %>% filter(DEG_Status == "DEG")

wb <- createWorkbook()
addWorksheet(wb, "Summary", gridLines = FALSE)
addWorksheet(wb, "Significant_DEGs", gridLines = FALSE)
addWorksheet(wb, "All_Tested_Genes", gridLines = FALSE)

title_style <- createStyle(
  fgFill = "#17365D", fontColour = "#FFFFFF",
  fontSize = 18, textDecoration = "bold"
)
header_style <- createStyle(
  fgFill = "#4472C4", fontColour = "#FFFFFF",
  textDecoration = "bold", halign = "center", valign = "center",
  wrapText = TRUE
)
label_style <- createStyle(
  fgFill = "#D9EAF7", fontColour = "#17365D", textDecoration = "bold"
)
sci_style <- createStyle(numFmt = "0.00E+00")

writeData(wb, "Summary", paste0(project_name, " differential-expression results"), startRow = 1)
mergeCells(wb, "Summary", cols = 1:8, rows = 1)
addStyle(wb, "Summary", title_style, rows = 1, cols = 1:8, gridExpand = TRUE)

summary_left <- data.frame(
  Item = c("Analysis", "Comparison", "Formal DEG rule", "Tested genes", "Formal DEGs", "Up", "Down"),
  Value = c(
    "DESeq2 count model with ashr log2FC shrinkage",
    paste(Sys.getenv("TREAT_GROUP"), "vs", Sys.getenv("CONTROL_GROUP")),
    paste0("Adjusted P < ", padj_cutoff, " and |ashr-shrunken log2FC| > ", lfc_cutoff),
    nrow(x), nrow(sig), sum(sig$Direction == "Up"), sum(sig$Direction == "Down")
  )
)
writeData(wb, "Summary", summary_left, startRow = 3, colNames = FALSE)
addStyle(wb, "Summary", label_style, rows = 3:9, cols = 1, gridExpand = TRUE)
writeData(
  wb, "Summary",
  "Exact numeric zeros are floating-point underflow for extremely small values. Use the display columns (<1e-300) for reporting; numeric columns are retained for computation.",
  startRow = 11
)
mergeCells(wb, "Summary", cols = 1:8, rows = 11)
setColWidths(wb, "Summary", cols = 1, widths = 22)
setColWidths(wb, "Summary", cols = 2:8, widths = 20)

writeDataTable(wb, "Significant_DEGs", sig, tableStyle = "TableStyleMedium2")
writeDataTable(wb, "All_Tested_Genes", x, tableStyle = "TableStyleMedium2")
for (sheet in c("Significant_DEGs", "All_Tested_Genes")) {
  freezePane(wb, sheet, firstActiveRow = 2, firstActiveCol = 4)
  setColWidths(wb, sheet, cols = 1:ncol(x), widths = "auto")
  setColWidths(wb, sheet, cols = c(2, 12, 13), widths = c(21, 16, 20))
  addStyle(wb, sheet, header_style, rows = 1, cols = 1:ncol(x), gridExpand = TRUE)
  addStyle(wb, sheet, sci_style, rows = 2:(nrow(if (sheet == "Significant_DEGs") sig else x) + 1),
           cols = c(10, 11), gridExpand = TRUE)
}

out <- file.path(diff_dir, paste0("00_", project_name, "_Differential_Expression.xlsx"))
saveWorkbook(wb, out, overwrite = TRUE)
message("DEG workbook written: ", out)
