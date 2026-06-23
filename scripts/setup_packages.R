cran_packages <- c(
  "data.table", "dplyr", "tibble", "tidyr", "stringr", "readxl",
  "ggplot2", "ggalluvial", "ggrepel", "Cairo", "pheatmap",
  "RColorBrewer", "msigdbr", "igraph", "ggraph", "patchwork",
  "purrr", "openxlsx", "ggpubr", "yaml", "jsonlite", "scales"
)
bioc_packages <- c(
  "limma", "DESeq2", "AnnotationDbi", "org.Hs.eg.db",
  "clusterProfiler", "DOSE", "enrichplot", "BiocParallel",
  "decoupleR", "GSVA", "STRINGdb"
)

local_library <- file.path(getwd(), ".r-library")
dir.create(local_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_library, .libPaths()))

cran_missing <- setdiff(cran_packages, rownames(installed.packages()))
if (length(cran_missing)) {
  install.packages(cran_missing, repos = "https://cloud.r-project.org", lib = local_library)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = local_library)
}
bioc_missing <- setdiff(bioc_packages, rownames(installed.packages()))
if (length(bioc_missing)) {
  BiocManager::install(bioc_missing, ask = FALSE, update = FALSE, lib = local_library)
}
cat("Dependencies installed.\n")
