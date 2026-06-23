#!/usr/bin/env Rscript

# ==============================================================================
# Shared workflow utility commands
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
command <- if (length(args) >= 1) args[[1]] else ""

usage <- function() {
  cat("Usage: Rscript scripts/support/workflow_utils.R <command>\n")
  cat("Commands:\n")
  cat("  check-dependencies     Check required and optional R packages\n")
  cat("  check-step01-config    Check current config against Step01 outputs\n")
  cat("  write-session-info     Write logs/session_info.txt for current mode\n")
}

check_dependencies <- function() {
  required_packages <- c(
    "data.table",
    "dplyr",
    "tibble",
    "tidyr",
    "stringr",
    "readxl",
    "ggplot2",
    "ggalluvial",
    "ggrepel",
    "Cairo",
    "pheatmap",
    "RColorBrewer",
    "limma",
    "DESeq2",
    "AnnotationDbi",
    "org.Hs.eg.db",
    "clusterProfiler",
    "DOSE",
    "enrichplot",
    "msigdbr",
    "BiocParallel",
    "decoupleR",
    "GSVA",
    "STRINGdb",
    "igraph",
    "ggraph",
    "patchwork",
    "purrr",
    "openxlsx",
    "ggpubr"
  )

  optional_packages <- c(
    "pathview"
  )

  installed <- rownames(installed.packages())
  missing_required <- setdiff(required_packages, installed)
  missing_optional <- setdiff(optional_packages, installed)

  if (length(missing_required) == 0) {
    cat("[OK] Required R packages are installed.\n")
  } else {
    cat("[MISSING] Required R packages:\n")
    cat(paste0("  - ", missing_required, collapse = "\n"), "\n", sep = "")
  }

  if (length(missing_optional) > 0) {
    cat("[OPTIONAL] Missing optional R packages:\n")
    cat(paste0("  - ", missing_optional, collapse = "\n"), "\n", sep = "")
    cat("Optional packages are only needed for optional steps such as KEGG pathview.\n")
  }

  if (length(missing_required) > 0) {
    quit(status = 1)
  }
}

check_step01_config <- function() {
  output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results/final")
  analysis_mode_tag <- Sys.getenv("ANALYSIS_MODE_TAG", unset = "nobatch")
  batch_column <- Sys.getenv("BATCH_COLUMN", unset = "")
  count_matrix_path <- Sys.getenv("COUNT_MATRIX_PATH", unset = "data/matrix_gene.count.xls")
  sample_metadata_path <- Sys.getenv("SAMPLE_METADATA_PATH", unset = "data/sample_metadata.csv")
  control_group <- Sys.getenv("CONTROL_GROUP", unset = "NC")
  treat_group <- Sys.getenv("TREAT_GROUP", unset = "Treatment")
  exclude_samples_raw <- Sys.getenv("EXCLUDE_SAMPLES", unset = "")
  exclude_samples <- trimws(unlist(strsplit(exclude_samples_raw, ",")))
  exclude_samples <- exclude_samples[!is.na(exclude_samples) & exclude_samples != ""]
  exclude_samples <- paste(exclude_samples, collapse = ",")

  config_file <- file.path(output_dir, "data_processed", "Analysis_Config_Used.csv")
  if (!file.exists(config_file)) {
    stop(paste0("Error: missing Step01 config file: ", config_file))
  }

  used_config <- read.csv(config_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(used_config) < 1) {
    stop(paste0("Error: empty Step01 config file: ", config_file))
  }
  used_config <- used_config[1, , drop = FALSE]

  current_config <- list(
    count_matrix_path = count_matrix_path,
    sample_metadata_path = sample_metadata_path,
    control_group = control_group,
    treat_group = treat_group,
    analysis_mode_tag = analysis_mode_tag,
    batch_column = batch_column,
    exclude_samples = exclude_samples
  )

  required_fields <- names(current_config)
  missing_fields <- setdiff(required_fields, colnames(used_config))
  if (length(missing_fields) > 0) {
    stop(paste0("Error: Step01 config is missing fields: ", paste(missing_fields, collapse = ", ")))
  }

  mismatch <- character(0)
  for (field in required_fields) {
    old_value <- as.character(used_config[[field]][1])
    new_value <- as.character(current_config[[field]])
    if (is.na(old_value)) {
      old_value <- ""
    }
    if (is.na(new_value)) {
      new_value <- ""
    }
    if (!identical(old_value, new_value)) {
      mismatch <- c(mismatch, field)
    }
  }

  if (length(mismatch) > 0) {
    cat("Current config differs from existing Step01 outputs.\n")
    cat("Mismatched fields:", paste(mismatch, collapse = ", "), "\n")
    cat("You changed BATCH_COLUMNS, EXCLUDE_SAMPLES, group contrast, count matrix, or metadata after QC.\n")
    cat("Please run: bash start.sh full\n")
    cat("Do not run downstream on stale Step01 outputs.\n")
    quit(status = 1)
  }

  cat("Step01 config matches current config.\n")
}

write_session_info <- function() {
  output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results/final")
  log_dir <- file.path(output_dir, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  session_file <- file.path(log_dir, "session_info.txt")

  sink(session_file)
  cat("R version:\n")
  print(R.version.string)
  cat("\nSession info:\n")
  print(sessionInfo())
  sink()

  cat("Session info written to:", session_file, "\n")
}

if (command == "check-dependencies") {
  check_dependencies()
} else if (command == "check-step01-config") {
  check_step01_config()
} else if (command == "write-session-info") {
  write_session_info()
} else {
  usage()
  quit(status = 1)
}

