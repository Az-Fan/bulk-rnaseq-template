# Execute exactly one public scientific module and persist its state for the DAG.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(4L, 5L)) {
  stop("Usage: run_module.R <initialize|module> <project> <run> <output_state> [input_state]")
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
module <- args[[1L]]
project <- args[[2L]]
run <- args[[3L]]
output_state <- args[[4L]]

source(file.path(root, "workflow", "functions.R"), local = globalenv())
module_files <- c("01_qc.R", "02_differential.R", "03_enrichment.R", "04_activity.R",
                  "05_network.R", "06_motif.R", "07_exploratory.R")
for (path in module_files) source(file.path(root, "workflow", path), local = globalenv())

if (identical(module, "initialize")) {
  context <- initialize_context(root, project, run)
} else {
  input_state <- args[[5L]]
  context <- readRDS(input_state)
  calls <- list(
    `01_QC` = run_qc,
    `02_Differential` = run_differential,
    `03_Enrichment` = run_enrichment,
    `04_Regulation` = run_activity,
    `05_Network` = run_network,
    `06_Motif` = run_motif,
    `07_Exploratory` = run_exploratory
  )
  if (is.null(calls[[module]])) stop("Unknown public module: ", module)
  context <- calls[[module]](context)
}

dir.create(dirname(output_state), recursive = TRUE, showWarnings = FALSE)
saveRDS(context, output_state, compress = "xz")
