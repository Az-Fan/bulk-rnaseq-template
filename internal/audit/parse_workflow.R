#!/usr/bin/env Rscript
files <- c("workflow/01_qc.R", "workflow/02_differential.R", "workflow/03_enrichment.R",
           "workflow/04_activity.R", "workflow/05_network.R", "workflow/06_motif.R",
           "workflow/07_exploratory.R", "workflow/functions.R", "workflow/run.R",
           "workflow/snakemake/run_module.R", "workflow/snakemake/finalize.R",
           "workflow/snakemake/qc_preview.R")
invisible(lapply(files, parse))
cat("R parse OK:", paste(files, collapse = ", "), "\n")
