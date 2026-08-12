# Public workflow orchestrator
#
# This file contains no statistical method. It constructs the confirmed project
# context, calls the seven public modules in order, validates rendered outputs,
# and seals a staging run only when every module has an explicit valid state.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: pixi run analyze-staging -- <project_dir> <staging_run_dir>")

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(file_argument)) stop("Cannot resolve workflow/run.R location")
script_path <- normalizePath(sub("^--file=", "", file_argument[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(dirname(dirname(script_path)), winslash = "/", mustWork = TRUE)
data.table::setDTthreads(1L)
set.seed(104729L)

module_files <- c("01_qc.R", "02_differential.R", "03_enrichment.R", "04_activity.R",
                  "05_network.R", "06_motif.R", "07_exploratory.R")
source(file.path(root, "workflow", "functions.R"), local = globalenv())
for (file in module_files) source(file.path(root, "workflow", file), local = globalenv())

context <- initialize_context(root, args[[1L]], args[[2L]])
module_calls <- list(
  `01_QC` = run_qc,
  `02_Differential` = run_differential,
  `03_Enrichment` = run_enrichment,
  `04_Regulation` = run_activity,
  `05_Network` = run_network,
  `06_Motif` = run_motif,
  `07_Exploratory` = run_exploratory
)

manifest_path <- file.path(context$run_dir, "run_manifest.json")
motif_phase <- identical(Sys.getenv("BULK_RNASEQ_MOTIF_PHASE"), "1")
write_manifest <- function(status, error = NULL) {
  module_status <- lapply(names(module_calls), function(name) context$module_status[[name]] %||%
                            list(status = "not_started", detail = ""))
  names(module_status) <- names(module_calls)
  write_json(list(
    schema_version = 2, project_id = context$config$project_id,
    run_id = basename(context$run_dir), status = status,
    started_at = context$started_at,
    updated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    counts_sha256 = context$count_hash,
    configuration_sha256 = list(
      project_yml = sha256_file(file.path(context$project_dir, "project.yml")),
      samples_tsv = sha256_file(context$paths$samples), contrasts_tsv = sha256_file(context$paths$contrasts),
      qc_approval_yml = sha256_file(context$paths$approval), source_manifest_yml = sha256_file(context$paths$source_manifest)
    ),
    workflow_sha256 = stats::setNames(vapply(c("functions.R", "run.R", module_files), function(file) {
      sha256_file(file.path(root, "workflow", file))
    }, character(1)), c("functions.R", "run.R", module_files)),
    offline_formal_run = TRUE, random_seed = as.integer(context$config$analysis$random_seed %||% 104729L),
    design_formula = context$config$design$formula,
    contrasts = as.character(context$contrasts$contrast_id), modules = module_status,
    error = error
  ), manifest_path)
}
if (motif_phase) {
  existing_manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  existing_state <- readRDS(file.path(context$run_dir, "Provenance", "workflow_context_state.rds"))
  context$state <- existing_state$state
  context$module_status <- existing_state$module_status
} else {
  write_manifest("running")
}

active_module <- NULL
tryCatch({
  modules_to_run <- if (motif_phase) "06_Motif" else names(module_calls)
  for (name in modules_to_run) {
    active_module <- name
    context <- module_calls[[name]](context)
    write_manifest("running")
  }
}, error = function(error) {
  detail <- conditionMessage(error)
  if (!is.null(active_module)) context <<- record_module_status(context, active_module, "failed_explicit", detail)
  write_manifest("failed", detail)
  stop(error)
})

if (!motif_phase && identical(context$module_status[["06_Motif"]]$status, "pending_motif_phase")) {
  saveRDS(list(state = list(annotation = context$state$annotation,
                            differential = context$state$differential),
               module_status = context$module_status),
          file.path(context$run_dir, "Provenance", "workflow_context_state.rds"), compress = "xz")
  write_manifest("awaiting_motif")
  message("Core phase complete; motif phase required: ", context$run_dir)
  quit(save = "no", status = 42L)
}
valid_states <- c("complete", "not_applicable", "skipped_by_user")
states <- vapply(context$module_status, `[[`, character(1), "status")
if (!identical(names(context$module_status), names(module_calls)) || any(!states %in% valid_states)) {
  write_manifest("failed", "One or more modules lack a valid terminal state")
  stop("One or more modules lack a valid terminal state")
}

# PDF is the formal default. Check each PDF signature and size; when PNG was
# explicitly requested, also decode it and verify a matching PDF companion.
pdf_files <- list.files(context$run_dir, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
pdf_rows <- lapply(pdf_files, function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  signature <- tryCatch(rawToChar(readBin(con, "raw", n = 4L)), error = function(e) "")
  data.frame(file = substring(path, nchar(context$run_dir) + 2L), format = "pdf",
             valid = identical(signature, "%PDF") && file.info(path)$size >= 1000,
             width_px = NA_integer_, height_px = NA_integer_,
             pdf_companion = TRUE, stringsAsFactors = FALSE)
})
png_files <- list.files(context$run_dir, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
png_rows <- lapply(png_files, function(path) {
  image <- tryCatch(png::readPNG(path, info = TRUE), error = identity)
  ok <- !inherits(image, "error")
  dimensions <- if (ok) dim(image) else c(NA_integer_, NA_integer_)
  stem <- sub("\\.png$", "", path)
  frozen_resource_raster <- grepl("/Pathview/", path, fixed = TRUE)
  data.frame(file = substring(path, nchar(context$run_dir) + 2L), format = "png", valid = ok,
             width_px = dimensions[[2L]], height_px = dimensions[[1L]],
             pdf_companion = file.exists(paste0(stem, ".pdf")),
             vector_exempt_reason = if (frozen_resource_raster) "frozen_database_raster" else "",
             stringsAsFactors = FALSE)
})
qa_rows <- c(pdf_rows, png_rows)
qa <- if (length(qa_rows)) data.table::rbindlist(qa_rows, fill = TRUE) else data.frame(
  file = character(), format = character(), valid = logical(), width_px = integer(),
  height_px = integer(), pdf_companion = logical(), vector_exempt_reason = character())
write_result_table(qa, file.path(context$run_dir, "Provenance", "visual_render_qa.tsv"))
png_qa <- qa[qa$format == "png", , drop = FALSE]
if ((nrow(qa) && any(!qa$valid)) ||
    (nrow(png_qa) && any(png_qa$width_px < 800 | png_qa$height_px < 600 |
                         (!png_qa$pdf_companion & png_qa$vector_exempt_reason == "")))) {
  write_manifest("failed", "Rendered figure QA failed")
  stop("Rendered figure QA failed")
}

configuration_dir <- file.path(context$run_dir, "Provenance", "Configuration")
dir.create(configuration_dir, recursive = TRUE, showWarnings = FALSE)
for (path in c(file.path(context$project_dir, "project.yml"), context$paths$samples, context$paths$contrasts,
               context$paths$approval, context$paths$source_manifest)) {
  file.copy(path, file.path(configuration_dir, basename(path)), overwrite = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(context$run_dir, "Provenance", "session_info.txt"))
methods <- c(
  "# Methods", "",
  paste0("The confirmed design formula was `", context$config$design$formula, "`."),
  "Raw integer counts were filtered using the declared minimum count and minimum sample rules.",
  "One DESeq2 model was fitted and reused for every confirmed contrast. ashr shrinkage was used for displayed LFC estimates.",
  paste0("DEG decision mode: `", context$config$thresholds$mode %||% "screening", "`."),
  "ORA used tested mapped genes as the universe and separated Up and Down sets. GSEA used the complete finite Wald-statistic ranking.",
  "All resources were read from the frozen SHA256-verified registry; formal analysis did not update resources or install software."
)
writeLines(methods, file.path(context$run_dir, "Provenance", "METHODS.md"), useBytes = TRUE)
limitations <- c(
  "# Limitations", "",
  paste0("Counts source method: `", context$config$counts_provenance$source_method, "`; annotation: `",
         context$config$counts_provenance$genome_annotation %||% "unknown", "`."),
  "Integer values are necessary but do not independently prove raw-count provenance.",
  "Expression-derived pathway/TF activity, motif enrichment and network centrality are associative or exploratory, not causal evidence.",
  "Any approved sample exclusion is a user-confirmed design decision; the workflow does not remove samples automatically."
)
writeLines(limitations, file.path(context$run_dir, "Provenance", "LIMITATIONS.md"), useBytes = TRUE)

# Provide a local-file-friendly index without duplicating scientific outputs.
module_cards <- vapply(names(module_calls), function(name) {
  status <- context$module_status[[name]]$status
  sprintf("<li><a href='%s/'>%s</a> — %s</li>", name, name, status)
}, character(1))
for (name in names(module_calls)) {
  directory <- file.path(context$run_dir, name)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(directory)
  page <- paste0("<!doctype html><html><head><meta charset='utf-8'><title>", name,
                 "</title></head><body><p><a href='../index.html'>← Results</a></p><h1>", name,
                 "</h1><p>Status: ", context$module_status[[name]]$status, "</p><ul>",
                 paste(sprintf("<li><a href='%s'>%s</a></li>", entries, entries), collapse = ""),
                 "</ul></body></html>")
  writeLines(page, file.path(directory, "index.html"), useBytes = TRUE)
}
index <- paste0("<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><title>",
                context$config$project_id, " results</title></head><body><h1>", context$config$project_id,
                " — validated staging result</h1><p>Open each module below. Comparisons: ",
                paste(context$contrasts$contrast_id, collapse = ", "), ".</p><ul>", paste(module_cards, collapse = ""),
                "</ul><p><a href='Provenance/'>Provenance</a></p></body></html>")
writeLines(index, file.path(context$run_dir, "index.html"), useBytes = TRUE)
write_manifest("complete")
state_cache <- file.path(context$run_dir, "Provenance", "workflow_context_state.rds")
if (file.exists(state_cache)) unlink(state_cache)
message("Validated staging analysis complete: ", context$run_dir)
