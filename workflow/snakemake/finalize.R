# Seal a Snakemake staging run. Statistical work is already complete; this
# script validates terminal states and creates navigation/provenance only.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: finalize.R <project> <run>")

project <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
run <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[[1L]])), "../.."), winslash = "/", mustWork = TRUE)
cfg <- yaml::read_yaml(file.path(project, "project.yml"))
modules <- c("01_QC", "02_Differential", "03_Enrichment", "04_Regulation", "05_Network", "06_Motif", "07_Exploratory")
valid <- c("complete", "not_applicable", "skipped_by_user")
status <- lapply(modules, function(module) {
  path <- file.path(run, "Provenance", paste0(module, "_status.tsv"))
  if (!file.exists(path)) stop("Missing module status: ", path)
  row <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(row) != 1L || !row$status[[1L]] %in% valid) stop("Invalid terminal module state: ", module)
  list(status = row$status[[1L]], detail = row$detail[[1L]])
})
names(status) <- modules

pdfs <- list.files(run, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
pdf_valid <- vapply(pdfs, function(path) {
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  identical(tryCatch(rawToChar(readBin(con, "raw", n = 4L)), error = function(e) ""), "%PDF") && file.info(path)$size >= 1000
}, logical(1))
qa <- data.frame(file = substring(pdfs, nchar(run) + 2L), format = "pdf", valid = pdf_valid, stringsAsFactors = FALSE)
utils::write.table(qa, file.path(run, "Provenance", "visual_render_qa.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
if (length(pdf_valid) && any(!pdf_valid)) stop("Rendered PDF QA failed")

config_dir <- file.path(run, "Provenance", "Configuration")
dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
for (name in c("project.yml", "samples.tsv", "contrasts.tsv", "qc_approval.yml")) {
  source <- file.path(project, name)
  if (file.exists(source)) file.copy(source, file.path(config_dir, name), overwrite = TRUE)
}
writeLines(capture.output(sessionInfo()), file.path(run, "Provenance", "session_info.txt"))
writeLines(c(
  "# Methods", "", paste0("Design: `", cfg$design$formula, "`."),
  "The seven public R modules preserve the declared legacy scientific contract; Snakemake manages dependencies, resources and resumability only.",
  "DESeq2 uses raw integer counts. ORA uses the tested mapped-gene universe and separate directions. GSEA uses the complete declared ranking.",
  "Frozen resources were SHA256 verified before formal offline analysis."
), file.path(run, "Provenance", "METHODS.md"))
writeLines(c(
  "# Limitations", "", "Expression-derived TF/pathway activity, motif enrichment and network centrality are associative, not causal or direct-binding evidence.",
  "Unknown upstream provenance remains explicitly reported; integer values alone do not prove raw-count origin."
), file.path(run, "Provenance", "LIMITATIONS.md"))

for (module in modules) {
  directory <- file.path(run, module); dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  entries <- setdiff(list.files(directory), "index.html")
  links <- paste(sprintf("<li><a href='%s'>%s</a></li>", entries, entries), collapse = "")
  page <- paste0("<!doctype html><html><head><meta charset='utf-8'><title>", module,
                 "</title></head><body><p><a href='../index.html'>Back to results</a></p><h1>", module,
                 "</h1><p>Status: ", status[[module]]$status, "</p><ul>", links, "</ul></body></html>")
  writeLines(page, file.path(directory, "index.html"), useBytes = TRUE)
}
cards <- paste(vapply(modules, function(module) sprintf("<li><a href='%s/index.html'>%s</a> — %s</li>", module, module, status[[module]]$status), character(1)), collapse = "")
writeLines(paste0("<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><title>", cfg$project_id,
                  " results</title></head><body><h1>", cfg$project_id, " — validated staging result</h1><ul>", cards,
                  "</ul><p><a href='Provenance/'>Provenance and validation</a></p></body></html>"), file.path(run, "index.html"), useBytes = TRUE)

manifest <- list(schema_version = 2L, project_id = cfg$project_id, run_id = basename(run), status = "complete",
                 updated_at = format(Sys.time(), tz = "UTC", usetz = TRUE), offline_formal_run = TRUE,
                 orchestrator = "Snakemake", modules = status, figures = list(pdf = length(pdfs), png = length(list.files(run, pattern = "\\.png$", recursive = TRUE))))
jsonlite::write_json(manifest, file.path(run, "run_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
