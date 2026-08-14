# Render pre-approval QC only. This stage cannot run DE or be published.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: qc_preview.R <project> <staging_run>")

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
source(file.path(root, "workflow", "functions.R"), local = globalenv())
source(file.path(root, "workflow", "01_qc.R"), local = globalenv())

data.table::setDTthreads(1L)
context <- initialize_context(root, args[[1L]], args[[2L]])
context <- run_qc(context)

pdfs <- list.files(context$run_dir, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
valid <- vapply(pdfs, function(path) {
  connection <- file(path, "rb")
  on.exit(close(connection), add = TRUE)
  identical(tryCatch(rawToChar(readBin(connection, "raw", n = 4L)), error = function(e) ""), "%PDF") &&
    file.info(path)$size >= 1000
}, logical(1))
if (!length(pdfs) || any(!valid)) stop("QC preview PDF validation failed")

write_result_table(data.frame(
  file = substring(pdfs, nchar(context$run_dir) + 2L), format = "pdf", valid = valid
), file.path(context$run_dir, "Provenance", "visual_render_qa.tsv"))
write_json(list(
  schema_version = 1L, project_id = context$config$project_id,
  run_id = basename(context$run_dir), status = "awaiting_user_qc_approval",
  counts_sha256 = context$count_hash, formal_analysis = FALSE,
  module = context$module_status[["01_QC"]],
  instruction = "Inspect QC, then explicitly approve exclusions and regenerate the formal run plan."
), file.path(context$run_dir, "run_manifest.json"))

writeLines(paste0(
  "<!doctype html><html><head><meta charset='utf-8'><title>", context$config$project_id,
  " QC preview</title></head><body><h1>", context$config$project_id,
  " — QC preview only</h1><p>This is not a formal analysis and cannot be published.</p>",
  "<p><a href='01_QC/index.html'>Open QC figures and tables</a></p>",
  "<p>After review, update samples.tsv and qc_approval.yml, then generate a new plan token.</p></body></html>"
), file.path(context$run_dir, "index.html"), useBytes = TRUE)

entries <- setdiff(list.files(file.path(context$run_dir, "01_QC")), "index.html")
links <- paste(sprintf("<li><a href='%s'>%s</a></li>", entries, entries), collapse = "")
writeLines(paste0(
  "<!doctype html><html><head><meta charset='utf-8'><title>QC preview</title></head><body>",
  "<p><a href='../index.html'>Back</a></p><h1>QC preview</h1><ul>", links, "</ul></body></html>"
), file.path(context$run_dir, "01_QC", "index.html"), useBytes = TRUE)
