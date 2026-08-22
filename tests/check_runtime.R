required <- c(
  "aPEAR" = "1.0.0",
  "org.Hs.eg.db" = "3.22.0",
  "GO.db" = "3.22.0",
  "reactome.db" = "1.95.0"
)

missing <- names(required)[!vapply(
  names(required),
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)]
if (length(missing)) {
  stop("Missing required runtime package(s): ", paste(missing, collapse = ", "))
}

actual <- vapply(
  names(required),
  function(package) as.character(utils::packageVersion(package)),
  character(1)
)
wrong <- names(required)[actual != required]
if (length(wrong)) {
  details <- paste0(wrong, "=", actual[wrong], " (expected ", required[wrong], ")")
  stop("Runtime package version mismatch: ", paste(details, collapse = "; "))
}

cat("Runtime packages verified:\n")
for (package in names(required)) {
  cat(sprintf("- %s %s\n", package, actual[[package]]))
}
