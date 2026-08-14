#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: inspect_msigdbr_resource.R RESOURCE_RDS")
x <- readRDS(args[[1]])
sets <- x$gene_sets
cat(sprintf("release=%s gene_sets=%d memberships=%d\n",
            x$msigdb_release, length(unique(sets$gs_name)), nrow(sets)))
for (key in c("H|", "C2|CP:KEGG", "C2|CP:REACTOME", "C3|TFT:GTRD",
              "C5|GO:BP", "C5|GO:CC", "C5|GO:MF")) {
  parts <- strsplit(paste0(key, " "), "|", fixed = TRUE)[[1]]
  parts[[2]] <- trimws(parts[[2]])
  keep <- sets$gs_cat == parts[[1]]
  if (nzchar(parts[[2]])) keep <- keep & sets$gs_subcat == parts[[2]]
  cat(sprintf("%s gene_sets=%d memberships=%d\n", key,
              length(unique(sets$gs_name[keep])), sum(keep)))
}
