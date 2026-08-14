#!/usr/bin/env Rscript

# Convert the archived msigdbr 7.5.1 package data into a standalone frozen
# resource. The package is never installed or loaded during formal analysis.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: freeze_msigdbr_7_5_1.R SOURCE_TARBALL OUTPUT_RDS")
}

source_tar <- normalizePath(args[[1]], mustWork = TRUE)
output_rds <- args[[2]]
dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)

scratch <- tempfile("msigdbr-7.5.1-")
dir.create(scratch)
on.exit(unlink(scratch, recursive = TRUE, force = TRUE), add = TRUE)
utils::untar(source_tar, files = "msigdbr/R/sysdata.rda", exdir = scratch)

data_env <- new.env(parent = emptyenv())
load(file.path(scratch, "msigdbr", "R", "sysdata.rda"), envir = data_env)
required <- c("msigdbr_genesets", "msigdbr_geneset_genes", "msigdbr_genes")
missing <- setdiff(required, ls(data_env, all.names = TRUE))
if (length(missing)) stop("Archived resource lacks objects: ", paste(missing, collapse = ", "))

flat <- merge(
  data_env$msigdbr_genesets,
  data_env$msigdbr_geneset_genes,
  by = "gs_id",
  all = FALSE,
  sort = FALSE
)
flat <- merge(flat, data_env$msigdbr_genes, by = "gene_id", all = FALSE, sort = FALSE)
flat$gene_symbol <- flat$human_gene_symbol
flat$entrez_gene <- flat$human_entrez_gene
flat$ensembl_gene <- flat$human_ensembl_gene

front <- c(
  "gs_cat", "gs_subcat", "gs_name", "gene_symbol", "entrez_gene",
  "ensembl_gene", "human_gene_symbol", "human_entrez_gene",
  "human_ensembl_gene"
)
flat <- flat[, c(front, setdiff(names(flat), front)), drop = FALSE]
flat <- unique(flat)
flat <- flat[order(flat$gs_name, flat$human_gene_symbol, flat$gene_symbol), , drop = FALSE]
rownames(flat) <- NULL

payload <- list(
  resource_id = "msigdb_7.5.1_human_gene_sets",
  species = "Homo sapiens",
  msigdb_release = "7.5.1",
  source_archive_sha256 = digest::digest(file = source_tar, algo = "sha256"),
  generated_by = "internal/resources/freeze_msigdbr_7_5_1.R",
  gene_sets = flat
)

if (file.exists(output_rds)) {
  existing <- readRDS(output_rds)
  if (!identical(existing$source_archive_sha256, payload$source_archive_sha256) ||
      !identical(existing$gene_sets, payload$gene_sets)) {
    stop("Refusing to overwrite a different frozen MSigDB resource: ", output_rds)
  }
} else {
  saveRDS(payload, output_rds, compress = "xz")
}

cat(sprintf(
  "resource=%s gene_sets=%d memberships=%d sha256=%s\n",
  output_rds,
  length(unique(flat$gs_name)),
  nrow(flat),
  digest::digest(file = output_rds, algo = "sha256")
))
