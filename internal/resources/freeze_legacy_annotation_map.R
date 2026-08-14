#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: freeze_legacy_annotation_map.R LEGACY_DE_CSV OUTPUT_RDS")
input <- normalizePath(args[[1]], mustWork = TRUE)
output <- args[[2]]
x <- utils::read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("gene_id", "gene", "ENTREZID")
if (length(setdiff(required, names(x)))) stop("Legacy DE table lacks annotation columns")
mapping <- unique(x[, required, drop = FALSE])
names(mapping) <- c("gene_id", "gene_symbol", "entrez_id")
mapping$gene_id <- sub("\\.[0-9]+$", "", as.character(mapping$gene_id))
mapping$gene_symbol <- as.character(mapping$gene_symbol)
mapping$entrez_id <- as.character(mapping$entrez_id)
if (anyDuplicated(mapping$gene_id)) stop("Legacy annotation map has duplicate gene IDs")
resource_id <- sub("\\.rds$", "", basename(output), ignore.case = TRUE)
payload <- list(
  resource_id = resource_id,
  source_table_sha256 = digest::digest(file = input, algo = "sha256"),
  mapping = mapping
)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
if (file.exists(output)) {
  if (!identical(readRDS(output), payload)) stop("Refusing to overwrite a different frozen mapping")
} else saveRDS(payload, output, compress = "xz")
cat(sprintf("rows=%d sha256=%s\n", nrow(mapping), digest::digest(file = output, algo = "sha256")))
