#!/usr/bin/env Rscript

# Prepare downloadable resources once, store under resources/, then reuse offline.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

project_root <- getwd()
resource_dir <- file.path(project_root, "resources")
cache_dir <- file.path(resource_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

count_matrix_path <- Sys.getenv("COUNT_MATRIX_PATH", unset = file.path(project_root, "data", "matrix_gene.count.xls"))
if (!file.exists(count_matrix_path)) {
  stop(paste0("Count matrix not found: ", count_matrix_path))
}

message(">>> resource cache dir: ", cache_dir)
message(">>> count matrix: ", count_matrix_path)

`%||%` <- function(x, y) if (!is.null(x)) x else y

extract_gene_symbols <- function(df) {
  symbol_col <- intersect(c("Symbol", "symbol", "gene_name", "gene"), colnames(df))[1] %||% colnames(df)[1]
  syms <- unique(as.character(df[[symbol_col]]))
  syms <- trimws(syms)
  syms <- syms[!is.na(syms) & syms != "" & syms != "-"]
  unique(syms)
}

prepare_collectri_progeny <- function() {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("Package decoupleR is required for collectri/progeny cache.")
  }

  collectri_file <- file.path(cache_dir, "collectri_human.rds")
  progeny_file <- file.path(cache_dir, "progeny_human_top500.rds")

  if (!file.exists(collectri_file)) {
    message(">>> downloading CollecTRI ...")
    net_collectri <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
    saveRDS(net_collectri, collectri_file)
    message(">>> saved: ", collectri_file, " (rows=", nrow(net_collectri), ")")
  } else {
    message(">>> exists: ", collectri_file)
  }

  if (!file.exists(progeny_file)) {
    message(">>> downloading PROGENy ...")
    net_progeny <- decoupleR::get_progeny(organism = "human", top = 500)
    saveRDS(net_progeny, progeny_file)
    message(">>> saved: ", progeny_file, " (rows=", nrow(net_progeny), ")")
  } else {
    message(">>> exists: ", progeny_file)
  }
}

prepare_stringdb_edges <- function(symbols) {
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("Package STRINGdb is required for stringdb cache.")
  }

  out_dir <- file.path(cache_dir, "stringdb")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  map_file <- file.path(out_dir, "string_map_project.csv")
  edges_file <- file.path(out_dir, "string_edges_project.csv")

  if (file.exists(map_file) && file.exists(edges_file)) {
    message(">>> exists: ", map_file)
    message(">>> exists: ", edges_file)
    return(invisible(NULL))
  }

  max_genes <- suppressWarnings(as.integer(Sys.getenv("STRING_PRELOAD_MAX_GENES", unset = "5000")))
  if (!is.na(max_genes) && max_genes > 0 && length(symbols) > max_genes) {
    symbols <- symbols[seq_len(max_genes)]
    message(">>> STRING preload gene symbols truncated to first ", max_genes)
  }

  message(">>> downloading STRING mappings/interactions ...")
  old_timeout <- getOption("timeout")
  options(timeout = max(1200, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  db <- STRINGdb::STRINGdb$new(
    version = "11.5",
    species = 9606,
    score_threshold = 400,
    input_directory = out_dir
  )

  sym_df <- data.frame(gene = symbols, stringsAsFactors = FALSE)
  mapped <- db$map(sym_df, "gene", removeUnmappedRows = TRUE)
  mapped <- mapped[!duplicated(mapped$STRING_id), , drop = FALSE]
  if (nrow(mapped) < 2) {
    stop("Too few mapped genes from STRINGdb.")
  }

  inter <- db$get_interactions(mapped$STRING_id)
  inter$from_symbol <- mapped$gene[match(inter$from, mapped$STRING_id)]
  inter$to_symbol <- mapped$gene[match(inter$to, mapped$STRING_id)]
  inter <- inter[!is.na(inter$from_symbol) & !is.na(inter$to_symbol), c("from", "to", "from_symbol", "to_symbol", "combined_score")]
  inter <- unique(inter)

  data.table::fwrite(mapped, map_file)
  data.table::fwrite(inter, edges_file)

  message(">>> saved: ", map_file, " (rows=", nrow(mapped), ")")
  message(">>> saved: ", edges_file, " (rows=", nrow(inter), ")")
}

prepare_kegg_cache <- function() {
  if (!requireNamespace("pathview", quietly = TRUE)) {
    message(">>> skip pathview cache: package pathview not installed.")
    return(invisible(NULL))
  }

  kegg_cache_dir <- file.path(cache_dir, "pathview_kegg")
  preplot_dir <- file.path(cache_dir, "pathview_preplot")
  dir.create(kegg_cache_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(preplot_dir, recursive = TRUE, showWarnings = FALSE)

  kegg_pathways <- Sys.getenv("KEGG_PATHWAYS", unset = "hsa04151,hsa04110")
  pathway_ids <- unique(trimws(unlist(strsplit(kegg_pathways, ","))))
  pathway_ids <- pathway_ids[pathway_ids != ""]
  if (length(pathway_ids) == 0) {
    return(invisible(NULL))
  }

  dummy_gene <- c(1, -1)
  names(dummy_gene) <- c("7157", "1956") # TP53 / EGFR (ENTREZ)

  old_wd <- getwd()
  setwd(preplot_dir)
  on.exit(setwd(old_wd), add = TRUE)

  for (pid in pathway_ids) {
    message(">>> preloading KEGG: ", pid)
    tryCatch(
      {
        pathview::pathview(
          gene.data = dummy_gene,
          pathway.id = pid,
          species = "hsa",
          kegg.dir = kegg_cache_dir,
          kegg.native = TRUE,
          out.suffix = "resource_cache"
        )
      },
      error = function(e) {
        message(">>> KEGG preload failed for ", pid, ": ", e$message)
      }
    )
  }
}

expr_raw <- data.table::fread(count_matrix_path, data.table = FALSE)
symbols <- extract_gene_symbols(expr_raw)
message(">>> detected symbols: ", length(symbols))

safe_run <- function(tag, fn) {
  message(">>> [", tag, "] start")
  tryCatch(
    {
      fn()
      message(">>> [", tag, "] done")
    },
    error = function(e) {
      msg <- paste0("[", tag, "] failed: ", e$message)
      message(">>> ", msg)
      writeLines(msg, con = file.path(cache_dir, paste0("SKIPPED_", tag, ".txt")))
    }
  )
}

safe_run("collectri_progeny", prepare_collectri_progeny)
prepare_stringdb <- Sys.getenv("PREPARE_STRINGDB", unset = "0")
if (prepare_stringdb %in% c("1", "true", "TRUE", "yes", "YES")) {
  safe_run("stringdb", function() prepare_stringdb_edges(symbols))
} else {
  msg <- "[stringdb] skipped by default (set PREPARE_STRINGDB=1 to enable full STRING preload)"
  message(">>> ", msg)
  writeLines(msg, con = file.path(cache_dir, "SKIPPED_stringdb.txt"))
}
prepare_kegg <- Sys.getenv("PREPARE_KEGG", unset = "0")
if (prepare_kegg %in% c("1", "true", "TRUE", "yes", "YES")) {
  safe_run("kegg", prepare_kegg_cache)
} else {
  msg <- "[kegg] skipped by default to keep resource free of image artifacts (set PREPARE_KEGG=1 to enable)"
  message(">>> ", msg)
  writeLines(msg, con = file.path(cache_dir, "SKIPPED_kegg.txt"))
}

message(">>> resource prepare finished. cache root: ", cache_dir)

