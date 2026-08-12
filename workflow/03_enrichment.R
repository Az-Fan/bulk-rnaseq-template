# Module 03 — over-representation and ranked-list enrichment
#
# ORA uses tested genes as its universe and is run separately for Up and Down
# DEG sets. GSEA uses every finite DESeq2 Wald statistic. All gene-set resources
# are frozen registry entries; no database is queried during analysis.

run_enrichment <- function(context) {
  decision <- module_decision(context$config, "enrichment")
  if (!identical(decision$status, "enabled")) {
    return(record_module_status(context, "03_Enrichment", decision$status, decision$reason %||% ""))
  }
  if (is.null(context$state$differential)) stop("Differential expression must complete before enrichment")

  resource_id <- context$config$resources$enrichment_gene_sets %||%
    if (identical(context$config$species, "human")) "msigdb.2025.1.Hs.rds" else NULL
  if (is.null(resource_id)) stop("An enabled enrichment module requires resources.enrichment_gene_sets")
  resource_file <- resource_path(context, resource_id)
  gene_sets <- readRDS(resource_file)
  required <- c("db_gene_symbol", "gs_name", "gs_collection", "gs_subcollection")
  if (length(setdiff(required, names(gene_sets)))) stop("Enrichment resource lacks required columns")
  target_code <- if (identical(context$config$species, "human")) "HS" else "MM"
  if ("db_target_species" %in% names(gene_sets)) gene_sets <- gene_sets[gene_sets$db_target_species == target_code, , drop = FALSE]
  collections <- context$config$enrichment$collections %||% c("H", "C2", "C5")
  gene_sets <- gene_sets[gene_sets$gs_collection %in% collections, , drop = FALSE]
  gene_sets <- gene_sets[!is.na(gene_sets$db_gene_symbol) & nzchar(gene_sets$db_gene_symbol), , drop = FALSE]
  gene_sets$database <- ifelse(gene_sets$gs_collection == "H", "Hallmark",
                               ifelse(nzchar(gene_sets$gs_subcollection), gene_sets$gs_subcollection, gene_sets$gs_collection))
  gene_sets <- unique(gene_sets[, c("gs_name", "db_gene_symbol", "database")])
  if (!nrow(gene_sets)) stop("No compatible gene sets remain after species/collection filtering")

  min_size <- as.integer(context$config$enrichment$min_size %||% 15L)
  max_size <- as.integer(context$config$enrichment$max_size %||% 500L)
  if (min_size < 10L) stop("Formal GSEA min_size must be at least 10; use exploratory small-set views instead")
  gsea_curves <- as.integer(context$config$enrichment$gsea_curves %||% 6L)
  profiles <- threshold_profiles(context$config)

  ora_one <- function(query, universe, term_table) {
    term_list <- split(term_table$db_gene_symbol, term_table$gs_name)
    database <- vapply(split(term_table$database, term_table$gs_name), function(x) x[[1]], character(1))
    rows <- lapply(names(term_list), function(term) {
      genes <- intersect(unique(term_list[[term]]), universe)
      overlap <- intersect(query, genes)
      query_size <- length(query); universe_size <- length(universe); term_size <- length(genes)
      if (!query_size || term_size < min_size || term_size > max_size || !length(overlap)) return(NULL)
      pvalue <- stats::phyper(length(overlap) - 1L, term_size, universe_size - term_size,
                             query_size, lower.tail = FALSE)
      data.frame(pathway = term, database = database[[term]], overlap_count = length(overlap),
                 query_size = query_size, term_size = term_size, universe_size = universe_size,
                 gene_ratio = length(overlap) / query_size, background_ratio = term_size / universe_size,
                 fold_enrichment = (length(overlap) / query_size) / (term_size / universe_size),
                 pvalue = pvalue, overlap_genes = paste(sort(overlap), collapse = "/"), stringsAsFactors = FALSE)
    })
    result <- do.call(rbind, rows)
    if (is.null(result)) return(data.frame())
    result$padj <- stats::p.adjust(result$pvalue, method = "BH")
    result[order(result$padj, -result$fold_enrichment), ]
  }

  results <- list()
  for (contrast_id in names(context$state$differential)) {
    differential <- context$state$differential[[contrast_id]]
    de <- differential$result
    root <- comparison_module_dir(context, contrast_id, "03_Enrichment")
    ora_root <- file.path(root, "ORA")
    gsea_root <- file.path(root, "GSEA")
    dir.create(file.path(ora_root, "Tables"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(ora_root, "Figures"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(gsea_root, "Tables"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(gsea_root, "Figures"), recursive = TRUE, showWarnings = FALSE)

    symbol_by_id <- stats::setNames(de$gene_symbol, de$gene_id)
    universe <- unique(symbol_by_id[!is.na(symbol_by_id) & nzchar(symbol_by_id)])
    profile_results <- list()
    for (profile_name in names(profiles)) {
      direction_column <- paste0("direction_", profile_name)
      ora_rows <- list()
      for (direction in c("Up", "Down")) {
        ids <- de$gene_id[de[[direction_column]] == direction]
        query <- unique(symbol_by_id[ids])
        query <- query[!is.na(query) & nzchar(query)]
        ora <- ora_one(query, universe, gene_sets)
        if (nrow(ora)) ora$direction <- direction
        ora_rows[[direction]] <- ora
      }
      ora <- do.call(rbind, ora_rows)
      if (is.null(ora)) ora <- data.frame()
      profile_results[[profile_name]] <- ora
      write_result_table(ora, file.path(ora_root, "Tables", paste0("ORA_", profile_name, "_Full.tsv")))
      if (nrow(ora)) {
        write_result_table(ora[!is.na(ora$padj) & ora$padj < profiles[[profile_name]]$padj, ],
                           file.path(ora_root, "Tables", paste0("ORA_", profile_name, "_Significant.tsv")))
        plot <- plot_ora_publication(ora, paste0(contrast_id, " — ", profile_name, " ORA"),
                                     top_n = as.integer(context$config$enrichment$ora_overview_terms %||% 12L))
        save_publication_figure(plot, file.path(ora_root, "Figures", paste0("ORA_", profile_name, "_Directional")),
                                183, max(120, 5.2 * min(24L, nrow(ora)) + 45))
      }
    }

    ranked <- de[is.finite(de$stat) & !is.na(de$gene_symbol) & nzchar(de$gene_symbol), c("gene_symbol", "stat")]
    ranked <- ranked[order(-ranked$stat, ranked$gene_symbol), ]
    ranked <- ranked[!duplicated(ranked$gene_symbol), ]
    ranks <- ranked$stat; names(ranks) <- ranked$gene_symbol
    write_result_table(data.frame(gene_symbol = names(ranks), wald_statistic = unname(ranks)),
                       file.path(gsea_root, "Tables", "GSEA_Ranked_Gene_List.tsv.gz"))
    pathway_list <- split(gene_sets$db_gene_symbol, gene_sets$gs_name)
    pathway_list <- lapply(pathway_list, unique)
    set.seed(as.integer(context$config$analysis$random_seed %||% 104729L))
    gsea <- fgsea::fgseaMultilevel(pathways = pathway_list, stats = ranks,
                                  minSize = min_size, maxSize = max_size,
                                  eps = 0, scoreType = "std")
    gsea <- as.data.frame(gsea)
    if (nrow(gsea)) {
      gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = "/", character(1))
      gsea$database <- vapply(gsea$pathway, function(x) gene_sets$database[match(x, gene_sets$gs_name)], character(1))
      gsea <- gsea[order(gsea$padj, -abs(gsea$NES)), ]
    }
    write_result_table(gsea, file.path(gsea_root, "Tables", "GSEA_Full.tsv"))
    if (nrow(gsea)) {
      overview <- plot_gsea_overview(gsea, paste0(contrast_id, " — GSEA overview"),
                                     top_n = as.integer(context$config$enrichment$overview_terms %||% 20L))
      save_publication_figure(overview, file.path(gsea_root, "Figures", "GSEA_Directional_Overview"), 178,
                              max(115, 5.3 * min(20L, nrow(gsea)) + 45))
      curve_rows <- head(gsea, min(gsea_curves, nrow(gsea)))
      curve_manifest <- list()
      for (j in seq_len(nrow(curve_rows))) {
        summary_row <- curve_rows[j, , drop = FALSE]
        curve <- plot_gsea_curve(pathway_list[[summary_row$pathway]], ranks, summary_row, contrast_id)
        stem <- sprintf("%02d_%s", j, safe_id(gsub("[^A-Za-z0-9_.-]+", "_", summary_row$pathway), "pathway"))
        save_publication_figure(curve$plot, file.path(gsea_root, "Figures", "Curves", stem), 178, 140)
        write_result_table(curve$running, file.path(gsea_root, "Tables", "Curve_Data", paste0(stem, "_Running_ES.tsv.gz")))
        write_result_table(curve$hits, file.path(gsea_root, "Tables", "Curve_Data", paste0(stem, "_Hits.tsv.gz")))
        curve_manifest[[j]] <- data.frame(pathway = summary_row$pathway, NES = summary_row$NES,
                                          padj = summary_row$padj, size = summary_row$size,
                                          figure_stem = stem, stringsAsFactors = FALSE)
      }
      write_result_table(do.call(rbind, curve_manifest), file.path(gsea_root, "Tables", "GSEA_Curve_Manifest.tsv"))
    }
    results[[contrast_id]] <- list(ora = profile_results, gsea = gsea, ranks = ranks,
                                   pathway_list = pathway_list, directory = root)
  }
  context$state$enrichment <- results
  record_module_status(context, "03_Enrichment", "complete",
                       sprintf("ORA and GSEA completed for %d contrasts using %s", length(results), resource_id))
}
