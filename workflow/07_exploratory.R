# Module 07 — declared exploratory and project-specific views

run_exploratory <- function(context) {
  keys <- c("custom_gene_sets", "pathview", "personalized", "wgcna")
  state <- canonical_module_state(context$config, keys)
  if (!identical(state$status, "enabled")) {
    return(record_module_status(context, "07_Exploratory", state$status, state$reason))
  }
  root <- file.path(context$run_dir, "07_Exploratory")
  tables <- file.path(root, "Tables"); figures <- file.path(root, "Figures")
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures, recursive = TRUE, showWarnings = FALSE)

  # Declared target-gene views reuse VST values and selected samples unchanged.
  target_genes <- unique(as.character(context$config$exploratory$target_genes %||% character()))
  target_rows <- list()
  if (length(target_genes)) {
    annotation <- context$state$annotation
    for (gene in target_genes) {
      gene_id <- annotation$gene_id[match(gene, annotation$gene_symbol)]
      if (is.na(gene_id) && gene %in% annotation$gene_id) gene_id <- gene
      if (is.na(gene_id) || !gene_id %in% rownames(context$state$vst)) next
      metadata <- context$state$samples[context$state$selected_samples, , drop = FALSE]
      data <- data.frame(sample = metadata$canonical_sample,
                         group = metadata[[context$state$display_factor]],
                         expression = as.numeric(context$state$vst[gene_id, ]), stringsAsFactors = FALSE)
      write_result_table(data, file.path(tables, paste0(safe_id(gene, "target gene"), "_Expression.tsv")))
      plot <- plot_expression_boxplot(data, paste0(gene, " expression"))
      save_publication_figure(plot, file.path(figures, paste0(safe_id(gene, "target gene"), "_Boxplot")), 95, 95)
      target_rows[[gene]] <- data.frame(target_gene = gene, gene_id = gene_id, status = "complete")
    }
  }
  if (!length(target_rows)) {
    target_rows[[1L]] <- data.frame(target_gene = NA_character_, gene_id = NA_character_,
                                    status = "no_target_genes_declared_or_mapped")
  }
  write_result_table(do.call(rbind, target_rows), file.path(tables, "Target_Gene_Status.tsv"))

  custom_status <- "not_enabled"
  if (module_enabled(context, "custom_gene_sets")) {
    custom_config <- context$config$exploratory$custom_gene_sets
    if (is.null(custom_config$resource_id)) stop("custom_gene_sets requires exploratory.custom_gene_sets.resource_id")
    formal_min_size <- as.integer(custom_config$formal_min_size %||% 3L)
    custom_max_size <- as.integer(custom_config$max_size %||% 500L)
    custom_seed <- as.integer(custom_config$random_seed %||% 123L)
    custom_rank_metric <- as.character(custom_config$rank_metric %||% "wald_statistic")
    custom_curve_mode <- as.character(custom_config$curve_mode %||% "top_reportable")
    custom_reporting_fdr <- as.numeric(custom_config$reporting_fdr %||% 0.25)
    options(bulk_rnaseq.gsea_reporting_fdr = custom_reporting_fdr)
    custom_curve_count <- as.integer(custom_config$curve_count %||% 10L)
    cluster_heatmap_rows <- isTRUE(custom_config$heatmap_cluster_rows)
    cluster_heatmap_columns <- isTRUE(custom_config$heatmap_cluster_columns)
    if (!custom_rank_metric %in% c("wald_statistic", "ashr_lfc")) stop("Unknown custom gene-set rank_metric")
    if (!custom_curve_mode %in% c("all_returned", "top_reportable")) stop("Unknown custom gene-set curve_mode")
    source <- resource_path(context, custom_config$resource_id)
    if (!grepl("\\.xlsx?$", source, ignore.case = TRUE)) stop("Custom gene-set resource must be an Excel workbook")
    raw_sets <- as.data.frame(readxl::read_excel(source, sheet = 1, col_names = FALSE), stringsAsFactors = FALSE)
    if (ncol(raw_sets) < 2L) stop("Custom gene-set workbook requires name and comma-delimited genes")
    raw_sets <- raw_sets[!is.na(raw_sets[[1]]) & !is.na(raw_sets[[2]]) & nzchar(raw_sets[[1]]) & nzchar(raw_sets[[2]]), ]
    sets <- stats::setNames(lapply(raw_sets[[2L]], function(x) {
      genes <- trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1L]])
      unique(genes[nzchar(genes)])
    }), as.character(raw_sets[[1L]]))
    sets <- sets[vapply(sets, length, integer(1)) >= 2L]
    sizes <- vapply(sets, length, integer(1))
    manifest <- data.frame(gene_set = names(sets), declared_size = sizes,
                           analysis_class = ifelse(sizes >= formal_min_size, "gsea_and_heatmap", "expression_heatmap_only"))
    write_result_table(manifest, file.path(tables, "Custom_Gene_Set_Manifest.tsv"))
    custom_root <- file.path(root, "Custom_Gene_Sets")
    dir.create(file.path(custom_root, "Tables"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(custom_root, "Figures", "Curves"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(custom_root, "Figures", "Heatmaps"), recursive = TRUE, showWarnings = FALSE)

    term2gene <- do.call(rbind, lapply(names(sets), function(name) {
      data.frame(term = name, gene = sets[[name]], stringsAsFactors = FALSE)
    }))
    write_result_table(term2gene, file.path(custom_root, "Tables", "Custom_Gene_Sets_TERM2GENE.tsv"))

    # Legacy heatmap matrix: symbol conversion followed by retention of the
    # duplicate symbol row with the highest mean VST expression.
    mapping <- context$state$annotation[rownames(context$state$vst), c("gene_id", "gene_symbol"), drop = FALSE]
    valid <- !is.na(mapping$gene_symbol) & nzchar(mapping$gene_symbol)
    heat_frame <- data.frame(gene_symbol = mapping$gene_symbol[valid],
                             context$state$vst[valid, , drop = FALSE], check.names = FALSE)
    heat_frame$average_expression <- rowMeans(heat_frame[, -1, drop = FALSE])
    heat_frame <- heat_frame[order(-heat_frame$average_expression), , drop = FALSE]
    heat_frame <- heat_frame[!duplicated(heat_frame$gene_symbol), , drop = FALSE]
    rownames(heat_frame) <- heat_frame$gene_symbol
    expression_by_symbol <- as.matrix(heat_frame[, setdiff(names(heat_frame), c("gene_symbol", "average_expression")), drop = FALSE])
    group_annotation <- data.frame(group = context$state$samples[context$state$selected_samples, context$state$display_factor])
    rownames(group_annotation) <- colnames(expression_by_symbol)
    for (name in names(sets)) {
      matched <- rownames(expression_by_symbol)[toupper(rownames(expression_by_symbol)) %in% toupper(sets[[name]])]
      if (length(matched) < 2L) next
      matrix <- expression_by_symbol[matched, , drop = FALSE]
      write_result_table(data.frame(gene_symbol = rownames(matrix), matrix, check.names = FALSE),
                         file.path(custom_root, "Tables", "Heatmap_Data", paste0(safe_id(name), ".tsv.gz")))
      heat <- pheatmap::pheatmap(
        matrix, cluster_rows = cluster_heatmap_rows, cluster_cols = cluster_heatmap_columns,
        annotation_col = group_annotation,
        show_rownames = TRUE, show_colnames = FALSE, scale = "row", border_color = NA,
        color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
        silent = TRUE, main = paste0("Custom set: ", name)
      )$gtable
      save_publication_grob(heat, file.path(custom_root, "Figures", "Heatmaps", paste0("Heatmap_", safe_id(name))),
                            165, max(105, min(260, 4 * length(matched) + 55)))
    }

    for (contrast_id in names(context$state$differential)) {
      de <- context$state$differential[[contrast_id]]$result
      de_input_order <- de[match(rownames(context$state$dds), de$gene_id), , drop = FALSE]
      rank_column <- if (identical(custom_rank_metric, "wald_statistic")) "stat" else "log2FoldChange_ashr"
      ranked <- de_input_order[is.finite(de_input_order[[rank_column]]) &
                                 !is.na(de_input_order$gene_symbol) & nzchar(de_input_order$gene_symbol), ]
      ranked <- ranked[!duplicated(ranked$gene_symbol), , drop = FALSE]
      ranks <- ranked[[rank_column]]; names(ranks) <- ranked$gene_symbol
      ranks <- sort(ranks, decreasing = TRUE)
      set.seed(custom_seed)
      object <- clusterProfiler::GSEA(
        geneList = ranks, TERM2GENE = term2gene, pvalueCutoff = 1,
        pAdjustMethod = "BH", minGSSize = formal_min_size, maxGSSize = custom_max_size,
        BPPARAM = BiocParallel::SerialParam(), seed = TRUE, by = "fgsea", verbose = FALSE
      )
      gsea <- as.data.frame(object)
      if (nrow(gsea)) {
        gsea$pathway <- gsea$ID; gsea$padj <- gsea$p.adjust
        gsea$size <- gsea$setSize; gsea$ES <- gsea$enrichmentScore
        gsea$leadingEdge <- gsea$core_enrichment
        gsea <- gsea[order(gsea$padj, -abs(gsea$NES)), , drop = FALSE]
        write_result_table(gsea, file.path(custom_root, "Tables", paste0(contrast_id, "_GSEA_CustomGeneSets_FullTable.tsv")))
        overview <- plot_gsea_overview(gsea, paste0(contrast_id, " — custom gene-set GSEA"), min(20L, nrow(gsea)))
        if (!is.null(overview)) {
          save_publication_figure(overview, file.path(custom_root, "Figures", paste0(contrast_id, "_Custom_GSEA_Overview")),
                                  165, max(110, 5.2 * min(20L, nrow(gsea)) + 40))
        }
        curve_rows <- if (identical(custom_curve_mode, "all_returned")) gsea else {
          reportable <- gsea[!is.na(gsea$padj) & gsea$padj <= custom_reporting_fdr, , drop = FALSE]
          head(reportable, custom_curve_count)
        }
        for (j in seq_len(nrow(curve_rows))) {
          row <- curve_rows[j, , drop = FALSE]
          curve <- plot_gsea_curve(sets[[row$ID]], ranks, row, contrast_id)
          stem <- paste0(safe_id(contrast_id), "_", safe_id(row$ID))
          save_publication_figure(curve$plot, file.path(custom_root, "Figures", "Curves", stem), 178, 140)
          write_result_table(curve$running, file.path(custom_root, "Tables", "Curve_Data", paste0(stem, "_Running_ES.tsv.gz")))
          write_result_table(curve$hits, file.path(custom_root, "Tables", "Curve_Data", paste0(stem, "_Hits.tsv.gz")))
        }
      }
    }
    custom_status <- sprintf("%d custom sets: GSEA rank=%s, minGSSize=%d plus declared heatmaps",
                             length(sets), custom_rank_metric, formal_min_size)
  }

  if (module_enabled(context, "pathview")) {
    pathview_config <- context$config$exploratory$pathview
    if (is.null(pathview_config$pathway_ids) || !length(pathview_config$pathway_ids)) stop("Enabled pathview requires explicit pathway_ids")
    pathview_root <- file.path(root, "Pathview")
    dir.create(pathview_root, recursive = TRUE, showWarnings = FALSE)
    for (pathway in pathview_config$pathway_ids) {
      resource_ids <- names(resource_registry(context)$entries)
      matches <- resource_ids[grepl(pathway, resource_ids, fixed = TRUE)]
      if (!length(matches)) stop("No frozen Pathview resource is registered for: ", pathway)
      for (resource_id in matches) file.copy(resource_path(context, resource_id), pathview_root, overwrite = FALSE)
    }
  }
  if (module_enabled(context, "personalized") && !isTRUE(context$config$exploratory$personalized$consumes_standard_outputs_only)) {
    stop("Personalized analysis must be declared as consuming standard outputs only")
  }
  record_module_status(context, "07_Exploratory", "complete",
                       sprintf("%d target-gene views; custom sets: %s; configured Pathview/personalized views emitted",
                               length(target_rows), custom_status))
}
