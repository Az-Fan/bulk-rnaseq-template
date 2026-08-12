# Module 07 — declared exploratory views

run_exploratory <- function(context) {
  keys <- c("custom_gene_sets", "pathview", "personalized", "wgcna")
  state <- canonical_module_state(context$config, keys)
  if (!identical(state$status, "enabled")) {
    return(record_module_status(context, "07_Exploratory", state$status, state$reason))
  }
  root <- file.path(context$run_dir, "07_Exploratory")
  tables <- file.path(root, "Tables"); figures <- file.path(root, "Figures")
  dir.create(tables, recursive = TRUE, showWarnings = FALSE); dir.create(figures, recursive = TRUE, showWarnings = FALSE)
  target_genes <- unique(as.character(context$config$exploratory$target_genes %||% character()))
  rows <- list()
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
      rows[[gene]] <- data.frame(target_gene = gene, gene_id = gene_id, status = "complete")
    }
  }
  if (!length(rows)) rows[[1L]] <- data.frame(target_gene = NA_character_, gene_id = NA_character_, status = "no_target_genes_declared_or_mapped")
  write_result_table(do.call(rbind, rows), file.path(tables, "Target_Gene_Status.tsv"))

  custom_status <- "not_enabled"
  if (module_enabled(context, "custom_gene_sets")) {
    custom_config <- context$config$exploratory$custom_gene_sets
    if (is.null(custom_config$resource_id)) stop("custom_gene_sets requires exploratory.custom_gene_sets.resource_id")
    source <- resource_path(context, custom_config$resource_id)
    if (!grepl("\\.xlsx?$", source, ignore.case = TRUE)) stop("Custom gene-set resource must be an Excel workbook or a supported GMT resource")
    raw_sets <- as.data.frame(readxl::read_excel(source, sheet = 1, col_names = FALSE), stringsAsFactors = FALSE)
    if (ncol(raw_sets) < 2L) stop("Custom gene-set workbook requires name and comma-delimited genes")
    sets <- stats::setNames(lapply(raw_sets[[2L]], function(x) unique(trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1L]]))),
                            as.character(raw_sets[[1L]]))
    minimum <- as.integer(custom_config$formal_min_size %||% 10L)
    sizes <- vapply(sets, length, integer(1))
    manifest <- data.frame(gene_set = names(sets), declared_size = sizes,
                           analysis_class = ifelse(sizes >= minimum, "formal_gsea", ifelse(sizes >= 3L, "small_set_exploratory", "expression_only")))
    write_result_table(manifest, file.path(tables, "Custom_Gene_Set_Manifest.tsv"))
    custom_root <- file.path(root, "Custom_Gene_Sets")
    dir.create(file.path(custom_root, "Tables"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(custom_root, "Figures"), recursive = TRUE, showWarnings = FALSE)
    for (contrast_id in names(context$state$differential)) {
      de <- context$state$differential[[contrast_id]]$result
      ranked <- de[is.finite(de$stat) & !is.na(de$gene_symbol) & nzchar(de$gene_symbol), ]
      ranked <- ranked[order(-ranked$stat), ]; ranked <- ranked[!duplicated(ranked$gene_symbol), ]
      ranks <- ranked$stat; names(ranks) <- ranked$gene_symbol
      formal <- sets[sizes >= minimum]
      if (length(formal)) {
        gsea <- as.data.frame(fgsea::fgseaMultilevel(formal, ranks, minSize = minimum, maxSize = 500L, eps = 0))
        if (nrow(gsea)) {
          gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = "/", character(1))
          gsea <- gsea[order(gsea$padj, -abs(gsea$NES)), ]
          write_result_table(gsea, file.path(custom_root, "Tables", paste0(contrast_id, "_Custom_GSEA.tsv")))
          overview <- plot_gsea_overview(gsea, paste0(contrast_id, " — custom gene-set GSEA"), min(20L, nrow(gsea)))
          save_publication_figure(overview, file.path(custom_root, "Figures", paste0(contrast_id, "_Custom_GSEA_Overview")),
                                  165, max(110, 5.2 * min(20L, nrow(gsea)) + 40))
          curve_rows <- head(gsea, min(6L, nrow(gsea)))
          for (j in seq_len(nrow(curve_rows))) {
            row <- curve_rows[j, , drop = FALSE]
            curve <- plot_gsea_curve(formal[[row$pathway]], ranks, row, contrast_id)
            save_publication_figure(curve$plot, file.path(custom_root, "Figures", "Curves",
                                                          paste0(safe_id(contrast_id), "_", safe_id(row$pathway))), 178, 140)
          }
        }
      }
      # Small sets are retained as labelled exploratory z-score heatmaps rather
      # than being misrepresented as pathway GSEA.
      for (name in names(sets)[sizes >= 3L & sizes < minimum]) {
        ids <- context$state$annotation$gene_id[match(sets[[name]], context$state$annotation$gene_symbol)]
        ids <- ids[!is.na(ids) & ids %in% rownames(context$state$vst)]
        if (length(ids) < 2L) next
        matrix <- t(scale(t(context$state$vst[ids, , drop = FALSE])))
        rownames(matrix) <- context$state$annotation$gene_symbol[match(ids, context$state$annotation$gene_id)]
        heat <- pheatmap::pheatmap(matrix, cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA,
                                  color = grDevices::colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(101),
                                  silent = TRUE, main = paste0(name, " — small-set exploratory expression"))$gtable
        save_publication_grob(heat, file.path(custom_root, "Figures", "Small_Sets", safe_id(name)), 145, max(80, 5 * length(ids) + 35))
      }
    }
    custom_status <- sprintf("%d formal; %d small-set exploratory; %d expression-only", sum(sizes >= minimum),
                             sum(sizes >= 3L & sizes < minimum), sum(sizes < 3L))
  }
  if (module_enabled(context, "pathview")) {
    pathview_config <- context$config$exploratory$pathview
    if (is.null(pathview_config$pathway_ids) || !length(pathview_config$pathway_ids)) stop("Enabled pathview requires explicit pathway_ids")
    pathview_root <- file.path(root, "Pathview"); dir.create(pathview_root, recursive = TRUE, showWarnings = FALSE)
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
                       sprintf("%d target-gene views; custom sets: %s; configured Pathview/personalized views emitted", length(rows), custom_status))
}
