# Module 02 — differential expression
#
# The confirmed design is fitted once. Every declared contrast reuses that
# model. Screening and formal effect-size tests are explicit, distinct modes.

run_differential <- function(context) {
  decision <- module_decision(context$config, "differential")
  if (!identical(decision$status, "enabled")) {
    return(record_module_status(context, "02_Differential", decision$status, decision$reason %||% ""))
  }
  if (is.null(context$state$counts)) stop("QC must complete before differential expression")

  metadata <- context$state$samples[context$state$selected_samples, , drop = FALSE]
  rownames(metadata) <- metadata$sample
  design <- stats::as.formula(context$config$design$formula)
  design_variables <- all.vars(design)
  missing_variables <- setdiff(design_variables, names(metadata))
  if (length(missing_variables)) stop("Design formula names missing samples.tsv columns: ", paste(missing_variables, collapse = ", "))
  model_matrix <- stats::model.matrix(design, metadata)
  if (nrow(model_matrix) <= ncol(model_matrix)) {
    stop("The design has no residual degrees of freedom; more biological replication or a simpler confirmed design is required")
  }
  if (qr(model_matrix)$rank < ncol(model_matrix)) stop("Confirmed design matrix is not full rank")
  write_result_table(data.frame(sample = rownames(model_matrix), model_matrix, check.names = FALSE),
                     file.path(context$run_dir, "02_Differential", "Tables", "Design_Matrix.tsv"))

  for (i in seq_len(nrow(context$contrasts))) {
    contrast <- context$contrasts[i, , drop = FALSE]
    factor_name <- as.character(contrast$factor[[1]])
    if (!factor_name %in% design_variables) stop("Contrast factor is absent from design formula: ", factor_name)
    observed <- unique(as.character(metadata[[factor_name]]))
    if (!all(c(contrast$numerator[[1]], contrast$denominator[[1]]) %in% observed)) {
      stop("Contrast levels are absent from approved samples: ", contrast$contrast_id[[1]])
    }
    level_counts <- table(as.character(metadata[[factor_name]]))
    if (any(level_counts[c(as.character(contrast$numerator[[1]]), as.character(contrast$denominator[[1]]))] < 2L)) {
      stop("Each contrasted level requires at least two approved biological samples: ", contrast$contrast_id[[1]])
    }
  }

  counts <- context$state$counts[context$state$keep_final, context$state$selected_samples, drop = FALSE]
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, metadata, design = design)
  dds <- DESeq2::DESeq(dds, quiet = TRUE, parallel = FALSE)
  normalized <- DESeq2::counts(dds, normalized = TRUE)
  vst <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = FALSE))
  diagnostic_root <- file.path(context$run_dir, "02_Differential", "Diagnostics")
  dir.create(file.path(diagnostic_root, "Figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(diagnostic_root, "Tables"), recursive = TRUE, showWarnings = FALSE)
  size_factor_data <- data.frame(sample = metadata$canonical_sample, size_factor = DESeq2::sizeFactors(dds),
                                 group = metadata[[context$config$design$display_factor]])
  write_result_table(size_factor_data, file.path(diagnostic_root, "Tables", "Size_Factors.tsv"))
  size_factor_plot <- ggplot2::ggplot(size_factor_data,
                                      ggplot2::aes(stats::reorder(sample, size_factor), size_factor, fill = group)) +
    ggplot2::geom_col(width = 0.72) + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = publication_palette(unique(size_factor_data$group))) +
    ggplot2::labs(title = "DESeq2 size factors", x = NULL, y = "Median-ratio size factor", fill = NULL) +
    theme_publication()
  save_publication_figure(size_factor_plot, file.path(diagnostic_root, "Figures", "Size_Factors"), 145, 105)
  dispersion_data <- data.frame(
    gene_id = rownames(dds), baseMean = SummarizedExperiment::rowData(dds)$baseMean,
    gene_estimate = SummarizedExperiment::rowData(dds)$dispGeneEst,
    fitted = SummarizedExperiment::rowData(dds)$dispFit,
    final = DESeq2::dispersions(dds), stringsAsFactors = FALSE
  )
  write_result_table(dispersion_data, file.path(diagnostic_root, "Tables", "Dispersion_Estimates.tsv.gz"))
  dispersion_plot_data <- dispersion_data[is.finite(dispersion_data$baseMean) & dispersion_data$baseMean > 0 &
                                            is.finite(dispersion_data$final) & dispersion_data$final > 0, ]
  dispersion_plot_data <- dispersion_plot_data[order(dispersion_plot_data$baseMean), ]
  dispersion_plot <- ggplot2::ggplot(dispersion_plot_data, ggplot2::aes(baseMean, final)) +
    ggplot2::geom_point(alpha = 0.28, size = 0.55, colour = "#4C78A8") +
    ggplot2::geom_line(ggplot2::aes(y = fitted), colour = "#D55E00", linewidth = 0.7, na.rm = TRUE) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::labs(title = "Dispersion estimates and fitted trend", x = "Mean normalized count", y = "Dispersion") +
    theme_publication()
  save_publication_figure(dispersion_plot, file.path(diagnostic_root, "Figures", "Dispersion_Trend"), 145, 110)
  saveRDS(dds, file.path(context$run_dir, "Provenance", "deseq2_model.rds"), compress = "xz")
  write_result_table(data.frame(gene_id = rownames(normalized), normalized, check.names = FALSE),
                     file.path(context$run_dir, "02_Differential", "Tables", "Normalized_Counts.tsv.gz"))

  profiles <- threshold_profiles(context$config)
  differential_config <- context$config$differential %||% list()
  primary_name <- differential_config$primary_profile %||% names(profiles)[[1L]]
  if (!primary_name %in% names(profiles)) stop("Unknown differential.primary_profile: ", primary_name)
  mode <- context$config$thresholds$mode %||% "screening"
  if (!mode %in% c("screening", "lfc_threshold_test")) stop("thresholds.mode must be screening or lfc_threshold_test")
  results_by_contrast <- list()
  for (i in seq_len(nrow(context$contrasts))) {
    contrast <- context$contrasts[i, , drop = FALSE]
    contrast_id <- as.character(contrast$contrast_id[[1]])
    contrast_vector <- c(as.character(contrast$factor[[1]]), as.character(contrast$numerator[[1]]),
                         as.character(contrast$denominator[[1]]))
    primary <- profiles[[primary_name]]
    if (identical(mode, "lfc_threshold_test")) {
      raw_result <- DESeq2::results(dds, contrast = contrast_vector, alpha = primary$padj,
                                    lfcThreshold = primary$abs_lfc, altHypothesis = "greaterAbs")
    } else {
      raw_result <- DESeq2::results(dds, contrast = contrast_vector, alpha = primary$padj)
    }
    shrunken <- DESeq2::lfcShrink(dds, contrast = contrast_vector, res = raw_result, type = "ashr")
    raw_df <- as.data.frame(raw_result)
    shrunken_df <- as.data.frame(shrunken)
    ids <- rownames(raw_df)
    result <- data.frame(
      gene_id = ids, gene_symbol = context$state$annotation[ids, "gene_symbol"],
      entrez_id = context$state$annotation[ids, "entrez_id"],
      baseMean = raw_df$baseMean, log2FoldChange_raw = raw_df$log2FoldChange,
      log2FoldChange_ashr = shrunken_df$log2FoldChange, lfcSE = shrunken_df$lfcSE,
      stat = raw_df$stat, pvalue = raw_df$pvalue, padj = raw_df$padj,
      stringsAsFactors = FALSE
    )
    result$gene_label <- ifelse(is.na(result$gene_symbol) | !nzchar(result$gene_symbol), result$gene_id, result$gene_symbol)
    for (profile_name in names(profiles)) {
      profile <- profiles[[profile_name]]
      if (identical(mode, "lfc_threshold_test") && !identical(profile_name, primary_name)) {
        stop("lfc_threshold_test currently supports exactly one declared threshold profile")
      }
      significant <- !is.na(result$padj) & result$padj < profile$padj
      if (identical(mode, "screening")) significant <- significant & abs(result$log2FoldChange_ashr) > profile$abs_lfc
      direction <- ifelse(significant & result$log2FoldChange_ashr > 0, "Up",
                          ifelse(significant & result$log2FoldChange_ashr < 0, "Down", "Not_significant"))
      result[[paste0("significant_", profile_name)]] <- significant
      result[[paste0("direction_", profile_name)]] <- direction
    }
    result$significant <- result[[paste0("significant_", primary_name)]]
    result$direction <- result[[paste0("direction_", primary_name)]]
    result <- result[order(result$padj, -abs(result$log2FoldChange_ashr), na.last = TRUE), ]

    root <- comparison_module_dir(context, contrast_id, "02_Differential")
    figures <- file.path(root, "Figures"); tables <- file.path(root, "Tables")
    dir.create(figures, recursive = TRUE, showWarnings = FALSE)
    dir.create(tables, recursive = TRUE, showWarnings = FALSE)
    write_result_table(result, file.path(tables, "Differential_Full.tsv"))
    pvalue_data <- result[is.finite(result$pvalue), c("gene_id", "pvalue", "padj")]
    write_result_table(data.frame(total_tested = nrow(result), finite_pvalues = nrow(pvalue_data),
                                  missing_pvalues = sum(is.na(result$pvalue)), missing_padj = sum(is.na(result$padj))),
                       file.path(tables, "PValue_Diagnostic_Summary.tsv"))
    pvalue_plot <- ggplot2::ggplot(pvalue_data, ggplot2::aes(pvalue)) +
      ggplot2::geom_histogram(binwidth = 0.05, boundary = 0, closed = "left",
                              fill = "#4C78A8", colour = "white", linewidth = 0.25) +
      ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      ggplot2::labs(title = paste0(contrast_id, " — raw P-value distribution"), x = "Wald-test P value", y = "Genes") +
      theme_publication()
    save_publication_figure(pvalue_plot, file.path(figures, "PValue_Histogram"), 125, 95)
    for (profile_name in names(profiles)) {
      significant_column <- paste0("significant_", profile_name)
      write_result_table(result[result[[significant_column]], ],
                         file.path(tables, paste0("Differential_", profile_name, "_Significant.tsv")))
    }

    ma_data <- data.frame(baseMean = raw_df$baseMean, Raw = raw_df$log2FoldChange,
                          `ashr-shrunken` = shrunken_df$log2FoldChange,
                          direction = result$direction[match(ids, result$gene_id)], check.names = FALSE)
    ma_long <- rbind(data.frame(baseMean = ma_data$baseMean, log2FoldChange = ma_data$Raw, estimate = "Raw", direction = ma_data$direction),
                     data.frame(baseMean = ma_data$baseMean, log2FoldChange = ma_data[["ashr-shrunken"]], estimate = "ashr-shrunken", direction = ma_data$direction))
    ma_plot <- ggplot2::ggplot(ma_long, ggplot2::aes(baseMean, log2FoldChange, colour = direction)) +
      ggplot2::geom_point(alpha = 0.48, size = 0.65) + ggplot2::scale_x_log10() +
      ggplot2::scale_colour_manual(values = direction_palette, drop = FALSE) +
      ggplot2::geom_hline(yintercept = 0, colour = "#343A40", linewidth = 0.4) +
      ggplot2::facet_wrap(~estimate, nrow = 1) +
      ggplot2::labs(title = paste0(contrast_id, " — MA plots"), subtitle = "Raw and ashr-shrunken estimates from the same fitted model",
                    x = "Mean normalized count", y = "log2 fold change", colour = NULL) + theme_publication() +
      ggplot2::theme(legend.position = "bottom")
    save_publication_figure(ma_plot, file.path(figures, "MA_Raw_and_Shrunk"), 178, 105)

    volcano <- plot_volcano_publication(
      result, contrast_id, primary$padj, primary$abs_lfc,
      label_count = as.integer(differential_config$volcano_label_count %||% 15L)
    )
    # Preserve the user's established single-panel volcano proportions.
    save_publication_figure(volcano, file.path(figures, "Volcano"), 178, 132)
    write_result_table(attr(volcano, "plot_data"), file.path(tables, "Volcano_Plot_Data.tsv.gz"))

    heatmap_per_direction <- as.integer(differential_config$heatmap_genes_per_direction %||% 25L)
    selected_genes <- unique(c(head(result$gene_id[result$direction == "Up"], heatmap_per_direction),
                               head(result$gene_id[result$direction == "Down"], heatmap_per_direction)))
    if (length(selected_genes) >= 2L) {
      matrix <- t(scale(t(vst[selected_genes, , drop = FALSE])))
      symbols <- context$state$annotation[selected_genes, "gene_symbol"]
      rownames(matrix) <- make.unique(ifelse(is.na(symbols) | !nzchar(symbols), selected_genes, symbols))
      colnames(matrix) <- metadata$canonical_sample
      annotation_col <- data.frame(group = metadata[[contrast$factor[[1]]]])
      rownames(annotation_col) <- metadata$canonical_sample
      heat <- pheatmap::pheatmap(matrix, annotation_col = annotation_col, border_color = NA,
                                color = grDevices::colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(101),
                                show_rownames = TRUE, fontsize_row = 6.5, silent = TRUE,
                                main = paste0(contrast_id, " — top DE genes"))$gtable
      save_publication_grob(heat, file.path(figures, "Top_DEG_Heatmap"), 178, max(125, 3.2 * nrow(matrix)))
      write_result_table(data.frame(gene = rownames(matrix), matrix, check.names = FALSE),
                         file.path(tables, "Top_DEG_Heatmap_Matrix.tsv.gz"))
    }

    cooks <- SummarizedExperiment::assays(dds)[["cooks"]]
    if (!is.null(cooks)) {
      cooks_table <- data.frame(gene_id = rownames(cooks), max_cooks = apply(cooks, 1, max, na.rm = TRUE),
                                sample_at_max = colnames(cooks)[max.col(cooks, ties.method = "first")])
      write_result_table(cooks_table[order(-cooks_table$max_cooks), ], file.path(tables, "Cooks_Distance_Gene_Diagnostics.tsv.gz"))
    }
    results_by_contrast[[contrast_id]] <- list(result = result, contrast = contrast, directory = root)
  }

  context$state$dds <- dds
  context$state$vst <- vst
  context$state$model_matrix <- model_matrix
  context$state$differential <- results_by_contrast
  saveRDS(list(vst = vst, metadata = metadata, annotation = context$state$annotation,
               tested_gene_ids = rownames(dds), contrasts = results_by_contrast,
               threshold_mode = mode, profiles = profiles),
          file.path(context$run_dir, "Provenance", "differential_state.rds"), compress = "xz")
  record_module_status(context, "02_Differential", "complete",
                       sprintf("One model fitted; %d confirmed contrasts evaluated", length(results_by_contrast)))
}
