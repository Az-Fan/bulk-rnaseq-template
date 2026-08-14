# Module 01 — count and sample quality control
#
# This module validates the frozen matrix, records filtering and sample metrics,
# and computes initial/final sample-level views. It never chooses exclusions.

run_qc <- function(context) {
  decision <- module_decision(context$config, "qc")
  if (!identical(decision$status, "enabled")) {
    return(record_module_status(context, "01_QC", decision$status, decision$reason %||% ""))
  }

  raw <- read_project_table(context$paths$counts)
  sample_columns <- as.character(context$samples$sample)
  missing_samples <- setdiff(sample_columns, names(raw))
  if (length(missing_samples)) stop("Count matrix is missing declared samples: ", paste(missing_samples, collapse = ", "))
  gene_id_column <- context$config$inputs$gene_id_column %||%
    if ("Gene_ID" %in% names(raw)) "Gene_ID" else names(raw)[[1L]]
  if (!gene_id_column %in% names(raw)) stop("Configured gene_id_column is absent: ", gene_id_column)
  gene_symbol_column <- context$config$inputs$gene_symbol_column %||% if ("Name" %in% names(raw)) "Name" else NULL
  entrez_column <- context$config$inputs$entrez_id_column %||%
    if ("Entrez_geneID" %in% names(raw)) "Entrez_geneID" else NULL

  count_frame <- raw[, sample_columns, drop = FALSE]
  count_frame[] <- lapply(count_frame, function(x) suppressWarnings(as.numeric(x)))
  count_matrix <- as.matrix(count_frame)
  if (any(!is.finite(count_matrix)) || any(count_matrix < 0) ||
      any(abs(count_matrix - round(count_matrix)) > 1e-8)) {
    stop("DESeq2 input must be finite, nonnegative integer raw counts")
  }
  storage.mode(count_matrix) <- "integer"
  strip_versions <- isTRUE(context$config$inputs$strip_gene_version %||% TRUE)
  gene_id <- as.character(raw[[gene_id_column]])
  if (strip_versions) gene_id <- sub("\\.[0-9]+$", "", gene_id)
  if (anyNA(gene_id) || any(!nzchar(gene_id))) stop("Gene IDs must be non-empty")
  if (anyDuplicated(gene_id)) stop("Duplicate gene IDs require an explicit upstream aggregation decision")
  rownames(count_matrix) <- gene_id
  input_symbol <- if (is.null(gene_symbol_column)) rep(NA_character_, length(gene_id)) else as.character(raw[[gene_symbol_column]])
  input_entrez <- if (is.null(entrez_column)) rep(NA_character_, length(gene_id)) else as.character(raw[[entrez_column]])
  annotation_resource <- context$config$resources$annotation_map %||% NULL
  if (!is.null(annotation_resource)) {
    frozen <- readRDS(resource_path(context, annotation_resource))$mapping
    idx <- match(gene_id, frozen$gene_id)
    # The legacy table contains the tested genes after count filtering. Genes
    # outside that set cannot enter DE/enrichment, so retain input annotation
    # for them while requiring complete coverage of the eventual tested set.
    gene_symbol <- input_symbol
    entrez_id <- input_entrez
    covered <- !is.na(idx)
    gene_symbol[covered] <- frozen$gene_symbol[idx[covered]]
    entrez_id[covered] <- frozen$entrez_id[idx[covered]]
  } else if (identical(context$config$species, "human") && all(grepl("^ENSG", head(gene_id, 10L)))) {
    # Legacy step 01 used org.Hs.eg.db mapIds(multiVals = "first") and only
    # fell back to the input symbol when SYMBOL was unmapped. Reusing input
    # annotation directly changes renamed genes, duplicate-symbol handling,
    # ORA universe membership and GSEA set sizes.
    gene_symbol <- unname(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db, keys = gene_id, column = "SYMBOL",
      keytype = "ENSEMBL", multiVals = "first"
    ))
    entrez_id <- unname(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db, keys = gene_id, column = "ENTREZID",
      keytype = "ENSEMBL", multiVals = "first"
    ))
    missing_symbol <- is.na(gene_symbol) | !nzchar(gene_symbol)
    gene_symbol[missing_symbol] <- input_symbol[missing_symbol]
    missing_symbol <- is.na(gene_symbol) | !nzchar(gene_symbol)
    gene_symbol[missing_symbol] <- gene_id[missing_symbol]
  } else {
    gene_symbol <- input_symbol
    entrez_id <- input_entrez
  }
  annotation <- data.frame(gene_id = gene_id, gene_symbol = gene_symbol,
                           entrez_id = entrez_id, stringsAsFactors = FALSE)
  rownames(annotation) <- gene_id

  samples <- context$samples
  rownames(samples) <- samples$sample
  display_factor <- context$config$design$display_factor %||% as.character(context$contrasts$factor[[1L]])
  if (!display_factor %in% names(samples)) stop("Display factor is absent from samples.tsv: ", display_factor)
  samples$display_group <- factor(samples[[display_factor]], levels = unique(as.character(samples[[display_factor]])))

  min_count <- as.integer(context$config$filtering$min_count %||% 10L)
  min_samples <- as.integer(context$config$filtering$min_samples %||% 3L)
  selected <- samples$sample[!samples$excluded]
  if (length(selected) < 2L) stop("Fewer than two approved samples remain")
  if (min_samples > length(selected)) stop("filtering.min_samples exceeds the approved sample count")
  keep_initial <- rowSums(count_matrix >= min_count) >= min_samples
  keep_final <- rowSums(count_matrix[, selected, drop = FALSE] >= min_count) >= min_samples
  if (sum(keep_final) < 2L) stop("Gene filtering retained fewer than two genes")

  root <- file.path(context$run_dir, "01_QC")
  figures <- file.path(root, "Figures"); tables <- file.path(root, "Tables")
  dir.create(figures, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)

  qc_state <- function(sample_ids, keep, label, stem) {
    metadata <- samples[sample_ids, , drop = FALSE]
    dds <- DESeq2::DESeqDataSetFromMatrix(count_matrix[keep, sample_ids, drop = FALSE], metadata, design = ~1)
    vst <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = TRUE))
    pca <- stats::prcomp(t(vst), scale. = FALSE)
    explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
    coordinates <- data.frame(
      sample = rownames(pca$x), canonical_sample = metadata[rownames(pca$x), "canonical_sample"],
      group = metadata[rownames(pca$x), "display_group"], PC1 = pca$x[, 1], PC2 = pca$x[, 2],
      stringsAsFactors = FALSE
    )
    write_result_table(coordinates, file.path(tables, paste0(stem, "_Coordinates.tsv")))
    plot <- plot_pca_publication(coordinates, explained[[1]], explained[[2]], label)
    save_publication_figure(plot, file.path(figures, stem), 145, 118)
    list(vst = vst, metadata = metadata, coordinates = coordinates, explained = explained)
  }
  initial <- qc_state(samples$sample, keep_initial, "Initial PCA — before approved exclusions", "Initial_PCA")
  final <- qc_state(selected, keep_final, "Final PCA — approved samples", "Final_PCA")

  metrics <- data.frame(
    sample = samples$sample, canonical_sample = samples$canonical_sample,
    group = as.character(samples$display_group), library_size = colSums(count_matrix),
    detected_genes = colSums(count_matrix > 0), zero_fraction = colMeans(count_matrix == 0),
    excluded = samples$excluded, exclusion_reason = samples$exclusion_reason,
    stringsAsFactors = FALSE
  )
  write_result_table(metrics, file.path(tables, "Sample_QC_Metrics.tsv"))
  palette <- publication_palette(levels(samples$display_group))
  library_plot <- ggplot2::ggplot(metrics, ggplot2::aes(stats::reorder(canonical_sample, library_size), library_size,
                                                        fill = group, alpha = !excluded)) +
    ggplot2::geom_col(width = 0.72) + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.38), guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::labs(title = "Library size and approved exclusions", x = NULL, y = "Total reads", fill = display_factor,
                  caption = "Faded bars are exclusions already approved in samples.tsv; this plot does not choose exclusions.") +
    theme_publication()
  save_publication_figure(library_plot, file.path(figures, "Library_Size"), 145, 105)

  correlation <- stats::cor(final$vst, method = "pearson")
  distance <- as.matrix(stats::dist(t(final$vst)))
  canonical <- final$metadata$canonical_sample
  rownames(correlation) <- colnames(correlation) <- canonical
  rownames(distance) <- colnames(distance) <- canonical
  write_result_table(data.frame(sample = rownames(correlation), correlation, check.names = FALSE),
                     file.path(tables, "Sample_Correlation_Pearson.tsv"))
  write_result_table(data.frame(sample = rownames(distance), distance, check.names = FALSE),
                     file.path(tables, "Sample_Distance_Euclidean.tsv"))
  annotation_col <- data.frame(group = final$metadata$display_group)
  rownames(annotation_col) <- canonical
  heat_colors <- grDevices::colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(101)
  cor_heat <- pheatmap::pheatmap(correlation, annotation_col = annotation_col, annotation_row = annotation_col,
                                color = heat_colors, border_color = NA, silent = TRUE,
                                main = "Sample correlation")$gtable
  dist_heat <- pheatmap::pheatmap(distance, annotation_col = annotation_col, annotation_row = annotation_col,
                                 color = rev(heat_colors), border_color = NA, silent = TRUE,
                                 main = "Sample distance")$gtable
  save_publication_grob(cor_heat, file.path(figures, "Sample_Correlation_Heatmap"), 145, 125)
  save_publication_grob(dist_heat, file.path(figures, "Sample_Distance_Heatmap"), 145, 125)

  clustering <- stats::hclust(stats::as.dist(distance), method = context$config$qc$clustering_method %||% "complete")
  write_result_table(data.frame(order = seq_along(clustering$order), canonical_sample = clustering$labels[clustering$order]),
                     file.path(tables, "Sample_Clustering_Order.tsv"))
  if ("png" %in% requested_figure_formats()) {
    grDevices::png(file.path(figures, "Sample_Hierarchical_Clustering.png"), width = 145 / 25.4,
                   height = 105 / 25.4, units = "in", res = 400, type = "cairo")
    graphics::par(mar = c(4, 4, 2, 1)); graphics::plot(clustering, main = "Sample hierarchical clustering", xlab = "", sub = "", hang = -1)
    grDevices::dev.off()
  }
  grDevices::cairo_pdf(file.path(figures, "Sample_Hierarchical_Clustering.pdf"), width = 145 / 25.4, height = 105 / 25.4)
  graphics::par(mar = c(4, 4, 2, 1)); graphics::plot(clustering, main = "Sample hierarchical clustering", xlab = "", sub = "", hang = -1)
  grDevices::dev.off()

  mean_correlation <- vapply(seq_len(ncol(correlation)), function(i) mean(correlation[-i, i]), numeric(1))
  correlation_qc <- data.frame(canonical_sample = canonical, mean_pairwise_correlation = mean_correlation,
                               group = as.character(final$metadata$display_group))
  write_result_table(correlation_qc, file.path(tables, "Sample_Correlation_QC.tsv"))
  correlation_plot <- ggplot2::ggplot(correlation_qc, ggplot2::aes(stats::reorder(canonical_sample, mean_pairwise_correlation),
                                                                   mean_pairwise_correlation, fill = group)) +
    ggplot2::geom_col(width = 0.72) + ggplot2::coord_flip() + ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(title = "Mean correlation with other approved samples", x = NULL,
                  y = "Mean Pearson correlation", fill = display_factor) + theme_publication()
  save_publication_figure(correlation_plot, file.path(figures, "Sample_Correlation_QC_Barplot"), 145, 105)
  write_result_table(data.frame(total_genes = nrow(count_matrix), retained_initial = sum(keep_initial),
                                retained_final = sum(keep_final), removed_final = sum(!keep_final),
                                min_count = min_count, min_samples = min_samples),
                     file.path(tables, "Gene_Filtering_Summary.tsv"))

  context$state$counts <- count_matrix
  context$state$annotation <- annotation
  context$state$samples <- samples
  context$state$selected_samples <- selected
  context$state$keep_final <- keep_final
  context$state$qc_vst <- final$vst
  context$state$display_factor <- display_factor
  saveRDS(list(vst = final$vst, metadata = final$metadata, filtering = list(min_count = min_count, min_samples = min_samples,
                                                                          keep_final = keep_final)),
          file.path(context$run_dir, "Provenance", "qc_state.rds"), compress = "xz")
  record_module_status(context, "01_QC", "complete",
                       sprintf("%d approved samples; %d/%d genes retained", length(selected), sum(keep_final), nrow(count_matrix)))
}
