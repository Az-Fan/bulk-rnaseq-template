# Module 04 — legacy-compatible regulation and activity analysis
#
# Migrates legacy steps 05 (GTRD), 06, 10, 11 and 13 without changing their
# statistical inputs. Snakemake schedules this module; it does not replace the
# scientific methods. All networks and gene sets are frozen resources.

run_activity <- function(context) {
  state <- canonical_module_state(context$config, c("regulation"))
  if (!identical(state$status, "enabled")) {
    return(record_module_status(context, "04_Regulation", state$status, state$reason))
  }
  if (is.null(context$state$vst) || is.null(context$state$differential)) {
    stop("Differential expression must complete before regulation analysis")
  }
  root <- file.path(context$run_dir, "04_Regulation")
  tf_root <- file.path(root, "TF_Activity")
  gsva_root <- file.path(root, "GSVA")
  progeny_root <- file.path(root, "PROGENy")
  tf_gsea_root <- file.path(root, "TF_Target_GSEA")
  for (path in c(file.path(tf_root, "Tables"), file.path(tf_root, "Figures"),
                 file.path(gsva_root, "Tables"), file.path(gsva_root, "Figures"),
                 file.path(progeny_root, "Tables"), file.path(progeny_root, "Figures"),
                 file.path(tf_gsea_root, "Tables"), file.path(tf_gsea_root, "Figures"))) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }

  ids <- unlist(context$config$resources$activity %||% list(), use.names = FALSE)
  collectri_id <- ids[grepl("collectri", ids, ignore.case = TRUE)][1]
  progeny_id <- ids[grepl("progeny", ids, ignore.case = TRUE)][1]
  if (is.na(collectri_id) || is.na(progeny_id)) {
    stop("Regulation requires frozen CollecTRI and PROGENy resources")
  }
  collectri <- readRDS(resource_path(context, collectri_id))
  progeny <- readRDS(resource_path(context, progeny_id))
  if (!all(c("source", "target", "mor") %in% names(collectri))) stop("CollecTRI resource lacks source/target/mor")
  if (!all(c("source", "target", "weight") %in% names(progeny))) stop("PROGENy resource lacks source/target/weight")
  collectri <- unique(collectri[, c("source", "target", "mor")])
  progeny <- unique(progeny[, c("source", "target", "weight")])

  annotation <- context$state$annotation
  expression <- context$state$vst
  map <- annotation[rownames(expression), c("gene_id", "gene_symbol"), drop = FALSE]
  valid <- !is.na(map$gene_symbol) & nzchar(map$gene_symbol)
  expr_frame <- data.frame(gene_symbol = map$gene_symbol[valid], expression[valid, , drop = FALSE], check.names = FALSE)
  expression_symbol <- stats::aggregate(. ~ gene_symbol, data = expr_frame, FUN = mean)
  rownames(expression_symbol) <- expression_symbol$gene_symbol
  expression_symbol$gene_symbol <- NULL
  expression_symbol <- as.matrix(expression_symbol)
  metadata <- context$state$samples[context$state$selected_samples, , drop = FALSE]
  rownames(metadata) <- metadata$sample
  regulation_config <- context$config$regulation %||% list()
  top_tfs_n <- as.integer(regulation_config$top_tfs %||% 30L)
  top_pathways_n <- as.integer(regulation_config$top_pathways %||% 30L)
  target_tfs <- as.character(regulation_config$target_tfs %||% character())
  focus_pathway <- as.character(regulation_config$focus_pathway %||% "MAPK")
  minimum_targets <- as.integer(regulation_config$minimum_targets %||% 5L)
  tf_gsea_fdr <- as.numeric(regulation_config$tf_gsea_fdr %||% 0.05)
  tf_gsea_curves <- as.integer(regulation_config$tf_gsea_curves %||% 5L)
  tf_gsea_min_size <- as.integer(regulation_config$tf_gsea_min_size %||% 10L)
  tf_gsea_max_size <- as.integer(regulation_config$tf_gsea_max_size %||% 500L)
  regulation_seed <- as.integer(regulation_config$random_seed %||% 123L)
  correlation_method <- as.character(regulation_config$correlation_method %||% "pearson")
  correlation_pathways_n <- as.integer(regulation_config$correlation_pathways %||% 50L)
  contrast_overview_tfs <- as.integer(regulation_config$contrast_overview_tfs %||% 25L)
  target_tf_genes_n <- as.integer(regulation_config$target_tf_genes %||% 30L)
  tf_expression_count <- as.integer(regulation_config$tf_expression_count %||% 20L)
  tf_activity_pvalue <- as.numeric(regulation_config$tf_activity_pvalue %||% 0.05)
  external_tf_gsea_rank <- as.character(regulation_config$external_tf_gsea_rank %||% "wald_statistic")
  missing_target_tfs <- setdiff(target_tfs, unique(collectri$source))
  if (length(missing_target_tfs)) {
    stop("Declared target_tfs are absent from the frozen CollecTRI resource: ", paste(missing_target_tfs, collapse = ", "))
  }
  if (!focus_pathway %in% unique(progeny$source)) {
    stop("Declared focus_pathway is absent from the frozen PROGENy resource: ", focus_pathway)
  }

  # Sample-level inference is the same run_ulm/run_mlm calculation as legacy.
  tf_sample <- decoupleR::run_ulm(
    mat = expression_symbol, net = collectri, .source = "source", .target = "target",
    .mor = "mor", minsize = minimum_targets
  )
  progeny_sample <- decoupleR::run_mlm(
    mat = expression_symbol, net = progeny, .source = "source", .target = "target",
    .mor = "weight", minsize = minimum_targets
  )
  tf_sample_df <- as.data.frame(tf_sample)
  progeny_sample_df <- as.data.frame(progeny_sample)
  write_result_table(tf_sample_df, file.path(tf_root, "Tables", "TF_Activity_Sample_Full.tsv"))
  write_result_table(progeny_sample_df, file.path(progeny_root, "Tables", "PROGENy_Sample_Activity.tsv"))

  tf_matrix_df <- tidyr::pivot_wider(tf_sample_df[, c("source", "condition", "score")],
                                     names_from = "condition", values_from = "score")
  tf_matrix <- as.matrix(tf_matrix_df[, -1, drop = FALSE]); rownames(tf_matrix) <- tf_matrix_df$source
  write_result_table(data.frame(TF = rownames(tf_matrix), tf_matrix, check.names = FALSE),
                     file.path(tf_root, "Tables", "TF_Activity_Sample_Matrix.tsv"))

  progeny_matrix_df <- tidyr::pivot_wider(progeny_sample_df[, c("condition", "source", "score")],
                                          names_from = "source", values_from = "score")
  progeny_matrix <- as.matrix(progeny_matrix_df[, -1, drop = FALSE]); rownames(progeny_matrix) <- progeny_matrix_df$condition

  # Hallmark GSVA exactly follows legacy gsvaParam(kcdf = "Gaussian").
  payload <- readRDS(resource_path(context, context$config$resources$enrichment_gene_sets))
  if (!identical(payload$msigdb_release, "7.5.1")) stop("Regulation regression requires frozen MSigDB 7.5.1")
  hallmark <- payload$gene_sets[payload$gene_sets$gs_cat == "H", c("gs_name", "gene_symbol")]
  hallmark_sets <- lapply(split(hallmark$gene_symbol, hallmark$gs_name), unique)
  gsva_par <- GSVA::gsvaParam(exprData = expression_symbol, geneSets = hallmark_sets, kcdf = "Gaussian")
  gsva_scores <- GSVA::gsva(gsva_par, BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
  write_result_table(data.frame(pathway = rownames(gsva_scores), gsva_scores, check.names = FALSE),
                     file.path(gsva_root, "Tables", "GSVA_Score_Matrix.tsv"))

  # Legacy correlation table uses the 30 most variable TFs (+ declared targets)
  # and at most the 50 most variable pathways. No samples are removed.
  shared_samples <- intersect(colnames(tf_matrix), colnames(gsva_scores))
  tf_sd <- apply(tf_matrix[, shared_samples, drop = FALSE], 1, stats::sd)
  selected_tfs <- head(names(sort(tf_sd, decreasing = TRUE)), min(top_tfs_n, length(tf_sd)))
  selected_tfs <- unique(c(selected_tfs, intersect(target_tfs, rownames(tf_matrix))))
  pathway_sd <- apply(gsva_scores[, shared_samples, drop = FALSE], 1, stats::sd)
  correlation_pathways <- head(names(sort(pathway_sd, decreasing = TRUE)),
                               min(correlation_pathways_n, nrow(gsva_scores)))
  correlation <- stats::cor(t(tf_matrix[selected_tfs, shared_samples, drop = FALSE]),
                            t(gsva_scores[correlation_pathways, shared_samples, drop = FALSE]),
                            method = correlation_method)
  write_result_table(data.frame(TF = rownames(correlation), correlation, check.names = FALSE),
                     file.path(gsva_root, "Tables", "TF_Pathway_Correlation_Matrix.tsv"))

  group_annotation <- data.frame(Group = metadata[colnames(gsva_scores), context$config$design$display_factor])
  rownames(group_annotation) <- colnames(gsva_scores)
  top_pathways <- head(names(sort(pathway_sd, decreasing = TRUE)), min(top_pathways_n, length(pathway_sd)))
  gsva_heatmap <- pheatmap::pheatmap(
    gsva_scores[top_pathways, , drop = FALSE], annotation_col = group_annotation,
    cluster_rows = TRUE, cluster_cols = TRUE, scale = "row", silent = TRUE,
    color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    border_color = NA, fontsize_row = 7, main = paste0("Top ", length(top_pathways), " variable Hallmark pathways")
  )
  save_publication_grob(gsva_heatmap$gtable, file.path(gsva_root, "Figures", "GSVA_Pathway_Heatmap"), 183, 155)
  cor_heatmap <- pheatmap::pheatmap(
    correlation, cluster_rows = TRUE, cluster_cols = TRUE, silent = TRUE,
    color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-1, 1, length.out = 101), border_color = "#F2F2F2",
    fontsize_row = 7, fontsize_col = 7, main = "TF activity–Hallmark pathway correlation"
  )
  save_publication_grob(cor_heatmap$gtable, file.path(gsva_root, "Figures", "TF_Pathway_Correlation_Heatmap"), 210, 175)

  # PROGENy sample heatmap uses the legacy across-sample z-score operation.
  progeny_z <- scale(progeny_matrix)
  progeny_heatmap <- pheatmap::pheatmap(
    progeny_z, cluster_rows = TRUE, cluster_cols = TRUE, silent = TRUE,
    color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    border_color = "white", main = "PROGENy pathway activity (z-score)"
  )
  save_publication_grob(progeny_heatmap$gtable, file.path(progeny_root, "Figures", "Pathway_Activity_Heatmap"), 178, 135)

  normalize_gsea <- function(x, database) {
    if (!nrow(x)) return(data.frame())
    x$database <- database; x$pathway <- x$ID; x$padj <- x$p.adjust
    x$size <- x$setSize; x$ES <- x$enrichmentScore; x$leadingEdge <- x$core_enrichment
    x
  }
  run_tf_gsea <- function(ranks, term2gene, database, contrast_id, output_dir) {
    set.seed(regulation_seed)
    object <- clusterProfiler::GSEA(
      ranks, TERM2GENE = term2gene, pvalueCutoff = 1, pAdjustMethod = "BH",
      minGSSize = tf_gsea_min_size, maxGSSize = tf_gsea_max_size, seed = TRUE, by = "fgsea",
      BPPARAM = BiocParallel::SerialParam(), verbose = FALSE
    )
    table <- normalize_gsea(as.data.frame(object), database)
    dir.create(file.path(output_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_dir, "Figures", "GSEA_Curve_Plots"), recursive = TRUE, showWarnings = FALSE)
    write_result_table(table, file.path(output_dir, "Tables", paste0("GSEA_Full_Table_", database, ".tsv")))
    if (nrow(table)) {
      overview <- plot_gsea_overview(table, paste0(contrast_id, " — TF-target GSEA: ", database), top_n = 12L)
      save_publication_figure(overview, file.path(output_dir, "Figures", paste0("Summary_Directional_", database)), 183, 135)
      reportable <- table[!is.na(table$p.adjust) & table$p.adjust < tf_gsea_fdr, , drop = FALSE]
      curve_rows <- head(reportable[order(reportable$p.adjust, -abs(reportable$NES)), , drop = FALSE],
                         tf_gsea_curves)
      pathways <- split(term2gene$gene, term2gene$term)
      for (i in seq_len(nrow(curve_rows))) {
        row <- curve_rows[i, , drop = FALSE]
        curve <- plot_gsea_curve(pathways[[row$ID]], ranks, row, contrast_id)
        stem <- sprintf("%02d_%s", i, safe_id(gsub("[^A-Za-z0-9_.-]+", "_", row$ID), "TF gene set"))
        save_publication_figure(curve$plot, file.path(output_dir, "Figures", "GSEA_Curve_Plots", stem), 178, 140)
        write_result_table(curve$running, file.path(output_dir, "Tables", "Curve_Data", paste0(stem, "_Running_ES.tsv.gz")))
        write_result_table(curve$hits, file.path(output_dir, "Tables", "Curve_Data", paste0(stem, "_Hits.tsv.gz")))
      }
    }
    list(object = object, table = table)
  }

  # Frozen external TF databases retain their original GMT content and hashes.
  external_tf_sets <- list()
  for (resource_id in context$config$resources$tf_databases %||% character()) {
    label <- if (grepl("Consensus", resource_id)) "ChEA_Consensus" else if (grepl("ENCODE_TF", resource_id)) "ENCODE_2015" else safe_id(resource_id)
    table <- clusterProfiler::read.gmt(resource_path(context, resource_id))
    names(table) <- c("term", "gene")
    external_tf_sets[[label]] <- unique(table)
  }
  frozen_sets <- payload$gene_sets
  gtrd <- unique(frozen_sets[frozen_sets$gs_cat == "C3" & frozen_sets$gs_subcat == "TFT:GTRD",
                             c("gs_name", "gene_symbol")])
  names(gtrd) <- c("term", "gene")

  contrasts <- list()
  for (contrast_id in names(context$state$differential)) {
    de <- context$state$differential[[contrast_id]]$result
    ranked_de <- de[is.finite(de$stat) & !is.na(de$gene_symbol) & nzchar(de$gene_symbol), c("gene_symbol", "stat")]
    ranked_de <- ranked_de[order(-abs(ranked_de$stat)), , drop = FALSE]
    ranked_de <- ranked_de[!duplicated(ranked_de$gene_symbol), , drop = FALSE]
    stat_matrix <- matrix(ranked_de$stat, ncol = 1, dimnames = list(ranked_de$gene_symbol, contrast_id))
    tf_contrast <- as.data.frame(decoupleR::run_ulm(
      mat = stat_matrix, net = collectri, .source = "source", .target = "target", .mor = "mor", minsize = minimum_targets
    ))
    progeny_contrast <- as.data.frame(decoupleR::run_mlm(
      mat = stat_matrix, net = progeny, .source = "source", .target = "target", .mor = "weight", minsize = minimum_targets
    ))
    tf_contrast <- tf_contrast[order(-tf_contrast$score), , drop = FALSE]
    write_result_table(tf_contrast, file.path(tf_root, "Tables", paste0("TF_Activity_Contrast_", contrast_id, ".tsv")))
    write_result_table(progeny_contrast, file.path(progeny_root, "Tables", paste0("PROGENy_Contrast_Activity_", contrast_id, ".tsv")))

    top <- head(tf_contrast[order(-abs(tf_contrast$score)), , drop = FALSE], contrast_overview_tfs)
    top$direction <- ifelse(top$score >= 0, "Up", "Down")
    top$source <- factor(top$source, levels = rev(top$source))
    tf_bar <- ggplot2::ggplot(top, ggplot2::aes(score, source, fill = direction)) +
      ggplot2::geom_col(width = 0.72) + ggplot2::geom_vline(xintercept = 0, colour = "#343A40", linewidth = 0.4) +
      ggplot2::scale_fill_manual(values = direction_palette[c("Down", "Up")]) +
      ggplot2::labs(title = paste0(contrast_id, " — CollecTRI ULM activity"),
                    subtitle = "Frozen network; complete DESeq2 Wald-statistic input",
                    x = "ULM activity score", y = NULL, fill = NULL,
                    caption = "Expression-derived TF activity is not direct binding or causal evidence.") + theme_publication()
    save_publication_figure(tf_bar, file.path(tf_root, "Figures", paste0("TF_Activity_", contrast_id, "_Overview")), 178, 145)

    progeny_contrast$direction <- ifelse(progeny_contrast$score >= 0, "Up", "Down")
    progeny_contrast$source <- factor(progeny_contrast$source,
                                      levels = progeny_contrast$source[order(progeny_contrast$score)])
    progeny_bar <- ggplot2::ggplot(progeny_contrast, ggplot2::aes(score, source, fill = direction)) +
      ggplot2::geom_col(width = 0.72) + ggplot2::geom_vline(xintercept = 0, colour = "#343A40", linewidth = 0.4) +
      ggplot2::scale_fill_manual(values = direction_palette[c("Down", "Up")]) +
      ggplot2::labs(title = paste0(contrast_id, " — PROGENy pathway activity"),
                    x = "MLM activity score", y = NULL, fill = NULL,
                    caption = "Expression-derived pathway activity is not causal evidence.") + theme_publication()
    save_publication_figure(progeny_bar, file.path(progeny_root, "Figures", paste0("Pathway_Activity_", contrast_id, "_Barplot")), 150, 115)

    focus <- focus_pathway
    focus_data <- progeny[progeny$source == focus & progeny$target %in% rownames(stat_matrix), , drop = FALSE]
    focus_data$stat <- stat_matrix[focus_data$target, 1]
    focus_data$concordance <- ifelse(sign(focus_data$weight) == sign(focus_data$stat), "Concordant", "Discordant")
    focus_plot <- ggplot2::ggplot(focus_data, ggplot2::aes(weight, stat, colour = concordance)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 3, colour = "#8C9298") +
      ggplot2::geom_vline(xintercept = 0, linetype = 3, colour = "#8C9298") +
      ggplot2::geom_point(size = 2, alpha = 0.75) +
      ggrepel::geom_text_repel(ggplot2::aes(label = target), max.overlaps = 20, seed = 104729, size = 2.5) +
      ggplot2::scale_colour_manual(values = c(Concordant = "#009E73", Discordant = "#CC79A7")) +
      ggplot2::labs(title = paste0(focus, " pathway target genes"), x = "PROGENy target weight",
                    y = "DESeq2 Wald statistic", colour = NULL) + theme_publication()
    save_publication_figure(focus_plot, file.path(progeny_root, "Figures", paste0(focus, "_Target_Genes_", contrast_id)), 155, 125)
    write_result_table(focus_data, file.path(progeny_root, "Tables", paste0(focus, "_Target_Genes_", contrast_id, ".tsv")))

    de_input_order <- de[match(rownames(context$state$dds), de$gene_id), , drop = FALSE]
    lfc_rank <- de_input_order[is.finite(de_input_order$log2FoldChange_ashr) & !is.na(de_input_order$gene_symbol) &
                                           nzchar(de_input_order$gene_symbol), c("gene_symbol", "log2FoldChange_ashr")]
    lfc_rank <- lfc_rank[!duplicated(lfc_rank$gene_symbol), , drop = FALSE]
    lfc_ranks <- lfc_rank$log2FoldChange_ashr; names(lfc_ranks) <- lfc_rank$gene_symbol
    lfc_ranks <- sort(lfc_ranks, decreasing = TRUE)
    wald_ranks <- context$state$enrichment[[contrast_id]]$ranks
    tf_gsea <- list(
      GTRD_TF = run_tf_gsea(wald_ranks, gtrd, "GTRD_TF", contrast_id, file.path(tf_gsea_root, "GTRD"))
    )
    for (label in names(external_tf_sets)) {
      folder <- if (grepl("ChEA", label)) "ChEA" else if (grepl("ENCODE", label)) "ENCODE" else label
      external_ranks <- if (identical(external_tf_gsea_rank, "legacy_mixed")) lfc_ranks else wald_ranks
      tf_gsea[[label]] <- run_tf_gsea(external_ranks, external_tf_sets[[label]], label, contrast_id,
                                      file.path(tf_gsea_root, folder))
    }

    # Legacy step 10 also reports agreement between CollecTRI ULM and its own
    # TF-target GSEA (both driven by the complete Wald-statistic vector).
    set.seed(regulation_seed)
    collectri_gsea <- clusterProfiler::GSEA(
      wald_ranks, TERM2GENE = unique(collectri[, c("source", "target")]),
      pvalueCutoff = 1, minGSSize = minimum_targets, maxGSSize = tf_gsea_max_size, eps = 0,
      nPermSimple = 10000, BPPARAM = BiocParallel::SerialParam(), verbose = FALSE
    )
    collectri_table <- as.data.frame(collectri_gsea)
    concordance <- merge(tf_contrast[, c("source", "score", "p_value")],
                         collectri_table[, c("ID", "NES", "p.adjust")],
                         by.x = "source", by.y = "ID", all = FALSE)
    concordance$classification <- ifelse(
      concordance$p_value < tf_activity_pvalue & concordance$p.adjust < tf_gsea_fdr &
        sign(concordance$score) == sign(concordance$NES), "Consistent", "Not_significant_or_discordant"
    )
    write_result_table(concordance, file.path(tf_root, "Tables", paste0("TF_ULM_GSEA_Concordance_", contrast_id, ".tsv")))
    concordance_plot <- ggplot2::ggplot(concordance, ggplot2::aes(score, NES, colour = classification)) +
      ggplot2::geom_hline(yintercept = 0, colour = "#D7DADE") + ggplot2::geom_vline(xintercept = 0, colour = "#D7DADE") +
      ggplot2::geom_point(alpha = 0.75, size = 1.7) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linetype = 2, colour = "#343A40", linewidth = 0.5) +
      ggplot2::scale_colour_manual(values = c(Consistent = "#D55E00", Not_significant_or_discordant = "#BDBDBD")) +
      ggplot2::labs(title = "TF inference method concordance", x = "CollecTRI ULM activity", y = "TF-target GSEA NES",
                    colour = NULL, caption = "Agreement is supportive consistency, not direct binding or causal evidence.") + theme_publication()
    save_publication_figure(concordance_plot, file.path(tf_root, "Figures", paste0("TF_ULM_GSEA_Concordance_", contrast_id)), 150, 120)

    # Declared TF-expression display family; legacy projects retain count 20.
    # selection priority: significant TF-GSEA hits followed by |ULM score|.
    gsea_ranked <- concordance$source[order(concordance$p.adjust, -abs(concordance$score))]
    gsea_ranked <- gsea_ranked[concordance$p.adjust[match(gsea_ranked, concordance$source)] < tf_gsea_fdr]
    activity_ranked <- tf_contrast$source[order(-abs(tf_contrast$score))]
    selected_tfs <- head(unique(c(gsea_ranked, activity_ranked)), tf_expression_count)
    selected_tfs <- selected_tfs[selected_tfs %in% rownames(expression_symbol)]
    tf_expression <- data.frame(
      TF = rep(selected_tfs, each = ncol(expression_symbol)),
      sample = rep(colnames(expression_symbol), times = length(selected_tfs)),
      expression = as.numeric(t(expression_symbol[selected_tfs, , drop = FALSE])),
      stringsAsFactors = FALSE
    )
    tf_expression$group <- metadata[tf_expression$sample, context$state$display_factor]
    tf_expression$TF <- factor(tf_expression$TF, levels = selected_tfs)
    write_result_table(tf_expression, file.path(tf_root, "Tables", paste0("Top20_TF_Expression_", contrast_id, ".tsv")))
    group_palette <- publication_palette(unique(tf_expression$group))
    expression_plot <- ggplot2::ggplot(tf_expression, ggplot2::aes(TF, expression, colour = group, fill = group)) +
      ggplot2::geom_boxplot(width = 0.42, alpha = 0.18, outlier.shape = NA,
                            position = ggplot2::position_dodge(width = 0.62), linewidth = 0.4) +
      ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.06, dodge.width = 0.62),
                          size = 1.6, alpha = 0.9) +
      ggplot2::scale_colour_manual(values = group_palette) + ggplot2::scale_fill_manual(values = group_palette) +
      ggplot2::labs(title = "Top TF expression",
                    subtitle = sprintf("Declared %d-TF selection rule; unchanged VST values", tf_expression_count),
                    x = NULL, y = "VST expression", colour = "Group", fill = "Group") + theme_publication() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "bottom")
    save_publication_figure(expression_plot, file.path(tf_root, "Figures", paste0("Top20_TF_Expression_", contrast_id)), 230, 125)

    # One multi-page report and one report per explicitly declared target TF.
    global_report <- file.path(tf_root, "Figures", paste0("Global_TF_Analysis_Report_", contrast_id, ".pdf"))
    grDevices::cairo_pdf(global_report, width = 8, height = 7, onefile = TRUE)
    print(tf_bar); print(concordance_plot); print(expression_plot)
    variable_tfs <- head(names(sort(apply(tf_matrix, 1, stats::sd), decreasing = TRUE)), min(top_tfs_n, nrow(tf_matrix)))
    print(pheatmap::pheatmap(tf_matrix[variable_tfs, , drop = FALSE], scale = "row", silent = TRUE,
                             color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                             main = "Top variable TF activities")$gtable)
    grDevices::dev.off()
    for (target_tf in target_tfs[target_tfs %in% tf_contrast$source]) {
      report_path <- file.path(tf_root, "Figures", "single_tf_plots", paste0(safe_id(target_tf), "_Analysis_Report.pdf"))
      dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
      subset <- head(tf_contrast, tf_expression_count)
      if (!target_tf %in% subset$source) subset <- rbind(subset, tf_contrast[tf_contrast$source == target_tf, ])
      subset$highlight <- ifelse(subset$source == target_tf, "Target", "Other")
      subset$source <- factor(subset$source, levels = subset$source[order(subset$score)])
      rank_plot <- ggplot2::ggplot(subset, ggplot2::aes(score, source, fill = highlight)) +
        ggplot2::geom_col() + ggplot2::scale_fill_manual(values = c(Target = "#E69F00", Other = "#BDBDBD")) +
        ggplot2::labs(title = paste0("Activity rank: ", target_tf), x = "ULM activity", y = NULL, fill = NULL) + theme_publication()
      target_genes <- head(collectri$target[collectri$source == target_tf][order(-abs(collectri$mor[collectri$source == target_tf]))],
                           target_tf_genes_n)
      target_genes <- intersect(target_genes, rownames(expression_symbol))
      grDevices::cairo_pdf(report_path, width = 8, height = 6, onefile = TRUE)
      print(rank_plot)
      if (target_tf %in% collectri_table$ID) {
        row <- normalize_gsea(collectri_table[collectri_table$ID == target_tf, , drop = FALSE], "CollecTRI")
        print(plot_gsea_curve(split(collectri$target, collectri$source)[[target_tf]], wald_ranks, row, contrast_id)$plot)
      }
      if (length(target_genes) >= 3L) {
        print(pheatmap::pheatmap(expression_symbol[target_genes, , drop = FALSE], scale = "row", silent = TRUE,
                                 annotation_col = group_annotation,
                                 color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                                 main = paste0("Top targets: ", target_tf))$gtable)
      }
      grDevices::dev.off()
    }

    workbook <- list(DESeq2_Wald_statistic = de, PROGENy_Sample_Activity = progeny_sample_df,
                     PROGENy_Contrast_Activity = progeny_contrast, Focus_Pathway_Targets = focus_data)
    openxlsx::write.xlsx(workbook, file.path(progeny_root, "Tables", paste0("All_Analysis_Results_", contrast_id, ".xlsx")),
                         rowNames = FALSE, overwrite = TRUE)
    contrasts[[contrast_id]] <- list(tf = tf_contrast, progeny = progeny_contrast,
                                     tf_gsea = tf_gsea, collectri_gsea = collectri_gsea,
                                     concordance = concordance)
  }

  write_result_table(data.frame(
    analysis = c("TF activity", "PROGENy", "GSVA", "TF-target GSEA"),
    method = c("decoupleR run_ulm", "decoupleR run_mlm", "GSVA gsvaParam Gaussian", "clusterProfiler GSEA"),
    input = c("VST samples / DESeq2 Wald statistic", "VST samples / DESeq2 Wald statistic",
              "VST expression", "Wald statistic for GTRD; ashr LFC for ChEA/ENCODE"),
    resource = c(collectri_id, progeny_id, context$config$resources$enrichment_gene_sets,
                 paste(c("MSigDB 7.5.1 GTRD", context$config$resources$tf_databases %||% character()), collapse = ";")),
    stringsAsFactors = FALSE
  ), file.path(root, "Regulation_Provenance.tsv"))
  context$state$activity <- list(tf_sample = tf_sample, progeny_sample = progeny_sample,
                                 gsva = gsva_scores, correlation = correlation, contrasts = contrasts)
  record_module_status(context, "04_Regulation", "complete",
                       sprintf("Legacy ULM/MLM/GSVA and GTRD/ChEA/ENCODE GSEA completed for %d comparisons", length(contrasts)))
}
