# Module 03 — legacy-compatible ORA and core GSEA
#
# Snakemake may schedule this module, but the scientific contract is inherited
# from legacy steps 04/04b/05/07/08. The frozen msigdbr 7.5.1 payload is data,
# not an installed R package. Formal analysis performs no network access.

run_enrichment <- function(context) {
  decision <- module_decision(context$config, "enrichment")
  if (!identical(decision$status, "enabled")) {
    return(record_module_status(context, "03_Enrichment", decision$status, decision$reason %||% ""))
  }
  if (is.null(context$state$differential)) stop("Differential expression must complete before enrichment")
  if (!identical(context$config$species, "human")) {
    stop("Legacy-compatible enrichment currently requires an explicitly frozen species resource")
  }

  resource_id <- context$config$resources$enrichment_gene_sets %||% NULL
  if (is.null(resource_id)) stop("Enabled enrichment requires resources.enrichment_gene_sets")
  payload <- readRDS(resource_path(context, resource_id))
  if (!is.list(payload) || !identical(payload$msigdb_release, "7.5.1") || is.null(payload$gene_sets)) {
    stop("Legacy-compatible enrichment requires the frozen MSigDB 7.5.1 human resource")
  }
  gene_sets <- payload$gene_sets
  required <- c("gs_cat", "gs_subcat", "gs_name", "gene_symbol", "entrez_gene")
  if (length(setdiff(required, names(gene_sets)))) stop("Frozen gene-set resource lacks required columns")
  gene_sets$entrez_gene <- as.character(gene_sets$entrez_gene)

  specifications <- list(
    GO_BP = c("C5", "GO:BP"), GO_CC = c("C5", "GO:CC"), GO_MF = c("C5", "GO:MF"),
    KEGG = c("C2", "CP:KEGG"), Reactome = c("C2", "CP:REACTOME"), Hallmark = c("H", "")
  )
  subset_db <- function(name, id_column) {
    spec <- specifications[[name]]
    keep <- gene_sets$gs_cat == spec[[1]]
    if (nzchar(spec[[2]])) keep <- keep & gene_sets$gs_subcat == spec[[2]]
    out <- unique(gene_sets[keep, c("gs_name", id_column), drop = FALSE])
    names(out) <- c("term", "gene")
    out <- out[!is.na(out$gene) & nzchar(out$gene), , drop = FALSE]
    out$gene <- as.character(out$gene)
    out
  }
  enrichment_config <- context$config$enrichment %||% list()
  ora_database_names <- as.character(enrichment_config$ora_databases %||% names(specifications))
  gsea_database_names <- as.character(enrichment_config$gsea_databases %||% c("Reactome", "Hallmark", "KEGG", "GO_BP"))
  unknown_databases <- setdiff(unique(c(ora_database_names, gsea_database_names)), names(specifications))
  if (length(unknown_databases)) stop("Unknown enrichment databases: ", paste(unknown_databases, collapse = ", "))
  ora_databases <- lapply(ora_database_names, subset_db, id_column = "entrez_gene")
  names(ora_databases) <- ora_database_names
  gsea_databases <- lapply(gsea_database_names, subset_db, id_column = "gene_symbol")
  names(gsea_databases) <- gsea_database_names

  min_size <- as.integer(enrichment_config$min_size %||% 10L)
  max_size <- as.integer(enrichment_config$max_size %||% 500L)
  profiles <- threshold_profiles(context$config)
  ora_profiles <- as.character(enrichment_config$ora_profiles %||% names(profiles)[[1L]])
  if (length(setdiff(ora_profiles, names(profiles)))) {
    stop("Unknown enrichment.ora_profiles: ", paste(setdiff(ora_profiles, names(profiles)), collapse = ", "))
  }
  ora_fdr <- as.numeric(enrichment_config$ora_fdr %||% 0.05)
  gsea_reporting_fdr <- as.numeric(enrichment_config$gsea_reporting_fdr %||% 0.25)
  redundancy_jaccard <- as.numeric(enrichment_config$redundancy_jaccard %||% 0.7)
  options(bulk_rnaseq.ora_fdr = ora_fdr, bulk_rnaseq.ora_jaccard = redundancy_jaccard)
  formal_seed <- as.integer(enrichment_config$random_seed %||% 123L)
  gsea_curves <- as.integer(enrichment_config$gsea_curves %||% 6L)
  results <- list()

  normalize_ora <- function(x, database) {
    if (!nrow(x)) return(data.frame())
    x$database <- database
    x$pathway <- x$ID
    x$direction <- as.character(x$Cluster)
    ratio_value <- function(value) vapply(strsplit(as.character(value), "/", fixed = TRUE), function(z) as.numeric(z[[1]]) / as.numeric(z[[2]]), numeric(1))
    x$gene_ratio <- ratio_value(x$GeneRatio)
    x$background_ratio <- ratio_value(x$BgRatio)
    x$fold_enrichment <- x$gene_ratio / x$background_ratio
    x$overlap_count <- x$Count
    x$query_size <- vapply(strsplit(as.character(x$GeneRatio), "/", fixed = TRUE), function(z) as.integer(z[[2]]), integer(1))
    x$term_size <- vapply(strsplit(as.character(x$BgRatio), "/", fixed = TRUE), function(z) as.integer(z[[1]]), integer(1))
    x$universe_size <- vapply(strsplit(as.character(x$BgRatio), "/", fixed = TRUE), function(z) as.integer(z[[2]]), integer(1))
    x$padj <- x$p.adjust
    x$overlap_genes <- x$geneID
    x
  }

  normalize_gsea <- function(x, database) {
    if (!nrow(x)) return(data.frame())
    x$database <- database
    x$pathway <- x$ID
    x$padj <- x$p.adjust
    x$size <- x$setSize
    x$ES <- x$enrichmentScore
    x$leadingEdge <- x$core_enrichment
    x
  }

  for (contrast_id in names(context$state$differential)) {
    de <- context$state$differential[[contrast_id]]$result
    root <- comparison_module_dir(context, contrast_id, "03_Enrichment")
    ora_root <- file.path(root, "ORA"); gsea_root <- file.path(root, "GSEA")
    for (path in c(file.path(ora_root, "Tables"), file.path(ora_root, "Figures"),
                   file.path(gsea_root, "Tables"), file.path(gsea_root, "Figures"))) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }

    universe <- unique(as.character(de$entrez_id[!is.na(de$entrez_id) & nzchar(de$entrez_id)]))
    write_result_table(data.frame(ENTREZID = universe), file.path(ora_root, "Tables", "ORA_Universe_Entrez.tsv"))
  profile_results <- list()
  # compareCluster resolves character function names in the caller search path.
  # Use a temporary environment on the search path so the legacy character
  # dispatch remains explicit without attaching the full package namespace.
  compare_env <- new.env(parent = emptyenv())
  compare_env$enricher <- clusterProfiler::enricher
  attach(compare_env, name = "bulk_rnaseq_clusterprofiler_dispatch")
  on.exit(detach("bulk_rnaseq_clusterprofiler_dispatch"), add = TRUE)
  for (profile_name in ora_profiles) {
      direction_column <- paste0("direction_", profile_name)
      up <- unique(as.character(de$entrez_id[de[[direction_column]] == "Up"]))
      down <- unique(as.character(de$entrez_id[de[[direction_column]] == "Down"]))
      up <- up[!is.na(up) & nzchar(up)]; down <- down[!is.na(down) & nzchar(down)]
      cluster <- list(Up = up, Down = down, All = unique(c(up, down)))
      database_results <- list()
      advanced_objects <- list()
      for (database in names(ora_databases)) {
        object <- clusterProfiler::compareCluster(
          cluster, fun = "enricher", TERM2GENE = ora_databases[[database]],
          universe = universe, pvalueCutoff = 1, qvalueCutoff = 1,
          minGSSize = min_size, maxGSSize = max_size,
          pAdjustMethod = "BH"
        )
        table <- normalize_ora(as.data.frame(object), database)
        database_results[[database]] <- table
        write_result_table(table, file.path(ora_root, "Tables", paste0("Enrichment_Full_", database, "_", profile_name, ".tsv")))
        if (database %in% c("GO_BP", "KEGG")) {
          advanced_objects[[database]] <- clusterProfiler::enricher(
            cluster$All, TERM2GENE = ora_databases[[database]], universe = universe,
            pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = min_size,
            maxGSSize = max_size, pAdjustMethod = "BH"
          )
        }
      }
      ora <- do.call(rbind, database_results)
      if (is.null(ora)) ora <- data.frame()
      rownames(ora) <- NULL
      profile_results[[profile_name]] <- ora
      write_result_table(ora, file.path(ora_root, "Tables", paste0("ORA_", profile_name, "_Full.tsv")))
      if (nrow(ora)) {
        plot <- plot_ora_publication(ora, paste0(contrast_id, " — ", profile_name, " ORA"),
                                     top_n = as.integer(enrichment_config$ora_overview_terms %||% 12L),
                                     fdr_cutoff = ora_fdr)
        save_publication_figure(plot, file.path(ora_root, "Figures", paste0("ORA_", profile_name, "_Directional")),
                                183, 125)
        sankey <- plot_ora_sankey_bubble(ora, paste0(contrast_id, " — ", profile_name, " ORA Sankey-bubble"))
        if (!is.null(sankey)) {
          save_publication_figure(sankey, file.path(ora_root, "Figures", paste0("ORA_", profile_name, "_Sankey_Bubble")),
                                  235, 205)
        }
        # Restore the legacy output families one-by-one. Each database gets
        # its own figures so a shared top-N quota can never hide a database.
        plot_top_n <- as.integer(enrichment_config$ora_overview_terms %||% 15L)
        for (database in names(database_results)) {
          database_table <- database_results[[database]]
          if (!nrow(database_table)) next
          figure_dir <- file.path(ora_root, "Figures", database)
          dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
          significant <- database_table[!is.na(database_table$padj) & database_table$padj < ora_fdr, , drop = FALSE]
          reduced <- reduce_ora_gene_overlap(significant, cutoff = redundancy_jaccard)
          write_result_table(reduced, file.path(ora_root, "Tables", paste0("Enrichment_Reduced_", database, "_", profile_name, ".tsv")))

          for (direction in intersect(c("All", "Down", "Up"), unique(database_table$direction))) {
            dot <- plot_ora_dotplot(database_table, database, direction,
                                    paste(database, direction, "ORA"), top_n = plot_top_n)
            if (!is.null(dot)) {
              save_publication_figure(dot, file.path(figure_dir, paste0(database, "_Dotplot_", direction)),
                                      183, max(105, 5.7 * plot_top_n + 45))
            }
          }
          for (direction in intersect(c("Down", "Up"), unique(database_table$direction))) {
            lollipop <- plot_ora_lollipop(database_table, database, direction,
                                          paste(database, direction, "ORA"), top_n = plot_top_n)
            if (!is.null(lollipop)) {
              save_publication_figure(lollipop, file.path(figure_dir, paste0(database, "_Lollipop_", direction)),
                                      183, max(105, 5.7 * plot_top_n + 45))
            }
          }
          diverging <- plot_ora_diverging(database_table, database,
                                          paste(database, "ORA direction comparison"), top_n = min(8L, plot_top_n))
          if (!is.null(diverging)) {
            save_publication_figure(diverging, file.path(figure_dir, paste0(database, "_Diverging_Bar")),
                                    183, 140)
          }
          database_sankey <- plot_ora_sankey_bubble(database_table,
                                                     paste(contrast_id, database, "ORA Sankey-bubble"), top_n = min(2L, plot_top_n))
          if (!is.null(database_sankey)) {
            save_publication_figure(database_sankey,
                                    file.path(ora_root, "Figures", "Sankey", paste0(database, "_Sankey_Bubble")),
                                    195, 145)
            write_result_table(database_sankey$data,
                               file.path(ora_root, "Tables", "Sankey", paste0(database, "_Sankey_Display_Data.tsv")))
            write_result_table(significant,
                               file.path(ora_root, "Tables", "Sankey", paste0(database, "_Filtered_Enrichment.tsv")))
          }
        }

        # cnetplot, heatplot and emapplot are distinct evidence views; aPEAR
        # must not stand in for any of them in file-level migration audits.
        fold_change <- de$log2FoldChange_ashr
        names(fold_change) <- as.character(de$entrez_id)
        fold_change <- fold_change[is.finite(fold_change) & !is.na(names(fold_change)) & nzchar(names(fold_change))]
        fold_change <- fold_change[!duplicated(names(fold_change))]
        fold_change_symbol <- de$log2FoldChange_ashr
        names(fold_change_symbol) <- as.character(de$gene_symbol)
        fold_change_symbol <- fold_change_symbol[
          is.finite(fold_change_symbol) & !is.na(names(fold_change_symbol)) & nzchar(names(fold_change_symbol))
        ]
        fold_change_symbol <- fold_change_symbol[!duplicated(names(fold_change_symbol))]
        for (database in names(advanced_objects)) {
          object <- advanced_objects[[database]]
          if (is.null(object)) next
          object@result <- object@result[!is.na(object@result$p.adjust) & object@result$p.adjust < ora_fdr, , drop = FALSE]
          n_terms <- nrow(object@result)
          if (!n_terms) next
          # Convert Entrez labels to symbols for display only. The enrichment
          # result, tested genes and statistics stay unchanged.
          readable_object <- tryCatch(
            clusterProfiler::setReadable(object, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID"),
            error = function(e) object
          )
          readable_object@result$Description <- clean_ora_term(readable_object@result$Description, width = 24L)
          advanced_dir <- file.path(ora_root, "Figures", "Advanced")
          set.seed(formal_seed)
          cnet <- enrichplot::cnetplot(
            readable_object, foldChange = fold_change_symbol,
            showCategory = min(5L, n_terms), node_label = "category",
            color_category = "#0072B2", color_item = "#BDBDBD",
            size_category = 1.15, size_item = 0.55, size_edge = 0.35
          ) +
            ggplot2::ggtitle(paste(database, "gene-concept network"))
          save_publication_figure(cnet, file.path(advanced_dir, paste0(database, "_Cnet")), 210, 175)
          heat <- enrichplot::heatplot(
            readable_object, foldChange = fold_change_symbol,
            showCategory = min(8L, n_terms)
          ) +
            ggplot2::ggtitle(paste(database, "gene-pathway heatmap")) +
            ggplot2::theme(
              axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, size = 5.2),
              axis.text.y = ggplot2::element_text(size = 7.5)
            )
          save_publication_figure(heat, file.path(advanced_dir, paste0(database, "_Heat")), 235, 150)
          if (n_terms > 1L) {
            similarity <- enrichplot::pairwise_termsim(readable_object, method = "JC")
            set.seed(formal_seed)
            emap <- enrichplot::emapplot(
              similarity, showCategory = min(12L, n_terms), layout = "kk",
              min_edge = 0.35, node_label_size = 3.2, label_format = 24,
              size_category = 0.8, size_edge = 0.35
            ) +
              ggplot2::ggtitle(paste(database, "enrichment map"))
            save_publication_figure(emap, file.path(advanced_dir, paste0(database, "_Emap")), 183, 170)
          }
        }

        # aPEAR is the pinned local Pixi dependency. It clusters only the
        # already-saved significant ORA rows; formal ORA results are unchanged.
        network_input <- ora[ora$direction %in% c("Up", "Down") & !is.na(ora$p.adjust) & ora$p.adjust < ora_fdr, , drop = FALSE]
        if (nrow(network_input) >= 3L) {
          network_input <- reduce_ora_gene_overlap(network_input, cutoff = redundancy_jaccard)
          network_input <- do.call(rbind, lapply(
            split(network_input, interaction(network_input$database, network_input$direction, drop = TRUE)),
            function(value) head(value[order(value$p.adjust, -value$fold_enrichment), , drop = FALSE], 3L)
          ))
          network_input$Description <- paste(
            network_input$database, network_input$direction,
            clean_ora_term(network_input$pathway, width = 24L), sep = ": "
          )
          write_result_table(network_input, file.path(ora_root, "Tables", "aPEAR_Display_Terms.tsv"))
          set.seed(formal_seed)
          apear_plot <- tryCatch(aPEAR::enrichmentNetwork(network_input), error = function(e) NULL)
          if (!is.null(apear_plot)) {
            save_publication_figure(apear_plot, file.path(ora_root, "Figures", paste0("ORA_", profile_name, "_aPEAR")), 235, 175)
          }
        }
      }
    }

    # Legacy step 05 removed duplicated symbols in the original count-matrix
    # order before sorting the statistic. Do not select the most significant
    # Ensembl row for a duplicated symbol: that changes the ranked list.
    de_input_order <- de[match(rownames(context$state$dds), de$gene_id), , drop = FALSE]
    de_input_order$ranking_metric <- ifelse(
      is.finite(de_input_order$stat), de_input_order$stat,
      ifelse(is.finite(de_input_order$log2FoldChange_raw) & is.finite(de_input_order$pvalue),
             sign(de_input_order$log2FoldChange_raw) *
               -log10(pmax(de_input_order$pvalue, .Machine$double.xmin)), NA_real_)
    )
    ranked <- de_input_order[is.finite(de_input_order$ranking_metric) &
                               !is.na(de_input_order$gene_symbol) & nzchar(de_input_order$gene_symbol),
                             c("gene_symbol", "ranking_metric")]
    names(ranked)[[2]] <- "stat"
    ranked <- ranked[!duplicated(ranked$gene_symbol), , drop = FALSE]
    ranked <- ranked[order(-ranked$stat), , drop = FALSE]
    ranks <- ranked$stat; names(ranks) <- ranked$gene_symbol
    write_result_table(data.frame(gene_symbol = names(ranks), wald_statistic = unname(ranks)),
                       file.path(gsea_root, "Tables", "GSEA_Ranked_Gene_List.tsv.gz"))

    database_gsea <- list(); pathway_lists <- list()
    for (database in names(gsea_databases)) {
      set.seed(formal_seed)
      object <- clusterProfiler::GSEA(
        ranks, TERM2GENE = gsea_databases[[database]], pvalueCutoff = 1,
        pAdjustMethod = "BH", minGSSize = min_size, maxGSSize = max_size,
        seed = TRUE, by = "fgsea", BPPARAM = BiocParallel::SerialParam(), verbose = FALSE
      )
      table <- normalize_gsea(as.data.frame(object), database)
      database_gsea[[database]] <- table
      pathway_lists[[database]] <- split(gsea_databases[[database]]$gene, gsea_databases[[database]]$term)
      write_result_table(table, file.path(gsea_root, "Tables", paste0("GSEA_Full_Table_", database, ".tsv")))
      if (nrow(table)) {
        overview <- plot_gsea_overview(table, paste0(contrast_id, " — ", database, " GSEA"),
                                       top_n = as.integer(enrichment_config$overview_terms %||% 20L),
                                       fdr_cutoff = gsea_reporting_fdr)
        if (!is.null(overview)) {
          save_publication_figure(overview, file.path(gsea_root, "Figures", database, paste0(database, "_NES_Overview")), 178, 130)
        }
        reportable <- table[!is.na(table$p.adjust) & table$p.adjust <= gsea_reporting_fdr, , drop = FALSE]
        curve_rows <- head(reportable[order(reportable$p.adjust, -abs(reportable$NES)), , drop = FALSE],
                           min(gsea_curves, nrow(reportable)))
        for (j in seq_len(nrow(curve_rows))) {
          summary_row <- curve_rows[j, , drop = FALSE]
          curve <- plot_gsea_curve(pathway_lists[[database]][[summary_row$ID]], ranks, summary_row, contrast_id)
          stem <- sprintf("%02d_%s", j, safe_id(gsub("[^A-Za-z0-9_.-]+", "_", summary_row$ID), "pathway"))
          save_publication_figure(curve$plot, file.path(gsea_root, "Figures", database, stem), 178, 140)
          write_result_table(curve$running, file.path(gsea_root, "Tables", "Curve_Data", database, paste0(stem, "_Running_ES.tsv.gz")))
          write_result_table(curve$hits, file.path(gsea_root, "Tables", "Curve_Data", database, paste0(stem, "_Hits.tsv.gz")))
        }
        # Preserve legacy step 07 as an automatic, reproducible core-gene
        # heatmap for the leading significant pathway of each database.
        heat_candidates <- table[!is.na(table$p.adjust) & table$p.adjust <= gsea_reporting_fdr, , drop = FALSE]
        heat_candidates <- heat_candidates[order(heat_candidates$p.adjust, -abs(heat_candidates$NES)), , drop = FALSE]
        if (nrow(heat_candidates)) {
          heat_row <- heat_candidates[1, , drop = FALSE]
          core_genes <- strsplit(heat_row$core_enrichment, "/", fixed = TRUE)[[1]]
          map <- context$state$annotation[rownames(context$state$vst), c("gene_id", "gene_symbol"), drop = FALSE]
          valid <- !is.na(map$gene_symbol) & nzchar(map$gene_symbol)
          frame <- data.frame(gene_symbol = map$gene_symbol[valid], context$state$vst[valid, , drop = FALSE], check.names = FALSE)
          frame$average <- rowMeans(frame[, -1, drop = FALSE]); frame <- frame[order(-frame$average), ]
          frame <- frame[!duplicated(frame$gene_symbol), ]; rownames(frame) <- frame$gene_symbol
          matched <- intersect(core_genes, rownames(frame))
          if (length(matched) >= 2L) {
            matrix <- as.matrix(frame[matched, context$state$selected_samples, drop = FALSE])
            annotation_col <- data.frame(group = context$state$samples[context$state$selected_samples, context$state$display_factor])
            rownames(annotation_col) <- context$state$selected_samples
            heat <- pheatmap::pheatmap(matrix, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
                                      annotation_col = annotation_col, show_colnames = FALSE, border_color = NA, silent = TRUE,
                                      color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                                      main = paste0(database, " core enrichment: ", heat_row$Description))$gtable
            save_publication_grob(heat, file.path(gsea_root, "Figures", database, "Core_Enrichment_Heatmap"),
                                  178, max(110, min(260, 4 * length(matched) + 55)))
            write_result_table(data.frame(gene_symbol = rownames(matrix), matrix, check.names = FALSE),
                               file.path(gsea_root, "Tables", "Heatmap_Data", paste0(database, "_Core_Enrichment.tsv.gz")))
          }
        }
      }
    }
    gsea <- do.call(rbind, database_gsea)
    if (is.null(gsea)) gsea <- data.frame()
    rownames(gsea) <- NULL
    write_result_table(gsea, file.path(gsea_root, "Tables", "GSEA_Full.tsv"))
    write_result_table(data.frame(
      analysis = "ORA/GSEA", resource_id = resource_id, msigdb_release = payload$msigdb_release,
      seed = formal_seed, gsea_rank = "DESeq2 Wald statistic", ora_id = "ENTREZID",
      ora_profiles = paste(ora_profiles, collapse = ","), ora_fdr = ora_fdr,
      gsea_reporting_fdr = gsea_reporting_fdr,
      ora_universe = "all tested mapped genes", stringsAsFactors = FALSE
    ), file.path(root, "Gene_Set_Provenance.tsv"))
    results[[contrast_id]] <- list(ora = profile_results, gsea = database_gsea,
                                   ranks = ranks, pathway_list = pathway_lists, directory = root)
  }

  context$state$enrichment <- results
  record_module_status(context, "03_Enrichment", "complete",
                       sprintf("Legacy-compatible ORA/GSEA completed for %d contrasts using %s", length(results), resource_id))
}
