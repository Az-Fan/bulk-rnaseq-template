# ==============================================================================
# Step04b ORA sankey-bubble plots
# Depends on Step01 DEG output and Step04 ORA tables.
# ==============================================================================

rm(list = ls())

output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results/final")
shared_cache_dir <- Sys.getenv("SHARED_CACHE_DIR", unset = "shared_cache")
local_r_lib <- file.path(shared_cache_dir, "R_libs")
dir.create(local_r_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(normalizePath(local_r_lib, mustWork = FALSE), .libPaths()))
if (Sys.getenv("XDG_DATA_HOME", unset = "") == "") {
  Sys.setenv(XDG_DATA_HOME = file.path(shared_cache_dir, "xdg_data"))
}
diff_analysis_dir <- file.path(output_dir, "diff_analysis")
functional_enrichment_dir <- file.path(output_dir, "functional_enrichment")
sankey_bubble_dir <- file.path(functional_enrichment_dir, "sankey_bubble")
dir.create(sankey_bubble_dir, showWarnings = FALSE, recursive = TRUE)

required_pkgs <- c(
  "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "dplyr", "tidyr",
  "purrr", "stringr", "ggplot2", "ggalluvial", "patchwork", "RColorBrewer",
  "Cairo", "DOSE"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
  library(patchwork)
  library(RColorBrewer)
  library(Cairo)
  library(DOSE)
})

has_reactome <- requireNamespace("ReactomePA", quietly = TRUE) &&
  requireNamespace("reactome.db", quietly = TRUE)

pvalue_filter <- as.numeric(Sys.getenv("ORA_SANKEY_PVALUE", unset = "0.05"))
padj_filter <- as.numeric(Sys.getenv("ORA_SANKEY_PADJ", unset = "0.05"))
fc_cutoff <- as.numeric(Sys.getenv("ORA_SANKEY_LOGFC", unset = "1"))
top_pathways <- as.integer(Sys.getenv("ORA_SANKEY_TOP_N", unset = "20"))

plot_width <- 13
plot_height <- 11
single_plot_min_width <- 15
single_plot_max_height <- 42
sankey_width_ratio <- 3
dot_width_ratio <- 2
plot_font <- "sans"
sankey_text_size <- 3.2
bubble_low_color <- "#2166ac"
bubble_high_color <- "#b2182b"
bubble_size_range <- c(2.5, 8)
pathway_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
)

load(file = file.path(diff_analysis_dir, "DEseq2_Diff_Annotated.Rdata"))
res_df <- as.data.frame(res_annotated)

if (!"ENTREZID" %in% colnames(res_df) && "entrez" %in% colnames(res_df)) {
  res_df$ENTREZID <- res_df$entrez
}

if (!"ENTREZID" %in% colnames(res_df)) {
  ids <- AnnotationDbi::select(
    x = org.Hs.eg.db,
    keys = unique(res_df$gene[!is.na(res_df$gene) & res_df$gene != ""]),
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "SYMBOL",
    skipValidKeysTest = TRUE
  )
  ids <- ids[!is.na(ids$ENTREZID) & ids$ENTREZID != "", c("SYMBOL", "ENTREZID")]
  ids <- ids[!duplicated(ids$SYMBOL), , drop = FALSE]
  res_df <- dplyr::left_join(res_df, ids, by = c("gene" = "SYMBOL"))
}

res_df$ENTREZID <- as.character(res_df$ENTREZID)
gene_mapping <- res_df %>%
  dplyr::filter(!is.na(.data$gene), .data$gene != "", !is.na(.data$ENTREZID), .data$ENTREZID != "") %>%
  dplyr::distinct(.data$ENTREZID, .keep_all = TRUE) %>%
  dplyr::transmute(entrez = as.character(.data$ENTREZID), symbol = as.character(.data$gene))

diff_df <- res_df %>%
  dplyr::filter(
    !is.na(.data$ENTREZID),
    .data$ENTREZID != "",
    .data$adj.P.Val < padj_filter,
    abs(.data$logFC) > fc_cutoff
  )
gene_all <- unique(as.character(diff_df$ENTREZID))
gene_all <- gene_all[grepl("^[0-9]+$", gene_all)]

cat("Sankey-bubble DEG Entrez IDs: ", length(gene_all), "\n", sep = "")
if (length(gene_all) == 0) {
  stop("No DEG Entrez IDs passed filters; cannot build sankey-bubble plots.")
}

map_gene_ids_to_symbols <- function(gene_id_string) {
  ids <- strsplit(as.character(gene_id_string), "/", fixed = TRUE)[[1]]
  mapped <- gene_mapping$symbol[match(ids, gene_mapping$entrez)]
  mapped[is.na(mapped) | mapped == ""] <- ids[is.na(mapped) | mapped == ""]
  paste(mapped, collapse = "/")
}

parse_ratio <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(parts) {
    if (length(parts) != 2) return(NA_real_)
    as.numeric(parts[1]) / as.numeric(parts[2])
  }, numeric(1))
}

read_step04_table <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if ("Cluster" %in% colnames(df) && any(df$Cluster == "All")) {
    df <- df[df$Cluster == "All", , drop = FALSE]
  }
  df
}

filter_enrich_df <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  needed <- c("Description", "GeneRatio", "pvalue", "p.adjust", "geneID", "Count")
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0) {
    warning("Skipping enrichment table with missing columns: ", paste(missing, collapse = ", "))
    return(NULL)
  }
  df <- df %>%
    dplyr::filter(!is.na(.data$pvalue), !is.na(.data$p.adjust)) %>%
    dplyr::filter(.data$pvalue < pvalue_filter, .data$p.adjust < padj_filter)
  if (nrow(df) == 0) return(NULL)
  df$geneID <- vapply(df$geneID, map_gene_ids_to_symbols, character(1))
  df
}

run_do_enrichment <- function() {
  obj <- tryCatch(
    DOSE::enrichDO(
      gene = gene_all,
      ont = "HDO",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ),
    error = function(e) {
      message("DO enrichment skipped: ", e$message)
      NULL
    }
  )
  if (is.null(obj)) {
    cached_do <- file.path(sankey_bubble_dir, "DO", "DO_filtered_enrichment.csv")
    if (file.exists(cached_do)) {
      message("DO enrichment fallback: using cached table ", cached_do)
      return(read.csv(cached_do, check.names = FALSE, stringsAsFactors = FALSE))
    }
    return(NULL)
  }
  as.data.frame(obj)
}

run_reactome_enrichment <- function() {
  if (!has_reactome) {
    message("Reactome skipped: ReactomePA/reactome.db is not installed.")
    return(NULL)
  }
  obj <- tryCatch(
    ReactomePA::enrichPathway(
      gene = gene_all,
      organism = "human",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    ),
    error = function(e) {
      message("Reactome enrichment skipped: ", e$message)
      NULL
    }
  )
  if (is.null(obj)) return(NULL)
  obj <- tryCatch(
    clusterProfiler::setReadable(obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
    error = function(e) obj
  )
  as.data.frame(obj)
}

make_palette <- function(values, base_colors) {
  values <- unique(as.character(values))
  if (length(values) == 0) return(character(0))
  colors <- if (length(values) <= length(base_colors)) {
    base_colors[seq_along(values)]
  } else {
    grDevices::colorRampPalette(base_colors)(length(values))
  }
  stats::setNames(colors, values)
}

create_sankey_bubble <- function(enrich_df, analysis_name, top_n = top_pathways) {
  enrich_df <- filter_enrich_df(enrich_df)
  if (is.null(enrich_df) || nrow(enrich_df) == 0) {
    cat(analysis_name, ": no significant result after filtering.\n", sep = "")
    return(NULL)
  }

  analysis_dir <- file.path(sankey_bubble_dir, analysis_name)
  dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(enrich_df, file.path(analysis_dir, paste0(analysis_name, "_filtered_enrichment.csv")), row.names = FALSE)

  enrich_top <- enrich_df %>%
    dplyr::arrange(.data$pvalue) %>%
    dplyr::slice_head(n = min(top_n, nrow(.))) %>%
    dplyr::mutate(
      GeneRatio_numeric = parse_ratio(.data$GeneRatio),
      Count = as.numeric(.data$Count),
      Description = stringr::str_wrap(as.character(.data$Description), width = 42),
      term_order = dplyr::row_number()
    )

  link_data <- enrich_top %>%
    tidyr::separate_rows(.data$geneID, sep = "/") %>%
    dplyr::rename(Gene = .data$geneID, Pathway = .data$Description) %>%
    dplyr::filter(!is.na(.data$Gene), .data$Gene != "", !is.na(.data$Pathway), .data$Pathway != "") %>%
    dplyr::mutate(link_id = dplyr::row_number(), y = 1)

  if (nrow(link_data) == 0) {
    cat(analysis_name, ": no gene-term links to plot.\n", sep = "")
    return(NULL)
  }

  lodes <- ggalluvial::to_lodes_form(
    link_data,
    axes = c("Gene", "Pathway"),
    key = "axis",
    value = "stratum",
    id = "link_id"
  ) %>%
    dplyr::left_join(link_data %>% dplyr::select("link_id", "Pathway"), by = "link_id")

  gene_values <- unique(link_data$Gene)
  pathway_values <- unique(enrich_top$Description)
  gene_colors <- make_palette(gene_values, RColorBrewer::brewer.pal(12, "Set3"))
  pathway_palette <- make_palette(pathway_values, pathway_colors)
  node_colors <- c(gene_colors, pathway_palette)

  base_theme <- ggplot2::theme(
    text = ggplot2::element_text(family = plot_font, face = "bold"),
    plot.margin = grid::unit(c(0.1, 0.1, 0.1, 0.1), "cm")
  )

  sankey_plot <- ggplot2::ggplot(
    lodes,
    ggplot2::aes(
      x = .data$axis,
      stratum = .data$stratum,
      alluvium = .data$link_id,
      y = .data$y
    )
  ) +
    ggalluvial::geom_flow(
      ggplot2::aes(fill = .data$Pathway),
      alpha = 0.28,
      width = 0.08,
      knot.pos = 0.35,
      color = "transparent"
    ) +
    ggalluvial::geom_stratum(
      ggplot2::aes(fill = .data$stratum),
      color = "grey35",
      linewidth = 0.15,
      width = 0.08
    ) +
    ggplot2::geom_text(
      stat = "stratum",
      ggplot2::aes(label = after_stat(stratum)),
      size = sankey_text_size,
      family = plot_font,
      fontface = "bold",
      hjust = 1,
      nudge_x = -0.055
    ) +
    ggplot2::scale_fill_manual(values = c(node_colors, pathway_palette), guide = "none") +
    ggplot2::scale_x_discrete(expand = c(0.28, 0.02)) +
    ggplot2::labs(x = paste0("Gene-Term relationship (", analysis_name, ")"), y = NULL) +
    ggplot2::theme_void() +
    base_theme +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8), size = 11),
      plot.margin = ggplot2::margin(5.5, 0, 5.5, 5.5)
    )

  sankey_build <- ggplot2::ggplot_build(sankey_plot)
  stratum_layer <- sankey_build$data[[2]]
  right_nodes <- stratum_layer %>%
    dplyr::filter(.data$x == max(.data$x)) %>%
    dplyr::mutate(
      node_name = as.character(.data$stratum),
      node_center_y = (.data$ymin + .data$ymax) / 2
    ) %>%
    dplyr::filter(.data$node_name %in% pathway_values) %>%
    dplyr::select("node_name", "node_center_y", "ymin", "ymax")

  bubble_data <- enrich_top %>%
    dplyr::left_join(right_nodes, by = c("Description" = "node_name"))

  bubble_plot <- ggplot2::ggplot(
    bubble_data,
    ggplot2::aes(x = .data$GeneRatio_numeric, y = .data$node_center_y, color = -log10(.data$pvalue))
  ) +
    ggplot2::geom_point(ggplot2::aes(size = .data$Count), alpha = 0.9, stroke = 0.4) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_color_gradient(low = bubble_low_color, high = bubble_high_color, name = "-log10(Pvalue)") +
    ggplot2::scale_radius(range = bubble_size_range, name = "Count") +
    ggplot2::labs(x = "GeneRatio", y = NULL) +
    ggplot2::theme_bw() +
    base_theme +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 4), size = 11),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8), size = 11),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      legend.position = "right",
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 0)
    )

  y_range <- sankey_build$layout$panel_params[[1]]$y.range
  combined_plot <- (
    sankey_plot + ggplot2::coord_cartesian(clip = "off", ylim = y_range)
  ) + (
    bubble_plot +
      ggplot2::annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = min(right_nodes$ymin),
        ymax = max(right_nodes$ymax),
        fill = NA,
        color = "black",
        linewidth = 0.3
      ) +
      ggplot2::coord_cartesian(ylim = y_range)
  ) +
    patchwork::plot_layout(widths = c(sankey_width_ratio, dot_width_ratio))

  pdf_path <- file.path(analysis_dir, paste0(analysis_name, "_sankey_bubble.pdf"))
  png_path <- file.path(analysis_dir, paste0(analysis_name, "_sankey_bubble.png"))
  plot_width_local <- max(plot_width, single_plot_min_width)
  plot_height_local <- max(
    plot_height,
    min(single_plot_max_height, length(gene_values) * 0.16 + 3),
    min(single_plot_max_height, length(pathway_values) * 0.65 + 4)
  )

  ggplot2::ggsave(
    pdf_path,
    combined_plot,
    width = plot_width_local,
    height = plot_height_local,
    device = grDevices::cairo_pdf,
    limitsize = FALSE
  )
  ggplot2::ggsave(
    png_path,
    combined_plot,
    width = plot_width_local,
    height = plot_height_local,
    dpi = 300,
    limitsize = FALSE
  )

  cat(analysis_name, " sankey-bubble saved: ", pdf_path, "\n", sep = "")
  list(plot = combined_plot, width = plot_width_local, height = plot_height_local)
}

tables <- list(
  KEGG = read_step04_table(file.path(functional_enrichment_dir, "kegg_analysis", "Enrichment_Full_KEGG.csv")),
  GO_BP = read_step04_table(file.path(functional_enrichment_dir, "go_analysis", "Enrichment_Full_GO_BP.csv")),
  GO_CC = read_step04_table(file.path(functional_enrichment_dir, "go_analysis", "Enrichment_Full_GO_CC.csv")),
  GO_MF = read_step04_table(file.path(functional_enrichment_dir, "go_analysis", "Enrichment_Full_GO_MF.csv")),
  DO = run_do_enrichment(),
  Reactome = run_reactome_enrichment()
)

plots <- purrr::imap(tables, function(df, nm) {
  tryCatch(
    create_sankey_bubble(df, nm),
    error = function(e) {
      message(nm, " sankey-bubble skipped: ", e$message)
      NULL
    }
  )
})
plots <- plots[!vapply(plots, is.null, logical(1))]

if (length(plots) >= 2) {
  plot_objects <- lapply(plots, `[[`, "plot")
  while (length(plot_objects) < 6) {
    plot_objects[[paste0("No_data_", length(plot_objects) + 1)]] <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = "No data",
        size = 8,
        color = "grey60",
        fontface = "bold"
      ) +
      ggplot2::theme_void()
  }
  plot_objects <- plot_objects[seq_len(6)]

  plot_tiles <- lapply(plot_objects, patchwork::wrap_elements)
  combined <- patchwork::wrap_plots(plot_tiles, ncol = 3) +
    patchwork::plot_annotation(
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.tag = ggplot2::element_text(size = 18, face = "bold", family = plot_font),
        plot.tag.position = c(0.01, 0.99)
      )
    )

  out_pdf <- file.path(sankey_bubble_dir, "Combined_sankey_bubble_2x3.pdf")
  out_png <- file.path(sankey_bubble_dir, "Combined_sankey_bubble_2x3.png")
  out_tiff <- file.path(sankey_bubble_dir, "Combined_sankey_bubble_2x3.tiff")
  ggplot2::ggsave(out_pdf, combined, width = 30, height = 20, device = grDevices::cairo_pdf, limitsize = FALSE)
  ggplot2::ggsave(out_png, combined, width = 30, height = 20, dpi = 150, limitsize = FALSE)
  ggplot2::ggsave(
    out_tiff,
    combined,
    width = 30,
    height = 20,
    dpi = 200,
    device = "tiff",
    compression = "lzw",
    limitsize = FALSE
  )
  cat("Combined sankey-bubble saved: ", out_pdf, "\n", sep = "")

  multipage_pdf <- file.path(sankey_bubble_dir, "Combined_sankey_bubble_multipage.pdf")
  grDevices::cairo_pdf(multipage_pdf, width = 18, height = 24, onefile = TRUE)
  for (nm in names(plots)) {
    page_plot <- plots[[nm]]$plot +
      patchwork::plot_annotation(
        title = nm,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5)
        )
      )
    print(page_plot)
  }
  grDevices::dev.off()
  cat("Multipage sankey-bubble saved: ", multipage_pdf, "\n", sep = "")
} else {
  cat("Fewer than two sankey-bubble plots were generated; combined 2x3 plot skipped.\n")
}

cat("Step04b sankey-bubble finished.\n")
