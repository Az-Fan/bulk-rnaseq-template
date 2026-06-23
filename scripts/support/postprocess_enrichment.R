options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})
source("scripts/support/project_plot.R")

output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results")
functional_dir <- file.path(output_dir, "functional_enrichment")
top_n <- as.integer(Sys.getenv("ORA_PLOT_TOP_N", unset = "10"))

ratio_numeric <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
    if (length(z) != 2) return(NA_real_)
    as.numeric(z[[1]]) / as.numeric(z[[2]])
  }, numeric(1))
}

reduce_gene_overlap <- function(df, cutoff = 0.7) {
  if (!"geneID" %in% names(df)) return(df)
  out <- list()
  for (cluster in unique(df$Cluster)) {
    x <- df[df$Cluster == cluster & !is.na(df$p.adjust), , drop = FALSE]
    x <- x[order(x$p.adjust), , drop = FALSE]
    kept <- integer(0)
    sets <- list()
    for (i in seq_len(nrow(x))) {
      genes <- unique(strsplit(as.character(x$geneID[i]), "/", fixed = TRUE)[[1]])
      overlaps <- vapply(sets, function(old) {
        length(intersect(genes, old)) / max(1, length(union(genes, old)))
      }, numeric(1))
      if (!length(overlaps) || max(overlaps) < cutoff) {
        kept <- c(kept, i)
        sets[[length(sets) + 1]] <- genes
      }
    }
    out[[cluster]] <- x[kept, , drop = FALSE]
  }
  bind_rows(out)
}

plot_ora_source <- function(path, source_name, database_label) {
  if (!file.exists(path)) return()
  df <- read.csv(path, check.names = FALSE)
  if (!nrow(df)) return()
  df$GeneRatioNumeric <- ratio_numeric(df$GeneRatio)
  df$Description <- sub("^Gobp\\s+", "", df$Description, ignore.case = TRUE)
  # Pre-filter: keep only significant terms to avoid O(n^2) slowdown
  df <- df[!is.na(df$p.adjust) & df$p.adjust < 0.05, , drop = FALSE]
  if (!nrow(df)) return()
  reduced <- reduce_gene_overlap(df, cutoff = 0.7)
  write.csv(
    reduced,
    file.path(dirname(path), paste0("Enrichment_Reduced_", source_name, ".csv")),
    row.names = FALSE
  )

  for (direction in intersect(c("Up", "Down"), unique(reduced$Cluster))) {
    x <- reduced %>%
      filter(Cluster == direction, !is.na(p.adjust), p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      slice_head(n = top_n) %>%
      mutate(
        score = -log10(pmax(p.adjust, .Machine$double.xmin)),
        Description = factor(Description, levels = rev(Description))
      )
    if (!nrow(x)) next

    if (plot_enabled("ora.lollipop")) {
      colors <- project_colors()
      p <- ggplot(x, aes(score, Description)) +
        geom_segment(aes(x = 0, xend = score, yend = Description), color = "grey75") +
        geom_point(
          aes(size = Count, color = GeneRatioNumeric),
          alpha = 0.95
        ) +
        scale_color_gradient(low = "#91BFDB", high = if (direction == "Up") colors[["up"]] else colors[["down"]]) +
        labs(
          title = paste(source_name, direction, "ORA"),
          subtitle = "Gene-overlap-reduced representative terms",
          x = expression(-log[10]("adjusted P")),
          y = NULL,
          color = "Gene ratio",
          size = "Genes"
        ) +
        project_theme()
      out <- file.path(dirname(path), paste0(source_name, "_Lollipop_", direction, ".png"))
      save_project_plot(p, out, width = 8, height = max(4.5, nrow(x) * 0.42 + 1.6))
      write_figure_note(
        out,
        paste(source_name, direction, "ORA lollipop"),
        basename(path),
        "Over-representation analysis; representative terms selected by gene-set Jaccard overlap < 0.7",
        c("adjusted P < 0.05", paste0("top ", top_n, " terms")),
        database_label
      )
    }
  }

  if (plot_enabled("ora.diverging_bar")) {
    x <- reduced %>%
      filter(Cluster %in% c("Up", "Down"), !is.na(p.adjust), p.adjust < 0.05) %>%
      group_by(Cluster) %>%
      arrange(p.adjust, .by_group = TRUE) %>%
      slice_head(n = min(8, top_n)) %>%
      ungroup() %>%
      mutate(
        signed_score = ifelse(Cluster == "Up", 1, -1) *
          -log10(pmax(p.adjust, .Machine$double.xmin)),
        label = paste0(ifelse(Cluster == "Up", "UP: ", "DOWN: "), Description),
        label = factor(label, levels = label[order(signed_score)])
      )
    if (nrow(x)) {
      colors <- project_colors()
      p <- ggplot(x, aes(signed_score, label, fill = Cluster)) +
        geom_col(width = 0.72) +
        geom_vline(xintercept = 0, color = "grey35") +
        scale_fill_manual(values = c(Up = colors[["up"]], Down = colors[["down"]])) +
        scale_x_continuous(labels = function(v) abs(v)) +
        labs(
          title = paste(source_name, "ORA direction comparison"),
          subtitle = "Bar length = -log10 adjusted P",
          x = expression(-log[10]("adjusted P")),
          y = NULL
        ) +
        project_theme() +
        theme(legend.position = "none")
      out <- file.path(dirname(path), paste0(source_name, "_Diverging_Bar.png"))
      save_project_plot(p, out, width = 9, height = max(5, nrow(x) * 0.34 + 1.6))
      write_figure_note(
        out,
        paste(source_name, "ORA direction comparison"),
        basename(path),
        "Direction-separated ORA with gene-overlap redundancy reduction",
        c("adjusted P < 0.05", "up and down analysed separately"),
        database_label
      )
    }
  }
}

plot_ora_source(
  file.path(functional_dir, "go_analysis", "Enrichment_Full_GO_BP.csv"),
  "GO_BP", "MSigDB C5 GO:BP"
)
plot_ora_source(
  file.path(functional_dir, "kegg_analysis", "Enrichment_Full_KEGG.csv"),
  "KEGG", "MSigDB C2 CP:KEGG"
)
plot_ora_source(
  file.path(functional_dir, "hallmark_analysis", "Enrichment_Full_Hallmark.csv"),
  "Hallmark", "MSigDB Hallmark"
)

if (plot_enabled("gsea.nes_overview")) {
  gsea_dir <- file.path(functional_dir, "gsea_analysis")
  table_dir <- file.path(gsea_dir, "tables")
  for (db in c("Hallmark", "KEGG", "Reactome", "GO_BP")) {
    path <- file.path(table_dir, paste0("GSEA_Full_Table_", db, ".csv"))
    if (!file.exists(path)) next
    x <- read.csv(path, check.names = FALSE) %>%
      filter(!is.na(p.adjust), p.adjust < 0.25) %>%
      arrange(p.adjust) %>%
      slice_head(n = 15) %>%
      mutate(Description = factor(Description, levels = Description[order(NES)]))
    if (!nrow(x)) next
    colors <- project_colors()
    p <- ggplot(x, aes(NES, Description)) +
      geom_vline(xintercept = 0, color = "grey45") +
      geom_segment(aes(x = 0, xend = NES, yend = Description), color = "grey75") +
      geom_point(aes(size = setSize, color = p.adjust), alpha = 0.95) +
      scale_color_gradient(low = colors[["positive"]], high = colors[["negative"]], trans = "reverse") +
      labs(
        title = paste(db, "GSEA overview"),
        x = "Normalized enrichment score",
        y = NULL,
        color = "FDR",
        size = "Gene-set size"
      ) +
      project_theme()
    out_dir <- file.path(gsea_dir, paste0(db, "_plots"))
    out <- file.path(out_dir, paste0(db, "_NES_Overview.png"))
    save_project_plot(p, out, width = 8.5, height = max(5, nrow(x) * 0.36 + 1.6))
    write_figure_note(
      out,
      paste(db, "GSEA NES overview"),
      basename(path),
      "Preranked GSEA using the DESeq2 Wald statistic",
      c("FDR < 0.25", "top 15 by FDR"),
      paste("MSigDB", db)
    )
  }
}


