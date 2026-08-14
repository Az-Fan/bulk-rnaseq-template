# Shared scientific helpers and publication figure layer
#
# Statistical decisions stay in the numbered modules. This file contains only
# reusable validation, I/O, plotting, resource lookup and export primitives.

`%||%` <- function(x, fallback) if (is.null(x) || length(x) == 0L) fallback else x

workflow_root <- function(candidate = getwd()) {
  candidate <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(candidate, "pixi.toml"))) {
    stop("The workflow root must contain pixi.toml: ", candidate)
  }
  candidate
}

as_flag <- function(x, field) {
  value <- tolower(trimws(as.character(x)))
  if (any(!value %in% c("true", "false"))) stop(field, " must contain only true/false")
  value == "true"
}

read_project_table <- function(path) {
  if (!file.exists(path)) stop("Required project table is missing: ", path)
  extension <- tolower(tools::file_ext(sub("\\.gz$", "", path)))
  signature <- if (!grepl("\\.gz$", path, ignore.case = TRUE)) readBin(path, "raw", n = 8L) else raw()
  is_ole_xls <- length(signature) >= 8L && identical(signature[1:8], as.raw(c(0xD0, 0xCF, 0x11, 0xE0, 0xA1, 0xB1, 0x1A, 0xE1)))
  is_zip_xlsx <- length(signature) >= 2L && identical(signature[1:2], charToRaw("PK"))
  if (extension %in% c("xlsx", "xls") && (is_ole_xls || is_zip_xlsx)) {
    return(as.data.frame(readxl::read_excel(path), check.names = FALSE, stringsAsFactors = FALSE))
  }
  separator <- if (extension == "csv") "," else "\t"
  data.table::fread(path, sep = separator, data.table = FALSE, check.names = FALSE)
}

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

write_result_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", na = "")
  invisible(path)
}

write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(path)
}

safe_id <- function(x, field = "identifier") {
  x <- as.character(x)
  if (any(!grepl("^[A-Za-z0-9_.-]+$", x))) stop(field, " contains an unsafe value")
  x
}

module_decision <- function(config, key) {
  decision <- config$analysis$modules[[key]]
  if (is.null(decision)) stop("project.yml is missing an explicit module decision: ", key)
  if (!isTRUE(decision$confirmed)) stop("Module decision is not confirmed: ", key)
  decision
}

module_enabled <- function(context, key) identical(module_decision(context$config, key)$status, "enabled")

module_confirmation_string <- function(config) {
  keys <- sort(names(config$analysis$modules))
  paste(vapply(keys, function(key) {
    decision <- module_decision(config, key)
    paste0(key, "=", decision$status)
  }, character(1)), collapse = ",")
}

plan_confirmation_string <- function(config, paths, project_file) {
  ordered <- c(project_file, paths$samples, paths$contrasts, paths$approval, paths$source_manifest)
  if (any(!file.exists(ordered))) stop("Cannot create a run-plan token while required decision inputs are missing")
  material <- c(config$project_id, vapply(ordered, sha256_file, character(1)))
  paste0("PLAN:", config$project_id, ":",
         substr(digest::digest(paste(material, collapse = "|"), algo = "sha256", serialize = FALSE), 1L, 24L))
}

canonical_module_state <- function(config, keys) {
  decisions <- lapply(keys, function(key) module_decision(config, key))
  enabled <- vapply(decisions, function(x) identical(x$status, "enabled"), logical(1))
  if (any(enabled)) return(list(status = "enabled", reason = ""))
  states <- vapply(decisions, `[[`, character(1), "status")
  reasons <- vapply(decisions, function(x) x$reason %||% "", character(1))
  list(status = if (all(states == "not_applicable")) "not_applicable" else "skipped_by_user",
       reason = paste(unique(reasons[nzchar(reasons)]), collapse = "; "))
}

publication_palette <- function(levels) {
  colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
              "#E69F00", "#56B4E9", "#7A5195", "#000000")
  levels <- unique(as.character(levels))
  if (length(levels) > length(colors)) {
    stop("More than eight display groups require a project-declared accessible palette")
  }
  stats::setNames(colors[seq_along(levels)], levels)
}

direction_palette <- c(Down = "#0072B2", Not_significant = "#BDBDBD", Up = "#D55E00")

theme_publication <- function(base_size = 9, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2.2,
                                         colour = "#20252B", margin = ggplot2::margin(b = 3)),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.5, colour = "#5D646C",
                                            margin = ggplot2::margin(b = 7)),
      plot.caption = ggplot2::element_text(size = base_size - 1.5, colour = "#5D646C", hjust = 0),
      axis.title = ggplot2::element_text(colour = "#20252B"),
      axis.text = ggplot2::element_text(colour = "#20252B"),
      axis.line = ggplot2::element_line(colour = "#20252B", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(colour = "#20252B", linewidth = 0.35),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.tag = ggplot2::element_text(face = "bold")
    )
}

requested_figure_formats <- function() {
  formats <- getOption("bulk_rnaseq.figure_formats", "pdf")
  unique(tolower(as.character(formats)))
}

save_publication_figure <- function(plot, stem, width_mm = 178, height_mm = 125,
                                    dpi = 400, formats = requested_figure_formats()) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  if ("pdf" %in% formats) {
    ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width_mm, height = height_mm,
                    units = "mm", device = grDevices::cairo_pdf, bg = "white")
  }
  if ("png" %in% formats) {
    ggplot2::ggsave(paste0(stem, ".png"), plot, width = width_mm, height = height_mm,
                    units = "mm", dpi = dpi, bg = "white", type = "cairo")
  }
  if ("svg" %in% formats) {
    ggplot2::ggsave(paste0(stem, ".svg"), plot, width = width_mm, height = height_mm,
                    units = "mm", device = svglite::svglite, bg = "white")
  }
  invisible(stem)
}

save_publication_grob <- function(grob, stem, width_mm = 178, height_mm = 140, dpi = 400,
                                  formats = requested_figure_formats()) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  if ("pdf" %in% formats) {
    grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width, height = height)
    grid::grid.newpage(); grid::grid.draw(grob); grDevices::dev.off()
  }
  if ("png" %in% formats) {
    grDevices::png(paste0(stem, ".png"), width = width, height = height, units = "in",
                   res = dpi, type = "cairo", bg = "white")
    grid::grid.newpage(); grid::grid.draw(grob); grDevices::dev.off()
  }
  if ("svg" %in% formats) {
    svglite::svglite(paste0(stem, ".svg"), width = width, height = height)
    grid::grid.newpage(); grid::grid.draw(grob); grDevices::dev.off()
  }
  invisible(stem)
}

plot_pca_publication <- function(data, pc1_percent, pc2_percent, title,
                                 subtitle = "All approved samples are shown; polygons are observed ranges, not confidence intervals") {
  groups <- unique(as.character(data$group))
  palette <- publication_palette(groups)
  hull <- do.call(rbind, lapply(groups, function(level) {
    x <- data[data$group == level, , drop = FALSE]
    if (nrow(x) < 3L) return(NULL)
    x[chull(x$PC1, x$PC2), , drop = FALSE]
  }))
  centers <- stats::aggregate(cbind(PC1, PC2) ~ group, data, mean)
  p <- ggplot2::ggplot(data, ggplot2::aes(PC1, PC2, colour = group, shape = group))
  if (!is.null(hull) && nrow(hull)) {
    p <- p + ggplot2::geom_polygon(data = hull, ggplot2::aes(PC1, PC2, fill = group, group = group),
                                  colour = NA, alpha = 0.09, inherit.aes = FALSE) +
      ggplot2::geom_polygon(data = hull, ggplot2::aes(PC1, PC2, colour = group, group = group),
                           fill = NA, linewidth = 0.55, linetype = 2, inherit.aes = FALSE)
  }
  p +
    ggplot2::geom_hline(yintercept = 0, colour = "#D7DADE", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, colour = "#D7DADE", linewidth = 0.35) +
    ggplot2::geom_point(size = 2.8, stroke = 0.8) +
    ggplot2::geom_point(data = centers, ggplot2::aes(PC1, PC2, colour = group),
                        shape = 4, size = 4, stroke = 1.1, inherit.aes = FALSE) +
    ggrepel::geom_text_repel(ggplot2::aes(label = canonical_sample), size = 2.7,
                             max.overlaps = Inf, show.legend = FALSE, seed = 104729) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = sprintf("PC1 (%.1f%%)", pc1_percent),
                  y = sprintf("PC2 (%.1f%%)", pc2_percent), colour = NULL, shape = NULL) +
    theme_publication()
}

plot_volcano_publication <- function(data, contrast_label, padj_cutoff, lfc_cutoff,
                                     label_count = 15L) {
  x <- data
  x$minus_log10_padj <- -log10(pmax(x$padj, .Machine$double.xmin))
  labels <- x[x$direction != "Not_significant" & !is.na(x$gene_label), , drop = FALSE]
  labels <- head(labels[order(labels$padj, -abs(labels$log2FoldChange_ashr)), , drop = FALSE], label_count)
  plot <- ggplot2::ggplot(x, ggplot2::aes(log2FoldChange_ashr, minus_log10_padj, colour = direction)) +
    ggplot2::geom_point(alpha = 0.65, size = 1) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2,
                        linewidth = 0.45, colour = "black") +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = 2,
                        linewidth = 0.45, colour = "black") +
    ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = gene_label),
                             size = 2.7, max.overlaps = Inf, seed = 104729,
                             show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = direction_palette, drop = FALSE) +
    ggplot2::labs(title = contrast_label,
                  x = "ashr-shrunken log2 fold change",
                  y = "-log10 adjusted P", colour = NULL) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title.position = "plot", legend.position = "right")
  attr(plot, "plot_data") <- x
  plot
}

plot_expression_boxplot <- function(data, title, y_label = "VST expression") {
  groups <- unique(as.character(data$group))
  palette <- publication_palette(groups)
  ggplot2::ggplot(data, ggplot2::aes(group, expression, fill = group, colour = group)) +
    ggplot2::geom_violin(width = 0.82, trim = FALSE, alpha = 0.16, linewidth = 0.45) +
    ggplot2::geom_boxplot(width = 0.36, outlier.shape = NA, alpha = 0.55, linewidth = 0.45) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.07, height = 0, seed = 104729),
                        shape = 21, fill = "white", size = 2.2, stroke = 0.65) +
    ggplot2::scale_fill_manual(values = palette) + ggplot2::scale_colour_manual(values = palette) +
    ggplot2::labs(title = title, x = NULL, y = y_label) + theme_publication() +
    ggplot2::theme(legend.position = "none")
}

clean_ora_term <- function(value, width = 44L) {
  value <- sub("^(GOBP|GOCC|GOMF|KEGG|REACTOME|HALLMARK)_", "", value)
  value <- gsub("_", " ", value, fixed = TRUE)
  vapply(value, function(item) paste(strwrap(item, width = width), collapse = "\n"), character(1))
}

reduce_ora_gene_overlap <- function(data, cutoff = getOption("bulk_rnaseq.ora_jaccard", 0.7)) {
  if (!nrow(data)) return(data)
  groups <- split(data, interaction(data$database, data$direction, drop = TRUE))
  reduced <- lapply(groups, function(x) {
    x <- x[order(x$padj, -x$fold_enrichment), , drop = FALSE]
    genes <- strsplit(as.character(x$overlap_genes), "/", fixed = TRUE)
    keep <- logical(nrow(x)); accepted <- list()
    for (i in seq_len(nrow(x))) {
      candidate <- unique(genes[[i]])
      similarity <- if (!length(accepted)) 0 else max(vapply(accepted, function(previous) {
        union_size <- length(union(candidate, previous))
        if (!union_size) 0 else length(intersect(candidate, previous)) / union_size
      }, numeric(1)))
      if (similarity < cutoff) {
        keep[[i]] <- TRUE
        accepted[[length(accepted) + 1L]] <- candidate
      }
    }
    x[keep, , drop = FALSE]
  })
  out <- do.call(rbind, reduced)
  rownames(out) <- NULL
  out
}

ora_plot_window <- function(data, database, direction, top_n = 15L, reduce_overlap = FALSE,
                            fdr_cutoff = getOption("bulk_rnaseq.ora_fdr", 0.05),
                            jaccard_cutoff = getOption("bulk_rnaseq.ora_jaccard", 0.7)) {
  x <- data[data$database == database & data$direction == direction &
              !is.na(data$padj) & data$padj < fdr_cutoff, , drop = FALSE]
  if (!nrow(x)) return(x)
  if (reduce_overlap) x <- reduce_ora_gene_overlap(x, cutoff = jaccard_cutoff)
  x <- head(x[order(x$padj, -x$fold_enrichment), , drop = FALSE], top_n)
  x$display_term <- clean_ora_term(x$pathway)
  x$term_label <- factor(x$display_term, levels = rev(unique(x$display_term)))
  x$minus_log10_padj <- -log10(pmax(x$padj, .Machine$double.xmin))
  x
}

plot_ora_publication <- function(data, title, top_n = 12L,
                                 fdr_cutoff = getOption("bulk_rnaseq.ora_fdr", 0.05)) {
  x <- data[data$direction %in% c("Up", "Down") & !is.na(data$padj) & data$padj < fdr_cutoff, , drop = FALSE]
  if (!nrow(x)) return(NULL)
  counts <- stats::aggregate(pathway ~ database + direction, x, length)
  names(counts)[[3L]] <- "significant_terms"
  ggplot2::ggplot(counts, ggplot2::aes(database, significant_terms, fill = direction)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.64) +
    ggplot2::geom_text(ggplot2::aes(label = significant_terms),
                       position = ggplot2::position_dodge(width = 0.72), vjust = -0.3, size = 2.7) +
    ggplot2::scale_fill_manual(values = direction_palette[c("Down", "Up")]) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(title = title,
                  subtitle = "Number of FDR-significant terms; detailed terms are shown in database-specific figures",
                  x = NULL, y = sprintf("Significant ORA terms (FDR < %.3g)", fdr_cutoff), fill = NULL,
                  caption = "No database competes for a shared top-N display quota; complete tested terms remain in the tables.") +
    theme_publication() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 28, hjust = 1))
}

plot_ora_dotplot <- function(data, database, direction, title, top_n = 15L,
                             fdr_cutoff = getOption("bulk_rnaseq.ora_fdr", 0.05)) {
  x <- ora_plot_window(data, database, direction, top_n, reduce_overlap = FALSE, fdr_cutoff = fdr_cutoff)
  if (!nrow(x)) return(NULL)
  ggplot2::ggplot(x, ggplot2::aes(fold_enrichment, term_label, size = overlap_count,
                                  colour = minus_log10_padj)) +
    ggplot2::geom_point(alpha = 0.92) +
    ggplot2::scale_colour_gradientn(colours = c("#3B4CC0", "#1FA187", "#FDE725")) +
    ggplot2::scale_size_continuous(range = c(2.2, 6.5)) +
    ggplot2::labs(title = title,
                  subtitle = sprintf("Top %d FDR-significant terms ordered by adjusted P", nrow(x)),
                  x = "Fold enrichment", y = NULL, size = "Overlap", colour = expression(-log[10](FDR)),
                  caption = "Dots show fold enrichment; every tested term remains in the tables.") +
    theme_publication() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(colour = "#E6E8EB", linewidth = 0.35))
}

plot_ora_lollipop <- function(data, database, direction, title, top_n = 15L,
                              fdr_cutoff = getOption("bulk_rnaseq.ora_fdr", 0.05),
                              jaccard_cutoff = getOption("bulk_rnaseq.ora_jaccard", 0.7)) {
  x <- ora_plot_window(data, database, direction, top_n, reduce_overlap = TRUE,
                       fdr_cutoff = fdr_cutoff, jaccard_cutoff = jaccard_cutoff)
  if (!nrow(x)) return(NULL)
  ggplot2::ggplot(x, ggplot2::aes(minus_log10_padj, term_label)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = minus_log10_padj, yend = term_label),
                          colour = "#CDD2D8", linewidth = 0.55) +
    ggplot2::geom_point(ggplot2::aes(size = overlap_count, colour = fold_enrichment), alpha = 0.95) +
    ggplot2::scale_colour_gradientn(colours = c("#2C3E50", "#2A9D8F", "#F4D35E")) +
    ggplot2::scale_size_continuous(range = c(2.2, 6.5)) +
    ggplot2::labs(title = title,
                  subtitle = sprintf("Gene-overlap-reduced representative terms (Jaccard < %.2f)", jaccard_cutoff),
                  x = expression(-log[10](adjusted~P)), y = NULL, size = "Overlap", colour = "Fold enrichment",
                  caption = "This reduced display does not remove rows from the complete ORA table.") +
    theme_publication() + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

plot_ora_diverging <- function(data, database, title, top_n = 8L,
                               fdr_cutoff = getOption("bulk_rnaseq.ora_fdr", 0.05),
                               jaccard_cutoff = getOption("bulk_rnaseq.ora_jaccard", 0.7)) {
  x <- data[data$database == database & data$direction %in% c("Up", "Down") &
              !is.na(data$padj) & data$padj < fdr_cutoff, , drop = FALSE]
  if (!nrow(x)) return(NULL)
  x <- reduce_ora_gene_overlap(x, cutoff = jaccard_cutoff)
  x <- do.call(rbind, lapply(split(x, x$direction), function(z) head(z[order(z$padj), , drop = FALSE], top_n)))
  x$signed_score <- ifelse(x$direction == "Up", 1, -1) * -log10(pmax(x$padj, .Machine$double.xmin))
  x$display_term <- paste0(ifelse(x$direction == "Up", "UP: ", "DOWN: "), clean_ora_term(x$pathway, 38L))
  x$term_label <- factor(x$display_term, levels = x$display_term[order(x$signed_score)])
  ggplot2::ggplot(x, ggplot2::aes(signed_score, term_label, fill = direction)) +
    ggplot2::geom_col(width = 0.68) +
    ggplot2::geom_vline(xintercept = 0, colour = "#596068", linewidth = 0.45) +
    ggplot2::scale_fill_manual(values = direction_palette[c("Down", "Up")]) +
    ggplot2::scale_x_continuous(labels = function(value) abs(value)) +
    ggplot2::labs(title = title,
                  subtitle = "Representative terms after within-direction gene-overlap reduction",
                  x = expression(-log[10](adjusted~P)), y = NULL,
                  caption = "Left = downregulated DEG subset; right = upregulated DEG subset.") +
    theme_publication() + ggplot2::theme(legend.position = "none", panel.grid.major.y = ggplot2::element_blank())
}

plot_ora_sankey_bubble <- function(data, title, top_n = 2L,
                                   fdr_cutoff = getOption("bulk_rnaseq.ora_fdr", 0.05)) {
  x <- data[data$direction %in% c("Up", "Down") & !is.na(data$padj) & data$padj < fdr_cutoff, , drop = FALSE]
  if (!nrow(x)) return(NULL)
  x <- do.call(rbind, lapply(split(x, interaction(x$database, x$direction, drop = TRUE)),
                            function(z) head(z[order(z$padj, -z$fold_enrichment), , drop = FALSE], top_n)))
  database_count <- length(unique(x$database))
  x$term_label <- gsub("^(GOBP|GOCC|GOMF|KEGG|REACTOME|HALLMARK)_", "", x$pathway)
  x$term_label <- gsub("_", " ", x$term_label, fixed = TRUE)
  label_width <- if (database_count > 1L) 23L else 30L
  x$term_label <- vapply(x$term_label, function(value) paste(strwrap(value, label_width), collapse = "\n"), character(1))
  ggplot2::ggplot(x, ggplot2::aes(axis1 = direction, axis2 = term_label,
                                  y = overlap_count, fill = direction)) +
    ggalluvial::geom_alluvium(width = 0.10, alpha = 0.58, knot.pos = 0.35) +
    ggalluvial::geom_stratum(width = 0.11, colour = "white", linewidth = 0.35) +
    ggplot2::geom_text(stat = ggalluvial::StatStratum,
                       ggplot2::aes(label = ggplot2::after_stat(stratum)),
                       size = if (database_count > 1L) 1.55 else 2.15, lineheight = 0.86) +
    ggplot2::scale_x_discrete(limits = c("Direction", "Enriched term"), expand = c(0.10, 0.10)) +
    ggplot2::scale_fill_manual(values = direction_palette[c("Down", "Up")]) +
    ggplot2::facet_wrap(~database, ncol = if (database_count > 1L) 2 else 1, scales = "free_y") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = title,
                  subtitle = "Flow width is the number of overlapping genes; terms and directions are unchanged",
                  x = NULL, y = "Overlapping genes", fill = NULL,
                  caption = "Only the display is reduced to the top terms; complete ORA tables retain every tested term.") +
    theme_publication() + ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_blank(),
                                        axis.ticks.y = ggplot2::element_blank(), panel.spacing = grid::unit(10, "pt"),
                                        plot.margin = ggplot2::margin(5.5, 35, 5.5, 5.5))
}

plot_gsea_overview <- function(data, title, top_n = 12L,
                               fdr_cutoff = getOption("bulk_rnaseq.gsea_reporting_fdr", 1)) {
  if (!nrow(data)) return(NULL)
  significant <- data[!is.na(data$padj) & data$padj <= fdr_cutoff, , drop = FALSE]
  if (!nrow(significant)) return(NULL)
  significant <- significant[order(significant$padj, -abs(significant$NES)), , drop = FALSE]
  x <- head(significant, top_n)
  x$direction <- ifelse(x$NES >= 0, "Up", "Down")
  x$pathway_label <- factor(x$pathway, levels = rev(x$pathway))
  ggplot2::ggplot(x, ggplot2::aes(NES, pathway_label, size = size, colour = direction)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#8C9298", linewidth = 0.45) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = pathway_label),
                          colour = "#D9DDE1", linewidth = 0.45) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_colour_manual(values = direction_palette[c("Down", "Up")]) +
    ggplot2::scale_size_continuous(range = c(2, 6)) +
    ggplot2::labs(title = title,
                  subtitle = sprintf("Top %d pathways ordered by FDR and |NES|; sign is preserved", nrow(x)),
                  x = "Normalized enrichment score (NES)", y = NULL,
                  size = "Set size", colour = "Rank direction",
                  caption = "NES sign follows the declared contrast; GSEA is not causal evidence.") +
    theme_publication() + ggplot2::theme(panel.grid.major.x = ggplot2::element_line(colour = "#E6E8EB", linewidth = 0.35))
}

plot_gsea_curve <- function(pathway_genes, ranks, summary_row, contrast_label) {
  curve <- fgsea::plotEnrichment(pathway_genes, ranks)
  running <- as.data.frame(curve$data)
  hits <- as.data.frame(curve$layers[[2L]]$data)
  names(hits) <- c("rank", "rank_metric")
  peak_index <- if (summary_row$ES[[1]] >= 0) which.max(running$ES) else which.min(running$ES)
  peak_rank <- running$rank[[peak_index]]
  peak_es <- running$ES[[peak_index]]
  color <- if (summary_row$NES[[1]] >= 0) direction_palette[["Up"]] else direction_palette[["Down"]]
  p_curve <- ggplot2::ggplot(running, ggplot2::aes(rank, ES)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#777D84", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = peak_rank, colour = color, linetype = 3, linewidth = 0.45) +
    ggplot2::geom_line(colour = color, linewidth = 0.8) +
    ggplot2::geom_point(data = running[peak_index, , drop = FALSE], shape = 21, fill = "white",
                        colour = color, size = 2.4, stroke = 0.8) +
    ggplot2::annotate("label", x = max(running$rank) * 0.98,
                      y = if (peak_es >= 0) max(running$ES) * 0.82 else min(running$ES) * 0.82,
                      label = sprintf("NES = %.2f\nFDR = %.2g\nSet size = %d",
                                      summary_row$NES, summary_row$padj, summary_row$size),
                      hjust = 1, size = 2.5, linewidth = 0.25) +
    ggplot2::labs(y = "Running enrichment score", x = NULL) + theme_publication() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  p_hits <- ggplot2::ggplot(hits, ggplot2::aes(rank, 0)) +
    ggplot2::geom_segment(ggplot2::aes(xend = rank, y = -0.7, yend = 0.7), linewidth = 0.28, colour = "#343A40") +
    ggplot2::geom_vline(xintercept = peak_rank, colour = color, linetype = 3, linewidth = 0.45) +
    ggplot2::coord_cartesian(ylim = c(-1, 1)) + ggplot2::labs(x = NULL, y = NULL) + theme_publication() +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  metric <- data.frame(rank = seq_along(ranks), metric = unname(ranks))
  p_metric <- ggplot2::ggplot(metric, ggplot2::aes(rank, metric)) +
    ggplot2::geom_area(data = metric[metric$metric >= 0, ], fill = direction_palette[["Up"]], alpha = 0.8) +
    ggplot2::geom_area(data = metric[metric$metric < 0, ], fill = direction_palette[["Down"]], alpha = 0.8) +
    ggplot2::geom_hline(yintercept = 0, colour = "#777D84", linewidth = 0.35) +
    ggplot2::labs(x = "Rank in complete DESeq2 Wald-statistic list", y = "Rank metric") + theme_publication()
  combined <- patchwork::wrap_plots(p_curve, p_hits, p_metric, ncol = 1, heights = c(3.5, 0.65, 1.4)) +
    patchwork::plot_annotation(
      title = summary_row$pathway[[1]],
      subtitle = paste0(contrast_label, " · running score, hit positions and complete ranked metric"),
      caption = "NES and FDR are taken from the saved GSEA result; enrichment does not establish pathway activation or causality.",
      theme = theme_publication()
    )
  list(plot = combined, running = running, hits = hits, metric = metric)
}

threshold_profiles <- function(config) {
  profiles <- config$thresholds$profiles
  if (is.null(profiles)) {
    profiles <- list(Primary = list(padj = config$thresholds$padj, abs_lfc = config$thresholds$abs_lfc))
  }
  for (name in names(profiles)) {
    if (is.null(profiles[[name]]$padj) || is.null(profiles[[name]]$abs_lfc)) {
      stop("Threshold profile requires padj and abs_lfc: ", name)
    }
  }
  profiles
}

comparison_module_dir <- function(context, contrast_id, module_name) {
  if (nrow(context$contrasts) > 1L) {
    file.path(context$run_dir, "Comparisons", contrast_id, module_name)
  } else {
    file.path(context$run_dir, module_name, contrast_id)
  }
}

resource_registry <- function(context) {
  registry_path <- file.path(context$root, "resources", "registry.yml")
  registry <- yaml::read_yaml(registry_path)
  resources <- registry$resources %||% list()
  names(resources) <- vapply(resources, `[[`, character(1), "resource_id")
  list(path = registry_path, entries = resources)
}

resource_path <- function(context, resource_id) {
  registry <- resource_registry(context)
  entry <- registry$entries[[resource_id]]
  if (is.null(entry)) stop("Frozen resource is not registered: ", resource_id)
  path <- normalizePath(file.path(dirname(registry$path), entry$local_path), winslash = "/", mustWork = TRUE)
  if (!identical(tolower(sha256_file(path)), tolower(entry$sha256))) stop("Frozen resource SHA256 mismatch: ", resource_id)
  if (!entry$species %in% c(context$config$species, "not_applicable")) {
    stop("Resource species does not match project: ", resource_id)
  }
  path
}

record_module_status <- function(context, module, status, detail = "") {
  allowed <- c("complete", "not_applicable", "skipped_by_user", "failed_explicit")
  if (!status %in% allowed) stop("Invalid module status: ", status)
  context$module_status[[module]] <- list(status = status, detail = detail)
  write_result_table(data.frame(module = module, status = status, detail = detail),
                     file.path(context$run_dir, "Provenance", paste0(module, "_status.tsv")))
  context
}

apply_factor_levels <- function(samples, config) {
  levels_config <- config$design$factor_levels %||% list()
  for (name in names(levels_config)) {
    if (!name %in% names(samples)) stop("factor_levels names a missing samples.tsv column: ", name)
    declared <- as.character(levels_config[[name]])
    observed <- unique(as.character(samples[[name]]))
    if (!setequal(declared, observed)) {
      stop("Declared factor levels do not equal observed levels for ", name,
           ": declared=", paste(declared, collapse = ","), "; observed=", paste(observed, collapse = ","))
    }
    samples[[name]] <- factor(samples[[name]], levels = declared)
  }
  samples
}

initialize_context <- function(root, project_dir, run_dir) {
  root <- workflow_root(root)
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = FALSE)
  expected_parent <- normalizePath(file.path(project_dir, "work", "staging"), winslash = "/", mustWork = FALSE)
  if (!startsWith(run_dir, paste0(expected_parent, "/"))) {
    stop("Formal analyses must be staged under project/work/staging/<run_id>")
  }
  motif_phase <- identical(Sys.getenv("BULK_RNASEQ_MOTIF_PHASE"), "1")
  if (!motif_phase && dir.exists(run_dir)) {
    existing <- list.files(run_dir, all.files = TRUE, no.. = TRUE)
    # Snakemake creates the declared log directory before the initialize rule.
    # An otherwise-empty Logs directory is not a prior analysis and is safe.
    log_entries <- if (identical(existing, "Logs")) {
      list.files(file.path(run_dir, "Logs"), all.files = TRUE, no.. = TRUE, full.names = TRUE)
    } else character()
    provenance_entries <- if ("Provenance" %in% existing) {
      list.files(file.path(run_dir, "Provenance"), all.files = TRUE, no.. = TRUE, full.names = TRUE)
    } else character()
    allowed_preinit_provenance <- c("upstream_contract_checked")
    permitted <- identical(Sys.getenv("BULK_RNASEQ_SNAKEMAKE"), "1") &&
      all(existing %in% c("Logs", "Provenance")) &&
      all(basename(provenance_entries) %in% allowed_preinit_provenance) &&
      (!length(log_entries) || setequal(basename(log_entries), "00_initialize.log"))
    if (length(existing) && !permitted) {
      stop("Staging run directory already exists and is not empty: ", run_dir,
           " [snakemake=", Sys.getenv("BULK_RNASEQ_SNAKEMAKE"),
           "; entries=", paste(existing, collapse = ","),
           "; logs=", paste(basename(log_entries), collapse = ","), "]")
    }
  }
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  config_path <- file.path(project_dir, "project.yml")
  config <- yaml::read_yaml(config_path)
  qc_preview_mode <- identical(Sys.getenv("BULK_RNASEQ_QC_PREVIEW"), "1")
  if (qc_preview_mode) {
    qc_decision <- config$analysis$modules$qc
    if (!isTRUE(qc_decision$confirmed) || !identical(qc_decision$status, "enabled")) {
      stop("QC preview requires an explicitly confirmed enabled QC module")
    }
  } else {
    expected_confirmation <- module_confirmation_string(config)
    actual_confirmation <- Sys.getenv("BULK_RNASEQ_MODULE_CONFIRMATION", unset = "")
    if (!identical(actual_confirmation, expected_confirmation)) {
      stop("Every analysis requires a fresh exact module confirmation. Expected: ", expected_confirmation)
    }
  }
  figure_formats <- config$export$formats %||% "pdf"
  if (!"pdf" %in% figure_formats) stop("PDF is the mandatory default figure format")
  options(bulk_rnaseq.figure_formats = figure_formats)
  inputs <- config$inputs %||% list()
  paths <- list(
    counts = file.path(project_dir, inputs$counts_file %||% "input/counts.tsv.gz"),
    samples = file.path(project_dir, inputs$samples_file %||% "samples.tsv"),
    contrasts = file.path(project_dir, inputs$contrasts_file %||% "contrasts.tsv"),
    approval = file.path(project_dir, inputs$qc_approval_file %||% "qc_approval.yml"),
    source_manifest = file.path(project_dir, "input", "source_manifest.yml")
  )
  missing <- names(paths)[!vapply(paths, file.exists, logical(1))]
  if (length(missing)) stop("Required project inputs are missing: ", paste(missing, collapse = ", "))
  expected_plan <- plan_confirmation_string(config, paths, config_path)
  actual_plan <- Sys.getenv("BULK_RNASEQ_PLAN_CONFIRMATION", unset = "")
  if (!identical(actual_plan, expected_plan)) {
    stop("Every analysis requires fresh approval of the complete scientific plan. Expected: ", expected_plan)
  }
  samples <- read_project_table(paths$samples)
  contrasts <- read_project_table(paths$contrasts)
  required_samples <- c("sample", "excluded", "exclusion_reason")
  if (length(setdiff(required_samples, names(samples)))) stop("samples.tsv lacks required columns")
  samples$canonical_sample <- samples$canonical_sample %||% samples$sample
  samples$excluded <- as_flag(samples$excluded, "samples.tsv excluded")
  samples <- apply_factor_levels(samples, config)
  required_contrasts <- c("contrast_id", "factor", "numerator", "denominator", "confirmed")
  if (length(setdiff(required_contrasts, names(contrasts)))) stop("contrasts.tsv lacks required columns")
  safe_id(contrasts$contrast_id, "contrast_id")
  if (anyDuplicated(contrasts$contrast_id)) stop("contrast_id values must be unique")
  if (!qc_preview_mode && any(!as_flag(contrasts$confirmed, "contrasts.tsv confirmed"))) {
    stop("Every contrast must be explicitly confirmed")
  }
  approval <- yaml::read_yaml(paths$approval)
  source_manifest <- yaml::read_yaml(paths$source_manifest)
  count_hash <- sha256_file(paths$counts)
  if (!identical(tolower(count_hash), tolower(source_manifest$sha256))) stop("Count matrix SHA256 differs from source_manifest.yml")
  if (!qc_preview_mode && (!isTRUE(approval$approved) ||
                           !identical(tolower(count_hash), tolower(approval$counts_sha256)))) {
    stop("QC approval is missing or stale for the current count matrix")
  }
  if (!identical(config$counts_provenance$normalization, "raw_counts")) stop("DESeq2 requires confirmed raw_counts")
  env <- new.env(parent = emptyenv())
  env$root <- root; env$project_dir <- project_dir; env$run_dir <- run_dir
  env$config <- config; env$paths <- paths; env$samples <- samples; env$contrasts <- contrasts
  env$approval <- approval; env$source_manifest <- source_manifest; env$count_hash <- count_hash
  env$module_status <- list(); env$state <- list(); env$started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  dir.create(file.path(run_dir, "Provenance"), recursive = TRUE, showWarnings = FALSE)
  env
}
