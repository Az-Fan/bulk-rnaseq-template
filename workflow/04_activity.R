# Module 04 — pathway and regulator activity
#
# Sample-level scores reuse the fitted model matrix and contrast directions.
# Weighted-regulator outputs are descriptive estimates; no uncalibrated normal
# approximation is reported as a confirmatory P value.

run_activity <- function(context) {
  state <- canonical_module_state(context$config, c("regulation"))
  if (!identical(state$status, "enabled")) {
    return(record_module_status(context, "04_Regulation", state$status, state$reason))
  }
  if (is.null(context$state$vst)) stop("Differential expression must complete before activity analysis")
  root <- file.path(context$run_dir, "04_Regulation")
  tables <- file.path(root, "Tables"); figures <- file.path(root, "Figures")
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures, recursive = TRUE, showWarnings = FALSE)

  resource_ids <- context$config$resources$activity %||% list()
  if (!length(resource_ids)) stop("An enabled regulation module requires resources.activity declarations")
  expression <- context$state$vst
  annotation <- context$state$annotation
  symbols <- annotation[rownames(expression), "gene_symbol"]
  valid <- !is.na(symbols) & nzchar(symbols)
  expression <- expression[valid, , drop = FALSE]
  rownames(expression) <- symbols[valid]
  expression <- expression[!duplicated(rownames(expression)), , drop = FALSE]
  metadata <- context$state$samples[context$state$selected_samples, , drop = FALSE]
  rownames(metadata) <- metadata$sample
  score_tables <- list()

  score_weighted_network <- function(network, weight_column, label) {
    required <- c("source", "target", weight_column)
    if (length(setdiff(required, names(network)))) stop(label, " resource lacks required columns")
    network <- network[network$target %in% rownames(expression), , drop = FALSE]
    rows <- lapply(split(network, network$source), function(x) {
      weights <- as.numeric(x[[weight_column]])
      denominator <- sqrt(sum(weights^2))
      if (!is.finite(denominator) || denominator == 0 || nrow(x) < 5L) return(NULL)
      values <- expression[x$target, , drop = FALSE]
      data.frame(source = x$source[[1]], sample = colnames(values),
                 score = as.numeric(crossprod(weights, values)) / denominator,
                 targets_used = nrow(x), stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }

  for (resource_id in unlist(resource_ids, use.names = FALSE)) {
    network <- readRDS(resource_path(context, resource_id))
    if (all(c("source", "target", "mor") %in% names(network))) {
      score <- score_weighted_network(network, "mor", "TF activity")
      score$method <- "CollecTRI weighted expression score"
      score$type <- "TF"
    } else if (all(c("source", "target", "weight") %in% names(network))) {
      score <- score_weighted_network(network, "weight", "PROGENy")
      score$method <- "PROGENy weighted expression score"
      score$type <- "Pathway"
    } else {
      stop("Unknown activity resource format: ", resource_id)
    }
    score$resource_id <- resource_id
    score_tables[[resource_id]] <- score
  }
  scores <- do.call(rbind, score_tables)
  write_result_table(scores, file.path(tables, "Activity_Sample_Scores.tsv.gz"))
  contrasts <- list()
  for (i in seq_len(nrow(context$contrasts))) {
    contrast <- context$contrasts[i, ]
    factor_values <- metadata[scores$sample, contrast$factor]
    selected <- factor_values %in% c(contrast$numerator, contrast$denominator)
    subset <- scores[selected, , drop = FALSE]
    subset$group <- factor_values[selected]
    estimates <- do.call(rbind, lapply(split(subset, interaction(subset$type, subset$source, drop = TRUE)), function(x) {
      numerator <- x$score[x$group == contrast$numerator]
      denominator <- x$score[x$group == contrast$denominator]
      data.frame(contrast_id = contrast$contrast_id, type = x$type[[1]], source = x$source[[1]],
                 estimate = mean(numerator) - mean(denominator),
                 numerator_mean = mean(numerator), denominator_mean = mean(denominator),
                 targets_used = x$targets_used[[1]], method = x$method[[1]], stringsAsFactors = FALSE)
    }))
    estimates$direction <- ifelse(estimates$estimate >= 0, "Up", "Down")
    contrasts[[contrast$contrast_id]] <- estimates
    write_result_table(estimates, file.path(tables, paste0("Activity_", contrast$contrast_id, "_Descriptive_Contrasts.tsv")))
    top <- head(estimates[order(-abs(estimates$estimate)), ], 25L)
    top$source <- factor(top$source, levels = rev(top$source))
    plot <- ggplot2::ggplot(top, ggplot2::aes(estimate, source, fill = direction)) +
      ggplot2::geom_col(width = 0.72) + ggplot2::geom_vline(xintercept = 0, colour = "#343A40", linewidth = 0.4) +
      ggplot2::facet_grid(type ~ ., scales = "free_y", space = "free_y") +
      ggplot2::scale_fill_manual(values = direction_palette[c("Down", "Up")]) +
      ggplot2::labs(title = paste0(contrast$contrast_id, " — expression-derived activity"),
                    subtitle = "Descriptive group-mean differences; no uncalibrated normal-approximation P values",
                    x = paste0("Mean score difference (", contrast$numerator, " − ", contrast$denominator, ")"),
                    y = NULL, fill = NULL,
                    caption = "Activity scores are expression-derived associations, not direct binding or causal evidence.") +
      theme_publication()
    save_publication_figure(plot, file.path(figures, paste0("Activity_", contrast$contrast_id, "_Overview")),
                            178, max(115, 5.2 * nrow(top) + 45))
  }
  context$state$activity <- list(scores = scores, contrasts = contrasts)
  record_module_status(context, "04_Regulation", "complete",
                       sprintf("%d frozen activity resources; descriptive contrasts for %d comparisons",
                               length(resource_ids), nrow(context$contrasts)))
}
