# Module 06 — promoter and optional peak-aware motif analysis
#
# Counts-only promoter motif is exploratory. Up and Down genes are analysed
# separately. Controls are drawn from tested non-DE genes after exact promoter-
# length matching and stratified GC/base-expression matching. Motif enrichment
# is not direct binding evidence.

run_motif <- function(context) {
  promoter <- module_decision(context$config, "motif_promoter")
  peaks <- module_decision(context$config, "motif_peaks")
  if (!identical(promoter$status, "enabled") && !identical(peaks$status, "enabled")) {
    reason <- paste(unique(c(promoter$reason %||% "", peaks$reason %||% "")), collapse = "; ")
    return(record_module_status(context, "06_Motif",
                                if (all(c(promoter$status, peaks$status) == "not_applicable")) "not_applicable" else "skipped_by_user",
                                reason))
  }
  if (identical(peaks$status, "enabled")) {
    stop("Peak-aware motif is enabled but no generic peak executor has been implemented")
  }
  if (!identical(Sys.getenv("BULK_RNASEQ_MOTIF_PHASE"), "1")) {
    write_result_table(data.frame(module = "06_Motif", status = "pending_motif_phase",
                                  detail = "The control plane will invoke the locked motif environment"),
                       file.path(context$run_dir, "Provenance", "06_Motif_status.tsv"))
    context$module_status[["06_Motif"]] <- list(status = "pending_motif_phase",
                                                 detail = "Awaiting locked motif environment")
    return(context)
  }
  executable <- vapply(c("streme", "ame", "bedtools", "samtools"), Sys.which, character(1))
  if (any(!nzchar(executable))) {
    stop("Motif is enabled but the active Pixi environment lacks: ", paste(names(executable)[!nzchar(executable)], collapse = ", "))
  }
  matching <- context$config$motif$background_matching %||% list()
  required <- c("gc", "promoter_length", "base_expression")
  if (!all(vapply(required, function(x) isTRUE(matching[[x]]), logical(1)))) {
    stop("Promoter motif requires confirmed matching on GC, promoter length and base expression")
  }
  fasta <- resource_path(context, context$config$resources$motif_fasta)
  gtf <- resource_path(context, context$config$resources$motif_gtf)
  motif_database <- resource_path(context, context$config$resources$motif_database)
  root <- file.path(context$run_dir, "06_Motif", "Promoter")
  tables <- file.path(root, "Tables"); figures <- file.path(root, "Figures")
  dir.create(tables, recursive = TRUE, showWarnings = FALSE); dir.create(figures, recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("motif-"); dir.create(temporary); on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  genome <- file.path(temporary, "genome.fa")
  if (grepl("\\.gz$", fasta, ignore.case = TRUE)) {
    input <- gzfile(fasta, "rb"); output <- file(genome, "wb")
    repeat { block <- readBin(input, "raw", 1024L * 1024L); if (!length(block)) break; writeBin(block, output) }
    close(input); close(output)
  } else file.copy(fasta, genome)
  status <- system2(executable[["samtools"]], c("faidx", genome))
  if (status != 0L) stop("samtools faidx failed")

  command <- if (grepl("\\.gz$", gtf, ignore.case = TRUE)) {
    sprintf("gzip -dc %s | awk -F '\\t' '$3 == \"gene\"'", shQuote(gtf))
  } else {
    sprintf("awk -F '\\t' '$3 == \"gene\"' %s", shQuote(gtf))
  }
  genes <- data.table::fread(cmd = command, sep = "\t", header = FALSE, quote = "", fill = TRUE,
                             select = c(1L, 4L, 5L, 7L, 9L), data.table = FALSE)
  promoters <- do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
    x <- genes[i, ]
    match <- regexec('gene_id "([^"]+)"', x[[5L]])
    id <- regmatches(x[[5L]], match)[[1L]]
    if (length(id) < 2L) return(NULL)
    id <- sub("\\.[0-9]+$", "", id[[2L]])
    start <- as.integer(x[[2L]]); end <- as.integer(x[[3L]])
    if (x[[4L]] == "+") c(x[[1L]], max(0L, start - 1001L), start + 100L, id, ".", x[[4L]])
    else c(x[[1L]], max(0L, end - 100L), end + 1001L, id, ".", x[[4L]])
  }))
  promoters <- as.data.frame(promoters, stringsAsFactors = FALSE)
  names(promoters) <- c("chrom", "start", "end", "gene_id", "score", "strand")
  promoters$start <- as.integer(promoters$start); promoters$end <- as.integer(promoters$end)
  promoters <- promoters[!duplicated(promoters$gene_id), ]
  write_result_table(promoters, file.path(tables, "All_Annotated_Promoters.tsv.gz"))

  write_fasta <- function(bed, stem) {
    bed_path <- file.path(temporary, paste0(stem, ".bed"))
    fasta_path <- file.path(temporary, paste0(stem, ".fa"))
    data.table::fwrite(bed[, c("chrom", "start", "end", "gene_id", "score", "strand")], bed_path,
                       sep = "\t", col.names = FALSE)
    status <- system2(executable[["bedtools"]], c("getfasta", "-s", "-name", "-fi", genome, "-bed", bed_path, "-fo", fasta_path))
    if (status != 0L) stop("bedtools getfasta failed for ", stem)
    fasta_path
  }
  fasta_gc <- function(path) {
    lines <- readLines(path, warn = FALSE)
    header <- startsWith(lines, ">")
    index <- cumsum(header)
    ids <- sub("^>([^:]+).*", "\\1", lines[header])
    sequence <- vapply(split(lines[!header], index[!header]), paste, collapse = "", character(1))
    gc <- vapply(sequence, function(x) {
      bases <- strsplit(toupper(x), "", fixed = TRUE)[[1L]]
      mean(bases %in% c("G", "C"))
    }, numeric(1))
    data.frame(gene_id = ids, gc_fraction = gc, sequence_length = nchar(sequence), stringsAsFactors = FALSE)
  }

  # Extract all tested promoter sequences once, then compute matching covariates.
  tested <- unique(unlist(lapply(context$state$differential, function(x) x$result$gene_id)))
  tested_promoters <- promoters[promoters$gene_id %in% tested, ]
  all_fasta <- write_fasta(tested_promoters, "tested_promoters")
  covariates <- fasta_gc(all_fasta)
  base_mean <- do.call(rbind, lapply(names(context$state$differential), function(id) {
    x <- context$state$differential[[id]]$result
    data.frame(contrast_id = id, gene_id = x$gene_id, baseMean = x$baseMean, direction = x$direction)
  }))
  summaries <- list(); known_results <- list()
  for (contrast_id in names(context$state$differential)) {
    de <- context$state$differential[[contrast_id]]$result
    background_candidates <- merge(de[de$direction == "Not_significant", c("gene_id", "baseMean")], covariates, by = "gene_id")
    for (direction in c("Up", "Down")) {
      foreground <- merge(de[de$direction == direction, c("gene_id", "baseMean")], covariates, by = "gene_id")
      if (nrow(foreground) < 10L) {
        summaries[[paste(contrast_id, direction)]] <- data.frame(contrast_id, direction, status = "not_applicable_too_few_foreground",
                                                                  foreground = nrow(foreground), background = 0L)
        next
      }
      pool <- background_candidates
      gc_breaks <- unique(stats::quantile(c(foreground$gc_fraction, pool$gc_fraction),
                                          probs = seq(0, 1, length.out = 6), na.rm = TRUE))
      pool$gc_bin <- cut(pool$gc_fraction, breaks = gc_breaks, include.lowest = TRUE)
      foreground$gc_bin <- cut(foreground$gc_fraction, breaks = gc_breaks, include.lowest = TRUE)
      expression_breaks <- unique(stats::quantile(log1p(c(foreground$baseMean, pool$baseMean)), probs = seq(0, 1, length.out = 6), na.rm = TRUE))
      pool$expression_bin <- cut(log1p(pool$baseMean), breaks = expression_breaks, include.lowest = TRUE)
      foreground$expression_bin <- cut(log1p(foreground$baseMean), breaks = expression_breaks, include.lowest = TRUE)
      set.seed(as.integer(context$config$analysis$random_seed) + match(direction, c("Up", "Down")))
      used <- character()
      background <- do.call(rbind, lapply(seq_len(nrow(foreground)), function(i) {
        available <- pool[!pool$gene_id %in% used, ]
        same_stratum <- available$gc_bin == foreground$gc_bin[[i]] &
          available$expression_bin == foreground$expression_bin[[i]]
        candidates <- available[!is.na(same_stratum) & same_stratum, ]
        if (!nrow(candidates)) candidates <- available
        gc_scale <- stats::sd(pool$gc_fraction) %||% 1
        expression_scale <- stats::sd(log1p(pool$baseMean)) %||% 1
        distance <- abs(candidates$gc_fraction - foreground$gc_fraction[[i]]) / max(gc_scale, 1e-6) +
          abs(log1p(candidates$baseMean) - log1p(foreground$baseMean[[i]])) / max(expression_scale, 1e-6) +
          abs(candidates$sequence_length - foreground$sequence_length[[i]])
        chosen <- candidates[which.min(distance), , drop = FALSE]
        used <<- c(used, chosen$gene_id)
        chosen
      }))
      if (nrow(background) < 10L) stop("Matched motif background retained fewer than ten unique genes")
      foreground_bed <- tested_promoters[match(foreground$gene_id, tested_promoters$gene_id), ]
      background_bed <- tested_promoters[match(background$gene_id, tested_promoters$gene_id), ]
      stem <- paste0(safe_id(contrast_id), "_", direction)
      fg_fasta <- write_fasta(foreground_bed, paste0(stem, "_foreground"))
      bg_fasta <- write_fasta(background_bed, paste0(stem, "_background"))
      out <- file.path(root, stem); dir.create(out, recursive = TRUE, showWarnings = FALSE)
      streme_log <- file.path(out, "STREME.log")
      streme <- system2(executable[["streme"]], c("--p", fg_fasta, "--n", bg_fasta, "--oc", file.path(out, "STREME"),
                                                   "--dna", "--minw", "6", "--maxw", "15", "--thresh", "0.05"),
                        stdout = streme_log, stderr = streme_log)
      if (streme != 0L) stop("STREME failed for ", stem, ": ", paste(tail(readLines(streme_log, warn = FALSE), 6L), collapse = " | "))
      ame_path <- file.path(out, "AME_JASPAR.tsv")
      ame_log <- file.path(out, "AME.log")
      ame <- system2(executable[["ame"]], c("--text", "--control", bg_fasta, fg_fasta, motif_database), stdout = ame_path, stderr = ame_log)
      if (ame != 0L) stop("AME failed for ", stem, ": ", paste(tail(readLines(ame_log, warn = FALSE), 6L), collapse = " | "))
      ame_table <- tryCatch(read.delim(ame_path, comment.char = "#", check.names = FALSE), error = function(e) data.frame())
      if (nrow(ame_table)) {
        ame_table$contrast_id <- contrast_id; ame_table$direction <- direction
        known_results[[stem]] <- ame_table
        p_column <- intersect(c("adj_p-value", "p-value"), names(ame_table))[[1L]]
        name_column <- intersect(c("motif_alt_ID", "motif_ID"), names(ame_table))[[1L]]
        top <- head(ame_table[order(ame_table[[p_column]]), ], 15L)
        top$motif <- factor(top[[name_column]], levels = rev(top[[name_column]]))
        top$minus_log10_p <- -log10(pmax(top[[p_column]], .Machine$double.xmin))
        plot <- ggplot2::ggplot(top, ggplot2::aes(minus_log10_p, motif)) +
          ggplot2::geom_segment(ggplot2::aes(x = 0, xend = minus_log10_p, yend = motif), colour = "#D9DDE1", linewidth = 0.5) +
          ggplot2::geom_point(colour = if (direction == "Up") direction_palette[["Up"]] else direction_palette[["Down"]], size = 2.8) +
          ggplot2::labs(title = paste0(contrast_id, " — ", direction, " promoter motifs"),
                        subtitle = "JASPAR known-motif enrichment against GC/expression/length-matched tested promoters",
                        x = expression(-log[10](adjusted~P)), y = NULL,
                        caption = "Counts-only promoter motif is exploratory and does not establish direct TF binding.") + theme_publication()
        save_publication_figure(plot, file.path(figures, paste0(stem, "_Known_Motif_Enrichment")), 145, 115)
      } else {
        status_plot <- ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.58, label = "No reportable AME known-motif enrichment",
                            size = 5.2, fontface = "bold", colour = "#20252B") +
          ggplot2::annotate("text", x = 0.5, y = 0.42,
                            label = "The complete AME output and logs are retained; no motif hit is invented or hidden.",
                            size = 3.2, colour = "#606871") +
          ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
          ggplot2::labs(title = paste0(contrast_id, " — ", direction, " promoter motifs"),
                        subtitle = "JASPAR known-motif enrichment status",
                        caption = "Counts-only promoter motif is exploratory and does not establish direct TF binding.") +
          theme_publication() +
          ggplot2::theme(axis.line = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
                         axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
        save_publication_figure(status_plot, file.path(figures, paste0(stem, "_Known_Motif_Status")), 145, 85)
      }
      balance <- rbind(data.frame(set = "foreground", foreground[, c("gene_id", "gc_fraction", "sequence_length", "baseMean")]),
                       data.frame(set = "background", background[, c("gene_id", "gc_fraction", "sequence_length", "baseMean")]))
      write_result_table(balance, file.path(tables, paste0(stem, "_Matching_Balance.tsv.gz")))
      balance_long <- rbind(
        data.frame(set = balance$set, metric = "GC fraction", value = balance$gc_fraction),
        data.frame(set = balance$set, metric = "Promoter length (bp)", value = balance$sequence_length),
        data.frame(set = balance$set, metric = "log1p(baseMean)", value = log1p(balance$baseMean))
      )
      balance_long$set <- factor(balance_long$set, levels = c("background", "foreground"))
      balance_plot <- ggplot2::ggplot(balance_long, ggplot2::aes(set, value, fill = set, colour = set)) +
        ggplot2::geom_violin(trim = FALSE, alpha = 0.16, linewidth = 0.4) +
        ggplot2::geom_boxplot(width = 0.28, outlier.shape = NA, alpha = 0.55, linewidth = 0.4) +
        ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.06, height = 0, seed = 104729),
                            size = 0.75, alpha = 0.55) +
        ggplot2::facet_wrap(~ metric, scales = "free_y", nrow = 1) +
        ggplot2::scale_fill_manual(values = c(background = "#999999", foreground = direction_palette[[direction]])) +
        ggplot2::scale_colour_manual(values = c(background = "#666666", foreground = direction_palette[[direction]])) +
        ggplot2::labs(title = paste0(contrast_id, " — ", direction, " motif background matching"),
                      subtitle = "Every foreground promoter is compared with one non-DE tested promoter",
                      x = NULL, y = NULL,
                      caption = "Matched covariates reduce background bias; motif results remain exploratory.") +
        theme_publication() + ggplot2::theme(legend.position = "none")
      save_publication_figure(balance_plot, file.path(figures, paste0(stem, "_Background_Matching")), 178, 95)
      summaries[[stem]] <- data.frame(contrast_id, direction, status = "complete", foreground = nrow(foreground), background = nrow(background),
                                      foreground_gc_mean = mean(foreground$gc_fraction), background_gc_mean = mean(background$gc_fraction),
                                      foreground_log_expression_mean = mean(log1p(foreground$baseMean)),
                                      background_log_expression_mean = mean(log1p(background$baseMean)))
    }
  }
  write_result_table(do.call(rbind, summaries), file.path(tables, "Motif_Run_Summary.tsv"))
  if (length(known_results)) write_result_table(do.call(rbind, known_results), file.path(tables, "AME_JASPAR_All.tsv.gz"))
  context$state$motif <- summaries
  record_module_status(context, "06_Motif", "complete", "Up/Down promoter motif with matched tested-gene backgrounds")
}
