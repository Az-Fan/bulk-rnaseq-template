# Presentation-only RNA-seq figures: consume saved result values and the full rank.
# No statistical tests, fitting, pathway-resource downloads, or data imputation.
enrichment_plot_modes <- function(config) {
  style <- config$enrichment$visualization %||% "compact"
  if (!style %in% c("compact", "legacy", "both")) stop("Unknown enrichment visualization")
  list(compact = style %in% c("compact", "both"),
       legacy = style %in% c("legacy", "both") || isTRUE(config$migration$requires_legacy_parity))
}

enrichment_members <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) return(character())
  unique(strsplit(as.character(x), "/", fixed = TRUE)[[1L]])
}

select_enrichment_display <- function(table, mode, cutoff, limit, jaccard = .7) {
  if (!nrow(table)) return(table)
  x <- table[is.finite(table$padj) & table$padj <= cutoff, , drop = FALSE]
  if (mode == "ORA") x <- x[x$direction %in% c("Up", "Down"), , drop = FALSE]
  if (!nrow(x)) return(x)
  if (mode == "GSEA") x$direction <- ifelse(x$NES >= 0, "Up", "Down")
  x <- x[order(x$padj, x$pathway, x$direction), , drop = FALSE]
  member_column <- if (mode == "GSEA") "leadingEdge" else "overlap_genes"
  selected <- integer()
  quota <- ceiling(limit / 2)
  for (direction in c("Up", "Down")) {
    accepted <- integer()
    for (i in which(x$direction == direction)) {
      current <- enrichment_members(x[[member_column]][i])
      redundant <- any(vapply(accepted, function(k) {
        previous <- enrichment_members(x[[member_column]][k])
        union_size <- length(union(current, previous))
        union_size > 0 && length(intersect(current, previous)) / union_size > jaccard
      }, logical(1)))
      if (!redundant) accepted <- c(accepted, i)
      if (length(accepted) >= quota) break
    }
    selected <- c(selected, accepted)
  }
  x[head(selected, limit), , drop = FALSE]
}

# Build traceable display data from the SAME identifier representative used by
# the complete GSEA ranking (first count-matrix row); never choose by significance.
write_compact_enrichment <- function(context, contrast_id, de, table, database,
                                    mode, profile = "", ranks = NULL, pathways = NULL) {
  cfg <- context$config$enrichment %||% list()
  cutoff <- if (mode == "GSEA") cfg$gsea_reporting_fdr %||% .25 else cfg$ora_fdr %||% .05
  root <- comparison_module_dir(context, contrast_id, "03_Enrichment")
  label <- paste(c(mode, database, profile[nzchar(profile)]), collapse = "_")
  tables <- file.path(root, mode, "Tables", "Display", label)
  figures <- file.path(root, mode, "Figures", database)
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures, recursive = TRUE, showWarnings = FALSE)
  limit <- if (mode == "GSEA") min(18L, cfg$overview_terms %||% 18L) else
    min(16L, 2L * (cfg$ora_overview_terms %||% 8L))
  selected <- select_enrichment_display(table, mode, cutoff, limit, cfg$redundancy_jaccard %||% .7)
  status <- if (nrow(selected)) "complete" else "no_reportable_pathways"
  write_result_table(data.frame(status = status, reporting_fdr = cutoff, selected_terms = nrow(selected)),
                     file.path(tables, "Display_Status.tsv"))
  if (!nrow(selected)) return(invisible(NULL))
  contrast <- context$contrasts[context$contrasts$contrast_id == contrast_id, , drop = FALSE]
  if (nrow(contrast) != 1L) stop("Display requires exactly one declared contrast")
  numerator <- as.character(contrast$numerator); denominator <- as.character(contrast$denominator)
  factor <- as.character(contrast$factor)
  member_column <- if (mode == "GSEA") "leadingEdge" else "overlap_genes"
  selected$core_gene_count <- vapply(selected[[member_column]], function(x) length(enrichment_members(x)), integer(1))
  # A safe numeric key avoids filesystem truncation/collisions in long pathway IDs.
  selected$term_id <- paste(database, selected$pathway, selected$direction, sep = "|")
  write_result_table(selected, file.path(tables, "Selected_Pathways.tsv"))
  jobs <- list()
  add_job <- function(kind, data, stem, extra = list()) {
    path <- file.path(tables, paste0(stem, "_Display.tsv"))
    write_result_table(data, path)
    jobs[[length(jobs) + 1L]] <<- c(list(kind = kind, data = path,
      pdf = file.path(figures, paste0(stem, ".pdf")), database = database,
      contrast_label = contrast_id, numerator = numerator, denominator = denominator,
      fdr_cutoff = cutoff, profile = profile), extra)
  }
  add_job(if (mode == "GSEA") "gsea_overview" else "ora_overview", selected,
          paste0(label, "_Compact_Overview"))
  id_column <- if (mode == "GSEA" || context$config$species == "mouse") "gene_symbol" else "entrez_id"
  input_ids <- rownames(context$state$vst)
  mapped <- de[match(input_ids, de$gene_id), , drop = FALSE]
  mapped$gene_key <- as.character(mapped[[id_column]])
  mapped <- mapped[!is.na(mapped$gene_key) & nzchar(mapped$gene_key), , drop = FALSE]
  write_result_table(mapped[, c("gene_key", "gene_id", "gene_symbol", "stat", "log2FoldChange_ashr", "padj")],
                     file.path(tables, "Identifier_Representatives_All_Candidates.tsv"))
  mapped <- mapped[!duplicated(mapped$gene_key), , drop = FALSE]
  write_result_table(mapped[, c("gene_key", "gene_id", "gene_symbol")],
                     file.path(tables, "Identifier_Representatives_Selected.tsv"))
  gene_info <- function(ids) {
    m <- mapped[match(ids, mapped$gene_key), , drop = FALSE]
    symbol <- as.character(m$gene_symbol)
    symbol[is.na(symbol) | !nzchar(symbol)] <- ids[is.na(symbol) | !nzchar(symbol)]
    data.frame(gene_key = ids, gene_id = m$gene_id, gene_symbol = symbol,
               rank_metric = if (mode == "GSEA") unname(ranks[ids]) else m$stat,
               log2FC = m$log2FoldChange_ashr, gene_padj = m$padj,
               mapping_status = ifelse(is.na(m$gene_id), "unmapped", "matched"), stringsAsFactors = FALSE)
  }
  samples <- context$state$selected_samples
  groups <- as.character(context$state$samples[samples, factor])
  sample_order <- order(match(groups, c(denominator, numerator)), seq_along(samples))
  sample_order <- sample_order[groups[sample_order] %in% c(denominator, numerator)]
  samples <- samples[sample_order]; groups <- groups[sample_order]
  if (!length(samples) || !all(c(numerator, denominator) %in% groups)) stop("Contrast groups missing in expression display")
  matrix <- context$state$vst[, samples, drop = FALSE]
  write_result_table(data.frame(sample = samples, group = groups, sample_order = seq_along(samples)),
                     file.path(tables, "Sample_Order.tsv"))
  if (mode == "GSEA") {
    rank_path <- file.path(tables, "Complete_Ranks.tsv")
    write_result_table(data.frame(gene_key = names(ranks), rank_metric = unname(ranks)), rank_path)
    curves <- cfg$gsea_curves %||% 6L
    curve_rows <- selected[order(selected$padj, -abs(selected$NES), selected$pathway), , drop = FALSE]
    for (i in seq_len(min(curves, nrow(curve_rows)))) {
      row <- curve_rows[i, , drop = FALSE]
      ids <- enrichment_members(row$leadingEdge)
      info <- gene_info(ids)
      info <- info[order(-sign(row$NES) * info$rank_metric, info$gene_key, na.last = TRUE), , drop = FALSE]
      stem <- sprintf("%s_%02d_Core_Combo", label, i)
      write_result_table(info, file.path(tables, paste0(stem, "_All_Core_Genes.tsv")))
      # Constant or missing values stay NA; no artificial zero or imputation.
      shown <- head(info[info$mapping_status == "matched", , drop = FALSE], 30L)
      if (!nrow(shown)) {
        write_result_table(data.frame(status = "no_mapped_core_genes"), file.path(tables, paste0(stem, "_Status.tsv")))
        next
      }
      values <- matrix[shown$gene_id, , drop = FALSE]
      z <- t(vapply(seq_len(nrow(values)), function(k) {
        v <- values[k, ]; sd <- stats::sd(v)
        if (!is.finite(sd) || sd == 0 || any(!is.finite(v))) rep(NA_real_, length(v)) else (v - mean(v)) / sd
      }, numeric(length(samples))))
      long <- shown[rep(seq_len(nrow(shown)), each = length(samples)), , drop = FALSE]
      long$sample <- rep(samples, nrow(shown)); long$group <- rep(groups, nrow(shown))
      long$sample_order <- rep(seq_along(samples), nrow(shown))
      long$gene_order <- rep(seq_len(nrow(shown)), each = length(samples))
      long$vst <- as.vector(t(values)); long$zscore <- as.vector(t(z))
      member_path <- file.path(tables, paste0(stem, "_All_Tested_Members.tsv"))
      write_result_table(data.frame(gene_key = intersect(pathways[[row$pathway]], names(ranks))), member_path)
      add_job("core_combo", long, stem, list(ranks = rank_path, members = member_path,
        pathway = row$pathway, NES = row$NES, pathway_padj = row$padj,
        tested_members = length(intersect(pathways[[row$pathway]], names(ranks))),
        total_genes = nrow(info), shown_genes = nrow(shown)))
    }
  }
  # Keep Up and Down as separate terms even when the pathway identifier repeats.
  terms <- head(selected[order(selected$padj, selected$pathway, selected$direction), , drop = FALSE], 6L)
  memberships <- lapply(terms[[member_column]], enrichment_members)
  ids <- unique(unlist(memberships))
  if (length(ids)) {
    info <- gene_info(ids)
    info$membership_count <- vapply(ids, function(id) sum(vapply(memberships, function(m) id %in% m, logical(1))), integer(1))
    info <- info[order(-info$membership_count, -abs(info$rank_metric), info$gene_key, na.last = TRUE), , drop = FALSE]
    reserved <- unique(unlist(lapply(memberships, function(m) {
      k <- info[info$gene_key %in% m, , drop = FALSE]
      head(k$gene_key[order(-abs(k$rank_metric), k$gene_key, na.last = TRUE)], 5L)
    })))
    chosen <- head(unique(c(reserved, info$gene_key)), 40L)
    shown_ids <- info$gene_key[info$gene_key %in% chosen]
    long <- do.call(rbind, lapply(seq_len(nrow(terms)), function(i) {
      d <- info; d$term_id <- terms$term_id[i]; d$term_order <- i
      d$pathway <- terms$pathway[i]; d$term_direction <- terms$direction[i]
      d$pathway_padj <- terms$padj[i]; d$is_member <- d$gene_key %in% memberships[[i]]
      d$log2FC[!d$is_member] <- NA_real_
      d$gene_order <- match(d$gene_key, info$gene_key)
      d
    }))
    write_result_table(long, file.path(tables, "Pathway_Gene_All_Memberships.tsv"))
    add_job("gene_pathway", long[long$gene_key %in% shown_ids, , drop = FALSE],
            paste0(label, "_Pathway_Genes"), list(mode = paste(mode, profile),
            total_genes = nrow(info), shown_genes = length(shown_ids)))
  }
  write_result_table(data.frame(expression = "DESeq2 VST blind=FALSE; contrast groups only",
    scaling = "per gene sample SD (n-1); constant/missing rows stay NA; no imputation",
    clustering = "none; distance and linkage not applicable",
    gene_order = "combo: rank toward enriched end; map: membership count then absolute rank",
    identifier_rule = "first count-matrix row; candidates and chosen gene IDs retained",
    effect = "ASHR-shrunken log2FC", gene_fdr = "saved DESeq2 BH padj",
    color_limits = "[-2,2]; squish display only", stringsAsFactors = FALSE),
    file.path(tables, "Visualization_Provenance.tsv"))
  portable_jobs <- lapply(jobs, function(job) {
    prefix <- paste0(normalizePath(context$run_dir, winslash = "/"), "/")
    for (key in intersect(c("data", "pdf", "members", "ranks"), names(job))) {
      path <- normalizePath(job[[key]], winslash = "/", mustWork = FALSE)
      if (!startsWith(path, prefix)) stop("Display artifacts must stay inside the run")
      job[[key]] <- substring(path, nchar(prefix) + 1L)
    }
    job
  })
  jsonlite::write_json(portable_jobs, file.path(tables, "Visualization_Jobs.json"), auto_unbox = TRUE, pretty = TRUE)
  render_enrichment_jobs(jobs, context$run_dir)
  invisible(jobs)
}

render_enrichment_jobs <- function(jobs, run_root) {
  library(ggplot2)
  library(patchwork)
read_table <- function(p) {
  columns <- names(read.delim(p,nrows=0,check.names=FALSE))
  types <- rep(NA_character_,length(columns)); types[columns %in% c('gene_key','gene_symbol','gene_id')] <- 'character'
  read.delim(p,check.names=FALSE,colClasses=types)
}
old_theme <- theme_set(theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(),
  plot.title=element_text(face='bold'),plot.margin=margin(8,10,8,10)))
on.exit(theme_set(old_theme),add=TRUE)

nice <- function(x) {
  x <- sub('^KEGG_(mmu|hsa)[0-9]+_','',x)
  x <- sub(' - (Mus musculus \\(house mouse\\)|Homo sapiens \\(human\\))$','',x)
  y <- sub('^(HALLMARK_|GOBP_|GOCC_|GOMF_|REACTOME_|KEGG_|WIKIPATHWAYS_|WP_)','',x)
  y <- tools::toTitleCase(tolower(gsub('_',' ',y)))
  replacements <- c('Tgf Beta'='TGF-beta','Mtorc1'='mTORC1','Tnfa'='TNF-alpha',
    'Nfkb'='NF-kappaB','E2f'='E2F','G2m'='G2/M','Il2'='IL2','Stat5'='STAT5',
    'Stat3'='STAT3','Il6'='IL6','Kras'='KRAS','Myc'='MYC','Pi3k'='PI3K',
    'Akt'='AKT','Mtor'='mTOR','Dna'='DNA','Rna'='RNA','Atp'='ATP','Gtp'='GTP',
    'Vegf'='VEGF','Notch'='Notch','Mapk'='MAPK')
  for (k in names(replacements)) y <- gsub(paste0('\\b',k,'\\b'),replacements[[k]],y,ignore.case=TRUE)
  y
}
wrap <- function(x,width=42,max_lines=Inf) vapply(x,function(s) {
  lines <- strwrap(s,width=width)
  if (length(lines)>max_lines) lines <- c(lines[seq_len(max_lines-1)],paste0(lines[max_lines],' ...'))
  paste(lines,collapse='\n')
},character(1))
qformat <- function(x) ifelse(x>=.001,sprintf('%.4f',x),formatC(x,format='e',digits=1))
export <- function(p,j,height=190,width=240) {
  save_publication_figure(p,sub('\\.pdf$','',j$pdf),width,height)
  cat(basename(j$pdf),'\n')
}
fc_scale <- function(name='log2FC') scale_fill_gradient2(low='#397DA8',mid='white',high='#C54D3D',
  midpoint=0,limits=c(-2,2),oob=scales::squish,na.value='#ECECEC',name=name)

for (j in jobs) {
  for (key in intersect(c('data','pdf','members','ranks'),names(j))) {
    if (!startsWith(j[[key]],'/')) j[[key]] <- file.path(run_root,j[[key]])
  }
  d <- read_table(j$data)
  if (j$kind=='gsea_overview') {
    d <- d[order(d$NES),]
    labels <- paste0(wrap(nice(d$pathway),max_lines=3),ifelse(d$padj>.05,' \u2020',''))
    label_lines <- vapply(strsplit(labels,'\n'),length,integer(1))
    d$y <- cumsum(label_lines)-label_lines/2
    d$significance <- -log10(pmax(d$padj,1e-50))
    p <- ggplot(d,aes(NES,y)) + geom_vline(xintercept=0,colour='grey65',linetype=2) +
      geom_point(aes(size=core_gene_count,colour=significance),alpha=.95) +
      scale_y_continuous(breaks=d$y,labels=labels,limits=c(0,sum(label_lines)),expand=c(0,0)) +
      scale_colour_gradient(low='#93C5D4',high='#943746',name='-log10(FDR)') +
      scale_size_continuous(range=c(3,7),name='Core genes') +
      labs(title=paste0(j$contrast_label,': ',j$database,' pathway overview'),
        subtitle=wrap(paste0('GSEA | ',nrow(d),' selected pathways | Positive NES: ',j$numerator,' end'),width=100),
        x='Normalized enrichment score (NES)',y=NULL,
        caption=paste0('Point size: leading-edge genes. Color: pathway FDR (-log10, capped at 50).\n',
          paste0('\u2020 Exploratory: FDR > 0.05 (reporting cutoff ',j$fdr_cutoff,'). Display selection only.\nFull labels and complete results remain in the accompanying tables.'))) +
      theme(axis.text.y=element_text(colour='#202020',size=10),legend.position='right',
        plot.title=element_text(size=14),plot.subtitle=element_text(size=10),plot.caption=element_text(hjust=0,size=9))
    lines <- sum(vapply(strsplit(labels,'\n'),length,integer(1)))
    export(p,j,height=max(140,65+4.8*lines))
  } else if (j$kind=='ora_overview') {
    d$direction <- factor(d$direction,levels=c('Up','Down'))
    d <- d[order(d$direction,d$gene_ratio),]
    labels <- wrap(nice(d$pathway),max_lines=3)
    duplicate <- duplicated(labels)|duplicated(labels,fromLast=TRUE)
    labels[duplicate] <- paste0(labels[duplicate], ' [', seq_len(nrow(d))[duplicate], ']')
    d$label <- factor(seq_len(nrow(d)),levels=seq_len(nrow(d)),labels=labels)
    d$significance <- -log10(pmax(d$padj,1e-50))
    p <- ggplot(d,aes(gene_ratio,label)) + geom_point(aes(size=overlap_count,colour=significance)) +
      facet_grid(direction~.,scales='free_y',space='free_y') +
      scale_colour_gradient(low='#93C5D4',high='#943746',name='-log10(FDR)') +
      scale_size_continuous(range=c(3,7),name='Overlap genes') +
      labs(title=paste0(j$contrast_label,': ',j$database,' ORA'),
        subtitle=paste0(j$profile,' profile | Up / Down foregrounds | Pathway FDR <= ',j$fdr_cutoff),
        x='GeneRatio (overlap / mapped foreground)',y=NULL,
        caption='GeneRatio is a proportion, not an effect size. Color: pathway FDR (-log10, capped at 50).\nDisplay selection only. Full labels, tested pathways and foreground/universe tables are retained.') +
      theme(axis.text.y=element_text(size=10,colour='black'),plot.title=element_text(size=14),
        plot.subtitle=element_text(size=10),plot.caption=element_text(hjust=0,size=9),strip.background=element_rect(fill='grey95'))
    max_lines <- max(vapply(strsplit(as.character(d$label),'\n'),length,integer(1)))
    export(p,j,height=max(95,65+4.2*max_lines*nrow(d)))
  } else if (j$kind=='core_combo') {
    genes <- unique(d$gene_symbol[order(d$gene_order)])
    # Original gene identifiers remain available, including symbol collisions.
    symbols <- setNames(as.character(d$gene_symbol),as.character(d$gene_key))
    ids <- unique(as.character(d$gene_key[order(d$gene_order)]))
    lab <- symbols[ids]
    duplicate <- duplicated(lab) | duplicated(lab,fromLast=TRUE)
    lab[duplicate] <- paste0(lab[duplicate],' (',ids[duplicate],')')
    d$gene <- factor(as.character(d$gene_key),levels=rev(ids))
    d$sample <- factor(d$sample,levels=unique(d$sample[order(d$sample_order)]))
    d$group <- factor(d$group,levels=c(j$denominator,j$numerator))
    hm <- ggplot(d,aes(sample,gene,fill=zscore)) + geom_tile(colour='white',linewidth=.3) +
      facet_grid(.~group,scales='free_x',space='free_x') +
      scale_y_discrete(labels=setNames(lab,ids)) +
      scale_fill_gradient2(low='#397DA8',mid='white',high='#C54D3D',midpoint=0,
        limits=c(-2,2),oob=scales::squish,na.value='#ECECEC',name='Row z-score') +
      labs(title=paste0('Core genes: ',j$shown_genes,' of ',j$total_genes,' shown'),x=NULL,y=NULL) +
      theme(panel.grid=element_blank(),axis.text.y=element_text(size=10,face='italic',colour='black'),
        axis.text.x=element_text(size=8.5,angle=45,hjust=1),legend.position='bottom',
        strip.background=element_rect(fill='grey95'),strip.text=element_text(size=9),plot.title=element_text(size=12))
    ef <- d[!duplicated(d$gene),]
    effect <- ggplot(ef,aes(x=1,y=gene)) + geom_tile(aes(fill=log2FC),width=.85,height=.94) +
      geom_text(aes(label=sprintf('%+.2f',log2FC)),size=3) +
      geom_text(aes(x=2,label=qformat(gene_padj)),size=2.8,colour='grey25') + fc_scale() +
      scale_x_continuous(breaks=c(1,2),labels=c('log2FC','Gene FDR'),limits=c(.5,2.6),expand=c(0,0)) +
      labs(title='Gene effects',x=NULL,y=NULL) +
      theme(panel.grid=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=9),legend.position='bottom',plot.title=element_text(size=12))
    r <- read_table(j$ranks); ranks <- setNames(r$rank_metric,as.character(r$gene_key))
    m <- read_table(j$members)
    curve <- fgsea::plotEnrichment(as.character(m$gene_key),ranks) +
      labs(title=wrap(paste0(nice(j$pathway),': enrichment signal and core genes'),width=75),
        subtitle=paste0(j$contrast_label,' | NES = ',sprintf('%.2f',j$NES),
          ' | Pathway FDR = ',qformat(j$pathway_padj),if(j$pathway_padj>.05)' (exploratory)' else '',
          ' | Members = ',j$tested_members),
        x=paste0('Gene rank (',j$numerator,' end on the left)'),y='Running enrichment score') +
      theme(plot.title=element_text(size=14),plot.subtitle=element_text(size=9))
    curve$layers[[1]]$aes_params$colour <- '#237D8A'
    caption <- paste0('Core genes = saved GSEA leading edge; order: strongest rank metric toward the enriched end. ',
      if(j$total_genes>j$shown_genes)'Display subset; all core genes in accompanying table.' else 'All core genes shown.',
      '\nExpression: per-gene z-scores of model-related VST (blind = FALSE), using only the contrast groups. No clustering or imputation.',
      '\nlog2FC: ASHR-shrunken numerator vs denominator; FDR: saved DESeq2 BH-adjusted P value. Pathway FDR and gene FDR are distinct; core membership does not imply gene significance.',
      '\nColor scales capped at +/-2 for readability. Grey expression tiles indicate unavailable values.')
    p <- (curve / (hm + effect + plot_layout(widths=c(if(length(levels(d$sample))>7)4.3 else 3.3,1.7)))) +
      plot_layout(heights=c(1,max(1.3,j$shown_genes/9))) +
      plot_annotation(caption=wrap(caption,width=145),theme=theme(plot.caption=element_text(hjust=0,size=8.5),plot.margin=margin(8,8,8,8)))
    export(p,j,height=max(220,150+j$shown_genes*4.3),width=250)
  } else if (j$kind=='gene_pathway') {
    ids <- unique(as.character(d$gene_key[order(d$gene_order)]))
    symbols <- setNames(as.character(d$gene_symbol),as.character(d$gene_key))
    lab <- symbols[ids]; duplicate <- duplicated(lab)|duplicated(lab,fromLast=TRUE)
    lab[duplicate] <- paste0(lab[duplicate],' (',ids[duplicate],')')
    terms <- unique(d$term_id[order(d$term_order)])
    d$gene <- factor(as.character(d$gene_key),levels=rev(ids))
    d$term <- factor(d$term_id,levels=terms,labels=sprintf('P%02d',seq_along(terms)))
    hm <- ggplot(d,aes(term,gene,fill=log2FC)) + geom_tile(colour='white',linewidth=.4) +
      scale_y_discrete(labels=setNames(lab,ids)) + fc_scale('log2FC in member genes') +
      labs(title=paste0('Genes: ',j$shown_genes,' of ',j$total_genes,' shown'),x='Pathway key',y=NULL) +
      theme(panel.grid=element_blank(),axis.text.y=element_text(size=9,face='italic',colour='black'),
        legend.position='bottom',plot.title=element_text(size=12))
    key <- data.frame(x=0,y=seq_along(terms),
      label=paste0(sprintf('P%02d',seq_along(terms)),': ',
        if(startsWith(j$mode,'ORA'))paste0(d$term_direction[match(terms,d$term_id)],' / ') else '',
        nice(d$pathway[match(terms,d$term_id)])),
      padj=vapply(terms,function(t) d$pathway_padj[match(t,d$term_id)],numeric(1)))
    key$label <- paste0(wrap(key$label,width=42),'\nPathway FDR: ',qformat(key$padj))
    keys <- ggplot(key,aes(x,y,label=label)) + geom_text(hjust=0,vjust=.5,size=3.2,lineheight=1.1) +
      scale_y_reverse(limits=c(length(terms)+.6,.4)) + scale_x_continuous(limits=c(0,1),expand=c(0,0)) +
      labs(title='Pathway names') + theme_void() + theme(plot.title=element_text(size=12,face='bold'),plot.margin=margin(10,10,25,5))
    p <- (hm + keys + plot_layout(widths=c(1.05,1.3))) +
      plot_annotation(title=paste0(j$contrast_label,': pathway-gene map'),subtitle=j$mode,
        caption=paste0('Colored cells: actual core/overlap membership. Grey: no membership or unavailable effect. Color: saved log2FC (capped at +/-2).',
          '\nSelection reserves up to 5 genes per pathway, then fills by membership count and absolute rank metric. No clustering.',
          '\nAll genes and memberships remain in tables; rows ordered by membership count, then absolute rank metric.',
          '\nCore/overlap membership does not imply individual gene significance. Pathway FDR belongs to the original testing family.'),
        theme=theme(plot.title=element_text(size=15,face='bold'),plot.subtitle=element_text(size=10),
          plot.caption=element_text(hjust=0,size=8.5),plot.margin=margin(8,8,8,8)))
    export(p,j,height=max(160,100+j$shown_genes*4.1),width=250)
  } else stop('Unsupported plot kind')
}

}
