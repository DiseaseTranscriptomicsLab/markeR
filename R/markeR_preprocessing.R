# =============================================================================
# markeR_preprocessing.R  -  Sequential Preprocessing Tab (v7)
# =============================================================================
#' @importFrom shiny moduleServer NS reactiveValues reactiveVal reactive observe
#'   observeEvent req isolate renderPlot renderUI uiOutput plotOutput tagList
#'   tags div p h5 numericInput radioButtons selectInput checkboxInput
#'   checkboxGroupInput actionButton showNotification withProgress incProgress
#'   downloadButton downloadHandler fluidRow column HTML hr helpText validate
#'   need updateNumericInput updateSelectInput showModal modalDialog modalButton
#'   removeModal textAreaInput icon
#' @importFrom bslib card card_header layout_sidebar sidebar tooltip nav_panel
#' @importFrom ggplot2 ggplot aes geom_line geom_density geom_boxplot geom_col
#'   geom_vline geom_hline geom_point annotate scale_fill_manual theme_bw theme
#'   element_blank element_line element_text xlab ylab ggtitle labs margin
#'   theme_void ggsave scale_color_manual scale_x_continuous scale_y_continuous
#'   coord_cartesian facet_wrap
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom DESeq2 DESeqDataSetFromMatrix vst
#' @importFrom SummarizedExperiment assay
#' @importFrom limma voom lmFit removeBatchEffect
#' @importFrom DT DTOutput renderDT datatable formatStyle styleEqual
#' @importFrom base64enc dataURI
#' @importFrom grDevices png dev.off col2rgb rgb colorRampPalette
#' @importFrom stats prcomp model.matrix reformulate setNames median quantile var sd
#' @importFrom utils write.csv head tail download.file
#' @importFrom jsonlite toJSON
#' @importFrom plotly ggplotly plotly_build plotlyOutput renderPlotly
# Note: AnnotationDbi is optional (Suggests); accessed via :: with requireNamespace guard

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L) return(b)
  # Safe NA check that handles data.frames (a[1] returns a 1-column df, not a scalar)
  na_flag <- tryCatch({
    scalar <- if (is.data.frame(a)) a[[1L]][[1L]] else a[[1L]]
    isTRUE(is.na(scalar))
  }, error = function(e) FALSE)
  if (na_flag) return(b)
  a
}

# =============================================================================
# INTERNAL HELPERS
# =============================================================================

.pp_complexity_byexprgenes <- function(data, samples_to_plot,
                                        reads_in_percentage = FALSE,
                                        progress_fn = NULL) {
  n_s     <- length(samples_to_plot)
  df_list <- lapply(seq_along(samples_to_plot), function(i) {
    samp <- samples_to_plot[i]
    if (!is.null(progress_fn))
      progress_fn(i / n_s, paste0("Processing sample ", i, "/", n_s, " - ", samp))
    vals <- sort(data[, samp], decreasing = TRUE)
    data.frame(Sample = samp, CumulativeReads = cumsum(vals),
               GeneIndex = seq_along(vals), stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, df_list)
  if (!reads_in_percentage) {
    plt <- ggplot2::ggplot(df,
        ggplot2::aes(x = GeneIndex, y = log10(CumulativeReads + 1),
                     group = Sample, colour = Sample, text = Sample)) +
      ggplot2::geom_line(alpha = 0.75, linewidth = 0.7) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.border     = ggplot2::element_blank(),
                     axis.line        = ggplot2::element_line(colour = "black")) +
      ggplot2::xlab("Gene index (desc.)") + ggplot2::ylab("Cumulative reads (log₁₀)")
  } else {
    totals <- tapply(df$CumulativeReads, df$Sample, max)
    df$Pct <- 100 * df$CumulativeReads / totals[df$Sample]
    plt <- ggplot2::ggplot(df,
        ggplot2::aes(x = GeneIndex, y = Pct, group = Sample, colour = Sample, text = Sample)) +
      ggplot2::geom_line(alpha = 0.75, linewidth = 0.7) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.border     = ggplot2::element_blank(),
                     axis.line        = ggplot2::element_line(colour = "black")) +
      ggplot2::xlab("Gene index (desc.)") + ggplot2::ylab("Cumulative reads (%)")
  }
  list(plot = plt, data_to_plot = df)
}

.pp_complexity_samplingreads <- function(data, NumberGenes, NumberIterations,
                                          stepReads, samples_to_plot,
                                          reads_in_percentage = FALSE,
                                          progress_fn = NULL) {
  NumberGenes <- min(as.integer(NumberGenes), 1000L, nrow(data))
  # Only sample from genes with non-zero total counts (avoids empty vectors)
  row_totals <- rowSums(data[, samples_to_plot, drop = FALSE])
  eligible   <- which(row_totals > 0)
  if (length(eligible) == 0L) return(NULL)
  gene_idx <- sample(eligible, min(NumberGenes, length(eligible)))
  n_s <- length(samples_to_plot)
  out <- data.frame()
  for (si in seq_along(samples_to_plot)) {
    samp <- samples_to_plot[si]
    if (!is.null(progress_fn))
      progress_fn(si / n_s, paste0("Sampling reads for sample ", si, "/", n_s, " - ", samp))
    sub     <- data.frame(gene = rownames(data)[gene_idx],
                           count = data[gene_idx, samp], stringsAsFactors = FALSE)
    # Remove genes with zero counts in this sample
    sub     <- sub[sub$count > 0, , drop = FALSE]
    vec     <- unlist(mapply(rep, sub$gene, ceiling(sub$count), SIMPLIFY = FALSE))
    libsize <- length(vec)
    if (libsize < 1L) next
    vec   <- sample(vec)
    steps <- unique(c(seq(1L, libsize, by = max(1L, as.integer(stepReads))), libsize))
    n_unique <- vapply(steps, function(n)
      mean(replicate(NumberIterations,
                     length(unique(sample(vec, n, replace = FALSE))))), numeric(1))
    out <- rbind(out, data.frame(Sample = samp, NumberReads = steps,
      PercentageReads = steps / libsize, NumberUniqueGenesDetected = n_unique,
      stringsAsFactors = FALSE))
  }
  if (nrow(out) == 0) return(NULL)
  if (reads_in_percentage) {
    plt <- ggplot2::ggplot(out,
      ggplot2::aes(x = PercentageReads, y = NumberUniqueGenesDetected,
                   colour = Sample, group = Sample, text = Sample)) +
      ggplot2::geom_line() + ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::xlab("Fraction of reads") + ggplot2::ylab("Unique genes detected")
  } else {
    plt <- ggplot2::ggplot(out,
      ggplot2::aes(x = NumberReads, y = NumberUniqueGenesDetected,
                   colour = Sample, group = Sample, text = Sample)) +
      ggplot2::geom_line() + ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::xlab("Number of reads") + ggplot2::ylab("Unique genes detected")
  }
  list(plot = plt, data_to_plot = out)
}

# ── Boxplot helper (shared by Step 3 reactive renders) ────────────────────────

# samps_plot: if provided, bypasses selection logic (ensures before/after use same samples)
.pp_make_boxplot <- function(log_d, all_s, display, n_box, color_by, meta, title_str,
                              samps_plot = NULL, hide_xlabels = FALSE,
                              y_label = "log₂(count + 1)") {
  if (is.null(samps_plot)) {
    samps_plot <- switch(display,
      random  = { set.seed(42L + n_box); sample(all_s, min(n_box, length(all_s))) },
      extreme = {
        meds <- apply(log_d, 2, stats::median)
        ord  <- order(meds); nb <- floor(n_box / 2L); nt <- n_box - nb
        unique(c(all_s[utils::head(ord, nb)], all_s[utils::tail(ord, nt)]))
      },
      all_s
    )
  }
  sub_d <- log_d[, samps_plot, drop = FALSE]
  df    <- do.call(rbind, lapply(samps_plot, function(s)
    data.frame(Sample = s, LogCount = sub_d[, s], stringsAsFactors = FALSE)))
  if (!is.null(color_by) && !is.null(meta) && color_by %in% colnames(meta)) {
    id_col <- if ("SampleID" %in% colnames(meta)) "SampleID" else colnames(meta)[1]
    cmap   <- stats::setNames(as.character(meta[[color_by]]),
                               as.character(meta[[id_col]]))
    df$Group <- cmap[df$Sample]
  } else { df$Group <- df$Sample; color_by <- "Sample" }
  # When samps_plot is provided by caller, preserve that order so before/after plots align.
  # Otherwise sort by median (ascending) for a meaningful default ordering.
  if (!is.null(samps_plot)) {
    df$Sample <- factor(df$Sample, levels = samps_plot)
  } else {
    med_lvls  <- names(sort(tapply(df$LogCount, df$Sample, stats::median, na.rm = TRUE)))
    df$Sample <- factor(df$Sample, levels = med_lvls)
  }
  x_text_theme <- if (hide_xlabels)
    ggplot2::element_blank()
  else
    ggplot2::element_text(angle = 90, hjust = 1, size = 8)
  ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = LogCount, fill = Group, text = Sample)) +
    ggplot2::geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3, lwd = 0.3) +
    ggplot2::scale_fill_manual(values = stats::setNames(
      rep(.pp_palette, length.out = length(unique(df$Group))),
      unique(df$Group))) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x  = x_text_theme,
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line    = ggplot2::element_line(colour = "black"),
      legend.position = "bottom") +
    ggplot2::xlab("") + ggplot2::ylab(y_label) +
    ggplot2::ggtitle(paste0(title_str, " - ", length(samps_plot), " samples")) +
    ggplot2::labs(fill = color_by)
}

# ── PCA: compute once, render many times ──────────────────────────────────────

.pp_run_pca <- function(log_counts, n_pcs = 10L, scale = FALSE, center = TRUE) {
  pca   <- stats::prcomp(t(as.matrix(log_counts)), scale = scale, center = center)
  n_use <- min(n_pcs, ncol(pca$x))
  list(
    df = as.data.frame(pca$x[, 1:n_use, drop = FALSE]),
    ev = pca$sdev^2
  )
}

.pp_pca_from_df <- function(pca_df, ev, meta, color_by = NULL, title = "",
                              pc_x = "PC1", pc_y = "PC2") {
  df <- pca_df
  cols_avail <- colnames(df)
  if (!pc_x %in% cols_avail) pc_x <- cols_avail[1]
  if (!pc_y %in% cols_avail) pc_y <- cols_avail[min(2L, length(cols_avail))]
  pc_x_i <- as.integer(sub("PC", "", pc_x))
  pc_y_i <- as.integer(sub("PC", "", pc_y))
  pc1pct <- if (pc_x_i <= length(ev)) round(100 * ev[pc_x_i] / sum(ev), 1) else 0
  pc2pct <- if (pc_y_i <= length(ev)) round(100 * ev[pc_y_i] / sum(ev), 1) else 0

  if (!is.null(meta) && !is.null(color_by) && color_by %in% colnames(meta)) {
    id_col   <- if ("SampleID" %in% colnames(meta)) "SampleID" else colnames(meta)[1]
    meta_ids <- as.character(meta[[id_col]])
    idx      <- match(rownames(df), meta_ids)
    df$ColorVar <- as.character(meta[[color_by]])[idx]
  } else {
    df$ColorVar <- rownames(df); color_by <- "Sample"
  }
  df$x_coord   <- df[[pc_x]]
  df$y_coord   <- df[[pc_y]]
  df$SampleLbl <- rownames(df)
  # Clamp axes to data range so zero-reference lines don't collapse the plot
  xl <- range(df$x_coord, na.rm = TRUE); px <- diff(xl) * 0.06
  yl <- range(df$y_coord, na.rm = TRUE); py <- diff(yl) * 0.06
  grps <- unique(df$ColorVar)
  ggplot2::ggplot(df, ggplot2::aes(x = x_coord, y = y_coord,
                                    colour = ColorVar, text = SampleLbl)) +
    ggplot2::geom_point(size = 3.5, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
    ggplot2::scale_colour_manual(values = stats::setNames(
      rep(.pp_palette, length.out = length(grps)), grps)) +
    ggplot2::coord_cartesian(xlim = c(xl[1] - px, xl[2] + px),
                              ylim = c(yl[1] - py, yl[2] + py)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                   legend.margin = ggplot2::margin(),
                   panel.border = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")) +
    ggplot2::xlab(paste0(pc_x, ": ", pc1pct, "%")) +
    ggplot2::ylab(paste0(pc_y, ": ", pc2pct, "%")) +
    ggplot2::ggtitle(title) + ggplot2::labs(colour = color_by)
}

.pp_scree_plot <- function(ev, n_show = 10L, pc_x = "PC1", pc_y = "PC2") {
  n   <- min(n_show, length(ev))
  pct <- round(100 * ev[1:n] / sum(ev), 1)
  pc_labels <- paste0("PC", 1:n)
  sel_nums  <- suppressWarnings(as.integer(sub("PC", "", c(pc_x, pc_y))))
  sel_nums  <- sel_nums[!is.na(sel_nums)]
  df <- data.frame(Num = 1:n, Pct = pct, stringsAsFactors = FALSE)
  df$Selected <- df$Num %in% sel_nums
  ggplot2::ggplot(df, ggplot2::aes(x = .data$Num, y = .data$Pct, fill = .data$Selected)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#d1d5db", "TRUE" = "#EBB43E"),
                                guide = "none") +
    ggplot2::scale_x_continuous(breaks = 1:n, labels = 1:n) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border     = ggplot2::element_blank(),
                   axis.line        = ggplot2::element_line(colour = "black"),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x      = ggplot2::element_text(size = 8)) +
    ggplot2::xlab("PC") + ggplot2::ylab("% variance")
}

# ── Interactive report helpers ────────────────────────────────────────────────

# 5-number summary (per column) for interactive boxplot report
.pp_box_stats_report <- function(mat) {
  mat <- as.matrix(mat)
  stats::setNames(lapply(seq_len(ncol(mat)), function(j) {
    x <- mat[, j]; x <- x[is.finite(x)]
    if (length(x) == 0L) return(list(q1=0, median=0, q3=0, lf=0, uf=0))
    qs  <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    iqr <- qs[[3L]] - qs[[1L]]
    list(q1     = round(qs[[1L]], 3),
         median = round(qs[[2L]], 3),
         q3     = round(qs[[3L]], 3),
         lf     = round(max(min(x), qs[[1L]] - 1.5 * iqr), 3),
         uf     = round(min(max(x), qs[[3L]] + 1.5 * iqr), 3))
  }), colnames(mat))
}

# Shared colour palette — used in both ggplot2 (Shiny) and Plotly (HTML report)
.pp_palette <- c(
  "#E63946", "#457B9D", "#2A9D8F", "#F4A261", "#6A4C93",
  "#1D8348", "#C77DFF", "#F72585", "#4CC9F0", "#F9C74F",
  "#90BE6D", "#577590", "#7B2D8B", "#F94144", "#43AA8B",
  "#FF6B6B", "#118AB2", "#06D6A0", "#FFD166", "#EF476F"
)
.pp_report_palette <- jsonlite::toJSON(.pp_palette, auto_unbox = FALSE)

# Align metadata to a vector of sample names → named list of character vectors per column
.pp_align_meta <- function(samples, meta, meta_cols) {
  if (is.null(meta) || length(meta_cols) == 0L) return(list())
  id_col   <- if ("SampleID" %in% colnames(meta)) "SampleID" else colnames(meta)[1L]
  meta_ids <- as.character(meta[[id_col]])
  idx      <- match(samples, meta_ids)
  out <- lapply(meta_cols, function(col) {
    if (all(!is.na(idx))) as.character(meta[[col]])[idx] else rep("NA", length(samples))
  })
  stats::setNames(out, meta_cols)
}

# Interactive Plotly boxplot div with colour-by dropdown (uses pre-computed box stats)
.report_box_div <- function(box_stats, meta, meta_cols, initial_color, div_id,
                             title = "", height = "460px", y_label = "log₂(count + 1)") {
  samples <- names(box_stats)
  if (length(samples) == 0L) return("<p style='color:#9ca3af;'>No data</p>")
  meta_aln <- .pp_align_meta(samples, meta, meta_cols)
  meta_aln[["_sample_"]] <- samples
  all_cols <- c("_sample_", meta_cols)
  opts <- paste(sapply(all_cols, function(col) {
    lbl <- if (col == "_sample_") "Sample" else col
    sel <- if (identical(col, initial_color %||% "_sample_")) " selected" else ""
    sprintf('<option value="%s"%s>%s</option>', col, sel, lbl)
  }), collapse = "\n      ")
  paste0(
    '<div style="margin-bottom:8px;">',
    '<label style="font-size:0.84em;font-weight:600;margin-right:6px;">Colour by:</label>',
    '<select id="', div_id, '_sel" onchange="boxRC_', div_id, '(this.value)"',
    ' style="font-size:0.84em;padding:2px 8px;border-radius:4px;border:1px solid #d1d5db;">',
    opts, '</select></div>',
    '<div id="', div_id, '" style="width:100%;height:', height, ';"></div>',
    '<script>(function(){',
    'var S=', jsonlite::toJSON(samples, auto_unbox=FALSE), ',',
    'BS=', jsonlite::toJSON(box_stats, auto_unbox=TRUE), ',',
    'M=', jsonlite::toJSON(meta_aln, auto_unbox=FALSE), ',',
    'pal=', .pp_report_palette, ';',
    'var initCol=', jsonlite::toJSON(initial_color %||% "_sample_", auto_unbox=TRUE), ';',
    'var ttl=', jsonlite::toJSON(title, auto_unbox=TRUE), ';',
    # One trace per sample, each with a single colour derived from its group.
    # legendgroup+showlegend:first gives one legend entry per group, no boxmode needed.
    'function tr(col){',
    '  var vals=M[col]||S.map(function(_,i){return "S"+(i+1);});',
    '  var uniq=[...new Set(vals)].sort();',
    '  var cm={}; uniq.forEach(function(g,i){cm[g]=pal[i%pal.length];});',
    '  var seen={};',
    '  return S.map(function(s,i){',
    '    var g=vals[i],st=BS[s]||{q1:0,median:0,q3:0,lf:0,uf:0};',
    '    var first=!seen[g]; seen[g]=true;',
    '    return{type:"box",name:g,x:[s],',
    '      q1:[st.q1],median:[st.median],q3:[st.q3],',
    '      lowerfence:[st.lf],upperfence:[st.uf],',
    '      marker:{color:cm[g],opacity:0.85},boxpoints:false,',
    '      legendgroup:g,showlegend:first};',
    '  });',
    '}',
    'var lay={title:{text:ttl,font:{size:13}},',
    '  xaxis:{tickangle:-45,tickfont:{size:9},automargin:true},',
    '  yaxis:{title:', jsonlite::toJSON(y_label, auto_unbox=TRUE), ',titlefont:{size:11}},',
    '  showlegend:true,legend:{orientation:"h",y:-0.25,font:{size:11}},',
    '  margin:{l:55,r:10,t:ttl?40:20,b:80},',
    '  plot_bgcolor:"white",paper_bgcolor:"white"};',
    'function __init(){',
    '  Plotly.newPlot("', div_id, '",tr(initCol),lay,{responsive:true,displayModeBar:true});',
    '  window["boxRC_', div_id, '"]=function(col){Plotly.react("', div_id, '",tr(col),lay);};',
    '}',
    'if(typeof Plotly!=="undefined"){__init();}',
    'else{window.addEventListener("load",__init);}',
    '})();</script>'
  )
}

# Interactive PCA scatter div with colour-by dropdown
# Convert a ggplot to an interactive Plotly div for the HTML report.
# Falls back to a static PNG if plotly is not available or conversion fails.
.report_ggplotly_div <- function(gg, div_id, height = "420px") {
  if (is.null(gg)) return("<p style='color:#9ca3af;font-size:0.85em;'>(not available)</p>")
  if (!requireNamespace("plotly", quietly = TRUE))
    return(plot_to_img(gg))
  tryCatch({
    pl      <- plotly::ggplotly(gg, tooltip = "text")
    built   <- plotly::plotly_build(pl)$x
    pl_data <- jsonlite::toJSON(built$data,   auto_unbox = TRUE, null = "null")
    pl_lay  <- jsonlite::toJSON(built$layout, auto_unbox = TRUE, null = "null")
    paste0(
      '<div id="', div_id, '" style="width:100%;height:', height, ';"></div>',
      '<script>(function(){',
      'var d=', pl_data, ',l=', pl_lay, ';',
      'function __init(){Plotly.newPlot("', div_id, '",d,l,{responsive:true,displayModeBar:true});}',
      'if(typeof Plotly!=="undefined"){__init();}',
      'else{window.addEventListener("load",__init);}',
      '})();</script>'
    )
  }, error = function(e) plot_to_img(gg))
}

.report_pca_div <- function(pca_df, ev, meta, meta_cols, initial_color, pc_x, pc_y,
                             div_id, height = "620px") {
  if (is.null(pca_df) || nrow(pca_df) == 0L)
    return("<p style='color:#9ca3af;'>PCA not available</p>")
  samples <- rownames(pca_df)
  x_vals  <- round(pca_df[[pc_x]], 4)
  y_vals  <- round(pca_df[[pc_y]], 4)
  pc_x_i  <- as.integer(sub("PC", "", pc_x))
  pc_y_i  <- as.integer(sub("PC", "", pc_y))
  pct_x   <- if (!is.null(ev) && pc_x_i <= length(ev)) round(100 * ev[pc_x_i] / sum(ev), 1) else 0
  pct_y   <- if (!is.null(ev) && pc_y_i <= length(ev)) round(100 * ev[pc_y_i] / sum(ev), 1) else 0
  meta_aln <- .pp_align_meta(samples, meta, meta_cols)
  meta_aln[["_sample_"]] <- samples
  all_cols <- c("_sample_", meta_cols)
  opts <- paste(sapply(all_cols, function(col) {
    lbl <- if (col == "_sample_") "Sample" else col
    sel <- if (identical(col, initial_color %||% "_sample_")) " selected" else ""
    sprintf('<option value="%s"%s>%s</option>', col, sel, lbl)
  }), collapse = "\n      ")
  xlbl <- paste0(pc_x, ": ", pct_x, "%")
  ylbl <- paste0(pc_y, ": ", pct_y, "%")
  paste0(
    '<div style="margin-bottom:8px;">',
    '<label style="font-size:0.84em;font-weight:600;margin-right:6px;">Colour by:</label>',
    '<select id="', div_id, '_sel" onchange="pcaRC_', div_id, '(this.value)"',
    ' style="font-size:0.84em;padding:2px 8px;border-radius:4px;border:1px solid #d1d5db;">',
    opts, '</select></div>',
    '<div id="', div_id, '" style="width:100%;height:', height, ';"></div>',
    '<script>(function(){',
    'var S=', jsonlite::toJSON(samples, auto_unbox=FALSE), ',',
    'X=', jsonlite::toJSON(x_vals, auto_unbox=FALSE), ',',
    'Y=', jsonlite::toJSON(y_vals, auto_unbox=FALSE), ',',
    'M=', jsonlite::toJSON(meta_aln, auto_unbox=FALSE), ',',
    'pal=', .pp_report_palette, ';',
    'var initCol=', jsonlite::toJSON(initial_color %||% "_sample_", auto_unbox=TRUE), ';',
    'var xlbl=', jsonlite::toJSON(xlbl, auto_unbox=TRUE), ',',
    'ylbl=', jsonlite::toJSON(ylbl, auto_unbox=TRUE), ';',
    'function tr(col){',
    '  var vals=M[col]||S.map(function(_,i){return "S"+(i+1);});',
    '  var uniq=[...new Set(vals)].sort();',
    '  var cm={}; uniq.forEach(function(g,i){cm[g]=pal[i%pal.length];});',
    '  var base=[{type:"scatter",mode:"markers",x:X,y:Y,text:S,',
    '    hovertemplate:"%{text}<br>"+xlbl+": %{x:.2f}<br>"+ylbl+": %{y:.2f}<extra></extra>",',
    '    marker:{size:10,opacity:0.85,color:vals.map(function(v){return cm[v];}),',
    '            line:{width:0.5,color:"rgba(255,255,255,0.7)"}},showlegend:false}];',
    '  return base.concat(uniq.map(function(g,i){',
    '    return{type:"scatter",mode:"markers",name:g,x:[null],y:[null],',
    '           marker:{size:10,color:pal[i%pal.length]},showlegend:true};',
    '  }));',
    '}',
    'var lay={xaxis:{title:xlbl,zeroline:true,zerolinecolor:"#e5e7eb"},',
    '  yaxis:{title:ylbl,zeroline:true,zerolinecolor:"#e5e7eb"},',
    '  legend:{orientation:"h",y:-0.2},margin:{l:60,r:10,t:20,b:80},',
    '  plot_bgcolor:"white",paper_bgcolor:"white"};',
    'function __init(){',
    '  Plotly.newPlot("', div_id, '",tr(initCol),lay,{responsive:true,displayModeBar:true});',
    '  window["pcaRC_', div_id, '"]=function(col){Plotly.react("', div_id, '",tr(col),lay);};',
    '}',
    'if(typeof Plotly!=="undefined"){__init();}',
    'else{window.addEventListener("load",__init);}',
    '})();</script>'
  )
}

# ── Gene ID detection ──────────────────────────────────────────────────────────

# Map of species key → list(label, pkg, biomart_dataset, ensembl_prefix)
# pkg:             org.*.eg.db fallback (optional — needs separate install)
# biomart_dataset: Ensembl BioMart dataset name (primary — needs only biomaRt)
# ensembl_prefix:  regex to auto-detect species from gene IDs
.pp_species_map <- list(
  human       = list(label="Human (Homo sapiens)",         pkg="org.Hs.eg.db",   biomart_dataset="hsapiens_gene_ensembl",      ensembl_prefix="^ENSG[0-9]"),
  mouse       = list(label="Mouse (Mus musculus)",         pkg="org.Mm.eg.db",   biomart_dataset="mmusculus_gene_ensembl",     ensembl_prefix="^ENSMUSG[0-9]"),
  rat         = list(label="Rat (Rattus norvegicus)",      pkg="org.Rn.eg.db",   biomart_dataset="rnorvegicus_gene_ensembl",   ensembl_prefix="^ENSRNOG[0-9]"),
  zebrafish   = list(label="Zebrafish (Danio rerio)",      pkg="org.Dr.eg.db",   biomart_dataset="drerio_gene_ensembl",        ensembl_prefix="^ENSDARG[0-9]"),
  fly         = list(label="Fruit fly (D. melanogaster)",  pkg="org.Dm.eg.db",   biomart_dataset="dmelanogaster_gene_ensembl", ensembl_prefix="^FBgn[0-9]"),
  worm        = list(label="C. elegans",                   pkg="org.Ce.eg.db",   biomart_dataset="celegans_gene_ensembl",      ensembl_prefix="^WBGene[0-9]"),
  yeast       = list(label="Yeast (S. cerevisiae)",        pkg="org.Sc.sgd.db",  biomart_dataset="scerevisiae_gene_ensembl",   ensembl_prefix=NULL),
  arabidopsis = list(label="Arabidopsis thaliana",         pkg="org.At.tair.db", biomart_dataset="athaliana_eg_gene",          ensembl_prefix="^AT[0-9]G[0-9]"),
  chicken     = list(label="Chicken (Gallus gallus)",      pkg="org.Gg.eg.db",   biomart_dataset="ggallus_gene_ensembl",       ensembl_prefix="^ENSGALG[0-9]"),
  pig         = list(label="Pig (Sus scrofa)",             pkg="org.Ss.eg.db",   biomart_dataset="sscrofa_gene_ensembl",       ensembl_prefix="^ENSSSCG[0-9]"),
  bovine      = list(label="Bovine (Bos taurus)",          pkg="org.Bt.eg.db",   biomart_dataset="btaurus_gene_ensembl",       ensembl_prefix="^ENSBTAG[0-9]"),
  macaque     = list(label="Macaque (M. mulatta)",         pkg="org.Mmu.eg.db",  biomart_dataset="mmulatta_gene_ensembl",      ensembl_prefix="^ENSMMUG[0-9]"),
  dog         = list(label="Dog (Canis lupus familiaris)", pkg="org.Cf.eg.db",   biomart_dataset="cfamiliaris_gene_ensembl",   ensembl_prefix="^ENSCAFG[0-9]"),
  horse       = list(label="Horse (Equus caballus)",       pkg=NULL,             biomart_dataset="ecaballus_gene_ensembl",     ensembl_prefix="^ENSECAG[0-9]"),
  sheep       = list(label="Sheep (Ovis aries)",           pkg=NULL,             biomart_dataset="oaries_gene_ensembl",        ensembl_prefix="^ENSOARG[0-9]"),
  cat         = list(label="Cat (Felis catus)",            pkg=NULL,             biomart_dataset="fcatus_gene_ensembl",        ensembl_prefix="^ENSFCAG[0-9]"),
  rabbit      = list(label="Rabbit (Oryctolagus cuniculus)",pkg=NULL,            biomart_dataset="ocuniculus_gene_ensembl",    ensembl_prefix="^ENSOCUG[0-9]")
)

# Convert gene IDs to symbols via biomaRt (queries Ensembl API — no per-species install needed)
.pp_convert_biomart <- function(ids_clean, keytype, biomart_dataset) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) return(NULL)
  bm_filter <- switch(keytype,
    ENSEMBL  = "ensembl_gene_id",
    ENTREZID = "entrezgene_id",
    REFSEQ   = "refseq_mrna",
    NULL)
  if (is.null(bm_filter) || is.null(biomart_dataset)) return(NULL)
  mart <- tryCatch(
    suppressMessages(biomaRt::useEnsembl("genes", dataset = biomart_dataset)),
    error = function(e) NULL)
  if (is.null(mart)) return(NULL)
  res <- tryCatch(
    suppressMessages(biomaRt::getBM(
      attributes = c(bm_filter, "external_gene_name"),
      filters    = bm_filter,
      values     = ids_clean,
      mart       = mart)),
    error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0L) return(NULL)
  # Return named vector: input_id → symbol (first hit per ID)
  res <- res[res$external_gene_name != "", , drop = FALSE]
  res <- res[!duplicated(res[[bm_filter]]), , drop = FALSE]
  stats::setNames(res$external_gene_name, res[[bm_filter]])[ids_clean]
}

# Convert via org.*.eg.db (requires per-species package installed)
.pp_convert_orgdb <- function(ids_clean, keytype, pkg) {
  if (is.null(pkg) || !requireNamespace(pkg, quietly = TRUE)) return(NULL)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) return(NULL)
  org_db <- tryCatch(get(pkg, envir = asNamespace(pkg)), error = function(e) NULL)
  if (is.null(org_db)) return(NULL)
  tryCatch(
    suppressMessages(AnnotationDbi::mapIds(org_db, keys = ids_clean,
                                            column = "SYMBOL", keytype = keytype,
                                            multiVals = "first")),
    error = function(e) NULL)
}

.pp_detect_gene_id_type <- function(gene_ids) {
  s <- utils::head(gene_ids[nchar(gene_ids) > 0], 30)
  if (length(s) == 0L) return("symbol")
  # Check species-specific Ensembl prefixes first
  for (sp in names(.pp_species_map)) {
    pfx <- .pp_species_map[[sp]]$ensembl_prefix
    if (!is.null(pfx) && mean(grepl(pfx, s)) > 0.5)
      return(paste0("ensembl_", sp))
  }
  if (mean(grepl("^ENS[A-Z]+[0-9]", s)) > 0.5) return("ensembl_other")
  if (mean(grepl("^(NM_|NR_|XM_|XR_)[0-9]", s)) > 0.5) return("refseq")
  if (mean(grepl("^[0-9]+$", s)) > 0.5) return("entrez")
  return("symbol")
}

# Infer default species key from detected id_type
.pp_default_species <- function(id_type) {
  if (grepl("^ensembl_", id_type)) {
    sp <- sub("^ensembl_", "", id_type)
    if (sp %in% names(.pp_species_map)) return(sp)
  }
  "human"
}

# ── Colour warning ─────────────────────────────────────────────────────────────

.pp_color_warn <- function(meta, color_by, n_thresh = 10L) {
  if (is.null(meta) || is.null(color_by) || !color_by %in% colnames(meta)) return(NULL)
  vals   <- meta[[color_by]]
  if (is.numeric(vals)) return(NULL)
  n_lvls <- length(unique(vals[!is.na(vals)]))
  if (n_lvls > n_thresh)
    return(shiny::div(class = "alert alert-warning",
      style = "font-size:0.82em;padding:6px 10px;margin:6px 0 0 0;",
      shiny::HTML(paste0("⚠️ <b>", color_by, "</b> has <b>", n_lvls,
        "</b> unique values. Choose a variable with ≤", n_thresh, " levels."))))
  NULL
}

# ── Done card ─────────────────────────────────────────────────────────────────

.pp_done_card <- function(step_num, step_name, summary_text, back_btn_id,
                          plot_ui = NULL, note = NULL) {
  shiny::tags$details(
    style = "border:1px solid #dee2e6;border-radius:6px;margin-bottom:12px;background:#fff;",
    shiny::tags$summary(
      style = paste0("display:flex;align-items:center;justify-content:space-between;",
        "padding:10px 14px;cursor:pointer;list-style:none;",
        "background:#f8f9fa;border-radius:6px;user-select:none;"),
      shiny::tags$div(style = "display:flex;align-items:center;gap:8px;",
        shiny::tags$span(class = "pp-num pp-num-done",
          style = "width:20px;height:20px;font-size:0.72em;flex-shrink:0;", "✓"),
        shiny::tags$strong(paste0("Step ", step_num, " - ", step_name)),
        shiny::tags$span(style = "font-size:0.82em;color:#28a745;margin-left:4px;", summary_text)
      ),
      shiny::tags$span(style = "font-size:0.78em;color:#6c757d;", "▼ expand to review")
    ),
    shiny::tags$div(style = "padding:10px 14px;border-top:1px solid #dee2e6;",
      if (!is.null(note) && nchar(trimws(note)) > 0L)
        shiny::div(style = "margin-bottom:8px;padding:7px 10px;background:#fffbeb;border:1px solid #fde68a;border-radius:4px;font-size:0.83em;color:#374151;",
          shiny::tags$span(style = "font-weight:600;color:#92400e;margin-right:4px;", "\U0001f4dd Note:"),
          shiny::HTML(gsub("\n", "<br>", .html_esc(note)))),
      if (!is.null(plot_ui)) shiny::div(style = "margin-bottom:10px;", plot_ui),
      shiny::actionButton(back_btn_id, paste0("↩ Edit Step ", step_num),
                           class = "btn-sm btn-outline-secondary")
    )
  )
}

# ── HTML escaping helper (shiny::htmlEscape is not exported) ──────────────────
.html_esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

# ── Remove suspicious samples section ─────────────────────────────────────────

.pp_remove_section_ui <- function(picker_id, btn_id) {
  shiny::tags$details(
    style = "margin-top:12px;border:1px dashed #d1d5db;border-radius:6px;padding:8px 12px;",
    shiny::tags$summary(
      style = "cursor:pointer;font-size:0.84em;color:#6b7280;user-select:none;",
      "⚠️ Remove suspicious samples (optional)"
    ),
    shiny::div(style = "margin-top:8px;",
      shiny::uiOutput(picker_id),
      shiny::actionButton(btn_id, "Remove selected & update data",
                           class = "btn-sm btn-outline-danger", style = "margin-top:6px;")
    )
  )
}

# ── Native listbox ─────────────────────────────────────────────────────────────

.pp_native_picker <- function(input_id, choices) {
  shiny::selectInput(inputId = input_id, label = "Select samples to remove:",
    choices = choices, selected = NULL, multiple = TRUE,
    selectize = FALSE, size = min(6L, max(3L, length(choices))), width = "100%")
}

# ── Expand button ──────────────────────────────────────────────────────────────

.pp_expand_btn <- function(btn_id, label = "⤢ Expand / Download") {
  shiny::div(style = "display:flex;align-items:center;justify-content:flex-end;gap:8px;margin-top:4px;",
    shiny::tags$small(style = "color:#9ca3af;font-style:italic;",
      "Full-size interactive view & download →"),
    shiny::actionButton(btn_id, label, class = "btn-outline-primary btn-sm",
      style = "font-size:0.78em;padding:2px 10px;"))
}

# ── Plot modal (uses plotly for interactive/resizable) ─────────────────────────

.pp_plot_modal <- function(plot_ns_id, dl_btn_id) {
  use_plotly <- requireNamespace("plotly", quietly = TRUE)
  shiny::modalDialog(
    size = "xl", easyClose = TRUE,
    if (use_plotly) plotly::plotlyOutput(plot_ns_id, height = "72vh")
    else            shiny::plotOutput(plot_ns_id,   height = "72vh"),
    footer = shiny::tagList(
      shiny::downloadButton(dl_btn_id, "⬇ Download PNG"),
      shiny::modalButton("Close")
    )
  )
}

# ── Slow-render warning ────────────────────────────────────────────────────────

.pp_slow_warn <- function(text = NULL) {
  msg <- text %||% "⏳ Rendering may take a moment for large datasets - results will appear when ready."
  shiny::div(class = "pp-slow-warn", msg)
}

# =============================================================================
# UI
# =============================================================================
#' @export
preprocessingUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::tags$style(shiny::HTML("
      .pp-num { display:inline-flex;align-items:center;justify-content:center;
        width:22px;height:22px;border-radius:50%;font-weight:700;font-size:0.75em;
        color:#fff;flex-shrink:0;position:relative;z-index:2; }
      .pp-num-active { background:#EBB43E; }
      .pp-num-done   { background:#28a745; }
      .pp-num-locked { background:#9ca3af; }
      .pp-locked { filter:blur(3px);opacity:0.38;pointer-events:none!important;user-select:none; }
      .pp-step-card  { margin-bottom:12px; }
      .pp-hint { font-size:0.82em;color:#6c757d;font-style:italic;margin:4px 0 8px 0; }
      .pp-run-btn    { margin:10px 0 4px 0; }
      .pp-proceed-btn { margin-top:14px; }
      .pp-info-note { font-size:0.8em;color:#EBB43E;margin:2px 0 6px 0; }
      .pp-slow-warn { font-size:0.81em;color:#6b7280;background:#f9fafb;border-radius:4px;
        padding:5px 10px;margin-bottom:8px;border-left:3px solid #EBB43E; }
      details summary::-webkit-details-marker { display:none; }
      details > summary { outline:none; }
      .pp-main-scroll { overflow-y:auto;max-height:calc(100vh - 80px);padding-right:4px; }
      /* Ensure the bslib sidebar panel scrolls and never hides behind the fixed log bar */
      .bslib-sidebar-layout > .sidebar,
      .bslib-sidebar-layout > .sidebar > .sidebar-content {
        overflow-y: auto !important;
        padding-bottom: 64px !important;
      }
      .pp-summary-card { background:#fff;border:2px solid #EBB43E;border-radius:8px;
        padding:16px 18px;margin-top:8px; }
      .pp-summary-row { display:flex;align-items:flex-start;gap:8px;padding:3px 0;
        font-size:0.84em;border-bottom:1px solid #f3f4f6; }
      .pp-summary-row:last-child { border-bottom:none; }
      .pp-summary-label { font-weight:600;color:#374151;min-width:140px;flex-shrink:0; }
      .pp-dl-btn { width:100%;margin-top:4px;font-size:0.78em; }
      .pp-geneid-banner { background:#eef2ff;border:1px solid #c7d2fe;border-radius:8px;
        padding:14px 16px;margin-bottom:12px; }
      /* Progress tracker */
      .pp-progress-item { display:flex;align-items:flex-start;gap:8px; }
      .pp-progress-dot-col { display:flex;flex-direction:column;align-items:center;flex-shrink:0; }
      .pp-progress-connector { width:2px;flex:1;min-height:14px;background:#dee2e6;margin:2px 0; }
      .pp-progress-connector.done { background:#28a745; }
      .pp-progress-label { padding-bottom:14px;width:100%; }
      /* Norm status badge */
      .pp-norm-badge { display:inline-block;font-size:0.75em;padding:2px 7px;
        border-radius:12px;background:#d1fae5;color:#065f46;font-weight:600;margin-left:6px; }
      /* Loading spinner overlay on recalculating plots */
      .shiny-plot-output.recalculating {
        position: relative; min-height: 60px;
      }
      .shiny-plot-output.recalculating::after {
        content: '';
        position: absolute; top: 0; left: 0; right: 0; bottom: 0;
        background: rgba(255,255,255,0.78);
        z-index: 10;
      }
      .shiny-plot-output.recalculating::before {
        content: 'Loading...';
        position: absolute; top: 50%; left: 50%;
        transform: translate(-50%, -50%);
        background: white; border: 1px solid #dee2e6; border-radius: 6px;
        padding: 8px 16px; font-size: 0.88em; color: #6b7280;
        z-index: 11; white-space: nowrap;
        box-shadow: 0 2px 8px rgba(0,0,0,0.08);
      }
    ")),

    bslib::layout_sidebar(fill = FALSE,
      sidebar = bslib::sidebar(width = 230, open = TRUE,
        style = "overflow-y: auto; max-height: 100vh; padding-bottom: 55px;",
        shiny::tags$h6("Progress", style = "font-weight:700;margin-bottom:6px;"),
        shiny::uiOutput(ns("sidebar_progress")),
        shiny::hr(style = "margin:8px 0;"),
        shiny::uiOutput(ns("sidebar_data_info")),
        shiny::hr(style = "margin:8px 0;"),
        shiny::actionButton(ns("reset_preprocessing"), "↺ Reset all",
          class = "btn-outline-danger btn-sm w-100", style = "margin-bottom:4px;",
          title = "Reset all preprocessing steps and restore original data"),
        shiny::actionButton(ns("skip_preprocessing"), "Skip (use data as-is)",
          class = "btn-outline-secondary btn-sm w-100",
          title = "Bypass preprocessing - imported data flows to all analyses"),
        shiny::uiOutput(ns("sidebar_downloads")),
        # Spacer so content is never obscured by the fixed global log bar
        shiny::div(style = "min-height: 60px; flex-shrink: 0;")
      ),
      shiny::div(class = "pp-main-scroll", style = "padding:10px;",
        shiny::uiOutput(ns("gene_id_banner")),
        shiny::uiOutput(ns("gene_id_map_ui")),
        shiny::uiOutput(ns("step1_card")),
        shiny::uiOutput(ns("step2_card")),
        shiny::uiOutput(ns("step3_card")),
        shiny::uiOutput(ns("step4_card")),
        shiny::uiOutput(ns("pp_summary_card")),
        shiny::uiOutput(ns("pp_data_tables"))
      )
    )
  )
}

# =============================================================================
# SERVER
# =============================================================================
#' @export
preprocessingServer <- function(id, get_expr, get_meta, log_fn = NULL, get_log = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns   <- session$ns
    glog <- function(msg) if (!is.null(log_fn)) log_fn(msg)

    # ── Startup: check gene ID conversion backends ────────────────────────────
    local({
      has_biomart  <- requireNamespace("biomaRt",      quietly = TRUE)
      has_annotdbi <- requireNamespace("AnnotationDbi", quietly = TRUE)
      if (!has_biomart && !has_annotdbi) {
        glog(paste0("[WARN] Neither biomaRt nor AnnotationDbi is installed. ",
                    "Gene ID conversion will not work. ",
                    "Recommended: BiocManager::install('biomaRt')"))
      } else if (!has_biomart) {
        glog(paste0("[INFO] biomaRt not installed — gene ID conversion will use local ",
                    "org.*.eg.db packages only. For full species support install: ",
                    "BiocManager::install('biomaRt')"))
      }
    })

    modal_plot        <- shiny::reactiveVal(NULL)
    .data_fingerprint <- shiny::reactiveVal(NULL)  # tracks identity of loaded data
    .finalized_fp     <- shiny::reactiveVal(NULL)   # fingerprint at finalization (fed-back output)

    rv <- shiny::reactiveValues(
      step              = 1L, active_step = 1L,
      data_s0           = NULL, data_s1 = NULL, data_s2 = NULL,
      norm_linear       = NULL, data_s3 = NULL,
      meta_s            = NULL, final_data = NULL, finalized = 0L,
      orig_n_genes      = NULL, orig_n_samples = NULL,
      s1_summary        = "", s2_summary = "", s3_summary = "",
      s1_run            = FALSE, s2_run = FALSE,
      s3_norm_done      = FALSE, s3_run = FALSE, s3_method_used = NULL,
      s4_run            = FALSE, s4_bc_run = FALSE,
      s4_pca_scale      = FALSE, s4_pca_center = TRUE,
      s1_plot_a         = NULL, s1_plot_b = NULL,
      s4_pca_df_before  = NULL, s4_pca_ev_before = NULL,
      s4_pca_df_after   = NULL, s4_pca_ev_after  = NULL,
      s4_bc_result      = NULL,
      bc_batch_cols_used  = NULL,
      bc_effect_cols_used = NULL,
      report_s2_before  = NULL, report_s2_after   = NULL,
      report_s3_before  = NULL, report_s3_after   = NULL,
      report_s4_before  = NULL, report_s4_after   = NULL,
      s2_thr            = 1.0,
      gene_id_type      = "symbol",
      gene_id_dismissed = FALSE,
      gene_id_map       = NULL,
      data_s0_orig      = NULL,   # original data before any gene ID conversion
      downloads_available = FALSE,
      pp_summary        = NULL,
      note_s1 = "", note_s2 = "", note_s3 = "", note_s4 = "",
      report_interactive_data = NULL
    )

    shiny::observe({
      e <- get_expr(); m <- get_meta(); shiny::req(e, m)
      # Fingerprint: nrow × ncol + first/last colname (fast, no hashing)
      fp <- paste(nrow(e), ncol(e), colnames(e)[1], colnames(e)[ncol(e)], sep="|")
      prev_fp <- shiny::isolate(.data_fingerprint())
      .data_fingerprint(fp)  # always keep fingerprint current

      # Guard: if already finalized and the fingerprint matches the finalization fingerprint,
      # this change is just our own output being fed back by the parent app. Skip reset.
      if (shiny::isolate(rv$finalized) > 0L) {
        fin_fp <- shiny::isolate(.finalized_fp())
        if (!is.null(fin_fp) && fp == fin_fp) return()
        # Fingerprint differs from finalized data - genuine new dataset loaded; fall through to reset.
      }

      d0 <- as.data.frame(e)
      if (!is.null(prev_fp) && fp != prev_fp) {
        # New dataset loaded - reset all preprocessing state to avoid stale-data crashes
        rv$step <- 1L; rv$active_step <- 1L
        rv$data_s2 <- NULL; rv$norm_linear <- NULL; rv$data_s3 <- NULL; rv$final_data <- NULL
        rv$s1_run <- FALSE; rv$s2_run <- FALSE; rv$s3_run <- FALSE; rv$s3_norm_done <- FALSE
        rv$s4_run <- FALSE; rv$s4_bc_run <- FALSE; rv$s2_thr <- 1.0
        rv$s1_plot_a <- NULL; rv$s1_plot_b <- NULL
        rv$s4_pca_df_before <- NULL; rv$s4_pca_ev_before <- NULL
        rv$s4_pca_df_after  <- NULL; rv$s4_pca_ev_after  <- NULL
        rv$s4_bc_result <- NULL; rv$pp_summary <- NULL; rv$finalized <- 0L
        rv$downloads_available <- FALSE; rv$gene_id_dismissed <- FALSE; rv$gene_id_map <- NULL
        rv$s1_summary <- ""; rv$s2_summary <- ""; rv$s3_summary <- ""
        .finalized_fp(NULL)
        # (no notification — silent reset to avoid noise during normal data-tab usage)
      }
      rv$data_s0 <- d0; rv$data_s1 <- d0; rv$meta_s <- m
      rv$data_s0_orig <- d0   # keep an unmodified copy for full reset
      rv$orig_n_genes <- nrow(d0); rv$orig_n_samples <- ncol(d0)
      rv$gene_id_type <- .pp_detect_gene_id_type(rownames(d0))
      if (is.null(prev_fp)) rv$gene_id_dismissed <- FALSE  # only clear on first load
    })

    # =========================================================================
    # GENE ID DETECTION BANNER
    # =========================================================================

    output$gene_id_banner <- shiny::renderUI({
      shiny::req(rv$data_s0)
      id_type <- rv$gene_id_type
      if (id_type == "symbol" || rv$gene_id_dismissed) return(NULL)
      type_label <- if (grepl("^ensembl_", id_type)) {
        sp <- sub("^ensembl_", "", id_type)
        lbl <- .pp_species_map[[sp]]$label %||% sp
        paste0("Ensembl (", lbl, ")")
      } else switch(id_type,
        ensembl_other = "Ensembl (other species)",
        refseq        = "RefSeq - NM_, NR_, ...",
        entrez        = "Entrez Gene IDs",
        "unknown format"
      )
      default_species <- .pp_default_species(id_type)
      sp_choices <- stats::setNames(
        names(.pp_species_map),
        vapply(.pp_species_map, `[[`, character(1), "label"))
      shiny::div(class = "pp-geneid-banner",
        shiny::tags$h6(style = "margin-bottom:6px;font-weight:700;color:#3730a3;",
          "\U0001f9ec Gene ID format detected: ", type_label),
        shiny::tags$p(style = "font-size:0.86em;margin-bottom:8px;",
          "Your gene identifiers appear to be in ", shiny::tags$b(type_label), " format. ",
          "Would you like to convert them to gene symbols?",
          shiny::tags$br(),
          shiny::tags$span(style = "color:#6b7280;font-size:0.9em;",
            "Version numbers (e.g. .2 in ENSG00000141736.2) will be stripped automatically.")),
        shiny::selectInput(ns("gene_id_species"), "Species:",
          choices = sp_choices, selected = default_species, width = "300px"),
        shiny::div(style = "display:flex;gap:8px;margin-top:8px;",
          shiny::actionButton(ns("convert_gene_ids"), "✓ Convert to gene symbols",
                               class = "btn-primary btn-sm"),
          shiny::actionButton(ns("skip_gene_id"), "✗ Keep original IDs",
                               class = "btn-outline-secondary btn-sm")
        )
      )
    })

    shiny::observeEvent(input$skip_gene_id, {
      rv$gene_id_dismissed <- TRUE
      glog("Gene ID conversion skipped - keeping original identifiers.")
    })

    shiny::observeEvent(input$convert_gene_ids, {
      shiny::req(rv$data_s0)
      id_type <- shiny::isolate(rv$gene_id_type)
      species <- shiny::isolate(input$gene_id_species %||% .pp_default_species(id_type))
      sp_info <- .pp_species_map[[species]] %||% list()
      shiny::withProgress(message = "\U0001f9ec Converting gene IDs to symbols…", value = 0, {
        gene_ids  <- rownames(rv$data_s0)
        ids_clean <- sub("\\.[0-9]+$", "", gene_ids)
        keytype   <- if (grepl("^ensembl", id_type)) "ENSEMBL"
                     else switch(id_type, refseq = "REFSEQ", entrez = "ENTREZID", "SYMBOL")

        # ── 1. Try biomaRt (no per-species install required) ─────────────────
        symbols <- NULL
        if (requireNamespace("biomaRt", quietly = TRUE)) {
          shiny::showModal(shiny::modalDialog(
            title = NULL, footer = NULL, easyClose = FALSE,
            shiny::tags$div(
              style = "text-align:center; padding:40px 20px;",
              shiny::tags$div(class = "spinner-border text-primary",
                style = "width:3.5rem; height:3.5rem;", role = "status",
                shiny::tags$span(class = "visually-hidden", "Loading…")),
              shiny::tags$h5(
                paste0("Querying Ensembl BioMart for ", sp_info$label %||% species, "…"),
                style = "margin-top:18px; color:#333;"),
              shiny::tags$p(
                "Mapping gene IDs to symbols via Ensembl — this may take 15–30 s.",
                style = "color:#777; font-size:0.9em; margin-top:6px;")
            )
          ))
          glog(paste0("Gene ID conversion: trying biomaRt (", sp_info$biomart_dataset, ")…"))
          symbols <- .pp_convert_biomart(ids_clean, keytype, sp_info$biomart_dataset)
          shiny::removeModal()
          if (!is.null(symbols))
            glog("Gene ID conversion: biomaRt query succeeded.")
          else
            glog("Gene ID conversion: biomaRt returned no results, trying local db…")
        }

        # ── 2. Fall back to org.*.eg.db if biomaRt unavailable/failed ────────
        if (is.null(symbols)) {
          pkg <- sp_info$pkg
          if (!is.null(pkg) && requireNamespace(pkg, quietly = TRUE) &&
              requireNamespace("AnnotationDbi", quietly = TRUE)) {
            shiny::incProgress(0.2, paste0("Using local annotation (", pkg, ")…"))
            glog(paste0("Gene ID conversion: using ", pkg, "."))
            symbols <- .pp_convert_orgdb(ids_clean, keytype, pkg)
          }
        }

        # ── 3. Neither method worked ──────────────────────────────────────────
        if (is.null(symbols)) {
          has_biomart <- requireNamespace("biomaRt", quietly = TRUE)
          install_hint <- if (!has_biomart)
            "BiocManager::install('biomaRt')  # works for all species, recommended"
          else
            paste0("BiocManager::install('", sp_info$pkg %||% "org.*.eg.db", "')")
          msg <- paste0("Gene ID conversion failed for species '", species, "'. ",
                        "Install biomaRt (recommended) or the species-specific package: ",
                        install_hint)
          glog(paste0("[ERROR] ", msg))
          shiny::showNotification(
            shiny::HTML(paste0(
              "<b>Conversion failed.</b> Install biomaRt (works for all species):<br>",
              "<code>BiocManager::install('biomaRt')</code>")),
            type = "error", duration = 30)
          return(NULL)
        }

        shiny::incProgress(0.5, "Updating gene identifiers…")
        new_ids <- ifelse(is.na(symbols) | symbols == "", ids_clean, as.character(symbols))
        dups <- duplicated(new_ids)
        if (any(dups)) {
          for (dv in unique(new_ids[dups])) {
            idx <- which(new_ids == dv); new_ids[idx] <- paste0(dv, "_dup", seq_along(idx))
          }
        }
        shiny::incProgress(0.5, "Updating gene identifiers…")
        mapped_flag <- !is.na(symbols) & symbols != ""
        n_mapped    <- sum(mapped_flag)
        n_dropped   <- sum(is.na(symbols) | symbols == "")
        # Store mapping table for display
        rv$gene_id_map <- data.frame(
          original_id = gene_ids,
          cleaned_id  = ids_clean,
          symbol      = ifelse(is.na(symbols), "(not found)", as.character(symbols)),
          new_id      = new_ids,   # actual rowname used after conversion
          status      = ifelse(mapped_flag, "mapped", "not found"),
          stringsAsFactors = FALSE
        )
        rownames(rv$data_s0) <- new_ids; rownames(rv$data_s1) <- new_ids
        rv$gene_id_dismissed <- TRUE
        glog(paste0("Gene ID Conversion: ", n_mapped, "/", length(gene_ids),
                     " IDs mapped from ", id_type, " to gene symbols (", species, "). ",
                     n_dropped, " not found."))
        shiny::showNotification(
          shiny::HTML(paste0("✅ ", format(n_mapped, big.mark=","), " IDs converted; ",
            n_dropped, " not found (kept as original ID).")),
          type = "default", duration = 8)
        shiny::incProgress(0.1, "Done.")
      })
    })

    output$gene_id_map_ui <- shiny::renderUI({
      shiny::req(!is.null(rv$gene_id_map))
      mp <- rv$gene_id_map
      n_mapped   <- sum(mp$status == "mapped")
      n_notfound <- sum(mp$status == "not found")
      n_removed  <- sum(mp$status == "removed")
      n_kept_orig <- sum(mp$status == "kept (original ID)")
      # Summary line: show what happened
      summary_parts <- paste0(format(n_mapped, big.mark=","), " mapped")
      if (n_notfound > 0)  summary_parts <- paste0(summary_parts, ", ", n_notfound, " not found")
      if (n_removed > 0)   summary_parts <- paste0(summary_parts, ", ", n_removed, " removed")
      if (n_kept_orig > 0) summary_parts <- paste0(summary_parts, ", ", n_kept_orig, " kept with original ID")
      shiny::div(style="margin-bottom:12px;",
        shiny::div(
          style="background:#f0f9ff;border:1px solid #bae6fd;border-radius:8px;padding:10px 14px;margin-bottom:8px;",
          shiny::tags$span(style="font-weight:600;color:#0369a1;",
            "\U0001f9ec Gene ID conversion complete: "),
          paste0(summary_parts, ".")
        ),
        if (n_notfound > 0) shiny::div(
          style="background:#fefce8;border:1px solid #fde68a;border-radius:8px;padding:10px 14px;margin-bottom:8px;",
          shiny::tags$span(style="font-weight:600;color:#854d0e;",
            paste0(n_notfound, " genes could not be mapped.")),
          shiny::tags$p(style="font-size:0.86em;margin:6px 0 8px;",
            "Choose what to do with them:"),
          shiny::div(style="display:flex;gap:8px;flex-wrap:wrap;",
            shiny::actionButton(ns("geneid_keep_unmatched"), "Keep (use original ID)",
              class="btn-sm btn-outline-secondary"),
            shiny::actionButton(ns("geneid_remove_unmatched"),
              paste0("Remove (drop ", n_notfound, " genes)"),
              class="btn-sm btn-outline-danger")
          )
        ),
        shiny::tags$details(
          shiny::tags$summary(
            style="cursor:pointer;font-size:0.86em;color:#6b7280;padding:4px 0;",
            "Show full mapping table"),
          DT::DTOutput(ns("tbl_gene_id_map"), width="100%")
        )
      )
    })

    shiny::observeEvent(input$geneid_keep_unmatched, {
      # Do nothing - unmatched genes already kept; just dismiss by marking map as handled
      shiny::showNotification(
        shiny::HTML(paste0("✅ Unmatched genes kept with original IDs.")),
        type="default", duration=4)
      # Update the map so the choice buttons disappear
      mp <- rv$gene_id_map
      if (!is.null(mp)) {
        mp$status <- ifelse(mp$status == "not found", "kept (original ID)", mp$status)
        rv$gene_id_map <- mp
      }
    })

    shiny::observeEvent(input$geneid_remove_unmatched, {
      shiny::req(rv$data_s0, rv$gene_id_map)
      mp <- rv$gene_id_map
      # Use new_id (the actual current rownames) not original_id (pre-conversion)
      keep_ids <- mp$new_id[mp$status == "mapped"]
      n_before <- nrow(rv$data_s0)
      n_drop   <- sum(mp$status == "not found")
      if (length(keep_ids) == 0L) {
        shiny::showNotification("No mapped genes to keep — aborting remove.", type="error"); return()
      }
      rv$data_s0 <- rv$data_s0[rownames(rv$data_s0) %in% keep_ids, , drop=FALSE]
      rv$data_s1 <- rv$data_s1[rownames(rv$data_s1) %in% keep_ids, , drop=FALSE]
      # Update map status
      mp$status[mp$status == "not found"] <- "removed"
      rv$gene_id_map <- mp
      glog(paste0("Gene ID: ", n_drop, " unmatched genes removed. ",
                  format(nrow(rv$data_s0), big.mark=","), " genes remaining."))
      shiny::showNotification(
        shiny::HTML(paste0("🗑️ ", n_drop, " unmatched genes removed. ",
                           format(nrow(rv$data_s0), big.mark=","), " genes remaining.")),
        type="warning", duration=6)
    })

    output$tbl_gene_id_map <- DT::renderDT({
      shiny::req(rv$gene_id_map)
      mp <- rv$gene_id_map
      all_statuses <- c("mapped", "not found", "kept (original ID)", "removed")
      mp$status <- factor(mp$status, levels = all_statuses)
      DT::datatable(mp,
        colnames = c("Original ID","Cleaned ID","Symbol","New ID","Status"),
        options  = list(scrollX=TRUE, pageLength=15, dom="ftipr"),
        class    = "compact stripe hover",
        rownames = FALSE,
        filter   = "top") |>
        DT::formatStyle("status",
          color = DT::styleEqual(
            all_statuses,
            c("#15803d", "#b91c1c", "#92400e", "#6b7280")
          ),
          fontWeight = "600")
    })

    # =========================================================================
    # MODAL PLOT - plotly if available, otherwise static
    # =========================================================================

    output$modal_plot_render <- (
      if (requireNamespace("plotly", quietly = TRUE))
        plotly::renderPlotly
      else
        shiny::renderPlot
    )({
      shiny::req(modal_plot())
      if (requireNamespace("plotly", quietly = TRUE)) {
        p <- modal_plot()
        # If plot has a 'text' aesthetic use it as the tooltip (sample name hover)
        has_text <- tryCatch("text" %in% names(p$mapping), error = function(e) FALSE)
        if (has_text) plotly::ggplotly(p, tooltip = c("text", "x", "y"))
        else plotly::ggplotly(p)
      } else {
        modal_plot()
      }
    })

    output$dl_plot_png <- shiny::downloadHandler(
      filename = function() paste0("markeR_plot_", Sys.Date(), ".png"),
      content  = function(file) {
        shiny::req(modal_plot())
        ggplot2::ggsave(file, modal_plot(), width = 12, height = 8, dpi = 150, bg = "white")
      }
    )

    # ── Expand button handlers ─────────────────────────────────────────────────
    shiny::observeEvent(input$expand_s1a, {
      shiny::req(rv$s1_plot_a)
      modal_plot(rv$s1_plot_a$plot)
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })
    shiny::observeEvent(input$expand_s1b, {
      shiny::req(rv$s1_plot_b)
      modal_plot(rv$s1_plot_b$plot)
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })
    shiny::observeEvent(input$expand_s2_all, {
      shiny::req(rv$data_s1)
      d <- rv$data_s1; thr <- shiny::isolate(rv$s2_thr)
      df <- data.frame(log_mean = log2(rowMeans(as.matrix(d)) + 1))
      modal_plot(ggplot2::ggplot(df, ggplot2::aes(x = log_mean)) +
        ggplot2::geom_density(fill = "#EBB43E", alpha = 0.45, colour = "#b8860b", linewidth = 0.8) +
        ggplot2::geom_vline(xintercept = thr, linetype = "dashed", colour = "firebrick", linewidth = 1.2) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), panel.border=ggplot2::element_blank(),
                       axis.line=ggplot2::element_line(colour="black")) +
        ggplot2::xlab("Mean log₂(count + 1)") + ggplot2::ylab("Density") +
        ggplot2::ggtitle("All genes - expression density"))
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })
    shiny::observeEvent(input$expand_s2_after, {
      shiny::req(rv$data_s1)
      d <- rv$data_s1; thr <- shiny::isolate(rv$s2_thr)
      lm <- log2(rowMeans(as.matrix(d)) + 1); kv <- lm[lm >= thr]
      if (length(kv) == 0) return(NULL)
      modal_plot(ggplot2::ggplot(data.frame(log_mean=kv), ggplot2::aes(x=log_mean)) +
        ggplot2::geom_density(fill="#28a745", alpha=0.4, colour="#155724", linewidth=0.8) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor=ggplot2::element_blank(), panel.border=ggplot2::element_blank(),
                       axis.line=ggplot2::element_line(colour="black")) +
        ggplot2::xlab("Mean log₂(count + 1)") + ggplot2::ylab("Density") +
        ggplot2::ggtitle(paste0("Retained genes (thr ≥ ", round(thr,2), ")")))
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })
    shiny::observeEvent(input$expand_s3_before, {
      shiny::req(rv$s3_norm_done, rv$data_s2)
      log_d      <- log2(as.matrix(rv$data_s2) + 1); all_s <- colnames(log_d)
      col_by     <- shiny::isolate(input$s3_color_by)
      hide_xlabs <- shiny::isolate(isTRUE(input$s3_hide_xlabels))
      samps      <- shiny::isolate(.s3_samps())
      modal_plot(.pp_make_boxplot(log_d, all_s, NULL, NULL, col_by, rv$meta_s,
                                   "Raw counts (before normalisation)",
                                   samps_plot = samps,
                                   hide_xlabels = hide_xlabs))
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })
    shiny::observeEvent(input$expand_s3_after, {
      shiny::req(rv$s3_norm_done, rv$norm_linear)
      log_d      <- log2(as.matrix(rv$norm_linear) + 1); all_s <- colnames(log_d)
      col_by     <- shiny::isolate(input$s3_color_by)
      hide_xlabs <- shiny::isolate(isTRUE(input$s3_hide_xlabels))
      samps      <- shiny::isolate(.s3_samps())
      method     <- rv$s3_method_used %||% "TMM"
      y_lbl      <- if (grepl("DESeq2", method, ignore.case = TRUE)) "VST"
                    else "log₂(CPM + 1)"
      modal_plot(.pp_make_boxplot(log_d, all_s, NULL, NULL, col_by, rv$meta_s,
                                   "After normalisation",
                                   samps_plot = samps,
                                   hide_xlabels = hide_xlabs,
                                   y_label = y_lbl))
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })
    shiny::observeEvent(input$expand_pca_before, {
      shiny::req(rv$s4_pca_df_before, rv$s4_pca_ev_before)
      plt <- .pp_pca_from_df(rv$s4_pca_df_before, rv$s4_pca_ev_before, rv$meta_s,
                              shiny::isolate(input$s4_pca_color),
                              if (rv$s4_bc_run) "PCA - before batch correction" else "PCA",
                              shiny::isolate(input$s4_pc_x %||% "PC1"),
                              shiny::isolate(input$s4_pc_y %||% "PC2"))
      modal_plot(plt)
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })
    shiny::observeEvent(input$expand_pca_after, {
      shiny::req(rv$s4_pca_df_after, rv$s4_pca_ev_after)
      plt <- .pp_pca_from_df(rv$s4_pca_df_after, rv$s4_pca_ev_after, rv$meta_s,
                              shiny::isolate(input$s4_pca_after_color %||% input$s4_pca_color),
                              "PCA - after batch correction",
                              shiny::isolate(input$s4_pc_after_x %||% "PC1"),
                              shiny::isolate(input$s4_pc_after_y %||% "PC2"))
      modal_plot(plt)
      shiny::showModal(.pp_plot_modal(ns("modal_plot_render"), ns("dl_plot_png")))
    })

    # =========================================================================
    # SIDEBAR - progress tracker with connecting line
    # =========================================================================

    step_info <- list(
      list(n=1L, label="Sample QC",
        tip="Check sample quality before filtering. Poor-quality or outlier samples caught here prevent downstream artefacts."),
      list(n=2L, label="Gene Filtering",
        tip="Lowly expressed genes add noise without signal. Removing them improves statistical power and speeds up analysis."),
      list(n=3L, label="Normalisation",
        tip="Samples differ in sequencing depth. Normalisation makes them comparable so biological differences aren't confounded by technical ones."),
      list(n=4L, label="PCA & Batch",
        tip="PCA reveals major sources of variation. If technical effects dominate, batch correction refocuses the data on biology.")
    )

    output$sidebar_progress <- shiny::renderUI({
      step <- rv$step
      n_steps <- length(step_info)
      shiny::div(
        lapply(seq_along(step_info), function(i) {
          x      <- step_info[[i]]
          is_last <- i == n_steps
          is_done <- x$n < step
          is_act  <- x$n == step
          if (is_done)       { dot_cls <- "pp-num pp-num-done";   txt <- "✓" }
          else if (is_act)   { dot_cls <- "pp-num pp-num-active"; txt <- as.character(x$n) }
          else               { dot_cls <- "pp-num pp-num-locked";  txt <- as.character(x$n) }
          shiny::div(class = "pp-progress-item",
            shiny::div(class = "pp-progress-dot-col",
              shiny::tags$span(class = dot_cls,
                style = "width:20px;height:20px;font-size:0.72em;", txt),
              if (!is_last)
                shiny::div(class = paste0("pp-progress-connector", if (is_done) " done" else ""))
            ),
            shiny::div(class = "pp-progress-label",
              shiny::div(style = "display:flex;align-items:center;justify-content:space-between;",
                shiny::tags$span(x$label, style = paste0("font-size:0.82em;",
                  if (is_act) "font-weight:700;" else "")),
                bslib::tooltip(
                  shiny::tags$span(style = "cursor:help;color:#9ca3af;font-size:0.8em;", "ⓘ"),
                  x$tip, placement = "right"
                )
              )
            )
          )
        })
      )
    })

    output$sidebar_data_info <- shiny::renderUI({
      e_curr <- if (!is.null(rv$data_s2)) rv$data_s2
                else if (!is.null(rv$data_s1)) rv$data_s1
                else get_expr()
      if (is.null(e_curr)) return(shiny::helpText("No data loaded.", style="font-size:0.82em;"))
      cg <- nrow(e_curr); cs <- ncol(e_curr)
      og <- rv$orig_n_genes; os <- rv$orig_n_samples
      filtered <- !is.null(og) && (cg != og || cs != os)
      shiny::tags$div(style = "font-size:0.82em;color:#555;",
        shiny::tags$b("Current data:"), shiny::tags$br(),
        paste0(format(cg, big.mark=","), " genes"), shiny::tags$br(),
        paste0(format(cs, big.mark=","), " samples"),
        if (filtered) shiny::tags$div(style="color:#9ca3af;font-size:0.82em;margin-top:3px;",
          paste0("(original: ", format(og,big.mark=","), "×", format(os,big.mark=","), ")")))
    })

    output$sidebar_downloads <- shiny::renderUI({
      if (!rv$downloads_available) return(NULL)
      # Determine analysis display settings for boxplot sample choice
      rid  <- rv$report_interactive_data
      s3d  <- if (!is.null(rid$s3_display_mode)) rid$s3_display_mode else "all"
      s3n  <- if (!is.null(rid$s3_n_box)) as.integer(rid$s3_n_box) else 20L
      has_subsamp <- !identical(s3d, "all")
      # Build the sample-choice radio (only shown when normalisation was done)
      show_report <- !isTRUE(rv$pp_summary$norm_method == "None (data used as-is)")
      report_ui <- if (show_report) {
        subsamp_label <- if (s3d == "random")
          paste0("Random ", s3n, " (as in analysis)")
        else if (s3d == "extreme")
          paste0("Extreme ", s3n, " (as in analysis)")
        else NULL
        sample_ctrl <- if (has_subsamp)
          shiny::radioButtons(
            ns("report_sample_choice"), label = NULL,
            choices  = c("All samples" = "all",
                         setNames(s3d, subsamp_label)),
            # Preserve the user's current selection on re-render; default to s3d
            selected = shiny::isolate(input$report_sample_choice) %||% s3d,
            inline   = FALSE
          )
        else
          shiny::tags$p(
            style = "font-size:0.75em;color:#6b7280;margin:2px 0 4px;",
            "Showing all samples (no subsampling used in analysis)."
          )
        shiny::tagList(
          shiny::tags$div(
            style = "margin-top:8px;margin-bottom:2px;",
            shiny::tags$label(
              style = "font-size:0.75em;font-weight:600;color:#6b7280;display:block;margin-bottom:3px;",
              "Samples in report plots:"
            ),
            sample_ctrl
          ),
          shiny::downloadButton(ns("dl_report"), "Summary report (.html)",
                                 class="btn-outline-secondary btn-sm pp-dl-btn")
        )
      } else NULL

      shiny::tagList(
        shiny::hr(style = "margin:8px 0;"),
        shiny::tags$div(style="font-size:0.8em;font-weight:600;color:#555;margin-bottom:2px;",
          "\U0001f4be Download results"),
        shiny::tags$p(style="font-size:0.75em;color:#9ca3af;margin:0 0 4px 0;",
          "Optional - data is already available to all analyses."),
        shiny::downloadButton(ns("dl_expr"), "Expression data (.csv)",
                               class="btn-outline-primary btn-sm pp-dl-btn"),
        shiny::downloadButton(ns("dl_meta"), "Metadata (.csv)",
                               class="btn-outline-secondary btn-sm pp-dl-btn"),
        report_ui,
        if (!is.null(get_log))
          shiny::downloadButton(ns("dl_log"), "Processing log (.txt)",
                                 class="btn-outline-secondary btn-sm pp-dl-btn")
      )
    })

    output$dl_expr <- shiny::downloadHandler(
      filename = function() paste0("markeR_expr_preprocessed_", Sys.Date(), ".csv"),
      content  = function(file) { shiny::req(rv$final_data); utils::write.csv(rv$final_data, file) }
    )
    output$dl_meta <- shiny::downloadHandler(
      filename = function() paste0("markeR_metadata_preprocessed_", Sys.Date(), ".csv"),
      content  = function(file) { shiny::req(rv$meta_s); utils::write.csv(rv$meta_s, file, row.names=FALSE) }
    )
    output$dl_log <- shiny::downloadHandler(
      filename = function() paste0("markeR_log_", Sys.Date(), ".txt"),
      content  = function(file) { writeLines(if (!is.null(get_log)) get_log() else character(0), file) }
    )

    output$dl_report <- shiny::downloadHandler(
      filename = function() paste0("markeR_preprocessing_report_", format(Sys.time(), "%Y-%m-%d_%H%M"), ".html"),
      content  = function(file) {
        shiny::req(rv$pp_summary, rv$final_data)
        s  <- rv$pp_summary
        fd <- rv$final_data
        mt <- rv$meta_s

        # ── Sample subsetting for plots ───────────────────────────────────────
        # Mirror the display mode chosen during preprocessing (step 3).
        samp_choice  <- isolate(input$report_sample_choice) %||% "all"
        rid_rep      <- rv$report_interactive_data
        s3d          <- if (!is.null(rid_rep$s3_display_mode)) rid_rep$s3_display_mode else "all"
        s3n_raw      <- if (!is.null(rid_rep$s3_n_box)) as.integer(rid_rep$s3_n_box) else 20L
        # Effective display mode: "all" if user chose "All samples", else mirror step 3
        display_mode <- if (identical(samp_choice, "all")) "all" else s3d

        # All available sample names across both matrices (consistent set for all plots)
        # Use data_s2 column order as canonical — matches .s3_samps() in the app
        # so random/extreme sample selection is identical between UI and report.
        all_samp_names <- if (!is.null(rv$data_s2)) colnames(rv$data_s2)
                          else if (!is.null(rv$norm_linear)) colnames(rv$norm_linear)
                          else character(0)
        # Cap s3n to at most floor(n_available/2) to avoid over-requesting samples
        s3n <- min(s3n_raw, max(2L, floor(length(all_samp_names) / 2L)))

        # Pick a CONSISTENT set using the same logic as .pp_make_boxplot
        plot_samp_names <- if (display_mode == "random") {
          n_pick <- min(s3n, length(all_samp_names))
          set.seed(42L + s3n)
          sample(all_samp_names, n_pick, replace = FALSE)
        } else if (display_mode == "extreme" &&
                   (!is.null(rv$norm_linear) || !is.null(rv$data_s2))) {
          # Use post-normalisation data for ranking (matches step 3 display)
          ref_mat <- if (!is.null(rv$norm_linear)) rv$norm_linear else rv$data_s2
          ref_d   <- log2(as.matrix(ref_mat) + 1)
          avail   <- intersect(all_samp_names, colnames(ref_d))
          meds   <- apply(ref_d[, avail, drop = FALSE], 2L, median)
          ord    <- order(meds)
          nb     <- floor(s3n / 2L); nt <- s3n - nb
          c(avail[ord[seq_len(nb)]],
            avail[ord[seq(length(avail) - nt + 1L, length(avail))]])
        } else {
          NULL  # "all" — no subsetting
        }

        # Helper to subset a matrix to the chosen sample names (consistent set)
        .subs <- function(m) {
          if (is.null(m) || is.null(plot_samp_names)) return(m)
          keep <- intersect(plot_samp_names, colnames(m))
          if (length(keep) == 0L) return(m)
          m[, keep, drop = FALSE]
        }
        samp_note <- if (!is.null(plot_samp_names)) {
          mode_label <- if (display_mode == "random") "random selection"
                        else "extreme values (N/2 lowest + N/2 highest median)"
          paste0("<p style='font-size:0.82em;color:#6b7280;font-style:italic;margin:4px 0 8px;'>",
                 "ℹ️ Plots show ", length(plot_samp_names), " samples (",
                 mode_label, ", same as analysis).</p>")
        } else ""

        # ── Helper: render a ggplot to an inline base64 img tag ─────────────
        plot_to_img <- function(p, w = 900, h = 420, res = 130) {
          if (is.null(p)) return("<p style='color:#9ca3af;font-size:0.85em;'>(plot not available)</p>")
          tmp <- tempfile(fileext = ".png")
          on.exit(tryCatch(unlink(tmp), error = function(e) NULL), add = TRUE)
          tryCatch({
            grDevices::png(tmp, width = w, height = h, res = res, bg = "white")
            print(p)
            grDevices::dev.off()
            uri <- base64enc::dataURI(file = tmp, mime = "image/png")
            paste0("<img src='", uri, "' style='max-width:100%;border:1px solid #e5e7eb;border-radius:6px;margin-top:8px;'>")
          }, error = function(e) {
            tryCatch(grDevices::dev.off(), error = function(e2) NULL)
            paste0("<p style='color:#ef4444;font-size:0.85em;'>Plot error: ", e$message, "</p>")
          })
        }

        # ── Helper: build a DataTable HTML snippet ────────────────────────────
        make_dt <- function(df, id, max_cols = 30, caption = "") {
          if (is.null(df) || nrow(df) == 0) return("<p>No data.</p>")
          cols <- utils::head(colnames(df), max_cols)
          has_rn <- !is.null(rownames(df)) && !identical(rownames(df), as.character(seq_len(nrow(df))))
          thead_cols <- if (has_rn) c("Gene", cols) else cols
          rows <- sapply(seq_len(nrow(df)), function(i) {
            cells <- if (has_rn)
              c(rownames(df)[i],
                sapply(cols, function(cc) {
                  v <- df[i, cc]
                  if (is.numeric(v)) round(v, 4) else as.character(v)
                }))
            else
              sapply(cols, function(cc) {
                v <- df[i, cc]
                if (is.numeric(v)) round(v, 4) else as.character(v)
              })
            paste0("<tr>", paste0("<td>", cells, "</td>", collapse=""), "</tr>")
          })
          trunc_note <- if (ncol(df) > max_cols)
            paste0("<p style='font-size:0.82em;color:#6b7280;'>Showing first ", max_cols,
                   " of ", ncol(df), " columns.</p>") else ""
          paste0(
            trunc_note,
            "<div style='overflow-x:auto;'>",
            "<table id='", id, "' class='display compact' style='width:100%;font-size:0.82em;'>",
            "<caption style='caption-side:top;font-weight:600;text-align:left;padding-bottom:4px;'>", caption, "</caption>",
            "<thead><tr>", paste0("<th>", thead_cols, "</th>", collapse=""), "</tr></thead>",
            "<tbody>", paste(rows, collapse=""), "</tbody>",
            "</table></div>"
          )
        }

        # ── Step plots ────────────────────────────────────────────────────────
        # Step 1: interactive Plotly (falls back to static PNG)
        s1a_div <- .report_ggplotly_div(
          if (!is.null(rv$s1_plot_a)) rv$s1_plot_a$plot else NULL, "s1a")
        s1b_div <- .report_ggplotly_div(
          if (!is.null(rv$s1_plot_b)) rv$s1_plot_b$plot else NULL, "s1b")

        # Step 2: before/after density plots (static PNG)
        s2_before_img <- if (!is.null(rv$report_s2_before))
          plot_to_img(rv$report_s2_before)
        else "<p style='color:#9ca3af;font-size:0.85em;'>(not available)</p>"
        s2_after_img <- if (!is.null(rv$report_s2_after))
          plot_to_img(rv$report_s2_after)
        else "<p style='color:#9ca3af;font-size:0.85em;'>(no retained genes)</p>"

        # Step 3: interactive boxplots with colour-by dropdown
        rid   <- rv$report_interactive_data
        mcols <- if (!is.null(rid)) rid$meta_cols else setdiff(colnames(rv$meta_s), "SampleID")
        s3c   <- if (!is.null(rid)) rid$s3_color_by else NULL
        # Subsample expression matrices for boxplots using a CONSISTENT sample set
        d_s2_plot   <- .subs(rv$data_s2)
        d_norm_plot <- .subs(rv$norm_linear)
        # Subset metadata to the SAME sample set (use whichever matrix is available)
        .subs_meta <- function(m, d) {
          if (is.null(m) || is.null(d)) return(m)
          keep <- intersect(rownames(m), colnames(d))
          if (length(keep) == 0L) return(m)
          m[keep, , drop = FALSE]
        }
        # Determine shared sample names across both matrices so metadata aligns both
        shared_samp <- if (!is.null(d_s2_plot) && !is.null(d_norm_plot))
          intersect(colnames(d_s2_plot), colnames(d_norm_plot))
        else if (!is.null(d_s2_plot)) colnames(d_s2_plot)
        else if (!is.null(d_norm_plot)) colnames(d_norm_plot)
        else character(0)
        mt_plot <- if (length(shared_samp) > 0L && !is.null(rv$meta_s)) {
          keep <- intersect(rownames(rv$meta_s), shared_samp)
          if (length(keep) > 0L) rv$meta_s[keep, , drop = FALSE] else rv$meta_s
        } else rv$meta_s
        norm_method_rep <- s$norm_method %||% "TMM"
        after_ylabel <- if (grepl("DESeq2", norm_method_rep, ignore.case = TRUE)) "VST"
                        else "log₂(CPM + 1)"
        s3_before_div <- tryCatch(
          if (!is.null(d_s2_plot))
            .report_box_div(.pp_box_stats_report(log2(as.matrix(d_s2_plot) + 1)),
                            mt_plot, mcols, s3c, "s3bef", "Before normalisation",
                            y_label = "log₂(count + 1)")
          else plot_to_img(rv$report_s3_before, h=420),
          error = function(e) plot_to_img(rv$report_s3_before, h=420)
        )
        s3_after_div  <- tryCatch(
          if (!is.null(d_norm_plot))
            .report_box_div(.pp_box_stats_report(log2(as.matrix(d_norm_plot) + 1)),
                            mt_plot, mcols, s3c, "s3aft", "After normalisation",
                            y_label = after_ylabel)
          else plot_to_img(rv$report_s3_after, h=420),
          error = function(e) plot_to_img(rv$report_s3_after, h=420)
        )

        # Step 4: interactive PCA with colour-by dropdown
        pc_x  <- if (!is.null(rid)) rid$pc_x else "PC1"
        pc_y  <- if (!is.null(rid)) rid$pc_y else "PC2"
        s4c   <- if (!is.null(rid)) rid$s4_pca_color else NULL
        s4_before_div <- tryCatch(
          .report_pca_div(rv$s4_pca_df_before, rv$s4_pca_ev_before,
                          rv$meta_s, mcols, s4c, pc_x, pc_y, "pca_bef"),
          error = function(e) plot_to_img(rv$report_s4_before)
        )
        s4_after_div  <- if (!is.null(rv$s4_pca_df_after))
          tryCatch(
            .report_pca_div(rv$s4_pca_df_after, rv$s4_pca_ev_after,
                            rv$meta_s, mcols, s4c, pc_x, pc_y, "pca_aft"),
            error = function(e) plot_to_img(rv$report_s4_after)
          )
        else NULL

        # ── Summary fields ────────────────────────────────────────────────────
        bc_txt <- if (isTRUE(s$bc_applied)) {
          bc_p  <- if (!is.null(s$bc_vars)   && length(s$bc_vars)   > 0L)
                     paste0("Batch removed: [", paste(s$bc_vars,   collapse=", "), "]") else NULL
          eff_p <- if (!is.null(s$bc_effect) && length(s$bc_effect) > 0L)
                     paste0("Retained: [", paste(s$bc_effect, collapse=", "), "]") else NULL
          paste0("Applied - ", paste(Filter(Negate(is.null), list(bc_p, eff_p)), collapse=" | "))
        } else "Not applied"

        kv <- function(label, val)
          paste0("<div class='kv'><div class='label'>", label,
                 "</div><div class='val'>", val, "</div></div>")

        expr_tbl  <- make_dt(fd, "tbl_expr", caption =
          paste0("Final expression data (", format(nrow(fd),big.mark=","),
                 " genes x ", ncol(fd), " samples)"))
        meta_tbl  <- if (!is.null(mt) && nrow(mt)>0)
          make_dt(mt, "tbl_meta", caption="Metadata") else "<p>No metadata.</p>"

        # ── Gene ID conversion block ──────────────────────────────────────────
        geneid_block <- local({
          gmap <- rv$gene_id_map
          if (is.null(gmap)) {
            "<p style='color:#9ca3af;font-size:0.86em;'>Gene ID conversion was not applied.</p>"
          } else {
            n_total      <- nrow(gmap)
            n_mapped     <- sum(gmap$status == "mapped")
            n_notfnd     <- sum(gmap$status == "not found")
            n_kept_orig  <- sum(gmap$status == "kept (original ID)")
            n_removed    <- sum(gmap$status == "removed")
            id_type      <- .html_esc(rv$gene_id_type %||% "symbol")
            pct          <- if (n_total > 0) round(100 * n_mapped / n_total, 1) else 0
            # Build a compact mapping table (show original, new, status columns)
            tbl_df <- gmap[, c("original_id", "new_id", "symbol", "status"), drop = FALSE]
            colnames(tbl_df) <- c("Original ID", "Converted ID", "Gene Symbol", "Status")
            map_tbl <- make_dt(tbl_df, "tbl_geneid",
                               caption = paste0("Gene ID mapping (", n_total, " genes)"))
            paste0(
              "<div class='kv-grid'>",
              kv("Conversion type", id_type),
              kv("Total genes",     format(n_total, big.mark = ",")),
              kv("Mapped to symbol", paste0(format(n_mapped, big.mark = ","), " (", pct, "%)")),
              if (n_kept_orig > 0)
                kv("Kept with original ID",
                   paste0(format(n_kept_orig, big.mark = ","),
                          " — kept in data using original identifiers")) else "",
              if (n_notfnd > 0)
                kv("Not found (unresolved)",
                   paste0(format(n_notfnd, big.mark = ","),
                          " — no symbol found; kept in data with cleaned ID")) else "",
              if (n_removed > 0)
                kv("Removed from data", format(n_removed, big.mark = ",")) else
                kv("Removed from data", "0 — all genes retained"),
              "</div>",
              map_tbl
            )
          }
        })

        # ── CSS + JS (DataTables from CDN) ────────────────────────────────────
        head_block <- paste0(
          "<meta charset='UTF-8'>",
          "<title>markeR Preprocessing Report</title>",
          "<script src='https://cdn.plot.ly/plotly-2.30.0.min.js'></script>",
          "<link rel='stylesheet' href='https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css'>",
          "<style>",
          "body{font-family:system-ui,sans-serif;max-width:1200px;margin:30px auto;color:#1f2937;padding:0 24px;}",
          "h1{color:#EBB43E;border-bottom:3px solid #EBB43E;padding-bottom:10px;margin-bottom:4px;}",
          "h2{color:#374151;margin-top:36px;font-size:1.05em;border-left:4px solid #EBB43E;padding-left:10px;}",
          "h3{color:#6b7280;font-size:0.92em;margin:16px 0 4px;}",
          ".kv-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(220px,1fr));gap:8px;margin:12px 0;}",
          ".kv{background:#f9fafb;border:1px solid #e5e7eb;border-radius:6px;padding:9px 13px;}",
          ".kv .label{font-size:0.78em;color:#6b7280;margin-bottom:2px;text-transform:uppercase;letter-spacing:.03em;}",
          ".kv .val{font-weight:600;font-size:0.95em;}",
          ".step-section{margin-top:28px;padding:16px;background:#fff;border:1px solid #e5e7eb;border-radius:8px;}",
          ".step-header{display:flex;align-items:center;gap:10px;margin-bottom:12px;}",
          ".step-num{display:inline-flex;align-items:center;justify-content:center;",
          "  width:26px;height:26px;border-radius:50%;background:#EBB43E;color:#fff;",
          "  font-weight:700;font-size:0.8em;flex-shrink:0;}",
          ".step-title{font-weight:700;font-size:1em;}",
          ".step-summary{font-size:0.86em;color:#374151;background:#f0fdf4;",
          "  border:1px solid #bbf7d0;border-radius:5px;padding:6px 10px;margin-bottom:8px;}",
          ".step-note{font-size:0.85em;color:#374151;background:#fffbeb;",
          "  border:1px solid #fde68a;border-radius:5px;padding:6px 10px;margin-bottom:10px;}",
          ".img-row{display:grid;grid-template-columns:1fr 1fr;gap:12px;}",
          "@media(max-width:700px){.img-row{grid-template-columns:1fr;}}",
          "table.dataTable thead th{background:#f3f4f6!important;color:#374151;}",
          "footer{margin-top:48px;color:#9ca3af;font-size:0.8em;border-top:1px solid #e5e7eb;padding-top:14px;}",
          "</style>"
        )
        script_block <- paste0(
          "<script src='https://code.jquery.com/jquery-3.7.0.min.js'></script>",
          "<script src='https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js'></script>",
          "<script>$(function(){",
          "$('#tbl_expr').DataTable({pageLength:15,scrollX:true,order:[]});",
          "$('#tbl_meta').DataTable({pageLength:10,scrollX:true,order:[]});",
          "if($('#tbl_geneid').length)$('#tbl_geneid').DataTable({pageLength:15,scrollX:true,order:[]});",
          "});</script>"
        )

        # ── Assemble HTML ─────────────────────────────────────────────────────
        html <- paste0(
          "<!DOCTYPE html><html><head>", head_block, "</head><body>",
          "<h1>markeR - Preprocessing Report</h1>",
          "<p style='color:#6b7280;font-size:0.86em;'>Generated: ",
          format(Sys.time(), "%Y-%m-%d %H:%M"), "</p>",

          # Gene ID conversion (before summary — it's the first step)
          "<h2>Gene ID Conversion</h2>",
          geneid_block,

          # Summary
          "<h2>Summary</h2>",
          "<div class='kv-grid'>",
          kv("Final dataset", paste0(format(s$genes_kept,big.mark=","), " genes x ", s$n_samples, " samples")),
          if (!is.null(s$genes_removed) && s$genes_removed > 0L)
            kv("Genes removed", paste0(format(s$genes_removed,big.mark=","), " (below threshold)")) else "",
          if (!is.null(s$samples_removed) && s$samples_removed > 0L)
            kv("Samples removed", as.character(s$samples_removed)) else "",
          kv("Normalisation", s$norm_method),
          kv("Gene filtering threshold", paste0("log2(count + 1) >= ", rv$s2_thr %||% "N/A")),
          kv("Batch correction", bc_txt),
          kv("Data scale", s$data_scale),
          "</div>",

          # Step 1
          "<div class='step-section'>",
          "<div class='step-header'><span class='step-num'>1</span><span class='step-title'>Sample Complexity QC</span></div>",
          if (nchar(rv$s1_summary) > 0)
            paste0("<div class='step-summary'>", rv$s1_summary, "</div>") else "",
          if (nchar(trimws(rv$note_s1 %||% "")) > 0L)
            paste0("<div class='step-note'>\U0001f4dd ", .html_esc(rv$note_s1), "</div>") else "",
          "<h3>Library complexity</h3>", s1a_div,
          "<h3>Read sampling</h3>", s1b_div,
          "</div>",

          # Step 2
          "<div class='step-section'>",
          "<div class='step-header'><span class='step-num'>2</span><span class='step-title'>Gene Filtering</span></div>",
          if (nchar(rv$s2_summary) > 0)
            paste0("<div class='step-summary'>", rv$s2_summary,
                   " | Threshold: log2(count + 1) >= ", rv$s2_thr %||% "N/A", "</div>") else "",
          if (nchar(trimws(rv$note_s2 %||% "")) > 0L)
            paste0("<div class='step-note'>\U0001f4dd ", .html_esc(rv$note_s2), "</div>") else "",
          "<div class='img-row'>",
          "<div><h3>Before filtering</h3>", s2_before_img, "</div>",
          "<div><h3>After filtering</h3>", s2_after_img, "</div>",
          "</div>",
          "</div>",

          # Step 3 — interactive boxplots
          "<div class='step-section'>",
          "<div class='step-header'><span class='step-num'>3</span><span class='step-title'>Normalisation (", s$norm_method, ")</span></div>",
          if (nchar(rv$s3_summary) > 0)
            paste0("<div class='step-summary'>", rv$s3_summary, "</div>") else "",
          if (nchar(trimws(rv$note_s3 %||% "")) > 0L)
            paste0("<div class='step-note'>\U0001f4dd ", .html_esc(rv$note_s3), "</div>") else "",
          samp_note,
          "<div class='img-row'>",
          "<div><h3>Before normalisation</h3>", s3_before_div, "</div>",
          "<div><h3>After normalisation</h3>",  s3_after_div,  "</div>",
          "</div>",
          "</div>",

          # Step 4 — interactive PCA
          "<div class='step-section'>",
          "<div class='step-header'><span class='step-num'>4</span><span class='step-title'>PCA & Batch Correction</span></div>",
          if (isTRUE(s$bc_applied))
            paste0("<div class='step-summary'>", bc_txt, "</div>") else "",
          if (nchar(trimws(rv$note_s4 %||% "")) > 0L)
            paste0("<div class='step-note'>\U0001f4dd ", .html_esc(rv$note_s4), "</div>") else "",
          samp_note,
          "<div class='img-row'>",
          "<div><h3>", if (isTRUE(s$bc_applied)) "Before batch correction" else "PCA", "</h3>", s4_before_div, "</div>",
          if (!is.null(s4_after_div))
            paste0("<div><h3>After batch correction</h3>", s4_after_div, "</div>") else "",
          "</div>",
          "</div>",

          # Data tables
          "<h2>Expression Data</h2>", expr_tbl,
          "<h2>Metadata</h2>", meta_tbl,

          "<footer>Generated by markeR - Disease Transcriptomics Lab</footer>",
          script_block,
          "</body></html>"
        )
        writeLines(html, file)
      }
    )

    output$pp_data_tables <- shiny::renderUI({
      shiny::req(rv$finalized > 0L, !is.null(rv$final_data))
      shiny::tagList(
        shiny::div(style="margin-top:18px;",
          shiny::tags$h6(style="color:#374151;font-weight:600;margin-bottom:6px;",
            "Processed expression data"),
          DT::DTOutput(ns("tbl_final_expr"), width="100%")
        ),
        if (!is.null(rv$meta_s))
          shiny::div(style="margin-top:18px;",
            shiny::tags$h6(style="color:#374151;font-weight:600;margin-bottom:6px;",
              "Metadata"),
            DT::DTOutput(ns("tbl_final_meta"), width="100%")
          )
      )
    })

    output$tbl_final_expr <- DT::renderDT({
      shiny::req(rv$finalized > 0L, rv$final_data)
      df <- rv$final_data
      df <- round(df, 3)
      DT::datatable(df,
        options = list(scrollX=TRUE, pageLength=10, dom="tip"),
        class   = "compact stripe hover",
        rownames= TRUE)
    })

    output$tbl_final_meta <- DT::renderDT({
      shiny::req(rv$finalized > 0L, rv$meta_s)
      DT::datatable(rv$meta_s,
        options = list(scrollX=TRUE, pageLength=5, dom="tip",
                       lengthMenu=list(c(5,10,25,-1),c("5","10","25","All"))),
        class   = "compact stripe hover",
        rownames= FALSE)
    })

    shiny::observeEvent(input$skip_preprocessing, {
      d <- rv$data_s0 %||% tryCatch(as.data.frame(get_expr()), error=function(e) NULL)
      if (is.null(d)) {
        shiny::showNotification("No data loaded - cannot skip.", type="warning"); return()
      }
      n_g <- nrow(d); n_s <- ncol(d)
      rv$final_data <- d
      rv$pp_summary <- list(
        genes_kept=n_g, genes_removed=0L,
        n_samples=n_s, samples_removed=0L,
        norm_method="None (data used as-is)",
        bc_applied=FALSE, bc_vars=NULL, bc_effect=NULL,
        data_scale="as imported (no normalisation applied)",
        bc_label="not applied"
      )
      rv$s1_summary <- "not applied"
      rv$s2_summary <- "not applied"
      rv$s3_summary <- "not applied"
      rv$downloads_available <- TRUE
      rv$step <- 5L; rv$active_step <- 5L
      rv$finalized <- rv$finalized + 1L      # triggers parent app to pick up final_data
      glog(paste0("Preprocessing skipped - using imported data as-is (",
        format(n_g,big.mark=","), " genes × ", n_s, " samples)."))
      shiny::showNotification(
        shiny::HTML(paste0("✅ Data passed through as-is - ",
          format(n_g,big.mark=","), " genes × ", n_s, " samples.")),
        type="default", duration=5)
    })

    shiny::observeEvent(input$reset_preprocessing, {
      # Use rv$data_s0_orig to restore state BEFORE any gene ID conversion.
      # rv$data_s0 may have converted rownames (symbols) after gene ID conversion;
      # rv$data_s0_orig always holds the data exactly as first imported.
      e <- rv$data_s0_orig %||% rv$data_s0
      m <- tryCatch(get_meta(), error = function(x) rv$meta_s)
      # Restore rv$data_s0 to pre-conversion state and reset gene ID tracking
      rv$data_s0           <- e
      rv$gene_id_map       <- NULL
      rv$gene_id_dismissed <- FALSE
      rv$gene_id_type      <- if (!is.null(e)) .pp_detect_gene_id_type(rownames(e)) else "symbol"
      rv$step <- 1L; rv$active_step <- 1L
      rv$data_s1 <- e; rv$data_s2 <- NULL
      rv$norm_linear <- NULL; rv$data_s3 <- NULL; rv$final_data <- NULL
      rv$meta_s <- m
      rv$orig_n_genes <- if (!is.null(e)) nrow(e) else rv$orig_n_genes
      rv$orig_n_samples <- if (!is.null(e)) ncol(e) else rv$orig_n_samples
      rv$s1_run <- FALSE; rv$s2_run <- FALSE; rv$s3_run <- FALSE; rv$s3_norm_done <- FALSE
      rv$s4_run <- FALSE; rv$s4_bc_run <- FALSE; rv$s2_thr <- 1.0; rv$s3_method_used <- NULL
      rv$s1_summary <- ""; rv$s2_summary <- ""; rv$s3_summary <- ""
      rv$s1_plot_a <- NULL; rv$s1_plot_b <- NULL
      rv$s4_pca_df_before <- NULL; rv$s4_pca_ev_before <- NULL
      rv$s4_pca_df_after  <- NULL; rv$s4_pca_ev_after  <- NULL
      rv$s4_bc_result <- NULL; rv$pp_summary <- NULL; rv$finalized <- 0L
      rv$downloads_available <- FALSE
      rv$report_s2_before <- NULL; rv$report_s2_after <- NULL
      rv$report_s3_before <- NULL; rv$report_s3_after <- NULL
      rv$report_s4_before <- NULL; rv$report_s4_after <- NULL
      rv$note_s1 <- ""; rv$note_s2 <- ""; rv$note_s3 <- ""; rv$note_s4 <- ""
      rv$report_interactive_data <- NULL
      .finalized_fp(NULL)
      # (silent reset — no notification needed)
    })

    # =========================================================================
    # STEP 1 - Sample Complexity QC
    # =========================================================================

    .s1_plot_samples_calc <- shiny::reactive({
      e <- rv$data_s0; shiny::req(e); all_s <- colnames(e)
      method <- input$s1_subset_method %||% "all"
      switch(method,
        all      = all_s,
        specific = { sel <- input$s1_specific_samples; if (length(sel)>0L) intersect(sel,all_s) else all_s },
        random   = sample(all_s, min(as.integer(input$s1_n_samp %||% 10L), length(all_s))),
        topbottom = {
          n <- as.integer(input$s1_n_samp %||% 10L); tots <- log2(colSums(e)+1); ord <- order(tots)
          nb <- floor(n/2L); nt <- n-nb
          unique(c(all_s[utils::head(ord,nb)], all_s[utils::tail(ord,nt)]))
        },
        all_s
      )
    })

    output$step1_card <- shiny::renderUI({
      e <- rv$data_s0
      if (is.null(e))
        return(bslib::card(class="pp-step-card",
          bslib::card_header(shiny::tags$div(style="display:flex;align-items:center;",
            shiny::tags$span(class="pp-num pp-num-active","1"),
            shiny::tags$strong("Step 1 - Sample Complexity QC (Optional)"))),
          shiny::div(style="padding:14px;", shiny::helpText("Load expression data first (Import tab)."))))
      if (rv$step > 1L && rv$active_step != 1L)
        return(.pp_done_card(1,"Sample Complexity QC",rv$s1_summary,ns("back_to_s1"),
          plot_ui=shiny::plotOutput(ns("plot_s1a"),height="360px"),
          note=rv$note_s1))
      controls <- shiny::fluidRow(
        shiny::column(6,
          shiny::tags$b("Plot type"),
          shiny::radioButtons(ns("s1_plot_type"),label=NULL,
            choices=c("Cumulative reads by gene index"="byexpr",
                      "Unique genes vs reads sampled (slow)"="sampling"),
            selected="byexpr"),
          shiny::checkboxInput(ns("s1_pct"),"Show reads as % of total",FALSE),
          shiny::conditionalPanel(sprintf("input['%s'] == 'sampling'",ns("s1_plot_type")),
            shiny::fluidRow(
              shiny::column(6, shiny::numericInput(ns("s1_n_genes"),"N genes (≤1000)",500,min=10,max=1000,step=50)),
              shiny::column(6, shiny::numericInput(ns("s1_n_iter"),"Iterations",5,min=1,max=20))),
            shiny::numericInput(ns("s1_step_reads"),"Read step size",500,min=10,max=10000,step=100),
            shiny::helpText("⚠️ Sampling is VERY slow as it is more accurate. Use fewer samples and genes to speed it up.",
                             style="color:#856404;font-size:0.8em;"))),
        shiny::column(6,
          shiny::tags$b("Samples to visualise"),
          shiny::div(class="pp-info-note",
            "⚠️ This selection is for the plot only - it does not filter your data."),
          shiny::radioButtons(ns("s1_subset_method"),label=NULL,
            choices=c("All samples"="all","Specific samples"="specific",
                      "Random N"="random","Top/Bottom N by lib-size"="topbottom"),
            selected="all"),
          shiny::conditionalPanel(sprintf("input['%s'] == 'specific'",ns("s1_subset_method")),
            shiny::uiOutput(ns("s1_specific_picker"))),
          shiny::conditionalPanel(sprintf("input['%s'] != 'all' && input['%s'] != 'specific'",
              ns("s1_subset_method"),ns("s1_subset_method")),
            {
              n_s1 <- if (!is.null(rv$data_s0)) ncol(rv$data_s0) else 20L
              def_n1 <- min(10L, max(2L, floor(n_s1 / 2L)))
              shiny::numericInput(ns("s1_n_samp"), "N samples",
                value = min(shiny::isolate(input$s1_n_samp) %||% def_n1, n_s1),
                min = 2L, max = n_s1, step = 1L)
            }))
      )
      bslib::card(class="pp-step-card",
        bslib::card_header(shiny::tags$div(style="display:flex;align-items:center;",
          shiny::tags$span(class="pp-num pp-num-active","1"),
          shiny::tags$strong("Step 1 - Sample Complexity QC"),
          shiny::tags$span(style="margin-left:8px;font-size:0.78em;color:#6c757d;font-style:italic;",
            "(optional)"))),
        shiny::div(style="padding:14px;", controls,
          shiny::hr(style="margin:8px 0;"),
          shiny::uiOutput(ns("s1_run_section"))))
    })

    output$s1_run_section <- shiny::renderUI({
      run <- rv$s1_run
      n_s <- if (!is.null(rv$data_s0)) ncol(rv$data_s0) else 0L
      if (!run)
        return(shiny::div(class="pp-run-btn",
          .pp_slow_warn(paste0("⏳ This processes each sample individually. ",
            "For large datasets (>50 samples), allow a minute.")),
          shiny::actionButton(ns("s1_run_btn"),"\U0001f4ca Show Complexity Plot",
            class="btn-primary btn-sm",style="margin-right:6px;"),
          shiny::actionButton(ns("skip_s1"),"Skip Step 1 →",class="btn-outline-secondary btn-sm")))
      n_sel <- tryCatch(length(shiny::isolate(.s1_plot_samples_calc())), error = function(e) n_s)
      shiny::tagList(
        shiny::div(style="margin-bottom:6px;",
          .pp_slow_warn(paste0("Plots show the selected ", n_sel, " sample(s). Click 'Reset step' to clear and start again.")),
          shiny::actionButton(ns("s1_reset_btn"),"\U0001f504 Reset step",
            class="btn-sm btn-outline-danger",style="margin-right:6px;"),
          shiny::actionButton(ns("skip_s1"),"Skip & proceed →",class="btn-sm btn-outline-secondary")),
        shiny::fluidRow(
          shiny::column(6,
            if (!is.null(rv$s1_plot_a))
              shiny::tagList(shiny::plotOutput(ns("plot_s1a"),height="360px"),
                             .pp_expand_btn(ns("expand_s1a"),
                               "⤢ Expand (hover to see sample name)"))
            else
              shiny::div(style="height:360px;display:flex;align-items:center;justify-content:center;color:#9ca3af;",
                "Click 'Show Complexity Plot' above.")),
          shiny::column(6,
            shiny::conditionalPanel(sprintf("input['%s'] == 'sampling'",ns("s1_plot_type")),
              shiny::tagList(
                shiny::actionButton(ns("run_sampling"),"▶ Run sampling plot",
                  class="btn-sm btn-outline-secondary",style="margin-bottom:6px;"),
                shiny::plotOutput(ns("plot_s1b"),height="360px"),
                if (!is.null(rv$s1_plot_b)) .pp_expand_btn(ns("expand_s1b"),
                  "⤢ Expand (hover to see sample name)"))),
            shiny::conditionalPanel(sprintf("input['%s'] != 'sampling'",ns("s1_plot_type")),
              shiny::div(style="display:flex;align-items:center;justify-content:center;height:360px;color:#9ca3af;font-size:0.88em;text-align:center;",
                shiny::HTML("Select <em>Unique genes vs reads sampled</em><br>above for a second plot."))))),
        .pp_remove_section_ui(ns("s1_remove_picker"),ns("s1_remove_btn")),
        shiny::div(style="margin-top:10px;",
          shiny::textAreaInput(ns("s1_notes"), "\U0001f4dd Step notes (optional — included in report):",
            value = shiny::isolate(rv$note_s1), rows = 2,
            placeholder = "e.g. Removed outlier samples based on low library complexity.",
            width = "100%")),
        shiny::div(class="pp-proceed-btn",
          shiny::actionButton(ns("proceed_s1"),"Proceed to Gene Filtering →",class="btn-warning btn-sm")))
    })

    output$s1_specific_picker <- shiny::renderUI({
      e <- rv$data_s0; shiny::req(e)
      shiny::selectInput(ns("s1_specific_samples"),label="Samples to show:",
        choices=colnames(e),selected=colnames(e),multiple=TRUE,selectize=TRUE,width="100%")
    })
    output$s1_remove_picker <- shiny::renderUI({
      e <- rv$data_s1; shiny::req(e, rv$s1_run)
      .pp_native_picker(ns("s1_remove_samples"), colnames(e))
    })
    output$plot_s1a <- shiny::renderPlot({
      shiny::req(rv$s1_run, !is.null(rv$s1_plot_a)); rv$s1_plot_a$plot
    })
    output$plot_s1b <- shiny::renderPlot({
      shiny::req(rv$s1_run); res <- rv$s1_plot_b
      shiny::validate(shiny::need(!is.null(res),"Click '▶ Run sampling plot' to generate.")); res$plot
    })

    shiny::observeEvent(input$s1_run_btn, {
      shiny::req(rv$data_s0)
      samps <- shiny::isolate(.s1_plot_samples_calc()); pct <- shiny::isolate(isTRUE(input$s1_pct))
      shiny::withProgress(message=paste0("\U0001f4ca Computing complexity - ",length(samps)," sample(s)…"),value=0,{
        res <- tryCatch(.pp_complexity_byexprgenes(rv$data_s0,samps,pct,
          progress_fn=function(frac,msg) shiny::incProgress(frac*0.9,msg)),
          error=function(e){shiny::showNotification(e$message,type="error");NULL})
        shiny::incProgress(0.1,"Rendering plot…"); rv$s1_plot_a <- res; rv$s1_run <- TRUE
      })
    })
    shiny::observeEvent(input$run_sampling, {
      shiny::req(rv$data_s0)
      samps <- shiny::isolate(.s1_plot_samples_calc())
      ng <- min(as.integer(shiny::isolate(input$s1_n_genes %||% 500L)),1000L)
      ni <- as.integer(shiny::isolate(input$s1_n_iter %||% 5L))
      sr <- as.integer(shiny::isolate(input$s1_step_reads %||% 500L))
      pct <- shiny::isolate(isTRUE(input$s1_pct))
      shiny::withProgress(message=paste0("⏳ Sampling reads - ",length(samps)," sample(s)…"),value=0,{
        res <- tryCatch(.pp_complexity_samplingreads(rv$data_s0,ng,ni,sr,samps,pct,
          progress_fn=function(frac,msg) shiny::incProgress(frac*0.9,msg)),
          error=function(e){shiny::showNotification(e$message,type="error",duration=10);NULL})
        shiny::incProgress(0.1,"Rendering…"); rv$s1_plot_b <- res
      })
    })
    shiny::observeEvent(input$s1_remove_btn, {
      to_rm <- input$s1_remove_samples; shiny::req(length(to_rm)>=1L,rv$data_s1)
      keep <- setdiff(colnames(rv$data_s1),to_rm)
      if (length(keep)<2L){shiny::showNotification("Cannot remove: fewer than 2 samples would remain.",type="error");return()}
      rv$data_s1 <- rv$data_s1[,keep,drop=FALSE]; .align_meta(keep)
      glog(paste0("Sample QC: removed ",length(to_rm)," sample(s): [",paste(to_rm,collapse=", "),"]."))
      shiny::showNotification(paste0(length(to_rm)," sample(s) removed. Re-run to refresh."),type="warning",duration=5)
    })
    # Reset step 1 - clear plots and allow user to re-run with different settings
    shiny::observeEvent(input$s1_reset_btn, {
      rv$s1_run    <- FALSE
      rv$s1_plot_a <- NULL
      rv$s1_plot_b <- NULL
    })
    shiny::observeEvent(input$skip_s1, { .go_to_s2() })
    shiny::observeEvent(input$proceed_s1, { .go_to_s2() })

    .go_to_s2 <- function() {
      shiny::req(rv$data_s1)
      rv$note_s1 <- shiny::isolate(input$s1_notes) %||% ""
      n_s <- ncol(rv$data_s1); removed <- (rv$orig_n_samples %||% n_s) - n_s
      rv$s1_summary <- paste0(n_s," samples",if(removed>0L) paste0(" (",removed," removed)") else "")
      glog(paste0("Sample Complexity QC: ",n_s," sample(s) passed",
        if(removed>0L) paste0(", ",removed," removed") else "","."))
      rv$step <- 2L; rv$active_step <- 2L; rv$s2_run <- TRUE; rv$s2_thr <- 1.0
    }
    shiny::observeEvent(input$back_to_s1, {
      rv$active_step <- 1L; rv$step <- 1L; rv$s1_run <- FALSE
      rv$data_s1 <- rv$data_s0; rv$meta_s <- get_meta()
      rv$data_s2 <- NULL; rv$norm_linear <- NULL; rv$data_s3 <- NULL
      rv$final_data <- NULL; rv$finalized <- 0L; rv$pp_summary <- NULL
      rv$downloads_available <- FALSE
      rv$report_s2_before <- NULL; rv$report_s2_after <- NULL
      rv$report_s3_before <- NULL; rv$report_s3_after <- NULL
      rv$report_s4_before <- NULL; rv$report_s4_after <- NULL
      rv$note_s2 <- ""; rv$note_s3 <- ""; rv$note_s4 <- ""
      rv$report_interactive_data <- NULL
      .finalized_fp(NULL)
    })

    # =========================================================================
    # STEP 2 - Gene Filtering
    # =========================================================================

    shiny::observeEvent(input$s2_plot_click, {
      shiny::req(input$s2_plot_click)
      val <- round(as.numeric(input$s2_plot_click$x),3)
      rv$s2_thr <- val; shiny::updateNumericInput(session,"s2_threshold_input",value=val)
    })
    shiny::observeEvent(input$s2_threshold_input, {
      val <- suppressWarnings(as.numeric(input$s2_threshold_input))
      if (!is.na(val)) rv$s2_thr <- val
    }, ignoreInit=TRUE)

    output$step2_card <- shiny::renderUI({
      if (rv$step < 2L)
        return(shiny::div(class="pp-locked",
          bslib::card(class="pp-step-card",
            bslib::card_header(shiny::tags$div(style="display:flex;align-items:center;",
              shiny::tags$span(class="pp-num pp-num-locked","2"),
              shiny::tags$strong("Step 2 - Filter Lowly Expressed Genes"))),
            shiny::div(style="padding:14px;color:#9ca3af;","Complete Step 1 (or skip it) to unlock."))))
      if (rv$step > 2L && rv$active_step != 2L)
        return(.pp_done_card(2,"Filter Lowly Expressed Genes",rv$s2_summary,ns("back_to_s2"),
          plot_ui=shiny::fluidRow(
            shiny::column(6, shiny::plotOutput(ns("plot_s2"),      height="360px")),
            shiny::column(6, shiny::plotOutput(ns("plot_s2_after"), height="360px"))),
          note=rv$note_s2))
      bslib::card(class="pp-step-card",
        bslib::card_header(shiny::tags$div(style="display:flex;align-items:center;",
          shiny::tags$span(class="pp-num pp-num-active","2"),
          shiny::tags$strong("Step 2 - Filter Lowly Expressed Genes"))),
        shiny::div(style="padding:14px;",
          .pp_slow_warn("⏳ Plots update as you adjust the threshold - allow a moment for large gene sets."),
          shiny::tags$p(class="pp-hint","\U0001f4a1 Click on the left plot to move the threshold, or type a value below."),
          shiny::fluidRow(
            shiny::column(4,
              shiny::numericInput(ns("s2_threshold_input"),"Threshold:",
                value=shiny::isolate(rv$s2_thr),min=0,step=0.1,width="150px"),
              shiny::uiOutput(ns("s2_gene_count")),
              shiny::div(style="margin-top:10px;",
                .pp_expand_btn(ns("expand_s2_all"),"⤢ Expand left"),
                .pp_expand_btn(ns("expand_s2_after"),"⤢ Expand right"))),
            shiny::column(8,
              shiny::fluidRow(
                shiny::column(6,
                  shiny::tags$b(style="font-size:0.85em;","All genes"),
                  shiny::plotOutput(ns("plot_s2"),click=ns("s2_plot_click"),height="250px")),
                shiny::column(6,
                  shiny::tags$b(style="font-size:0.85em;","Genes passing filter"),
                  shiny::plotOutput(ns("plot_s2_after"),height="250px"))))),
          shiny::div(style="margin-top:10px;",
            shiny::textAreaInput(ns("s2_notes"), "\U0001f4dd Step notes (optional — included in report):",
              value = shiny::isolate(rv$note_s2), rows = 2,
              placeholder = "e.g. Threshold chosen to remove genes with near-zero counts across all samples.",
              width = "100%")),
          shiny::div(class="pp-proceed-btn",
            shiny::actionButton(ns("proceed_s2"),
              "Apply filter & proceed to Normalisation →",class="btn-warning btn-sm"))))
    })

    output$s2_gene_count <- shiny::renderUI({
      shiny::req(rv$step>=2L,rv$data_s1)
      d <- rv$data_s1; thr <- rv$s2_thr; n_total <- nrow(d)
      n_keep <- sum(log2(rowMeans(as.matrix(d))+1) >= thr)
      shiny::tags$div(style="font-size:0.83em;color:#555;margin-top:8px;line-height:1.7;",
        shiny::HTML(paste0("<b>",format(n_keep,big.mark=","),"</b> genes kept<br>",
          "<b>",format(n_total-n_keep,big.mark=","),"</b> removed (",
          round(100*(n_total-n_keep)/n_total,1),"%)")))
    })
    output$plot_s2 <- shiny::renderPlot({
      shiny::req(rv$s2_run,rv$step>=2L,rv$data_s1)
      d <- rv$data_s1; thr <- rv$s2_thr
      df <- data.frame(log_mean = log2(rowMeans(as.matrix(d))+1))
      # Only show the 'click to set threshold' hint while step 2 is the active step
      step2_active <- shiny::isolate(rv$active_step) == 2L
      title_txt <- if (step2_active) "All genes - click to set threshold" else "All genes"
      ggplot2::ggplot(df,ggplot2::aes(x=log_mean)) +
        ggplot2::geom_density(fill="#EBB43E",alpha=0.45,colour="#b8860b",linewidth=0.8) +
        ggplot2::geom_vline(xintercept=as.numeric(thr),linetype="dashed",colour="firebrick",linewidth=1.2) +
        ggplot2::annotate("text",x=as.numeric(thr),y=Inf,vjust=1.8,hjust=-0.1,
          label=paste0("thr = ",round(thr,3)),colour="firebrick",size=3.5) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),panel.border=ggplot2::element_blank(),
          axis.line=ggplot2::element_line(colour="black")) +
        ggplot2::xlab("Mean log₂(count + 1)") + ggplot2::ylab("Density") +
        ggplot2::ggtitle(title_txt)
    })
    output$plot_s2_after <- shiny::renderPlot({
      shiny::req(rv$s2_run,rv$step>=2L,rv$data_s1)
      d <- rv$data_s1; thr <- rv$s2_thr
      kv <- log2(rowMeans(as.matrix(d))+1); kv <- kv[kv >= thr]
      if (length(kv)==0) return(ggplot2::ggplot() +
        ggplot2::annotate("text",x=0.5,y=0.5,label="No genes pass this threshold",size=4.5,colour="#9ca3af")+
        ggplot2::theme_void())
      ggplot2::ggplot(data.frame(log_mean=kv),ggplot2::aes(x=log_mean)) +
        ggplot2::geom_density(fill="#28a745",alpha=0.4,colour="#155724",linewidth=0.8) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),panel.border=ggplot2::element_blank(),
          axis.line=ggplot2::element_line(colour="black")) +
        ggplot2::xlab("Mean log₂(count + 1)") + ggplot2::ylab("Density") +
        ggplot2::ggtitle(paste0("Retained genes (thr ≥ ",round(thr,2),")"))
    })
    shiny::observeEvent(input$proceed_s2, {
      shiny::req(rv$step>=2L,rv$data_s1)
      rv$note_s2 <- shiny::isolate(input$s2_notes) %||% ""
      d <- rv$data_s1; thr <- as.numeric(rv$s2_thr)
      keep <- rownames(d)[log2(rowMeans(as.matrix(d))+1) >= thr]
      if (length(keep)<10L){shiny::showNotification("Too few genes retained - lower the threshold.",type="error");return()}
      rv$data_s2 <- d[keep,,drop=FALSE]
      n_rm <- nrow(d) - length(keep)
      rv$s2_summary <- paste0(format(length(keep),big.mark=",")," genes kept, ",
        format(n_rm,big.mark=",")," removed (thr=",thr,")")
      glog(paste0("Gene Filtering: threshold=",thr," - ",format(length(keep),big.mark=","),
        " genes retained, ",format(n_rm,big.mark=",")," removed."))
      rv$step <- 3L; rv$active_step <- 3L; rv$s3_run <- FALSE; rv$s3_norm_done <- FALSE
    })
    shiny::observeEvent(input$back_to_s2, {
      rv$active_step <- 2L; rv$step <- 2L; rv$s2_run <- TRUE
      rv$data_s2 <- NULL; rv$norm_linear <- NULL; rv$data_s3 <- NULL
      rv$s3_norm_done <- FALSE; rv$s3_run <- FALSE
      rv$final_data <- NULL; rv$finalized <- 0L; rv$pp_summary <- NULL
      rv$downloads_available <- FALSE
      rv$report_s3_before <- NULL; rv$report_s3_after <- NULL
      rv$report_s4_before <- NULL; rv$report_s4_after <- NULL
      .finalized_fp(NULL)
    })

    # =========================================================================
    # STEP 3 - Normalisation
    # Normalisation runs ONCE on button click (heavy).
    # Boxplots are reactive to display params - they auto-update without re-running TMM.
    # =========================================================================

    output$step3_card <- shiny::renderUI({
      if (rv$step < 3L)
        return(shiny::div(class="pp-locked",
          bslib::card(class="pp-step-card",
            bslib::card_header(shiny::tags$div(style="display:flex;align-items:center;",
              shiny::tags$span(class="pp-num pp-num-locked","3"),
              shiny::tags$strong("Step 3 - Normalisation"))),
            shiny::div(style="padding:14px;color:#9ca3af;","Complete Steps 1–2 to unlock."))))
      if (rv$step > 3L && rv$active_step != 3L)
        return(.pp_done_card(3,"Normalisation",rv$s3_summary,ns("back_to_s3"),
          plot_ui=shiny::fluidRow(
            shiny::column(6, shiny::plotOutput(ns("plot_s3_before"),height="360px")),
            shiny::column(6, shiny::plotOutput(ns("plot_s3_after"), height="360px"))),
          note=rv$note_s3))
      meta_cols <- setdiff(colnames(rv$meta_s),"SampleID")
      bslib::card(class="pp-step-card",
        bslib::card_header(shiny::tags$div(style="display:flex;align-items:center;",
          shiny::tags$span(class="pp-num pp-num-active","3"),
          shiny::tags$strong("Step 3 - Normalisation & Quality Distribution"),
          if (rv$s3_norm_done) shiny::tags$span(class="pp-norm-badge",
            paste0(rv$s3_method_used %||% "Norm", " done")))),
        shiny::div(style="padding:14px;",
          # Normalisation method (always visible)
          shiny::fluidRow(
            shiny::column(4,
              shiny::tags$b("Normalisation method"),
              shiny::radioButtons(ns("s3_method"),label=NULL,
                choices=c("TMM (edgeR)"="TMM","VST (DESeq2)"="DESeq2"),selected="TMM"),
              shiny::uiOutput(ns("s3_run_btns"))
            ),
            shiny::column(8,
              if (rv$s3_norm_done)
                shiny::tagList(
                  shiny::tags$p(style="font-size:0.82em;color:#6b7280;margin-bottom:4px;",
                    "Display options - changing these updates plots immediately without re-normalising:"),
                  shiny::fluidRow(
                    shiny::column(5,
                      shiny::tags$b(style="font-size:0.85em;","Samples in boxplot"),
                      bslib::tooltip(
                        shiny::tags$span(style="cursor:help;color:#9ca3af;font-size:0.8em;margin-left:4px;","ⓘ"),
                        "Normalisation runs on ALL samples - this selection is for display only.",
                        placement="top"),
                      shiny::radioButtons(ns("s3_sample_display"),label=NULL,
                        choiceValues=c("all","random","extreme"),
                        choiceNames=list("All","Random N",
                          shiny::tags$span("Extreme N",
                            bslib::tooltip(shiny::tags$span(style="cursor:help;color:#9ca3af;font-size:0.8em;margin-left:3px;","ⓘ"),
                              "N/2 highest + N/2 lowest median samples - useful to spot outliers.",placement="top"))),
                        selected = {
                          prev <- shiny::isolate(input$s3_sample_display)
                          if (!is.null(prev)) prev
                          else if (!is.null(rv$data_s2) && ncol(rv$data_s2) > 20L) "extreme"
                          else "all"
                        }),
                      shiny::conditionalPanel(sprintf("input['%s'] != 'all'",ns("s3_sample_display")),
                        {
                          n_s3 <- if (!is.null(rv$data_s2)) ncol(rv$data_s2) else 40L
                          def_n3 <- min(20L, max(2L, floor(n_s3 / 2L)))
                          shiny::numericInput(ns("s3_n_box"), "N",
                            value = min(shiny::isolate(input$s3_n_box) %||% def_n3, n_s3),
                            min = 2L, max = n_s3, step = 1L)
                        })),
                    shiny::column(7,
                      shiny::tags$b(style="font-size:0.85em;","Colour by"),
                      if (length(meta_cols)>0L)
                        shiny::selectInput(ns("s3_color_by"),label=NULL,choices=meta_cols,selected=meta_cols[1])
                      else shiny::helpText("No metadata columns."),
                      shiny::uiOutput(ns("s3_color_warn")),
                      shiny::checkboxInput(ns("s3_hide_xlabels"), "Hide x-axis labels", value = FALSE))
                  )
                )
              else
                shiny::div(style="color:#9ca3af;font-size:0.85em;padding:14px 0;",
                  "Run normalisation first to unlock display options.")
            )
          ),
          shiny::uiOutput(ns("s3_plots_section"))))
    })

    output$s3_run_btns <- shiny::renderUI({
      n_samp <- if (!is.null(rv$data_s2)) ncol(rv$data_s2) else "?"
      method <- shiny::isolate(input$s3_method %||% "TMM")
      if (!rv$s3_norm_done)
        shiny::div(class="pp-run-btn",
          .pp_slow_warn(paste0("Normalisation runs on all ", n_samp,
            " samples once - boxplots update automatically afterwards.")),
          shiny::actionButton(ns("s3_run_btn"),"\U0001f4ca Run Normalisation",class="btn-primary btn-sm"))
      else
        shiny::div(style="margin-top:8px;",
          shiny::actionButton(ns("s3_run_btn"),"\U0001f504 Re-run Normalisation",
            class="btn-sm btn-outline-secondary",style="margin-bottom:6px;"),
          shiny::tags$p(style="font-size:0.78em;color:#6b7280;margin:0;",
            paste0("Re-runs normalisation from scratch (use if method changed).")))
    })

    output$s3_color_warn <- shiny::renderUI({ .pp_color_warn(rv$meta_s, input$s3_color_by) })

    output$s3_plots_section <- shiny::renderUI({
      if (!rv$s3_norm_done) return(NULL)
      shiny::tagList(
        shiny::hr(style="margin:10px 0;"),
        shiny::fluidRow(
          shiny::column(6,
            shiny::tags$b(style="font-size:0.85em;","Before normalisation"),
            shiny::plotOutput(ns("plot_s3_before"),height="300px"),
            .pp_expand_btn(ns("expand_s3_before"))),
          shiny::column(6,
            shiny::tags$b(style="font-size:0.85em;","After normalisation"),
            shiny::plotOutput(ns("plot_s3_after"),height="300px"),
            .pp_expand_btn(ns("expand_s3_after")))
        ),
        .pp_remove_section_ui(ns("s3_remove_picker"),ns("s3_remove_btn")),
        shiny::div(style="margin-top:10px;",
          shiny::textAreaInput(ns("s3_notes"), "\U0001f4dd Step notes (optional — included in report):",
            value = shiny::isolate(rv$note_s3), rows = 2,
            placeholder = "e.g. TMM selected; VST requires integer counts.",
            width = "100%")),
        shiny::div(class="pp-proceed-btn",
          shiny::actionButton(ns("proceed_s3"),
            "Confirm normalisation & proceed to PCA →",class="btn-warning btn-sm"))
      )
    })

    # TMM runs ONLY here - boxplots are handled by reactive renderPlots below
    shiny::observeEvent(input$s3_run_btn, {
      shiny::req(rv$step >= 3L, rv$data_s2)
      d <- as.matrix(rv$data_s2)
      method <- shiny::isolate(input$s3_method %||% "TMM")
      shiny::withProgress(
        message = paste0("⚙️ Running ", method, " on all ", ncol(d), " samples × ",
                          format(nrow(d),big.mark=","), " genes…"),
        value = 0, {
          shiny::incProgress(0.2, paste0("Applying ", method, "…"))
          res <- if (method == "TMM") {
            dge <- edgeR::DGEList(counts=d); dge <- edgeR::calcNormFactors(dge)
            nl  <- edgeR::cpm(dge, log=FALSE)
            list(norm_linear=nl, log_norm=log2(nl+1), method="TMM")
          } else {
            # DESeq2 requires non-negative raw integer counts
            if (any(d < 0)) {
              msg <- "DESeq2/VST requires non-negative counts. Data contains negative values — it may be log-transformed or pre-normalised. Use TMM instead, or supply raw counts."
              shiny::showNotification(msg, type = "error", duration = NULL)
              glog(paste0("[DESeq2 ERROR] ", msg))
              return(NULL)
            }
            # Warn if data looks pre-normalised (many fractional values or very low max)
            frac_pct <- mean(abs(d - round(d)) > 0.01) * 100
            if (frac_pct > 30) {
              msg <- paste0("DESeq2 WARNING: ", round(frac_pct), "% of values are non-integer. DESeq2 expects raw integer read counts. Values rounded — results may be unreliable if data is pre-normalised (TPM, FPKM, CPM). Consider using TMM instead.")
              shiny::showNotification(shiny::HTML(paste0("<b>DESeq2 warning:</b> ", msg)),
                type = "warning", duration = 15)
              glog(msg)
            } else if (frac_pct > 0) {
              msg <- paste0("DESeq2: ", round(frac_pct, 1), "% of values were non-integer and have been rounded to the nearest integer.")
              shiny::showNotification(msg, type = "message", duration = 8)
              glog(msg)
            }
            int_d <- round(d); storage.mode(int_d) <- "integer"
            col_df <- if(!is.null(rv$meta_s)) rv$meta_s else data.frame(row.names=colnames(int_d))
            dds <- tryCatch(DESeq2::DESeqDataSetFromMatrix(int_d,colData=col_df,design=~1),
              error=function(e){shiny::showNotification(e$message,type="error");NULL})
            if (is.null(dds)) return(NULL)
            vst_obj <- tryCatch(DESeq2::vst(dds,blind=TRUE),
              error=function(e){shiny::showNotification(e$message,type="error");NULL})
            if (is.null(vst_obj)) return(NULL)
            lg <- SummarizedExperiment::assay(vst_obj)
            # VST output is in a log-like scale; store as-is for PCA/downstream.
            # norm_linear holds a pseudo-linear back-transform used for batch correction input.
            list(norm_linear=2^lg-1, log_norm=lg, method="DESeq2/VST")
          }
          if (is.null(res)) return(NULL)
          rv$norm_linear     <- res$norm_linear
          rv$data_s3         <- as.data.frame(res$log_norm)
          rv$s3_method_used  <- res$method   # store actual method for badge / report
          shiny::incProgress(0.8, "Done - boxplots updating…")
          rv$s3_norm_done <- TRUE; rv$s3_run <- TRUE
        })
    })

    # Shared sample selection - computed ONCE from display params, reused for both plots
    # so "random N" and "extreme N" show the SAME samples before and after normalisation.
    .s3_samps <- shiny::reactive({
      shiny::req(rv$s3_norm_done, rv$data_s2, rv$norm_linear)
      all_s   <- colnames(rv$data_s2)
      # Use same default threshold as the radioButtons selected= logic
      display <- input$s3_sample_display %||%
        if (length(all_s) > 20L) "extreme" else "all"
      n_box   <- as.integer(input$s3_n_box %||% min(20L, max(2L, floor(length(all_s) / 2L))))
      switch(display,
        random  = { set.seed(42L + n_box); sample(all_s, min(n_box, length(all_s))) },
        extreme = {
          # Rank by post-normalisation medians so "extreme" reflects the corrected data
          log_d <- log2(as.matrix(rv$norm_linear) + 1)
          meds  <- apply(log_d, 2, stats::median)
          ord   <- order(meds); nb <- floor(n_box / 2L); nt <- n_box - nb
          unique(c(all_s[utils::head(ord, nb)], all_s[utils::tail(ord, nt)]))
        },
        all_s
      )
    })

    # Boxplots are fully reactive - update automatically when display params change.
    # Both use the SAME sample set so comparisons are valid.
    output$plot_s3_before <- shiny::renderPlot({
      shiny::req(rv$s3_norm_done, rv$data_s2)
      samps       <- .s3_samps()
      log_d       <- log2(as.matrix(rv$data_s2) + 1)
      color_by    <- input$s3_color_by
      hide_xlabs  <- isTRUE(input$s3_hide_xlabels)
      .pp_make_boxplot(log_d, colnames(log_d),
                        input$s3_sample_display %||% "all",
                        as.integer(input$s3_n_box %||% 20L),
                        color_by, rv$meta_s, "Before normalisation",
                        samps_plot = samps, hide_xlabels = hide_xlabs,
                        y_label = "log₂(count + 1)")
    })
    output$plot_s3_after <- shiny::renderPlot({
      shiny::req(rv$s3_norm_done, rv$norm_linear)
      samps       <- .s3_samps()
      method      <- rv$s3_method_used %||% "TMM"
      y_lbl       <- if (grepl("DESeq2", method, ignore.case = TRUE)) "VST"
                     else "log₂(CPM + 1)"
      log_d       <- log2(as.matrix(rv$norm_linear) + 1)
      color_by    <- input$s3_color_by
      hide_xlabs  <- isTRUE(input$s3_hide_xlabels)
      .pp_make_boxplot(log_d, colnames(log_d),
                        input$s3_sample_display %||% "all",
                        as.integer(input$s3_n_box %||% 20L),
                        color_by, rv$meta_s, "After normalisation",
                        samps_plot = samps, hide_xlabels = hide_xlabs,
                        y_label = y_lbl)
    })

    output$s3_remove_picker <- shiny::renderUI({
      e <- rv$data_s2; shiny::req(e, rv$s3_run)
      .pp_native_picker(ns("s3_remove_samples"), colnames(e))
    })
    shiny::observeEvent(input$s3_remove_btn, {
      to_rm <- input$s3_remove_samples; shiny::req(length(to_rm)>=1L,rv$data_s2)
      keep <- setdiff(colnames(rv$data_s2),to_rm)
      if (length(keep)<2L){shiny::showNotification("Cannot remove: fewer than 2 samples would remain.",type="error");return()}
      rv$data_s2 <- rv$data_s2[,keep,drop=FALSE]
      if (!is.null(rv$data_s1)) rv$data_s1 <- rv$data_s1[,keep,drop=FALSE]
      .align_meta(keep)
      glog(paste0("Normalisation step: removed ",length(to_rm)," sample(s): [",paste(to_rm,collapse=", "),"]."))
      shiny::showNotification(paste0(length(to_rm)," sample(s) removed. Re-run to refresh."),type="warning",duration=5)
    })
    shiny::observeEvent(input$proceed_s3, {
      shiny::req(rv$step>=3L, rv$norm_linear, rv$data_s3)
      rv$note_s3 <- shiny::isolate(input$s3_notes) %||% ""
      # Use stored method (radioButtons may have reset to default "TMM" after card re-render)
      method <- rv$s3_method_used %||% shiny::isolate(input$s3_method %||% "TMM")
      rv$s3_summary <- paste0(method," - ",format(nrow(rv$data_s3),big.mark=","),
        " genes × ",ncol(rv$data_s3)," samples")
      glog(paste0("Normalisation: ",method," applied - ",format(nrow(rv$data_s3),big.mark=","),
        " genes × ",ncol(rv$data_s3)," samples. Final data will be in non-log (CPM) scale."))
      rv$step <- 4L; rv$active_step <- 4L; rv$s4_run <- FALSE; rv$s4_bc_run <- FALSE
    })
    shiny::observeEvent(input$back_to_s3, {
      rv$active_step <- 3L; rv$step <- 3L; rv$s3_run <- FALSE; rv$s3_norm_done <- FALSE
      rv$norm_linear <- NULL; rv$data_s3 <- NULL
      rv$final_data <- NULL; rv$finalized <- 0L; rv$pp_summary <- NULL
      rv$downloads_available <- FALSE
      rv$report_s4_before <- NULL; rv$report_s4_after <- NULL
      .finalized_fp(NULL)
    })

    # =========================================================================
    # STEP 4 - PCA & Optional Batch Correction
    # PCA is computed ONCE on button click; colour/PC changes only redraw.
    # Scree plot shows per-PC variance; user picks which PCs to display.
    # =========================================================================

    output$step4_card <- shiny::renderUI({
      if (rv$step < 4L)
        return(shiny::div(class="pp-locked",
          bslib::card(class="pp-step-card",
            bslib::card_header(shiny::tags$div(style="display:flex;align-items:center;",
              shiny::tags$span(class="pp-num pp-num-locked","4"),
              shiny::tags$strong("Step 4 - PCA & Batch Correction (Optional)"))),
            shiny::div(style="padding:14px;color:#9ca3af;","Complete Steps 1–3 to unlock."))))
      if (rv$step > 4L && rv$active_step != 4L)
        return(.pp_done_card(4,"PCA & Batch Correction",
          if (!is.null(rv$pp_summary)) rv$pp_summary$bc_label else "done",
          ns("back_to_s4"),
          plot_ui=shiny::fluidRow(
            shiny::column(6, shiny::plotOutput(ns("pca_before"),height="380px")),
            shiny::column(6, if (rv$s4_bc_run) shiny::plotOutput(ns("pca_after"),height="380px"))),
          note=rv$note_s4))
      meta_cols <- setdiff(colnames(rv$meta_s),"SampleID")
      bslib::card(class="pp-step-card",
        bslib::card_header(shiny::tags$div(style="display:flex;align-items:center;",
          shiny::tags$span(class="pp-num pp-num-active","4"),
          shiny::tags$strong("Step 4 - PCA & Batch Correction (Optional)"))),
        shiny::div(style="padding:14px;",
          shiny::uiOutput(ns("s4_run_section"))))
    })

    output$s4_run_section <- shiny::renderUI({
      meta_cols <- setdiff(colnames(rv$meta_s), "SampleID")
      n_g <- if (!is.null(rv$data_s3)) nrow(rv$data_s3) else "?"
      n_s <- if (!is.null(rv$data_s3)) ncol(rv$data_s3) else "?"

      if (!rv$s4_run)
        return(shiny::div(class = "pp-run-btn",
          .pp_slow_warn(paste0("⏳ PCA decomposes the full normalised matrix (", n_g, " genes × ", n_s,
            " samples). This runs once — changing colour or PCs only redraws.")),
          shiny::div(style = "margin: 10px 0 8px 0;",
            shiny::tags$p(style = "font-size:0.82em;font-weight:600;margin-bottom:4px;", "PCA options:"),
            shiny::checkboxInput(ns("s4_pca_center"), "Center genes (subtract mean per gene — recommended)",
              value = shiny::isolate(input$s4_pca_center) %||% TRUE),
            shiny::checkboxInput(ns("s4_pca_scale"),
              "Scale genes (divide by SD — equalises variance; use with caution for RNA-seq)",
              value = shiny::isolate(input$s4_pca_scale) %||% FALSE)
          ),
          shiny::actionButton(ns("s4_run_btn"), "\U0001f4ca Compute PCA", class = "btn-primary btn-sm")))

      n_pcs_avail <- if (!is.null(rv$s4_pca_df_before)) ncol(rv$s4_pca_df_before) else 2L
      pc_choices  <- paste0("PC", 1:n_pcs_avail)
      # bc_done depends only on rv state - not checkbox - so layout survives checkbox toggles
      bc_done     <- rv$s4_bc_run && !is.null(rv$s4_bc_result)
      # ── Finalize ─────────────────────────────────────────────────────────────
      finalize_ui <- if (bc_done) {
        shiny::tagList(
          shiny::tags$p(style = "font-size:0.84em;color:#555;margin-top:10px;font-weight:600;",
            "Batch correction complete. Choose which data to use downstream:"),
          shiny::div(style = "display:flex;gap:8px;flex-wrap:wrap;margin-top:4px;",
            shiny::actionButton(ns("finalize_bc"), "✓ Use batch-corrected data",
              class = "btn-warning btn-sm"),
            shiny::actionButton(ns("finalize_norm"), "✓ Use normalised only (no batch correction)",
              class = "btn-outline-secondary btn-sm"))
        )
      } else {
        shiny::div(style = "display:flex;align-items:center;gap:12px;margin-top:12px;",
          shiny::actionButton(ns("finalize_btn"), "✓ Finalise & apply preprocessing",
            class = "btn-success btn-sm"),
          shiny::tags$span(style = "font-size:0.82em;color:#6c757d;",
            "Updates expression data used by all downstream analyses."))
      }

      # ── PCA options panel: rendered separately so checkbox changes react live ──
      pca_top_ui <- shiny::uiOutput(ns("s4_pca_top_ui"))

      # ── BC controls ───────────────────────────────────────────────────────────
      bc_controls_ui <- shiny::div(
        shiny::checkboxInput(ns("s4_do_bc"), "Apply batch correction (voom/lmFit)",
          value = shiny::isolate(input$s4_do_bc) %||% FALSE),
        shiny::uiOutput(ns("s4_bc_section"))
      )

      # ── Column 1: Before PCA (always visible after PCA computed) ─────────────
      col1_title <- if (bc_done) "Before batch correction" else "PCA"
      col1 <- shiny::div(
        shiny::tags$p(style = "font-weight:700;font-size:0.9em;text-align:center;margin-bottom:10px;
          border-bottom:2px solid #dee2e6;padding-bottom:6px;", col1_title),
        shiny::tags$b(style = "font-size:0.82em;", "Colour by"),
        if (length(meta_cols) > 0L)
          shiny::selectInput(ns("s4_pca_color"), label = NULL, choices = meta_cols,
            selected = shiny::isolate(input$s4_pca_color) %||% meta_cols[1], width = "100%")
        else shiny::helpText("No metadata variables."),
        shiny::uiOutput(ns("s4_color_warn")),
        shiny::fluidRow(
          shiny::column(6,
            shiny::tags$b(style = "font-size:0.82em;", "X axis"),
            shiny::selectInput(ns("s4_pc_x"), label = NULL, choices = pc_choices,
              selected = shiny::isolate(input$s4_pc_x) %||% "PC1", width = "100%")),
          shiny::column(6,
            shiny::tags$b(style = "font-size:0.82em;", "Y axis"),
            shiny::selectInput(ns("s4_pc_y"), label = NULL, choices = pc_choices,
              selected = shiny::isolate(input$s4_pc_y) %||%
                if (n_pcs_avail >= 2L) "PC2" else "PC1", width = "100%"))
        ),
        shiny::tags$p(style = "font-size:0.75em;color:#9ca3af;margin:0 0 6px 0;",
          "Colour/PC changes redraw instantly — no re-PCA."),
        shiny::plotOutput(ns("pca_before"), height = "320px"),
        .pp_expand_btn(ns("expand_pca_before")),
        shiny::tags$p(style = "font-size:0.8em;font-weight:600;color:#6b7280;
          margin:12px 0 4px 0;text-align:center;", "Scree plot"),
        shiny::plotOutput(ns("scree_before"), height = "160px"),
        shiny::tags$p(style = "font-size:0.75em;color:#9ca3af;margin:2px 0 0 0;",
          "Highlighted bars = selected PCs.")
      )

      # ── Column 2: After PCA (or placeholder) ──────────────────────────────────
      col2 <- if (!isTRUE(input$s4_do_bc)) {
        shiny::div(
          style = "display:flex;align-items:center;justify-content:center;
            min-height:420px;border:1px dashed #dee2e6;border-radius:8px;",
          shiny::div(style = "text-align:center;color:#9ca3af;font-size:0.85em;padding:24px;",
            shiny::HTML("Check <b>Apply batch correction</b> above<br>to compare PCA before and after.")))
      } else if (!bc_done) {
        shiny::div(
          style = "display:flex;align-items:center;justify-content:center;
            min-height:420px;border:1px dashed #dee2e6;border-radius:8px;",
          shiny::div(style = "text-align:center;color:#9ca3af;font-size:0.85em;padding:24px;",
            shiny::HTML("Select batch and biological variables above,<br>
              then click <b>Run Batch Correction</b>.")))
      } else {
        n_pcs_bc <- if (!is.null(rv$s4_pca_df_after)) ncol(rv$s4_pca_df_after) else 2L
        pc_ch_bc <- paste0("PC", 1:n_pcs_bc)
        shiny::div(
          shiny::tags$p(style = "font-weight:700;font-size:0.9em;text-align:center;margin-bottom:10px;
            border-bottom:2px solid #dee2e6;padding-bottom:6px;", "After batch correction"),
          shiny::tags$b(style = "font-size:0.82em;", "Colour by"),
          if (length(meta_cols) > 0L)
            shiny::selectInput(ns("s4_pca_after_color"), label = NULL, choices = meta_cols,
              selected = shiny::isolate(input$s4_pca_after_color) %||% meta_cols[1], width = "100%")
          else shiny::helpText("No metadata."),
          shiny::fluidRow(
            shiny::column(6,
              shiny::tags$b(style = "font-size:0.82em;", "X axis"),
              shiny::selectInput(ns("s4_pc_after_x"), label = NULL, choices = pc_ch_bc,
                selected = shiny::isolate(input$s4_pc_after_x) %||% "PC1", width = "100%")),
            shiny::column(6,
              shiny::tags$b(style = "font-size:0.82em;", "Y axis"),
              shiny::selectInput(ns("s4_pc_after_y"), label = NULL, choices = pc_ch_bc,
                selected = shiny::isolate(input$s4_pc_after_y) %||%
                  if (n_pcs_bc >= 2L) "PC2" else "PC1", width = "100%"))
          ),
          shiny::tags$p(style = "font-size:0.75em;color:#9ca3af;margin:0 0 6px 0;",
            "Colour/PC changes redraw instantly - no re-PCA."),
          shiny::plotOutput(ns("pca_after"), height = "320px"),
          .pp_expand_btn(ns("expand_pca_after")),
          shiny::tags$p(style = "font-size:0.8em;font-weight:600;color:#6b7280;
            margin:12px 0 4px 0;text-align:center;", "Scree plot"),
          shiny::plotOutput(ns("scree_after"), height = "160px"),
          shiny::tags$p(style = "font-size:0.75em;color:#9ca3af;margin:2px 0 0 0;",
            "Highlighted bars = selected PCs.")
        )
      }

      shiny::tagList(
        pca_top_ui,
        bc_controls_ui,
        shiny::hr(style = "margin:10px 0;"),
        shiny::fluidRow(
          shiny::column(6, col1),
          shiny::column(6, col2)
        ),
        shiny::hr(style = "margin:10px 0;"),
        shiny::textAreaInput(ns("s4_notes"), "\U0001f4dd Step notes (optional — included in report):",
          value = shiny::isolate(rv$note_s4), rows = 2,
          placeholder = "e.g. Batch corrected for sequencing run; retained condition effect.",
          width = "100%"),
        finalize_ui
      )
    })

    output$s4_color_warn <- shiny::renderUI({ .pp_color_warn(rv$meta_s, input$s4_pca_color) })

    output$s4_bc_section <- shiny::renderUI({
      if (!isTRUE(input$s4_do_bc)) return(NULL)
      meta_cols <- setdiff(colnames(rv$meta_s), "SampleID")
      # Use persisted selections after BC has been run (survive renderUI rebuilds)
      sel_batch  <- rv$bc_batch_cols_used  %||% shiny::isolate(input$s4_batch_cols)
      sel_effect <- rv$bc_effect_cols_used %||% shiny::isolate(input$s4_effect_cols)
      shiny::tags$details(open = TRUE,
        style = "border:1px solid #dee2e6;border-radius:6px;padding:8px 12px;margin-top:4px;",
        shiny::tags$summary(style = "cursor:pointer;font-weight:600;font-size:0.88em;user-select:none;",
          "Batch correction options"),
        shiny::div(style = "margin-top:8px;",
          shiny::fluidRow(
            shiny::column(6,
              shiny::tags$b(style = "font-size:0.88em;", "Variable(s) to remove (batch)"),
              shiny::helpText("Technical variation to eliminate (e.g. DatasetID, Batch).",
                style = "font-size:0.78em;margin:0;"),
              shiny::div(style = "max-height:110px;overflow-y:auto;border:1px solid #dee2e6;border-radius:4px;padding:4px 8px;",
                shiny::checkboxGroupInput(ns("s4_batch_cols"), label = NULL, choices = meta_cols,
                  selected = sel_batch))),
            shiny::column(6,
              shiny::tags$b(style = "font-size:0.88em;", "Variable(s) to retain (signal)"),
              shiny::helpText("Meaningful signal to protect from removal (e.g. Condition, Age). Optional.",
                style = "font-size:0.78em;margin:0;"),
              shiny::div(style = "max-height:110px;overflow-y:auto;border:1px solid #dee2e6;border-radius:4px;padding:4px 8px;",
                shiny::checkboxGroupInput(ns("s4_effect_cols"), label = NULL, choices = meta_cols,
                  selected = sel_effect)))
          ),
          shiny::div(style="margin-top:8px;font-size:0.78em;color:#6b7280;",
            "Batch correction follows the approach described in ",
            shiny::tags$a(href="https://elifesciences.org/articles/88623", target="_blank",
              style="color:#3b82f6;", "voyAGEr (eLife 2024)"),
            " and ",
            shiny::tags$a(href="https://academic.oup.com/nargab/article/8/2/lqag057/8703711#565027660",
              target="_blank", style="color:#3b82f6;", "markeR (NAR Genomics and Bioinformatics 2026)"),
            "."),
          shiny::div(style = "margin-top:10px;",
            shiny::actionButton(ns("bc_run_btn"), "\U0001f4ca Run Batch Correction",
              class = "btn-sm btn-primary"))
        )
      )
    })

    # PCA COMPUTATION (once, or re-run when settings change)
    shiny::observeEvent(input$s4_run_btn, {
      shiny::req(rv$step>=4L, rv$data_s3)
      n_g <- nrow(rv$data_s3); n_s <- ncol(rv$data_s3)
      do_scale  <- isTRUE(shiny::isolate(input$s4_pca_scale)  %||% FALSE)
      do_center <- isTRUE(shiny::isolate(input$s4_pca_center) %||% TRUE)
      shiny::withProgress(
        message=paste0("\U0001f4ca Computing PCA - ",format(n_g,big.mark=",")," genes × ",n_s," samples…"),
        value=0, {
          shiny::incProgress(0.4,"Decomposing variance…")
          result <- tryCatch(
            .pp_run_pca(rv$data_s3, n_pcs=min(10L,n_s-1L,n_g-1L),
                        scale=do_scale, center=do_center),
            error=function(e){shiny::showNotification(paste("PCA error:",e$message),type="error");NULL})
          if (!is.null(result)) {
            rv$s4_pca_df_before <- result$df
            rv$s4_pca_ev_before <- result$ev
            rv$s4_pca_scale     <- do_scale
            rv$s4_pca_center    <- do_center
            n_pcs <- ncol(result$df); pc_ch <- paste0("PC",1:n_pcs)
            shiny::updateSelectInput(session,"s4_pc_x",choices=pc_ch,selected="PC1")
            shiny::updateSelectInput(session,"s4_pc_y",choices=pc_ch,
              selected=if(n_pcs>=2L)"PC2" else "PC1")
          }
          shiny::incProgress(0.6,"Done.")
          rv$s4_run <- TRUE
        })
    })

    # PCA OPTIONS PANEL — separate renderUI so checkbox changes react without
    # re-rendering the whole step-4 layout (plots, BC section, etc.)
    output$s4_pca_top_ui <- shiny::renderUI({
      if (!rv$s4_run) return(NULL)   # pre-run panel handles options itself
      # Read inputs WITHOUT isolate so changes propagate
      cur_scale  <- input$s4_pca_scale  %||% rv$s4_pca_scale
      cur_center <- input$s4_pca_center %||% rv$s4_pca_center
      pca_stale  <- !isTRUE(cur_scale  == rv$s4_pca_scale) ||
                    !isTRUE(cur_center == rv$s4_pca_center)
      shiny::div(
        style = "background:#f9fafb;border:1px solid #e5e7eb;border-radius:6px;padding:10px 14px;margin-bottom:10px;",
        shiny::div(style = "display:flex;align-items:flex-start;gap:24px;flex-wrap:wrap;",
          shiny::div(
            shiny::tags$p(style = "font-size:0.82em;font-weight:600;margin:0 0 4px 0;", "PCA options:"),
            shiny::checkboxInput(ns("s4_pca_center"), "Center (subtract mean per gene)", value = cur_center),
            shiny::checkboxInput(ns("s4_pca_scale"),  "Scale (divide by SD per gene)",  value = cur_scale)
          ),
          shiny::div(style = "display:flex;align-items:center;gap:8px;margin-top:22px;",
            if (pca_stale)
              shiny::div(
                shiny::tags$span(style = "font-size:0.8em;color:#b45309;font-weight:600;margin-right:4px;",
                  "⚠ Settings changed —"),
                shiny::actionButton(ns("s4_reset_pca"), "↺ Reset & Recompute",
                  class = "btn-warning btn-sm", style = "padding:3px 10px;font-size:0.8em;"))
            else
              shiny::tagList(
                shiny::tags$span(style = "font-size:0.75em;color:#9ca3af;",
                  paste0(if (cur_center) "centered" else "not centered",
                         ", ", if (cur_scale) "scaled" else "not scaled")),
                shiny::actionButton(ns("s4_reset_pca"), "↺ Reset PCA",
                  class = "btn-outline-secondary btn-sm",
                  style = "font-size:0.75em;padding:3px 9px;"))
          )
        )
      )
    })

    # RESET PCA — clears results so the user can change options and re-run
    shiny::observeEvent(input$s4_reset_pca, {
      rv$s4_run           <- FALSE
      rv$s4_bc_run        <- FALSE
      rv$s4_pca_df_before <- NULL; rv$s4_pca_ev_before <- NULL
      rv$s4_pca_df_after  <- NULL; rv$s4_pca_ev_after  <- NULL
      rv$s4_bc_result     <- NULL
    })

    # PCA RENDERING - reacts to colour + PC selectors (no recomputation)
    output$pca_before <- shiny::renderPlot({
      shiny::req(rv$s4_run, rv$s4_pca_df_before, rv$s4_pca_ev_before)
      col  <- input$s4_pca_color
      pc_x <- input$s4_pc_x %||% "PC1"
      pc_y <- input$s4_pc_y %||% "PC2"
      title <- if (isTRUE(input$s4_do_bc) && rv$s4_bc_run) "PCA - before batch correction" else "PCA"
      .pp_pca_from_df(rv$s4_pca_df_before, rv$s4_pca_ev_before, rv$meta_s, col, title, pc_x, pc_y)
    })

    output$scree_before <- shiny::renderPlot({
      shiny::req(rv$s4_run, rv$s4_pca_ev_before)
      pc_x <- input$s4_pc_x %||% "PC1"
      pc_y <- input$s4_pc_y %||% "PC2"
      .pp_scree_plot(rv$s4_pca_ev_before, n_show=min(10L,length(rv$s4_pca_ev_before)), pc_x, pc_y)
    })

    # (s4_after_pc_selectors removed - PC selectors are now embedded in the two-column layout)

    # BATCH CORRECTION
    shiny::observeEvent(input$bc_run_btn, {
      shiny::req(rv$step>=4L, rv$norm_linear, rv$meta_s)
      batch_cols  <- shiny::isolate(input$s4_batch_cols)
      effect_cols <- shiny::isolate(input$s4_effect_cols) %||% character(0)
      if (length(batch_cols)<1L){shiny::showNotification("Select at least one batch variable.",type="warning");return()}
      # Persist the column selections so .do_finalize can read them even after UI rebuilds
      rv$bc_batch_cols_used  <- batch_cols
      rv$bc_effect_cols_used <- effect_cols
      norm_counts <- as.matrix(rv$norm_linear); meta <- rv$meta_s; n_s <- ncol(norm_counts)
      shiny::withProgress(
        message=paste0("⚙️ Batch correction - batch: [",paste(batch_cols,collapse=", "),"] - ",n_s," samples…"),
        value=0, {
          shiny::incProgress(0.15,"Building design matrices…")
          corrdf <- tryCatch({
            bc_form <- stats::reformulate(paste(batch_cols,collapse=" + "),intercept=FALSE)
            mmbatch <- stats::model.matrix(bc_form,data=meta)
            # Protect biological variables if specified; otherwise use intercept-only
            if (length(effect_cols) > 0L) {
              eff_form <- stats::reformulate(paste(effect_cols,collapse=" + "),intercept=FALSE)
              mmkeep   <- stats::model.matrix(eff_form,data=meta)
            } else {
              mmkeep <- matrix(1, nrow=nrow(meta), ncol=1,
                               dimnames=list(NULL,"intercept"))
            }
            mm <- cbind(mmkeep,mmbatch)
            shiny::incProgress(0.30,"Running voom + lmFit…")
            D0  <- edgeR::DGEList(norm_counts); D0 <- edgeR::calcNormFactors(D0)
            y   <- limma::voom(D0,mm,plot=FALSE)
            fit <- limma::lmFit(y,mm)
            n_keep <- ncol(mmkeep)
            beta   <- fit$coefficients[,-(1:n_keep),drop=FALSE]; beta[is.na(beta)] <- 0
            corrcounts <- as.matrix(y$E) - beta %*% t(mmbatch)
            offset <- apply(corrcounts,1,min) - apply(log2(norm_counts+1),1,min)
            as.data.frame(2^(corrcounts - offset))
          }, error=function(e){
            shiny::showNotification(paste("Batch correction failed:",e$message),type="error",duration=12); NULL})
          if (!is.null(corrdf)) {
            shiny::incProgress(0.25,"Computing post-correction PCA…")
            log_bc    <- log2(as.matrix(corrdf)+1)
            pca_after <- tryCatch(
              .pp_run_pca(log_bc, n_pcs=min(10L,n_s-1L),
                          scale=isTRUE(rv$s4_pca_scale), center=isTRUE(rv$s4_pca_center)),
              error=function(e){shiny::showNotification(paste("Post-BC PCA:",e$message),type="error");NULL})
            rv$s4_bc_result    <- corrdf
            rv$s4_pca_df_after <- if (!is.null(pca_after)) pca_after$df else NULL
            rv$s4_pca_ev_after <- if (!is.null(pca_after)) pca_after$ev else NULL
            rv$s4_bc_run       <- TRUE
            if (!is.null(pca_after)) {
              pc_ch_bc <- paste0("PC",1:ncol(pca_after$df))
              shiny::updateSelectInput(session,"s4_pc_after_x",choices=pc_ch_bc,selected="PC1")
              shiny::updateSelectInput(session,"s4_pc_after_y",choices=pc_ch_bc,
                selected=if(length(pc_ch_bc)>=2L)"PC2" else "PC1")
            }
            glog(paste0("Batch Correction: removed [",paste(batch_cols,collapse=", "),
              "], retained [",paste(effect_cols,collapse=", "),"] - ",n_s," samples."))
          }
          shiny::incProgress(0.30,"Done.")
        })
    })

    # AFTER-BC PCA RENDERING - reacts to colour + PC selectors
    output$pca_after <- shiny::renderPlot({
      shiny::req(rv$s4_run, rv$s4_bc_run, rv$s4_pca_df_after, rv$s4_pca_ev_after)
      col  <- input$s4_pca_after_color %||% input$s4_pca_color
      pc_x <- input$s4_pc_after_x %||% "PC1"
      pc_y <- input$s4_pc_after_y %||% "PC2"
      .pp_pca_from_df(rv$s4_pca_df_after, rv$s4_pca_ev_after, rv$meta_s, col,
                       "PCA - after batch correction", pc_x, pc_y)
    })
    output$scree_after <- shiny::renderPlot({
      shiny::req(rv$s4_bc_run, rv$s4_pca_ev_after)
      pc_x <- input$s4_pc_after_x %||% "PC1"
      pc_y <- input$s4_pc_after_y %||% "PC2"
      .pp_scree_plot(rv$s4_pca_ev_after, n_show=min(10L,length(rv$s4_pca_ev_after)), pc_x, pc_y)
    })

    # ── Finalise ──────────────────────────────────────────────────────────────
    .do_finalize <- function(use_bc) {
      shiny::req(rv$step >= 4L)
      rv$note_s4 <- shiny::isolate(input$s4_notes) %||% ""
      if (use_bc) {
        bc <- rv$s4_bc_result
        if (is.null(bc)){shiny::showNotification("Batch correction not completed.",type="error");return()}
        rv$final_data <- bc  # linear scale (2^corrected)
      } else {
        # Use TMM-CPM (non-log, linear scale) - NOT log2 data
        rv$final_data <- as.data.frame(rv$norm_linear)
      }
      n_g  <- nrow(rv$final_data); n_s <- ncol(rv$final_data)
      orig_g <- rv$orig_n_genes %||% n_g; orig_s <- rv$orig_n_samples %||% n_s
      genes_removed   <- orig_g - n_g; genes_kept <- n_g; samples_removed <- orig_s - n_s
      norm_method <- rv$s3_method_used %||% shiny::isolate(input$s3_method %||% "TMM")
      bc_vars     <- if (use_bc) rv$bc_batch_cols_used else NULL

      rv$pp_summary <- list(
        genes_kept=genes_kept, genes_removed=genes_removed,
        n_samples=n_s, samples_removed=samples_removed,
        norm_method=norm_method, bc_applied=use_bc,
        bc_vars=bc_vars,
        bc_effect=if(use_bc) rv$bc_effect_cols_used else NULL,
        data_scale=if(use_bc) {
          paste0("non-log (linear, batch-corrected ", if(norm_method=="DESeq2") "VST→CPM" else "CPM", ")")
        } else if (norm_method == "DESeq2") {
          "non-log (VST→linear scale)"
        } else {
          "non-log (TMM-CPM)"
        },
        bc_label=if(use_bc)
          paste0("batch-corrected ([",paste(bc_vars,collapse=", "),"])")
          else "normalised (no batch correction)"
      )
      glog(paste0("Preprocessing finalised",
        if(use_bc) paste0(" WITH batch correction [",paste(bc_vars,collapse=", "),"]")
        else " (no batch correction)",
        " - ",format(genes_kept,big.mark=",")," genes kept",
        if(genes_removed>0L) paste0(", ",format(genes_removed,big.mark=",")," removed") else "",
        " - ",n_s," samples. Final data: ",rv$pp_summary$data_scale,"."))
      # ── Snapshot plots for the HTML report ───────────────────────────────────
      # Normalisation boxplots: use the user's current display settings
      tryCatch({
        disp     <- shiny::isolate(input$s3_sample_display %||% "all")
        n_box    <- as.integer(shiny::isolate(input$s3_n_box %||% 20L))
        col_by   <- shiny::isolate(input$s3_color_by)
        hide_xl  <- shiny::isolate(isTRUE(input$s3_hide_xlabels))
        if (!is.null(rv$data_s2)) {
          log_d <- log2(as.matrix(rv$data_s2) + 1)
          rv$report_s3_before <- .pp_make_boxplot(log_d, colnames(log_d), disp, n_box,
                                                   col_by, rv$meta_s,
                                                   "Before normalisation", hide_xlabels = hide_xl)
        }
        if (!is.null(rv$norm_linear)) {
          log_n <- log2(as.matrix(rv$norm_linear) + 1)
          rv$report_s3_after  <- .pp_make_boxplot(log_n, colnames(log_n), disp, n_box,
                                                   col_by, rv$meta_s,
                                                   "After normalisation", hide_xlabels = hide_xl)
        }
      }, error = function(e) glog(paste0("Report snapshot (norm): ", e$message)))
      # PCA plots: use the user's chosen colour variable
      tryCatch({
        pca_col  <- shiny::isolate(input$s4_pca_color)
        pc_x     <- shiny::isolate(input$s4_pc_x %||% "PC1")
        pc_y     <- shiny::isolate(input$s4_pc_y %||% "PC2")
        if (rv$s4_run && !is.null(rv$s4_pca_df_before)) {
          title_b <- if (use_bc) "PCA - before batch correction" else "PCA"
          rv$report_s4_before <- .pp_pca_from_df(rv$s4_pca_df_before, rv$s4_pca_ev_before,
                                                   rv$meta_s, pca_col, title_b, pc_x, pc_y)
        }
        if (use_bc && rv$s4_bc_run && !is.null(rv$s4_pca_df_after)) {
          pca_col_after <- shiny::isolate(input$s4_pca_after_color %||% input$s4_pca_color)
          pc_ax <- shiny::isolate(input$s4_pc_after_x %||% "PC1")
          pc_ay <- shiny::isolate(input$s4_pc_after_y %||% "PC2")
          rv$report_s4_after  <- .pp_pca_from_df(rv$s4_pca_df_after, rv$s4_pca_ev_after,
                                                   rv$meta_s, pca_col_after,
                                                   "PCA - after batch correction", pc_ax, pc_ay)
        }
      }, error = function(e) glog(paste0("Report snapshot (PCA): ", e$message)))
      # Filter density plots (step 2)
      tryCatch({
        if (!is.null(rv$data_s1) && rv$s2_run) {
          thr <- rv$s2_thr %||% 1.0; d <- rv$data_s1
          lm  <- log2(rowMeans(as.matrix(d)) + 1)
          rv$report_s2_before <- ggplot2::ggplot(data.frame(log_mean=lm), ggplot2::aes(x=log_mean)) +
            ggplot2::geom_density(fill="#EBB43E",alpha=0.45,colour="#b8860b",linewidth=0.8) +
            ggplot2::geom_vline(xintercept=thr,linetype="dashed",colour="firebrick",linewidth=1.2) +
            ggplot2::annotate("text",x=thr,y=Inf,vjust=1.8,hjust=-0.1,
              label=paste0("threshold = ",round(thr,3)),colour="firebrick",size=3.5) +
            ggplot2::theme_bw() +
            ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
                           panel.border=ggplot2::element_blank(),
                           axis.line=ggplot2::element_line(colour="black")) +
            ggplot2::xlab("Mean log2(count+1)") + ggplot2::ylab("Density") +
            ggplot2::ggtitle(paste0("All genes (threshold >= ",round(thr,3),")"))
          kv2 <- lm[lm >= thr]
          rv$report_s2_after <- if (length(kv2) > 0)
            ggplot2::ggplot(data.frame(log_mean=kv2), ggplot2::aes(x=log_mean)) +
              ggplot2::geom_density(fill="#28a745",alpha=0.4,colour="#155724",linewidth=0.8) +
              ggplot2::theme_bw() +
              ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),
                             panel.border=ggplot2::element_blank(),
                             axis.line=ggplot2::element_line(colour="black")) +
              ggplot2::xlab("Mean log2(count+1)") + ggplot2::ylab("Density") +
              ggplot2::ggtitle(paste0("Retained genes (",format(length(kv2),big.mark=","),")"))
          else NULL
        }
      }, error = function(e) glog(paste0("Report snapshot (filter): ", e$message)))
      # Store data needed for interactive report (box stats + PCA coords + color prefs)
      tryCatch({
        meta_cols_ri <- setdiff(colnames(rv$meta_s), "SampleID")
        rv$report_interactive_data <- list(
          s3_color_by     = shiny::isolate(input$s3_color_by),
          s4_pca_color    = shiny::isolate(input$s4_pca_color),
          pc_x            = shiny::isolate(input$s4_pc_x %||% "PC1"),
          pc_y            = shiny::isolate(input$s4_pc_y %||% "PC2"),
          meta_cols       = meta_cols_ri,
          s3_display_mode = shiny::isolate(input$s3_sample_display %||% "all"),
          s3_n_box        = as.integer(shiny::isolate(input$s3_n_box %||% 20L))
        )
      }, error = function(e) glog(paste0("Interactive report data: ", e$message)))
      # Record the fingerprint of the final data so we can detect genuinely new data later
      fd <- rv$final_data
      if (!is.null(fd) && ncol(fd) > 0L) {
        .finalized_fp(paste(nrow(fd), ncol(fd), colnames(fd)[1], colnames(fd)[ncol(fd)], sep="|"))
      }
      rv$downloads_available <- TRUE; rv$finalized <- rv$finalized+1L
      rv$step <- 5L; rv$active_step <- 5L
      shiny::showNotification(
        shiny::HTML(paste0("✅ Preprocessing complete - ",format(genes_kept,big.mark=","),
          " genes × ",n_s," samples (",rv$pp_summary$data_scale,").")),
        type="default", duration=8)
    }

    shiny::observeEvent(input$finalize_btn,  { .do_finalize(FALSE) })
    shiny::observeEvent(input$finalize_bc,   { .do_finalize(TRUE)  })
    shiny::observeEvent(input$finalize_norm, { .do_finalize(FALSE) })

    shiny::observeEvent(input$back_to_s4, {
      rv$active_step <- 4L; rv$step <- 4L; rv$s4_run <- FALSE; rv$s4_bc_run <- FALSE
      rv$final_data <- NULL; rv$finalized <- 0L; rv$pp_summary <- NULL
      rv$downloads_available <- FALSE
      rv$report_s4_before <- NULL; rv$report_s4_after <- NULL
      .finalized_fp(NULL)
    })

    # =========================================================================
    # SUMMARY CARD
    # =========================================================================

    output$pp_summary_card <- shiny::renderUI({
      shiny::req(rv$finalized>0L, !is.null(rv$pp_summary))
      s <- rv$pp_summary
      .row <- function(label, value)
        shiny::tags$div(class="pp-summary-row",
          shiny::tags$span(class="pp-summary-label", label),
          shiny::tags$span(value))
      shiny::div(class="pp-summary-card",
        shiny::tags$h6(style="color:#EBB43E;font-weight:700;margin-bottom:12px;",
          "✅ Preprocessing Summary"),
        .row("Final dataset:", paste0(format(s$genes_kept,big.mark=",")," genes × ",s$n_samples," samples")),
        .row("Genes kept:",    format(s$genes_kept,big.mark=",")),
        if (s$genes_removed>0L)
          .row("Genes removed:", paste0(format(s$genes_removed,big.mark=",")," (below expression threshold)")),
        if (s$samples_removed>0L)
          .row("Samples removed:", as.character(s$samples_removed)),
        .row("Normalisation:", s$norm_method),
        .row("Batch correction:",
          if (s$bc_applied) {
            bc_part  <- if (!is.null(s$bc_vars)    && length(s$bc_vars)    > 0L)
                          paste0("batch: [",   paste(s$bc_vars,   collapse=", "), "]") else NULL
            eff_part <- if (!is.null(s$bc_effect) && length(s$bc_effect) > 0L)
                          paste0("signal retained: [",paste(s$bc_effect, collapse=", "), "]") else NULL
            parts <- Filter(Negate(is.null), list(bc_part, eff_part))
            paste0("✓ Applied - ", paste(parts, collapse=" | "))
          } else "Not applied"),
        .row("Data scale:", shiny::tags$span(style="font-weight:600;color:#065f46;", s$data_scale)),
        .row("Data used downstream:",
          if (s$bc_applied) "Batch-corrected expression data" else "Normalised expression data"),
        shiny::hr(style="margin:12px 0;"),
        shiny::div(
          style="background:#f0fdf4;border:1px solid #bbf7d0;border-radius:6px;padding:12px 16px;",
          shiny::tags$p(style="font-weight:600;color:#15803d;margin-bottom:4px;",
            "✅ Your preprocessed data is ready for analysis."),
          shiny::tags$p(style="font-size:0.86em;color:#374151;margin-bottom:0;",
            "Use the tabs at the top of the page to navigate to any analysis - the preprocessed data is automatically available everywhere.")
        ),
        shiny::div(style="margin-top:10px;font-size:0.82em;color:#9ca3af;",
          "Downloads are available in the sidebar (optional).")
      )
    })

    # =========================================================================
    # HELPERS
    # =========================================================================

    .align_meta <- function(keep) {
      meta <- rv$meta_s; if (is.null(meta)) return()
      id_col <- if ("SampleID" %in% colnames(meta)) "SampleID" else colnames(meta)[1]
      idx    <- match(keep, as.character(meta[[id_col]])); valid <- !is.na(idx)
      if (sum(valid)>=1L) rv$meta_s <- meta[idx[valid],,drop=FALSE]
    }

    list(final_data = shiny::reactive(rv$final_data),
         finalized  = shiny::reactive(rv$finalized))
  })
}
