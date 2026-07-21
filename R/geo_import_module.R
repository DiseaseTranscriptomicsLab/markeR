# =============================================================================
# geo_import_module.R
# Adaptive GEO Data Import Pipeline for RNA-Seq - markeR Shiny Module
#
# Phase 1 Analysis Summary (20 GEO accessions)
# ─────────────────────────────────────────────
# FORMAT HETEROGENEITY
#   XLSX  : GSE334803, GSE331505, GSE299445
#   XLS   : GSE299859, GSE334767
#   TSV   : GSE335198
#   TXT   : GSE335122, GSE303976, GSE300974, GSE318564 (TAR of TXTs), GSE320441, GSE314271
#   CSV   : GSE299456, GSE296187, GSE314271, GSE309671 (Salmon .sf)
#   SF    : GSE309671  → Salmon quant.sf per-sample files (gene-level tximport needed)
#   TAR   : GSE318564  → scRNA-seq 10x Genomics barcodes/features/matrix (EXCLUDE / WARN)
#
# EXCLUDE (non-RNA-Seq bulk)
#   GSE334799 : BED/BW  → ChIP-seq / ATAC-seq chromatin accessibility
#   GSE300943 : Hi-C    → chromatin conformation, no expression matrix
#   GSE305083 : CUT&RUN → chromatin profiling, no expression matrix
#
# METADATA MAPPING CHALLENGES
#   • GSE335177 (SuperSeries) : multiple GPL platforms → must iterate sub-GSEs
#   • GSE309671 (Salmon SF)   : sample-per-file naming (GSM######_quant.sf)
#   • GSE318564 (scRNA/TAR)   : 10x barcodes are not bulk samples
#   • In many non-standard cases, Series Matrix "Sample_title" ≠ column header
#     in supplementary file → fuzzy matching required
#
# EXCLUSION LOGIC
#   Keyword signals in platform annotation / GSE title used to flag & skip:
#   Hi-C, ChIP, ATAC, CUT&RUN, CUT&TAG, IDAT, BED, bigWig, methylation, bisulfite
#
# =============================================================================

#' @importFrom shiny moduleServer NS reactive reactiveVal observeEvent req
#'   withProgress incProgress showNotification renderUI uiOutput tagList tags
#'   selectInput actionButton helpText HTML p strong hr h5 updateSelectInput
#'   updateRadioButtons
#' @importFrom bslib card card_header
#' @importFrom DT renderDT DTOutput datatable
#' @importFrom GEOquery getGEO getGEOSuppFiles
#' @importFrom Biobase exprs pData annotation experimentData
#' @importFrom readxl read_xls read_xlsx excel_sheets
#' @importFrom data.table fread
#' @importFrom R.utils gunzip
#' @importFrom shinyWidgets pickerInput pickerOptions
# Note: stringdist and tximport are optional (Suggests); accessed via :: with requireNamespace guards

# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────

.GEO_EXCLUDE_KEYWORDS <- c(
  "Hi-C", "HiC", "ChIP", "ATAC", "CUT&RUN", "CUT&TAG",
  "CUT.RUN", "CUT.TAG", "IDAT", "bisulfite", "methylation",
  "RRBS", "WGBS", "bigwig", "bigWig", "BED", "bedgraph",
  "scRNA", "single.cell", "single cell", "10[Xx]",
  # DNA-repair / chromatin assays (batch 2 findings)
  "SAR-seq", "SARseq", "SAR seq", "END-seq", "ENDseq", "END seq",
  # Small/non-coding RNA assays
  "miRNA", "mirna", "small.rna", "smallrna", "small RNA",
  "non-coding rna", "noncoding rna", "ncRNA profiling"
)

.EXPR_KEYWORDS <- c(
  "count", "counts", "tpm", "fpkm", "rpkm", "cpm",
  "expression", "expr", "raw", "normalized", "norm",
  "readcount", "read_count", "gene_count", "featurecount",
  # Additional quantification tools (batch 2 findings)
  "rsem", "tmm", "htseq", "kallisto", "salmon",
  "matrix", "profile"
)

.NON_EXPR_KEYWORDS <- c(
  "peak", "bigwig", "bw", "bed", "idat", "bam", "bai",
  "metadata", "phenotype", "samplesheet", "readme", "annotation",
  "sra", "fastq", "barcode",
  # NOTE: "feature" removed - false-positive on "featureCounts" expression files
  # Scoped 10x file detection is handled in classify_supp_file via barcodes/features.tsv
  # Analysis results (not raw counts): should not be used as expression matrix
  "gsea", "enrichment", "pathway", "volcano", "gene_set",
  "diffexpr", "diff_expr", "deseq_result", "edger_result",
  "limma_result", "gseaResult", "ranked_gene"
)

.SALMON_KEYWORDS <- c("quant.sf", "quant.genes.sf", "\\.sf$")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1 - ASSAY-TYPE CLASSIFICATION
# ─────────────────────────────────────────────────────────────────────────────

#' Classify a GEO object's assay type
#'
#' Returns one of: "bulk_rnaseq", "microarray", "scrna", "excluded", "unknown"
#' @param obj A GEO ExpressionSet or GSE object
#' @param gse_accession Character. Used only for error messages.
classify_geo_assay <- function(obj, gse_accession = "") {

  platform_title <- tryCatch(tolower(Biobase::annotation(obj)), error = function(e) "")
  if (length(platform_title)   == 0L) platform_title   <- ""

  experiment_title <- tryCatch(tolower(Biobase::experimentData(obj)@title), error = function(e) "")
  if (length(experiment_title) == 0L) experiment_title <- ""

  # Official GEO experiment type field - the most authoritative classification signal
  experiment_type <- tryCatch({
    other <- Biobase::experimentData(obj)@other
    if (is.list(other) && !is.null(other[["type"]])) {
      tolower(paste(other[["type"]], collapse = " "))
    } else ""
  }, error = function(e) "")
  if (length(experiment_type) == 0L) experiment_type <- ""

  # ── Priority 1: scRNA / single-cell ──────────────────────────────────────────
  # Check title + type: "[scRNA-seq]" often appears in the title even when
  # experiment_type says "Expression profiling by HTS".
  scrna_text <- paste(experiment_title, experiment_type, platform_title)
  if (grepl("10x|scrna|single.cell|single cell|cellranger|seurat|anndata|\\[sc ",
            scrna_text, ignore.case = TRUE)) {
    return(list(type   = "scrna",
                reason = paste0("Single-cell / 10x signals for ", gse_accession)))
  }

  # ── Priority 2: Non-coding RNA - check experiment_type only (authoritative) ──
  if (grepl("non-coding rna profiling|ncrna profiling|mirna sequencing|small rna sequencing",
            experiment_type, ignore.case = TRUE)) {
    return(list(type   = "excluded",
                reason = paste0("Non-coding/miRNA experiment type for ", gse_accession,
                                ": '", experiment_type, "'")))
  }

  # ── Priority 3: Expression profiling → confirmed RNA-seq or microarray ───────
  if (grepl("expression profiling by high throughput sequencing",
            experiment_type, ignore.case = TRUE)) {
    return(list(type   = "bulk_rnaseq",
                reason = "GEO experiment type confirmed as expression profiling by HTS"))
  }
  if (grepl("expression profiling by array", experiment_type, ignore.case = TRUE)) {
    return(list(type   = "microarray",
                reason = "GEO experiment type confirmed as microarray"))
  }

  # ── Priority 4: Unknown / "Other" experiment type ────────────────────────────
  # Do NOT exclude based on title keywords - titles of multi-assay studies (e.g.
  # "Hi-C and RNA-seq of X") would cause false exclusions.  Instead, pass through
  # and let the supplementary-file scorer be the gatekeeper.  If no expression
  # matrix is found in the files, import_geo_data returns status="warning".
  # Only the platform annotation (GPL string) is checked - it's more reliable.
  if (grepl("affymetrix|illumina.*bead|agilent", platform_title, ignore.case = TRUE)) {
    expr_mat <- tryCatch(Biobase::exprs(obj), error = function(e) NULL)
    if (!is.null(expr_mat) && nrow(expr_mat) > 0)
      return(list(type = "microarray", reason = "Microarray platform with expression matrix"))
  }

  return(list(type   = "bulk_rnaseq",
              reason = paste0("Experiment type '", experiment_type,
                              "' - proceeding to file-level classification")))
}


#' Classify a supplementary file by assay type
#'
#' Returns "rnaseq", "salmon", "scrna", "excluded", or "unknown"
classify_supp_file <- function(fname) {
  fl <- tolower(fname)

  # Salmon .sf
  if (any(sapply(.SALMON_KEYWORDS, function(p) grepl(p, fl, ignore.case = TRUE)))) {
    return("salmon")
  }

  # 10x / scRNA
  if (grepl("barcodes\\.tsv|features\\.tsv|matrix\\.mtx|filtered_feature_bc", fl)) {
    return("scrna")
  }

  # Excluded assay types
  for (kw in .GEO_EXCLUDE_KEYWORDS) {
    if (grepl(kw, fl, ignore.case = TRUE)) return("excluded")
  }

  # Non-expression files
  for (kw in .NON_EXPR_KEYWORDS) {
    if (grepl(kw, fl, ignore.case = TRUE)) return("non_expr")
  }

  "rnaseq"
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2 - SMART SUPPLEMENTARY FILE SCORER
# ─────────────────────────────────────────────────────────────────────────────

#' Score candidate supplementary files by expression-matrix likelihood
#'
#' @param file_names Character vector of basenames or full paths.
#' @return Named numeric vector of scores (0–1), highest = most likely expression.
score_supp_files <- function(file_names) {
  scores <- numeric(length(file_names))
  names(scores) <- file_names

  for (i in seq_along(file_names)) {
    fn <- file_names[i]
    if (is.na(fn) || !nzchar(fn)) { scores[i] <- 0; next }
    fname <- tolower(basename(fn))
    ext   <- tools::file_ext(fname)
    if (length(ext) == 0L || is.na(ext)) { scores[i] <- 0.1; next }

    # Base score from extension
    ext_score <- switch(ext,
      "gz"   = 0.5,
      "txt"  = 0.5,
      "tsv"  = 0.55,
      "csv"  = 0.55,
      "xls"  = 0.6,
      "xlsx" = 0.6,
      "sf"   = 0.7,   # Salmon - handle separately
      "tar"  = 0.2,   # Archive - needs inspection
      "zip"  = 0.2,
      0.1
    )

    # For .gz files, peek at the inner extension to refine the base score
    if (ext == "gz") {
      inner_ext <- tolower(tools::file_ext(sub("\\.gz$", "", fname, ignore.case = TRUE)))
      if (inner_ext %in% c("csv", "tsv")) ext_score <- 0.55
      if (inner_ext == "txt")              ext_score <- 0.50
      if (inner_ext %in% c("xls", "xlsx")) ext_score <- 0.60
    }

    # Expression keyword bonus
    expr_bonus <- sum(sapply(.EXPR_KEYWORDS,
      function(kw) grepl(kw, fname, ignore.case = TRUE))) * 0.1

    # Non-expression keyword penalty
    non_expr_penalty <- sum(sapply(.NON_EXPR_KEYWORDS,
      function(kw) grepl(kw, fname, ignore.case = TRUE))) * 0.2

    # Classify assay type
    assay_type <- classify_supp_file(fname)
    assay_penalty <- switch(assay_type,
      "excluded" = 1.0,
      "non_expr" = 0.5,
      "scrna"    = 0.8,
      0
    )

    scores[i] <- max(0, min(1, ext_score + expr_bonus - non_expr_penalty - assay_penalty))
  }

  scores
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3 - SWISS ARMY KNIFE FILE PARSER
# ─────────────────────────────────────────────────────────────────────────────

#' Unified expression file reader
#'
#' Handles .xlsx, .xls, .csv, .tsv, .txt, .gz (of any text format).
#' Returns a numeric matrix: rows = genes, cols = samples.
#' Uses data.table::fread for text formats and readxl for Excel.
#'
#' @param path      Full path to file (may be compressed).
#' @param orig_name Original filename for extension detection.
#' @param sheet     For Excel: sheet number or name (default = 1).
#' @return List with elements: `matrix` (numeric matrix), `warnings` (character).
parse_expression_file <- function(path, orig_name = NULL, sheet = 1) {

  fname <- if (!is.null(orig_name)) orig_name else basename(path)
  warnings_out <- character(0)

  # ── Step 1: Decompress if needed ──────────────────────────────────────────
  actual_path <- path
  inner_name  <- fname

  if (grepl("\\.gz$", fname, ignore.case = TRUE)) {
    inner_name  <- sub("\\.gz$", "", fname, ignore.case = TRUE)
    tmp_ungz    <- tempfile(fileext = paste0(".", tools::file_ext(inner_name)))
    tryCatch(
      R.utils::gunzip(path, destname = tmp_ungz, overwrite = TRUE, remove = FALSE),
      error = function(e) stop("Failed to decompress .gz: ", e$message)
    )
    actual_path <- tmp_ungz
  }

  if (grepl("\\.tar\\.gz$|\\.tar$", fname, ignore.case = TRUE)) {
    return(.parse_tar_archive(path, fname))
  }

  ext <- tolower(tools::file_ext(inner_name))

  # ── Step 2: Read by extension ──────────────────────────────────────────────
  df <- tryCatch({

    if (ext %in% c("xlsx", "xls")) {

      # Try to find the right sheet
      all_sheets <- tryCatch(
        readxl::excel_sheets(actual_path),
        error = function(e) as.character(sheet)
      )

      # Score sheets by name
      best_sheet <- .pick_best_excel_sheet(all_sheets)
      raw <- readxl::read_excel(actual_path, sheet = best_sheet, col_names = TRUE)
      as.data.frame(raw, stringsAsFactors = FALSE)

    } else if (ext %in% c("sf")) {

      # Salmon quant.sf - return tximport-ready path, not matrix
      return(list(matrix = NULL, type = "salmon_sf", path = actual_path,
                  warnings = "Salmon .sf file detected - use tximport for gene-level aggregation."))

    } else {

      # Generic text: use fread for speed and auto-detection
      dt <- data.table::fread(
        actual_path,
        header        = TRUE,
        check.names   = FALSE,
        stringsAsFactors = FALSE,
        data.table    = FALSE,
        fill          = TRUE,
        sep           = "auto"
      )
      dt

    }
  }, error = function(e) {
    stop("File read failed for '", fname, "': ", e$message)
  })

  if (is.null(df) || nrow(df) == 0) {
    stop("Empty file or failed to parse: '", fname, "'")
  }

  # ── Step 3: Identify gene column and expression columns ───────────────────
  result <- .extract_expr_matrix(df, fname)
  result$warnings <- c(warnings_out, result$warnings)
  result
}


#' Internal: pick best Excel sheet by name scoring
.pick_best_excel_sheet <- function(sheets) {
  if (length(sheets) == 1) return(sheets[1])

  scores <- sapply(sheets, function(s) {
    sl <- tolower(s)
    sum(sapply(.EXPR_KEYWORDS, function(k) grepl(k, sl))) * 2 -
    sum(sapply(.NON_EXPR_KEYWORDS, function(k) grepl(k, sl)))
  })

  sheets[which.max(scores)]
}


#' Internal: extract numeric expression matrix from a data.frame
.extract_expr_matrix <- function(df, fname = "") {

  warnings_out <- character(0)

  # Remove fully NA or fully empty columns
  df <- df[, sapply(df, function(col) sum(!is.na(col) & col != "") > 0), drop = FALSE]

  if (ncol(df) == 0 || nrow(df) == 0) {
    stop("Data frame is empty after removing blank columns in '", fname, "'")
  }

  # Classify each column as numeric or character
  num_fraction <- sapply(df, function(col) {
    suppressWarnings(mean(!is.na(as.numeric(as.character(col)))))
  })

  numeric_cols  <- which(num_fraction > 0.8)
  char_cols     <- which(num_fraction < 0.2)

  # Fallback: try transposed layout (samples as rows)
  if (length(numeric_cols) < 2 && ncol(df) > nrow(df)) {
    warnings_out <- c(warnings_out,
      "Matrix appears transposed - attempting automatic transpose.")
    df <- as.data.frame(t(df), stringsAsFactors = FALSE)
    # Re-assess
    num_fraction  <- sapply(df, function(col) suppressWarnings(
      mean(!is.na(as.numeric(as.character(col))))
    ))
    numeric_cols  <- which(num_fraction > 0.8)
    char_cols     <- which(num_fraction < 0.2)
  }

  if (length(numeric_cols) < 1) {
    stop("No numeric columns found in '", fname, "'")
  }

  # Gene column: prefer explicit "gene", "symbol", "id" header, else first char col
  gene_col_idx <- NA_integer_

  gene_headers <- c("gene", "gene_id", "gene_name", "gene_symbol",
                    "symbol", "feature", "id", "ensembl", "entrez", "name")

  for (gc in gene_headers) {
    hit <- which(tolower(colnames(df)) == gc)
    if (length(hit) > 0) { gene_col_idx <- hit[1]; break }
  }

  if (is.na(gene_col_idx) && length(char_cols) > 0) {
    gene_col_idx <- char_cols[1]
  }

  # Extract genes
  if (!is.na(gene_col_idx)) {
    genes <- trimws(as.character(df[[gene_col_idx]]))
  } else if (!is.null(rownames(df)) && !all(rownames(df) == as.character(seq_len(nrow(df))))) {
    genes <- rownames(df)
  } else {
    stop("Cannot determine gene identifier column in '", fname, "'")
  }

  # Build numeric matrix
  expr_cols <- setdiff(numeric_cols, gene_col_idx)
  expr      <- df[, expr_cols, drop = FALSE]

  # Coerce to numeric (handle comma decimals)
  expr[] <- lapply(expr, function(col) {
    x <- as.character(col)
    x <- gsub(",", ".", x)
    suppressWarnings(as.numeric(x))
  })

  mat <- as.matrix(expr)
  storage.mode(mat) <- "numeric"

  # Remove rows with all-NA
  valid <- !apply(mat, 1, function(r) all(is.na(r))) & !is.na(genes) & genes != ""
  mat   <- mat[valid, , drop = FALSE]
  genes <- genes[valid]

  # Remove duplicate gene names (keep first)
  dup <- duplicated(genes)
  if (any(dup)) {
    warnings_out <- c(warnings_out,
      paste0(sum(dup), " duplicate gene names removed (kept first occurrence)."))
    mat   <- mat[!dup, , drop = FALSE]
    genes <- genes[!dup]
  }

  rownames(mat) <- genes

  if (nrow(mat) < 2 || ncol(mat) < 1) {
    stop("Matrix too small after parsing '", fname, "': ",
         nrow(mat), " genes × ", ncol(mat), " samples")
  }

  list(matrix = mat, type = "matrix", warnings = warnings_out)
}


#' Internal: parse a TAR or TAR.GZ archive
.parse_tar_archive <- function(path, fname) {

  extract_dir <- tempfile(pattern = "geo_tar_")
  dir.create(extract_dir)
  tryCatch(
    utils::untar(path, exdir = extract_dir),
    error = function(e) stop("Failed to untar '", fname, "': ", e$message)
  )

  all_files <- list.files(extract_dir, full.names = TRUE, recursive = TRUE)

  # Detect 10x / scRNA content → warn and return NULL
  is_10x <- any(grepl("barcodes\\.tsv|features\\.tsv|matrix\\.mtx", all_files))
  if (is_10x) {
    return(list(
      matrix   = NULL,
      type     = "scrna_10x",
      warnings = paste0(
        "⚠️ TAR archive contains 10x Genomics / scRNA-seq files. ",
        "markeR is designed for bulk RNA-seq. ",
        "For scRNA-seq data, pre-aggregate to pseudo-bulk and upload manually."
      )
    ))
  }

  # Detect BigWig / BW-only archives (e.g. RAW.tar of BW from ChIP/ATAC/SAR/END-seq)
  is_bw_only <- length(all_files) > 0 &&
    all(grepl("\\.(bw|bigwig|bedgraph|bg|bam|bai|bed)(\\.gz)?$",
              all_files, ignore.case = TRUE))
  if (is_bw_only) {
    return(list(
      matrix   = NULL,
      type     = "bw_tar",
      warnings = paste0(
        "⚠️ TAR archive contains only BigWig/BW/BED files - ",
        "not a bulk RNA-seq expression matrix. Dataset excluded."
      )
    ))
  }

  # Filter for text-like expression files
  text_files <- all_files[grepl("\\.(txt|tsv|csv|xlsx?)(|\\.gz)$",
                                 all_files, ignore.case = TRUE)]

  if (length(text_files) == 0) {
    return(list(
      matrix   = NULL,
      type     = "unknown_tar",
      warnings = "⚠️ TAR archive contains no recognisable expression files."
    ))
  }

  # Score and pick best
  scores <- score_supp_files(text_files)
  ranked <- text_files[order(-scores)]

  matrices <- list()
  all_warns <- character(0)

  for (f in ranked) {
    res <- tryCatch(
      parse_expression_file(f, basename(f)),
      error = function(e) list(matrix = NULL, warnings = e$message)
    )
    if (!is.null(res$matrix)) {
      matrices[[basename(f)]] <- res$matrix
    }
    all_warns <- c(all_warns, res$warnings)
  }

  if (length(matrices) == 0) {
    return(list(matrix = NULL, type = "tar_failed", warnings = all_warns))
  }

  list(matrix = matrices, type = "tar_multi", warnings = all_warns)
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4 - TIDY METADATA EXTRACTOR
# ─────────────────────────────────────────────────────────────────────────────

#' Extract tidy metadata from a GEO ExpressionSet
#'
#' @param obj GEO ExpressionSet
#' @return data.frame: rows = samples, columns = attributes (SampleID always first)
extract_geo_metadata <- function(obj) {

  meta <- tryCatch(Biobase::pData(obj), error = function(e) NULL)
  if (is.null(meta) || nrow(meta) == 0) {
    stop("No metadata found in GEO object.")
  }

  # Ensure SampleID column
  if (!"geo_accession" %in% colnames(meta)) {
    meta$geo_accession <- rownames(meta)
  }

  # Move geo_accession to first column
  meta <- meta[, c("geo_accession",
                   setdiff(colnames(meta), "geo_accession")), drop = FALSE]

  # Clean up column names (remove trailing ":ch1" GEO suffixes)
  colnames(meta) <- gsub(":ch1$", "", colnames(meta))
  colnames(meta) <- trimws(colnames(meta))

  meta
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5 - EXPRESSION–METADATA BRIDGE (FUZZY MATCHING)
# ─────────────────────────────────────────────────────────────────────────────

#' Internal: token-overlap match between expression column names and metadata values
#'
#' Tokenises both sides by splitting on non-alphanumeric characters, then scores
#' each (expr_col, meta_val) pair by the summed length of shared tokens (≥2 chars).
#' Designed to handle cases like expr col "counts.sampA" vs meta value "library: sampA".
#'
#' @param expr_cols  Character vector of expression matrix column names.
#' @param meta_vals  Character vector of values from one metadata column.
#' @param min_score  Minimum token-length score to accept a match (default 3).
#' @return List: idx (best meta index per expr col), scores, accepted (logical), n_matched.
.try_token_match <- function(expr_cols, meta_vals, min_score = 2L) {

  tokenize <- function(x) {
    toks <- unlist(strsplit(as.character(x), "[^a-zA-Z0-9]+"))
    toks <- tolower(toks)
    # Keep tokens of length >= 2 (e.g. "doxo", "control") AND purely numeric
    # tokens of any length. Numeric tokens are frequently the ONLY thing that
    # distinguishes replicates (e.g. "doxo_1" vs "doxo_2") - dropping single
    # digits here would make every replicate of a condition score identically
    # against each other, causing ties in the matcher below.
    keep <- nchar(toks) >= 2L | grepl("^[0-9]+$", toks)
    toks[keep]
  }

  tok_expr <- lapply(expr_cols,               tokenize)
  tok_meta <- lapply(as.character(meta_vals), tokenize)

  n_e <- length(expr_cols)
  n_m <- length(meta_vals)
  score_mat <- matrix(0L, nrow = n_e, ncol = n_m)

  for (i in seq_len(n_e)) {
    te <- tok_expr[[i]]
    if (length(te) == 0L) next
    for (j in seq_len(n_m)) {
      shared <- intersect(te, tok_meta[[j]])
      if (length(shared) > 0L)
        score_mat[i, j] <- sum(nchar(shared))
    }
  }

  # Greedy one-to-one assignment, highest-scoring pairs first. A plain
  # per-row which.max() (the previous approach) lets two different
  # expression columns both claim the SAME metadata row when their scores
  # tie - which silently duplicates that row's values (e.g. "title") across
  # two distinct samples instead of reporting a failed/partial match. Once a
  # metadata row is claimed here, no other expression column may take it.
  idx      <- rep(NA_integer_, n_e)
  accepted <- rep(FALSE, n_e)
  best_scores <- apply(score_mat, 1L, max)

  pairs <- which(score_mat >= min_score, arr.ind = TRUE)
  if (nrow(pairs) > 0L) {
    ord   <- order(score_mat[pairs], decreasing = TRUE)
    pairs <- pairs[ord, , drop = FALSE]
    used_m <- logical(n_m)
    for (k in seq_len(nrow(pairs))) {
      i <- pairs[k, 1L]; j <- pairs[k, 2L]
      if (!is.na(idx[i]) || used_m[j]) next
      idx[i] <- j; used_m[j] <- TRUE; accepted[i] <- TRUE
    }
  }
  # Any expression column left without a unique match (score below
  # min_score, or its best partner(s) were already claimed by a
  # higher-scoring pair) still gets its single best-scoring row as a
  # fallback so downstream code always has an index to work with, but stays
  # flagged as NOT accepted so callers can treat it as an unmatched sample.
  unmatched <- which(is.na(idx))
  for (i in unmatched) {
    idx[i] <- if (any(score_mat[i, ] > 0)) which.max(score_mat[i, ]) else 1L
  }

  list(
    idx       = idx,
    scores    = best_scores,
    accepted  = accepted,
    n_matched = sum(accepted)
  )
}


#' Align expression matrix columns to metadata rows
#'
#' Tries exact match first, then fuzzy string-distance matching.
#' Also tries: strip common prefixes/suffixes, lowercase, alphanumeric-only.
#'
#' @param expr_mat  Numeric matrix: rows = genes, cols = samples.
#' @param meta_df   data.frame: rows = samples, with SampleID column.
#' @param id_col    Column in meta_df containing sample identifiers.
#' @param max_dist  Maximum string distance for fuzzy match (default 5).
#' @return List: `expr` (reordered/renamed matrix), `meta` (aligned metadata),
#'               `match_type` ("exact", "fuzzy", "partial"), `warnings`.
align_expr_to_meta <- function(expr_mat, meta_df,
                                id_col    = "geo_accession",
                                max_dist  = 5,
                                match_col = NULL) {

  warnings_out <- character(0)

  expr_cols  <- colnames(expr_mat)
  meta_ids   <- as.character(meta_df[[id_col]])

  # ── 1. Exact match ─────────────────────────────────────────────────────────
  if (all(expr_cols %in% meta_ids)) {
    meta_aligned <- meta_df[match(expr_cols, meta_ids), , drop = FALSE]
    rownames(meta_aligned) <- expr_cols
    return(list(expr = expr_mat, meta = meta_aligned,
                match_type = "exact", warnings = warnings_out))
  }

  # ── 2. Case-insensitive match ──────────────────────────────────────────────
  expr_lower <- tolower(expr_cols)
  meta_lower <- tolower(meta_ids)
  if (all(expr_lower %in% meta_lower)) {
    idx <- match(expr_lower, meta_lower)
    meta_aligned <- meta_df[idx, , drop = FALSE]
    rownames(meta_aligned) <- expr_cols
    warnings_out <- c(warnings_out,
      "Sample IDs matched case-insensitively between expression matrix and metadata.")
    return(list(expr = expr_mat, meta = meta_aligned,
                match_type = "case_insensitive", warnings = warnings_out))
  }

  # ── 3. Normalised match (alphanumeric only, lowercase) ────────────────────
  normalise <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  expr_norm <- normalise(expr_cols)
  meta_norm <- normalise(meta_ids)

  if (all(expr_norm %in% meta_norm)) {
    idx <- match(expr_norm, meta_norm)
    meta_aligned <- meta_df[idx, , drop = FALSE]
    rownames(meta_aligned) <- expr_cols
    warnings_out <- c(warnings_out,
      "Sample IDs matched after removing non-alphanumeric characters.")
    return(list(expr = expr_mat, meta = meta_aligned,
                match_type = "normalised", warnings = warnings_out))
  }

  # ── 4. Fuzzy (string-distance) match ──────────────────────────────────────
  if (requireNamespace("stringdist", quietly = TRUE)) {

    dist_mat <- stringdist::stringdistmatrix(
      expr_norm, meta_norm,
      method = "lv"   # Levenshtein
    )

    # For each expression column, find the closest metadata ID
    min_dists <- apply(dist_mat, 1, min)
    best_idx  <- apply(dist_mat, 1, which.min)

    # Accept only matches within max_dist
    accepted <- min_dists <= max_dist

    if (sum(accepted) == length(expr_cols)) {
      meta_aligned <- meta_df[best_idx, , drop = FALSE]
      rownames(meta_aligned) <- expr_cols
      warnings_out <- c(warnings_out,
        paste0("Sample IDs fuzzy-matched (Levenshtein, max_dist=", max_dist,
               "). Max distance used: ", max(min_dists), ". ",
               "Please verify the alignment is correct."))
      return(list(expr = expr_mat, meta = meta_aligned,
                  match_type = "fuzzy", warnings = warnings_out))
    }

    # Partial alignment: keep only matched samples
    if (sum(accepted) > 0) {
      expr_sub  <- expr_mat[, accepted, drop = FALSE]
      meta_aligned <- meta_df[best_idx[accepted], , drop = FALSE]
      rownames(meta_aligned) <- expr_cols[accepted]
      warnings_out <- c(warnings_out,
        paste0("Only ", sum(accepted), "/", length(expr_cols),
               " samples could be matched between expression matrix and metadata. ",
               "Unmatched samples removed."))
      return(list(expr = expr_sub, meta = meta_aligned,
                  match_type = "partial_fuzzy", warnings = warnings_out))
    }
  }

  # ── 5. No match found ─────────────────────────────────────────────────────

  # ── 6. User-specified column - token-overlap matching ──────────────────────
  # Invoked when the user selects a specific metadata column via the Shiny UI.
  # Handles mismatches like expr col "counts.sampA" vs meta value "library: sampA".
  if (!is.null(match_col) && match_col %in% colnames(meta_df)) {

    meta_vals <- as.character(meta_df[[match_col]])
    tmatch    <- .try_token_match(expr_cols, meta_vals)

    if (tmatch$n_matched == length(expr_cols)) {
      meta_aligned <- meta_df[tmatch$idx, , drop = FALSE]
      rownames(meta_aligned) <- expr_cols
      warnings_out <- c(warnings_out, paste0(
        "Sample IDs matched via token overlap against metadata column '",
        match_col, "'. Verify the alignment in the preview."))
      return(list(expr       = expr_mat,
                  meta       = meta_aligned,
                  match_type = "token_match",
                  warnings   = warnings_out))
    }

    if (tmatch$n_matched > 0L) {
      accepted     <- tmatch$accepted
      expr_sub     <- expr_mat[, accepted, drop = FALSE]
      meta_aligned <- meta_df[tmatch$idx[accepted], , drop = FALSE]
      rownames(meta_aligned) <- expr_cols[accepted]
      warnings_out <- c(warnings_out, paste0(
        "Partial token match via column '", match_col, "': ",
        tmatch$n_matched, "/", length(expr_cols), " samples matched. ",
        "Unmatched samples were dropped."))
      return(list(expr       = expr_sub,
                  meta       = meta_aligned,
                  match_type = "partial_token_match",
                  warnings   = warnings_out))
    }

    warnings_out <- c(warnings_out,
      paste0("Token match via column '", match_col,
             "' also failed - no shared tokens found."))
  }

  warnings_out <- c(warnings_out,
    "⚠️ Could not align sample IDs between expression matrix and metadata. ",
    "Returning expression matrix with metadata unaligned.")
  list(expr = expr_mat, meta = meta_df,
       match_type = "unmatched", warnings = warnings_out)
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6 - SUPERSERIES HANDLER
# ─────────────────────────────────────────────────────────────────────────────

#' Detect if a GEO accession is a SuperSeries
#'
#' Returns a character vector of sub-series GSE IDs, or character(0) if not.
detect_superseries <- function(gse_accession) {
  tryCatch({
    gse <- GEOquery::getGEO(gse_accession, GSEMatrix = FALSE)
    relations <- GEOquery::Meta(gse)$relation
    subseries <- relations[grepl("SuperSeries of", relations, ignore.case = TRUE)]
    gsm_ids <- gsub(".*?(GSE[0-9]+).*", "\\1", subseries)
    gsm_ids[grepl("^GSE[0-9]+$", gsm_ids)]
  }, error = function(e) character(0))
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7 - SALMON (.sf) HANDLER via tximport
# ─────────────────────────────────────────────────────────────────────────────

#' Import Salmon .sf files using tximport
#'
#' @param sf_files  Named character vector: names = sample IDs, values = paths.
#' @param tx2gene   Optional data.frame with columns (tx_id, gene_id).
#'                  If NULL, transcript-level counts are returned as-is.
#' @return Numeric matrix: rows = genes/transcripts, cols = samples.
import_salmon_sf <- function(sf_files, tx2gene = NULL) {

  if (!requireNamespace("tximport", quietly = TRUE)) {
    stop("Package 'tximport' is required for Salmon .sf files. ",
         "Install with: BiocManager::install('tximport')")
  }

  missing <- !file.exists(sf_files)
  if (any(missing)) {
    stop("Missing Salmon .sf files: ",
         paste(names(sf_files)[missing], collapse = ", "))
  }

  txi <- if (!is.null(tx2gene)) {
    tximport::tximport(sf_files, type = "salmon",
                       tx2gene = tx2gene, ignoreTxVersion = TRUE)
  } else {
    tximport::tximport(sf_files, type = "salmon", txOut = TRUE)
  }

  counts <- round(txi$counts)
  storage.mode(counts) <- "numeric"
  counts
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8 - MAIN ORCHESTRATOR
# ─────────────────────────────────────────────────────────────────────────────

#' Full GEO import pipeline
#'
#' Downloads GEO data, detects assay type, parses expression, extracts metadata,
#' and aligns them. Returns a structured result list.
#'
#' @param gse_accession Character. GEO Series accession (e.g., "GSE334803").
#' @param progress_fn   Optional function(value, detail) for Shiny progress.
#' @return List with elements:
#'   \describe{
#'     \item{status}{"ok", "warning", "excluded", "error"}
#'     \item{expr}{Numeric matrix (genes × samples), or NULL}
#'     \item{meta}{data.frame (samples × attributes), or NULL}
#'     \item{supp_files}{data.frame from getGEOSuppFiles(), or NULL}
#'     \item{geo_objects}{Raw GEO list}
#'     \item{assay_type}{"bulk_rnaseq", "microarray", "scrna", "excluded", "salmon_sf"}
#'     \item{messages}{Character vector of info/warnings/errors}
#'     \item{candidates}{Named list of candidate matrices (if multiple found)}
#'   }
import_geo_data <- function(gse_accession,
                             progress_fn = NULL,
                             tmp_dir     = tempdir()) {

  msgs    <- character(0)
  .prog   <- function(val, detail = "") {
    if (!is.null(progress_fn)) progress_fn(val, detail)
  }

  result <- list(
    status      = "ok",
    expr        = NULL,
    meta        = NULL,
    supp_files  = NULL,
    geo_objects = NULL,
    assay_type  = "unknown",
    messages    = msgs,
    candidates  = NULL
  )

  # ── 1. Download GEO Series Matrix ─────────────────────────────────────────
  .prog(0.05, "Fetching GEO Series Matrix…")
  geo_objs <- tryCatch(
    GEOquery::getGEO(gse_accession, GSEMatrix = TRUE),
    error = function(e) stop("getGEO failed: ", e$message)
  )
  result$geo_objects <- geo_objs
  .prog(0.20, "Classifying assay type…")

  # ── 2. Classify assay type of first object ─────────────────────────────────
  if (length(geo_objs) == 0) stop("No GEO objects returned for ", gse_accession)

  assay_info <- classify_geo_assay(geo_objs[[1]], gse_accession)
  result$assay_type <- assay_info$type

  if (assay_info$type == "excluded") {
    result$status   <- "excluded"
    result$messages <- c(msgs, assay_info$reason)
    return(result)
  }

  # ── 3. Try built-in expression matrix (Series Matrix) ─────────────────────
  .prog(0.30, "Checking Series Matrix expression data…")
  for (obj in geo_objs) {
    expr <- tryCatch(Biobase::exprs(obj), error = function(e) NULL)
    if (!is.null(expr) && nrow(expr) > 0 && ncol(expr) > 0) {
      result$expr <- expr
      result$meta <- extract_geo_metadata(obj)
      break
    }
  }

  # ── 4. Download supplementary files ───────────────────────────────────────
  .prog(0.40, "Downloading supplementary files…")
  supp <- tryCatch(
    GEOquery::getGEOSuppFiles(gse_accession,
                               makeDirectory = FALSE,
                               baseDir       = tmp_dir),
    error = function(e) {
      msgs <<- c(msgs, paste("getGEOSuppFiles warning:", e$message))
      NULL
    }
  )
  result$supp_files <- supp

  # If expression matrix already found, we are done with parsing
  if (!is.null(result$expr)) {
    result$messages <- msgs
    return(result)
  }

  # ── 5. Parse supplementary files ──────────────────────────────────────────
  if (is.null(supp) || nrow(supp) == 0) {
    result$status   <- "warning"
    result$messages <- c(msgs, "No supplementary files found and Series Matrix has no expression data.")
    return(result)
  }

  supp_paths <- rownames(supp)
  scores     <- score_supp_files(supp_paths)
  ranked     <- supp_paths[order(-scores)]

  .prog(0.55, "Parsing supplementary expression files…")

  candidates  <- list()
  parse_warns <- character(0)
  salmon_files <- character(0)

  for (fp in ranked) {
    bname      <- basename(fp)
    file_type  <- classify_supp_file(bname)

    if (file_type == "excluded") {
      msgs <- c(msgs, paste0("Skipped (non-RNA-seq): ", bname))
      next
    }

    if (file_type == "salmon") {
      salmon_files <- c(salmon_files, fp)
      next
    }

    if (file_type == "scrna") {
      msgs <- c(msgs, paste0(
        "⚠️ Skipped (single-cell file): ", bname, ". ",
        "markeR requires bulk RNA-seq. Pre-aggregate to pseudo-bulk before uploading."))
      next
    }

    res <- tryCatch(
      parse_expression_file(fp, bname),
      error = function(e) list(matrix = NULL, type = "parse_error", warnings = e$message)
    )

    parse_warns <- c(parse_warns, res$warnings)

    if (!is.null(res$matrix)) {
      if (is.list(res$matrix)) {
        # Multiple matrices from archive
        candidates <- c(candidates, res$matrix)
      } else {
        candidates[[bname]] <- res$matrix
      }
    }

    res_type <- if (is.null(res$type) || length(res$type) == 0L) "unknown" else res$type
    if (res_type %in% c("scrna_10x", "salmon_sf")) {
      msgs <- c(msgs, res$warnings)
    }
  }

  # ── 6. Salmon .sf handling ─────────────────────────────────────────────────
  if (length(salmon_files) > 0 && length(candidates) == 0) {
    result$assay_type <- "salmon_sf"
    # Name files by GSM ID extracted from filename
    sf_names <- gsub(".*?(GSM[0-9]+).*", "\\1", basename(salmon_files))
    names(salmon_files) <- ifelse(grepl("^GSM", sf_names), sf_names, basename(salmon_files))
    result$messages <- c(msgs, parse_warns,
      paste0("Salmon .sf files detected (", length(salmon_files), " samples). ",
             "Use import_salmon_sf() with a tx2gene table, or provide a pre-aggregated counts matrix."))
    result$status <- "warning"
    return(result)
  }

  if (length(candidates) == 0) {
    result$status   <- "warning"
    result$messages <- c(msgs, parse_warns,
      "No readable expression matrix found in supplementary files.")
    return(result)
  }

  # ── 7. Select best candidate ───────────────────────────────────────────────
  .prog(0.75, "Selecting best expression matrix…")

  if (length(candidates) == 1) {
    result$expr <- candidates[[1]]
  } else {
    # Score candidate matrices by size (more genes × samples = better)
    sizes <- sapply(candidates, function(m) nrow(m) * ncol(m))
    result$expr       <- candidates[[which.max(sizes)]]
    result$candidates <- candidates
    msgs <- c(msgs,
      paste0("Multiple candidate matrices found (", length(candidates),
             "). Largest selected automatically. Others available in `candidates`."))
  }

  # ── 8. Extract metadata if not yet done ────────────────────────────────────
  if (is.null(result$meta) && length(geo_objs) > 0) {
    result$meta <- tryCatch(
      extract_geo_metadata(geo_objs[[1]]),
      error = function(e) { msgs <<- c(msgs, paste("Metadata extraction failed:", e$message)); NULL }
    )
  }

  # ── 9. Align expression ↔ metadata ────────────────────────────────────────
  .prog(0.88, "Aligning expression matrix to metadata…")

  if (!is.null(result$expr) && !is.null(result$meta)) {
    aligned <- tryCatch(
      align_expr_to_meta(result$expr, result$meta),
      error = function(e) {
        msgs <<- c(msgs, paste("Alignment failed:", e$message))
        list(expr = result$expr, meta = result$meta,
             match_type = "failed", warnings = e$message)
      }
    )
    result$expr <- aligned$expr
    result$meta <- .ensure_sampleid_col(aligned$meta)
    msgs <- c(msgs, aligned$warnings,
              paste0("Sample ID alignment: ", aligned$match_type))
  }

  .prog(1.0, "Done.")
  result$messages <- c(msgs, parse_warns)
  result
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8b - SAMPLEID COLUMN HELPER
# ─────────────────────────────────────────────────────────────────────────────

#' Ensure metadata has a SampleID first column derived from rownames
#'
#' After any alignment step the matched sample names live in rownames(meta).
#' This helper promotes them into an explicit first column called "SampleID",
#' preserving every other column (including the original geo_accession column).
#'
#' @param meta data.frame returned by align_expr_to_meta or extract_geo_metadata.
#' @return data.frame with SampleID as the first column.
.ensure_sampleid_col <- function(meta) {
  if (is.null(meta) || nrow(meta) == 0L) return(meta)

  rn <- rownames(meta)

  # If rownames are just 1,2,3,... (no meaningful names) → nothing to promote
  if (is.null(rn) || all(rn == as.character(seq_len(nrow(meta))))) return(meta)

  # Drop any pre-existing SampleID column (will be rebuilt from rownames)
  other_cols <- setdiff(colnames(meta), "SampleID")
  meta_clean <- meta[, other_cols, drop = FALSE]

  # Use data.frame() - cbind() does not accept stringsAsFactors and would
  # create a spurious column named "stringsAsFactors" instead of ignoring it.
  result <- data.frame(SampleID = rn, meta_clean,
                       check.names = FALSE, stringsAsFactors = FALSE)
  rownames(result) <- rn
  result
}


# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9 - SHINY MODULE
# ─────────────────────────────────────────────────────────────────────────────

#' GEO Import Module - UI
#'
#' @param id Shiny module namespace ID.
#' @export
geoImportUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(

    shiny::div(
      class = "alert alert-warning", style = "font-size:0.85em; padding:8px 12px;",
      shiny::icon("flask"),
      shiny::tags$b(" GEO import is in beta."),
      " GEO submissions aren't standardised, so some accessions may not parse or",
      " import correctly. We apologise for any inconvenience.",
      " If an accession doesn't work as expected, please ",
      shiny::tags$a(href = "https://github.com/DiseaseTranscriptomicsLab/markeR/issues",
                    target = "_blank", rel = "noopener noreferrer",
                    "open an issue on markeR's GitHub page"),
      " and we'll look into it."
    ),

    shiny::textInput(
      ns("accession"),
      "Enter GEO Accession (e.g., GSE334803)",
      placeholder = "GSE334803"
    ),

    shiny::actionButton(
      ns("fetch"),
      label = shiny::tagList(
        shiny::icon("download"), "Fetch from GEO"
      ),
      class = "btn-primary"
    ),

    shiny::hr(),

    # GEO dataset info card (title, organism, abstract popover)
    shiny::uiOutput(ns("geo_info_ui")),

    # Status / messages area
    shiny::uiOutput(ns("status_ui")),

    # GEO object selector (multi-platform series)
    shiny::uiOutput(ns("object_selector")),

    # SuperSeries sub-series selector
    shiny::uiOutput(ns("superseries_selector")),

    # Supplementary file selector
    shiny::uiOutput(ns("supp_selector")),

    # Candidate matrix selector (after auto-parse)
    shiny::uiOutput(ns("candidate_selector")),

    # Metadata column configuration
    shiny::uiOutput(ns("meta_col_selector")),

    # Sample alignment panel - mismatch resolution + optional rename
    shiny::uiOutput(ns("sample_alignment_ui")),

    shiny::hr()
  )
}


#' GEO Import Module - Server
#'
#' @param id         Shiny module namespace ID.
#' @param session    Parent Shiny session (for notifications).
#' @return A list of reactives: `expr_data`, `meta_data`, `status`.
#' @export
geoImportServer <- function(id, log_fn = NULL) {

  shiny::moduleServer(id, function(input, output, session) {

    ns <- session$ns

    # Internal logging helper - calls log_fn (app-level) if provided.
    # isolate() prevents taking reactive dependencies inside observers.
    glog <- function(msg) {
      if (!is.null(log_fn)) log_fn(msg)
    }

    # ── Reactive storage ────────────────────────────────────────────────────
    rv <- shiny::reactiveValues(
      geo_result       = NULL,   # full import_geo_data() result
      expr             = NULL,   # final numeric matrix
      meta             = NULL,   # final metadata data.frame
      messages         = character(0),
      status           = "idle", # "idle","loading","ok","warning","excluded","error"
      subseries        = character(0),  # detected sub-series for SuperSeries
      match_col_needed = FALSE,  # TRUE when alignment failed → show column-picker UI
      example_expr_col = "",     # one unmatched expr col name to show in the hint
      geo_info         = NULL    # list(accession, title, abstract, organism, n_samples)
    )

    # Internal helper: check alignment after expr/meta update and set flag
    .check_alignment <- function() {
      e <- rv$expr
      m <- rv$meta
      if (is.null(e) || is.null(m)) { rv$match_col_needed <- FALSE; return() }
      # Aligned means every expr column appears as a rowname in meta
      all_matched <- all(colnames(e) %in% rownames(m))
      rv$match_col_needed <- !all_matched
      if (!all_matched) {
        unmatched <- setdiff(colnames(e), rownames(m))
        rv$example_expr_col <- if (length(unmatched) > 0) unmatched[1] else colnames(e)[1]
      }
    }

    # ── 1. Fetch button ─────────────────────────────────────────────────────
    shiny::observeEvent(input$fetch, {
      shiny::req(input$accession)
      accession <- trimws(toupper(input$accession))

      if (!grepl("^GSE[0-9]+$", accession)) {
        shiny::showNotification(
          "Invalid accession format. Expected 'GSExxxxxx'.",
          type = "error", duration = 5
        )
        glog(paste0("Invalid accession entered: '", accession, "' - expected GSExxxxxx format."))
        return()
      }

      glog(paste0("=== GEO import started: ", accession, " ==="))
      rv$status   <- "loading"
      rv$messages <- character(0)
      rv$expr     <- NULL
      rv$meta     <- NULL
      rv$geo_result <- NULL

      # Centred full-page loading overlay
      shiny::showModal(shiny::modalDialog(
        title     = NULL,
        footer    = NULL,
        easyClose = FALSE,
        shiny::tags$div(
          style = "text-align:center; padding:40px 20px;",
          shiny::tags$div(
            class = "spinner-border text-primary",
            style = "width:3.5rem; height:3.5rem;",
            role  = "status",
            shiny::tags$span(class = "visually-hidden", "Loading…")
          ),
          shiny::tags$h5(
            paste("Fetching", accession, "from GEO…"),
            style = "margin-top:18px; color:#333;"
          ),
          shiny::tags$p(
            "Downloading metadata & supplementary files - this usually takes 30–90 s.",
            style = "color:#777; font-size:0.9em; margin-top:6px;"
          )
        )
      ))

      shiny::withProgress(
        message = paste("Importing", accession, "from GEO…"),
        value   = 0,
        {
          result <- tryCatch(
            import_geo_data(
              gse_accession = accession,
              progress_fn   = function(v, d) shiny::incProgress(v, detail = d),
              tmp_dir       = tempdir()
            ),
            error = function(e) {
              list(status = "error", expr = NULL, meta = NULL,
                   messages = paste("Import error:", e$message),
                   geo_objects = NULL, supp_files = NULL,
                   assay_type = "unknown", candidates = NULL)
            }
          )

          rv$geo_result <- result
          rv$status     <- result$status
          rv$messages   <- result$messages
          rv$expr       <- result$expr
          rv$meta       <- .ensure_sampleid_col(result$meta)
          .check_alignment()

          # ── Pipe all pipeline messages to the app log ──────────────────────
          for (m in result$messages) glog(m)
          if (!is.null(result$expr))
            glog(paste0("Expression matrix: ", nrow(result$expr), " genes × ",
                        ncol(result$expr), " samples."))
          if (!is.null(result$meta))
            glog(paste0("Metadata: ", nrow(result$meta), " rows × ",
                        ncol(result$meta), " variables."))
          glog(paste0("Import status: ", result$status))

          # ── Extract dataset summary for the info card ──────────────────
          if (!is.null(result$geo_objects) && length(result$geo_objects) > 0) {
            obj0 <- result$geo_objects[[1]]
            rv$geo_info <- list(
              accession = accession,
              title     = tryCatch(
                Biobase::experimentData(obj0)@title, error = function(e) ""),
              abstract  = tryCatch(
                Biobase::experimentData(obj0)@abstract, error = function(e) ""),
              organism  = tryCatch(
                paste(unique(na.omit(Biobase::pData(obj0)[["organism_ch1"]])),
                      collapse = ", "),
                error = function(e) ""),
              n_samples = tryCatch(nrow(Biobase::pData(obj0)), error = function(e) 0L)
            )
          }

          # Detect SuperSeries
          if (!is.null(result$geo_objects) && length(result$geo_objects) > 0) {
            subs <- tryCatch(
              detect_superseries(accession),
              error = function(e) character(0)
            )
            rv$subseries <- subs
          }
        }
      )

      shiny::removeModal()

      # Build a clean match-summary string from the aligned data.
      # Internal pipeline messages (like "Returning expression matrix with
      # metadata unaligned") are intentionally not shown to the user.
      .fmt_match <- function() {
        e <- rv$expr; m <- rv$meta
        if (is.null(e) || is.null(m)) return(NULL)
        n_expr  <- ncol(e)
        # Check against SampleID column, first column, and rownames - same
        # logic as align_meta_to_expr in the main app
        id_col   <- if ("SampleID" %in% colnames(m)) "SampleID" else colnames(m)[1]
        meta_ids <- unique(c(as.character(m[[id_col]]), rownames(m)))
        matched  <- sum(colnames(e) %in% meta_ids)
        dropped  <- n_expr - matched
        msg <- paste0(matched, "/", n_expr, " sample",
                      if (n_expr != 1) "s" else "", " matched")
        if (dropped > 0)
          msg <- paste0(msg, " - ", dropped, " dropped (no metadata)")
        msg
      }
      match_summary <- .fmt_match()
      # Log import summary (counts only - no alignment performed yet)
      if (!is.null(rv$expr) && !is.null(rv$meta)) {
        glog(paste0("GEO import complete: ",
                    ncol(rv$expr), " expression samples, ",
                    nrow(rv$meta), " metadata rows loaded. ",
                    "Use the Data tab to review and align before proceeding."))
      }

      # Show result notification
      import_counts <- if (!is.null(rv$expr) && !is.null(rv$meta))
        paste0(ncol(rv$expr), " expression samples, ", nrow(rv$meta), " metadata rows")
      else match_summary
      if (rv$status == "ok") {
        shiny::showNotification(
          paste0(accession, " imported - ", import_counts),
          type = "default", duration = 6
        )
      } else if (rv$status == "excluded") {
        shiny::showNotification(
          shiny::HTML(paste0(
            "⛔ <b>", accession, " excluded:</b><br>",
            paste(rv$messages, collapse = "<br>")
          )),
          type = "error", duration = NULL
        )
      } else if (rv$status == "warning") {
        # Show match summary + only user-relevant warning lines
        user_warns <- grep(
          "format|file|parse|GSM|column|TAR|gz|supplementary",
          rv$messages, value = TRUE, ignore.case = TRUE
        )
        body <- if (!is.null(import_counts)) import_counts else ""
        if (length(user_warns) > 0)
          body <- paste(c(body, user_warns), collapse = "<br>")
        shiny::showNotification(
          shiny::HTML(paste0(
            "⚠️ <b>", accession, " loaded with warnings:</b><br>", body
          )),
          type = "warning", duration = 12
        )
      } else if (rv$status == "error") {
        shiny::showNotification(
          shiny::HTML(paste0(
            "❌ <b>Import failed for ", accession, ":</b><br>",
            paste(rv$messages, collapse = "<br>")
          )),
          type = "error", duration = NULL
        )
      }
    })

    # ── 2. Status UI ────────────────────────────────────────────────────────
    output$status_ui <- shiny::renderUI({
      if (rv$status == "idle")    return(NULL)
      if (rv$status == "loading") return(shiny::p("⏳ Loading…"))

      msgs <- rv$messages
      if (length(msgs) == 0) return(NULL)

      colour <- switch(rv$status,
        "ok"       = "#155724",
        "warning"  = "#856404",
        "excluded" = "#721c24",
        "error"    = "#721c24",
        "#333"
      )
      bg <- switch(rv$status,
        "ok"       = "#d4edda",
        "warning"  = "#fff3cd",
        "excluded" = "#f8d7da",
        "error"    = "#f8d7da",
        "#f8f9fa"
      )

      shiny::tags$div(
        style = paste0("background:", bg, "; color:", colour,
                       "; border-radius:4px; padding:10px; margin-bottom:10px;"),
        shiny::tags$ul(
          lapply(msgs, function(m) shiny::tags$li(m))
        )
      )
    })

    # ── 3. GEO object (platform) selector ───────────────────────────────────
    output$object_selector <- shiny::renderUI({
      res <- rv$geo_result
      if (is.null(res) || is.null(res$geo_objects)) return(NULL)
      objs <- res$geo_objects
      if (length(objs) <= 1) return(NULL)

      labels <- sapply(objs, function(obj) {
        title    <- tryCatch(Biobase::experimentData(obj)@title, error = function(e) "")
        platform <- tryCatch(Biobase::annotation(obj), error = function(e) "")
        paste0(title, " [", platform, "]")
      })

      shiny::selectInput(
        session$ns("selected_object"),
        "Select GEO platform/object:",
        choices  = stats::setNames(seq_along(objs), labels),
        selected = 1
      )
    })

    # Reload expr + meta when user picks a different GEO object
    shiny::observeEvent(input$selected_object, {
      res <- rv$geo_result
      shiny::req(res, res$geo_objects)
      idx <- as.integer(input$selected_object)
      obj <- res$geo_objects[[idx]]

      expr <- tryCatch(Biobase::exprs(obj), error = function(e) NULL)
      meta <- tryCatch(extract_geo_metadata(obj), error = function(e) NULL)

      if (!is.null(expr) && nrow(expr) > 0) rv$expr <- expr
      if (!is.null(meta)) rv$meta <- .ensure_sampleid_col(meta)
      .check_alignment()

      glog(paste0("Platform/object selected: #", idx,
                  if (!is.null(expr) && nrow(expr) > 0)
                    paste0(" - ", nrow(expr), " genes × ", ncol(expr), " samples")
                  else " - no expression matrix"))
    })

    # ── 4. SuperSeries selector ─────────────────────────────────────────────
    output$superseries_selector <- shiny::renderUI({
      if (length(rv$subseries) == 0) return(NULL)
      shiny::tagList(
        shiny::helpText(shiny::HTML(paste0(
          "⚠️ <b>", input$accession, "</b> is a SuperSeries containing ",
          length(rv$subseries), " sub-series. ",
          "Select one to import, or import the full series above."
        ))),
        shiny::selectInput(
          session$ns("selected_subseries"),
          "Import a sub-series instead:",
          choices = c("(Use full series)" = "", rv$subseries)
        )
      )
    })

    shiny::observeEvent(input$selected_subseries, {
      shiny::req(input$selected_subseries != "")
      # Update text input and trigger fetch
      shiny::updateTextInput(session, "accession",
                              value = input$selected_subseries)
      shiny::showNotification(
        paste("Click 'Fetch from GEO' to import", input$selected_subseries),
        type = "message", duration = 5
      )
    })

    # ── 5. Supplementary file selector ──────────────────────────────────────
    output$supp_selector <- shiny::renderUI({
      res <- rv$geo_result
      if (is.null(res) || is.null(res$supp_files)) return(NULL)
      supp  <- res$supp_files
      if (nrow(supp) == 0) return(NULL)

      paths  <- rownames(supp)
      scores <- score_supp_files(paths)
      labels <- paste0(basename(paths),
                       " (score: ", round(scores, 2), ")")

      score_uid <- paste0("score_info_", gsub("[^a-zA-Z0-9]", "", paste0(sample(letters, 6), collapse="")))

      shiny::tagList(
        shiny::h5("Manual supplementary file selection"),

        # Score help row: text + ⓘ popover
        shiny::tags$div(
          style = "display:flex;align-items:center;gap:6px;margin-bottom:6px;",
          shiny::tags$span(
            style = "font-size:0.85em;color:#6c757d;",
            "Files are ranked by expression-likelihood score (0–1)."
          ),
          shiny::tags$span(
            id        = score_uid,
            style     = "cursor:pointer;color:#4f46e5;font-size:1em;",
            tabindex  = "0",
            `data-bs-toggle`    = "popover",
            `data-bs-trigger`   = "hover focus",
            `data-bs-placement` = "right",
            `data-bs-html`      = "true",
            `data-bs-title`     = "How scores are calculated",
            `data-bs-content`   = paste0(
              "<b>Base score</b> from file extension:<br>",
              "&nbsp;• .tsv / .csv → 0.55<br>",
              "&nbsp;• .txt → 0.50<br>",
              "&nbsp;• .xls / .xlsx → 0.60<br>",
              "&nbsp;• .gz → refined from inner extension<br>",
              "&nbsp;• .tar / .zip → 0.20<br>",
              "<br><b>+Bonus</b> for expression keywords in filename ",
              "(e.g. <i>count, tpm, fpkm, expr, normalized</i>): +0.10 each<br>",
              "<br><b>−Penalty</b> for non-expression keywords ",
              "(e.g. <i>clinical, metadata, manifest</i>): −0.20 each<br>",
              "<br><b>−Penalty</b> for other-assay signals ",
              "(e.g. <i>chip, atac, hic, methylat</i>): −0.40–1.00<br>",
              "<br>Final score is clamped to [0, 1]. ",
              "Higher → more likely to be a bulk RNA-seq expression matrix."
            ),
            shiny::icon("circle-info")
          )
        ),

        shiny::tags$script(shiny::HTML(paste0(
          "(function(){\n",
          "  var tmpl='<div class=\"popover\" role=\"tooltip\">",
          "<div class=\"popover-arrow\"></div>",
          "<h3 class=\"popover-header\"></h3>",
          "<div class=\"popover-body\" style=\"max-height:320px;overflow-y:auto;",
          "font-size:0.82em;line-height:1.5;\"></div></div>';\n",
          "  function init(){\n",
          "    var el=document.getElementById('", score_uid, "');\n",
          "    if(el&&window.bootstrap&&window.bootstrap.Popover){\n",
          "      new window.bootstrap.Popover(el,{sanitize:false,container:'body',template:tmpl});\n",
          "    }\n",
          "  }\n",
          "  setTimeout(init,300);\n",
          "  setTimeout(init,900);\n",
          "})();"
        ))),

        shiny::selectInput(
          session$ns("selected_supp"),
          "Choose a supplementary file:",
          choices  = c("(Auto-selected)" = "",
                       stats::setNames(paths, labels)[order(-scores)]),
          selected = ""
        )
      )
    })

    shiny::observeEvent(input$selected_supp, {
      shiny::req(input$selected_supp != "")
      fp  <- input$selected_supp

      glog(paste0("Manual supplementary file selected: ", basename(fp)))

      shiny::withProgress(message = "Parsing supplementary file…", value = 0, {
        shiny::incProgress(0.3, detail = "Reading file…")
        res <- tryCatch(
          parse_expression_file(fp, basename(fp)),
          error = function(e) list(matrix = NULL, warnings = e$message)
        )
        shiny::incProgress(0.7, detail = "Validating…")

        if (!is.null(res$matrix)) {
          mat <- if (is.list(res$matrix)) res$matrix[[1]] else res$matrix
          rv$expr <- mat
          glog(paste0("Supplementary file parsed: ", nrow(mat), " genes × ", ncol(mat), " samples."))
          shiny::showNotification("Expression data loaded from selected file.",
                                   type = "default", duration = 4)
        } else {
          glog(paste0("Supplementary file parse failed: ",
                      paste(res$warnings, collapse = "; ")))
          shiny::showNotification(
            shiny::HTML(paste0("❌ Could not parse file:<br>",
                               paste(res$warnings, collapse = "<br>"))),
            type = "error", duration = NULL
          )
        }
      })
    })

    # ── 6. Candidate matrix selector ────────────────────────────────────────
    output$candidate_selector <- shiny::renderUI({
      res <- rv$geo_result
      if (is.null(res) || is.null(res$candidates)) return(NULL)
      cands <- res$candidates

      sizes  <- sapply(cands, function(m) paste0(nrow(m), " genes × ", ncol(m), " samples"))
      labels <- paste0(names(cands), "  [", sizes, "]")

      shiny::tagList(
        shiny::h5("Multiple expression matrices found"),
        shiny::helpText("The largest was selected automatically. Choose another if needed."),
        shiny::selectInput(
          session$ns("selected_candidate"),
          "Select matrix:",
          choices = stats::setNames(seq_along(cands), labels),
          selected = which.max(sapply(cands, function(m) nrow(m) * ncol(m)))
        )
      )
    })

    shiny::observeEvent(input$selected_candidate, {
      res <- rv$geo_result
      shiny::req(res, res$candidates)
      idx <- as.integer(input$selected_candidate)
      rv$expr <- res$candidates[[idx]]
      nm  <- names(res$candidates)[idx]
      m   <- res$candidates[[idx]]
      glog(paste0("Expression matrix selected: '", nm, "' - ",
                  nrow(m), " genes × ", ncol(m), " samples."))
    })

    # ── 7. GEO dataset info card ────────────────────────────────────────────
    output$geo_info_ui <- shiny::renderUI({
      info <- rv$geo_info
      if (is.null(info)) return(NULL)

      uid <- paste0("geo_abs_", gsub("[^a-zA-Z0-9]", "", info$accession))

      # Escape content for Bootstrap data attribute (double-quotes break it)
      esc <- function(x) {
        x <- gsub("&",  "&amp;",  as.character(x))
        x <- gsub("<",  "&lt;",   x)
        x <- gsub(">",  "&gt;",   x)
        x <- gsub('"',  "&quot;", x)
        x <- gsub("'",  "&#39;",  x)
        x
      }

      abs_html <- esc(info$abstract)
      title_s  <- esc(info$title)

      shiny::tagList(

        # Popover CSS - scrollable body, reasonable max-width
        shiny::tags$style(shiny::HTML(
          ".geo-abs-pop .popover { max-width:420px; }
           .geo-abs-pop .popover-body {
             max-height:280px; overflow-y:auto;
             font-size:0.83em; line-height:1.5; white-space:normal;
           }"
        )),

        shiny::tags$div(
          style = paste0(
            "background:#eef2ff;border:1px solid #c7d2fe;border-radius:6px;",
            "padding:9px 12px;margin-bottom:8px;"
          ),

          # Top row: accession + external link + abstract popover icon
          shiny::tags$div(
            style = "display:flex;align-items:center;gap:8px;margin-bottom:4px;",
            shiny::tags$strong(
              style = "font-size:0.9em;",
              info$accession
            ),
            shiny::tags$a(
              href   = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
                              info$accession),
              target = "_blank",
              style  = "font-size:0.8em;color:#4f46e5;",
              shiny::icon("arrow-up-right-from-square")
            ),
            # Abstract hover icon - Bootstrap 5 popover
            shiny::tags$span(
              id                = uid,
              style             = "cursor:pointer;color:#4f46e5;margin-left:auto;",
              tabindex          = "0",
              `data-bs-toggle`  = "popover",
              `data-bs-trigger` = "hover focus",
              `data-bs-placement` = "bottom",
              `data-bs-html`    = "true",
              `data-bs-title`   = "Abstract",
              `data-bs-content` = abs_html,
              `data-bs-custom-class` = "geo-abs-pop",
              shiny::icon("circle-info")
            )
          ),

          # Title (truncated)
          shiny::tags$div(
            style = paste0(
              "font-size:0.82em;color:#1e1b4b;",
              "display:-webkit-box;-webkit-line-clamp:2;",
              "-webkit-box-orient:vertical;overflow:hidden;"
            ),
            info$title
          ),

          # Organism + sample count chips
          shiny::tags$div(
            style = "margin-top:5px;display:flex;flex-wrap:wrap;gap:5px;",
            if (nchar(info$organism) > 0)
              shiny::tags$span(
                style = paste0(
                  "font-size:0.75em;background:#c7d2fe;color:#3730a3;",
                  "border-radius:99px;padding:1px 8px;"
                ),
                info$organism
              ),
            shiny::tags$span(
              style = paste0(
                "font-size:0.75em;background:#ddd6fe;color:#5b21b6;",
                "border-radius:99px;padding:1px 8px;"
              ),
              paste0(info$n_samples, " samples")
            )
          )
        ),

        # Init Bootstrap popover via JS (safe: retries until element exists).
        # Custom template puts overflow-y:auto directly on .popover-body so it
        # cannot be overridden by Bootstrap's own stylesheet rules.
        shiny::tags$script(shiny::HTML(paste0(
          "(function() {",
          "  var tmpl = '<div class=\"popover\" role=\"tooltip\">",
          "<div class=\"popover-arrow\"></div>",
          "<h3 class=\"popover-header\"></h3>",
          "<div class=\"popover-body\" style=\"max-height:260px;overflow-y:auto;",
          "font-size:0.83em;line-height:1.5;white-space:normal;\"></div></div>';",
          "  function init() {",
          "    var el = document.getElementById('", uid, "');",
          "    if (el && window.bootstrap && window.bootstrap.Popover) {",
          "      new window.bootstrap.Popover(el, {",
          "        sanitize: false, container: 'body', template: tmpl",
          "      });",
          "    }",
          "  }",
          "  setTimeout(init, 200);",
          "  setTimeout(init, 800);",
          "})();"
        )))
      )
    })


    # ── 8. Metadata column selector ─────────────────────────────────────────
    output$meta_col_selector <- shiny::renderUI({
      meta <- rv$meta
      if (is.null(meta)) return(NULL)

      all_cols <- colnames(meta)
      # SampleID is always kept; let the user pick which additional cols to include
      pick_cols <- setdiff(all_cols, "SampleID")

      shiny::tagList(
        shiny::h5("Metadata columns to keep"),
        shinyWidgets::pickerInput(
          inputId  = session$ns("keep_cols"),
          label    = NULL,
          choices  = pick_cols,
          selected = pick_cols,
          multiple = TRUE,
          options  = list(`actions-box` = TRUE)
        ),
        shiny::actionButton(session$ns("apply_meta"), "Apply",
                             class = "btn-sm btn-secondary")
      )
    })

    shiny::observeEvent(input$apply_meta, {
      shiny::req(rv$meta)
      meta <- rv$meta

      # Always preserve SampleID (first col) + user-selected additional cols
      keep <- unique(c("SampleID", input$keep_cols))
      keep <- intersect(keep, colnames(meta))
      rv$meta <- meta[, keep, drop = FALSE]
      glog(paste0("Metadata columns kept (", length(keep), "): ",
                  paste(keep, collapse = ", "), "."))

      # Re-align if expression also present
      if (!is.null(rv$expr)) {
        aligned <- tryCatch(
          align_expr_to_meta(rv$expr, rv$meta, id_col = "SampleID"),
          error = function(e) list(expr = rv$expr, meta = rv$meta,
                                   match_type = "failed", warnings = e$message)
        )
        rv$expr <- aligned$expr
        rv$meta <- .ensure_sampleid_col(aligned$meta)
        if (length(aligned$warnings) > 0) {
          rv$messages <- c(rv$messages, aligned$warnings)
          for (w in aligned$warnings) glog(w)
        }
        glog(paste0("Re-alignment after column filter: ", aligned$match_type, "."))
      }
      .check_alignment()
    })

    # ── 8. Unified sample-alignment panel ────────────────────────────────────
    #
    # Two modes:
    #  A) Mismatch: amber header + all three fix options shown with descriptions
    #  B) Aligned:  compact "Edit sample names" section for optional renaming
    #
    output$sample_alignment_ui <- shiny::renderUI({
      expr <- rv$expr
      meta <- rv$meta
      if (is.null(expr) && is.null(meta)) return(NULL)

      n_expr <- if (!is.null(expr)) ncol(expr) else 0L
      n_meta <- if (!is.null(meta)) nrow(meta) else 0L
      expr_nms <- if (!is.null(expr)) colnames(expr) else character(0)
      meta_ids  <- if (!is.null(meta)) {
        if ("SampleID" %in% colnames(meta)) meta$SampleID else rownames(meta)
      } else character(0)
      if (length(meta_ids) == 0) meta_ids <- as.character(seq_len(n_meta))
      all_meta_cols <- if (!is.null(meta)) colnames(meta) else character(0)

      # ── helper: build a textarea block ──────────────────────────────────────
      make_textarea <- function(input_id, vals, n) {
        shiny::textAreaInput(
          session$ns(input_id),
          label  = NULL,
          value  = paste(vals, collapse = "\n"),
          rows   = min(6, n),
          resize = "vertical",
          width  = "100%"
        )
      }

      # ── helper: build the remove-samples pickers ─────────────────────────────
      # Separate pickers for expression columns and metadata rows so the user
      # can remove independently from each (e.g., drop a QC sample from expr
      # only, or drop an extra metadata row with no matching expression data).
      remove_samples_ui <- function() {
        if (n_expr == 0 && n_meta == 0) return(NULL)

        make_picker <- function(input_id, names, label_text) {
          if (length(names) == 0) return(NULL)
          shiny::tagList(
            shiny::tags$p(
              style = "font-size:0.82em;font-weight:600;margin:0 0 2px;",
              label_text
            ),
            shiny::tags$p(
              style = "font-size:0.78em;color:#666;margin:0 0 5px;",
              shiny::HTML(paste0(
                "All selected = kept. <b>Deselect</b> to remove."
              ))
            ),
            shinyWidgets::pickerInput(
              session$ns(input_id),
              label    = NULL,
              choices  = names,
              selected = names,
              multiple = TRUE,
              width    = "100%",
              options  = shinyWidgets::pickerOptions(
                actionsBox          = TRUE,
                liveSearch          = TRUE,
                selectedTextFormat  = "count > 3",
                countSelectedText   = "{0} of {1} kept"
              )
            )
          )
        }

        shiny::tagList(
          make_picker("keep_expr_samples",
                      expr_nms,
                      paste0("Expression columns (", n_expr, " samples)")),
          if (n_expr > 0 && n_meta > 0)
            shiny::tags$div(style = "height:6px;"),
          make_picker("keep_meta_samples",
                      meta_ids,
                      paste0("Metadata rows (", n_meta, " samples)")),
          shiny::tags$div(
            style = "margin-top:8px;",
            shiny::actionButton(
              session$ns("apply_remove_samples"),
              shiny::tagList(shiny::icon("trash-can"), "Remove deselected"),
              class = "btn-sm btn-danger"
            )
          )
        )
      }

      # ── OPTION CARD helper ───────────────────────────────────────────────────
      # Creates a bordered card with a numbered badge, title, subtitle, and body.
      # NOTE: overflow must NOT be hidden - Bootstrap dropdown menus need to
      # escape the card boundaries (they use position:absolute).
      opt_card <- function(num, title, subtitle, body, border_color = "#dee2e6") {
        shiny::tags$div(
          style = paste0(
            "border:1px solid ", border_color, ";border-radius:6px;",
            "margin-bottom:10px;"
          ),
          # Header row
          shiny::tags$div(
            style = paste0(
              "display:flex;align-items:center;gap:8px;",
              "padding:8px 10px;background:#f8f9fa;",
              "border-bottom:1px solid ", border_color, ";"
            ),
            shiny::tags$span(
              style = paste0(
                "display:inline-flex;align-items:center;justify-content:center;",
                "width:22px;height:22px;border-radius:50%;",
                "background:#495057;color:#fff;font-size:0.75em;font-weight:700;",
                "flex-shrink:0;"
              ),
              num
            ),
            shiny::tags$div(
              shiny::tags$div(
                style = "font-weight:600;font-size:0.9em;line-height:1.2;",
                title
              ),
              shiny::tags$div(
                style = "font-size:0.8em;color:#666;margin-top:2px;",
                subtitle
              )
            )
          ),
          # Body
          shiny::tags$div(
            style = "padding:10px;",
            body
          )
        )
      }

      # ────────────────────────────────────────────────────────────────────────
      # MODE A: mismatch - show all three options
      # ────────────────────────────────────────────────────────────────────────
      if (rv$match_col_needed) {
        ex_col <- rv$example_expr_col

        shiny::tagList(
          # Amber header
          shiny::tags$div(
            style = paste0(
              "background:#fff3cd;color:#856404;border-radius:6px;",
              "padding:10px 12px;margin-bottom:10px;"
            ),
            shiny::tags$strong("⚠️ Sample names could not be matched automatically"),
            shiny::tags$p(
              style = "margin:5px 0 0;font-size:0.85em;",
              shiny::HTML(paste0(
                "The expression matrix columns and metadata rows use different identifiers.",
                if (nchar(ex_col) > 0)
                  paste0(" Example unmatched column: <code>", ex_col, "</code>.")
                else "",
                "<br>Choose one of the options below to resolve this."
              ))
            )
          ),

          # Option 1 - Remove samples
          opt_card(
            "1",
            "Remove samples",
            paste0("Deselect the samples you want to drop. ",
                   "Expression columns and metadata rows are listed separately - ",
                   "you can remove from one or both independently."),
            remove_samples_ui(),
            border_color = "#f5c6cb"
          ),

          # Option 2 - Token / keyword match
          opt_card(
            "2",
            "Match by keyword in a metadata column",
            paste0("Use when a metadata column contains text that partially overlaps ",
                   "with the expression column names ",
                   "(e.g. 'Library: sampleA' vs 'counts.sampleA')."),
            shiny::tagList(
              shiny::selectInput(
                session$ns("match_col"),
                "Select metadata column to search:",
                choices  = all_meta_cols,
                selected = all_meta_cols[1]
              ),
              shiny::actionButton(
                session$ns("try_match_col"),
                shiny::tagList(shiny::icon("magnifying-glass"), "Try keyword match"),
                class = "btn-warning btn-sm"
              )
            ),
            border_color = "#ffc107"
          ),

          # Option 3 - Rename expression columns
          opt_card(
            "3",
            "Rename expression columns to match metadata",
            paste0("Use when you know the correct sample IDs. ",
                   "Edit the names below so they match the metadata exactly."),
            shiny::tagList(
              make_textarea("rename_expr_names", expr_nms, n_expr),
              shiny::actionButton(
                session$ns("apply_rename_expr"),
                shiny::tagList(shiny::icon("check"), "Apply & retry"),
                class = "btn-sm btn-secondary"
              )
            )
          ),

          # Option 4 - Set metadata sample IDs
          opt_card(
            "4",
            "Assign sample names to metadata rows",
            paste0("Use when the metadata has no usable sample IDs. ",
                   "Type the correct name for each row in order."),
            shiny::tagList(
              make_textarea("override_meta_names", meta_ids, n_meta),
              shiny::actionButton(
                session$ns("apply_override_meta"),
                shiny::tagList(shiny::icon("check"), "Apply & retry"),
                class = "btn-sm btn-secondary"
              )
            )
          )
        )

      # ────────────────────────────────────────────────────────────────────────
      # MODE B: aligned (or data not fully loaded) - compact edit section
      # ────────────────────────────────────────────────────────────────────────
      } else {
        shiny::tags$details(
          style = "margin-top:8px;",
          shiny::tags$summary(
            style = paste0(
              "cursor:pointer;font-size:0.88em;font-weight:600;",
              "color:#6c757d;padding:4px 0;user-select:none;"
            ),
            shiny::icon("pen-to-square"),
            shiny::HTML(" &nbsp;Edit sample names")
          ),
          shiny::tags$div(
            style = "margin-top:8px;",
            if (n_expr > 0)
              shiny::tagList(
                shiny::tags$p(
                  style = "font-size:0.82em;font-weight:600;margin:0 0 3px;",
                  paste0("Expression columns (", n_expr, ")")
                ),
                shiny::tags$p(
                  style = "font-size:0.78em;color:#888;margin:0 0 5px;",
                  "Edit names in the expression matrix."
                ),
                make_textarea("rename_expr_names", expr_nms, n_expr),
                shiny::actionButton(
                  session$ns("apply_rename_expr"),
                  shiny::tagList(shiny::icon("check"), "Apply"),
                  class = "btn-xs btn-secondary"
                ),
                shiny::tags$hr(style = "margin:10px 0;")
              ),
            if (n_meta > 0)
              shiny::tagList(
                shiny::tags$p(
                  style = "font-size:0.82em;font-weight:600;margin:0 0 3px;",
                  paste0("Metadata sample IDs (", n_meta, " rows)")
                ),
                shiny::tags$p(
                  style = "font-size:0.78em;color:#888;margin:0 0 5px;",
                  "Edit the SampleID column values."
                ),
                make_textarea("override_meta_names", meta_ids, n_meta),
                shiny::actionButton(
                  session$ns("apply_override_meta"),
                  shiny::tagList(shiny::icon("check"), "Apply"),
                  class = "btn-xs btn-secondary"
                )
              ),
            if (n_expr > 0 || n_meta > 0)
              shiny::tagList(
                shiny::tags$hr(style = "margin:10px 0;"),
                shiny::tags$p(
                  style = "font-size:0.82em;font-weight:600;margin:0 0 3px;",
                  "Remove samples"
                ),
                remove_samples_ui()
              )
          )
        )
      }
    })


    # ── Observer: keyword/token match (Option 1) ─────────────────────────────
    shiny::observeEvent(input$try_match_col, {
      shiny::req(rv$expr, rv$meta, input$match_col)
      mc <- input$match_col

      glog(paste0("Trying column-based match using metadata column '", mc, "'..."))

      aligned <- tryCatch(
        align_expr_to_meta(rv$expr, rv$meta,
                           id_col    = colnames(rv$meta)[1],
                           match_col = mc),
        error = function(e) list(expr       = rv$expr,
                                 meta       = rv$meta,
                                 match_type = "failed",
                                 warnings   = e$message)
      )

      rv$expr <- aligned$expr
      rv$meta <- .ensure_sampleid_col(aligned$meta)

      for (w in aligned$warnings) glog(w)

      if (aligned$match_type %in% c("token_match", "partial_token_match")) {
        rv$messages <- c(
          rv$messages[!grepl("Could not align|unmatched|Sample ID alignment: unmatched",
                             rv$messages, ignore.case = TRUE)],
          paste0("✅ Keyword match via '", mc, "' (", aligned$match_type, ")"),
          aligned$warnings
        )
        glog(paste0("Column match succeeded via '", mc, "' (", aligned$match_type, ")."))
        shiny::showNotification(
          shiny::HTML(paste0("✅ Matched via <b>", mc, "</b> (", aligned$match_type, ").")),
          type = "default", duration = 8
        )
      } else {
        rv$messages <- c(rv$messages, aligned$warnings)
        glog(paste0("Column match via '", mc, "' failed (", aligned$match_type, ")."))
        shiny::showNotification(
          shiny::HTML(paste0("❌ Keyword match via <b>", mc, "</b> failed.")),
          type = "error", duration = NULL
        )
      }
      .check_alignment()
    })


    # ── Observer: rename expression columns (Option 2) ───────────────────────
    shiny::observeEvent(input$apply_rename_expr, {
      shiny::req(rv$expr, input$rename_expr_names)
      new_names <- trimws(strsplit(input$rename_expr_names, "\n")[[1]])
      new_names <- new_names[nchar(new_names) > 0]
      n_cols    <- ncol(rv$expr)

      if (length(new_names) != n_cols) {
        glog(paste0("Expression rename failed: expected ", n_cols,
                    " names, got ", length(new_names), "."))
        shiny::showNotification(
          paste0("❌ Expected ", n_cols, " names, got ", length(new_names), "."),
          type = "error", duration = 8
        )
        return()
      }

      glog(paste0("Expression columns renamed (", n_cols, " samples)."))
      colnames(rv$expr) <- new_names

      if (!is.null(rv$meta)) {
        aligned <- tryCatch(
          align_expr_to_meta(rv$expr, rv$meta),
          error = function(e) list(expr       = rv$expr,
                                   meta       = rv$meta,
                                   match_type = "failed",
                                   warnings   = e$message)
        )
        rv$expr <- aligned$expr
        rv$meta <- .ensure_sampleid_col(aligned$meta)
        for (w in aligned$warnings) glog(w)
        rv$messages <- c(
          rv$messages[!grepl("Could not align|unmatched|Sample ID alignment: unmatched",
                             rv$messages, ignore.case = TRUE)],
          paste0("✅ Expression columns renamed; alignment: ", aligned$match_type)
        )
        glog(paste0("Post-rename alignment: ", aligned$match_type, "."))
        .check_alignment()
        shiny::showNotification(
          if (!rv$match_col_needed) "✅ Columns renamed and aligned."
          else "⚠️ Columns renamed - alignment still unresolved.",
          type = if (!rv$match_col_needed) "default" else "warning",
          duration = 6
        )
      } else {
        rv$messages <- c(rv$messages,
          paste0("✅ Expression matrix columns renamed (", n_cols, " samples)."))
        shiny::showNotification("✅ Expression columns renamed.", type = "default", duration = 5)
      }
    })


    # ── Observer: set metadata sample IDs (Option 3) ─────────────────────────
    shiny::observeEvent(input$apply_override_meta, {
      shiny::req(rv$meta, input$override_meta_names)
      new_names <- trimws(strsplit(input$override_meta_names, "\n")[[1]])
      new_names <- new_names[nchar(new_names) > 0]
      n_rows    <- nrow(rv$meta)

      if (length(new_names) != n_rows) {
        glog(paste0("Metadata ID override failed: expected ", n_rows,
                    " names, got ", length(new_names), "."))
        shiny::showNotification(
          paste0("❌ Expected ", n_rows, " names, got ", length(new_names), "."),
          type = "error", duration = 8
        )
        return()
      }

      glog(paste0("Metadata sample IDs set manually (", n_rows, " rows)."))
      meta       <- rv$meta
      rownames(meta) <- new_names
      other_cols <- setdiff(colnames(meta), "SampleID")
      rv$meta    <- cbind(SampleID = new_names,
                          meta[, other_cols, drop = FALSE],
                          stringsAsFactors = FALSE)

      if (!is.null(rv$expr)) {
        aligned <- tryCatch(
          align_expr_to_meta(rv$expr, rv$meta),
          error = function(e) list(expr       = rv$expr,
                                   meta       = rv$meta,
                                   match_type = "failed",
                                   warnings   = e$message)
        )
        rv$expr <- aligned$expr
        rv$meta <- .ensure_sampleid_col(aligned$meta)
        for (w in aligned$warnings) glog(w)
        rv$messages <- c(
          rv$messages[!grepl("Could not align|unmatched|Sample ID alignment: unmatched",
                             rv$messages, ignore.case = TRUE)],
          paste0("✅ Metadata sample IDs set; alignment: ", aligned$match_type)
        )
        glog(paste0("Post-override alignment: ", aligned$match_type, "."))
        .check_alignment()
        shiny::showNotification(
          if (!rv$match_col_needed) "✅ Sample IDs set and aligned."
          else "⚠️ Sample IDs set - alignment still unresolved.",
          type = if (!rv$match_col_needed) "default" else "warning",
          duration = 6
        )
      } else {
        rv$messages <- c(rv$messages,
          paste0("✅ Metadata sample IDs set manually (", n_rows, " rows)."))
        shiny::showNotification("✅ Metadata sample IDs updated.", type = "default", duration = 5)
      }
    })


    # ── Observer: remove samples (Option 1 / Mode B) ─────────────────────────
    # keep_expr_samples / keep_meta_samples hold what the user wants to KEEP.
    shiny::observeEvent(input$apply_remove_samples, {

      n_before_expr <- if (!is.null(rv$expr)) ncol(rv$expr) else 0L
      n_before_meta <- if (!is.null(rv$meta)) nrow(rv$meta) else 0L

      # ── Expression columns ────────────────────────────────────────────────
      if (!is.null(rv$expr)) {
        keep_e <- input$keep_expr_samples
        if (is.null(keep_e) || length(keep_e) == 0) {
          shiny::showNotification(
            "❌ Cannot remove all expression columns.",
            type = "error", duration = 6
          )
          return()
        }
        valid_keep <- intersect(keep_e, colnames(rv$expr))
        rv$expr <- rv$expr[, valid_keep, drop = FALSE]
      }

      # ── Metadata rows ─────────────────────────────────────────────────────
      if (!is.null(rv$meta)) {
        keep_m <- input$keep_meta_samples
        if (is.null(keep_m) || length(keep_m) == 0) {
          shiny::showNotification(
            "❌ Cannot remove all metadata rows.",
            type = "error", duration = 6
          )
          return()
        }
        meta_ids_cur <- if ("SampleID" %in% colnames(rv$meta))
          as.character(rv$meta$SampleID)
        else
          rownames(rv$meta)
        keep_rows <- meta_ids_cur %in% keep_m
        if (any(keep_rows)) {
          rv$meta <- rv$meta[keep_rows, , drop = FALSE]
        }
      }

      n_removed_expr <- n_before_expr - if (!is.null(rv$expr)) ncol(rv$expr) else 0L
      n_removed_meta <- n_before_meta - if (!is.null(rv$meta)) nrow(rv$meta) else 0L

      parts <- character(0)
      if (n_removed_expr > 0)
        parts <- c(parts, paste0(n_removed_expr, " expression column(s)"))
      if (n_removed_meta > 0)
        parts <- c(parts, paste0(n_removed_meta, " metadata row(s)"))

      if (length(parts) == 0) {
        shiny::showNotification("No samples were removed.", type = "message", duration = 4)
        return()
      }

      msg <- paste0("✅ Removed: ", paste(parts, collapse = " and "), ".")
      rv$messages <- c(rv$messages, msg)
      if (n_removed_expr > 0)
        glog(paste0("Samples removed from expression: ", n_removed_expr,
                    " removed, ", ncol(rv$expr), " remaining."))
      if (n_removed_meta > 0)
        glog(paste0("Rows removed from metadata: ", n_removed_meta,
                    " removed, ", nrow(rv$meta), " remaining."))
      shiny::showNotification(msg, type = "default", duration = 6)
      .check_alignment()
    })


    # ── 9. Return public reactives ───────────────────────────────────────────
    list(
      expr_data  = shiny::reactive(rv$expr),
      meta_data  = shiny::reactive(rv$meta),
      status     = shiny::reactive(rv$status),
      messages   = shiny::reactive(rv$messages),
      geo_result = shiny::reactive(rv$geo_result)
    )
  })
}
