#!/usr/bin/env Rscript
# =============================================================================
# Cell composition correction for age-senescence associations in GTEx
# Addresses reviewer comment 3.1
#
# Workflow per tissue:
#   1. Run xCell2 to estimate cell type enrichment scores per sample
#   2. Select the most variable cell types as covariates (data-driven)
#   3. Build baseline (~ AGE) and corrected (~ AGE + cell types) design matrices
#   4. Run calculateDE() with each design -> limma moderated t-statistics for AGE
#   5. Run runGSEA() on both ranked gene lists
#   6. Collect nominal p-values; apply BH correction across all tissues at end
#   7. Visualise NES with and without correction across all 49 tissues
#
# Usage:
#   Rscript composition_correction_job.R
#   nohup Rscript composition_correction_job.R > logs/composition_correction.log 2>&1 &
#
# Output files (written to OUT_DIR):
#   composition_correction_results.rds   — full results data frame
#   gsea_composition_correction.pdf      — visualisation
# =============================================================================

t_start <- Sys.time()
message("=============================================================")
message("  Cell composition correction analysis — all GTEx tissues")
message("  Started: ", format(t_start, "%Y-%m-%d %H:%M:%S"))
message("=============================================================\n")

# -----------------------------------------------------------------------------
# HOW xCell2 WORKS
# -----------------------------------------------------------------------------
# xCell2 estimates cell type enrichment scores from bulk RNA-seq data using a
# reference panel of known pure cell type expression profiles. It operates in
# two stages:
#
# (1) TRAINING: For each cell type, xCell2 identifies genes that are
#     specifically and consistently highly expressed in that cell type relative
#     to all others in the reference panel (not just highly expressed, but
#     *relatively* highly expressed). Multiple signatures per cell type are
#     built using different parameter combinations (stringency thresholds,
#     reference sample subsets), producing an ensemble that is more robust than
#     a single signature. This is why internal signature names contain suffixes
#     like "#_0.333_1_37_0.97" encoding those parameters.
#
# (2) SCORING: For each sample, xCell2 ranks all genes by expression and asks
#     whether the signature genes for a given cell type are collectively ranked
#     higher than expected by chance — conceptually similar to ssGSEA. A high
#     score means that cell type's transcriptional signal is enriched in the
#     sample; a score near zero means no detectable enrichment. A spillover
#     correction is then applied to reduce cross-talk between related cell
#     types that share marker genes.
#
# Unlike CIBERSORTx, which deconvolves bulk samples into cell type proportions
# that sum to 1 by solving a mixture model, xCell2 scores each cell type
# independently. Scores are not proportions, not bounded, and not comparable
# across cell types in absolute terms. Within a cell type however, higher
# scores consistently reflect greater presence of that cell type's
# transcriptional signal across samples — which is sufficient for our purpose
# of using them as covariates to control for age-related compositional
# variation in GTEx tissues.
#
# ON MULTIPLE SIGNATURES PER CELL TYPE:
# Internally, xCell2 stores multiple gene set signatures per cell type,
# identified by names such as:
#   "CD4-positive, alpha-beta T cell#_0.333_1_37_0.97"
# The suffix encodes training parameters (quantile threshold, number of genes,
# performance score). This ensemble approach makes the enrichment scores more
# robust. xCell2Analysis() aggregates across all signatures per cell type
# and returns a single score per cell type per sample — the xcell_scores matrix
# contains one row per cell type (43 total). The internal naming is irrelevant
# for downstream use.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(xCell2)
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

# Load markeR functions directly from source
fun_dir <- "../functions_Sept18_2025/"
files   <- list.files(fun_dir, pattern = "\\.R$", full.names = TRUE)
invisible(lapply(files, source))

# -----------------------------------------------------------------------------
# 1. Setup
# -----------------------------------------------------------------------------

# Number of cell type covariates to include in the corrected model.
# Using 6 covariates provides meaningful compositional correction while
# keeping the model parsimonious relative to typical GTEx tissue sample sizes
# (~100-500 samples). Too many covariates risk overfitting and collinearity.
N_COVARIATES <- 6

# Tissues with significant age-senescence signal in the paper (Figure 7),
# used to annotate the visualisation.
TISSUE_SIGNAL <- c(
  "Artery - Aorta",
  "Breast - Mammary Tissue",
  "Cells - Cultured fibroblasts",
  "Thyroid"
)

# Output directory
OUT_DIR <- "../Figures/Comment3_1"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# -----------------------------------------------------------------------------
# 2. Load data
# -----------------------------------------------------------------------------

message("Loading data...")

# GTEx v8 expression matrix: normalised NON-log-transformed counts (genes x samples)
# Pre-processed and batch-corrected as per Schneider et al. 2024 (voyAGEr)
GTEx_alltissues <- readRDS("../data_aux/GTExV8_voyagercorrected.rds")

# Sample metadata: must contain columns SMTSD (tissue) and AGE (numeric, donor age)
metadata_GTEx_alltissues <- readRDS("../data_aux/GTExV8_metadata.rds")

# Senescence gene sets with directionality (data frames with cols: gene, enrichment)
signatures_bidirectional <- readRDS("../data_aux/SenescenceSigntures_Bidirectional.rds")

# We analyse HernandezSegura only — the most robust senescence signature from
# the paper's benchmarking, derived from transcriptomic data across multiple
# senescence models, with explicit up/downregulation directionality
gene_sets <- list(
  HernandezSegura = signatures_bidirectional$HernandezSegura
)

# Load xCell2 BlueprintEncode pre-trained reference.
# This reference covers 43 cell types derived from Blueprint and ENCODE
# bulk RNA-seq datasets — a broad, well-characterised panel suitable for
# diverse human tissue types.
message("  Loading xCell2 BlueprintEncode reference...")
data("BlueprintEncode.xCell2Ref")

# Get all tissues present in the metadata
TISSUES <- unique(metadata_GTEx_alltissues$SMTSD)
message("  Tissues to analyse: ", length(TISSUES))
message("  Data loaded.\n")

# -----------------------------------------------------------------------------
# 3. Helper functions
# -----------------------------------------------------------------------------

# Select the N most variable cell types from xCell2 output as covariates.
#
# Rationale: rather than specifying tissue-specific cell type lists (which
# requires prior biological knowledge and introduces subjectivity), we select
# the cell types that vary most across samples within each tissue. These are
# the cell types most likely to correlate with age and therefore confound the
# age-senescence association. Purely invariant cell types carry no information
# as covariates and would only consume degrees of freedom.
select_xcell_covariates <- function(xcell_mat, n_covars = N_COVARIATES) {
  vars     <- apply(xcell_mat, 1, var)
  # sort descending by variance and take top n
  selected <- names(sort(vars, decreasing = TRUE))[seq_len(n_covars)]
  selected
}

# Build baseline design matrix: intercept + AGE.
# This replicates the model used in the paper for Figure 7.
build_design_baseline <- function(metadata) {
  model.matrix(~ AGE, data = metadata)
}

# Build corrected design matrix: intercept + AGE + scaled xCell2 scores.
#
# Cell type scores are scaled (z-scored) before inclusion to:
#   (a) improve numerical stability in limma's matrix algebra
#   (b) make coefficients comparable across cell types with different score ranges
#
# The AGE coefficient from this model represents the age effect on gene
# expression after accounting for variation in cell type composition —
# i.e., the age signal that cannot be explained by compositional shifts alone.
build_design_corrected <- function(metadata, xcell_mat, covariates) {
  samples         <- rownames(metadata)
  ct_df           <- as.data.frame(t(xcell_mat[covariates, samples, drop = FALSE]))
  colnames(ct_df) <- make.names(colnames(ct_df))
  ct_df_scaled    <- as.data.frame(scale(ct_df))
  df <- data.frame(AGE = metadata[samples, "AGE"],
                   ct_df_scaled,
                   row.names = samples)
  model.matrix(~ AGE + ., data = df)
}

# Determine the minimum shared genes threshold dynamically.
#
# xCell2 requires that the expression matrix shares a minimum fraction of
# genes with its reference. Our GTEx matrices are heavily filtered (only
# the most highly expressed genes are retained per tissue after the voyAGEr
# preprocessing pipeline), so the observed overlap (~0.45-0.65) is lower than
# xCell2's default threshold of 0.9. This is not an ID mismatch — gene symbols
# are consistent — but a consequence of expression-based filtering.
# We set the threshold to 10% below the observed overlap to ensure the analysis
# runs while remaining as stringent as possible, with a hard floor of 0.3 to
# avoid running with genuinely insufficient overlap.
compute_min_shared_genes <- function(expr_mat, xcell_ref, margin = 0.1, floor = 0.3) {
  overlap <- length(intersect(rownames(expr_mat), xcell_ref@genes_used)) /
    length(xcell_ref@genes_used)
  max(floor, overlap - margin)
}

# Extract GSEA results for a single signature and model from runGSEA output.
# Saves nominal p-value only — BH correction is applied globally at the end
# to ensure fair comparison across all tissues (see Step 6).
extract_gsea_result <- function(gsea_out, sig_name, model_label, tissue,
                                covariates) {
  df <- gsea_out[["AGE"]]
  df <- df[df$pathway == sig_name, ]
  data.frame(
    tissue      = tissue,
    signature   = sig_name,
    model       = model_label,
    NES         = df$NES,
    pval        = df$pval,
    # store which cell types were used as covariates for transparency
    covariates  = paste(covariates, collapse = "; "),
    n_covars    = length(covariates)
  )
}

# -----------------------------------------------------------------------------
# 4. Main loop — all tissues
# -----------------------------------------------------------------------------

results_all  <- list()
overlap_info <- list()

for (tissue in TISSUES) {
  
  t_tissue <- Sys.time()
  message("----------------------------------------------------------")
  message("  Tissue: ", tissue, " (", which(TISSUES == tissue), "/",
          length(TISSUES), ")")
  message("----------------------------------------------------------")
  
  # Subset expression matrix and metadata to this tissue
  meta_t           <- metadata_GTEx_alltissues[metadata_GTEx_alltissues$SMTSD == tissue, ]
  rownames(meta_t) <- meta_t$SAMPID
  expr_t           <- GTEx_alltissues[, meta_t$SAMPID]
  
  # --- Step 1: Compute gene overlap and set dynamic threshold ----------------
  # Record overlap for reporting in the methods section
  overlap_pct <- length(intersect(rownames(expr_t),
                                  BlueprintEncode.xCell2Ref@genes_used)) /
    length(BlueprintEncode.xCell2Ref@genes_used)
  overlap_info[[tissue]] <- round(overlap_pct * 100, 1)
  message("  Gene overlap with xCell2 reference: ", overlap_info[[tissue]], "%")
  
  min_shared <- compute_min_shared_genes(expr_t, BlueprintEncode.xCell2Ref)
  message("  minSharedGenes set to: ", round(min_shared, 3))
  
  # --- Step 2: Run xCell2 deconvolution --------------------------------------
  # Estimate cell type enrichment scores for all 43 cell types across samples
  message("  Running xCell2...")
  xcell_scores <- tryCatch(
    xCell2Analysis(
      mix            = expr_t,
      xcell2object   = BlueprintEncode.xCell2Ref,
      minSharedGenes = min_shared
    ),
    error = function(e) {
      message("  WARNING: xCell2 failed for ", tissue, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(xcell_scores)) {
    message("  Skipping tissue due to xCell2 failure.\n")
    next
  }
  
  # --- Step 3: Select covariates ---------------------------------------------
  # Select the N most variable cell types — these are most likely to correlate
  # with age and confound the senescence signal
  covariates <- select_xcell_covariates(xcell_scores, n_covars = N_COVARIATES)
  message("  Selected covariates (top ", N_COVARIATES, " by variance): ",
          paste(covariates, collapse = ", "))
  
  # --- Step 4: Build design matrices -----------------------------------------
  design_base <- build_design_baseline(meta_t)
  design_corr <- build_design_corrected(meta_t, xcell_scores, covariates)
  
  # --- Step 5: Differential expression ---------------------------------------
  # calculateDE fits a limma linear model and returns moderated t-statistics
  # for each gene. The AGE coefficient captures the per-gene age association.
  # Note: calculateDE log2-transforms internally, so raw counts are passed.
  message("  Running calculateDE (baseline)...")
  de_base <- tryCatch(
    calculateDE(data = expr_t, metadata = meta_t, modelmat = design_base),
    error = function(e) {
      message("  WARNING: calculateDE (baseline) failed: ", e$message); NULL
    }
  )
  
  message("  Running calculateDE (corrected)...")
  de_corr <- tryCatch(
    calculateDE(data = expr_t, metadata = meta_t, modelmat = design_corr),
    error = function(e) {
      message("  WARNING: calculateDE (corrected) failed: ", e$message); NULL
    }
  )
  
  if (is.null(de_base) || is.null(de_corr)) {
    message("  Skipping tissue due to calculateDE failure.\n")
    next
  }
  
  # --- Step 6: GSEA ----------------------------------------------------------
  # Genes are ranked by their moderated t-statistic for AGE (handled internally
  # by runGSEA for bidirectional gene sets). The resulting NES summarises the
  # coordinated age-associated expression change of the signature genes.
  message("  Running GSEA...")
  gsea_base <- tryCatch(
    runGSEA(DEGList = de_base, gene_sets = gene_sets),
    error = function(e) {
      message("  WARNING: GSEA (baseline) failed: ", e$message); NULL
    }
  )
  
  gsea_corr <- tryCatch(
    runGSEA(DEGList = de_corr, gene_sets = gene_sets),
    error = function(e) {
      message("  WARNING: GSEA (corrected) failed: ", e$message); NULL
    }
  )
  
  if (is.null(gsea_base) || is.null(gsea_corr)) {
    message("  Skipping tissue due to GSEA failure.\n")
    next
  }
  
  # --- Collect results -------------------------------------------------------
  for (sig_name in names(gene_sets)) {
    results_all[[paste(tissue, sig_name, "Baseline",  sep = "|")]] <-
      extract_gsea_result(gsea_base, sig_name, "Baseline",  tissue, covariates)
    results_all[[paste(tissue, sig_name, "Corrected", sep = "|")]] <-
      extract_gsea_result(gsea_corr, sig_name, "Corrected", tissue, covariates)
  }
  
  t_elapsed <- round(difftime(Sys.time(), t_tissue, units = "mins"), 1)
  message("  Done in ", t_elapsed, " min\n")
}

# -----------------------------------------------------------------------------
# 5. Compile results
# -----------------------------------------------------------------------------

results_df <- do.call(rbind, results_all)
rownames(results_df) <- NULL

message("\nGene overlap per tissue (% of xCell2 reference genes):")
print(as.data.frame(overlap_info))

# -----------------------------------------------------------------------------
# 6. Multiple testing correction
# -----------------------------------------------------------------------------
# BH correction is applied separately for baseline and corrected models,
# across all tissues and signatures tested. This mirrors the approach used
# in the original paper (Figure 7) and ensures fair comparison: each model
# is corrected for the same number of tests (one per tissue per signature).
# Correcting jointly across both models would be inappropriate as they are
# not independent tests.

results_df <- results_df %>%
  group_by(model, signature) %>%
  mutate(padj = p.adjust(pval, method = "BH"),
         sig  = padj < 0.05) %>%
  ungroup()

message("\nResults (", nrow(results_df), " rows):")
print(as.data.frame(results_df))

# Save full results
out_rds <- file.path(OUT_DIR, "composition_correction_results.rds")
saveRDS(results_df, out_rds)
message("\nResults saved to: ", out_rds)

# -----------------------------------------------------------------------------
# 7. Visualisation — all 49 tissues
# -----------------------------------------------------------------------------
# With 49 tissues, a standard dot/line plot becomes unreadable on a single
# panel. We use a heatmap-style layout showing NES for baseline and corrected
# side by side, with significance overlaid as asterisks. Tissues are grouped
# by whether they showed signal in the paper. This allows the reviewer to see
# at a glance: (a) which tissues retain signal after correction, (b) whether
# any new signal emerges, and (c) how much the NES changes overall.
# 
# model_colours <- c("Baseline"  = "#2171b5",
#                    "Corrected" = "#cb181d")
# 
# plot_df <- results_df %>%
#   mutate(
#     tissue_group = ifelse(tissue %in% TISSUE_SIGNAL,
#                           "Signal in paper", "No signal in paper"),
#     tissue_group = factor(tissue_group,
#                           levels = c("Signal in paper", "No signal in paper")),
#     model        = factor(model, levels = c("Baseline", "Corrected")),
#     sig_label    = ifelse(sig, "*", "")
#   )
# 
# # Ordered by baseline NES within each group for easier reading
# tissue_order <- plot_df %>%
#   filter(model == "Baseline") %>%
#   arrange(tissue_group, NES) %>%
#   pull(tissue) %>%
#   unique()
# 
# plot_df <- plot_df %>%
#   mutate(tissue = factor(tissue, levels = tissue_order))
# 
# # Paired dot plot with connecting line, faceted by tissue group
# # Each tissue appears as two dots (baseline blue, corrected red) connected by
# # a grey line. The direction of the line shows how correction changes the NES.
# # Filled circle = FDR < 0.05; open circle = not significant.
# plot_all <- ggplot(plot_df,
#                    aes(x = NES, y = tissue, colour = model)) +
#   geom_vline(xintercept = 0, linetype = "dashed",
#              colour = "grey50", linewidth = 0.4) +
#   # grey line connects baseline to corrected for each tissue
#   geom_line(aes(group = tissue), colour = "grey75", linewidth = 0.5) +
#   # filled = significant, open = not significant
#   geom_point(aes(shape = sig), size = 3, stroke = 1.2, fill = "white") +
#   geom_text(aes(label = sig_label, x = NES + sign(NES) * 0.12),
#             size = 3.5, show.legend = FALSE) +
#   scale_shape_manual(
#     values = c("FALSE" = 21, "TRUE"  = 19),
#     labels = c("FALSE" = "FDR \u2265 0.05", "TRUE" = "FDR < 0.05"),
#     name   = NULL
#   ) +
#   scale_colour_manual(values = model_colours, name = "Model") +
#   facet_grid(tissue_group ~ ., scales = "free_y", space = "free_y") +
#   labs(
#     x        = "Normalised Enrichment Score (NES)",
#     y        = NULL,
#     title    = "HernandezSegura age-senescence GSEA across 49 GTEx tissues",
#     subtitle = "Baseline (~ AGE) vs. cell composition-corrected (~ AGE + xCell2)\nFilled = FDR < 0.05 | * marks significance"
#   ) +
#   theme_cowplot(11) +
#   theme(
#     strip.background  = element_rect(fill = "grey92"),
#     strip.text.y      = element_text(angle = 0),
#     legend.position   = "bottom",
#     axis.text.y       = element_text(size = 8)
#   )
# 
# out_pdf <- file.path(OUT_DIR, "gsea_composition_correction_all_tissues.pdf")
# ggsave(out_pdf, plot = plot_all, width = 8, height = 18)
# message("Plot saved to: ", out_pdf)
# 
# # Also save a focused version showing only the 8 selected tissues for the
# # main figure in the revision (signal + expected-no-signal)
# TISSUES_FOCUSED <- c(
#   "Artery - Aorta", "Breast - Mammary Tissue",
#   "Cells - Cultured fibroblasts", "Thyroid",
#   "Artery - Tibial", "Skin - Sun Exposed (Lower leg)",
#   "Lung", "Pancreas"
# )
# 
# plot_focused <- plot_all %+%
#   filter(plot_df, tissue %in% TISSUES_FOCUSED)
# 
# out_pdf_focused <- file.path(OUT_DIR, "gsea_composition_correction_focused.pdf")
# ggsave(out_pdf_focused, plot = plot_focused, width = 8, height = 6)
# message("Focused plot saved to: ", out_pdf_focused)

# -----------------------------------------------------------------------------
# 8. Done
# -----------------------------------------------------------------------------

t_total <- round(difftime(Sys.time(), t_start, units = "mins"), 1)
message("\n=============================================================")
message("  Analysis complete!")
message("  Total time: ", t_total, " min")
message("  Tissues analysed: ", length(unique(results_df$tissue)))
message("  Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("=============================================================")