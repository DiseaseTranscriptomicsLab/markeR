# =============================================================================
# Null distribution of NES values via age label permutation
# Per tissue, per model (baseline and corrected)
# xCell2 is run once per tissue; only AGE is permuted
# =============================================================================

t_start <- Sys.time()
message("Started: ", format(t_start, "%Y-%m-%d %H:%M:%S"))

suppressPackageStartupMessages({
  library(xCell2)
  library(limma)
  library(dplyr)
})

fun_dir <- "../functions_Sept18_2025/"
invisible(lapply(list.files(fun_dir, pattern = "\\.R$", full.names = TRUE), source))

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
N_PERM       <- 2   # increase to 500-1000 if time allows
N_COVARIATES <- 6
OUT_DIR      <- "../Figures/Comment3_1"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
message("Loading data...")
GTEx_alltissues          <- readRDS("../data_aux/GTExV8_voyagercorrected.rds")
metadata_GTEx_alltissues <- readRDS("../data_aux/GTExV8_metadata.rds")
signatures_bidirectional <- readRDS("../data_aux/SenescenceSigntures_Bidirectional.rds")

gene_sets <- list(HernandezSegura = signatures_bidirectional$HernandezSegura)

data("BlueprintEncode.xCell2Ref")
TISSUES <- unique(metadata_GTEx_alltissues$SMTSD)
message("  Tissues: ", length(TISSUES), " | Permutations per model: ", N_PERM)

# -----------------------------------------------------------------------------
# Helpers (minimal)
# -----------------------------------------------------------------------------
compute_min_shared_genes <- function(expr_mat, xcell_ref,
                                     margin = 0.1, floor = 0.3) {
  overlap <- length(intersect(rownames(expr_mat), xcell_ref@genes_used)) /
    length(xcell_ref@genes_used)
  max(floor, overlap - margin)
}

select_xcell_covariates <- function(xcell_mat, n = N_COVARIATES) {
  names(sort(apply(xcell_mat, 1, var), decreasing = TRUE))[seq_len(n)]
}

build_design_baseline <- function(meta) model.matrix(~ AGE, data = meta)

build_design_corrected <- function(meta, xcell_mat, covariates) {
  samples  <- rownames(meta)
  ct_df    <- as.data.frame(t(xcell_mat[covariates, samples, drop = FALSE]))
  colnames(ct_df) <- make.names(colnames(ct_df))
  df <- data.frame(AGE = meta[samples, "AGE"], ct_df, row.names = samples)
  model.matrix(~ AGE + ., data = df)
}

run_one <- function(expr_t, meta_t, design) {
  # returns NES for HernandezSegura, or NA on failure
  de <- tryCatch(
    calculateDE(data = expr_t, metadata = meta_t, modelmat = design),
    error = function(e) NULL
  )
  if (is.null(de)) return(NA_real_)
  gs <- tryCatch(
    runGSEA(DEGList = de, gene_sets = gene_sets),
    error = function(e) NULL
  )
  if (is.null(gs)) return(NA_real_)
  gs[["AGE"]]$NES[gs[["AGE"]]$pathway == "HernandezSegura"]
}

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------
null_results <- list()

for (tissue in TISSUES) {
  
  message("\n[", which(TISSUES == tissue), "/", length(TISSUES), "] ", tissue)
  
  meta_t           <- metadata_GTEx_alltissues[
    metadata_GTEx_alltissues$SMTSD == tissue, ]
  rownames(meta_t) <- meta_t$SAMPID
  expr_t           <- GTEx_alltissues[, meta_t$SAMPID]
  
  # xCell2 runs ONCE per tissue — reused across all permutations
  min_shared   <- compute_min_shared_genes(expr_t, BlueprintEncode.xCell2Ref)
  xcell_scores <- tryCatch(
    xCell2Analysis(mix          = expr_t,
                   xcell2object = BlueprintEncode.xCell2Ref,
                   minSharedGenes = min_shared),
    error = function(e) { message("  xCell2 failed: ", e$message); NULL }
  )
  if (is.null(xcell_scores)) next
  
  covariates <- select_xcell_covariates(xcell_scores)
  message("  Covariates: ", paste(covariates, collapse = ", "))
  message("  Running ", N_PERM, " permutations...")
  
  perm_rows <- vector("list", N_PERM * 2)
  idx       <- 1L
  
  for (i in seq_len(N_PERM)) {
    
    # Permute age labels — breaks age-expression association, preserves structure
    meta_perm     <- meta_t
    meta_perm$AGE <- sample(meta_t$AGE)
    
    # Baseline model
    design_base <- build_design_baseline(meta_perm)
    nes_base    <- run_one(expr_t, meta_perm, design_base)
    
    # Corrected model (same permuted AGE, same covariates)
    design_corr <- build_design_corrected(meta_perm, xcell_scores, covariates)
    nes_corr    <- run_one(expr_t, meta_perm, design_corr)
    
    perm_rows[[idx]]     <- data.frame(tissue  = tissue,
                                       model   = "Non-corrected",
                                       perm_id = i,
                                       NES     = nes_base)
    perm_rows[[idx + 1]] <- data.frame(tissue  = tissue,
                                       model   = "Corrected",
                                       perm_id = i,
                                       NES     = nes_corr)
    idx <- idx + 2L
    
    if (i %% 50 == 0) message("    ", i, "/", N_PERM, " done")
  }
  
  null_results[[tissue]] <- do.call(rbind, perm_rows)
}

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
null_df  <- do.call(rbind, null_results)
out_file <- file.path(OUT_DIR, "null_NES_permutations.rds")
saveRDS(null_df, out_file)
message("\nSaved to: ", out_file)
message("Total time: ",
        round(difftime(Sys.time(), t_start, units = "mins"), 1), " min")