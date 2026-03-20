# =============================================================================
# Null distribution of NES values via age label permutation
# Per tissue, per model (baseline and corrected)
# xCell2 is run once per tissue; only AGE is permuted
# =============================================================================
#!/usr/bin/env Rscript

t_start <- Sys.time()
message("Started: ", format(t_start, "%Y-%m-%d %H:%M:%S"))

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(parallel)
})

fun_dir <- "../functions_Sept18_2025/"
invisible(lapply(list.files(fun_dir, pattern = "\\.R$", full.names = TRUE), source))

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------
N_PERM   <- 1000
N_CORES  <- min(detectCores() - 1, 49)
#OUT_DIR  <- "../Figures/Comment3_1"
#if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

message("Using ", N_CORES, " cores")

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
message("Loading data...")
#GTEx_alltissues          <- readRDS("../data_aux/GTExV8_voyagercorrected.rds")
GTEx_alltissues <- readRDS("../data_aux/data_gtex.rds")
metadata_GTEx_alltissues <- readRDS("../data_aux/GTExV8_metadata.rds")
signatures_bidirectional <- readRDS("../data_aux/SenescenceSigntures_Bidirectional.rds")
gene_sets <- list(HernandezSegura = signatures_bidirectional$HernandezSegura)

TISSUES <- unique(metadata_GTEx_alltissues$SMTSD)
message("Tissues: ", length(TISSUES), " | Permutations: ", N_PERM,
        " | Cores: ", N_CORES)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
run_one <- function(expr_t, meta_t, design) {
  de <- tryCatch(
    calculateDE(data = expr_t, metadata = meta_t, modelmat = design),
    error = function(e) NULL)
  if (is.null(de)) return(NA_real_)
  gs <- tryCatch(
    runGSEA(DEGList = de, gene_sets = gene_sets),
    error = function(e) NULL)
  if (is.null(gs)) return(NA_real_)
  gs[["AGE"]]$NES[gs[["AGE"]]$pathway == "HernandezSegura"]
}

# -----------------------------------------------------------------------------
# Function for one tissue
# -----------------------------------------------------------------------------
process_tissue <- function(tissue) {
  
  t_tiss <- Sys.time()
  idx    <- which(TISSUES == tissue)
  pfx    <- sprintf("[%02d/%02d] %s", idx, length(TISSUES), tissue)
  
  message(pfx, " — starting")
  
  meta_t           <- metadata_GTEx_alltissues[
    metadata_GTEx_alltissues$SMTSD == tissue, ]
  rownames(meta_t) <- meta_t$SAMPID
  #expr_t           <- GTEx_alltissues[, meta_t$SAMPID]
  expr_t           <- 2^GTEx_alltissues[[tissue]]
  expr_t           <- expr_t[, meta_t$SAMPID]
  
  message(pfx, " — running ", N_PERM, " permutations")
  
  rows <- vector("list", N_PERM)
  for (i in seq_len(N_PERM)) {
    meta_perm     <- meta_t
    meta_perm$AGE <- sample(meta_t$AGE)
    nes           <- run_one(expr_t, meta_perm,
                             model.matrix(~ AGE, data = meta_perm))
    rows[[i]] <- data.frame(tissue  = tissue,
                            perm_id = i,
                            NES     = nes)
    
    if (i %% 50 == 0)
      message(pfx, " — ", i, "/", N_PERM, " permutations done")
  }
  
  t_elapsed <- round(difftime(Sys.time(), t_tiss, units = "mins"), 1)
  message(pfx, " — DONE (", t_elapsed, " min)")
  
  do.call(rbind, rows)
}

# -----------------------------------------------------------------------------
# Run in parallel
# -----------------------------------------------------------------------------
null_results <- mclapply(
  TISSUES,
  process_tissue,
  mc.cores       = N_CORES,
  mc.preschedule = FALSE
)

# -----------------------------------------------------------------------------
# Compile and save
# -----------------------------------------------------------------------------
null_df  <- do.call(rbind, Filter(is.data.frame, null_results))
#out_file <- file.path(OUT_DIR, "null_NES_permutations_simpler_1000perms_fullGTExdata.rds")
saveRDS(null_df, "../data_aux/null_NES_permutations_simpler_1000perms_fullGTExdata.rds")

message("Saved.")
message("Total time: ",
        round(difftime(Sys.time(), t_start, units = "mins"), 1), " min")