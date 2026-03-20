library(parallel) 
library(ggpubr)
fun_dir <- "../functions_Sept18_2025/"
invisible(lapply(list.files(fun_dir, pattern = "\\.R$", full.names = TRUE), source))

GTEx_alltissues          <- readRDS("../data_aux/GTExV8_voyagercorrected.rds")
metadata_GTEx_alltissues <- readRDS("../data_aux/GTExV8_metadata.rds")
signatures_bidirectional <- readRDS("../data_aux/SenescenceSigntures_Bidirectional.rds")

nperm    <- 1000
methods  <- c("logmedian", "ranking", "ssGSEA")
gene_set <- list(
  HernandezSegura = signatures_bidirectional$HernandezSegura,
  SAUL_SEN_MAYO   = signatures_bidirectional$SAUL_SEN_MAYO
)
tissues  <- unique(metadata_GTEx_alltissues$SMTSD)
N_CORES  <- min(detectCores() - 1, length(tissues))

t_start  <- Sys.time()
message(sprintf("[%s] Starting — %d tissues, %d cores, %d permutations",
                format(t_start, "%H:%M:%S"), length(tissues), N_CORES, nperm))

# -----------------------------------------------------------------------------
# Function for one tissue
# -----------------------------------------------------------------------------
process_tissue <- function(tissue) {
  
  t0  <- Sys.time()
  idx <- which(tissues == tissue)
  message(sprintf("[%s] START  [%02d/%02d] %s",
                  format(t0, "%H:%M:%S"), idx, length(tissues), tissue))
  
  subset_metadata      <- metadata_GTEx_alltissues[
    metadata_GTEx_alltissues$SMTSD == tissue, ]
  subset_metadata$AGE  <- subset_metadata$AGE - mean(subset_metadata$AGE)
  subset_data          <- GTEx_alltissues[, subset_metadata$SAMPID]
  
  tissue_score_rows <- data.frame(NULL)
  tissue_gsea_rows  <- data.frame(NULL)
  
  for (sig in names(gene_set)) {
    
    signature        <- list(gene_set[[sig]])
    names(signature) <- sig
    
    # ---- GSEA ----
    data_varassoc_gsea <- suppressWarnings(suppressMessages(
      GSEA_VariableAssociation(
        data     = subset_data,
        metadata = subset_metadata,
        cols     = "AGE",
        mode     = "simple",
        gene_set = signature
      )$data
    ))
    data_varassoc_gsea           <- data_varassoc_gsea[, c("NES", "pval", "Contrast")]
    data_varassoc_gsea$signature <- sig
    data_varassoc_gsea$method    <- "GSEA"
    data_varassoc_gsea$tissue    <- tissue
    tissue_gsea_rows             <- rbind(tissue_gsea_rows, data_varassoc_gsea)
    
    # ---- Score-based + permutation FPR ----
    for (method in methods) {
      
      data_varassoc_score <- suppressWarnings(suppressMessages(
        Score_VariableAssociation(
          data     = subset_data,
          metadata = subset_metadata,
          cols     = "AGE",
          method   = method,
          gene_set = signature,
          mode     = "simple",
          printplt = FALSE
        )$Overall
      ))
      
      df_permutations <- data.frame(NULL)
      for (j in seq_len(nperm)) {
        set.seed(j)
        metadata_shuffled     <- subset_metadata
        metadata_shuffled$AGE <- sample(metadata_shuffled$AGE)
        
        cohend_shuffle <- suppressWarnings(suppressMessages(
          Score_VariableAssociation(
            data     = subset_data[, metadata_shuffled$SAMPID],
            metadata = metadata_shuffled,
            cols     = "AGE",
            method   = method,
            gene_set = signature,
            mode     = "simple",
            printplt = FALSE
          )$Overall
        ))
        cohend_shuffle$signature <- sig
        cohend_shuffle$method    <- method
        cohend_shuffle$tissue    <- tissue
        df_permutations          <- rbind(df_permutations, cohend_shuffle)
      }
      
      fpr <- sum(df_permutations$Cohen_f > data_varassoc_score$Cohen_f) /
        nrow(df_permutations)
      
      data_varassoc_score$signature <- sig
      data_varassoc_score$method    <- method
      data_varassoc_score$tissue    <- tissue
      data_varassoc_score$fpr       <- fpr
      tissue_score_rows             <- rbind(tissue_score_rows, data_varassoc_score)
    }
  }
  
  elapsed <- round(difftime(Sys.time(), t0, units = "mins"), 1)
  message(sprintf("[%s] DONE   [%02d/%02d] %s (%s min)",
                  format(Sys.time(), "%H:%M:%S"), idx, length(tissues), tissue, elapsed))
  
  list(score = tissue_score_rows, gsea = tissue_gsea_rows)
}

# -----------------------------------------------------------------------------
# Run in parallel
# -----------------------------------------------------------------------------
results_list <- mclapply(
  tissues,
  process_tissue,
  mc.cores       = N_CORES,
  mc.preschedule = FALSE
)

# -----------------------------------------------------------------------------
# Combine and save
# -----------------------------------------------------------------------------
results_df_score <- do.call(rbind, lapply(results_list, `[[`, "score"))
results_df_gsea  <- do.call(rbind, lapply(results_list, `[[`, "gsea"))

saveRDS(results_df_gsea,  "../data_aux/results_df_gsea_GTEx_FPR_agecentered.rds")
saveRDS(results_df_score, "../data_aux/results_df_score_GTEx_FPR_agecentered.rds")

message(sprintf("[%s] All done — total time: %s min",
                format(Sys.time(), "%H:%M:%S"),
                round(difftime(Sys.time(), t_start, units = "mins"), 1)))