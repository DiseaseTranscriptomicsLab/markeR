# Last run: 20 Nov 2025
set.seed("20112025")

library("ggplot2")
library("colorspace")
library("scales")
library("scater") 
library("reshape2")
#library("markeR")
library("ggbreak")
library("ggnewscale")
library("dplyr")
library("tidyr")  
library(patchwork)
library("ggpubr") 
library(purrr)
library(edgeR)
library(RColorBrewer)
library(pROC) 


# Your Cohen's d function
compute_cohens_d <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  if(n1 < 2 || n2 < 2) return(NA)
  m1 <- mean(x)
  m2 <- mean(y)
  s1 <- sd(x)
  s2 <- sd(y)
  pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (pooled_sd == 0) return(NA)
  d <- abs((m1 - m2) / pooled_sd)
  return(d)
}

# Path to your folder with functions
fun_dir <- "functions_Sept18_2025"

# Get all .R files in the folder
files <- list.files(paste0("../",fun_dir), pattern = "\\.R$", full.names = TRUE)

# Source them one by one
invisible(lapply(files, source))




metadata <- readRDS("../data_aux/metadata.rds")
corrcounts <- readRDS("../data_aux/corrcounts.rds")
signatures_bidirectional <- readRDS("../data_aux/SenescenceSigntures_Bidirectional.rds") # Divided by direction


# Expression data: genes (rows) × samples (columns)
df_expr <- corrcounts

# Class vector: named vector or factor
group_vec <-  metadata$Condition
names(group_vec) <- metadata$sampleID

# Signature gene sets
gene_sets <- list(
  CellAge = signatures_bidirectional$CellAge,
  HernandezSegura = signatures_bidirectional$HernandezSegura
)

# Percentages to evaluate
sample_pcts <- seq(1, 100, by = 1)

# Start at 5 genes, increment by 5 up to just below the smallest signature; 
#include the exact max of the small signature; add a few small steps above it; 
#then medium jumps to roughly half the largest signature; 
#finally, larger jumps up to the largest signature, including its maximum.


# Number of genes detected in the expression matrix
n_detected <- sapply(gene_sets, function(sig) {
  if (is.data.frame(sig)) {
    sum(sig[,1] %in% rownames(df_expr))
  } else {
    sum(sig %in% rownames(df_expr))
  }
})
 
len_HS <- n_detected["HernandezSegura"]
len_CA <- n_detected["CellAge"]

# 1. Small steps of 5 up to just below HS
vec1 <- seq(1, len_HS - 1, by=1)       # 5,10,…,50

# 2. Include the exact HS max
vec2 <- len_HS                           # 54

# 3. Slightly above HS
vec3 <- seq(len_HS + 1, len_HS + 15, by=1)  # 55, 59, etc.

# 4. Medium jumps to ~half of CellAge
vec4 <- seq(max(vec3)+10, round(len_CA/2), by=10)

# 5. Bigger jumps to max CellAge
vec5 <- seq(max(vec4)+50, len_CA, by=30)
vec5 <- vec5[vec5 <= len_CA]

# Combine
gene_numbers <- unique(c(vec1, vec2, vec3, vec4, vec5, len_CA))
gene_numbers

# testing
#gene_numbers <- c(5,54,500, 1206)
#sample_pcts <-  c(10,50,100)

# Result container
results <- data.frame()

# Loop through gene sets
for (sig_name in names(gene_sets)) {
  
  siggenes_original <- gene_sets[[sig_name]]
  
  # Subset to genes that exist in your expression matrix
  if (is.data.frame(siggenes_original)){
    full_genes <- siggenes_original[siggenes_original[,1] %in% row.names(df_expr),]
    n_full <- nrow(full_genes)
  } else {
    full_genes <- intersect(siggenes_original, rownames(df_expr))
    n_full <- length(full_genes)
  }
  
  # Define the number of genes to sample for this signature
  if (sig_name == "HernandezSegura") {
    n_genes_vec <- gene_numbers[gene_numbers <= n_full]  # don't exceed signature length
  } else {
    n_genes_vec <- gene_numbers[gene_numbers <= n_full]
  }
  
  # Loop over absolute gene numbers instead of percentages
  for (n_genes_to_sample in n_genes_vec) {
    
    set.seed(12345)
    if (is.data.frame(full_genes)){
      sel_genes <- full_genes[sample(1:nrow(full_genes), n_genes_to_sample),]
    } else {
      sel_genes <- sample(full_genes, n_genes_to_sample)
    }
    
    sigs_list <- list(sel_genes)
    names(sigs_list) <- sig_name
    
    # Loop over sample percentages as before
    for (s_pct in sample_pcts) {
      
      print(paste0("Percentage Samples: ", s_pct,"% | Genes sampled: ", n_genes_to_sample))
      
      # Subset samples
      group1 <- names(group_vec[group_vec == "Senescent"])
      group2 <- names(group_vec[group_vec == "Proliferative"])
      
      set.seed(12345)
      sel_group1 <- sample(group1, max(2, round(length(group1) * s_pct / 100)))
      set.seed(12345)
      sel_group2 <- sample(group2, max(2, round(length(group2) * s_pct / 100)))
      
      sel_samples <- c(sel_group1, sel_group2)
      
      df_subset <- df_expr[, sel_samples]
      metadata_subset <- metadata[metadata$sampleID %in% sel_samples, c("sampleID","Condition")]
      
      # SCORE approach
      df_Scores <- CalculateScores(data = df_subset[, metadata_subset$sampleID],
                                   metadata = metadata_subset,
                                   method = "logmedian",
                                   gene_sets = sigs_list)
      df_Scores <- as.data.frame(df_Scores[[sig_name]])
      cohen_d_score <- compute_cohens_d(df_Scores[df_Scores$Condition == "Senescent","score"], 
                                        df_Scores[df_Scores$Condition == "Proliferative","score"])
      
      p_score <- tryCatch({
        t.test(df_Scores$score ~ df_Scores$Condition)$p.value
      }, error = function(e) NA)
      
      # ENRICHMENT approach
      DEGs <- calculateDE(data = df_subset[, metadata_subset$sampleID],
                          metadata = metadata_subset,
                          variables = "Condition",
                          contrasts = c("Senescent - Proliferative"))
      
      res <- runGSEA(DEGList = DEGs, gene_sets = sigs_list)
      res <- res$`Senescent-Proliferative`
      
      NES <- res$NES
      p_enrich <- res$pval
      
      # Save
      results <- rbind(results, data.frame(
        signature = sig_name,
        gene_n = n_genes_to_sample,  # note change: store number of genes instead of percentage
        sample_pct = s_pct,
        cohend = cohen_d_score,
        p_score = p_score,
        NES = NES,
        p_enrich = p_enrich
      ))
      
      
    }
  }
}


saveRDS(results, "../data_aux/results_tradeoffs_absolutenbgenes_biggerres.rds")







# CODE FOR RUNNING WITH % GENES INSTEAD OF ABSOLUTE NUMBER
# metadata <- readRDS("metadata.rds")
# corrcounts <- readRDS("corrcounts.rds")
# signatures_bidirectional <- readRDS("SenescenceSigntures_Bidirectional.rds") # Divided by direction
# 
# 
# # Expression data: genes (rows) × samples (columns)
# df_expr <- corrcounts
# 
# # Class vector: named vector or factor
# group_vec <-  metadata$Condition
# names(group_vec) <- metadata$sampleID
# 
# # Signature gene sets
# gene_sets <- list(
#   CellAge = signatures_bidirectional$CellAge,
#   HernandezSegura = signatures_bidirectional$HernandezSegura
# )
# 
# # Percentages to evaluate
# gene_pcts <- seq(10, 100, by = 10)
# sample_pcts <- seq(10, 100, by = 10) 
# # Result container
# results <- data.frame()
# 
# # Loop through gene sets
# for (sig_name in names(gene_sets)) {
#   
#   siggenes_original <- gene_sets[[sig_name]]
#   
#   if (is.data.frame(siggenes_original)){
#     full_genes <- siggenes_original[siggenes_original[,1] %in% row.names(df_expr),]
#   } else {
#     full_genes <- intersect(siggenes_original, rownames(df_expr))
#   }
#   
#   for (g_pct in gene_pcts) {
#     
#     
#     if (is.data.frame(siggenes_original)){
#       ngenes_to_sample <-  max(1, round(nrow(full_genes) * g_pct / 100))
#       set.seed("12345")
#       sel_genes <- full_genes[sample(1:nrow(full_genes), ngenes_to_sample),]
#     } else {
#       ngenes_to_sample <-  max(1, round(length(full_genes) * g_pct / 100))
#       set.seed("12345")
#       sel_genes <- full_genes[sample(1:length(full_genes), ngenes_to_sample)]
#     }
#     sigs_list <- list(sel_genes)
#     names(sigs_list) <- sig_name
#     
#     for (s_pct in sample_pcts) {
#       
#       print(paste0("Percentage Samples: ", s_pct,"% | Percentage Genes: ", g_pct, "%"))
#       
#       # Subset samples
#       group1 <- names(group_vec[group_vec == "Senescent"])
#       group2 <- names(group_vec[group_vec == "Proliferative"])
#       
#       set.seed("12345")
#       sel_group1 <- sample(group1, max(2, round(length(group1) * s_pct / 100)))
#       set.seed("12345")
#       sel_group2 <- sample(group2, max(2, round(length(group2) * s_pct / 100)))
#       
#       sel_samples <- c(sel_group1, sel_group2)
#       #sel_expr <- expr_mat[sel_genes, sel_samples]
#       
#       df_subset <- df_expr[,sel_samples]
#       metadata_subset <- metadata[metadata$sampleID %in% sel_samples,c("sampleID","Condition")]
#       
#       # SCORE approach: mean expression score per sample
#       df_Scores <- CalculateScores(data = df_subset[,metadata_subset$sampleID],
#                                    metadata = metadata_subset,
#                                    method = "logmedian",
#                                    gene_sets =  sigs_list)
#       df_Scores <- as.data.frame(df_Scores[[sig_name]])
#       cohen_d_score <- compute_cohens_d(df_Scores[df_Scores$Condition == "Senescent","score"], df_Scores[df_Scores$Condition == "Proliferative","score"])
#       
#       # Score p-value (t-test)
#       p_score <- tryCatch({
#         t.test(df_Scores$score ~ df_Scores$Condition)$p.value
#       }, error = function(e) NA)
#       
#       # ENRICHMENT approach: gene ranking for fgsea
#       DEGs <- calculateDE(data = df_subset[,metadata_subset$sampleID],
#                           metadata = metadata_subset,
#                           variables = "Condition",
#                           contrasts = c("Senescent - Proliferative"))
#       #DEGs <- DEGs$`Senescent-Proliferative`
#       
#       res <- runGSEA(DEGList = DEGs,
#                      gene_sets = sigs_list )
#       res <- res$`Senescent-Proliferative`
#       
#       NES <-  res$NES
#       p_enrich <-res$pval
#       
#       # Save
#       results <- rbind(results, data.frame(
#         signature = sig_name,
#         gene_pct = g_pct,
#         sample_pct = s_pct,
#         cohend = cohen_d_score,
#         p_score = p_score,
#         NES = NES,
#         p_enrich = p_enrich
#       ))
#     }
#   }
# }
#
# # results$p_score <- p.adjust(results$p_score, method = "BH")
# # results$p_enrich <- p.adjust(results$p_enrich, method = "BH")



# SAVE RESULTS!!!!!!






