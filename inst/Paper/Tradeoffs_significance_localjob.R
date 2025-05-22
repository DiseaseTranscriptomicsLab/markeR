library("ggplot2")
library("colorspace")
library("scales")
library("scater")
library("reshape2")
library("markeR")
library("ggbreak")
library("ggnewscale")

corrcounts <- readRDS("data/corrcounts.rds")
metadata <- readRDS("data/metadata.rds")
signatures_bidirectional <- readRDS("data/SenescenceSigntures_Bidirectional.rds") # Divided by direction


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
gene_pcts <-  seq(2,100,by=1)
sample_pcts <- seq(2,100,by=1)


# Result container
results <- data.frame()

# Loop through gene sets
for (sig_name in names(gene_sets)) {

  siggenes_original <- gene_sets[[sig_name]]

  if (is.data.frame(siggenes_original)){
    full_genes <- siggenes_original[siggenes_original[,1] %in% row.names(df_expr),]
  } else {
    full_genes <- intersect(siggenes_original, rownames(df_expr))
  }

  for (g_pct in gene_pcts) {


    if (is.data.frame(siggenes_original)){
      ngenes_to_sample <-  max(1, round(nrow(full_genes) * g_pct / 100))
      set.seed("12345")
      sel_genes <- full_genes[sample(1:nrow(full_genes), ngenes_to_sample),]
    } else {
      ngenes_to_sample <-  max(1, round(length(full_genes) * g_pct / 100))
      set.seed("12345")
      sel_genes <- full_genes[sample(1:length(full_genes), ngenes_to_sample)]
    }
    sigs_list <- list(sel_genes)
    names(sigs_list) <- sig_name

    for (s_pct in sample_pcts) {

      print(paste0("Percentage Samples: ", s_pct,"% | Percentage Genes: ", g_pct, "%"))

      # Subset samples
      group1 <- names(group_vec[group_vec == "Senescent"])
      group2 <- names(group_vec[group_vec == "Proliferative"])

      set.seed("12345")
      sel_group1 <- sample(group1, max(2, round(length(group1) * s_pct / 100)))
      set.seed("12345")
      sel_group2 <- sample(group2, max(2, round(length(group2) * s_pct / 100)))

      sel_samples <- c(sel_group1, sel_group2)
      #sel_expr <- expr_mat[sel_genes, sel_samples]

      df_subset <- df_expr[,sel_samples]
      metadata_subset <- metadata[metadata$sampleID %in% sel_samples,c("sampleID","Condition")]

      # SCORE approach: mean expression score per sample
      df_Scores <- CalculateScores(data = df_subset[,metadata_subset$sampleID],
                                   metadata = metadata_subset,
                                   method = "logmedian",
                                   gene_sets =  sigs_list)
      df_Scores <- as.data.frame(df_Scores[[sig_name]])
      cohen_d_score <- compute_cohens_d(df_Scores[df_Scores$Condition == "Senescent","score"], df_Scores[df_Scores$Condition == "Proliferative","score"])

      # Score p-value (t-test)
      p_score <- tryCatch({
        t.test(df_Scores$score ~ df_Scores$Condition)$p.value
      }, error = function(e) NA)

      # ENRICHMENT approach: gene ranking for fgsea
      DEGs <- calculateDE(data = df_subset[,metadata_subset$sampleID],
                          metadata = metadata_subset,
                          variables = "Condition",
                          contrasts = c("Senescent - Proliferative"))
      #DEGs <- DEGs$`Senescent-Proliferative`

      res <- runGSEA(DEGList = DEGs,
                     gene_sets = sigs_list )
      res <- res$`Senescent-Proliferative`

      NES <-  res$NES
      p_enrich <-res$pval

      # Save
      results <- rbind(results, data.frame(
        signature = sig_name,
        gene_pct = g_pct,
        sample_pct = s_pct,
        cohend = cohen_d_score,
        p_score = p_score,
        NES = NES,
        p_enrich = p_enrich
      ))
    }
  }
}

# results$p_score <- p.adjust(results$p_score, method = "BH")
# results$p_enrich <- p.adjust(results$p_enrich, method = "BH")

#saveRDS(results,"data/Tradeoffs_v1.rds") # correct for BH for all p-values
#saveRDS(results,"data/Tradeoffs_v2.rds") # no correction, correction done after based on the type of plot
saveRDS(results,"data/Tradeoffs_v3.rds") # Fix order of samples id data subset
