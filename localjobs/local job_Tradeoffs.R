# Last run: 20 Nov 2025
set.seed("20112025")

metadata <- readRDS("metadata.rds")
corrcounts <- readRDS("corrcounts.rds")
signatures_bidirectional <- readRDS("SenescenceSigntures_Bidirectional.rds") # Divided by direction


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
gene_pcts <- seq(10, 100, by = 10)
sample_pcts <- seq(10, 100, by = 10)


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
