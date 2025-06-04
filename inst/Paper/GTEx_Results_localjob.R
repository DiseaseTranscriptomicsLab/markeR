
devtools::load_all()
setwd("data")
nperm=1000

signatures_bidirectional <- readRDS("SenescenceSigntures_Bidirectional.rds") # Divided by direction
GTEx_alltissues <- readRDS("GTExV8_voyagercorrected.rds") # https://github.com/DiseaseTranscriptomicsLab/voyAGEr/tree/main/Corrected_Counts
metadata_GTEx_alltissues <- readRDS("GTExV8_metadata.rds") # Restricted access



methods <- c("logmedian","ranking","ssGSEA")
gene_set <- list(HernandezSegura=signatures_bidirectional$HernandezSegura,
                 SAUL_SEN_MAYO=signatures_bidirectional$SAUL_SEN_MAYO)
tissues <- unique(metadata_GTEx_alltissues$SMTSD)

results_df_score <- data.frame(NULL)
results_df_gsea <- data.frame(NULL)

# Initialize progress bar
pb <- txtProgressBar(min = 0, max = length(tissues), style = 3)

for (i in seq_along(tissues)) {

  tissue <- tissues[i]

  subset_metadata <- metadata_GTEx_alltissues[metadata_GTEx_alltissues$SMTSD == tissue,]
  subset_data <- GTEx_alltissues[,subset_metadata$SAMPID]


  # Update progress bar
  setTxtProgressBar(pb, i)

  for (sig in names(gene_set)){

    signature <- list(gene_set[[sig]])
    names(signature) <- sig

    data_varassoc_gsea <- suppressWarnings(suppressMessages(GSEA_VariableAssociation(data=subset_data,
                                                                                     metadata=subset_metadata,
                                                                                     cols=c("AGE"),
                                                                                     mode="simple",
                                                                                     gene_set=signature)$data))

    data_varassoc_gsea <- data_varassoc_gsea[,c("NES","pval","Contrast")]
    data_varassoc_gsea$signature <- sig
    data_varassoc_gsea$method <- "GSEA"
    data_varassoc_gsea$tissue <- tissue

    results_df_gsea <- rbind(results_df_gsea, data_varassoc_gsea)

    for (method in methods){

      data_varassoc_score <- suppressWarnings(suppressMessages(Score_VariableAssociation(data=subset_data,
                                                                                         metadata=subset_metadata,
                                                                                         cols=c("AGE"), # SMRIN was a variable for batch correction
                                                                                         method=method,
                                                                                         gene_set = signature,
                                                                                         mode="simple", printplt = F)$Overall))


      df_permutations <- data.frame(NULL)
      for (j in 1:nperm){
        set.seed(j)
        metadata_subset_shuffleAGE <- subset_metadata
        metadata_subset_shuffleAGE$AGE <- sample(metadata_subset_shuffleAGE$AGE)

        cohend_shuffle <- suppressWarnings(suppressMessages(Score_VariableAssociation(data=subset_data[,metadata_subset_shuffleAGE$SAMPID],
                                                                                           metadata=metadata_subset_shuffleAGE,
                                                                                           cols=c("AGE"), # SMRIN was a variable for batch correction
                                                                                           method=method,
                                                                                           gene_set = signature,
                                                                                           mode="simple", printplt = F)$Overall))


        cohend_shuffle$signature <- sig
        cohend_shuffle$method <- method
        cohend_shuffle$tissue <- tissue
        df_permutations <- rbind(df_permutations, cohend_shuffle)

      }

      # calculateFPR

      df_permutations$Cohen_f
      fpr <- sum(df_permutations$Cohen_f > data_varassoc_score$Cohen_f)/length(df_permutations$Cohen_f)


      data_varassoc_score$signature <- sig
      data_varassoc_score$method <- method
      data_varassoc_score$tissue <- tissue
      data_varassoc_score$fpr <- fpr
      results_df_score <- rbind(results_df_score, data_varassoc_score)




    }

  }

}

# Close the progress bar
close(pb)


saveRDS( results_df_gsea,"results_df_gsea_GTEx_FPR.rds")
saveRDS( results_df_score,"results_df_score_GTEx_FPR.rds")
