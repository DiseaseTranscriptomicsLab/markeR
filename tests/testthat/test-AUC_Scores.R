library(devtools)

test_that("AUC_Scores returns a ggpubr/ggplot object and AUCs are in [0.5, 1]", {
  set.seed(44)
  expr <- as.data.frame(matrix(rexp(60, rate = 0.2), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5)
  )
  gene_sets <- list(
    Signature1 = c("Gene1", "Gene2"),
    Signature2 = c("Gene4", "Gene5")
  )
  p <- AUC_Scores(expr, metadata, gene_sets, method = "all", variable = "Group", nrow = 1)
  # Should be a ggpubr or gg object
  expect_true(any(c("ggpubr", "gg") %in% class(p)))
  # Extract heatmap data and check AUC range
  auc_list <- ROCAUC_Scores_Calculate(expr, metadata, gene_sets, method = "all", variable = "Group")
  for (method_name in names(auc_list)) {
    for (sig in names(auc_list[[method_name]])) {
      for (con in names(auc_list[[method_name]][[sig]])) {
        auc <- as.numeric(auc_list[[method_name]][[sig]][[con]]$AUC)
        expect_true(auc >= 0.5 && auc <= 1)
      }
    }
  }
})
