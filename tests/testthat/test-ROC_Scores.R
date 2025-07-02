library(devtools)

test_that("ROCAUC_Scores_Calculate returns correct structure and AUC for two-group variable", {
  set.seed(42)
  expr <- as.data.frame(matrix(rexp(40, rate = 0.2), nrow = 4, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:4)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5)
  )
  gene_sets <- list(Signature1 = c("Gene1", "Gene2"))

  # Run the ROC calculation
  res <- ROCAUC_Scores_Calculate(expr, metadata, gene_sets, method = "logmedian", variable = "Group")
  expect_true(is.list(res))
  expect_true("logmedian" %in% names(res))
  expect_true("Signature1" %in% names(res$logmedian))
  # Only one contrast expected for 2-group (A - B)
  contr <- names(res$logmedian$Signature1)
  expect_true(length(contr) == 1)
  auc <- res$logmedian$Signature1[[contr]]$AUC
  # Compute manually with pROC
  scores <- CalculateScores(expr, metadata, gene_sets, method = "logmedian")$Signature1
  scores$Group <- as.factor(scores$Group)
  manual_roc <- pROC::roc(Group ~ score, data = scores, quiet = TRUE)
  manual_auc <- as.numeric(pROC::auc(manual_roc))
  # Allow for flip if AUC < 0.5
  if (manual_auc < 0.5) manual_auc <- 1 - manual_auc
  expect_equal(as.numeric(auc), manual_auc, tolerance = 1e-5)
})

test_that("ROC_Scores returns a ggplot object (or list) and has correct titles", {
  set.seed(43)
  expr <- as.data.frame(matrix(rexp(40, rate = 0.2), nrow = 4, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:4)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5)
  )
  gene_sets <- list(Signature1 = c("Gene1", "Gene2"))
  plt <- ROC_Scores(expr, metadata, gene_sets, method = "logmedian", variable = "Group", grid = FALSE)
  expect_true(is.list(plt))
  # Should be a list of ggplot objects, each with correct title and subtitle
  for (p in plt) {
    expect_true("gg" %in% class(p))
    expect_true(grepl("Signature1", p$labels$title))
    expect_true(any(grepl("A", p$labels$subtitle), grepl("B", p$labels$subtitle)))
  }
})
