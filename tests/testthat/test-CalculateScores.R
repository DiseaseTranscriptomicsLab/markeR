library(devtools)

test_that("CalculateScores logmedian matches manual calculation for unidirectional gene set", {
  set.seed(123)
  # Simulate positive expression data (genes as rows, samples as columns)
  expr <- as.data.frame(matrix(rexp(60, rate = 0.2), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)

  # Define a unidirectional gene set
  gene_sets <- list(Signature1 = c("Gene1", "Gene2", "Gene3"))

  # Manual calculation for Sample1
  genes <- gene_sets$Signature1
  data_subset <- log2(as.matrix(expr[genes, , drop=FALSE]) + 1)
  data_subset <- data_subset - apply(data_subset, 1, median)
  expected <- sum(data_subset[, "Sample1"]) / length(genes)

  # Calculate via function
  scores_list <- CalculateScores(expr, NULL, gene_sets, method="logmedian")
  score <- scores_list$Signature1$score[scores_list$Signature1$sample == "Sample1"]

  expect_equal(as.numeric(score), as.numeric(expected), tolerance=1e-10)
})

test_that("CalculateScores ranking method matches manual calculation for unidirectional gene set", {
  set.seed(123)
  expr <- as.data.frame(matrix(runif(30, 1, 100), nrow = 6, ncol = 5))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:5)
  gene_sets <- list(Signature1 = c("Gene1", "Gene3", "Gene5"))

  # Manual calculation for Sample2
  sample <- "Sample2"
  expr_sample <- expr[, sample]
  names(expr_sample) <- rownames(expr)
  expr_sample <- expr_sample[order(expr_sample, decreasing = FALSE)]
  manual_ranks <- match(gene_sets$Signature1, names(expr_sample))
  manual_score <- sum(manual_ranks, na.rm = TRUE) / nrow(expr)

  scores_list <- CalculateScores(expr, NULL, gene_sets, method="ranking")
  score <- scores_list$Signature1$score[scores_list$Signature1$sample == sample]

  expect_equal(as.numeric(score), as.numeric(manual_score), tolerance=1e-10)
})

test_that("CalculateScores returns all methods when method='all'", {
  set.seed(1)
  expr <- as.data.frame(matrix(rexp(18, rate = 0.2), nrow = 6, ncol = 3))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:3)
  gene_sets <- list(Signature1 = c("Gene1", "Gene2", "Gene3"))
  result <- CalculateScores(expr, NULL, gene_sets, method = "all")
  expect_true(all(c("ssGSEA","logmedian","ranking") %in% names(result)))
})

test_that("CalculateScores handles bidirectional gene set for logmedian", {
  set.seed(1)
  expr <- as.data.frame(matrix(abs(rnorm(60, 5, 2)), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  sig <- data.frame(Gene = c("Gene1", "Gene2", "Gene3", "Gene4"), Direction = c(1,1,-1,-1))
  gene_sets <- list(BiSig = sig)
  scores <- CalculateScores(expr, NULL, gene_sets, method = "logmedian")
  expect_true("score" %in% colnames(scores$BiSig))
  expect_equal(nrow(scores$BiSig), 10)
})


test_that("CalculateScores handles bidirectional gene set for ssGSEA", {
  set.seed(1)
  expr <- as.data.frame(matrix(abs(rnorm(60, 5, 2)), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  sig <- data.frame(Gene = c("Gene1", "Gene2", "Gene3", "Gene4"), Direction = c(1,1,-1,-1))
  gene_sets <- list(BiSig = sig)
  scores <- CalculateScores(expr, NULL, gene_sets, method = "ssGSEA")
  expect_true("score" %in% colnames(scores$BiSig))
  expect_equal(nrow(scores$BiSig), 10)
})

test_that("CalculateScores handles bidirectional gene set for ranking", {
  set.seed(1)
  expr <- as.data.frame(matrix(abs(rnorm(60, 5, 2)), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  sig <- data.frame(Gene = c("Gene1", "Gene2", "Gene3", "Gene4"), Direction = c(1,1,-1,-1))
  gene_sets <- list(BiSig = sig)
  scores <- CalculateScores(expr, NULL, gene_sets, method = "ranking")
  expect_true("score" %in% colnames(scores$BiSig))
  expect_equal(nrow(scores$BiSig), 10)
})
