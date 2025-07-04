library(mockery)
library(devtools)

test_that("PlotScores returns a ggplot object for categorical variable", {
  set.seed(42)
  expr <- as.data.frame(matrix(rexp(60, rate = 0.2), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5),
    Age = seq(30, 75, length.out = 10)
  )
  gene_sets <- list(
    Signature1 = c("Gene1", "Gene2", "Gene3"),
    Signature2 = c("Gene4", "Gene5", "Gene6")
  )
  plt <- PlotScores(
    data = expr,
    metadata = metadata,
    gene_sets = gene_sets,
    method = "logmedian",
    Variable = "Group"
  )
  expect_true(inherits(plt, "gg"))
})

test_that("PlotScores returns a ggplot object for numeric variable", {
  set.seed(43)
  expr <- as.data.frame(matrix(rexp(60, rate = 0.2), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5),
    Age = seq(30, 75, length.out = 10)
  )
  gene_sets <- list(
    Signature1 = c("Gene1", "Gene2", "Gene3"),
    Signature2 = c("Gene4", "Gene5", "Gene6")
  )
  plt <- PlotScores(
    data = expr,
    metadata = metadata,
    gene_sets = gene_sets,
    method = "logmedian",
    Variable = "Age"
  )
  expect_true(inherits(plt, "gg"))
})

test_that("PlotScores returns a ggplot object for density plot when Variable is NULL", {
  set.seed(44)
  expr <- as.data.frame(matrix(rexp(60, rate = 0.2), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5),
    Age = seq(30, 75, length.out = 10)
  )
  gene_sets <- list(
    Signature1 = c("Gene1", "Gene2", "Gene3"),
    Signature2 = c("Gene4", "Gene5", "Gene6")
  )
  plt <- PlotScores(
    data = expr,
    metadata = metadata,
    gene_sets = gene_sets,
    method = "logmedian"
  )
  expect_true(inherits(plt, "gg"))
})

test_that("PlotScores returns a list with heatmap and volcano for method='all'", {
  set.seed(45)
  expr <- as.data.frame(matrix(rexp(60, rate = 0.2), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5),
    Age = seq(30, 75, length.out = 10)
  )
  gene_sets <- list(
    Signature1 = c("Gene1", "Gene2", "Gene3"),
    Signature2 = c("Gene4", "Gene5", "Gene6")
  )
  res <- PlotScores(
    data = expr,
    metadata = metadata,
    gene_sets = gene_sets,
    method = "all",
    Variable = "Group"
  )
  expect_true(is.list(res))
  expect_true(all(c("heatmap", "volcano") %in% names(res)))
  expect_true(inherits(res$heatmap, "gg"))
  expect_true(inherits(res$volcano, "gg"))
})



test_that("PlotScores calls PlotScores_Categorical for categorical variable", {
  # Setup mock
  mock_cat <- mock(list(dummy="categorical"))
  stub(PlotScores, "PlotScores_Categorical", mock_cat)
  set.seed(1)
  expr <- as.data.frame(matrix(rexp(20), nrow=4))
  rownames(expr) <- paste0("Gene", 1:4)
  colnames(expr) <- paste0("Sample", 1:5)
  metadata <- data.frame(sample = colnames(expr), Group = factor(rep(c("A", "B"), length.out=5)))
  gene_sets <- list(Signature1 = c("Gene1", "Gene2"))

  res <- PlotScores(expr, metadata, gene_sets, method="logmedian", Variable="Group")
  expect_equal(res$dummy, "categorical")
  expect_called(mock_cat, 1)
})

test_that("PlotScores calls PlotScores_Numeric for numeric variable", {
  mock_num <- mock(list(dummy="numeric"))
  stub(PlotScores, "PlotScores_Numeric", mock_num)
  set.seed(2)
  expr <- as.data.frame(matrix(rexp(20), nrow=4))
  rownames(expr) <- paste0("Gene", 1:4)
  colnames(expr) <- paste0("Sample", 1:5)
  metadata <- data.frame(sample = colnames(expr), Age = seq(20, 60, length.out=5))
  gene_sets <- list(Signature1 = c("Gene1", "Gene2"))

  res <- PlotScores(expr, metadata, gene_sets, method="logmedian", Variable="Age")
  expect_equal(res$dummy, "numeric")
  expect_called(mock_num, 1)
})

test_that("PlotScores calls PlotScores_Categorical for NULL variable", {
  mock_cat <- mock(list(dummy="density"))
  stub(PlotScores, "PlotScores_Categorical", mock_cat)
  set.seed(3)
  expr <- as.data.frame(matrix(rexp(20), nrow=4))
  rownames(expr) <- paste0("Gene", 1:4)
  colnames(expr) <- paste0("Sample", 1:5)
  metadata <- data.frame(sample = colnames(expr), Group = factor(rep(c("A", "B"), length.out=5)))
  gene_sets <- list(Signature1 = c("Gene1", "Gene2"))

  res <- PlotScores(expr, metadata, gene_sets, method="logmedian")
  expect_equal(res$dummy, "density")
  expect_called(mock_cat, 1)
})


# Helper for null coalesce (avoids error if no $plots)
`%||%` <- function(a, b) if (!is.null(a)) a else b


test_that("cohen_d returns the correct value for two groups", {
  set.seed(2025)
  x <- c(1.1, 2.3, 2.5, 3.1, 1.8)
  y <- c(5.2, 5.7, 6.1, 5.9, 6.3)
  manual_d <- (mean(x) - mean(y)) / sqrt(((var(x) * (length(x)-1)) + (var(y) * (length(y)-1))) / (length(x) + length(y) - 2))
  expect_equal(cohen_d(x, y), manual_d, tolerance = 1e-10)
})

test_that("compute_cohen_d returns expected values and contrasts", {
  # Simulate some scores with a grouping variable
  set.seed(2025)
  df <- data.frame(
    score = c(1,2,3,4,5,10,11,12,13,14),
    group = rep(c("A", "B"), each = 5)
  )
  results <- compute_cohen_d(df, variable = "group", quantitative_var = "score", mode = "simple")
  # Manual Cohen's d
  d_manual <- cohen_d(df$score[df$group == "A"], df$score[df$group == "B"])
  expect_equal(results$CohenD[1], d_manual, tolerance = 1e-10)
  expect_equal(results$Group1[1], "A")
  expect_equal(results$Group2[1], "B")
})

test_that("CohenD_allConditions returns correct structure and value", {
  set.seed(2025)
  expr <- as.data.frame(matrix(runif(40, 1, 10), nrow = 4, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:4)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5)
  )
  gene_sets <- list(Signature1 = c("Gene1", "Gene2"))
  res <- CohenD_allConditions(expr, metadata, gene_sets, variable = "Group", mode = "simple")
  expect_true(is.list(res))
  expect_true("Signature1" %in% names(res))
  # Check for correct structure
  expect_true(all(c("CohenD","PValue","padj") %in% names(res$Signature1)))
  # Check value is as expected for logmedian method
  dval <- res$Signature1$CohenD["logmedian", 1]
  # Calculate manually
  scores_list <- CalculateScores(expr, metadata, gene_sets, method = "logmedian")
  df <- scores_list$Signature1
  d_manual <- cohen_d(df$score[df$Group == "A"], df$score[df$Group == "B"])
  expect_equal(as.numeric(dval), d_manual, tolerance = 1e-10)
})

test_that("CohenF_allConditions returns correct structure and value", {
  set.seed(2025)
  expr <- as.data.frame(matrix(runif(40, 1, 10), nrow = 4, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:4)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Age = seq(20, 65, length.out = 10)
  )
  gene_sets <- list(Signature1 = c("Gene1", "Gene2"))
  res <- CohenF_allConditions(expr, metadata, gene_sets, variable = "Age")
  expect_true(is.list(res))
  expect_true("Signature1" %in% names(res))
  expect_true(all(c("CohenF","PValue","padj") %in% names(res$Signature1)))
  # Check value for logmedian method
  fval <- res$Signature1$CohenF["logmedian", 1]
  # Manual calculation
  scores_list <- CalculateScores(expr, metadata, gene_sets, method = "logmedian")
  df <- scores_list$Signature1
  model <- lm(score ~ Age, data = df)
  # Compute Cohen's f squared: f2 = R2/(1-R2); f = sqrt(f2)
  r2 <- summary(model)$r.squared
  manual_f <- sqrt(r2 / (1 - r2))
  expect_equal(as.numeric(fval), manual_f, tolerance = 1e-10)
})

test_that("PlotScores calculates p-values and cohen's d for categorical variables", {
  set.seed(42)
  expr <- as.data.frame(matrix(rexp(60, rate = 0.2), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5),
    Age = seq(30, 75, length.out = 10)
  )
  gene_sets <- list(
    Signature1 = c("Gene1", "Gene2", "Gene3"),
    Signature2 = c("Gene4", "Gene5", "Gene6")
  )
  plt <- PlotScores(
    data = expr,
    metadata = metadata,
    gene_sets = gene_sets,
    method = "logmedian",
    Variable = "Group",
    pvalcalc = TRUE,
    compute_cohen =TRUE
  )
  expect_true(inherits(plt, "gg"))
})


test_that("PlotScores calculates p-values and cohen's f for numeric variables", {
  set.seed(43)
  expr <- as.data.frame(matrix(rexp(60, rate = 0.2), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5),
    Age = seq(30, 75, length.out = 10)
  )
  gene_sets <- list(
    Signature1 = c("Gene1", "Gene2", "Gene3"),
    Signature2 = c("Gene4", "Gene5", "Gene6")
  )
  plt <- PlotScores(
    data = expr,
    metadata = metadata,
    gene_sets = gene_sets,
    method = "logmedian",
    Variable = "Age",
    pvalcalc = TRUE,
    compute_cohen =TRUE
  )
  expect_true(inherits(plt, "gg"))
})
