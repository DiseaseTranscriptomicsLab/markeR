test_that("identify_variable_type classifies variable types correctly", {
  df <- data.frame(
    num = 1:20,
    bin = rep(c("A", "B"), 10),
    multi = rep(letters[1:4], 5),
    stringsAsFactors = FALSE
  )
  res <- identify_variable_type(df, cols = c("num", "bin", "multi"))
  expect_equal(as.character(res["num"]), "Numeric")
  expect_equal(as.character(res["bin"]), "Categorical Bin")
  expect_equal(as.character(res["multi"]), "Categorical Multi")
})

test_that("compute_stat_tests returns correct structure and methods", {
  set.seed(1701)
  df <- data.frame(
    score = rnorm(30),
    age = sample(20:50, 30, replace = TRUE),
    bin = rep(c("A", "B"), 15),
    multi = rep(letters[1:3], 10),
    stringsAsFactors = FALSE
  )
  res <- compute_stat_tests(df, target_var = "score")
  expect_true(is.list(res))
  expect_true(all(c("age", "bin", "multi") %in% names(res)))
  expect_true(all(c("metric", "p_value") %in% colnames(res$age)))
})

test_that("Score_VariableAssociation returns all expected outputs", {
  set.seed(1702)
  expr <- as.data.frame(matrix(rnorm(200), nrow = 20, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:20)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sampleID = colnames(expr),
    Group = rep(c("X", "Y"), each = 5),
    Age = sample(20:60, 10),
    stringsAsFactors = FALSE
  )
  gene_set <- list(Sig1 = paste0("Gene", 1:5))
  res <- Score_VariableAssociation(
    data = expr,
    metadata = metadata,
    cols = c("Group", "Age"),
    method = "logmedian",
    gene_set = gene_set,
    mode = "simple",
    printplt = FALSE
  )
  expect_true(is.list(res))
  expect_true(all(c("Overall", "Contrasts", "plot", "plot_contrasts", "plot_overall", "plot_distributions") %in% names(res)))
  expect_true(is.data.frame(res$Overall))
})

test_that("GSEA_VariableAssociation returns plot and data", {
  set.seed(1703)
  expr <- as.data.frame(matrix(rnorm(200), nrow = 20, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:20)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sampleID = colnames(expr),
    Group = rep(c("X", "Y"), each = 5),
    stringsAsFactors = FALSE
  )
  gene_set <- list(SigA = paste0("Gene", 1:5))
  res <- GSEA_VariableAssociation(
    data = expr,
    metadata = metadata,
    cols = "Group",
    gene_set = gene_set,
    printplt = FALSE
  )
  expect_true(is.list(res))
  expect_true(all(c("plot", "data") %in% names(res)))
  expect_true(is.data.frame(res$data))
  expect_true(inherits(res$plot, "ggplot"))
})

test_that("VariableAssociation runs for both score and GSEA methods", {
  set.seed(1704)
  expr <- as.data.frame(matrix(rnorm(120), nrow = 12, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:12)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sampleID = colnames(expr),
    Group = rep(c("X", "Y"), each = 5),
    Age = sample(20:60, 10),
    stringsAsFactors = FALSE
  )
  gene_set <- list(Sig1 = paste0("Gene", 1:4))
  res_score <- VariableAssociation(
    method = "ssGSEA",
    data = expr,
    metadata = metadata,
    cols = c("Group", "Age"),
    gene_set = gene_set,
    printplt = FALSE
  )
  expect_true(is.list(res_score))
  expect_true("Overall" %in% names(res_score))

  res_gsea <- VariableAssociation(
    method = "GSEA",
    data = expr,
    metadata = metadata,
    cols = "Group",
    gene_set = gene_set,
    printplt = FALSE
  )
  expect_true(is.list(res_gsea))
  expect_true(all(c("plot", "data") %in% names(res_gsea)))
})
