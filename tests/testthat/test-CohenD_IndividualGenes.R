test_that("CohenD_IndividualGenes returns a barplot and data for single group", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(2201)
  expr <- matrix(abs(rnorm(40)), nrow = 4, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:4)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    Sample = colnames(expr),
    Condition = rep(c("Case", "Control"), each = 5),
    stringsAsFactors = FALSE
  )
  res <- CohenD_IndividualGenes(
    data = expr,
    metadata = metadata,
    genes = rownames(expr),
    condition_var = "Condition",
    class = "Case",
    group_var = NULL,
    title = "Cohen's d Barplot",
    params = list(colors = c("#CCCCCC"))
  )
  expect_true(is.list(res))
  expect_true("plot" %in% names(res))
  expect_true("data" %in% names(res))
  expect_true(is.data.frame(res$data))
})

test_that("CohenD_IndividualGenes returns heatmap for multiple groups", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(2202)
  expr <- matrix(abs(rnorm(40)), nrow = 4, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:4)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    Sample = colnames(expr),
    Condition = rep(c("Case", "Control"), each = 5),
    Group = rep(c("A", "B"), times = 5),
    stringsAsFactors = FALSE
  )
  res <- CohenD_IndividualGenes(
    data = expr,
    metadata = metadata,
    genes = rownames(expr),
    condition_var = "Condition",
    class = "Case",
    group_var = "Group",
    title = "Cohen's d Heatmap",
    params = list(limits = c(0, 2))
  )
  expect_true(is.list(res))
  expect_true("plot" %in% names(res))
  expect_true("data" %in% names(res))
})

test_that("CohenD_IndividualGenes errors if samples in data are not in metadata", {
  set.seed(2203)
  expr <- matrix(abs(rnorm(12)), nrow = 3, ncol = 4)
  rownames(expr) <- paste0("Gene", 1:3)
  colnames(expr) <- paste0("Sample", 1:4)
  metadata <- data.frame(
    Sample = c("SampleA", "SampleB", "SampleC", "SampleD"),
    Condition = rep(c("Case", "Control"), each = 2),
    stringsAsFactors = FALSE
  )
  expect_error(
    CohenD_IndividualGenes(
      data = expr,
      metadata = metadata[-1, , drop = FALSE], # Remove one sample to break alignment
      genes = rownames(expr),
      condition_var = "Condition",
      class = "Case"
    ),
    "Not all samples in the data are described in the metadata"
  )
})

test_that("CohenD_IndividualGenes errors if group_var not in metadata", {
  set.seed(2204)
  expr <- matrix(abs(rnorm(12)), nrow = 3, ncol = 4)
  rownames(expr) <- paste0("Gene", 1:3)
  colnames(expr) <- paste0("Sample", 1:4)
  metadata <- data.frame(
    Sample = colnames(expr),
    Condition = rep(c("Case", "Control"), each = 2),
    stringsAsFactors = FALSE
  )
  expect_error(
    CohenD_IndividualGenes(
      data = expr,
      metadata = metadata,
      genes = rownames(expr),
      condition_var = "Condition",
      class = "Case",
      group_var = "NOT_A_COLUMN"
    ),
    "not found in the metadata"
  )
})

test_that("CohenD_IndividualGenes warns for too many genes", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(2205)
  expr <- matrix(abs(rnorm(320)), nrow = 32, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:32)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    Sample = colnames(expr),
    Condition = rep(c("Case", "Control"), each = 5),
    stringsAsFactors = FALSE
  )
  expect_warning(
    CohenD_IndividualGenes(
      data = expr,
      metadata = metadata,
      genes = rownames(expr),
      condition_var = "Condition",
      class = "Case"
    ),
    "Too many genes selected"
  )
})
