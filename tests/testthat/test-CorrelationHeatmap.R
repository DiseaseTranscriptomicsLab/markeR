test_that("CorrelationHeatmap returns a single heatmap for selected genes", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(1901)
  expr <- as.data.frame(matrix(runif(60, min = 0, max = 10), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  res <- CorrelationHeatmap(
    data = expr,
    genes = rownames(expr),
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )
  expect_true(is.list(res))
  expect_true("data" %in% names(res))
  expect_true("plot" %in% names(res))
  expect_true(is.matrix(res$data))
})

test_that("CorrelationHeatmap returns heatmaps separated by group", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(1902)
  expr <- as.data.frame(matrix(runif(60, min = 0, max = 10), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    SampleID = colnames(expr),
    Condition = rep(c("A", "B"), each = 5)
  )
  res <- CorrelationHeatmap(
    data = expr,
    metadata = metadata,
    genes = rownames(expr),
    separate.by = "Condition",
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )
  expect_true(is.list(res))
  expect_true("data" %in% names(res))
  expect_true("plot" %in% names(res))
  expect_true(is.list(res$data))
  expect_true(all(c("A", "B") %in% names(res$data)))
})

test_that("CorrelationHeatmap allows custom colorlist and limits_colorscale", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(1903)
  expr <- as.data.frame(matrix(runif(60, min = 0, max = 10), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  colorlist <- list(low = "black", mid = "yellow", high = "green")
  res <- CorrelationHeatmap(
    data = expr,
    genes = rownames(expr),
    colorlist = colorlist,
    limits_colorscale = c(-1, 0, 1),
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )
  expect_true(is.list(res))
  expect_true(is.matrix(res$data))
})

test_that("CorrelationHeatmap returns aux when detailedresults=TRUE", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(1904)
  expr <- as.data.frame(matrix(runif(60, min = 0, max = 10), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  res <- CorrelationHeatmap(
    data = expr,
    genes = rownames(expr),
    detailedresults = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE
  )
  expect_true("aux" %in% names(res))
  expect_true(is.list(res$aux))
})

test_that("CorrelationHeatmap errors if separate.by is given but metadata is missing", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(1905)
  expr <- as.data.frame(matrix(runif(60, min = 0, max = 10), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  expect_error(
    CorrelationHeatmap(
      data = expr,
      genes = rownames(expr),
      separate.by = "Condition",
      cluster_rows = FALSE,
      cluster_columns = FALSE
    ),
    "separate.by is not NULL but metadata is missing"
  )
})
