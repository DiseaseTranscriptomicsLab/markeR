test_that("ExpressionHeatmap returns scaled data and a ComplexHeatmap object", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(2001)
  expr <- matrix(rnorm(25), nrow = 5, ncol = 5)
  rownames(expr) <- paste0("Gene", 1:5)
  colnames(expr) <- paste0("Sample", 1:5)
  metadata <- data.frame(
    Sample = colnames(expr),
    Condition = rep(c("A", "B"), length.out = 5),
    Batch = rep(c("X", "Y"), length.out = 5),
    stringsAsFactors = FALSE
  )
  res <- ExpressionHeatmap(
    data = expr,
    metadata = metadata,
    genes = rownames(expr),
    annotate.by = c("Condition", "Batch"),
    cluster_columns = FALSE,
    title = "Demo Expression Heatmap",
    scale_position = "right",
    legend_position = "top",
    titlesize = 14
  )
  expect_true(is.list(res))
  expect_true(all(c("data", "plot") %in% names(res)))
  expect_true(is.matrix(res$data))
})

test_that("ExpressionHeatmap uses annotation_colors if provided", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(2002)
  expr <- matrix(rnorm(25), nrow = 5, ncol = 5)
  rownames(expr) <- paste0("Gene", 1:5)
  colnames(expr) <- paste0("Sample", 1:5)
  metadata <- data.frame(
    Sample = colnames(expr),
    Condition = rep(c("A", "B"), length.out = 5),
    Batch = rep(c("X", "Y"), length.out = 5),
    stringsAsFactors = FALSE
  )
  annotation_colors <- list(
    Condition = c(A = "orange", B = "purple"),
    Batch = c(X = "green", Y = "blue")
  )
  res <- ExpressionHeatmap(
    data = expr,
    metadata = metadata,
    genes = rownames(expr),
    annotate.by = c("Condition", "Batch"),
    annotation_colors = annotation_colors,
    cluster_columns = TRUE,
    scale_position = "right",
    legend_position = "top"
  )
  expect_true(is.list(res))
  expect_true(all(c("data", "plot") %in% names(res)))
})

test_that("ExpressionHeatmap errors if annotate.by is specified but metadata is NULL", {
  skip_if_not_installed("ComplexHeatmap")
  set.seed(2003)
  expr <- matrix(rnorm(16), nrow = 4, ncol = 4)
  rownames(expr) <- paste0("Gene", 1:4)
  colnames(expr) <- paste0("Sample", 1:4)
  expect_error(
    ExpressionHeatmap(
      data = expr,
      metadata = NULL,
      genes = rownames(expr),
      annotate.by = "Condition"
    ),
    "annotate.by is specified but metadata is NULL"
  )
})

test_that("ExpressionHeatmap works with different legend and scale positions", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(2004)
  expr <- matrix(rnorm(25), nrow = 5, ncol = 5)
  rownames(expr) <- paste0("Gene", 1:5)
  colnames(expr) <- paste0("Sample", 1:5)
  metadata <- data.frame(
    Sample = colnames(expr),
    Condition = rep(c("A", "B"), length.out = 5),
    Batch = rep(c("X", "Y"), length.out = 5),
    stringsAsFactors = FALSE
  )
  # Try all valid combinations
  for (scale_pos in c("right", "top", "bottom")) {
    for (legend_pos in c("top", "right", "bottom")) {
      res <- ExpressionHeatmap(
        data = expr,
        metadata = metadata,
        genes = rownames(expr),
        annotate.by = c("Condition", "Batch"),
        scale_position = scale_pos,
        legend_position = legend_pos
      )
      expect_true(is.list(res))
      expect_true(all(c("data", "plot") %in% names(res)))
    }
  }
})

test_that("ExpressionHeatmap show_row_names and show_column_names arguments work", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(2005)
  expr <- matrix(rnorm(25), nrow = 5, ncol = 5)
  rownames(expr) <- paste0("Gene", 1:5)
  colnames(expr) <- paste0("Sample", 1:5)
  metadata <- data.frame(
    Sample = colnames(expr),
    Condition = rep(c("A", "B"), length.out = 5),
    Batch = rep(c("X", "Y"), length.out = 5),
    stringsAsFactors = FALSE
  )
  res <- ExpressionHeatmap(
    data = expr,
    metadata = metadata,
    genes = rownames(expr),
    annotate.by = c("Condition", "Batch"),
    show_row_names = FALSE,
    show_column_names = TRUE
  )
  expect_true(is.list(res))
  expect_true(all(c("data", "plot") %in% names(res)))
})
