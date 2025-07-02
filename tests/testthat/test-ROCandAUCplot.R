test_that("ROCandAUCplot returns ROC plot and AUC values for standard input", {
  skip_if_not_installed("pROC")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(2101)
  expr <- matrix(rexp(30, rate = 0.5), nrow = 3, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:3)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    SampleID = colnames(expr),
    Condition = rep(c("A", "B"), each = 5),
    Group = rep(c("G1", "G2"), times = 5),
    stringsAsFactors = FALSE
  )
  res <- ROCandAUCplot(
    data = expr,
    metadata = metadata,
    genes = rownames(expr),
    condition_var = "Condition",
    class = "A",
    plot_type = "all",
    roc_params = list(nrow=1),
    title = "Test ROC"
  )
  expect_true(is.list(res))
  expect_true("roc_plot" %in% names(res))
  expect_true("auc_plot" %in% names(res))
  expect_true(is.data.frame(res$auc_values))
})

test_that("ROCandAUCplot returns auc_plot for plot_type='auc'", {
  skip_if_not_installed("pROC")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(2102)
  expr <- matrix(rexp(24, rate = 0.5), nrow = 3, ncol = 8)
  rownames(expr) <- paste0("Gene", 1:3)
  colnames(expr) <- paste0("Sample", 1:8)
  metadata <- data.frame(
    SampleID = colnames(expr),
    Condition = rep(c("A", "B"), each = 4),
    stringsAsFactors = FALSE
  )
  res <- ROCandAUCplot(
    data = expr,
    metadata = metadata,
    genes = rownames(expr),
    condition_var = "Condition",
    class = "A",
    plot_type = "auc",
    auc_params = list(colors = "#3B415B")
  )
  expect_true(is.list(res))
  expect_true("auc_plot" %in% names(res))
  expect_true("auc_values" %in% names(res))
  expect_false("roc_plot" %in% names(res))
})


test_that("ROCandAUCplot errors if samples in data are not in metadata", {
  set.seed(2104)
  expr <- matrix(rexp(12, rate = 0.5), nrow = 3, ncol = 4)
  rownames(expr) <- paste0("Gene", 1:3)
  colnames(expr) <- paste0("Sample", 1:4)
  metadata <- data.frame(
    SampleID = c("SampleA", "SampleB", "SampleC", "SampleD"),
    Condition = rep(c("A", "B"), each = 2),
    stringsAsFactors = FALSE
  )
  expect_error(
    ROCandAUCplot(
      data = expr,
      metadata = metadata[-1, , drop = FALSE], # Remove one sample to break alignment
      genes = rownames(expr),
      condition_var = "Condition",
      class = "A"
    ),
    "Not all samples in the data are described in the metadata"
  )
})

test_that("ROCandAUCplot errors if group_var not in metadata", {
  set.seed(2105)
  expr <- matrix(rexp(12, rate = 0.5), nrow = 3, ncol = 4)
  rownames(expr) <- paste0("Gene", 1:3)
  colnames(expr) <- paste0("Sample", 1:4)
  metadata <- data.frame(
    SampleID = colnames(expr),
    Condition = rep(c("A", "B"), each = 2),
    stringsAsFactors = FALSE
  )
  expect_error(
    ROCandAUCplot(
      data = expr,
      metadata = metadata,
      genes = rownames(expr),
      condition_var = "Condition",
      class = "A",
      group_var = "NOT_A_COLUMN"
    ),
    "not found in the metadata"
  )
})

test_that("ROCandAUCplot warns for too many genes", {
  skip_if_not_installed("pROC")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  set.seed(2106)
  expr <- matrix(rexp(320, rate = 0.5), nrow = 32, ncol = 10)
  rownames(expr) <- paste0("Gene", 1:32)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    SampleID = colnames(expr),
    Condition = rep(c("X", "Y"), each = 5),
    stringsAsFactors = FALSE
  )
  expect_warning(
    ROCandAUCplot(
      data = expr,
      metadata = metadata,
      genes = rownames(expr),
      condition_var = "Condition",
      class = "X"
    ),
    "Too many genes selected"
  )
})
