test_that("calculateDE returns expected structure with variables and contrast", {
  set.seed(1001)
  expr <- matrix(rpois(100, lambda = 20), nrow = 10, ncol = 10)
  rownames(expr) <- paste0("gene", 1:10)
  colnames(expr) <- paste0("sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )

  res <- calculateDE(
    data = expr,
    metadata = metadata,
    variables = "Group",
    contrasts = "A-B"
  )
  expect_true(is.list(res))
  expect_true("A-B" %in% names(res))
  expect_true(all(c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val","B") %in% names(res[["A-B"]])))
  expect_equal(nrow(res[["A-B"]]), 10)
})

test_that("calculateDE works with manually supplied design matrix", {
  set.seed(1002)
  expr <- matrix(rpois(100, lambda = 20), nrow = 10, ncol = 10)
  rownames(expr) <- paste0("gene", 1:10)
  colnames(expr) <- paste0("sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )
  design <- model.matrix(~0 + Group, data = metadata)
  colnames(design) <- c("A", "B")
  res <- calculateDE(
    data = expr,
    metadata = metadata,
    modelmat = design,
    contrasts = "A-B"
  )
  expect_true(is.list(res))
  expect_true("A-B" %in% names(res))
})

test_that("calculateDE returns all coefficients when no contrast is given", {
  set.seed(1003)
  expr <- matrix(rpois(100, lambda = 20), nrow = 10, ncol = 10)
  rownames(expr) <- paste0("gene", 1:10)
  colnames(expr) <- paste0("sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )
  res <- calculateDE(
    data = expr,
    metadata = metadata,
    variables = "Group"
  )
  expect_true(is.list(res))
  expect_true("A" %in% names(res))
  expect_true("B" %in% names(res))
})

test_that("calculateDE errors with mismatched design matrix", {
  set.seed(1004)
  expr <- matrix(rpois(100, lambda = 20), nrow = 10, ncol = 10)
  rownames(expr) <- paste0("gene", 1:10)
  colnames(expr) <- paste0("sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )
  design_bad <- matrix(1, nrow = 9, ncol = 2)
  expect_error(
    calculateDE(
      data = expr,
      metadata = metadata,
      modelmat = design_bad,
      contrasts = "A-B"
    ),
    "Rows in 'modelmat' must match the number of samples"
  )
})

test_that("calculateDE handles NAs with ignore_NAs", {
  set.seed(1005)
  expr <- matrix(rpois(100, lambda = 20), nrow = 10, ncol = 10)
  rownames(expr) <- paste0("gene", 1:10)
  colnames(expr) <- paste0("sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )
  metadata$Group[1] <- NA
  res <- calculateDE(
    data = expr,
    metadata = metadata,
    variables = "Group",
    contrasts = "A-B",
    ignore_NAs = TRUE
  )
  expect_true(is.list(res))
  expect_true("A-B" %in% names(res))
  expect_true(nrow(res[["A-B"]]) == 10)
})

test_that("calculateDE output structure is correct", {
  set.seed(1006)
  expr <- matrix(rpois(100, lambda = 20), nrow = 10, ncol = 10)
  rownames(expr) <- paste0("gene", 1:10)
  colnames(expr) <- paste0("sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )
  res <- calculateDE(
    data = expr,
    metadata = metadata,
    variables = "Group",
    contrasts = "A-B"
  )
  tab <- res[[1]]
  expect_true(is.data.frame(tab))
  expect_equal(ncol(tab), 6) # logFC, AveExpr, t, P.Value, adj.P.Val, B
  expect_equal(nrow(tab), 10)
})
