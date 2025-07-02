test_that("plotPCA returns ggplot or ggarrange object and data for basic use", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("edgeR")
  set.seed(2301)
  data <- abs(matrix(rnorm(1000), nrow=50, ncol=20))
  colnames(data) <- paste0("Sample", 1:20)
  rownames(data) <- paste0("Gene", 1:50)
  metadata <- data.frame(Sample = colnames(data),
                         Group = rep(c("A", "B"), each = 10),
                         stringsAsFactors = FALSE)
  res <- plotPCA(data, metadata, ColorVariable = "Group", pointSize = 3)
  expect_true(is.list(res))
  expect_true(all(c("plt", "data") %in% names(res)))
})

test_that("plotPCA supports custom gene selection and multiple PCs", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("edgeR")
  set.seed(2302)
  data <- abs(matrix(rnorm(400), nrow=20, ncol=20))
  colnames(data) <- paste0("Sample", 1:20)
  rownames(data) <- paste0("Gene", 1:20)
  metadata <- data.frame(Sample = colnames(data),
                         Group = rep(c("A", "B"), each = 10),
                         stringsAsFactors = FALSE)
  res <- plotPCA(data, metadata, genes = rownames(data)[1:10], PCs = list(c(1,2), c(2,3)), ColorVariable = "Group", pointSize = 6, legend_nrow = 1, legend_position = "bottom")
  expect_true(is.list(res))
  expect_true(all(c("plt", "data") %in% names(res)))
})

test_that("plotPCA errors with too few genes for chosen PCs", {
  skip_if_not_installed("edgeR")
  set.seed(2303)
  data <- abs(matrix(rnorm(6), nrow=2, ncol=3))
  colnames(data) <- paste0("Sample", 1:3)
  rownames(data) <- paste0("Gene", 1:2)
  metadata <- data.frame(Sample = colnames(data), Group = rep("A", 3), stringsAsFactors = FALSE)
  expect_error(
    plotPCA(data, metadata, genes = rownames(data)[1:2], PCs = list(c(1,2), c(2,3)), ColorVariable = "Group"),
    "Number of genes too low for number of chosen PCs"
  )
})

test_that("plotPCA errors with only one gene", {
  skip_if_not_installed("edgeR")
  set.seed(2304)
  data <- abs(matrix(rnorm(3), nrow=1, ncol=3))
  colnames(data) <- paste0("Sample", 1:3)
  rownames(data) <- "Gene1"
  metadata <- data.frame(Sample = colnames(data), Group = rep("A", 3), stringsAsFactors = FALSE)
  expect_error(
    plotPCA(data, metadata, genes = "Gene1", ColorVariable = "Group"),
    "Error: Number of genes should be >1"
  )
})

test_that("plotPCA supports custom colors and legend settings", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggpubr")
  skip_if_not_installed("edgeR")
  set.seed(2305)
  data <- abs(matrix(rnorm(100), nrow=10, ncol=10))
  colnames(data) <- paste0("Sample", 1:10)
  rownames(data) <- paste0("Gene", 1:10)
  metadata <- data.frame(Sample = colnames(data),
                         Group = rep(c("A", "B"), each = 5),
                         stringsAsFactors = FALSE)
  ColorValues <- c(A = "#440154FF", B = "#21908CFF")
  res <- plotPCA(data, metadata, ColorVariable = "Group", ColorValues = ColorValues, legend_nrow = 1, legend_position = "bottom", pointSize = 6)
  expect_true(is.list(res))
  expect_true(all(c("plt", "data") %in% names(res)))
})
