test_that("IndividualGenes_Violins returns a plot and data frame for basic use", {
  set.seed(1801)
  data <- data.frame(
    A = c(10, 20, 30),
    B = c(5, 15, 25),
    C = c(2, 12, 22)
  )
  rownames(data) <- c("Gene1", "Gene2", "Gene3")
  metadata <- data.frame(
    sample = c("A", "B", "C"),
    Group = c("Control", "Treatment", "Control"),
    stringsAsFactors = FALSE
  )
  genes <- c("Gene1", "Gene2")
  res <- IndividualGenes_Violins(data, metadata, genes, "Group", plot = FALSE)
  expect_true(is.list(res))
  expect_true(all(c("plot", "data") %in% names(res)))
  expect_true(inherits(res$plot, "ggplot"))
  expect_true(is.data.frame(res$data))
  expect_true(all(res$data$gene %in% genes))
})

test_that("IndividualGenes_Violins facets by gene and custom divide variable", {
  set.seed(1802)
  data <- data.frame(
    A = c(10, 20, 30, 40),
    B = c(5, 15, 25, 35),
    C = c(2, 12, 22, 32)
  )
  rownames(data) <- c("Gene1", "Gene2", "Gene3", "Gene4")
  metadata <- data.frame(
    sample = c("A", "B", "C"),
    Group = c("Control", "Treatment", "Control"),
    Block = c("X", "Y", "X"),
    stringsAsFactors = FALSE
  )
  genes <- c("Gene1", "Gene2")
  res <- IndividualGenes_Violins(
    data, metadata, genes, "Group",
    plot = FALSE, divide = "Block", invert_divide = TRUE
  )
  expect_true(is.list(res))
  expect_true(inherits(res$plot, "ggplot"))
})

