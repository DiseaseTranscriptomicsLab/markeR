library(testthat)
library(markeR)  # Replace with your package name if different

# A dummy test to ensure testing works
test_that("Basic math works", {
  expect_equal(1 + 1, 2)
})

# Test for calculateScore_logmedian_unidirectional function
test_that("calculateScore_logmedian_unidirectional returns a data frame with expected columns", {
  # Create a small example expression matrix with 2 genes and 3 samples
  data <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
  rownames(data) <- c("gene1", "gene2")
  colnames(data) <- c("sample1", "sample2", "sample3")

  # Define a gene signature using one gene name
  signature <- "gene1"

  # Call the function from your package
  result <- calculateScore_logmedian_unidirectional(data, signature)

  # Verify that the result is a data frame with the expected columns
  expect_true(is.data.frame(result))
  expect_true(all(c("sample", "score") %in% colnames(result)))
})
