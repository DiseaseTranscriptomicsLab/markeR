test_that("plotNESlollipop returns list of ggplot objects for single contrast", {
  gsea_mock <- list(
    "A" = data.frame(
      pathway = c("Pathway 1", "Pathway 2"),
      NES = c(1.2, -1.5),
      padj = c(0.01, 0.2),
      stat_used = c("t", "B")
    )
  )
  result <- plotNESlollipop(gsea_mock)
  expect_type(result, "list")
  expect_true(all(vapply(result, function(x) inherits(x, "gg"), logical(1))))
})

test_that("plotNESlollipop does not fail when all padj > sig_threshold", {
  gsea_mock <- list(
    "A" = data.frame(
      pathway = c("Pathway 1", "Pathway 2"),
      NES = c(1, -1),
      padj = c(0.2, 0.4),
      stat_used = c("t", "B")
    )
  )
  expect_silent(plotNESlollipop(gsea_mock))
})

test_that("plotNESlollipop handles empty input gracefully", {
  gsea_mock <- list()
  expect_silent(result <- plotNESlollipop(gsea_mock))
  expect_length(result, 0)
})
