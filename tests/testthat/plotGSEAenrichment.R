
test_that("plotNESlollipop respects padj coloring and sig_threshold", {
  set.seed(1503)
  padjs <- c(0.001, 0.02, 0.1, 0.6)
  GSEA_results <- list(
    ContrastX = data.frame(
      NES = seq(-2,1,length.out=4),
      padj = padjs,
      pathway = paste("Set", 1:4),
      stat_used = c("B", "B", "B", "B")
    )
  )
  # Should work for padj below and above threshold, and with default colors
  res <- plotNESlollipop(GSEA_results, sig_threshold = 0.05)
  expect_true(is.list(res))
  expect_true(all(vapply(res, inherits, logical(1), "ggplot")))
})


test_that("plotNESlollipop handles a single pathway", {
  set.seed(1505)
  GSEA_results <- list(
    One = data.frame(
      NES = 2.1,
      padj = 0.005,
      pathway = "SinglePathway",
      stat_used = "t"
    )
  )
  res <- plotNESlollipop(GSEA_results)
  expect_true(is.list(res))
  expect_equal(length(res), 1)
  expect_true(all(vapply(res, inherits, logical(1), "ggplot")))
})
