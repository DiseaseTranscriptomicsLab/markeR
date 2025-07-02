test_that("plotCombinedGSEA returns a ggplot for standard input", {
  set.seed(1601)
  GSEA_results <- list(
    Contrast1 = data.frame(
      NES = rnorm(3),
      padj = c(0.01, 0.2, 0.04),
      pathway = paste("Pathway", 1:3),
      stat_used = c("t", "B", "B")
    ),
    Contrast2 = data.frame(
      NES = rnorm(3),
      padj = c(0.03, 0.06, 0.5),
      pathway = paste("Pathway", 4:6),
      stat_used = c("t", "B", "B")
    )
  )
  p <- plotCombinedGSEA(GSEA_results)
  expect_true(is_ggplot(p))
})


test_that("plotCombinedGSEA handles all pathways below and above threshold", {
  set.seed(1604)
  GSEA_results <- list(
    C1 = data.frame(
      NES = rnorm(3),
      padj = c(0.001, 0.01, 0.02),
      pathway = paste("P", 1:3),
      stat_used = rep("t", 3)
    ),
    C2 = data.frame(
      NES = rnorm(3),
      padj = c(0.9, 0.8, 0.7),
      pathway = paste("Q", 1:3),
      stat_used = rep("t", 3)
    )
  )
  p <- plotCombinedGSEA(GSEA_results, sig_threshold = 0.05)
  expect_true(is_ggplot(p))
})

