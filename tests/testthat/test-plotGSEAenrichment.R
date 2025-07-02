test_that("plotGSEAenrichment works with bidirectional gene set", {
  set.seed(1402)
  GSEA_results <- list(
    Contrast1 = data.frame(
      NES = -2.1,
      padj = 0.01,
      pathway = "SigB",
      stat_used = "t"
    )
  )
  DEGList <- list(
    Contrast1 = data.frame(
      t = rnorm(40),
      B = rnorm(40),
      row.names = paste0("Gene", 1:40)
    )
  )
  gene_sets <- list(
    SigB = data.frame(Gene = paste0("Gene", 1:10), Direction = sample(c(1, -1), 10, replace = TRUE))
  )
  res <- plotGSEAenrichment(GSEA_results, DEGList, gene_sets)
  expect_true(is.list(res))
  expect_true(all(vapply(res, inherits, logical(1), "ggplot")))
})

test_that("plotGSEAenrichment arranges plots in a grid when grid=TRUE", {
  set.seed(1403)
  GSEA_results <- list(
    Contrast1 = data.frame(
      NES = rnorm(2),
      padj = runif(2),
      pathway = c("S1", "S2"),
      stat_used = rep("t", 2)
    ),
    Contrast2 = data.frame(
      NES = rnorm(2),
      padj = runif(2),
      pathway = c("S1", "S2"),
      stat_used = rep("t", 2)
    )
  )
  DEGList <- list(
    Contrast1 = data.frame(
      t = rnorm(50),
      B = rnorm(50),
      row.names = paste0("Gene", 1:50)
    ),
    Contrast2 = data.frame(
      t = rnorm(50),
      B = rnorm(50),
      row.names = paste0("Gene", 1:50)
    )
  )
  gene_sets <- list(
    S1 = paste0("Gene", 1:20),
    S2 = paste0("Gene", 21:40)
  )
  res <- plotGSEAenrichment(GSEA_results, DEGList, gene_sets, grid = TRUE)
  expect_true(inherits(res, "gtable") || inherits(res, "ggarrange"))
})

test_that("plotGSEAenrichment skips missing gene sets", {
  set.seed(1404)
  GSEA_results <- list(
    Contrast1 = data.frame(
      NES = 1.3,
      padj = 0.04,
      pathway = "NotInGeneSets",
      stat_used = "t"
    )
  )
  DEGList <- list(
    Contrast1 = data.frame(
      t = rnorm(20),
      B = rnorm(20),
      row.names = paste0("Gene", 1:20)
    )
  )
  gene_sets <- list(
    SomeOtherSet = paste0("Gene", 1:10)
  )
  res <- plotGSEAenrichment(GSEA_results, DEGList, gene_sets)
  expect_true(is.list(res))
  expect_equal(length(res), 0)
})

test_that("plotGSEAenrichment works with multiple contrasts and signatures", {
  set.seed(1405)
  GSEA_results <- list(
    Contrast1 = data.frame(
      NES = rnorm(2),
      padj = runif(2),
      pathway = c("A", "B"),
      stat_used = rep("t", 2)
    ),
    Contrast2 = data.frame(
      NES = rnorm(2),
      padj = runif(2),
      pathway = c("A", "B"),
      stat_used = rep("t", 2)
    )
  )
  DEGList <- list(
    Contrast1 = data.frame(
      t = rnorm(30),
      B = rnorm(30),
      row.names = paste0("Gene", 1:30)
    ),
    Contrast2 = data.frame(
      t = rnorm(30),
      B = rnorm(30),
      row.names = paste0("Gene", 1:30)
    )
  )
  gene_sets <- list(
    A = paste0("Gene", 1:15),
    B = paste0("Gene", 16:30)
  )
  res <- plotGSEAenrichment(GSEA_results, DEGList, gene_sets)
  expect_true(is.list(res))
  expect_equal(length(res), 4)
  expect_true(all(vapply(res, inherits, logical(1), "ggplot")))
})
