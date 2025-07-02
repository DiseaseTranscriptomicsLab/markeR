test_that("runGSEA runs and returns expected structure for vector gene set", {
  set.seed(1301)
  deg <- data.frame(
    t = rnorm(50),
    B = rnorm(50),
    row.names = paste0("Gene", 1:50)
  )
  DEGList <- list(Contrast1 = deg)
  gene_sets <- list(Set1 = c("Gene1", "Gene20", "Gene30"))
  res <- runGSEA(DEGList, gene_sets, nPermSimple = 100)
  expect_true(is.list(res))
  expect_true("Contrast1" %in% names(res))
  expect_true(is.data.frame(res$Contrast1))
  expect_true("padj" %in% names(res$Contrast1))
  expect_true("stat_used" %in% names(res$Contrast1))
})

test_that("runGSEA works with bidirectional data.frame gene set", {
  set.seed(1302)
  deg <- data.frame(
    t = rnorm(50),
    B = rnorm(50),
    row.names = paste0("Gene", 1:50)
  )
  DEGList <- list(ContrastA = deg)
  gs_df <- data.frame(Gene = c("Gene2", "Gene10", "Gene15"),
                      Direction = c(1, -1, 1))
  gene_sets <- list(Bidir = gs_df)
  res <- runGSEA(DEGList, gene_sets, nPermSimple = 100)
  expect_true(is.data.frame(res$ContrastA))
  expect_equal(unique(res$ContrastA$stat_used), "t")
})

test_that("runGSEA applies ContrastCorrection for multiple contrasts", {
  set.seed(1303)
  deg1 <- data.frame(t = rnorm(50), B = rnorm(50), row.names = paste0("Gene", 1:50))
  deg2 <- data.frame(t = rnorm(50), B = rnorm(50), row.names = paste0("Gene", 1:50))
  DEGList <- list(C1 = deg1, C2 = deg2)
  gene_sets <- list(Set1 = c("Gene1", "Gene20", "Gene30"))
  res1 <- runGSEA(DEGList, gene_sets, nPermSimple = 100, ContrastCorrection = FALSE)
  res2 <- runGSEA(DEGList, gene_sets, nPermSimple = 100, ContrastCorrection = TRUE)
  # padj should be different if correction is applied across contrasts
  expect_true(any(res1$C1$padj != res2$C1$padj))
  expect_true(any(res1$C2$padj != res2$C2$padj))
})

test_that("runGSEA uses user-specified stat argument", {
  set.seed(1304)
  deg <- data.frame(
    t = rnorm(50),
    B = rnorm(50),
    row.names = paste0("Gene", 1:50)
  )
  DEGList <- list(Contr = deg)
  gene_sets <- list(GS = c("Gene1", "Gene20", "Gene30"))
  res <- runGSEA(DEGList, gene_sets, stat = "B", nPermSimple = 100)
  expect_equal(unique(res$Contr$stat_used), "B")
})
