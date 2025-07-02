test_that("plotVolcano creates a basic volcano plot without errors", {
  set.seed(101)
  expr <- matrix(rpois(100, lambda = 20), nrow = 10, ncol = 10)
  rownames(expr) <- paste0("gene", 1:10)
  colnames(expr) <- paste0("sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )
  de_res <- calculateDE(
    data = expr,
    metadata = metadata,
    variables = "Group",
    contrasts = "A-B"
  )
  p <- plotVolcano(
    DEResultsList = de_res,
    genes = NULL,
    x = "logFC",
    y = "-log10(adj.P.Val)",
    pointSize = 2,
    title = "Test Volcano"
  )
  expect_true(is_ggplot(p) || inherits(p, "gtable"))
})


test_that("plotVolcano supports multiple contrasts and signatures", {
  set.seed(104)
  expr <- matrix(rpois(200, lambda = 20), nrow = 20, ncol = 10)
  rownames(expr) <- paste0("gene", 1:20)
  colnames(expr) <- paste0("sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Group = rep(c("A", "B"), each = 5),
    Batch = rep(c("X", "Y"), times = 5),
    stringsAsFactors = FALSE
  )
  de_res <- calculateDE(
    data = expr,
    metadata = metadata,
    variables = "Group",
    contrasts = c("A-B","B-A")
  )
  sig1 <- rownames(de_res[["A-B"]])[1:2]
  sig2 <- rownames(de_res[["A-B"]])[3:4]
  p <- plotVolcano(
    DEResultsList = de_res,
    genes = list(Sig1 = sig1, Sig2 = sig2),
    x = "logFC",
    y = "-log10(adj.P.Val)",
    title = "Grid Test"
  )
  expect_true(is_ggplot(p) || inherits(p, "gtable"))
})
