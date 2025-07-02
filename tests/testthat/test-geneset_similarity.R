library(devtools)

test_that("geneset_similarity computes expected Jaccard index", {
  set1 <- c("GeneA", "GeneB", "GeneC")
  set2 <- c("GeneB", "GeneC", "GeneD")
  sim <-geneset_similarity(signatures = list(set1=set1), other_user_signatures = list(set2=set2), metric="jaccard")
  sim <- sim$data$Score
  expect_gt(sim, 0)
  expect_lt(sim, 1)
})

test_that("geneset_similarity returns jaccard 1 for identical sets", {
  gs <- c("GeneA", "GeneB")
  sim <-geneset_similarity(signatures = list(set1=gs), other_user_signatures = list(set2=gs), metric="jaccard")
  sim <- sim$data$Score
  expect_equal(sim, 1)
})

test_that("geneset_similarity pval_threshold filters labels as expected for odds_ratio", {
  # Fake simple signatures to ensure overlap
  sig1 <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5")
  sig2 <- c("GENE1", "GENE2", "GENE6", "GENE7", "GENE8")
  signatures <- list(A = sig1)
  others <- list(B = sig2)
  # Large universe to ensure non-perfect overlap
  universe <- toupper(paste0("GENE", 1:50))

  # Run with a permissive pval_threshold (should label if pval < threshold)
  res <- geneset_similarity(
    signatures = signatures,
    other_user_signatures = others,
    metric = "odds_ratio",
    universe = universe,
    pval_threshold = 1 # allow all p-values
  )
  d <- res$data
  expect_true(any(d$Label != "")) # Some label should be shown

  # Run with a restrictive pval_threshold (should blank out most labels)
  res2 <- geneset_similarity(
    signatures = signatures,
    other_user_signatures = others,
    metric = "odds_ratio",
    universe = universe,
    pval_threshold = 1e-10 # extremely small, almost nothing will label
  )
  d2 <- res2$data
  expect_true(all(d2$Label == "")) # All labels should be blank
})


test_that("geneset_similarity returns expected odds ratio values", {
  # Simple signatures with known overlap
  sig1 <- c("GENE1", "GENE2", "GENE3", "GENE4")
  sig2 <- c("GENE1", "GENE2", "GENE5", "GENE6")
  signatures <- list(A = sig1)
  others <- list(B = sig2)
  # Define a universe that is the union of all genes plus some extras
  universe <- toupper(c(sig1, sig2, "GENE7", "GENE8", "GENE9", "GENE10"))

  # Compute expected contingency table
  # a = overlap (GENE1, GENE2): 2
  # b = sig1 only (GENE3, GENE4): 2
  # c = sig2 only (GENE5, GENE6): 2
  # d = universe - union: (GENE7, GENE8, GENE9, GENE10): 4

  cont_tbl <- matrix(c(2, 2, 2, 4), nrow = 2)
  fisher_res <- fisher.test(cont_tbl)
  expected_or <- as.numeric(fisher_res$estimate)
  expected_log10or <- log10(expected_or)

  # Run the function
  res <- geneset_similarity(
    signatures = signatures,
    other_user_signatures = others,
    metric = "odds_ratio",
    universe = universe
  )
  d <- res$data

  # Check the actual odds ratio (on log10 scale) is close to expected
  expect_equal(d$Score, expected_log10or, tolerance = 1e-6)
})

test_that("geneset_similarity returns expected Jaccard index values with H collection", {
  # Use a minimal set that will be in H hallmark for testing
  sig1 <- c("TP53", "BRCA1", "MYC", "EGFR", "CDK2")         # from example
  signatures <- list(A = sig1)

  # Get hallmark (H) gene sets using msigdbr
  gs <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  hallmark_sets <- split(toupper(gs$gene_symbol), gs$gs_name)

  # Pick a hallmark set and compute expected Jaccard index
  # Let's use "HALLMARK_MYC_TARGETS_V1" if present
  if (!"HALLMARK_MYC_TARGETS_V1" %in% names(hallmark_sets)) {
    skip("HALLMARK_MYC_TARGETS_V1 not present in msigdbr::msigdbr() output.")
  }
  h_set <- hallmark_sets[["HALLMARK_MYC_TARGETS_V1"]]
  expected_jaccard <- length(intersect(sig1, h_set)) / length(union(sig1, h_set))

  # Run the function, using only this hallmark as msig_subset to keep the test fast
  res <- geneset_similarity(
    signatures = signatures,
    metric = "jaccard",
    collection = "H",
    msig_subset = "HALLMARK_MYC_TARGETS_V1"
  )
  d <- res$data
  # Find the row for this comparison
  idx <- which(d$Compared_Signature == "HALLMARK_MYC_TARGETS_V1")
  expect_true(length(idx) == 1)
  actual_jaccard <- d$Score[idx]
  expect_equal(actual_jaccard, expected_jaccard, tolerance = 1e-10)
})
