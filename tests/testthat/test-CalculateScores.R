test_that("CalculateScores runs without error", {
  # Very simple test
  data(metadata_example)
  data(counts_example)

  expect_silent(CalculateScores(data = counts_example,
                                metadata = metadata_example,
                                method = "logmedian",
                                gene_sets = list(Senescence=SimpleSenescenceSignature)))


})
