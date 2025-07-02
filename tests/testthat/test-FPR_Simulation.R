library(devtools)

test_that("FPR_Simulation returns FPR=1 when simulated and original signatures are identical", {
  set.seed(42)
  expr <- as.data.frame(matrix(abs(rnorm(40)), nrow = 4, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:4)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Condition = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )
  sigs <- list(Sig1 = c("Gene1", "Gene2", "Gene3"))
  # gene_list is exactly the signature, so all simulations sample the same genes as the original
  result <- FPR_Simulation(
    data = expr,
    metadata = metadata,
    original_signatures = sigs,
    Variable = "Condition",
    gene_list = c("Gene1", "Gene2", "Gene3"),
    number_of_sims = 10,
    title = "Identical Sims"
  )
  # Extract the per-signature data frame
  df <- result$data$Sig1
  # FPR is only computed for Original rows
  fpr_vals <- df$FPR[df$type == "Original"]
  expect_true(all(abs(fpr_vals - 1) < 1e-8))
})





test_that("FPR_Simulation: FPR is between 0 and 1 for random gene sets", {
  set.seed(324)
  expr <- as.data.frame(matrix(abs(rnorm(40)), nrow = 4, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:4)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Condition = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )
  sigs <- list(SigA = c("Gene1", "Gene2", "Gene3"))
  result <- FPR_Simulation(
    data = expr,
    metadata = metadata,
    original_signatures = sigs,
    Variable = "Condition",
    number_of_sims = 12
  )
  df <- result$data$SigA
  fpr_vals <- df$FPR[df$type == "Original"]
  expect_true(all(fpr_vals >= 0 & fpr_vals <= 1))
})

test_that("FPR_Simulation handles multiple signatures and returns one plot per signature", {
  set.seed(325)
  expr <- as.data.frame(matrix(abs(rnorm(60)), nrow = 6, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Condition = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )
  sigs <- list(
    Sig1 = c("Gene1", "Gene2"),
    Sig2 = c("Gene3", "Gene4")
  )
  result <- FPR_Simulation(
    data = expr,
    metadata = metadata,
    original_signatures = sigs,
    Variable = "Condition",
    number_of_sims = 9
  )
  expect_true(length(result$data) == length(sigs))
})

test_that("FPR_Simulation: Simulated rows count matches number_of_sims times methods times contrasts", {
  set.seed(326)
  expr <- as.data.frame(matrix(abs(rnorm(40)), nrow = 4, ncol = 10))
  rownames(expr) <- paste0("Gene", 1:4)
  colnames(expr) <- paste0("Sample", 1:10)
  metadata <- data.frame(
    sample = colnames(expr),
    Condition = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )
  sigs <- list(SigX = c("Gene1", "Gene2"))
  n_sims <- 7
  result <- FPR_Simulation(
    data = expr,
    metadata = metadata,
    original_signatures = sigs,
    Variable = "Condition",
    number_of_sims = n_sims
  )
  df <- result$data$SigX
  n_methods <- length(unique(df$method))
  n_contrasts <- length(unique(df$contrast))
  n_sim <- sum(df$type == "Simulated")
  expect_equal(n_sim, n_sims * n_methods * n_contrasts)
})

