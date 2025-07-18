test_that("VisualiseIndividualGenes dispatches correctly for each type", {
  set.seed(123)
  expr_data <- matrix(rexp(1000, rate = 1), nrow = 50, ncol = 20)
  rownames(expr_data) <- paste0("Gene", 1:50)
  colnames(expr_data) <- paste0("Sample", 1:20)

  sample_info <- data.frame(
    SampleID = colnames(expr_data),
    Condition = rep(c("A", "B"), each = 10),
    Diagnosis = rep(c("Disease", "Control"), each = 10),
    stringsAsFactors = FALSE
  )
  rownames(sample_info) <- sample_info$SampleID
  selected_genes <- rownames(expr_data)[1:5]

  # Minimal mocks (replace with real tests or mocks using mockery or testthat::with_mock)
  mock_return <- ggplot2::ggplot() + ggplot2::geom_blank()
  local_mocked_bindings(
    IndividualGenes_Violins = function(...) mock_return,
    CorrelationHeatmap = function(...) mock_return,
    ExpressionHeatmap = function(...) mock_return,
    ROCandAUCplot = function(...) mock_return,
    CohenD_IndividualGenes = function(...) mock_return,
    plotPCA = function(...) mock_return
  )

  expect_s3_class(
    VisualiseIndividualGenes(
      type = "violin",
      data = expr_data,
      metadata = sample_info,
      genes = selected_genes,
      GroupingVariable = "Condition"
    ),
    "ggplot"
  )

  expect_s3_class(
    VisualiseIndividualGenes(
      type = "correlation",
      data = expr_data,
      genes = selected_genes
    ),
    "ggplot"
  )

  expect_s3_class(
    VisualiseIndividualGenes(
      type = "expression",
      data = expr_data,
      genes = selected_genes
    ),
    "ggplot"
  )

  expect_s3_class(
    VisualiseIndividualGenes(
      type = "pca",
      data = expr_data,
      genes = selected_genes
    ),
    "ggplot"
  )

  for (ptype in c("roc", "auc", "rocauc", "cohend")) {
    expect_s3_class(
      VisualiseIndividualGenes(
        type = ptype,
        data = expr_data,
        metadata = sample_info,
        genes = selected_genes,
        condition_var = "Diagnosis",
        class = "Disease"
      ),
      "ggplot"
    )
  }
})

test_that("VisualiseIndividualGenes throws errors for missing required arguments", {
  expr_data <- matrix(1, nrow = 5, ncol = 2)
  colnames(expr_data) <- c("S1", "S2")
  rownames(expr_data) <- paste0("Gene", 1:5)

  metadata <- data.frame(
    SampleID = c("S1", "S2"),
    Condition = c("A", "B"),
    Diagnosis = c("Disease", "Control"),
    row.names = c("S1", "S2"),
    stringsAsFactors = FALSE
  )


  # Missing genes
  expect_error(VisualiseIndividualGenes(type = "violin", data = expr_data),
               "Argument 'genes' is required")

  # Missing metadata for violin
  expect_error(VisualiseIndividualGenes(type = "violin", data = expr_data, genes = c("Gene1")),
               "metadata")

  # Missing GroupingVariable
  expect_error(
    VisualiseIndividualGenes(
      type = "violin",
      data = expr_data,
      metadata = metadata,
      genes = c("Gene1")
    ),
    "GroupingVariable"
  )

  # Missing condition_var or class for roc
  expect_error(
    VisualiseIndividualGenes(
      type = "roc",
      data = expr_data,
      metadata = metadata,
      genes = c("Gene1"),
      class = "Disease"
    ),
    "condition_var"
  )

  expect_error(
    VisualiseIndividualGenes(
      type = "roc",
      data = expr_data,
      metadata = metadata,
      genes = c("Gene1"),
      condition_var = "Diagnosis"
    ),
    "class"
  )

  # Invalid plot type
  expect_error(
    VisualiseIndividualGenes(
      type = "banana",
      data = expr_data,
      genes = c("Gene1")
    ),
    "Invalid type"
  )
})

