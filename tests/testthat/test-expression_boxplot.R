# tests/testthat/test-expressionBoxplot.R

test_that("expressionBoxplot produces correct plot structure", {
  skip_if_not_installed("ggplot2")

  # Test basic plot creation
  p <- expressionBoxplot("HLA-DQA1")

  # Test plot is ggplot object
  expect_s3_class(p, "ggplot")

  # Test required aesthetic mappings
  expect_true("x" %in% names(p$mapping))
  expect_true("y" %in% names(p$mapping))
  expect_true("fill" %in% names(p$mapping))

  # Test plot labels
  expect_true(grepl("HLA-DQA1", p$labels$title))
  expect_equal(p$labels$x, NULL)  # x label should be removed
  expect_equal(p$labels$y, "Expression level")
})

test_that("expressionBoxplot data filtering works correctly", {
  skip_if_not_installed("ggplot2")

  data("example_expression", package = "CeDExplorer")

  # Test that only the specified gene is plotted
  p <- expressionBoxplot("HLA-DQA1", dataset = example_expression)
  plotted_genes <- unique(p$data$gene)
  expect_equal(plotted_genes, "HLA-DQA1")
  expect_equal(length(plotted_genes), 1)
})

test_that("expressionBoxplot plot types work", {
  skip_if_not_installed("ggplot2")

  # Test boxplot (default)
  p_box <- expressionBoxplot("HLA-DQA1", violin = FALSE)
  expect_s3_class(p_box, "ggplot")

  # Test violin plot
  p_violin <- expressionBoxplot("HLA-DQA1", violin = TRUE)
  expect_s3_class(p_violin, "ggplot")

  # Both should have the same data
  expect_equal(nrow(p_box$data), nrow(p_violin$data))
})

test_that("expressionBoxplot handles different conditions", {
  skip_if_not_installed("ggplot2")

  data("example_expression", package = "CeDExplorer")

  p <- expressionBoxplot("HLA-DQA1", dataset = example_expression)

  # Check that both conditions are present
  conditions <- unique(p$data$condition)
  expect_true("CeD" %in% conditions)
  expect_true("Control" %in% conditions)
  expect_equal(length(conditions), 2)
})

test_that("expressionBoxplot handles edge cases", {
  skip_if_not_installed("ggplot2")

  data("example_expression", package = "CeDExplorer")

  # Test with non-existent gene
  expect_warning(
    result <- expressionBoxplot("NONEXISTENT_GENE", dataset = example_expression),
    "not found in dataset"
  )
  expect_null(result)

  # Test error for missing gene argument
  expect_error(
    expressionBoxplot(),
    "Please provide a gene symbol"
  )
})
