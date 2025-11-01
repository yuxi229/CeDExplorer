# tests/testthat/test-dashboard.R

test_that("geneEvidenceView works with valid inputs", {
  skip_if_not_installed("ggplot2")

  # Test basic functionality
  result <- geneEvidenceView("HLA-DQA1")

  # Could be a combined plot or list
  if (requireNamespace("patchwork", quietly = TRUE)) {
    expect_s3_class(result, "gg")
  } else {
    expect_type(result, "list")
  }
})

test_that("geneEvidenceView input validation", {
  # Test error for missing gene
  expect_error(geneEvidenceView(), "Please provide a gene symbol")

  # Test error for empty gene
  expect_error(geneEvidenceView(character(0)), "Please provide a gene symbol")
})

test_that("geneEvidenceView handles different evidence types", {
  skip_if_not_installed("ggplot2")

  # Test with only expression
  result1 <- geneEvidenceView("HLA-DQA1", show_gwas = FALSE, show_network = FALSE, show_hla = FALSE)

  # Test with only network
  result2 <- geneEvidenceView("HLA-DQB1", show_expression = FALSE, show_gwas = FALSE, show_hla = FALSE)

  # Should handle these without crashing
  expect_true(TRUE)
})
