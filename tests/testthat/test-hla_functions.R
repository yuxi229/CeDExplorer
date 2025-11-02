# tests/testthat/test-hla_functions.R

test_that("hlaInteractionPlot works with valid inputs", {
  skip_if_not_installed("ggplot2")

  # Test basic functionality
  p <- hlaInteractionPlot("HLA-DQA1")
  expect_s3_class(p, "ggplot")

  # Test with different parameters
  p2 <- hlaInteractionPlot("HLA-DQB1", point_size = 3, add_trend = TRUE)
  expect_s3_class(p2, "ggplot")
})

test_that("hlaInteractionPlot input validation", {
  # Test error for missing gene
  expect_error(hlaInteractionPlot(), "Please provide a gene symbol")

  # Test error for invalid gene (now checking for the correct error message)
  expect_error(
    hlaInteractionPlot("INVALID_GENE"),
    "Dataset missing required columns: INVALID_GENE"
  )

  # Test with a gene that exists in dataset but not in allele frequencies
  # We need to create a custom test case for this
  skip_if_not_installed("ggplot2")

  # Create custom dataset with an extra gene
  custom_data <- data.frame(
    sample = c("S1", "S2"),
    condition = c("CeD", "Control"),
    `EXTRA_GENE` = c(5, 6),
    check.names = FALSE
  )

  # This should fail because EXTRA_GENE is not in allele frequency data
  expect_error(
    hlaInteractionPlot("EXTRA_GENE", dataset = custom_data),
    "not found in allele frequency data"
  )
})

test_that("hlaInteractionPlot handles edge cases", {
  skip_if_not_installed("ggplot2")

  # Test with empty dataset
  empty_data <- data.frame(sample = character(0), condition = character(0))
  expect_error(
    hlaInteractionPlot("HLA-DQA1", dataset = empty_data),
    "Dataset missing required columns"
  )

  # Test with invalid allele frequency data
  invalid_freq <- data.frame(wrong_col = "test")
  expect_error(
    hlaInteractionPlot("HLA-DQA1", allele_freq = invalid_freq),
    "allele_freq must contain columns: allele, frequency, condition"
  )
})
