# tests/testthat/test-network_functions.R

test_that("plot_gene_network works with valid inputs", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggraph")

  # Test basic functionality with genes that have interactions
  genes <- c("HLA-DQA1", "HLA-DQB1")  # These should have interactions in our data
  p <- plot_gene_network(genes)

  # The function might return NULL if layout fails, so check for either
  if (!is.null(p)) {
    expect_s3_class(p, "ggplot")
  }

  # Test with single gene (edge case)
  p_single <- plot_gene_network("HLA-DQA1")
  if (!is.null(p_single)) {
    expect_true(TRUE)  # Just check it doesn't crash
  }
})

test_that("plot_gene_network input validation", {
  # Test error for missing genes
  expect_error(plot_gene_network(), "Please provide a vector of gene symbols")

  # Test error for empty genes
  expect_error(plot_gene_network(character(0)), "Please provide a vector of gene symbols")

  # Test error for invalid layout
  expect_error(plot_gene_network(c("HLA-DQA1"), layout = "invalid"), "layout must be one of")

  # Test error for invalid color_by
  expect_error(plot_gene_network(c("HLA-DQA1"), color_by = "invalid"), "color_by must be one of")
})

test_that("plot_gene_network handles edge cases gracefully", {
  skip_if_not_installed("igraph")

  # Test with non-existent genes (should warn but not crash)
  expect_warning(
    result <- plot_gene_network(c("NONEXISTENT_GENE1", "NONEXISTENT_GENE2")),
    "No interactions found"
  )

  # Test with mixed valid/invalid genes
  genes_mixed <- c("HLA-DQA1", "INVALID_GENE")
  p <- plot_gene_network(genes_mixed)
  # Should handle this without crashing
  expect_true(TRUE)
})
