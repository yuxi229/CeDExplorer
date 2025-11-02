# Test network functions
test_that('plotGeneNetwork works with valid inputs', {
  # Use genes that are likely to have interactions
  genes <- c('HLA-DQA1', 'HLA-DQB1', 'IL2', 'STAT1')
  
  # Test that function runs without error
  expect_error(plot <- plotGeneNetwork(genes), NA)
  
  # If plot is created, check it's a ggplot
  if(!is.null(plot)) {
    expect_s3_class(plot, 'ggplot')
  }
})

test_that('plotGeneNetwork handles single genes', {
  # Single gene should handle gracefully
  expect_error(plot <- plotGeneNetwork('HLA-DQA1'), NA)
})

test_that('plotGeneNetwork validates input correctly', {
  # Empty input should give informative error
  expect_error(
    plotGeneNetwork(character(0)),
    "Please provide a vector of gene symbols"
  )
  
  # Non-existent genes should handle gracefully (no error)
  expect_error(plotGeneNetwork(c('FAKE_GENE_1', 'FAKE_GENE_2')), NA)
})

# [END]
