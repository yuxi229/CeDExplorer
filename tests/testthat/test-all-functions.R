# tests/testthat/test-all-functions.R

test_that("all exported functions work", {
  skip_on_cran()
  
  # Test data loading
  expect_s3_class(get_celiac_gwas(fetch = FALSE), "data.frame")
  
  # Test plotting functions return ggplot objects
  gwas_data <- get_celiac_gwas(fetch = FALSE)
  p1 <- plot_manhattan(gwas_data)
  if (!is.null(p1)) expect_s3_class(p1, "ggplot")
  
  p2 <- expression_boxplot("HLA-DQA1")
  if (!is.null(p2)) expect_s3_class(p2, "ggplot")
  
  # Test network function
  p3 <- plot_gene_network(c("HLA-DQA1", "HLA-DQB1"))
  if (!is.null(p3)) expect_s3_class(p3, "ggplot")
  
  # Test HLA function
  p4 <- hla_interaction_plot("HLA-DQA1")
  expect_s3_class(p4, "ggplot")
  
  # Test dashboard function
  result <- geneEvidenceView("HLA-DQA1")
  expect_true(!is.null(result))
})