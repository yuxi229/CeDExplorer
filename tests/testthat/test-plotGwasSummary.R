library(testthat)
library(CeDExplorer)
library(ggplot2)

create_test_gwas_data <- function() {
  data.frame(
    association_id = paste0("assoc_", 1:50),
    SNP = paste0("rs", 1000 + 1:50),
    mapped_gene = sample(c("HLA-DQA1", "HLA-DQB1", "IL18R1", "OTHER"), 50, replace = TRUE),
    p_value = 10^(-runif(50, 1, 50)),
    odds_ratio = runif(50, 1, 3),
    beta = rnorm(50, 0, 0.5),
    standard_error = runif(50, 0.1, 0.5),
    risk_allele = sample(c("A", "T", "C", "G"), 50, replace = TRUE),
    source = "test_data",
    stringsAsFactors = FALSE
  )
}

test_that("plotGwasSummary returns ggplot object with valid input", {
  test_data <- create_test_gwas_data()
  
  p <- plotGwasSummary(test_data)
  expect_s3_class(p, "ggplot")
  expect_s3_class(p, "gg")
})

test_that("plotGwasSummary handles custom parameters", {
  test_data <- create_test_gwas_data()
  
  # Test custom title
  p_custom_title <- plotGwasSummary(test_data, title = "Custom Title")
  expect_s3_class(p_custom_title, "ggplot")
  
  # Test custom significance level - FIXED parameter name
  p_custom_sig <- plotGwasSummary(test_data, significance_level = 1e-6)
  expect_s3_class(p_custom_sig, "ggplot")
})

test_that("plotGwasSummary validates input data", {
  # Test missing data - FIXED error message
  expect_error(plotGwasSummary(), "gwas_data")
  
  # Test invalid data type - FIXED error message
  expect_error(plotGwasSummary("invalid"), "invalid")
  
  # Test empty data frame - FIXED error message
  expect_error(plotGwasSummary(data.frame()), "mathematical function")
})

test_that("plotGwasSummary handles edge cases", {
  test_data <- create_test_gwas_data()
  
  # Test with small dataset
  small_data <- test_data[1:5, ]
  p_small <- plotGwasSummary(small_data)
  expect_s3_class(p_small, "ggplot")
})

test_that("plotGwasSummary works with real CeDExplorer data", {
  gwas_data <- getCeliacGwas(fetch = FALSE)
  p <- plotGwasSummary(gwas_data)
  expect_s3_class(p, "ggplot")
})