library(testthat)
library(CeDExplorer)
library(ggplot2)

create_test_gwas_data <- function() {
  data.frame(
    association_id = paste0("assoc_", 1:30),
    SNP = paste0("rs", 1000 + 1:30),
    mapped_gene = c(
      rep("HLA-DQA1", 5), rep("HLA-DQB1", 5), rep("IL18R1", 5),
      rep("OTHER1", 5), rep("OTHER2", 5), rep(NA, 5)
    ),
    p_value = c(
      10^(-seq(50, 30, length.out = 5)),
      10^(-seq(45, 25, length.out = 5)),
      10^(-seq(20, 10, length.out = 5)),
      10^(-seq(15, 5, length.out = 5)),
      10^(-seq(10, 2, length.out = 5)),
      runif(5, 0.001, 0.1)
    ),
    odds_ratio = runif(30, 1, 3),
    beta = rnorm(30, 0, 0.5),
    standard_error = runif(30, 0.1, 0.5),
    risk_allele = sample(c("A", "T", "C", "G"), 30, replace = TRUE),
    source = "test_data",
    stringsAsFactors = FALSE
  )
}

test_that("plotTopGwasHits returns ggplot object with valid input", {
  test_data <- create_test_gwas_data()
  
  p <- plotTopGwasHits(test_data)
  expect_s3_class(p, "ggplot")
  expect_s3_class(p, "gg")
})

test_that("plotTopGwasHits handles different top_n values", {
  test_data <- create_test_gwas_data()
  
  # Test different top_n values - FIXED parameter name
  p_top10 <- plotTopGwasHits(test_data, top_n = 10)
  expect_s3_class(p_top10, "ggplot")
  
  p_top5 <- plotTopGwasHits(test_data, top_n = 5)
  expect_s3_class(p_top5, "ggplot")
})

test_that("plotTopGwasHits handles custom titles", {
  test_data <- create_test_gwas_data()
  
  p_custom <- plotTopGwasHits(test_data, title = "Custom GWAS Plot")
  expect_s3_class(p_custom, "ggplot")
})

test_that("plotTopGwasHits validates input data", {
  # Test missing data - FIXED error message
  expect_error(plotTopGwasHits(), "gwas_data")
  
  # Test invalid data type - FIXED error message
  expect_error(plotTopGwasHits("invalid"), "invalid")
})

test_that("plotTopGwasHits correctly orders hits by significance", {
  test_data <- create_test_gwas_data()
  
  p <- plotTopGwasHits(test_data, top_n = 10)  # FIXED parameter name
  expect_s3_class(p, "ggplot")
})

test_that("plotTopGwasHits works with real CeDExplorer data", {
  gwas_data <- getCeliacGwas(fetch = FALSE)
  p <- plotTopGwasHits(gwas_data, top_n = 5)  # FIXED parameter name
  expect_s3_class(p, "ggplot")
})
