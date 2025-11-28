# tests/testthat/test_gwas.R
context("GWAS Functions")

test_that("getCeliacGwas returns data frame with expected columns", {
  # Test with bundled data (always works)
  data <- getCeliacGwas(fetch = FALSE)
  
  expect_s3_class(data, "data.frame")
  expect_true(nrow(data) > 0)
  
  # Check for essential columns
  expected_cols <- c("association_id", "SNP", "mapped_gene", "p_value", "source")
  expect_true(all(expected_cols %in% names(data)))
})

test_that("getCeliacGwas handles fetch = TRUE gracefully", {
  # Test that it doesn't crash with fetch = TRUE
  # It should either return real data or fall back to bundled data
  data <- getCeliacGwas(fetch = TRUE)
  
  expect_s3_class(data, "data.frame")
  expect_true(nrow(data) > 0)
  expect_true("p_value" %in% names(data))
})

test_that("plotTopGwasHits creates ggplot object", {
  data <- getCeliacGwas(fetch = FALSE)
  
  # Test basic plot
  plot <- plotTopGwasHits(data, top_n = 5)
  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot, "gg")
  
  # Test with different top_n values
  plot_10 <- plotTopGwasHits(data, top_n = 10)
  expect_s3_class(plot_10, "ggplot")
})

test_that("plotTopGwasHits handles edge cases", {
  data <- getCeliacGwas(fetch = FALSE)
  
  # Test with very small top_n
  plot <- plotTopGwasHits(data, top_n = 2)
  expect_s3_class(plot, "ggplot")
  
  # Test with top_n larger than data
  plot <- plotTopGwasHits(data, top_n = 100)
  expect_s3_class(plot, "ggplot")
})

test_that("plotGwasSummary creates ggplot object", {
  data <- getCeliacGwas(fetch = FALSE)
  
  plot <- plotGwasSummary(data)
  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot, "gg")
})

test_that("GWAS data has reasonable p-values", {
  data <- getCeliacGwas(fetch = FALSE)
  
  # Check p-values are numeric and in reasonable range
  expect_type(data$p_value, "double")
  expect_true(all(data$p_value >= 0 & data$p_value <= 1))
  
  # Check -log10 p-values can be calculated
  log_p <- -log10(data$p_value)
  expect_true(all(is.finite(log_p) | is.na(log_p)))
})

test_that("GWAS functions work with real data structure", {
  # Simulate the structure we get from GWAS Catalog
  mock_gwas_data <- data.frame(
    association_id = paste0("test_", 1:10),
    SNP = paste0("rs", 1000 + 1:10),
    mapped_gene = c("HLA-DQA1", "HLA-DQB1", "IL2", "IL21", "SH2B3", 
                    "CCR1", "CCR3", "TAGAP", "PTPN2", "CTLA4"),
    p_value = 10^(-sample(5:50, 10)),
    odds_ratio = runif(10, 1.0, 3.0),
    source = "test_data",
    stringsAsFactors = FALSE
  )
  
  # Test that plotting functions work with this structure
  plot1 <- plotTopGwasHits(mock_gwas_data)
  plot2 <- plotGwasSummary(mock_gwas_data)
  
  expect_s3_class(plot1, "ggplot")
  expect_s3_class(plot2, "ggplot")
})
