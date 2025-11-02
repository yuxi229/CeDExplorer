# tests/testthat/test_gwas.R
context("GWAS Functions")

test_that("get_celiac_gwas returns data frame with expected columns", {
  # Test with bundled data (always works)
  data <- get_celiac_gwas(fetch = FALSE)
  
  expect_s3_class(data, "data.frame")
  expect_true(nrow(data) > 0)
  
  # Check for essential columns
  expected_cols <- c("association_id", "SNP", "mapped_gene", "p_value", "source")
  expect_true(all(expected_cols %in% names(data)))
})

test_that("get_celiac_gwas handles fetch = TRUE gracefully", {
  # Test that it doesn't crash with fetch = TRUE
  # It should either return real data or fall back to bundled data
  data <- get_celiac_gwas(fetch = TRUE)
  
  expect_s3_class(data, "data.frame")
  expect_true(nrow(data) > 0)
  expect_true("p_value" %in% names(data))
})

test_that("plot_top_gwas_hits creates ggplot object", {
  data <- get_celiac_gwas(fetch = FALSE)
  
  # Test basic plot
  plot <- plot_top_gwas_hits(data, top_n = 5)
  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot, "gg")
  
  # Test with different top_n values
  plot_10 <- plot_top_gwas_hits(data, top_n = 10)
  expect_s3_class(plot_10, "ggplot")
})

test_that("plot_top_gwas_hits handles edge cases", {
  data <- get_celiac_gwas(fetch = FALSE)
  
  # Test with very small top_n
  plot <- plot_top_gwas_hits(data, top_n = 2)
  expect_s3_class(plot, "ggplot")
  
  # Test with top_n larger than data
  plot <- plot_top_gwas_hits(data, top_n = 100)
  expect_s3_class(plot, "ggplot")
})

test_that("plot_gwas_summary creates ggplot object", {
  data <- get_celiac_gwas(fetch = FALSE)
  
  plot <- plot_gwas_summary(data)
  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot, "gg")
})

test_that("GWAS data has reasonable p-values", {
  data <- get_celiac_gwas(fetch = FALSE)
  
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
  plot1 <- plot_top_gwas_hits(mock_gwas_data)
  plot2 <- plot_gwas_summary(mock_gwas_data)
  
  expect_s3_class(plot1, "ggplot")
  expect_s3_class(plot2, "ggplot")
})