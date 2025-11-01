# In test-plot_manhattan.R
test_that("plot_manhattan runs and returns a plot object", {
  # Test with ggplot2 backend (consistent return type)
  p <- plot_manhattan(use_qqman = FALSE)

  # Should return a ggplot object
  expect_s3_class(p, "ggplot")

  # Should have the expected title
  expect_match(p$labels$title, "Celiac disease GWAS")

  # Test that it works with custom data
  data <- get_celiac_gwas(fetch = FALSE)
  p_custom <- plot_manhattan(data = data, title = "Test Manhattan Plot", use_qqman = FALSE)
  expect_s3_class(p_custom, "ggplot")
  expect_match(p_custom$labels$title, "Test Manhattan Plot")
})
