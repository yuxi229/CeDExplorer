test_that("scoreInflammationSignature returns expected structure", {
  data("example_expression_matrix", package = "CeDExplorer")
  data("example_metadata", package = "CeDExplorer")

  res <- scoreInflammationSignature(
    dataset = example_expression_matrix,
    signature = "IFN_gamma",
    plot = TRUE
  )

  expect_type(res, "list")
  expect_true(all(c("scores", "plot") %in% names(res)))
  expect_true(is.numeric(res$scores))
  expect_s3_class(res$plot, "ggplot")
})
