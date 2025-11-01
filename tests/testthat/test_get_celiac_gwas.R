test_that("get_celiac_gwas returns data.frame with expected columns", {
  df <- get_celiac_gwas(fetch = FALSE)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("SNP", "mapped_gene", "chromosome", "position", "p_value") %in% names(df)))
  expect_true(nrow(df) > 0)
})
