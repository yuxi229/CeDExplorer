# tests/testthat/test-plot_functions.R

test_that("plot_manhattan function works correctly", {
  # Skip if required packages not available
  skip_if_not_installed("ggplot2")

  # Load example data
  data("celiac_gwas_example", package = "CeDExplorer")

  # Test 1: Basic functionality with default data
  test_that("plot_manhattan works with default data", {
    p <- plot_manhattan()
    expect_s3_class(p, "ggplot")
    expect_true("position" %in% names(p$data))
    expect_true("logp" %in% names(p$data))
  })

  # Test 2: Works with custom data
  test_that("plot_manhattan works with custom data", {
    custom_data <- celiac_gwas_example
    p <- plot_manhattan(data = custom_data, title = "Test Plot")
    expect_s3_class(p, "ggplot")
    expect_match(p$labels$title, "Test Plot")
  })

  # Test 3: Handles highlighting
  test_that("plot_manhattan handles highlighting", {
    # Test with a few SNPs to highlight
    highlight_snps <- celiac_gwas_example$SNP[1:3]
    p <- plot_manhattan(data = celiac_gwas_example, highlight = highlight_snps)
    expect_s3_class(p, "ggplot")
  })

  # Test 4: Custom genome-wide line
  test_that("plot_manhattan accepts custom genome-wide line", {
    p <- plot_manhattan(data = celiac_gwas_example, genomewide_line = 1e-6)
    expect_s3_class(p, "ggplot")
    # Check that the line is plotted (should be at -log10(1e-6) = 6)
    expect_true(6 %in% p$layers[[2]]$data)  # geom_hline yintercept
  })

  # Test 5: Error handling for missing columns
  test_that("plot_manhattan errors on missing columns", {
    bad_data <- celiac_gwas_example[, c("SNP", "mapped_gene")]  # Missing required cols
    expect_error(
      plot_manhattan(data = bad_data),
      "Input data must contain columns: chromosome, position, p_value"
    )
  })

  # Test 6: Works with minimal valid data
  test_that("plot_manhattan works with minimal valid data", {
    minimal_data <- data.frame(
      chromosome = c("1", "1", "2"),
      position = c(1000, 2000, 1500),
      p_value = c(0.001, 0.0001, 0.01)
    )
    p <- plot_manhattan(data = minimal_data)
    expect_s3_class(p, "ggplot")
  })
})

test_that("expression_boxplot function works correctly", {
  skip_if_not_installed("ggplot2")

  # Load example data
  data("example_expression", package = "CeDExplorer")

  test_that("expression_boxplot works with valid gene", {
    p <- expression_boxplot("HLA-DQA1")
    expect_s3_class(p, "ggplot")
    expect_match(p$labels$title, "HLA-DQA1")
    expect_true("condition" %in% names(p$data))
    expect_true("expression" %in% names(p$data))
  })

  test_that("expression_boxplot creates violin plots", {
    p <- expression_boxplot("HLA-DQA1", violin = TRUE)
    expect_s3_class(p, "ggplot")
  })

  test_that("expression_boxplot works with custom dataset", {
    custom_data <- example_expression
    p <- expression_boxplot("HLA-DQB1", dataset = custom_data)
    expect_s3_class(p, "ggplot")
    expect_match(p$labels$title, "HLA-DQB1")
  })

  test_that("expression_boxplot handles non-existent genes", {
    expect_warning(
      p <- expression_boxplot("NONEXISTENT_GENE"),
      "not found in dataset"
    )
    expect_null(p)
  })

  test_that("expression_boxplot errors on invalid dataset", {
    bad_data <- example_expression[, c("gene", "sample")]  # Missing condition, expression
    expect_error(
      expression_boxplot("HLA-DQA1", dataset = bad_data),
      "Dataset must contain columns: gene, sample, condition, expression"
    )
  })

  test_that("expression_boxplot filters correct gene", {
    p <- expression_boxplot("IL2RA")
    expect_s3_class(p, "ggplot")
    # Check that only IL2RA data is in the plot
    expect_true(all(p$data$gene == "IL2RA"))
  })
})

# Integration tests
test_that("plot functions integrate well with data functions", {
  skip_if_not_installed("ggplot2")

  # Test 1: plot_manhattan with fetched data
  test_that("plot_manhattan works with get_celiac_gwas data", {
    gwas_data <- get_celiac_gwas(fetch = FALSE)  # Use bundled data
    p <- plot_manhattan(data = gwas_data)
    expect_s3_class(p, "ggplot")
    expect_true(nrow(p$data) > 0)
  })

  # Test 2: expression_boxplot with all available genes
  test_that("expression_boxplot works with all genes in example data", {
    data("example_expression", package = "CeDExplorer")
    genes <- unique(example_expression$gene)

    for (gene in genes) {
      p <- expression_boxplot(gene)
      expect_s3_class(p, "ggplot")
      expect_match(p$labels$title, gene)
    }
  })
})

# Edge case tests
test_that("plot functions handle edge cases", {
  skip_if_not_installed("ggplot2")

  # Test 1: plot_manhattan with very small dataset
  test_that("plot_manhattan works with single chromosome", {
    single_chr_data <- celiac_gwas_example[celiac_gwas_example$chromosome == "6", ]
    p <- plot_manhattan(data = single_chr_data)
    expect_s3_class(p, "ggplot")
  })

  # Test 2: expression_boxplot with NA values
  test_that("expression_boxplot handles NA values gracefully", {
    data_with_na <- example_expression
    data_with_na$expression[1:5] <- NA
    p <- expression_boxplot("HLA-DQA1", dataset = data_with_na)
    expect_s3_class(p, "ggplot")
  })

  # Test 3: plot_manhattan with custom colors (if implemented)
  test_that("plot_manhattan accepts visual customization", {
    p <- plot_manhattan(data = celiac_gwas_example, title = "Custom Title")
    expect_s3_class(p, "ggplot")
    expect_match(p$labels$title, "Custom Title")
  })
})
