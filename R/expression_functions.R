# Purpose: Create boxplot or violin plot comparing gene expression between CeD and control samples
# Author: Yuxi Zhang
# Date: Nov 27, 2025
# Version: 1.0
# Bugs and Issues: None known

#' Expression boxplot for a given gene
#'
#' @name expressionBoxplot
#' 
#' @description
#' Creates a boxplot comparing expression levels between CeD and control samples.
#'
#' @param gene Character; gene symbol to plot (e.g., "HLA-DQA1").
#' @param dataset Data frame; optional dataset with expression data.
#' @param violin Logical; whether to draw violin plots instead of boxplots.
#' 
#' @references
#' Hintze, J. L., & Nelson, R. D. (1998). Violin plots: A box plot-density 
#' trace synergism. The American Statistician, 52(2), 181-184. 
#' https://doi.org/10.1080/00031305.1998.10480559
#'
#' Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold 
#' change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550. 
#' https://doi.org/10.1186/s13059-014-0550-8
#'
#' McCarthy, D. J., Chen, Y., & Smyth, G. K. (2012). Differential expression 
#' analysis of multifactor RNA-Seq experiments with respect to biological 
#' variation. Nucleic Acids Research, 40(10), 4288-4297. 
#' https://doi.org/10.1093/nar/gks042
#'
#' Wilcoxon, F. (1945). Individual comparisons by ranking methods. 
#' Biometrics Bulletin, 1(6), 80-83. 
#' https://doi.org/10.2307/3001968
#' 
#' @return A ggplot2 object showing expression distribution between groups.
#' @export
#' @importFrom utils data
#' @importFrom stats wilcox.test
#' @importFrom grDevices colorRampPalette
#' @importFrom utils head
#' @examples
#' \donttest{
#' # Basic boxplot for HLA-DQA1
#' expressionBoxplot("HLA-DQA1")
#'
#' # Violin plot for HLA-DQB1
#' expressionBoxplot("HLA-DQB1", violin = TRUE)
#'
#' # Using custom dataset
#' custom_data <- data.frame(
#'   gene = rep(c("HLA-DQA1", "HLA-DQB1"), each = 20),
#'   sample = paste0("S", 1:40),
#'   condition = rep(rep(c("CeD", "Control"), each = 10), 2),
#'   expression = c(rnorm(10, 8, 1), rnorm(10, 5, 1), 
#'                  rnorm(10, 7, 1), rnorm(10, 4, 1))
#' )
#' expressionBoxplot("HLA-DQA1", dataset = custom_data)
#' }
expressionBoxplot <- function(gene, dataset = NULL, violin = FALSE) {
  # ---- Argument checks ----
  if (missing(gene) || is.null(gene) || !is.character(gene) || length(gene) == 0) {
    stop("Please provide a gene symbol.")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install 'ggplot2'.")
  }

  # ---- Load bundled dataset if none provided ----
  if (is.null(dataset)) {
    data("example_expression", package = "CeDExplorer", envir = environment())
    dataset <- example_expression
  }

  # ---- Validate structure ----
  required_cols <- c("gene", "sample", "condition", "expression")
  if (!all(required_cols %in% colnames(dataset))) {
    stop("Dataset must contain columns: gene, sample, condition, expression.")
  }

  # ---- Filter for gene ----
  # FIX: Use base R filtering instead of !! operator
  df <- dataset[dataset$gene == gene, , drop = FALSE]

  if (nrow(df) == 0) {
    warning(paste("Gene", gene, "not found in dataset. Available genes:",
                  paste(utils::head(unique(dataset$gene), 5), collapse = ", "), "..."))
    return(NULL)
  }

  # ---- Plot ----
  p <- ggplot2::ggplot(df, ggplot2::aes(x = condition, y = expression, fill = condition))
  
  # Use violin plot if requested (Hintze & Nelson, 1998)
  if (isTRUE(violin)) {
    p <- p + ggplot2::geom_violin(trim = FALSE, alpha = 0.5)
  } else {
    p <- p + ggplot2::geom_boxplot(alpha = 0.6)
  }

  p <- p +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = paste("Expression of", gene, "in CeD vs Control"),
      x = NULL, y = "Expression level"
    ) +
    ggplot2::theme(legend.position = "none")

  # ---- Optional Wilcoxon test annotation ----
  # Wilcoxon rank-sum test (Wilcoxon, 1945)
  if (length(unique(df$condition)) == 2) {
    pval <- suppressWarnings(stats::wilcox.test(expression ~ condition, data = df)$p.value)
    p <- p + ggplot2::annotate(
      "text", x = 1.5, y = max(df$expression, na.rm = TRUE),
      label = paste0("p = ", signif(pval, 3)), size = 4
    )
  }

  return(p)
}

# [END] 


