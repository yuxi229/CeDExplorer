# Purpose: Visualize relationship between HLA allele frequencies and gene expression
# Author: Yuxi Zhang
# Date: Nov 27, 2025
# Version: 1.0
# Bugs and Issues: Uses simulated data for demonstration
#'
#' Plot HLA allele frequency vs gene expression
#'
#' @name hlaInteractionPlot
#' 
#' @description
#' Visualizes the relationship between HLA allele frequencies and gene expression levels
#' across different conditions (e.g., CeD vs Control). Uses simulated data for demonstration.
#'
#' @param gene Character; gene symbol to plot (e.g., "HLA-DQA1").
#' @param dataset Expression data.frame with columns: sample, condition, and gene columns.
#'        If NULL, uses example data.
#' @param allele_freq Data.frame with HLA allele frequencies. If NULL, uses example data.
#' @param condition_col Character; column name in dataset indicating condition (default: "condition").
#' @param point_size Numeric; size of points in scatter plot (default: 2).
#' @param add_trend Logical; whether to add trend lines (default: TRUE).
#' 
#' @references
#' Sollid, L. M., & Jabri, B. (2013). Triggers and drivers of autoimmunity: 
#' Lessons from coeliac disease. Nature Reviews Immunology, 13(4), 294-302.
#' https://doi.org/10.1038/nri3407
#'
#' Trynka, G., Hunt, K. A., Bockett, N. A., Romanos, J., Mistry, V., 
#' Szperl, A., Bakker, S. F., Bardella, M. T., Bhaw-Rosun, L., Castillejo, G., 
#' de la Concha, E. G., de Almeida, R. C., Dias, K. R., van Diemen, C. C., 
#' Dubois, P. C., Duerr, R. H., Edkins, S., Franke, L., Fransen, K., ... 
#' van Heel, D. A. (2011). Dense genotyping identifies and localizes multiple 
#' common and rare variant association signals in celiac disease. 
#' Nature Genetics, 43(12), 1193-1201. https://doi.org/10.1038/ng.998
#' 
#' @return A ggplot2 object showing HLA allele frequency vs gene expression.
#' @export
#' @importFrom utils data
#' @importFrom stats cor.test
#' @importFrom stats rnorm
#' @examples
#' \dontrun{
#' # Basic plot with HLA-DQA1
#' hlaInteractionPlot("HLA-DQA1")
#'
#' # Plot with HLA-DQB1 and custom styling
#' hlaInteractionPlot("HLA-DQB1", point_size = 3, add_trend = TRUE)
#' }
hlaInteractionPlot <- function(gene,
                               dataset = NULL,
                               allele_freq = NULL,
                               condition_col = "condition", 
                               point_size = 2,
                               add_trend = TRUE) {
  
  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2: install.packages('ggplot2')")
  }
  
  # Validate inputs
  if (missing(gene) || !is.character(gene) || length(gene) == 0) {
    stop("Please provide a gene symbol.")
  }
  
  # Load example data if not provided
  if (is.null(dataset)) {
    dataset <- getHlaExpressionData()
  }
  
  if (is.null(allele_freq)) {
    allele_freq <- getHlaAlleleFrequencies()
  }
  
  # Validate dataset structure
  required_cols <- c("sample", condition_col, gene)
  missing_cols <- setdiff(required_cols, colnames(dataset))
  if (length(missing_cols) > 0) {
    stop("Dataset missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Validate allele frequency data
  if (!all(c("allele", "frequency", "condition") %in% colnames(allele_freq))) {
    stop("allele_freq must contain columns: allele, frequency, condition")
  }
  
  # Check if gene exists in allele frequency data
  if (!gene %in% allele_freq$allele) {
    stop("Gene '", gene, "' not found in allele frequency data. Available: ",
         paste(unique(allele_freq$allele), collapse = ", "))
  }
  
  # Filter allele frequency data for the specific gene
  gene_freq <- allele_freq[allele_freq$allele == gene, ]
  
  # Merge expression and frequency data
  plot_data <- merge(dataset[, c("sample", condition_col, gene)],
                     gene_freq[, c("condition", "frequency")],
                     by.x = condition_col,
                     by.y = "condition",
                     all.x = TRUE)
  
  # Rename columns for clarity
  colnames(plot_data)[colnames(plot_data) == gene] <- "expression"
  colnames(plot_data)[colnames(plot_data) == condition_col] <- "condition"
  
  # Remove any NA values
  plot_data <- plot_data[!is.na(plot_data$expression) & !is.na(plot_data$frequency), ]
  
  if (nrow(plot_data) == 0) {
    stop("No data available for plotting after merging expression and frequency data.")
  }
  
  # Create the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = frequency, y = expression, color = condition)) +
    ggplot2::geom_point(size = point_size, alpha = 0.7) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = paste("HLA Allele Frequency vs", gene, "Expression"),
      x = paste(gene, "Allele Frequency"),
      y = paste(gene, "Expression Level"),
      color = "Condition"
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  # Add trend lines if requested
  if (add_trend && length(unique(plot_data$condition)) > 1) {
    p <- p + ggplot2::geom_smooth(ggplot2::aes(group = condition), method = "lm", se = TRUE, alpha = 0.2)
  } else if (add_trend) {
    p <- p + ggplot2::geom_smooth(method = "lm", se = TRUE, alpha = 0.2, color = "blue")
  }
  
  # Add correlation annotation
  corr_text <- calculateCorrelation(plot_data$frequency, plot_data$expression)
  p <- p + ggplot2::annotate("text",
                             x = max(plot_data$frequency, na.rm = TRUE),
                             y = max(plot_data$expression, na.rm = TRUE),
                             label = corr_text,
                             hjust = 1, vjust = 1,
                             size = 4, color = "darkred")
  
  return(p)
}

#' Calculate correlation text for annotation
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @return Formatted correlation string
#' @keywords internal
#' @importFrom stats cor.test
calculateCorrelation <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) {
    return("Insufficient data")
  }
  
  cor_test <- cor.test(x, y, method = "pearson")
  r <- round(cor_test$estimate, 3)
  p <- round(cor_test$p.value, 4)
  
  if (p < 0.001) {
    p_text <- "p < 0.001"
  } else {
    p_text <- paste("p =", p)
  }
  
  return(paste0("r = ", r, "\n", p_text))
}

#' Get example HLA expression data
#'
#' @return Data.frame with HLA expression data
#' @keywords internal
#' @importFrom stats rnorm
getHlaExpressionData <- function() {
  set.seed(123)
  samples <- paste0("S", 1:100)
  conditions <- rep(c("CeD", "Control"), each = 50)
  
  # Simulate HLA expression data with some biological patterns
  data.frame(
    sample = samples,
    condition = conditions,
    `HLA-DQA1` = ifelse(conditions == "CeD",
                        rnorm(50, mean = 8, sd = 1.5),
                        rnorm(50, mean = 5, sd = 1.2)),
    `HLA-DQB1` = ifelse(conditions == "CeD",
                        rnorm(50, mean = 7.5, sd = 1.3),
                        rnorm(50, mean = 4.8, sd = 1.1)),
    `HLA-DRA` = ifelse(conditions == "CeD",
                       rnorm(50, mean = 6.8, sd = 1.1),
                       rnorm(50, mean = 6.5, sd = 1.0)),
    `HLA-DRB1` = ifelse(conditions == "CeD",
                        rnorm(50, mean = 7.2, sd = 1.4),
                        rnorm(50, mean = 5.0, sd = 1.3)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Get example HLA allele frequency data
#'
#' @return Data.frame with HLA allele frequencies
#' @keywords internal
getHlaAlleleFrequencies <- function() {
  # Simulate allele frequency data based on known CeD associations
  data.frame(
    allele = rep(c("HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1"), each = 2),
    condition = rep(c("CeD", "Control"), 4),
    frequency = c(
      # HLA-DQA1: Higher in CeD
      0.15, 0.05,
      # HLA-DQB1: Higher in CeD
      0.12, 0.04,
      # HLA-DRA: Similar frequencies
      0.08, 0.07,
      # HLA-DRB1: Higher in CeD
      0.18, 0.06
    ),
    stringsAsFactors = FALSE
  )
}

# [END] 
