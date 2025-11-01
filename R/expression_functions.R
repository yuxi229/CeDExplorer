#' Expression boxplot for a given gene
#'
#' @description
#' Creates a boxplot comparing expression levels between CeD and control samples.
#'
#' @param gene Character; gene symbol to plot (e.g., "HLA-DQA1").
#' @param dataset Data frame; optional dataset with expression data.
#' @param violin Logical; whether to draw violin plots instead of boxplots.
#' @return A ggplot2 object showing expression distribution between groups.
#' @export
#' @importFrom utils data
#' @importFrom stats wilcox.test
#' @importFrom grDevices colorRampPalette
#' @importFrom utils head
expression_boxplot <- function(gene, dataset = NULL, violin = FALSE) {
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
  if (length(unique(df$condition)) == 2) {
    pval <- suppressWarnings(stats::wilcox.test(expression ~ condition, data = df)$p.value)
    p <- p + ggplot2::annotate(
      "text", x = 1.5, y = max(df$expression, na.rm = TRUE),
      label = paste0("p = ", signif(pval, 3)), size = 4
    )
  }

  return(p)
}



#' Heatmap of selected genes across samples
#'
#' @description
#' Creates a clustered heatmap of selected genes across samples,
#' using ComplexHeatmap or pheatmap.
#'
#' @param genes Character vector; list of genes to plot.
#' @param dataset Expression matrix or data frame (rows = genes, columns = samples).
#'        If NULL, uses bundled example dataset.
#' @param scale_rows Logical; whether to scale expression values by gene (default TRUE).
#' @param cluster Logical; whether to cluster genes/samples (default TRUE).
#'
#' @return A heatmap object (ComplexHeatmap or pheatmap).
#' @export
#' @examples
#' \dontrun{
#' celiac_heatmap(c("HLA-DQA1", "IL2", "TNF"))
#' }
celiac_heatmap <- function(genes, dataset = NULL, scale_rows = TRUE, cluster = TRUE) {
  if (missing(genes) || length(genes) == 0) {
    stop("Please provide at least one gene name.")
  }

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Please install 'ComplexHeatmap'.")
  }

  if (is.null(dataset)) {
    data("example_expression_matrix", package = "CeDExplorer", envir = environment())
    dataset <- example_expression_matrix
  }

  # Ensure matrix format (genes as rownames)
  if (is.data.frame(dataset)) {
    dataset <- as.matrix(dataset)
  }

  rownames(dataset) <- gsub("\\s+", "", rownames(dataset))
  genes <- genes[genes %in% rownames(dataset)]

  if (length(genes) == 0) {
    warning("None of the provided genes found in dataset.")
    return(NULL)
  }

  mat <- dataset[genes, , drop = FALSE]

  # Optional scaling per gene
  if (scale_rows) {
    mat <- t(scale(t(mat)))
  }

  # Define annotation colors
  hm <- ComplexHeatmap::Heatmap(
    mat,
    name = "Expression (z-score)",
    col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    cluster_rows = cluster,
    cluster_columns = cluster,
    row_names_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 8),
    column_title = "Samples",
    row_title = "Genes"
  )

  ComplexHeatmap::draw(hm)
  invisible(hm)
}

