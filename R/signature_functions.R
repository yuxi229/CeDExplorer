#' Compute inflammation signature scores for CeD samples
#'
#' @description
#' Calculates per-sample signature scores from expression data for a given
#' inflammation-related gene set (e.g., IFN-gamma or TNF-alpha response).
#' Supports simple mean expression or GSVA scoring.
#'
#' @param dataset Matrix or data.frame of gene expression
#'        (rows = genes, columns = samples). If NULL, uses example data.
#' @param signature Character string; name of the gene signature
#'        (e.g., "IFN_gamma", "TNF_alpha").
#' @param method Character; scoring method ("mean" or "gsva"). Default "mean".
#' @param plot Logical; whether to return a ggplot2 boxplot comparing CeD vs Control.
#'
#' @return A list with elements:
#' \describe{
#'   \item{scores}{Named numeric vector of per-sample signature scores.}
#'   \item{plot}{ggplot2 object (if plot = TRUE).}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' res <- score_inflammation_signature(signature = "IFN_gamma")
#' res$plot
#' }
score_inflammation_signature <- function(dataset = NULL,
                                         signature = "IFN_gamma",
                                         method = "mean",
                                         plot = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install 'ggplot2'.")
  }

  # Optionally use GSVA for more advanced scoring
  if (tolower(method) == "gsva" && !requireNamespace("GSVA", quietly = TRUE)) {
    stop("Please install 'GSVA' package for GSVA-based scoring.")
  }

  # Load example data if none provided
  if (is.null(dataset)) {
    data("example_expression_matrix", package = "CeDExplorer", envir = environment())
    dataset <- example_expression_matrix
  }

  # Load example metadata (CeD vs Control)
  if (!exists("example_metadata", envir = environment())) {
    data("example_metadata", package = "CeDExplorer", envir = environment())
  }
  meta <- get("example_metadata", envir = environment())

  # Example built-in gene signatures
  sig_db <- list(
    IFN_gamma = c("IFNG", "STAT1", "CXCL9", "CXCL10", "IRF1"),
    TNF_alpha = c("TNF", "NFKB1", "IL1B", "RELA", "CCL2")
  )

  if (!signature %in% names(sig_db)) {
    stop(paste("Signature not recognized. Choose from:",
               paste(names(sig_db), collapse = ", ")))
  }

  genes <- sig_db[[signature]]
  genes_present <- genes[genes %in% rownames(dataset)]

  if (length(genes_present) == 0) {
    warning("None of the signature genes found in dataset.")
    return(list(scores = numeric(0), plot = NULL))
  }

  # Compute signature score per sample
  if (tolower(method) == "mean") {
    scores <- colMeans(dataset[genes_present, , drop = FALSE], na.rm = TRUE)
  } else if (tolower(method) == "gsva") {
    gsva_res <- GSVA::gsva(expr = as.matrix(dataset),
                           gset.idx.list = list(signature = genes_present),
                           method = "ssgsea",
                           verbose = FALSE)
    scores <- gsva_res[1, ]
  } else {
    stop("Unknown scoring method. Use 'mean' or 'gsva'.")
  }

  out <- list(scores = scores)

  # Optional boxplot comparing CeD vs Control samples
  if (plot && !is.null(meta)) {
    df <- data.frame(
      sample = names(scores),
      score = scores,
      condition = meta$condition[match(names(scores), meta$sample)]
    )
    df <- subset(df, !is.na(condition))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = condition, y = score, fill = condition)) +
      ggplot2::geom_boxplot(alpha = 0.6) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = paste("Inflammation Signature Score:", signature),
        x = NULL, y = "Signature score"
      ) +
      ggplot2::theme(legend.position = "none")

    # Add Wilcoxon p-value if two groups exist
    if (length(unique(df$condition)) == 2) {
      pval <- suppressWarnings(wilcox.test(score ~ condition, data = df)$p.value)
      p <- p + ggplot2::annotate("text", x = 1.5, y = max(df$score, na.rm = TRUE),
                                 label = paste0("p = ", signif(pval, 3)), size = 4)
    }

    out$plot <- p
  }

  return(out)
}

