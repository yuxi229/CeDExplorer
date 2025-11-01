#' Pathway enrichment analysis for a given gene list
#'
#' @description
#' Runs pathway enrichment (GO, KEGG, Reactome, etc.) using g:Profiler
#' for a user-supplied list of genes, returning an enrichment results table
#' and a simple ggplot2 barplot of top terms.
#'
#' @param genes Character vector of gene symbols (e.g., c("HLA-DQA1", "IL2", "TNF")).
#' @param source Character; database(s) to use for enrichment (default: "GO:BP", "KEGG", "REAC").
#' @param organism Character; organism code for g:Profiler (default: "hsapiens").
#' @param top Integer; number of top terms to display in the barplot (default: 10).
#'
#' @return A list with elements:
#' \describe{
#'   \item{results}{data.frame of enriched terms and p-values.}
#'   \item{plot}{ggplot object showing the top enriched pathways.}
#' }
#' @export
#' @examples
#' \dontrun{
#' data("celiac_gwas_example", package = "CeDExplorer")
#' genes <- unique(na.omit(celiac_gwas_example$mapped_gene))
#' res <- map_to_pathways(genes)
#' res$plot
#' }
map_to_pathways <- function(genes,
                            source = c("GO:BP", "KEGG", "REAC"),
                            organism = "hsapiens",
                            top = 10) {
  if (missing(genes) || length(genes) == 0) {
    stop("Please provide a vector of gene symbols.")
  }

  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop("Please install the 'gprofiler2' package.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install the 'ggplot2' package.")
  }

  # Run g:Profiler enrichment
  gp <- gprofiler2::gost(
    query = genes,
    organism = organism,
    sources = source,
    correction_method = "fdr"
  )

  if (is.null(gp$result) || nrow(gp$result) == 0) {
    warning("No enriched pathways found for the given genes.")
    return(list(results = data.frame(), plot = NULL))
  }

  # Prepare results
  res <- gp$result[, c("term_name", "source", "p_value", "term_size", "intersection_size")]
  res <- res[order(res$p_value), ]
  res$term_name <- factor(res$term_name, levels = rev(res$term_name))

  # Select top pathways
  top_res <- head(res, top)

  # Barplot using ggplot2
  p <- ggplot2::ggplot(top_res, ggplot2::aes(x = term_name, y = -log10(p_value), fill = source)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      x = NULL,
      y = expression(-log[10]("FDR-adjusted p-value")),
      title = "Top Enriched Pathways"
    )

  return(list(results = res, plot = p))
}
