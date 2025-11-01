#' Plot gene interaction network using curated PPI data
#'
#' @description
#' Visualizes protein-protein interaction networks for celiac-related genes
#' using curated interaction data included in the package.
#'
#' @param genes Character vector of gene symbols to include in the network.
#' @param min_score Numeric; minimum interaction confidence score (0-1, default: 0.4).
#' @param layout Character; layout algorithm: "fr" (default), "circle", "kk", "dh".
#' @param color_by Character; node color by: "type" (default) or "importance".
#' @param show_labels Logical; whether to show gene labels (default: TRUE).
#'
#' @return A ggplot2 object showing the gene interaction network.
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic network with HLA genes
#' genes <- c("HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1")
#' plot_gene_network(genes)
#'
#' # Larger network with cytokine signaling
#' genes <- c("HLA-DQA1", "IL2", "IL2RA", "STAT1", "IFNG", "TNF")
#' plot_gene_network(genes, min_score = 0.5, layout = "circle")
#' }
plot_gene_network <- function(genes, min_score = 0.4, layout = "fr",
                              color_by = "type", show_labels = TRUE) {

  # Check required packages
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Please install igraph: install.packages('igraph')")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2: install.packages('ggplot2')")
  }

  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("Please install ggraph: install.packages('ggraph')")
  }

  # Validate inputs
  if (missing(genes) || length(genes) == 0) {
    stop("Please provide a vector of gene symbols.")
  }

  valid_layouts <- c("fr", "circle", "kk", "dh")
  if (!layout %in% valid_layouts) {
    stop("layout must be one of: ", paste(valid_layouts, collapse = ", "))
  }

  valid_color_by <- c("type", "importance")
  if (!color_by %in% valid_color_by) {
    stop("color_by must be one of: ", paste(valid_color_by, collapse = ", "))
  }

  # Load curated PPI data
  data("ppi_data", package = "CeDExplorer", envir = environment())

  # Filter interactions for our genes
  network_edges <- ppi_data$edges
  network_nodes <- ppi_data$nodes

  # Get interactions between our genes of interest
  selected_edges <- network_edges[
    network_edges$gene1 %in% genes &
      network_edges$gene2 %in% genes &
      network_edges$score >= min_score,
  ]

  # Get node information for our genes
  selected_nodes <- network_nodes[network_nodes$gene %in% genes, ]

  if (nrow(selected_edges) == 0) {
    warning("No interactions found between the provided genes with score >= ", min_score)

    # Create graph with isolated nodes only - SIMPLIFIED APPROACH
    g <- igraph::make_empty_graph(directed = FALSE)
    g <- igraph::add_vertices(g, nv = nrow(selected_nodes),
                              name = selected_nodes$gene,
                              type = selected_nodes$type,
                              importance = selected_nodes$importance)
  } else {
    # Create graph from edges - SIMPLIFIED APPROACH
    # First create all nodes, then add edges
    all_genes <- unique(c(selected_edges$gene1, selected_edges$gene2, selected_nodes$gene))
    node_attrs <- network_nodes[network_nodes$gene %in% all_genes, ]

    g <- igraph::graph_from_data_frame(
      d = selected_edges[, c("gene1", "gene2")],
      vertices = node_attrs,
      directed = FALSE
    )

    # Add edge weights
    igraph::E(g)$weight <- selected_edges$score
  }

  # Add degree centrality
  igraph::V(g)$degree <- igraph::degree(g)
  igraph::V(g)$label <- igraph::V(g)$name

  # Choose layout function
  layout_func <- switch(layout,
                        "fr" = igraph::layout_with_fr,
                        "circle" = igraph::layout_in_circle,
                        "kk" = igraph::layout_with_kk,
                        "dh" = igraph::layout_with_dh
  )

  # Create the plot with explicit layout calculation
  tryCatch({
    # Calculate layout explicitly
    coords <- layout_func(g)

    p <- ggraph::ggraph(g, layout = "manual", x = coords[,1], y = coords[,2])

    # Add edges if they exist
    if (igraph::ecount(g) > 0) {
      p <- p + ggraph::geom_edge_link(
        ggplot2::aes(alpha = weight),
        color = "gray50",
        show.legend = FALSE
      )
    }

    # Add nodes with coloring
    if (color_by == "type") {
      p <- p + ggraph::geom_node_point(
        ggplot2::aes(size = degree, color = type),
        alpha = 0.8
      )
    } else {
      p <- p + ggraph::geom_node_point(
        ggplot2::aes(size = degree, color = importance),
        alpha = 0.8
      ) + ggplot2::scale_color_viridis_c()
    }

    # Add labels if requested
    if (show_labels) {
      p <- p + ggraph::geom_node_text(
        ggplot2::aes(label = label),
        repel = TRUE,
        size = 3,
        max.overlaps = 20
      )
    }

    # Final styling
    p <- p +
      ggplot2::theme_void() +
      ggplot2::labs(
        title = "Gene Interaction Network",
        subtitle = paste("Nodes:", igraph::vcount(g), "| Edges:", igraph::ecount(g)),
        color = ifelse(color_by == "type", "Gene Type", "Importance")
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5),
        legend.position = "bottom"
      )

    return(p)

  }, error = function(e) {
    # Fallback: simple plot without ggraph
    warning("Advanced layout failed, creating simple network plot")

    if (igraph::vcount(g) == 0) {
      stop("No genes to plot")
    }

    # Simple base plot
    plot(g,
         vertex.label = igraph::V(g)$name,
         vertex.size = igraph::degree(g) * 3 + 5,
         vertex.color = "lightblue",
         main = paste("Gene Network:", paste(genes, collapse = ", ")))

    # Return NULL since we can't return a ggplot object in this case
    return(invisible(NULL))
  })
}
