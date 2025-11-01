#' Composite evidence view for a gene
#'
#' @description
#' Creates a comprehensive dashboard showing multiple evidence types for a gene
#' including expression, GWAS associations, network interactions, and HLA associations.
#'
#' @param gene Character; gene symbol to analyze (e.g., "HLA-DQA1").
#' @param show_expression Logical; whether to show expression boxplot (default: TRUE).
#' @param show_gwas Logical; whether to show GWAS Manhattan plot highlight (default: TRUE).
#' @param show_network Logical; whether to show network interactions (default: TRUE).
#' @param show_hla Logical; whether to show HLA associations (default: TRUE).
#' @param ncol Number of columns in the layout (default: 2).
#'
#' @return A composite ggplot2 object or list of plots.
#' @export
#'
#' @examples
#' \dontrun{
#' # Comprehensive view for HLA-DQA1
#' geneEvidenceView("HLA-DQA1")
#'
#' # Focus on specific evidence types
#' geneEvidenceView("HLA-DQB1", show_network = FALSE, show_hla = FALSE)
#' }
geneEvidenceView <- function(gene,
                             show_expression = TRUE,
                             show_gwas = TRUE,
                             show_network = TRUE,
                             show_hla = TRUE,
                             ncol = 2) {

  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2: install.packages('ggplot2')")
  }

  # Validate inputs
  if (missing(gene) || !is.character(gene) || length(gene) == 0) {
    stop("Please provide a gene symbol.")
  }

  # Initialize plot list
  plot_list <- list()

  # 1. Expression Boxplot
  if (show_expression) {
    if (requireNamespace("CeDExplorer", quietly = TRUE)) {
      tryCatch({
        expr_plot <- expression_boxplot(gene)
        if (!is.null(expr_plot)) {
          expr_plot <- expr_plot +
            ggplot2::labs(title = paste("Expression:", gene)) +
            ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
          plot_list[["expression"]] <- expr_plot
        }
      }, error = function(e) {
        message("Could not create expression plot: ", e$message)
      })
    }
  }

  # 2. GWAS Manhattan Plot Highlight
  if (show_gwas) {
    tryCatch({
      # Load GWAS data and create a focused plot around the gene
      gwas_plot <- create_gwas_gene_highlight(gene)
      if (!is.null(gwas_plot)) {
        plot_list[["gwas"]] <- gwas_plot
      }
    }, error = function(e) {
      message("Could not create GWAS plot: ", e$message)
    })
  }

  # 3. Network Interactions
  if (show_network) {
    if (requireNamespace("igraph", quietly = TRUE)) {
      tryCatch({
        # Get genes that interact with our target gene
        network_genes <- get_interacting_genes(gene)
        if (length(network_genes) > 0) {
          network_plot <- plot_gene_network(network_genes) +
            ggplot2::labs(title = paste("Network:", gene)) +
            ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
          plot_list[["network"]] <- network_plot
        }
      }, error = function(e) {
        message("Could not create network plot: ", e$message)
      })
    }
  }

  # 4. HLA Associations
  if (show_hla) {
    tryCatch({
      hla_plot <- create_hla_association_plot(gene)
      if (!is.null(hla_plot)) {
        plot_list[["hla"]] <- hla_plot
      }
    }, error = function(e) {
      message("Could not create HLA plot: ", e$message)
    })
  }

  # Check if we have any plots
  if (length(plot_list) == 0) {
    stop("No plots could be generated for gene: ", gene)
  }

  # Combine plots
  if (requireNamespace("patchwork", quietly = TRUE)) {
    # Use patchwork for elegant combination
    combined_plot <- patchwork::wrap_plots(plot_list, ncol = ncol) +
      patchwork::plot_annotation(
        title = paste("Evidence Dashboard:", gene),
        theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14))
      )
    return(combined_plot)
  } else {
    # Return list if patchwork not available
    message("Install 'patchwork' for combined layout. Returning list of plots.")
    return(plot_list)
  }
}

#' Create GWAS gene highlight plot
#' @param gene Gene symbol
#' @return ggplot2 object or NULL
#' @keywords internal
create_gwas_gene_highlight <- function(gene) {
  tryCatch({
    # Load GWAS data
    data("celiac_gwas_example", package = "CeDExplorer", envir = environment())

    # Filter for genes near our target (simple approach)
    gwas_data <- celiac_gwas_example
    target_genes <- unique(gwas_data$mapped_gene)

    # Check if gene is in GWAS data
    if (gene %in% target_genes) {
      # Create a focused Manhattan plot around this gene
      gene_data <- gwas_data[gwas_data$mapped_gene == gene, ]

      if (nrow(gene_data) > 0) {
        p <- ggplot2::ggplot(gene_data, ggplot2::aes(x = as.numeric(position), y = -log10(p_value))) +
          ggplot2::geom_point(color = "red", size = 2, alpha = 0.7) +
          ggplot2::geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "blue") +
          ggplot2::theme_minimal() +
          ggplot2::labs(
            title = paste("GWAS:", gene),
            x = "Genomic Position",
            y = "-log10(p-value)"
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = 10),
            axis.text = ggplot2::element_text(size = 8)
          )
        return(p)
      }
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}

#' Get genes that interact with target gene
#' @param gene Gene symbol
#' @return Vector of interacting genes
#' @keywords internal
get_interacting_genes <- function(gene) {
  tryCatch({
    data("ppi_data", package = "CeDExplorer", envir = environment())

    # Find genes that interact with our target
    interactions <- ppi_data$edges
    interacting_genes <- unique(c(
      interactions$gene1[interactions$gene2 == gene],
      interactions$gene2[interactions$gene1 == gene]
    ))

    # Include the target gene itself
    all_genes <- unique(c(gene, interacting_genes))

    # Limit to reasonable number for plotting
    if (length(all_genes) > 10) {
      all_genes <- c(gene, head(interacting_genes, 9))
    }

    return(all_genes)
  }, error = function(e) {
    return(character(0))
  })
}

#' Create HLA association plot
#' @param gene Gene symbol
#' @return ggplot2 object or NULL
#' @keywords internal
create_hla_association_plot <- function(gene) {
  tryCatch({
    # Check if it's an HLA gene
    if (grepl("^HLA-", gene)) {
      # Create a simple HLA association plot
      hla_data <- data.frame(
        condition = c("CeD", "Control"),
        frequency = c(0.15, 0.05),  # Example frequencies
        expression = c(8, 5)        # Example expression
      )

      p <- ggplot2::ggplot(hla_data, ggplot2::aes(x = frequency, y = expression, color = condition)) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_segment(ggplot2::aes(xend = frequency, yend = expression), size = 1) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste("HLA:", gene),
          x = "Allele Frequency",
          y = "Expression",
          color = "Condition"
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 10),
          legend.position = "bottom"
        )
      return(p)
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}

#' Create summary statistics for gene
#' @param gene Gene symbol
#' @return Data.frame with summary stats
#' @keywords internal
get_gene_summary <- function(gene) {
  # This would typically query various databases
  # For now, return example data
  data.frame(
    metric = c("GWAS Associations", "Network Partners", "Expression Fold-Change", "HLA Risk"),
    value = c("5 SNPs", "8 genes", "1.6x", "High"),
    stringsAsFactors = FALSE
  )
}
