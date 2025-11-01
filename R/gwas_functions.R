#' Get curated celiac disease GWAS summary table
#'
#' @description
#' Loads a bundled curated CeD GWAS dataset, or optionally fetches the latest
#' GWAS associations from the GWAS Catalog via the gwasrapidd API.
#'
#' @param fetch Logical; if TRUE, retrieves celiac disease data from the GWAS Catalog.
#' @param save Logical; if TRUE, saves fetched results to "celiac_gwas_fetched.tsv".
#'
#' @return A data.frame with SNP, mapped_gene, chromosome, position, p_value, odds_ratio, source.
#' @export
get_celiac_gwas <- function(fetch = FALSE, save = FALSE) {
  if (!fetch) {
    data("celiac_gwas_example", package = "CeDExplorer", envir = environment())
    return(celiac_gwas_example)
  }

  if (!requireNamespace("gwasrapidd", quietly = TRUE)) {
    stop("Please install the 'gwasrapidd' package to fetch live data.")
  }

  message("Fetching celiac disease GWAS data from GWAS Catalog...")
  studies <- gwasrapidd::get_studies(efo_trait = "celiac disease")
  study_ids <- studies@studies$study_id

  if (length(study_ids) == 0) {
    stop("No celiac disease studies found in GWAS Catalog.")
  }

  assoc_list <- lapply(study_ids, function(id) {
    x <- gwasrapidd::get_associations(study_id = id)
    df <- x@associations

    # Skip if no associations
    if (nrow(df) == 0) return(NULL)

    # Safely extract columns (use NA if missing)
    data.frame(
      SNP = if ("variant_id" %in% colnames(df)) df$variant_id else NA,
      mapped_gene = if ("mapped_gene" %in% colnames(df)) df$mapped_gene else NA,
      chromosome = if ("chromosome_name" %in% colnames(df)) df$chromosome_name else NA,
      position = if ("chromosome_position" %in% colnames(df)) df$chromosome_position else NA,
      p_value = if ("pvalue" %in% colnames(df)) df$pvalue else NA,
      odds_ratio = if ("or" %in% colnames(df)) df$or else NA,
      source = "GWAS Catalog",
      stringsAsFactors = FALSE
    )
  })

  # Remove NULL elements (empty studies)
  assoc_list <- assoc_list[!vapply(assoc_list, is.null, logical(1))]

  if (length(assoc_list) == 0) {
    stop("No association data available for celiac disease.")
  }

  all_assocs <- do.call(rbind, assoc_list)

  # Clean up numeric columns
  all_assocs$p_value <- as.numeric(all_assocs$p_value)
  all_assocs$position <- as.numeric(all_assocs$position)
  all_assocs$chromosome <- as.character(all_assocs$chromosome)

  # Optional save
  if (save) {
    utils::write.table(all_assocs, "celiac_gwas_fetched.tsv",
                       sep = "\t", row.names = FALSE, quote = FALSE)
  }

  message("Fetched ", nrow(all_assocs), " total associations.")
  return(all_assocs)
}


#' Fetch GWAS Catalog entries for a given trait
#'
#' @description
#' Programmatically retrieves studies and associations for a user-specified trait
#' from the GWAS Catalog via the \pkg{gwasrapidd} API.
#'
#' @param trait Character; name of the trait (e.g., "celiac disease", "type 2 diabetes").
#' @param save Logical; if TRUE, saves the results to "<trait>_gwas.tsv".
#'
#' @return A `data.frame` with SNP-level associations including mapped genes and p-values.
#' @examples
#' \dontrun{
#' fetch_gwas_catalog("celiac disease")
#' fetch_gwas_catalog("Crohn's disease")
#' }
#' @export
fetch_gwas_catalog <- function(trait, save = FALSE) {
  if (missing(trait)) stop("Please provide a trait name, e.g., fetch_gwas_catalog('celiac disease').")

  if (!requireNamespace("gwasrapidd", quietly = TRUE)) {
    stop("Package 'gwasrapidd' is required. Install via install.packages('gwasrapidd').")
  }

  message("Fetching studies for trait: ", trait)
  studies <- gwasrapidd::get_studies(efo_trait = trait)
  study_ids <- studies@studies$study_id

  if (length(study_ids) == 0) {
    warning("No studies found for trait: ", trait)
    return(NULL)
  }

  assoc_list <- lapply(study_ids, function(id) {
    gwasrapidd::get_associations(study_id = id)
  })

  results <- do.call(rbind, lapply(assoc_list, function(x) {
    data.frame(
      SNP = x@associations$variant_id,
      mapped_gene = x@associations$mapped_gene,
      chromosome = x@associations$chromosome_name,
      position = x@associations$chromosome_position,
      p_value = x@associations$pvalue,
      odds_ratio = x@associations$or,
      source = "GWAS Catalog",
      stringsAsFactors = FALSE
    )
  }))

  if (save) {
    filename <- paste0(gsub(" ", "_", trait), "_gwas.tsv")
    write.table(results, filename, sep = "\t", row.names = FALSE, quote = FALSE)
  }

  return(results)
}


#' Get curated celiac disease GWAS summary table
#'
#' @description
#' Loads a bundled curated CeD GWAS dataset, or optionally fetches the latest
#' GWAS associations from the GWAS Catalog via the gwasrapidd API.
#' If fetched data lack variant coordinates, a warning is shown and
#' the bundled example dataset is returned instead.
#'
#' @param fetch Logical; if TRUE, retrieves CeD data from the GWAS Catalog.
#' @param save Logical; if TRUE, saves fetched results to "celiac_gwas_fetched.tsv".
#'
#' @return A data.frame with SNP, mapped_gene, chromosome, position, p_value, odds_ratio, source.
#' @export
get_celiac_gwas <- function(fetch = FALSE, save = FALSE) {
  if (!fetch) {
    data("celiac_gwas_example", package = "CeDExplorer", envir = environment())
    return(celiac_gwas_example)
  }

  if (!requireNamespace("gwasrapidd", quietly = TRUE)) {
    stop("Please install the 'gwasrapidd' package to fetch live data.")
  }

  message("Fetching celiac disease GWAS data from GWAS Catalog...")
  studies <- gwasrapidd::get_studies(efo_trait = "celiac disease")
  study_ids <- studies@studies$study_id

  if (length(study_ids) == 0) {
    warning("No CeD studies found in GWAS Catalog; returning bundled example data.")
    data("celiac_gwas_example", package = "CeDExplorer", envir = environment())
    return(celiac_gwas_example)
  }

  assoc_list <- lapply(study_ids, function(id) {
    x <- gwasrapidd::get_associations(study_id = id)
    df <- x@associations
    if (nrow(df) == 0) return(NULL)
    data.frame(
      SNP = if ("variant_id" %in% names(df)) df$variant_id else NA,
      mapped_gene = if ("mapped_gene" %in% names(df)) df$mapped_gene else NA,
      chromosome = if ("chromosome_name" %in% names(df)) df$chromosome_name else NA,
      position = if ("chromosome_position" %in% names(df)) df$chromosome_position else NA,
      p_value = if ("pvalue" %in% names(df)) df$pvalue else NA,
      odds_ratio = if ("or" %in% names(df)) df$or else NA,
      source = "GWAS Catalog",
      stringsAsFactors = FALSE
    )
  })

  assoc_list <- assoc_list[!vapply(assoc_list, is.null, logical(1))]
  if (length(assoc_list) == 0) {
    warning("No association data returned; using bundled example dataset instead.")
    data("celiac_gwas_example", package = "CeDExplorer", envir = environment())
    return(celiac_gwas_example)
  }

  all_assocs <- do.call(rbind, assoc_list)

  # Check for missing coordinates
  if (all(is.na(all_assocs$chromosome)) || all(is.na(all_assocs$position))) {
    warning("Fetched GWAS data have no chromosome or position information.\n",
            "Returning bundled example dataset for reproducibility.")
    data("celiac_gwas_example", package = "CeDExplorer", envir = environment())
    return(celiac_gwas_example)
  }

  all_assocs$p_value <- as.numeric(all_assocs$p_value)
  all_assocs$position <- as.numeric(all_assocs$position)
  all_assocs$chromosome <- as.character(all_assocs$chromosome)

  if (save) {
    utils::write.table(all_assocs, "celiac_gwas_fetched.tsv",
                       sep = "\t", row.names = FALSE, quote = FALSE)
  }

  message("Fetched ", nrow(all_assocs), " association records.")
  return(all_assocs)
}

#' Manhattan plot for celiac disease GWAS results
#'
#' @description
#' Generates a Manhattan plot showing genome-wide SNP-level associations
#' (-log10 p-values) for celiac disease (CeD) or a user-supplied GWAS dataset.
#' The function can either use the bundled example dataset or take a custom data.frame
#' with at least `chromosome`, `position`, and `p_value` columns.
#'
#' @param data A `data.frame` containing GWAS summary statistics. If `NULL`, uses
#' the bundled example dataset `celiac_gwas_example`.
#' @param highlight Optional vector of SNP IDs (rsIDs) or gene names to highlight.
#' @param title Character, plot title.
#' @param genomewide_line Numeric p-value threshold for the genome-wide significance line.
#' Default = 5e-8.
#' @param use_qqman Logical; if TRUE, uses qqman package for plotting. If FALSE or
#' qqman not available, uses ggplot2.
#' @return A ggplot2 object (or base R plot if use_qqman = TRUE)
#' @examples
#' \dontrun{
#' # Default with bundled data
#' plot_manhattan()
#'
#' # With custom GWAS table
#' gwas <- get_celiac_gwas(fetch = TRUE)
#' plot_manhattan(gwas, highlight = c("HLA-DQA1", "HLA-DQB1"))
#' }
#' @export
plot_manhattan <- function(data = NULL, highlight = NULL, title = "Celiac disease GWAS",
                           genomewide_line = 5e-8, use_qqman = FALSE) {
  # Load bundled dataset if user didn't provide one
  if (is.null(data)) {
    data("celiac_gwas_example", package = "CeDExplorer", envir = environment())
    data <- celiac_gwas_example
  }

  # Check input validity
  required_cols <- c("chromosome", "position", "p_value")
  if (!all(required_cols %in% colnames(data))) {
    stop("Input data must contain columns: chromosome, position, p_value.")
  }

  # Ensure we have SNP column for highlighting
  if (!"SNP" %in% colnames(data)) {
    data$SNP <- paste0("rs", seq_len(nrow(data)))
  }

  # Try qqman only if explicitly requested and available
  if (use_qqman && requireNamespace("qqman", quietly = TRUE)) {
    message("Using qqman for Manhattan plot...")
    # Prepare qqman-format data
    qq_df <- data.frame(
      SNP = data$SNP,
      CHR = as.numeric(data$chromosome),
      BP = as.numeric(data$position),
      P = as.numeric(data$p_value)
    )

    # Create the plot and return the data frame (qqman doesn't return ggplot)
    qqman::manhattan(qq_df,
                     main = title,
                     genomewideline = -log10(genomewide_line),
                     highlight = highlight,
                     cex = 0.6, col = c("skyblue3", "steelblue4"))
    return(invisible(qq_df))
  } else {
    # Use ggplot2 approach (consistent return type)
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Please install ggplot2 for Manhattan plot functionality.")
    }

    message("Using ggplot2 for Manhattan plot...")

    # Clean data for plotting
    plot_data <- data
    plot_data$chromosome <- as.numeric(plot_data$chromosome)
    plot_data$position <- as.numeric(plot_data$position)
    plot_data$p_value <- as.numeric(plot_data$p_value)

    # Remove NA values
    plot_data <- plot_data[!is.na(plot_data$chromosome) &
                             !is.na(plot_data$position) &
                             !is.na(plot_data$p_value), ]

    # Calculate -log10(p)
    plot_data$logp <- -log10(plot_data$p_value)

    # Create chromosome factor for alternating colors
    plot_data$chr_factor <- as.factor(plot_data$chromosome)

    # Prepare highlight data if provided
    highlight_data <- NULL
    if (!is.null(highlight)) {
      highlight_data <- plot_data[plot_data$SNP %in% highlight |
                                    plot_data$mapped_gene %in% highlight, ]
    }

    # Create ggplot
    plt <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$position, y = .data$logp,
                                                   color = .data$chr_factor)) +
      ggplot2::geom_point(alpha = 0.7, size = 1) +
      ggplot2::geom_hline(yintercept = -log10(genomewide_line),
                          linetype = "dashed", color = "red", linewidth = 0.8) +
      ggplot2::scale_color_manual(values = rep(c("skyblue3", "steelblue4"),
                                               length.out = length(unique(plot_data$chromosome)))) +
      ggplot2::facet_wrap(~chromosome, scales = "free_x", nrow = 1) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(
        title = title,
        x = "Genomic position",
        y = expression(-log[10](p-value))
      ) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.spacing = ggplot2::unit(0.1, "lines"),
        strip.text = ggplot2::element_text(size = 8)
      )

    # Add highlights if provided
    if (!is.null(highlight_data) && nrow(highlight_data) > 0) {
      plt <- plt +
        ggplot2::geom_point(data = highlight_data, color = "red", size = 2, alpha = 0.9) +
        ggplot2::geom_text(data = highlight_data,
                           ggplot2::aes(label = .data$SNP),
                           color = "red", size = 3, vjust = -0.5, check_overlap = TRUE)
    }

    return(plt)
  }
}

