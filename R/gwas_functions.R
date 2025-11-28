# Purpose: Retrieve celiac disease GWAS data from GWAS Catalog or use bundled data
# Author: Yuxi Zhang
# Date: Nov 27, 2025
# Version: 1.0
# Bugs and Issues: Requires internet connection for live data fetching

#' Get celiac disease GWAS data from GWAS Catalog
#'
#' @name getCeliacGwas
#' 
#' @description
#' Retrieve celiac disease GWAS data from GWAS Catalog or use bundled data
#' 
#' @param fetch Logical, whether to fetch from GWAS Catalog or use bundled data
#' @param force_fetch Logical, whether to force fetching even if bundled data exists
#' @return Data frame with GWAS data
#' 
#' @references
#' MacDonald, J. R., Zuber, V., & Ong, J. S. (2020). gwasrapidd: An R 
#' package to query, download and wrangle GWAS Catalog data. 
#' Bioinformatics, 36(2), 649-650. 
#' https://doi.org/10.1093/bioinformatics/btz605
#'
#' Trynka, G., Hunt, K. A., Bockett, N. A., Romanos, J., Mistry, V., 
#' Szperl, A., Bakker, S. F., Bardella, M. T., Bhaw-Rosun, L., Castillejo, G., 
#' de la Concha, E. G., de Almeida, R. C., Dias, K. R., van Diemen, C. C., 
#' Dubois, P. C., Duerr, R. H., Edkins, S., Franke, L., Fransen, K., ... 
#' van Heel, D. A. (2011). Dense genotyping identifies and localizes multiple 
#' common and rare variant association signals in celiac disease. 
#' Nature Genetics, 43(12), 1193-1201. https://doi.org/10.1038/ng.998
#' 
#' @importFrom gwasrapidd get_studies get_associations
#' @importFrom dplyr left_join distinct
#' @export
#' @examples
#' \donttest{
#' # Get bundled example data (default)
#' gwas_data <- getCeliacGwas()
#' head(gwas_data)
#'
#' # Try to fetch live data (requires internet and gwasrapidd package)
#' \dontrun{
#' live_data <- getCeliacGwas(fetch = TRUE)
#' }
#'
#' # Force fetch even if bundled data exists
#' \dontrun{
#' force_data <- getCeliacGwas(forceFetch = TRUE)
#' }
#' }
getCeliacGwas <- function(fetch = FALSE, force_fetch = FALSE) {
  
  if(fetch || force_fetch) {
    message("Fetching celiac disease GWAS data from GWAS Catalog...")
    
    tryCatch({
      # Ensure gwasrapidd is available
      if(!requireNamespace("gwasrapidd", quietly = TRUE)) {
        stop("gwasrapidd package required for fetching GWAS data")
      }
      
      # Get associations for celiac disease
      associations <- gwasrapidd::get_associations(efo_trait = "celiac disease")
      
      if(nrow(associations@associations) == 0) {
        stop("No celiac disease associations found in GWAS Catalog")
      }
      
      cat("DEBUG: Raw associations retrieved:", nrow(associations@associations), "\n")
      
      # Start with basic association data
      gwas_data <- data.frame(
        association_id = associations@associations$association_id,
        p_value = associations@associations$pvalue,
        odds_ratio = associations@associations$or_per_copy_number,
        beta = associations@associations$beta_number,
        standard_error = associations@associations$standard_error,
        source = "GWAS Catalog",
        stringsAsFactors = FALSE
      )
      
      # Extract SNP information from risk_alleles
      if(nrow(associations@risk_alleles) > 0) {
        risk_alleles <- as.data.frame(associations@risk_alleles)
        # Get unique variant information
        snp_info <- risk_alleles[, c("association_id", "variant_id", "risk_allele")]
        names(snp_info)[names(snp_info) == "variant_id"] <- "SNP"
        
        # Merge with main data
        # Using gwasrapidd package (MacDonald et al., 2020)
        gwas_data <- merge(gwas_data, snp_info, by = "association_id", all.x = TRUE)
      } else {
        gwas_data$SNP <- NA_character_
        gwas_data$risk_allele <- NA_character_
      }
      
      # Extract gene information
      if(nrow(associations@genes) > 0) {
        genes <- as.data.frame(associations@genes)
        # For associations with multiple genes, take the first one
        gene_info <- genes[!duplicated(genes$association_id), 
                           c("association_id", "gene_name")]
        names(gene_info)[names(gene_info) == "gene_name"] <- "mapped_gene"
        
        # Merge with main data
        gwas_data <- merge(gwas_data, gene_info, by = "association_id", all.x = TRUE)
      } else {
        gwas_data$mapped_gene <- NA_character_
      }
      
      # Remove rows with missing p-values
      gwas_data <- gwas_data[!is.na(gwas_data$p_value), ]
      
      cat("DEBUG: Final dataset has", nrow(gwas_data), "associations\n")
      cat("DEBUG: Sample data:\n")
      print(head(gwas_data[, c("association_id", "SNP", "mapped_gene", "p_value")]))
      
      if(nrow(gwas_data) == 0) {
        stop("No complete GWAS records found after filtering")
      }
      
      message(sprintf("Successfully retrieved %d associations from GWAS Catalog", nrow(gwas_data)))
      
      return(gwas_data)
      
    }, error = function(e) {
      warning("Failed to fetch from GWAS Catalog: ", e$message, 
              "\nReturning bundled example dataset.")
      return(getCeliacGwas(fetch = FALSE))
    })
    
  } else {
    # Return bundled example data
    curated_data <- data.frame(
      association_id = c("example_1", "example_2", "example_3"),
      SNP = c("rs2187668", "rs7454108", "rs9272346"),
      mapped_gene = c("HLA-DQA1", "HLA-DQB1", "IL18R1"),
      p_value = c(1e-50, 1e-45, 3e-12),
      odds_ratio = c(2.8, 3.1, 1.4),
      beta = c(NA_real_, NA_real_, NA_real_),
      standard_error = c(NA_real_, NA_real_, NA_real_),
      risk_allele = c(NA_character_, NA_character_, NA_character_),
      source = rep("Bundled curated data", 3),
      stringsAsFactors = FALSE
    )
    return(curated_data)
  }
}


#' Plot GWAS summary (volcano plot) since we lack chromosome data for Manhattan plot
#' @param gwas_data GWAS data from getCeliacGwas()
#' @param title Plot title
#' @param significance_level Genome-wide significance level
#' 
#' @references
#' Dudbridge, F., & Gusnanto, A. (2008). Estimation of significance 
#' thresholds for genomewide association scans. *Genetic Epidemiology*, *32*(3), 
#' 227-234. https://doi.org/10.1002/gepi.20297
#' 
#' @import ggplot2
#' @importFrom dplyr arrange desc
#' @export
#' 
#' @examples
#' \donttest{
#' # Get GWAS data and create summary plot
#' gwas_data <- getCeliacGwas()
#' plotGwasSummary(gwas_data)
#'
#' # Custom title and significance level
#' plotGwasSummary(gwas_data, 
#'                 title = "Celiac Disease Association Results",
#'                 significanceLevel = 1e-6)
#' }
plotGwasSummary <- function(gwas_data, title = "Celiac Disease GWAS Summary", 
                              significance_level = 5e-8) {
  
  # Calculate -log10 p-values
  gwas_data$log_p <- -log10(gwas_data$p_value)
  
  # Remove any infinite or missing values
  gwas_data <- gwas_data[is.finite(gwas_data$log_p) & !is.na(gwas_data$log_p), ]
  
  # Create a simple point plot (similar to volcano but without effect size)
  p <- ggplot(gwas_data, aes(x = 1:nrow(gwas_data), y = log_p)) +
    geom_point(alpha = 0.6, size = 1, color = "blue") +
    # Genome-wide significance threshold (Dudbridge & Gusnanto, 2008)
    geom_hline(yintercept = -log10(significance_level), 
               color = "red", linetype = "dashed", linewidth = 0.8) +
    labs(
      title = title,
      x = "SNP Index",
      y = "-log10(p-value)"
    ) +
    theme_minimal()
  
  return(p)
}

#' Plot top GWAS hits with gene labels
#' @param gwas_data GWAS data from getCeliacGwas()
#' @param top_n Number of top hits to display
#' @param title Plot title
#' 
#' @references
#' Sollid, L. M., & Jabri, B. (2013). Triggers and drivers of autoimmunity: 
#' Lessons from coeliac disease. Nature Reviews Immunology, 13(4), 294-302.
#' https://doi.org/10.1038/nri3407
#' 
#' @import ggplot2
#' @export
#' 
#' @examples
#' \donttest{
#' # Plot top 20 hits
#' gwas_data <- getCeliacGwas()
#' plotTopGwasHits(gwas_data)
#'
#' # Plot top 10 hits with custom title
#' plotTopGwasHits(gwas_data, topN = 10, 
#'                 title = "Top 10 Celiac Disease GWAS Hits")
#' }
plotTopGwasHits <- function(gwas_data, top_n = 20, title = "Top Celiac Disease GWAS Hits") {
  
  # Sort by p-value and take top hits using base R
  gwas_data_sorted <- gwas_data[order(gwas_data$p_value), ]
  top_hits <- head(gwas_data_sorted, top_n)
  
  # Create meaningful labels
  top_hits$label <- ifelse(!is.na(top_hits$mapped_gene) & top_hits$mapped_gene != "",
                           paste0(top_hits$mapped_gene, "\n(", top_hits$SNP, ")"),
                           top_hits$SNP)
  
  # Order by significance for plotting (reverse order for bar plot)
  top_hits <- top_hits[order(-top_hits$p_value), ]  # Reverse for plot
  
  p <- ggplot(top_hits, aes(x = reorder(label, -log10(p_value)), 
                            y = -log10(p_value))) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = -log10(5e-8), color = "red", 
               linetype = "dashed", linewidth = 0.8) +
    labs(
      title = title,
      x = NULL,
      y = "-log10(p-value)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major.x = element_blank()
    )
  
  return(p)
}
