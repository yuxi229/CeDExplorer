#' Example celiac disease GWAS dataset
#'
#' A curated dataset containing GWAS summary statistics for celiac disease.
#'
#' @format A data frame with 3 rows and 7 variables:
#' \describe{
#'   \item{SNP}{Variant identifier}
#'   \item{mapped_gene}{Mapped gene}
#'   \item{chromosome}{Chromosome}
#'   \item{position}{Genomic position}
#'   \item{p_value}{Association p-value}
#'   \item{odds_ratio}{Odds ratio}
#'   \item{source}{Data source}
#' }
#' @examples
#' data(celiac_gwas_example)
#' head(celiac_gwas_example)
"celiac_gwas_example"

#' Example expression dataset
#'
#' A simulated dataset containing gene expression data for celiac disease and control samples.
#'
#' @format A data frame with 240 rows and 4 variables:
#' \describe{
#'   \item{gene}{Gene symbol}
#'   \item{sample}{Sample identifier}
#'   \item{condition}{Condition, either 'CeD' or 'Control'}
#'   \item{expression}{Simulated expression value}
#' }
#' @examples
#' data(example_expression)
#' head(example_expression)
"example_expression"

#' Example expression matrix
#'
#' A simulated gene expression matrix for demonstration purposes.
#'
#' @format A matrix with 12 rows (genes) and 10 columns (samples)
#' @examples
#' data(example_expression_matrix)
#' head(example_expression_matrix)
"example_expression_matrix"

#' Example sample metadata
#'
#' Sample annotation data matching the example expression datasets.
#'
#' @format A data frame with 10 rows and 2 variables:
#' \describe{
#'   \item{sample}{Sample identifier}
#'   \item{condition}{Condition, either 'CeD' or 'Control'}
#' }
#' @examples
#' data(example_metadata)
#' head(example_metadata)
"example_metadata"
#' Curated Protein-Protein Interaction Data
#'
#' A curated dataset containing protein-protein interactions for
#' celiac disease-related genes.
#'
#' @format A list with two components:
#' \describe{
#'   \item{edges}{A data frame with interaction data:
#'     \describe{
#'       \item{gene1}{First interacting gene}
#'       \item{gene2}{Second interacting gene}
#'       \item{score}{Interaction confidence score (0-1)}
#'       \item{evidence}{Type of evidence: complex, signaling, costimulation, predicted}
#'     }
#'   }
#'   \item{nodes}{A data frame with gene information:
#'     \describe{
#'       \item{gene}{Gene symbol}
#'       \item{type}{Gene type: MHC, Cytokine/Chemokine, Transcription Factor, Cell Surface, Other}
#'       \item{importance}{Relative importance score (0-1)}
#'     }
#'   }
#' }
#' @examples
#' data(ppi_data)
#' head(ppi_data$edges)
#' head(ppi_data$nodes)
"ppi_data"