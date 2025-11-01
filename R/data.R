# Data documentation for CeDExplorer package
#' Curated Protein-Protein Interaction Data for CeDExplorer
#'
#' A curated dataset containing protein-protein interactions for
#' celiac disease-related genes. Includes both known biological interactions
#' and predicted associations.
#'
#' @format A list with two components:
#' \describe{
#'   \item{edges}{A data frame with 25 rows and 4 variables:
#'     \describe{
#'       \item{gene1}{First interacting gene}
#'       \item{gene2}{Second interacting gene}
#'       \item{score}{Interaction confidence score (0-1)}
#'       \item{evidence}{Type of evidence: complex, signaling, costimulation, predicted}
#'     }
#'   }
#'   \item{nodes}{A data frame with 28 rows and 3 variables:
#'     \describe{
#'       \item{gene}{Gene symbol}
#'       \item{type}{Gene type: MHC, Cytokine/Chemokine, Transcription Factor, Cell Surface, Other}
#'       \item{importance}{Relative importance score (0-1)}
#'     }
#'   }
#' }
#' @source Curated from known biological pathways and simulated data
#' @examples
#' data(ppi_data)
#' head(ppi_data$edges)
#' head(ppi_data$nodes)
"ppi_data"
