#' Example celiac disease GWAS dataset
#'
#' @description
#' A small curated example dataset of celiac disease-associated SNPs and mapped genes,
#' used to demonstrate CeDExplorer functions.
#'
#' @format A data frame with 7 columns:
#' \describe{
#'   \item{SNP}{Variant identifier (rsID)}
#'   \item{mapped_gene}{Nearest or implicated gene symbol}
#'   \item{chromosome}{Chromosome number}
#'   \item{position}{Base pair position (GRCh37)}
#'   \item{p_value}{Association p-value}
#'   \item{odds_ratio}{Reported odds ratio}
#'   \item{source}{Data source description}
#' }
#'
#' @source Manually curated subset of celiac disease GWAS summary data from the GWAS Catalog.
#' @usage data(celiac_gwas_example)
"celiac_gwas_example"
