#' Launch CeDExplorer Shiny Application
#'
#' Interactive web application for exploring gene expression in celiac disease.
#' Users can upload CSV files with expression data or use built-in demo data.
#'
#' @details
#' The Shiny app provides:
#' - Gene expression visualization using boxplots or violin plots
#' - Support for user-uploaded CSV data
#' - Built-in demo data for testing
#' - Data summary statistics
#'
#' File format requirements for upload:
#' - CSV format with samples as rows
#' - Required columns: 'sample', 'condition' 
#' - Additional columns: gene expression values
#' - Condition values: 'CeD' or 'Control'
#'
#' Demo data location: inst/extdata/demo_expression_data.csv
#'
#' @export
#' @examples
#' \dontrun{
#' # Launch the Shiny app
#' runCeDExplorerApp()
#' }
runCeDExplorerApp <- function() {
  app_dir <- system.file("shiny-scripts", package = "CeDExplorer")
  if (app_dir == "") {
    stop("Could not find Shiny app directory. Try re-installing CeDExplorer.")
  }
  
  shiny::runApp(app_dir, display.mode = "normal")
}