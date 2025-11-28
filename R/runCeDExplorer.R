# Launch CeDExplorer Shiny Application
#
# Starts the interactive Shiny web application for exploring
# celiac disease gene expression data.

#' Run CeDExplorer Shiny Application
#'
#' @name runCeDExplorerApp
#'
#' @details
#' The Shiny application provides:
#' - Interactive visualization of gene expression differences between 
#'   celiac disease (CeD) and control samples
#' - Support for user-uploaded CSV data in long format
#' - Built-in demo data for testing
#' - Dynamic gene selection based on available data
#' - Option for boxplots or violin plots
#'
#' @section Demo Data:
#' The app includes built-in demo data (`example_expression`) containing
#' expression values for immune-related genes. Users can also upload their
#' own CSV files with the following required columns:
#' - `gene`: Gene symbol (e.g., "HLA-DQA1")
#' - `sample`: Sample identifier  
#' - `condition`: "CeD" or "Control"
#' - `expression`: Numeric expression values
#'
#' Demo data location: Use built-in `example_expression` dataset or
#' create CSV files following the required format.
#'
#' @return Shiny application object
#'
#' @export
#' @importFrom shiny runApp
#' @examples
#' \dontrun{
#' # Launch the Shiny app
#' runCeDExplorerApp()
#' 
#' # The app will open in your default web browser
#' # Use the demo data or upload your own CSV files
#' }
runCeDExplorerApp <- function() {
  app_dir <- system.file("shiny-scripts", package = "CeDExplorer")
  if (app_dir == "") {
    stop("Could not find Shiny app directory. Try re-installing CeDExplorer.")
  }
  
  shiny::runApp(app_dir, display.mode = "normal")
}