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
#' @references
#' Chang, W., Cheng, J., Allaire, J., Sievert, C., Schloerke, B., Xie, Y., 
#' Allen, J., McPherson, J., Dipert, A., & Borges, B. (2023). shiny: Web 
#' Application Framework for R. R package version 1.7.5. 
#' https://CRAN.R-project.org/package=shiny
#'
#' Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. 
#' Springer-Verlag New York. ISBN 978-3-319-24277-4.
#' 
#' Hintze, J. L., & Nelson, R. D. (1998). Violin plots: A box plot-density 
#' trace synergism. The American Statistician, 52(2), 181-184. 
#' https://doi.org/10.1080/00031305.1998.10480559
#'
#' Wilcoxon, F. (1945). Individual comparisons by ranking methods. 
#' Biometrics Bulletin, 1(6), 80-83. 
#' https://doi.org/10.2307/3001968
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
  
  return (app)
}

# [END] 