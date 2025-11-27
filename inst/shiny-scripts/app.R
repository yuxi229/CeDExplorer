#'
#' CeDExplorer Shiny Application
#' 
#' Interactive web application for exploring gene expression in celiac disease.
#' Users can upload their own expression data or use demo data.
#'
library(shiny)
library(CeDExplorer)
library(ggplot2)

# Define UI
ui <- fluidPage(
  titlePanel("CeDExplorer: Gene Expression Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # Data Source Selection
      radioButtons("data_source", "Data Source:",
                   choices = c("Use Demo Data" = "demo", 
                               "Upload Your Data" = "upload"),
                   selected = "demo"),
      
      # File Upload (only shown when upload is selected)
      conditionalPanel(
        condition = "input.data_source == 'upload'",
        fileInput("file_upload", "Upload Expression Data (CSV):",
                  accept = c(".csv"),
                  placeholder = "Choose CSV file"),
        helpText("File format: CSV with samples as rows, genes as columns"),
        helpText("Required columns: 'sample', 'condition', plus gene columns")
      ),
      
      # Gene Selection
      selectInput("gene", "Select Gene to Visualize:",
                  choices = c("HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1"),
                  selected = "HLA-DQA1"),
      
      # Plot Options
      checkboxInput("violin", "Use Violin Plot", value = TRUE),
      
      # Action Button
      actionButton("update_plot", "Generate Plot", 
                   class = "btn-primary"),
      
      # Demo Data Info
      conditionalPanel(
        condition = "input.data_source == 'demo'",
        hr(),
        helpText("Demo data: example_expression dataset in long format with 6 HLA genes")
      )
    ),
    
    mainPanel(
      width = 9,
      
      # Plot Output
      plotOutput("expression_plot", height = "500px"),
      
      # Data Summary
      verbatimTextOutput("data_summary"),
      
      # File Format Instructions
      # In the UI, update the file format instructions:
      conditionalPanel(
        condition = "input.data_source == 'upload'",
        wellPanel(
          h4("File Format Requirements:"),
          tags$ul(
            tags$li("File type: CSV (comma-separated values)"),
            tags$li("Required columns: 'gene', 'sample', 'condition', 'expression'"),
            tags$li("Gene: Gene symbol (e.g., HLA-DQA1)"),
            tags$li("Sample: Sample identifier"),
            tags$li("Condition: 'CeD' or 'Control'"),
            tags$li("Expression: Numeric expression values")
          ),
          h5("Example structure (LONG format):"),
          tags$pre(
            "gene,sample,condition,expression
HLA-DQA1,S1,CeD,8.2
HLA-DQA1,S2,Control,5.1
HLA-DQB1,S1,CeD,7.5
HLA-DQB1,S2,Control,4.8"
          )
        )
      )
    )
  )
)

# Define server logic
# Define server logic
server <- function(input, output, session) {
  
  # Reactive expression for data
  plot_data <- reactive({
    if (input$data_source == "demo") {
      # Use example_expression directly - it's already in the correct LONG format
      cat("=== DEBUG: Using demo data in long format ===\n")
      return(CeDExplorer::example_expression)
    } else {
      # Use uploaded data
      req(input$file_upload)
      
      # Read uploaded file
      user_data <- read.csv(input$file_upload$datapath)
      
      # Validate required columns for LONG format
      required_cols <- c("gene", "sample", "condition", "expression")
      missing_cols <- setdiff(required_cols, colnames(user_data))
      
      if (length(missing_cols) > 0) {
        showNotification(
          paste("Missing required columns:", paste(missing_cols, collapse = ", ")),
          type = "error",
          duration = 10
        )
        return(NULL)
      }
      
      return(user_data)
    }
  })
  
  # Generate plot automatically when inputs change
  output$expression_plot <- renderPlot({
    # Require that we have data and a gene selected
    data <- plot_data()
    req(data, input$gene)
    
    cat("=== DEBUG: Starting plot generation ===\n")
    cat("Selected gene:", input$gene, "\n")
    cat("Data dimensions:", dim(data), "\n")
    cat("Data columns:", paste(colnames(data), collapse = ", "), "\n")
    
    tryCatch({
      # Generate the plot - the function will filter for the selected gene
      p <- CeDExplorer::expressionBoxplot(
        gene = input$gene,
        dataset = data,
        violin = input$violin
      )
      cat("=== DEBUG: Plot generated successfully ===\n")
      return(p)
    }, error = function(e) {
      cat("=== DEBUG: Plot error:", e$message, "===\n")
      showNotification(
        paste("Error generating plot:", e$message),
        type = "error",
        duration = 10
      )
      # Return a simple error plot
      ggplot2::ggplot() + 
        ggplot2::geom_text(ggplot2::aes(x = 0.5, y = 0.5, 
                                        label = paste("Error:", e$message))) +
        ggplot2::theme_void()
    })
  })
  
  # Data summary
  output$data_summary <- renderPrint({
    data <- plot_data()
    req(data)
    
    cat("Data Summary:\n")
    cat("=============\n")
    cat("Total observations:", nrow(data), "\n")
    cat("Unique samples:", length(unique(data$sample)), "\n")
    cat("Conditions:", paste(unique(data$condition), collapse = ", "), "\n")
    cat("Available genes:", paste(unique(data$gene), collapse = ", "), "\n")
    cat("Samples per condition:\n")
    print(table(data$condition))
  })
}

# Create Shiny app
shinyApp(ui = ui, server = server)