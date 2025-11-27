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
        helpText("Required columns: 'gene', 'sample', 'condition', 'expression'")
      ),
      
      # Gene Selection (initially empty, will be populated by server)
      selectInput("gene", "Select Gene to Visualize:",
                  choices = NULL,  # Will be populated dynamically
                  selected = NULL),
      
      # Plot Options
      checkboxInput("violin", "Use Violin Plot", value = TRUE),
      
      # Demo Data Info
      conditionalPanel(
        condition = "input.data_source == 'demo'",
        hr(),
        helpText("Demo data: example_expression with immune-related genes")
      )
    ),
    
    mainPanel(
      width = 9,
      
      # Plot Output
      plotOutput("expression_plot", height = "500px"),
      
      # Data Summary
      verbatimTextOutput("data_summary"),
      
      # File Format Instructions
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
Gene1,S1,CeD,8.2
Gene1,S2,Control,5.1
Gene2,S1,CeD,7.5
Gene2,S2,Control,4.8"
          )
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive expression for data
  plot_data <- reactive({
    if (input$data_source == "demo") {
      # Use demo data
      data <- CeDExplorer::example_expression
      return(data)
    } else {
      # Use uploaded data
      req(input$file_upload)
      
      # Read uploaded file
      user_data <- read.csv(input$file_upload$datapath)
      
      # Validate required columns
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
  
  # Update gene selection based on available data
  observe({
    data <- plot_data()
    req(data)
    
    available_genes <- sort(unique(data$gene))
    
    # Update the select input
    updateSelectInput(session, "gene", 
                      choices = available_genes,
                      selected = available_genes[1])
  })
  
  # Generate plot automatically when inputs change
  output$expression_plot <- renderPlot({
    data <- plot_data()
    req(data, input$gene)
    
    # Additional validation: check if selected gene exists
    if (!input$gene %in% data$gene) {
      return(ggplot2::ggplot() + 
               ggplot2::geom_text(ggplot2::aes(x = 0.5, y = 0.5, 
                                               label = paste("Gene", input$gene, "not found in data"))) +
               ggplot2::theme_void())
    }
    
    tryCatch({
      p <- CeDExplorer::expressionBoxplot(
        gene = input$gene,
        dataset = data,
        violin = input$violin
      )
      return(p)
    }, error = function(e) {
      showNotification(
        paste("Error generating plot:", e$message),
        type = "error",
        duration = 10
      )
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
    cat("Available genes:", paste(sort(unique(data$gene)), collapse = ", "), "\n")
    cat("Samples per condition:\n")
    print(table(data$condition))
  })
}

# Create Shiny app
shinyApp(ui = ui, server = server)