# Load necessary packages
library(shiny)
library(g3viz)
library(dplyr)
library(DT)
library(shinydashboard)
library(shinyjs)

# Define the file path for the variants dataset
variants_file_path <- "mutations.tsv"

# Define the user interface
ui <- dashboardPage(
  dashboardHeader(title = "Protein Changes"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Variants", tabName = "variants_tab", icon = icon("dna")),
      menuItem("Single-Cell", tabName = "singlecell_tab", icon = icon("microscope"))
    )
  ),
  dashboardBody(
    useShinyjs(),
    tabItems(
      # Variants Tab
      tabItem(
        tabName = "variants_tab",
        h2("Variants Visualization"),
        fluidRow(
          column(4,
                 selectizeInput("gene", "Select Gene:", choices = NULL, options = list(
                   placeholder = 'Please select a gene',
                   onInitialize = I('function() { this.setValue(""); }')
                 )),
                 checkboxGroupInput("disease", "Select Disease:", 
                                    choices = c("All", "NK", "T"),
                                    selected = "All"),
                 selectInput("mutation_type", "Mutation Type:", 
                             choices = c("All", "Silent", "Non-silent"), 
                             selected = "All"),
                 sliderInput("aa_range", "Amino Acid Position Range:", min = 0, max = 1000, value = c(0, 1000)),
                 downloadButton("download_data", "Download Data")
          ),
          column(8,
                 tabsetPanel(
                   tabPanel("Lollipop Plot", g3LollipopOutput("lollipopPlot", height = "400px")),
                   tabPanel("Data Table", DTOutput("data_table")),
                   tabPanel("Summary", verbatimTextOutput("summary_stats"))
                 )
          )
        )
      ),
      
      # Single-Cell Tab
      tabItem(
        tabName = "singlecell_tab",
        h2("Single-Cell Data Placeholder")
      )
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  # Variants Logic
  # Define the blacklist vector
  blacklisted_genes <- c("C1orf87")  # Add more genes to this list as needed
  
  # Use reactiveFileReader to monitor the TSV file
  data <- reactiveFileReader(
    intervalMillis = 1000,
    session = session,
    filePath = variants_file_path,
    readFunc = function(filePath) {
      tryCatch({
        test <- read.delim(filePath, sep = "\t", stringsAsFactors = FALSE)
        
        # Filter out blacklisted genes
        test <- test %>% filter(!Hugo_Symbol %in% blacklisted_genes)
        
        test$AA_Position <- sapply(test$Protein_Change, function(x) {
          pos <- as.numeric(gsub("p\\.([A-Z*])([0-9]+)([A-Z*]|fs.*|del.*|ins.*)", "\\2", x))
          if (is.na(pos)) {
            return(0)  # Assign 0 for missing positions
          } else {
            return(pos)
          }
        })
        # Combine GD with T
        test$disease[test$disease == "GD"] <- "T (GD)"
        # Group mutations into Silent and Non-silent
        test$mutation_group <- ifelse(test$Variant_Classification == "Silent", "Silent", "Non-silent")
        test
      }, error = function(e) {
        showNotification(paste("Error reading file:", e$message), type = "error")
        return(NULL)
      })
    }
  )
  
  # Update gene choices
  observe({
    req(data())
    gene_choices <- unique(data()$Hugo_Symbol)
    updateSelectizeInput(session, "gene", choices = gene_choices, server = TRUE)
  })
  
  # Reactive expression to filter data based on selected gene, mutation type, and diseases
  filtered_data <- reactive({
    req(input$gene, input$disease, data())
    if (input$gene %in% blacklisted_genes) {
      showNotification(paste("Gene", input$gene, "is blacklisted and cannot be visualized."), type = "warning")
      return(NULL)
    }
    
    df <- data() %>% filter(Hugo_Symbol == input$gene)
    
    # Filter by mutation type
    if (input$mutation_type != "All") {
      df <- df %>% filter(mutation_group == input$mutation_type)
    }
    
    # Filter by disease
    if ("All" %in% input$disease) {
      df
    } else {
      selected_diseases <- input$disease
      if ("T" %in% selected_diseases) {
        selected_diseases <- c(selected_diseases, "T (GD)")
      }
      df %>% filter(disease %in% selected_diseases)
    }
  })
  
  # Filtered data with AA range
  filtered_data_with_range <- reactive({
    req(filtered_data(), input$aa_range)
    filtered_data() %>%
      filter(AA_Position >= input$aa_range[1] & AA_Position <= input$aa_range[2])
  })
  
  # Update AA range slider based on selected gene and mutation type
  observe({
    req(filtered_data())
    if (nrow(filtered_data()) > 0) {
      max_aa <- max(filtered_data()$AA_Position, na.rm = TRUE)
      updateSliderInput(session, "aa_range", min = 0, max = max_aa, value = c(0, max_aa))
      enable("aa_range")  # Enable the slider
    } else {
      updateSliderInput(session, "aa_range", min = 0, max = 1000, value = c(0, 1000))
      disable("aa_range")  # Disable the slider using shinyjs
    }
  })
  
  # Render the lollipop plot
  output$lollipopPlot <- renderG3Lollipop({
    req(input$gene)  # Ensure a gene is selected
    filtered_data <- filtered_data_with_range()
    
    # If no data exists for the selected filters, return NULL and show a notification
    if (is.null(filtered_data) || nrow(filtered_data) == 0) {
      return(NULL)
    }
    
    # Render the lollipop plot with available data
    g3Lollipop(
      filtered_data,
      gene.symbol = input$gene,
      gene.symbol.col = "Hugo_Symbol",
      aa.pos.col = "AA_Position",
      protein.change.col = "Protein_Change",
      factor.col = "mutation_group",
      plot.options = g3Lollipop.theme(theme.name = "cbioportal")
    )
  })
  
  # Render the data table
  output$data_table <- renderDT({
    datatable(filtered_data_with_range(), options = list(pageLength = 10))
  })
  
  # Render summary statistics
  output$summary_stats <- renderPrint({
    req(filtered_data_with_range())
    cat("Summary Statistics:\n\n")
    cat("Total Mutations:", nrow(filtered_data_with_range()), "\n")
    cat("Unique Protein Changes:", length(unique(filtered_data_with_range()$Protein_Change)), "\n")
    cat("\nMutations by Type:\n")
    print(table(filtered_data_with_range()$mutation_group))
    cat("\nMutations by Disease:\n")
    print(table(filtered_data_with_range()$disease))
  })
  
  # Download handler for filtered data
  output$download_data <- downloadHandler(
    filename = function() {
      paste(input$gene, "_protein_changes_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_data_with_range(), file, row.names = FALSE)
    }
  )
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

