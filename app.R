# Load necessary packages
library(shiny)
library(g3viz)
library(dplyr)
library(DT)
library(shinydashboard)
library(shinyjs)
library(biomaRt)
library(memoise)

# Define the file path
file_path <- "mutations.tsv"

# Define the user interface
ui <- dashboardPage(
  dashboardHeader(title = "Protein Changes Visualization"),
  dashboardSidebar(
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
    div(style = "display: flex; justify-content: space-between; padding: 10px;",
        downloadButton("download_data", "Download Data")
    ),
    hr(),
    sliderInput("aa_range", "Amino Acid Position Range:", min = 0, max = 1000, value = c(0, 1000))
  ),
  dashboardBody(
    useShinyjs(),
    tabsetPanel(
      tabPanel("Lollipop Plot", g3LollipopOutput("lollipopPlot", height = "600px")),
      tabPanel("Data Table", DTOutput("data_table")),
      tabPanel("Summary", verbatimTextOutput("summary_stats"))
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  # Define the blacklist vector
  blacklisted_genes <- c("C1orf87")  # Add more genes to this list as needed
  
  # Memoize biomaRt functions for caching
  memoised_useMart <- memoise(useMart)
  memoised_getBM <- memoise(getBM)
  
  # Use reactiveFileReader to monitor the TSV file
  data <- reactiveFileReader(
    intervalMillis = 1000,
    session = session,
    filePath = file_path,
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
    
    df <- data() %>% 
      filter(Hugo_Symbol == input$gene)
    
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
  
  # Fetch domain information with caching and error handling
  get_domain_info <- reactive({
    req(input$gene)
    
    tryCatch({
      ensembl <- memoised_useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      
      # First, get the Ensembl gene ID
      gene_info <- memoised_getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
                                  filters = "hgnc_symbol",
                                  values = input$gene,
                                  mart = ensembl)
      
      if(nrow(gene_info) == 0) {
        showNotification(paste("No information found for gene:", input$gene), type = "warning")
        return(NULL)
      }
      
      # Then, use the Ensembl gene ID to fetch domain information
      domains <- memoised_getBM(attributes = c("ensembl_gene_id", "interpro", "interpro_start", "interpro_end"),
                                filters = "ensembl_gene_id",
                                values = gene_info$ensembl_gene_id,
                                mart = ensembl)
      
      # Rename columns to match g3Lollipop expectations
      colnames(domains) <- c("ensembl_gene_id", "protein_domain_id", "protein_domain_start", "protein_domain_end")
      
      return(domains)
      
    }, error = function(e) {
      showNotification(paste("Error fetching domain information:", e$message), type = "error")
      return(NULL)
    })
  })
  
  # Update AA range slider based on selected gene and mutation type
  observe({
    req(filtered_data())
    
    # Check if there are any rows in filtered data
    if (nrow(filtered_data()) > 0) {
      max_aa <- max(filtered_data()$AA_Position, na.rm = TRUE)
      # Update slider with valid range
      updateSliderInput(session, "aa_range", max = max_aa, value = c(0, max_aa))
    } else {
      # Set slider to default range or disable it if no data is available
      updateSliderInput(session, "aa_range", max = 1000, value = c(0, 1000))
    }
  })
  
  # Filtered data with AA range
  filtered_data_with_range <- reactive({
    req(filtered_data(), input$aa_range)
    filtered_data() %>%
      filter(AA_Position >= input$aa_range[1] & AA_Position <= input$aa_range[2])
  })
  
  # Render the lollipop plot with cBioPortal theme
  output$lollipopPlot <- renderG3Lollipop({
    req(filtered_data_with_range())
    
    domains <- get_domain_info()
    
    if (is.null(filtered_data_with_range()) || nrow(filtered_data_with_range()) == 0) {
      showNotification("No data available for the selected criteria.", type = "warning")
      return(NULL)
    }
    
    # Prepare domain data for g3Lollipop
    domain_data <- NULL
    
    if (!is.null(domains) && nrow(domains) > 0) {
      domain_data <- list(
        pfam = lapply(1:nrow(domains), function(i) {
          list(
            name = domains$protein_domain_id[i],
            start = domains$protein_domain_start[i],
            end = domains$protein_domain_end[i]
          )
        })
      )
    }
    
    plot.options <- g3Lollipop.theme(
      theme.name = "cbioportal",
      title.text = paste(input$gene, "Protein Changes"),
      y.axis.label = "# of Mutations"
    )
    
    plot.options$chart$width <- 1200
    plot.options$chart$height <- 800
    
    if (!is.null(domain_data)) { 
      plot.options$chart$domains <- domain_data 
    }
    
    g3Lollipop(
      filtered_data_with_range(),
      gene.symbol = input$gene,
      gene.symbol.col = "Hugo_Symbol",
      aa.pos.col = "AA_Position",
      protein.change.col = "Protein_Change",
      factor.col = "mutation_group",
      plot.options = plot.options
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
  
  # Download handler for data
  output$download_data <- downloadHandler(
    filename = function() {
      paste(input$gene, "_protein_changes_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_data_with_range(), file, row.names = FALSE)
    }
  )

}

# Run Shiny app
shinyApp(ui = ui, server = server)
