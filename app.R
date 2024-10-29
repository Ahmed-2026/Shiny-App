# Load necessary packages
library(shiny)
library(g3viz)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(shinydashboard)
library(shinyjs)

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
    conditionalPanel(
      condition = "input.disease.includes('T')",
      checkboxInput("include_gd", "Include GD in T", value = TRUE)
    ),
    selectInput("mutation_type", "Mutation Type:", 
                choices = c("All", "Missense_Mutation", "Nonsense_Mutation", "Silent"), 
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
      tabPanel("Interactive Plot", plotlyOutput("interactivePlot", height = "800px")),  # Increased height
      tabPanel("Data Table", DTOutput("data_table")),
      tabPanel("Summary", verbatimTextOutput("summary_stats"))
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  # Use reactiveFileReader to monitor the TSV file
  data <- reactiveFileReader(
    intervalMillis = 1000, # Check every second
    session = session,
    filePath = file_path,
    readFunc = function(filePath) {
      tryCatch({
        test <- read.delim(filePath, sep = "\t", stringsAsFactors = FALSE)
        test$AA_Position <- sapply(test$Protein_Change, function(x) {
          pos <- as.numeric(gsub("p\\.([A-Z*])(\\d+)([A-Z*])", "\\2", x))
          if (is.na(pos)) {
            return(0)  # Assign 0 for missing positions
          } else {
            return(pos)
          }
        })
        # Combine GD with T
        test$disease[test$disease == "GD"] <- "T (GD)"
        print(paste("Successfully read", nrow(test), "rows from", filePath))
        test
      }, error = function(e) {
        print(paste("Error reading file:", e$message))
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
    df <- data() %>% 
      filter(Hugo_Symbol == input$gene) %>%
      mutate(mutation_type = Variant_Classification)  # Ensure this column is used for mutation type
    
    # Filter by mutation type
    if (input$mutation_type != "All") {
      df <- df %>% filter(Variant_Classification == input$mutation_type)
    }
    
    # Filter by disease
    if ("All" %in% input$disease) {
      df
    } else {
      selected_diseases <- input$disease
      if ("T" %in% selected_diseases && !input$include_gd) {
        selected_diseases <- selected_diseases[selected_diseases != "T (GD)"]
      }
      df %>% filter(disease %in% selected_diseases)
    }
  })
  
  # Check if the selected gene has a stop codon or lacks UniProt mapping
  check_gene_status <- reactive({
    req(filtered_data())
    nonsense_mutations <- any(grepl("\\*", filtered_data()$Protein_Change))
    no_uniprot_mapping <- input$gene %in% c("C1ORF87")  # Example, add more genes if needed
    list(nonsense = nonsense_mutations, no_mapping = no_uniprot_mapping)
  })
  
  # Update AA range slider based on selected gene
  observe({
    req(filtered_data())
    max_aa <- max(filtered_data()$AA_Position)
    updateSliderInput(session, "aa_range", max = max_aa, value = c(0, max_aa))
  })
  
  # Filtered data with AA range
  filtered_data_with_range <- reactive({
    req(filtered_data(), input$aa_range)
    filtered_data() %>%
      filter(AA_Position >= input$aa_range[1] & AA_Position <= input$aa_range[2])
  })
  
  # Automatically switch to interactive plot if no UniProt mapping or stop codon
  observe({
    status <- check_gene_status()
    if (status$no_mapping || status$nonsense) {
      updateTabsetPanel(session, "tabs", selected = "Interactive Plot")
    }
  })
  
  # Perform clustering on the filtered data
  clustered_data <- reactive({
    req(filtered_data_with_range())
    df <- filtered_data_with_range()
    
    # Ensure enough rows for clustering
    if (nrow(df) < 2) {
      return(df)
    }
    
    # Perform clustering (k-means clustering for demonstration)
    set.seed(123)  # For reproducibility
    clusters <- kmeans(df$AA_Position, centers = 3)  # Adjust the number of clusters as needed
    df$Cluster <- as.factor(clusters$cluster)  # Assign clusters to data
    return(df)
  })
  
  # Render the lollipop plot with the cBioPortal theme
  output$lollipopPlot <- renderG3Lollipop({
    req(filtered_data_with_range())
    if (is.null(filtered_data_with_range()) || nrow(filtered_data_with_range()) == 0) {
      return(NULL)
    }
    
    plot.options <- g3Lollipop.theme(
      theme.name = "cbioportal", 
      title.text = paste(input$gene, "Protein Changes (cbioportal theme)"),
      y.axis.label = "# of Mutations"
    )
    
    g3Lollipop(filtered_data_with_range(), 
               gene.symbol = input$gene, 
               gene.symbol.col = "Hugo_Symbol", 
               aa.pos.col = "AA_Position", 
               protein.change.col = "Protein_Change", 
               factor.col = "Variant_Classification",  # Mutations will be colored by their type
               plot.options = plot.options)
  })
  
  # Render the interactive plot with clustering
  output$interactivePlot <- renderPlotly({
    req(clustered_data())
    if (is.null(clustered_data()) || nrow(clustered_data()) == 0) {
      return(NULL)
    }
    
    # Abbreviate disease labels for better visualization
    clustered_data()$disease <- gsub("NK", "N", clustered_data()$disease)
    clustered_data()$disease <- gsub("T", "T", clustered_data()$disease)  # You can customize abbreviations
    
    p <- ggplot(clustered_data(), aes(x = AA_Position, y = disease, color = Cluster)) +  # Color by cluster
      geom_point(size = 3) +
      theme_minimal() +
      labs(title = paste(input$gene, "Protein Changes (Clustered)"),
           x = "Amino Acid Position",
           y = "Disease") +
      theme(plot.title = element_text(size = 18, face = "bold"))
    
    ggplotly(p, tooltip = c("x", "y", "color", "Protein_Change"))
  })
  
  # Render the data table
  output$data_table <- renderDT({
    datatable(filtered_data_with_range(), options = list(pageLength = 10))
  })
  
  # Add a summary statistics output
  output$summary_stats <- renderPrint({
    req(filtered_data_with_range())
    cat("Summary Statistics:\n\n")
    cat("Total Mutations:", nrow(filtered_data_with_range()), "\n")
    cat("Unique Protein Changes:", length(unique(filtered_data_with_range()$Protein_Change)), "\n")
    cat("\nMutations by Type:\n")
    print(table(filtered_data_with_range()$Variant_Classification))
    cat("\nMutations by Disease:\n")
    print(table(filtered_data_with_range()$disease))
  })
  
  # Download handler for the data
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
