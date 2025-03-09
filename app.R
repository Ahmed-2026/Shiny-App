# Load necessary packages
library(shiny)
library(g3viz)
library(dplyr)
library(DT)
library(shinydashboard)
library(shinyjs)
library(Seurat)
library(SCpubr)
library(ggplot2)

# Define the file path for the variants dataset
variants_file_path <- "mutations.tsv"

# Define the file path for the single-cell dataset
single_cell_file_path <- "seurat_object.rds"

# Define the user interface
ui <- dashboardPage(
  dashboardHeader(title = "Protein Changes and Single-Cell Analysis"),
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
        h2("Single-Cell Data Visualization"),
        fluidRow(
          column(4,
                 selectInput("reduction", "Select Reduction:",
                             choices = c("umap", "pca"), selected = "umap"),
                 selectInput("group_by", "Group By:",
                             choices = c("Seurat Clusters" = "seurat_clusters", "Original Identity" = "orig.ident"), selected = "seurat_clusters"),
                 checkboxInput("label", "Add Labels", FALSE),
                 checkboxInput("plot_axes", "Show Axes", TRUE),
                 numericInput("legend_ncol", "Legend Columns", 2, min = 1, max = 5),
                 selectInput("assay", "Select Assay:", choices = NULL),
                 selectInput("plot_type", "Plot Type:",
                             choices = c("DimPlot", "NebulosaPlot", "ViolinPlot", "DotPlot", "GeyserPlot"),
                             selected = "DimPlot"),
                 uiOutput("feature_select") # Add UI for feature selection
          ),
          column(8,
                 plotOutput("plot", height = "600px")
          )
        )
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
  
  # Single-cell Logic
  # Load single-cell data
  seurat_object <- reactive({
    req(single_cell_file_path)
    tryCatch({
      sample <- readRDS(single_cell_file_path)
      return(sample)
    }, error = function(e) {
      showNotification(paste("Error loading Seurat object:", e$message), type = "error")
      return(NULL)
    })
  })
  
  # Update assay choices
  observe({
    req(seurat_object())
    sample <- seurat_object()
    assays <- names(sample@assays)
    updateSelectInput(session, "assay", choices = assays, selected = "RNA")
  })
  
  # Create UI for feature selection based on plot type
  output$feature_select <- renderUI({
    req(seurat_object(), input$plot_type, input$assay)
    sample <- seurat_object()
    DefaultAssay(sample) <- input$assay  # Switch to the selected assay
    
    if (input$plot_type %in% c("ViolinPlot", "DotPlot", "GeyserPlot", "NebulosaPlot")) {
      # Access gene names using rownames from the selected assay
      feature_choices <- rownames(sample@assays[[input$assay]])
      
      selectizeInput("features", "Select Features:",
                     choices = feature_choices,
                     multiple = TRUE,  # Allow multiple selections
                     selected = if(input$plot_type == "NebulosaPlot") feature_choices[1] else NULL, # Select first feature as default for Nebulosa
                     options = list(placeholder = "Please select features",
                                    onInitialize = I('function() { this.setValue(""); }')))
    } else {
      return(NULL) # No feature selection for DimPlot
    }
  })
  
  output$plot <- renderPlot({
    req(seurat_object(), input$plot_type, input$assay)
    
    sample <- seurat_object()
    DefaultAssay(sample) <- input$assay  # Switch to the selected assay
    reduction <- input$reduction
    group.by <- input$group_by
    label <- input$label
    plot.axes <- input$plot_axes
    legend.position <- input$legend_ncol
    legend.ncol <- input$legend_ncol
    
    # Use input$features only when it's available and for relevant plot types
    if (!is.null(input$features) && input$plot_type %in% c("ViolinPlot", "DotPlot", "GeyserPlot", "NebulosaPlot")) {
      features <- input$features
    } else {
      features <- VariableFeatures(sample)[1:min(10, length(VariableFeatures(sample)))] # Default features
    }
    
    # Handling different plot types
    if (input$plot_type == "DimPlot") {
      p <- SCpubr::do_DimPlot(sample = sample, reduction = reduction, group.by = group.by,
                              label = label, plot.axes = plot.axes,
                              legend.ncol = legend.ncol)  # Remove legend.position
      print(p)
    } else if (input$plot_type == "NebulosaPlot") {
      req(input$features) #Require some feature to be selected
      
      # Ensure the features are present in the Seurat object
      valid_features <- input$features[input$features %in% rownames(sample)]
      
      if (length(valid_features) > 0) {
        
        # Correctly extract reduction embeddings
        reduction_embeddings <- Embeddings(sample[[input$reduction]])
        
        # Access normalized data for the selected features
        expression_matrix <- GetAssayData(sample, assay = "RNA", slot = "data")[valid_features, , drop = FALSE]
        
        # Create a data frame for ggplot2
        df <- data.frame(reduction_embeddings, t(as.matrix(expression_matrix)))
        colnames(df)[1:2] <- c("x", "y")  # Rename columns for easier reference
        colnames(df)[3] <- "expression" #rename expression of feature for plotting
        
        # Generate the Nebulosa plot
        p <- ggplot(df, aes(x = x, y = y, color = expression)) +  # Use the 'expression' column
          geom_point(size = 1) +
          scale_color_viridis_c(option = "magma") +
          labs(x = paste0(toupper(input$reduction), "1"),
               y = paste0(toupper(input$reduction), "2")) +
          theme_minimal()
        
        # Show or hide axes based on user input
        if (!input$plot_axes) {
          p <- p + theme(axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         panel.border = element_blank(),
                         panel.grid = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank())
        }
        print(p)
        return(p)
        
      } else {
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Selected features not found in Seurat object")
        return(NULL)
      }
      
    } else if (input$plot_type == "ViolinPlot") {
      if (!is.null(features)) {
        p <- SCpubr::do_ViolinPlot(sample = sample, features = features)
        print(p)
      } else {
        top_features <- VariableFeatures(sample)[1:min(10, length(VariableFeatures(sample)))]
        p <- SCpubr::do_ViolinPlot(sample = sample, features = top_features)
        print(p)
      }
    } else if (input$plot_type == "DotPlot") {
      if (!is.null(features)) {
        p <- SCpubr::do_DotPlot(sample = sample, features = features)
        print(p)
      } else {
        top_features <- VariableFeatures(sample)[1:min(10, length(VariableFeatures(sample)))]
        p <- SCpubr::do_DotPlot(sample = sample, features = top_features)
        print(p)
      }
    } else if (input$plot_type == "GeyserPlot") {
      if (!is.null(features) && length(features) > 0) {
        p <- SCpubr::do_GeyserPlot(sample = sample, features = features[1])
        print(p)
      } else {
        top_feature <- VariableFeatures(sample)[1]
        p <- SCpubr::do_GeyserPlot(sample = sample, features = top_feature)
        print(p)
      }
    }
    return(p)
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
