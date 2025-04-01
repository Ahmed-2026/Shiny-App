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

# Define file paths
variants_file_path <- "mutationsF.tsv"
single_cell_file_path <- "seurat_object.rds"

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "Protein Changes"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Variants", tabName = "variants_tab", icon = icon("dna")),
      menuItem("Patient Mutations", tabName = "patient_tab", icon = icon("user")),
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
                 selectizeInput("gene", "Select Gene:", choices = NULL, 
                                options = list(placeholder = 'Select a gene')),
                 selectInput("mutation_type", "Mutation Type:",
                             choices = c("All", "Silent", "Non-silent"),
                             selected = "All"),
                 sliderInput("aa_range", "Amino Acid Position Range:", 
                             min = 0, max = 1000, value = c(0, 1000)),
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
      
      # Patient Mutations Tab
      tabItem(
        tabName = "patient_tab",
        h2("Patient Mutations"),
        fluidRow(
          column(4,
                 selectizeInput("patient", "Select Patient:", choices = NULL,
                                options = list(placeholder = 'Select a patient')),
                 downloadButton("download_patient_data", "Download Patient Data")
          ),
          column(8,
                 tabsetPanel(
                   tabPanel("Mutations Table", DTOutput("patient_mutations_table")),
                   tabPanel("Summary", verbatimTextOutput("patient_summary_stats"))
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
                             choices = c("Seurat Clusters" = "seurat_clusters",
                                         "Original Identity" = "orig.ident")),
                 checkboxInput("label", "Add Labels", FALSE),
                 checkboxInput("plot_axes", "Show Axes", TRUE),
                 numericInput("legend_ncol", "Legend Columns", 2, min = 1, max = 5),
                 selectInput("assay", "Select Assay:", choices = NULL),
                 selectInput("plot_type", "Plot Type:",
                             choices = c("DimPlot", "NebulosaPlot", "ViolinPlot",
                                         "DotPlot", "GeyserPlot")),
                 uiOutput("feature_select")
          ),
          column(8,
                 plotOutput("plot", height = "600px")
          )
        )
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  # Variants Logic
  blacklisted_genes <- c("C1orf87")
  
  variants_data <- reactiveFileReader(
    intervalMillis = 1000,
    session = session,
    filePath = variants_file_path,
    readFunc = function(filePath) {
      read.delim(filePath, sep = "\t", stringsAsFactors = FALSE) %>%
        filter(!Hugo_Symbol %in% blacklisted_genes) %>%
        mutate(
          AA_Position = sapply(Protein_Change, function(x) {
            pos <- as.numeric(gsub(".*?([0-9]+).*", "\\1", x))
            ifelse(is.na(pos), 0, pos)
          }),
          mutation_group = case_when(
            Variant_Classification == "Silent" ~ "Silent",
            TRUE ~ "Non-silent"
          )
        )
    }
  )
  
  observe({
    updateSelectizeInput(session, "gene", 
                         choices = unique(variants_data()$Hugo_Symbol),
                         server = TRUE)
  })
  
  filtered_variants <- reactive({
    req(input$gene)
    df <- variants_data() %>%
      filter(Hugo_Symbol == input$gene) %>%
      {if(input$mutation_type != "All") filter(., mutation_group == input$mutation_type) else .}
    
    validate(need(nrow(df) > 0, "No mutations found for selected filters"))
    df
  })
  
  observe({
    req(filtered_variants())
    updateSliderInput(session, "aa_range",
                      min = min(filtered_variants()$AA_Position),
                      max = max(filtered_variants()$AA_Position))
  })
  
  output$lollipopPlot <- renderG3Lollipop({
    req(filtered_variants())
    g3Lollipop(
      filtered_variants(),
      gene.symbol = input$gene,
      gene.symbol.col = "Hugo_Symbol",
      aa.pos.col = "AA_Position",
      protein.change.col = "Protein_Change",
      factor.col = "mutation_group",
      plot.options = g3Lollipop.theme(theme.name = "cbioportal")
    )
  })
  
  output$data_table <- renderDT({
    datatable(filtered_variants(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$summary_stats <- renderPrint({
    cat("Total mutations:", nrow(filtered_variants()), "\n")
    cat("Unique protein changes:", length(unique(filtered_variants()$Protein_Change)), "\n")
    cat("\nMutation types:\n")
    print(table(filtered_variants()$mutation_group))
  })
  
  # Patient Mutations Logic
  observe({
    updateSelectizeInput(session, "patient",
                         choices = unique(variants_data()$Tumor_Sample_Barcode),
                         server = TRUE)
  })
  
  patient_data <- reactive({
    req(input$patient)
    variants_data() %>% filter(Tumor_Sample_Barcode == input$patient)
  })
  
  output$patient_mutations_table <- renderDT({
    datatable(patient_data(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$patient_summary_stats <- renderPrint({
    cat("Patient:", input$patient, "\n")
    cat("Total mutations:", nrow(patient_data()), "\n")
    cat("Unique genes:", length(unique(patient_data()$Hugo_Symbol)), "\n")
    cat("\nMutation types:\n")
    print(table(patient_data()$Variant_Classification))
  })
  
  # Single-cell Logic
  seurat_obj <- reactive({
    req(file.exists(single_cell_file_path))
    readRDS(single_cell_file_path) %>% UpdateSeuratObject()
  })
  
  observe({
    req(seurat_obj())
    updateSelectInput(session, "assay", 
                      choices = names(seurat_obj()@assays),
                      selected = DefaultAssay(seurat_obj()))
  })
  
  output$feature_select <- renderUI({
    req(seurat_obj(), input$plot_type, input$assay)
    if(input$plot_type != "DimPlot") {
      features <- rownames(seurat_obj()[[input$assay]])
      selectizeInput("features", "Select Features:", 
                     choices = features, 
                     multiple = TRUE,
                     options = list(placeholder = 'Select features'))
    }
  })
  
  output$plot <- renderPlot({
    req(seurat_obj(), input$plot_type)
    
    base_params <- list(
      sample = seurat_obj(),
      reduction = input$reduction,
      group.by = input$group_by,
      label = input$label,
      legend.ncol = input$legend_ncol
    )
    
    if(input$plot_type == "DimPlot") {
      do.call(do_DimPlot, base_params)
    }
    else if(input$plot_type == "NebulosaPlot") {
      req(input$features)
      p <- SCpubr::do_NebulosaPlot(
        sample = seurat_obj(),
        features = input$features,
        reduction = input$reduction,
        joint = length(input$features) > 1
      )
      
      if(!input$plot_axes) {
        p <- p & theme(
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()
        )
      }
      print(p)
    }
    else if(input$plot_type == "ViolinPlot") {
      req(input$features)
      do.call(do_ViolinPlot, c(base_params, list(features = input$features)))
    }
    else if(input$plot_type == "DotPlot") {
      req(input$features)
      do.call(do_DotPlot, c(base_params, list(features = input$features)))
    }
    else if(input$plot_type == "GeyserPlot") {
      req(input$features)
      do.call(do_GeyserPlot, c(base_params, list(features = input$features[1])))
    }
  })
}

# Run application
shinyApp(ui, server)
