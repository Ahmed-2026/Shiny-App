# Load necessary packages
library(shiny)
library(g3viz)
library(dplyr)
library(DT)
library(shinydashboard)
library(shinyjs)
library(ggplot2)

# Define file path
variants_file_path <- "mutationsF.tsv"

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "Protein Changes"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Variants", tabName = "variants_tab", icon = icon("dna")),
      menuItem("Patient Mutations", tabName = "patient_tab", icon = icon("user"))
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
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  blacklisted_genes <- c("C1orf87")
  
  variants_data <- reactiveFileReader(
    intervalMillis = 2000,
    session = session,
    filePath = variants_file_path,
    readFunc = function(filePath) {
      # Skip 0 lines if there are no metadata lines; adjust skip if needed!
      df <- read.delim(filePath, sep = "\t", skip = 0, stringsAsFactors = FALSE)
      validate(need("Hugo_Symbol" %in% colnames(df), "CRITICAL ERROR: Missing Hugo_Symbol column"))
      df %>%
        filter(!Hugo_Symbol %in% blacklisted_genes) %>%
        mutate(
          AA_Position = sapply(Protein_Change, function(x) {
            if(is.na(x)) return(0)
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
  
  # Gene selector: unique, in file order!
  observe({
    genes <- variants_data()$Hugo_Symbol
    genes <- genes[!duplicated(genes)]  # Keep only unique, in file order
    updateSelectizeInput(session, "gene", 
                         choices = genes, 
                         selected = genes[1],
                         server = TRUE)
  })
  
  # Variants filter
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
    patients <- variants_data()$Tumor_Sample_Barcode
    patients <- patients[!is.na(patients) & 
                           patients != "" & 
                           !grepl("^[><=]", patients)]
    updateSelectizeInput(session, "patient",
                         choices = sort(unique(patients)),
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
}

# Run application
shinyApp(ui, server)
