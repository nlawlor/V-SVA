options(repos = BiocInstaller::biocinstallRepos())
getOption("repos")
library(shiny)
library(rsconnect)
library(iasva)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(DelayedArray)
library(irlba)
library(Rtsne)
library(RColorBrewer)
library(pheatmap)
library(corrplot)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(clusterProfiler)
library(udunits2)
library(DOSE)
library(edgeR)
library(scran)
library(preprocessCore)
library(corrplot)
library(ggplot2)
library(plotly)
library(DT)
library(shinyjs)
library(shinythemes)
library(rintrojs)
library(htmlwidgets)

################################################################################################################
# Define UI
################################################################################################################
ui <- shinyUI(fluidPage(theme = shinytheme("cerulean"),
                        # Application title
                        titlePanel(div("Welcome to the Shiny App for Iteratively Adjusted Surrogate Variable Analysis (IA-SVA)", 
                                       img(height = 88, width = 217, 
                                           src = "jax.logo.gif", 
                                           class = "pull-right"))),
                        # side bar where user first interacts 
                        sidebarPanel(
                          # call function to use package rintrojs                                
                          introjsUI(),
                          # add first button for guided tutorial
                          introBox(
                            actionButton("help", "Click here for a guided tutorial of the app!", icon = icon("list-ol")),
                            data.step = 1,
                            data.intro = "This is the start of the tutorial. Click 'Next Step' to continue."
                          ),
                          
                          # Or provide download button for manual
                          introBox(
                            downloadButton("Download_Manual", "Download Tutorial Document"),
                            data.step = 2,
                            data.intro = "Click this button to download an in-depth tutorial/manual for this app."
                          ),
                          
                          # download test datasets
                          introBox(
                            downloadButton("Download_Test_Exp", "Download Test Expression Data"),
                            data.step = 3,
                            data.intro = "Click this button to download an example gene expression file. This file may be used to test out the app."
                          ),
                          
                          introBox(
                            downloadButton("Download_Test_Meta", "Download Test Metadata"),
                            data.step = 4,
                            data.intro = "Click this button to download an example sample metadata file. This file may be used to test out the app."
                          ),
                      
                          # Load user data
                          introBox(
                            h4("1. Load Data"),
                            fileInput(
                              inputId = "input_exprs",
                              "Choose gene/ADT file to upload",
                              accept = c(
                                ".Rds",
                                ".csv",
                                ".txt"
                              )
                            ),
                            fileInput(
                              inputId = "input_meta",
                              "Choose sample metadata file to upload",
                              accept = c(
                                ".Rds",
                                ".csv",
                                ".txt"
                              )
                            ),
                            data.step = 5, data.intro = "Click the Browse buttons to provide your gene/ADT expression data and sample metadata.
                            Gene expression data should be in the format of: gene symbol names in rows, sample names in columns. Sample meta data should
                            have sample names in rows and meta data variables in columns.",
                            data.hint = "The following input file formats are accepted: .Rds, tab-delimited text, .csv"),
                          introBox(
                            h4("2. Data Preprocessing"),
                            div(style="display: inline-block;vertical-align:top; width: 100px;",
                                numericInput(inputId = "Cell_num", label = "# of Cells",
                                            min = 0, value = 5)),
                            div(style="display: inline-block;vertical-align:top; width: 100px;",
                                numericInput(inputId = "Count_num", label = "# of Counts",
                                            min = 0, value = 5)),
                            div(style="display: inline-block;vertical-align:top; width: 200px;",
                                selectInput(inputId = "norm_method", label = "Normalization Method",
                                            choices = c("CPM", "Quantile", "scran", "None"), selected = "CPM")),
                            data.step = 6, data.intro = "Specify the criteria to filter genes. E.g., genes with 3 or more counts in at least 5 cells. From the drop-down menu, select your desired normalization method.",
                            data.hint = "CPM = Counts per million from edgeR package, Quantile = quantile normalization from preprcessCore R package,
                            scran = method to deconvolute size factors from cell pools from scran R package."
                            ),
                          introBox(
                            # after the user loads the data, update the select input options to reflect metadata table column names
                            h4("3. Specify Known Factors to Adjust For"),
                            div(style="display: inline-block;vertical-align:top; width: 200px;",
                                selectInput(inputId = "known_factors", label = "Known Factors",
                                            choices = c(""), selected = NULL, multiple = TRUE)),
                            data.step = 7, data.intro = "From the drop-down menu, select the known factor(s) you would like to adjust for. You may choose a single or multiple factors",
                            data.hint = "If you have not uploaded a metadata file, no choices will be available yet."
                          ),
                          introBox(
                            # after the user loads the data, update the select input options to reflect metadata table column names
                            h4("4. IA-SVA Analysis"),
                            # options to customize parameters
                            div(style="display: inline-block;vertical-align:top; width: 200px;",
                                numericInput(inputId = "pct_cutt", label = "% Threshold for SV retention", value = 1,
                                             min = 1, max = 99)),
                            br(),
                            shiny::actionButton("do", "Run Analysis", icon = icon("paper-plane"), class = "btn-primary"),
                            data.step = 8, data.intro = "Click this button to run IA-SVA and identify unknown sources (surrogate variables) of variation within your data",
                            data.hint = "Percentage threshold for SV retention. IA-SVA computes the percentage of unmodeled variance explained by the putative hidden factor and compare it with the user-defined threshold. If the percentage is greater than the threshold, SV is retained."
                          ),
                          
                          br(),
                          h4("Author: "),
                          introBox(
                            tags$div(class = "header", checked = NA,  
                                     tags$i("Nathan Lawlor (nathan.lawlor03@gmail.com)"),
                                     tags$p(""),
                                     tags$a(href = "https://www.jax.org/research-and-faculty/research-labs/the-ucar-lab",
                                            "Visit the Ucar Lab Here!", target = "_blank")),
                            data.step = 9, data.intro = "Please contact me with any questions about the app. I'm happy to help! Also, checkout the resources below!"),
                          h4("Other Resources: "),
                          tags$div(class = "header", checked = NA,  
                                   tags$i(""),
                                   tags$p(""),
                                   tags$a(href = "https://www.bioconductor.org/packages/devel/bioc/html/iasva.html",
                                          "Check out IASVA on Bioconductor", target = "_blank")),
                          tags$div(class = "header", checked = NA,  
                                   tags$i(""),
                                   tags$p(""),
                                   tags$a(href = "https://www.biorxiv.org/content/early/2018/04/24/151217",
                                          "Check out the IASVA manuscript on bioRxiv", target = "_blank"))
                        ),
                        
                        # main panel with multiple tabs: detect SVs, find marker genes, annotate genes, visualize data with tsne
                        introBox(
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Data/QC",
                                       verbatimTextOutput("load_data"),
                                       br(),
                                       # add two plots below this text output with qc stats
                                       fluidRow(column(width = 6,
                                                       plotOutput("Detect", height = 500, width = 500)
                                        ), 
                                        column(width = 6, plotOutput("Total_Genes", height = 500, width = 500)
                                        )
                                      )
                              ),
                              tabPanel("Surrogate Variables", downloadButton("Download_SVs", "Download All Surrogate Variables Plot"),
                                       downloadButton("Download_Correlation", "Download Correlation Plot"),
                                       fluidRow(
                                         tags$div(
                                                tags$h5("Graphing Options")),
                                                # specify how to color points in all sv plot
                                                div(style="display: inline-block;vertical-align:left; width: 200px;",
                                                    selectInput(inputId = "All_SV_Color", label = "Color Points",
                                                                choices = c(""), selected = NULL, multiple = FALSE)),
                                                div(style="display: inline-block;vertical-align:left; width: 200px;",
                                                    selectInput(inputId = "All_SV_Num", label = "SV's to Plot",
                                                                choices = c(""), selected = NULL, multiple = TRUE)),
                                                br(),
                                                shiny::actionButton("updateAllSV", "Update Plots", icon = icon("refresh")),
                                                align="left"),
                                         fluidRow(column(width = 6,
                                           plotOutput("corrplot", height = 500, width = 500)
                                        ), column(width = 6,
                                           plotOutput("SVplot", height = 500, width = 500)
                                       )
                                       )),
                              tabPanel("Pairwise Surrogate Variable", 
                                       downloadButton("Download_Pairwise", "Download Pairwise Plot"),
                                       fluidRow(
                                         column(width = 3,
                                                tags$div(
                                                  tags$h5("")),
                                                # specify how to color points in pairwise plot, which svs to visualize
                                                div(style="display: inline-block;vertical-align:top; width: 200px;",
                                                    selectInput(inputId = "Pair_SV_Color", label = "Color Points",
                                                                choices = c(""), selected = NULL, multiple = FALSE)),
                                                div(style="display: inline-block;vertical-align:top; width: 100px;",
                                                    selectInput(inputId = "sv_x", label = "SV X-axis",
                                                                choices = c(""), selected = NULL, multiple = FALSE)),
                                                div(style="display: inline-block;vertical-align:top; width: 100px;",
                                                    selectInput(inputId = "sv_y", label = "SV Y-axis",
                                                                choices = c(""), selected = NULL, multiple = FALSE)),
                                                shiny::actionButton("updatePairSV", "Update Pairwise SV Plot", icon = icon("refresh")),
                                                align="left"),
                                         column(width = 9, plotlyOutput("PairSV", height = 600, width = 600)))),
                              tabPanel("Identify Marker Genes", downloadButton("Download_Heatmap", "Download Heatmap"),
                                       br(),
                                       # Select svs to choose for marker gene identification
                                       selectInput(inputId = "SV_marks", label = "Choose SVs", 
                                                   choices = "", multiple = TRUE, selected = NULL),
                                       
                                       # specify mult testing correction
                                       selectInput(inputId = "mark_sig", label = "P-Value Adjustment", 
                                                   choices = c("BH", "bonferroni", "none"),
                                                   selected = "BH"),
                                       # specify how to color points in all sv plot
                                       div(style="display: inline-block;vertical-align:left; width: 200px;",
                                       # specify sig cutof
                                       numericInput(inputId = "mark_cutoff", label = "Adjusted P-value cutoff", 
                                                   min = 0, max = 1, value = 0.05, step = 0.05)
                                       ),
                                       div(style="display: inline-block;vertical-align:left; width: 200px;",
                                       # specify r2 cutoff
                                       numericInput(inputId = "rsqcutoff", label = "R-squared cutoff",
                                                    min = 0, max = 1, value = 0.3, step = 0.1)
                                       ),
                                       br(),
                                       # identify known factors to include in heatmap
                                       selectInput(inputId = "heatmap_known_factors", label = "Known Factors to Include in Heatmap",
                                                   choices = c(""), selected = NULL, multiple = TRUE),
                                       br(),
                                       shiny::actionButton("FindGenes", "Identify Marker Genes", icon = icon("search"), class = "btn-primary"),
                                       plotOutput("MarkerHeatmap", height = 800, width = 800)),
                              tabPanel("Marker Genes Table", downloadButton("Download_Markers", "Download Marker Genes Table"), DT::dataTableOutput("genes_table")),
                              tabPanel("Gene Enrichment Analysis",
                                       downloadButton("Download_path_plot", "Download Enrichment Analysis Plot"),
                                       # Select svs to choose for marker gene identification
                                       selectInput(inputId = "Path_Type", label = "Enrichment Analysis Type", 
                                                   choices = c("Gene Ontology Biological Process",
                                                               "Gene Ontology Cellular Component",
                                                               "Gene Ontology Molecular Function",
                                                               "KEGG", "Homo sapiens Immune Modules",
                                                               "Homo sapiens PBMC Cell Specific Modules",
                                                               "Custom"),
                                                   multiple = FALSE,
                                                   selected = "Gene Ontology Biological Process"),
                                       
                                       # add conditional panel for custom gene lists
                                       conditionalPanel(
                                         condition = "input.Path_Type == 'Custom'",
                                         fileInput(
                                           inputId = "input_gene_mod",
                                           "Upload custom gene and module file",
                                           accept = ".txt"
                                          )
                                       ),
                                       conditionalPanel(
                                         condition = "input.Path_Type == 'Custom'",
                                         fileInput(
                                           inputId = "input_mod_names",
                                           "Upload custom module identifier file",
                                           accept = ".txt"
                                         )
                                       ),
                                       # indicate species type
                                       selectInput(inputId = "Species_Type", label = "Species", 
                                                   choices = c("Homo sapiens", "Mus musculus",
                                                               "Rattus norvegicus"), multiple = FALSE,
                                                   selected = "Homo sapiens"),
                                       div(style="display: inline-block;vertical-align:left; width: 200px;",
                                       selectInput(inputId = "pvalue_correct", label = "P-value Adjustment",
                                                   choices = c("BH", "bonferroni", "none"),
                                                   selected = "BH")
                                       ),
                                       # specify sig cutof
                                       div(style="display: inline-block;vertical-align:left; width: 200px;",
                                       numericInput(inputId = "path_cutoff", label = "Adjusted P-value cutoff", 
                                                    min = 0, max = 1, value = 0.2, step = 0.05)
                                       ),
                                       br(),
                                       # indicate max num of results to display
                                       div(style="display: inline-block;vertical-align:left; width: 200px;",
                                           numericInput(inputId = "path_viz_num", label = "Max # of Results to Visualize", 
                                                        min = 0, value = 10, step = 1)
                                       ),
                                       br(),
                                       shiny::actionButton("Path_Analysis", "Run Enrichment Analysis", icon = icon("paper-plane"), class = "btn-primary"),
                                       plotOutput("Enrich_Plot", height = 800, width = 800)
                                       ),
                              tabPanel("Gene Enrichment Analysis Table", 
                                       downloadButton("Download_path_table", "Download Enrichment Analysis Results Table"), DT::dataTableOutput("enrich_table")),
                              tabPanel("Dimension Reduction and Visualization", downloadButton("download_dim_all", "Download All Genes Plot"),
                                       downloadButton("download_dim_marker", "Download Marker Genes Plot"),
                                       # choose dimension reduction method
                                       selectInput(inputId = "Dim_Type", label = "Dimension Reduction", 
                                                   choices = c("PCA", "t-SNE"), multiple = FALSE,
                                                   selected = "PCA"),
                                       # choose variable to color points
                                       selectInput(inputId = "Dim_Color", label = "Color Points", 
                                                   choices = "", multiple = FALSE,
                                                   selected = NULL),
                                       # this button compute the pca/tsne
                                       shiny::actionButton("Dim_Analysis", "Run Dimension Reduction Analysis", icon = icon("paper-plane"), class = "btn-primary"),
                                       # this button can be used to just change plot coloring (after computing tsne/pca)
                                       shiny::actionButton("Update_Dim", "Update Plot Point Colors", icon = icon("refresh")),
                                       br(),
                                       br(),
                                       fluidRow(
                                         column(width = 6, align = "left",
                                                plotlyOutput("Dim_Plot_Orig", height = 600, width = 600)
                                                ),
                                         column(width = 6, align = "right",
                                                plotlyOutput("Dim_Plot_Markers", height = 600, width = 600))
                                      )
                                    )
                                  )
                                )
                              )
  )
)



################################################################################################################
# Define server logic required to plot data
################################################################################################################
server <- shinyServer(function(input, output, session) {
  set.seed(1)
  # file upload size limit (set to 2000 MB)
  options(shiny.maxRequestSize=2000*1024^2) 
  # color function
  color.vec <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
                 "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "black")
  
  # reactive values to store data
  dataTables <- reactiveValues(
    exp_mat = NULL,
    exp_filt = NULL,
    meta_df = NULL,
    exp_norm = NULL,
    iasva.res = NULL,
    summ_exp = NULL,
    markers = NULL,
    iasva_vars = NULL,
    markers_formatted = NULL,
    gene.df = NULL,
    species = NULL,
    species_kegg = NULL,
    univ_gene = NULL,
    enrich_res = NULL,
    pre_dim_mark = NULL,
    pre_dim_orig = NULL,
    dim_mark = NULL,
    dim_orig = NULL
  )
  
  # load in expression data, check file extension
  observeEvent(input$input_exprs, {
    if (grepl(x = input$input_exprs$datapath, pattern = ".csv")) {
      dataTables$exp_mat <- read.csv(file = input$input_exprs$datapath,
                                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    } else if (grepl(x = input$input_exprs$datapath, pattern = ".txt")) {
      dataTables$exp_mat <- read.delim(file = input$input_exprs$datapath,
                                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    } else if (grepl(x = input$input_exprs$datapath, pattern = ".Rds")) {
      dataTables$exp_mat <- readRDS(file = input$input_exprs$datapath)  
    }
    # change file to matrix
    dataTables$exp_mat <- as.matrix(dataTables$exp_mat)
  })
    
  observeEvent(input$input_meta, {
    # load in metadata
    if (grepl(x = input$input_meta$datapath, pattern = ".csv")) {
      dataTables$meta_df <- read.csv(file = input$input_meta$datapath,
                                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    } else if (grepl(x = input$input_meta$datapath, pattern = ".txt")) {
      dataTables$meta_df <- read.delim(file = input$input_meta$datapath,
                                       header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    } else if (grepl(x = input$input_meta$datapath, pattern = ".Rds")) {
      dataTables$meta_df <- readRDS(file = input$input_meta$datapath)
    }
    dataTables$meta_df <- as.data.frame(dataTables$meta_df)
    # change all variables to factors
    for (fac in 1:ncol(dataTables$meta_df)) {
      dataTables$meta_df[, fac] <- as.factor(dataTables$meta_df[, fac])
    }
    
    # identify variables in meta data 
    updateSelectInput(session = session, inputId = "known_factors", label = "Known Factors",
                      choices = colnames(dataTables$meta_df))
    updateSelectInput(session = session, inputId = "heatmap_known_factors", label = "Known Factors to Include in Heatmap",
                      choices = colnames(dataTables$meta_df))
    
  })
  
  
  # shiny app error validation for loading data, and plots, tables
    # validate functions must go inside the plot/table rendering functions
      # option could be to have a data loading panel in which we provide qc and make sure the data is loaded properly
      # this first panel would also conduct all iasva analyses
  
  valid_load <- function() {
    shiny::validate(
      # adding an isolate() will make the error messages only appear/disappear when click run analysis button
      # before can run ia-sva need to have loaded data
      need(isolate(dataTables$exp_mat != ""), "Please load a gene/ADT file before running IA-SVA"),
      need(isolate(dataTables$meta_df != ""), "Please load a metadata file before running IA-SVA"),
      # make sure specify known variables
      need(isolate(input$known_factors != ""), "Please specify known factor(s) to adjust for before running IA-SVA"),
      # maybe have option to specify no variables?
      need(isolate(input$pct_cutt >= 1), "Percent threshold must be greater than or equal to 1"),
      need(isolate(input$pct_cutt < 100), "Percent threshold must be less than 100"),
      # cell filtering
      need(isolate(input$Cell_num >= 0), "Number of cells must be greather than or equal to 0"),
      need(isolate(input$Count_num >= 0), "Number of gene/ADT counts must be greather than or equal to 0"),
      # check for numeric inputs
      need(isolate(is.numeric(input$pct_cutt)), "Percent threshold must be numeric"),
      need(isolate(is.numeric(input$Cell_num)), "Number of cells must be numeric"),
      need(isolate(is.numeric(input$Count_num)), "Number of gene/ADT counts must be numeric")
      
    )
  }
  
  valid_iasva <- function() {
    shiny::validate(
      need(isolate(!is.null(dataTables$iasva.res)), "No surrogate variables obtained... Please adjust input parameters")
    )
  }
  
  # validation for surrogate variable pair plots
  valid_pair_plot <- function() {
    shiny::validate(
      need(isolate(length(as.numeric(input$All_SV_Num)) >= 2), "Error: Please select 2 or more SVs in order to visualize this plot")
    )
  }
  
  # validation for marker genes
  valid_markers <- function() {
    shiny::validate(
      need(isolate(length(as.numeric(input$SV_marks)) >= 1), "Error: Please select 1 or more SVs to identify markers for"),
      need(isolate(input$heatmap_known_factors != ""), "Please specify known factor(s) to include in heatmap visualization"),
      need(isolate(input$mark_cutoff > 0), "Adjusted p-value cutoff must be greater than 0"),
      need(isolate(input$mark_cutoff <= 1), "Adjusted p-value cutoff must be less than or equal to 1"),
      need(isolate(input$rsqcutoff > 0), "R-squared cutoff must be greater than 0"),
      need(isolate(input$rsqcutoff <= 1), "R-squared cutoff must be less than or equal to 1")
    )
  }
  
  # validation for enrichment analysis
  valid_enrich <- function() {
    shiny::validate(
      need(isolate(length(as.numeric(input$SV_marks)) >= 1), "Error: Please select 1 or more SVs to identify markers for"),
      need(isolate(!is.null(dataTables$markers)), "Error: No marker genes were identified... Please adjust r-squared and significance cutoffs and re-calculate before proceeding with enrichment analysis"),
      need(isolate(input$path_cutoff > 0), "Enrichment analysis adjusted p-value cutoff must be greater than 0"),
      need(isolate(input$path_cutoff <= 1), "Enrichment analysis adjusted p-value cutoff must be less than or equal to 1")
    )
  }
  
  # shiny app error validation for dim plots
  valid_dim <- function() {
    shiny::validate(
      need(isolate(!is.null(dataTables$markers)), "Error: No marker genes were identified... Please adjust r-squared and significance cutoffs and re-calculate before proceeding with dimension reduction analysis")
    )
  }
  
  
  # extract known factors to adjust for, create model matrix
  observeEvent(input$do, {
    output$load_data <- renderText({
      # error handling
      valid_load()
      # by filter genes with low counts
      withProgress(expr = filter <- apply(dataTables$exp_mat, 1, function(x) length(x[x>isolate(input$Count_num)])>=isolate(input$Cell_num)),
                   message = "Preprocessing data, please wait")
      dataTables$exp_filt <- dataTables$exp_mat[filter,]
      
      # normalize the data
      if (isolate(input$norm_method) == "CPM") {
        dataTables$exp_norm <- edgeR::cpm(dataTables$exp_filt)
      } else if (isolate(input$norm_method) == "Quantile") {
        dataTables$exp_norm <- normalize.quantiles(dataTables$exp_filt)
      } else if (isolate(input$norm_method) == "scran") {
        sce <- SingleCellExperiment(list(counts=dataTables$exp_filt))
        sce <- computeSumFactors(sce)
        sce <- normalize(sce)
        dataTables$exp_norm <- exprs(sce)
      } else if (isolate(input$norm_method) == "None") {
        dataTables$exp_norm <- dataTables$exp_filt
      }
      
      # ensure that the gene and sample matrices are ordered the same
      if (all(table(colnames(dataTables$exp_norm) == rownames(dataTables$meta_df))) == FALSE) {
        id_samp <- NULL
        withProgress(expr =
                       for (s_num in 1:nrow(dataTables$meta_df)) {
                         id_s <- which(colnames(dataTables$exp_norm) == rownames(dataTables$meta_df)[s_num])
                         id_samp <- c(id_samp, id_s)
                       }
                     , message = "Gene/ADT and sample matrices have different name orders, re-ordering gene/ADT matrix now...")
        dataTables$exp_norm <- dataTables$exp_norm[, id_samp]
      } else {}
      
      # extract known factors
      id_mod <- which(colnames(dataTables$meta_df) %in% isolate(input$known_factors))
      if (length(id_mod) > 1) {
        formdf1 <- as.formula(paste("~", colnames(dataTables$meta_df)[id_mod][1], "+", paste(colnames(dataTables$meta_df)[id_mod[2:length(id_mod)]],collapse="+"), sep = ""))
        mod <- model.matrix(formdf1, data = dataTables$meta_df)
      } else {
        varf1 <- as.factor(dataTables$meta_df[, id_mod])
        mod <- model.matrix(~varf1, data = dataTables$meta_df)
      }
      # IA-SVA analysis
      summ_exp <- SummarizedExperiment(assays = as.matrix(dataTables$exp_norm))
      dataTables$summ_exp <- summ_exp
      dataTables$iasva.res <- withProgress(expr = fast_iasva(summ_exp, mod[,-1, drop = F], verbose=FALSE,
                                                  pct.cutoff = isolate(input$pct_cutt), num.sv = NULL),
                                message = "IA-SVA analysis in progress, please wait")
      # if no SV's are calculated inform user
      valid_iasva()
      
      # display SV's
      iasva.sv <- as.data.frame(dataTables$iasva.res$sv)
      rownames(iasva.sv) <- colnames(dataTables$exp_norm)
      
      # update coloring scheme for total plot
      updateSelectInput(session = session, inputId = "All_SV_Color", label = "Color Points",
                        choices = colnames(dataTables$meta_df), selected = colnames(dataTables$meta_df[1]))
      # update number of svs to visualize
      updateSelectInput(session = session, inputId = "All_SV_Num", label = "SV Number",
                        choices = 1:ncol(iasva.sv), selected = 1:ncol(iasva.sv))
      
      # update coloring scheme for pairwise
      updateSelectInput(session = session, inputId = "Pair_SV_Color", label = "Color Points",
                        choices = colnames(dataTables$meta_df), selected = colnames(dataTables$meta_df)[1])
      # update coloring scheme for dimension reduction plots
      updateSelectInput(session = session, inputId = "Dim_Color", label = "Color Points",
                        choices = colnames(dataTables$meta_df), selected = colnames(dataTables$meta_df)[1])
      # update sv x/y selections
      updateSelectInput(session = session, inputId = "sv_x", label = "SV X-axis",
                        choices = colnames(iasva.sv), selected = "SV1")
      updateSelectInput(session = session, inputId = "sv_y", label = "SV Y-axis",
                        choices = colnames(iasva.sv), selected = "SV2")
      # update sv selections in marker panel
      updateSelectInput(session = session, inputId = "SV_marks", label = "Choose SVs",
                        choices = colnames(iasva.sv), selected = "SV1")
      
      # output messages
      print(paste("Data loaded and IA-SVA analysis finished! \n",
                  "Gene number: ", nrow(dataTables$exp_norm), " \n", 
                  "Sample number: ", ncol(dataTables$exp_norm), " \n",
                  "Number of SV's identified: ", ncol(dataTables$iasva.res$sv), sep = ""))
      
    })
    
    # make a histogram of detected genes
    output$Detect <- renderPlot({
      valid_load()
      valid_iasva()
      # binarize data
      bin_data <- dataTables$exp_norm
      bin_data[bin_data < 1] <- 0
      bin_data[bin_data >= 1] <- 1
      num.exp <- apply(bin_data,2,sum)
      summ <- summary(num.exp)
      hist(num.exp, col = "dodgerblue", main="Features detected in each sample", 
           ylab = "Samples (n)", xlab = "Number of features detected")
      legend("topright", legend = paste(names(summ), round(summ, digits = 2), sep = " "), title = "Summary of Features Detected")
    })
    
    # make a violin plot of total log2 scaled counts
    output$Total_Genes <- renderPlot({
      valid_load()
      valid_iasva()
      csum <- colSums(dataTables$exp_filt)
      boxplot(csum, main = "Total features in each sample (total raw counts after preprocessing)",
              ylab = "Total Raw Counts", xlab = "")
      csum_summ <- summary(csum)
      legend("topright", legend = paste(names(csum_summ), round(csum_summ, digits = 2), sep = " "), title = "Summary of Feature Totals")
    })
    
    # make a grid of all svs by default and no coloring
    output$SVplot <- renderPlot({
      valid_load()
      valid_iasva()
      par(mar = c(5, 5, 4, 2), xpd = TRUE)
      iasva.sv <- as.data.frame(dataTables$iasva.res$sv)
      rownames(iasva.sv) <- colnames(dataTables$exp_norm)
      pairs(iasva.sv, main="", pch=20, cex=0.5, lower.panel = NULL)
    })
    
    # make a corrplot of the svs and known factors (that they choose to adjust for?)
    output$corrplot <- renderPlot({
      valid_load()
      valid_iasva()
      # change factors to numeric for correlation
      meta_sel <- dataTables$meta_df
      for (jcol in 1:ncol(meta_sel)) {
        meta_sel[,jcol] <- as.numeric(as.factor(meta_sel[,jcol]))
      }
      iasva_vars <- cbind(dataTables$iasva.res$sv, meta_sel)
      # need to append column names to matrix
      colnames(iasva_vars) <- c(paste("SV", 1:ncol(dataTables$iasva.res$sv), sep = ""),
                                colnames(dataTables$meta_df))
      dataTables$iasva_vars <- iasva_vars
      par(mar = c(5, 5, 4, 2), xpd = TRUE)
      corrplot(cor(iasva_vars), type = "upper")
    })
    
    # make interactive plots
    output$PairSV <- renderPlotly({
      valid_load()
      valid_iasva()
      iasva.sv <- as.data.frame(dataTables$iasva.res$sv)
      rownames(iasva.sv) <- colnames(dataTables$exp_norm)
      plot_ly(iasva.sv, x = ~SV1, y = ~SV2, type = "scatter",
              mode = "markers", text = paste("Cell ID: ", rownames(iasva.sv), sep = ""))
    })
    
  })
  
  # if user specifies an option for pairwise SV plots, then update
    # also update the correlation plot
    # code gives error, only one column in argument to pairs
  observeEvent(input$updateAllSV, {
    output$SVplot <- renderPlot({
      valid_load()
      valid_iasva()
      valid_pair_plot()
      isolate({
        # paired SV
        par(mar = c(5, 5, 4, 2), xpd = TRUE)
        id_fac <- which(colnames(dataTables$meta_df) == input$All_SV_Color)
        fac_int <- as.factor(dataTables$meta_df[, id_fac])
        pairs(dataTables$iasva.res$sv[, c(as.numeric(input$All_SV_Num))], main="", pch=20, cex=0.5,
              col = color.vec[fac_int], bg = color.vec[fac_int], lower.panel = NULL)
        legend("bottomleft", levels(fac_int), fill=color.vec,
               title = as.character(colnames(dataTables$meta_df)[id_fac]))
      })
    })
    output$corrplot <- renderPlot({
      valid_load()
      valid_iasva()
      valid_pair_plot()
      # correlation plot
      isolate({
        # change factors to numeric for correlation
        meta_sel <- dataTables$meta_df
        for (jcol in 1:ncol(meta_sel)) {
          meta_sel[,jcol] <- as.numeric(as.factor(meta_sel[,jcol]))
        }
        iasva_vars <- cbind(dataTables$iasva.res$sv[, c(as.numeric(input$All_SV_Num))], meta_sel)
        # need to append column names to matrix
        colnames(iasva_vars) <- c(paste("SV", c(as.numeric(input$All_SV_Num)), sep = ""),
                                  colnames(dataTables$meta_df))
        corrplot(cor(iasva_vars), type = "upper")
      })
    })
  })
  
  # generate updated plotly plot
  observeEvent(input$updatePairSV, {
    output$PairSV <- renderPlotly({
      valid_load()
      valid_iasva()
      isolate ({
        # grouping color variable
        id_fac <- which(colnames(dataTables$meta_df) == input$Pair_SV_Color)
        fac_int <- as.factor(dataTables$meta_df[, id_fac])
        # sv axis labels
        x_ax <- list(
          title = as.character(input$sv_x)
        )
        y_ax <- list(
          title = as.character(input$sv_y)
        )
        # extract chosen SVs
        id_svx <- which(colnames(dataTables$iasva.res$sv) == input$sv_x)
        id_svy <- which(colnames(dataTables$iasva.res$sv) == input$sv_y)
        
        plot_ly(as.data.frame(dataTables$iasva.res$sv), x = dataTables$iasva.res$sv[,id_svx],
                y = dataTables$iasva.res$sv[,id_svy], type = "scatter",
                mode = "markers", text = paste("Cell ID: ", rownames(dataTables$meta_df), sep = ""),
                color = ~fac_int) %>%
          layout(xaxis = x_ax, yaxis = y_ax)
      })
    })
  })
  
  observeEvent(input$FindGenes, {
    output$MarkerHeatmap <- renderPlot({
      valid_load()
      valid_iasva()
      valid_markers()
      # identify marker genes (plot all by default)
      # identify which svs were chosen
      id_sv_mark <- which(colnames(dataTables$iasva.res$sv) %in% isolate(input$SV_marks))
      marker_genes <- withProgress(expr = iasva::find_markers(Y = dataTables$summ_exp,
                                                              iasva.sv = as.matrix(dataTables$iasva.res$sv[, id_sv_mark, drop=FALSE]),
                                                              rsq.cutoff = isolate(input$rsqcutoff), method = isolate(input$mark_sig), sig.cutoff = isolate(input$mark_cutoff)),
                                   message = paste("Identifying marker genes for ", paste(colnames(dataTables$iasva.res$sv)[id_sv_mark], collapse = ","), " please wait", sep = ""))
      dataTables$markers <- marker_genes
      # specify which metadata to plot on heatmap
      id_mod <- which(colnames(dataTables$meta_df) %in% isolate(input$heatmap_known_factors))
      anno.col <- as.data.frame(dataTables$meta_df[, id_mod, drop = F])
      # subset expression by markers
      all_marks <- marker_genes$All_Unique_Markers
      log_mat <- log(as.matrix(dataTables$exp_norm[all_marks,])+1)
      log_mat <- log_mat[complete.cases(log_mat),]
      
      withProgress(expr = pheatmap(log_mat, show_colnames = FALSE,
               show_rownames = TRUE,
               clustering_method = "ward.D2",
               annotation_col = anno.col),
               message = "Visualizing marker genes, please wait")
    })
    
    # render a table of the marker genes
    output$genes_table <- DT::renderDataTable({
      valid_load()
      valid_iasva()
      valid_markers()
      # create a data frame of marker genes
      marker_out <- matrix(data = NA, nrow = length(dataTables$markers$All_Unique_Markers), ncol = 2)
      colnames(marker_out) <- c("Gene", "SV")
      marker_out[,1] <- dataTables$markers$All_Unique_Markers
      am_mat <- do.call("cbind", dataTables$markers)
      # determine which SVs the gene is in
      SV_nams <- NULL
      if (nrow(marker_out) > 2) {
        for (sv_gen in 1:nrow(marker_out)) {
          sv_vec <- NULL
          for (sv_col in 1:(ncol(am_mat)-1)) {
            if (as.character(marker_out[sv_gen, 1]) %in% am_mat[, sv_col]) {
              sv_vec <- c(sv_vec, colnames(am_mat)[sv_col])
            } else {}
            # concatenate names
            if (length(sv_vec) > 1) {
              sv_fin <- paste(sv_vec, collapse = ", ")
            } else {
              sv_fin <- sv_vec
            }
          }
          SV_nams <- c(SV_nams, sv_fin)
        }
      } else if (nrow(marker_out) == 2) {
        SV_nams <- rep(names(dataTables$markers)[1], nrow(marker_out))
      }
      # append to matrix
      marker_out[,2] <- SV_nams
      marker_df <- as.data.frame(marker_out)
      dataTables$markers_formatted <- marker_df
      DT::datatable(marker_df)
    })
  })
  
  # gene enrichment analysis (using the markers_formatted data frame)
  # retreive gene information first based on species selection
  observeEvent(input$Path_Analysis, {
    # plot results
    output$Enrich_Plot <- renderPlot({
      valid_load()
      valid_iasva()
      valid_enrich()
      gene <- dataTables$markers_formatted[,1]
      if (isolate(input$Species_Type) == "Homo sapiens") {
        # convert from gene symbol to entrez id
        gene.df <- bitr(gene, fromType = "SYMBOL",
                        toType = c("ENSEMBL", "ENTREZID"),
                        OrgDb = org.Hs.eg.db)
        dataTables$species <- org.Hs.eg.db
        dataTables$species_kegg <- "hsa"
        # universe genes
        univ.df <- bitr(rownames(dataTables$exp_norm), fromType = "SYMBOL",
                        toType = c("ENSEMBL", "ENTREZID"),
                        OrgDb = org.Hs.eg.db)
        dataTables$univ_gene <- univ.df
        
      } else if (isolate(input$Species_Type) == "Mus musculus") {
        gene.df <- bitr(gene, fromType = "SYMBOL",
                        toType = c("ENSEMBL", "ENTREZID"),
                        OrgDb = org.Mm.eg.db)
        dataTables$species <- org.Mm.eg.db
        dataTables$species_kegg <- "mmu"
        # universe genes
        univ.df <- bitr(rownames(dataTables$exp_norm), fromType = "SYMBOL",
                        toType = c("ENSEMBL", "ENTREZID"),
                        OrgDb = org.Mm.eg.db)
        dataTables$univ_gene <- univ.df
        
      } else if (isolate(input$Species_Type == "Rattus norvegicus")) {
        gene.df <- bitr(gene, fromType = "SYMBOL",
                        toType = c("ENSEMBL", "ENTREZID"),
                        OrgDb = org.Rn.eg.db)
        dataTables$species <- org.Rn.eg.db
        dataTables$species_kegg <- "rno"
        # universe genes
        univ.df <- bitr(rownames(dataTables$exp_norm), fromType = "SYMBOL",
                        toType = c("ENSEMBL", "ENTREZID"),
                        OrgDb = org.Rn.eg.db)
        dataTables$univ_gene <- univ.df
      }
      # save to reactive value
      dataTables$gene.df <- gene.df
      
      # what type of analysis chosen
      if (isolate(input$Path_Type) == "Gene Ontology Biological Process") {
        ego <- withProgress(expr = enrichGO(gene = dataTables$gene.df$ENTREZID,
                                            OrgDb = dataTables$species,
                                            keyType = "ENTREZID",
                                            ont = "BP",
                                            pvalueCutoff = 0.05, pAdjustMethod = isolate(input$pvalue_correct),
                                            qvalueCutoff = isolate(input$path_cutoff),
                                            minGSSize = 5,
                                            readable = TRUE),
                            message = "Performing GO Biological Process Enrichment Analysis, please wait")
      } else if (isolate(input$Path_Type) == "Gene Ontology Cellular Component") {
        ego <- withProgress(expr = enrichGO(gene = dataTables$gene.df$ENTREZID,
                                            OrgDb = dataTables$species,
                                            keyType = "ENTREZID",
                                            ont = "CC",
                                            pvalueCutoff = 0.05, pAdjustMethod = isolate(input$pvalue_correct),
                                            qvalueCutoff = isolate(input$path_cutoff),
                                            minGSSize = 5,
                                            readable = TRUE),
                            message = "Performing GO Cellular Component Enrichment Analysis, please wait")
      } else if (isolate(input$Path_Type) == "Gene Ontology Molecular Function") {
        ego <- withProgress(expr = enrichGO(gene = dataTables$gene.df$ENTREZID,
                                            OrgDb = dataTables$species,
                                            keyType = "ENTREZID",
                                            ont = "MF",
                                            pvalueCutoff = 0.05, pAdjustMethod = isolate(input$pvalue_correct),
                                            qvalueCutoff = isolate(input$path_cutoff),
                                            minGSSize = 5,
                                            readable = TRUE),
                            message = "Performing GO Molecular Function Enrichment Analysis, please wait")
      } else if (isolate(input$Path_Type) == "KEGG") {
        ego <- withProgress(expr = enrichKEGG(gene = dataTables$gene.df$ENTREZID,
                                              organism = dataTables$species_kegg,
                                              keyType = "kegg",
                                              pvalueCutoff = 0.05, pAdjustMethod = isolate(input$pvalue_correct),
                                              qvalueCutoff = isolate(input$path_cutoff),
                                              minGSSize = 5),
                            message = "Performing KEGG Enrichment Analysis, please wait")
      } else if (isolate(input$Path_Type) == "Homo sapiens Immune Modules") {
        # load in necessary files
        gen_df <- read.delim("Data/Tmod.gene.to.immune.module.term.txt", header = F, check.names = F, stringsAsFactors = F)
        mod_df <- read.delim("Data/Tmod.immune.module.term.and.name.txt", header = T, check.names = F, stringsAsFactors = F)
        # enrichment analysis with gene symbols
        ego <- withProgress(expr = enricher(gene = dataTables$markers_formatted[,1],
                                              pvalueCutoff = 0.05, pAdjustMethod = isolate(input$pvalue_correct),
                                              qvalueCutoff = isolate(input$path_cutoff),
                                              minGSSize = 5, TERM2GENE = gen_df, TERM2NAME = mod_df),
                            message = "Performing Homo sapiens Immune Modules Enrichment Analysis, please wait")
      } else if (isolate(input$Path_Type) == "Homo sapiens PBMC Cell Specific Modules") {
        # load in necessary files
        gen_df <- read.delim("Data/Human.PBMC.Modules.and.genes.txt", header = F, check.names = F, stringsAsFactors = F)
        mod_df <- read.delim("Data/Human.PBMC.Modules.term.names.txt", header = F, check.names = F, stringsAsFactors = F)
        # enrichment analysis with gene symbols
        ego <- withProgress(expr = enricher(gene = dataTables$markers_formatted[,1],
                                            pvalueCutoff = 0.05, pAdjustMethod = isolate(input$pvalue_correct),
                                            qvalueCutoff = isolate(input$path_cutoff),
                                            minGSSize = 5, TERM2GENE = gen_df, TERM2NAME = mod_df),
                            message = "Performing Homo sapiens PBMC Cell Specific Modules Enrichment Analysis, please wait")
        
      } else if (isolate(input$Path_Type) == "Custom") {
        # load in gene and module file
        gen_df <- read.delim(file = isolate(input$input_gene_mod$datapath),
                                           header = F, check.names = F, stringsAsFactors = F)
        mod_df <- read.delim(file = isolate(input$input_mod_names$datapath),
                                                    header = F, check.names = F, stringsAsFactors = F)
        # enrichment analysis with gene symbols
        ego <- withProgress(expr = enricher(gene = dataTables$markers_formatted[,1],
                                            pvalueCutoff = 0.05, pAdjustMethod = isolate(input$pvalue_correct),
                                            qvalueCutoff = isolate(input$path_cutoff),
                                            minGSSize = 5, TERM2GENE = isolate(gen_df), TERM2NAME = isolate(mod_df)),
                            message = "Performing Custom Enrichment Analysis, please wait")
        
      }
      # save results to reactive value
      dataTables$enrich_res <- ego
      dp <- clusterProfiler::dotplot(object = dataTables$enrich_res, showCategory = isolate(input$path_viz_num)) + ggtitle(isolate(input$Path_Type))
      withProgress(expr = plot(dp),
                   message = "Visualizing gene enrichment results, please wait")
    
    })
    
    # render table of results
    output$enrich_table <- DT::renderDataTable({
      valid_load()
      valid_iasva()
      valid_enrich()
      DT::datatable(as.data.frame(dataTables$enrich_res))
    })
  })
  
  ## dimension reduction analysis
  observeEvent(input$Dim_Analysis, {
    if (isolate(input$Dim_Type) == "PCA") {
      # make interactive plots
      output$Dim_Plot_Orig <- renderPlotly({
        valid_load()
        valid_iasva()
        valid_dim()
        # original matrix plot
        # transpose matrix
        set.seed(1)
        trans_orig <- t(dataTables$exp_norm)
        # remove any zeros
        dataTables$pre_dim_orig <- trans_orig[, apply(trans_orig, 2, var, na.rm = TRUE) != 0]
        withProgress(expr = dim_orig <- prcomp(x = dataTables$pre_dim_orig, center = TRUE, scale. = TRUE),
                     message = "Performing PCA reduction of all genes, please wait")
        dim_orig_mat <- dim_orig$x
        rownames(dim_orig_mat) <- colnames(dataTables$exp_norm)
        
        # marker matrix plot
        # transpose matrix
        trans_mark <- t(dataTables$exp_norm[dataTables$markers_formatted[,1],])
        # remove any zeros
        dataTables$pre_dim_mark <- trans_mark[, apply(trans_mark, 2, var, na.rm = TRUE) != 0]
        withProgress(expr = dim_mark <- prcomp(x = dataTables$pre_dim_mark, center = TRUE, scale. = TRUE),
                     message = "Performing PCA reduction of IASVA selected genes, please wait")
        dim_mark_mat <- dim_mark$x
        rownames(dim_mark_mat) <- colnames(dataTables$exp_norm)
        dataTables$dim_mark <- as.data.frame(dim_mark_mat)
        dataTables$dim_orig <- as.data.frame(dim_orig_mat)
        
        # factor for coloring
        id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
        fac_int <- as.factor(dataTables$meta_df[, id_fac])
        plot_ly(dataTables$dim_orig, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = ""),
                color = ~fac_int) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = ""))
      })
      
      output$Dim_Plot_Markers <- renderPlotly({
        valid_load()
        valid_iasva()
        valid_dim()
        id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
        fac_int <- as.factor(dataTables$meta_df[, id_fac])
        plot_ly(dataTables$dim_mark, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = ""),
                color = ~fac_int) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), ")", sep = ""))
      })
      
    } else if (isolate(input$Dim_Type) == "t-SNE") {
      # make interactive plots
      output$Dim_Plot_Orig <- renderPlotly({
        valid_load()
        valid_iasva()
        valid_dim()
        # original matrix plot
        set.seed(1)
        trans_orig <- t(dataTables$exp_norm)
        # remove any zeros
        dataTables$pre_dim_orig <- trans_orig[, apply(trans_orig, 2, var, na.rm = TRUE) != 0]
        withProgress(expr = dim_orig <- Rtsne(X = dataTables$pre_dim_orig, dims = 3),
                     message = "Performing t-SNE reduction of all genes, please wait")
        dim_orig_mat <- dim_orig$Y
        rownames(dim_orig_mat) <- colnames(dataTables$exp_norm)
        colnames(dim_orig_mat) <- c("tSNE1", "tSNE2", "tSNE3")
        
        # marker matrix plot
        # transpose matrix
        trans_mark <- t(dataTables$exp_norm[dataTables$markers_formatted[,1],])
        # remove any zeros
        dataTables$pre_dim_mark <- trans_mark[, apply(trans_mark, 2, var, na.rm = TRUE) != 0]
        withProgress(expr = dim_mark <- Rtsne(X = dataTables$pre_dim_mark, dims = 3),
                     message = "Performing t-SNE reduction of IASVA selected genes, please wait")
        dim_mark_mat <- dim_mark$Y
        rownames(dim_mark_mat) <- colnames(dataTables$exp_norm)
        colnames(dim_mark_mat) <- c("tSNE1", "tSNE2", "tSNE3")
        # remove any zeros
        dim_mark_mat <- dim_mark_mat[, apply(dim_mark_mat, 2, var, na.rm = TRUE) != 0]
        
        dataTables$dim_mark <- as.data.frame(dim_mark_mat)
        dataTables$dim_orig <- as.data.frame(dim_orig_mat)
        
        # factor for coloring
        id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
        fac_int <- as.factor(dataTables$meta_df[, id_fac])
        plot_ly(dataTables$dim_orig, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = ""),
                color = ~fac_int) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = ""))
      })
      
      output$Dim_Plot_Markers <- renderPlotly({
        valid_load()
        valid_iasva()
        valid_dim()
        id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
        fac_int <- as.factor(dataTables$meta_df[, id_fac])
        plot_ly(dataTables$dim_mark, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = ""),
                color = ~fac_int) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), ")", sep = ""))
      })
    }
  })
  
  # update interactive plots
  observeEvent(input$Update_Dim, {
    if (isolate(input$Dim_Type) == "PCA") {
      # make interactive plots
      output$Dim_Plot_Orig <- renderPlotly({
        valid_load()
        valid_iasva()
        valid_dim()
        # factor for coloring
        id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
        fac_int <- as.factor(dataTables$meta_df[, id_fac])
        plot_ly(dataTables$dim_orig, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = ""),
                color = ~fac_int) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = ""))
      })
      
      output$Dim_Plot_Markers <- renderPlotly({
        valid_load()
        valid_iasva()
        valid_dim()
        id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
        fac_int <- as.factor(dataTables$meta_df[, id_fac])
        plot_ly(dataTables$dim_mark, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = ""),
                color = ~fac_int) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), ")", sep = ""))
      })
    } else if (isolate(input$Dim_Type) == "t-SNE") {
      # make interactive plots
      output$Dim_Plot_Orig <- renderPlotly({
        valid_load()
        valid_iasva()
        valid_dim()
        # factor for coloring
        id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
        fac_int <- as.factor(dataTables$meta_df[, id_fac])
        plot_ly(dataTables$dim_orig, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = ""),
                color = ~fac_int) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = ""))
      })
      
      output$Dim_Plot_Markers <- renderPlotly({
        valid_load()
        valid_iasva()
        valid_dim()
        # factor for coloring
        id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
        fac_int <- as.factor(dataTables$meta_df[, id_fac])
        plot_ly(dataTables$dim_mark, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = ""),
                color = ~fac_int) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), ")", sep = ""))
      })
    }
    
  })
  
  ########################################################################################
  # guided tutorial
  ########################################################################################
  # initiate hints on startup with custom button and event
  hintjs(session, options = list("hintButtonLabel"="Close Hint"),
         events = list())
  observeEvent(input$help,
               introjs(session, options = list("nextLabel"="Next Step",
                                               "prevLabel"="Previous Step",
                                               "skipLabel"="Quit"),
                       events = list("oncomplete"=I('alert("End of Tutorial")')))
  )
  
  ######################################################################################################
  # Download figures and tables
  ######################################################################################################
  # download figure panel
  output$Download_SVs <- downloadHandler(
    filename = function() { paste("IASVA", input$norm_method, "normalized", input$known_factors, "adjusted.SV.plot.pdf", sep=".") },
    content = function(file) {
      pdf(file)
      par(mar = c(5, 5, 4, 2), xpd = TRUE)
      valid_load()
      valid_iasva()
      valid_pair_plot()
      id_fac <- which(colnames(dataTables$meta_df) == input$All_SV_Color)
      fac_int <- as.factor(dataTables$meta_df[, id_fac])
      pairs(dataTables$iasva.res$sv[, c(as.numeric(input$All_SV_Num))], main="", pch=20, cex=0.5,
            col = color.vec[fac_int], bg = color.vec[fac_int], lower.panel = NULL)
      legend("bottomleft", levels(fac_int), fill=color.vec,
             title = as.character(colnames(dataTables$meta_df)[id_fac]))
      dev.off()
    })
  
  # download correlation plot
  output$Download_Correlation <- downloadHandler(
    filename = function() { paste("IASVA", input$norm_method, "normalized", input$known_factors, "adjusted.correlation.plot.pdf", sep=".") },
    content = function(file) {
      pdf(file)
      valid_load()
      valid_iasva()
      valid_pair_plot()
      # change factors to numeric for correlation
      meta_sel <- dataTables$meta_df
      for (jcol in 1:ncol(meta_sel)) {
        meta_sel[,jcol] <- as.numeric(as.factor(meta_sel[,jcol]))
      }
      iasva_vars <- cbind(dataTables$iasva.res$sv[, c(as.numeric(input$All_SV_Num))], meta_sel)
      # need to append column names to matrix
      colnames(iasva_vars) <- c(paste("SV", c(as.numeric(input$All_SV_Num)), sep = ""),
                                colnames(dataTables$meta_df))
      corrplot(cor(iasva_vars), type = "upper")
      dev.off()
    })
  
  # download pairwise plots
  output$Download_Pairwise <- downloadHandler(
    filename = function() { paste("IASVA", input$norm_method, "normalized", input$known_factors, "adjusted.pair.svs.html", sep=".") },
    content = function(file) {
      valid_load()
      valid_iasva()
      id_fac <- which(colnames(dataTables$meta_df) == input$Pair_SV_Color)
      fac_int <- as.factor(dataTables$meta_df[, id_fac])
      # sv axis labels
      x_ax <- list(
        title = as.character(input$sv_x)
      )
      y_ax <- list(
        title = as.character(input$sv_y)
      )
      # extract chosen SVs
      id_svx <- which(colnames(dataTables$iasva.res$sv) == input$sv_x)
      id_svy <- which(colnames(dataTables$iasva.res$sv) == input$sv_y)
      
      py <- plot_ly(as.data.frame(dataTables$iasva.res$sv), x = dataTables$iasva.res$sv[,id_svx],
              y = dataTables$iasva.res$sv[,id_svy], type = "scatter",
              mode = "markers", text = paste("Cell ID: ", rownames(dataTables$meta_df), sep = ""),
              color = ~fac_int) %>%
        layout(xaxis = x_ax, yaxis = y_ax)
      htmlwidgets::saveWidget(py, file)
  })
  
  # download heatmap of markers
  output$Download_Heatmap <- downloadHandler(
    filename = function() { paste("IASVA", input$norm_method, "normalized", input$known_factors, "adjusted.heatmap.plot.pdf", sep=".") },
    content = function(file) {
      valid_load()
      valid_iasva()
      valid_markers()
      pdf(file)
      id_mod <- which(colnames(dataTables$meta_df) %in% isolate(input$heatmap_known_factors))
      anno.col <- as.data.frame(dataTables$meta_df[, id_mod, drop = F])
      # subset expression by markers
      all_marks <- dataTables$markers$All_Unique_Markers
      log_mat <- log(as.matrix(dataTables$exp_norm[all_marks,])+1)
      log_mat <- log_mat[complete.cases(log_mat),]
      pheatmap(log_mat, show_colnames = FALSE,
               show_rownames = TRUE,
               clustering_method = "ward.D2",
               annotation_col = anno.col)
      dev.off()
    })
  
  # download marker genes in table
  output$Download_Markers <- downloadHandler(
    filename = function() {
      paste("IASVA", input$norm_method, "normalized", input$known_factors, "adjusted.marker.genes.csv", sep=".")
    },
    content = function(file) {
      valid_load()
      valid_iasva()
      valid_markers()
      write.csv(dataTables$markers_formatted, file, row.names = FALSE)
    }
  )
  
  # download enrichment analysis plot
  output$Download_path_plot <- downloadHandler(
    filename = function() { 
      paste("IASVA", input$norm_method, "normalized", input$known_factors, "adjusted.enrichment.plot.pdf", sep = ".")
    },
    content = function(file) {
      valid_load()
      valid_iasva()
      valid_markers()
      valid_enrich()
      pdf(file)
      dp <- clusterProfiler::dotplot(object = dataTables$enrich_res, showCategory = isolate(input$path_viz_num)) + ggtitle(isolate(input$Path_Type)) + 
        theme(text = element_text(size=8), axis.text.y = element_text(size = 8))
      plot(dp)
      dev.off()
    }
  )
  
  # download enrichment analysis table
  output$Download_path_table <- downloadHandler(
    filename = function() {
      paste("IASVA", input$norm_method, "normalized", input$known_factors, 
            "adjusted.enrichment.results.csv", sep=".")
    },
    content = function(file) {
      valid_load()
      valid_iasva()
      valid_markers()
      valid_enrich()
      write.csv(as.data.frame(dataTables$enrich_res), file, row.names = FALSE)
    }
  )
  
  # download dimension reduction plots
  output$download_dim_all <- downloadHandler(
    filename = function() { paste("IASVA", input$norm_method, "normalized", input$known_factors, "adjusted", input$Dim_Type, "all.genes.html", sep=".") },
    content = function(file) {
      valid_load()
      valid_iasva()
      valid_dim()
      if (input$Dim_Type == "PCA") {
        # factor for coloring
        id_fac <- which(colnames(dataTables$meta_df) == input$Dim_Color)
        fac_int <- as.factor(dataTables$meta_df[, id_fac])
        # make interactive plots
        py <- plot_ly(dataTables$dim_orig, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = ""),
                color = ~fac_int) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = ""))
        htmlwidgets::saveWidget(py, file)
      } else if (input$Dim_Type == "t-SNE") {
        # factor for coloring
        id_fac <- which(colnames(dataTables$meta_df) == input$Dim_Color)
        fac_int <- as.factor(dataTables$meta_df[, id_fac])
        # make interactive plots
        py <- plot_ly(dataTables$dim_orig, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = ""),
                  color = ~fac_int) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = ""))
        htmlwidgets::saveWidget(py, file)
      }
    }
  )
  output$download_dim_marker <- downloadHandler(
    filename = function() { paste("IASVA", input$norm_method, "normalized", input$known_factors, "adjusted", input$Dim_Type, "marker.genes.html", sep=".") },
        content = function(file) {
          valid_load()
          valid_iasva()
          valid_dim()
          if (input$Dim_Type == "PCA") {
            # factor for coloring
            id_fac <- which(colnames(dataTables$meta_df) == input$Dim_Color)
            fac_int <- as.factor(dataTables$meta_df[, id_fac])
            py <- plot_ly(dataTables$dim_mark, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                          mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = ""),
                          color = ~fac_int) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), ")", sep = ""))
            htmlwidgets::saveWidget(py, file)
          } else if (input$Dim_Type == "t-SNE") {
            id_fac <- which(colnames(dataTables$meta_df) == input$Dim_Color)
            fac_int <- as.factor(dataTables$meta_df[, id_fac])
            py <- plot_ly(dataTables$dim_mark, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                          mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = ""),
                          color = ~fac_int) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), ")", sep = ""))
            htmlwidgets::saveWidget(py, file)
          }
        }
    )
  
  
  # download manual/instructions
  output$Download_Manual <- downloadHandler(
    filename = function() {
      paste("Data/IASVA_Shiny_App_Manual.docx", sep = "")
    },
    content = function(file) {
      file.copy("Data/IASVA_Shiny_App_Manual.docx", file)
    }
  )
  
  # download test data sets (csv files)
  output$Download_Test_Exp <- downloadHandler(
    filename = function() {
      paste("Data/test.exp.csv", sep = "")
    },
    content = function(file) {
      file.copy("Data/test.exp.csv", file)
    }
  )
  
  output$Download_Test_Meta <- downloadHandler(
    filename = function() {
      paste("Data/test.metadata.csv", sep = "")
    },
    content = function(file) {
      file.copy("Data/test.metadata.csv", file)
    }
  )
  
})

##################################################################################################
# Run the app 
##################################################################################################
shinyApp(ui = ui, server = server)
