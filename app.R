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
library(shinythemes)
library(rintrojs)
library(htmlwidgets)

################################################################################################################
# Define UI
################################################################################################################
ui <- shinyUI(fluidPage(theme = shinytheme("cerulean"),
                        
                        tags$head(includeScript("google-analytics.js")),
                        
                        # Application title
                        titlePanel(div("VIASVA: A Shiny app for Visual Iteratively Adjusted Surrogate Variable Analysis", 
                                       img(height = 88, width = 217, 
                                           src = "jax.logo.gif", 
                                           class = "pull-right"))),
                        # side bar where user first interacts 
                        sidebarPanel(
                          # call function to use package rintrojs                                
                          introjsUI(),
                          # add first button for guided tutorial
                          div(style="display: inline-block;vertical-align:top; width: 200px;",
                          introBox(
                            actionButton("help", "Guided Tutorial!", icon = icon("list-ol")),
                            data.step = 1,
                            data.intro = "This is the start of the tutorial. Click 'Next Step' to continue."
                          )
                          ),
                          
                          # Or provide download button for manual
                          div(style="display: inline-block;vertical-align:top; width: 200px;",
                          introBox(
                            downloadButton("Download_Manual", "Download Tutorial"),
                            data.step = 2,
                            data.intro = "Click this button to download an in-depth tutorial/manual for this app."
                          )
                          ),
                          
                          # download test datasets
                          div(style="display: inline-block;vertical-align:top; width: 200px;",
                          introBox(
                            downloadButton("Download_Test_Exp", "Download Test Expression"),
                            data.step = 3,
                            data.intro = "Click this button to download an example gene expression file. This file may be used to test out the app."
                            )
                          ),
                          
                          div(style="display: inline-block;vertical-align:top; width: 200px;",
                          introBox(
                            downloadButton("Download_Test_Meta", "Download Test Metadata"),
                            data.step = 4,
                            data.intro = "Click this button to download an example sample metadata file. This file may be used to test out the app."
                            )
                          ),
                      
                          # Load user data
                          introBox(
                            h4("1. Load Data"),
                            fileInput(
                              inputId = "input_exprs",
                              "(Required) Choose gene/ADT file to upload",
                              accept = c(
                                ".Rds",
                                ".csv",
                                ".txt"
                              )
                            ),
                            
                            # conditional panel if they have metadata
                            div(style="display: inline-block;vertical-align:top; width: 200px;",
                                selectInput(inputId = "have_meta", label = "Do You Have Sample Metadata?",
                                            choices = c("Yes", "No"), 
                                            selected = "No")
                            ),
                            conditionalPanel(
                              condition = "input.have_meta == 'Yes'",
                              fileInput(
                                inputId = "input_meta",
                                "(Optional) Choose sample metadata file to upload",
                                accept = c(
                                  ".Rds",
                                  ".csv",
                                  ".txt"
                                )
                              )
                            ),
                            data.step = 5, data.intro = "Click the Browse buttons to provide your gene/ADT expression data.
                            Gene expression data should be in the format of: gene symbol names in rows, sample names in columns.
                            Providing sample metadata is optional. If you do choose to provide metadata, the file should be formatted such that
                            sample names are in the rows and meta data variables are in columns. Step 1: The following input file formats are accepted: .Rds, tab-delimited text, .csv."),
                          
                          # optional step to down-sample large datasets
                          introBox(
                            h5("2a. Cell/Sample Filtering (Optional)"),
                            div(style="display: inline-block;vertical-align:top; width: 200px;",
                                selectInput(inputId = "cellfilt", label = "Filter Cells/Samples?",
                                            choices = c("Yes", "No"), 
                                            selected = "No")
                            ),
                            conditionalPanel(
                              condition = "input.cellfilt == 'Yes'",
                              div(style="display: inline-block;vertical-align:top; width: 200px;",
                                  numericInput(inputId = "cellfilt_num", label = "Minimum Genes Detected in each Cell",
                                               min = 0, value = 500)),
                              br(),
                              shiny::actionButton("Run_cellfilt", "Run Cell Filtering", icon = icon("paper-plane"))
                            ),
                            data.step = 6, data.intro = "Remove low-quality cells based on the number of genes detected (at least one read count) in each cell.
                            Please select yes or no from the drop down menu to indicate if you would like to filter cells.
                            If you choose yes, specify a minimum number of genes detected in each cell."
                          ),
                          
                          introBox(
                            h5("2b. Down-sample Cells (Optional)"),
                            div(style="display: inline-block;vertical-align:top; width: 200px;",
                                selectInput(inputId = "downsample", label = "Down-sample?",
                                            choices = c("Yes", "No"), 
                                            selected = "No")
                            ),
                            conditionalPanel(
                              condition = "input.downsample == 'Yes'",
                              div(style="display: inline-block;vertical-align:top; width: 200px;",
                                  numericInput(inputId = "downsample_num", label = "Minimum number of cells",
                                               min = 0, value = 100)),
                              br(),
                              shiny::actionButton("Run_downsample", "Run Cell Down-sampling", icon = icon("paper-plane"))
                            ),
                            data.step = 7, data.intro = "For large datasets, users may choose to randomly down-sample to a specified number of samples to decrease
                            computational time for downstream analyses. Note down-sampling may reduce power of analyses. Please select yes or no from the drop down menu to indicate if you would like to down-sample your data.
                            If you choose yes, specify a minimum number of samples to randomly down-sample to."
                          ),
                          introBox(
                            h5("2c. Gene Filtering and Normalization (Optional)"),
                            div(style="display: inline-block;vertical-align:top; width: 200px;",
                                selectInput(inputId = "gene_filter", label = "Filter Genes?",
                                            choices = c("Yes", "No"), 
                                            selected = "No")
                            ),
                            conditionalPanel(
                              condition = "input.gene_filter == 'Yes'",
                              div(style="display: inline-block;vertical-align:top; width: 200px;",
                                  numericInput(inputId = "Count_num", label = "# of Counts in each Cell",
                                               min = 0, value = 5)),
                              div(style="display: inline-block;vertical-align:top; width: 200px;",
                                  numericInput(inputId = "Cell_num", label = "Detected in # of Cells",
                                               min = 0, value = 5)),
                              div(style="display: inline-block;vertical-align:top; width: 200px;",
                                  selectInput(inputId = "norm_method", label = "Normalization Method",
                                              choices = c("CPM", "Quantile", "scran", "None"), selected = "CPM")),
                              br(),
                              shiny::actionButton("Run_gene_filter", "Run Gene Filtering", icon = icon("paper-plane"))
                            ),
                            data.step = 8, data.intro = "Users have the option to remove lowly expressed genes from their data. If 'yes' is chosen, please specify the criteria to filter genes. E.g., genes with 3 or more counts in at least 5 cells. Next, from the drop-down menu, select your desired normalization method. If no is chosen, move on to the next step.
                            CPM = Counts per million from edgeR package, Quantile = quantile normalization from preprcessCore R package,
                            scran = method to deconvolute size factors from cell pools from scran R package." 
                          ),
                          introBox(
                            # after the user loads the data, update the select input options to reflect metadata table column names
                            h4("3. Specify Known Factors to Adjust For"),
                            div(style="display: inline-block;vertical-align:top; width: 200px;",
                                selectInput(inputId = "known_factors", label = "Known Factors",
                                            choices = c(""), selected = NULL, multiple = TRUE)),
                            data.step = 9, data.intro = "Click in the box to select the known factor(s) from your metadata file you would like to adjust for. You may choose a single or multiple factors.
                            If no metadata file was provided, two factors Genes_Detected (number of genes detected in each sample) and Log_Total_Counts (log transformed total read counts in each sample) will be calculated for you.
                            Please note, at least one known factor must be selected prior to IA-SVA analysis."
                          ),
                          introBox(
                            # after the user loads the data, update the select input options to reflect metadata table column names
                            h4("4. IA-SVA Analysis"),
                            
                            # either choose pct cutoff or num.sv parameter
                            div(style="display: inline-block;vertical-align:top; width: 200px;",
                                selectInput(inputId = "iasva_param", label = "Parameter Choice",
                                            choices = c("Percentage Threshold", "Number of SVs"), 
                                            selected = "Percentage Threshold")),
                            
                            # conditional panel to choose parameter for IA-SVA analysis (percent threshold)
                            conditionalPanel(
                              condition = "input.iasva_param == 'Percentage Threshold'",
                              # options to customize parameters
                              div(style="display: inline-block;vertical-align:top; width: 200px;",
                                  numericInput(inputId = "pct_cutt", label = "% Threshold for SV retention", value = 1,
                                               min = 1, max = 99))
                            ),
                            
                            # conditional panel to choose parameter for IA-SVA analysis (number of sv's)
                            conditionalPanel(
                              condition = "input.iasva_param == 'Number of SVs'",
                              # options to customize parameters
                              div(style="display: inline-block;vertical-align:top; width: 200px;",
                                  numericInput(inputId = "num_of_svs", label = "Number of SVs to Estimate", value = 5,
                                               min = 1, max = 99))
                            ),
                            
                            br(),
                            shiny::actionButton("do", "Run Analysis", icon = icon("paper-plane"), class = "btn-primary"),
                            data.step = 10, data.intro = "Click this button to run IA-SVA and identify unknown sources (surrogate variables) of variation within your data.
                            Step 4: 'Percentage threshold for SV retention': IA-SVA computes the percentage of unmodeled variance explained by the putative hidden factor and compare it with the user-defined threshold. If the percentage is greater than the threshold, SV is retained.
                            Alternatively, 'Number of SVs' forces IA-SVA to estimate a specified number of SVs."
                          ),
                          
                          br(),
                          h4("Author/Contact: "),
                          introBox(
                            tags$div(class = "header", checked = NA,  
                                     tags$i("Nathan Lawlor (nathan.lawlor03@gmail.com)")
                                     ),
                          h4("Questions/Issues?: "),
                          tags$div(class = "header", checked = NA,
                                   tags$i(""),
                                   tags$p(""),
                                   tags$a(href = "https://github.com/nlawlor/iasva_shiny",
                                          "Visit the App Github Page Here!", target = "_blank")),
                          h4("Other Resources: "),
                          tags$div(class = "header", checked = NA,
                            tags$i(""),
                            tags$p(""),
                            tags$a(href = "https://www.jax.org/research-and-faculty/research-labs/the-ucar-lab",
                                 "Visit the Ucar Lab Here!", target = "_blank")),
                          tags$div(class = "header", checked = NA,  
                                   tags$i(""),
                                   tags$p(""),
                                   tags$a(href = "https://www.bioconductor.org/packages/devel/bioc/html/iasva.html",
                                          "Check out the R package on Bioconductor", target = "_blank")),
                          tags$div(class = "header", checked = NA,  
                                   tags$i(""),
                                   tags$p(""),
                                   tags$a(href = "https://www.biorxiv.org/content/early/2018/04/24/151217",
                                          "Check out the IASVA manuscript on bioRxiv", target = "_blank")),
                          data.step = 11, data.intro = "Please contact me or visit the Github page with any questions about the app. Also, checkout the resources below!")
                        ),
                        # main panel with multiple tabs: detect SVs, find marker genes, annotate genes, visualize data with tsne
                        introBox(
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Data/QC",
                                       verbatimTextOutput("exp_load"),
                                       br(),
                                       # add two plots below this text output with qc stats
                                       downloadButton("Download_QC", "Download Histogram"),
                                       plotOutput("Detect", height = 600, width = 600),
                                       br(),
                                       verbatimTextOutput("meta_load"),
                                       verbatimTextOutput("cell_filtering_summ"),
                                       verbatimTextOutput("cell_downsampling_summ"),
                                       verbatimTextOutput("gene_filtering_summ"),
                                       verbatimTextOutput("load_data")
                              ),
                              tabPanel("Surrogate Variables", 
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
                                       div(style="display: inline-block;vertical-align:left; width: 800px;",
                                           h4("Correlation Plot"),
                                           div(style="display: inline-block;vertical-align:left; width: 200px;",
                                               downloadButton("Download_Correlation", "Download Correlation Plot")
                                           ),
                                           div(style="display: inline-block;vertical-align:left; width: 200px;",
                                               downloadButton("Download_Correlation_Table", "Download Correlation Table")
                                           ),
                                           plotOutput("corrplot", width = 800, height = 800)
                                       ),
                                       div(style="display: inline-block;vertical-align:right; width: 800px;",
                                           h4("Surrogate Variables Plot"),
                                           div(style="display: inline-block;vertical-align:left; width: 200px;",
                                               downloadButton("Download_SVs", "Download SV Plot")
                                           ),
                                           div(style="display: inline-block;vertical-align:left; width: 200px;",
                                               downloadButton("Download_SVs_Table", "Download SV Table")
                                           ),
                                           plotOutput("SVplot", width = 800, height = 800)
                                           
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
                                                shiny::actionButton("updatePairSV", "Update Pairwise SV Plot", icon = icon("refresh"), class = "btn-primary"),
                                                align="left"),
                                         column(width = 9, plotlyOutput("PairSV", height = 600, width = 600)))),
                              tabPanel("Identify Marker Genes", 
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
                                       downloadButton("Download_Heatmap", "Download Heatmap"),
                                       plotOutput("MarkerHeatmap", height = 800, width = 800),
                                       
                                       # add marker genes after heatmap
                                       br(),
                                       h4("Marker Genes Table: "),
                                       downloadButton("Download_Markers", "Download Marker Genes Table"), DT::dataTableOutput("genes_table")
                              ),
                              
                              # dimension reduction and clustering
                              tabPanel("Visualization with Marker Genes",
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
                                       div(style="display: inline-block;vertical-align:left; width: 200px;",
                                           downloadButton("download_dim_all", "Download All Genes Plot")
                                       ),
                                       div(style="display: inline-block;vertical-align:left; width: 200px;",
                                           downloadButton("download_all_coordinates", "Download All Genes Coordinates")   
                                       ),
                                       br(),
                                       # plot with all genes
                                       div(style="display: inline-block;vertical-align:left; width: 1000px;",
                                           plotlyOutput("Dim_Plot_Orig", height = 800, width = 1000)),
                                       br(),
                                       # plot with ia-sva marker genes
                                       div(style="display: inline-block;vertical-align:left; width: 250px;",
                                           downloadButton("download_dim_marker", "Download Marker Genes Plot")
                                       ),
                                       div(style="display: inline-block;vertical-align:left; width: 250px;",
                                           downloadButton("download_marker_coordinates", "Download Marker Genes Coordinates")
                                       ),
                                       div(style="display: inline-block;vertical-align:left; width: 1000px;",
                                           plotlyOutput("Dim_Plot_Markers", height = 800, width = 1000))
                              ),
                              
                              # gene enrichment analysis figure and table
                              tabPanel("Gene Enrichment Analysis",
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
                                       downloadButton("Download_path_plot", "Download Enrichment Analysis Plot"),
                                       plotOutput("Enrich_Plot", height = 800, width = 800),
                                       
                                       # tabPanel("Gene Enrichment Analysis Table", 
                                       # include table of pathway results
                                       br(),
                                       h4("Gene Enrichment Analysis Table"),
                                       downloadButton("Download_path_table", "Download Enrichment Analysis Results Table"), 
                                       DT::dataTableOutput("enrich_table")
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
    meta_df = NULL,
    cell_filter_choice = FALSE,
    cell_downsample_choice = FALSE,
    gene_filter_choice = FALSE,
    detect_num = NULL,
    exp_norm = NULL,
    iasva.res = NULL,
    summ_exp = NULL,
    markers = NULL,
    iasva_vars = NULL,
    markers_formatted = NULL,
    chosen_svs = NULL,
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
      withProgress(expr = dataTables$exp_norm <- read.csv(file = input$input_exprs$datapath,
                                     header = T, check.names = F, stringsAsFactors = F, row.names = 1),
                   message = "Loading expression file, please wait")
    } else if (grepl(x = input$input_exprs$datapath, pattern = ".txt")) {
      withProgress(expr = dataTables$exp_norm <- read.delim(file = input$input_exprs$datapath,
                                       header = T, check.names = F, stringsAsFactors = F, row.names = 1),
                   message = "Loading expression file, please wait")
    } else if (grepl(x = input$input_exprs$datapath, pattern = ".Rds")) {
      withProgress(expr = dataTables$exp_norm <- readRDS(file = input$input_exprs$datapath),
                   message = "Loading expression file, please wait")
    }
    # change file to matrix
    dataTables$exp_norm <- as.matrix(dataTables$exp_norm)
    
    # message to display
    output$exp_load <- renderText({
      print(paste("Expression file loaded! \n", 
                  "Gene number: ", isolate(nrow(dataTables$exp_norm)), " \n",
                  "Cell number: ", isolate(ncol(dataTables$exp_norm)), " \n",
                  sep = ""))
    })
    
    # display a histogram of genes detected to help with filtering of cells
    output$Detect <- renderPlot({
      # validate
      shiny::validate(
        need(isolate(input$input_exprs$datapath != ""), "Error: No gene/ADT file loaded")
      )
      isolate({
        withProgress(message = "Plotting the number of detected genes in each sample, please wait", {
          # binarize data
          bin_data <- dataTables$exp_norm
          bin_data[bin_data < 1] <- 0
          bin_data[bin_data >= 1] <- 1
          num.exp <- apply(bin_data,2,sum)
          dataTables$detect_num <- num.exp
          summ <- summary(num.exp)
          hist(num.exp, col = "dodgerblue", main="", 
               ylab = "Samples (n)", xlab = "Number of genes detected in each sample")
          legend("topright", legend = paste(names(summ), round(summ, digits = 2), sep = " "), title = "Summary of genes detected")
      
          # create a meta data object by default
          dataTables$meta_df <- data.frame(Genes_Detected = num.exp, Log_Total_Counts = log(colSums(dataTables$exp_norm)+1))
          # change all variables to factors
          for (fac in 1:ncol(dataTables$meta_df)) {
            dataTables$meta_df[, fac] <- as.factor(dataTables$meta_df[, fac])
          }
          # identify variables in meta data 
          updateSelectInput(session = session, inputId = "known_factors", label = "Known Factors",
                            choices = colnames(dataTables$meta_df))
          updateSelectInput(session = session, inputId = "heatmap_known_factors", label = "Known Factors to Include in Heatmap",
                            choices = c(colnames(dataTables$meta_df), "None"))
      
          # make progress bar seen
          for (i in 1:5) {
            incProgress(1/5)
            Sys.sleep(0.25)
          }
        })
      })
    })
  })
  
  observeEvent(input$input_meta, {
    # load in metadata
    if (grepl(x = input$input_meta$datapath, pattern = ".csv")) {
      withProgress(expr = dataTables$meta_df <- read.csv(file = input$input_meta$datapath,
                                     header = T, check.names = F, stringsAsFactors = F, row.names = 1),
                   message = "Loading metadata file, please wait")
    } else if (grepl(x = input$input_meta$datapath, pattern = ".txt")) {
      withProgress(expr = dataTables$meta_df <- read.delim(file = input$input_meta$datapath,
                                       header = T, check.names = F, stringsAsFactors = F, row.names = 1),
                   message = "Loading metadata file, please wait")
    } else if (grepl(x = input$input_meta$datapath, pattern = ".Rds")) {
      withProgress(expr = dataTables$meta_df <- readRDS(file = input$input_meta$datapath),
                   message = "Loading metadata file, please wait")
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
                      choices = c(colnames(dataTables$meta_df), "None"))
    
    # render summary of metadata file
    output$meta_load <- renderText({
      # if metadata and exp files are different
      # the expression and metadata tables need to have same number of rows and columns
      shiny::validate(
        need(isolate(ncol(dataTables$exp_norm)) == isolate(nrow(dataTables$meta_df)), "Error: Sample metadata file row number is not equal to Expression file column number!")
      )
      
      print(paste("Metadata file loaded! \n", 
                  "Cell number: ", isolate(nrow(dataTables$meta_df)), " \n",
                  "# of metadata variables: ", isolate(ncol(dataTables$meta_df)), " \n",
                  sep = ""))
    })
  })
  
  
  # shiny app error validation for loading data, and plots, tables
    # validate functions must go inside the plot/table rendering functions
      # option could be to have a data loading panel in which we provide qc and make sure the data is loaded properly
      # this first panel would also conduct all iasva analyses
  
  valid_load <- function() {
    shiny::validate(
      # adding an isolate() will make the error messages only appear/disappear when click run analysis button
      # before can run ia-sva need to have loaded data
      need(isolate(dataTables$exp_norm != ""), "Error: Please load a gene/ADT file before running IA-SVA"),
      need(isolate(dataTables$meta_df != ""), "Error: Please load a metadata file before running IA-SVA"),
      # make sure specify known variables
      need(isolate(input$known_factors != ""), "Error: Please specify known factor(s) to adjust for before running IA-SVA"),
      # maybe have option to specify no variables?
      need(isolate(input$pct_cutt >= 1), "Error: Percent threshold must be greater than or equal to 1"),
      need(isolate(input$num_of_svs >= 1), "Error: Number of SVs to estimate must be greater than or equal to 1"),
      need(isolate(input$pct_cutt < 100), "Error: Percent threshold must be less than 100"),
      # cell filtering
      need(isolate(input$Cell_num >= 0), "Error: Number of cells must be greather than or equal to 0"),
      need(isolate(input$Count_num >= 0), "Error: Number of gene/ADT counts must be greather than or equal to 0"),
      # check for numeric inputs
      need(isolate(is.numeric(input$pct_cutt)), "Error: Percent threshold must be numeric"),
      need(isolate(is.numeric(input$num_of_svs)), "Error: Number of SVs to estimate must be a numeric value"),
      need(isolate(is.numeric(input$Cell_num)), "Error: Number of cells must be numeric"),
      need(isolate(is.numeric(input$Count_num)), "Error: Number of gene/ADT counts must be numeric"),
      # the expression and metadata tables need to have same number of rows and columns
      need(ncol(dataTables$exp_norm) == nrow(dataTables$meta_df), "Error: Sample metadata file row number is not equal to Expression file column number!")
    )
  }
  
  valid_downsample <- function() {
    shiny::validate(
      # if user chooses a downsampling number, make sure the number is greater than the input expression size
      need(isolate(dataTables$exp_norm != ""), "Error: Please load a gene/ADT file before running IA-SVA"),
      need(isolate(dataTables$meta_df != ""), "Error: Please load a metadata file before running IA-SVA"),
      need(ncol(dataTables$exp_norm) == nrow(dataTables$meta_df), "Error: Sample metadata file row number is not equal to Expression file column number!"),
      need(isolate(is.numeric(input$downsample_num)), "Error: Down-sampling value must be numeric"),
      need(isolate(ncol(dataTables$exp_norm) >= input$downsample_num),  "Error: Down-sampling value too large, must be less than input sample size")
    )
  }
  
  valid_iasva <- function() {
    shiny::validate(
      need(isolate(!is.null(dataTables$iasva.res)), "Error: No surrogate variables obtained... Please adjust input parameters")
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
      need(isolate(input$heatmap_known_factors != ""), "Error: Please specify known factor(s) to include in heatmap visualization"),
      need(isolate(input$mark_cutoff > 0), "Error: Adjusted p-value cutoff must be greater than 0"),
      need(isolate(input$mark_cutoff <= 1), "Error: Adjusted p-value cutoff must be less than or equal to 1"),
      need(isolate(input$rsqcutoff > 0), "Error: R-squared cutoff must be greater than 0"),
      need(isolate(input$rsqcutoff <= 1), "Error: R-squared cutoff must be less than or equal to 1")
    )
  }
  
  # validation for enrichment analysis
  valid_enrich <- function() {
    shiny::validate(
      need(isolate(length(as.numeric(input$SV_marks)) >= 1), "Error: Please select 1 or more SVs to identify markers for"),
      need(isolate(!is.null(dataTables$markers)), "Error: No marker genes were identified... Please adjust r-squared and significance cutoffs and re-calculate before proceeding with enrichment analysis"),
      need(isolate(input$path_cutoff > 0), "Error: Enrichment analysis adjusted p-value cutoff must be greater than 0"),
      need(isolate(input$path_cutoff <= 1), "Error: Enrichment analysis adjusted p-value cutoff must be less than or equal to 1")
    )
  }
  
  # shiny app error validation for dim plots
  valid_dim <- function() {
    shiny::validate(
      need(isolate(!is.null(dataTables$markers)), "Error: No marker genes were identified... Please adjust r-squared and significance cutoffs and re-calculate before proceeding with dimension reduction analysis")
    )
  }
  
  # observe for filtering of cells
  observeEvent(input$Run_cellfilt, {
    # filter cells based on numbers of genes detected (at least one count)
    withProgress(message = paste("Removing cells that have less than ", isolate(input$cellfilt_num),
                                 " genes detected, please wait", sep = ""), {
                                   # use reactive value
                                   num.sel <- dataTables$detect_num[dataTables$detect_num >= isolate(input$cellfilt_num)]
                                   # make progress bar seen
                                   for (i in 1:5) {
                                     incProgress(1/5)
                                     Sys.sleep(0.25)
                                   }
                                   # subset data
                                   dataTables$exp_norm <- dataTables$exp_norm[, names(num.sel)]
                                   dataTables$meta_df <- dataTables$meta_df[names(num.sel), ]
                                   # change reactive value
                                   dataTables$cell_filter_choice <- TRUE
                                 })
  })
  
  
  # observe for down-sampling of data
    observeEvent(input$Run_downsample, {
      # ensure down-sample value is correct
        valid_downsample()
        set.seed(1)
        dw_samp <- base::sample(x = isolate(1:ncol(dataTables$exp_norm)), size = isolate(input$downsample_num))
        # subset data
        withProgress(message = paste("Down-sampling data to ", isolate(input$downsample_num), " samples, please wait", sep = ""), {
          dataTables$exp_norm <- dataTables$exp_norm[, dw_samp]
          # make progress bar seen
          for (i in 1:5) {
            incProgress(1/5)
            Sys.sleep(0.25)
          }
        })
        dataTables$meta_df <- dataTables$meta_df[dw_samp, ]
        # change reactive value
        dataTables$cell_downsample_choice <- TRUE
  })
  
  # filtering of genes
  observeEvent(input$Run_gene_filter, {
      # by filter genes with low counts
      withProgress(expr = filter <- apply(dataTables$exp_norm, 1, function(x) length(x[x>isolate(input$Count_num)])>=isolate(input$Cell_num)),
                   message = "Filtering genes, please wait")
      dataTables$exp_norm <- dataTables$exp_norm[filter,]
      
      # normalize the data
      if (isolate(input$norm_method) == "CPM") {
        dataTables$exp_norm <- edgeR::cpm(dataTables$exp_norm)
      } else if (isolate(input$norm_method) == "Quantile") {
        dataTables$exp_norm <- normalize.quantiles(dataTables$exp_norm)
      } else if (isolate(input$norm_method) == "scran") {
        sce <- SingleCellExperiment(list(counts=dataTables$exp_norm))
        sce <- computeSumFactors(sce)
        sce <- normalize(sce)
        dataTables$exp_norm <- exprs(sce)
      } else if (isolate(input$norm_method) == "None") {
        dataTables$exp_norm <- dataTables$exp_norm
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
      # update reactive value
      dataTables$gene_filter_choice <- TRUE
  })
  
  # output QC messages
  output$cell_filtering_summ <- renderText({
    if (dataTables$cell_filter_choice) {
      print(paste("Cell filtering finished! \n",
                  "Current cell number: ", isolate(ncol(dataTables$exp_norm)), " \n",
                  sep = ""))
    }
  })
  
  output$cell_downsampling_summ <- renderText({
    if (dataTables$cell_downsample_choice) {
      valid_downsample()
      print(paste("Down-sampling finished! \n", 
                  "Current cell number: ", isolate(ncol(dataTables$exp_norm)), " \n", 
                  sep = ""))
    }
  })
  
  output$gene_filtering_summ <- renderText({
    if (dataTables$gene_filter_choice) {
      print(paste("Gene filtering and normalization finished! \n", 
                  "Current gene number: ", isolate(nrow(dataTables$exp_norm)), " \n", 
                  sep = ""))
    }
  })
  
  # extract known factors to adjust for, create model matrix
  observeEvent(input$do, {
    output$load_data <- renderText({
      # error handling
      valid_load()
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
      
      # depending on which ia-sva parameters were chosen, evaluate
      if (isolate(input$iasva_param == "Percentage Threshold")) {
        dataTables$iasva.res <- withProgress(expr = fast_iasva(summ_exp, mod[,-1, drop = F], verbose=FALSE,
                                                               pct.cutoff = isolate(input$pct_cutt), num.sv = NULL),
                                             message = "IA-SVA analysis in progress, please wait")
      } else if (isolate(input$iasva_param == "Number of SVs")) {
        dataTables$iasva.res <- withProgress(expr = fast_iasva(summ_exp, mod[,-1, drop = F], verbose=FALSE,
                                                               pct.cutoff = isolate(input$pct_cutt), num.sv = isolate(input$num_of_svs)),
                                             message = "IA-SVA analysis in progress, please wait")
      }
      
      # if no SV's are calculated inform user
      valid_iasva()
      
      # display SV's
      iasva.sv <- as.data.frame(dataTables$iasva.res$sv)
      rownames(iasva.sv) <- colnames(dataTables$exp_norm)
      
      # update coloring scheme for total plot
      updateSelectInput(session = session, inputId = "All_SV_Color", label = "Color Points",
                        choices = c(colnames(dataTables$meta_df), "No Color"), selected = "No Color")
      # update number of svs to visualize
      updateSelectInput(session = session, inputId = "All_SV_Num", label = "SV Number",
                        choices = 1:ncol(iasva.sv), selected = 1:ncol(iasva.sv))
      
      # update coloring scheme for pairwise
      updateSelectInput(session = session, inputId = "Pair_SV_Color", label = "Color Points",
                        choices = c(colnames(dataTables$meta_df), "No Color"), selected = "No Color")
      # update coloring scheme for dimension reduction plots
      updateSelectInput(session = session, inputId = "Dim_Color", label = "Color Points",
                        choices = c(colnames(dataTables$meta_df), "No Color"), selected = "No Color")
      # update sv x/y selections
      updateSelectInput(session = session, inputId = "sv_x", label = "SV X-axis",
                        choices = colnames(iasva.sv), selected = "SV1")
      updateSelectInput(session = session, inputId = "sv_y", label = "SV Y-axis",
                        choices = colnames(iasva.sv), selected = "SV2")
      # update sv selections in marker panel
      updateSelectInput(session = session, inputId = "SV_marks", label = "Choose SVs",
                        choices = colnames(iasva.sv), selected = "SV1")
      
      # output messages
      print(paste("IA-SVA analysis finished! \n",
                  "Number of SV's identified: ", isolate(ncol(dataTables$iasva.res$sv)), sep = ""))
      
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
      col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
      corrplot(abs(cor(iasva_vars)), type = "upper", method = "color",
               col = col(200), number.cex = 1,
               addCoef.col = "black",
               tl.col = "black", tl.srt = 90, diag = FALSE)
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
        if (isolate(input$All_SV_Color == "No Color")) {
          pairs(dataTables$iasva.res$sv[, c(as.numeric(input$All_SV_Num))], main="", pch=20, cex=0.5,
               lower.panel = NULL)
        } else {
          id_fac <- which(colnames(dataTables$meta_df) == input$All_SV_Color)
          fac_int <- as.factor(dataTables$meta_df[, id_fac])
          pairs(dataTables$iasva.res$sv[, c(as.numeric(input$All_SV_Num))], main="", pch=20, cex=0.5,
                col = color.vec[fac_int], bg = color.vec[fac_int], lower.panel = NULL)
          legend("bottomleft", levels(fac_int), fill=color.vec,
                 title = as.character(colnames(dataTables$meta_df)[id_fac]))
        }
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
        col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
        corrplot(abs(cor(iasva_vars)), type = "upper", method = "color",
                 col = col(200), number.cex = 1,
                 addCoef.col = "black",
                 tl.col = "black", tl.srt = 90, diag = FALSE)
      })
    })
  })
  
  # generate updated plotly plot
  observeEvent(input$updatePairSV, {
    output$PairSV <- renderPlotly({
      valid_load()
      valid_iasva()
      isolate ({
        if (isolate(input$Pair_SV_Color == "No Color")) {
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
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$meta_df), sep = "")) %>%
            layout(xaxis = x_ax, yaxis = y_ax)
        } else {
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
        }
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
      dataTables$chosen_svs <- as.character(paste(colnames(dataTables$iasva.res$sv)[id_sv_mark], collapse = ","))
      
      # if no coloring on heatmap for metadata
      if (isolate(input$heatmap_known_factors == "None")) {
        # subset expression by markers
        all_marks <- marker_genes$All_Unique_Markers
        log_mat <- log(as.matrix(dataTables$exp_norm[all_marks,])+1)
        log_mat <- log_mat[complete.cases(log_mat),]
        
        withProgress(expr = pheatmap(log_mat, show_colnames = FALSE,
                                     show_rownames = TRUE,
                                     clustering_method = "ward.D2"),
                     message = "Visualizing marker genes, please wait")
      } else {
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
      }
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
        
        # if no coloring
        if (isolate(input$Dim_Color == "No Color")) {
          plot_ly(dataTables$dim_orig, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = "")) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = ""))
        } else {
          # factor for coloring
          id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
          fac_int <- as.factor(dataTables$meta_df[, id_fac])
          plot_ly(dataTables$dim_orig, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = ""),
                  color = ~fac_int) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = "")) 
        }
      })
      
      output$Dim_Plot_Markers <- renderPlotly({
        valid_load()
        valid_iasva()
        valid_dim()
        
        if (isolate(input$Dim_Color == "No Color")) {
          plot_ly(dataTables$dim_mark, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = "")) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), "; ", dataTables$chosen_svs, ")", sep = ""))
        } else {
          id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
          fac_int <- as.factor(dataTables$meta_df[, id_fac])
          plot_ly(dataTables$dim_mark, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = ""),
                  color = ~fac_int) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), "; ", dataTables$chosen_svs, ")", sep = ""))
        }
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
        
        # if no coloring 
        if (isolate(input$Dim_Color == "No Color")) {
          plot_ly(dataTables$dim_orig, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = "")) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = ""))  
        } else {
          # factor for coloring
          id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
          fac_int <- as.factor(dataTables$meta_df[, id_fac])
          plot_ly(dataTables$dim_orig, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = ""),
                  color = ~fac_int) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = "")) 
        }
      })
      
      output$Dim_Plot_Markers <- renderPlotly({
        valid_load()
        valid_iasva()
        valid_dim()
        
        if (isolate(input$Dim_Color == "No Color")) {
          plot_ly(dataTables$dim_mark, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = "")) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), "; ", dataTables$chosen_svs, ")", sep = ""))  
        } else {
          id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
          fac_int <- as.factor(dataTables$meta_df[, id_fac])
          plot_ly(dataTables$dim_mark, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = ""),
                  color = ~fac_int) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), "; ", dataTables$chosen_svs, ")", sep = "")) 
        }
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
        
        # if no color
        if (isolate(input$Dim_Color == "No Color")) {
          plot_ly(dataTables$dim_orig, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = "")) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = ""))
        } else {
          # factor for coloring
          id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
          fac_int <- as.factor(dataTables$meta_df[, id_fac])
          plot_ly(dataTables$dim_orig, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = ""),
                  color = ~fac_int) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = ""))
        }
      })
      
      output$Dim_Plot_Markers <- renderPlotly({
        valid_load()
        valid_iasva()
        valid_dim()
        
        # if no color
        if (isolate(input$Dim_Color == "No Color")) {
          plot_ly(dataTables$dim_mark, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = "")) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), "; ", dataTables$chosen_svs, ")", sep = ""))
        } else {
          id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
          fac_int <- as.factor(dataTables$meta_df[, id_fac])
          plot_ly(dataTables$dim_mark, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = ""),
                  color = ~fac_int) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), "; ", dataTables$chosen_svs, ")", sep = "")) 
        }
      })
    } else if (isolate(input$Dim_Type) == "t-SNE") {
      # make interactive plots
      output$Dim_Plot_Orig <- renderPlotly({
        valid_load()
        valid_iasva()
        valid_dim()
        
        # if no coloring
        if (isolate(input$Dim_Color == "No Color")) {
          plot_ly(dataTables$dim_orig, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = "")) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = ""))
        } else {
          # factor for coloring
          id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
          fac_int <- as.factor(dataTables$meta_df[, id_fac])
          plot_ly(dataTables$dim_orig, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = ""),
                  color = ~fac_int) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = "")) 
        }
      })
      
      output$Dim_Plot_Markers <- renderPlotly({
        valid_load()
        valid_iasva()
        valid_dim()
        
        if (isolate(input$Dim_Color == "No Color")) {
          plot_ly(dataTables$dim_mark, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = "")) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), "; ", dataTables$chosen_svs, ")", sep = ""))
        } else {
          # factor for coloring
          id_fac <- which(colnames(dataTables$meta_df) == isolate(input$Dim_Color))
          fac_int <- as.factor(dataTables$meta_df[, id_fac])
          plot_ly(dataTables$dim_mark, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                  mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = ""),
                  color = ~fac_int) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), "; ", dataTables$chosen_svs, ")", sep = "")) 
        }
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
  
  # download qc histogram
  output$Download_QC <- downloadHandler(
    filename = function() { paste("IASVA", input$norm_method, "normalized", input$known_factors, "adjusted.QC.histogram.pdf", sep=".") },
    content = function(file) {
      pdf(file)
      shiny::validate(
        need(isolate(dataTables$exp_norm != ""), "Error: No gene/ADT file loaded")
      )
        # call reactive data
        summ <- summary(dataTables$detect_num)
        hist(dataTables$detect_num, col = "dodgerblue", main="", 
             ylab = "Samples (n)", xlab = "Number of genes detected in each sample")
        legend("topright", legend = paste(names(summ), round(summ, digits = 2), sep = " "), title = "Summary of genes detected")
      dev.off()
    })
  
  
  # download figure panel
  output$Download_SVs <- downloadHandler(
    filename = function() { paste("IASVA", input$norm_method, "normalized", input$known_factors, "adjusted.SV.plot.pdf", sep=".") },
    content = function(file) {
      pdf(file)
      par(mar = c(5, 5, 4, 2), xpd = TRUE)
      valid_load()
      valid_iasva()
      valid_pair_plot()
    
      if (isolate(input$All_SV_Color == "No Color")) {
        pairs(dataTables$iasva.res$sv[, c(as.numeric(input$All_SV_Num))], main="", pch=20, cex=0.5,
              lower.panel = NULL)
      } else {
        id_fac <- which(colnames(dataTables$meta_df) == input$All_SV_Color)
        fac_int <- as.factor(dataTables$meta_df[, id_fac])
        pairs(dataTables$iasva.res$sv[, c(as.numeric(input$All_SV_Num))], main="", pch=20, cex=0.5,
              col = color.vec[fac_int], bg = color.vec[fac_int], lower.panel = NULL)
        legend("bottomleft", levels(fac_int), fill=color.vec,
               title = as.character(colnames(dataTables$meta_df)[id_fac]))
      }
      dev.off()
    })
  
  # download SV coordinates table for samples
  output$Download_SVs_Table <- downloadHandler(
    filename = function() {
      paste("IASVA", input$norm_method, "normalized", input$known_factors, "SV.coordinates.csv", sep=".")
    },
    content = function(file) {
      valid_load()
      valid_iasva()
      meta_sel <- dataTables$meta_df
      for (jcol in 1:ncol(meta_sel)) {
        meta_sel[,jcol] <- as.numeric(as.factor(meta_sel[,jcol]))
      }
      sv_vals <- dataTables$iasva.res$sv[, c(as.numeric(input$All_SV_Num))]
      # need to append column names to matrix
      colnames(sv_vals) <- paste("SV", c(as.numeric(input$All_SV_Num)), sep = "")
      rownames(sv_vals) <- rownames(meta_sel)
      write.csv(sv_vals, file, row.names = TRUE)
    }
  )
  
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
      col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
      corrplot(abs(cor(iasva_vars)), type = "upper", method = "color",
               col = col(200), number.cex = 1,
               addCoef.col = "black",
               tl.col = "black", tl.srt = 90, diag = FALSE)
      dev.off()
    })
  
  # download correlation table values
  output$Download_Correlation_Table <- downloadHandler(
    filename = function() {
      paste("IASVA", input$norm_method, "normalized", input$known_factors, "correlation.values.csv", sep=".")
    },
    content = function(file) {
      valid_load()
      valid_iasva()
      # change factors to numeric for correlation
      meta_sel <- dataTables$meta_df
      for (jcol in 1:ncol(meta_sel)) {
        meta_sel[,jcol] <- as.numeric(as.factor(meta_sel[,jcol]))
      }
      iasva_vars <- cbind(dataTables$iasva.res$sv[, c(as.numeric(input$All_SV_Num))], meta_sel)
      # need to append column names to matrix
      colnames(iasva_vars) <- c(paste("SV", c(as.numeric(input$All_SV_Num)), sep = ""),
                                colnames(dataTables$meta_df))
      write.csv(abs(cor(iasva_vars)), file, row.names = FALSE)
    }
  )
  
  # download pairwise plots
  output$Download_Pairwise <- downloadHandler(
    filename = function() { paste("IASVA", input$norm_method, "normalized", input$known_factors, "adjusted.pair.svs.html", sep=".") },
    content = function(file) {
      valid_load()
      valid_iasva()
      
      if (isolate(input$Pair_SV_Color == "No Color")) {
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
                mode = "markers", text = paste("Cell ID: ", rownames(dataTables$meta_df), sep = "")) %>%
          layout(xaxis = x_ax, yaxis = y_ax)
      } else {
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
        
        py <- plot_ly(as.data.frame(dataTables$iasva.res$sv), x = dataTables$iasva.res$sv[,id_svx],
                y = dataTables$iasva.res$sv[,id_svy], type = "scatter",
                mode = "markers", text = paste("Cell ID: ", rownames(dataTables$meta_df), sep = ""),
                color = ~fac_int) %>%
          layout(xaxis = x_ax, yaxis = y_ax)
      }
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
      
      # if no coloring on heatmap for metadata
      if (isolate(input$heatmap_known_factors == "None")) {
        # subset expression by markers
        all_marks <- dataTables$markers$All_Unique_Markers
        log_mat <- log(as.matrix(dataTables$exp_norm[all_marks,])+1)
        log_mat <- log_mat[complete.cases(log_mat),]
        pheatmap(log_mat, show_colnames = FALSE,
                                     show_rownames = TRUE,
                                     clustering_method = "ward.D2")
      } else {
        # specify which metadata to plot on heatmap
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
      }
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
        
        # if no coloring
        if (isolate(input$Dim_Color == "No Color")) {
          # make interactive plots
          py <- plot_ly(dataTables$dim_orig, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                        mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = "")) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = ""))
          htmlwidgets::saveWidget(py, file)  
        } else {
          # factor for coloring
          id_fac <- which(colnames(dataTables$meta_df) == input$Dim_Color)
          fac_int <- as.factor(dataTables$meta_df[, id_fac])
          # make interactive plots
          py <- plot_ly(dataTables$dim_orig, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                        mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = ""),
                        color = ~fac_int) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = ""))
          htmlwidgets::saveWidget(py, file)  
        }
        
      } else if (input$Dim_Type == "t-SNE") {
        
        #if no coloring
        if (isolate(input$Dim_Color == "No Color")) {
          py <- plot_ly(dataTables$dim_orig, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                        mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_orig), sep = "")) %>% layout(title = paste("All Genes (n = ", nrow(dataTables$exp_norm), ")", sep = ""))
          htmlwidgets::saveWidget(py, file)
        } else {
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
    }
  )
  output$download_dim_marker <- downloadHandler(
    filename = function() { paste("IASVA", input$norm_method, "normalized", input$known_factors, "adjusted", input$Dim_Type, "marker.genes.html", sep=".") },
        content = function(file) {
          valid_load()
          valid_iasva()
          valid_dim()
          if (input$Dim_Type == "PCA") {
            # if no color
            if (isolate(input$Dim_Color == "No Color")) {
              py <- plot_ly(dataTables$dim_mark, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                            mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = "")) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), "; ", dataTables$chosen_svs, ")", sep = ""))
              htmlwidgets::saveWidget(py, file)
            } else {
              # factor for coloring
              id_fac <- which(colnames(dataTables$meta_df) == input$Dim_Color)
              fac_int <- as.factor(dataTables$meta_df[, id_fac])
              py <- plot_ly(dataTables$dim_mark, x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d",
                            mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = ""),
                            color = ~fac_int) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), "; ", dataTables$chosen_svs, ")", sep = ""))
              htmlwidgets::saveWidget(py, file)
            }
            
          } else if (input$Dim_Type == "t-SNE") {
            if (isolate(input$Dim_Color == "No Color")) {
              py <- plot_ly(dataTables$dim_mark, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                            mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = "")) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), "; ", dataTables$chosen_svs, ")", sep = ""))
              htmlwidgets::saveWidget(py, file)
            } else {
              id_fac <- which(colnames(dataTables$meta_df) == input$Dim_Color)
              fac_int <- as.factor(dataTables$meta_df[, id_fac])
              py <- plot_ly(dataTables$dim_mark, x = ~tSNE1, y = ~tSNE2, z = ~tSNE3, type = "scatter3d",
                            mode = "markers", text = paste("Cell ID: ", rownames(dataTables$dim_mark), sep = ""),
                            color = ~fac_int) %>% layout(title = paste("IA-SVA Genes (n = ", nrow(dataTables$exp_norm[dataTables$markers_formatted[,1],]), "; ", dataTables$chosen_svs, ")", sep = ""))
              htmlwidgets::saveWidget(py, file)
            }
          }
        }
    )
  
  # download dimension reduction coordinates
  output$download_all_coordinates <- downloadHandler(
    filename = function() { paste("IASVA", input$norm_method, "normalized", input$known_factors, "adjusted", input$Dim_Type, "all.genes.coordinates.csv", sep=".") },
    content = function(file) {
      valid_load()
      valid_iasva()
      valid_dim()
      write.csv(as.data.frame(dataTables$dim_orig), file, row.names = TRUE)
    }
  )
  
  output$download_marker_coordinates <- downloadHandler(
    filename = function() { paste("IASVA", input$norm_method, "normalized", input$known_factors, "adjusted", input$Dim_Type, "marker.genes.coordinates.csv", sep=".") },
    content = function(file) {
      valid_load()
      valid_iasva()
      valid_dim()
      write.csv(as.data.frame(dataTables$dim_mark), file, row.names = TRUE)
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
