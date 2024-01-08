#' Run STEGO application.
#' @name STEGO
#' @importFrom stats chisq.test
#' @import ClusTCR2
#' @export



runSTEGO <- function(){
  require(sysfonts)
  # font ------
  font <- as.data.frame(font_families())
  font
  names(font) <- "Fonts"

  ### packages ------
  suppressWarnings(require("DescTools"))
  suppressMessages(require("bslib"))
  suppressMessages(require("circlize")) # colorRamp2
  suppressMessages(suppressWarnings(require("ClusTCR2")))
  suppressMessages(require("colourpicker")) # select visual colour
  suppressMessages(require("ComplexHeatmap"))
  suppressMessages(require("corrplot"))
  suppressMessages(require("doParallel"))
  suppressMessages(require("dplyr"))
  suppressMessages(require("dior"))
  suppressMessages(require("DT"))
  suppressMessages(require("forcats"))
  suppressMessages(require("foreach"))
  suppressMessages(require("fpc")) #
  suppressMessages(require("GGally"))
  suppressMessages(require("ggnet"))
  suppressMessages(require("ggpattern"))
  suppressMessages(require("ggplot2"))
  suppressMessages(require("ggrepel"))
  suppressMessages(require("ggridges")) # ridges plot
  suppressMessages(require("ggseqlogo")) # create logo figure
  suppressMessages(require("gridExtra"))
  suppressMessages(require("harmony"))
  suppressMessages(require("igraph"))
  suppressMessages(require("iterators"))
  suppressMessages(require("linkcomm")) # create the network graphs.
  suppressMessages(require("lubridate")) # create the network graphs.
  suppressMessages(require("Matrix"))
  suppressMessages(require("network"))
  suppressMessages(require("plyr"))
  suppressMessages(require("purrr"))
  suppressMessages(require("randomcoloR"))
  suppressMessages(require("RColorBrewer"))
  suppressMessages(require("Rcpp"))
  suppressMessages(require("readr"))
  suppressMessages(require("reshape2")) # acast function
  suppressMessages(require("scales")) # acast function
  suppressMessages(require("scGate"))
  suppressMessages(require("Seurat"))
  suppressMessages(require("SeuratObject"))
  suppressMessages(require("SeuratData"))
  suppressMessages(require("SeuratDisk"))
  suppressMessages(require("shiny"))
  suppressMessages(require("shinyBS"))
  suppressMessages(require("shinybusy"))
  suppressMessages(require("shinyWidgets"))
  suppressMessages(require("showtext"))
  suppressMessages(require("showtextdb"))
  suppressMessages(require("stringr"))
  suppressMessages(require("sysfonts"))
  suppressMessages(require("tibble"))
  suppressMessages(require("tidyr"))
  suppressMessages(require("VLF"))
  require("chisq.posthoc.test")

  suppressWarnings(source(system.file("scGATE","custom_df_scGATE.R",package = "STEGO.R")))
  suppressWarnings(source(system.file("scGATE","custom_annotation_models.R",package = "STEGO.R")))

  col_markers <- c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu","YlOrBr", "YlOrRd")
  options(shiny.maxRequestSize = 100000*1024^2)

  # error messages------
  error_message_val1 <- "Upload call tags file"
  error_message_val2 <- "Upload TCR file"
  error_message_val3 <- "Upload count file"
  error_message_val4 <- "Check for missing uploaded files"
  error_message_val_10x_barcode <- "Upload 10x Barcode file"
  error_message_val_10x_features <- "Upload 10x Features file"
  error_message_val_sc <- "Upload raw processed count matrix for Seurat"
  error_message_val_UMAP <- "Upload .h5Seurat file"

  # cellTypist -----
  cellTypistModels <- c("Immune_All_Low.pkl",   #immune sub-populations combined from 20 tissues of 18 studies
                        "Immune_All_High.pkl" , # immune populations combined from 20 tissues of 18 studies
                        "Adult_Mouse_Gut.pkl"  , #cell types in the adult mouse gut combined from eight datasets
                        "Autopsy_COVID19_Lung.pkl" , # cell types from the lungs of 16 SARS-CoV-2 infected COVID-19 autopsy adult donors
                        "COVID19_Immune_Landscape.pkl" ,  #immune subtypes from lung and blood of COVID-19 patients and healthy controls
                        "Cells_Fetal_Lung.pkl"  ,# cell types from human embryonic and fetal lungs
                        "Cells_Intestinal_Tract.pkl",  # intestinal cells from fetal, pediatric (healthy and Crohn's disease) and adult human gut
                        "Cells_Lung_Airway.pkl" ,  #cell populations from scRNA-seq of five locations of the human lungs and airways
                        "Developing_Human_Brain.pkl"  ,# cell types from the first-trimester developing human brain
                        "Developing_Human_Thymus.pkl" ,# cell populations in embryonic, fetal, pediatric, and adult stages of the human thymus
                        "Developing_Mouse_Brain.pkl"  ,# cell types from the embryonic mouse brain between gastrulation and birth
                        "Healthy_COVID19_PBMC.pkl" ,#  peripheral blood mononuclear cell types from healthy and COVID-19 individuals
                        "Human_IPF_Lung.pkl" ,#  cell types from idiopathic pulmonary fibrosis, chronic obstructive pulmonary disease, and healthy lungs of adult humans
                        "Human_Lung_Atlas.pkl" ,#  integrated Human Lung Cell Atlas (HLCA) combining 46 datasets of the human respiratory system
                        "Human_PF_Lung.pkl"  ,# cell types from different forms of pulmonary fibrosis lungs of adult humans
                        "Lethal_COVID19_Lung.pkl" ,#  cell types from the lungs of individuals who died of COVID-19 and control individuals
                        "Nuclei_Lung_Airway.pkl",#   cell populations from snRNA-seq of five locations of the human lungs and airways
                        "Pan_Fetal_Human.pkl"  # stromal and immune populations from the human fetus
  )

  # functions -----
  test_fun <- function() {
    for (i in 1:15) {
      incProgress(1/15)
      sum(runif(1000000,0,1))
    }
  }

  gg_fill_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  PercentageFeatureSet2 <- function (sc = sc, assay = 'RNA', pattern = "Mt", col.name = "mtDNA") {
    features <- grep(pattern = pattern, x = rownames(sc[[assay]][]), value = TRUE)
    features
    counts_of_selected <- colSums(x = GetAssayData(object = sc,
                                                   assay = assay, slot = "counts")[features, , drop = FALSE])
    counts_of_selected
    total_of_object <- colSums(x = GetAssayData(object = sc,
                                                assay = assay, slot = "counts"))
    percent_of_set <- counts_of_selected/total_of_object*100
    percent_of_set
    col.name2 = paste0(col.name,"_percent")
    col.name2
    sc <- AddMetaData(object = sc, metadata = percent_of_set,
                      col.name = col.name2)
    #
    sc
  }

  TCR_Expanded_fun <- function (sc,Samp_col,V_gene_sc) {
    sc <- sc
    req(V_gene_sc, Samp_col)

    df3.meta <- sc@meta.data
    df3.meta2 <- df3.meta[,names(df3.meta) %in% c(Samp_col,V_gene_sc)]
    names(df3.meta2)[names(df3.meta2) %in% Samp_col] <- "ID_Column"
    names(df3.meta2)[names(df3.meta2) %in% V_gene_sc] <- "v_gene_selected"
    # }
    df3.meta3 <- df3.meta2
    df3.meta3$v_gene_selected <- ifelse(df3.meta3$v_gene_selected=="_._","Unknown",df3.meta3$v_gene_selected)
    df3.meta3$v_gene_selected <- ifelse(df3.meta3$v_gene_selected=="NA_NA & NA_NA","Unknown",df3.meta3$v_gene_selected)
    df3.meta3$v_gene_selected <- ifelse(df3.meta3$v_gene_selected=="","Unknown",df3.meta3$v_gene_selected)
    df3.meta3$v_gene_selected <- ifelse(df3.meta3$v_gene_selected=="NA","Unknown",df3.meta3$v_gene_selected)

    if (nrow(df3.meta3[-c(grep("Unknown",df3.meta3$v_gene_selected )),])>0) {
      df3.meta3 <- df3.meta3[-c(grep("Unknown",df3.meta3$v_gene_selected )),]
    }

    meta2.names <- names(df3.meta3)
    df3.meta3$samp.count <- 1
    total.condition <- as.data.frame(ddply(df3.meta3,"ID_Column",numcolwise(sum)))
    emtpy <- matrix(nrow =dim(df3.meta3)[1],ncol=dim(total.condition)[1])

    for (i in 1:dim(df3.meta3)[1]) {
      emtpy[i,] <- ifelse(df3.meta3$ID_Column[i]==total.condition$ID_Column[1:dim(total.condition)[1]],
                          total.condition[total.condition$ID_Column==total.condition$ID_Column[1:dim(total.condition)[1]],2],F)
    }
    as.data.frame(emtpy)
    #
    df3.meta3$frequency <- 1/rowSums(emtpy)
    df3.meta3$percent <- 1/rowSums(emtpy)*100

    df3 <- as.data.frame(ddply(df3.meta3,meta2.names,numcolwise(sum)))
    df3 <- df3[order(df3$samp.count,decreasing = T),]

    df4 <- df3 %>%
      mutate(Frequency_expanded = case_when(
        frequency <=1e-4 ~ "1. Rare (0 > X < 1e-4)",
        frequency <= 0.001 ~ "2. Small (1e-4 > X <= 0.001)",
        frequency <= 0.01 ~ "3. Medium (0.001 > X <= 0.01)",
        frequency <= 0.10 ~ "4. Large (0.01 > X <= 0.1)",
        frequency <= 0.50 ~ "5. Gigantic (0.1 > X <= 0.5)",
        frequency <= 1 ~ "6. Hyperexpanded (0.5 > X <= 1)",
        TRUE ~ "Other"))

    df4 <- df4 %>%
      mutate(Number_expanded = case_when(
        samp.count <=1 ~ "1. Single (0 < X <= 1)",
        samp.count <=5 ~ "2. Small (1 < X <= 5)",
        samp.count <=20 ~ "3. Medium (5 < X <= 20)",
        samp.count <=100 ~ "4. Large (20 < X <= 100)",
        samp.count <= 500 ~ "5. Hyperexpanded (100 < X <= 500)",
        TRUE ~ "6. Hyperexpanded (>500)"))
    df4
  }
  mat_sum <- function(sc,Samp_col,V_gene_sc) {
    sc <- sc
    df <- sc@meta.data
    df <- as.data.frame(df)
    names(df)
    unique.df <- unique(df[,names(df) %in% c(Samp_col,V_gene_sc) ])
    names(unique.df) <- c("group","chain")
    unique.df <- subset(unique.df,unique.df$chain != "NA")
    unique.df <- subset(unique.df,unique.df$group != "NA")
    unique.df$cloneCount <- 1
    mat <- acast(unique.df, chain~group, value.var="cloneCount")
    mat[is.na(mat)] <- 0
    Count_data <- as.data.frame(rowSums(mat))
    names(Count_data) <- "V1"
    unique.df <- (df[,names(df) %in% c(Samp_col,V_gene_sc) ])
    names(unique.df) <- c("group","chain")
    unique.df <- subset(unique.df,unique.df$chain != "NA")
    unique.df <- subset(unique.df,unique.df$group != "NA")
    unique.df$cloneCount <- 1
    mat <- acast(unique.df, chain~group, value.var="cloneCount")
    mat[is.na(mat)] <- 0
    sum_data <- as.data.frame(rowSums(mat))
    names(sum_data) <- "V1"
    mat <- as.data.frame(mat)
    mat$No.TimePoints <-Count_data$V1
    mat$CloneTotal <-sum_data$V1
    mat
  }




  # UI page -----
  ui <- fluidPage(
    # shinythemes::themeSelector(),
    tags$style(HTML("
                    .tabbable > .nav > li > a {background-color: #808080;  color:white}
                    ")
    ),
    navbarPage(title = "STEGO.R",
               # theme = "cerulean",  # <--- To use a theme, uncomment this
               # "shinythemes",
               # theme=bs_theme(version = 5, bootswatch = "default"),

               navbarMenu("Quality control",
                          # tabPanel("Convert_to_AIRR_format",
                          #
                          #
                          #          ),
                          # 10x_Genomics ----

                          tabPanel("10x_Genomics",

                                   sidebarLayout(
                                     sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                                  # selectInput("dataset_10x", "Choose a dataset:", choices = c("test_data_10x", "own_data_10x")),
                                                  selectInput("Source_type_10x","Input types",choices = c("Raw",".h5")),
                                                  conditionalPanel(condition="input.Source_type_10x=='Raw'",
                                                                   fileInput('file_calls_10x', 'Barcode file (.tsv.gz or .tsv)'),
                                                                   fileInput('file_features_10x', 'Features file (.tsv.gz or .tsv)'),
                                                                   fileInput('file_matrix_10x', 'Matrix file (.mtx.gz or .mtx)'),

                                                  ),
                                                  selectInput("csv_contig_file","format of the contig file",choices  = c("csv/csv.gz","tsv")),
                                                  fileInput('file_TCR_10x', 'filtered contig annotations'),
                                                  # textInput("sample_name_10x","Add file and sample name","Treatment_group"),
                                                  h6("File name"),
                                                  fluidRow(
                                                    column(6,textInput("group10x","Add Treatment/Group name","Group")),
                                                    column(6,textInput("Indiv10x","Add Individual name","Indiv"))
                                                  ),
                                                  selectInput("BCR_TCR_10x","Type of data",choices = c("TCR only","BCR only")),
                                     ),
                                     # 10x main panel -----
                                     mainPanel(

                                       tabsetPanel(id = "panel_10x",
                                                   tabPanel("Uploaded data",value = 1,
                                                            conditionalPanel(condition="input.Source_type_10x=='Raw'",
                                                                             div(DT::dataTableOutput("test.files.10x1")),
                                                                             div(DT::dataTableOutput("test.files.10x2")),
                                                            ),

                                                            # div(DT::dataTableOutput("test.files.10x3")),
                                                            div(DT::dataTableOutput("test.files.10x4")),
                                                   ),
                                                   tabPanel("TCRex",
                                                            add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                            div(DT::dataTableOutput("tb_TCRex_10x_df")),
                                                            downloadButton('downloaddf_TCRex_10x','Download table')

                                                   ),
                                                   tabPanel("Seurat QC",value = 2,
                                                            tags$head(tags$style("#tb_10x_matrix2  {white-space: nowrap;  }")),
                                                            add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                            div(DT::dataTableOutput("tb_10x_matrix2")),
                                                            tags$head(tags$style("#sum_tb_10x1  {white-space: nowrap;  }")),
                                                            div(DT::dataTableOutput("sum_tb_10x1")),
                                                            div(DT::dataTableOutput("tb_10x_meta1")),
                                                            fluidRow(column(3,downloadButton('downloadtb_10x_matrix2','Download matrix')),
                                                                     column(3,downloadButton('downloadtb_10x_metadata2','Download metadata'))
                                                            )

                                                   ),
                                                   tabPanel("ClusTCR",value = 3,
                                                            tags$head(tags$style("#tb_10x_contigues1  {white-space: nowrap;  }")),
                                                            div(DT::dataTableOutput("tb_10x_contigues1")),
                                                            selectInput("chain_clusTCR2_10x","Select to download",choices = c("AG","BD","IgH","IgLK")),
                                                            downloadButton('downloadtb_10x_contigues1','Download clusTCR')
                                                   ),
                                                   tabPanel("TCR_Explore",value = 4,
                                                            tags$head(tags$style("#dt_TCR_Explore_10x  {white-space: nowrap;  }")),
                                                            div(DT::dataTableOutput("dt_TCR_Explore_10x")),
                                                            downloadButton('downloaddt_TCR_Explore_10x','Download TCR_Explore')

                                                   ),
                                       ),
                                     )
                                   )
                          ),
                          # 10x_Genomics end -----
                          # BD Rhapsody  ------
                          tabPanel("BD rhapsody data",
                                   sidebarLayout(
                                     sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 1000px; position:relative;", width=3,
                                                  # UPLOAD the three files...
                                                  # selectInput("dataset_BD", "Choose a dataset:", choices = c("test_data_BD", "own_data_BD")),
                                                  textInput("name.BD",h5("Add File Name"),value = ""),
                                                  fluidRow(
                                                    column(6,  numericInput("no_lines_skip_Tags","If needed, skip first 7 lines of Sample Tag",value = 7,min=0,max=20,step=7),),
                                                    column(6, fileInput('file_calls_BD', 'Sample Tag Calls (.csv)',
                                                                        accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),)
                                                  ),


                                                  h5("Upload the matrix files"),
                                                  selectInput("Format_bd","Format type",choices = c("cellXgene","Barcode_features_matrix"),selected = "Barcode_features_matrix"),
                                                  conditionalPanel(condition="input.Format_bd=='Barcode_features_matrix'",
                                                                   fileInput('file_barcode_bd', 'Barcode file (.tsv.gz or .tsv)'),
                                                                   fileInput('file_matrix_bd', 'Matrix file'),
                                                                   fileInput('file_features_bd', 'Features file (.tsv.gz or .tsv)'),

                                                  ),
                                                  conditionalPanel(condition="input.Format_bd=='cellXgene'",
                                                                   fluidRow(
                                                                     column(6, numericInput("no_lines_skip_counts","If needed, skip first 7 lines of Count Matrix",value = 7,min=0,max=10,step=7),),
                                                                     column(6,  fileInput('file_counts_BD', 'Counts (.csv)',
                                                                                          accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),),
                                                                   ),
                                                  ),
                                                  h5("Upload the TCR Contig file"),
                                                  selectInput("filtered_list","Contig Format",choices = c("Paired","Dominant","Unfiltered")),
                                                  conditionalPanel(condition="input.filtered_list=='Paired'",
                                                                   fluidRow(
                                                                     column(6, numericInput("no_lines_skip_TCR","If needed, skip first 7 lines of VDJ Contig file",value = 7,min=0,max=10,step=7),),
                                                                     column(6,fileInput('file_TCR_BD', 'TCR file (.csv)',
                                                                                        accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))),),
                                                  ),
                                                  conditionalPanel(condition="input.filtered_list=='Dominant' || input.filtered_list=='Unfiltered'",
                                                                   fileInput('file_TCR_bd2', 'Contig AIRR file (.tsv)',
                                                                             accept=c('.tsv','tsv'))
                                                  ),
                                                  ### filter out non-function TCR and un-paired TCR
                                                  conditionalPanel(condition="input.filtered_list=='Dominant' || input.filtered_list=='Unfiltered'",
                                                                   column(6,selectInput("locus_column",h5("Chain (e.g. locus)"),"")),
                                                  ),
                                                  conditionalPanel(condition="input.filtered_list=='Paired'",
                                                                   fluidRow(
                                                                     column(6, selectInput("V_gene_AG_BDrhap",h6("Alpha/Gamma V gene"),""),),
                                                                     column(6, selectInput("Junction_AG_BDrhap",h6("Alpha/Gamma junction"),""),),
                                                                     column(6, selectInput("V_gene_BD_BDrhap",h6("Beta/Delta V gene"),"")),
                                                                     column(6, selectInput("Junction_BD_BDrhap",h6("Beta/Delta junction"),"") )
                                                                   ),
                                                  ),
                                                  fluidRow(
                                                    column(6, checkboxInput("filtering_TCR", "paired chains?", value = FALSE, width = NULL),),
                                                    column(6, checkboxInput("BCR_present", "BCR present?", value = FALSE, width = NULL),),
                                                  ),
                                     ),
                                     # main panel ------
                                     mainPanel(
                                       tabsetPanel(
                                         tabPanel("Imported data",
                                                  add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                  div(DT::dataTableOutput("test.files")),
                                                  conditionalPanel(condition="input.Format_bd=='cellXgene'",
                                                                   add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                   div(DT::dataTableOutput("test.files3"))
                                                  ),
                                                  conditionalPanel(condition="input.Format_bd=='Barcode_features_matrix'",
                                                                   add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                   div(DT::dataTableOutput("test.files.bd1")),
                                                                   add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                   div(DT::dataTableOutput("test.files.bd2")),
                                                                   add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                   div(DT::dataTableOutput("test.files.bd3")),
                                                  ),
                                                  conditionalPanel(condition="input.filtered_list=='Paired'",
                                                                   add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                   div(DT::dataTableOutput("test.files2"))
                                                  ),
                                                  conditionalPanel(condition="input.filtered_list=='Dominant' || input.filtered_list=='Unfiltered'",
                                                                   add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                   div(DT::dataTableOutput("test.files.bd4")),
                                                  ),
                                         ),
                                         # tabPanel("Checking Merge",
                                         #
                                         #          div(DT::dataTableOutput("Check_table")),
                                         #
                                         #          ),
                                         tabPanel("clusTCR2",
                                                  tags$head(tags$style("#tb_clusTCR  {white-space: nowrap;  }")),
                                                  add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                  div(DT::dataTableOutput("tb_clusTCR")),
                                                  selectInput("chain_clusTCR2_bd","Select to download",choices = c("AG","BD","IgH","IgLK")),
                                                  downloadButton('downloaddf_clusTCR','Download table')
                                         ),
                                         tabPanel("TCRex",
                                                  add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                  div(DT::dataTableOutput("tb_TCRex_BDrap_df")),
                                                  downloadButton('downloaddf_TCRex_BDrap','Download table')
                                         ),
                                         tabPanel("For Seurat",
                                                  tags$head(tags$style("#tb_count_matrix  {white-space: nowrap;  }")),
                                                  add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                  div(DT::dataTableOutput("tb_count_matrix")),
                                                  add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                  div(DT::dataTableOutput("tb_metadata_sc")),
                                                  fluidRow(
                                                    column(3,downloadButton('downloadtb_count_matrix','Download count table')),
                                                    column(3),
                                                    column(3,downloadButton('downloadtb_metadata_sc','Download meta.data table')),
                                                  ),
                                         ),
                                         tabPanel("TCR_Explore",
                                                  tags$head(tags$style("#tb_TCR_Explore  {white-space: nowrap;  }")),
                                                  add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                  div(DT::dataTableOutput("tb_TCR_Explore")),
                                                  downloadButton('downloadtb_TCR_Explore','Download table')
                                         ),
                                         tabPanel("Multi-TCR",
                                                  div(DT::dataTableOutput("tb_multiTCR")),
                                                  downloadButton('downloadtb_multiTCR','Download Multi-TCR table')
                                         ),
                                         tabPanel("Create Sample Tags file",
                                                  tags$head(tags$style("#tb_sample_tags_created  {white-space: nowrap;  }")),
                                                  textInput("sample_tags_name","Name of sample",value = "BD EA splenocyte"),
                                                  div(DT::dataTableOutput("tb_sample_tags_created")),
                                                  downloadButton('downloadtb_sample_tags','Download Tags')
                                         ),
                                       )
                                     )

                                   ),
                          ),
                          # BD rhapsody end ----
                          # Convert from .h5Seurat to .rds --------
                          tabPanel("Convert format",
                                   sidebarLayout(
                                     sidebarPanel(
                                       fileInput("file1_h5Seurat.file",
                                                 "Choose .h5Seurat files from directory",
                                                 multiple = TRUE,
                                                 accept=c('.h5Seurat','h5Seurat')),
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       textInput("project_name5","Name of Project",value = ""),
                                       downloadButton('downloaddf_SeruatObj_Convert_to_RDS','Download converted .rds Seurat Obj')


                                     ),
                                     mainPanel(
                                       tabsetPanel(id = "Converting_formatting",
                                                   tabPanel("Converting", value = "converting_PA",
                                                            p("Convert .h5Seurat (V4 Seurat) to .rds"),
                                                            p(""),
                                                            p(""),

                                                            verbatimTextOutput('Convert_to_RDS_out')

                                                   )
                                       )
                                     )
                                   )
                          ),
                          # array format -----
                          tabPanel("Array",
                                   sidebarLayout(
                                     sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                                  fileInput('file_calls_Array', 'Matrix file (.txt.gz)',
                                                            accept=c('text/csv', 'text/comma-separated-values,text/plain', 'txt.gz')),
                                                  fileInput('file_contig_Array', 'TCR file (.txt.gz)',
                                                            accept=c('text/csv', 'text/comma-separated-values,text/plain', 'txt.gz')),
                                                  checkboxInput("pairing_TCR_Array","Paired only?"),
                                                  textInput("sample_name_Array","Add sample name","Treatment_group"),
                                                  textInput("name.array","Add file name","FileName")
                                     ),
                                     mainPanel(
                                       tabsetPanel(
                                         tabPanel("Check files",
                                                  add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                  div(DT::dataTableOutput("test.files_array_Matrix")),
                                                  add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                  div(DT::dataTableOutput("test.files_array_contig")),

                                         ),
                                         tabPanel("Filtering",
                                                  add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                  div(DT::dataTableOutput("test.files_array_contig_Filtered")),
                                         ),
                                         tabPanel("clusTCR",
                                                  add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                  div(DT::dataTableOutput("test.files_ClusTCR2_array")),
                                                  downloadButton("download_ClusTCR2_labs_array"),
                                         ),
                                         tabPanel("TCRex"),
                                         tabPanel("Seurat",
                                                  add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                  div(DT::dataTableOutput("tb_array_matrix2")),
                                                  div(DT::dataTableOutput("tb_Array_meta1")),
                                                  fluidRow(column(3,downloadButton('downloadtb_array_matrix2','Download matrix')),
                                                           column(3,downloadButton('downloadtb_array_metadata2','Download metadata'))
                                                  ),





                                         ),
                                         tabPanel("TCR_Explore")


                                       )
                                     )





                                   )

                          ),

                          # array data -----

               ), # NavbarMenu
               # TCRex merge files ----
               tabPanel("TCRex merge",
                        sidebarLayout(
                          sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                       # selectInput("dataset2", "Choose a dataset:", choices = c("test_data_clusTCR2","own_data_clusTCR2")),
                                       fileInput('file2_TCRexMerge', 'Select files to merge',
                                                 multiple = TRUE,
                                                 accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                       downloadButton('downloaddf_TCRexFiltered','Download table')
                          ),
                          mainPanel(
                            tabsetPanel(id = "TCRex_tabs",
                                        tabPanel("Merge Multiple Files",value = "merge",
                                                 div(DT::dataTableOutput("DEx_TCRexFiltered")),
                                                 # downloadButton('downloaddf_multiple_ClusTCR2','Download table')
                                        ),


                            )
                          )
                        )
               ),

               # TCR clustering with ClusTCR2 -----
               tabPanel("ClusTCR2",
                        sidebarLayout(
                          sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                       # selectInput("dataset2", "Choose a dataset:", choices = c("test_data_clusTCR2","own_data_clusTCR2")),

                                       selectInput("Clust_lab_tab_output","Add prefix to file",choices = c("AG","BD","IgH","IgKL")),

                                       conditionalPanel(condition="input.clusTCR2_tabs=='merge'",
                                                        fileInput('file1_ClusTCR2_multiple', 'Select files to merge',
                                                                  multiple = TRUE,
                                                                  accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))

                                       ),

                                       conditionalPanel(condition="input.clusTCR2_tabs=='processing1' || input.clusTCR2_tabs=='processing2'",
                                                        fileInput('file2_ClusTCR2', 'Select file for single samples',
                                                                  accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),

                                                        actionButton("run_ClusTCR2","Update Clustering"),
                                                        fluidRow(
                                                          column(6,selectInput( "clusTCR2_names",label = h5("CDR3"),"")),
                                                          column(6,selectInput( "clusTCR2_Vgene",label = h5("V gene"),"")),
                                                        ),
                                                        fluidRow(
                                                          column(6,checkboxInput("allele_ClusTCR2","Remove allele *00?", value = T)),
                                                          # column(6, numericInput("cores_ClusTCR2","Number of cores to parallel",value=1))
                                                        ),
                                       ),


                                       conditionalPanel(condition="input.clusTCR2_tabs=='processing2'",

                                                        downloadButton('download_ClusTCR_labels','Download Cluster table')

                                       ),

                          ),
                          mainPanel(
                            tabsetPanel(id = "clusTCR2_tabs",
                                        tabPanel("Merge Multiple Files",value = "merge",
                                                 add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                 div(DT::dataTableOutput("DEx_multiple_ClusTCR2")),
                                                 downloadButton('downloaddf_multiple_ClusTCR2','Download table')
                                        ),
                                        tabPanel("Clustering inputs",value = "processing1",
                                                 add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                 div(DT::dataTableOutput("clust_dt2_table")),
                                        ),

                                        tabPanel("Outputs",value = "processing2",
                                                 tabsetPanel(
                                                   tabPanel("Processing",
                                                            add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                            verbatimTextOutput('ClusTCR2_Time'),
                                                            verbatimTextOutput('verbatum_ClusTCR2')
                                                   ),


                                                   # div(DT::dataTableOutput("")),
                                                   tabPanel("Table for analysis",
                                                            add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                            div(DT::dataTableOutput("ClusTCR2_lab")),

                                                            p(" "),

                                                            # verbatimTextOutput('ClusTCR2_lab'),
                                                            # div(DT::dataTableOutput("")),
                                                   ),
                                                   tabPanel("Figures",
                                                            # cluster number
                                                            fluidRow(
                                                              column(3,numericInput("selected_Cluster","Selected cluster",value = 1)),
                                                              column(3,numericInput("filter_connections","Keep connections >",value = 1)),
                                                              #Name (CDR3_V_gene_Cluster), cluster, CDR3, V_gene, Len (length of CDR3 sequence),CDR3_selected,Name_selected,cluster_selected, (_selected only prints names of the chosen cluster), None
                                                              column(3,selectInput('lab_clust_by',"Label cluster by:",choices = c("Name", "cluster", "CDR3", "V_gene", "Len","CDR3_selected","Name_selected", "cluster_selected","None"),selected = "cluster")),
                                                              column(3,selectInput('Clust_size_order',"Order of cluster",choices = c("cluster", "Original_cluster", "Clust_size_order"),selected = "Clust_size_order")),
                                                              column(3,selectInput('colour_ClusTCR2',"Type of colouring",choices = c("color_all", "color_test"),selected = "color_test")),
                                                              column(3,numericInput("text_size1","Size of selected cluster",value = 4)),
                                                              column(3,numericInput("text_size2","Size of non-selected cluster",value = 2)),
                                                            ),

                                                            fluidRow(
                                                              conditionalPanel(condition="input.colour_ClusTCR2 == 'color_all'",
                                                                               column(3,selectInput("colour_ClusTCR2_types","Colour panel",choices = c("default","rainbow","random","heat.colors","terrain.colors","topo.colors","hcl.colors"))),

                                                              ),
                                                              conditionalPanel(condition="input.colour_ClusTCR2 == 'color_test'",
                                                                               fluidRow(
                                                                                 column(3,colourInput("sel_colour_netplot","Selected colour","purple")),
                                                                                 column(3,colourInput("nonsel_colour_netplot","Non-selected colour","grey80")),
                                                                               ))
                                                            ),
                                                            fluidRow(
                                                              column(8,plotOutput("NP_ClusTCR",height="600px")),
                                                              column(4,plotOutput("MP_ClusTCR",height="200px")),
                                                            ),
                                                            fluidRow(
                                                              column(1,numericInput("width_Network_plot2", "Width of PDF", value=10)),
                                                              column(1,numericInput("height_Network_plot2", "Height of PDF", value=8)),
                                                              column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_Network_plot2','Download Network PDF')),
                                                              column(2,numericInput("width_png_Network_plot2","Width of PNG", value = 1200)),
                                                              column(2,numericInput("height_png_Network_plot2","Height of PNG", value = 1000)),
                                                              column(2,numericInput("resolution_PNG_Network_plot2","Resolution of PNG", value = 144)),
                                                              column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_Network_plot2','Download Network PNG')),

                                                            ),
                                                            fluidRow(
                                                              column(1,numericInput("width_Motif_plot2", "Width of PDF", value=10)),
                                                              column(1,numericInput("height_Motif_plot2", "Height of PDF", value=3.5)),
                                                              column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_Motif_plot2','Download Motif PDF')),
                                                              column(2,numericInput("width_png_Motif_plot2","Width of PNG", value = 1200)),
                                                              column(2,numericInput("height_png_Motif_plot2","Height of PNG", value = 600)),
                                                              column(2,numericInput("resolution_PNG_Motif_plot2","Resolution of PNG", value = 144)),
                                                              column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_Motif_plot2','Download Motif PNG'))),

                                                   ),
                                                   # tabPanel("Time",

                                                   # )
                                                 ))

                            )
                          )
                        )
               ),
               # end TCR clustering ------
               # Quality control side bar panel -----
               tabPanel("Seurat QC",
                        sidebarLayout(
                          sidebarPanel(shinyjs::useShinyjs(),
                                       id = "side-panel",
                                       style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,

                                       # selectInput("dataset_sc", "Choose a dataset:", choices = c("test_data_sc", "own_data_sc")),
                                       # upload the file
                                       fileInput('file_SC', 'Load csv (BDrhap), csv.gz (10x), .h5 (10x)',

                                       ),
                                       fileInput('file_SC_meta', 'Upload file meta.data file (.csv.gz or .csv)',

                                       ),

                                       textInput("project_name","Name of sample",value = ""),
                                       # selectInput("species","Species",choices = c("human","mouse","other")),
                                       selectInput("df_seruatobj_type","Data type", choices = c("10x_Genomics (raw)","10x_Genomics (.h5)","BD Rhapsody (Human Immune panel)","BD Rhapsody (Mouse)","Array")),
                                       selectInput("stored_in_expression","Does the .h5 object has multiple part?",choices = c("no","yes")),
                                       uiOutput("feature_input"),
                                       actionButton("run","Update Violin plot"),


                                       conditionalPanel(condition = ("input.run != 0"),
                                                        fluidRow(
                                                          column(6,numericInput("dimension_heatmap.min","View heatmap dimensions.min", value = 1)),
                                                          column(6,numericInput("dimension_heatmap.max","View heatmap dimensions.max", value = 10)),
                                                          column(6,numericInput("numberofcells","Number of cells to use for heatmap", value = 500)),
                                                        ),
                                                        selectInput("method_Seurat","Method",choices = c("vst","dispersion","mvp"),selected = "vst"),

                                                        fluidRow(
                                                          column(6,numericInput("dimension_sc","Max dimensions for clustering", value = 15)),
                                                          column(6,numericInput("resolution","Resolution of clusters", value = 1)),
                                                        ),
                                                        actionButton("run_reduction","Run clustering"),
                                                        p(""),
                                                        actionButton("run_metadata","Impute metadata after clustering"),
                                       ),


                                       conditionalPanel(condition = ("input.run_metadata != 0"),

                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        downloadButton('downloaddf_SeruatObj','Download Seurat')


                                       ),
                                       tags$hr(),
                                       # actionButton("reset_input", "Reset inputs")

                          ),
                          # QC main panel -----
                          mainPanel(



                            tabsetPanel(
                              # tabPanel("Instructions"
                              #          ),
                              tabPanel("Header check",
                                       div(DT::dataTableOutput("DEx_header_name_check.dt")),
                              ),
                              tabPanel("Violin and correlation",
                                       tabsetPanel(
                                         tabPanel("Before",
                                                  add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                  plotOutput("before_plot_sc", height = "600px"),
                                                  fluidRow(
                                                    column(1,numericInput("width_before_plot_sc", "Width of PDF", value=10)),
                                                    column(1,numericInput("height_before_plot_sc", "Height of PDF", value=8)),
                                                    column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_before_plot_sc','Download Network PDF')),
                                                    column(2,numericInput("width_png_before_plot_sc","Width of PNG", value = 1200)),
                                                    column(2,numericInput("height_png_before_plot_sc","Height of PNG", value = 1000)),
                                                    column(2,numericInput("resolution_PNG_before_plot_sc","Resolution of PNG", value = 144)),
                                                    column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_before_plot_sc','Download Network PNG')),
                                                  ),
                                         ),
                                         tabPanel("After",
                                                  p("hit 'update Violin plot' to check cut-offs"),
                                                  add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                  plotOutput("after_plot_sc", height = "600px"),

                                                  fluidRow(
                                                    column(1,numericInput("width_after_plot_sc", "Width of PDF", value=10)),
                                                    column(1,numericInput("height_after_plot_sc", "Height of PDF", value=8)),
                                                    column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_after_plot_sc','Download Network PDF')),
                                                    column(2,numericInput("width_png_after_plot_sc","Width of PNG", value = 1200)),
                                                    column(2,numericInput("height_png_after_plot_sc","Height of PNG", value = 1000)),
                                                    column(2,numericInput("resolution_PNG_after_plot_sc","Resolution of PNG", value = 144)),
                                                    column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_after_plot_sc','Download Network PNG')),
                                                  ),
                                         )
                                       ),
                              ),
                              # Variable features -----
                              tabPanel("Top variable features",
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       plotOutput("plot_10_features_sc", height = "600px"),
                                       fluidRow(
                                         column(1,numericInput("width_plot_10_features_sc", "Width of PDF", value=10)),
                                         column(1,numericInput("height_plot_10_features_sc", "Height of PDF", value=8)),
                                         column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_plot_10_features_sc','Download Network PDF')),
                                         column(2,numericInput("width_png_plot_10_features_sc","Width of PNG", value = 1200)),
                                         column(2,numericInput("height_png_plot_10_features_sc","Height of PNG", value = 1000)),
                                         column(2,numericInput("resolution_PNG_plot_10_features_sc","Resolution of PNG", value = 144)),
                                         column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_plot_10_features_sc','Download PNG')),
                                       ),
                              ),
                              # Elbow and heatmap  -----
                              tabPanel("Elbow Plot",
                                       plotOutput("create_elbowPlot_sc", height = "600px"),
                                       fluidRow(
                                         column(1,numericInput("width_create_elbowPlot_sc", "Width of PDF", value=10)),
                                         column(1,numericInput("height_create_elbowPlot_sc", "Height of PDF", value=8)),
                                         column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_create_elbowPlot_sc','Download PDF')),
                                         column(2,numericInput("width_png_create_elbowPlot_sc","Width of PNG", value = 1200)),
                                         column(2,numericInput("height_png_create_elbowPlot_sc","Height of PNG", value = 1000)),
                                         column(2,numericInput("resolution_PNG_create_elbowPlot_sc","Resolution of PNG", value = 144)),
                                         column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_create_elbowPlot_sc','Download PNG')),
                                       ),

                              ),
                              tabPanel("DimHeatmap",
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),

                                       plotOutput("create_PCA_heatmap_sc", height = "1000px"),

                                       fluidRow(
                                         column(1,numericInput("width_heatmap_sc", "Width of PDF", value=10)),
                                         column(1,numericInput("height_heatmap_sc", "Height of PDF", value=8)),
                                         column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_heatmap_sc','Download Network PDF')),
                                         column(2,numericInput("width_png_heatmap_sc","Width of PNG", value = 1200)),
                                         column(2,numericInput("height_png_heatmap_sc","Height of PNG", value = 1000)),
                                         column(2,numericInput("resolution_PNG_heatmap_sc","Resolution of PNG", value = 144)),
                                         column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_heatmap_sc','Download Network PNG')),
                                       ),

                              ),
                              # tabPanel("Resolution plot"),
                              # UMAP  -----
                              tabPanel("UMAP",
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       plotOutput("create_UMAP_sc", height = "600px")
                              ), # Export a table with meta.data and expression.
                              tabPanel("Add metadata",
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       div(DT::dataTableOutput("DEx_view.meta.dt")),
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       div(DT::dataTableOutput("DEx_table_meta.data")),
                                       # selectInput("","comaprison 1")
                              ),


                            ),


                          ),
                          # end of QC -----


                        ),
               ),
               ###################
               # Merge multiple Seurat objects -----
               tabPanel("Merge SC",
                        sidebarLayout(
                          # Sidebar with a slider input
                          sidebarPanel(id = "tPanel5",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                       selectInput("sample.type.source","Species",choices = c("hs","mm")),
                                       fileInput("file1_rds.file",
                                                 "Choose .rds files from directory",
                                                 multiple = TRUE,
                                                 accept=c("rds",".rds")),
                                       textInput("project_name2","Name of Project",value = "Pro"),
                                       downloadButton('downloaddf_SeruatObj_merged','Download Merged Seurat')
                          ),

                          # Show a plot of the generated distribution
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Upload files",
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       verbatimTextOutput("testing_mult"),
                                       verbatimTextOutput("testing_mult2")
                              ),
                              tabPanel("Variable data",
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       actionButton("run_var","Run VariableFeatures"),
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       verbatimTextOutput("var_harmony_verbrose")
                              ),
                              tabPanel("Scale data",
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       actionButton("run_scale","Run Scale"),
                                       div(DT::dataTableOutput("Tb_scaling_features_for_annotation")),
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       verbatimTextOutput("scale_harmony_verbrose")
                              ),
                              tabPanel("PCA",
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       actionButton("run_PCA","Run PCA"),
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       verbatimTextOutput("PCA_harmony_verbrose")
                              ),
                              tabPanel("harmony",
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       actionButton("run_harmony","Run Harmony"),
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       verbatimTextOutput("harmony_verbrose"),
                              ),
                              tabPanel("Dimentional reduction",
                                       fluidRow(
                                         column(3,numericInput("dimension_Merged","Max number of dimensions", value = 30)),
                                         column(6,numericInput("res_merged","Resolution of clusters", value = 0.5)),
                                       ),
                                       actionButton("run_reduction_harmony","Run Dimentional reduction"),
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       verbatimTextOutput("testing_mult3"),
                                       # actionButton("download_reduction_harmony","download"),


                              ),
                              tabPanel("UMAP",
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       actionButton("download_reduction_harmony","download"),
                                       plotOutput("create_UMAP_merged",height="600px"),

                                       fluidRow(
                                         column(1,numericInput("width_sc_merged", "Width of PDF", value=10)),
                                         column(1,numericInput("height_sc_merged", "Height of PDF", value=8)),
                                         column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_sc_merged','Download Network PDF')),
                                         column(2,numericInput("width_png_sc_merged","Width of PNG", value = 1200)),
                                         column(2,numericInput("height_png_sc_merged","Height of PNG", value = 1000)),
                                         column(2,numericInput("resolution_PNG_sc_merged","Resolution of PNG", value = 144)),
                                         column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_sc_merged','Download Network PNG')),

                                       ),

                              ),

                            ))
                        )

               ),
               ###################
               # remove cells based on one factor -----
               tabPanel("Remove Samps",
                        sidebarLayout(
                          sidebarPanel(id = "tPanelSamps",style = "max-height: 800px; position:relative;", width=4,
                                       fileInput("file1_rds.fileSampsRemove",
                                                 "Upload .rds file",
                                                 multiple = F,
                                                 accept=c('.rds','rds')),
                                       textInput("project_name4","Name of Project",value = ""),
                                       selectInput("Samp_col_SampToRemove","Column name",choices = ""),
                                       downloadButton('downloaddf_SeruatObj_annotated_SampToKeep','Download .rds'),
                          ),
                          mainPanel(
                            h5("Before Samples are removed"),
                            verbatimTextOutput("Preliminary_samp_to_remove"),

                            selectInput("ID_Column_factor_SampToRemove","Order of graph",choices = "", multiple = T, width = "1400px"),

                            h5("After Samples are removed"),
                            verbatimTextOutput("Filtered_samp_to_remove"),

                          )
                        )),
               # Add annotations -----
               tabPanel("Annotations",
                        sidebarLayout(
                          sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                       uiOutput("Detect_version"),


                                       selectInput("Data_types","Source",choices = c("10x_HS","BD_HS.Immune.Panel","BD_HS.Full.Panel","10x_MM","BD_MM_Full.Panel","BD_MM_Immune.Panel",
                                                                                     "TCR-seq")),
                                       selectInput("sample.type.source.markers","Species",choices = c("hs","mm")),
                                       fileInput("file1_rds.file2",
                                                 "Choose merged or single .rds files from directory",
                                                 multiple = TRUE,
                                                 accept=c('.rds','rds')),
                                       textInput("project_name3","Name of Project",value = ""),
                                       selectInput("Require_custom_geneset","Require custom genes?", choices = c("no","yes")),
                                       uiOutput("scGate_cutoffs"),

                                       downloadButton('downloaddf_SeruatObj_annotated','Download Annotated Seurat'),


                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Upload",
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       verbatimTextOutput("testing_mult_anno")
                              ),

                              tabPanel("scGATE",

                                       # custom annotations databases -----
                                       conditionalPanel(condition="input.Require_custom_geneset == 'yes'",
                                                        fluidRow(
                                                          column(2,checkboxInput("GeneSet1_scGate","GeneSet1", value = F)),
                                                          column(2,textInput("geneset1_name","Name",value = "GeneSet1")),

                                                          column(2,checkboxInput("GeneSet2_scGate","GeneSet2", value = F)),
                                                          column(2,textInput("geneset2_name","Name",value = "GeneSet2")),

                                                          column(2,checkboxInput("GeneSet3_scGate","GeneSet3", value = F)),
                                                          column(2,textInput("geneset3_name","Name",value = "GeneSet3")),
                                                          column(2,checkboxInput("GeneSet4_scGate","GeneSet4", value = F)),
                                                          column(2,textInput("geneset4_name","Name",value = "GeneSet4")),
                                                          column(2,checkboxInput("GeneSet5_scGate","GeneSet5", value = F)),
                                                          column(2,textInput("geneset5_name","Name",value = "GeneSet5")),
                                                          column(2,checkboxInput("GeneSet6_scGate","GeneSet6", value = F)),
                                                          column(2,textInput("geneset6_name","Name",value = "GeneSet6")),
                                                          column(2,checkboxInput("GeneSet7_scGate","GeneSet7", value = F)),
                                                          column(2,textInput("geneset7_name","Name",value = "GeneSet7")),
                                                          column(2,checkboxInput("GeneSet8_scGate","GeneSet8", value = F)),
                                                          column(2,textInput("geneset8_name","Name",value = "GeneSet8")),
                                                          column(2,checkboxInput("GeneSet9_scGate","GeneSet9", value = F)),
                                                          column(2,textInput("geneset9_name","Name",value = "GeneSet9")),
                                                        )
                                       ),


                                       # human 10x annotations -----
                                       conditionalPanel(condition="input.Data_types == '10x_HS' || input.Data_types == 'BD_HS.Full.Panel'",
                                                        selectInput("GenericID_scGATE","Generic To include",choices = c("Bcell","CD4T","CD8T","CD8TIL" ,"Erythrocyte" ,"Megakaryocyte" , "MoMacDC","Myeloid","NK","PanBcell","panDC","PlasmaCell","Tcell","Tcell.alphabeta"), selected = c("Bcell","MoMacDC","NK","CD8T","CD4T","PlasmaCell"), multiple = T,width = "1200px"),
                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        fluidRow(
                                                          add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                          column(3,checkboxInput("generic_scGATE","Generic (Human)", value = F)),
                                                          column(3,checkboxInput("CD4_scGATE","CD4 T cell (Human)", value = F)),
                                                          column(3,checkboxInput("CD8_scGATE","CD8 T cell (Human)", value = F)),
                                                          column(3,checkboxInput("GeneralMarkers_scGATE","ESGA general markers (Human; Lit)", value = F)),
                                                          column(3,checkboxInput("Function_scGATE","T cell Function types (Human; Lit)", value = F)),
                                                          column(3,checkboxInput("Ex.Sen_scGATE","Exhausted/Senescence (Human)", value = F)),
                                                          column(3,checkboxInput("COVID_scGATE","COVID markers (Human)", value = F)),
                                                          column(3,checkboxInput("Activation_scGATE","Activation markers (Human)", value = F)),
                                                          column(3,checkboxInput("IFNgTNFa_scGATE","IFNg and TNFa (Human)", value = F)),
                                                          column(3,checkboxInput("GNLY.PFR1.GZMB_scGATE","GNLY.PFR1.GZMB markers (Human)", value = F)),
                                                          column(3,checkboxInput("Interlukin_scGATE","Interleukins markers (Human)", value = F))),
                                       ),

                                       # BD rhapsody human immune panel -----
                                       conditionalPanel("input.Data_types == 'BD_HS.Immune.Panel'",
                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        fluidRow(
                                                          add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                          column(3,checkboxInput("BDrhapsody_scGATE_Tcells","Main T cell markers (CD4 and CD8)", value = F)),
                                                          add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                          column(3,checkboxInput("BDrhapsody_scGATE_GNLY.PRF1.GZMB","GNLY.PRF1.GZMB", value = F)),
                                                          add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                          column(3,checkboxInput("BDrhapsody_scGATE_TNF.IFNG","TNF.IFNG", value = F)),
                                                          add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                          column(3,checkboxInput("BDrhapsody_scGATE_Effector_CD8","T_cell_types", value = F)),
                                                          add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                          column(3,checkboxInput("BDrhapsody_scGATE_Lung.Residence.markers","BD Rhapsody (Human)", value = F)),
                                                          add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                          column(3,checkboxInput("BDrhapsody_scGATE_NK_receptors","BD Rhapsody (Human)", value = F)),
                                                          add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                          column(3,checkboxInput("BDrhapsody_scGATE_Activation","BD Rhapsody (Human)", value = F)),
                                                          add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                          column(3,checkboxInput("BDrhapsody_scGATE_Other","BD Rhapsody (Human)", value = F)),
                                                        ),

                                       ),

                                       # BD rhapsody MM full panel ----
                                       conditionalPanel("input.Data_types == 'BD_MM_Full.Panel' || input.Data_types =='10x_MM'",
                                                        h5("Under development"),
                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        h6("Curated from Sharland lab"),
                                                        fluidRow(
                                                          column(3,checkboxInput("BDrhapsody_scGATE.MM.Tcell","Major T cell popualtions", value = F)),
                                                          column(3,checkboxInput("BDrhapsody_scGATE.MM.Memory","Memory", value = F)),
                                                          column(3,checkboxInput("BDrhapsody_scGATE.MM.signatures","signatures", value = F)),
                                                          column(3,checkboxInput("BDrhapsody_scGATE.MM.Innate.NK","Innate & NK", value = F)),
                                                        ),
                                                        h6("From the human BD Rhapsody searches"),
                                                        fluidRow(
                                                          column(3,checkboxInput("BDrhapsody_scGATE.MM.TNF.IFNg","TNF.IFNg", value = F)),
                                                          column(3,checkboxInput("BDrhapsody_scGATE.MM.subtypes","Subtypes", value = F)),
                                                          column(3,checkboxInput("BDrhapsody_scGATE.MM.other","other", value = F)),

                                                        ),

                                       ),
                                       # BD rhapsody MM immune panel ----
                                       conditionalPanel("input.Data_types == 'BD_MM_Immune.Panel'",
                                                        h5("Under development")

                                       ),
                                       conditionalPanel(condition="input.Require_custom_geneset == 'yes'",
                                                        verbatimTextOutput("scGATE_verbatum_GeneSet1"),
                                                        verbatimTextOutput("scGATE_verbatum_GeneSet2"),
                                                        verbatimTextOutput("scGATE_verbatum_GeneSet3"),
                                                        verbatimTextOutput("scGATE_verbatum_GeneSet4"),
                                                        verbatimTextOutput("scGATE_verbatum_GeneSet5"),
                                                        verbatimTextOutput("scGATE_verbatum_GeneSet6"),
                                                        verbatimTextOutput("scGATE_verbatum_GeneSet7"),
                                                        verbatimTextOutput("scGATE_verbatum_GeneSet8"),
                                                        verbatimTextOutput("scGATE_verbatum_GeneSet9"),

                                       ),
                                       # human 10x annotations Verbatium -----
                                       conditionalPanel(condition="input.Data_types == '10x_HS' || input.Data_types == 'BD_HS.Full.Panel'",
                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        verbatimTextOutput("scGATE_verbatum_Generic"),
                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        verbatimTextOutput("scGATE_verbatum_CD4"),
                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        verbatimTextOutput("scGATE_verbatum_CD8"),
                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        verbatimTextOutput("scGATE_verbatum_GeneralMarkers"),
                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        verbatimTextOutput("scGATE_verbatum_Function"),
                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        verbatimTextOutput("scGATE_verbatum_Ex.Sen"),
                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        verbatimTextOutput("scGATE_verbatum_COVID"),
                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        verbatimTextOutput("scGATE_verbatum_Activation"),
                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        verbatimTextOutput("scGATE_verbatum_IFNgTNFa"),
                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        verbatimTextOutput("scGATE_verbatum_GNLY.PFR1.GZMB"),
                                                        add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                        verbatimTextOutput("scGATE_verbatum_Interlukin")
                                       ),
                                       conditionalPanel("input.Data_types == 'BD_HS.Immune.Panel'",
                                                        # verbatimTextOutput("scGATE_verbatum_BDrhapsody_scGATE"),
                                       ),



                                       conditionalPanel("input.Data_types == 'BD_MM_Full.Panel' || input.Data_types =='10x_MM'",
                                                        verbatimTextOutput("scGATE_verbatum_BDrhapsody_MM.FP.Tcell"),
                                                        verbatimTextOutput("scGATE_verbatum_BDrhapsody_MM.FP.Memory"),
                                                        verbatimTextOutput("scGATE_verbatum_BDrhapsody_MM.FP.signatures"),
                                                        verbatimTextOutput("scGATE_verbatum_BDrhapsody_MM.FP.Innate.NK"),
                                                        verbatimTextOutput("scGATE_verbatum_BDrhapsody_MM.FP.TNF.IFNg"),
                                                        verbatimTextOutput("scGATE_verbatum_BDrhapsody_MM.FP.subtypes"),
                                                        verbatimTextOutput("scGATE_verbatum_BDrhapsody_MM.FP.other")
                                       ),

                              ),
                              tabPanel("UMAP check",
                                       fluidRow(
                                         column(4, plotOutput("create_custom_1", height = "600px")),
                                         column(4,plotOutput("create_custom_2", height = "600px"),),
                                         column(4, plotOutput("create_custom_3", height = "600px")),
                                         column(4,plotOutput("create_custom_4", height = "600px"),),
                                         column(4, plotOutput("create_custom_5", height = "600px")),
                                         column(4,plotOutput("create_custom_6", height = "600px"),),
                                         column(4, plotOutput("create_custom_7", height = "600px")),
                                         column(4,plotOutput("create_custom_8", height = "600px"),),
                                         column(4, plotOutput("create_custom_9", height = "600px")),
                                       )
                                       # plotOutput("create_custom_1", height = "600px"),


                              ),

                              # classification based on TCR_seq -----

                              tabPanel("TCR-seq",
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "purple"),
                                       div(DT::dataTableOutput("TCR_seq_classification_df")),
                              ),
                              # display metadata -----
                              tabPanel("Marker check",
                                       conditionalPanel(condition="input.Panel_DEX==5",
                                                        fluidRow(
                                                          column(3,numericInput("min.ptc.sc","minimum point",value = 0.25)),
                                                          column(3,numericInput("logfc.ptc.sc","Log fold change",value = 0.25)),
                                                          column(3,selectInput("normalN","Type of Differnetial expression",
                                                                               choices = c("wilcox","bimod","roc","t","negbinom","poisson","LR","MAST","DESeq2")
                                                          )),
                                                          column(3,style = "margin-top: 25px;",actionButton("run_differental.exp","run differental expression"),)
                                                        ),

                                       ),

                                       conditionalPanel(condition="input.Panel_DEX==55 || input.Panel_DEX==5",

                                                        fluidRow(
                                                          column(4,selectInput("multiple_group_sc","Include group comparison",choices=c("no","yes"))),
                                                          column(4,selectInput("meta_data_sc_clust","Cluster by",choices = "")),
                                                          column(4,selectInput("meta_data_sc_","Add group",choices = "")),

                                                        ),

                                       ),

                                       tabsetPanel(id = "Panel_DEX",
                                                   # Cluster table -----
                                                   # tabPanel("Checking files",
                                                   #          div(DT::dataTableOutput("list_of_genes")),
                                                   #
                                                   #          # verbatimTextOutput("checking_files_markers_featurePlot_sc"),
                                                   #          ),
                                                   tabPanel("Cluster differences (Feature plot)",
                                                            actionButton("run_string.data3","View Feature plot"),
                                                            fluidRow(column(12, selectInput("string.data3","column names for summary","",multiple = T, width = "1200px") )),
                                                            fluidRow(
                                                              column(2,checkboxInput("label_is_true_features","Add plot lables",value=T)),
                                                              column(2,selectInput("norm_expression_for_all","Set Maximum",choices=c("no","yes"))),
                                                              column(2,numericInput("max_norm_FP","Set maximum scale value",value = 10, step = 1, min=1)),
                                                              column(2,colourInput("lower_col_FP","Min (Colour)",value = "grey90")),
                                                              column(2,colourInput("upper_col_FP","Max (colour)",value = "Darkblue"))
                                                            ),
                                                            plotOutput("markers_featurePlot_sc", height = "600px"),

                                                            fluidRow(
                                                              column(3,numericInput("width_markers_featurePlot_sc", "Width of PDF", value=10)),
                                                              column(3,numericInput("height_markers_featurePlot_sc", "Height of PDF", value=8)),
                                                              column(3),
                                                              column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_markers_featurePlot_sc','Download PDF'))),
                                                            fluidRow(
                                                              column(3,numericInput("width_png_markers_featurePlot_sc","Width of PNG", value = 1200)),
                                                              column(3,numericInput("height_png_markers_featurePlot_sc","Height of PNG", value = 1000)),
                                                              column(3,numericInput("resolution_PNG_markers_featurePlot_sc","Resolution of PNG", value = 144)),
                                                              column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_markers_featurePlot_sc','Download PNG'))
                                                            )

                                                   ),
                                                   # differential expression within clusters ----
                                                   tabPanel("Treatment differences within clusters",value = 55,

                                                            actionButton("run_update_clust","Update comparisons"),
                                                            fluidRow(
                                                              column(4,selectInput("unique.Idents1","comaprison 1",choices = "")),
                                                              column(4,selectInput("unique.Idents2","comaprison 2",choices = "")),
                                                            ),
                                                            tabsetPanel(
                                                              tabPanel("Table",
                                                                       div(DT::dataTableOutput("DEx_table_comparison")),
                                                                       downloadButton('downloaddf_DEx_sc','Download Table (.csv)'),
                                                                       downloadButton('downloaddf_DEx_sc_ggVolcanoR','Download ggVolcanoR compatible table (.csv)')
                                                              ),
                                                              tabPanel("Plot",
                                                                       plotOutput("volc_plot_cluster", height = "600px"))
                                                            ),
                                                   ),
                                                   tabPanel("Cluster differences (All markers)", value = 5,
                                                            add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                            div(DT::dataTableOutput("DEx_table_clusters")),
                                                            downloadButton('downloaddf_DEx_table_clusters','Download Table (.csv)')
                                                   ),


                                       )

                              ),

                              # meta data table ------
                              tabPanel("Meta data table",
                                       fluidRow(
                                         # column(3,checkboxInput("add.kmeans","Add K-means classification", value = F)),
                                         # column(3,checkboxInput("add.scGATE","Add scGATE classifications", value = T))
                                       ),
                                       add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                       div(DT::dataTableOutput("DEx_table_TcellClass_scGATE")),
                              )
                            ),
                          )

                        )
               ),


               ###################
               # Analysis (UI side panel) ---------
               tabPanel("Analysis",
                        sidebarLayout(
                          sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 1000px; position:relative;", width=3,
                                       selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR"), selected = "BD_Rhapsody_AIRR"),
                                       selectInput("species_analysis","Species",choices = c("hs","mm")),
                                       selectInput("SeuratVersion","Seurat Version", choices = c("Version 4","Version 5"), selected = "Version 4"),
                                       # selectInput("STEGO_R_pro","QC processed",choices = c("STEGO_R (.rds)")), #,"Seurat (.rds)"
                                       # textInput("name.file_clust","Name added to files",value = ""),
                                       conditionalPanel(condition="input.check_up_files== 'up'",
                                                        fileInput('file_SC_pro', 'Upload seurat file',
                                                                  accept=c("rds",'.rds','rds')),

                                                        selectInput("add_additional_lables", "add additional labels", choices = c("no","yes")),
                                                        conditionalPanel(condition="input.add_additional_lables== 'yes'",
                                                                         p("The .csv file first column should be label 'ID' and match the selected column"),
                                                                         selectInput("Samp_col2","Sample column name",choices = ""),

                                                                         fileInput("file_Labels_to_add","Upload other identifiers (.csv)",
                                                                                   accept=c('.csv','csv')
                                                                         ),
                                                        ),
                                                        selectInput("Type_of_receptor","Type of receptor",choices = c("TCR","BCR"),selected = "TCR"),
                                                        conditionalPanel(condition="input.Type_of_receptor== 'TCR'",

                                                                         fileInput('file_cluster_file_AG', 'Upload AG clusTCR2 file (.csv)',
                                                                                   accept=c('.csv','csv')),
                                                                         fileInput('file_cluster_file_BD', 'Upload BD clusTCR2 file (.csv)',
                                                                                   accept=c('.csv','csv')),

                                                        ),

                                                        conditionalPanel(condition="input.Type_of_receptor== 'BCR'",

                                                                         fileInput('file_cluster_file_IgH', 'Upload IgH clusTCR2 file (.csv)',
                                                                                   accept=c('.csv','csv')),
                                                                         fileInput('file_cluster_file_IgKL', 'Upload IgKL clusTCR2 file (.csv)',
                                                                                   accept=c('.csv','csv')),
                                                        ),
                                                        numericInput("skip_TCRex_up","Skip # of lines for TCRex file",value = 7),
                                                        fileInput('upload_TCRex_file', 'Upload TCRex (.tsv)',
                                                                  accept=c('tsv','.tsv'))

                                       ),
                                       selectInput("Samp_col","Column name",choices = ""),
                                       selectInput("V_gene_sc","V gene with/without CDR3",choices = ""),
                                       selectInput("colourtype","Type of colouring",choices = c("default","rainbow","random","heat.colors","terrain.colors","topo.colors","hcl.colors","one colour")),
                                       conditionalPanel( condition="input.check_up_files != 'up'",
                                                         conditionalPanel( condition="input.Panel_TCRUMAP != 'Expanded'",
                                                                           conditionalPanel( condition="input.check_up_files != 'up2'",
                                                                                             uiOutput("classification_to_add")
                                                                           )

                                                         )
                                       ),
                                       # Need to check the colouring by, may need to reduce to 1?
                                       conditionalPanel( condition="input.check_up_files == 'up2'",
                                                         uiOutput("classification_to_add_overview")
                                       ),

                                       # expansion priority UI -------
                                       conditionalPanel( condition="input.check_up_files == 'Prior' || input.check_up_files == 'TCR_and_GEX_tb' ",
                                                         fluidRow(
                                                           column(6,selectInput("Split_group_by_","Split graph by:",choices ="")),
                                                           column(6,numericInput("cutoff.expanded","Cut off greater than", value = 0.5, step = 0.01, min = 0,max = 0.99)),
                                                           column(6,uiOutput("cut.off_expanded2")),
                                                           column(6,uiOutput("classification_to_add2")),
                                                         ),

                                                         conditionalPanel( condition="input.PriorTBMods == 'PriorClustTB' || input.Panel_TCRUMAP == 'ClusTCR2'",
                                                                           selectInput("Clusters_to_dis_PIE","Clusters to display",choices = "",multiple = F)
                                                         ),
                                                         fluidRow(

                                                           column(6,numericInput("cut.off_percent_rep","Percent of Repertoire", value = 1, step = 1, min = 1,max = 100)),
                                                           column(6,numericInput("size.dot.umap","size of UMAP dot's",value =2,step = 1, min = 1))
                                                         )
                                       ),

                                       conditionalPanel( condition="input.PriorTBMods == 'PriorRepertoireTB' || input.Panel_TCRUMAP == 'Expanded'",
                                                         fluidRow(
                                                           column(12, uiOutput("Expanded.dotplot.cutoffs")),


                                                           column(12,selectizeInput("selected_Indiv_Ex_1","Samp 1",choices = "", multiple = T)),
                                                           column(12,selectizeInput("selected_Indiv_Ex_2","Samp 2",choices = "", multiple = T)),


                                                         )
                                       ),

                                       conditionalPanel(condition="input.Panel_TCRUMAP=='ClusTCR2'",

                                                        selectInput("chain_TCR","Chains included",choices = c("TRAG","TRBD","IgH","IgKL")),
                                                        # column(4,selectInput("V_call_clust_sc","V gene",choices = "")),
                                                        # column(4,selectInput("junction_clust_sc","Junction",choices = "")
                                                        # ),
                                       ),
                                       #


                                       # Expanded stat cut-offs -----

                                       conditionalPanel( condition="input.PriorTBMods == 'PriorRepertoireTB' || input.check_up_files == 'TCR_and_GEX_tb' ",
                                                         fluidRow(
                                                           column(4,numericInput("min_point_","Min point cut off", value = 0.25)),
                                                           column(4,numericInput("LogFC_","Min LogFC cut off", value = 0.25)),
                                                           column(4,numericInput("pval.ex.filter","adj.p-val cut-off", value = 0.1)),
                                                         )),
                                       conditionalPanel(condition = "input.Panel_TCRUMAP == 'Marker'",

                                                        column(12,selectInput("col_marker_scale","Colour scale",choices = col_markers,selected = col_markers[1])),
                                       ),

                                       conditionalPanel(condition="input.check_up_files != 'up' ",

                                                        fluidRow(
                                                          column(12,selectInput("colourtype","Type of colouring",choices = c("default","rainbow","random","heat.colors","terrain.colors","topo.colors","hcl.colors","one colour"))),
                                                          column(6,colourInput("one.colour.default","One colour","grey50")),
                                                          column(6,colourInput("NA_col_analysis","NA colour","grey90"),)

                                                        ),

                                                        conditionalPanel(condition="input.check_up_files == 'up2' || input.check_up_files == 'Prior'",
                                                                         fluidRow(column(12,selectInput("Graph_type_bar","Type of graph",choices = c("Number_expanded","Frequency_expanded","Top_clonotypes")))),
                                                                         fluidRow(
                                                                           # conditionalPanel(condition="input.Graph_type_bar =='Top_clonotypes'",
                                                                           uiOutput("Top_clone_number")
                                                                           # column(6,numericInput("top_no_clonotypes","Top clonotypes per group",value = 1,step = 1, min = 0, max = 20)),
                                                                           # ),
                                                                         ),
                                                        ),
                                                        fluidRow( column(6, numericInput("wrap_row",h5("Wrap rows"), value = 3))),

                                                        h4("What individuals to include"),
                                                        fluidRow(
                                                          column(6,selectInput("Split_by_group","Include group comparison",choices=c("no","yes"))),
                                                          # column(6,conditionalPanel(condition="input.STEGO_R_pro == 'STEGO_R (.rds)'",

                                                          # conditionalPanel(condition = "input.Split_by_group == 'yes'",

                                                          # ),

                                                        ),
                                                        fluidRow(
                                                          column(6,selectInput("by_indiv_pie_epi","Display one individual?",choices = c("no","yes"))),

                                                          column(6,selectInput("selected_Indiv","Display one individual",choices = ""),)
                                                        ),

                                                        selectInput("font_type",label = h4("Type of font"),choices = font,selected = "serif"),
                                                        fluidRow(
                                                          column(6,numericInput("text_size","Size of #",value=16)),
                                                          column(6,numericInput("title.text.sizer2","Axis text size",value=30)),

                                                          column(6,numericInput("Bar_legend_size","Legend text size",value=12)),
                                                          column(6, selectInput("legend_position","Legend location",choices = c("top","bottom","left","right","none"),selected = "right")),
                                                        ),
                                                        fluidRow(
                                                          column(6, numericInput("Strip_text_size","Strip text size", value = 16)),
                                                          column(6, numericInput("anno_text_size","Annotation text size", value = 6)),


                                                        ),



                                       )
                          ),

                          # add in clustering  (why did I add this comment?) -----
                          mainPanel(

                            tabsetPanel(id = "check_up_files",
                                        tabPanel("Check files uploaded",value = 'up',
                                                 add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),


                                                 fluidRow(
                                                   column(12,  div(DT::dataTableOutput("meta.data_check_upload"))),
                                                   column(12,  div(DT::dataTableOutput("Sample_names_merging_sc"))),
                                                   # column(12,  div(DT::dataTableOutput("Tb_TCR_clonotypes.Umap"))),
                                                   column(12,   div(DT::dataTableOutput("Tb_ClusTCR_test"))),
                                                   column(12,   div(DT::dataTableOutput("Tb_tcrex_test")))

                                                 ),

                                        ),
                                        #prioritization strategy ------

                                        tabPanel("Priortisation",value = "Prior",
                                                 # sidebarLayout(
                                                 # sidebarPanel(
                                                 h5("under construction"),
                                                 # ),
                                                 # mainPanel(
                                                 tabsetPanel(id = "PriorTBMods",
                                                             tabPanel("Analysis steps", value = "ModstoConsid",

                                                                      p("There are three section"),
                                                                      p("1. Clonotype"),
                                                                      p("2. Cluster"),
                                                                      p("3. Epitope"),
                                                                      p("   \t"),
                                                                      # p("Preliminary rules for decision tree"),
                                                                      # p("1. dataset size, one or more conditions, TCR-seq coverage"),
                                                                      # p("2. Total clonal expansion phenotype vs immunodominant clones?"),
                                                                      # p("3. TCR clustering single chain"),
                                                                      # p("4. TCRex predicted TCR - distinct TCR seuqences with same predicted epitope"),
                                                                      # p("5. Time point interrogations"),
                                                                      # p("6. Further prioritize based on function: relative to whole data set, relative to specific time points,relative to a specific group and or sample"),
                                                                      # p("7. Need to create a small ontology based on genes present: FindMarker genes,Expressed genes (e.g., CD8, CD4)"),
                                                                      # p("8. Multiple chains per cell for testing  BAA (CD8 ab T cells) or ABB (NKT?) - 2 different member of a dual beta (Dex binding)")
                                                             ),


                                                             # modules of priority ------
                                                             # tabPanel("things to consider", value = "PriorRepTB",
                                                             #
                                                             #          # actionButton("Run_TCR_check_group_check","TCR/Group check"),
                                                             #
                                                             #            ),
                                                             tabPanel("Clonotype",value = "PriorRepertoireTB",
                                                                      add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "purple"),
                                                                      verbatimTextOutput("Simple_workflow_step1"),
                                                                      verbatimTextOutput("Number_of_clonotypes_to_"),
                                                                      add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "purple"),
                                                                      uiOutput("Module_case_statements"),
                                                                      div(DT::dataTableOutput("Test_table_1")),


                                                             ),
                                                             tabPanel("Cluster",value = "PriorClustTB",
                                                                      fluidRow(
                                                                        column(3,uiOutput("Default_priority_cutoffAG")),
                                                                        column(3,uiOutput("Default_priority_cutoffBD")),
                                                                      ),
                                                                      verbatimTextOutput("number_clusters_to_analyse_AG"),
                                                                      verbatimTextOutput("number_clusters_to_analyse_BD"),

                                                                      add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "purple"),
                                                                      div(DT::dataTableOutput("PriorClustTB_Tab")),
                                                                      add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "purple"),
                                                                      uiOutput("Cluster_dowload_button_prior"),
                                                                      # div(DT::dataTableOutput("colors.top_dt")),

                                                             ),
                                                             tabPanel("Epitope",value = "PriorEpiTB"),
                                                             # add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "purple")
                                                 )
                                        ),

                                        # UMAP -> TCR -----
                                        tabPanel("Overview", value = 'up2',
                                                 fluidRow(
                                                   add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                   column(12,selectInput("ID_Column_factor","Order of graph",choices = "", multiple = T, width = "1200px")),
                                                 ),
                                                 fluidRow(
                                                   column(3,selectInput("Split_group_by_overview","Split graph by:",choices ="")),
                                                 ),
                                                 tabsetPanel(id = "QC_panel",
                                                             # TCR overview ontop of UMAP -----
                                                             tabPanel("TCR",
                                                                      tabsetPanel(id = "TCR",
                                                                                  tabPanel("Overlap",
                                                                                           tabsetPanel(
                                                                                             tabPanel("Table",
                                                                                                      div(DT::dataTableOutput("Upset_plot_overlap_Tb")),
                                                                                                      downloadButton('downloaddf_Upset_plot_overlap_Tb','Download table')
                                                                                             ),
                                                                                             tabPanel("Upset Plot",
                                                                                                      plotOutput("Upset_plot_overlap"),
                                                                                                      fluidRow(
                                                                                                        column(1,numericInput("width_Upset_plot_overlap", "Width of PDF", value=10)),
                                                                                                        column(1,numericInput("height_Upset_plot_overlap", "Height of PDF", value=8)),
                                                                                                        column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_Upset_plot_overlap','Download PDF')),
                                                                                                        column(2,numericInput("width_png_Upset_plot_overlap","Width of PNG", value = 1200)),
                                                                                                        column(2,numericInput("height_png_Upset_plot_overlap","Height of PNG", value = 1000)),
                                                                                                        column(2,numericInput("resolution_PNG_Upset_plot_overlap","Resolution of PNG", value = 144)),
                                                                                                        column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_Upset_plot_overlap','Download PNG'))),
                                                                                             ),
                                                                                           )
                                                                                  ),

                                                                                  tabPanel("Clonal expansion plots",value = 2,
                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           fluidRow(
                                                                                             # uiOutput('myPanel_top_clonal_plot'),
                                                                                             conditionalPanel(condition="input.Graph_type_bar=='Number_expanded' || input.Graph_type_bar=='Frequency_expanded'",
                                                                                                              fluidRow(
                                                                                                                column(3,
                                                                                                                       wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                                                 uiOutput('myPanel_clonal_plot'))),
                                                                                                                column(9, plotOutput("clonality.bar.graph",height="600px"))
                                                                                                              )
                                                                                             ),
                                                                                             conditionalPanel(condition="input.Graph_type_bar=='Top_clonotypes'",
                                                                                                              fluidRow(
                                                                                                                column(3,
                                                                                                                       wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                                                 uiOutput('myPanel_top_clonal_plot'))),
                                                                                                                column(9, plotOutput("clonality.bar.graph2",height="600px"))

                                                                                                              ),
                                                                                             ),
                                                                                           ),

                                                                                           fluidRow(
                                                                                             column(1,numericInput("width_clonality.bar.graph", "Width of PDF", value=10)),
                                                                                             column(1,numericInput("height_clonality.bar.graph", "Height of PDF", value=8)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_clonaity.bar.graph','Download PDF')),
                                                                                             column(2,numericInput("width_png_clonality.bar.graph","Width of PNG", value = 1200)),
                                                                                             column(2,numericInput("height_png_clonality.bar.graph","Height of PNG", value = 1000)),
                                                                                             column(2,numericInput("resolution_PNG_clonality.bar.graph","Resolution of PNG", value = 144)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_clonaity.bar.graph','Download PNG'))
                                                                                           ),
                                                                                  ),
                                                                                  # UMAP clonality -> TCR -----
                                                                                  tabPanel("UMAP with clonality (counts)", value = 3,
                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           fluidRow(
                                                                                             column(3,selectInput("filter_umap_expand","Filter plot",choices = c("no","yes"))),


                                                                                           ),
                                                                                           conditionalPanel(condition="input.filter_umap_expand == 'yes'",
                                                                                                            fluidRow(
                                                                                                              column(3,numericInput("UMAP_1x","UMAP_1 <",value = 10)),
                                                                                                              column(3,numericInput("UMAP_1y","UMAP_1 >",value = -10)),
                                                                                                              column(3,numericInput("UMAP_2x","UMAP_2 <",value = 10)),
                                                                                                              column(3,numericInput("UMAP_2y","UMAP_2 >",value = -10))
                                                                                                            ),

                                                                                           ),
                                                                                           conditionalPanel(condition="input.Graph_type_bar=='Number_expanded' || input.Graph_type_bar=='Frequency_expanded'",
                                                                                                            fluidRow(
                                                                                                              column(3,
                                                                                                                     wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                                               uiOutput('cols_UMAP_clonal_plot'))),
                                                                                                              column(9, plotOutput("clonality.TCR.UMAP",height="600px"))
                                                                                                            ),

                                                                                                            fluidRow(
                                                                                                              column(1,numericInput("width_TCR.UMAP", "Width of PDF", value=10)),
                                                                                                              column(1,numericInput("height_TCR.UMAP", "Height of PDF", value=8)),
                                                                                                              column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_TCR.UMAP','Download PDF')),
                                                                                                              column(2,numericInput("width_png_TCR.UMAP","Width of PNG", value = 1600)),
                                                                                                              column(2,numericInput("height_png_TCR.UMAP","Height of PNG", value = 1000)),
                                                                                                              column(2,numericInput("resolution_PNG_TCR.UMAP","Resolution of PNG", value = 144)),
                                                                                                              column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_TCR.UMAP','Download PNG'))
                                                                                                            ),


                                                                                           ),
                                                                                           conditionalPanel(condition="input.Graph_type_bar=='Top_clonotypes'",
                                                                                                            # column(12,   div(DT::dataTableOutput("Tb_For_colouring_check"))),

                                                                                                            fluidRow(
                                                                                                              # column(4,selectInput("Split_by_group","Include group comparison",choices=c("no","yes"))),
                                                                                                              column(4,selectInput("display_all_samps","Display all sample",choices=c("yes","no"))),

                                                                                                            ),

                                                                                                            column(12,selectInput("ID_Column_metadata","Select to display",choices = "", multiple = T, width = "1200px")),

                                                                                                            add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                                            fluidRow(
                                                                                                              column(3,
                                                                                                                     wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                                               uiOutput('cols_UMAP_Topclonotypes'))),
                                                                                                              column(9, plotOutput("clonality.TCR.UMAP.top",height="600px")),
                                                                                                              fluidRow(
                                                                                                                column(1,numericInput("width_TCR.UMAP_top", "Width of PDF", value=10)),
                                                                                                                column(1,numericInput("height_TCR.UMAP_top", "Height of PDF", value=8)),
                                                                                                                column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_TCR.UMAP_top','Download PDF')),
                                                                                                                column(2,numericInput("width_png_TCR.UMAP_top","Width of PNG", value = 1800)),
                                                                                                                column(2,numericInput("height_png_TCR.UMAP_top","Height of PNG", value = 1000)),
                                                                                                                column(2,numericInput("resolution_PNG_TCR.UMAP_top","Resolution of PNG", value = 144)),
                                                                                                                column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_TCR.UMAP_top','Download PNG'))),

                                                                                                            ),
                                                                                           ),

                                                                                  ),
                                                                      ),
                                                             ),


                                                             # T cell classification ------
                                                             tabPanel("GEX", value = "GEX_panel",
                                                                      tabsetPanel(id = "Panel_class",
                                                                                  tabPanel("Percentage",value = 16,
                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           verbatimTextOutput("Percent_tab"),
                                                                                           downloadButton('downloaddf_Percent_tab','Download table')
                                                                                  ),
                                                                                  tabPanel("UMAP plot",value = 14,
                                                                                           fluidRow(column(3,selectInput("show_selected","Show all labels?",choices=c("All","Selected_list"))),
                                                                                                    column(9,uiOutput("SiteNumInput",width = "900px")),
                                                                                           ),

                                                                                           fluidRow(
                                                                                             column(2,numericInput("Filter_lower_UMAP1_marker_GEX","UMAP_1 >",value = -20)),
                                                                                             column(2,numericInput("Filter_lower_UMAP1_marker2_GEX","UMAP_1 <",value = 20)),
                                                                                             column(2,numericInput("Filter_lower_UMAP2_marker_GEX","UMAP_2 >",value = -20)),
                                                                                             column(2,numericInput("Filter_lower_UMAP2_marker2_GEX","UMAP_2 <",value = 20)),
                                                                                           ),

                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           fluidRow(
                                                                                             column(3,
                                                                                                    wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                              uiOutput('cols_UMAP_all_classification'))),
                                                                                             column(9, plotOutput("UMAP_all_classification2",height="600px"))
                                                                                           ),
                                                                                           fluidRow(
                                                                                             column(1,numericInput("width_UMAP_all_classification", "Width of PDF", value=10)),
                                                                                             column(1,numericInput("height_UMAP_all_classification", "Height of PDF", value=8)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_UMAP_all_classification','Download PDF')),
                                                                                             column(2,numericInput("width_png_UMAP_all_classification","Width of PNG", value = 1200)),
                                                                                             column(2,numericInput("height_png_UMAP_all_classification","Height of PNG", value = 1000)),
                                                                                             column(2,numericInput("resolution_PNG_UMAP_all_classification","Resolution of PNG", value = 144)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_UMAP_all_classification','Download PNG'))),
                                                                                  ),

                                                                                  tabPanel("Pie chart",value = 15,
                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           fluidRow(
                                                                                             # column(3,numericInput("strip_size","Size of header",value = 10)),
                                                                                             column(9,p("Fill = Select function type; Group = Select cluster"))
                                                                                           ),
                                                                                           fluidRow(column(3,
                                                                                                           wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                                     uiOutput('myPanel_pie'))), # pie chart
                                                                                                    column(9, plotOutput("Classification_clonotype_pie",height="600px"))
                                                                                           ),
                                                                                           # fluidRow(
                                                                                           #   div(DT::dataTableOutput("table_pie")),
                                                                                           # ),

                                                                                           fluidRow(
                                                                                             column(1,numericInput("width_Classification_clonotype_pie", "Width of PDF", value=10)),
                                                                                             column(1,numericInput("height_Classification_clonotype_pie", "Height of PDF", value=8)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_Classification_clonotype_pie','Download PDF')),
                                                                                             column(2,numericInput("width_png_Classification_clonotype_pie","Width of PNG", value = 1200)),
                                                                                             column(2,numericInput("height_png_Classification_clonotype_pie","Height of PNG", value = 1000)),
                                                                                             column(2,numericInput("resolution_PNG_Classification_clonotype_pie","Resolution of PNG", value = 144)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_Classification_clonotype_pie','Download PNG'))),

                                                                                  ),
                                                                                  # tabPanel("Longitudinal (χ2) ",
                                                                                  #          h5("balloonplot of χ2"),
                                                                                  #          p("The plots include the pearsons residuals for each cell (aka standardized residuals), for the proportion of cell contribution. The Contrib = contribution in percentage (%)"),p("See:",tags$a(href="http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r", " for Information")),
                                                                                  #          selectInput("type_res","Type of comparison",choices = c("Residuals","Contrib")),
                                                                                  #          conditionalPanel(condition="input.type_res=='Residuals'",
                                                                                  #                           fluidRow(
                                                                                  #                             column(3,colourInput("lower_col_chi","Lower colour:","purple")),
                                                                                  #                             column(3,colourInput("mid_col_chi","Mid colour:","White")),
                                                                                  #                             column(3,colourInput("high_col_chi","Higher colour:","Gold")),
                                                                                  #                           )),
                                                                                  #          conditionalPanel(condition="input.type_res=='Contrib'",
                                                                                  #                           fluidRow(
                                                                                  #                             column(3,colourInput("lower_col_chi2","Lower colour:","purple")),
                                                                                  #                             # column(3,colourInput("mid_col_chi","Mid colour:","White")),
                                                                                  #                             column(3,colourInput("high_col_chi2","Higher colour:","Gold")),
                                                                                  #                           )),
                                                                                  #
                                                                                  #          tabsetPanel(
                                                                                  #            tabPanel("Chi-square test",
                                                                                  #                     # div(DT::dataTableOutput()),
                                                                                  #                     verbatimTextOutput("Chi_tab_before"),
                                                                                  #                     h5("Post hoc test"),
                                                                                  #                     div(DT::dataTableOutput("Post_hoc_chi")),
                                                                                  #                     # verbatimTextOutput(),
                                                                                  #            ),
                                                                                  #            tabPanel("balloonplot",
                                                                                  #                     plotOutput("Chi_square_plot",height="600px"),
                                                                                  #                     fluidRow(
                                                                                  #                       column(1,numericInput("width_Chi_square_plot", "Width of PDF", value=10)),
                                                                                  #                       column(1,numericInput("height_Chi_square_plot", "Height of PDF", value=8)),
                                                                                  #                       column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_Chi_square_plot','Download PDF')),
                                                                                  #                       column(2,numericInput("width_png_Chi_square_plot","Width of PNG", value = 1200)),
                                                                                  #                       column(2,numericInput("height_png_Chi_square_plot","Height of PNG", value = 1000)),
                                                                                  #                       column(2,numericInput("resolution_PNG_Chi_square_plot","Resolution of PNG", value = 144)),
                                                                                  #                       column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_Chi_square_plot','Download PNG'))),
                                                                                  #            )
                                                                                  #          ),
                                                                                  # )

                                                                      ),
                                                             ),
                                                             # tabPanel("Clonotype overlap per cluster (upset plot)",
                                                             #          # add in upset plot per cluster
                                                             #          )
                                                 )
                                        ),
                                        # end of differential expression -----
                                        # TCR and GEX analysis section-----
                                        tabPanel("TCR and GEX", value = "TCR_and_GEX_tb",
                                                 #### User interface above the TCR with GEx

                                                 fluidRow(column(3, conditionalPanel( condition="input.Panel_TCRUMAP == 'ClusTCR2'",
                                                                                      selectInput("Samp_col_cluster","Column name",choices = ""))
                                                 ),


                                                 ),

                                                 # column(3, selectInput("V_gene_sc","Select V with/without CDR3",choices ="")),
                                                 fluidRow(
                                                   column(3,
                                                          conditionalPanel(
                                                            condition="input.Panel_TCRUMAP=='top_clone'",
                                                            actionButton("run_top","Update selected clonotype")),
                                                   ),
                                                   column(9,
                                                          conditionalPanel(
                                                            condition="input.Panel_TCRUMAP=='top_clone'",
                                                            selectInput("Selected_clonotype","Select clonotype:",choices ="",width = "900px") ),

                                                   ),
                                                 ),
                                                 fluidRow(


                                                   # column(2,conditionalPanel(condition="input.Panel_TCRUMAP=='top_clone' || input.Panel_TCRUMAP=='Epitope'",
                                                   # selectInput("Split_by_group2","Split by Function?",choices=c("no","yes")))),

                                                   column(9,
                                                          conditionalPanel(condition="input.Panel_TCRUMAP=='top_clone' || input.Panel_TCRUMAP=='Epitope'",
                                                                           selectInput("Graph_split_order","Order of split by:",choices = "", multiple = T, width = "900px"))
                                                   ),
                                                 ),



                                                 conditionalPanel(condition="input.Panel_TCRUMAP=='Epitope'",
                                                                  fluidRow(
                                                                    column(3,uiOutput("classification_to_add_epitope")),
                                                                    column(3,uiOutput("classification_to_add_epitope2")),
                                                                    column(3,actionButton("Update_epi","Impute Epitope list")),
                                                                    column(3,selectInput("Epi_of_interest","Epitope of interest",""))

                                                                  ),
                                                 ),
                                                 # Classification to include ------
                                                 tabsetPanel(id = "Panel_TCRUMAP",
                                                             # top clonotypes plot -----
                                                             tabPanel("Top clonotypes", value = "top_clone",
                                                                      tabsetPanel(
                                                                        # tabPanel("Table",
                                                                        #          add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                        #          div(DT::dataTableOutput("test.table_ridge")),
                                                                        #
                                                                        # ),

                                                                        tabPanel("Summary table",
                                                                                 add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                 div(DT::dataTableOutput("Top_clonotype_sum")),
                                                                                 downloadButton('download_Top_clonotype_sum','Download table')



                                                                        ),
                                                                        tabPanel("Bar graph",
                                                                                 add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),

                                                                                 fluidRow(
                                                                                   column(3,
                                                                                          wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                    uiOutput('myPanel_Top_bar_clonotype'))),
                                                                                   column(9, plotOutput("top_clonotype",height="600px"))
                                                                                 ),
                                                                                 fluidRow(
                                                                                   column(1,numericInput("width_top_clonotype", "Width of PDF", value=10)),
                                                                                   column(1,numericInput("height_top_clonotype", "Height of PDF", value=8)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_top_clonotype','Download PDF')),
                                                                                   column(2,numericInput("width_png_top_clonotype","Width of PNG", value = 1200)),
                                                                                   column(2,numericInput("height_png_top_clonotype","Height of PNG", value = 1000)),
                                                                                   column(2,numericInput("resolution_PNG_top_clonotype","Resolution of PNG", value = 144)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_top_clonotype','Download PNG'))),

                                                                        ),
                                                                        # tabPanel("Pie labs",
                                                                        #          add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                        #          div(DT::dataTableOutput("Top_clonotype_Labs")),
                                                                        #
                                                                        # ),
                                                                        tabPanel("Pie/UMAP chart",

                                                                                 add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                 fluidRow(
                                                                                   column(3,selectInput("Plot_type_selected","Plot",choices = c("pie","UMAP"))),
                                                                                   # column(3,numericInput("strip_size_pie","Size of header",value = 10)),
                                                                                   column(3,numericInput("size_selected_top","Size of Point",value = 2)),
                                                                                 ),
                                                                                 fluidRow(
                                                                                   column(3,
                                                                                          wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                    uiOutput('myPanel_Top_pie_clonotype'))),
                                                                                   column(9, plotOutput("top_clonotype_pie",height="600px")),
                                                                                 ),
                                                                                 fluidRow(
                                                                                   column(1,numericInput("width_top_clonotype_pie", "Width of PDF", value=10)),
                                                                                   column(1,numericInput("height_top_clonotype_pie", "Height of PDF", value=8)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_top_clonotype_pie','Download PDF')),
                                                                                   column(2,numericInput("width_png_top_clonotype_pie","Width of PNG", value = 1200)),
                                                                                   column(2,numericInput("height_png_top_clonotype_pie","Height of PNG", value = 1000)),
                                                                                   column(2,numericInput("resolution_PNG_top_clonotype_pie","Resolution of PNG", value = 144)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_top_clonotype_pie','Download PNG'))),
                                                                        ),
                                                                        # add in find marker for comparing population to other for top clonotype
                                                                        # tabPanel("FindMarker"),

                                                                        # tabPanel("Expression",



                                                                        # tabsetPanel(
                                                                        tabPanel("Ridge/Violin plots",

                                                                                 actionButton("run_string.data_Exp_top","View Ridge plot"),
                                                                                 fluidRow(column(12, selectInput("string.data_Exp_top","column names for summary","",multiple = F, width = "1200px") )),
                                                                                 fluidRow(
                                                                                   column(3, checkboxInput("restric_ex","Restrict to above a threshold?", value = F )),
                                                                                   column(3, numericInput("Gre_ex","Expression above:", value = 0 )),
                                                                                   column(3, selectInput("plot_type_ridgvi","Plot type", choices = c("Ridge (selected clonotype)","Ridge (compare)","Violin (selected clonotype)", "Violin (compare)"))),
                                                                                 ),

                                                                                 add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                 fluidRow(
                                                                                   column(6, div(DT::dataTableOutput("Ridge_chart_alpha_gamma_stat"))),
                                                                                   column(6, plotOutput("Ridge_chart_alpha_gamma_plot_out",height="600px"))

                                                                                 ),
                                                                                 column(6,downloadButton('downloaddf_clusTCR_GEx','Download stats')),
                                                                                 fluidRow(
                                                                                   column(1,numericInput("width_Ridge_chart_alpha_gamma_plot_out", "Width of PDF", value=10)),
                                                                                   column(1,numericInput("height_Ridge_chart_alpha_gamma_plot_out", "Height of PDF", value=8)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_Ridge_chart_alpha_gamma_plot_out','Download PDF')),
                                                                                   column(2,numericInput("width_png_Ridge_chart_alpha_gamma_plot_out","Width of PNG", value = 1200)),
                                                                                   column(2,numericInput("height_png_Ridge_chart_alpha_gamma_plot_out","Height of PNG", value = 1000)),
                                                                                   column(2,numericInput("resolution_PNG_Ridge_chart_alpha_gamma_plot_out","Resolution of PNG", value = 144)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_Ridge_chart_alpha_gamma_plot_out','Download PNG'))),

                                                                        ),

                                                                        tabPanel("Stats",
                                                                                 add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                 div(DT::dataTableOutput("Ridge_chart_alpha_gamma_stat_comp")),
                                                                                 downloadButton('downloaddf_FindMarker_Top','Download stat (Right)')
                                                                        ),
                                                                        # dotplot top-----
                                                                        tabPanel("Dotplot",
                                                                                 fluidRow(
                                                                                   column(2,colourInput("low.dotplot","Lower color:","darkblue")),
                                                                                   column(2,colourInput("middle.dotplot","Middle color:","white")),
                                                                                   column(2,colourInput("high.dotplot","High color:","darkred")),
                                                                                   column(2,checkboxInput("restict_no_points","Restrict Label",value = F)),
                                                                                   column(2,numericInput("pval.ex.top_genes","Top genes to display", value = 40)),
                                                                                 ),
                                                                                 add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),

                                                                                 plotOutput("all_expression_dotplot_top",height="400px"),

                                                                                 textInput("name_clonotype_selected","Name of clone","clone 1"),
                                                                                 fluidRow(
                                                                                   column(1,numericInput("width_all_expression_dotplot_top", "Width of PDF", value=20)),
                                                                                   column(1,numericInput("height_all_expression_dotplot_top", "Height of PDF", value=4)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_all_expression_dotplot_top','Download PDF')),
                                                                                   column(2,numericInput("width_png_all_expression_dotplot_top","Width of PNG", value = 2400)),
                                                                                   column(2,numericInput("height_png_all_expression_dotplot_top","Height of PNG", value = 700)),
                                                                                   column(2,numericInput("resolution_PNG_all_expression_dotplot_top","Resolution of PNG", value = 144)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_all_expression_dotplot_top','Download PNG'))),

                                                                        ),
                                                                        tabPanel("Over-representation",
                                                                                 fluidRow(
                                                                                   column(3,numericInput("in.geneset.cutoff_top","Min number of genes in GeneSet",value = 1, min = 0,step = 1, max = 60000)),
                                                                                   column(3,numericInput("p.val_cutoff_top","p-val cut-off",value = 0.05, min = 0, max = 1)),
                                                                                   # column(3,numericInput("adjust_cutoff_top","BH cut-off",value = 1, min = 0, max = 1)),
                                                                                 ),
                                                                                 add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                 div(DT::dataTableOutput("Over_rep_Top_clones_Tab")),
                                                                                 downloadButton('downloadtb_over.rep.Top_Ex','Download table')
                                                                                 #
                                                                        ),
                                                                        # ),
                                                                        # ),
                                                                      ),
                                                             ),

                                                             # expanded phenotype -----
                                                             tabPanel("Expanded", value = "Expanded",
                                                                      fluidRow(
                                                                        column(3,selectInput("Samp_col_expanded","Sample column name",choices = "")),
                                                                        column(9,selectInput("ID_Column_factor_expanded","ID's to include",choices = "", multiple = T, width = "1200px"))
                                                                      ),

                                                                      fluidRow(
                                                                        # column(3,actionButton("update_Expansion_run","Upate expansion")),


                                                                      ),
                                                                      tabsetPanel(id = "ExPan",
                                                                                  tabPanel("Table",
                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           div(DT::dataTableOutput("Expansion_check")),
                                                                                  ),
                                                                                  tabPanel("ex.UMAP",value = "ExPan_UMAP",
                                                                                           fluidRow(
                                                                                             add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                             column(3,
                                                                                                    wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                              uiOutput('cols_UMAP_Expanded'))),
                                                                                             column(9, plotOutput("UMAP_Expanded",height="600px"))
                                                                                           ),
                                                                                           fluidRow(
                                                                                             column(1,numericInput("width_UMAP_Expanded", "Width of PDF", value=10)),
                                                                                             column(1,numericInput("height_UMAP_Expanded", "Height of PDF", value=8)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_UMAP_Expanded','Download PDF')),
                                                                                             column(2,numericInput("width_png_UMAP_Expanded","Width of PNG", value = 1200)),
                                                                                             column(2,numericInput("height_png_UMAP_Expanded","Height of PNG", value = 1000)),
                                                                                             column(2,numericInput("resolution_PNG_UMAP_Expanded","Resolution of PNG", value = 144)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_UMAP_Expanded','Download PNG'))
                                                                                           ),
                                                                                  ),
                                                                                  tabPanel("Check expansion",
                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           div(DT::dataTableOutput("Expansion_check2")),
                                                                                  ),
                                                                                  tabPanel("Stats",value = "ExPan_stat",
                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           div(DT::dataTableOutput("compare.stat_Ex")),
                                                                                           downloadButton('downloadtb_compare.stat_Ex','Download table')
                                                                                  ),
                                                                                  tabPanel("Dotplot",value = "ExPan_dot",
                                                                                           fluidRow(
                                                                                             column(2,colourInput("low.dotplot.ex","Lower color:","darkblue")),
                                                                                             column(2,colourInput("middle.dotplot.ex","Middle color:","white")),
                                                                                             column(2,colourInput("high.dotplot.ex","High color:","darkred")),

                                                                                           ),
                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           plotOutput("relative_expression_dotplot_ex",height="600px"),
                                                                                           fluidRow(
                                                                                             column(1,numericInput("width_all_expression_dotplot_ex", "Width of PDF", value=20)),
                                                                                             column(1,numericInput("height_all_expression_dotplot_ex", "Height of PDF", value=4)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_all_expression_dotplot_ex','Download PDF')),
                                                                                             column(2,numericInput("width_png_all_expression_dotplot_ex","Width of PNG", value = 2400)),
                                                                                             column(2,numericInput("height_png_all_expression_dotplot_ex","Height of PNG", value = 600)),
                                                                                             column(2,numericInput("resolution_PNG_all_expression_dotplot_ex","Resolution of PNG", value = 144)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_all_expression_dotplot_ex','Download PNG'))
                                                                                           ),
                                                                                  ),
                                                                                  tabPanel("Over-representation",value="ExPan_OvRep",
                                                                                           fluidRow(
                                                                                             column(3,numericInput("in.geneset.cutoff_Exp","Min number of genes in GeneSet",value = 1, min = 0,step = 1, max = 60000)),
                                                                                             column(3,numericInput("p.val_cutoff_Exp","p-val cut-off",value = 0.05, min = 0, max = 1)),
                                                                                             # column(3,numericInput("adjust_cutoff_Exp","BH cut-off",value = 1, min = 0, max = 1)),
                                                                                             # column(3,numericInput("adjust_cutoff_Exp","BH cut-off",value = 0.05, min = 0, max = 1)),

                                                                                           ),

                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           div(DT::dataTableOutput("Over_rep_Exp_Tab")),
                                                                                           downloadButton('downloadtb_over.rep_Exp','Download table')
                                                                                           #
                                                                                  )

                                                                      )

                                                             ),

                                                             # epitope analysis -----
                                                             tabPanel("Epitope", value = "Epitope",
                                                                      tabsetPanel(id = "EpitipeTabs",
                                                                                  # tabPanel("Uploaded Epitope file",
                                                                                  #          add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                  #          div(DT::dataTableOutput("MainTcell_Check")),
                                                                                  # ),
                                                                                  tabPanel("Summary Table",
                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           div(DT::dataTableOutput("Pie_Epitope_dt")),
                                                                                           downloadButton('downloaddf_Pie_Epitope_dt','Download table')
                                                                                  ),

                                                                                  tabPanel("Heatmap",
                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           plotOutput("Heatmap_epi_plot",height="600px"),

                                                                                           fluidRow(
                                                                                             column(1,numericInput("width_Heatmap_epi_plot", "Width of PDF", value=10)),
                                                                                             column(1,numericInput("height_Heatmap_epi_plot", "Height of PDF", value=8)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_Heatmap_epi_plot','Download PDF')),
                                                                                             column(2,numericInput("width_png_Heatmap_epi_plot","Width of PNG", value = 1200)),
                                                                                             column(2,numericInput("height_png_Heatmap_epi_plot","Height of PNG", value = 1000)),
                                                                                             column(2,numericInput("resolution_PNG_Heatmap_epi_plot","Resolution of PNG", value = 144)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_Heatmap_epi_plot','Download PNG'))),


                                                                                  ),
                                                                                  tabPanel("UMAP",

                                                                                           numericInput("value_size_epi_umap","Size of epitope dots",value = 2),

                                                                                           # column(3,selectInput("epitope_umap_selected","Select",choices = c("beta","epitope","pathology"),selected = "pathology")),
                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           fluidRow(
                                                                                             column(3,
                                                                                                    wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                              uiOutput('myPanel_cols_epitope'))),
                                                                                             column(9, plotOutput("UMAP_Epitope_plot",height="600px"))
                                                                                           ),

                                                                                           fluidRow(
                                                                                             column(1,numericInput("width_UMAP_Epitope", "Width of PDF", value=10)),
                                                                                             column(1,numericInput("height_UMAP_Epitope", "Height of PDF", value=8)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_UMAP_Epitope','Download PDF')),
                                                                                             column(2,numericInput("width_png_UMAP_Epitope","Width of PNG", value = 1200)),
                                                                                             column(2,numericInput("height_png_UMAP_Epitope","Height of PNG", value = 1000)),
                                                                                             column(2,numericInput("resolution_PNG_UMAP_Epitope","Resolution of PNG", value = 144)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_UMAP_Epitope','Download PNG'))),

                                                                                  ),
                                                                                  tabPanel("Pie (Expression)",

                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           fluidRow(
                                                                                             column(3,
                                                                                                    wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                              uiOutput('myPanel_cols_epitope_pie'))),
                                                                                             column(9, plotOutput("Pie_Epitope_plot",height="600px"))
                                                                                           ),

                                                                                           fluidRow(
                                                                                             column(1,numericInput("width_Pie_Epitope", "Width of PDF", value=10)),
                                                                                             column(1,numericInput("height_Pie_Epitope", "Height of PDF", value=8)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_Pie_Epitope','Download PDF')),
                                                                                             column(2,numericInput("width_png_Pie_Epitope","Width of PNG", value = 1200)),
                                                                                             column(2,numericInput("height_png_Pie_Epitope","Height of PNG", value = 1000)),
                                                                                             column(2,numericInput("resolution_PNG_Pie_Epitope","Resolution of PNG", value = 144)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_Pie_Epitope','Download PNG'))),
                                                                                  ),

                                                                                  tabPanel("Stats",value = "EpiPan_stat",
                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           div(DT::dataTableOutput("Epi_of_interest_DF")),
                                                                                           div(DT::dataTableOutput("compare.stat_Epi_DT")),
                                                                                           downloadButton('downloadtb_compare.stat_Epi','Download table')
                                                                                  ),
                                                                                  tabPanel("Dotplot",value = "EpiPan_dot",
                                                                                           fluidRow(
                                                                                             column(2,colourInput("low.dotplot.epi","Lower color:","darkblue")),
                                                                                             column(2,colourInput("middle.dotplot.epi","Middle color:","white")),
                                                                                             column(2,colourInput("high.dotplot.epi","High color:","darkred")),
                                                                                             column(2,checkboxInput("restrict.dotpot.epi","Restrict to top list",value = F)),
                                                                                             column(2,numericInput("restrict.dotpot.num.epi","Total genes to display:", value = 10))
                                                                                           ),
                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),


                                                                                           # div(DT::dataTableOutput("checking_epi_dot_issue")),

                                                                                           plotOutput("all_expression_dotplot_epi",height="400px"),

                                                                                           fluidRow(
                                                                                             column(1,numericInput("width_all_expression_dotplot_epi", "Width of PDF", value=20)),
                                                                                             column(1,numericInput("height_all_expression_dotplot_epi", "Height of PDF", value=4)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_all_expression_dotplot_epi','Download PDF')),
                                                                                             column(2,numericInput("width_png_all_expression_dotplot_epi","Width of PNG", value = 2400)),
                                                                                             column(2,numericInput("height_png_all_expression_dotplot_epi","Height of PNG", value = 700)),
                                                                                             column(2,numericInput("resolution_PNG_all_expression_dotplot_epi","Resolution of PNG", value = 144)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_all_expression_dotplot_epi','Download PNG'))
                                                                                           ),



                                                                                  ),
                                                                                  tabPanel("Over-representation",value="EpiPan_OvRep",
                                                                                           fluidRow(
                                                                                             column(3,numericInput("in.geneset.cutoff_Epi","Min number of genes in GeneSet",value = 1, min = 0,step = 1, max = 60000)),
                                                                                             column(3,numericInput("p.val_cutoff_Epi","p-val cut-off",value = 0.05, min = 0, max = 1)),
                                                                                             # column(3,numericInput("adjust_cutoff_Epi","BH cut-off",value = 1, min = 0, max = 1)),
                                                                                             # column(3,numericInput("adjust_cutoff_Exp","BH cut-off",value = 0.05, min = 0, max = 1)),

                                                                                           ),

                                                                                           add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                           div(DT::dataTableOutput("Over_rep_Epi_Tab")),
                                                                                           downloadButton('downloadtb_over.rep.Epi','Download table')
                                                                                           #
                                                                                  )

                                                                      ),

                                                             ),
                                                             # ClusTCR2 Analysis -----
                                                             tabPanel("ClusTCR2",value = "ClusTCR2",
                                                                      tabsetPanel(
                                                                        tabPanel("Table.Clust",
                                                                                 add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                 div(DT::dataTableOutput("Tb_ClusTCR_selected")),
                                                                                 downloadButton('downloadtb_Tb_ClusTCR_selected','Download table')
                                                                        ),
                                                                        tabPanel("UMAP",
                                                                                 fluidRow(
                                                                                   column(3,
                                                                                          wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                    uiOutput('myPanel_cols_clust_UMAP'))),
                                                                                   column(9, plotOutput("UMAP_ClusTCR2_plot",height="600px"))
                                                                                 ),
                                                                                 fluidRow(
                                                                                   column(1,numericInput("width_UMAP_ClusTCR2_plot", "Width of PDF", value=10)),
                                                                                   column(1,numericInput("height_UMAP_ClusTCR2_plot", "Height of PDF", value=8)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_UMAP_ClusTCR2_plot','Download PDF')),
                                                                                   column(2,numericInput("width_png_UMAP_ClusTCR2_plot","Width of PNG", value = 1200)),
                                                                                   column(2,numericInput("height_png_UMAP_ClusTCR2_plot","Height of PNG", value = 1000)),
                                                                                   column(2,numericInput("resolution_PNG_UMAP_ClusTCR2_plot","Resolution of PNG", value = 144)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_UMAP_ClusTCR2_plot','Download PNG'))),

                                                                        ),
                                                                        tabPanel("motif",
                                                                                 plotOutput("Motif_ClusTCR2_cluster",height="300px"),
                                                                                 verbatimTextOutput("print_unique_cases"),
                                                                                 # div(DT::dataTableOutput("Tb_motif_cluster")),

                                                                                 fluidRow(
                                                                                   column(1,numericInput("width_Motif_ClusTCR2_cluster", "Width of PDF", value=10)),
                                                                                   column(1,numericInput("height_Motif_ClusTCR2_cluster", "Height of PDF", value=4)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_Motif_ClusTCR2_cluster','Download PDF')),
                                                                                   column(2,numericInput("width_png_Motif_ClusTCR2_cluster","Width of PNG", value = 2400)),
                                                                                   column(2,numericInput("height_png_Motif_ClusTCR2_cluster","Height of PNG", value = 600)),
                                                                                   column(2,numericInput("resolution_PNG_Motif_ClusTCR2_cluster","Resolution of PNG", value = 144)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_Motif_ClusTCR2_cluster','Download PNG'))),

                                                                        ),
                                                                        tabPanel("Pie (Expression)", value = "ClusPie",
                                                                                 fluidRow(
                                                                                   column(3,
                                                                                          wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                    uiOutput('myPanel_cols_clusTCR2_pie')
                                                                                          )
                                                                                   ),
                                                                                   column(9, plotOutput("Pie_ClusTCR2_plot",height="600px"))
                                                                                 ),
                                                                                 fluidRow(
                                                                                   column(1,numericInput("width_Pie_ClusTCR2_plot", "Width of PDF", value=10)),
                                                                                   column(1,numericInput("height_Pie_ClusTCR2_plot", "Height of PDF", value=8)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_Pie_ClusTCR2_plot','Download PDF')),
                                                                                   column(2,numericInput("width_png_Pie_ClusTCR2_plot","Width of PNG", value = 1200)),
                                                                                   column(2,numericInput("height_png_Pie_ClusTCR2_plot","Height of PNG", value = 1000)),
                                                                                   column(2,numericInput("resolution_PNG_Pie_ClusTCR2_plot","Resolution of PNG", value = 144)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_Pie_ClusTCR2_plot','Download PNG'))
                                                                                 ),
                                                                        ),

                                                                        tabPanel("Stats",value = "ClusPan_stat",
                                                                                 add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),

                                                                                 div(DT::dataTableOutput("compare.stat_Cluster_DT")),
                                                                                 downloadButton('downloadtb_compare.stat_Cluster','Download table')
                                                                        ),

                                                                        # dotplot Cluster ------
                                                                        tabPanel("Dotplot",value = "ClusPan_dot",
                                                                                 fluidRow(
                                                                                   column(2,colourInput("low.dotplot.clust","Lower color:","darkblue")),
                                                                                   column(2,colourInput("middle.dotplot.clust","Middle color:","white")),
                                                                                   column(2,colourInput("high.dotplot.clust","High color:","darkred")),
                                                                                   column(2,checkboxInput("restrict.dotpot.clust","Restrict to top list",value = F)),
                                                                                   column(2,numericInput("restrict.dotpot.num.clust","Total genes to display:", value = 10))
                                                                                 ),
                                                                                 add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),

                                                                                 plotOutput("all_expression_dotplot_cluster",height="400px"),

                                                                                 fluidRow(
                                                                                   column(1,numericInput("width_all_expression_dotplot_clust", "Width of PDF", value=20)),
                                                                                   column(1,numericInput("height_all_expression_dotplot_clust", "Height of PDF", value=8)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_all_expression_dotplot_clust','Download PDF')),
                                                                                   column(2,numericInput("width_png_all_expression_dotplot_clust","Width of PNG", value = 2400)),
                                                                                   column(2,numericInput("height_png_all_expression_dotplot_clust","Height of PNG", value = 700)),
                                                                                   column(2,numericInput("resolution_PNG_all_expression_dotplot_clust","Resolution of PNG", value = 144)),
                                                                                   column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_all_expression_dotplot_clust','Download PNG'))
                                                                                 ),
                                                                        ),
                                                                        tabPanel("Over-representation",value="ClusPan_OvRep",
                                                                                 fluidRow(
                                                                                   column(3,numericInput("in.geneset.cutoff_Clust","Min number of genes in GeneSet",value = 1, min = 0,step = 1, max = 60000)),
                                                                                   column(3,numericInput("p.val_cutoff_Clust","p-val cut-off",value = 0.05, min = 0, max = 1)),
                                                                                   # column(3,numericInput("adjust_cutoff_Clust","BH cut-off",value = 1, min = 0, max = 1)),
                                                                                 ),
                                                                                 add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),

                                                                                 div(DT::dataTableOutput("Over_rep_Cluster_Tab")),
                                                                                 downloadButton('downloadtb_over.rep.cluster','Download table')
                                                                        )
                                                                      )
                                                             ),
                                                             # marker specific TCR analysis --------
                                                             tabPanel("Marker",value = "Marker",
                                                                      actionButton("load_marker_genes","Load genes"),
                                                                      # actionButton("Marker_analysis", "Imput markers"),

                                                                      tabsetPanel(id = "Marker_Panel",
                                                                                  tabPanel("Single marker",

                                                                                           fluidRow(
                                                                                             column(4,selectizeInput("Var_to_col_marker","Marker col","")),
                                                                                             column(2,numericInput("Filter_lower_UMAP1_marker","UMAP_1 >",value = -20)),
                                                                                             column(2,numericInput("Filter_lower_UMAP1_marker2","UMAP_1 <",value = 20)),
                                                                                             column(2,numericInput("Filter_lower_UMAP2_marker","UMAP_2 >",value = -20)),
                                                                                             column(2,numericInput("Filter_lower_UMAP2_marker2","UMAP_2 <",value = 20)),
                                                                                           ),
                                                                                           # conditionalPanel(condition="input.Marker_Panel =='Marker_Panel_plot_stats'",
                                                                                           fluidRow(
                                                                                             column(3,selectInput("Analysis_marker_stats_type","Analysis type",choices = c("Population","Expanded"))),
                                                                                             column(3,selectInput("Help_marer","Definitions",choices = c("No","Yes"))),
                                                                                             column(3, numericInput("cutoff_marker_gt","Marker +ve cut-off (>)",value = 0,step = 0.1)),
                                                                                             column(3,numericInput("pos_expanded_cut_off","clone count (>) for +ve",value = 2, min=1)),
                                                                                           ),
                                                                                           conditionalPanel(condition = "input.Help_marer == 'Yes'",
                                                                                                            # p(strong("1. Total."), "Refers to the marker of interest compared to all negative cells"),
                                                                                                            p(strong("1. Population."), "Refers to the the restricted UMAP selected population and compared the +ve to -ve marker population."),
                                                                                                            p(strong("2. Expanded."), "Compares the clonal expaned (min 3) vs non-expanded population for all cells +ve for the marker in the restricted UMAP space."),
                                                                                           ),

                                                                                           tabsetPanel(
                                                                                             tabPanel("Table", value = "Marker_Panel_table",
                                                                                                      div(DT::dataTableOutput("marker_selected_tab")),
                                                                                             ),
                                                                                             tabPanel("UMAP plot",value = "Marker_Panel_plot_UMAP",
                                                                                                      h4("Select area of the plot to keep for the specific marker"),
                                                                                                      h6("Recommended to filter to broad populations based on UMAP e.g., CD4, CD8 or other"),
                                                                                                      add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                                      numericInput("max_scale","MAX scale",value = ""),

                                                                                                      plotOutput("marker_selected_UMAP_plot",height="600px"),

                                                                                                      fluidRow(
                                                                                                        column(3,numericInput("width_marker_selected_UMAP_plot", "Width of PDF", value=10)),
                                                                                                        column(3,numericInput("height_marker_selected_UMAP_plot", "Height of PDF", value=8)),
                                                                                                        column(3),
                                                                                                        column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_marker_selected_UMAP_plot','Download PDF'))),
                                                                                                      fluidRow(
                                                                                                        column(3,numericInput("width_png_marker_selected_UMAP_plot","Width of PNG", value = 1200)),
                                                                                                        column(3,numericInput("height_png_marker_selected_UMAP_plot","Height of PNG", value = 1000)),
                                                                                                        column(3,numericInput("resolution_PNG_marker_selected_UMAP_plot","Resolution of PNG", value = 144)),
                                                                                                        column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_marker_selected_UMAP_plot','Download PNG'))
                                                                                                      )


                                                                                             ),

                                                                                             tabPanel("Violin/Ridge plot",value = "Marker_Panel_plot_VR",
                                                                                                      h4("Filter marker of interest based on threshold"),
                                                                                                      fluidRow(
                                                                                                        column(3,selectInput("select_plot_vio.ridge","Plot type",choices = c("Violin","Ridge")))
                                                                                                      ),
                                                                                                      add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                                      plotOutput("marker_selected_VioRidge_plot",height="600px")

                                                                                             ),
                                                                                             # tabPanel("UMAP"),
                                                                                             tabPanel("TCR/BCR mapped", value = "Marker_Panel_plot_TCR",
                                                                                                      h4("TCR and/or BCR seqeunces that are positive for that marker"),
                                                                                                      add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                                      fluidRow(
                                                                                                        column(6,div(DT::dataTableOutput("TCR_marker_positive_count"))),
                                                                                                        column(6,div(DT::dataTableOutput("TCR_marker_neg_count")))
                                                                                                      ),

                                                                                                      div(DT::dataTableOutput("merged_marker_hist_table")),
                                                                                                      downloadButton('downloaddf_clonotype_distribution','Download table'),

                                                                                                      add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                                      plotOutput("marker_selected_histogram_plot",height="600px")



                                                                                             ),

                                                                                             tabPanel("Stats",value = "MP_plot_stats",
                                                                                                      add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                                      div(DT::dataTableOutput("Compare.stat_marker")),
                                                                                                      downloadButton('downloaddf_Marker_stats','Download table')

                                                                                                      # filtered on marker of interest, broad population
                                                                                                      # determine if
                                                                                                      # segregated into expanded vs non-expaned  expanded clonotypes?


                                                                                             ),
                                                                                             # tabPanel("Dotplot",
                                                                                             #
                                                                                             #
                                                                                             # ),
                                                                                             # tabPanel("Over-representation",
                                                                                             #
                                                                                             # ),

                                                                                           ),
                                                                                           # cut-off 1,
                                                                                           # umap filtering?


                                                                                  ),
                                                                                  # dual marker analysis ------
                                                                                  tabPanel("Dual marker analysis",
                                                                                           # marker 1 # cut-off 1
                                                                                           # marker 2 @ cut-off 2 (negative control)
                                                                                           fluidRow(
                                                                                             column(4,selectizeInput("Var_to_col_marker2","X-axis Marker","")),
                                                                                             column(4,selectizeInput("Var_to_col_marker3","Y-axis Marker",""))
                                                                                           ),
                                                                                           fluidRow(
                                                                                             column(2,numericInput("Filter_dual_UMAP1_marker","UMAP_1 >",value = -20)),
                                                                                             column(2,numericInput("Filter_dual_UMAP1_marker2","UMAP_1 <",value = 20)),
                                                                                             column(2,numericInput("Filter_dual_UMAP2_marker","UMAP_2 >",value = -20)),
                                                                                             column(2,numericInput("Filter_dual_UMAP2_marker2","UMAP_2 <",value = 20)),
                                                                                             # column(2),
                                                                                             column(2,numericInput("X_axis_dot_dual","X-axis line",value = 0, step = 0.1)),
                                                                                             column(2,numericInput("Y_axis_dot_dual","Y-axis line",value = 0, step = 0.1)),

                                                                                           ),

                                                                                           tabsetPanel(
                                                                                             tabPanel("table",
                                                                                                      add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                                      div(DT::dataTableOutput("meta_data_for_features_scale2_df")),

                                                                                             ),
                                                                                             tabPanel("Feature plots",
                                                                                                      add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                                      fluidRow(
                                                                                                        column(3,numericInput("max_scale2","MAX scale (left)",value = "")),
                                                                                                        column(3,numericInput("max_scale3","MAX scale (right)",value = "")),
                                                                                                      ),

                                                                                                      add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                                      fluidRow(
                                                                                                        column(6,plotOutput("marker_selected_UMAP_plot2",height="600px")),
                                                                                                        column(6,plotOutput("marker_selected_UMAP_plot3",height="600px")),
                                                                                                      )
                                                                                             ),
                                                                                             tabPanel("Dot-plot",

                                                                                                      add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "150px",width = "150px", color = "blue"),
                                                                                                      plotOutput("df_dotplot_marker_plot",height="600px"),

                                                                                                      fluidRow(
                                                                                                        column(1,numericInput("width_df_dotplot_marker_plot", "Width of PDF", value=10)),
                                                                                                        column(1,numericInput("height_df_dotplot_marker_plot", "Height of PDF", value=8)),
                                                                                                        column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_df_dotplot_marker_plot','Download Network PDF')),
                                                                                                        column(2,numericInput("width_png_df_dotplot_marker_plot","Width of PNG", value = 1200)),
                                                                                                        column(2,numericInput("height_png_df_dotplot_marker_plot","Height of PNG", value = 1000)),
                                                                                                        column(2,numericInput("resolution_PNG_df_dotplot_marker_plot","Resolution of PNG", value = 144)),
                                                                                                        column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_df_dotplot_marker_plot','Download Network PNG')),

                                                                                                      ),


                                                                                             ),
                                                                                             # tabPanel("Violin plot"),# will include splitting by T cell/B cell markers? Ig?
                                                                                             # tabPanel("UMAP"),
                                                                                             tabPanel("TCR quadrant table",
                                                                                                      div(DT::dataTableOutput("dual_maker_TCR_Sum_DT")),

                                                                                                      downloadButton('Dule_marker_TCRsummary_DT','Download table')

                                                                                             )
                                                                                           ))
                                                                      )

                                                             ),

                                                             # Overlap ------

                                                             # tabPanel("Clonotypes per cluster (Pie/bar plot)"),
                                                             # tabPanel("Upset plot")
                                                 ),
                                        ),

                            )
                          )
                        )
               ),
               # end of Integration -----
               tabPanel("Post Analysis",
                        sidebarLayout(
                          sidebarPanel(),
                          mainPanel(
                            tabsetPanel(id = "Other_post_analysis",
                                        tabPanel("Contig design"), # upload the unfiltered AIRR file and clonotypes for designing TCR contigs
                                        tabPanel("Converting", value = "converting_PA",
                                                 p("Convert .h5Seurat to .rds")

                                        )


                            )
                          )
                        )

               ),

               # tabPanel("Epitope")
    ) # nav page
  )

  ########
  # server ------
  server <- function(input, output,session) {

    # convert ------
    Convert_to_RDS <- reactive({
      inFile_sc_pro2 <- input$file1_h5Seurat.file
      if (is.null(inFile_sc_pro2)) return(NULL)
      else {
        dataframe = LoadH5Seurat(inFile_sc_pro2$datapath)

      }

    })

    output$Convert_to_RDS_out <- renderPrint({
      sc <- input$file1_h5Seurat.file
      validate(
        need(nrow(sc)>0,
             "Upload files")
      )
      df <- Convert_to_RDS()
      print(df)

    })

    output$downloaddf_SeruatObj_Convert_to_RDS <- downloadHandler(
      filename = function(){
        x = today()
        # paste(input$project_name,"_SC.obj_",x,".h5Seurat", sep = "")
        paste(input$project_name5,"_V4_Seurat",x,".rds", sep = "")
      },
      content = function(file){
        SaveSeuratRds(Convert_to_RDS(), file)
        # SaveH5Seurat(vals_meta.sc$metadata_SCobj,file)
      })

    # add UI ------
    # observeEvent(input$reset_input, {
    #   shinyjs::reset("side-panel")
    # })

    output$scGate_cutoffs <- renderUI({
      if (input$Data_types == '10x_HS' || input$Data_types == 'BD_HS.Full.Panel' || input$Data_types == '10x_MM'|| input$Data_types == 'BD_MM_Full.Panel') {
        numericInput("threshold_scGate","scGate threshold",value = 0.2)
      }

      else {
        numericInput("threshold_scGate","scGate threshold",value = 0.5)
      }
      #c("10x_HS","BD_HS.Immune.Panel","BD_HS.Full.Panel","10x_MM","BD_MM_Full.Panel","BD_MM_Immune.Panel"

    })



    output$classification_to_add <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      df3.meta <- c(names(sc@meta.data))
      # df3.meta <- df3.meta[!grepl("RNA",df3.meta) & !grepl("BCR",df3.meta) & !grepl("TCR",df3.meta)& !grepl("_gene",df3.meta) & !grepl("allele",df3.meta) & !grepl("percent",df3.meta) & !grepl("cdr3",df3.meta)]

      # tags$style("#Colour_By_this {background-color:#e5e5e5;}")

      selectInput("Colour_By_this","Colour by: ",choices = df3.meta,selected="seurat_clusters")

    })
    ## overview colouring
    output$classification_to_add_overview <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      df3.meta <- c(names(sc@meta.data))
      # df3.meta <- df3.meta[!grepl("RNA",df3.meta) & !grepl("BCR",df3.meta) & !grepl("TCR",df3.meta)& !grepl("_gene",df3.meta) & !grepl("allele",df3.meta) & !grepl("percent",df3.meta) & !grepl("cdr3",df3.meta)]
      # tags$style("#Colour_By_this_overview {background-color:#e5e5e5;}")
      selectInput("Colour_By_this_overview","Colour by: ",choices = df3.meta,selected="seurat_clusters")

    })
    output$classification_to_add_epitope <- renderUI({
      sc <- UMAP_metadata_with_labs()

      # summarising the clonality
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )
      df3.meta <- sc@meta.data
      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

      }

      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }

      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta,epi,by="CDR3_beta")
      selectInput("epitope_umap_selected","Colour Pie by (hm = y-axis):",choices = names(df3.meta),selected="beta")
    })
    output$classification_to_add_epitope2 <- renderUI({
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )
      df3.meta <- sc@meta.data
      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

      }
      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }
      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta,epi,by="CDR3_beta")

      selectInput("epitope_umap_selected2","Split Pie by (hm = x-axis):",choices = names(df3.meta),selected="epitope")

    })

    # user interface parameters-----
    output$feature_input <- renderUI({
      if (input$df_seruatobj_type=="10x_Genomics (raw)") {
        fluidRow(
          column(6,numericInput("features.min","minimum features (>)", value = 200)),
          column(6,numericInput("features.max","Maximum features (<)", value = 6000)),
          column(6,numericInput("percent.mt","Mitochondrial DNA cut-off (<)", value = 20)),
          column(6,numericInput("percent.rb","Ribosomal RNA cut-off (>)", value = 5)),
        )
      }

      else if (input$df_seruatobj_type=="10x_Genomics (.h5)") {
        fluidRow(
          column(6,numericInput("features.min","minimum features (>)", value = 200)),
          column(6,numericInput("features.max","Maximum features (<)", value = 6000)),
          column(6,numericInput("percent.mt","Mitochondrial DNA cut-off (<)", value = 20)),
          column(6,numericInput("percent.rb","Ribosomal RNA cut-off (>)", value = 5)),
        )
      }

      else if (input$df_seruatobj_type=="BD Rhapsody (Mouse)") {
        fluidRow(
          column(6,numericInput("features.min","minimum features (>)", value = 200)),
          column(6,numericInput("features.max","Maximum features (<)", value = 6000)),
          column(6,numericInput("percent.mt","Mitochondrial DNA cut-off (<)", value = 20)),
          column(6,numericInput("percent.rb","Ribosomal RNA cut-off (>)", value = 0)),
        )
      }

      else if (input$df_seruatobj_type=="BD Rhapsody (Human Immune panel)") {
        fluidRow(
          column(6,numericInput("features.min","minimum features (>)", value = 45)),
          column(6,numericInput("features.max","Maximum features (<)", value = 160)),
          column(6,numericInput("percent.mt","Mitochondrial DNA cut-off (<)", value = 0)),
          column(6,numericInput("percent.rb","Ribosomal RNA cut-off (>)", value = 0)),
        )
      }

      else {
        fluidRow(
          column(6,numericInput("features.min","minimum features (>)", value = 200)),
          column(6,numericInput("features.max","Maximum features (<)", value = 6000)),
          column(6,numericInput("percent.mt","Mitochondrial DNA cut-off (<)", value = 20)),
          column(6,numericInput("percent.rb","Ribosomal RNA cut-off (>)", value = 0)),
        )
      }

    })
    # human BD rhapsody data -----
    # three files required for BD data: Sample Tag calls, TCR file and count ----
    input.data.calls.bd <- reactive({
      inFile12 <- input$file_calls_BD
      if (is.null(inFile12)) return(NULL)
      else {
        dataframe <- read.csv(
          inFile12$datapath, skip=input$no_lines_skip_Tags, header = T)}
    })

    output$test.files <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- input.data.calls.bd()
      validate(
        need(nrow(calls)>0,
             error_message_val1)
      )
      calls
    })

    input.data.TCR.BD <- reactive({
      inFileTCRBD <- input$file_TCR_BD
      if (is.null(inFileTCRBD)) return(NULL)

      else {
        dataframe <- read.csv(
          inFileTCRBD$datapath,skip=input$no_lines_skip_TCR, header = T)}

    })
    output$test.files2 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- input.data.TCR.BD()
      validate(
        need(nrow(calls)>0,
             error_message_val2)
      )
      calls
    })
    input.data.count.BD <- reactive({
      inFilecountBD <- input$file_counts_BD
      if (is.null(inFilecountBD)) return(NULL)

      else {
        dataframe <- read.csv(inFilecountBD$datapath,skip=input$no_lines_skip_counts, header = T)}

    })
    output$test.files3 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- as.data.frame(input.data.count.BD())
      validate(
        need(nrow(calls)>0,
             error_message_val3)
      )
      calls[1:6,1:6]
    })

    # 10x format is now in BD rhapsody ----
    # barcode file -----
    input.data.barcode.bd <- reactive({
      inFile_bd_barcode <- input$file_barcode_bd
      if (is.null(inFile_bd_barcode)) return(NULL)

      else {
        dataframe <- read.table(
          inFile_bd_barcode$datapath)}
    })

    output$test.files.bd1 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- input.data.barcode.bd()
      validate(
        need(nrow(calls)>0,
             "Upload file")
      )
      calls
    })

    # features file -----
    input.data.features.bd2 <- reactive({
      inFile_bd_features <- input$file_features_bd
      if (is.null(inFile_bd_features)) return(NULL)

      else {
        dataframe <- read.table(
          inFile_bd_features$datapath)}
    })

    output$test.files.bd2 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- input.data.features.bd2()
      validate(
        need(nrow(calls)>0,
             "upload file")
      )
      calls
    })

    # Matrix file ----
    input.data.matrix.bd <- reactive({
      inFile_bd_matrix <- input$file_matrix_bd
      if (is.null(inFile_bd_matrix)) return(NULL)

      else {
        dataframe <- Matrix::readMM(inFile_bd_matrix$datapath)}
    })
    output$test.files.bd3 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- as.data.frame(input.data.matrix.bd())
      validate(
        need(nrow(calls)>0,
             "Upload Matrix")
      )
      calls[1:10,1:10]
    })



    # contig files ----
    input.data.TCR.bd2 <- reactive({
      inFile_bd2_TCR <- input$file_TCR_bd2
      if (is.null(inFile_bd2_TCR)) return(NULL)

      else {
        dataframe <- read.table(inFile_bd2_TCR$datapath,sep = "\t", header = T)}
    })
    output$test.files.bd4 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- input.data.TCR.bd2()
      validate(
        need(nrow(calls)>0,
             "Upload AIRR Contigs (Dominant)")
      )
      calls
    })

    ## BD combining Cell x genes ----
    TCR_Filtering_Paired <- function() {

      TCR <- as.data.frame(input.data.TCR.BD())
      sample_tags <- input.data.calls.bd()

      validate(
        need(nrow(TCR)>0 && nrow(sample_tags)>0,
             "Upload files")
      )

      calls_TCR <- merge(sample_tags,TCR, by ="Cell_Index")
      calls_TCR_count <- calls_TCR

      #filtering to paired TCR
      calls_TCR_count$TRAG <-ifelse(grepl("TR",calls_TCR_count$TCR_Alpha_Gamma_V_gene_Dominant),"TRAG","Missing TRAG gene")
      calls_TCR_count$TRBD <-ifelse(grepl("TR",calls_TCR_count$TCR_Beta_Delta_V_gene_Dominant),"TRBD","Missing TRBD gene")
      calls_TCR_count$TRAG_fun <-ifelse(grepl("[*]",calls_TCR_count$TCR_Alpha_Gamma_CDR3_Translation_Dominant),"Non-functional",
                                        ifelse(grepl("Missing",calls_TCR_count$TRAG),"Missing TCR",

                                               ifelse(grepl("TRGV10",calls_TCR_count$TCR_Alpha_Gamma_V_gene_Dominant),"pseudogene","productive")))

      calls_TCR_count$TRBD_fun <- ifelse(grepl("[*]",calls_TCR_count$TCR_Beta_Delta_CDR3_Translation_Dominant),"Non-functional",
                                         ifelse(grepl("Missing",calls_TCR_count$TRBD),"Missing TCR","productive"))


      calls_TCR_count$paired_TCR <- ifelse(calls_TCR_count$TRAG=="TRAG" & calls_TCR_count$TRBD=="TRBD","paired TCR",
                                           ifelse(calls_TCR_count$TRAG=="TRAG" | calls_TCR_count$TRBD=="TRBD", "Unpaired TCR","No TCR"
                                           ))
      calls_TCR_count$Productive_TCR <- ifelse(calls_TCR_count$TRBD_fun=="productive" & calls_TCR_count$TRAG_fun=="productive","productive TCR",
                                               ifelse(calls_TCR_count$TRBD_fun=="productive" & calls_TCR_count$TRAG_fun=="Missing TCR","productive TCR",
                                                      ifelse(calls_TCR_count$TRBD_fun=="Missing TCR" & calls_TCR_count$TRAG_fun=="productive","productive TCR",
                                                             "unproductive")))

      if (input$BCR_present ==T) {
        calls_TCR_count$v_gene_IgL <- ifelse(grepl("IG",calls_TCR_count$BCR_Light_V_gene_Dominant),"v_gene_IgL","Missing IgL gene")
        calls_TCR_count$v_gene_IgH <- ifelse(grepl("IG",calls_TCR_count$BCR_Heavy_V_gene_Dominant),"v_gene_IgH","Missing IgL gene")
        calls_TCR_count$v_gene_IgL_fun <- ifelse(grepl("[*]",calls_TCR_count$BCR_Light_CDR3_Translation_Dominant),"Non-functional",
                                                 ifelse(grepl("Missing",calls_TCR_count$v_gene_IgL),"Missing BCR","productive"))
        calls_TCR_count$v_gene_IgH_fun <- ifelse(grepl("[*]",calls_TCR_count$BCR_Heavy_CDR3_Translation_Dominant),"Non-functional",
                                                 ifelse(grepl("Missing",calls_TCR_count$v_gene_IgH),"Missing BCR","productive"))
        #
        calls_TCR_count$paired_BCR <- ifelse(calls_TCR_count$v_gene_IgL=="v_gene_IgL" & calls_TCR_count$v_gene_IgH=="v_gene_IgH","paired BCR",
                                             ifelse(calls_TCR_count$v_gene_IgL=="v_gene_IgL" | calls_TCR_count$v_gene_IgH=="v_gene_IgH", "Unpaired BCR","No BCR"
                                             ))


        calls_TCR_count$Productive_BCR <- ifelse(calls_TCR_count$v_gene_IgL_fun=="productive" & calls_TCR_count$v_gene_IgH_fun=="productive","productive BCR",
                                                 ifelse(calls_TCR_count$v_gene_IgL_fun=="productive" & calls_TCR_count$v_gene_IgH_fun=="Missing BCR","productive BCR",
                                                        ifelse(calls_TCR_count$v_gene_IgL_fun=="Missing BCR" & calls_TCR_count$v_gene_IgH_fun=="productive","productive BCR",
                                                               "unproductive")))

        calls_TCR_count$BCR_TCR <- ifelse(calls_TCR_count$Productive_BCR=="productive BCR" & calls_TCR_count$Productive_TCR=="productive TCR","Both TCR and BCR present",NA)
      }

      else {
        calls_TCR_count <- calls_TCR_count
      }

      # filter out non functional TCR
      if (input$BCR_present ==T) {
        productive <- calls_TCR_count[calls_TCR_count$Productive_TCR %in% "productive TCR" | calls_TCR_count$Productive_BCR %in% "productive BCR",] }

      else {
        productive <- calls_TCR_count[calls_TCR_count$Productive_TCR %in% "productive TCR",]
      }
      # filter out un-paired TCR's
      if (input$filtering_TCR ==T) {

        if (input$BCR_present ==T) {
          paired <- productive[productive$paired_TCR %in% "paired TCR" | productive$paired_BCR %in% "paired BCR",]
        }
        else {
          paired <- productive[productive$paired_TCR %in% "paired TCR",]
        }
      }
      else {
        paired <- productive
      }
      paired

    }
    ## BD barcode contig and features ----
    observe({
      if (input$filtered_list == "Paired") {
        updateSelectInput(
          session,
          "locus_column",
          choices=names(input.data.TCR.BD()),
          selected = "chain") }

      else {
        updateSelectInput(
          session,
          "locus_column",
          choices=names(input.data.TCR.bd2()),
          selected = "locus") }
    })


    filtering_required <- reactive({
      contigs <- input.data.TCR.bd2()

      validate(
        need(nrow(contigs)>0,
             "Upload files")
      )

      TCR_unfiltered <- contigs
      TCR_unfiltered$seq_issue <- ifelse(grepl("[*]",TCR_unfiltered$sequence_alignment_aa),"stop-codon","productive")
      TCR_unfiltered_prod <- subset(TCR_unfiltered,TCR_unfiltered$seq_issue=="productive")
      TCR_unfiltered_prod <- subset(TCR_unfiltered_prod,TCR_unfiltered_prod$cdr3_length<30)
      TCR_unfiltered_prod <- subset(TCR_unfiltered_prod,TCR_unfiltered_prod$cdr3_length>5)
      df_unique <- TCR_unfiltered_prod[,names(TCR_unfiltered_prod) %in% c("locus","junction_aa")]
      df_unique$Clonal_Expanded <- 1
      df_unique_sum <- ddply(df_unique,names(df_unique)[-c(3)] ,numcolwise(sum))
      TCR_unfiltered_prod <- merge(TCR_unfiltered_prod,df_unique_sum,by = c("locus","junction_aa"))
      TCR_unfiltered_prod$juct_issue <- ifelse(grepl("CTCSAR",TCR_unfiltered_prod$junction_aa),"CTCSAR",
                                               ifelse(grepl(".TCSAR",TCR_unfiltered_prod$junction_aa),".TCSAR",
                                                      ifelse(grepl('^C', TCR_unfiltered_prod$junction_aa) & grepl('F$', TCR_unfiltered_prod$junction_aa),"C___F",
                                                             ifelse(grepl('^C', TCR_unfiltered_prod$junction_aa),"C__",
                                                                    ifelse(TCR_unfiltered_prod$junction_aa == '',"blank",
                                                                           "other")))))

      TCR_unfiltered_prod_junc <- subset(TCR_unfiltered_prod,TCR_unfiltered_prod$juct_issue!="blank")

      count.chain <- as.data.frame(table(TCR_unfiltered_prod_junc$cell_id))
      names(count.chain) <- c("cell_id","seq_identified")

      TCR_unfiltered_prod_junc_seq <- merge(TCR_unfiltered_prod_junc,count.chain, by = "cell_id")
      sum_tab <- TCR_unfiltered_prod_junc_seq[,names(TCR_unfiltered_prod_junc_seq) %in% c("locus","cell_id")]
      sum_tab$count <- 1
      sum_tab_2 <- as.data.frame(ddply(sum_tab,names(sum_tab)[1:2],numcolwise(sum)))
      count.chain <- as.data.frame(table(TCR_unfiltered_prod_junc_seq$cell_id,TCR_unfiltered_prod_junc_seq$locus))
      names(count.chain) <- c("cell_id","chain","seq_identified")
      TRA <- subset(count.chain,count.chain$chain=="TRA")
      names(TRA)[2:3] <- paste0(names(TRA),"_TRA")[2:3]
      TRB <- subset(count.chain,count.chain$chain=="TRB")
      names(TRB)[2:3] <- paste0(names(TRB),"_TRB")[2:3]
      TRG <- subset(count.chain,count.chain$chain=="TRG")
      names(TRG)[2:3] <- paste0(names(TRG),"_TRG")[2:3]
      TRD <- subset(count.chain,count.chain$chain=="TRD")
      names(TRD)[2:3] <- paste0(names(TRD),"_TRD")[2:3]

      TRAB <-  merge(TRA,TRB,by = "cell_id")
      TRGD <- merge(TRG,TRD,by = "cell_id")
      TRAB_TRGD <- merge(TRAB,TRGD,by = "cell_id")
      names(TRAB_TRGD)
      head(TRAB_TRGD)
      TRAB_TRGD_2 <- TRAB_TRGD %>%
        mutate(Clonality = case_when(
          seq_identified_TRA>0 & seq_identified_TRB>0 & seq_identified_TRG>0 & seq_identified_TRD >0 ~ "AB GD", # all
          # three chains called
          seq_identified_TRA>0 & seq_identified_TRB>0 & seq_identified_TRG>0 & seq_identified_TRD ==0 ~ "AB G",
          seq_identified_TRA>0 & seq_identified_TRB>0 & seq_identified_TRG==0 & seq_identified_TRD >0 ~ "AB D",
          seq_identified_TRA>0 & seq_identified_TRB==0 & seq_identified_TRG>0 & seq_identified_TRD >0 ~ "A GD",
          seq_identified_TRA==0 & seq_identified_TRB>0 & seq_identified_TRG>0 & seq_identified_TRD >0 ~ "B GD",
          # two chains called
          seq_identified_TRA>0 & seq_identified_TRB>0 & seq_identified_TRG==0 & seq_identified_TRD ==0 ~ "AB",
          seq_identified_TRA==0 & seq_identified_TRB==0 & seq_identified_TRG>0 & seq_identified_TRD >0 ~ "GD",

          seq_identified_TRA==0 & seq_identified_TRB>0 & seq_identified_TRG>0 & seq_identified_TRD ==0 ~ "B G",
          seq_identified_TRA==0 & seq_identified_TRB>0 & seq_identified_TRG==0 & seq_identified_TRD >0 ~ "B D",

          seq_identified_TRA>0 & seq_identified_TRB==0 & seq_identified_TRG==0 & seq_identified_TRD >0 ~ "A D",
          seq_identified_TRA>0 & seq_identified_TRB==0 & seq_identified_TRG>0  & seq_identified_TRD==0  ~ "A G",
          # one chains called
          seq_identified_TRA>0 & seq_identified_TRB==0 & seq_identified_TRG==0 & seq_identified_TRD ==0 ~ "A",
          seq_identified_TRA==0 & seq_identified_TRB>0 & seq_identified_TRG==0 & seq_identified_TRD ==0 ~ "B",
          seq_identified_TRA==0 & seq_identified_TRB==0 & seq_identified_TRG>0 & seq_identified_TRD ==0 ~ "G",
          seq_identified_TRA==0 & seq_identified_TRB==0 & seq_identified_TRG==0 & seq_identified_TRD >0 ~ "D",
          TRUE ~ "Other"))
      TRAB_TRGD_paired <- TRAB_TRGD_2
      TRAB_TRGD_paired <- TRAB_TRGD_paired[,!grepl("chain_",names(TRAB_TRGD_paired))]
      TRAB_TRGD_paired$sum <- rowSums(TRAB_TRGD_paired[2:5])
      names(TRAB_TRGD_paired)
      head(TRAB_TRGD_paired)
      TRAB_TRGD_paired$pairing_type <- ifelse(TRAB_TRGD_paired$Clonality =="AB" & TRAB_TRGD_paired$sum==2,"standard AB",
                                              ifelse(TRAB_TRGD_paired$Clonality =="GD" & TRAB_TRGD_paired$sum==2,"standard GD",
                                                     ifelse(TRAB_TRGD_paired$sum==1,"Unpaired",
                                                            "non-standard")))

      df <- merge(TCR_unfiltered_prod_junc_seq,TRAB_TRGD_paired,by = "cell_id")
      df_2 <- df
      Standard_AB <- subset(df_2,df_2$pairing_type == "standard AB")
      Standard_GD <- subset(df_2,df_2$pairing_type == "standard GD")
      non_standard <- subset(df_2,df_2$pairing_type == "non-standard")
      non_standard <- non_standard[order(non_standard$consensus_count, decreasing = T),]
      non_standard$filter <- paste(non_standard$cell_id,non_standard$locus)
      non_standard2 = non_standard[!duplicated(non_standard$filter),]
      non_standard2 <- non_standard2[order(non_standard2$cell_id, decreasing = F),]
      unique(non_standard2$Clonality)
      non_standard2$Standard_to_keep <- ifelse(non_standard2$Clonality=="AB G" & non_standard2$locus == "TRG", "remove_G",
                                               ifelse(non_standard2$Clonality=="AB D" & non_standard2$locus == "TRD", "remove_D",
                                                      ifelse(non_standard2$Clonality== "B GD" & non_standard2$locus == "TRB", "remove_B",
                                                             ifelse(non_standard2$Clonality== "A GD" & non_standard2$locus == "TRA", "remove_A",
                                                                    ifelse(non_standard2$Clonality== "AB GD" & c(non_standard2$locus == "TRG" | non_standard2$locus == "TRD"), "remove_GD","keep")))))

      non_standard <- subset(non_standard2,non_standard2$Standard_to_keep=="keep")
      non_standard <- non_standard[,!names(non_standard) %in% c("filter","Standard_to_keep")]
      dim(non_standard)


      Filtered_TCR_list <- rbind(Standard_AB,Standard_GD,non_standard)
      Filtered_TCR_list <- Filtered_TCR_list[order(Filtered_TCR_list$cell_id, decreasing = F),]
      contigs <- Filtered_TCR_list
      contigs

    })

    tb_bd_meta.data_TCR <- function () {
      contigs <- input.data.TCR.bd2()
      sample_tags <- input.data.calls.bd()
      validate(
        need(nrow(contigs)>0 && nrow(sample_tags)>0,
             "Upload files")
      )

      if (input$filtered_list == "Unfiltered") {

        contigs <- filtering_required()

      } else {
        contigs <- contigs
      }
      names(contigs)[names(contigs) %in% "cell_id"] <- "Cell_Index"
      contigs_merge <- merge(contigs,sample_tags,by="Cell_Index")

      # remove non-functional sequences
      if (nrow(contigs_merge[-c(grep("[*]",contigs_merge$junction_aa)),]>0)) {
        contigs_merge <- contigs_merge[-c(grep("[*]",contigs_merge$junction_aa)),]
      }

      # removed undefined sequences

      contigs_merge <- subset(contigs_merge,contigs_merge$Sample_Name != "Multiplet" )
      contigs_merge <- subset(contigs_merge,contigs_merge$Sample_Name != "Undetermined")

      if (input$filtered_list == "Unfiltered") {
        contigs_lim <- contigs_merge[!names(contigs_merge) %in% c("consensus_count","sequence_id","duplicate_count","germline_alignment","reads","length","cdr3","rev_comp","complete_vdj",names(contigs_merge[grep("fwr",names(contigs_merge))]),"cdr1","cdr2","productive",names(contigs_merge[grep("sequence",names(contigs_merge))]),names(contigs_merge[grep("cigar",names(contigs_merge))]),names(contigs_merge[grep("support",names(contigs_merge))]),"Sample_Tag","sum","dominant","putative_cell","juct_issue","seq_issue","c_call"
                                                                  # "sum","dominant","putative_cell","seq_issue","juct_issue"
        ) ]
      } else {
        contigs_lim <- contigs_merge[!names(contigs_merge) %in% c("consensus_count","sequence_id","duplicate_count","germline_alignment","reads","length","cdr3","rev_comp","complete_vdj",names(contigs_merge[grep("fwr",names(contigs_merge))]),"cdr1","cdr2","productive",names(contigs_merge[grep("sequence",names(contigs_merge))]),names(contigs_merge[grep("cigar",names(contigs_merge))]),names(contigs_merge[grep("support",names(contigs_merge))]),"Sample_Tag"
                                                                  # "sum","dominant","putative_cell","seq_issue","juct_issue"
        ) ]

      }

      contigs_lim$v_gene <- gsub("[*]0.","",contigs_lim$v_call)
      contigs_lim$j_gene <- gsub("[*]0.","",contigs_lim$j_call)
      contigs_lim$d_gene <- gsub("[*]0.","",contigs_lim$d_call)
      # contigs_lim$v_gene_BD <- gsub("[*]0.","",calls_TCR_paired.fun$TCR_Beta_Delta_V_gene_Dominant)
      # names(contigs_lim)   <-  gsub("_call","_gene",names(contigs_lim))
      contigs_lim

      names(contigs_lim)[names(contigs_lim) %in% input$locus_column] <- "chain"

      contig_AG <- subset(contigs_lim,contigs_lim$chain=="TRA" | contigs_lim$chain=="TRG")

      if (input$filtered_list == "Unfiltered") {
        contig_AG <- contig_AG[!duplicated(contig_AG$Cell_Index),]
        name.list <- names(contig_AG[c(names(contig_AG[grep("gene",names(contig_AG))]),
                                       names(contig_AG[grep("call",names(contig_AG))]),
                                       names(contig_AG[grep("cdr3",names(contig_AG))]),
                                       names(contig_AG[grep("cdr2",names(contig_AG))]),
                                       names(contig_AG[grep("cdr1",names(contig_AG))]),
                                       names(contig_AG[grep("junction",names(contig_AG))]),
                                       "chain",
                                       "cell_type_experimental",
                                       "Clonal_Expanded"
        )])
      } else {
        name.list <- names(contig_AG[c(names(contig_AG[grep("gene",names(contig_AG))]),
                                       names(contig_AG[grep("call",names(contig_AG))]),
                                       names(contig_AG[grep("cdr3",names(contig_AG))]),
                                       names(contig_AG[grep("cdr2",names(contig_AG))]),
                                       names(contig_AG[grep("cdr1",names(contig_AG))]),
                                       names(contig_AG[grep("junction",names(contig_AG))]),
                                       "chain",
                                       "cell_type_experimental"
        )])
      }

      contig_AG <- contig_AG %>%
        select(all_of(name.list), everything())
      names(contig_AG)[1:summary(name.list)[1]] <-paste(names(contig_AG[names(contig_AG) %in% name.list]),"_AG",sep="")
      # contig_AG

      if (input$filtered_list == "Unfiltered") {
        if (length(contigs_lim$Clonality[contigs_lim$Clonality == "B D"])>0) {


          contig_D <- subset(contigs_lim,contigs_lim$chain=="TRD" & contigs_lim$Clonality == "B D")
          name.list <- names(contig_D[c(names(contig_D[grep("gene",names(contig_D))]),
                                        names(contig_D[grep("call",names(contig_D))]),
                                        names(contig_D[grep("cdr3",names(contig_D))]),
                                        names(contig_D[grep("cdr2",names(contig_D))]),
                                        names(contig_D[grep("cdr1",names(contig_D))]),
                                        names(contig_D[grep("junction",names(contig_D))]),
                                        "chain",
                                        "cell_type_experimental",

                                        "Clonal_Expanded"
          )])

          contig_D <- contig_D %>%
            select(all_of(name.list), everything())
          names(contig_D)[1:summary(name.list)[1]] <-paste(names(contig_D[names(contig_D) %in% name.list]),"_AG",sep="")

          contig_AG <- rbind(contig_AG,contig_D)
        } else {
          contig_AG <- contig_AG
        }
      }

      contig_BD <- subset(contigs_lim,contigs_lim$chain=="TRB" | contigs_lim$chain=="TRD")

      if (input$filtered_list == "Unfiltered") {
        name.list <- names(contig_BD[c(names(contig_BD[grep("gene",names(contig_BD))]),
                                       names(contig_BD[grep("call",names(contig_BD))]),
                                       names(contig_BD[grep("cdr3",names(contig_BD))]),
                                       names(contig_BD[grep("cdr2",names(contig_BD))]),
                                       names(contig_BD[grep("cdr1",names(contig_BD))]),
                                       names(contig_BD[grep("junction",names(contig_BD))]),
                                       "chain",
                                       "cell_type_experimental",
                                       "Clonal_Expanded"
        )])
      } else {

        name.list <- names(contig_BD[c(names(contig_BD[grep("gene",names(contig_BD))]),
                                       names(contig_BD[grep("call",names(contig_BD))]),
                                       names(contig_BD[grep("cdr3",names(contig_BD))]),
                                       names(contig_BD[grep("cdr2",names(contig_BD))]),
                                       names(contig_BD[grep("cdr1",names(contig_BD))]),
                                       names(contig_BD[grep("junction",names(contig_BD))]),
                                       "chain",
                                       "cell_type_experimental"
        )])
      }

      contig_BD <- contig_BD %>%
        select(all_of(name.list), everything())
      names(contig_BD)[1:summary(name.list)[1]] <-paste(names(contig_BD[names(contig_BD) %in% name.list]),"_BD",sep="")

      if (input$filtered_list == "Unfiltered") {
        if (length(contigs_lim$Clonality[contigs_lim$Clonality == "B D"])>0) {

          contig_B <- subset(contigs_lim,contigs_lim$chain=="TRB" & contigs_lim$Clonality == "B D")

          name.list <- names(contig_B[c(names(contig_B[grep("gene",names(contig_B))]),
                                        names(contig_B[grep("call",names(contig_B))]),
                                        names(contig_B[grep("cdr3",names(contig_B))]),
                                        names(contig_B[grep("cdr2",names(contig_B))]),
                                        names(contig_B[grep("cdr1",names(contig_B))]),
                                        names(contig_B[grep("junction",names(contig_B))]),
                                        "chain",
                                        "cell_type_experimental",
                                        "Clonal_Expanded"
          )])

          contig_B <- contig_B %>%
            select(all_of(name.list), everything())
          names(contig_B)[1:summary(name.list)[1]] <-paste(names(contig_B[names(contig_B) %in% name.list]),"_BD",sep="")

          head(contig_B)
          contig_BD$other <- ifelse(contig_BD$Cell_Index %in% contig_B$Cell_Index & contig_BD$chain_BD=="TRD","duplicate","other")
          contig_BD <- subset(contig_BD,contig_BD$other=="other")
          contig_BD <- contig_BD[,!names(contig_BD) %in% "other"]
          contig_BD <- rbind(contig_BD,contig_B)

        }}


      contig_BD

      if (input$filtered_list == "Unfiltered") {
        contig_paired <- merge(contig_AG,contig_BD, by=c("Cell_Index","Sample_Name",names(contig_BD[grep("seq_",names(contig_BD))]),"Clonality", "pairing_type"),all = T)
      } else {
        contig_paired <- merge(contig_AG,contig_BD, by=c("cell_id","Sample_Name"),all = T)
      }
      contig_paired

      # if (input$filtered_list == "Unfiltered") {
      #   if (length(contigs_lim$Clonality[contigs_lim$Clonality == "B D"])>0) {
      contig_paired$pairing <- ifelse(contig_paired$Clonality=="AB","abTCR Paired",
                                      ifelse(contig_paired$Clonality=="GD","gdTCR Paired",
                                             ifelse(contig_paired$Clonality=="B D", "bdTCR Paired",
                                                    ifelse(contig_paired$Clonality=="B G", "bgTCR Paired",
                                                           ifelse(contig_paired$Clonality=="A D", "adTCR Paired",
                                                                  ifelse(contig_paired$Clonality == "B","orphan B",
                                                                         ifelse(contig_paired$Clonality == "A","orphan A",
                                                                                ifelse(contig_paired$Clonality == "G","orphan G",
                                                                                       ifelse(contig_paired$Clonality == "D","orphan D",NA)))))))))

      # }} else {
      #   contig_paired$pairing <- ifelse(contig_paired$chain_BD=="TRB" & contig_paired$chain_AG=="TRA","abTCR Paired",
      #                                   ifelse(contig_paired$chain_BD=="TRD" & contig_paired$chain_AG=="TRG","gdTCR Paired",NA
      #                                   ))
      #   }




      contig_paired$pairing[is.na(contig_paired$pairing)] <- "OTHER"

      # name.list2 <- names(contig_paired)[!names(contig_paired)%in% c("d_gene_AG", "c_gene_AG","c_gene_BD")]
      # contig_paired <- contig_paired[,names(contig_paired) %in% name.list2]
      contig_paired

      contig_paired_only <- contig_paired
      contig_paired_only[is.na(contig_paired_only)] <- "None"
      contig_paired_only[contig_paired_only == ''] <- "None"
      contig_paired_only$d_gene_BD[contig_paired_only$d_gene_BD == 'None'] <- "_"
      # contig_paired_only$d_gene_AG[contig_paired_only$d_gene_AG == 'None'] <- "_"
      contig_paired_only
      if (input$filtering_TCR ==T) {
        contig_paired_only <- subset(contig_paired_only,contig_paired_only$junction_AG!="None")
        contig_paired_only <- subset(contig_paired_only,contig_paired_only$junction_BD!="None")
        contig_paired_only
      }
      contig_paired_only
      contig_paired_only$vj_gene_AG <- paste(contig_paired_only$v_gene_AG,contig_paired_only$j_gene_AG,sep = ".")
      contig_paired_only$vj_gene_AG <- gsub("None.None","",contig_paired_only$vj_gene_AG)
      contig_paired_only
      contig_paired_only$vj_gene_BD <- paste(contig_paired_only$v_gene_BD,contig_paired_only$j_gene_BD,sep = ".")
      contig_paired_only$vj_gene_BD <- gsub(".NA.",".",contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub("[.]None[.]",".",contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub("None.None","",contig_paired_only$vj_gene_BD)

      contig_paired_only$vdj_gene_BD <- paste(contig_paired_only$v_gene_BD,contig_paired_only$d_gene_BD,contig_paired_only$j_gene_BD,sep = ".")
      contig_paired_only$vdj_gene_BD <- gsub("[.]_[.]",".",contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("None.None",".",contig_paired_only$vdj_gene_BD)
      #
      contig_paired_only$vj_gene_cdr3_AG <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$junction_aa_AG,sep = "_")
      contig_paired_only$vj_gene_cdr3_AG <- gsub("_None","",contig_paired_only$vj_gene_cdr3_AG)
      #
      contig_paired_only$vj_gene_cdr3_BD <- paste(contig_paired_only$vj_gene_BD,contig_paired_only$junction_aa_BD,sep = "_")
      contig_paired_only$vj_gene_cdr3_BD <- gsub("_None","",contig_paired_only$vj_gene_cdr3_BD)
      #
      contig_paired_only$vdj_gene_cdr3_BD <- paste(contig_paired_only$vdj_gene_BD,contig_paired_only$junction_aa_BD,sep = "_")
      contig_paired_only$vdj_gene_cdr3_BD <- gsub("_None","",contig_paired_only$vdj_gene_cdr3_BD)
      #
      contig_paired_only$vj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$vj_gene_BD,sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$vdj_gene_BD,sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- gsub("^ & ","",contig_paired_only$vdj_gene_AG_BD)
      contig_paired_only$vdj_gene_AG_BD <- gsub(" & $","",contig_paired_only$vdj_gene_AG_BD)
      #
      # #updating names to be consistant....
      contig_paired_only$vj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vj_gene_cdr3_BD,sep = " & ")
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub("^ & ","",contig_paired_only$vj_gene_cdr3_AG_BD)
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub(" & $","",contig_paired_only$vj_gene_cdr3_AG_BD)

      contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vdj_gene_cdr3_BD,sep = " & ")
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub("^ & ","",contig_paired_only$vdj_gene_cdr3_AG_BD)
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub(" & $","",contig_paired_only$vdj_gene_cdr3_AG_BD)
      # contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vdj_gene_cdr3_BD,sep = " & ")
      names(contig_paired_only)[names(contig_paired_only) %in% "cell_id"] <- "Cell_Index"
      contig_paired_only <- contig_paired_only[!duplicated(contig_paired_only$Cell_Index),] # remove duplicates


      contig_paired_only <- merge(contig_paired_only,sample_tags,by=c("Cell_Index","Sample_Name"),all = T)

      if (input$filtered_list == "Unfiltered") {
        contig_paired_only <- contig_paired_only %>%
          select(all_of(c("Cell_Index","Sample_Name","Sample_Tag","pairing","pairing_type","Clonality","seq_identified","seq_identified_TRA","seq_identified_TRB","seq_identified_TRG","seq_identified_TRD","cell_type_experimental_AG","cell_type_experimental_BD","chain_AG","chain_BD",
                          names(contig_paired_only[grep("call",names(contig_paired_only))]))),
                 everything())
      }

      else {
        contig_paired_only <- contig_paired_only %>%
          select(all_of(c("Cell_Index","Sample_Name","Sample_Tag","pairing","cell_type_experimental_AG","cell_type_experimental_BD","chain_AG","chain_BD",
                          names(contig_paired_only[grep("call",names(contig_paired_only))]))),
                 everything())

      }


      # names(contig_paired_only) <- gsub("_aa","",names(contig_paired_only) )
      contig_paired_only
      # merge(sample_tags,contig_paired_only,by.x="cell_id",by.y="Cell_Index",all=T)
    }


    ## create Sample_tags file =====

    samp.tags <- reactive({
      if (input$filtered_list == "Paired") {
        TCR <- as.data.frame(TCR_Filtering_Paired())
      }

      else {
        TCR <- input.data.TCR.bd2()
      }


      validate(
        need(nrow(TCR)>0,
             error_message_val2)
      )
      if (input$filtered_list == "Paired") {
        Sample_Tags <- as.data.frame(TCR$Cell_Index)

      }

      else {
        Sample_Tags <- as.data.frame(unique(TCR$cell_id))
      }

      names(Sample_Tags) <- "Cell_Index"
      Sample_Tags$Sample_Tag <- "SampleTag01"
      #
      Sample_Tags$Sample_Name <- input$sample_tags_name

      Sample_Tags <- Sample_Tags[order(Sample_Tags$Cell_Index),]
      # Sample_Tags
      Sample_Tags

    })

    output$tb_sample_tags_created <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      samp.tags()
    })

    output$downloadtb_sample_tags <- downloadHandler(
      filename = function(){
        x = today()
        paste(input$name.BD,"_Sample_Tags_",x,".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(samp.tags())
        write.csv(df,file, row.names = F)
      } )

    ### testing table merge ------
    observe({
      updateSelectInput(
        session,
        "V_gene_AG_BDrhap",
        choices=names(input.data.TCR.BD()),
        selected = "TCR_Alpha_Gamma_V_gene_Dominant") }) # junction sequence
    observe({
      updateSelectInput(
        session,
        "Junction_AG_BDrhap",
        choices=names(input.data.TCR.BD()),
        selected = "TCR_Alpha_Gamma_CDR3_Translation_Dominant") }) # junction sequence
    observe({
      updateSelectInput(
        session,
        "V_gene_BD_BDrhap",
        choices=names(input.data.TCR.BD()),
        selected = "TCR_Beta_Delta_V_gene_Dominant") }) # junction sequence
    observe({
      updateSelectInput(
        session,
        "Junction_BD_BDrhap",
        choices=names(input.data.TCR.BD()),
        selected = "TCR_Beta_Delta_CDR3_Translation_Dominant") }) # junction sequence

    output$Check_table <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      df_nocouts <- TCR_Filtering_Paired()

      if (nrow( df_nocouts[df_nocouts$Total_VDJ_Read_Count != 0, ]>0)) {
        df_nocouts2 <- df_nocouts[df_nocouts$Total_VDJ_Read_Count != 0, ]
      }
      df_nocouts2

      df_nocouts2_AG <- df_nocouts2[,names(df_nocouts2) %in% c("Sample_Name",input$V_gene_AG_BDrhap,input$Junction_AG_BDrhap)]
      names(df_nocouts2_AG) <- c("Sample_Name","v_call","junction_aa")
      df_nocouts2_AG

      df_nocouts2_BD <- df_nocouts2[,names(df_nocouts2) %in% c("Sample_Name",input$V_gene_BD_BDrhap,input$Junction_BD_BDrhap ,"","")]
      names(df_nocouts2_BD) <- c("Sample_Name","v_call","junction_aa")

      df_nocouts3 <-as.data.frame(rbind(df_nocouts2_AG,df_nocouts2_BD))
      df_nocouts3$v_call <- gsub("^$",NA, df_nocouts3$v_call)
      df_nocouts3 <- df_nocouts3[complete.cases(df_nocouts3$v_call), ]

      df_nocouts3$junction_aa <- gsub("^$",NA, df_nocouts3$junction_aa)
      df_nocouts3 <- df_nocouts3[complete.cases(df_nocouts3$junction_aa), ]

      if (nrow(df_nocouts3[-c(grep("[*]",df_nocouts3$junction_aa)),]>0)) {
        df_nocouts3 <- df_nocouts3[-c(grep("[*]",df_nocouts3$junction_aa)),]
      }


      df_nocouts3 <- df_nocouts3[!duplicated(df_nocouts3$junction_aa), ]
      df_nocouts3

    })

    ## for the clusTCR ----
    df_clusTCR <- function () {
      df_nocouts <- TCR_Filtering_Paired()

      if (nrow( df_nocouts[df_nocouts$Total_VDJ_Read_Count != 0, ]>0)) {
        df_nocouts2 <- df_nocouts[df_nocouts$Total_VDJ_Read_Count != 0, ]
      }

      #
      df_nocouts2_AG <- df_nocouts2[,names(df_nocouts2) %in% c("Sample_Name",input$V_gene_AG_BDrhap,input$Junction_AG_BDrhap)]
      names(df_nocouts2_AG) <- c("Sample_Name","v_call","junction_aa")
      df_nocouts2_AG

      df_nocouts2_BD <- df_nocouts2[,names(df_nocouts2) %in% c("Sample_Name",input$V_gene_BD_BDrhap,input$Junction_BD_BDrhap)]
      names(df_nocouts2_BD) <- c("Sample_Name","v_call","junction_aa")

      if (input$BCR_present ==T) {
        df_nocouts2_IgL <- df_nocouts2[,names(df_nocouts2) %in% c("Sample_Name","BCR_Light_CDR3_Translation_Dominant","BCR_Light_V_gene_Dominant")]
        df_nocouts2_IgH <- df_nocouts2[,names(df_nocouts2) %in% c("Sample_Name","BCR_Heavy_CDR3_Translation_Dominant","BCR_Heavy_V_gene_Dominant")]
        names(df_nocouts2_IgL) <- c("Sample_Name","v_call","junction_aa")
        names(df_nocouts2_IgH) <- c("Sample_Name","v_call","junction_aa")
        df_nocouts3 <-as.data.frame(rbind(df_nocouts2_AG,df_nocouts2_BD,df_nocouts2_IgL,df_nocouts2_IgH))
      }

      else {
        df_nocouts3 <-as.data.frame(rbind(df_nocouts2_AG,df_nocouts2_BD))
      }

      df_nocouts3$v_call <- gsub("^$",NA, df_nocouts3$v_call)
      df_nocouts3 <- df_nocouts3[complete.cases(df_nocouts3$v_call), ]

      df_nocouts3$junction_aa <- gsub("^$",NA, df_nocouts3$junction_aa)
      df_nocouts3 <- df_nocouts3[complete.cases(df_nocouts3$junction_aa), ]

      if (nrow(df_nocouts3[-c(grep("[*]",df_nocouts3$junction_aa)),]>0)) {
        df_nocouts3 <- df_nocouts3[-c(grep("[*]",df_nocouts3$junction_aa)),]
      }

      df_nocouts3 <- df_nocouts3[!duplicated(df_nocouts3$junction_aa), ]
      df_nocouts3

    }

    tb_bd_contigues_contig <- reactive({
      contigs <- input.data.TCR.bd2()
      sample_tags <- input.data.calls.bd()
      validate(
        need(nrow(contigs)>0 && nrow(sample_tags)>0,
             "Upload files")
      )

      if (input$filtered_list == "Unfiltered") {

        df <- filtering_required()
        contigs <- contigs[contigs$cell_id %in% unique(df$cell_id),]

      }
      else {
        contigs <- contigs
      }

      contigs_merge <- merge(contigs,sample_tags,by.x="cell_id",by.y="Cell_Index")

      contigs_merge <- subset(contigs_merge,contigs_merge$Sample_Name != "Multiplet" )
      contigs_merge <- subset(contigs_merge,contigs_merge$Sample_Name != "Undetermined")

      contigs2 <- contigs[names(contigs) %in% c("v_call","junction_aa")]
      contigs2$v_call <- gsub("[*]0.","",contigs2$v_call)

      # contigs2$Sample_Name <- sample_tags$Sample_Name
      contigs2[is.na(contigs2)] <- "None"
      contigs2[contigs2 == ''] <- "None"

      if (nrow(contigs2[-c(grep("[*]",contigs2$junction_aa)),]>0)) {
        contigs2 <- contigs2[-c(grep("[*]",contigs2$junction_aa)),]
      }

      contigs2 <- subset(contigs2,contigs2$junction_aa!= "None")
      contigs2[!duplicated(contigs2[,c('v_call','junction_aa')]),]
    })



    df_clusTCR_Filtered <- reactive({
      if (input$filtered_list == "Paired") {

        df <- df_clusTCR()
      }
      else {
        df <- tb_bd_contigues_contig()
        df <- df[grep("^C",df$junction_aa),]
      }

      if (input$chain_clusTCR2_bd == "AG") {
        (rbind(df[grep("TRAV",df$v_call),],df[grep("TRGV",df$v_call),]))

      }


      else if (input$chain_clusTCR2_bd == "IgH") {
        rbind(df[grep("IGH",df$v_call),])

      }

      else if (input$chain_clusTCR2_bd == "IgLK") {
        (rbind(df[grep("IGL",df$v_call),],df[grep("IGK",df$v_call),]))

      }


      else {
        (rbind(df[grep("TRBV",df$v_call),],df[grep("TRDV",df$v_call),]))
      }

    })

    output$tb_clusTCR <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      if (input$filtered_list == "Paired") {
        df_clusTCR()
      }
      else {
        df_clusTCR_Filtered()
      }
    })
    output$downloaddf_clusTCR <- downloadHandler(
      filename = function(){
        x = today()
        if (input$chain_clusTCR2_bd == "AG") {
          paste("AG_",input$name.BD,"","_ClusTCR_",x,".csv", sep = "")

        }

        else if (input$chain_clusTCR2_bd == "IgH") {
          paste("IgH_",input$name.BD,"_ClusTCR_",x,".csv", sep = "")
        }

        else if (input$chain_clusTCR2_bd == "IgLK") {
          paste("IgLK_",input$name.BD,"_ClusTCR_",x,".csv", sep = "")

        }

        else {
          paste("BD_",input$name.BD,"_ClusTCR_",x,".csv", sep = "")
        }

      },
      content = function(file){

        df <- as.data.frame(df_clusTCR_Filtered())
        write.csv(df,file, row.names = F)

      } )

    ## count table for Seurat----
    tb_bd_matrix <- function () {
      barcode <- input.data.barcode.bd()
      features <- input.data.features.bd2()

      mmMat <- as.data.frame(input.data.matrix.bd())

      validate(
        need(nrow(barcode)>0 & nrow(features)>0 & nrow(mmMat)>0,
             "Upload files")
      )

      rownames(mmMat) <- make.unique(features$V2) # select which feature to label genes...
      names(mmMat) <- barcode$V1
      mmMat$Gene_Name <- rownames(mmMat)
      mmMat <- mmMat %>%
        select(all_of("Gene_Name"), everything())
      mmMat

    }

    df_count.matrix_bd <- reactive( {
      mat <- as.data.frame(input.data.count.BD())

      validate(
        need(nrow(mat)>0,
             "Upload file")
      )
      mat <- as.data.frame(mat)
      head(mat)[1:6]
      rownames(mat) <- mat$Cell_Index

      mat <- subset(mat, select = -Cell_Index)
      head(mat)[1:6]

      mat <- as.data.frame(t(mat))
      head(mat)[1:6]

      rownames(mat) <- gsub("^X","",rownames(mat))
      head(mat)[1:6]
      names(mat) <- gsub("^X","",names(mat))
      mat
    })

    output$tb_count_matrix <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      if (input$Format_bd=='cellXgene'){
        (df_count.matrix_bd())[1:6,1:6]
      }
      else {
        tb_bd_matrix()[1:6,1:6]

      }


    })
    output$downloadtb_count_matrix <- downloadHandler(
      filename = function(){

        if (input$Format_bd=='cellXgene'){
          paste(input$name.BD,"_BD_Count_Matrix_",gsub("-", ".", Sys.Date()),".csv", sep = "")
        }
        else {
          paste(input$name.BD,"_count-matrix_bd_",gsub("-", ".", Sys.Date()),".csv.gz", sep = "")

        }
      },
      content = function(file){
        if (input$Format_bd=='cellXgene'){

          df <- as.data.frame(df_count.matrix_bd())
          # write.table(,file, row.names = T)
          write.csv(df, file, row.names = T)
        }
        else {
          df <- as.data.frame(tb_bd_matrix())
          # write.table(,file, row.names = T)
          write_csv(df,gzfile(file))

        }

      })

    ### meta.data for seruat ----
    meta.data_for_Seuratobj <- function () {

      if (input$filtered_list == "Paired") {
        contig_paired_only <- TCR_Filtering_Paired()
        contig_paired_only
        contig_paired_only$v_gene_AG <- gsub("[*]0.","",contig_paired_only$TCR_Alpha_Gamma_V_gene_Dominant)
        contig_paired_only$j_gene_AG <- gsub("[*]0.","",contig_paired_only$TCR_Alpha_Gamma_J_gene_Dominant)
        contig_paired_only$cdr3_AG <- gsub("[*]0.","",contig_paired_only[,names(contig_paired_only) %in% input$Junction_AG_BDrhap])
        contig_paired_only$v_gene_BD <- gsub("[*]0.","",contig_paired_only$TCR_Beta_Delta_V_gene_Dominant)
        contig_paired_only$j_gene_BD <- gsub("[*]0.","",contig_paired_only$TCR_Beta_Delta_J_gene_Dominant)
        contig_paired_only$d_gene_BD <- gsub("[*]0.","",contig_paired_only$TCR_Beta_Delta_D_gene_Dominant)
        contig_paired_only$cdr3_BD <- gsub("[*]0.","",contig_paired_only[,names(contig_paired_only) %in% input$Junction_BD_BDrhap])

        contig_paired_only$d_gene_BD <- gsub("^$","NA", contig_paired_only$d_gene_BD)
        #
        contig_paired_only$vj_gene_AG <- paste(contig_paired_only$v_gene_AG,contig_paired_only$j_gene_AG,sep = ".")
        contig_paired_only$vj_gene_AG <- gsub("NA.NA","",contig_paired_only$vj_gene_AG)
        #
        contig_paired_only$vj_gene_BD <- paste(contig_paired_only$v_gene_BD,contig_paired_only$j_gene_BD,sep = ".")
        contig_paired_only$vj_gene_BD <- gsub(".NA.",".",contig_paired_only$vj_gene_BD)
        contig_paired_only$vj_gene_BD <- gsub(".None.",".",contig_paired_only$vj_gene_BD)
        contig_paired_only$vj_gene_BD <- gsub("NA.NA","",contig_paired_only$vj_gene_BD)
        #
        contig_paired_only$vdj_gene_BD <- paste(contig_paired_only$v_gene_BD,contig_paired_only$d_gene_BD,contig_paired_only$j_gene_BD,sep = ".")
        contig_paired_only$vdj_gene_BD <- gsub(".NA.",".",contig_paired_only$vdj_gene_BD)
        contig_paired_only$vdj_gene_BD <- gsub("[.]NA[.]",".",contig_paired_only$vdj_gene_BD)
        contig_paired_only$vdj_gene_BD <- gsub("NA.NA","",contig_paired_only$vdj_gene_BD)
        #
        contig_paired_only$vj_gene_cdr3_AG <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$cdr3_AG,sep = "_")
        contig_paired_only$vj_gene_cdr3_AG <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_AG)
        #
        contig_paired_only$vj_gene_cdr3_BD <- paste(contig_paired_only$vj_gene_BD,contig_paired_only$cdr3_BD,sep = "_")
        contig_paired_only$vj_gene_cdr3_BD <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_BD)
        #
        contig_paired_only$vdj_gene_cdr3_BD <- paste(contig_paired_only$vdj_gene_BD,contig_paired_only$cdr3_BD,sep = "_")
        contig_paired_only$vdj_gene_cdr3_BD <- gsub("_NA","",contig_paired_only$vdj_gene_cdr3_BD)
        #
        contig_paired_only$vj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$vj_gene_BD,sep = " & ")
        contig_paired_only$vdj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$vdj_gene_BD,sep = " & ")
        contig_paired_only$vdj_gene_AG_BD <- gsub("^ & ","",contig_paired_only$vdj_gene_AG_BD)
        contig_paired_only$vdj_gene_AG_BD <- gsub(" & $","",contig_paired_only$vdj_gene_AG_BD)
        #
        # #updating names to be consistant....
        contig_paired_only$vj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vj_gene_cdr3_BD,sep = " & ")
        contig_paired_only$vj_gene_cdr3_AG_BD <- gsub("^ & ","",contig_paired_only$vj_gene_cdr3_AG_BD)
        contig_paired_only$vj_gene_cdr3_AG_BD <- gsub(" & $","",contig_paired_only$vj_gene_cdr3_AG_BD)

        contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vdj_gene_cdr3_BD,sep = " & ")
        contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub("^ & ","",contig_paired_only$vdj_gene_cdr3_AG_BD)
        contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub(" & $","",contig_paired_only$vdj_gene_cdr3_AG_BD)


        contig_paired_only <- contig_paired_only[!names(contig_paired_only) %in% c("Total_VDJ_Read_Count","Total_VDJ_Molecule_Count","TCR_Alpha_Gamma_Read_Count","TCR_Alpha_Gamma_Molecule_Count","TCR_Beta_Delta_Read_Count","TCR_Beta_Delta_Molecule_Count","Sample_Tag","TCR_Alpha_Gamma_C_gene_Dominant","TCR_Beta_Delta_C_gene_Dominant","TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant","TCR_Beta_Delta_CDR3_Nucleotide_Dominant") ]
        contig_paired_only
        if (input$BCR_present ==T) {

          contig_paired_only$v_gene_IgL <- gsub("[*]0.","",contig_paired_only$BCR_Light_V_gene_Dominant)
          contig_paired_only$j_gene_IgL <- gsub("[*]0.","",contig_paired_only$BCR_Light_J_gene_Dominant)
          contig_paired_only$v_gene_IgH <- gsub("[*]0.","",contig_paired_only$BCR_Heavy_V_gene_Dominant)
          contig_paired_only$d_gene_IgH <- gsub("[*]0.","",contig_paired_only$BCR_Heavy_D_gene_Dominant)
          contig_paired_only$j_gene_IgH <- gsub("[*]0.","",contig_paired_only$BCR_Heavy_J_gene_Dominant)

          contig_paired_only$v_gene_IgH_L <- paste(contig_paired_only$v_gene_IgH,contig_paired_only$v_gene_IgL,sep = " & ")

          contig_paired_only$v_gene_cdr3_IgH <- paste(contig_paired_only$v_gene_IgH,contig_paired_only$BCR_Heavy_CDR3_Translation_Dominant,sep = "_")
          contig_paired_only$v_gene_cdr3_IgL <- paste(contig_paired_only$v_gene_IgL,contig_paired_only$BCR_Light_CDR3_Translation_Dominant,sep = "_")
          contig_paired_only$v_gene_cdr3_IgH_L <- paste(contig_paired_only$v_gene_cdr3_IgH,contig_paired_only$v_gene_cdr3_IgL,sep = ".")

          names(contig_paired_only)[names(contig_paired_only) %in% c("BCR_Light_CDR3_Translation_Dominant")] <- "cdr3_IgL"
          names(contig_paired_only)[names(contig_paired_only) %in% c("BCR_Light_V_gene_Dominant")] <- "v_allele_IgL"
          names(contig_paired_only)[names(contig_paired_only) %in% c("BCR_Light_J_gene_Dominant")] <- "j_allele_IgL"
          names(contig_paired_only)[names(contig_paired_only) %in% c("BCR_Heavy_CDR3_Translation_Dominant")] <- "cdr3_IgH"
          names(contig_paired_only)[names(contig_paired_only) %in% c("BCR_Heavy_V_gene_Dominant")] <- "v_allele_IgH"
          names(contig_paired_only)[names(contig_paired_only) %in% c("BCR_Heavy_J_gene_Dominant")] <- "j_allele_IgH"
          names(contig_paired_only)[names(contig_paired_only) %in% c("BCR_Heavy_D_gene_Dominant")] <- "d_allele_IgH"
          contig_paired_only

        }

        else {
          contig_paired_only
        }
      }

      else {
        tb_bd_meta.data_TCR()
      }

    }

    output$tb_metadata_sc <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      meta.data_for_Seuratobj()
    })

    output$downloadtb_metadata_sc <- downloadHandler(
      filename = function(){
        x = today()
        paste(input$name.BD,"_Meta.data_",x,".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(meta.data_for_Seuratobj())
        write.csv(df,file, row.names = F)
      } )


    # for TCRex output -----
    TCRex_BDrap_df <- reactive({

      calls_TCR_paired.fun <- TCR_Filtering_Paired()


      if (nrow(calls_TCR_paired.fun[-c(grep("[*]",calls_TCR_paired.fun$TCR_Beta_Delta_CDR3_Translation_Dominant)),]>0)) {
        calls_TCR_paired.fun <- calls_TCR_paired.fun[-c(grep("[*]",calls_TCR_paired.fun$TCR_Beta_Delta_CDR3_Translation_Dominant)),] # remove st
      }

      calls_TCR_paired.fun$TRBV_gene <- gsub("[*]0.","",calls_TCR_paired.fun$TCR_Beta_Delta_V_gene_Dominant)
      calls_TCR_paired.fun$CDR3_beta <- paste("C",calls_TCR_paired.fun$TCR_Beta_Delta_CDR3_Translation_Dominant,"F",sep="")
      calls_TCR_paired.fun$TRBJ_gene<- gsub("[*]0.","",calls_TCR_paired.fun$TCR_Beta_Delta_J_gene_Dominant)
      calls_TCR_paired.fun$cloneCount <- 1
      calls_TCR_paired.fun2 <- calls_TCR_paired.fun[,names(calls_TCR_paired.fun) %in% c("TRBV_gene","CDR3_beta","TRBJ_gene","cloneCount")]

      calls_TCR_paired.fun3 <- ddply(calls_TCR_paired.fun2,names(calls_TCR_paired.fun2)[-c(4)] ,numcolwise(sum))
      calls_TCR_paired.fun3 <- subset(calls_TCR_paired.fun3,calls_TCR_paired.fun3$CDR3_beta != "CF")

      calls_TCR_paired.fun3 <- calls_TCR_paired.fun3[!grepl("TRD", calls_TCR_paired.fun3$TRBV_gene),]
      calls_TCR_paired.fun3

    })

    TCRex_Bd_df <- reactive({
      contigs <- input.data.TCR.bd2()
      sample_tags <- input.data.calls.bd()
      validate(
        need(nrow(contigs)>0 && nrow(sample_tags)>0,
             "Upload files")
      )
      if (input$filtered_list == "Unfiltered") {
        df <- filtering_required()
        contigs <- contigs[contigs$cell_id %in% unique(df$cell_id),]

      }
      else {
        contigs <- contigs
      }
      contigs_merge <- merge(contigs,sample_tags,by.x="cell_id",by.y="Cell_Index")

      contigs_merge <- subset(contigs_merge,contigs_merge$Sample_Name != "Multiplet" )
      contigs_merge <- subset(contigs_merge,contigs_merge$Sample_Name != "Undetermined")

      contigs2 <- contigs_merge[,names(contigs_merge) %in% c("v_call","j_call","junction_aa")]
      contigs2[is.na(contigs2)] <- "None"
      contigs2[contigs2 == ''] <- "None"
      contigs2
      names(contigs2)[names(contigs2) %in% c("junction_aa")] <- "CDR3_beta"
      if (nrow(contigs2[-c(grep("[*]",contigs2$CDR3_beta)),]>0)) {
        contigs2 <- contigs2[-c(grep("[*]",contigs2$CDR3_beta)),]
      }
      names(contigs2)[names(contigs2) %in% c("v_call")] <- "TRBV_gene"
      names(contigs2)[names(contigs2) %in% c("j_call")] <- "TRBJ_gene"
      contigs2 <- contigs2 %>%
        select(TRBV_gene,CDR3_beta,everything())
      contigs2 <- subset(contigs2,contigs2$CDR3_beta!= "None")
      contigs2 <- contigs2
      contigs2$TRBV_gene <- gsub("[*]0.","",contigs2$TRBV_gene)
      contigs2$TRBJ_gene <- gsub("[*]0.","",contigs2$TRBJ_gene)
      contigs2$cloneCount <- 1
      calls_TCR_paired.fun3 <- ddply(contigs2,names(contigs2)[-c(4)] ,numcolwise(sum))
      calls_TCR_paired.fun3 <- calls_TCR_paired.fun3[grepl("TRB", calls_TCR_paired.fun3$TRBV_gene),]
      calls_TCR_paired.fun3
    })

    tb_TCRex_BDrap_df <- reactive ({
      if (input$filtered_list == "Paired") {
        TCRex_BDrap_df()
      }
      else {
        TCRex_Bd_df()
      }

    })

    output$tb_TCRex_BDrap_df <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      tb_TCRex_BDrap_df()
    })


    output$downloaddf_TCRex_BDrap <- downloadHandler(
      filename = function(){
        x = today()
        paste(input$name.BD,"_TCRex_",x,".tsv", sep = "")

      },
      content = function(file){

        df <- as.data.frame(tb_TCRex_BDrap_df())
        write.table(df,file, row.names = F,sep="\t", quote = F)


      } )

    ## TCR_Explore compatible -----
    df_TCR_Explore <- function () {
      calls_TCR_paired.fun <- meta.data_for_Seuratobj()
      calls_TCR_paired.fun$cloneCount <- 1 # added for TCR_Explore (total number of cells with that gene and sequence)

      calls_TCR_paired.fun <- calls_TCR_paired.fun[!names(calls_TCR_paired.fun) %in% c("Total_VDJ_Read_Count","Total_VDJ_Molecule_Count","TCR_Alpha_Gamma_Read_Count","TCR_Alpha_Gamma_Molecule_Count","TCR_Beta_Delta_Read_Count","TCR_Beta_Delta_Molecule_Count","Cell_Type_Experimental") ]


      calls_TCR_paired.fun[order(calls_TCR_paired.fun$vdj_gene_AG_BD),]

      calls_TCR_paired.fun$vdj_gene_AG_BD <- gsub("^$",NA, calls_TCR_paired.fun$vdj_gene_AG_BD)
      calls_TCR_paired.fun <- calls_TCR_paired.fun[complete.cases(calls_TCR_paired.fun$vdj_gene_AG_BD), ]

      calls_TCR_paired.fun <- calls_TCR_paired.fun %>%
        select(cloneCount, everything())

      calls_TCR_paired.fun
    }

    output$tb_TCR_Explore <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      df_TCR_Explore()
    })

    output$downloadtb_TCR_Explore <- downloadHandler(
      filename = function(){
        paste(input$name.BD,"TCR_Explore.compatible",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(df_TCR_Explore())
        write.csv(df,file, row.names = F)
      } )


    # multi-TCR files -----

    multi_TCR_BDrhap <- reactive({
      contigs <- input.data.TCR.bd2()
      sample_tags <- input.data.calls.bd()
      validate(
        need(nrow(contigs)>0 && nrow(sample_tags)>0,
             "Upload files")
      )

      if (input$filtered_list == "Unfiltered") {
        sample_tags$order <- 1:dim(sample_tags)[1]
        head(sample_tags)
        names(contigs)
        contigs_merge <- merge(contigs,sample_tags,by.x="cell_id",by.y="Cell_Index") # need this for BD rhap to remove multiplets/undetermined

        # remove non-functional sequences
        if (nrow(contigs_merge[-c(grep("[*]",contigs_merge$junction_aa)),]>0)) {
          contigs_merge <- contigs_merge[-c(grep("[*]",contigs_merge$junction_aa)),]
        }
        names(contigs_merge)

        contigs_merge <- contigs_merge[contigs_merge[,names(contigs_merge) %in% "productive"]=="True",]

        contigs_merge <- contigs_merge[grep("^C",contigs_merge$junction_aa),]
        contigs_merge <- subset(contigs_merge,contigs_merge$Sample_Name != "Multiplet" )
        contigs_merge <- subset(contigs_merge,contigs_merge$Sample_Name != "Undetermined")
        contigs_lim <- contigs_merge[!names(contigs_merge) %in% c("consensus_count","sequence_id","duplicate_count","germline_alignment","reads","length","cdr3","rev_comp","complete_vdj",names(contigs_merge[grep("fwr",names(contigs_merge))]),names(contigs_merge[grep("cdr1",names(contigs_merge))]),"junction",names(contigs_merge[grep("length",names(contigs_merge))]),"productive",names(contigs_merge[grep("cdr2",names(contigs_merge))]),"productive",names(contigs_merge[grep("sequence",names(contigs_merge))]),names(contigs_merge[grep("cigar",names(contigs_merge))]),names(contigs_merge[grep("support",names(contigs_merge))]),"Sample_Tag","sum","dominant","putative_cell","juct_issue","seq_issue","c_call","cell_type_experimental")]
        contigs_lim$v_gene <- gsub("[*]0.","",contigs_lim$v_call)
        contigs_lim$j_gene <- gsub("[*]0.","",contigs_lim$j_call)
        contigs_lim$d_gene <- gsub("[*]0.","",contigs_lim$d_call)
        contigs_lim[grep("call",names(contigs_lim))]

        names(contigs_lim) %in% names(contigs_lim[grep("call",names(contigs_lim))])

        contigs_lim <- contigs_lim[,!names(contigs_lim) %in% names(contigs_lim[grep("call",names(contigs_lim))])]
        head(contigs_lim)
        contigs_lim$cdr3_aa <- gsub("$^","none",contigs_lim$cdr3_aa)
        contigs_lim <- subset(contigs_lim,contigs_lim$cdr3_aa!= "none")

        # contigs_lim$v_gene_BD <- gsub("[*]0.","",calls_TCR_paired.fun$TCR_Beta_Delta_V_gene_Dominant)
        # names(contigs_lim)   <-  gsub("_call","_gene",names(contigs_lim))
        names(contigs_lim)

        names(contigs_lim)[names(contigs_lim) %in% "locus"] <- "chain" #input$locus_column

        contig_A <- subset(contigs_lim,contigs_lim$chain=="TRA")
        head(contig_A)
        name.list <- names(contig_A[c(names(contig_A[grep("gene",names(contig_A))]),
                                      names(contig_A[grep("call",names(contig_A))]),
                                      names(contig_A[grep("cdr3",names(contig_A))]),
                                      names(contig_A[grep("junction",names(contig_A))]),
                                      "chain"
        )])

        contig_A <- contig_A %>%
          select(all_of(name.list), everything())
        names(contig_A)[1:summary(name.list)[1]] <-paste(names(contig_A[names(contig_A) %in% name.list]),"_A",sep="")

        names(contig_A) %in% c("cell_id","Sample_Name")

        contig_G <- subset(contigs_lim,contigs_lim$chain=="TRG")
        head(contig_G)
        name.list <- names(contig_G[c(names(contig_G[grep("gene",names(contig_G))]),
                                      names(contig_G[grep("call",names(contig_G))]),
                                      names(contig_G[grep("cdr3",names(contig_G))]),
                                      names(contig_G[grep("junction",names(contig_G))]),
                                      "chain"
        )])
        head(contig_A)
        contig_G <- contig_G %>%
          select(all_of(name.list), everything())
        names(contig_G)[1:summary(name.list)[1]] <-paste(names(contig_G[names(contig_G) %in% name.list]),"_G",sep="")

        names(contig_G) %in% c("cell_id","Sample_Name")


        # contig_A

        contig_B <- subset(contigs_lim,contigs_lim$chain=="TRB")
        head(contig_B)
        name.list <- names(contig_B[c(names(contig_B[grep("gene",names(contig_B))]),
                                      names(contig_B[grep("call",names(contig_B))]),
                                      names(contig_B[grep("cdr3",names(contig_B))]),
                                      names(contig_B[grep("junction",names(contig_B))]),
                                      "chain"
        )])

        contig_B <- contig_B %>%
          select(all_of(name.list), everything())
        names(contig_B)[1:summary(name.list)[1]] <-paste(names(contig_B[names(contig_B) %in% name.list]),"_B",sep="")

        names(contig_B) %in% c("cell_id","Sample_Name")


        contig_D <- subset(contigs_lim,contigs_lim$chain=="TRD")
        head(contig_D)
        name.list <- names(contig_D[c(names(contig_D[grep("gene",names(contig_D))]),
                                      names(contig_D[grep("call",names(contig_D))]),
                                      names(contig_D[grep("cdr3",names(contig_D))]),
                                      names(contig_D[grep("junction",names(contig_D))]),
                                      "chain"
        )])

        contig_D <- contig_D %>%
          select(all_of(name.list), everything())
        names(contig_D)[1:summary(name.list)[1]] <-paste(names(contig_D[names(contig_D) %in% name.list]),"_D",sep="")

        names(contig_D) %in% c("cell_id","Sample_Name")

        contig_paired_AB <- merge(contig_A,contig_B, by=c("cell_id","Sample_Name","order"),all = T)
        contig_paired_GD <- merge(contig_G,contig_D, by=c("cell_id","Sample_Name","order"),all = T)
        contig_paired <- merge(contig_paired_AB,contig_paired_GD, by=c("cell_id","Sample_Name","order"),all = T)

        contig_paired <- contig_paired[!duplicated(contig_paired[names(contig_paired) %in% c("cell_id","Sample_Name","order","junction_aa_A","junction_aa_B","junction_aa_G","junction_aa_D")]),]
        dim(contig_paired)

        contig_paired[is.na(contig_paired)] <- "-"
        contig_paired <- as.data.frame(contig_paired)

        head(contig_paired)

        contig_paired$pairing <-  ifelse(contig_paired$chain_B=="TRB" & contig_paired$chain_A=="TRA" & contig_paired$chain_D=="TRD" & contig_paired$chain_G=="TRG","AB GD",
                                         ifelse(contig_paired$chain_B=="TRB" & contig_paired$chain_A=="TRA" & contig_paired$chain_D!="TRD" & contig_paired$chain_G=="TRG","AB G",
                                                ifelse(contig_paired$chain_B=="TRB" & contig_paired$chain_A=="TRA" & contig_paired$chain_D=="TRD" & contig_paired$chain_G!="TRG", "AB D",
                                                       ifelse(contig_paired$chain_B!="TRB" & contig_paired$chain_A=="TRA" & contig_paired$chain_D=="TRD" & contig_paired$chain_G=="TRG","AG D",
                                                              ifelse(contig_paired$chain_B=="TRB" & contig_paired$chain_A!="TRA" & contig_paired$chain_D=="TRD" & contig_paired$chain_G=="TRG","B GD",
                                                                     ifelse(contig_paired$chain_B=="TRB" & contig_paired$chain_A=="TRA" & contig_paired$chain_D!="TRD" & contig_paired$chain_G!="TRG","AB",
                                                                            ifelse(contig_paired$chain_B!="TRB" & contig_paired$chain_A!="TRA" & contig_paired$chain_D=="TRD" & contig_paired$chain_G=="TRG","GD",

                                                                                   ifelse(contig_paired$chain_B=="TRB" & contig_paired$chain_A!="TRA" & contig_paired$chain_D!="TRD" & contig_paired$chain_G=="TRG","B G",
                                                                                          ifelse(contig_paired$chain_B!="TRB" & contig_paired$chain_A=="TRA" & contig_paired$chain_D=="TRD" & contig_paired$chain_G!="TRG","A D",
                                                                                                 "other"
                                                                                          )))))))))





        unique(contig_paired$pairing)
        contig_paired

        contig_paired_only <- contig_paired
        head(contig_paired_only)
        # contig_paired_only <- subset(contig_paired_only,contig_paired_only$pairing!="single")

        AB_G <- contig_paired_only[contig_paired_only$pairing %in% c("AB GD","AB G","B GD","B G"),]

        grep("junction_aa",(names(AB_G)))

        junction_only <- contig_paired_only[,c(names(contig_paired_only)[grep("junction_aa",(names(contig_paired_only)))],names(contig_paired_only)[grep("v_gene",(names(contig_paired_only)))], names(contig_paired_only)[grep("chain_",(names(contig_paired_only)))]),]
        names.list <- names(junction_only)
        junction_only$cloneCount <- 1

        grep("^C",junction_only$junction_aa_A)

        junction_only_sum <- ddply(junction_only,names.list ,numcolwise(sum))

        junction_only_sum$pairing <-  ifelse(junction_only_sum$chain_B=="TRB" & junction_only_sum$chain_A=="TRA" & junction_only_sum$chain_D=="TRD" & junction_only_sum$chain_G=="TRG","AB GD",
                                             ifelse(junction_only_sum$chain_B=="TRB" & junction_only_sum$chain_A=="TRA" & junction_only_sum$chain_D!="TRD" & junction_only_sum$chain_G=="TRG","AB G",
                                                    ifelse(junction_only_sum$chain_B=="TRB" & junction_only_sum$chain_A=="TRA" & junction_only_sum$chain_D=="TRD" & junction_only_sum$chain_G!="TRG", "AB D",
                                                           ifelse(junction_only_sum$chain_B!="TRB" & junction_only_sum$chain_A=="TRA" & junction_only_sum$chain_D=="TRD" & junction_only_sum$chain_G=="TRG","AG D",
                                                                  ifelse(junction_only_sum$chain_B=="TRB" & junction_only_sum$chain_A!="TRA" & junction_only_sum$chain_D=="TRD" & junction_only_sum$chain_G=="TRG","B GD",
                                                                         ifelse(junction_only_sum$chain_B=="TRB" & junction_only_sum$chain_A=="TRA" & junction_only_sum$chain_D!="TRD" & junction_only_sum$chain_G!="TRG","AB",
                                                                                ifelse(junction_only_sum$chain_B!="TRB" & junction_only_sum$chain_A!="TRA" & junction_only_sum$chain_D=="TRD" & junction_only_sum$chain_G=="TRG","GD",
                                                                                       ifelse(junction_only_sum$chain_B=="TRB" & junction_only_sum$chain_A!="TRA" & junction_only_sum$chain_D!="TRD" & junction_only_sum$chain_G=="TRG","B G",
                                                                                              ifelse(junction_only_sum$chain_B=="TRB" & junction_only_sum$chain_A!="TRA" & junction_only_sum$chain_D=="TRD" & junction_only_sum$chain_G!="TRG","B D",
                                                                                                     ifelse(junction_only_sum$chain_B!="TRB" & junction_only_sum$chain_A=="TRA" & junction_only_sum$chain_D=="TRD" & junction_only_sum$chain_G!="TRG","A D",
                                                                                                            ifelse(junction_only_sum$chain_B!="TRB" & junction_only_sum$chain_A=="TRA" & junction_only_sum$chain_D!="TRD" & junction_only_sum$chain_G=="TRG","A G",
                                                                                                                   ifelse(junction_only_sum$chain_B!="TRB" & junction_only_sum$chain_A=="TRA" & junction_only_sum$chain_D!="TRD" & junction_only_sum$chain_G!="TRG","A",
                                                                                                                          ifelse(junction_only_sum$chain_B=="TRB" & junction_only_sum$chain_A!="TRA" & junction_only_sum$chain_D!="TRD" & junction_only_sum$chain_G!="TRG","B",
                                                                                                                                 ifelse(junction_only_sum$chain_B!="TRB" & junction_only_sum$chain_A!="TRA" & junction_only_sum$chain_D!="TRD" & junction_only_sum$chain_G=="TRG","G",
                                                                                                                                        ifelse(junction_only_sum$chain_B!="TRB" & junction_only_sum$chain_A!="TRA" & junction_only_sum$chain_D=="TRD" & junction_only_sum$chain_G!="TRG","D",
                                                                                                                                               "other"
                                                                                                                                        )))))))))))))))

        head(junction_only_sum)
        junction_only_sum_selected <- subset(junction_only_sum,junction_only_sum$junction_aa_A=="CVLSASSSFSKLVF")
        head(junction_only_sum_selected)
        merged_junction_contig <-  merge(contig_paired_only,junction_only_sum, by = c("junction_aa_A","junction_aa_B","junction_aa_G","junction_aa_D","v_gene_A","v_gene_B","v_gene_G","v_gene_D","chain_A","chain_B","chain_G","chain_D","pairing")) # "v_gene_B","v_gene_D","v_gene_G"
        name.list_all <- c("cell_id","Sample_Name","pairing","order")
        merged_junction_contig <- merged_junction_contig %>%
          select(all_of(name.list_all), everything())

        # merged_junction_contig <- merged_junction_contig[!names(merged_junction_contig) %in% c(names(merged_junction_contig[grep("v_gene",names(merged_junction_contig))]))]
        merged_junction_contig <- merged_junction_contig[order(merged_junction_contig$cell_id),]
        merged_junction_contig
      }

      else {
        as.data.frame("Upload unfilted AIRR file")
      }

    })

    output$tb_multiTCR <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      multi_TCR_BDrhap()
    })

    output$downloadtb_multiTCR <- downloadHandler(
      filename = function(){
        paste(input$name.BD,"_Multi_TCR_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(multi_TCR_BDrhap())
        write.csv(df,file, row.names = F)
      } )

    # 10x_Genomics data -----
    ## barcode file -----
    input.data.barcode.10x <- reactive({
      inFile_10x_barcode <- input$file_calls_10x
      if (is.null(inFile_10x_barcode)) return(NULL)

      else {
        dataframe <- read.table(
          inFile_10x_barcode$datapath)}
    })

    output$test.files.10x1 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- input.data.barcode.10x()
      validate(
        need(nrow(calls)>0,
             error_message_val_10x_barcode)
      )
      calls
    })
    ## features file -----
    input.data.features.10x <- reactive({
      inFile_10x_features <- input$file_features_10x
      if (is.null(inFile_10x_features)) return(NULL)

      else {
        dataframe <- read.table(
          inFile_10x_features$datapath)}
    })

    output$test.files.10x2 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- input.data.features.10x()
      validate(
        need(nrow(calls)>0,
             error_message_val_10x_features)
      )
      calls
    })

    # Matrix file
    input.data.matrix.10x <- reactive({
      inFile_10x_matrix <- input$file_matrix_10x
      if (is.null(inFile_10x_matrix)) return(NULL)

      else {
        dataframe <- Matrix::readMM(inFile_10x_matrix$datapath)}
    })
    output$test.files.10x3 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- as.data.frame(input.data.matrix.10x())
      validate(
        need(nrow(calls)>0,
             error_message_val_10x_features)
      )
      calls[1:10,1:10]
    })

    ## contig files ----
    input.data.TCR.10x <- reactive({

      inFile_10x_TCR <- input$file_TCR_10x
      if (is.null(inFile_10x_TCR)) return(NULL)

      else {
        if (input$csv_contig_file == "csv/csv.gz") {
          dataframe <- read.csv(inFile_10x_TCR$datapath, na.strings = c("", "NA"))
        }
        else {
          dataframe <- read.table(inFile_10x_TCR$datapath,sep = "\t",header =T)
        }
        dataframe
      }
    })
    output$test.files.10x4 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      if (input$BCR_TCR_10x == "TCR only") {
        calls <- input.data.TCR.10x()
        validate(
          need(nrow(calls)>0,
               error_message_val_10x_features)
        )
        tb_10x_meta.data_TCR()
      }

      else {
        calls <- input.data.TCR.10x()
        validate(
          need(nrow(calls)>0,
               error_message_val_10x_features)
        )
        tb_10x_meta.data_BCR()
      }

    })
    ## meta.data for seurat ----
    tb_10x_for.tcr.exploreAIRR <- reactive({
      contigs <- input.data.TCR.10x()
      validate(
        need(nrow(contigs)>0,
             "Upload file")
      )
      contigs
      names(contigs) <- gsub("call","gene",names(contigs))

      names(contigs) <- gsub("junction_aa","cdr3",names(contigs))
      contigs
      contigs_lim <- contigs[!names(contigs) %in% c("is_cell","contig_id","high_confidence","raw_consensus_id","exact_subclonotype_id","umis","reads","length","cdr3_nt","germline_alignment","cdr3_nt_id","cdr3_nt_alignment","rev_comp",names(contigs[grep("_end",names(contigs))]),names(contigs[grep("cigar",names(contigs))]),names(contigs[grep("length",names(contigs))]),names(contigs[grep("count",names(contigs))]),names(contigs[grep("sequence",names(contigs))]),names(contigs[grep("fwr",names(contigs))]),names(contigs[grep("cdr1",names(contigs))]),names(contigs[grep("cdr2",names(contigs))])
      )]


      contigs_lim$chain <- ifelse(grepl("TRA",contigs_lim$v_gene),"TRA",
                                  ifelse(grepl("TRB",contigs_lim$v_gene),"TRB",
                                         ifelse(grepl("TRG",contigs_lim$v_gene),"TRG",
                                                ifelse(grepl("TRD",contigs_lim$v_gene),"TRD",
                                                       ""
                                                ))))

      contigs_lim
      names(contigs_lim) <- gsub("junction","cdr3_nt",names(contigs_lim))

      contig_AG <- subset(contigs_lim,contigs_lim$chain=="TRA" | contigs_lim$chain=="TRG")
      head(contig_AG)
      name.list <- names(contig_AG[c(names(contig_AG[grep("gene",names(contig_AG))]),
                                     names(contig_AG[grep("cdr3",names(contig_AG))]),

                                     "chain")])
      name.list

      contig_AG <- contig_AG %>%
        select(all_of(name.list), everything())
      names(contig_AG)[1:summary(name.list)[1]] <-paste(names(contig_AG[names(contig_AG) %in% name.list]),"_AG",sep="")
      contig_AG

      contig_BD <- subset(contigs_lim,contigs_lim$chain=="TRB" | contigs_lim$chain=="TRD")
      name.list <- names(contig_BD[c(names(contig_BD[grep("gene",names(contig_BD))]),
                                     names(contig_BD[grep("cdr3",names(contig_BD))]),
                                     "chain")])
      contig_BD <- contig_BD %>%
        select(all_of(name.list), everything())


      names(contig_BD)[1:summary(name.list)[1]] <-paste(names(contig_BD[names(contig_BD) %in% name.list]),"_BD",sep="")
      contig_BD
      # contig_paired <- merge(contig_AG,contig_BD, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      # contig_paired <- merge(contig_AG,contig_BD, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      contig_paired <- merge(contig_AG,contig_BD, by=c("cell_id", "clone_id"),all = T)

      contig_paired$pairing <- ifelse(contig_paired$chain_BD=="TRB" & contig_paired$chain_AG=="TRA","abTCR Paired",
                                      ifelse(contig_paired$chain_BD=="TRD" & contig_paired$chain_AG=="TRG","gdTCR Paired",NA
                                      ))
      contig_paired
      contig_paired$pairing[is.na(contig_paired$pairing)] <- "unpaired"
      contig_paired <- contig_paired[!names(contig_paired) %in% c("d_gene_AG")]
      contig_paired_only <- contig_paired
      contig_paired_only <- subset(contig_paired_only,contig_paired_only$cdr3_BD!="None")
      contig_paired_only <- subset(contig_paired_only,contig_paired_only$cdr3_AG!="None")
      contig_paired_only$d_gene_BD <- gsub("^$","NA", contig_paired_only$d_gene_BD)
      #
      contig_paired_only$vj_gene_AG <- paste(contig_paired_only$v_gene_AG,contig_paired_only$j_gene_AG,sep = ".")
      contig_paired_only$vj_gene_AG <- gsub("NA.NA","",contig_paired_only$vj_gene_AG)
      #
      contig_paired_only$vj_gene_BD <- paste(contig_paired_only$v_gene_BD,contig_paired_only$j_gene_BD,sep = ".")
      contig_paired_only$vj_gene_BD <- gsub(".NA.",".",contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub(".None.",".",contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub("NA.NA","",contig_paired_only$vj_gene_BD)
      #
      contig_paired_only$vdj_gene_BD <- paste(contig_paired_only$v_gene_BD,contig_paired_only$d_gene_BD,contig_paired_only$j_gene_BD,sep = ".")
      contig_paired_only$vdj_gene_BD <- gsub(".NA.",".",contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("[.]NA[.]",".",contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("NA.NA","",contig_paired_only$vdj_gene_BD)
      #
      contig_paired_only$vj_gene_cdr3_AG <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$cdr3_AG,sep = "_")
      contig_paired_only$vj_gene_cdr3_AG <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_AG)
      #
      contig_paired_only$vj_gene_cdr3_BD <- paste(contig_paired_only$vj_gene_BD,contig_paired_only$cdr3_BD,sep = "_")
      contig_paired_only$vj_gene_cdr3_BD <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_BD)
      #
      contig_paired_only$vdj_gene_cdr3_BD <- paste(contig_paired_only$vdj_gene_BD,contig_paired_only$cdr3_BD,sep = "_")
      contig_paired_only$vdj_gene_cdr3_BD <- gsub("_NA","",contig_paired_only$vdj_gene_cdr3_BD)
      #
      contig_paired_only$vj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$vj_gene_BD,sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$vdj_gene_BD,sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- gsub("^ & ","",contig_paired_only$vdj_gene_AG_BD)
      contig_paired_only$vdj_gene_AG_BD <- gsub(" & $","",contig_paired_only$vdj_gene_AG_BD)
      #
      # #updating names to be consistant....
      contig_paired_only$vj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vj_gene_cdr3_BD,sep = " & ")
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub("^ & ","",contig_paired_only$vj_gene_cdr3_AG_BD)
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub(" & $","",contig_paired_only$vj_gene_cdr3_AG_BD)

      contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vdj_gene_cdr3_BD,sep = " & ")
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub("^ & ","",contig_paired_only$vdj_gene_cdr3_AG_BD)
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub(" & $","",contig_paired_only$vdj_gene_cdr3_AG_BD)


      names(contig_paired_only)[names(contig_paired_only) %in% "cell_id"] <- "Cell_Index"

      contig_paired_only$Sample_Name <- paste(input$group10x,"_",input$Indiv10x,sep ="")

      contig_paired_only <- contig_paired_only %>%
        select(all_of(c("Cell_Index","Sample_Name")), everything())
      contig_paired_only

    })
    tb_10x_meta.data_TCR_CSV <- function () {
      contigs <- input.data.TCR.10x()
      validate(
        need(nrow(contigs)>0,
             "Upload file")
      )
      contigs_lim <- contigs[!names(contigs) %in% c("is_cell","contig_id","high_confidence","raw_consensus_id","exact_subclonotype_id","umis","reads","length","cdr3_nt",names(contigs[grep("fwr",names(contigs))]),names(contigs[grep("cdr1",names(contigs))]),names(contigs[grep("cdr2",names(contigs))])
      )]
      contigs_lim
      contig_AG <- subset(contigs_lim,contigs_lim$chain=="TRA" | contigs_lim$chain=="TRG")
      name.list <- names(contig_AG[c(names(contig_AG[grep("gene",names(contig_AG))]),
                                     names(contig_AG[grep("cdr3",names(contig_AG))]),
                                     "chain")])
      contig_AG <- contig_AG %>%
        select(all_of(name.list), everything())
      names(contig_AG)[1:summary(name.list)[1]] <-paste(names(contig_AG[names(contig_AG) %in% name.list]),"_AG",sep="")
      contig_AG

      contig_BD <- subset(contigs_lim,contigs_lim$chain=="TRB" | contigs_lim$chain=="TRD")
      name.list <- names(contig_BD[c(names(contig_BD[grep("gene",names(contig_BD))]),
                                     names(contig_BD[grep("cdr3",names(contig_BD))]),
                                     "chain")])
      contig_BD <- contig_BD %>%
        select(all_of(name.list), everything())


      names(contig_BD)[1:summary(name.list)[1]] <-paste(names(contig_BD[names(contig_BD) %in% name.list]),"_BD",sep="")
      contig_BD
      # contig_paired <- merge(contig_AG,contig_BD, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      # contig_paired <- merge(contig_AG,contig_BD, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      contig_paired <- merge(contig_AG,contig_BD, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)

      contig_paired$pairing <- ifelse(contig_paired$chain_BD=="TRB" & contig_paired$chain_AG=="TRA","abTCR Paired",
                                      ifelse(contig_paired$chain_BD=="TRD" & contig_paired$chain_AG=="TRG","gdTCR Paired",NA
                                      ))
      contig_paired
      contig_paired$pairing[is.na(contig_paired$pairing)] <- "unpaired"
      contig_paired <- contig_paired[!names(contig_paired) %in% c("d_gene_AG")]
      contig_paired_only <- contig_paired
      contig_paired_only <- subset(contig_paired_only,contig_paired_only$cdr3_BD!="None")
      contig_paired_only <- subset(contig_paired_only,contig_paired_only$cdr3_AG!="None")
      contig_paired_only$d_gene_BD <- gsub("^$","NA", contig_paired_only$d_gene_BD)
      #
      contig_paired_only$vj_gene_AG <- paste(contig_paired_only$v_gene_AG,contig_paired_only$j_gene_AG,sep = ".")
      contig_paired_only$vj_gene_AG <- gsub("NA.NA","",contig_paired_only$vj_gene_AG)
      #
      contig_paired_only$vj_gene_BD <- paste(contig_paired_only$v_gene_BD,contig_paired_only$j_gene_BD,sep = ".")
      contig_paired_only$vj_gene_BD <- gsub(".NA.",".",contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub(".None.",".",contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub("NA.NA","",contig_paired_only$vj_gene_BD)
      #
      contig_paired_only$vdj_gene_BD <- paste(contig_paired_only$v_gene_BD,contig_paired_only$d_gene_BD,contig_paired_only$j_gene_BD,sep = ".")
      contig_paired_only$vdj_gene_BD <- gsub(".NA.",".",contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("[.]NA[.]",".",contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("NA.NA","",contig_paired_only$vdj_gene_BD)
      #
      contig_paired_only$vj_gene_cdr3_AG <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$cdr3_AG,sep = "_")
      contig_paired_only$vj_gene_cdr3_AG <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_AG)
      #
      contig_paired_only$vj_gene_cdr3_BD <- paste(contig_paired_only$vj_gene_BD,contig_paired_only$cdr3_BD,sep = "_")
      contig_paired_only$vj_gene_cdr3_BD <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_BD)
      #
      contig_paired_only$vdj_gene_cdr3_BD <- paste(contig_paired_only$vdj_gene_BD,contig_paired_only$cdr3_BD,sep = "_")
      contig_paired_only$vdj_gene_cdr3_BD <- gsub("_NA","",contig_paired_only$vdj_gene_cdr3_BD)
      #
      contig_paired_only$vj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$vj_gene_BD,sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$vdj_gene_BD,sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- gsub("^ & ","",contig_paired_only$vdj_gene_AG_BD)
      contig_paired_only$vdj_gene_AG_BD <- gsub(" & $","",contig_paired_only$vdj_gene_AG_BD)
      #
      # #updating names to be consistant....
      contig_paired_only$vj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vj_gene_cdr3_BD,sep = " & ")
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub("^ & ","",contig_paired_only$vj_gene_cdr3_AG_BD)
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub(" & $","",contig_paired_only$vj_gene_cdr3_AG_BD)

      contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vdj_gene_cdr3_BD,sep = " & ")
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub("^ & ","",contig_paired_only$vdj_gene_cdr3_AG_BD)
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub(" & $","",contig_paired_only$vdj_gene_cdr3_AG_BD)
      # contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vdj_gene_cdr3_BD,sep = " & ")
      names(contig_paired_only)[names(contig_paired_only) %in% "barcode"] <- "Cell_Index"
      contig_paired_only <- contig_paired_only[!duplicated(contig_paired_only$Cell_Index),] # remove duplicates
      contig_paired_only$Sample_Name <- paste(input$group10x,"_",input$Indiv10x,sep ="")

      contig_paired_only <- contig_paired_only %>%
        select(all_of(c("Cell_Index","Sample_Name")), everything())

      contig_paired_only
    }
    tb_10x_meta.data_TCR <- reactive({
      if (input$csv_contig_file == "csv/csv.gz") {
        dataframe <- tb_10x_meta.data_TCR_CSV()
      }
      else {
        dataframe <- tb_10x_for.tcr.exploreAIRR()
      }
      dataframe
    })


    tb_10x_meta.data_BCR <- function () {
      contigs <- input.data.TCR.10x()
      validate(
        need(nrow(contigs)>0,
             "Upload file")
      )
      contigs <- contigs[order(contigs$umis,decreasing = T),]

      contigs_lim <- contigs[!names(contigs) %in% c("is_cell","contig_id","high_confidence","raw_consensus_id","exact_subclonotype_id","reads","length","cdr3_nt",names(contigs[grep("fwr",names(contigs))]),names(contigs[grep("cdr1",names(contigs))]),names(contigs[grep("cdr2",names(contigs))])
      )]
      contigs_lim
      contig_LK <- subset(contigs_lim,contigs_lim$chain=="IGL" | contigs_lim$chain=="IGK")
      name.list <- names(contig_LK[c(names(contig_LK[grep("gene",names(contig_LK))]),
                                     names(contig_LK[grep("cdr3",names(contig_LK))]),
                                     "chain")])
      name.list
      contig_LK <- contig_LK %>%
        select(all_of(name.list), everything())
      names(contig_LK)[1:summary(name.list)[1]] <-paste(names(contig_LK[names(contig_LK) %in% name.list]),"_IgL",sep="")
      head(contig_LK)

      contig_H <- subset(contigs_lim,contigs_lim$chain=="IGH")

      name.list <- names(contig_H[c(names(contig_H[grep("gene",names(contig_H))]),
                                    names(contig_H[grep("cdr3",names(contig_H))]),
                                    "chain")])
      contig_H <- contig_H %>%
        select(all_of(name.list), everything())


      names(contig_H)[1:summary(name.list)[1]] <-paste(names(contig_H[names(contig_H) %in% name.list]),"_IgH",sep="")
      head(contig_H)
      # contig_paired <- merge(contig_LK,contig_H, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      # contig_paired <- merge(contig_LK,contig_H, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      contig_paired <- merge(contig_LK,contig_H, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)

      contig_paired$pairing <- ifelse(contig_paired$chain_IgH=="IGH" & contig_paired$chain_IgL=="IGK","IGK Paired",
                                      ifelse(contig_paired$chain_IgH=="IGH" & contig_paired$chain_IgL=="IGL","IGL Paired",NA
                                      ))

      contig_paired
      contig_paired$pairing[is.na(contig_paired$pairing)] <- "unpaired"
      contig_paired <- contig_paired[!names(contig_paired) %in% c("d_gene_IgL")]
      contig_paired_only <- contig_paired
      contig_paired_only <- subset(contig_paired_only,contig_paired_only$cdr3_IgH!="None")
      contig_paired_only <- subset(contig_paired_only,contig_paired_only$cdr3_IgL!="None")
      dim(contig_paired_only)

      contig_paired_only$d_gene_IgH <- sub("^$","NA", contig_paired_only$d_gene_IgH)
      #
      contig_paired_only$vj_gene_LK <- paste(contig_paired_only$v_gene_IgL,contig_paired_only$j_gene_IgL,sep = ".")
      contig_paired_only$vj_gene_LK <- gsub("NA.NA","",contig_paired_only$vj_gene_LK)
      #
      contig_paired_only$vj_gene_H <- paste(contig_paired_only$v_gene_IgH,contig_paired_only$j_gene_IgH,sep = ".")
      contig_paired_only$vj_gene_H <- gsub(".NA.",".",contig_paired_only$vj_gene_H)
      contig_paired_only$vj_gene_H <- gsub(".None.",".",contig_paired_only$vj_gene_H)
      contig_paired_only$vj_gene_H <- gsub("NA.NA","",contig_paired_only$vj_gene_H)
      #
      contig_paired_only$vdj_gene_H <- paste(contig_paired_only$v_gene_IgH,contig_paired_only$d_gene_IgH,contig_paired_only$j_gene_IgH,sep = ".")
      contig_paired_only$vdj_gene_H <- gsub(".NA.",".",contig_paired_only$vdj_gene_H)
      contig_paired_only$vdj_gene_H <- gsub(".None.",".",contig_paired_only$vdj_gene_H)
      contig_paired_only$vdj_gene_H <- gsub("NA.NA","",contig_paired_only$vdj_gene_H)
      #
      contig_paired_only$vj_gene_cdr3_LK <- paste(contig_paired_only$vj_gene_LK,contig_paired_only$cdr3_IgL,sep = "_")
      contig_paired_only$vj_gene_cdr3_LK <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_LK)
      #
      contig_paired_only$vj_gene_cdr3_H <- paste(contig_paired_only$vj_gene_H,contig_paired_only$cdr3_IgH,sep = "_")
      contig_paired_only$vj_gene_cdr3_H <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_H)
      #
      contig_paired_only$vdj_gene_cdr3_H <- paste(contig_paired_only$vdj_gene_H,contig_paired_only$cdr3_IgH,sep = "_")
      contig_paired_only$vdj_gene_cdr3_H <- gsub("_NA","",contig_paired_only$vdj_gene_cdr3_H)
      #
      contig_paired_only$vj_gene_LK_H <- paste(contig_paired_only$vj_gene_LK,contig_paired_only$vj_gene_H,sep = " & ")
      contig_paired_only$vdj_gene_LK_H <- paste(contig_paired_only$vj_gene_LK,contig_paired_only$vdj_gene_H,sep = " & ")
      contig_paired_only$vdj_gene_LK_H <- gsub("^ & ","",contig_paired_only$vdj_gene_LK_H)
      contig_paired_only$vdj_gene_LK_H <- gsub(" & $","",contig_paired_only$vdj_gene_LK_H)
      #
      # #updating names to be consistant....
      contig_paired_only$vj_gene_cdr3_LK_H <- paste(contig_paired_only$vj_gene_cdr3_LK,contig_paired_only$vj_gene_cdr3_H,sep = " & ")
      contig_paired_only$vj_gene_cdr3_LK_H <- gsub("^ & ","",contig_paired_only$vj_gene_cdr3_LK_H)
      contig_paired_only$vj_gene_cdr3_LK_H <- gsub(" & $","",contig_paired_only$vj_gene_cdr3_LK_H)

      contig_paired_only$vdj_gene_cdr3_LK_H <- paste(contig_paired_only$vj_gene_cdr3_LK,contig_paired_only$vdj_gene_cdr3_H,sep = " & ")
      contig_paired_only$vdj_gene_cdr3_LK_H <- gsub("^ & ","",contig_paired_only$vdj_gene_cdr3_LK_H)
      contig_paired_only$vdj_gene_cdr3_LK_H <- gsub(" & $","",contig_paired_only$vdj_gene_cdr3_LK_H)
      # contig_paired_only$vdj_gene_cdr3_LK_H <- paste(contig_paired_only$vj_gene_cdr3_LK,contig_paired_only$vdj_gene_cdr3_H,sep = " & ")
      names(contig_paired_only)[names(contig_paired_only) %in% "barcode"] <- "Cell_Index"
      dup <- contig_paired_only[duplicated(contig_paired_only$Cell_Index),]
      contig_paired_only <- contig_paired_only[order(contig_paired_only$Cell_Index, contig_paired_only$umis.x,contig_paired_only$umis.y,decreasing = T),]

      contig_paired_only_dup <- contig_paired_only[!duplicated(contig_paired_only$Cell_Index),] # remove duplicate barcodes.
      names(contig_paired_only_dup)
      contig_paired_only_dup <- contig_paired_only_dup[!names(contig_paired_only_dup) %in% c("umis.x","umis.y")]
      contig_paired_only_dup$Sample_Name <- paste(input$group10x,"_",input$Indiv10x,sep ="")

      contig_paired_only_dup <- contig_paired_only_dup %>%
        select(all_of(c("Cell_Index","Sample_Name")), everything())

      contig_paired_only_dup

      contig_paired_only_dup
    }

    output$tb_10x_meta1 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      if (input$BCR_TCR_10x=="TCR only") {
        tb_10x_meta.data_TCR()
      }
      else if (input$BCR_TCR_10x=="BCR only") {
        tb_10x_meta.data_BCR()
      }
      else {
        # TCR <- tb_10x_meta.data_TCR()
        # BCR <- tb_10x_meta.data_BCR()
        # contig_paired <- names(BCR)[!grepl("_IgH",names(BCR)) & !grepl("_IgL",names(BCR))]
        # contig_paired <- merge(TCR,BCR,by = merge.names, all=T)
        # contig_paired
      }

    })
    output$downloadtb_10x_metadata2 <- downloadHandler(
      filename = function(){
        paste(input$group10x,"_",input$Indiv10x,"_metadata_10x_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){

        if (input$BCR_TCR_10x=="TCR only") {
          df <-  tb_10x_meta.data_TCR()

        }

        else if (input$BCR_TCR_10x=="BCR only") {
          df <-  tb_10x_meta.data_BCR()
        }
        else {

        }
        # write.table(,file, row.names = T)
        write.csv(df, file,row.names = F)

      } )

    # meta-data summary ----
    sum_tb_10x <- function () {
      contigs <- input.data.TCR.10x()
      validate(
        need(nrow(contigs)>0,
             "Upload file")
      )

      if (input$BCR_TCR_10x=="TCR only") {
        contig_paired <-  tb_10x_meta.data_TCR()

      }

      else if (input$BCR_TCR_10x=="BCR only") {
        contig_paired <-  tb_10x_meta.data_BCR()
      }
      else {


        TCR <- tb_10x_meta.data_TCR()
        BCR <- tb_10x_meta.data_BCR()
        contig_paired <- names(BCR)[!grepl("_IgH",names(BCR)) & !grepl("_IgL",names(BCR))]
        contig_paired <- merge(TCR,BCR,by = merge.names, all=T)

      }

      # contig_paired <- as.data.frame(tb_10x_meta.data())

      count.df <- contig_paired[names(contig_paired) %in% c(names(contig_paired)[grepl("chain",names(contig_paired))],"pairing")]
      count.df$count <- 1
      # ddply(count.df,names(count.df) ,numcolwise(sum))
      df1 <- ddply(count.df,names(count.df)[c(-4)] ,numcolwise(sum))
      df1

    }
    output$sum_tb_10x1 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      sum_tb_10x()
    })

    ## ClusTCR T cell receptor sequences  -----
    tb_10x_contigues_contig <- reactive({
      contigs <- input.data.TCR.10x()
      validate(
        need(nrow(contigs)>0,
             error_message_val_10x_features)
      )
      if (input$csv_contig_file == "tsv") {
        names(contigs) <- gsub("call","gene",names(contigs))
        names(contigs) <- gsub("junction_aa","cdr3",names(contigs))
      }
      else {

      }

      contigs2 <- contigs[names(contigs) %in% c("v_gene","cdr3")]
      contigs2$Sample_Name <-  paste(input$group10x,"_",input$Indiv10x, sep = "")
      names(contigs2)[names(contigs2) %in% c("cdr3")] <- "junction_aa"
      names(contigs2)[names(contigs2) %in% c("v_gene")] <- "v_call"
      if (nrow(contigs2[-c(grep("[*]",contigs2$junction_aa)),]>0)) {
        contigs2 <- contigs2[-c(grep("[*]",contigs2$junction_aa)),]
      }

      contigs2$junction_aa <- gsub("^$","None", contigs2$junction_aa)

      if (nrow(contigs2[-c(grep("None",contigs2$junction_aa)),]>0)) {
        contigs2 <- contigs2[-c(grep("None",contigs2$junction_aa)),]
      }

      df <- contigs2[!duplicated(contigs2[,c('v_call','junction_aa')]),]

      if (input$chain_clusTCR2_10x == "AG") {
        (rbind(df[grep("TRAV",df$v_call),],df[grep("TRGV",df$v_call),]))

      }

      else if (input$chain_clusTCR2_10x == "IgH") {
        rbind(df[grep("IGH",df$v_call),])

      }

      else if (input$chain_clusTCR2_10x == "IgLK") {
        (rbind(df[grep("IGL",df$v_call),],df[grep("IGK",df$v_call),]))

      }

      else {
        (rbind(df[grep("TRBV",df$v_call),],df[grep("TRDV",df$v_call),]))
      }

    })

    output$tb_10x_contigues1 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      tb_10x_contigues_contig()
    })

    output$downloadtb_10x_contigues1 <- downloadHandler(
      filename = function(){
        x = today()
        paste(input$chain_clusTCR2_10x,"_",input$group10x,"_",input$Indiv10x,"_clusTCR_10x_",x,".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(tb_10x_contigues_contig())
        write.csv(df,file, row.names = F)
      } )


    ## TCR explore 10x -----
    TCR_Explore_10x <- function () {
      if (input$BCR_TCR_10x=="TCR only") {
        if (input$csv_contig_file == "csv/csv.gz") {

          contig_paired_only <-  tb_10x_meta.data_TCR()


        }
        else {
          contig_paired_only <-  tb_10x_for.tcr.exploreAIRR()
        }


      }

      else if (input$BCR_TCR_10x=="BCR only") {
        contig_paired_only <-  tb_10x_meta.data_BCR()
      }
      else {


        TCR <- tb_10x_meta.data_TCR()
        BCR <- tb_10x_meta.data_BCR()
        contig_paired_only <- names(BCR)[!grepl("_IgH",names(BCR)) & !grepl("_IgL",names(BCR))]
        contig_paired_only <- merge(TCR,BCR,by = merge.names, all=T)

      }

      contig_paired_only$cloneCount <- 1
      contig_paired_only$group <- input$group10x
      contig_paired_only$Indiv <- input$Indiv10x
      contig_paired_only <-  contig_paired_only %>%
        select(cloneCount,group, Indiv, everything())
      contig_paired_only

    }

    output$dt_TCR_Explore_10x <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      TCR_Explore_10x()
    })

    output$downloaddt_TCR_Explore_10x <- downloadHandler(
      filename = function(){
        paste(input$group10x,"_",input$Indiv10x,"_TCR_Explore_10x_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(TCR_Explore_10x())
        write.csv(df,file, row.names = F)
      } )

    ## count table ----
    tb_10x_matrix <- function () {
      barcode <- input.data.barcode.10x()
      features <- input.data.features.10x()
      mmMat <- as.data.frame(input.data.matrix.10x())

      validate(
        need(nrow(barcode)>0 & nrow(features)>0 & nrow(mmMat)>0,
             "Upload files")
      )

      rownames(mmMat) <- make.unique(features$V2) # select which feature to label genes...
      names(mmMat) <- barcode$V1
      mmMat$Gene_Name <- rownames(mmMat)
      mmMat <- mmMat %>%
        select(all_of("Gene_Name"), everything())
      mmMat

    }

    tb_10x_matrix_h5 <- function () {
      mmMat <- as.data.frame(input.data.h5.10x())

      validate(
        need(nrow(mmMat)>0,
             "Upload h5 files")
      )
      mmMat
    }

    output$tb_10x_matrix2 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      tb_10x_matrix()[1:6,1:6]
    })

    output$downloadtb_10x_matrix2 <- downloadHandler(
      filename = function(){
        x = today()
        paste(input$group10x,"_",input$Indiv10x,"_count-matrix_10x_",x,".csv.gz", sep = "")
      },
      content = function(file){
        df <- as.data.frame(tb_10x_matrix())
        # write.table(,file, row.names = T)
        write_csv(df,gzfile(file))
      })

    # for TCRex output -----
    TCRex_10x_df <- function () {
      contigs <- input.data.TCR.10x()
      validate(
        need(nrow(contigs)>0,
             error_message_val_10x_features)
      )

      if (input$csv_contig_file == "tsv") {
        names(contigs) <- gsub("call","gene",names(contigs))
        names(contigs) <- gsub("junction_aa","cdr3",names(contigs))
      }
      else {

      }

      contigs2 <- contigs[,names(contigs) %in% c("v_gene","j_gene","cdr3")]

      names(contigs2)[names(contigs2) %in% c("cdr3")] <- "CDR3_beta"
      if (nrow(contigs2[-c(grep("[*]",contigs2$CDR3_beta)),]>0)) {
        contigs2 <- contigs2[-c(grep("[*]",contigs2$CDR3_beta)),]
      }
      #

      names(contigs2)[names(contigs2) %in% c("v_gene")] <- "TRBV_gene"
      names(contigs2)[names(contigs2) %in% c("j_gene")] <- "TRBJ_gene"
      contigs2 <- contigs2 %>%
        select(TRBV_gene,CDR3_beta,everything())


      contigs2$CDR3_beta <- gsub("^$","None", contigs2$CDR3_beta)

      if (nrow(contigs2[-c(grep("None",contigs2$CDR3_beta)),]>0)) {
        contigs2 <- subset(contigs2,contigs2$CDR3_beta!= "None")
      }

      contigs2$TRBJ_gene <- gsub("^$","None", contigs2$TRBJ_gene)

      if (nrow(contigs2[-c(grep("None",contigs2$TRBJ_gene)),]>0)) {
        contigs2 <- subset(contigs2,contigs2$TRBJ_gene!= "None")
      }


      contigs2 <- contigs2
      contigs2$TRBV_gene <- gsub("[*]0.","",contigs2$TRBV_gene)
      contigs2$TRBJ_gene <- gsub("[*]0.","",contigs2$TRBJ_gene)
      contigs2$cloneCount <- 1
      calls_TCR_paired.fun3 <- ddply(contigs2,names(contigs2)[-c(4)] ,numcolwise(sum))
      calls_TCR_paired.fun3 <- calls_TCR_paired.fun3[grepl("TRB", calls_TCR_paired.fun3$TRBV_gene),]
      calls_TCR_paired.fun3
    }

    output$tb_TCRex_10x_df <- DT::renderDataTable(filter = list(position = 'top', clear = FALSE), escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      TCRex_10x_df()
    })

    output$downloaddf_TCRex_10x <- downloadHandler(
      filename = function(){
        x = today()
        paste( input$group10x,"_",input$Indiv10x,"_TCRex_",x,".tsv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(TCRex_10x_df())
        df <- df[,names(df) %in% c("TRBV_gene","CDR3_beta","TRBJ_gene")]
        write.table(df,file, row.names = F,sep="\t", quote = F)
      } )




    # Array ------
    input.data.calls.array <- reactive({
      inFile_arrayM <- input$file_calls_Array
      if (is.null(inFile_arrayM)) return(NULL)
      else {
        dataframe <- read.table(inFile_arrayM$datapath)
      }
    })
    output$test.files_array_Matrix <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- input.data.calls.array()
      validate(
        need(nrow(calls)>0,
             "Upload Matrix file")
      )
      head(calls)[1:6]
    })

    input.data.calls.array_contig <- reactive({
      inFile_arrayC <- input$file_contig_Array
      if (is.null(inFile_arrayC)) return(NULL)
      else {
        dataframe <- read.table(inFile_arrayC$datapath,sep ="\t",header = T, row.names = 1)
      }

    })
    output$test.files_array_contig <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- input.data.calls.array_contig()
      validate(
        need(nrow(calls)>0,
             "Upload Contig/TCR file")
      )
      calls
    })


    Filtering_TCR_Array <- reactive({
      df <- input.data.calls.array_contig()
      df$Cell_Index <- rownames(df)
      df_filtering <- subset(df,df$clonotype!="None")

      if(input$pairing_TCR_Array==T) {
        # "Single pair"     "None"            "Orphan beta"     "Extra alpha"     "Extra beta"      "Orphan alpha"    "Two full chains"
        df_filtering <- subset(df_filtering,df_filtering$chain_pairing!="Orphan beta")
        df_filtering <- subset(df_filtering,df_filtering$chain_pairing!="Orphan alpha")
      }

      else {
        df_filtering
      }
      df_filtering <-  df_filtering[,!grepl("_2_",names(df_filtering))]
      df_filtering <-  df_filtering[,!grepl("_junction_ins",names(df_filtering))]
      df_filtering <-  df_filtering[,!grepl("_expr",names(df_filtering))]
      df_filtering <-  df_filtering[,!grepl("clon",names(df_filtering))]
      df_filtering <-  df_filtering[,!grepl("multi_chain",names(df_filtering))]
      df_filtering <-  df_filtering[,!grepl("TRA_1_d_gene",names(df_filtering))]
      # df_filtering <-  df_filtering[,c(!grepl("multi_chain",names(df_filtering),"" ))]
      calls_TCR_paired.fun <- df_filtering
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TRA_1_cdr3")] <- "cdr3_AG"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TRA_1_j_gene")] <- "j_gene_AG"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TRA_1_v_gene")] <- "v_gene_AG"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TRA_1_c_gene")] <- "c_gene_AG"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TRB_1_cdr3")] <- "cdr3_BD"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TRB_1_d_gene")] <- "d_gene_BD"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TRB_1_j_gene")] <- "j_gene_BD"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TRB_1_v_gene")] <- "v_gene_BD"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TRB_1_c_gene")] <- "c_gene_BD"
      contig_paired_only <- calls_TCR_paired.fun

      contig_paired_only$vj_gene_AG <- paste(contig_paired_only$v_gene_AG,contig_paired_only$j_gene_AG,sep = ".")
      contig_paired_only$vj_gene_AG <- gsub("NA.NA","",contig_paired_only$vj_gene_AG)
      contig_paired_only$vj_gene_AG <- gsub("None.None","",contig_paired_only$vj_gene_AG)
      #
      contig_paired_only$vj_gene_BD <- paste(contig_paired_only$v_gene_BD,contig_paired_only$j_gene_BD,sep = ".")
      contig_paired_only$vj_gene_BD <- gsub(".NA.",".",contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub(".None.",".",contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub("NA.NA","",contig_paired_only$vj_gene_BD)
      #
      contig_paired_only$vdj_gene_BD <- paste(contig_paired_only$v_gene_BD,contig_paired_only$d_gene_BD,contig_paired_only$j_gene_BD,sep = ".")
      contig_paired_only$vdj_gene_BD <- gsub(".NA.",".",contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub(".None.",".",contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("NA.NA","",contig_paired_only$vdj_gene_BD)
      #
      contig_paired_only$vj_gene_cdr3_AG <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$cdr3_AG,sep = "_")
      contig_paired_only$vj_gene_cdr3_AG <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_AG)
      contig_paired_only$vj_gene_cdr3_AG <- gsub("_None","",contig_paired_only$vj_gene_cdr3_AG)
      #
      contig_paired_only$vj_gene_cdr3_BD <- paste(contig_paired_only$vj_gene_BD,contig_paired_only$cdr3_BD,sep = "_")
      contig_paired_only$vj_gene_cdr3_BD <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_BD)
      #
      contig_paired_only$vdj_gene_cdr3_BD <- paste(contig_paired_only$vdj_gene_BD,contig_paired_only$cdr3_BD,sep = "_")
      contig_paired_only$vdj_gene_cdr3_BD <- gsub("_NA","",contig_paired_only$vdj_gene_cdr3_BD)
      #
      contig_paired_only$vj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$vj_gene_BD,sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG,contig_paired_only$vdj_gene_BD,sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- gsub("^ & ","",contig_paired_only$vdj_gene_AG_BD)
      contig_paired_only$vdj_gene_AG_BD <- gsub(" & $","",contig_paired_only$vdj_gene_AG_BD)
      #
      # #updating names to be consistant....
      contig_paired_only$vj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vj_gene_cdr3_BD,sep = " & ")
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub("^& ","",contig_paired_only$vj_gene_cdr3_AG_BD)
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub(" & $","",contig_paired_only$vj_gene_cdr3_AG_BD)

      contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vdj_gene_cdr3_BD,sep = " & ")
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub("^ & ","",contig_paired_only$vdj_gene_cdr3_AG_BD)
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub(" & $","",contig_paired_only$vdj_gene_cdr3_AG_BD)
      # contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vdj_gene_cdr3_BD,sep = " & ")
      names(contig_paired_only)[names(contig_paired_only) %in% "barcode"] <- "Cell_Index"
      contig_paired_only$Sample_Name <- input$sample_name_Array

      contig_paired_only <- contig_paired_only %>%
        select(all_of(c("Cell_Index","Sample_Name")), everything())

      contig_paired_only

    })


    output$test.files_array_contig_Filtered <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      calls <- Filtering_TCR_Array()
      validate(
        need(nrow(calls)>0,
             "Upload Contig TCR file")
      )
      calls
    })

    # clusTCR array -----
    ClusTCR2_array <- reactive({
      df <- Filtering_TCR_Array()
      df_nocouts2_AG <- df[,names(df) %in% c("Sample_Name","cdr3_AG","v_gene_AG")]
      df_nocouts2_BD <- df[,names(df) %in% c("Sample_Name","cdr3_BD","v_gene_BD")]
      names(df_nocouts2_AG) <- c("Sample_Name","junction_aa","v_call")
      names(df_nocouts2_BD) <- c("Sample_Name","junction_aa","v_call")



      df_nocouts3 <-as.data.frame(rbind(df_nocouts2_AG,df_nocouts2_BD))
      if (nrow(df_nocouts3[-c(grep("[*]",df_nocouts3$junction_aa)),]>0)) {
        df_nocouts3 <- df_nocouts3[-c(grep("[*]",df_nocouts3$junction_aa)),]
      }
      df_nocouts3 <- subset(df_nocouts3,df_nocouts3$v_call!="None")
      df_nocouts3

    })

    output$test.files_ClusTCR2_array <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      calls <- ClusTCR2_array()
      calls
    })
    output$download_ClusTCR2_labs_array <- downloadHandler(
      filename = function(){
        paste("ClusTCR2_output_Array",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(ClusTCR2_array())
        write.csv(df,file, row.names = F)
      } )
    # TCRex array  -----

    ## seurat count table ----
    tb_Array_matrix <- function () {
      mmMat <- as.data.frame(input.data.calls.array())

      validate(
        need(nrow(mmMat)>0,
             "Upload files")
      )
      mmMat
    }

    output$tb_array_matrix2 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      tb_Array_matrix()[1:6,1:6]
    })

    output$downloadtb_array_matrix2 <- downloadHandler(
      filename = function(){
        x = today()
        paste(input$name.array,"_count-matrix_array_",x,".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(t(tb_Array_matrix()))
        # write.table(,file, row.names = T)
        write.csv(df,file)
      })

    ## meta.data for seurat ----
    tb_Array_meta.data_TCR <- function () {
      contigs <- Filtering_TCR_Array()
      validate(
        need(nrow(contigs)>0,
             "Upload file")
      )
      contigs$Cell_Index <- rownames(contigs)
      contigs$Cell_Index <- gsub("[:]","-",contigs$Cell_Index)
      contigs$Cell_Index <- gsub("[.]","-",contigs$Cell_Index)
      contigs$Sample_Name <- input$sample_name_Array

      contig_paired_only <- contigs %>%
        select(all_of(c("Cell_Index","Sample_Name")), everything())

      contig_paired_only
    }

    output$tb_Array_meta1 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      tb_Array_meta.data_TCR()
    })
    output$downloadtb_array_metadata2 <- downloadHandler(
      filename = function(){
        paste(input$name.array,"_metadata_array_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){

        df <-  tb_Array_meta.data_TCR()
        write.csv(df, file, row.names = F)

      } )

    # TCRex Merge  ------
    input.data_TCRexMerge <- reactive({
      inFile2_TCRexMerge <- input$file2_TCRexMerge
      validate(
        need(nrow(inFile2_TCRexMerge)>0,
             "Upload mutliple TCRex files to merge")
      )
      num <- dim(inFile2_TCRexMerge)[1]
      samples_list <- vector("list", length = num)
      samples_list
      for (i in 1:num) {
        sc <- read.table(input$file2_TCRexMerge[[i, 'datapath']],sep="\t",header = T)
        samples_list[[i]] <- sc

      }

      samples_list
    })

    merged_TCRexFiltered <- reactive({
      myfiles <- input.data_TCRexMerge()
      df <- rbind(myfiles[[1]])
      for (i in 2:length(myfiles)) {
        df <- rbind(df,myfiles[[i]])

      }

      if (nrow(df[-c(grep("[*]",df$CDR3_beta)),]>0)) {
        df <- df[-c(grep("[*]",df$CDR3_beta)),]
      }

      df$CDR3_beta <- gsub("^$","None", df$CDR3_beta)

      if (nrow(df[-c(grep("None",df$CDR3_beta)),]>0)) {
        df <- df[-c(grep("None",df$CDR3_beta)),]
      }

      df$TRBJ_gene <- gsub("^$","None", df$TRBJ_gene)

      if (nrow(df[-c(grep("None",df$TRBJ_gene)),]>0)) {
        df <- subset(df,df$TRBJ_gene!= "None")
      }

      df[!duplicated(df$CDR3_beta),]
    })

    output$DEx_TCRexFiltered <-  DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      df <- merged_TCRexFiltered()

      df
    })


    output$downloaddf_TCRexFiltered <- downloadHandler(
      filename = function(){
        x = today()
        paste("TCRex_merged_",x,".tsv", sep = "")
      },
      content = function(file){
        df <- merged_TCRexFiltered()
        write.table(df,file, row.names = F,sep = "\t")
      })

    # observeEvent(input$downloaddf_TCRexFiltered,{
    #   write.table(merged_TCRexFiltered(),paste("TCRex_merged_",today(),".tsv", sep = ""), row.names = F,sep = "\t")
    # }
    #
    #             )



    # ClusTCR2 ------
    input.data_ClusTCR2  <- reactive({
      inFile2_ClusTCR2 <- input$file2_ClusTCR2
      if (is.null(inFile2_ClusTCR2)) return(NULL)
      else {
        dataframe <- read.csv(
          inFile2_ClusTCR2$datapath,header=TRUE)
      }
    })

    input.data_ClusTCR2_multiple <- reactive({
      inFile.seq <- input$file1_ClusTCR2_multiple
      num <- dim(inFile.seq)[1]
      samples_list <- vector("list", length = num)
      samples_list
      for (i in 1:num) {
        sc <- read.csv(input$file1_ClusTCR2_multiple[[i, 'datapath']])
        samples_list[[i]] <- sc

      }
      samples_list
    })

    merged_Clust_filtered <- reactive({
      inFile2_ClusTCR2 <- input$file1_ClusTCR2_multiple

      validate(
        need(nrow(inFile2_ClusTCR2)>0,
             "Upload mutliple CLusTCR2 files to merge")
      )

      myfiles <- input.data_ClusTCR2_multiple()
      df <- rbind(myfiles[[1]])
      for (i in 2:length(myfiles)) {
        df <- rbind(df,myfiles[[i]])

      }

      if (nrow(df[-c(grep("[*]",df$junction_aa)),]>0)) {
        df <- df[-c(grep("[*]",df$junction_aa)),]
      }

      df$CDR3_beta <- gsub("^$","None", df$junction_aa)

      if (nrow(df[-c(grep("None",df$CDR3_beta)),]>0)) {
        df <- df[-c(grep("None",df$CDR3_beta)),]
      }


      df[!duplicated(df$junction_aa),]
    })

    output$DEx_multiple_ClusTCR2 <-  DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      df <- merged_Clust_filtered()

      df
    })

    output$downloaddf_multiple_ClusTCR2 <- downloadHandler(
      filename = function(){
        x = today()
        paste(input$Clust_lab_tab_output,"_Multi_ClusTCR_",x,".csv", sep = "")
      },
      content = function(file){
        df <- merged_Clust_filtered()
        write.csv(df,file, row.names = F)
      })

    ## clustering images ------

    observe({
      updateSelectInput(
        session,
        "clusTCR2_names",
        choices=names(input.data_ClusTCR2()),
        selected = "junction_aa") }) # junction sequence
    observe({
      updateSelectInput(
        session,
        "clusTCR2_Vgene",
        choices=names(input.data_ClusTCR2()),
        selected = "v_call")}) # V_gene

    vals_ClusTCR2_old <- reactiveValues(output_dt_old=NULL)
    clust_dt2 <- reactive({
      df1 <- input.data_ClusTCR2()
      validate(
        need(nrow(df1)>0,
             "Upload ClusTCR file")
      )
      df1
    })
    output$clust_dt2_table <- DT::renderDataTable({
      df1 <- input.data_ClusTCR2()
      validate(
        need(nrow(df1)>0,
             "Upload ClusTCR file")
      )
      df1
    })

    # ## run clustering ------

    vals_ClusTCR2 <- reactiveValues(output_dt2=NULL)

    observeEvent(input$run_ClusTCR2,{

      df1 <- input.data_ClusTCR2()
      validate(
        need(nrow(df1)>0,
             "Upload ClusTCR file")
      )

      req(input$clusTCR2_names,input$clusTCR2_Vgene)

      clust_dt_DATA_5 <- df1[,names(df1) %in% c(input$clusTCR2_names,input$clusTCR2_Vgene)]
      ptm <- proc.time()

      df_cluster <- ClusTCR(clust_dt_DATA_5, allele =input$allele_ClusTCR2, v_gene = input$clusTCR2_Vgene)

      if (dim(df_cluster)[1]>1) {
        cluster_lab <- mcl_cluster(df_cluster,expansion = 1,inflation = 1)
        end <- proc.time() - ptm
        cluster_lab[[3]] <- end
        vals_ClusTCR2$output_dt2 <- cluster_lab
        vals_ClusTCR2$output_dt2
      }
    })

    output$verbatum_ClusTCR2 <- renderPrint({
      df1 <- input.data_ClusTCR2()
      validate(
        need(nrow(df1)>0,
             "Upload ClusTCR file")
      )
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")

      # message("Clusters exist")
      if ( is.null(vals_ClusTCR2$output_dt2) ) {
        message("No clusters == 1 edit distance and therefore MCL not performed")
      }

      else {

      }

      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })



    ClusTCR2_lab_df <- reactive({
      df1 <- input.data_ClusTCR2()
      validate(
        need(nrow(df1)>0,
             "Upload ClusTCR file")
      )
      req(vals_ClusTCR2$output_dt2)
      output_dt <- vals_ClusTCR2$output_dt2
      output_dt[[1]]
    })


    output$ClusTCR2_lab <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      df1 <- input.data_ClusTCR2()
      validate(
        need(nrow(df1)>0,
             "Upload ClusTCR file")
      )
      ClusTCR2_lab_df()
    })

    output$download_ClusTCR_labels <- downloadHandler(
      filename = function(){
        x = today()
        paste(input$Clust_lab_tab_output,"_ClusTCR2_output_",x,".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(ClusTCR2_lab_df())
        write.csv(df,file, row.names = F)
      } )

    output$ClusTCR2_Time <- renderPrint({

      df1 <- input.data_ClusTCR2()
      validate(
        need(nrow(df1)>0,
             "Upload ClusTCR file")
      )
      req(vals_ClusTCR2$output_dt2)

      output_dt <- vals_ClusTCR2$output_dt2

      if (is.null(output_dt)){return(NULL)}
      output_dt[[3]]
    })


    # create plots
    Network_plot_clusTCR2 <- reactive({
      df1 <- input.data_ClusTCR2()
      validate(
        need(nrow(df1)>0,
             "Upload ClusTCR file")
      )
      req(vals_ClusTCR2$output_dt2,input$filter_connections,input$selected_Cluster,input$lab_clust_by,input$Clust_size_order, input$text_size1,input$text_size2)
      Network_df <- vals_ClusTCR2$output_dt2
      set.seed(123)
      # ?netplot
      netplot_ClusTCR2(Network_df,
                       filter_plot = input$filter_connections,
                       Clust_selected = input$selected_Cluster,
                       label=input$lab_clust_by,
                       Clust_column_name = input$Clust_size_order,
                       colour = input$colour_ClusTCR2,
                       selected_text_size = input$text_size1,
                       non_selected_text_size = input$text_size2,
                       alpha_selected = 1,alpha_non_selected = 1,
                       selected_text_col = "black",
                       non_selected_text_col = "black",
                       all.colour = input$colour_ClusTCR2_types,
                       selected_col = input$sel_colour_netplot,
                       non_selected_col = input$nonsel_colour_netplot

      )

    })
    output$NP_ClusTCR <- renderPlot({
      df1 <- input.data_ClusTCR2()
      validate(
        need(nrow(df1)>0,
             "Upload ClusTCR file")
      )
      Network_plot_clusTCR2()
    })

    motif_plot_clusTCR2 <- reactive({
      df1 <- input.data_ClusTCR2()
      validate(
        need(nrow(df1)>0,
             "Upload ClusTCR file")
      )
      req(vals_ClusTCR2$output_dt2)
      Network_df <- vals_ClusTCR2$output_dt2
      set.seed(123)
      motif_plot(Network_df,Clust_selected = input$selected_Cluster,Clust_column_name = input$Clust_size_order)

    })

    output$MP_ClusTCR <- renderPlot({
      df1 <- input.data_ClusTCR2()
      validate(
        need(nrow(df1)>0,
             "Upload ClusTCR file")
      )

      motif_plot_clusTCR2()
    })

    ## downloading plot -----
    output$downloadPlot_Network_plot2 <- downloadHandler(
      filename = function() {
        x <- today()
        paste("clusTCR2_Network_plot_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Network_plot2,height=input$height_Network_plot2, onefile = FALSE) # open the pdf device
        plot(Network_plot_clusTCR2())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Network_plot2 <- downloadHandler(
      filename = function() {
        x <- today()
        paste("clusTCR2_Network_plot_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Network_plot2,
            height = input$height_png_Network_plot2,
            res = input$resolution_PNG_Network_plot2)
        plot(Network_plot_clusTCR2())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    # downlaod Motif
    output$downloadPlot_Motif_plot2 <- downloadHandler(
      filename = function() {
        x <- today()
        paste("clusTCR2_Motif_plot"," Cluster ",input$selected_Cluster," ",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Motif_plot2,height=input$height_Motif_plot2, onefile = FALSE) # open the pdf device
        plot(motif_plot_clusTCR2())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Motif_plot2 <- downloadHandler(
      filename = function() {
        x <- today()
        paste("clusTCR2_Motif_plot"," Cluster ",input$selected_Cluster," ", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Motif_plot2, height = input$height_png_Motif_plot2, res = input$resolution_PNG_Motif_plot2)
        plot(motif_plot_clusTCR2())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )


    # Seurat -----
    ## uploading the raw files ----
    input.data_sc <- reactive({
      inFile_sc <- input$file_SC
      if (is.null(inFile_sc)) return(NULL)
      else {
        if (input$df_seruatobj_type =="10x_Genomics (raw)") {
          dataframe = read.csv(inFile_sc$datapath)
        }

        else if (input$df_seruatobj_type =="10x_Genomics (.h5)") {
          dataframe <- suppressMessages(Read10X_h5(inFile_sc$datapath, use.names = TRUE, unique.features = TRUE))
          # if (input$stored_in_expression=="yes") {
          #   rownames(df.test$`Gene Expression`) <- gsub("GRCh38___","",rownames(df.test$`Gene Expression`))
          #   sc <- CreateSeuratObject(counts = df.test$`Gene Expression`, project = input$project_name)
          #
          # }
          #
          # else {
          #   sc <- CreateSeuratObject(counts = df.test, project = input$project_name)
          # }

        }

        else if (input$df_seruatobj_type =="BD Rhapsody (Mouse)") {
          dataframe = read.csv(inFile_sc$datapath,row.names = 1)
        }

        else if (input$df_seruatobj_type =="BD Rhapsody (Human Immune panel)") {
          dataframe = read.csv(inFile_sc$datapath,row.names = 1)
        }

        else {
          dataframe = read.csv(inFile_sc$datapath,row.names = 1)
        }
      }
    })

    output$DEx_header_name_check.dt <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      df.test <- input.data_sc()
      validate(
        need(nrow(df.test)>0,
             error_message_val_sc)
      )
      if (input$df_seruatobj_type =="10x_Genomics (raw)") {
        names(df.test) <- gsub("[.]1","-1",names(df.test) )
        rownames(df.test) <- make.unique(df.test$Gene_Name)
        df.test2 <- df.test[,!names(df.test) %in% c("Gene_Name")]
      }

      else if (input$df_seruatobj_type =="10x_Genomics (.h5)") {
        df.test2 <- as.data.frame("file possibly too big and will not be rendered")
      }
      else if (input$df_seruatobj_type =="BD Rhapsody (Mouse)") {
        names(df.test) <- gsub("X","",names(df.test))
        df.test2 <- df.test[!rownames(df.test) %in% c("Cell_Index"),]
      }

      else if (input$df_seruatobj_type =="BD Rhapsody (Human Immune panel)") {
        names(df.test) <- gsub("X","",names(df.test))
        df.test2 <- df.test[!rownames(df.test) %in% c("Cell_Index"),]
      }

      else {
        names(df.test) <- gsub("[.]","-",names(df.test))
        rownames(df.test) <- gsub("[.]","-",rownames(df.test))
        df.test2 <- df.test[!rownames(df.test) %in% c("Cell_Index"),]

      }
      head(df.test2)[1:6]
    })

    ## reading in 10x and BD data ----
    df_seruatobj <- reactive({
      df.test <- input.data_sc()
      validate(
        need(length(df.test)>0,
             error_message_val_sc)
      )
      if (input$df_seruatobj_type =="10x_Genomics (raw)") {

        names(df.test) <- gsub("[.]1","-1",names(df.test) )
        rownames(df.test) <- make.unique(df.test$Gene_Name)
        df.test2 <- df.test[,!names(df.test) %in% c("Gene_Name")]
        sc <- CreateSeuratObject(counts = df.test2, assay = "RNA", project = input$project_name)
        sc <- PercentageFeatureSet(sc, pattern = "^MT-",col.name = "mtDNA")
        sc <- PercentageFeatureSet(sc, pattern = "^RP[SL]",col.name = "rRNA")
        sc
      }

      else if (input$df_seruatobj_type=="10x_Genomics (.h5)") {
        # rownames(df.test$`Gene Expression`) <- gsub("GRCh38___","",rownames(df.test$`Gene Expression`))
        # sc <- CreateSeuratObject(counts = df.test[[1]], project = input$project_name)
        sc <- CreateSeuratObject(counts = df.test, project = input$project_name)
        sc <- PercentageFeatureSet(sc, pattern = "^MT-",col.name = "mtDNA")
        sc <- PercentageFeatureSet(sc, pattern = "^RP[SL]",col.name = "rRNA")
        sc
      }
      else if (input$df_seruatobj_type =="BD Rhapsody (Mouse)") {
        names(df.test) <- as.character(gsub("X","",names(df.test)))

        # rownames(df.test) <- make.unique(df.test$Gene_Name)
        # df.test2 <- df.test[,!names(df.test) %in% c("Gene_Name")]

        sc <- CreateSeuratObject(counts = df.test, assay = "RNA",project = input$project_name)
        sc <- PercentageFeatureSet(sc, pattern = "^Mt",col.name = "mtDNA")
        sc <- PercentageFeatureSet(sc, pattern = "Rp[sl]",col.name = "rRNA")
        sc
      }

      else if (input$df_seruatobj_type =="BD Rhapsody (Human Immune panel)") {
        names(df.test) <- as.character(gsub("X","",names(df.test)))

        # rownames(df.test) <- make.unique(df.test$Gene_Name)
        # df.test2 <- df.test[,!names(df.test) %in% c("Gene_Name")]

        sc <- CreateSeuratObject(counts = df.test, assay = "RNA",project = input$project_name)
        sc <- PercentageFeatureSet(sc, pattern = "^MT-",col.name = "mtDNA")
        sc <- PercentageFeatureSet(sc, pattern = "^RP[SL]",col.name = "rRNA")
        sc
      }

      else {
        names(df.test) <- gsub("[.]","-",names(df.test))
        rownames(df.test) <- gsub("[.]","-",rownames(df.test))
        sc <- CreateSeuratObject(counts = df.test, assay = "RNA", project = input$project_name)
        sc <- PercentageFeatureSet(sc, pattern = "^MT-",col.name = "mtDNA")
        sc <- PercentageFeatureSet(sc, pattern = "^RP[SL]",col.name = "rRNA")

      }
    })
    before_plot <- reactive({
      sc <- df_seruatobj()
      validate(
        need(nrow(sc)>0,
             "Upload files")
      )
      VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "mtDNA","rRNA"), ncol = 2)
    })

    output$before_plot_sc <- renderPlot({
      withProgress(message = 'Figure is being generated...',
                   detail = '', value = 0, {
                     test_fun()
                   })
      before_plot()
    })


    output$downloadPlot_before_plot_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name2,"_before_plot_sc_",x,".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_before_plot_sc,height=input$height_before_plot_sc, onefile = FALSE) # open the pdf device
        plot(before_plot())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_before_plot_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name2,"_before_plot_sc_",x,".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_before_plot_sc, height = input$height_png_before_plot_sc, res = input$resolution_PNG_before_plot_sc)
        plot(before_plot())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    ##### after filtering plot -------


    vals2 <- reactiveValues(after_violin_plot=NULL)
    observeEvent(input$run,{
      sc <- df_seruatobj()
      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )
      vals2$after_violin_plot <- subset(sc, subset = nFeature_RNA >= input$features.min & nFeature_RNA <= input$features.max & mtDNA <= input$percent.mt & rRNA >= input$percent.rb)
    })

    output$after_plot_sc <- renderPlot({
      sc <- vals2$after_violin_plot
      validate(
        need(nrow(sc)>0,
             "Run Filtering")
      )
      VlnPlot(vals2$after_violin_plot, features = c("nFeature_RNA", "nCount_RNA", "mtDNA","rRNA"), ncol = 2)
    })

    output$downloadPlot_after_plot_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name2,"_after_plot_sc_",x,".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_after_plot_sc,height=input$height_after_plot_sc, onefile = FALSE) # open the pdf device
        plot_after <- VlnPlot(vals2$after_violin_plot, features = c("nFeature_RNA", "nCount_RNA", "mtDNA","rRNA"), ncol = 2)
        plot(plot_after)
        dev.off()}, contentType = "application/pdf")

    output$downloadPlotPNG_after_plot_sc <- downloadHandler(
      filename = function() {
        x <- today()

        paste(input$project_name2,"_after_plot_sc_",x,".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_after_plot_sc, height = input$height_png_after_plot_sc, res = input$resolution_PNG_after_plot_sc)
        plot_after <- VlnPlot(vals2$after_violin_plot, features = c("nFeature_RNA", "nCount_RNA", "mtDNA","rRNA"), ncol = 2)
        plot(plot_after)
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    ### normalisationa and feature plot ------
    feature_serartobj <- reactive({
      sc <- vals2$after_violin_plot
      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )
      sc <- NormalizeData(sc)
      sc <- FindVariableFeatures(sc, selection.method = input$method_Seurat)
    })

    plot_10_features <- reactive({
      sc <- feature_serartobj()
      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )
      # Identify the 10 most highly variable genes
      top10 <- head(VariableFeatures(sc), 10)
      # plot variable features with and without labels
      plot1 <- VariableFeaturePlot(sc)
      plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

      plot2
    })
    output$plot_10_features_sc<- renderPlot({
      plot_10_features()
    })

    ###### download variable features (QC) -----
    output$downloadPlot_plot_10_features_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name,"_10_features_sc_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_plot_10_features_sc,height=input$height_plot_10_features_sc, onefile = FALSE) # open the pdf device
        plot(plot_10_features())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_plot_10_features_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name,"_10_features_sc_",x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_plot_10_features_sc, height = input$height_png_plot_10_features_sc, res = input$resolution_PNG_plot_10_features_sc)
        plot(plot_10_features())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )


    ## PCA and choosing # of dimensions to reduce ----
    create_PCA <- reactive({
      sc <- feature_serartobj()
      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )
      all.genes <- rownames(sc)
      sc <- ScaleData(sc, features = all.genes)
      sc <- RunPCA(sc, features = VariableFeatures(object = sc))
      sc
    })
    create_PCA_heatmap <- reactive({
      DimHeatmap(create_PCA(), dims = input$dimension_heatmap.min:input$dimension_heatmap.max, cells = input$numberofcells, balanced = TRUE)
    })
    output$create_PCA_heatmap_sc<- renderPlot({
      create_PCA_heatmap()
    })


    output$downloadPlot_heatmap_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name2,"_heatmap_sc_",x,".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_heatmap_sc,height=input$height_heatmap_sc, onefile = FALSE) # open the pdf device
        plot(create_PCA_heatmap())
        dev.off()}, contentType = "application/pdf")

    output$downloadPlotPNG_after_plot_sc <- downloadHandler(
      filename = function() {
        x <- today()

        paste(input$project_name2,"_heatmap_sc_",x,".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_heatmap_sc, height = input$height_png_heatmap_sc, res = input$resolution_PNG_heatmap_sc)
        plot(create_PCA_heatmap())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    ##### create elbow plot ----
    create_elbowPlot <- reactive({
      ElbowPlot(create_PCA())
    })
    output$create_elbowPlot_sc<- renderPlot({
      create_elbowPlot()
    })

    #######  download Elbow plot -----
    output$downloadPlot_create_elbowPlot_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name,"_elbowPlot_sc_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_create_elbowPlot_sc,height=input$height_create_elbowPlot_sc, onefile = FALSE) # open the pdf device
        plot(create_elbowPlot())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_create_elbowPlot_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name,"_elbowPlot_sc_",x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_create_elbowPlot_sc, height = input$height_png_create_elbowPlot_sc, res = input$resolution_PNG_create_elbowPlot_sc)
        plot(create_elbowPlot())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )


    ######  UMAP plot -----
    vals_clust <- reactiveValues(sc_clustering=NULL)

    observeEvent(input$run_reduction,{
      sc <- create_PCA()

      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )

      sc<- FindNeighbors(sc, dims = 1:input$dimension_sc)
      sc <- FindClusters(sc, resolution = input$resolution)
      sc <- RunUMAP(sc, dims = 1:input$dimension_sc)
      vals_clust$sc_clustering <- sc
      vals_clust$sc_clustering
    })

    create_UMAP <- reactive({
      sc <- vals_clust$sc_clustering

      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )

      DimPlot(sc, reduction = "umap")
    })
    output$create_UMAP_sc<- renderPlot({
      create_UMAP()
    })

    ## Differential expression with two conditions -----

    input.data_sc_meta <- reactive({
      inFile_sc_meta <- input$file_SC_meta
      if (is.null(inFile_sc_meta)) return(NULL)
      else {
        dataframe = read.csv(inFile_sc_meta$datapath)
      }
    })

    output$DEx_view.meta.dt <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      sc <- input.data_sc_meta()

      validate(
        need(nrow(sc)>0,
             "Upload metadata")
      )
      sc
    })
    vals_meta.sc <- reactiveValues(metadata_SCobj=NULL)
    observeEvent(input$run_metadata,{
      sc <- vals_clust$sc_clustering
      validate(
        need(nrow(sc)>0,
             "Run clustering or add in metadata")
      )

      meta.data.import <- input.data_sc_meta()
      validate(
        need(nrow(meta.data.import)>0,
             "Impute Metadata")
      )
      sc@meta.data$Cell_Index <- rownames(sc@meta.data)

      sc@meta.data$order <- 1:length(sc@meta.data$Cell_Index)
      scMeta.data <- sc@meta.data
      meta.data2 <- merge(scMeta.data,meta.data.import,by="Cell_Index",all.x=T)
      sc@meta.data <- meta.data2
      rownames(sc@meta.data) <- sc@meta.data$Cell_Index
      sc@meta.data <- sc@meta.data[order((sc@meta.data$order)),]
      vals_meta.sc$metadata_SCobj <- sc
      vals_meta.sc$metadata_SCobj
    })

    output$DEx_table_meta.data <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- vals_meta.sc$metadata_SCobj
      validate(
        need(nrow(calls)>0,
             "Impute metadata")
      )
      calls@meta.data$Cell_Index <- gsub("[.]","-",calls@meta.data$Cell_Index)
      calls@meta.data
    })

    # save Seurat object -----
    output$downloaddf_SeruatObj <- downloadHandler(
      filename = function(){
        x = today()
        # paste(input$project_name,"_SC.obj_",x,".h5Seurat", sep = "")
        paste(input$project_name,"_SC.obj_",x,".rds", sep = "")
      },
      content = function(file){
        SaveSeuratRds(vals_meta.sc$metadata_SCobj, file)
        # SaveH5Seurat(vals_meta.sc$metadata_SCobj,file)
      })

    # merging multiple Seurat Obj -----
    getData <- reactive({
      inFile.seq <- input$file1_rds.file
      num <- dim(inFile.seq)[1]
      seurat_object_list <- vector("list", length = num)
      seurat_object_list
      for (i in 1:num) {
        sc <- LoadSeuratRds(input$file1_rds.file[[i, 'datapath']])
        sc <-  DietSeurat(
          sc,
          scale.data = NULL
        )
        sc@meta.data$Cell_Index_old <- sc@meta.data$Cell_Index
        sc@meta.data$Noval.ID <- paste(sc@meta.data$orig.ident,sc@meta.data$Cell_Index,i,sep = "_")
        sc@meta.data$Noval.ID <- ifelse(grepl("NA[.]",sc@meta.data$Noval.ID),"",sc@meta.data$Noval.ID)
        seurat_object_list[[i]] <- sc
      }
      seurat_object_list
    })

    output$testing_mult <- renderPrint({
      sc <- input$file1_rds.file
      validate(
        need(nrow(sc)>0,
             "Upload files")
      )
      df <- getData()
      print(df)
    })
    merging_sc <- reactive({
      sc <- input$file1_rds.file
      validate(
        need(nrow(sc)>0,
             "Upload files")
      )
      samples_list <- getData()
      num <- length(samples_list)
      pbmc.normalized <- (samples_list[[1]])
      for (i in 2:num ) {
        pbmc.normalized <- merge(pbmc.normalized, y = samples_list[[i]],merge.data = TRUE)
      }

      pbmc.normalized
    })

    output$testing_mult2 <- renderPrint({

      merging_sc()

    })

    Vals_norm <- reactiveValues(Norm1=NULL)
    Vals_norm1 <- reactiveValues(Norm1=NULL)
    Vals_norm2 <- reactiveValues(Norm1=NULL)
    Vals_norm3 <- reactiveValues(Norm1=NULL)
    Vals_norm4 <- reactiveValues(Norm1=NULL)
    Vals_normk <- reactiveValues(anno=NULL)
    # Vals_normk <- reactiveValues(anno=NULL)

    # find variable features
    observeEvent(input$run_var,{
      sc <- merging_sc()
      validate(
        need(nrow(sc)>0,
             "Run Variable")
      )
      all.genes <- rownames(sc)

      if(length(all.genes)<2000) {
        sc <- FindVariableFeatures(sc, selection.method = "vst")
      } else {
        sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
      }
      Vals_norm4$Norm1 <- sc
    })

    output$var_harmony_verbrose <- renderPrint({
      sc <- Vals_norm4$Norm1
      validate(
        need(nrow(sc)>0,
             "Run Variable")
      )
      sc
    })
    # scale data
    observeEvent(input$run_scale,{
      sc <- Vals_norm4$Norm1
      validate(
        need(nrow(sc)>0,
             "Run Scale")
      )
      req(Vals_norm4$Norm1)
      kmeans <- read.csv(system.file("Kmean","Kmeans.requires.annotation.csv",package = "STEGO.R"))


      if (length(sc@assays$RNA@meta.data$var.features)>0) {
        var.genes <- as.data.frame(sc@assays$RNA@meta.data$var.features)

      } else {
        var.genes <- as.data.frame(sc@assays$RNA@var.features)
      }

      names(var.genes) <- "V1"
      if (input$sample.type.source == 'hs') {
        kmeans2 <- as.data.frame(kmeans$Human)
        names(kmeans2) <- "V1"
      }
      else {

        kmeans2 <- as.data.frame(kmeans$Mouse)
        names(kmeans2) <- "V1"
      }
      Vals_normk$anno <- unique(rbind(var.genes,kmeans2))
    })

    output$Tb_scaling_features_for_annotation <-  DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      sc <- Vals_norm4$Norm1
      validate(
        need(nrow(sc)>0,
             "Run Scale")
      )
      Vals_normk$anno
    })


    observeEvent(input$run_scale,{
      sc <- Vals_norm4$Norm1
      validate(
        need(nrow(sc)>0,
             "Run Scale")
      )
      all.genes <- Vals_normk$anno
      sc <- ScaleData(sc, features = all.genes$V1)
      Vals_norm1$Norm1 <- sc

    })

    output$scale_harmony_verbrose <- renderPrint({
      sc <- Vals_norm1$Norm1

      validate(
        need(nrow(sc)>0,
             "Run Scale")
      )
      sc

    })

    observeEvent(input$run_PCA,{
      sc <- Vals_norm1$Norm1
      validate(
        need(nrow(sc)>0,
             "Run PCA")
      )
      sc <- RunPCA(sc)
      Vals_norm2$Norm1 <- sc
    })

    output$PCA_harmony_verbrose <- renderPrint({
      sc <- Vals_norm2$Norm1
      validate(
        need(nrow(sc)>0,
             "Run PCA")
      )
      head(sc@meta.data)
    })

    observeEvent(input$run_harmony,{
      sc <- Vals_norm2$Norm1
      validate(
        need(nrow(sc)>0,
             "Run Harmony")
      )
      sc <- RunHarmony(sc, "orig.ident", plot_convergence = TRUE)
      Vals_norm3$Norm1 <- sc
    })


    output$harmony_verbrose <- renderPrint({
      sc <- Vals_norm3$Norm1
      validate(
        need(nrow(sc)>0,
             "Run Harmony")
      )
      sc
    })

    observeEvent(input$run_reduction_harmony,{
      sc <- Vals_norm3$Norm1
      sc <- sc %>%
        RunUMAP(reduction = "harmony", dims = 1:input$dimension_Merged) %>%
        FindNeighbors(reduction = "harmony", dims = 1:input$dimension_Merged) %>%
        FindClusters(resolution = input$res_merged) %>%
        identity()
      sc@meta.data$Cell_Index_2 <- sc@meta.data$Cell_Index
      sc@meta.data$Cell_Index <- rownames(sc@meta.data)
      Vals_norm$Norm1 <- sc


    })
    output$testing_mult3 <- renderPrint({
      df <- Vals_norm$Norm1
      validate(
        need(nrow(df)>0,
             "Run reduction")
      )
      df
    })
    output$create_UMAP_merged <- renderPlot({
      sc <- Vals_norm$Norm1
      validate(
        need(nrow(sc)>0,
             "Run reduction")
      )
      DimPlot(sc, reduction = "umap", group.by = "orig.ident", pt.size = 1)
    })

    observeEvent(input$download_reduction_harmony,{
      req(Vals_norm$Norm1)
      x = today()
      as.h5Seurat(Vals_norm$Norm1,paste(input$project_name2,"_merged_",x,".h5Seurat", sep = ""))
    })

    # download Harmony merged ----
    output$downloadPlot_sc_merged <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name2,"_SC_Merged_UMAP",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_sc_merged,height=input$height_sc_merged, onefile = FALSE) # open the pdf device
        plot(DimPlot(Vals_norm$Norm1, reduction = "umap", group.by = "orig.ident", pt.size = 1))
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_sc_merged <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name2,"_SC_Merged_UMAP", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_sc_merged, height = input$height_png_sc_merged, res = input$resolution_PNG_sc_merged)
        plot(DimPlot(Vals_norm$Norm1, reduction = "umap", group.by = "orig.ident", pt.size = 1))
        dev.off()},   contentType = "application/png" # MIME type of the image
    )



    # output$downloaddf_SeruatObj_merged <- downloadHandler(
    #   filename = function(){
    #     x = today()
    #     paste(input$project_name2,"_merged_",x,".h5Seurat", sep = "")
    #   },
    #   content = function(file){
    #     as.h5Seurat(Vals_norm$Norm1,file)
    #   } )


    output$downloaddf_SeruatObj_merged <- downloadHandler(
      filename = function(){
        x = today()
        paste(input$project_name2,"_merged_",x,".rds", sep = "")
      },
      content = function(file){
        SaveSeuratRds(Vals_norm$Norm1, file)
      })

    # remove unwanted cells from scOBJ ------
    getData_SampRemove <- reactive({
      inFile_sc_SampRemove <- input$file1_rds.fileSampsRemove
      if (is.null(inFile_sc_SampRemove))
        return(NULL)

      LoadSeuratRds(inFile_sc_SampRemove$datapath)

    })

    output$Preliminary_samp_to_remove <- renderPrint({
      inFile_sc_SampRemove <- input$file1_rds.fileSampsRemove
      if (is.null(inFile_sc_SampRemove))
        return(NULL)

      LoadSeuratRds(inFile_sc_SampRemove$datapath)
    })


    observe({
      sc <- getData_SampRemove()
      validate(
        need(nrow(sc)>0,
             "Upload .h5Seurat object")
      )

      df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "Samp_col_SampToRemove",
        choices=names(df3.meta),
        selected = "Sample_Name")
    })



    select_group_metadata_SampToRemove <- reactive({
      sc <- getData_SampRemove()

      validate(
        need(nrow(sc)>0,
             "upload file")
      )
      df <- sc@meta.data
      df2 <- as.data.frame(unique(df[names(df) %in% input$Samp_col_SampToRemove]))
      df2 <- as.data.frame(df2)
      df2

    })

    observe({
      df2 <- select_group_metadata_SampToRemove()
      validate(
        need(nrow(df2)>0,
             error_message_val1)
      )
      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- as.data.frame(df2[order(df2$V1),])
      names(df2) <- "V1"
      df2
      # df2 <- subset(df2,df2$V1 != "NA")

      df3 <- subset(df2,df2$V1 != "NA")

      updateSelectInput(
        session,
        "ID_Column_factor_SampToRemove",
        choices=df2$V1,
        selected = df3$V1
      )
    })

    Filtered_samp_to_remove_process <- reactive({
      sc <- input$file1_rds.fileSampsRemove
      validate(
        need(nrow(sc)>0,
             "Upload files")
      )
      sc <- getData_SampRemove()

      sc@meta.data$selected <- sc@meta.data[,names(sc@meta.data) %in% input$Samp_col_SampToRemove]
      sc@meta.data$keep <- ifelse(sc@meta.data$selected %in% c(input$ID_Column_factor_SampToRemove),"keep","NS")
      sc <- subset(x = sc, subset = keep == "keep")
      sc@meta.data <- sc@meta.data[,!names(sc@meta.data) %in% c("selected","keep")]
      sc

    })

    output$Filtered_samp_to_remove <- renderPrint({
      sc <- input$file1_rds.fileSampsRemove
      validate(
        need(nrow(sc)>0,
             "Upload files")
      )
      sc2 <- Filtered_samp_to_remove_process()
      print(sc2)

    })

    # output$downloaddf_SeruatObj_annotated_SampToKeep <- downloadHandler(
    #   filename = function(){
    #     x = today()
    #     paste(input$project_name4,"_keep_",x,".h5Seurat", sep = "")
    #   },
    #   content = function(file){
    #       sc <- Filtered_samp_to_remove_process()
    #       as.h5Seurat(sc,file)
    #
    #   } )

    output$downloaddf_SeruatObj_annotated_SampToKeep <- downloadHandler(
      filename = function(){
        x = today()
        paste(input$project_name4,"_keep_",x,".rds", sep = "")
      },
      content = function(file){
        sc <- Filtered_samp_to_remove_process()
        SaveSeuratRds(sc, file)
      })

    # Add cell annotations to merged Seurat object -----
    getData_2 <- reactive({
      inFile_sc_pro2 <- input$file1_rds.file2
      if (is.null(inFile_sc_pro2)) return(NULL)
      else {
        dataframe = LoadSeuratRds(inFile_sc_pro2$datapath)
      }

    })

    output$testing_mult_anno <- renderPrint({
      sc <- input$file1_rds.file2
      validate(
        need(nrow(sc)>0,
             "Upload files")
      )
      df <- getData_2()

      print(df)

    })

    output$Detect_version <- renderUI({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload files")
      )
      VersionControl <- sc@version
      if(length(grep("4.",VersionControl))>0) {
        selectInput("SeuratVersion2","Seurat Version", choices = c("Version 4","Version 5"), selected = "Version 4")
      } else if (length(grep("5.",VersionControl))>0) {
        selectInput("SeuratVersion2","Seurat Version", choices = c("Version 4","Version 5"), selected = "Version 5")
      } else {
        selectInput("SeuratVersion2","Seurat Version", choices = c("Version 4","Version 5"), selected = "")
      }

    })


    observe({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      df3.meta <- sc@meta.data

      updateSelectInput(
        session,
        "V_gene_Class_2",
        choices=names(df3.meta),
        selected = "vdj_gene_cdr3_AG_BD")
    })

    Vals_norm2 <- reactiveValues(Norm2=NULL)
    Annotation <- reactiveValues(LengthofAnno=NULL)
    ## add classification based on TCR-seq ----
    TCR_seq_classification <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload files")
      )

      if (input$sample.type.source.markers == "hs") {
        sc@meta.data$unconventional <- ifelse(sc@meta.data$vj_gene_AG == "TRAV1-2.TRAJ33","MAIT",
                                              ifelse(sc@meta.data$vj_gene_AG == "TRAV1-2.TRAJ12","MAIT",
                                                     ifelse(sc@meta.data$vj_gene_AG == "TRAV1-2.TRAJ23","MAIT",
                                                            ifelse(sc@meta.data$vj_gene_AG == "TRAV10.TRAJ18","iNKT",
                                                                   ifelse(sc@meta.data$v_gene_BD == "TRBV4-1","possible CD1c-restricted",
                                                                          ifelse(sc@meta.data$chain_AG == 'TRG' & sc@meta.data$chain_BD == 'TRB',"gb T cell",
                                                                                 ifelse(sc@meta.data$chain_AG == 'TRG' & sc@meta.data$chain_BD == 'TRD',"gd T cell",NA

                                                                                 )))))))
      }

      else {

      }


      sc

    })

    output$TCR_seq_classification_df <-  DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      sc <- TCR_seq_classification()
      sc@meta.data
    })



    # scGATE annotations HS 10x -------


    scGATE_anno_generic <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      req(input$threshold_scGate)

      len <- length(rownames(sc@assays$RNA$scale.data))

      if (len>3000){
        len = 3000
      }

      if (input$generic_scGATE==T) {
        scGate_models_DB <- custom_db_scGATE(system.file("scGATE","human/generic",package = "STEGO.R"))
        models.list <- scGate_models_DB[c(input$GenericID_scGATE)]
        obj <- scGate(sc, model = models.list,
                      pos.thr = input$threshold_scGate,
                      neg.thr = input$threshold_scGate,
                      nfeatures = len,
                      ncores = 8 )
        sc@meta.data$generic <- obj@meta.data$scGate_multi
        sc
      }
      else {
        sc
      }
      sc
    })
    scGATE_anno_CD4 <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      req(input$threshold_scGate)

      len <- length(rownames(sc@assays$RNA$scale.data))

      if (len>3000){
        len = 3000
      }


      if (input$CD4_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/CD4_TIL",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$CD4_TIL <- obj@meta.data$scGate_multi
      }
      else {
        sc
      }

      sc

    })
    scGATE_anno_CD8 <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>3000){
        len = 3000
      }
      if (input$CD8_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/CD8_TIL",package = "STEGO.R")))
        models.list <- scGate_models_DB
        models.list

        obj <- obj <- obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$CD8_TIL <- obj@meta.data$scGate_multi
      }
      else {
        sc
      }
      sc
    })
    scGATE_anno_GeneralMarkers <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>3000){
        len = 3000
      }
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])

      names(df) <- "Cell_Index"
      if (input$GeneralMarkers_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/ECSA",package = "STEGO.R")))
        models.list <- scGate_models_DB
        models.list
        obj <- obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$GeneralMarkers <- obj@meta.data$scGate_multi
      }
      else {
        sc
      }
      sc
    })
    scGATE_anno_Function <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }
      if (input$Function_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/Function",package = "STEGO.R")))
        models.list <- scGate_models_DB
        models.list
        obj <- obj <- obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$Function <- obj@meta.data$scGate_multi
      }
      else {
        sc
      }
      sc
    })
    scGATE_anno_Ex.Sen <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }
      if (input$Ex.Sen_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/Exhausted_Senescence",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$Exhausted_Senescence <- obj@meta.data$scGate_multi
      }
      else {
        sc
      }
      sc
    })
    scGATE_anno_COVID <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }
      if (input$COVID_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/COVID",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$COVID <- obj@meta.data$scGate_multi
      }
      else {
        sc
      }
      sc
    })
    scGATE_anno_Activation <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }
      if (input$Activation_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/Activation",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$Activation <- obj@meta.data$scGate_multi
      }
      else {
        sc
      }
      sc
    })
    scGATE_anno_IFNgTNFa <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }
      if (input$IFNgTNFa_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/IFNgTNFa",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$IFNgTNFa <- obj@meta.data$scGate_multi
      }
      else {
        sc
      }
      sc
    })
    scGATE_anno_GNLY.PFR1.GZMB <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }
      if (input$GNLY.PFR1.GZMB_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/GNLY.PFR1.GZMB",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$GNLY.PFR1.GZMB <- obj@meta.data$scGate_multi
      }
      else {
        sc
      }
      sc
    })
    scGATE_anno_Interlukin <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }
      if (input$Interlukin_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/IL",package = "STEGO.R")))
        models.list <- scGate_models_DB
        models.list
        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$Interleukins <- obj@meta.data$scGate_multi
      }
      else {
        sc
      }
      sc
    })
    scGATE_anno_NK_common <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }
      if (input$Interlukin_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/IL",package = "STEGO.R")))
        models.list <- scGate_models_DB
        models.list
        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$Interleukins <- obj@meta.data$scGate_multi
      }
      else {
        sc
      }
      sc
    }) # need to do...
    # BD rhapsody gates HS immune panel -----
    scGATE_anno_BD_rhapsody <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      sc<- JoinLayers(sc)
      if (input$Data_types=="BD_HS.Immune.Panel") {

        req(input$threshold_scGate)
        # obtained the list of markers from the following article based on Mouse gut: # https://www.frontiersin.org/articles/10.3389/fimmu.2020.563414/full

        # Main T cell markers ------

        if (input$BDrhapsody_scGATE_Tcells==T ) {
          models.list <- list()
          my_scGate_model <- gating_model(name = "CD8B",level = 1, signature = c("CD8B","CD8A"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD4neg",level = 3, signature = c("CD4"), positive = F)
          models.list$CD8 <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD4",level = 1, signature = c("CD4"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD8",level = 2, signature = c("CD8B","CD8A"), positive = F) # positive for
          models.list$CD4 <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD8ANeg",level = 1, signature = c("CD8A","CD8B"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD4neg",level = 2, signature = c("CD4"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD3E",level = 3, signature = c("CD3E"), positive =T)
          models.list$DN <- my_scGate_model
          obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate)
          sc@meta.data$T_cells <- obj@meta.data$scGate_multi

        }


        # GNLY.PRF1.GZMB -----
        if (input$BDrhapsody_scGATE_GNLY.PRF1.GZMB==T ) {
          models.list <- list()
          my_scGate_model <- gating_model(name = "GZMB",level = 1, signature = c("GZMB"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "PRF1",level = 2, signature = c("PRF1"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "GNLY",level = 3, signature = c("GNLY"), positive = T)

          models.list$GNLY.PRF1.GZMB <- my_scGate_model
          my_scGate_model <- gating_model(name = "GZMB",level = 1, signature = c("GZMB"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "PRF1",level = 2, signature = c("PRF1"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "GNLY",level = 3, signature = c("GNLY"), positive = F)
          models.list$PRF1.GZMB <- my_scGate_model

          my_scGate_model <- gating_model(name = "GZMB",level = 1, signature = c("GZMB"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "PRF1",level = 2, signature = c("PRF1"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "GNLY",level = 3, signature = c("GNLY"), positive = T)
          models.list$PRF1.GNLY <- my_scGate_model

          my_scGate_model <- gating_model(name = "GZMB",level = 1, signature = c("GZMB"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "PRF1",level = 2, signature = c("PRF1"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "GNLY",level = 3, signature = c("GNLY"), positive = T)

          obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate)
          sc@meta.data$GNLY.PRF1.GZMB <- obj@meta.data$scGate_multi

        }

        # TNF and IFNg----
        if (input$BDrhapsody_scGATE_TNF.IFNG==T ) {
          models.list <- list()
          models.list

          my_scGate_model <- gating_model(name = "IFNG",level = 1, signature = c("IFNG"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "TNF",level = 2, signature = c("TNF"), positive = T)
          models.list$IFNG.TNF <- my_scGate_model

          my_scGate_model <- gating_model(name = "IFNG",level = 1, signature = c("IFNG"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "TNF",level = 2, signature = c("TNF"), positive = F)
          models.list$IFNG <- my_scGate_model

          my_scGate_model <- gating_model(name = "IFNG",level = 1, signature = c("IFNG"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "TNF",level = 2, signature = c("TNF"), positive = T)
          models.list$TNF <- my_scGate_model

          obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate)

          sc@meta.data$TNF.IFNG <- obj@meta.data$scGate_multi
          # FeaturePlot(sc,features = "IFNG") | FeaturePlot(sc,features = "TNF")|DimPlot(sc, group.by = "TNF.IFNG") + scale_colour_manual(values = rainbow(3), na.value = "grey90")
        }
        # T cell types (CD4_based)-----
        if (input$BDrhapsody_scGATE_Effector_CD8==T ) {
          models.list <- list()
          my_scGate_model <- gating_model(name = "CCR4",level = 1, signature = c("CCR4"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "FOXP3",level = 2, signature = c("FOXP3"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD8",level = 3, signature = c("CD8A","CD8B"), positive = F)
          models.list$Th2_like <- my_scGate_model
          my_scGate_model <- gating_model(name = "CXCR3",level = 1, signature = c("CXCR3"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "TBX21",level = 2, signature = c("TBX21"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD8",level = 3, signature = c("CD8A","CD8B"), positive = F)
          models.list$Th1_like <- my_scGate_model
          my_scGate_model <- gating_model(name = "RORC",level = 1, signature = c("RORC"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "FOXP3",level = 2, signature = c("FOXP3"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD8",level = 3, signature = c("CD8A","CD8B"), positive = F)
          models.list$Th17_like <- my_scGate_model

          my_scGate_model <- gating_model(name = "CXCR5",level = 1, signature = c("CXCR5"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "FOXP3",level = 2, signature = c("FOXP3"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD8",level = 3, signature = c("CD8A","CD8B"), positive = F)
          models.list$Tfh_like <- my_scGate_model

          my_scGate_model <- gating_model(name = "FOXP3",level = 1, signature = c("FOXP3"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD8",level = 2, signature = c("CD8A","CD8B"), positive = F)
          models.list$Treg_like <- my_scGate_model

          my_scGate_model <- gating_model(name = "KLRK1",level = 1, signature = c("KLRK1"), positive = T) # https://pubmed.ncbi.nlm.nih.gov/36705564/
          my_scGate_model <- gating_model(my_scGate_model,name = "IL7R",level = 2, signature = c("IL7R"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "CXCR3",level = 3, signature = c("CXCR3"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "TBX21",level = 3, signature = c("TBX21"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD8A",level =4, signature = c("CD8A"), positive = T)
          models.list$KILR_CD8_effector <- my_scGate_model

          my_scGate_model <- gating_model(name = "GZMA",level = 1, signature = c("GZMA"), positive = T) # https://pubmed.ncbi.nlm.nih.gov/36705564/
          my_scGate_model <- gating_model(my_scGate_model,name = "GZMK",level = 2, signature = c("GZMK"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD8",level = 3, signature = c("CD8A","CD8B"), positive = T)
          models.list$Effector_CD8 <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD3E",level = 1, signature = c("CD3E"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD14",level = 2, signature = c("CD14"), positive = T)
          models.list$Mono.Tcell.Complex <- my_scGate_model

          obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate)

          sc@meta.data$T_cell_types <- obj@meta.data$scGate_multi
        }
        # DimPlot(sc, group.by = "T_cell_types") + scale_colour_manual(values = rainbow(10), na.value = "grey90")

        #### lung resident T cells ------
        if (input$BDrhapsody_scGATE_Lung.Residence.markers ==T ) {
          models.list <- list()
          my_scGate_model <- gating_model(name = "CD69",level = 1, signature = c("CD69"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "ITGAE",level = 2, signature = c("ITGAE"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "ITGA4",level = 3, signature = c("ITGA4"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRK1",level = 4, signature = c("KLRK1"), positive = T)
          models.list$CD69.CD103.CD49d.NKG2D <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD69",level = 1, signature = c("CD69"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "ITGAE",level = 2, signature = c("ITGAE"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "ITGA4",level = 3, signature = c("ITGA4"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRK1",level = 4, signature = c("KLRK1"), positive = T)
          models.list$CD69.CD103.NKG2D <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD69",level = 1, signature = c("CD69"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "ITGAE",level = 2, signature = c("ITGAE"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "ITGA4",level = 3, signature = c("ITGA4"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRK1",level = 4, signature = c("KLRK1"), positive = T)
          models.list$CD103.CD49d.NKG2D <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD69",level = 1, signature = c("CD69"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "ITGAE",level = 2, signature = c("ITGAE"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "ITGA4",level = 3, signature = c("ITGA4"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRK1",level = 4, signature = c("KLRK1"), positive = T)
          models.list$CD69.CD49d.NKG2D <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD69",level = 1, signature = c("CD69"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "ITGAE",level = 2, signature = c("ITGAE"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "ITGA4",level = 3, signature = c("ITGA4"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRK1",level = 4, signature = c("KLRK1"), positive = F)
          models.list$CD69.CD103 <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD69",level = 1, signature = c("CD69"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "ITGAE",level = 2, signature = c("ITGAE"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "ITGA4",level = 3, signature = c("ITGA4"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRK1",level = 4, signature = c("KLRK1"), positive = F)
          models.list$CD103.CD49d <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD69",level = 1, signature = c("CD69"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "ITGAE",level = 2, signature = c("ITGAE"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "ITGA4",level = 3, signature = c("ITGA4"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRK1",level = 4, signature = c("KLRK1"), positive = F)
          models.list$CD69.CD49d <- my_scGate_model


          # models.list <- list()
          # my_scGate_model <- gating_model(name = "ITGAE",level = 1, signature = c("ITGAE"), positive = T)
          # models.list$CD103 <- my_scGate_model
          obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate)

          sc@meta.data$Lung.Residence.markers <- obj@meta.data$scGate_multi
          sc@meta.data$Lung.Residence.markers <- ifelse( sc@meta.data$Lung.Residence.markers=="NA",NA, sc@meta.data$Lung.Residence.markers)
        }

        # NK-Like phenotypes ------
        if (input$BDrhapsody_scGATE_NK_receptors==T ) {
          models.list <- list()

          my_scGate_model <- gating_model(name = "CD161",level = 1, signature = c("KLRB1"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRC1",level = 2, signature = c("KLRC1"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2E",level = 3, signature = c("KLRC3"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2F",level = 4, signature = c("KLRC4"), positive = T) #NKG2F
          models.list$NKG2A.NKG2E.NKG2F.CD161_NK1.1 <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD161",level = 1, signature = c("KLRB1"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRC1",level = 2, signature = c("KLRC1"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2E",level = 3, signature = c("KLRC3"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2F",level = 4, signature = c("KLRC4"), positive = F) #NKG2F
          models.list$NKG2A.NKG2E.CD161_NK1.1 <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD161",level = 1, signature = c("KLRB1"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRC1",level = 2, signature = c("KLRC1"), positive = F) # and HLA-E
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2E",level = 3, signature = c("KLRC3"), positive = T) #NKG2E and HLA-E
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2F",level = 4, signature = c("KLRC4"), positive = F) #NKG2F and HLA-E
          models.list$NKG2E.CD161_NK1.1 <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD161",level = 1, signature = c("KLRB1"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRC1",level = 2, signature = c("KLRC1"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2E",level = 3, signature = c("KLRC3"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2F",level = 4, signature = c("KLRC4"), positive = F) #NKG2F
          models.list$CD161_NK1.1 <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD161",level = 1, signature = c("KLRB1"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRC1",level = 2, signature = c("KLRC1"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2E",level = 3, signature = c("KLRC3"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2F",level = 4, signature = c("KLRC4"), positive = T) #NKG2F
          models.list$NKG2A.NKG2E.NKG2F <- my_scGate_model


          my_scGate_model <- gating_model(name = "CD161",level = 1, signature = c("KLRB1"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRC1",level = 2, signature = c("KLRC1"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2E",level = 3, signature = c("KLRC3"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2F",level = 4, signature = c("KLRC4"), positive = T) #NKG2F
          models.list$NKG2E.NKG2F <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD161",level = 1, signature = c("KLRB1"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRC1",level = 2, signature = c("KLRC1"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2E",level = 3, signature = c("KLRC3"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2F",level = 4, signature = c("KLRC4"), positive = T) #NKG2F
          models.list$NKG2A.NKG2F <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD161",level = 1, signature = c("KLRB1"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRC1",level = 2, signature = c("KLRC1"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2E",level = 3, signature = c("KLRC3"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2F",level = 4, signature = c("KLRC4"), positive = F) #NKG2F
          models.list$NKG2A.NKG2E <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD161",level = 1, signature = c("KLRB1"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRC1",level = 2, signature = c("KLRC1"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2E",level = 3, signature = c("KLRC3"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2F",level = 4, signature = c("KLRC4"), positive = F) #NKG2F
          models.list$NKG2A <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD161",level = 1, signature = c("KLRB1"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRC1",level = 2, signature = c("KLRC1"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2E",level = 3, signature = c("KLRC3"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2F",level = 4, signature = c("KLRC4"), positive = F) #NKG2F
          models.list$NKG2E <- my_scGate_model

          my_scGate_model <- gating_model(name = "CD161",level = 1, signature = c("KLRB1"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "KLRC1",level = 2, signature = c("KLRC1"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2E",level = 3, signature = c("KLRC3"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "NKG2F",level = 4, signature = c("KLRC4"), positive = T) #NKG2F
          models.list$NKG2F <- my_scGate_model

          obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate)
          # my_scGate_model <- gating_model(my_scGate_model,name = "NKp80",level = 3, signature = c("KLRF1"), positive = T)
          # FeaturePlot(sc,features = "KLRF1")

          sc@meta.data$NK_receptors <- obj@meta.data$scGate_multi

          # DimPlot(sc, group.by = "NK_like_cells") + scale_colour_manual(values = rainbow(5), na.value = "grey90")
        }
        # Activated -----
        if (input$BDrhapsody_scGATE_Other==T ) {
          models.list <- list()
          my_scGate_model <- gating_model(name = "HLA",level = 1, signature = c("HLA-DRA"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "IL2RA",level = 2, signature = c("IL2RA"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD38",level = 3, signature = c("CD38"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD69",level = 4, signature = c("CD69"), positive = T)
          models.list$CD69 <- my_scGate_model
          my_scGate_model <- gating_model(name = "HLA",level = 1, signature = c("HLA-DRA"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "IL2RA",level = 2, signature = c("IL2RA"), positive = T)
          models.list$CD25_HLA.DRA <- my_scGate_model
          my_scGate_model <- gating_model(name = "HLA",level = 1, signature = c("HLA-DRA"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "IL2RA",level = 2, signature = c("IL2RA"), positive = F)
          models.list$HLA.DRA <- my_scGate_model
          my_scGate_model <- gating_model(name = "HLA",level = 1, signature = c("HLA-DRA"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "IL2RA",level = 2, signature = c("IL2RA"), positive = T)
          models.list$CD25 <- my_scGate_model
          my_scGate_model <- gating_model(name = "HLA",level = 1, signature = c("HLA-DRA"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "IL2RA",level = 2, signature = c("IL2RA"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "CD38",level = 3, signature = c("CD38"), positive = T)
          # my_scGate_model <- gating_model(my_scGate_model,name = "CD69",level = 4, signature = c("CD69"), positive = F)
          models.list$CD38 <- my_scGate_model
          models.list

          obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate)

          sc@meta.data$Activation <- obj@meta.data$scGate_multi
        }
        # DimPlot(sc, group.by = "Activation") + scale_colour_manual(values = rainbow(4), na.value = "grey90")
        # FeaturePlot(sc,features = "CD69") | FeaturePlot(sc,features = "CD38")

        # other -----
        if (input$BDrhapsody_scGATE_Activation==T ) {
          models.list <- list()
          my_scGate_model <- gating_model(name = "PDCD1",level = 1, signature = c("PDCD1","TIGIT"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "B3GAT1",level = 2, signature = c("B3GAT1","KLRG1"), positive = T)
          models.list$exhausted_and_senescent <- my_scGate_model

          my_scGate_model <- gating_model(name = "PDCD1",level = 1, signature = c("PDCD1","TIGIT"), positive = T)
          my_scGate_model <- gating_model(my_scGate_model,name = "B3GAT1",level = 2, signature = c("B3GAT1","KLRG1"), positive = F)
          models.list$exhausted <- my_scGate_model

          my_scGate_model <- gating_model(name = "PDCD1",level = 1, signature = c("PDCD1","TIGIT"), positive = F)
          my_scGate_model <- gating_model(my_scGate_model,name = "B3GAT1",level = 2, signature = c("B3GAT1","KLRG1"), positive = T)
          models.list$senescent <- my_scGate_model
          my_scGate_model <- gating_model(name = "MKI67",level = 1, signature = c("MKI67","TOP2A"), positive = T)
          # my_scGate_model <- gating_model(my_scGate_model,name = "TOP2A",level = 2, signature = c("TOP2A"), positive = T)
          models.list$CellCycling <- my_scGate_model

          obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate)
          sc@meta.data$Other <- obj@meta.data$scGate_multi
          # DimPlot(sc, group.by = "Other") + scale_colour_manual(values = rainbow(6), na.value = "grey90") | FeaturePlot(sc,features = "IL2RA")
        }
      }
      else {
      }

      sc
    })

    # BD rhapsody gates Mouse Full panel -----

    scGATE_anno_BD_MM.FP_T.cell <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )

      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }

      if (input$BDrhapsody_scGATE.MM.Tcell==T)  {
        models.list <- list()
        my_scGate_model <- gating_model(name = "Cd8b1",level = 1, signature = c("Cd8b1","Cd8a"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "CD4neg",level = 2, signature = c("Cd4"), positive = F)
        models.list$CD8 <- my_scGate_model

        my_scGate_model <- gating_model(name = "Cd8b1",level = 1, signature = c("Cd8b1","Cd8a"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "CD4neg",level = 2, signature = c("Cd4"), positive = T)
        models.list$CD4 <- my_scGate_model

        my_scGate_model <- gating_model(name = "Cd8b1",level = 1, signature = c("Cd8b1","Cd8a"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "CD4neg",level = 2, signature = c("Cd4"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "CD3E",level = 3, signature = c("Cd3e"), positive =T)
        models.list$DN <- my_scGate_model
        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )

        sc@meta.data$T_cells <- obj@meta.data$scGate_multi
        sc
      }
      else {
        sc
      }

    })
    output$scGATE_verbatum_BDrhapsody_MM.FP.Tcell <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$BDrhapsody_scGATE.MM.Tcell==T) {
        scGATE_anno_BD_MM.FP_T.cell()
      }
      else{
        print("MM T cell not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })

    scGATE_anno_BD_MM.FP_Memory <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )

      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }
      if (input$BDrhapsody_scGATE.MM.Memory==T)  {

        models.list <- list()
        my_scGate_model <- gating_model(name = "Cd44",level = 1, signature = c("Cd44"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "Sell",level = 2, signature = c("Sell"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Ly6c2",level = 3, signature = c("Ly6c2"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "Itga4",level = 4, signature = c("Itga4"), positive = F)
        models.list$Naive <- my_scGate_model

        my_scGate_model <- gating_model(name = "Cd44",level = 1, signature = c("Cd44"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Sell",level = 2, signature = c("Sell"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Ly6c2",level = 3, signature = c("Ly6c2"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Itga4",level = 4, signature = c("Itga4"), positive = F)
        models.list$VM <- my_scGate_model

        my_scGate_model <- gating_model(name = "Cd44",level = 1, signature = c("Cd44"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Sell",level = 2, signature = c("Sell"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Ly6c2",level = 3, signature = c("Ly6c2"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Itga4",level = 4, signature = c("Itga4"), positive = T)
        models.list$CM <- my_scGate_model

        my_scGate_model <- gating_model(name = "Cd44",level = 1, signature = c("Cd44"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Sell",level = 2, signature = c("Sell"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "Itga4",level = 3, signature = c("Itga4"), positive = T)
        models.list$Eff.or.RM <- my_scGate_model

        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )

        sc@meta.data$Memory <- obj@meta.data$scGate_multi
      }
      else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_BDrhapsody_MM.FP.Memory <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$BDrhapsody_scGATE.MM.Memory==T) {
        scGATE_anno_BD_MM.FP_Memory()
      }
      else{
        print("MM memory not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })

    scGATE_anno_BD_MM.FP_signatures <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )

      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }

      if (input$BDrhapsody_scGATE.MM.signatures==T)  {
        models.list <- list()
        my_scGate_model <- gating_model(name = "Gzma",level = 1, signature = c("Gzma"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Gzmb",level = 2, signature = c("Gzmb"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Prf1",level = 3, signature = c("Prf1"), positive = T)
        models.list$Eff <- my_scGate_model

        my_scGate_model <- gating_model(name = "Sell",level = 1, signature = c("Sell"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Ccr7",level = 2, signature = c("Ccr7"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Cxcr6",level = 3, signature = c("Cxcr6"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "Itga1",level = 4, signature = c("Itga1"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "Prdm1",level = 5, signature = c("Prdm1"), positive = F)
        models.list$Recirculation <- my_scGate_model

        my_scGate_model <- gating_model(name = "Sell",level = 1, signature = c("Sell"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "Ccr7",level = 2, signature = c("Ccr7"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "Cxcr6",level = 3, signature = c("Cxcr6"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Itga1",level = 4, signature = c("Itga1"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Prdm1",level = 5, signature = c("Prdm1"), positive = T)
        models.list$Tiss_res <- my_scGate_model

        my_scGate_model <- gating_model(name = "Il2ra",level = 2, signature = c("Il2ra"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Foxp3",level = 2, signature = c("Foxp3"), positive = T)
        models.list$Treg <- my_scGate_model

        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$Signature <- obj@meta.data$scGate_multi
        sc@meta.data$Signature <- ifelse( sc@meta.data$Signature=="NA",NA, sc@meta.data$Signature)
      }
      else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_BDrhapsody_MM.FP.signatures <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$BDrhapsody_scGATE.MM.signatures==T) {
        scGATE_anno_BD_MM.FP_signatures()
      }
      else{
        print("MM signatures not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })

    scGATE_anno_BD_MM.FP_Innate.NK <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )

      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }

      if (input$BDrhapsody_scGATE.MM.Innate.NK==T)  {
        models.list <- list()
        my_scGate_model <- gating_model(name = "Gzma",level = 1, signature = c("Gzma"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Gzmb",level = 2, signature = c("Gzmb"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Prf1",level = 3, signature = c("Prf1"), positive = T)
        models.list$Eff <- my_scGate_model

        my_scGate_model <- gating_model(name = "Sell",level = 1, signature = c("Sell"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Ccr7",level = 2, signature = c("Ccr7"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Cxcr6",level = 3, signature = c("Cxcr6"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "Itga1",level = 4, signature = c("Itga1"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "Prdm1",level = 5, signature = c("Prdm1"), positive = F)
        models.list$Recirculation <- my_scGate_model

        my_scGate_model <- gating_model(name = "Sell",level = 1, signature = c("Sell"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "Ccr7",level = 2, signature = c("Ccr7"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "Cxcr6",level = 3, signature = c("Cxcr6"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Itga1",level = 4, signature = c("Itga1"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Prdm1",level = 5, signature = c("Prdm1"), positive = T)
        models.list$Tiss_res <- my_scGate_model

        my_scGate_model <- gating_model(name = "Il2ra",level = 2, signature = c("Il2ra"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Foxp3",level = 2, signature = c("Foxp3"), positive = T)
        models.list$Treg <- my_scGate_model

        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$Innate.NK <- obj@meta.data$scGate_multi
        sc@meta.data$Innate.NK <- ifelse( sc@meta.data$Innate.NK=="NA",NA, sc@meta.data$Innate.NK)
      }
      else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_BDrhapsody_MM.FP.Innate.NK <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$BDrhapsody_scGATE.MM.Innate.NK==T) {
        scGATE_anno_BD_MM.FP_Innate.NK()
      }
      else{
        print("MM Innate.NK not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })

    scGATE_anno_BD_MM.FP_TNF.IFNg <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )

      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }

      if (input$BDrhapsody_scGATE.MM.TNF.IFNg==T)  {
        models.list <- list()
        models.list

        my_scGate_model <- gating_model(name = "Ifng",level = 1, signature = c("Ifng"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Tnf",level = 2, signature = c("Tnf"), positive = T)
        models.list$IFNG.TNF <- my_scGate_model

        my_scGate_model <- gating_model(name = "Ifng",level = 1, signature = c("Ifng"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Tnf",level = 2, signature = c("Tnf"), positive = F)
        models.list$IFNG <- my_scGate_model

        my_scGate_model <- gating_model(name = "Ifng",level = 1, signature = c("Ifng"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "Tnf",level = 2, signature = c("Tnf"), positive = T)
        models.list$TNF <- my_scGate_model


        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )

        sc@meta.data$TNF.IFNg <- obj@meta.data$scGate_multi
        sc@meta.data$TNF.IFNg <- ifelse( sc@meta.data$TNF.IFNg=="NA",NA, sc@meta.data$TNF.IFNg)
      }
      else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_BDrhapsody_MM.FP.TNF.IFNg <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$BDrhapsody_scGATE.MM.TNF.IFNg==T) {
        scGATE_anno_BD_MM.FP_TNF.IFNg()
      }
      else{
        print("MM TNF.IFNg not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })

    scGATE_anno_BD_MM.FP_subtypes  <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )

      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }

      if (input$BDrhapsody_scGATE.MM.subtypes==T)  {
        models.list <- list()
        my_scGate_model <- gating_model(name = "Ccr4",level = 1, signature = c("Ccr4"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Foxp3",level = 2, signature = c("Foxp3"), positive = F)
        models.list$Th2_like <- my_scGate_model

        my_scGate_model <- gating_model(name = "Cxcr3",level = 1, signature = c("Cxcr3"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Tbx21",level = 2, signature = c("Tbx21"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "CD8",level = 3, signature = c("Cd8a","Cd8b1"), positive = F)
        models.list$Th1_like <- my_scGate_model
        my_scGate_model <- gating_model(name = "RORC",level = 1, signature = c("Rorc"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Foxp3",level = 2, signature = c("Foxp3"), positive = F)
        models.list$Th17_like <- my_scGate_model
        my_scGate_model <- gating_model(name = "Cxcr5",level = 1, signature = c("Cxcr5"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Foxp3",level = 2, signature = c("Foxp3"), positive = F)
        models.list$Tfh_like <- my_scGate_model
        my_scGate_model <- gating_model(name = "Foxp3",level = 1, signature = c("Foxp3"), positive = T)
        models.list$Treg_like <- my_scGate_model

        my_scGate_model <- gating_model(name = "Klrk1",level = 1, signature = c("Klrk1"), positive = T) # https://pubmed.ncbi.nlm.nih.gov/36705564/
        my_scGate_model <- gating_model(my_scGate_model,name = "Il7r",level = 2, signature = c("Il7r"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Cxcr3",level = 3, signature = c("Cxcr3"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "Tbx21",level = 3, signature = c("Tbx21"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "Cd8a",level =4, signature = c("Cd8a"), positive = T)
        models.list$KILR_CD8_effector <- my_scGate_model

        my_scGate_model <- gating_model(name = "Gzma",level = 1, signature = c("Gzma"), positive = T) # https://pubmed.ncbi.nlm.nih.gov/36705564/
        my_scGate_model <- gating_model(my_scGate_model,name = "Gzmk",level = 2, signature = c("Gzmk"), positive = T)
        models.list$Effector_CD8 <- my_scGate_model

        my_scGate_model <- gating_model(name = "CD3e",level = 1, signature = c("Cd3e"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "Cd14",level = 2, signature = c("Cd14"), positive = T)
        models.list$Mono.Tcell.Complex <- my_scGate_model

        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$subtypes <- obj@meta.data$scGate_multi
        sc@meta.data$subtypes <- ifelse( sc@meta.data$subtypes=="NA",NA, sc@meta.data$subtypes)
      }
      else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_BDrhapsody_MM.FP.subtypes <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$BDrhapsody_scGATE.MM.subtypes==T) {
        scGATE_anno_BD_MM.FP_subtypes()
      }
      else{
        print("MM subtypes not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })

    scGATE_anno_BD_MM.FP_other  <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )

      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len>2000){
        len = 2000
      }

      if (input$BDrhapsody_scGATE.MM.other==T)  {
        models.list <- list()
        my_scGate_model <- gating_model(name = "Pdcd1",level = 1, signature = c("Pdcd1","Tigit"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "B3GAT1",level = 2, signature = c("B3gat1","Klrg1"), positive = T)
        models.list$exhausted_and_senescent <- my_scGate_model

        my_scGate_model <- gating_model(name = "PDCD1",level = 1, signature = c("Pdcd1","Tigit"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model,name = "B3GAT1",level = 2, signature = c("B3gat1","Klrg1"), positive = F)
        models.list$exhausted <- my_scGate_model

        my_scGate_model <- gating_model(name = "PDCD1",level = 1, signature = c("Pdcd1","Tigit"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model,name = "B3GAT1",level = 2, signature = c("B3gat1","Klrg1"), positive = T)
        models.list$senescent <- my_scGate_model
        my_scGate_model <- gating_model(name = "MKI67",level = 1, signature = c("Mki67","Top2a"), positive = T)
        models.list$CellCycling <- my_scGate_model

        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$other <- obj@meta.data$scGate_multi
        sc@meta.data$other <- ifelse( sc@meta.data$other=="NA",NA, sc@meta.data$other)
      }
      else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_BDrhapsody_MM.FP.other <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$BDrhapsody_scGATE.MM.other==T) {
        scGATE_anno_BD_MM.FP_other()
      }
      else{
        print("MM other not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    # BD rhapsody gates Mouse immune panel -----

    # verbatium outputs -----
    output$scGATE_verbatum_Generic <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$generic_scGATE==T) {
        scGATE_anno_generic()
      }

      else{
        print("Generic not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")

    })
    output$scGATE_verbatum_CD4 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$CD4_scGATE==T) {
        scGATE_anno_CD4()}
      else{
        print("CD4 not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    output$scGATE_verbatum_CD8 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$CD8_scGATE==T) {
        scGATE_anno_CD8()}
      else{
        print("CD8 not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    output$scGATE_verbatum_GeneralMarkers <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$GeneralMarkers_scGATE==T) {
        scGATE_anno_GeneralMarkers()
      }
      else{
        print("General Markers not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    output$scGATE_verbatum_Function <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$Function_scGATE==T) {
        scGATE_anno_Function()
      }
      else{
        print("Function (e.g. Th1) not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    output$scGATE_verbatum_Ex.Sen <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$Ex.Sen_scGATE==T ) {
        scGATE_anno_Ex.Sen()
      }
      else{
        print("Exhausted/Senescence not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    output$scGATE_verbatum_COVID <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$COVID_scGATE==T) {
        scGATE_anno_COVID()
      }
      else{
        print("COVID not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    output$scGATE_verbatum_Activation <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$Activation_scGATE==T) {
        scGATE_anno_Activation()
      }
      else{
        print("Activation not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    output$scGATE_verbatum_IFNgTNFa <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$IFNgTNFa_scGATE==T) {
        scGATE_anno_IFNgTNFa()
      }
      else{
        print("IFNg/TNFa not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    output$scGATE_verbatum_GNLY.PFR1.GZMB <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$GNLY.PFR1.GZMB_scGATE==T) {
        scGATE_anno_GNLY.PFR1.GZMB()
      }
      else{
        print("GNLY.PFR1.GZMB not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    output$scGATE_verbatum_Interlukin <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$Interlukin_scGATE==T) {
        scGATE_anno_Interlukin()
      }
      else{
        print("Interlukins not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    # output$scGATE_verbatum_BDrhapsody_scGATE <- renderPrint({
    #   FN <- tempfile()
    #   zz <- file(FN, open = "wt")
    #   sink(zz ,type = "output")
    #   sink(zz, type = "message")
    #   if (input$BDrhapsody_scGATE==T) {
    #     scGATE_anno_BD_rhapsody()
    #   }
    #   else{
    #     print("BD Rhapsody (homo Sapian) not run")
    #   }
    #   sink(type = "message")
    #   sink(type = "output")
    #   cat(readLines(FN), sep="\n")
    # })

    # scGATE user_Custom -------
    # geneset1
    scGate_anno_GeneSet1 <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      if (input$GeneSet1_scGate==T) {
        req(input$threshold_scGate)
        req(scGate_models_DB_geneset1)
        len <- length(rownames(sc@assays$RNA$scale.data))
        if (len>2000){
          len = 2000
        }
        models.list <- scGate_models_DB_geneset1
        obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$geneSet1 <- obj@meta.data$scGate_multi

        sc
      }
      else {
        sc
      }
      sc
    })

    create_UMAP_custom_1 <- reactive({
      sc <- scGate_anno_GeneSet1()

      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet1", pt.size = 1)
    })
    output$create_custom_1<- renderPlot({
      create_UMAP_custom_1()
    })


    output$scGATE_verbatum_GeneSet1 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$GeneSet1_scGate==T) {
        scGate_anno_GeneSet1()
      }
      else{
        print("Geneset1 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    # geneset2
    scGate_anno_GeneSet2 <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      if (input$GeneSet2_scGate==T) {
        req(input$threshold_scGate)
        req(scGate_models_DB_geneset2)
        len <- length(rownames(sc@assays$RNA$scale.data))
        if (len>2000){
          len = 2000
        }
        models.list <- scGate_models_DB_geneset2
        obj <- obj <- obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$geneSet2 <- obj@meta.data$scGate_multi

        sc
      }
      else {
        sc
      }
      sc
    })

    create_UMAP_custom_2 <- reactive({
      sc <- scGate_anno_GeneSet2()

      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet2", pt.size = 1)
    })
    output$create_custom_2<- renderPlot({
      create_UMAP_custom_2()
    })

    output$scGATE_verbatum_GeneSet2 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$GeneSet2_scGate==T) {
        scGate_anno_GeneSet2()
      }
      else{
        print("Geneset2 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    # geneset3
    scGate_anno_GeneSet3 <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      if (input$GeneSet3_scGate==T) {
        req(input$threshold_scGate)
        req(scGate_models_DB_geneset3)
        len <- length(rownames(sc@assays$RNA$scale.data))
        if (len>2000){
          len = 2000
        }
        models.list <- scGate_models_DB_geneset3
        obj <- obj <- obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$geneSet3 <- obj@meta.data$scGate_multi

        sc
      }
      else {
        sc
      }
      sc
    })

    create_UMAP_custom_3 <- reactive({
      sc <- scGate_anno_GeneSet3()

      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet3", pt.size = 1)
    })
    output$create_custom_3<- renderPlot({
      create_UMAP_custom_3()
    })

    output$scGATE_verbatum_GeneSet3 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$GeneSet3_scGate==T) {
        scGate_anno_GeneSet3()
      }
      else{
        print("Geneset3 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    # geneset4
    scGate_anno_GeneSet4 <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      if (input$GeneSet4_scGate==T) {
        req(input$threshold_scGate)
        len <- length(rownames(sc@assays$RNA$scale.data))
        if (len>2000){
          len = 2000
        }
        models.list <- scGate_models_DB_geneset4
        obj <- obj <- obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$geneSet4 <- obj@meta.data$scGate_multi

        sc
      }
      else {
        sc
      }
      sc
    })

    create_UMAP_custom_4 <- reactive({
      sc <- scGate_anno_GeneSet4()

      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet4", pt.size = 1)
    })
    output$create_custom_4<- renderPlot({
      create_UMAP_custom_4()
    })

    output$scGATE_verbatum_GeneSet4 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$GeneSet4_scGate==T) {
        scGate_anno_GeneSet4()
      }
      else{
        print("Geneset4 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    # geneset5
    scGate_anno_GeneSet5 <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      if (input$GeneSet5_scGate==T) {
        req(input$threshold_scGate)
        len <- length(rownames(sc@assays$RNA$scale.data))
        if (len>2000){
          len = 2000
        }
        models.list <- scGate_models_DB_geneset5
        obj <- obj <- obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$geneSet5 <- obj@meta.data$scGate_multi

        sc
      }
      else {
        sc
      }
      sc
    })

    create_UMAP_custom_5 <- reactive({
      sc <- scGate_anno_GeneSet5()

      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet5", pt.size = 1)
    })
    output$create_custom_5<- renderPlot({
      create_UMAP_custom_5()
    })

    output$scGATE_verbatum_GeneSet5 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$GeneSet5_scGate==T) {
        scGate_anno_GeneSet5()
      }
      else{
        print("Geneset5 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    # geneset6
    scGate_anno_GeneSet6 <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      if (input$GeneSet6_scGate==T) {
        req(input$threshold_scGate)
        len <- length(rownames(sc@assays$RNA$scale.data))
        if (len>2000){
          len = 2000
        }
        models.list <- scGate_models_DB_geneset6
        obj <- obj <- obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$geneSet6 <- obj@meta.data$scGate_multi

        sc
      }
      else {
        sc
      }
      sc
    })

    create_UMAP_custom_6 <- reactive({
      sc <- scGate_anno_GeneSet6()

      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet6", pt.size = 1)
    })
    output$create_custom_6<- renderPlot({
      create_UMAP_custom_6()
    })

    output$scGATE_verbatum_GeneSet6 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$GeneSet6_scGate==T) {
        scGate_anno_GeneSet6()
      }
      else{
        print("Geneset6 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    # geneset7
    scGate_anno_GeneSet7 <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      if (input$GeneSet7_scGate==T) {
        req(input$threshold_scGate)
        len <- length(rownames(sc@assays$RNA$scale.data))
        if (len>2000){
          len = 2000
        }
        models.list <- scGate_models_DB_geneset7
        obj <- obj <- obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$geneSet7 <- obj@meta.data$scGate_multi

        sc
      }
      else {
        sc
      }
      sc
    })

    create_UMAP_custom_7 <- reactive({
      sc <- scGate_anno_GeneSet7()

      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet7", pt.size = 1)
    })
    output$create_custom_7<- renderPlot({
      create_UMAP_custom_7()
    })

    output$scGATE_verbatum_GeneSet7 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$GeneSet7_scGate==T) {
        scGate_anno_GeneSet7()
      }
      else{
        print("Geneset7 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    # geneset8
    scGate_anno_GeneSet8 <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      if (input$GeneSet8_scGate==T) {
        req(input$threshold_scGate)
        len <- length(rownames(sc@assays$RNA$scale.data))
        if (len>2000){
          len = 2000
        }
        models.list <- scGate_models_DB_geneset8
        obj <- obj <- obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$geneSet8 <- obj@meta.data$scGate_multi

        sc
      }
      else {
        sc
      }
      sc
    })

    create_UMAP_custom_8 <- reactive({
      sc <- scGate_anno_GeneSet8()

      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet8", pt.size = 1)
    })
    output$create_custom_8<- renderPlot({
      create_UMAP_custom_8()
    })
    output$scGATE_verbatum_GeneSet8 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$GeneSet8_scGate==T) {
        scGate_anno_GeneSet8()
      }
      else{
        print("Geneset8 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    # geneset8
    scGate_anno_GeneSet9 <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      if (input$GeneSet9_scGate==T) {
        req(input$threshold_scGate)
        len <- length(rownames(sc@assays$RNA$scale.data))
        if (len>2000){
          len = 2000
        }
        models.list <- scGate_models_DB_geneset9
        obj <- obj <- obj <- scGate(sc, model = models.list,pos.thr = input$threshold_scGate,neg.thr = input$threshold_scGate,nfeatures = len,ncores = 8 )
        sc@meta.data$geneSet9 <- obj@meta.data$scGate_multi

        sc
      }
      else {
        sc
      }
      sc
    })

    create_UMAP_custom_9 <- reactive({
      sc <- scGate_anno_GeneSet9()

      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet9", pt.size = 1)
    })
    output$create_custom_9<- renderPlot({
      create_UMAP_custom_9()
    })

    output$scGATE_verbatum_GeneSet9 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$GeneSet9_scGate==T) {
        scGate_anno_GeneSet9()
      }
      else{
        print("Geneset9 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })


    # creating the final object ------
    scGATE_anno_HS.10x <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      if (input$generic_scGATE==T) {
        obj <- scGATE_anno_generic()
        sc@meta.data$generic <- obj@meta.data$generic
      }
      if (input$CD4_scGATE==T) {
        obj <- scGATE_anno_CD4()
        sc@meta.data$CD4_TIL <- obj@meta.data$CD4_TIL
      }
      if (input$CD8_scGATE==T) {
        obj <- scGATE_anno_CD8()
        sc@meta.data$CD8_TIL <- obj@meta.data$CD8_TIL
      }
      if (input$GeneralMarkers_scGATE==T) {
        obj <- scGATE_anno_GeneralMarkers()
        sc@meta.data$GeneralMarkers <- obj@meta.data$GeneralMarkers
      }
      if (input$Function_scGATE==T) {
        obj <- scGATE_anno_Function()
        sc@meta.data$Function <- obj@meta.data$Function
      }
      if (input$Ex.Sen_scGATE==T) {
        obj <- scGATE_anno_Ex.Sen()
        sc@meta.data$Exhausted_Senescence <- obj@meta.data$Exhausted_Senescence
      }
      if (input$COVID_scGATE==T) {
        obj <- scGATE_anno_COVID()
        sc@meta.data$COVID <- obj@meta.data$COVID
      }
      if (input$Activation_scGATE==T) {
        obj <- scGATE_anno_Activation()
        sc@meta.data$Activation <- obj@meta.data$Activation
      }
      if (input$IFNgTNFa_scGATE==T) {
        obj <- scGATE_anno_IFNgTNFa()
        sc@meta.data$IFNgTNFa <- obj@meta.data$IFNgTNFa
      }
      if (input$GNLY.PFR1.GZMB_scGATE==T) {
        obj <- scGATE_anno_GNLY.PFR1.GZMB()
        sc@meta.data$GNLY.PFR1.GZMB <- obj@meta.data$GNLY.PFR1.GZMB
      }
      if (input$Interlukin_scGATE==T) {
        obj <- scGATE_anno_Interlukin()
        sc@meta.data$Interleukins <- obj@meta.data$Interleukins
      }
      if (input$GeneSet1_scGate==T) {
        obj <- scGate_anno_GeneSet1()
        sc@meta.data$geneSet1 <- obj@meta.data$geneSet1
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet1"] <- input$geneset1_name
      }
      if (input$GeneSet2_scGate==T) {
        obj <- scGate_anno_GeneSet2()
        sc@meta.data$geneSet2 <- obj@meta.data$geneSet2
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet2"] <- input$geneset2_name
      }
      if (input$GeneSet3_scGate==T) {
        obj <- scGate_anno_GeneSet3()
        sc@meta.data$geneSet3 <- obj@meta.data$geneSet3
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet3"] <- input$geneset3_name
      }
      if (input$GeneSet4_scGate==T) {
        obj <- scGate_anno_GeneSet4()
        sc@meta.data$geneSet4 <- obj@meta.data$geneSet4
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet4"] <- input$geneset4_name
      }
      if (input$GeneSet5_scGate==T) {
        obj <- scGate_anno_GeneSet5()
        sc@meta.data$geneSet5 <- obj@meta.data$geneSet5
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet5"] <- input$geneset5_name
      }

      if (input$GeneSet6_scGate==T) {
        obj <- scGate_anno_GeneSet6()
        sc@meta.data$geneSet6 <- obj@meta.data$geneSet6
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet6"] <- input$geneset6_name
      }
      if (input$GeneSet7_scGate==T) {
        obj <- scGate_anno_GeneSet7()
        sc@meta.data$geneSet7 <- obj@meta.data$geneSet7
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet7"] <- input$geneset7_name
      }
      if (input$GeneSet8_scGate==T) {
        obj <- scGate_anno_GeneSet8()
        sc@meta.data$geneSet8 <- obj@meta.data$geneSet8
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet8"] <- input$geneset8_name
      }
      if (input$GeneSet9_scGate==T) {
        obj <- scGate_anno_GeneSet9()
        sc@meta.data$geneSet9 <- obj@meta.data$geneSet9
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet9"] <- input$geneset9_name
      }
      else {
        sc
      }

      sc

    })
    scGATE_anno_BD_HS.IP <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )

      sc <- scGATE_anno_BD_rhapsody()
      if (input$GeneSet1_scGate==T) {
        obj <- scGate_anno_GeneSet1()
        sc@meta.data$geneSet1 <- obj@meta.data$geneSet1
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet1"] <- input$geneset1_name
      }
      if (input$GeneSet2_scGate==T) {
        obj <- scGate_anno_GeneSet2()
        sc@meta.data$geneSet2 <- obj@meta.data$geneSet2
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet2"] <- input$geneset2_name
      }
      if (input$GeneSet3_scGate==T) {
        obj <- scGate_anno_GeneSet3()
        sc@meta.data$geneSet3 <- obj@meta.data$geneSet3
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet3"] <- input$geneset3_name
      }
      if (input$GeneSet4_scGate==T) {
        obj <- scGate_anno_GeneSet4()
        sc@meta.data$geneSet4 <- obj@meta.data$geneSet4
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet4"] <- input$geneset4_name
      }
      if (input$GeneSet5_scGate==T) {
        obj <- scGate_anno_GeneSet5()
        sc@meta.data$geneSet5 <- obj@meta.data$geneSet5
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet5"] <- input$geneset5_name
      }
      if (input$GeneSet6_scGate==T) {
        obj <- scGate_anno_GeneSet6()
        sc@meta.data$geneSet6 <- obj@meta.data$geneSet6
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet6"] <- input$geneset6_name
      }
      if (input$GeneSet7_scGate==T) {
        obj <- scGate_anno_GeneSet7()
        sc@meta.data$geneSet7 <- obj@meta.data$geneSet7
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet7"] <- input$geneset7_name
      }
      if (input$GeneSet8_scGate==T) {
        obj <- scGate_anno_GeneSet8()
        sc@meta.data$geneSet8 <- obj@meta.data$geneSet8
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet8"] <- input$geneset8_name
      }
      if (input$GeneSet9_scGate==T) {
        obj <- scGate_anno_GeneSet9()
        sc@meta.data$geneSet9 <- obj@meta.data$geneSet9
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet9"] <- input$geneset9_name
      }

      sc

    })
    scGATE_anno_BD_HS.FP <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )

      if (input$GeneSet1_scGate==T) {
        obj <- scGate_anno_GeneSet1()
        sc@meta.data$geneSet1 <- obj@meta.data$geneSet1
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet1"] <- input$geneset1_name
      }

      if (input$GeneSet2_scGate==T) {
        obj <- scGate_anno_GeneSet2()
        sc@meta.data$geneSet2 <- obj@meta.data$geneSet2
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet2"] <- input$geneset2_name
      }

      if (input$GeneSet3_scGate==T) {
        obj <- scGate_anno_GeneSet3()
        sc@meta.data$geneSet3 <- obj@meta.data$geneSet3
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet3"] <- input$geneset3_name
      }
      if (input$GeneSet4_scGate==T) {
        obj <- scGate_anno_GeneSet4()
        sc@meta.data$geneSet4 <- obj@meta.data$geneSet4
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet4"] <- input$geneset4_name
      }

      if (input$GeneSet5_scGate==T) {
        obj <- scGate_anno_GeneSet5()
        sc@meta.data$geneSet5 <- obj@meta.data$geneSet5
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet5"] <- input$geneset5_name
      }

      if (input$GeneSet6_scGate==T) {
        obj <- scGate_anno_GeneSet6()
        sc@meta.data$geneSet6 <- obj@meta.data$geneSet6
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet6"] <- input$geneset6_name
      }
      if (input$GeneSet7_scGate==T) {
        obj <- scGate_anno_GeneSet7()
        sc@meta.data$geneSet7 <- obj@meta.data$geneSet7
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet7"] <- input$geneset7_name
      }
      if (input$GeneSet8_scGate==T) {
        obj <- scGate_anno_GeneSet8()
        sc@meta.data$geneSet8 <- obj@meta.data$geneSet8
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet8"] <- input$geneset8_name
      }
      if (input$GeneSet9_scGate==T) {
        obj <- scGate_anno_GeneSet9()
        sc@meta.data$geneSet9 <- obj@meta.data$geneSet9
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet9"] <- input$geneset9_name
      }

      sc

    })
    scGATE_anno_BD_MM.FP <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      if (input$BDrhapsody_scGATE.MM.Tcell==T) {
        obj <- scGATE_anno_BD_MM.FP_T.cell()
        sc@meta.data$T_cells <- obj@meta.data$T_cells
      }
      if (input$BDrhapsody_scGATE.MM.Memory==T) {
        obj <- scGATE_anno_BD_MM.FP_Memory()
        sc@meta.data$Memory <- obj@meta.data$Memory
      }
      if (input$BDrhapsody_scGATE.MM.signatures ==T) {
        obj <- scGATE_anno_BD_MM.FP_signatures()
        sc@meta.data$signatures <- obj@meta.data$signatures
      }
      if (input$BDrhapsody_scGATE.MM.Innate.NK ==T) {
        obj <- scGATE_anno_BD_MM.FP_Innate.NK()
        sc@meta.data$Innate.NK <- obj@meta.data$Innate.NK
      }
      if (input$BDrhapsody_scGATE.MM.TNF.IFNg ==T) {
        obj <- scGATE_anno_BD_MM.FP_TNF.IFNg()
        sc@meta.data$TNF.IFNg <- obj@meta.data$TNF.IFNg
      }
      if (input$BDrhapsody_scGATE.MM.subtypes ==T) {
        obj <- scGATE_anno_BD_MM.FP_subtypes()
        sc@meta.data$subtypes <- obj@meta.data$subtypes
      }
      if (input$BDrhapsody_scGATE.MM.other ==T) {
        obj <- scGATE_anno_BD_MM.FP_other()
        sc@meta.data$other <- obj@meta.data$other
      }

      if (input$GeneSet1_scGate==T) {
        obj <- scGate_anno_GeneSet1()
        sc@meta.data$geneSet1 <- obj@meta.data$geneSet1
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet1"] <- input$geneset1_name
      }
      if (input$GeneSet2_scGate==T) {
        obj <- scGate_anno_GeneSet2()
        sc@meta.data$geneSet2 <- obj@meta.data$geneSet2
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet2"] <- input$geneset2_name
      }
      if (input$GeneSet3_scGate==T) {
        obj <- scGate_anno_GeneSet3()
        sc@meta.data$geneSet3 <- obj@meta.data$geneSet3
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet3"] <- input$geneset3_name
      }
      if (input$GeneSet4_scGate==T) {
        obj <- scGate_anno_GeneSet4()
        sc@meta.data$geneSet4 <- obj@meta.data$geneSet4
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet4"] <- input$geneset4_name
      }
      if (input$GeneSet5_scGate==T) {
        obj <- scGate_anno_GeneSet5()
        sc@meta.data$geneSet5 <- obj@meta.data$geneSet5
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet5"] <- input$geneset5_name
      }
      if (input$GeneSet6_scGate==T) {
        obj <- scGate_anno_GeneSet6()
        sc@meta.data$geneSet6 <- obj@meta.data$geneSet6
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet6"] <- input$geneset6_name
      }
      if (input$GeneSet7_scGate==T) {
        obj <- scGate_anno_GeneSet7()
        sc@meta.data$geneSet7 <- obj@meta.data$geneSet7
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet7"] <- input$geneset7_name
      }
      if (input$GeneSet8_scGate==T) {
        obj <- scGate_anno_GeneSet8()
        sc@meta.data$geneSet8 <- obj@meta.data$geneSet8
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet8"] <- input$geneset8_name
      }
      if (input$GeneSet9_scGate==T) {
        obj <- scGate_anno_GeneSet9()
        sc@meta.data$geneSet9 <- obj@meta.data$geneSet9
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet9"] <- input$geneset9_name
      }

      else {
        sc
      }

      sc
    })
    scGATE_anno_BD_MM.IP <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )

      if (input$GeneSet1_scGate==T) {
        obj <- scGate_anno_GeneSet1()
        sc@meta.data$geneSet1 <- obj@meta.data$geneSet1
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet1"] <- input$geneset1_name
      }
      if (input$GeneSet2_scGate==T) {
        obj <- scGate_anno_GeneSet2()
        sc@meta.data$geneSet2 <- obj@meta.data$geneSet2
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet2"] <- input$geneset2_name
      }
      if (input$GeneSet3_scGate==T) {
        obj <- scGate_anno_GeneSet3()
        sc@meta.data$geneSet3 <- obj@meta.data$geneSet3
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet3"] <- input$geneset3_name
      }
      if (input$GeneSet4_scGate==T) {
        obj <- scGate_anno_GeneSet4()
        sc@meta.data$geneSet4 <- obj@meta.data$geneSet4
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet4"] <- input$geneset4_name
      }
      if (input$GeneSet5_scGate==T) {
        obj <- scGate_anno_GeneSet5()
        sc@meta.data$geneSet5 <- obj@meta.data$geneSet5
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet5"] <- input$geneset5_name
      }
      if (input$GeneSet6_scGate==T) {
        obj <- scGate_anno_GeneSet6()
        sc@meta.data$geneSet6 <- obj@meta.data$geneSet6
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet6"] <- input$geneset6_name
      }
      if (input$GeneSet7_scGate==T) {
        obj <- scGate_anno_GeneSet7()
        sc@meta.data$geneSet7 <- obj@meta.data$geneSet7
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet7"] <- input$geneset7_name
      }
      if (input$GeneSet8_scGate==T) {
        obj <- scGate_anno_GeneSet8()
        sc@meta.data$geneSet8 <- obj@meta.data$geneSet8
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet8"] <- input$geneset8_name
      }
      if (input$GeneSet9_scGate==T) {
        obj <- scGate_anno_GeneSet9()
        sc@meta.data$geneSet9 <- obj@meta.data$geneSet9
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet9"] <- input$geneset9_name
      }

      else {
        sc
      }

      sc
    })

    # all.annotations added -----
    output$DEx_table_TcellClass_scGATE <-  DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),
                                                               options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100),
                                                                              pageLength = 5, scrollX = TRUE),{
                                                                                # c("10x_HS","BD_HS.Immune.Panel","BD_HS.Full.Panel","10x_MM","BD_MM_Full.Panel","BD_MM_Immune.Panel",)
                                                                                sc <- getData_2()
                                                                                validate(
                                                                                  need(nrow(sc)>0,
                                                                                       "Upload file for annotation")
                                                                                )
                                                                                if (input$Data_types=="10x_HS") {
                                                                                  sc <- scGATE_anno_HS.10x()
                                                                                  as.data.frame(sc@meta.data)
                                                                                }
                                                                                else if (input$Data_types=="BD_HS.Immune.Panel") {
                                                                                  sc <- scGATE_anno_BD_HS.IP()
                                                                                  as.data.frame(sc@meta.data)
                                                                                }
                                                                                else if (input$Data_types=="BD_HS.Full.Panel") {
                                                                                  sc <- scGATE_anno_HS.10x()
                                                                                  as.data.frame(sc@meta.data)
                                                                                }

                                                                                else if (input$Data_types=="10x_MM") {
                                                                                  sc <- scGATE_anno_BD_MM.FP()
                                                                                  as.data.frame(sc@meta.data)
                                                                                }

                                                                                else if (input$Data_types=="BD_MM_Full.Panel") {
                                                                                  sc <- scGATE_anno_BD_MM.FP()
                                                                                  as.data.frame(sc@meta.data)
                                                                                }

                                                                                else if (input$Data_types=="BD_MM_Immune.Panel") {
                                                                                  sc <- scGATE_anno_BD_MM.IP()
                                                                                  as.data.frame(sc@meta.data)
                                                                                }

                                                                                else if (input$Data_types=="TCR-seq") {
                                                                                  sc <- TCR_seq_classification()
                                                                                  as.data.frame(sc@meta.data)

                                                                                }


                                                                                else {
                                                                                  df <- as.data.frame("Other panels are under development")
                                                                                  names(df) <- "V1"
                                                                                  df
                                                                                }

                                                                              })

    output$downloaddf_SeruatObj_annotated <- downloadHandler(
      filename = function(){
        x = today()
        paste(input$project_name3,"_annotated_",x,".rds", sep = "")
      },
      content = function(file){
        if (input$Data_types=="10x_HS") {
          sc <- scGATE_anno_HS.10x()
          SaveSeuratRds(sc,file)
        }
        else if (input$Data_types=="BD_HS.Immune.Panel") {
          sc <- scGATE_anno_BD_HS.IP()
          SaveSeuratRds(sc,file)
        }

        else if (input$Data_types=="BD_HS.Full.Panel") {
          sc <- scGATE_anno_BD_HS.FP()
          SaveSeuratRds(sc,file)
        }
        # mouse panels
        else if (input$Data_types=="10x_MM") {
          sc <- scGATE_anno_BD_MM.FP()
          SaveSeuratRds(sc,file)
        }
        else if (input$Data_types=="BD_MM_Full.Panel") {
          sc <- scGATE_anno_BD_MM.FP()
          SaveSeuratRds(sc,file)
        }
        else if (input$Data_types=="BD_MM_Immune.Panel") {
          sc <- scGATE_anno_BD_MM.IP()
          SaveSeuratRds(sc,file)
        }

        else if (input$Data_types=="TCR-seq") {
          sc <- TCR_seq_classification()
          SaveSeuratRds(sc,file)
        }

        else {

        }
      })

    # output$downloaddf_SeruatObj_annotated <- downloadHandler(
    #   filename = function(){
    #     x = today
    #     paste(input$project_name3,"_annotated_",gsub("-", ".", Sys.Date()),".h5Seurat", sep = "")
    #   },
    #   content = function(file){
    #     if (input$Data_types=="10x_HS") {
    #       sc <- scGATE_anno_HS.10x()
    #       as.h5Seurat(sc,file)
    #     }
    #     else if (input$Data_types=="BD_HS.Immune.Panel") {
    #       sc <- scGATE_anno_BD_HS.IP()
    #       as.h5Seurat(sc,file)
    #     }
    #
    #     else if (input$Data_types=="BD_HS.Full.Panel") {
    #       sc <- scGATE_anno_BD_HS.FP()
    #       as.h5Seurat(sc,file)
    #     }
    #     # mouse panels
    #     else if (input$Data_types=="10x_MM") {
    #       sc <- scGATE_anno_BD_MM.FP()
    #       as.h5Seurat(sc,file)
    #     }
    #     else if (input$Data_types=="BD_MM_Full.Panel") {
    #       sc <- scGATE_anno_BD_MM.FP()
    #       as.h5Seurat(sc,file)
    #     }
    #     else if (input$Data_types=="BD_MM_Immune.Panel") {
    #       sc <- scGATE_anno_BD_MM.IP()
    #       as.h5Seurat(sc,file)
    #     }
    #
    #     else if (input$Data_types=="TCR-seq") {
    #       sc <- TCR_seq_classification()
    #       as.h5Seurat(sc,file)
    #     }
    #
    #     else {
    #
    #     }
    #   } )

    #

    ## differential expression -----
    observe({
      df.test <- getData_2()
      validate(
        need(nrow(df.test)>0,
             error_message_val_sc)
      )
      updateSelectInput(
        session,
        "meta_data_sc_clust",
        choices=names(df.test@meta.data),
        selected = c("seurat_clusters"))
    })
    observe({
      df.test <- getData_2()
      validate(
        need(nrow(df.test)>0,
             error_message_val_sc)
      )
      updateSelectInput(
        session,
        "meta_data_sc_",
        choices=names(df.test@meta.data),
        selected = c("orig.ident"))
    })
    df_sc_clust <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )

      if (input$multiple_group_sc =="yes") {
        Idents(sc) <- paste(sc@meta.data[[input$meta_data_sc_clust]],
                            sc@meta.data[[input$meta_data_sc_]] )
      }
      else {
        Idents(sc) <- paste(sc@meta.data[[input$meta_data_sc_clust]])
      }

      sc

    })


    list_of_genes <- reactive({
      df.test <- getData_2()
      validate(
        need(nrow(df.test)>0,
             error_message_val_sc)
      )

      kmeans <- read.csv(system.file("Kmean","Kmeans.requires.annotation.csv",package = "STEGO.R"))
      var.genes <- as.data.frame(unique(df.test@assays$RNA@var.features))
      names(var.genes) <- "V1"

      if (input$sample.type.source.markers == 'hs') {

        kmeans2 <- as.data.frame(kmeans$Human)
        names(kmeans2) <- "V1"
      }
      else {

        kmeans2 <- as.data.frame(kmeans$Mouse)
        names(kmeans2) <- "V1"
      }
      kmeans2
      var.genes

      name.df <- unique(c(kmeans$V1,var.genes$V1))
      name.df <- as.data.frame(name.df)
      names(name.df) <- "Gene_Name"
      name.df <- as.data.frame(name.df[order(name.df$Gene_Name),])
      names(name.df) <- "Gene_Name"

      name.df

    })




    output$list_of_genes_df <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      df.test <- getData_2()
      validate(
        need(nrow(df.test)>0,
             error_message_val_sc)
      )
      list_of_genes()
    })


    output$meta_data_comaprison_check <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      calls <- sc@meta.data
      calls
    })

    vals_clust_markers <- reactiveValues(markers_for_table=NULL)
    observeEvent(input$run_differental.exp,{
      sc <- df_sc_clust()

      validate(
        need(nrow(sc)>0,
             "Run DEx")
      )


      sc.markers <- FindAllMarkers(df_sc_clust(),
                                   only.pos = TRUE,
                                   min.pct = input$min.ptc.sc,
                                   logfc.threshold = input$logfc.ptc.sc,
                                   test.use = input$normalN
      )
      vals_clust_markers$markers_for_table <- sc.markers
    })

    output$DEx_table_clusters <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- vals_clust_markers$markers_for_table
      calls
    })

    output$downloaddf_DEx_table_clusters <- downloadHandler(
      filename = function(){
        paste("DEx_all_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(vals_clust_markers$markers_for_table)
        write.csv(df,file, row.names = T)
      } )

    #
    # markers_heatmap <- reactive({
    #   sc.markers %>%
    #     group_by(cluster) %>%
    #     top_n(n = 5, wt = avg_log2FC) -> top10
    #   DoHeatmap(sc, features = top10$gene) + NoLegend()
    # })
    #

    ## View featurePlot markers ----
    observeEvent(input$run_string.data3,{
      df.test <- getData_2()
      validate(
        need(nrow(df.test)>0,
             error_message_val_sc)
      )
      name.df <- as.data.frame(list_of_genes())
      names(name.df) <- "Gene_Name"
      updateSelectInput(
        session,
        "string.data3",
        choices=name.df$Gene_Name,
        selected = c("GZMB","CD4","CD8A"))
    })

    #
    markers_featurePlot <- reactive({
      sc <- df_sc_clust()
      Feture_plots <- list()
      feature_name <- c(input$string.data3)
      req(input$string.data3)
      x <- length(feature_name)

      if (input$norm_expression_for_all=="yes") {

        for (i in 1:x) {
          Feture_plots[[i]] <- FeaturePlot(sc,features = feature_name[i],raster=FALSE,label=input$label_is_true_features) +
            scale_color_gradientn(colours = c(input$lower_col_FP,input$upper_col_FP),limits=c(0,input$max_norm_FP)) +
            theme(plot.background = element_rect(fill="white", color = NA)) +
            theme(plot.title = element_text(size=10))

        }
      }

      else {
        for (i in 1:x) {
          Feture_plots[[i]] <- FeaturePlot(sc,features = feature_name[i],raster=FALSE,label=input$label_is_true_features) +
            scale_color_gradientn(colours = c(input$lower_col_FP,input$upper_col_FP)) +
            theme(plot.background = element_rect(fill="white", color = NA)) +
            theme(plot.title = element_text(size=10))

        }
      }

      n <- length(Feture_plots)
      nCol <- round(sqrt(n),0)
      do.call("grid.arrange", c(Feture_plots, ncol=nCol))


    })
    output$markers_featurePlot_sc <- renderPlot({
      markers_featurePlot()
    })

    # download markers_featurePlot_sc
    output$downloadPlot_markers_featurePlot_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste("TCR_Explore_markers_featurePlot_sc_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_markers_featurePlot_sc,height=input$height_markers_featurePlot_sc, onefile = FALSE) # open the pdf device

        plot(markers_featurePlot())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_markers_featurePlot_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste("TCR_Explore_markers_featurePlot_sc_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_markers_featurePlot_sc,height = input$height_png_markers_featurePlot_sc,res = input$resolution_PNG_markers_featurePlot_sc)
        plot(markers_featurePlot())
        dev.off()},   contentType = "application/png")



    # Analysis -----
    observe({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      meta.data <- sc@meta.data
      updateSelectInput(
        session,
        "Split_group_by_overview",
        choices=names(meta.data),
        selected = "Sample_Name")
    })


    ## uploading seruat obj ----
    input.data_sc_pro <- reactive({
      inFile_sc_pro <- input$file_SC_pro
      if (is.null(inFile_sc_pro)) return(NULL)
      else {
        dataframe = LoadSeuratRds(inFile_sc_pro$datapath)
      }

    })
    input.data_sc_clusTCR_AG <- reactive({
      inFile_cluster_fileAG <- input$file_cluster_file_AG
      if (is.null(inFile_cluster_fileAG)) return(NULL)
      else {
        dataframe = read.csv(inFile_cluster_fileAG$datapath)
      }
    })


    input.data_sc_clusTCR_BD <- reactive({
      inFile_cluster_fileBD <- input$file_cluster_file_BD
      if (is.null(inFile_cluster_fileBD)) return(NULL)
      else {
        dataframe = read.csv(inFile_cluster_fileBD$datapath)
      }
    })


    input.data_sc_clusTCR_IgH <- reactive({
      inFile_cluster_fileIgH <- input$file_cluster_file_IgH
      if (is.null(inFile_cluster_fileIgH)) return(NULL)
      else {
        dataframe = read.csv(inFile_cluster_fileIgH$datapath)
      }
    })

    input.data_sc_clusTCR_IgKL <- reactive({
      inFile_cluster_fileIgKL <- input$file_cluster_file_IgKL
      if (is.null(inFile_cluster_fileIgKL)) return(NULL)
      else {
        dataframe = read.csv(inFile_cluster_fileIgH$datapath)
      }
    })

    ## Add in additional sample identifers - matched to Sample_Name or Orig.indent
    input.data_sc_Labels_to_add <- reactive({
      inFile_Labels_to_add <- input$file_Labels_to_add
      if (is.null(inFile_Labels_to_add)) return(NULL)
      else {
        dataframe = read.csv(inFile_Labels_to_add$datapath)
      }
    })
    ## add additional sample labels to file. -----
    UMAP_metadata_with_labs <- reactive({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      req(input$Samp_col2)
      sc@meta.data$order <- 1:dim(sc@meta.data)[1]
      reduction <- (sc@reductions$umap)
      UMAP <- as.data.frame(reduction@cell.embeddings)
      names(UMAP)[1:2] <- toupper(names(UMAP)[1:2])
      UMAP$Cell_Index <- rownames(UMAP)
      meta.data <- as.data.frame(sc@meta.data)
      meta.data <- meta.data
      umap.meta <- merge(UMAP,meta.data,by="Cell_Index",sort = F)

      if(input$add_additional_lables == "yes") {
        labs <- input.data_sc_Labels_to_add()
        validate(
          need(nrow(labs)>0,
               "Upload lab file")
        )
        names(umap.meta)[names(umap.meta) %in% input$Samp_col2] <- "ID"
        umap.meta <- merge(labs,umap.meta,by="ID", sort = F)
        names(umap.meta)[names(umap.meta) %in% "ID"] <- input$Samp_col2
      }
      rownames(umap.meta) <- umap.meta$Cell_Index
      umap.meta <- umap.meta[order(umap.meta$order),]
      sc@meta.data <- umap.meta
      sc

    })

    output$Sample_names_merging_sc <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      umap.meta <- sc@meta.data
      names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
      names(umap.meta)[names(umap.meta) %in% input$V_gene_sc] <- "v_gene_selected"
      if(nrow(umap.meta)>0) {
        umap.meta
      }

      else {
        as.data.frame("Processing")

      }
    })


    # checking issue with analysis ------
    output$meta.data_check_upload <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )

      umap.meta <- sc@meta.data
      names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
      names(umap.meta)[names(umap.meta) %in% input$V_gene_sc] <- "v_gene_selected"
      umap.meta
      req(TCR_Expanded())
      sc_merged <- merge(umap.meta,TCR_Expanded(),by=c("v_gene_selected","ID_Column"),all.x=T)
      sc_merged

      if(nrow(sc_merged)>0) {
        sc_merged
      }

      else {
        as.data.frame("Processing")

      }

    })


    # upload TCRex file ----
    input.data_sc_TCRex <- reactive({
      inupload_TCRex_file <- input$upload_TCRex_file
      if (is.null(inupload_TCRex_file)) return(NULL)
      else {
        dataframe = read.table(inupload_TCRex_file$datapath,skip = input$skip_TCRex_up,header = T,sep="\t")
      }
    })

    ## uploaded clusTCR table -----
    output$Tb_ClusTCR_test <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      AG_calls <- input.data_sc_clusTCR_AG()
      BD_calls <- input.data_sc_clusTCR_BD()
      validate(
        need(nrow(AG_calls)>0 & nrow(BD_calls)>0,
             "Upload clusTCR table, which is needed for TCR -> UMAP section")
      )
      rbind(AG_calls,BD_calls)
    })

    output$Tb_tcrex_test <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- input.data_sc_TCRex()
      validate(
        need(nrow(calls)>0,
             "Upload TCRex table")
      )
      calls
    })

    ### observe ----
    observe({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )

      df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "Samp_col",
        choices=names(df3.meta),
        selected = "Sample_Name")
    })


    #
    observe({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "V_gene_sc",
        choices=names(df3.meta),
        selected = "vdj_gene_cdr3_AG_BD")

    })

    observe({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "Samp_col2",
        choices=names(df3.meta),
        selected = "Sample_Name")
    })


    output$Tb_TCR_clonotypes.Umap <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      calls <- select_group_metadata()
      validate(
        need(nrow(calls)>0,
             error_message_val_UMAP)
      )
      calls
      UMAP.wt.clonality2 <- For_col_top()
      colorblind_vector <-as.data.frame(unlist(colors_UMAP_Topclonotypes()))

      if (dim(colorblind_vector)[1]==0) {
        num <- as.data.frame(unique(UMAP.wt.clonality2$topclones))

        if (input$colourtype == "default") {
          colorblind_vector <- c(gg_fill_hue(dim(num)[1]))
        } else if (input$colourtype == "hcl.colors") {
          colorblind_vector <- c(hcl.colors(dim(num)[1], palette = "viridis"))
        } else if (input$colourtype == "topo.colors") {
          colorblind_vector <- c(topo.colors(dim(num)[1]))
        } else if (input$colourtype == "heat.colors") {
          colorblind_vector <- c(heat.colors(dim(num)[1]))
        } else if (input$colourtype == "terrain.colors") {
          colorblind_vector <- c(terrain.colors(dim(num)[1]))
        } else if (input$colourtype == "rainbow") {
          colorblind_vector <- c(rainbow(dim(num)[1]))
        } else if (input$colourtype == "random") {
          colorblind_vector <- distinctColorPalette(dim(num)[1])

        }  else {

        }

      }
      colorblind_vector <- as.data.frame(colorblind_vector)


      topclones_col <- as.data.frame(unique(UMAP.wt.clonality2$topclones))
      names(topclones_col) <- "topclones"
      topclones_col$col <- colorblind_vector
      topclones_col

    })

    output$Tb_For_colouring_check <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      df2 <- as.data.frame(select_group_metadata())
      validate(
        need(nrow(df2)>0,
             error_message_val_UMAP)
      )
      names(df2) <- "V1"
      df2 <- as.data.frame(df2[order(df2$V1),])
      names(df2) <- "V1"
      df2
    })

    ## TCR expansion ----
    TCR_Expanded <- reactive({
      sc <- UMAP_metadata_with_labs()

      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )

      req(input$V_gene_sc, input$Samp_col)

      df3.meta <- sc@meta.data
      df3.meta2 <- df3.meta[,names(df3.meta) %in% c(input$Samp_col,input$V_gene_sc)]
      names(df3.meta2)[names(df3.meta2) %in% input$Samp_col] <- "ID_Column"
      names(df3.meta2)[names(df3.meta2) %in% input$V_gene_sc] <- "v_gene_selected"
      # }
      df3.meta3 <- df3.meta2
      df3.meta3$v_gene_selected <- ifelse(df3.meta3$v_gene_selected=="_._","Unknown",df3.meta3$v_gene_selected)
      df3.meta3$v_gene_selected <- ifelse(df3.meta3$v_gene_selected=="NA_NA & NA_NA","Unknown",df3.meta3$v_gene_selected)
      df3.meta3$v_gene_selected <- ifelse(df3.meta3$v_gene_selected=="","Unknown",df3.meta3$v_gene_selected)
      df3.meta3$v_gene_selected <- ifelse(df3.meta3$v_gene_selected=="NA","Unknown",df3.meta3$v_gene_selected)
      df3.meta3$v_gene_selected[is.na(df3.meta3$v_gene_selected)] <- "Unknown"
      df3.meta3
      if (nrow(df3.meta3[-c(grep("Unknown",df3.meta3$v_gene_selected )),]>0)) {
        df3.meta3 <- df3.meta3[-c(grep("Unknown",df3.meta3$v_gene_selected )),]
      }
      df3.meta3
      meta2.names <- names(df3.meta3)
      df3.meta3$samp.count <- 1
      total.condition <- as.data.frame(ddply(df3.meta3,"ID_Column",numcolwise(sum)))
      total.condition
      emtpy <- matrix(nrow =dim(df3.meta3)[1],ncol=dim(total.condition)[1])

      for (i in 1:dim(df3.meta3)[1]) {

        emtpy[i,] <- ifelse(df3.meta3$ID_Column[i]==total.condition$ID_Column[1:dim(total.condition)[1]],
                            total.condition[total.condition$ID_Column==total.condition$ID_Column[1:dim(total.condition)[1]],2],F)
      }
      as.data.frame(emtpy)
      #
      df3.meta3$frequency <- 1/rowSums(emtpy)
      df3.meta3$percent <- 1/rowSums(emtpy)*100

      df3 <- as.data.frame(ddply(df3.meta3,meta2.names,numcolwise(sum)))
      df3 <- df3[order(df3$samp.count,decreasing = T),]

      df4 <- df3 %>%
        mutate(Frequency_expanded = case_when(
          frequency <=1e-4 ~ "1. Rare (0 > X < 1e-4)",
          frequency <= 0.001 ~ "2. Small (1e-4 > X <= 0.001)",
          frequency <= 0.01 ~ "3. Medium (0.001 > X <= 0.01)",
          frequency <= 0.10 ~ "4. Large (0.01 > X <= 0.1)",
          frequency <= 0.50 ~ "5. Gigantic (0.1 > X <= 0.5)",
          frequency <= 1 ~ "6. Hyperexpanded (0.5 > X <= 1)",
          TRUE ~ "Other"))

      df4 <- df4 %>%
        mutate(Number_expanded = case_when(
          samp.count <=1 ~ "1. Single (0 < X <= 1)",
          samp.count <=5 ~ "2. Small (1 < X <= 5)",
          samp.count <=20 ~ "3. Medium (5 < X <= 20)",
          samp.count <=100 ~ "4. Large (20 < X <= 100)",
          samp.count <= 500 ~ "5. Hyperexpanded (100 < X <= 500)",
          TRUE ~ "6. Hyperexpanded (>500)"))
      df4
    })

    output$Tb_TCR_clonotypes.table <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      umap.meta <- sc@meta.data
      names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
      names(umap.meta)[names(umap.meta) %in% input$V_gene_sc] <- "v_gene_selected"
      umap.meta
      FILE_MERGED <- merge(umap.meta,TCR_Expanded(),by=c("v_gene_selected","ID_Column"),all.x=T)
      FILE_MERGED$ID_Column[FILE_MERGED$v_gene_selected == "NA"] <- "unknown"
      subset(FILE_MERGED,FILE_MERGED$v_gene_selected != "unknown")
    })

    ## Umap -----
    create_UMAP2 <- reactive({
      sc <- sc()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      umap.meta <- sc@meta.data
      names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
      names(umap.meta)[names(umap.meta) %in% input$V_gene_sc] <- "v_gene_selected"
      umap.meta


      sc <- merge(umap.meta,TCR_Expanded(),by=c("v_gene_selected","ID_Column"),all.x=T)
      ggplot(sc,aes(x=UMAP_1,UMAP_2,colour=seurat_clusters))+
        geom_point()+
        scale_color_manual(values = rainbow(length(unique(sc$seurat_clusters))) , na.value=input$NA_col_analysis)+
        theme_bw() +
        theme(
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=0),
          axis.title.x = element_text(colour="black",angle=0,vjust=.5,face="plain",family=input$font_type,size = input$title.text.sizer2),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )
    })

    output$create_UMAP_sc2 <- renderPlot({
      create_UMAP2()
    })

    output$downloadPlot_UMAP2 <- downloadHandler(
      filename = function() {
        x <- today()
        paste("TCR_Explore_UMAP_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_UMAP2,height=input$height_UMAP2, onefile = FALSE) # open the pdf device
        plot(create_UMAP2())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_UMAP2 <- downloadHandler(
      filename = function() {
        x <- today()
        paste("TCR_Explore_UMAP_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_UMAP2,
            height = input$height_png_UMAP2,
            res = input$resolution_PNG_UMAP2)
        plot(create_UMAP2())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    # clonotypes Plot -----
    cols_clonal_plot <- reactive({
      df4 <- TCR_Expanded()
      df4 <- df4[order(df4[,names(df4) %in% input$Graph_type_bar]),]
      df_col <- unique(df4[,names(df4) %in% input$Graph_type_bar])

      num <- as.data.frame(unique(df4[,names(df4) %in% input$Graph_type_bar]))
      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep="_"), paste(num[i,]), col.gg[i])
        })
      }
      else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep="_"), paste(num[i,]), col.hcl[i])
        })
      }
      else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep="_"), paste(num[i,]), col.topo[i])
        })
      }
      else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep="_"), paste(num[i,]), heat_col[i])
        })
      }
      else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep="_"), paste(num[i,]), col.terrain[i])
        })
      }

      else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep="_"), paste(num[i,]), palette_rainbow[i])
        }) }
      else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep="_"), paste(num[i,]), palette1[i])
        })
      }
      else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep="_"), paste(num[i,]), input$one.colour.default)
        })


      } # one colour

    })

    output$myPanel_clonal_plot <- renderUI({cols_clonal_plot()})

    colors_clonal_plot <- reactive({
      df4 <- TCR_Expanded()
      df4 <- df4[order(df4[,names(df4) %in% input$Graph_type_bar]),]
      df_col <- unique(df4[,names(df4) %in% input$Graph_type_bar])

      num <- as.data.frame(unique(df4[,names(df4) %in% input$Graph_type_bar]))
      # num <- unique(top_BD_cluster$Selected_function)

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.clonotype", i, sep="_")]]
      })
    })
    clonal_plot <- reactive({

      df4 <- TCR_Expanded()
      names(df4)[names(df4) %in% input$Samp_col] <- "ID_Column"
      df4 <- df4[df4$ID_Column %in% input$ID_Column_factor,]
      df4$ID_Column <- as.character(df4$ID_Column)
      df4$ID_Column <- factor(df4$ID_Column,levels = input$ID_Column_factor)
      df.col.1 <- unlist(colors_clonal_plot())

      # ggplot(df4, aes(y = frequency, x = Sample_Name, fill=Clonality, label = Clonality)) +
      #   # ggplot(df4,aes(x=Sample_Name,y=frequency,fill=Clonality))+
      #   geom_bar(stat="identity")+
      #   theme_bw() +
      #   scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 20),values=alpha(heat.colors(5), 1))

      ggplot(df4,aes(x=ID_Column,y=frequency,fill=get(input$Graph_type_bar),colour= get(input$Graph_type_bar),label=get(input$Graph_type_bar)))+
        geom_bar(stat="identity")+
        theme_bw() +
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 50),values=alpha(df.col.1, 1),na.value = input$NA_col_analysis) +
        scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 50),values = alpha(df.col.1, 1)) +
        theme(
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=90),
          axis.title.x = element_blank(),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )

    })

    # , na.value=input$NA_col_analysis
    ### top clonotypes -----
    vals <- reactiveValues(top10=NULL)
    cols_top_clonal_plot <- reactive({
      df4 <- TCR_Expanded()
      top10 <- df4 %>%
        group_by(ID_Column) %>%
        top_n(n = input$top_no_clonotypes, wt = frequency)

      top10 <- top10[order(top10$frequency,decreasing = F),]
      unique.top <- unique(top10$v_gene_selected)
      top10$v_gene_selected <- factor(top10$v_gene_selected,
                                      levels = unique.top)

      num <- as.data.frame(unique(top10$v_gene_selected))
      # num <- unique(top_BD_cluster$Selected_function)

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")


      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep="_"), paste(num[i,]), col.gg[i])
        })
      }
      else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep="_"), paste(num[i,]), col.hcl[i])
        })
      }
      else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep="_"), paste(num[i,]), col.topo[i])
        })
      }
      else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep="_"), paste(num[i,]), heat_col[i])
        })
      }
      else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep="_"), paste(num[i,]), col.terrain[i])
        })
      }

      else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep="_"), paste(num[i,]), palette_rainbow[i])
        }) }
      else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep="_"), paste(num[i,]), palette1[i])
        })

      }

      else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep="_"), paste(num[i,]), input$one.colour.default)
        })


      }

    })
    output$myPanel_top_clonal_plot <- renderUI({cols_top_clonal_plot()})
    colors_top_clonal_plot <- reactive({
      df4 <- TCR_Expanded()
      df4
      top10 <- df4 %>%
        group_by(ID_Column) %>%
        top_n(n = input$top_no_clonotypes, wt = frequency)

      top10 <- top10[order(top10$frequency,decreasing = F),]
      unique.top <- unique(top10$v_gene_selected)
      top10$v_gene_selected <- factor(top10$v_gene_selected,
                                      levels = unique.top)

      num <- as.data.frame(unique(top10$v_gene_selected))

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.clonotype_top", i, sep="_")]]
      })
    })

    top_clonal_plot <- reactive({
      df4 <- TCR_Expanded()
      df4
      top10 <- df4 %>%
        group_by(ID_Column) %>%
        top_n(n = input$top_no_clonotypes, wt = frequency)

      top10 <- top10[order(top10$frequency,decreasing = F),]
      unique.top <- unique(top10$v_gene_selected)
      top10$v_gene_selected <- factor(top10$v_gene_selected,
                                      levels = unique.top)

      df.col.1 <- unlist(colors_top_clonal_plot())
      unique.top$col <- as.list(df.col.1)
      ggplot(top10,aes(x=ID_Column,y=frequency,fill=v_gene_selected,label=v_gene_selected,colour="black"))+
        geom_bar(stat="identity")+
        theme_bw() +
        scale_fill_manual(values=alpha(unique.top$col, 1),na.value = input$NA_col_analysis,labels = ~ stringr::str_wrap(.x, width = 20)) +
        scale_colour_manual(values = "black",guide = "none")+
        theme(
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=90),
          axis.title.x = element_blank(),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )

      # ggplot(top10,aes(x=Sample_Name,y=frequency,fill=v_gene_cdr3_AB_GD, colour="Black"))+
      #   geom_bar()+
      #   scale_fill_manual(values=alpha(df.col.1, 1)) +
      #   scale_colour_manual(values = "black",guide = "none")
    })

    output$clonality.bar.graph <- renderPlot({
      if (input$Graph_type_bar=="Frequency_expanded") {
        clonal_plot()
      }
      else if (input$Graph_type_bar=="Number_expanded")  {
        clonal_plot()
      }
      else {
        top_clonal_plot()
      }

    })

    output$clonality.bar.graph2 <- renderPlot({
      if (input$Graph_type_bar=="Frequency_expanded") {
        clonal_plot()
      }
      else if (input$Graph_type_bar=="Number_expanded")  {
        clonal_plot()
      }
      else {
        top_clonal_plot()
      }

    })
    # Downloading the bar plot -------
    output$downloadPlot_clonaity.bar.graph <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$V_gene_sc,"_clonal_expanded_",x, ".pdf", sep = "")

      },
      content = function(file) {
        pdf(file, width=input$width_clonality.bar.graph,height=input$height_clonality.bar.graph, onefile = FALSE) # open the pdf device
        if (input$Graph_type_bar=="Frequency_expanded") {
          plot(clonal_plot())
        }
        else if (input$Graph_type_bar=="Number_expanded")  {
          plot(clonal_plot())
        }
        else {
          plot(top_clonal_plot())
        }
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_clonaity.bar.graph <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$V_gene_sc,"_clonal_expanded_",x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_clonality.bar.graph,
            height = input$height_png_clonality.bar.graph,
            res = input$resolution_PNG_clonality.bar.graph)
        if (input$Graph_type_bar=="Clonality") {
          plot(clonal_plot())
        }
        else if (input$Graph_type_bar=="Number_expanded")  {
          plot(clonal_plot())
        }
        else {
          plot(top_clonal_plot())
        }
        dev.off()},   contentType = "application/png" # MIME type of the image
    )


    # UMAP clonotype -> TCR -----
    cols_UMAP_clonal_plot <- reactive({
      UMAP.wt.clonality <- UMAP.TCRclonalit()
      UMAP.wt.clonality <-  UMAP.wt.clonality[order( UMAP.wt.clonality$TYPE.clonality),]
      UMAP.wt.clonality$TYPE.clonality <- factor(UMAP.wt.clonality$TYPE.clonality,levels = unique(UMAP.wt.clonality$TYPE.clonality))
      num <- as.data.frame(unique(UMAP.wt.clonality$TYPE.clonality))
      num <- as.data.frame(num[complete.cases(num)==T,])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep="_"), paste(num[i,]), col.gg[i])
        })
      }
      else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep="_"), paste(num[i,]), col.hcl[i])
        })
      }
      else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep="_"), paste(num[i,]), col.topo[i])
        })
      }
      else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep="_"), paste(num[i,]), heat_col[i])
        })
      }
      else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep="_"), paste(num[i,]), col.terrain[i])
        })
      }

      else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep="_"), paste(num[i,]), palette_rainbow[i])
        }) }
      else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep="_"), paste(num[i,]), palette1[i])
        })

      }
      else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep="_"), paste(num[i,]), input$one.colour.default)
        })


      } # one colour

    })
    output$cols_UMAP_clonal_plot <- renderUI({cols_UMAP_clonal_plot()})
    colors_UMAP_clonal_plot <- reactive({
      UMAP.wt.clonality <- UMAP.TCRclonalit()
      UMAP.wt.clonality <-  UMAP.wt.clonality[order( UMAP.wt.clonality$TYPE.clonality),]
      UMAP.wt.clonality$TYPE.clonality <- factor(UMAP.wt.clonality$TYPE.clonality,levels = unique(UMAP.wt.clonality$TYPE.clonality))
      # UMAP.wt.clonality <- UMAP.wt.clonality[order(UMAP.wt.clonality$TYPE.clonality),]
      num <- as.data.frame(unique(UMAP.wt.clonality$TYPE.clonality))
      num <- as.data.frame(num[complete.cases(num)==T,])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.UMAP_clonotype", i, sep="_")]]
      })
    })

    observe({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )

      df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "RDS_V_gene_A",
        choices=names(df3.meta),
        selected = "Av_gene")

      updateSelectInput(
        session,
        "RDS_V_gene_B",
        choices=names(df3.meta),
        selected = "Bv_gene")

      updateSelectInput(
        session,
        "RDS_cdr3_A",
        choices=names(df3.meta),
        selected = "Acdr3")
      updateSelectInput(
        session,
        "RDS_cdr3_B",
        choices=names(df3.meta),
        selected = "Bcdr3")

    })

    UMAP.TCRclonalit <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      umap.meta <- sc@meta.data
      names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
      names(umap.meta)[names(umap.meta) %in% input$V_gene_sc] <- "v_gene_selected"
      umap.meta
      UMAP.wt.clonality <- merge(umap.meta,TCR_Expanded(),by=c("v_gene_selected","ID_Column"),all.x=T)
      if (input$Graph_type_bar== "Number_expanded") {
        UMAP.wt.clonality$TYPE.clonality <- paste(UMAP.wt.clonality$Number_expanded)
      }
      else {
        UMAP.wt.clonality$TYPE.clonality <- paste(UMAP.wt.clonality$Frequency_expanded)

      }
      # UMAP.wt.clonality$TYPE.clonality<- ifelse(grepl("NA", UMAP.wt.clonality$TYPE.clonality),NA,UMAP.wt.clonality$TYPE.clonality)
      UMAP.wt.clonality <- subset(UMAP.wt.clonality,UMAP.wt.clonality$TYPE.clonality != "NA")
      UMAP.wt.clonality

    })

    UMAP.TCRclonalit2 <- reactive({
      UMAP.wt.clonality <- UMAP.TCRclonalit()
      UMAP.wt.clonality <-  UMAP.wt.clonality[order(UMAP.wt.clonality$TYPE.clonality),]

      colorblind_vector <-as.data.frame(unlist(colors_UMAP_clonal_plot()))

      if (dim(colorblind_vector)[1]==0) {
        num <- as.data.frame(unique(UMAP.wt.clonality$TYPE.clonality))

        if (input$colourtype == "default") {
          colorblind_vector <- c(gg_fill_hue(dim(num)[1]))
        } else if (input$colourtype == "hcl.colors") {
          colorblind_vector <- c(hcl.colors(dim(num)[1], palette = "viridis"))
        } else if (input$colourtype == "topo.colors") {
          colorblind_vector <- c(topo.colors(dim(num)[1]))
        } else if (input$colourtype == "heat.colors") {
          colorblind_vector <- c(heat.colors(dim(num)[1]))
        } else if (input$colourtype == "terrain.colors") {
          colorblind_vector <- c(terrain.colors(dim(num)[1]))
        } else if (input$colourtype == "rainbow") {
          colorblind_vector <- c(rainbow(dim(num)[1]))
        } else if (input$colourtype == "random") {
          colorblind_vector <- distinctColorPalette(dim(num)[1])

        }  else {

        }

      }
      colorblind_vector <- as.data.frame(colorblind_vector)
      names(colorblind_vector) <- "cols"

      names(UMAP.wt.clonality)[names(UMAP.wt.clonality) %in% input$Samp_col] <- "ID_Column"

      if (input$Split_by_group=="yes") {
        UMAP.wt.clonality <- UMAP.wt.clonality[UMAP.wt.clonality$ID_Column %in% input$ID_Column_factor,]
      }

      UMAP.wt.clonality$ID_Column <- factor(UMAP.wt.clonality$ID_Column,levels = input$ID_Column_factor)

      if (input$filter_umap_expand == 'yes') {
        UMAP.wt.clonality <- subset(UMAP.wt.clonality,UMAP.wt.clonality$UMAP_1 < input$UMAP_1x & UMAP.wt.clonality$UMAP_1 > input$UMAP_1y)
        UMAP.wt.clonality <- subset(UMAP.wt.clonality,UMAP.wt.clonality$UMAP_2 < input$UMAP_2x & UMAP.wt.clonality$UMAP_2 > input$UMAP_2y)
      }

      plot <- ggplot(UMAP.wt.clonality,aes(x=UMAP_1,UMAP_2,colour=TYPE.clonality,alpha = TYPE.clonality,label=TYPE.clonality))+
        geom_point()+
        scale_color_manual(values = colorblind_vector$cols, na.value=input$NA_col_analysis,labels = ~ stringr::str_wrap(.x, width = 40))+
        scale_alpha_manual(values = rep(1,length(unique(UMAP.wt.clonality$TYPE.clonality))), na.value=0.5,labels = ~ stringr::str_wrap(.x, width = 40))+
        theme_bw() +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=0),
          axis.title.x = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )
      if (input$Split_by_group=="no") {
        plot
      }
      else {
        plot+ facet_wrap(~ID_Column, nrow = input$wrap_row)
      }

    })
    output$clonality.TCR.UMAP <- renderPlot({
      UMAP.TCRclonalit2()
    })

    output$downloadPlot_TCR.UMAP <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$V_gene_sc,"_TCR.UMAP_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_TCR.UMAP,height=input$height_TCR.UMAP, onefile = FALSE) # open the pdf device

        plot(UMAP.TCRclonalit2())
        dev.off()
      }, contentType = "application/pdf" )

    output$downloadPlotPNG_TCR.UMAP <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$V_gene_sc,"_TCR.UMAP_",x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_TCR.UMAP,height = input$height_png_TCR.UMAP,res = input$resolution_PNG_TCR.UMAP)
        plot(UMAP.TCRclonalit2())
        dev.off()
      },   contentType = "application/png" # MIME type of the image
    )
    # Freq ------
    observe({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "Samp_col_UMAP",
        choices=names(df3.meta),
        selected = "ID_Column")
    }) # ID_Column

    For_col_top <- reactive({
      df4 <- TCR_Expanded()
      df4 <- subset(df4, df4$samp.count!=1)
      top10 <- df4 %>%
        group_by(ID_Column) %>%
        top_n(n = input$top_no_clonotypes, wt = frequency)

      top10 <- top10[order(top10$frequency,decreasing = F),]
      unique.top <- unique(top10$v_gene_selected)
      top10$v_gene_selected <- factor(top10$v_gene_selected,
                                      levels = unique.top)

      top10.2 <- as.data.frame(top10$v_gene_selected)
      names(top10.2) <- "v_gene_selected"
      top10.2$ID_Column <- top10$ID_Column
      top10.2$topclones <- top10.2$v_gene_selected

      UMAP.wt.clonality <- UMAP.TCRclonalit()
      UMAP.wt.clonality2 <- merge(UMAP.wt.clonality,top10.2,by=c("v_gene_selected","ID_Column"))
      UMAP.wt.clonality2 <- UMAP.wt.clonality2[order(UMAP.wt.clonality2$Cell_Index),]
      UMAP.wt.clonality2 %>% mutate(across(where(is.factor), as.character)) -> UMAP.wt.clonality2
      UMAP.wt.clonality2$topclones[is.na(UMAP.wt.clonality2$topclones)]<- "not selected"
      UMAP.wt.clonality2$topclones <- factor(UMAP.wt.clonality2$topclones,levels = c("not selected",unique(as.character(top10.2$v_gene_selected))))
      UMAP.wt.clonality2 <- UMAP.wt.clonality2[order(UMAP.wt.clonality2$topclones),]

    })

    cols_UMAP_Topclonotypes <- reactive({
      UMAP.wt.clonality2 <- For_col_top()
      num <- as.data.frame(unique(UMAP.wt.clonality2$topclones))
      col.gg <- c(gg_fill_hue(dim(num)[1]))
      palette_rainbow <- c(rainbow(dim(num)[1]))
      heat_col <- c(heat.colors(dim(num)[1]))
      col.terrain <- c(terrain.colors(dim(num)[1]))
      col.topo <- c(topo.colors(dim(num)[1]))
      col.hcl <- c(hcl.colors(dim(num)[1], palette = "viridis"))


      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype_top", i, sep="_"), paste(num[i,]), col.gg[i])
        })
      }
      else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype_top", i, sep="_"), paste(num[i,]), col.hcl[i])
        })
      }
      else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype_top", i, sep="_"), paste(num[i,]), col.topo[i])
        })
      }
      else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype_top", i, sep="_"), paste(num[i,]), heat_col[i])
        })
      }
      else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype_top", i, sep="_"), paste(num[i,]), col.terrain[i])
        })
      }

      else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype_top", i, sep="_"), paste(num[i,]), palette_rainbow[i])
        }) }
      else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype_top", i, sep="_"), paste(num[i,]), palette1[i])
        })

      }
      else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep="_"), paste(num[i,]), input$one.colour.default)
        })


      } # one colour

    })

    output$cols_UMAP_Topclonotypes <- renderUI({cols_UMAP_Topclonotypes()})

    colors_UMAP_Topclonotypes <- reactive({
      UMAP.wt.clonality2 <- For_col_top()

      num <- as.data.frame(unique(UMAP.wt.clonality2$topclones))

      col.gg <- c(gg_fill_hue(dim(num)[1]))
      palette_rainbow <- c(rainbow(dim(num)[1]))
      heat_col <- c(heat.colors(dim(num)[1]))
      col.terrain <- c(terrain.colors(dim(num)[1]))
      col.topo <- c(topo.colors(dim(num)[1]))
      col.hcl <- c(hcl.colors(dim(num)[1], palette = "viridis"))

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.UMAP_clonotype_top", i, sep="_")]]
      })
    })

    select_group_metadata <- reactive ({
      sc <- UMAP_metadata_with_labs()

      validate(
        need(nrow(sc)>0,
             "upload file")
      )
      df <- sc@meta.data
      df2 <- as.data.frame(unique(df[names(df) %in% input$Samp_col]))
      df2 <- as.data.frame(df2)
      df2

    })

    observe({

      df2 <- select_group_metadata()

      validate(
        need(nrow(df2)>0,
             error_message_val1)
      )

      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- as.data.frame(df2[order(df2$V1),])
      names(df2) <- "V1"
      df2

      updateSelectInput(
        session,
        "ID_Column_metadata",
        choices=df2$V1,
        selected = df2$V1[1]
      )


    })
    observe({
      df2 <- select_group_metadata()
      validate(
        need(nrow(df2)>0,
             error_message_val1)
      )
      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- as.data.frame(df2[order(df2$V1),])
      names(df2) <- "V1"
      df2
      # df2 <- subset(df2,df2$V1 != "NA")

      df3 <- subset(df2,df2$V1 != "NA")

      updateSelectInput(
        session,
        "ID_Column_factor",
        choices=df2$V1,
        selected = df3$V1
      )
    })
    observe({
      df2 <- select_group_metadata()
      validate(
        need(nrow(df2)>0,
             error_message_val1)
      )
      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- as.data.frame(df2[order(df2$V1),])
      names(df2) <- "V1"
      df2
      df2 <- subset(df2,df2$V1 != "NA")

      updateSelectInput(
        session,
        "selected_Indiv",
        choices=df2$V1,
        selected = df2$V1[1]
      )
    })

    Topclonotypes <- reactive({
      UMAP.wt.clonality2 <- For_col_top()
      UMAP.wt.clonality2 <- UMAP.wt.clonality2[UMAP.wt.clonality2$ID_Column %in% input$ID_Column_factor,]
      UMAP.wt.clonality2$ID_Column <- factor(UMAP.wt.clonality2$ID_Column, levels = input$ID_Column_factor)

      colorblind_vector <-as.data.frame(unlist(colors_UMAP_Topclonotypes()))

      if (dim(colorblind_vector)[1]==0) {
        num <- as.data.frame(unique(UMAP.wt.clonality2$topclones))

        if (input$colourtype == "default") {
          colorblind_vector <- c(gg_fill_hue(dim(num)[1]))
        } else if (input$colourtype == "hcl.colors") {
          colorblind_vector <- c(hcl.colors(dim(num)[1], palette = "viridis"))
        } else if (input$colourtype == "topo.colors") {
          colorblind_vector <- c(topo.colors(dim(num)[1]))
        } else if (input$colourtype == "heat.colors") {
          colorblind_vector <- c(heat.colors(dim(num)[1]))
        } else if (input$colourtype == "terrain.colors") {
          colorblind_vector <- c(terrain.colors(dim(num)[1]))
        } else if (input$colourtype == "rainbow") {
          colorblind_vector <- c(rainbow(dim(num)[1]))
        } else if (input$colourtype == "random") {
          colorblind_vector <- distinctColorPalette(dim(num)[1])

        }  else {

        }

      }
      colorblind_vector <- as.data.frame(colorblind_vector)


      names(colorblind_vector) <- "col"
      topclones_col <- as.data.frame(unique(UMAP.wt.clonality2$topclones))
      names(topclones_col) <- "topclones"
      topclones_col$col <- colorblind_vector$col
      topclones_col



      if (input$display_all_samps=="yes" & input$Split_by_group=="no") {
        topclones_col
      } else if (input$display_all_samps=="yes" & input$Split_by_group=="yes") {
        topclones_col
      } else if (input$display_all_samps=="no" & input$Split_by_group=="yes") {
        UMAP.wt.clonality2 <- UMAP.wt.clonality2[UMAP.wt.clonality2$ID_Column %in% input$ID_Column_metadata,]
        topclones_col <- topclones_col[topclones_col$topclones %in% unique(UMAP.wt.clonality2$topclones),]
        topclones_col
      } else {
        UMAP.wt.clonality2 <- UMAP.wt.clonality2[UMAP.wt.clonality2$ID_Column %in% input$ID_Column_metadata,]
        topclones_col <- topclones_col[topclones_col$topclones %in% unique(UMAP.wt.clonality2$topclones),]
        topclones_col
      }

      UMAP.wt.clonality2$topclones2 <- gsub(" & "," ",UMAP.wt.clonality2$topclones)

      names(UMAP.wt.clonality2)[names(UMAP.wt.clonality2) %in% input$Samp_col] <- "ID_Column"
      UMAP.wt.clonality2 <- UMAP.wt.clonality2[UMAP.wt.clonality2$ID_Column %in% input$ID_Column_factor,]
      UMAP.wt.clonality2$ID_Column <- factor(UMAP.wt.clonality2$ID_Column,levels = input$ID_Column_factor)



      plot <- ggplot(data=UMAP.wt.clonality2,aes(x=UMAP_1,UMAP_2,colour=topclones2,label =topclones2 ))+
        geom_point()+
        scale_color_manual(values = topclones_col$col, na.value=input$NA_col_analysis,labels = ~ stringr::str_wrap(.x, width = 50))+
        # scale_color_manual(values = colorblind_vector)+
        scale_alpha_manual(values = 1, na.value = 0.1,labels = ~ stringr::str_wrap(.x, width = 50))+
        theme_bw() +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=0),
          axis.title.x = element_text(colour="black",angle=0,vjust=.5,face="plain",family=input$font_type,size = input$title.text.sizer2),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )

      if (input$Split_by_group=="no") {
        plot
      }
      else {
        # UMAP.wt.clonality2$ID.names <- UMAP.wt.clonality2[,names(UMAP.wt.clonality2) %in% input$meta_data_sc_]
        # UMAP.wt.clonality2$ID.names <- UMAP.wt.clonality2[,names(UMAP.wt.clonality2) %in% input$Samp_col_UMAP]

        plot+ facet_wrap(~ID_Column, nrow = input$wrap_row)
      }
    })

    output$clonality.TCR.UMAP.top <- renderPlot({
      Topclonotypes()
    })

    output$downloadPlot_TCR.UMAP_top <- downloadHandler(
      filename = function() {
        x <- today()
        paste("TCR_Explore_TCR.UMAP_top_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_TCR.UMAP_top,height=input$height_TCR.UMAP_top, onefile = FALSE) # open the pdf device

        plot(Topclonotypes())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_TCR.UMAP_top <- downloadHandler(
      filename = function() {
        x <- today()
        paste("TCR_Explore_TCR.UMAP_top_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_TCR.UMAP_top,height = input$height_png_TCR.UMAP_top,res = input$resolution_PNG_TCR.UMAP_top)
        plot(Topclonotypes())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )


    ### check Idents ----

    vals_meta.sc2 <- reactiveValues(AddIndents_SCobj=NULL)
    vals_meta.sc3 <- reactiveValues(AddIndents_SCobj2=NULL)

    observeEvent(input$run_update_clust,{
      sc <- df_sc_clust()
      df.to.order <- as.data.frame(unique(Idents(df_sc_clust())))
      names(df.to.order) <- "V1"
      df.to.order2 <- as.data.frame(do.call(rbind, strsplit(as.character(df.to.order$V1), " ")))
      df.to.order2$uniqueID <- df.to.order$V1

      if (input$multiple_group_sc =="yes") {
        df.to.order2 <- df.to.order2 %>%
          arrange(V1,V2)
      }
      else {
        df.to.order2 <- df.to.order2 %>%
          arrange(V1)
      }
      vals_meta.sc3$AddIndents_SCobj2 <- df.to.order2
      vals_meta.sc3$AddIndents_SCobj2
    })
    observe({
      unique.Idents <- vals_meta.sc3$AddIndents_SCobj2
      updateSelectInput(
        session,
        "unique.Idents1",
        choices=unique.Idents$uniqueID,
        selected = unique.Idents$uniqueID[1])
    }) # ident.1
    observe({
      unique.Idents <- vals_meta.sc3$AddIndents_SCobj2
      updateSelectInput(
        session,
        "unique.Idents2",
        choices=unique.Idents$uniqueID,
        selected = unique.Idents$uniqueID[2])
    }) # ident.2

    DEx_sc <- reactive({
      b.interferon.response <- FindMarkers(df_sc_clust(),
                                           ident.1 = input$unique.Idents1,
                                           ident.2 = input$unique.Idents2,
                                           logfc.threshold = 0.25, # Give me ALL results
                                           min.pct = 0.25)
      b.interferon.response
    })

    output$DEx_table_comparison <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      DEx_sc()
    })

    output$downloaddf_DEx_sc <- downloadHandler(
      filename = function(){
        paste(" Differential expresion of ",input$unique.Idents1, " vs. ", input$unique.Idents2,gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(DEx_sc())
        write.csv(df,file, row.names = T)
      } )


    output$downloaddf_DEx_sc_ggVolcanoR <- downloadHandler(
      filename = function(){
        paste(" Differential expresion of ",input$unique.Idents1, " vs. ", input$unique.Idents2,gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(DEx_sc())
        df1 <- df
        df1$ID <- rownames(df1)
        df1$Pvalue <- df1$p_val_adj
        df1$logFC <- df1$avg_log2FC
        df2 <- df1[,names(df1) %in% c("ID","logFC","Pvalue")]
        write.csv(df2,file, row.names = F)
      } )

    ggplot_DEx_sc <- reactive({

      dat <- DEx_sc()
      dat <- as.data.frame(dat)
      dat <- dat[order(dat$p_val_adj),]
      dat$ID <- rownames(dat)
      neg <- -1*0.58
      pos <- 0.58

      maximum <- input$max

      mutateddf <- mutate(dat, sig=ifelse(dat$p_val_adj<0.05, "Pvalue<0.05", "Not Sig"))
      sig <- subset(dat, dat$p_val_adj<0.05 & abs(dat$avg_log2FC)>0.58)
      top <- sig[(1:30),]

      gene_list <- top$ID

      mutateddf.gene <- mutate(mutateddf, top=ifelse(mutateddf$ID %in% gene_list, "top", "other"))
      mutateddf.gene

      # no labels -----
      sub.mutateddf.gene <- mutate(mutateddf.gene,
                                   colour=ifelse(mutateddf.gene$p_val_adj<0.05 & mutateddf.gene$avg_log2FC>pos,"sig_up",
                                                 ifelse(mutateddf.gene$p_val_adj<0.05& mutateddf.gene$avg_log2FC<neg,"sig_down","NS")),
                                   alpha=ifelse(mutateddf.gene$p_val_adj<0.05& mutateddf.gene$avg_log2FC>pos,0.25,
                                                ifelse(mutateddf.gene$p_val_adj<0.05& mutateddf.gene$avg_log2FC<neg,0.25,0.1)),
                                   shape=ifelse(mutateddf.gene$p_val_adj<0.05& mutateddf.gene$avg_log2FC>pos,19,
                                                ifelse(mutateddf.gene$p_val_adj<0.05& mutateddf.gene$avg_log2FC<neg,19,1)),
                                   size=ifelse(mutateddf.gene$p_val_adj<0.05& mutateddf.gene$avg_log2FC>pos,3,
                                               ifelse(mutateddf.gene$p_val_adj<0.05& mutateddf.gene$avg_log2FC<neg,3,1)))
      # range of genes -----
      sub.mutateddf.gene2 <- mutate(mutateddf.gene,
                                    colour=ifelse(mutateddf.gene$ID %in% gene_list & mutateddf.gene$avg_log2FC>pos & mutateddf.gene$p_val_adj<0.05, "Labelled_up",
                                                  ifelse(mutateddf.gene$ID %in% gene_list & mutateddf.gene$avg_log2FC<neg & mutateddf.gene$p_val_adj<0.05, "Labelled_down",                                                                                           ifelse(mutateddf.gene$p_val_adj<0.05& mutateddf.gene$avg_log2FC>pos,"Significant-up",
                                                                                                                                                                                                                                                                             ifelse(mutateddf.gene$p_val_adj<0.05& mutateddf.gene$avg_log2FC<neg,"Significant-down","Non-significant")))),

                                    size=ifelse(mutateddf.gene$ID %in% gene_list & mutateddf.gene$avg_log2FC>pos & mutateddf.gene$p_val_adj<0.05, 3,
                                                ifelse(mutateddf.gene$ID %in% gene_list & mutateddf.gene$avg_log2FC<neg & mutateddf.gene$p_val_adj<0.05, 3,                                                                                           ifelse(mutateddf.gene$p_val_adj<0.05& mutateddf.gene$avg_log2FC>pos,3,
                                                                                                                                                                                                                                                             ifelse(mutateddf.gene$p_val_adj<0.05& mutateddf.gene$avg_log2FC<neg,3,1)))),

                                    shape=ifelse(mutateddf.gene$ID %in% gene_list & mutateddf.gene$avg_log2FC>pos & mutateddf.gene$p_val_adj<0.05, 19,
                                                 ifelse(mutateddf.gene$ID %in% gene_list & mutateddf.gene$avg_log2FC<neg & mutateddf.gene$p_val_adj<0.05, 19,                                                                                           ifelse(mutateddf.gene$p_val_adj<0.05& mutateddf.gene$avg_log2FC>pos,15,
                                                                                                                                                                                                                                                               ifelse(mutateddf.gene$p_val_adj<0.05& mutateddf.gene$avg_log2FC<neg,15,1)))),

      )

      y_lable1 <- bquote("-"~Log[10]~(.("p_val_adj")))
      y_lable1

      colour_class3 <- c("Significant-down","Significant-up","Labelled_down","Labelled_up","labelled-Non-significant","Non-significant")
      colour.df3 <- as.data.frame(c("Significant-down","Significant-up","Labelled_down","Labelled_up","labelled-Non-significant","Non-significant"))
      names(colour.df3) <- "label"
      # colour.df3$V1 <- c(input$down,input$up,input$col_lab1,input$col_lab2,input$col_lab3,input$NS)
      # colour.df3$shape <- c(input$shape2,input$shape1.1,input$shape1,input$shape1,input$shape1,input$shape3)
      # colour.df3$size <- c(input$size2,input$size1.1,input$size1,input$size1,input$size1,input$size3)
      # colour.df3$alpha <- c(input$alpha2,input$alpha2,input$alpha1,input$alpha1,input$alpha1,input$alpha3)
      # colour.class4 <- colour.df3[colour.df3$label %in% unique(sub.mutateddf.gene2$colour),]
      # sub.mutateddf.gene2$colour <- factor(sub.mutateddf.gene2$colour, levels = colour.class4$label)

      ggplot() +
        geom_point(aes(x=sub.mutateddf.gene2$avg_log2FC, y=-log10(sub.mutateddf.gene2$p_val_adj),
                       col=sub.mutateddf.gene2$colour,
                       # shape=sub.mutateddf.gene2$shape,
                       # alpha=sub.mutateddf.gene2$colour,
                       size=sub.mutateddf.gene2$size
        ),
        ) +
        # scale_color_manual(name="legend",values=colour.class4$V1, labels = colour.class4$label) +
        # scale_shape_manual(name="legend",values=colour.class4$shape, labels=colour.class4$label)+
        # scale_size_manual(name="legend",values=colour.class4$size, labels=colour.class4$label)+
        # scale_alpha_manual(name="legend",values=colour.class4$alpha, labels=colour.class4$label) +

        geom_text_repel(data=sub.mutateddf.gene2[sub.mutateddf.gene2$ID %in% gene_list,]
                        ,aes(x=sub.mutateddf.gene2$avg_log2FC[sub.mutateddf.gene2$ID %in% gene_list],
                             y=  -log10(sub.mutateddf.gene2$p_val_adj)[sub.mutateddf.gene2$ID %in% gene_list],
                             label= sub.mutateddf.gene2$ID[sub.mutateddf.gene2$ID %in% gene_list]),
                        size=6,
                        family='serif',
                        segment.alpha = 0.5,
                        show.legend = F,box.padding = unit(1, 'lines'),
                        max.overlaps = Inf) +
        guides(shape = guide_legend(override.aes = list(size = 5))) +
        theme_bw(base_size = 18)+
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        geom_vline(xintercept=pos, linetype="dashed", color = "darkorange") +
        geom_vline(xintercept=neg, linetype="dashed", color = "darkorange") +
        geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "darkorange") +
        # theme(text=element_text(size=20,family=input$font),
        #       axis.title = element_text(colour="black", size=input$axis,family=input$font),
        #       axis.text.x = element_text(colour="black",size=input$axis_text,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font),
        #       axis.text.y = element_text(colour="black",size=input$axis_text,angle=0,hjust=1,vjust=0,face="plain",family=input$font),
        #       axis.title.x=element_text(colour="black",size=input$axis,angle=0,hjust=.5,vjust=.5,face="plain",family=input$font),
        #       axis.title.y = element_text(colour="black",size=input$axis,angle=90,hjust=.5,vjust=.5,face="plain",family=input$font),
        #       legend.title  =element_blank(),
        #       legend.text = element_text(size=input$legend_size),
        #       legend.position = input$legend_location,
        #       legend.box="vertical",
        #       legend.margin=margin(),
      #       legend.justification = "top")+
      labs(y=y_lable1,
           x=expression(Log[2]~Fold~Change),
           title="")
      # guides(size="none", col = guide_legend(ncol=input$col))+
      # scale_y_continuous(limits = c(0, input$yhigh) ,breaks = seq(0, input$yhigh, by = input$ybreaks))+
      # scale_x_continuous(limits = c(input$xlow, input$xhigh), breaks = seq(input$xlow, input$xhigh, by = input$xbreaks))
    })

    output$volc_plot_cluster <- renderPlot({

      ggplot_DEx_sc()
    })



    # TCR -> UMAP ------
    ## upset plot ----

    ## TCR -> UMAP ------
    Percent_tab_df <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      meta.data <- sc@meta.data
      totals <- meta.data[,names(meta.data) %in% c(input$Samp_col,input$Colour_By_this_overview)]
      names(totals) <- c("groups","Function")
      tab <- table(totals$Function,totals$groups)
      as.matrix(round(t(tab)/colSums(tab)*100,2))

    })

    # output$Percent_tab <- renderPrint(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
    output$Percent_tab <-renderPrint({
      Percent_tab_df()
    })
    output$downloaddf_Percent_tab <- downloadHandler(
      filename = function(){
        x = today()
        paste("GEx_percent_",input$Colour_By_this_overview,"_",x,".tsv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(Percent_tab_df())
        write.table(df,file, row.names = F,sep="\t", quote = F)
      } )

    output$SiteNumInput <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP))

      df_class <- sc@meta.data
      df_class$selected <- df_class[,names(df_class) %in% input$Colour_By_this_overview]
      selectInput("SiteNumInput", "Show on graph", choices = unique(df_class$selected), selected = NULL, multiple = TRUE)
    })

    # Umap classification plot -------
    cols_UMAP_all_classification <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )

      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this_overview]
      top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function),NA,top_BD_cluster$Selected_function)
      top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num)==T,])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep="_"), paste(num[i,]), col.gg[i])
        })
      }
      else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep="_"), paste(num[i,]), col.hcl[i])
        })
      }
      else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep="_"), paste(num[i,]), col.topo[i])
        })
      }
      else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep="_"), paste(num[i,]), heat_col[i])
        })
      }
      else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep="_"), paste(num[i,]), col.terrain[i])
        })
      }

      else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep="_"), paste(num[i,]), palette_rainbow[i])
        }) }
      else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep="_"), paste(num[i,]), palette1[i])
        })

      }
      else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep="_"), paste(num[i,]), input$one.colour.default)
        })


      } # one colour

    })

    output$cols_UMAP_all_classification <- renderUI({cols_UMAP_all_classification()})

    colors_UMAP_all_classification <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this_overview]
      top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function),NA,top_BD_cluster$Selected_function)
      top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num)==T,])

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.UMAP_all_classification", i, sep="_")]]
      })
    })

    UMAP_all_classification <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this_overview]
      top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function),NA,top_BD_cluster$Selected_function)
      # top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
      col.file <- as.data.frame(unique(top_BD_cluster$Selected_function))
      col.file <- as.data.frame(col.file[complete.cases(col.file)==T,])
      names(col.file) <- "V1"
      if (input$show_selected=="Selected_list") {
        col.file$col <- unlist(colors_UMAP_all_classification())
        col.file <- col.file[col.file$V1 %in% input$SiteNumInput,]
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = col.file$V1)
      }

      else {
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))

        col.file$col <- unlist(colors_UMAP_all_classification())
        col.file
      }

      names(top_BD_cluster)[names(top_BD_cluster) %in% input$Samp_col] <- "ID_Column"

      if (input$by_indiv_pie_epi == "yes") {
        top_BD_cluster <- top_BD_cluster[top_BD_cluster$ID_Column %in% input$selected_Indiv,]
      }

      # top_BD_cluster$alpha_val <- ifelse(is.na(top_BD_cluster$Selected_function)==T,0.1,1)


      top_BD_cluster <- top_BD_cluster[top_BD_cluster$ID_Column %in% input$ID_Column_factor,]
      top_BD_cluster$ID_Column <- factor(top_BD_cluster$ID_Column,levels = input$ID_Column_factor)
      md <- top_BD_cluster
      md <- subset(md,md$UMAP_1 > input$Filter_lower_UMAP1_marker_GEX)
      md <- subset(md,md$UMAP_1 < input$Filter_lower_UMAP1_marker2_GEX)

      md <- subset(md,md$UMAP_2 > input$Filter_lower_UMAP2_marker_GEX)
      md <- subset(md,md$UMAP_2 < input$Filter_lower_UMAP2_marker2_GEX)
      top_BD_cluster <- md

      df <- ggplot(top_BD_cluster,aes(x=UMAP_1,UMAP_2,colour=Selected_function,alpha = Selected_function,label = Selected_function))+
        geom_point()+
        scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 20),values = col.file$col, na.value=input$NA_col_analysis)+
        scale_alpha_manual(labels = ~ stringr::str_wrap(.x, width = 20),values = rep(1,length(unique(top_BD_cluster$Selected_function))), na.value=0.1) +
        theme_bw() +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=0),
          axis.title.x = element_text(colour="black",angle=0,vjust=.5,face="plain",family=input$font_type,size = input$title.text.sizer2),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )

      if (input$Split_by_group=="no") {
        df <- df
      }
      else {
        df <- df + facet_wrap(~ID_Column, nrow = input$wrap_row)
      }
      df
    })

    output$UMAP_all_classification2 <- renderPlot({
      UMAP_all_classification()
    })

    output$downloadPlot_UMAP_all_classification  <- downloadHandler(
      filename = function() {
        x <- today()
        paste("UMAP_GEx_",input$Colour_By_this_overview,"_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_UMAP_all_classification ,height=input$height_UMAP_all_classification , onefile = FALSE) # open the pdf device
        plot(UMAP_all_classification())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_UMAP_all_classification  <- downloadHandler(
      filename = function() {
        x <- today()
        paste("UMAP_GEx_",input$Colour_By_this_overview,"_",x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_UMAP_all_classification ,
            height = input$height_png_UMAP_all_classification ,
            res = input$resolution_PNG_UMAP_all_classification )
        plot(UMAP_all_classification())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    #### Pie chart -----

    observe({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      meta.data <- sc@meta.data

      # if(names(meta.data) %in% "T_cells") {
      #   updateSelectInput(
      #     session,
      #     "Split_group_by_",
      #     choices=names(meta.data),
      #     selected = "T_cells")
      # } else {

      updateSelectInput(
        session,
        "Split_group_by_",
        choices=names(meta.data),
        selected = "orig.ident")
      # }


    })
    # pie colouring  ----
    cols_pie <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this_overview]
      top_BD_cluster$Selected_function <- ifelse(top_BD_cluster$Selected_function=="NA","-",top_BD_cluster$Selected_function)
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]

      # top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num)==T,])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep="_"), paste(num[i,]), col.gg[i])
        })
      }
      else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep="_"), paste(num[i,]), col.hcl[i])
        })
      }
      else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep="_"), paste(num[i,]), col.topo[i])
        })
      }
      else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep="_"), paste(num[i,]), heat_col[i])
        })
      }
      else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep="_"), paste(num[i,]), col.terrain[i])
        })
      }

      else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep="_"), paste(num[i,]), palette_rainbow[i])
        }) }
      else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep="_"), paste(num[i,]), palette1[i])
        })

      }
      else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep="_"), paste(num[i,]), input$one.colour.default)
        })


      } # one colour

    })
    output$myPanel_pie <- renderUI({cols_pie()})
    colors_pie <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this_overview]
      top_BD_cluster$Selected_function <- ifelse(top_BD_cluster$Selected_function=="NA","-",top_BD_cluster$Selected_function)
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num)==T,])

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.pie", i, sep="_")]]
      })
    })

    output$table_pie <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this_overview]
      top_BD_cluster$Selected_function <- ifelse(top_BD_cluster$Selected_function=="NA","-",top_BD_cluster$Selected_function)
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
      top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% input$Split_group_by_]
      top_BD_cluster$Selected_group <- factor(top_BD_cluster$Selected_group,level= input$Graph_split_order)
      top_BD_cluster$Selected_indiv <- top_BD_cluster[,names(top_BD_cluster) %in% input$Samp_col]
      top_BD_cluster
      df.col <- unlist(colors_pie())

      if (input$by_indiv_pie_epi == "yes") {
        top_BD_cluster <- top_BD_cluster[top_BD_cluster$Selected_indiv %in% input$selected_Indiv,]
      }

      top_BD_cluster
    })

    Pie_chart_Class <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this_overview]
      top_BD_cluster$Selected_function <- ifelse(top_BD_cluster$Selected_function=="NA","-",top_BD_cluster$Selected_function)
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
      top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% input$Split_group_by_overview]

      df.col <- unlist(colors_pie())
      names(top_BD_cluster)[names(top_BD_cluster) %in% input$Samp_col] <- "ID_Column"

      if (input$by_indiv_pie_epi == "yes") {
        top_BD_cluster <- top_BD_cluster[top_BD_cluster$ID_Column %in% input$selected_Indiv,]
      }

      df3.meta3 <-  as.data.frame(table(top_BD_cluster$Selected_group,top_BD_cluster$Selected_function))
      total.condition <- as.data.frame(ddply(df3.meta3,"Var1",numcolwise(sum)))
      emtpy <- matrix(nrow =dim(df3.meta3)[1],ncol=dim(total.condition)[1])

      for (j in 1:dim(total.condition)[1]){
        for (i in 1:dim(df3.meta3)[1]){
          emtpy[i,j] <- ifelse(df3.meta3$Var1[i]==total.condition$Var1[j],total.condition$Freq[j],F)

        }
      }

      df3.meta3$n <- df3.meta3$Freq/rowSums(emtpy)

      ggplot(df3.meta3,aes(x="", y=n, fill=Var2, group = Var1)) +
        geom_bar(stat="identity", width=1)+
        coord_polar("y", start=0)  +
        theme_void(20) +
        facet_wrap(~Var1, nrow = input$wrap_row) +
        theme(
          legend.key.size = unit(1, 'cm'))+
        scale_fill_manual(values = df.col, na.value = input$NA_col_analysis) +
        theme(strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
              legend.text = element_text(size = input$Bar_legend_size, family = input$font_type),
              legend.title = element_blank()

        )

      # theme(strip.background =element_rect(fill=input$strip.colour.tree))+
      # theme(strip.text = element_text(colour = input$strip.text.colour.tree))
      # scale_pattern_fill_manual(values = class_col.t)
    })

    output$Classification_clonotype_pie <- renderPlot({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      Pie_chart_Class()
    })
    output$downloadPlot_Classification_clonotype_pie <- downloadHandler(
      filename = function() {
        x <- today()
        paste("GEx_pie_",input$Colour_By_this_overview,"_",x,".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Classification_clonotype_pie,height=input$height_Classification_clonotype_pie, onefile = FALSE) # open the pdf device
        plot(Pie_chart_Class())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Classification_clonotype_pie <- downloadHandler(
      filename = function() {
        x <- today()
        paste("GEx_pie_",input$Colour_By_this_overview,"_",x,".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Classification_clonotype_pie,
            height = input$height_png_Classification_clonotype_pie,
            res = input$resolution_PNG_Classification_clonotype_pie)
        plot(Pie_chart_Class())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    # chi-square of expression -----
    chi_squ <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      meta.data <- sc@meta.data
      totals <- meta.data[,names(meta.data) %in% c(input$Samp_col,input$Split_group_by_overview)]

      totals %>%
        select(all_of(input$Samp_col), everything())

      names(totals) <- c("samps","split")
      totals$split <- as.character(totals$split)
      # totals <- totals[!totals$Function %in% c("NA"),]
      # totals <- totals[!totals$groups %in% c("NA"),]
      # totals <- totals[totals$groups %in% input$ID_Column_factor,]
      # totals$samps <- factor(totals$samps,levels = c(input$ID_Column_factor))
      totals
    })

    output$Chi_tab_before <-renderPrint({
      chi_squ()
      totals <- chi_squ()
      tb_totals <- table(totals$samps,totals$split)
      as.data.frame(tb_totals)
      chisq <- chisq.test(tb_totals)
      chisq
    })

    output$Post_hoc_chi <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      totals <- chi_squ()
      tb_totals <- table(totals$samps,totals$split)
      df <- tb_totals
      posthoc <- chisq.posthoc.test(df)
      resid <- subset(posthoc,posthoc$Value=="Residuals")
      pval <- subset(posthoc,posthoc$Value=="p values")
      as.data.frame(pval)
    })

    Chi_square_plot2 <- reactive({
      totals <- as.data.frame(chi_squ())
      totals <- totals[totals$samps %in% input$ID_Column_factor,]
      tb_totals <- table(totals$samps,totals$split)
      chisq <- chisq.test(tb_totals)
      chisq
      res <- chisq$residuals
      res <- setNames(melt(res), c('y', 'x', 'residuals'))


      plot_chi <- ggplot(res, aes(x = x, y = y,size = residuals, fill = residuals))

      plot_chi <- plot_chi +
        geom_point(shape = 21, colour = "black")+
        # scale_size_area(max_size = 10) +
        scale_fill_gradient2(
          low = input$lower_col_chi,
          mid = input$mid_col_chi,
          high = input$high_col_chi,
          space = "Lab",
          na.value = "grey90",
          guide = "colourbar",
          aesthetics = "fill")+
        theme_bw()+
        theme(
          axis.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          axis.text.x = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type,angle = 90),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          axis.line.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "right",
        )+
        guides(size = "none")
      plot_chi
    })

    output$Chi_square_plot <- renderPlot({
      Chi_square_plot2()
    })


    output$downloadPlot_Chi_square_plot <- downloadHandler(
      filename = function() {
        paste("_Chi_square_plot", ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Chi_square_plot,height=input$height_Chi_square_plot, onefile = FALSE) # open the pdf device
        plot(Chi_square_plot2())
        dev.off()
      }, contentType = "application/pdf")

    output$downloadPlotPNG_Chi_square_plot <- downloadHandler(
      filename = function() {
        paste("_Chi_square_plot_",".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Chi_square_plot,height = input$height_png_Chi_square_plot,res = input$resolution_PNG_Chi_square_plot)
        plot(Chi_square_plot2())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )


    # Select top clonotype ----
    Top_clonotype_df2 <- reactive({
      sc <- UMAP_metadata_with_labs()
      req(input$selected_Indiv,input$Samp_col,input$V_gene_sc)
      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      dataframe_one<- sc@meta.data
      names(dataframe_one)[names(dataframe_one) %in% input$Samp_col] <- "ID_Column"

      if (input$by_indiv_pie_epi == "yes"){
        df3.meta <- dataframe_one[dataframe_one$ID_Column %in% input$selected_Indiv,]
      }

      else {
        df3.meta <- dataframe_one
      }

      # summarising the clonality
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$V_gene_sc]
      df3.meta$cloneCount <- 1

      BD <- df3.meta[,names(df3.meta) %in% c("cluster_name","cloneCount")]
      BD
      BD_sum <- ddply(BD,names(BD)[-c(2)] ,numcolwise(sum))
      as.data.frame(BD_sum)
      BD_sum <- BD_sum[!BD_sum$cluster_name=="_",]
      BD_sum <- BD_sum[!BD_sum$cluster_name %in% "NA",]
      names(BD_sum)[2] <- "Total_count"
      BD_sum$frequency <- BD_sum$Total_count/sum(BD_sum[,c("Total_count")], na.rm=T)

      # BD_sum$frequency <- BD_sum$Total_count/sum(BD_sum$Total_count)
      BD_sum<- BD_sum[order(BD_sum$Total_count,decreasing = T),]

      BD_sum
    })

    output$Top_clonotype_sum <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength =10, scrollX = TRUE),{
      Top_clonotype_df2()

    })


    output$download_Top_clonotype_sum <- downloadHandler(
      filename = function(){
        x = today()
        paste("clonotype_summary_table_",x,".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(Top_clonotype_df2())
        write.csv(df,file, row.names = F)
      } )

    # top clonotypes observe events -----
    observe({
      BD_sum <- Top_clonotype_df2()
      validate(
        need(nrow(BD_sum)>0,
             error_message_val1)
      )

      # if (input$analysis_set.type)
      BD_sum <- subset(BD_sum,BD_sum$Total_count>1)
      if(dim(BD_sum)[1]>50) {

        BD_sum <- BD_sum[1:50,]
      }



      updateSelectInput(
        session,
        "Selected_clonotype",
        choices=BD_sum$cluster_name,
        selected = BD_sum$cluster_name[1])

    }

    )





    select_split_by <- reactive ({
      sc <- UMAP_metadata_with_labs()
      # df <- top_clonotype_bar_code()
      validate(
        need(nrow(sc)>0,
             "upload file")
      )
      req(input$Split_group_by_)
      df <- sc@meta.data
      df2 <- as.data.frame(unique(df[names(df) %in% input$Split_group_by_]))
      df2 <- as.data.frame(df2)
      df2

    })




    observe({
      df2 <- select_split_by()
      validate(
        need(nrow(df2)>0,
             error_message_val1)
      )
      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- (df2[order(df2$V1),])

      updateSelectInput(
        session,
        "Graph_split_order",
        choices=df2,
        selected=df2
      )
    })

    # top clonptype bar graph -----
    top_clonotype_bar_code <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload")
      )
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$V_gene_sc]
      top_BD_clonotype <- df3.meta[df3.meta$cluster_name %in% input$Selected_clonotype,]
      top_BD_clonotype

    })



    ##### bar graph of top clonotypes -----
    cols_Top_bar_clonotype <- reactive({
      dtop_clonotype_bar_code <- top_clonotype_bar_code()
      dtop_clonotype_bar_code$Selected_chain <- dtop_clonotype_bar_code[,names(dtop_clonotype_bar_code) %in% input$V_gene_sc]

      # top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(unique(dtop_clonotype_bar_code$Selected_chain))
      num <- as.data.frame(num[complete.cases(num)==T,])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_bar_clonotype", i, sep="_"), paste(num[i,]), col.gg[i])
        })
      }
      else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_bar_clonotype", i, sep="_"), paste(num[i,]), col.hcl[i])
        })
      }
      else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_bar_clonotype", i, sep="_"), paste(num[i,]), col.topo[i])
        })
      }
      else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_bar_clonotype", i, sep="_"), paste(num[i,]), heat_col[i])
        })
      }
      else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_bar_clonotype", i, sep="_"), paste(num[i,]), col.terrain[i])
        })
      }
      else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_bar_clonotype", i, sep="_"), paste(num[i,]), palette_rainbow[i])
        }) }
      else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_bar_clonotype", i, sep="_"), paste(num[i,]), palette1[i])
        })

      }
      else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep="_"), paste(num[i,]), input$one.colour.default)
        })


      } # one colour

    })
    output$myPanel_Top_bar_clonotype <- renderUI({cols_Top_bar_clonotype()})
    colors_cols_Top_bar_clonotype <- reactive({
      dtop_clonotype_bar_code <- top_clonotype_bar_code()
      dtop_clonotype_bar_code$Selected_chain <- dtop_clonotype_bar_code[,names(dtop_clonotype_bar_code) %in% input$V_gene_sc]
      # top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(unique(dtop_clonotype_bar_code$Selected_chain))
      num <- as.data.frame(num[complete.cases(num)==T,])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.Top_bar_clonotype", i, sep="_")]]
      })
    })

    ggplot_top_BD_clonotype_vs_SC <- reactive({
      dtop_clonotype_bar_code <- top_clonotype_bar_code()

      req(input$Graph_split_order)

      dtop_clonotype_bar_code$Selected_group <- dtop_clonotype_bar_code[,names(dtop_clonotype_bar_code) %in% input$Split_group_by_]
      colorblind_vector <- as.data.frame(unlist(colors_cols_Top_bar_clonotype()))

      if (dim(colorblind_vector)[1]==0) {
        num <- as.data.frame(unique(dtop_clonotype_bar_code$Selected_chain))
        # num <- as.data.frame(num[complete.cases(num)==T,])

        if (input$colourtype == "default") {
          colorblind_vector <- c(gg_fill_hue(dim(num)[1]))
        } else if (input$colourtype == "hcl.colors") {
          colorblind_vector <- c(hcl.colors(dim(num)[1], palette = "viridis"))
        } else if (input$colourtype == "topo.colors") {
          colorblind_vector <- c(topo.colors(dim(num)[1]))
        } else if (input$colourtype == "heat.colors") {
          colorblind_vector <- c(heat.colors(dim(num)[1]))
        } else if (input$colourtype == "terrain.colors") {
          colorblind_vector <- c(terrain.colors(dim(num)[1]))
        } else if (input$colourtype == "rainbow") {
          colorblind_vector <- c(rainbow(dim(num)[1]))
        } else if (input$colourtype == "random") {
          colorblind_vector <- distinctColorPalette(dim(num)[1])

        }  else {

        }

      }
      colorblind_vector <- as.data.frame(colorblind_vector)
      names(colorblind_vector) <- "cols"

      dtop_clonotype_bar_code$Selected_chain2 <- dtop_clonotype_bar_code[,names(dtop_clonotype_bar_code) %in% input$V_gene_sc]
      dtop_clonotype_bar_code$Selected_chain3 <- gsub("_"," ",dtop_clonotype_bar_code$Selected_chain2)
      dtop_clonotype_bar_code$Selected_chain3 <- gsub("[.]"," ",dtop_clonotype_bar_code$Selected_chain3)

      dtop_clonotype_bar_code <- dtop_clonotype_bar_code[dtop_clonotype_bar_code$Selected_group %in% input$Graph_split_order,]
      dtop_clonotype_bar_code$Selected_group <- factor(dtop_clonotype_bar_code$Selected_group,levels = input$Graph_split_order)

      ggplot(dtop_clonotype_bar_code, aes(x=Selected_group, fill=Selected_chain3,colour = Selected_chain3, label = Selected_chain3)) +
        geom_bar() +
        theme_bw()+
        scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector$cols, na.value=input$NA_col_analysis)+
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector$cols, na.value=input$NA_col_analysis)+
        # scale_alpha_manual(values = rep(1,length(unique(dtop_clonotype_bar_code$Selected_chain))), na.value=0.5)+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        theme(
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=90),
          axis.title.x = element_blank(),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.position = input$legend_position,
          legend.title = element_blank()
        ) +
        guides(color = "none", size = "none")

    })

    output$top_clonotype <- renderPlot({
      ggplot_top_BD_clonotype_vs_SC()
    })

    output$downloadPlot_top_clonotype <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Selected_clonotype,"top_clonotype_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_top_clonotype,height=input$height_top_clonotype, onefile = FALSE) # open the pdf device
        plot(ggplot_top_BD_clonotype_vs_SC())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_top_clonotype <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Selected_clonotype,"top_clonotype_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_top_clonotype,
            height = input$height_png_top_clonotype,
            res = input$resolution_PNG_top_clonotype)
        plot(ggplot_top_BD_clonotype_vs_SC())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    #### top clonotype pie -----
    cols_Top_pie_clonotype <- reactive({
      top_BD_cluster <-  top_clonotype_bar_code()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this]
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this]
      top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% input$Split_group_by_]
      # top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
      num <- (unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[!num == "NA"])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep="_"), paste(num[i,]), col.gg[i])
        })
      }
      else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep="_"), paste(num[i,]), col.hcl[i])
        })
      }
      else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep="_"), paste(num[i,]), col.topo[i])
        })
      }
      else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep="_"), paste(num[i,]), heat_col[i])
        })
      }
      else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep="_"), paste(num[i,]), col.terrain[i])
        })
      }
      else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep="_"), paste(num[i,]), palette_rainbow[i])
        }) }
      else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep="_"), paste(num[i,]), palette1[i])
        })

      }
      else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep="_"), paste(num[i,]), input$one.colour.default)
        })


      } # one colour

    })
    output$myPanel_Top_pie_clonotype <- renderUI({cols_Top_pie_clonotype()})
    colors_cols_Top_pie_clonotype <- reactive({
      top_BD_cluster <-  top_clonotype_bar_code()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this]
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this]
      top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% input$Split_group_by_]
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[!num == "NA"])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.Top_pie_clonotype", i, sep="_")]]
      })
    })

    pie_ag_bd_table <- reactive({
      top_BD_cluster <-  top_clonotype_bar_code()
      top_BD_cluster
      # req(input$Colour_By_this,input$Split_group_by_,input$Graph_split_order)

      top_BD_cluster$Selected_function <-  as.character((top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this]))
      top_BD_cluster$Selected_group <-  as.character((top_BD_cluster[,names(top_BD_cluster) %in% input$Split_group_by_]))
      top_BD_cluster <- top_BD_cluster[top_BD_cluster$Selected_group %in% input$Graph_split_order,]
      top_BD_cluster$Selected_group <- factor(top_BD_cluster$Selected_group,levels = input$Graph_split_order)

      df3.meta3 <- top_BD_cluster[,names(top_BD_cluster) %in% c("Selected_function","Selected_group")]
      df3.meta3$count <- 1
      df3.meta3
      total.condition <- ddply(df3.meta3,"Selected_group",numcolwise(sum))
      total.condition
      #
      total.group.condition <- ddply(df3.meta3,c("Selected_group","Selected_function"),numcolwise(sum))
      total.group.condition
      #
      nrow_count <- dim(total.group.condition)[1]
      ncol_count <- dim(total.condition)[1]
      #
      emtpy_mat <- matrix(nrow =nrow_count,ncol=ncol_count)

      for (i in 1:ncol_count){
        emtpy_mat[,i] <- (ifelse(total.group.condition$Selected_group==total.condition$Selected_group[i], total.condition$count[i],F))

      }
      as.data.frame(emtpy_mat)
      #
      total.group.condition$n <- total.group.condition$count/rowSums(emtpy_mat)
      total.group.condition$Selected_function <- ifelse(total.group.condition$Selected_function=="NA",NA,total.group.condition$Selected_function)
      total.group.condition
    })

    output$Top_clonotype_Labs <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength =10, scrollX = TRUE),{
      pie_ag_bd_table()

    })

    Pie_chart_alpha_gamma <- reactive({
      df.col <- unlist(colors_cols_Top_pie_clonotype())

      total.group.condition <- pie_ag_bd_table()


      total.group.condition$Selected_function <- gsub("_"," ",total.group.condition$Selected_function)
      total.group.condition$Selected_function <- gsub("[.]"," ",total.group.condition$Selected_function)

      # ggplot(dtop_clonotype_bar_code, aes(x=Selected_group, fill=Selected_chain3,colour = Selected_chain3, label = Selected_chain3)) +
      #   geom_bar() +
      #   theme_bw()+
      #   scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector, na.value=input$NA_col_analysis)+
      #   scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector, na.value=input$NA_col_analysis)+

      ggplot(total.group.condition,aes(x="", y=n, fill=as.character(Selected_function), group =as.character(Selected_group), label = Selected_function)) +
        geom_bar(stat="identity", width=1)+
        coord_polar("y", start=0)  +
        theme_void(20) +
        facet_wrap(~Selected_group, nrow = input$wrap_row) +
        theme(
          legend.key.size = unit(1, 'cm'),
          legend.title = element_blank()) +
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 10),values = df.col, na.value = input$NA_col_analysis) +
        theme(strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
              legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.position = input$legend_position,
              legend.title = element_blank()
        )

    })

    output$top_clonotype_pie <- renderPlot({
      if (input$Plot_type_selected=="pie") {
        Pie_chart_alpha_gamma()
      }
      else {
        UMAP_chart_alpha_gamma()


      }
    })

    #### top  UMAP -----
    UMAP_chart_alpha_gamma <- reactive({
      top_BD_cluster <-  top_clonotype_bar_code()
      req(input$Colour_By_this,input$Split_group_by_,input$Graph_split_order)
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this]
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]

      top_BD_cluster$Selected_function <-  as.character(top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this])

      top_BD_cluster$Selected_group <-  as.character(top_BD_cluster[,names(top_BD_cluster) %in% input$Split_group_by_])
      #
      top_BD_cluster <- top_BD_cluster[top_BD_cluster$Selected_group %in% input$Graph_split_order,]
      top_BD_cluster$Selected_group <- factor(top_BD_cluster$Selected_group,levels = input$Graph_split_order)
      top_BD_cluster$Selected_function <- ifelse(top_BD_cluster$Selected_function=="NA",NA,top_BD_cluster$Selected_function)

      df.col <- unlist(colors_cols_Top_pie_clonotype())

      plot_pie <- ggplot(data = top_BD_cluster,aes(x=UMAP_1,UMAP_2,colour=Selected_function))+
        geom_point()+
        scale_color_manual(values = df.col, na.value=input$NA_col_analysis) +
        # scale_size_manual(values = rep(input$size_selected_top,dim(top_BD_cluster)[1]),na.value = 1) +
        theme_bw()+
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )+
        scale_size(guide = 'none')

      # if (input$Split_by_group2=="no") {
      #   plot_pie
      # }
      # else {

      plot_pie <- plot_pie + facet_wrap(~top_BD_cluster$Selected_group, nrow = input$wrap_row)
      #

      plot_pie

    })

    # plot pie and UMAP ----

    output$downloadPlot_top_clonotype_pie <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Selected_clonotype,"_top_clonotype_pie_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_top_clonotype_pie ,height=input$height_top_clonotype_pie , onefile = FALSE) # open the pdf device
        if (input$Plot_type_selected=="pie") {
          plot(Pie_chart_alpha_gamma())
        }
        else {
          plot(UMAP_chart_alpha_gamma())
        }
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_top_clonotype_pie  <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Selected_clonotype,"_top_clonotype_pie_",x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_top_clonotype_pie ,
            height = input$height_png_top_clonotype_pie ,
            res = input$resolution_PNG_top_clonotype_pie )

        if (input$Plot_type_selected=="pie") {
          plot(Pie_chart_alpha_gamma())
        }
        else {
          plot(UMAP_chart_alpha_gamma())
        }


        dev.off()},   contentType = "application/png" # MIME type of the image
    )



    ### Add in the TCR clustering to the seurat object ----

    #### Ridge plot (expression per T cell) ----
    vals_Ridge_top <- reactiveValues(output_stats=NULL)

    compare.stat <- reactive({
      sc <-input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )

      df <- sc@meta.data
      df
      unique.df <- (df[,names(df) %in% c(input$Samp_col,input$V_gene_sc) ])
      names(unique.df) <- c("group","chain")

      # unique.df <- unique.df[unique.df$group %in% "LTR35-12",]

      unique.df <- subset(unique.df,unique.df$chain != "NA")
      unique.df <- subset(unique.df,unique.df$group != "NA")
      unique.df$cloneCount <- 1
      unique.df
      df_unique_sum <- ddply(unique.df,names(unique.df)[-c(3)] ,numcolwise(sum))
      df_unique_sum <- df_unique_sum[order(df_unique_sum$cloneCount, decreasing = T),]

      sc@meta.data$Vgene <- sc@meta.data[,names(sc@meta.data) %in% input$V_gene_sc]

      name.clone <- input$Selected_clonotype
      sc@meta.data$Gene_select <- ifelse(sc@meta.data$Vgene %in% name.clone,name.clone,"unselected")
      sc@meta.data
      unique(sc@meta.data$Gene_select)
      Idents(object = sc) <- sc@meta.data$Gene_select
      as.data.frame(Idents(object = sc))
      min.pct.expression<- input$min_point_ #standard setting: 0.25§
      min.logfc<-  input$LogFC_ #0.25 is standard
      p.val.cutoff <-  input$pval_top #(1/10^3) is standard, use (1/10^0) to ignore

      cluster.names <- unique(Idents(sc))[order(unique(Idents(sc)))]
      as.data.frame(cluster.names)
      # # print(paste0("calculating markers for cluster ",name.clone,". Total: ",length(cluster.names)," clusters"))
      markers.fm.list <- FindMarkers(sc, ident.1 = name.clone, min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
      markers.fm.list
      # markers.fm.list2 <- subset(markers.fm.list,markers.fm.list$p_val_adj < input$pval.ex.filter)
      # as.data.frame(markers.fm.list2)


    })
    output$Ridge_chart_alpha_gamma_stat_comp <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength =10, scrollX = TRUE),{

      sc <-input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )
      compare.stat()

    })

    output$downloaddf_FindMarker_Top <- downloadHandler(
      filename = function(){
        paste("Stats_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(compare.stat())
        write.csv(df,file, row.names = T)
      } )


    observeEvent(input$run_string.data_Exp_top,{
      # df <- compare.stat()

      sc <- input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )

      if(input$SeuratVersion == "Version 4") {
        df <- as.data.frame(sc@assays$RNA$scale.data)
      } else {
        df <- as.data.frame(sc@assays$RNA$scale.data)
      }


      updateSelectInput(
        session,
        "string.data_Exp_top",
        choices=df,
        selected = df[1])

    })


    output$test.table_ridge <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength =10, scrollX = TRUE),{
      sc <- input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )

      if(input$SeuratVersion == "Version 4") {
        df <- as.data.frame(sc@assays$RNA$scale.data)
      } else {
        df <- rownames(sc@assays$RNA$scale.data)
      }
      req(df)
      as.data.frame(df)[1:6]

    })

    #ridge plot ------
    Ridge_chart_alpha_gamma_df <- reactive({
      sc <-input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )

      if(input$SeuratVersion == "Version 4") {
        df <- as.data.frame(sc@assays$RNA$scale.data)
      } else {
        df <- as.data.frame(sc@assays$RNA$scale.data)
      }

      # if(input$SeuratVersion == "Version 4") {
      #   df= as.data.frame(sc@assays$RNA$scale.data)
      # } else {
      #
      #   df = as.data.frame(sc[["RNA"]]@layers$scale.data)
      # }

      req(input$string.data_Exp_top)

      MainTcell <- as.data.frame(t(df))
      meta.data <- as.data.frame(sc@meta.data)
      # names(MainTcell) <- ifelse(grepl("pAbO",names(MainTcell)),names(MainTcell),toupper(names(MainTcell)))
      MainTcell$Cell_Index <- rownames(MainTcell)
      # names(MainTcell) <- gsub("CD8B.","CD8B",names(MainTcell))

      gene_df <- MainTcell[,names(MainTcell) %in% c("Cell_Index",input$string.data_Exp_top)]

      top_BD_cluster <-  top_clonotype_bar_code()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Split_group_by_]
      top_BD_cluster <- top_BD_cluster[top_BD_cluster$Selected_function %in% input$Graph_split_order,]
      top_BD_cluster$Selected_function <-  as.character(top_BD_cluster[,names(top_BD_cluster) %in% input$Split_group_by_])
      top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = input$Graph_split_order)

      merge(top_BD_cluster,gene_df,by = "Cell_Index")
    })

    Ridge_chart_alpha_gamma_plot <- reactive({
      df <- Ridge_chart_alpha_gamma_df()
      df$expressed <- ifelse(df[,names(df) %in% input$string.data_Exp_top] >input$Gre_ex,"expressed","Not expressed")

      if(input$restric_ex == T) {
        df2 <- subset(df,df$expressed=="expressed")
      }
      else(
        df2 <- df
      )

      df2$Selected_function <- factor(df2$Selected_function,levels = input$Graph_split_order)

      ggplot(df2, aes(x = get(input$string.data_Exp_top), y = Selected_function, fill = Selected_function)) +
        geom_density_ridges() +
        theme_ridges() +
        theme(legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              # legend.title = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")+
        ggtitle(input$string.data_Exp_top)

    })

    Ridge_chart_alpha_gamma_stat_table <- reactive({
      df <- Ridge_chart_alpha_gamma_df()
      df$expressed <- ifelse(df[,names(df) %in% input$string.data_Exp_top] >input$Gre_ex,"expressed","Not expressed")

      if(input$restric_ex == T) {
        df2 <- subset(df,df$expressed=="expressed")
      }
      else(
        df2 <- df
      )
      df2[is.na(df2)] <- "-"
      at <- TukeyHSD(aov(get(input$string.data_Exp_top)~ get(input$Split_group_by_),data = df2))
      tab <- as.data.frame(at[1])
      names(tab) <- c("diff" ,"lwr","upr","p.adj")
      tab$stat <- ifelse(tab$p.adj<0.0001,"****",
                         ifelse(tab$p.adj<0.001,"***",
                                ifelse(tab$p.adj<0.01,"**",
                                       ifelse(tab$p.adj<0.05,"*","NS"
                                       ))))
      tab
    })

    output$Ridge_chart_alpha_gamma_stat <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength =10, scrollX = TRUE),{
      Ridge_chart_alpha_gamma_stat_table()

    })

    output$downloaddf_clusTCR_GEx <- downloadHandler(
      filename = function(){
        paste("Stats_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(Ridge_chart_alpha_gamma_stat_table())
        write.csv(df,file, row.names = T)
      } )
    ### Violin plots top clonotype -----
    Violin_plot_table <- reactive({
      sc <-input.data_sc_pro()
      need(nrow(sc)>0,
           error_message_val_sc)

      if(input$SeuratVersion == "Version 4") {
        df <- as.data.frame(sc@assays$RNA$scale.data)
      } else {
        df <- as.data.frame(sc@assays$RNA$scale.data)
      }

      MainTcell <- as.data.frame(t(df))
      meta.data <- as.data.frame(sc@meta.data)
      # names(MainTcell) <- ifelse(grepl("pAbO",names(MainTcell)),names(MainTcell),toupper(names(MainTcell)))
      MainTcell$Cell_Index <- rownames(MainTcell)
      # names(MainTcell) <- gsub("CD8B.","CD8B",names(MainTcell))

      gene_df <- MainTcell[,names(MainTcell) %in% c("Cell_Index",input$string.data_Exp_top)]

      top_BD_cluster <-  top_clonotype_bar_code()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Split_group_by_]
      top_BD_cluster <- top_BD_cluster[top_BD_cluster$Selected_function %in% input$Graph_split_order,]
      top_BD_cluster$Selected_function <-  as.character(top_BD_cluster[,names(top_BD_cluster) %in% input$Split_group_by_])
      top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = input$Graph_split_order)

      merge(top_BD_cluster,gene_df,by = "Cell_Index")


    })


    Violin_chart_alpha_gamma_plot <- reactive({
      df <- Ridge_chart_alpha_gamma_df()
      # df <- df[df$Selected_group %in% input$Graph_split_order,]
      # top_BD_cluster$Selected_group <- factor(top_BD_cluster$Selected_group,levels = input$Graph_split_order)
      df$expressed <- ifelse(df[,names(df) %in% input$string.data_Exp_top] >input$Gre_ex,"expressed","Not expressed")

      if(input$restric_ex == T) {
        df2 <- subset(df,df$expressed=="expressed")
      }
      else(
        df2 <- df
      )
      df2$Selected_function <- factor(df2$Selected_function,levels = input$Graph_split_order)
      ggplot(df2, aes(y = get(input$string.data_Exp_top), x = Selected_function, fill = Selected_function)) +
        geom_violin() +
        geom_jitter(height = 0, width = 0.1)+
        theme(legend.position = "none",

        )+
        theme_bw() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=90),
          axis.title.x = element_blank(),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_text(colour="black",size=20,family=input$font_type),
          legend.position = "none",
        )+
        ggtitle(input$string.data_Exp_top)

    })

    # vs non-selected clonotypes -----

    Ridge_chart_alpha_gamma_df_all <- reactive({
      sc <-UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )
      if(input$SeuratVersion == "Version 4") {
        df <- as.data.frame(sc@assays$RNA$scale.data)
      } else {
        df <- as.data.frame(sc@assays$RNA$scale.data)
      }

      MainTcell <- as.data.frame(t(df))
      meta.data <- as.data.frame(sc@meta.data)
      MainTcell$Cell_Index <- rownames(MainTcell)
      gene_df <- MainTcell[,names(MainTcell) %in% c("Cell_Index",input$string.data_Exp_top)]
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$V_gene_sc]
      df3.meta$selected_top_clonotype <- ifelse(df3.meta$cluster_name==input$Selected_clonotype,input$Selected_clonotype,"other")
      top_BD_clonotype2 <- merge(df3.meta,gene_df,by="Cell_Index")
      top_BD_clonotype2

    })
    Ridge_chart_alpha_gamma_plot_comp <- reactive({
      df <- Ridge_chart_alpha_gamma_df_all()
      df$expressed <- ifelse(df[,names(df) %in% input$string.data_Exp_top] >input$Gre_ex,"expressed","Not expressed")

      if(input$restric_ex == T) {
        df2 <- subset(df,df$expressed=="expressed")
      }
      else(
        df2 <- df
      )

      ggplot(df2, aes(x = get(input$string.data_Exp_top), y = selected_top_clonotype, fill = selected_top_clonotype)) +
        geom_density_ridges() +
        theme_ridges() +
        theme(legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.title = element_blank(),
              axis.text.x = element_blank(),
              # axis.text.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")

    })

    Violin_chart_alpha_gamma_plot_comp <- reactive({
      df <- Ridge_chart_alpha_gamma_df_all()
      df$expressed <- ifelse(df[,names(df) %in% input$string.data_Exp_top] >input$Gre_ex,"expressed","Not expressed")

      if(input$restric_ex == T) {
        df2 <- subset(df,df$expressed=="expressed")
      }
      else(
        df2 <- df
      )

      ggplot(df2, aes(y = get(input$string.data_Exp_top), x = selected_top_clonotype, fill = selected_top_clonotype)) +
        geom_violin() +
        geom_jitter(height = 0, width = 0.1)+
        theme(legend.position = "none",

        )+
        theme_bw() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )

    })
    output$Ridge_chart_alpha_gamma_plot_out <- renderPlot({
      if (input$plot_type_ridgvi =="Ridge (selected clonotype)") {
        Ridge_chart_alpha_gamma_plot()
      }
      else if (input$plot_type_ridgvi =="Ridge (compare)") {
        Ridge_chart_alpha_gamma_plot_comp()

      }

      else if (input$plot_type_ridgvi =="Violin (compare)") {

        Violin_chart_alpha_gamma_plot_comp()

      }
      else{
        Violin_chart_alpha_gamma_plot()
      }
    })

    output$downloadPlot_Ridge_chart_alpha_gamma_plot_out <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$string.data_Exp_top,"_",input$Selected_clonotype,"_",input$plot_type_ridgvi,"_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Ridge_chart_alpha_gamma_plot_out,height=input$height_Ridge_chart_alpha_gamma_plot_out, onefile = FALSE) # open the pdf device
        if (input$plot_type_ridgvi =="Ridge (selected clonotype)") {
          df <- Ridge_chart_alpha_gamma_plot()
        }
        else if (input$plot_type_ridgvi =="Ridge (compare)") {
          df <- Ridge_chart_alpha_gamma_plot_comp()
        }
        else if (input$plot_type_ridgvi =="Violin (compare)") {
          df <- Violin_chart_alpha_gamma_plot_comp()
        }
        else{
          df <- Violin_chart_alpha_gamma_plot()
        }
        plot(df)
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Ridge_chart_alpha_gamma_plot_out  <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$string.data_Exp_top,"_",input$Selected_clonotype,"_",input$plot_type_ridgvi,"_",x, ".png", sep = "")
        # paste("_",input$epitope_umap_selected,"_",input$epitope_umap_selected2,"_Heatmap_", x, "", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Ridge_chart_alpha_gamma_plot_out,
            height = input$height_png_Ridge_chart_alpha_gamma_plot_out,
            res = input$resolution_PNG_Ridge_chart_alpha_gamma_plot_out)
        if (input$plot_type_ridgvi =="Ridge (selected clonotype)") {
          df <- Ridge_chart_alpha_gamma_plot()
        }
        else if (input$plot_type_ridgvi =="Ridge (compare)") {
          df <- Ridge_chart_alpha_gamma_plot_comp()

        }

        else if (input$plot_type_ridgvi =="Violin (compare)") {

          df <- Violin_chart_alpha_gamma_plot_comp()

        }
        else{
          df <- Violin_chart_alpha_gamma_plot()
        }
        plot(df)
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    ### top dot plot -----
    all_expression_plot_top <- reactive({
      sc <-input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_sc))

      df <- sc@meta.data
      unique.df <- (df[,names(df) %in% c(input$Samp_col,input$V_gene_sc) ])
      names(unique.df) <- c("group","chain")

      # unique.df <- unique.df[unique.df$group %in% "LTR35-12",]

      unique.df <- subset(unique.df,unique.df$chain != "NA")
      unique.df <- subset(unique.df,unique.df$group != "NA")
      unique.df$cloneCount <- 1

      df_unique_sum <- ddply(unique.df,names(unique.df)[-c(3)] ,numcolwise(sum))
      df_unique_sum <- df_unique_sum[order(df_unique_sum$cloneCount, decreasing = T),]

      sc@meta.data$Vgene <- sc@meta.data[,names(sc@meta.data) %in% input$V_gene_sc]

      name.clone <- input$Selected_clonotype

      sc@meta.data$Gene_select <- ifelse(sc@meta.data$Vgene %in% name.clone,input$name_clonotype_selected,"unselected")
      unique(sc@meta.data$Gene_select)
      Idents(object = sc) <- sc@meta.data$Gene_select
      if (input$restict_no_points == F ) {
        list.names <- rownames(compare.stat())
      }

      else {
        list.names <- rownames(compare.stat())
        list.names <- list.names[1:input$pval.ex.top_genes]
      }

      size_legend = input$Bar_legend_size-2

      DotPlot(sc, features = list.names) +
        RotatedAxis() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size, angle = 90),
          axis.title.x = element_blank(),
          legend.title = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.text = element_text(colour="black", size=size_legend,family=input$font_type),
          legend.position = input$legend_position,
        ) +
        scale_colour_gradient2(low = input$low.dotplot, mid = input$middle.dotplot, high = input$high.dotplot)+
        scale_x_discrete(labels = label_wrap(20)) +
        scale_y_discrete(labels = label_wrap(20))


    })

    output$all_expression_dotplot_top <- renderPlot({
      all_expression_plot_top()
    })

    output$downloadPlot_all_expression_dotplot_top <- downloadHandler(
      filename = function() {
        paste(input$Selected_clonotype,"_dotplot","_",today(), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_all_expression_dotplot_top,height=input$height_all_expression_dotplot_top, onefile = FALSE) # open the pdf device
        df <- all_expression_plot_top()
        plot(df)
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_all_expression_dotplot_top  <- downloadHandler(
      filename = function() {
        paste(input$Selected_clonotype,"_dotplot","_",today(), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_all_expression_dotplot_top,
            height = input$height_png_all_expression_dotplot_top,
            res = input$resolution_PNG_all_expression_dotplot_top)
        df <- all_expression_plot_top()

        plot(df)
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    # Over representation analysis for top clonotypes -----

    Over_rep_Top_clones_old <- reactive({
      sc <-input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )

      df <- sc@meta.data
      # require()

      geneSet <- read.csv(system.file("OverRep","GeneSets.csv",package = "STEGO.R"),header = T)

      background.genes.name <- as.data.frame(rownames(sc@assays$RNA$scale.data))
      names(background.genes.name) <- "V1"
      background.genes <- length(rownames(sc@assays$RNA$scale.data))

      geneSet$background.genes <- background.genes

      markers.fm.list <- compare.stat()

      DEx.genes <- as.data.frame(rownames(markers.fm.list))
      names(DEx.genes) <- "V1"

      total.sig <- length(DEx.genes$V1)
      geneSet$total.sig <- length(DEx.genes$V1)

      geneSet$background.geneset <- NA
      geneSet$background.geneset.name <- NA
      geneSet$in.geneset <- NA
      geneSet$in.geneset.name <- NA

      if(input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        geneSet$GeneSet <- gsub("-",".",geneSet$GeneSet)
      }

      if(input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        require(stringr)
        geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
      }

      for (i in 1:dim(geneSet)[1]) {
        # listed GeneSet
        message(paste("GeneSet: ", i))
        Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet,";")[i])
        names(Gene.set.testing) <- "V1"
        Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
        names(Gene.set.testing2) <- "V1"
        # message(paste(dim(Gene.set.testing2)[1],"GeneSet total"))
        # background genes
        background.overlap <- merge(Gene.set.testing2,background.genes.name,by= "V1")
        # message(paste(dim(background.overlap)[1],"in Background"))
        geneSet$background.geneset[i] <- length(background.overlap$V1)
        geneSet$background.geneset.name[i] <- as.character(paste(unlist(background.overlap[1]), collapse=';'))
        # in sig gene list
        overlap <- merge(background.overlap,DEx.genes,by= "V1")
        # message(paste(dim(overlap)[1],"# Sig genes"))
        geneSet$in.geneset[i] <- length(overlap$V1)
        geneSet$in.geneset.name[i] <- as.character(paste(unlist(overlap[1]), collapse=';'))

      }

      geneSet2 <- subset(geneSet,geneSet$in.geneset>0)

      for (i in 1:dim(geneSet2)[1]) {
        tota.gene.set <- geneSet2$background.geneset[i] # genes that are identified in background
        tota.gene.set
        in.geneset <-  geneSet2$in.geneset[i]# DEx in geneset

        background.genes
        not.in.total <- background.genes - tota.gene.set
        not.in.geneset.sig <- total.sig - in.geneset
        d <- data.frame( gene.in.interest=c( in.geneset, not.in.geneset.sig),gene.not.interest=c( tota.gene.set, not.in.total))
        row.names(d) <- c("In_category", "not_in_category")

        if (in.geneset>0) {
          geneSet2$p.val[i] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
          geneSet2$lowerCI[i] <-  unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
          geneSet2$upperCI[i] <-unlist(fisher.test(d)$conf.int)[2]
          geneSet2$OR[i] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1],3)
        }

        else {
          geneSet2$p.value[i] <- "-"
          geneSet2$lowerCI[i] <-  "-"
          geneSet2$upperCI[i] <- "-"
          geneSet2$OR[i] <- "-"
        }
        # message(print(d))
        # message(print(round(prop.table(d),3)))
      }

      geneSet2

    })


    Over_rep_Top_clones<- reactive({
      geneSet2 <- Over_rep_Top_clones_old()

      validate(
        need(nrow(geneSet2)>0,
             error_message_val_sc)
      )
      geneSet2 <- geneSet2[order(geneSet2$p.val,decreasing = F),]
      geneSet2 <- subset(geneSet2,geneSet2$in.geneset>=input$in.geneset.cutoff_top)
      geneSet2 <- subset(geneSet2,geneSet2$p.val<=input$p.val_cutoff_top)
      geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
      geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
      # name.list <- c("Geneset_ID","p.val","FDR","Bonferroni","OR","lowerCI","upperCI","in.geneset.name","in.geneset","background.geneset","total.sig","background.genes","background.geneset.name")
      # geneSet2 <- geneSet2 %>%
      #     select(all_of(name.list), everything())

      geneSet2
    })

    output$Over_rep_Top_clones_Tab <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1,2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      Over_rep_Top_clones()
    })

    output$downloadtb_over.rep.Top_Ex <- downloadHandler(
      filename = function(){

        paste0(input$Selected_clonotype,"_",today(),"_over_rep.csv")
      },
      content = function(file){
        df <- as.data.frame(Over_rep_Top_clones())
        write.csv(df,file, row.names = F)
      } )

    # Expanded TCR interrogation (regarless of TCR sequence ) -------
    observe({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "Samp_col_expanded",
        choices=names(df3.meta),
        selected = "Sample_Name")
    })

    select_group_metadata_expanded <- reactive ({
      sc <- UMAP_metadata_with_labs()


      validate(
        need(nrow(sc)>0,
             "upload file")
      )
      req(input$Samp_col_expanded)
      df <- sc@meta.data
      df2 <- as.data.frame(unique(df[names(df) %in% input$Samp_col_expanded]))
      df2 <- as.data.frame(df2)
      df2

    })

    observe({
      df2 <- select_group_metadata_expanded()
      validate(
        need(nrow(df2)>0,
             error_message_val1)
      )

      req(input$Samp_col2)

      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- as.data.frame(df2[order(df2$V1),])
      names(df2) <- "V1"
      df2
      df2 <- subset(df2,df2$V1 != "NA")

      updateSelectInput(
        session,
        "ID_Column_factor_expanded",
        choices=df2$V1,
        selected = df2$V1[1]
      )
    })


    #   observe({
    #     sc <- UMAP_metadata_with_labs()
    #     validate(
    #       need(nrow(sc)>0,
    #            error_message_val_UMAP)
    #     )
    #     df3.meta <- sc@meta.data
    #       updateSelectInput(
    #         session,
    #         "Samp_col_expanded",
    #         choices=names(df3.meta),
    #         selected = "orig.ident")
    # })

    select_group_metadata_ex <- reactive ({
      df <- Expansion_check_table()

      validate(
        need(nrow(df)>0,
             error_message_val1)
      )
      req(df)
      # req(input$Samp_col_expanded,input$ID_Column_factor_expanded)
      df <- Expansion_check_table()
      df2 <- df$expansion.status
      df2 <- as.data.frame(df2)
      df2

    })

    observe({
      df2 <- select_group_metadata_ex()

      validate(
        need(nrow(df2)>0,
             error_message_val1)
      )
      req(df2)

      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- as.data.frame(df2[order(df2$V1),])
      names(df2) <- "V1"
      df2
      df2 <- subset(df2,df2$V1 != "NA")
      df3 <- unique(df2$V1)
      updateSelectizeInput(
        session,
        "selected_Indiv_Ex_1",
        choices=df3,
        selected = df3[1]
      )
    })

    observe({
      df2 <- select_group_metadata_ex()
      validate(
        need(nrow(df2)>0,
             error_message_val1)
      )
      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- as.data.frame(df2[order(df2$V1),])
      names(df2) <- "V1"
      df2
      df2 <- subset(df2,df2$V1 != "NA")
      df3 <- unique(df2$V1)
      updateSelectizeInput(
        session,
        "selected_Indiv_Ex_2",
        choices=df3,
        selected = df3[2]
      )
    })

    ### table check ------

    Vals_Expans_table <- reactiveValues(output_expTab=NULL)

    observeEvent(input$update_Expansion_run,{

    })

    Expansion_check_table <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )

      req(input$Samp_col_expanded,input$ID_Column_factor_expanded,input$V_gene_sc,input$cut.off_expanded, sc)

      message(paste("Calculating expansion"))

      Expansion.meta.data <- sc@meta.data
      names(Expansion.meta.data)[names(Expansion.meta.data) %in% input$Samp_col_expanded] <- "Samp"
      Expansion.meta.data <- Expansion.meta.data[Expansion.meta.data$Samp %in% c(input$ID_Column_factor_expanded),]
      names(Expansion.meta.data)[names(Expansion.meta.data) %in% "Samp"] <- input$Samp_col_expanded
      Expansion.meta.data
      CDR3_Vgene_all <- as.data.frame(Expansion.meta.data[,names(Expansion.meta.data) %in% input$V_gene_sc])
      names(CDR3_Vgene_all) <- "V1"
      CDR3_Vgene_all$count <- 1
      total.condition <- as.data.frame(ddply(CDR3_Vgene_all,"V1",numcolwise(sum)))
      total.condition <- total.condition[total.condition$V1 != "NA",]
      CDR3_Vgene_all2 <- as.data.frame(Expansion.meta.data[,names(Expansion.meta.data) %in% c(input$V_gene_sc,input$Samp_col_expanded)])
      names(CDR3_Vgene_all2)[names(CDR3_Vgene_all2) %in% input$V_gene_sc] <- "V1"
      names(CDR3_Vgene_all2)[names(CDR3_Vgene_all2) %in% input$Samp_col_expanded] <- "V2"
      CDR3_Vgene_all2$count.group <- 1
      total.condition.group <- as.data.frame(ddply(CDR3_Vgene_all2,c("V1","V2"),numcolwise(sum)))
      total.condition.group <- total.condition.group[total.condition.group$V1 != "NA",]
      total.condition.group.counts <- merge(total.condition.group,total.condition,by = "V1")
      total.condition.group.counts <- total.condition.group.counts[order(total.condition.group.counts$count.group, decreasing = T),]
      total.condition.group.counts$expand.singlets <- ifelse(total.condition.group.counts$count.group>input$cut.off_expanded,"Ex","NEx")
      names(total.condition.group.counts)[names(total.condition.group.counts) %in% "V1"] <-input$V_gene_sc
      names(total.condition.group.counts)[names(total.condition.group.counts) %in% "V2"] <-input$Samp_col_expanded
      total.condition.group.counts
      umap.meta2 <- merge(Expansion.meta.data, total.condition.group.counts, by = c(input$Samp_col_expanded,input$V_gene_sc),sort = F,all.x =T)
      rownames(umap.meta2) <- umap.meta2$Cell_Index
      message(paste("Finished calculating expansion"))
      names(umap.meta2)[names(umap.meta2) %in% input$Samp_col_expanded] <- "Selected_Status"
      umap.meta2$expansion.status <- paste0(umap.meta2$Selected_Status,"_",umap.meta2$expand.singlets)
      names(umap.meta2)[names(umap.meta2) %in% "Selected_Status"] <- input$Samp_col_expanded
      umap.meta2
    })


    output$Expansion_check2 <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      Expansion_check_table()
    })
    ### add in colouring specific to Expanded

    output$classification_to_add2 <- renderUI({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      req(Expansion_check_table())
      df3.meta <- c(names(Expansion_check_table()))

      # df3.meta <- df3.meta[!grepl("RNA",df3.meta) & !grepl("BCR",df3.meta) & !grepl("TCR",df3.meta)& !grepl("_gene",df3.meta) & !grepl("allele",df3.meta) & !grepl("percent",df3.meta) & !grepl("cdr3",df3.meta)]


      if(length(df3.meta)>0) {
        selectInput("Colour_By_this_Expanded","Expanded Colour by: ",choices = df3.meta,selected="expand.singlets")
      }

      else {
        selectInput("Colour_By_this_Expanded","Expanded Colour by: ",choices = "expand.singlets",selected="expand.singlets")
      }


    })


    ## Expansion UMAP plot -----
    cols_UMAP_Expanded <- reactive({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      req(input$Colour_By_this_Expanded)

      top_BD_cluster <-  Expansion_check_table()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this_Expanded]

      if (input$Colour_By_this_Expanded == "expand.singlets") {
        top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function),NA,top_BD_cluster$Selected_function)
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = c("Ex","NEx"))
        top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function,decreasing = F),]
      }
      else {
        top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function),NA,top_BD_cluster$Selected_function)
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))

      }
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num)==T,])
      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$Colour_By_this_Expanded == "expand.singlets") {
        col.gg[1] <- "#00BFC4"
        col.gg[2] <- "grey90"

        palette_rainbow[2] <- "grey90"
        heat_col[2] <- "grey90"
        col.terrain[2] <- "grey90"
        col.topo[2] <- "grey90"
        col.hcl[2] <- "grey90"
      }

      else {

      }


      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep="_"), paste(num[i,]), col.gg[i])
        })
      }
      else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep="_"), paste(num[i,]), col.hcl[i])
        })
      }
      else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep="_"), paste(num[i,]), col.topo[i])
        })
      }
      else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep="_"), paste(num[i,]), heat_col[i])
        })
      }
      else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep="_"), paste(num[i,]), col.terrain[i])
        })
      }

      else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep="_"), paste(num[i,]), palette_rainbow[i])
        }) }
      else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep="_"), paste(num[i,]), palette1[i])
        })

      }
      else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep="_"), paste(num[i,]), input$one.colour.default)
        })
      } # one colour

    })

    output$cols_UMAP_Expanded <- renderUI({cols_UMAP_Expanded()})

    colors_Expanded <- reactive({
      sc <- input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      top_BD_cluster <-  Expansion_check_table()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this_Expanded]

      if (input$Colour_By_this_Expanded == "expand.singlets") {
        top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function),NA,top_BD_cluster$Selected_function)
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = c("NEx","Ex"))
      }
      else {
        top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function),NA,top_BD_cluster$Selected_function)
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))

      }

      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num)==T,])

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.UMAP_Expanded", i, sep="_")]]
      })
    })


    UMAP_Expanded_plot <- reactive({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      req(input$Colour_By_this_Expanded,input$wrap_row)
      top_BD_cluster <-  Expansion_check_table()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$Colour_By_this_Expanded]

      if (input$Colour_By_this_Expanded == "expand.singlets") {
        top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function),NA,top_BD_cluster$Selected_function)
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = c("Ex","NEx"))
        # Ex <- subset(top_BD_cluster,top_BD_cluster$Selected_function=="Ex")
        # NEx <- subset(top_BD_cluster,top_BD_cluster$Selected_function=="NEx")
      }
      else {
        top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function),NA,top_BD_cluster$Selected_function)
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))

      }

      col.file <- as.data.frame(unique(top_BD_cluster$Selected_function))
      col.file <- as.data.frame(col.file[complete.cases(col.file)==T,])
      names(col.file) <- "V1"

      col.file$col <- unlist(colors_Expanded())

      colorblind_vector <-as.data.frame(unlist(colors_Expanded()))

      if (dim(colorblind_vector)[1]==0) {
        num <- as.data.frame(unique(top_BD_cluster$Selected_function))

        if (input$colourtype == "default") {
          colorblind_vector <- c(gg_fill_hue(dim(num)[1]))
        } else if (input$colourtype == "hcl.colors") {
          colorblind_vector <- c(hcl.colors(dim(num)[1], palette = "viridis"))
        } else if (input$colourtype == "topo.colors") {
          colorblind_vector <- c(topo.colors(dim(num)[1]))
        } else if (input$colourtype == "heat.colors") {
          colorblind_vector <- c(heat.colors(dim(num)[1]))
        } else if (input$colourtype == "terrain.colors") {
          colorblind_vector <- c(terrain.colors(dim(num)[1]))
        } else if (input$colourtype == "rainbow") {
          colorblind_vector <- c(rainbow(dim(num)[1]))
        } else if (input$colourtype == "random") {
          colorblind_vector <- distinctColorPalette(dim(num)[1])

        }  else {

        }

      }

      if (input$Colour_By_this_Expanded == "expand.singlets") {
        colorblind_vector[2] <- "grey90"
      }


      colorblind_vector <- as.data.frame(colorblind_vector)

      names(colorblind_vector) <- "cols"

      names(top_BD_cluster)[names(top_BD_cluster) %in% input$Samp_col_expanded] <- "ID_Column"
      top_BD_cluster <- top_BD_cluster[top_BD_cluster$ID_Column %in% input$ID_Column_factor_expanded,]
      top_BD_cluster$ID_Column <- factor(top_BD_cluster$ID_Column,levels = input$ID_Column_factor_expanded)

      df <- ggplot(top_BD_cluster,aes(x=UMAP_1,UMAP_2,colour=Selected_function,alpha = Selected_function,label = Selected_function))

      # if (input$Colour_By_this_Expanded == "expand.singlets") {
      #
      #   df <- df + geom_point(data= NEx, aes(x=UMAP_1,UMAP_2,colour=Selected_function,alpha = Selected_function))
      #   df <- df + geom_point(data= Ex, aes(x=UMAP_1,UMAP_2,colour=Selected_function,alpha = Selected_function))
      # }

      # else {
      df <-   df + geom_point()
      # }


      df <- df + scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 20),values = colorblind_vector$col, na.value=input$NA_col_analysis)+
        scale_alpha_manual(labels = ~ stringr::str_wrap(.x, width = 20),values = rep(1,length(unique(top_BD_cluster$Selected_function))), na.value=0.1) +
        theme_bw() +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=0),
          axis.title.x = element_text(colour="black",angle=0,vjust=.5,face="plain",family=input$font_type,size = input$title.text.sizer2),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )

      if (input$Split_by_group=="no") {
        df <- df
      }
      else {
        df <- df + facet_wrap(~ID_Column, nrow = input$wrap_row)
      }
      df
    })

    output$UMAP_Expanded <- renderPlot({
      UMAP_Expanded_plot()
    })



    # download UMAP expanded ------
    output$downloadPlot_UMAP_Expanded  <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_Expanded_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_UMAP_Expanded ,height=input$height_UMAP_Expanded , onefile = FALSE) # open the pdf device
        plot(UMAP_Expanded_plot())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_UMAP_Expanded  <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_Expanded_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_UMAP_Expanded,
            height = input$height_png_UMAP_Expanded ,
            res = input$resolution_PNG_UMAP_Expanded )
        plot(UMAP_Expanded_plot())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )


    ## Expanded Violin plot ------

    compare.stat_Violin_sc <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      sc@meta.data$Samp <- sc@meta.data[,names(sc@meta.data) %in% input$Samp_col_expanded]
      sc@meta.data$Samp.filter <- ifelse(sc@meta.data$Samp %in% c(input$ID_Column_factor_expanded), sc@meta.data$Samp,"Abs")
      sc <- subset(sc, Samp.filter != "Abs")
      sc@meta.data$filter <- sc@meta.data[,names(sc@meta.data) %in% input$V_gene_sc]
      sc@meta.data$filter <- gsub("NA","Unknown",sc@meta.data$filter )
      sc <- subset(sc, filter != "Unknown")
      sc@meta.data$filter <- gsub("Unknown","NA",sc@meta.data$filter )
      req(input$Samp_col_expanded,input$ID_Column_factor_expanded,input$V_gene_sc,input$cut.off_expanded)
      meta.data <- sc@meta.data
      expanded.meta.data <-  Expansion_check_table()
      expanded.meta.data <- expanded.meta.data[,names(expanded.meta.data) %in% c(input$V_gene_sc,"Cell_Index",input$Samp_col_expanded,"expand.singlets","expansion.status")]
      expanded.meta.data <- expanded.meta.data %>%
        select("Cell_Index", everything())

      meta.data_exp <- merge(meta.data,expanded.meta.data,by = c(input$V_gene_sc,"Cell_Index",input$Samp_col_expanded),all.x=T,sort = F)
      meta.data_exp <- meta.data_exp[order(meta.data_exp$order,decreasing = F),]
      rownames(meta.data_exp) <- meta.data_exp$Cell_Index
      sc@meta.data <- meta.data_exp
      sc

    })
    # Vals_expanded.stats <- reactiveValues(output_ex1=NULL)


    output$Expansion_check <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      sc <- compare.stat_Violin_sc()
      sc@meta.data
    })

    Vals_expanded.stats <- reactive({

      sc <-input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )

      sc <- compare.stat_Violin_sc()
      req(input$selected_Indiv_Ex_1,input$selected_Indiv_Ex_2,input$pval.ex.filter)

      message("Updating Ident")
      Idents(object = sc) <- sc@meta.data$expansion.status

      min.pct.expression<- input$min_point_ #standard setting: 0.25
      min.logfc <-  input$LogFC_ #0.25 is standard
      message(paste0(" Calculating markers for cluster ",input$selected_Indiv_Ex_1," vs ",c(input$selected_Indiv_Ex_2)))

      markers.fm.list <- FindMarkers(sc, ident.1 = input$selected_Indiv_Ex_1, ident.2 = c(input$selected_Indiv_Ex_2), min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
      markers.fm.list2 <- subset(markers.fm.list,markers.fm.list$p_val_adj < input$pval.ex.filter)
      Vals_expanded.stats <- as.data.frame(markers.fm.list2)
      Vals_expanded.stats

    })

    output$compare.stat_Ex <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength =10, scrollX = TRUE),{

      sc <-input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )
      as.data.frame(Vals_expanded.stats())
    })

    output$downloadtb_compare.stat_Ex <- downloadHandler(
      filename = function(){
        paste("compare.stat_Expanded_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(Vals_expanded.stats())
        write.csv(df,file, row.names = T)
      } )
    ## dot plot
    relative_expression_plot_ex <- reactive({
      sc <- compare.stat_Violin_sc()
      validate(
        need(nrow(sc)>0,
             error_message_val_sc))

      if(input$restrict.dotpot == "yes") {
        list.names <- rownames(Vals_expanded.stats())[1:input$restrict.dotpot.num]
      }
      else {
        list.names <- rownames(Vals_expanded.stats())
      }

      size_legend = input$Bar_legend_size-2
      Idents(object = sc) <- sc@meta.data$expansion.status
      DotPlot(sc,features = list.names) +
        RotatedAxis() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size, angle = 90),
          axis.title.x = element_blank(),
          legend.title = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.text = element_text(colour="black", size=size_legend,family=input$font_type),
          legend.position = input$legend_position,
        ) +
        scale_colour_gradient2(low = input$low.dotplot.ex, mid = input$middle.dotplot.ex, high = input$high.dotplot.ex)+
        scale_x_discrete(labels = label_wrap(20)) +
        scale_y_discrete(labels = label_wrap(20))
    })

    output$relative_expression_dotplot_ex <- renderPlot({
      relative_expression_plot_ex()
    })


    output$downloadPlot_all_expression_dotplot_ex <- downloadHandler(
      filename = function() {
        paste("Expanded_dotplot","_",today(), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_all_expression_dotplot_ex,height=input$height_all_expression_dotplot_ex, onefile = FALSE) # open the pdf device
        df <- relative_expression_plot_ex()
        plot(df)
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_all_expression_dotplot_ex  <- downloadHandler(
      filename = function() {
        paste("Expanded_dotplot","_",today(), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_all_expression_dotplot_ex,
            height = input$height_png_all_expression_dotplot_ex,
            res = input$resolution_PNG_all_expression_dotplot_ex)
        df <- relative_expression_plot_ex()

        plot(df)
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    # expression for one gene -------

    Violin_chart_Expanded <- reactive({
      df <- Ridge_chart_alpha_gamma_df()

      df2$Selected_function <- factor(df2$Selected_function,levels = input$Graph_split_order)
      ggplot(df2, aes(y = get(input$string.data_Exp_top), x = Selected_function, fill = Selected_function)) +
        geom_violin() +
        geom_jitter(height = 0, width = 0.1)+
        theme(legend.position = "none",

        )+
        theme_bw() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=90),
          axis.title.x = element_blank(),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_text(colour="black",size=20,family=input$font_type),
          legend.position = "none",
        )+
        ggtitle(input$string.data_Exp_top)

    })




    # Over representation analysis for Expanded  -----

    Over_rep_Exp_old <- reactive({
      sc <-input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )

      req(input$datasource,input$in.geneset.cutoff_Exp)

      geneSet <- read.csv(system.file("OverRep","GeneSets.csv",package = "STEGO.R"),header = T)

      background.genes.name <- as.data.frame(rownames(sc@assays$RNA$scale.data))
      names(background.genes.name) <- "V1"
      background.genes <- length(rownames(sc@assays$RNA$scale.data))


      geneSet$background.genes <- background.genes
      markers.fm.list <- Vals_expanded.stats() # needs to match the Exp stat
      DEx.genes <- as.data.frame(rownames(markers.fm.list))
      names(DEx.genes) <- "V1"
      total.sig <- length(DEx.genes$V1)
      geneSet$total.sig <- length(DEx.genes$V1)

      geneSet$background.geneset <- NA
      geneSet$background.geneset.name <- NA
      geneSet$in.geneset <- NA
      geneSet$in.geneset.name <- NA

      if(input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        geneSet$GeneSet <- gsub("-",".",geneSet$GeneSet)
      }

      if(input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        require(stringr)
        geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
      }

      for (i in 1:dim(geneSet)[1]) {
        # listed GeneSet
        message(paste("GeneSet: ", i))
        Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet,";")[i])
        names(Gene.set.testing) <- "V1"
        Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
        names(Gene.set.testing2) <- "V1"
        background.overlap <- merge(Gene.set.testing2,background.genes.name,by= "V1")
        # message(paste(dim(background.overlap)[1],"in Background"))
        geneSet$background.geneset[i] <- length(background.overlap$V1)
        geneSet$background.geneset.name[i] <- as.character(paste(unlist(background.overlap[1]), collapse=';'))
        # in sig gene list
        overlap <- merge(background.overlap,DEx.genes,by= "V1")
        # message(paste(dim(overlap)[1],"# Sig genes"))
        geneSet$in.geneset[i] <- length(overlap$V1)
        geneSet$in.geneset.name[i] <- as.character(paste(unlist(overlap[1]), collapse=';'))

      }

      geneSet2 <- subset(geneSet,geneSet$in.geneset>0)

      for (i in 1:dim(geneSet2)[1]) {
        tota.gene.set <- geneSet2$background.geneset[i] # genes that are identified in background
        tota.gene.set
        in.geneset <-  geneSet2$in.geneset[i]# DEx in geneset

        background.genes
        not.in.total <- background.genes - tota.gene.set
        not.in.geneset.sig <- total.sig - in.geneset
        d <- data.frame( gene.in.interest=c( in.geneset, not.in.geneset.sig),gene.not.interest=c( tota.gene.set, not.in.total))
        row.names(d) <- c("In_category", "not_in_category")

        if (in.geneset>0) {
          geneSet2$p.val[i] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
          geneSet2$lowerCI[i] <-  unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
          geneSet2$upperCI[i] <-unlist(fisher.test(d)$conf.int)[2]
          geneSet2$OR[i] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1],3)
        }

        else {
          geneSet2$p.value[i] <- "-"
          geneSet2$lowerCI[i] <-  "-"
          geneSet2$upperCI[i] <- "-"
          geneSet2$OR[i] <- "-"
        }
        # message(print(d))
        # message(print(round(prop.table(d),3)))
      }

      geneSet2
    })

    Over_rep_Exp <- reactive({
      geneSet2 <-Over_rep_Exp_old()

      validate(
        need(nrow(geneSet2)>0,
             error_message_val_sc)
      )

      req(input$datasource,input$in.geneset.cutoff_Exp)

      geneSet2 <- geneSet2[order(geneSet2$p.val,decreasing = F),]
      geneSet2 <- subset(geneSet2,geneSet2$in.geneset>=input$in.geneset.cutoff_Exp)

      # geneSet2 <- subset(geneSet2,geneSet2$pval.BH.adj<=input$adjust_cutoff_Exp)

      geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
      geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
      geneSet2 <- subset(geneSet2,geneSet2$p.val<=input$p.val_cutoff_Exp)
      # name.list <- c("Geneset_ID","p.val","FDR","Bonferroni","OR","lowerCI","upperCI","in.geneset.name","in.geneset","background.geneset","total.sig","background.genes","background.geneset.name")
      # geneSet2 <- geneSet2 %>%
      #   select(all_of(name.list), everything())
      geneSet2
    })
    output$Over_rep_Exp_Tab <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1,2,5,10,20,50,100), pageLength = 20, scrollX = TRUE),{
      Over_rep_Exp()
    })

    output$downloadtb_over.rep_Exp <- downloadHandler(
      filename = function(){

        paste0(input$Selected_clonotype,"_",today(),"_over_rep.csv")
      },
      content = function(file){
        df <- as.data.frame(Over_rep_Exp())
        write.csv(df,file, row.names = F)
      } )
    #### Epitope upload -----
    df_tcrex <- reactive({
      epi <- input.data_sc_TCRex()
      epi[!(duplicated(epi$CDR3_beta)),]
    })

    output$MainTcell_Check <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{

      df_tcrex()

    })
    #### Epitope heatmap -----
    heatmap_epitope <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )

      df3.meta <- sc@meta.data

      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

      }
      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }
      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta,epi,by="CDR3_beta")
      df3.meta$selected <- df3.meta[,names(df3.meta) %in% input$epitope_umap_selected]
      df3.meta <- df3.meta[order(df3.meta$selected,decreasing = F),]
      df3.meta$selected <- factor(df3.meta$selected,levels = unique(df3.meta$selected))
      df3.meta$Selected_group <- df3.meta[,names(df3.meta) %in% input$epitope_umap_selected2]
      df3.meta$cloneCount <- 1

      # df <- as.data.frame(ddply(df3.meta,(c("selected","Selected_group","cloneCount")),numcolwise(sum)))

      df.1 <- acast(df3.meta, selected~Selected_group, value.var="cloneCount")
      df.1
      head(df.1)
      df.1[is.na(df.1)] <- 0
      dim(df.1)
      min.FC <- min(df.1)
      med.FC <- max(df.1)/2
      max.FC <- max(df.1)

      ht <- Heatmap(df.1,
                    heatmap_legend_param = list(title = "Unique\nTCR's",
                                                title_gp = gpar(fontsize = 10,
                                                                fontface = "bold",fontfamily='serif'),
                                                labels_gp = gpar(fontsize = 10,fontfamily='serif')),
                    col = colorRamp2(c(min.FC,med.FC,max.FC), c("white","gold","purple")),
                    row_names_gp = grid::gpar(fontsize = 10,fontfamily='serif'),
                    column_names_gp = grid::gpar(fontsize = 10,fontfamily='serif'),

      )

      draw(ht, padding = unit(c(10, 10, 10, 10), "mm"))
    })
    output$Heatmap_epi_plot <- renderPlot({
      heatmap_epitope()

    })

    output$downloadPlot_Heatmap_epi_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_",input$epitope_umap_selected,"_",input$epitope_umap_selected2,"_Heatmap_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Heatmap_epi_plot ,height=input$height_Heatmap_epi_plot , onefile = FALSE) # open the pdf device
        plot(heatmap_epitope())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Heatmap_epi_plot  <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_",input$epitope_umap_selected,"_",input$epitope_umap_selected2,"_Heatmap_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Heatmap_epi_plot,
            height = input$height_png_Heatmap_epi_plot,
            res = input$resolution_PNG_Heatmap_epi_plot)
        plot(heatmap_epitope())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    #### epitope UMAP interrogation -----
    cols_epitope <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )
      df3.meta <- sc@meta.data

      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

      }

      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }

      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta,epi,by="CDR3_beta",all.x=T)
      df3.meta$selected <- df3.meta[,names(df3.meta) %in% input$epitope_umap_selected]
      df3.meta <- df3.meta[order(df3.meta$selected,decreasing = F),]
      df3.meta$selected <- factor(df3.meta$selected,levels = unique(df3.meta$selected))

      num <- as.data.frame(unique(df3.meta$selected))
      num <- as.data.frame(num[complete.cases(num)==T,])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep="_"), paste(num[i,]), col.gg[i])
        })
      }
      else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep="_"), paste(num[i,]), col.hcl[i])
        })
      }
      else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep="_"), paste(num[i,]), col.topo[i])
        })
      }
      else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep="_"), paste(num[i,]), heat_col[i])
        })
      }
      else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep="_"), paste(num[i,]), col.terrain[i])
        })
      }

      else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep="_"), paste(num[i,]), palette_rainbow[i])
        }) }
      else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep="_"), paste(num[i,]), palette1[i])
        })

      }
      else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep="_"), paste(num[i,]), input$one.colour.default)
        })


      } # one colour

    })

    output$myPanel_cols_epitope <- renderUI({cols_epitope()})
    colors_cols_epitope <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )
      df3.meta <- sc@meta.data

      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

      }

      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }

      else {

        df3.meta$CDR3_beta <- df3.meta$cdr3_BD

      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta,epi,by="CDR3_beta",all.x=T)
      df3.meta$selected <- df3.meta[,names(df3.meta) %in% input$epitope_umap_selected]


      df3.meta <- df3.meta[order(df3.meta$selected,decreasing = F),]
      df3.meta$selected <- factor(df3.meta$selected,levels = unique(df3.meta$selected))

      num <- as.data.frame(unique(df3.meta$selected))
      num <- as.data.frame(num[complete.cases(num)==T,])

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.cols_epitope", i, sep="_")]]
      })
    })

    UMAP_Epitope <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )
      df3.meta <- sc@meta.data

      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")
      }

      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }

      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta,epi,by="CDR3_beta")
      df3.meta$selected <- df3.meta[,names(df3.meta) %in% input$epitope_umap_selected]
      df3.meta <- df3.meta[order(df3.meta$selected,decreasing = F),]
      df3.meta$selected <- factor(df3.meta$selected,levels = unique(df3.meta$selected))


      names(df3.meta)[names(df3.meta) %in% input$Split_group_by_] <- "ID_Column"
      # df3.meta$ID_Column<- df3.meta[df3.meta$ID_Column %in% input$ID_Column_factor,]
      df3.meta$ID_Column <- factor(df3.meta$ID_Column,levels = input$Graph_split_order)

      num <- as.data.frame(unique(df3.meta$selected))
      num <- as.data.frame(num[complete.cases(num)==T,])
      palette_rainbow <- unlist(colors_cols_epitope())

      df <- ggplot(data=df3.meta,aes(x=UMAP_1,UMAP_2,colour=selected,size= selected,alpha = selected))+
        geom_point()+
        scale_color_manual(na.value=input$NA_col_analysis, values = palette_rainbow)+
        scale_size_manual(na.value=0.25,values = rep(input$value_size_epi_umap,dim(num)[1]))+
        scale_alpha_manual(na.value=0.25,values = rep(1,dim(num)[1]))+
        theme(
          legend.text = element_text(colour="black", size=12,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )+
        theme_bw()
      if (input$Split_by_group=="no") {
        df <- df
      }
      else {
        df <- df + facet_wrap(~ID_Column, nrow = input$wrap_row)
      }
      df
    })
    output$UMAP_Epitope_plot <- renderPlot({
      UMAP_Epitope()
    })

    output$downloadPlot_UMAP_Epitope <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_",input$epitope_umap_selected,"_UMAP_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_UMAP_Epitope,height=input$height_UMAP_Epitope, onefile = FALSE) # open the pdf device
        plot(UMAP_Epitope())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_UMAP_Epitope <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_",input$epitope_umap_selected,"_UMAP_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_UMAP_Epitope,
            height = input$height_png_UMAP_Epitope,
            res = input$resolution_PNG_UMAP_Epitope)
        plot(UMAP_Epitope())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    # Epitope pie chart function -----

    cols_epitope_pie <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )
      df3.meta <- sc@meta.data
      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

      }
      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }
      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta,epi,by="CDR3_beta")
      df3.meta$selected <- df3.meta[,names(df3.meta) %in% input$epitope_umap_selected]


      df3.meta <- df3.meta[order(df3.meta$selected,decreasing = F),]
      df3.meta$selected <- factor(df3.meta$selected,levels = unique(df3.meta$selected))

      num <- as.data.frame(unique(df3.meta$selected))
      num <- as.data.frame(num[complete.cases(num)==T,])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep="_"), paste(num[i,]), col.gg[i])
        })
      }
      else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep="_"), paste(num[i,]), col.hcl[i])
        })
      }
      else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep="_"), paste(num[i,]), col.topo[i])
        })
      }
      else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep="_"), paste(num[i,]), heat_col[i])
        })
      }
      else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep="_"), paste(num[i,]), col.terrain[i])
        })
      }

      else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep="_"), paste(num[i,]), palette_rainbow[i])
        }) }
      else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep="_"), paste(num[i,]), palette1[i])
        })

      }
      else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep="_"), paste(num[i,]), input$one.colour.default)
        })


      } # one colour

    })

    output$myPanel_cols_epitope_pie <- renderUI({cols_epitope_pie()})

    colors_cols_epitope_pie <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )
      df3.meta <- sc@meta.data

      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

      }
      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }
      else {

        df3.meta$CDR3_beta <- df3.meta$cdr3_BD

      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta,epi,by="CDR3_beta")
      df3.meta$selected <- df3.meta[,names(df3.meta) %in% input$epitope_umap_selected]

      df3.meta <- df3.meta[order(df3.meta$selected,decreasing = F),]
      df3.meta$selected <- factor(df3.meta$selected,levels = unique(df3.meta$selected))

      num <- as.data.frame(unique(df3.meta$selected))
      num <- as.data.frame(num[complete.cases(num)==T,])

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.cols_epitope_pie", i, sep="_")]]
      })
    })

    Pie_chart_Epitope <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )
      df3.meta <- sc@meta.data
      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")
      }
      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }
      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta,epi,by="CDR3_beta")
      df3.meta$Selected_function <- df3.meta[,names(df3.meta) %in% input$epitope_umap_selected]
      df3.meta <- df3.meta[order(df3.meta$Selected_function,decreasing = F),]
      df3.meta$Selected_function <- factor(df3.meta$Selected_function,levels = unique(df3.meta$Selected_function))

      top_BD_cluster <-  df3.meta
      top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% input$epitope_umap_selected2]
      top_BD_cluster$Selected_indiv <- top_BD_cluster[,names(top_BD_cluster) %in% input$Samp_col]

      if (input$by_indiv_pie_epi == "yes") {
        top_BD_cluster <- top_BD_cluster[top_BD_cluster$Selected_indiv %in% input$selected_Indiv,]
      }

      df.col <- unlist(colors_cols_epitope_pie())

      df3.meta3 <-  as.data.frame(table(top_BD_cluster$Selected_group,top_BD_cluster$Selected_function))
      total.condition <- as.data.frame(ddply(df3.meta3,"Var1",numcolwise(sum)))
      dim(total.condition)[1]
      dim(df3.meta3)[1]
      emtpy <- matrix(nrow =dim(df3.meta3)[1],ncol=dim(total.condition)[1])

      for (i in 1:dim(df3.meta3)[1]) {

        emtpy[i,] <- ifelse(df3.meta3$Var1[i]==total.condition$Var1[1:dim(total.condition)[1]],
                            total.condition[total.condition$Var1==total.condition$Var1[1:dim(total.condition)[1]],2],F)
      }
      df3.meta3$n <- df3.meta3$Freq/rowSums(emtpy)



      ggplot(df3.meta3,aes(x="", y=n, fill=Var2, group = Var1)) +
        geom_bar(stat="identity", width=1)+
        coord_polar("y", start=0)  +
        theme_void(20) +
        facet_wrap(~Var1, nrow = input$wrap_row) +
        theme(
          legend.key.size = unit(1, 'cm'),
          legend.title = element_blank()) +
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 20),values = df.col, na.value = input$NA_col_analysis) +
        theme(strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
              legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.position = input$legend_position,
              legend.title = element_blank()
        )

    })

    output$Pie_Epitope_plot <- renderPlot({
      Pie_chart_Epitope()
    })

    Pie_Epitope_dt_process <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )
      df3.meta <- sc@meta.data

      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")
      }

      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }

      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }

      df3.meta
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta,epi,by="CDR3_beta")
      df3.meta
      df3.meta$Selected_function <- df3.meta[,names(df3.meta) %in% input$epitope_umap_selected]
      top_BD_cluster <-  df3.meta
      top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% input$epitope_umap_selected2]
      top_BD_cluster$Selected_indiv <- top_BD_cluster[,names(top_BD_cluster) %in% input$Samp_col]
      top_BD_cluster$cloneCount <- 1
      top_BD_cluster2 <-  top_BD_cluster[,names(top_BD_cluster) %in% c("cloneCount","Selected_function","Selected_group","Selected_indiv")]
      top_BD_cluster2 <- top_BD_cluster2 %>%
        select(cloneCount, everything())
      df2 <- as.data.frame(ddply(top_BD_cluster2,names(top_BD_cluster2)[-c(1)],numcolwise(sum)))
      df2$fraction <- df2$cloneCount/sum(df2$cloneCount)
      df2$Percent <- round(df2$cloneCount/sum(df2$cloneCount)*100,2)
      df2$Selected_function <- ifelse(grepl("NA", df2$Selected_function),"-",df2$Selected_function)
      df2[order(df2$Percent,decreasing = T),]
    })

    output$Pie_Epitope_dt <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{

      Pie_Epitope_dt_process()

    })




    # Epitope stat -----
    Epitope_of_interest <- reactive ({
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )

      df3.meta <- sc@meta.data
      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

      }
      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }
      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }

      checking1 <- merge(df3.meta,epi,by="CDR3_beta")

      checking2 <- checking1[,names(checking1) %in% c("epitope","pathology")]
      names.check <- names(checking2)
      checking2$cloneCount <- 1
      sum.check2 <- ddply(checking2,names.check ,numcolwise(sum))
      sum.check2_Morethan2 <- subset(sum.check2,sum.check2$cloneCount>1)
      sum.check2_Morethan2 <- sum.check2_Morethan2[order(sum.check2_Morethan2$cloneCount,decreasing = T),]
      sum.check2_Morethan2

    })

    output$Epi_of_interest_DF <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )
      as.data.frame(Epitope_of_interest())

    })

    observeEvent(input$Update_epi,{
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )

      sum.check2_Morethan2 <- as.data.frame(Epitope_of_interest())

      updateSelectInput(
        session,
        "Epi_of_interest",
        choices=sum.check2_Morethan2$epitope,
        selected=sum.check2_Morethan2$epitope[1]
      )
    })

    compare.stat_Epi <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )
      req(input$Epi_of_interest)
      df3.meta <- sc@meta.data
      as.data.frame(Epitope_of_interest())
      epi <- epi[epi$epitope %in% input$Epi_of_interest,]
      epi <- epi[names(epi) %in% c("CDR3_beta","epitope","pathology")]
      epi <- unique(epi)

      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

      }
      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }
      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      df3.meta
      checking <- merge(df3.meta,epi,by="CDR3_beta",all.x=T)
      rownames(checking) <- checking$Cell_Index
      name.check.epi <- unlist(unique(epi$epitope))
      checking$epi_selected <- ifelse(as.character(checking$epitope) == name.check.epi,name.check.epi,"not selected")
      rownames(checking) <- (checking$Cell_Index)
      checking <- checking[order(checking$order,decreasing = F),]
      checking
      sc@meta.data <- checking
      Idents(object = sc) <- sc@meta.data$epi_selected

      min.pct.expression<- input$min_point_ #standard setting: 0.25
      min.logfc<-  input$LogFC_ #0.25 is standard

      cluster.names <- unique(Idents(sc))[order(unique(Idents(sc)))]

      markers.fm.list <- FindMarkers(sc, ident.1 = name.check.epi, min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
      markers.fm.list2 <- subset(markers.fm.list,markers.fm.list$p_val_adj < input$pval.ex.filter)
      as.data.frame(markers.fm.list2)

    })

    output$compare.stat_Epi_DT <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength =10, scrollX = TRUE),{
      sc <-input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )
      compare.stat_Epi()
    })

    #Epitope dotplot ----
    all_expression_plot_epi <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )
      req(input$Epi_of_interest)
      df3.meta <- sc@meta.data
      as.data.frame(Epitope_of_interest())
      epi <- epi[epi$epitope %in% input$Epi_of_interest,]
      epi <- epi[names(epi) %in% c("CDR3_beta","epitope","pathology")]
      epi <- unique(epi)

      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

      }
      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }
      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      df3.meta
      checking <- merge(df3.meta,epi,by="CDR3_beta",all.x=T)

      checking

      rownames(checking) <- checking$Cell_Index
      name.check.epi <- as.character(unlist(unique(epi$epitope)))
      checking$epitope <- as.character(checking$epitope)

      checking$epi_selected <- ifelse(as.character(checking$epitope) != name.check.epi,"NS",name.check.epi)
      checking$epi_selected[is.na(checking$epi_selected)] <- "NS"

      checking <- checking[order(checking$order,decreasing = F),]
      checking
      sc@meta.data <- checking
      Idents(object = sc) <- sc@meta.data$epi_selected

      if (input$restrict.dotpot.epi == F ) {
        list.names <- rownames(compare.stat_Epi())
      }

      else {
        list.names <- rownames(compare.stat_Epi())
        list.names <- list.names[1:input$restrict.dotpot.num.epi]
      }

      size_legend = input$Bar_legend_size-2

      DotPlot(sc, features = list.names) +
        RotatedAxis() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size, angle = 90),
          axis.title.x = element_blank(),
          legend.title = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.text = element_text(colour="black", size=size_legend,family=input$font_type),
          legend.position = input$legend_position,
        ) +
        scale_colour_gradient2(low = input$low.dotplot.epi, mid = input$middle.dotplot.epi, high = input$high.dotplot.epi)+
        scale_x_discrete(labels = label_wrap(20)) +
        scale_y_discrete(labels = label_wrap(20))
    })

    output$checking_epi_dot_issue <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1,2,5,10,20,50,100), pageLength = 20, scrollX = TRUE),{
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )
      req(input$Epi_of_interest)
      df3.meta <- sc@meta.data
      as.data.frame(Epitope_of_interest())
      epi <- epi[epi$epitope %in% input$Epi_of_interest,]
      epi <- epi[names(epi) %in% c("CDR3_beta","epitope","pathology")]
      epi <- unique(epi)

      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

      }
      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }
      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      df3.meta
      checking <- merge(df3.meta,epi,by="CDR3_beta",all.x=T)

      checking

      rownames(checking) <- checking$Cell_Index
      name.check.epi <- as.character(unlist(unique(epi$epitope)))
      checking$epitope <- as.character(checking$epitope)

      checking$epi_selected <- ifelse(as.character(checking$epitope) != name.check.epi,"NS",name.check.epi)
      checking$epi_selected[is.na(checking$epi_selected)] <- "NS"

      checking <- checking[order(checking$order,decreasing = F),]
      checking
      sc@meta.data <- checking
      Idents(object = sc) <- sc@meta.data$epi_selected
      as.data.frame(unique(Idents(object = sc)))
      # if (input$restict_no_points == F ) {
      #   list.names <- rownames(compare.stat_Epi())
      # }
      #
      # else {
      #   list.names <- rownames(compare.stat_Epi())
      #   list.names <- list.names[1:input$pval.ex.top_genes]
      # }


    })

    # add in percentage - default to 10%

    output$all_expression_dotplot_epi <- renderPlot({
      all_expression_plot_epi()
    })

    output$downloadPlot_all_expression_dotplot_epi <- downloadHandler(
      filename = function() {
        paste(input$Epi_of_interest,"_dotplot","_",today(), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_all_expression_dotplot_epi,height=input$height_all_expression_dotplot_epi, onefile = FALSE) # open the pdf device
        df <- all_expression_plot_epi()
        plot(df)
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_all_expression_dotplot_epi  <- downloadHandler(
      filename = function() {
        paste(input$Epi_of_interest,"_dotplot","_",today(), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_all_expression_dotplot_epi,
            height = input$height_png_all_expression_dotplot_epi,
            res = input$resolution_PNG_all_expression_dotplot_epi)
        df <- all_expression_plot_ex()

        plot(df)
        dev.off()},   contentType = "application/png" # MIME type of the image
    )
    # Over representation analysis for Epitope  -----
    Over_rep_Epi_old <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(sc)>0 & nrow(epi)>0,
             "Upload Files")
      )
      req(input$Epi_of_interest)
      df3.meta <- sc@meta.data
      as.data.frame(Epitope_of_interest())
      epi <- epi[epi$epitope %in% input$Epi_of_interest,]
      epi <- epi[names(epi) %in% c("CDR3_beta","epitope","pathology")]
      epi <- unique(epi)

      if(input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

      }
      else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      }
      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      df3.meta
      checking <- merge(df3.meta,epi,by="CDR3_beta",all.x=T)

      checking

      rownames(checking) <- checking$Cell_Index
      name.check.epi <- as.character(unlist(unique(epi$epitope)))
      checking$epitope <- as.character(checking$epitope)

      checking$epi_selected <- ifelse(as.character(checking$epitope) != name.check.epi,"NS",name.check.epi)
      checking$epi_selected[is.na(checking$epi_selected)] <- "NS"

      checking <- checking[order(checking$order,decreasing = F),]
      sc@meta.data <- checking
      Idents(object = sc) <- sc@meta.data$epi_selected

      geneSet <- read.csv(system.file("OverRep","GeneSets.csv",package = "STEGO.R"),header = T)

      background.genes.name <- as.data.frame(rownames(sc@assays$RNA$scale.data))
      names(background.genes.name) <- "V1"
      background.genes <- length(rownames(sc@assays$RNA$scale.data))

      geneSet$background.genes <- background.genes

      markers.fm.list <- compare.stat_Epi() # need to change to epi

      DEx.genes <- as.data.frame(rownames(markers.fm.list))
      names(DEx.genes) <- "V1"
      total.sig <- length(DEx.genes$V1)
      geneSet$total.sig <- length(DEx.genes$V1)

      geneSet$background.geneset <- NA
      geneSet$background.geneset.name <- NA
      geneSet$in.geneset <- NA
      geneSet$in.geneset.name <- NA

      if(input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        geneSet$GeneSet <- gsub("-",".",geneSet$GeneSet)
      }

      if(input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        require(stringr)
        geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
      }

      for (i in 1:dim(geneSet)[1]) {
        # listed GeneSet
        message(paste("GeneSet: ", i))
        Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet,";")[i])
        names(Gene.set.testing) <- "V1"
        Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
        names(Gene.set.testing2) <- "V1"
        # message(paste(dim(Gene.set.testing2)[1],"GeneSet total"))
        # background genes
        background.overlap <- merge(Gene.set.testing2,background.genes.name,by= "V1")
        # message(paste(dim(background.overlap)[1],"in Background"))
        geneSet$background.geneset[i] <- length(background.overlap$V1)
        geneSet$background.geneset.name[i] <- as.character(paste(unlist(background.overlap[1]), collapse=';'))
        # in sig gene list
        overlap <- merge(background.overlap,DEx.genes,by= "V1")
        # message(paste(dim(overlap)[1],"# Sig genes"))
        geneSet$in.geneset[i] <- length(overlap$V1)
        geneSet$in.geneset.name[i] <- as.character(paste(unlist(overlap[1]), collapse=';'))

      }

      geneSet2 <- subset(geneSet,geneSet$in.geneset>0)

      for (i in 1:dim(geneSet2)[1]) {
        tota.gene.set <- geneSet2$background.geneset[i] # genes that are identified in background
        tota.gene.set
        in.geneset <-  geneSet2$in.geneset[i]# DEx in geneset

        background.genes
        not.in.total <- background.genes - tota.gene.set
        not.in.geneset.sig <- total.sig - in.geneset
        d <- data.frame( gene.in.interest=c( in.geneset, not.in.geneset.sig),gene.not.interest=c( tota.gene.set, not.in.total))
        row.names(d) <- c("In_category", "not_in_category")

        if (in.geneset>0) {
          geneSet2$p.val[i] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
          geneSet2$lowerCI[i] <-  unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
          geneSet2$upperCI[i] <-unlist(fisher.test(d)$conf.int)[2]
          geneSet2$OR[i] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1],3)
        }

        else {
          geneSet2$p.value[i] <- "-"
          geneSet2$lowerCI[i] <-  "-"
          geneSet2$upperCI[i] <- "-"
          geneSet2$OR[i] <- "-"
        }
        # message(print(d))
        # message(print(round(prop.table(d),3)))
      }
      geneSet2
    })
    Over_rep_Epi <- reactive({
      geneSet2 <- Over_rep_Epi_old()

      validate(
        need(nrow(geneSet2)>0,
             "Upload Files")
      )
      geneSet2 <- geneSet2[order(geneSet2$p.val,decreasing = F),]
      geneSet2 <- subset(geneSet2,geneSet2$in.geneset>=input$in.geneset.cutoff_Epi)
      geneSet2 <- subset(geneSet2,geneSet2$p.val<=input$p.val_cutoff_Epi)

      geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
      geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
      name.list <- c("Geneset_ID","p.val","FDR","Bonferroni","OR","lowerCI","upperCI","in.geneset.name","in.geneset","background.geneset","total.sig","background.genes","background.geneset.name")
      geneSet2 <- geneSet2 %>%
        select(all_of(name.list), everything())
    })

    output$Over_rep_Epi_Tab <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1,2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      Over_rep_Epi()
    })

    output$downloadtb_over.rep.Epi <- downloadHandler(
      filename = function(){

        paste0("Epitope","_",today(),"_over_rep.csv")
      },
      content = function(file){
        df <- as.data.frame(Over_rep_Epi())
        write.csv(df,file, row.names = F)
      } )


    #### download Epitope files ------
    output$downloadPlot_Pie_Epitope <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$epitope_Pie_Epitope,"_Pie_epi_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Pie_Epitope,height=input$height_Pie_Epitope, onefile = FALSE) # open the pdf device
        plot(Pie_chart_Epitope())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Pie_Epitope <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$epitope_umap_selected,"_Pie_epi_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Pie_Epitope,
            height = input$height_png_Pie_Epitope,
            res = input$resolution_PNG_Pie_Epitope)
        plot(Pie_chart_Epitope())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )


    output$downloadPlot_UMAP_Epitope <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$epitope_umap_selected,"_UMAP_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_UMAP_Epitope,height=input$height_UMAP_Epitope, onefile = FALSE) # open the pdf device
        plot(UMAP_Epitope())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_UMAP_Epitope <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$epitope_umap_selected,"_UMAP_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_UMAP_Epitope,
            height = input$height_png_UMAP_Epitope,
            res = input$resolution_PNG_UMAP_Epitope)
        plot(UMAP_Epitope())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    output$downloaddf_Pie_Epitope_dt <- downloadHandler(
      filename = function(){
        paste("Summary_table_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(Pie_Epitope_dt_process())
        write.csv(df,file, row.names = F)
      } )



    #### clusTCR2 seting up the dataframe -----
    AG_cluster <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )
      md <- sc@meta.data
      x = today()

      if (length(input.data_sc_clusTCR_AG())>0) {
        clust <- input.data_sc_clusTCR_AG()
        req(clust,md)
        if(input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") {
          names(md)[names(md) %in% "v_gene_AG"] <- "Selected_V_AG"
          names(md)[names(md) %in% "junction_aa_AG"] <- "AminoAcid_AG"
        } else {
          names(md)[names(md) %in% "v_gene_AG"] <- "Selected_V_AG"
          names(md)[names(md) %in% "cdr3_AG"] <- "AminoAcid_AG"
        }

        req(md$Selected_V_AG,md$AminoAcid_AG)

        md$CDR3_Vgene <- paste(md$AminoAcid_AG,md$Selected_V_AG,sep="_")
        df <- merge(md,clust,by = "CDR3_Vgene")

        # Total clone Count -----
        df2 <- as.data.frame(df$Clust_size_order)
        names(df2) <- "Clust_size_order"
        df2$Total_cloneCount <- 1
        df3 <- as.data.frame(ddply(df2,"Clust_size_order",numcolwise(sum)))
        df3 <- subset(df3,df3$Total_cloneCount>2)
        df4 <- merge(df3,df,by = "Clust_size_order")
        df4
        # updated cluster count
        clusterCount <- df4[,names(df4) %in% c("Clust_size_order","CDR3_Vgene")]
        clusterCount <- clusterCount[!duplicated(clusterCount), ]
        clusterCount$ClusterCount <- 1
        clusterCount2 <- as.data.frame(ddply(clusterCount,c("Clust_size_order"),numcolwise(sum)))
        clusterCount2 <- subset(clusterCount2,clusterCount2$ClusterCount>1)
        df4_clusterCount2 <- merge(df4,clusterCount2,by = c("Clust_size_order"))
        df4_clusterCount2
        # Sample Count
        SampCount <- df4_clusterCount2[,names(df4_clusterCount2) %in% c("Clust_size_order",input$Samp_col)]
        SampCount <- SampCount[!duplicated(SampCount), ]
        SampCount$Sample_count <- 1
        df9 <- as.data.frame(ddply(SampCount,c("Clust_size_order"),numcolwise(sum)))
        df10 <- merge(df9,df4_clusterCount2,by = "Clust_size_order")

        # Calculating the priority
        df10$priority <- 1/(df10$Total_cloneCount * df10$ClusterCount * df10$Sample_count)
        df10 <- df10[order(df10$priority,decreasing = F),]
        df10
        df10$priority[is.na(df10$priority)] <- 0
        df10 <- subset(df10,df10$priority>0)
        df10
        # updated order
        updatedOrder <- df10[,names(df10) %in% c("priority","Clust_size_order")]
        updatedOrder <- updatedOrder[!duplicated(updatedOrder), ]
        updatedOrder <- updatedOrder[order(updatedOrder$priority,decreasing = F),]
        updatedOrder$Updated_order <- 1:dim(updatedOrder)[1]
        updatedOrder

        # final data frame
        df7 <- merge(updatedOrder, df10, by = c("Clust_size_order","priority"))
        df7 <- df7[order(df7$priority,decreasing = F),]


        clusterAG <- df7 %>%
          select(all_of(c(input$Samp_col,"Sample_count","Total_cloneCount","ClusterCount","priority","Updated_order")), everything())
      } else {

      }
    })
    BD_cluster <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )
      md <- sc@meta.data
      x = today()
      # req(input$priority_cutoffBD)
      if (length(input.data_sc_clusTCR_BD())>0) {
        sc <- UMAP_metadata_with_labs()
        validate(
          need(nrow(sc)>0,
               "Upload File")
        )
        md <- sc@meta.data
        x = today()
        # req(input$priority_cutoffBD)

        clust <- input.data_sc_clusTCR_BD()
        clust
        if(input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") {
          names(md)[names(md) %in% "v_gene_BD"] <- "Selected_V_BD"
          names(md)[names(md) %in% "junction_aa_BD"] <- "AminoAcid_BD"
        } else {
          names(md)[names(md) %in% "v_gene_BD"] <- "Selected_V_BD"
          names(md)[names(md) %in% "cdr3_BD"] <- "AminoAcid_BD"
        }
        md <- md[,names(md) %in% c("Selected_V_BD","AminoAcid_BD","UMAP_1","UMAP_2","Cell_Index","Sample_Name","orig.ident",input$Split_group_by_,input$Colour_By_this,input$Samp_col)]
        req(md$Selected_V_BD,md$AminoAcid_BD)
        md$CDR3_Vgene <- paste(md$AminoAcid_BD,md$Selected_V_BD,sep="_")
        df <- merge(md,clust,by = "CDR3_Vgene")

        # Total clone Count -----
        df2 <- as.data.frame(df$Clust_size_order)
        names(df2) <- "Clust_size_order"
        df2$Total_cloneCount <- 1
        df3 <- as.data.frame(ddply(df2,"Clust_size_order",numcolwise(sum)))
        df3 <- subset(df3,df3$Total_cloneCount>2)
        df4 <- merge(df3,df,by = "Clust_size_order")
        df4
        # updated cluster count
        clusterCount <- df4[,names(df4) %in% c("Clust_size_order","CDR3_Vgene")]
        clusterCount <- clusterCount[!duplicated(clusterCount), ]
        clusterCount$ClusterCount <- 1
        clusterCount2 <- as.data.frame(ddply(clusterCount,c("Clust_size_order"),numcolwise(sum)))
        clusterCount2 <- subset(clusterCount2,clusterCount2$ClusterCount>1)
        df4_clusterCount2 <- merge(df4,clusterCount2,by = c("Clust_size_order"))
        df4_clusterCount2
        # Sample Count
        SampCount <- df4_clusterCount2[,names(df4_clusterCount2) %in% c("Clust_size_order",input$Samp_col)]
        SampCount <- SampCount[!duplicated(SampCount), ]
        SampCount$Sample_count <- 1
        df9 <- as.data.frame(ddply(SampCount,c("Clust_size_order"),numcolwise(sum)))
        df10 <- merge(df9,df4_clusterCount2,by = "Clust_size_order")

        # Calculating the priority
        df10$priority <- 1/(df10$Total_cloneCount * df10$ClusterCount * df10$Sample_count)
        df10 <- df10[order(df10$priority,decreasing = F),]
        df10
        df10$priority[is.na(df10$priority)] <- 0
        df10 <- subset(df10,df10$priority>0)
        df10
        # updated order
        updatedOrder <- df10[,names(df10) %in% c("priority","Clust_size_order")]
        updatedOrder <- updatedOrder[!duplicated(updatedOrder), ]
        updatedOrder <- updatedOrder[order(updatedOrder$priority,decreasing = F),]
        updatedOrder$Updated_order <- 1:dim(updatedOrder)[1]
        updatedOrder

        # final data frame
        df7 <- merge(updatedOrder, df10, by = c("Clust_size_order","priority"))
        df7 <- df7[order(df7$priority,decreasing = F),]

        clusterBD <- df7 %>%
          select(all_of(c(input$Samp_col,"Sample_count","Total_cloneCount","ClusterCount","priority","Updated_order")), everything())
        clusterBD
      } else {

      }
    })

    clusTCR2_df <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )

      if (input$chain_TCR == "TRAG") {
        if (length(input.data_sc_clusTCR_AG())>0) {
          AG_cluster()
        }
      }
      else if (input$chain_TCR == "TRBD") {
        if (length(input.data_sc_clusTCR_BD())>0) {
          BD_cluster()
        }
      }

      else {

      }

    })

    output$Tb_ClusTCR_selected <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      cluster <- clusTCR2_df()
      validate(
        need(nrow(cluster)>0,
             "Upload clusTCR table, which is needed for TCR -> UMAP section")
      )
      cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]
      Network_df <- cluster[order(cluster$Updated_order),]
      Network_df <- Network_df[Network_df$Updated_order  %in% input$Clusters_to_dis_PIE,]
      Network_df$Selected <- Network_df[,names(Network_df) %in% input$Colour_By_this]
      Network_df


      df3.meta3 <-  as.data.frame(table(Network_df$ID_Column,Network_df$Selected))
      total.condition <- as.data.frame(ddply(df3.meta3,"Var1",numcolwise(sum)))
      dim(total.condition)[1]
      dim(df3.meta3)[1]
      emtpy <- matrix(nrow =dim(df3.meta3)[1],ncol=dim(total.condition)[1])

      for (i in 1:dim(df3.meta3)[1]) {

        emtpy[i,] <- ifelse(df3.meta3$Var1[i]==total.condition$Var1[1:dim(total.condition)[1]],
                            total.condition[total.condition$Var1==total.condition$Var1[1:dim(total.condition)[1]],2],F)
      }
      df3.meta3$n <- df3.meta3$Freq/rowSums(emtpy)
      df3.meta3

    })

    output$downloadtb_Tb_ClusTCR_selected <- downloadHandler(
      filename = function(){
        paste0(input$Clusters_to_dis_PIE,"_Cluster_",today(),"_over_rep.csv")
      },
      content = function(file){
        df <- as.data.frame(clusTCR2_df())
        write.csv(df,file, row.names = F)
      })

    # cluster to display ------
    observe({
      cluster <- clusTCR2_df()

      validate(
        need(nrow(cluster)>0,
             "upload clustering")
      )
      cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]

      # cluster <- cluster[cluster$Clust_size_order %in% input$lower_cluster:input$upper_cluster,]
      updateSelectInput(
        session,
        "Clusters_to_dis_PIE",
        choices=unique(cluster$Updated_order),
        selected = unique(cluster$Updated_order)[1]
      )
    }) # cluster to display


    # cluster UMAP (1 to display) ClusTCR -----
    cols_clust_UMAP <- reactive({
      cluster <- clusTCR2_df()
      validate(
        need(nrow(cluster)>0,
             "Upload clusTCR table")
      )
      req(cluster,input$Clusters_to_dis_PIE,input$Colour_By_this, input$Samp_col)
      cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]
      # names(cluster)[names(cluster) %in% input$Samp_col_cluster] <- "ID_Column"

      cluster <- cluster[cluster$Updated_order %in% input$Clusters_to_dis_PIE,]
      cluster$colour <- cluster[,names(cluster) %in% input$Colour_By_this]
      cluster$colour <- gsub("_"," ",cluster$colour)
      cluster$colour <- factor(cluster$colour, levels = unique(cluster$colour))
      cluster$colour <- gsub("NA",NA,cluster$colour)

      num <- as.data.frame(unique(cluster$colour))
      num <- as.data.frame(num[complete.cases(num)==T,])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep="_"), paste(num[i,]), col.gg[i])
        })
      }
      else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep="_"), paste(num[i,]), col.hcl[i])
        })
      }
      else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep="_"), paste(num[i,]), col.topo[i])
        })
      }
      else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep="_"), paste(num[i,]), heat_col[i])
        })
      }
      else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep="_"), paste(num[i,]), col.terrain[i])
        })
      }

      else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep="_"), paste(num[i,]), palette_rainbow[i])
        }) }
      else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep="_"), paste(num[i,]), palette1[i])
        })

      }
      else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep="_"), paste(num[i,]), input$one.colour.default)
        })


      } # one colour

    })
    output$myPanel_cols_clust_UMAP <- renderUI({cols_clust_UMAP()})

    colors_cols_cols_clust_UMAP <- reactive({
      cluster <- clusTCR2_df()
      validate(
        need(nrow(cluster)>0,
             "Upload clusTCR table")
      )
      req(cluster,input$Clusters_to_dis_PIE,input$Colour_By_this)
      cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]

      cluster <- cluster[cluster$Updated_order %in% input$Clusters_to_dis_PIE,]
      cluster$colour <- cluster[,names(cluster) %in% input$Colour_By_this]
      cluster$colour <- gsub("_"," ",cluster$colour)
      cluster$colour <- factor(cluster$colour, levels = unique(cluster$colour))
      cluster$colour <- gsub("NA",NA,cluster$colour)

      num <- as.data.frame(unique(cluster$colour))
      num <- as.data.frame(num[complete.cases(num)==T,])

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.cols_clust_UMAP", i, sep="_")]]
      })
    })

    UMAP_ClusTCR2 <- reactive({
      cluster <- clusTCR2_df()
      validate(
        need(nrow(cluster)>0,
             "Upload clusTCR table, which is needed for TCR -> UMAP section")
      )
      req(cluster,input$Clusters_to_dis_PIE,input$Colour_By_this)
      cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]

      cluster <- cluster[cluster$Updated_order %in% input$Clusters_to_dis_PIE,]
      cluster$colour <- cluster[,names(cluster) %in% input$Colour_By_this]
      cluster$colour <- gsub("_"," ",cluster$colour)
      cluster$colour <- factor(cluster$colour, levels = unique(cluster$colour))
      cluster$colour <- gsub("NA",NA,cluster$colour)

      len.colour <- length(unique(cluster$colour))
      col.df <- as.data.frame(unique(cluster$colour))
      col.df <- as.data.frame(col.df[complete.cases(col.df)==T,])
      col.df$col <- unlist(colors_cols_cols_clust_UMAP())

      figure <- ggplot(data=cluster,aes(x=UMAP_1,UMAP_2,colour=colour))+
        geom_point(size = input$size.dot.umap)+
        scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 20),values = col.df$col,na.value=input$NA_col_analysis) +
        theme_bw()+
        theme(
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
          # strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=0),
          axis.title.x = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),

        )

      if (input$Split_by_group=="no") {
        figure
      }
      else {
        figure+ facet_wrap(~ID_Column, nrow = input$wrap_row)
      }

    })

    output$UMAP_ClusTCR2_plot <- renderPlot({
      UMAP_ClusTCR2()
    })

    output$downloadPlot_UMAP_ClusTCR2_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_UMAP_ClusTCR2_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_UMAP_ClusTCR2_plot,height=input$height_UMAP_ClusTCR2_plot, onefile = FALSE) # open the pdf device
        plot(UMAP_ClusTCR2())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_UMAP_ClusTCR2_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_UMAP_ClusTCR2_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_UMAP_ClusTCR2_plot,
            height = input$height_png_UMAP_ClusTCR2_plot,
            res = input$resolution_PNG_UMAP_ClusTCR2_plot)
        plot(UMAP_ClusTCR2())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    ##### Clustering motif plot ----

    motif_plot_sc <- reactive({
      cluster <- clusTCR2_df()
      validate(
        need(nrow(cluster)>0,
             "Upload clusTCR table")
      )
      cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]
      Network_df <- cluster[order(cluster$Updated_order),]
      # Network_df2 <-
      Network_df %>% distinct(CDR3_Vgene, .keep_all = TRUE) # make Unique
      Network_df <- Network_df %>% distinct(CDR3_Vgene, .keep_all = TRUE) # make Unique
      Motif_from_cluster_file(Network_df,Clust_selected = input$Clusters_to_dis_PIE,selected_cluster_column = "Updated_order")
      # ?Motif_from_cluster_file
    })
    output$Motif_ClusTCR2_cluster <- renderPlot({
      motif_plot_sc()
    })

    # render which cases were contributing to the cluster
    output$print_unique_cases <- renderPrint({

      cluster <- clusTCR2_df()
      validate(
        need(nrow(cluster)>0,
             "Upload clusTCR table, which is needed for TCR -> UMAP section")
      )
      cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]

      df <- cluster[cluster$Updated_order  %in% input$Clusters_to_dis_PIE,]


      if (input$chain_TCR == "TRAG") {
        cat("This Motif if from: ",noquote(unique(df$ID_Column)),"; the TCR is:",noquote(unique(df$Selected_V_AG)),noquote(unique(df$j_gene_AG)))
      }
      else if (input$chain_TCR == "TRBD") {
        cat("This Motif if from:",noquote(unique(df$ID_Column)),"; the TCR is:",noquote(unique(df$Selected_V_BD)),noquote(unique(df$d_gene_BD)),noquote(unique(df$j_gene_BD)))
      }

      else { # BCR repertoire

      }

    })
    #download the motif plot -----
    output$downloadPlot_Motif_ClusTCR2_cluster <- downloadHandler(
      filename = function() {
        x <- today()
        paste("Motif_ClusTCR2_cluster_plot_",input$chain_TCR,"_",input$Clusters_to_dis_PIE,"_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Motif_ClusTCR2_cluster,height=input$height_Motif_ClusTCR2_cluster, onefile = FALSE) # open the pdf device
        plot(motif_plot_sc())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Motif_ClusTCR2_cluster <- downloadHandler(
      filename = function() {
        x <- today()
        paste("Motif_ClusTCR2_cluster_plot_",input$chain_TCR,"_",input$Clusters_to_dis_PIE,"_",x,  ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Motif_ClusTCR2_cluster,
            height = input$height_png_Motif_ClusTCR2_cluster,
            res = input$resolution_PNG_Motif_ClusTCR2_cluster)
        plot(motif_plot_sc())
        dev.off()},   contentType = "application/png" # MIME type of the image

    )

    #### cluster pie chart function -----
    cols_clusTCR2_pie <- reactive({
      cluster <- clusTCR2_df()
      validate(
        need(nrow(cluster)>0,
             "Upload clusTCR table, which is needed for TCR -> UMAP section")
      )
      cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]
      Network_df <- cluster[order(cluster$Updated_order),]
      Network_df <- Network_df[Network_df$Updated_order  %in% input$Clusters_to_dis_PIE,]
      Network_df$Selected <- Network_df[,names(Network_df) %in% input$Colour_By_this]
      Network_df$Selected <- gsub("NA",NA,Network_df$Selected)

      num <- as.data.frame(unique(Network_df$Selected))
      num <- as.data.frame(num[complete.cases(num)==T,])


      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), col.gg[i])
        })
      }
      else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), col.hcl[i])
        })
      }
      else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), col.topo[i])
        })
      }
      else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), heat_col[i])
        })
      }
      else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), col.terrain[i])
        })
      }

      else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), palette_rainbow[i])
        }) }
      else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), palette1[i])
        })

      }
      else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), input$one.colour.default)
        })


      } # one colour

    })

    output$myPanel_cols_clusTCR2_pie <- renderUI({cols_clusTCR2_pie()})

    colors_cols_clusTCR2_pie <- reactive({
      cluster <- clusTCR2_df()
      validate(
        need(nrow(cluster)>0,
             "Upload clusTCR table, which is needed for TCR -> UMAP section")
      )
      cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]
      Network_df <- cluster[order(cluster$Updated_order),]
      Network_df <- Network_df[Network_df$Updated_order  %in% input$Clusters_to_dis_PIE,]
      Network_df$Selected <- Network_df[,names(Network_df) %in% input$Colour_By_this]
      Network_df$Selected <- gsub("NA",NA,Network_df$Selected)

      num <- as.data.frame(unique(Network_df$Selected))
      num <- as.data.frame(num[complete.cases(num)==T,])


      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.cols_clusTCR2_pie", i, sep="_")]]
      })
    })

    Pie_chart_ClusTCR2 <- reactive({
      cluster <- clusTCR2_df()
      validate(
        need(nrow(cluster)>0,
             "Upload clusTCR table, which is needed for TCR -> UMAP section")
      )
      cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]
      Network_df <- cluster[order(cluster$Updated_order),]
      Network_df <- Network_df[Network_df$Updated_order  %in% input$Clusters_to_dis_PIE,]
      Network_df$Selected <- Network_df[,names(Network_df) %in% input$Colour_By_this]
      Network_df$Selected <- gsub("NA",NA,Network_df$Selected)

      df3.meta3 <-  as.data.frame(table(Network_df$ID_Column,Network_df$Selected))
      total.condition <- as.data.frame(ddply(df3.meta3,"Var1",numcolwise(sum)))
      dim(total.condition)[1]
      dim(df3.meta3)[1]
      emtpy <- matrix(nrow =dim(df3.meta3)[1],ncol=dim(total.condition)[1])

      for (i in 1:dim(df3.meta3)[1]) {

        emtpy[i,] <- ifelse(df3.meta3$Var1[i]==total.condition$Var1[1:dim(total.condition)[1]],
                            total.condition[total.condition$Var1==total.condition$Var1[1:dim(total.condition)[1]],2],F)
      }
      df3.meta3$n <- df3.meta3$Freq/rowSums(emtpy)
      # df3.meta3
      df.col <- unlist(colors_cols_clusTCR2_pie())

      ggplot(df3.meta3,aes(x="", y=n, fill=Var2, group = Var1)) +
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0)  +
        theme_void(20) +
        facet_wrap(~Var1, nrow = input$wrap_row) +
        theme(
          legend.key.size = unit(1, 'cm'),
          legend.title = element_blank()) +
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 20),values = df.col, na.value = input$NA_col_analysis) +
        theme(strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
              legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.position = input$legend_position,
              legend.title = element_blank()
        )

    })

    output$Pie_ClusTCR2_plot <- renderPlot({
      Pie_chart_ClusTCR2()
    })

    output$downloadPlot_Pie_ClusTCR2_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("Pie_ClusTCR2_plot_",input$chain_TCR,"_",input$Clusters_to_dis_PIE,"_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Pie_ClusTCR2_plot,height=input$height_Pie_ClusTCR2_plot, onefile = FALSE) # open the pdf device
        plot(Pie_chart_ClusTCR2())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Pie_ClusTCR2_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("Pie_ClusTCR2_plot_",input$chain_TCR,"_",input$Clusters_to_dis_PIE,"_",x,  ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Pie_ClusTCR2_plot,
            height = input$height_png_Pie_ClusTCR2_plot,
            res = input$resolution_PNG_Pie_ClusTCR2_plot)
        plot(Pie_chart_ClusTCR2())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )


    # Cluster stats   -----

    compare.stat_Cluster <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )

      # req(input$V_call_clust_sc,input$junction_clust_sc)
      md <- sc@meta.data

      cluster <- clusTCR2_df()

      cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]
      cluster <- cluster[order(cluster$Updated_order),]

      rownames(cluster) <- cluster$Cell_Index

      checking <- cluster[,names(cluster) %in% c("Updated_order","Cell_Index")]
      md.checking <- merge(md,checking,by="Cell_Index",all.x=T)
      md.checking <- md.checking[order(md.checking$order),]
      rownames(md.checking) <- md.checking$Cell_Index

      md.checking$Clust_selected <- ifelse(md.checking$Updated_order == input$Clusters_to_dis_PIE,input$Clusters_to_dis_PIE,"NS")
      md.checking$Clust_selected[is.na(md.checking$Clust_selected)] <- "NS"
      md.checking <- md.checking[order(md.checking$order),]
      md.checking
      sc@meta.data <- md.checking
      Idents(object = sc) <- sc@meta.data$Clust_selected

      name.check.clust <- input$Clusters_to_dis_PIE
      min.pct.expression<- input$min_point_ #standard setting: 0.25
      min.logfc<-  input$LogFC_ #0.25 is standard

      markers.fm.list <- FindMarkers(sc, ident.1 = name.check.clust, min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
      markers.fm.list2 <- subset(markers.fm.list,markers.fm.list$p_val_adj < input$pval.ex.filter)
      as.data.frame(markers.fm.list2)

    })

    output$compare.stat_Cluster_DT <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength =10, scrollX = TRUE),{
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )

      # req(input$V_call_clust_sc,input$junction_clust_sc)
      compare.stat_Cluster()
    })

    output$downloaddf_FindMarker_Cluster <- downloadHandler(
      filename = function(){
        x= today()
        paste(input$Clusters_to_dis_PIE,"_Cluster_Stats_",x,".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(compare.stat_Cluster())
        write.csv(df,file, row.names = T)
      } )

    # Cluster dot plot -----

    output$Cluster_of_interest_DF <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )

      # req(input$V_call_clust_sc,input$junction_clust_sc)
      md <- sc@meta.data

      cluster <- clusTCR2_df()

      cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]
      cluster <- cluster[order(cluster$Updated_order),]

      rownames(cluster) <- cluster$Cell_Index

      checking <- cluster[,names(cluster) %in% c("Updated_order","Cell_Index")]
      md.checking <- merge(md,checking,by="Cell_Index",all.x=T)
      md.checking <- md.checking[order(md.checking$order),]
      rownames(md.checking) <- md.checking$Cell_Index

      md.checking$Clust_selected <- ifelse(md.checking$Updated_order == input$Clusters_to_dis_PIE,input$Clusters_to_dis_PIE,"NS")
      md.checking$Clust_selected[is.na(md.checking$Clust_selected)] <- "NS"
      md.checking <- md.checking[order(md.checking$order),]
      md.checking


    })

    all_expression_plot_cluster <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )

      # req(input$V_call_clust_sc,input$junction_clust_sc)
      md <- sc@meta.data

      cluster <- clusTCR2_df()

      cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]
      cluster <- cluster[order(cluster$Updated_order),]

      rownames(cluster) <- cluster$Cell_Index

      checking <- cluster[,names(cluster) %in% c("Updated_order","Cell_Index")]
      md.checking <- merge(md,checking,by="Cell_Index",all.x=T)
      md.checking <- md.checking[order(md.checking$order),]
      rownames(md.checking) <- md.checking$Cell_Index

      md.checking$Clust_selected <- ifelse(md.checking$Updated_order == input$Clusters_to_dis_PIE,input$Clusters_to_dis_PIE,"NS")
      md.checking$Clust_selected[is.na(md.checking$Clust_selected)] <- "NS"
      md.checking <- md.checking[order(md.checking$order),]
      md.checking

      sc@meta.data <- md.checking
      Idents(object = sc) <- sc@meta.data$Clust_selected

      if (input$restrict.dotpot.clust == F ) {
        list.names <- rownames(compare.stat_Cluster())
      }

      else {
        list.names <- rownames(compare.stat_Cluster())
        list.names <- list.names[1:input$restrict.dotpot.num.clust]
      }

      size_legend = input$Bar_legend_size-2

      DotPlot(sc, features = list.names) +
        RotatedAxis() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size, angle = 90),
          axis.title.x = element_blank(),
          legend.title = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.text = element_text(colour="black", size=size_legend,family=input$font_type),
          legend.position = input$legend_position,
        ) +
        scale_colour_gradient2(low = input$low.dotplot.clust, mid = input$middle.dotplot.clust, high = input$high.dotplot.clust)+
        scale_x_discrete(labels = label_wrap(20)) +
        scale_y_discrete(labels = label_wrap(20))


    })

    output$all_expression_dotplot_cluster <- renderPlot({
      all_expression_plot_cluster()
    })

    output$downloadPlot_all_expression_dotplot_clust <- downloadHandler(
      filename = function() {
        paste(input$Clusters_to_dis_PIE,"_cluster_dotplot","_",today(), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_all_expression_dotplot_clust,
            height=input$height_all_expression_dotplot_clust, onefile = FALSE) # open the pdf device
        df <- all_expression_plot_cluster()
        plot(df)
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_all_expression_dotplot_clust  <- downloadHandler(
      filename = function() {
        paste(input$Clusters_to_dis_PIE,"_cluster_dotplot","_",today(), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_all_expression_dotplot_clust,
            height = input$height_png_all_expression_dotplot_clust,
            res = input$resolution_PNG_all_expression_dotplot_clust)
        df <- all_expression_plot_cluster()

        plot(df)
        dev.off()},   contentType = "application/png" # MIME type of the image
    )




    # Over representation analysis for Cluster  -----
    Over_rep_cluster_old <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )

      # req(input$V_call_clust_sc,input$junction_clust_sc)
      md <- sc@meta.data

      cluster <- clusTCR2_df()

      cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]
      cluster <- cluster[order(cluster$Updated_order),]

      rownames(cluster) <- cluster$Cell_Index

      checking <- cluster[,names(cluster) %in% c("Updated_order","Cell_Index")]
      md.checking <- merge(md,checking,by="Cell_Index",all.x=T)
      md.checking <- md.checking[order(md.checking$order),]
      rownames(md.checking) <- md.checking$Cell_Index

      md.checking$Clust_selected <- ifelse(md.checking$Updated_order == input$Clusters_to_dis_PIE,input$Clusters_to_dis_PIE,"NS")
      md.checking$Clust_selected[is.na(md.checking$Clust_selected)] <- "NS"
      md.checking <- md.checking[order(md.checking$order),]
      md.checking

      sc@meta.data <- md.checking
      Idents(object = sc) <- sc@meta.data$Clust_selected
      sc@meta.data <- checking
      # require()

      geneSet <- read.csv(system.file("OverRep","GeneSets.csv",package = "STEGO.R"),header = T)

      background.genes.name <- as.data.frame(rownames(sc@assays$RNA$scale.data))
      names(background.genes.name) <- "V1"
      background.genes <- length(rownames(sc@assays$RNA$scale.data))

      geneSet$background.genes <- background.genes
      head(geneSet)
      markers.fm.list <- compare.stat_Cluster()
      DEx.genes <- as.data.frame(rownames(markers.fm.list))
      names(DEx.genes) <- "V1"
      total.sig <- length(DEx.genes$V1)
      geneSet$total.sig <- length(DEx.genes$V1)

      geneSet$background.geneset <- NA
      geneSet$background.geneset.name <- NA
      geneSet$in.geneset <- NA
      geneSet$in.geneset.name <- NA

      if(input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        geneSet$GeneSet <- gsub("-",".",geneSet$GeneSet)
      }

      if(input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        require(stringr)
        geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
      }

      for (i in 1:dim(geneSet)[1]) {
        # listed GeneSet
        message(paste("GeneSet: ", i))
        Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet,";")[i])
        names(Gene.set.testing) <- "V1"
        Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
        names(Gene.set.testing2) <- "V1"
        # message(paste(dim(Gene.set.testing2)[1],"GeneSet total"))
        # background genes
        background.overlap <- merge(Gene.set.testing2,background.genes.name,by= "V1")
        # message(paste(dim(background.overlap)[1],"in Background"))
        geneSet$background.geneset[i] <- length(background.overlap$V1)
        geneSet$background.geneset.name[i] <- as.character(paste(unlist(background.overlap[1]), collapse=';'))
        # in sig gene list
        overlap <- merge(background.overlap,DEx.genes,by= "V1")
        # message(paste(dim(overlap)[1],"# Sig genes"))
        geneSet$in.geneset[i] <- length(overlap$V1)
        geneSet$in.geneset.name[i] <- as.character(paste(unlist(overlap[1]), collapse=';'))

      }

      geneSet2 <- subset(geneSet,geneSet$in.geneset>0)

      for (i in 1:dim(geneSet2)[1]) {
        tota.gene.set <- geneSet2$background.geneset[i] # genes that are identified in background
        tota.gene.set
        in.geneset <-  geneSet2$in.geneset[i]# DEx in geneset

        background.genes
        not.in.total <- background.genes - tota.gene.set
        not.in.geneset.sig <- total.sig - in.geneset
        d <- data.frame( gene.in.interest=c( in.geneset, not.in.geneset.sig),gene.not.interest=c( tota.gene.set, not.in.total))
        row.names(d) <- c("In_category", "not_in_category")

        if (in.geneset>0) {
          geneSet2$p.val[i] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
          geneSet2$lowerCI[i] <-  unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
          geneSet2$upperCI[i] <-unlist(fisher.test(d)$conf.int)[2]
          geneSet2$OR[i] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1],3)
        }

        else {
          geneSet2$p.value[i] <- "-"
          geneSet2$lowerCI[i] <-  "-"
          geneSet2$upperCI[i] <- "-"
          geneSet2$OR[i] <- "-"
        }
        # message(print(d))
        # message(print(round(prop.table(d),3)))
      }
      geneSet2
    })

    Over_rep_cluster <- reactive({
      geneSet2 <- Over_rep_cluster_old()
      validate(
        need(nrow(geneSet2)>0,
             "Upload Files")
      )

      geneSet2 <- geneSet2[order(geneSet2$p.val,decreasing = F),]
      geneSet2 <- subset(geneSet2,geneSet2$in.geneset>=input$in.geneset.cutoff_Clust)
      geneSet2 <- subset(geneSet2,geneSet2$p.val<=input$p.val_cutoff_Clust)
      geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
      geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
      geneSet2

    })

    output$Over_rep_Cluster_Tab <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1,2,5,10,20,50,100), pageLength = 20, scrollX = TRUE),{
      Over_rep_cluster()
    })

    output$downloadtb_over.rep.cluster <- downloadHandler(
      filename = function(){
        paste0(input$Clusters_to_dis_PIE,"_Cluster_",today(),"_over_rep.csv")
      },
      content = function(file){
        df <- as.data.frame(Over_rep_cluster())
        write.csv(df,file, row.names = F)
      })


    ##### Overlap -----
    #### upset plot -----

    Upset_plot <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "upload file")
      )
      df <- sc@meta.data
      df <- as.data.frame(df)
      unique.df <- unique(df[,names(df) %in% c(input$Samp_col,input$V_gene_sc) ])
      names(unique.df) <- c("group","chain")
      unique.df <- unique.df[unique.df$group %in% input$ID_Column_factor,]

      unique.df <- subset(unique.df,unique.df$chain != "NA")
      unique.df <- subset(unique.df,unique.df$group != "NA")
      unique.df$cloneCount <- 1
      mat <- acast(unique.df, chain~group, value.var="cloneCount")
      mat[is.na(mat)] <- 0
      mat <- as.data.frame(mat)
      df.x <- make_comb_mat(mat)
      list.names <- as.character(input$ID_Column_factor)

      ht = draw(UpSet(df.x,
                      pt_size = unit(5, "mm"),
                      lwd = 1,
                      row_names_gp =  gpar(fontfamily = input$font_type, fontsize = 12),
                      column_names_gp = gpar(fontfamily = input$font_type,fontsize =12),
                      top_annotation = upset_top_annotation(df.x,
                                                            add_numbers = T,
                                                            numbers_gp = gpar(fontfamily = input$font_type,fontsize =12),
                                                            annotation_name_gp = gpar(fontfamily = input$font_type, fontsize = 12),
                                                            gp = gpar(fill = "black"),
                      ),
                      right_annotation = upset_right_annotation(df.x,
                                                                add_numbers = T,
                                                                numbers_gp = gpar(fontfamily = input$font_type,fontsize =12),
                                                                annotation_name_gp = gpar(fontfamily = input$font_type,fontsize=12),
                                                                gp = gpar(fill = "black"),
                      ),
                      set_order  =  list.names

      ), padding = unit(c(20, 20, 20, 20), "mm"))
      ht
    })
    output$Upset_plot_overlap <- renderPlot({
      Upset_plot()
    })

    output$downloadPlot_Upset_plot_overlap <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Samp_col,"_",input$V_gene_sc,"_overlap_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Upset_plot_overlap,height=input$height_Upset_plot_overlap, onefile = FALSE) # open the pdf device
        plot(Upset_plot())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Upset_plot_overlap <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Samp_col,"_",input$V_gene_sc,"_overlap_",x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Upset_plot_overlap,height = input$height_png_Upset_plot_overlap,res = input$resolution_PNG_Upset_plot_overlap)
        plot(Upset_plot())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    # upset plot table -----
    Upset_plot_overlap <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "upload file")
      )

      df <- sc@meta.data
      df <- as.data.frame(df)
      unique.df <- unique(df[,names(df) %in% c(input$Samp_col,input$V_gene_sc) ])
      names(unique.df) <- c("group","chain")
      unique.df <- subset(unique.df,unique.df$chain != "NA")
      unique.df <- subset(unique.df,unique.df$group != "NA")
      unique.df$cloneCount <- 1
      mat <- acast(unique.df, chain~group, value.var="cloneCount")
      mat[is.na(mat)] <- 0
      Count_data <- as.data.frame(rowSums(mat))
      names(Count_data) <- "V1"
      unique.df <- (df[,names(df) %in% c(input$Samp_col,input$V_gene_sc) ])
      names(unique.df) <- c("group","chain")
      unique.df <- subset(unique.df,unique.df$chain != "NA")
      unique.df <- subset(unique.df,unique.df$group != "NA")
      unique.df$cloneCount <- 1
      mat <- acast(unique.df, chain~group, value.var="cloneCount")
      mat[is.na(mat)] <- 0
      sum_data <- as.data.frame(rowSums(mat))
      names(sum_data) <- "V1"
      mat <- as.data.frame(mat)
      mat$TotalSamps <-Count_data$V1
      mat$CloneTotal <-sum_data$V1
      mat
    })

    output$Upset_plot_overlap_Tb <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      Upset_plot_overlap()
    })

    output$downloaddf_Upset_plot_overlap_Tb <- downloadHandler(
      filename = function(){
        paste(input$name.BD,"_Overlap_",gsub("-", ".", Sys.Date()),".tsv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(Upset_plot_overlap())
        write.csv(df,file)
      } )


    # overlap table with UMAP and expression -----
    overlap_table <- reactive({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      reduction <- (sc@reductions$umap)
      UMAP <- as.data.frame(reduction@cell.embeddings)
      names(UMAP)[1:2] <- toupper(names(UMAP)[1:2])
      UMAP$Cell_Index <- rownames(UMAP)
      meta.data <- as.data.frame(sc@meta.data)
      umap.meta <- merge(UMAP,meta.data,by="Cell_Index")
      names(umap.meta)[names(umap.meta) %in% input$V_gene_sc] <- "chain"
      names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
      names(umap.meta)[names(umap.meta) %in% input$Colour_By_this] <- "Selected_function"
      names(umap.meta)[names(umap.meta) %in% input$Split_group_by_overview] <- "Selected_group"
      # sc <- merge(umap.meta,TCR_Expanded(),by=c("v_gene_selected","ID_Column"),all.x=T)

      Upset_plot_overlap <- Upset_plot_overlap()
      Upset_plot_overlap_top <- subset(Upset_plot_overlap,Upset_plot_overlap$sum>1)
      Upset_plot_overlap_top$chain <- rownames(Upset_plot_overlap_top)
      umap.meta.overlap <- merge(Upset_plot_overlap_top,umap.meta,by="chain")
      umap.meta.overlap
    })

    # overlap UMAP plot
    create_UMAP_overlap <- reactive({
      umap.meta.overlap <- overlap_table()
      validate(
        need(nrow(umap.meta.overlap)>0,
             error_message_val_UMAP)
      )

      ggplot(umap.meta.overlap,aes(x=UMAP_1,UMAP_2,colour=Selected_function))+
        geom_point()+
        scale_color_manual(values = rainbow(length(unique(umap.meta.overlap$Selected_function))), na.value=input$NA_col_analysis)+
        theme_bw() +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=0),
          axis.title.x = element_text(colour="black",angle=0,vjust=.5,face="plain",family=input$font_type,size = input$title.text.sizer2),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )


    })
    output$UMAP_overlap_classification <- renderPlot({
      create_UMAP_overlap()
    })


    # Umap Overlap plot -------

    Overlap_Pie_chart_Class <- reactive({
      top_BD_cluster <- overlap_table()

      # df.col <- unlist(colors_pie())

      df3.meta3 <-  as.data.frame(table(top_BD_cluster$Selected_group,top_BD_cluster$Selected_function))
      total.condition <- as.data.frame(ddply(df3.meta3,"Var1",numcolwise(sum)))
      emtpy <- matrix(nrow =dim(df3.meta3)[1],ncol=dim(total.condition)[1])
      for (i in 1:dim(df3.meta3)[1]) {
        emtpy[i,] <- ifelse(df3.meta3$Var1[i]==total.condition$Var1[1:dim(total.condition)[1]],
                            total.condition[total.condition$Var1==total.condition$Var1[1:dim(total.condition)[1]],2],F)
      }

      df3.meta3$n <- df3.meta3$Freq/rowSums(emtpy)


      ggplot(df3.meta3,aes(x="", y=n, fill=Var2, group = Var1)) +
        geom_bar(stat="identity", width=1)+
        coord_polar("y", start=0)  +
        theme_void(20) +
        facet_wrap(~Var1, nrow = input$wrap_row) +
        theme(
          legend.key.size = unit(1, 'cm'))+
        # scale_fill_manual(values = df.col, na.value = input$NA_col_analysis) +
        theme(strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
              legend.text = element_text(size = input$Bar_legend_size, family = input$font_type),
              legend.title = element_blank()
        )

    })

    output$Classification_Overlap_pie <- renderPlot({
      Overlap_Pie_chart_Class()
    })

    # Over representation analysis for Cluster  -----

    Over_rep_Overlap <- reactive({
      sc <-input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )

      df <- sc@meta.data
      # require()

      geneSet <- read.csv(system.file("OverRep","GeneSets.csv",package = "STEGO.R"),header = T)

      background.genes.name <- as.data.frame(rownames(sc@assays$RNA$scale.data))
      names(background.genes.name) <- "V1"
      background.genes <- length(rownames(sc@assays$RNA$scale.data))

      geneSet$background.genes <- background.genes
      head(geneSet)
      markers.fm.list <- compare.stat()
      DEx.genes <- as.data.frame(rownames(markers.fm.list))
      names(DEx.genes) <- "V1"
      total.sig <- length(DEx.genes$V1)
      geneSet$total.sig <- length(DEx.genes$V1)

      geneSet$background.geneset <- NA
      geneSet$background.geneset.name <- NA
      geneSet$in.geneset <- NA
      geneSet$in.geneset.name <- NA

      if(input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        geneSet$GeneSet <- gsub("-",".",geneSet$GeneSet)
      }

      if(input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        require(stringr)
        geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
      }

      for (i in 1:dim(geneSet)[1]) {
        # listed GeneSet
        message(paste("GeneSet: ", i))
        Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet,";")[i])
        names(Gene.set.testing) <- "V1"
        Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
        names(Gene.set.testing2) <- "V1"
        # message(paste(dim(Gene.set.testing2)[1],"GeneSet total"))
        # background genes
        background.overlap <- merge(Gene.set.testing2,background.genes.name,by= "V1")
        # message(paste(dim(background.overlap)[1],"in Background"))
        geneSet$background.geneset[i] <- length(background.overlap$V1)
        geneSet$background.geneset.name[i] <- as.character(paste(unlist(background.overlap[1]), collapse=';'))
        # in sig gene list
        overlap <- merge(background.overlap,DEx.genes,by= "V1")
        # message(paste(dim(overlap)[1],"# Sig genes"))
        geneSet$in.geneset[i] <- length(overlap$V1)
        geneSet$in.geneset.name[i] <- as.character(paste(unlist(overlap[1]), collapse=';'))

      }

      geneSet2 <- subset(geneSet,geneSet$in.geneset>0)

      for (i in 1:dim(geneSet2)[1]) {
        tota.gene.set <- geneSet2$background.geneset[i] # genes that are identified in background
        tota.gene.set
        in.geneset <-  geneSet2$in.geneset[i]# DEx in geneset

        background.genes
        not.in.total <- background.genes - tota.gene.set
        not.in.geneset.sig <- total.sig - in.geneset
        d <- data.frame( gene.in.interest=c( in.geneset, not.in.geneset.sig),gene.not.interest=c( tota.gene.set, not.in.total))
        row.names(d) <- c("In_category", "not_in_category")

        if (in.geneset>0) {
          geneSet2$p.val[i] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
          geneSet2$lowerCI[i] <-  unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
          geneSet2$upperCI[i] <-unlist(fisher.test(d)$conf.int)[2]
          geneSet2$OR[i] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1],3)
        }

        else {
          geneSet2$p.value[i] <- "-"
          geneSet2$lowerCI[i] <-  "-"
          geneSet2$upperCI[i] <- "-"
          geneSet2$OR[i] <- "-"
        }
        # message(print(d))
        # message(print(round(prop.table(d),3)))
      }
      geneSet2 <- geneSet2[order(geneSet2$p.val,decreasing = F),]
      geneSet2 <- geneSet2[order(geneSet2$p.val,decreasing = F),]
      geneSet2 <- subset(geneSet2,geneSet2$in.geneset>=input$in.geneset.cutoff_Exp)
      geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
      geneSet2$pval.bonferroni.adj <- p.adjust(geneSet2$p.val, method = "bonferroni")
      geneSet2 <- subset(geneSet2,geneSet2$p.val<=input$p.val_cutoff_Exp)
      # geneSet2 <- subset(geneSet2,geneSet2$pval.BH.adj<=input$adjust_cutoff_top)

      geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
      geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
      name.list <- c("Geneset_ID","p.val","FDR","Bonferroni","OR","lowerCI","upperCI","in.geneset.name","in.geneset","background.geneset","total.sig","background.genes","background.geneset.name")
      geneSet2 <- geneSet2 %>%
        select(all_of(name.list), everything())
    })

    output$Over_rep_overlap_Tab <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1,2,5,10,20,50,100), pageLength = 20, scrollX = TRUE),{
      Over_rep_Overlap()
    })

    output$downloadtb_over.rep.overlap <- downloadHandler(
      filename = function(){

        paste0("Overlap","_",today(),"_over_rep.csv")
      },
      content = function(file){
        df <- as.data.frame(Over_rep_Overlap())
        write.csv(df,file, row.names = F)
      } )

    # marker analysis -----

    ### single marker feature plot -----

    MainTcell_counts <- reactive({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )

      # req(input$V_call_clust_sc,input$junction_clust_sc)
      md <- sc@meta.data

      if(input$SeuratVersion == "Version 4") {
        df1 <- as.data.frame(sc@assays$RNA$counts[rownames(sc@assays$RNA$counts) %in% c(input$Var_to_col_marker,input$Var_to_col_marker2,input$Var_to_col_marker3,"CD4","CD8A","Cd4","Cd8","CD3E","Cd3e"),])
        names(df1) <- colnames(sc@assays$RNA$counts)
        df1 <- as.data.frame(t(df1))
        df1$Cell_Index <- rownames(df1)
        df1
      } else {
        df1 <- as.data.frame(sc@assays$RNA$counts[rownames(sc@assays$RNA$counts) %in% c(input$Var_to_col_marker,input$Var_to_col_marker2,input$Var_to_col_marker3,"CD4","CD8A","Cd4","Cd8","CD3E","Cd3e"),])
        names(df1) <- colnames(sc@assays$RNA$counts)
        df1 <- as.data.frame(t(df1))
        df1$Cell_Index <- rownames(df1)
        df1
      }




    })

    MainTcell_counts_names <- reactive({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )
      req(sc)

      df <- rownames(sc@assays$RNA$scale.data)[rowSums(sc@assays$RNA$scale.data) !=0]

      df <-as.data.frame(df)
      names(df) <- "V1"
      df
    })

    observeEvent(input$load_marker_genes,{
      sc <- (MainTcell_counts_names())
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )
      req(sc)

      updateSelectizeInput(
        session,
        "Var_to_col_marker",
        choices = sc$V1,
        select = "CD8A"
      )
    })

    selected_gene <- reactive({
      MainTcell <- MainTcell_counts()
      validate(
        need(nrow(MainTcell)>0,
             "Upload File")
      )

      MainTcell$Cell_Index <- rownames(MainTcell)
      MainTcell

    })

    meta_data_for_features <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )

      md <- sc@meta.data
      md <- merge(md,selected_gene(),by = "Cell_Index")
      md$log2.count <- log2(as.numeric(md[,names(md) %in% input$Var_to_col_marker]))
      md$log2.count[is.na(md$log2.count)] <- "-Inf"
      md
      md <- subset(md,md$UMAP_1 > input$Filter_lower_UMAP1_marker)
      md <- subset(md,md$UMAP_1 < input$Filter_lower_UMAP1_marker2)

      md <- subset(md,md$UMAP_2 > input$Filter_lower_UMAP2_marker)
      md <- subset(md,md$UMAP_2 < input$Filter_lower_UMAP2_marker2)

      subset(md,md$v_gene_AG != "NA")


    })



    # subset based on UMAP location, subset bashed on norm value for the marker of interest threshold as well. -----

    MainTcell_scale <- reactive({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )
      md <- sc@meta.data
      df1 <- as.data.frame(sc@assays$RNA$scale.data[rownames(sc@assays$RNA$scale.data) %in% c(input$Var_to_col_marker,input$Var_to_col_marker2,input$Var_to_col_marker3,"CD4","CD8A","Cd4","Cd8","CD3E","Cd3e","JCHAIN","Jchain","CD19","Cd19","MS4A1","Ms4a1"),])

      names(df1) <- colnames(sc@assays$RNA$scale.data)
      df1 <- as.data.frame(t(df1))
      df1$Cell_Index <- rownames(df1)
      df1
    })

    selected_scale <- reactive({
      MainTcell <- MainTcell_scale()
      validate(
        need(nrow(MainTcell)>0,
             "Upload File")
      )

      MainTcell$Cell_Index <- rownames(MainTcell)
      MainTcell

    })

    meta_data_for_features_scale <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )

      md <- sc@meta.data
      md <- merge(md,selected_scale(),by = "Cell_Index")
      md$scale <- as.numeric(md[,names(md) %in% input$Var_to_col_marker])

      md <- subset(md,md$UMAP_1 > input$Filter_lower_UMAP1_marker)
      md <- subset(md,md$UMAP_1 < input$Filter_lower_UMAP1_marker2)

      md <- subset(md,md$UMAP_2 > input$Filter_lower_UMAP2_marker)
      md <- subset(md,md$UMAP_2 < input$Filter_lower_UMAP2_marker2)
      subset(md,md$v_gene_AG != "NA")
    })


    # umap scale -------
    marker_selected_UMAP <- reactive({
      umap.meta <- meta_data_for_features_scale()
      validate(
        need(nrow(umap.meta)>0,
             "Upload File")
      )

      umap.meta <- meta_data_for_features_scale()
      umap.meta_pos <- subset(umap.meta,umap.meta$scale != "-Inf")

      umap.meta$scale <- ifelse(umap.meta$scale == "-Inf",NA,umap.meta$scale)
      umap.meta$scale <- as.numeric(umap.meta$scale)
      umap.meta_pos$scale <- as.numeric(umap.meta_pos$scale)

      ggplot(umap.meta,aes(x=UMAP_1,y=UMAP_2,color = scale))+
        geom_point(size = 1) +
        geom_point(data = umap.meta_pos,aes(x=UMAP_1,y=UMAP_2,color = scale),size = 1) +
        # scale_color_continuous(type = "viridis",na.value="grey75") +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
        ) +
        scale_color_distiller(direction = 1, palette = input$col_marker_scale,na.value = "grey85", limits = c(min(0),max(input$max_scale))) +
        theme(
          plot.title = element_text(colour="black",family=input$font_type,size = 24,face = "bold",vjust=.5),
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.title.x = element_text(colour="black",angle=0,vjust=.5,face="plain",family=input$font_type,size = input$title.text.sizer2),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        ) +
        labs(title = input$Var_to_col_marker)

    })

    output$marker_selected_UMAP_plot <- renderPlot({
      umap.meta <- MainTcell_counts_names()
      validate(
        need(nrow(umap.meta)>0,
             "Upload File")
      )
      marker_selected_UMAP()
    })


    output$downloadPlot_marker_selected_UMAP_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Var_to_col_marker,"_marker_selected_UMAP_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_marker_selected_UMAP_plot,height=input$height_marker_selected_UMAP_plot, onefile = FALSE) # open the pdf device

        plot(marker_selected_UMAP())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_marker_selected_UMAP_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Var_to_col_marker,"_marker_selected_UMAP_",x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_marker_selected_UMAP_plot,height = input$height_png_marker_selected_UMAP_plot,res = input$resolution_PNG_marker_selected_UMAP_plot)
        plot(marker_selected_UMAP())
        dev.off()},   contentType = "application/png")

    # scale violin plot ----
    marker_selected_Violine.ridge <- reactive({
      # umap.meta <- meta_data_for_features_scale()
      umap.meta <- meta_data_for_features()
      validate(
        need(nrow(umap.meta)>0,
             "Upload File")
      )

      umap.meta <- meta_data_for_features_scale()

      if (input$select_plot_vio.ridge == "Violin") {
        ggplot(umap.meta,aes(y=scale,x=Sample_Name))+
          geom_violin() +
          geom_jitter()+
          theme(
            legend.position = "none",
          )+
          theme_bw() +
          geom_hline(yintercept = input$cutoff_marker_gt,color = "orange", size = 3)
        # geom_vline(xintercept=neg, linetype="dashed", color = "darkorange")

      }

      else {
        ggplot(umap.meta,aes(x=scale,y=Sample_Name))+
          geom_density_ridges() +
          theme_ridges() +
          geom_vline(xintercept = input$cutoff_marker_gt,color = "orange",size = 3)
      }


    })


    output$marker_selected_VioRidge_plot <- renderPlot({
      umap.meta <- MainTcell_counts_names()
      validate(
        need(nrow(umap.meta)>0,
             "Upload File")
      )
      marker_selected_Violine.ridge()
    })


    filtered_positive_marker_TCRsum <- reactive({
      # umap.meta <- meta_data_for_features_scale()
      umap.meta <- meta_data_for_features()
      validate(
        need(nrow(umap.meta)>0,
             "Upload File")
      )

      umap.meta$log2.count <- ifelse(umap.meta$log2.count == "-Inf",-2,umap.meta$log2.count)
      umap.meta$log2.count <- as.numeric(umap.meta$log2.count)

      umap.meta_marker_pos <- subset(umap.meta,umap.meta$log2.count > input$cutoff_marker_gt )
      umap.meta_marker_pos$cloneCount <- 1
      umap.meta_marker_pos_TCR <- as.data.frame(umap.meta_marker_pos[,names(umap.meta_marker_pos) %in% c(input$V_gene_sc,"cloneCount")])
      file.names <- input$V_gene_sc

      df_unique_sum <- ddply(umap.meta_marker_pos_TCR,file.names ,numcolwise(sum))
      df_unique_sum$type <- "Pos"
      df_unique_sum[order(df_unique_sum$cloneCount, decreasing = T),]
    })

    filtered_negative_marker_TCRsum <- reactive({
      umap.meta <- meta_data_for_features()
      validate(
        need(nrow(umap.meta)>0,
             "Upload File")
      )

      umap.meta$log2.count <- ifelse(umap.meta$log2.count == "-Inf",-2,umap.meta$log2.count)
      umap.meta$log2.count <- as.numeric(umap.meta$log2.count)

      umap.meta_marker_neg <- subset(umap.meta,umap.meta$log2.count < input$cutoff_marker_gt )
      umap.meta_marker_neg$cloneCount <- 1
      umap.meta_marker_neg_TCR <- as.data.frame(umap.meta_marker_neg[,names(umap.meta_marker_neg) %in% c(input$V_gene_sc,"cloneCount")])
      file.names <- input$V_gene_sc

      df_unique_sum <- ddply(umap.meta_marker_neg_TCR,file.names ,numcolwise(sum))
      df_unique_sum$type <- "Neg"
      df_unique_sum[order(df_unique_sum$cloneCount, decreasing = T),]
    })

    output$TCR_marker_positive_count <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1,2,5,10,20,50,100), pageLength = 1, scrollX = TRUE),{
      umap.meta <- filtered_positive_marker_TCRsum()
      umap.meta

    })

    output$TCR_marker_neg_count <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1,2,5,10,20,50,100), pageLength = 1, scrollX = TRUE),{
      umap.meta <- filtered_negative_marker_TCRsum()
      umap.meta

    })

    # histogram of counts marker and TCR -----

    merged_marker_hist <- reactive({
      a <- filtered_negative_marker_TCRsum()
      a$expand <- ifelse(a$cloneCount>=2,"Ex","NEx")
      names(a)[2:4] <- paste(names(a)[2:4],"_neg",sep = "")

      b <- filtered_positive_marker_TCRsum()
      b$expand <- ifelse(b$cloneCount>=2,"Ex","NEx")
      names(b)[2:4] <- paste(names(b)[2:4],"_pos",sep = "")

      by_ab <- names(a)[1]
      ab <-  merge(b,a, by =by_ab,all = T)



      ab[is.na(ab)] <- 0

      ab$expand_pos <- ifelse(ab$expand_pos!=0,ab$expand_pos,"Abs")
      ab$expand_neg <- ifelse(ab$expand_neg!=0,ab$expand_neg,"Abs")
      ab <- ab[,!names(ab) %in% c("type_neg","type_pos")]
      ab$diff <- ab$cloneCount_pos -  ab$cloneCount_neg
      ab$ratio <- ab$cloneCount_pos/ab$cloneCount_neg
      ab[order(ab$diff,decreasing = T),]

    })

    output$merged_marker_hist_table <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1,2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      umap.meta <- merged_marker_hist()
      umap.meta
    })

    marker_histogram <- reactive({
      a <- filtered_negative_marker_TCRsum()
      b <- filtered_positive_marker_TCRsum()
      ab <- rbind(a,b)

      # ab <- subset(ab,ab$cloneCount>2)

      ab$log10_ab <- log10(ab$cloneCount)
      ggplot(ab, aes(x=cloneCount, color=type)) +
        geom_histogram(binwidth=0.1,fill="white", alpha=0.5, position="identity") +
        theme_classic() +
        # geom_density(alpha=.2) +
        theme(legend.position="top") +
        facet_wrap(~type)
    })

    output$marker_selected_histogram_plot <- renderPlot({
      umap.meta <- MainTcell_counts_names()
      validate(
        need(nrow(umap.meta)>0,
             "Upload File")
      )
      marker_histogram()
    })

    stats_marker_one <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )

      md <- sc@meta.data
      # md <- merge(md,selected_gene(),by = "Cell_Index",all.x = T, sort = F)
      md <- merge(md,selected_scale(),by = "Cell_Index",all.x = T, sort = F)
      rownames(md) <- md$Cell_Index
      sc@meta.data <- md


      sc@meta.data$scale <- (as.numeric(md[,names(md) %in% input$Var_to_col_marker]))

      sc <- subset(sc,UMAP_1 > input$Filter_lower_UMAP1_marker)
      sc <- subset(sc,UMAP_1 < input$Filter_lower_UMAP1_marker2)

      sc <- subset(sc,UMAP_2 > input$Filter_lower_UMAP2_marker)
      sc <- subset(sc,UMAP_2 < input$Filter_lower_UMAP2_marker2)
      sc <- subset(sc,v_gene_AG != "NA")

      sc@meta.data$pos_neg <- ifelse(sc@meta.data$scale > input$cutoff_marker_gt, "pos","neg")
      df <- filtered_positive_marker_TCRsum()

      df2 <- subset(df,df$cloneCount>input$pos_expanded_cut_off)

      list.gene <- unique(df2[,names(df2) %in% input$V_gene_sc])
      sc@meta.data$Selected_Gene <-  sc@meta.data[,names(sc@meta.data) %in% input$V_gene_sc]
      # sc@meta.data$in.list <- ifelse(sc@meta.data$Selected_Gene %in% list.gene & sc@meta.data$scale > input$cutoff_marker_gt,"In.Pos.list","other")
      sc@meta.data$Pos_Exp <- ifelse(sc@meta.data$Selected_Gene %in%  list.gene & sc@meta.data$pos_neg == "pos", "Ex_pos","other")
      # ifelse(sc@meta.data$Selected_Gene %in%  list.gene & sc@meta.data$pos_neg == "neg", "Ex_neg","other"
      # ))
      sc

    })

    output$marker_selected_tab <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1,2,5,10,20,50,100), pageLength = 20, scrollX = TRUE),{
      umap.meta <- meta_data_for_features()
      umap.meta
      # stats_marker_one()
      # meta_data_for_features_scale()
    })


    compare.stat_marker <- reactive({
      sc <-stats_marker_one()

      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )


      if (input$Analysis_marker_stats_type == "Expanded") {
        Idents(object = sc) <- sc@meta.data$Pos_Exp
      }

      else {
        Idents(object = sc) <- sc@meta.data$pos_neg

      }

      min.pct.expression<- input$min_point_ #standard setting: 0.25
      min.logfc<-  input$LogFC_Marker #0.25 is standard
      # p.val.cutoff <-  input$pval_top #(1/10^3) is standard, use (1/10^0) to ignore

      cluster.names <- unique(Idents(sc))[order(unique(Idents(sc)))]
      # print(paste0("calculating markers for cluster ",name.clone,". Total: ",length(cluster.names)," clusters"))

      if (input$Analysis_marker_stats_type == "Expanded") {
        markers.fm.list <- FindMarkers(sc, ident.1 = "Ex_pos", min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
      }

      else {
        markers.fm.list <- FindMarkers(sc, ident.1 = "pos", min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
      }

      markers.fm.list2 <- subset(markers.fm.list,markers.fm.list$p_val_adj < input$pval.ex.filter)
      as.data.frame(markers.fm.list2)

    })

    output$Compare.stat_marker <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength =10, scrollX = TRUE),{

      sc <-input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )
      compare.stat_marker()

    })

    # download tables ----

    output$downloaddf_clonotype_distribution <- downloadHandler(
      filename = function(){
        x = today()
        paste(input$Var_to_col_marker,"_clonotypes_", input$pos_expanded_cut_off,"_stats",x,".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(merged_marker_hist())
        write.csv(df,file, row.names = T)
      } )

    output$downloaddf_Marker_stats <- downloadHandler(
      filename = function(){
        x = today()
        paste(input$Var_to_col_marker,"_Exp_cuttoff", input$pos_expanded_cut_off,"_stats",x,".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(compare.stat_marker())
        write.csv(df,file, row.names = T)
      } )

    # display T cells assoicated with that marker. ----

    # Dual marker analysis -----

    observeEvent(input$load_marker_genes,{
      sc <- (MainTcell_counts_names())
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )
      req(sc)

      updateSelectizeInput(
        session,
        "Var_to_col_marker2",
        choices = sc$V1,
        select = "CD8B"
      )
    })

    observeEvent(input$load_marker_genes,{
      sc <- (MainTcell_counts_names())
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )
      req(sc)

      updateSelectizeInput(
        session,
        "Var_to_col_marker3",
        choices = sc$V1,
        select = "CD4"
      )
    })



    meta_data_for_features_scale2 <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )

      md <- sc@meta.data
      md <- merge(md,selected_scale(),by = "Cell_Index")
      md$scale2 <- as.numeric(md[,names(md) %in% input$Var_to_col_marker2])
      md$scale3 <- as.numeric(md[,names(md) %in% input$Var_to_col_marker3])

      md <- subset(md,md$UMAP_1 > input$Filter_dual_UMAP1_marker)
      md <- subset(md,md$UMAP_1 < input$Filter_dual_UMAP1_marker2)

      md <- subset(md,md$UMAP_2 > input$Filter_dual_UMAP2_marker)
      md <- subset(md,md$UMAP_2 < input$Filter_dual_UMAP2_marker2)
      subset(md,md$v_gene_AG != "NA")
    })


    output$meta_data_for_features_scale2_df <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1,2,5,10,20,50,100), pageLength = 1, scrollX = TRUE),{
      umap.meta <- meta_data_for_features_scale2()
      umap.meta

    })


    # umap dual scale -------
    marker_selected_UMAP_scale2 <- reactive({
      umap.meta <- meta_data_for_features_scale2()
      validate(
        need(nrow(umap.meta)>0,
             "Upload File")
      )
      umap.meta <- meta_data_for_features_scale2()
      umap.meta_pos <- subset(umap.meta,umap.meta$scale2 != "-Inf")

      umap.meta$scale2 <- ifelse(umap.meta$scale2 == "-Inf",NA,umap.meta$scale2)
      umap.meta$scale2 <- as.numeric(umap.meta$scale2)
      umap.meta_pos$scale2 <- as.numeric(umap.meta_pos$scale2)

      ggplot(umap.meta,aes(x=UMAP_1,y=UMAP_2,color = scale2))+
        geom_point(size = 1) +
        geom_point(data = umap.meta_pos,aes(x=UMAP_1,y=UMAP_2,color = scale2),size = 1) +
        # scale_color_continuous(type = "viridis",na.value="grey75") +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
        ) +
        scale_color_distiller(direction = 1, palette = input$col_marker_scale,na.value = "grey85", limits = c(min(0),max(input$max_scale2))) +
        theme(
          plot.title = element_text(colour="black",family=input$font_type,size = 24,face = "bold",vjust=.5),
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.title.x = element_text(colour="black",angle=0,vjust=.5,face="plain",family=input$font_type,size = input$title.text.sizer2),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        ) +
        labs(title = input$Var_to_col_marker2)
    })

    marker_selected_UMAP_scale3 <- reactive({
      umap.meta <- meta_data_for_features_scale2()
      validate(
        need(nrow(umap.meta)>0,
             "Upload File")
      )

      umap.meta <- meta_data_for_features_scale2()
      umap.meta_pos <- subset(umap.meta,umap.meta$scale3 != "-Inf")

      umap.meta$scale3 <- ifelse(umap.meta$scale3 == "-Inf",NA,umap.meta$scale2)
      umap.meta$scale3 <- as.numeric(umap.meta$scale3)
      umap.meta_pos$scale3 <- as.numeric(umap.meta_pos$scale3)

      ggplot(umap.meta,aes(x=UMAP_1,y=UMAP_2,color = scale3))+
        geom_point(size = 1) +
        geom_point(data = umap.meta_pos,aes(x=UMAP_1,y=UMAP_2,color = scale3),size = 1) +
        # scale_color_continuous(type = "viridis",na.value="grey75") +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
        ) +
        scale_color_distiller(direction = 1, palette = input$col_marker_scale,na.value = "grey85", limits = c(min(0),max(input$max_scale3))) +
        theme(
          plot.title = element_text(colour="black",family=input$font_type,size = 24,face = "bold",vjust=.5),
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.title.x = element_text(colour="black",angle=0,vjust=.5,face="plain",family=input$font_type,size = input$title.text.sizer2),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        ) +
        labs(title = input$Var_to_col_marker3)

    })

    output$marker_selected_UMAP_plot2 <- renderPlot({
      umap.meta <- MainTcell_counts_names()
      validate(
        need(nrow(umap.meta)>0,
             "Upload File")
      )
      marker_selected_UMAP_scale2()
    })
    output$marker_selected_UMAP_plot3 <- renderPlot({
      umap.meta <- MainTcell_counts_names()
      validate(
        need(nrow(umap.meta)>0,
             "Upload File")
      )
      marker_selected_UMAP_scale3()
    })

    # umap dual scale -------
    df_dotplot_marker <- reactive({

      umap.meta <- meta_data_for_features_scale2()

      validate(
        need(nrow(umap.meta)>0,
             "Upload File")
      )

      selected.col <- umap.meta[,names(umap.meta) %in% input$Colour_By_this]

      umap.meta$selected <- umap.meta[,names(umap.meta) %in% input$Colour_By_this]
      umap.meta$selected <- as.character(umap.meta$selected)
      umap.meta$selected <- gsub("NA",NA,umap.meta$selected)

      len.colour <- length(unique(selected.col))
      colour_markersby <- rainbow(len.colour)

      umap.meta$markerX <- umap.meta[,names(umap.meta) %in% input$Var_to_col_marker2]
      umap.meta$markerY <- umap.meta[,names(umap.meta) %in% input$Var_to_col_marker3]

      Q1 <- umap.meta[umap.meta$markerX<input$X_axis_dot_dual & umap.meta$markerY> input$Y_axis_dot_dual,]
      Q2 <- umap.meta[umap.meta$markerX>input$X_axis_dot_dual & umap.meta$markerY> input$Y_axis_dot_dual,]
      Q3 <- umap.meta[umap.meta$markerX<input$X_axis_dot_dual & umap.meta$markerY< input$Y_axis_dot_dual,]
      Q4 <- umap.meta[umap.meta$markerX>input$X_axis_dot_dual & umap.meta$markerY< input$Y_axis_dot_dual,]

      Q1.per <- dim(Q1)[1]/dim(umap.meta)[1]*100
      Q2.per <- dim(Q2)[1]/dim(umap.meta)[1]*100
      Q3.per <- dim(Q3)[1]/dim(umap.meta)[1]*100
      Q4.per <- dim(Q4)[1]/dim(umap.meta)[1]*100

      annotations <- data.frame(
        xpos = c(-Inf,-Inf,Inf,Inf),
        ypos =  c(-Inf, Inf,-Inf,Inf),
        annotateText = c(paste0("Q3: ",round(Q3.per,2),"%"),paste0("Q1: ",round(Q1.per,2),"%"),
                         paste0("Q4: ",round(Q4.per,2),"%"),paste0("Q2: ",round(Q2.per,2),"%")),
        hjustvar = c(-0.1,-0.1,1.1,1.2) ,
        vjustvar = c(-0.4,1.4,-0.4,1.4)) #<- adjust

      plot <-  ggplot() +
        geom_point(data=umap.meta,aes(x=get(input$Var_to_col_marker2),y=get(input$Var_to_col_marker3), colour = selected)) +
        theme_bw() +
        geom_hline(yintercept=input$Y_axis_dot_dual) +
        geom_vline(xintercept =input$X_axis_dot_dual) +
        xlab(input$Var_to_col_marker2) +
        ylab(input$Var_to_col_marker3) +
        scale_colour_manual(values = colour_markersby, na.value = "grey85") +
        geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), colour="black",family=input$font_type,size = input$anno_text_size) +
        theme(
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=0),
          axis.title.x = element_text(colour="black",angle=0,vjust=.5,face="plain",family=input$font_type,size = input$title.text.sizer2),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )

      plot

      # ggExtra::ggMarginal(plot,groupColour = TRUE, groupFill = TRUE)

      # ggplot(umap.meta.df,aes(x=scale2,y=scale3)) +
      #   geom_point() +
      #   theme_bw() +
      #   geom_hline(yintercept=-1) +
      #   geom_vline(xintercept = 0)

    })

    output$df_dotplot_marker_plot <- renderPlot({
      umap.meta <- MainTcell_counts_names()
      validate(
        need(nrow(umap.meta)>0,
             "Upload File")
      )
      df_dotplot_marker()
    })

    output$downloadPlot_df_dotplot_marker_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("Dual_marker_plot_",input$Var_to_col_marker2,"_",input$Var_to_col_marker3,"_",x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_df_dotplot_marker_plot,height=input$height_df_dotplot_marker_plot, onefile = FALSE) # open the pdf device
        plot(df_dotplot_marker())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_df_dotplot_marker_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("Dual_marker_plot_",input$Var_to_col_marker2,"_",input$Var_to_col_marker3,"_",x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_df_dotplot_marker_plot,
            height = input$height_png_df_dotplot_marker_plot,
            res = input$resolution_PNG_df_dotplot_marker_plot)
        plot(df_dotplot_marker())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    # summary of TCR per quad -----
    dual_maker_TCR_Sum <- reactive({

      umap.meta <- meta_data_for_features_scale2()

      validate(
        need(nrow(umap.meta)>0,
             "Upload File")
      )

      umap.meta$markerX <- umap.meta[,names(umap.meta) %in% input$Var_to_col_marker2]
      umap.meta$markerY <- umap.meta[,names(umap.meta) %in% input$Var_to_col_marker3]

      blank <- as.data.frame(names(umap.meta))
      names(blank) <- "V1"
      rownames(blank) <- blank$V1
      blank$blank1 <- "blank"
      blank$blank2 <- "blank"
      blank$blank3 <- "blank"
      blank$blank4 <- "blank"
      blank <- blank[,-c(1)]
      blank <- as.data.frame(t(blank))
      head(blank)

      list1 <- c(input$Y_axis_dot_dual-1,input$Y_axis_dot_dual-1,input$Y_axis_dot_dual+1,input$Y_axis_dot_dual+1)
      blank$markerY <- list1

      list2 <- c(input$X_axis_dot_dual-1,input$X_axis_dot_dual+1,input$X_axis_dot_dual-1,input$X_axis_dot_dual+1)
      blank$markerX <- list2
      blank

      umap.meta <- rbind(umap.meta,blank)
      umap.meta

      Q1 <- umap.meta[umap.meta$markerX<input$X_axis_dot_dual & umap.meta$markerY> input$Y_axis_dot_dual,]
      Q1$Q1 <- 1
      Q2 <- umap.meta[umap.meta$markerX>input$X_axis_dot_dual & umap.meta$markerY> input$Y_axis_dot_dual,]
      Q2$Q2 <- 1
      Q3 <- umap.meta[umap.meta$markerX<input$X_axis_dot_dual & umap.meta$markerY< input$Y_axis_dot_dual,]
      Q3$Q3 <- 1
      Q4 <- umap.meta[umap.meta$markerX>input$X_axis_dot_dual & umap.meta$markerY< input$Y_axis_dot_dual,]
      Q4$Q4 <- 1

      # as.data.frame(is.data.frame(subset(umap.meta, markerX>input$X_axis_dot_dual & markerY> input$Y_axis_dot_dual)))


      Q1_quad <- Q1[,names(Q1) %in% c(input$V_gene_sc,"Q1")]
      Q1_quad <- ddply(Q1_quad,input$V_gene_sc ,numcolwise(sum))

      Q2_quad <- Q2[,names(Q2) %in% c(input$V_gene_sc,"Q2")]
      Q2_quad <- ddply(Q2_quad,input$V_gene_sc ,numcolwise(sum))

      Q3_quad <- Q3[,names(Q3) %in% c(input$V_gene_sc,"Q3")]
      Q3_quad <- ddply(Q3_quad,input$V_gene_sc ,numcolwise(sum))

      Q4_quad <- Q4[,names(Q4) %in% c(input$V_gene_sc,"Q4")]
      Q4_quad <- ddply(Q4_quad,input$V_gene_sc ,numcolwise(sum))

      Q1_Q2 <- merge(Q1_quad,Q2_quad,by = input$V_gene_sc,all=T)
      Q3_Q4 <- merge(Q3_quad,Q4_quad,by = input$V_gene_sc,all=T)
      umeta.TCR <-  merge(Q1_Q2,Q3_Q4,by = input$V_gene_sc,all=T)

      umeta.TCR[is.na(umeta.TCR)] <- 0
      umeta.TCR_needed.names <- names(umeta.TCR)
      umeta.TCR$total <- rowSums(umeta.TCR[2:5])
      rownames(umeta.TCR) <- umeta.TCR[,c(1)]
      umeta.TCR <- umeta.TCR[!(row.names(umeta.TCR) %in% c("blank")), ]
      umeta.TCR <- umeta.TCR[,-c(1)]
      names(umeta.TCR) <- c(paste0(input$Var_to_col_marker3,"pos",input$Var_to_col_marker2,"neg"),
                            paste0(input$Var_to_col_marker3,"pos",input$Var_to_col_marker2,"pos"),
                            paste0(input$Var_to_col_marker3,"neg",input$Var_to_col_marker2,"neg"),
                            paste0(input$Var_to_col_marker3,"neg",input$Var_to_col_marker2,"pos"),
                            "total"
      )

      umeta.TCR
    })


    output$dual_maker_TCR_Sum_DT <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1,2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      umeta.TCR <- dual_maker_TCR_Sum()
      umeta.TCR

    })

    output$Dule_marker_TCRsummary_DT <- downloadHandler(
      filename = function(){
        x = today()
        paste("Dule_marker_TCRsummary_",input$Var_to_col_marker2,"_",input$Var_to_col_marker3,"_",x,".csv", sep = "")
      },
      content = function(file){
        df <- dual_maker_TCR_Sum()
        write.csv(df,file)
      } )


    # Add ridge plot for the distribution...

    # updating the UI for prioritisation -------
    output$Top_clone_number <- renderUI({
      sc <- UMAP_metadata_with_labs()

      validate(
        need(nrow(sc)>0,
             "upload file")
      )
      TCR_Expanded_Df <- TCR_Expanded_fun(sc,(input$Samp_col),(input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      TCR_Expanded_Df2 <-  subset(TCR_Expanded_Df,TCR_Expanded_Df$percent>input$cut.off_percent_rep)
      observations <- sum(TCR_Expanded_Df2$obs)

      if (length(observations)>0){
        column(6,numericInput("top_no_clonotypes","Top clonotypes per group",value = observations,step = 1, min = 0, max = 20))
      } else {
        # message("No immunodom present")
        observations <- 1
        column(6,numericInput("top_no_clonotypes","Top clonotypes per group",value = observations,step = 1, min = 0, max = 20))
      }

    })
    output$cut.off_expanded2 <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )

      req(input$Samp_col,input$V_gene_sc)
      TCR_Expanded_Df <- TCR_Expanded_fun(sc,(input$Samp_col),(input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1

      total.condition <- ddply(TCR_Expanded_Df,c("samp.count"),numcolwise(sum))
      total.condition$count_total <- total.condition$samp.count * total.condition$obs
      total.condition$frequency <- total.condition$count_total/sum(total.condition$count_total)
      total.condition <- total.condition %>% mutate(cum_freq = cumsum(frequency))
      total.condition <- subset(total.condition,total.condition$cum_freq>input$cutoff.expanded)
      num.cut.off <- total.condition$samp.count[1]

      # df3.meta <- df3.meta[!grepl("RNA",df3.meta) & !grepl("BCR",df3.meta) & !grepl("TCR",df3.meta)& !grepl("_gene",df3.meta) & !grepl("allele",df3.meta) & !grepl("percent",df3.meta) & !grepl("cdr3",df3.meta)]
      numericInput("cut.off_expanded","Cut off greater than", value = num.cut.off, step = 1, min = 1)

    })
    output$Expanded.dotplot.cutoffs <- renderUI({

      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "upload file")
      )

      TCR_Expanded_Df <- TCR_Expanded_fun(sc,(input$Samp_col),(input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      TCR_Expanded_Df2 <-  subset(TCR_Expanded_Df,TCR_Expanded_Df$percent>input$cut.off_percent_rep)
      observations <- sum(TCR_Expanded_Df2$obs)
      percentage <- sum(TCR_Expanded_Df2$percent)
      length.samp.ID <- length(unique(TCR_Expanded_Df$ID_Column))
      # print(length.samp.ID)
      mat <- mat_sum(sc,input$Samp_col,input$V_gene_sc)

      if(observations>0) {
        fluidRow(
          column(6,selectInput("restrict.dotpot","Restrict to top list",choices = c("no","yes"),selected = "no")),
          column(6,numericInput("restrict.dotpot.num","Total genes to display:", value = 40)),

        )
      } else {

        sc <- UMAP_metadata_with_labs()
        validate(
          need(nrow(sc)>0,
               error_message_val_UMAP)
        )

        markers.fm.list2 <- Vals_expanded.stats()
        req(markers.fm.list2)
        if (length(rownames(markers.fm.list2))>40) {
          selected.restrict.dotpot <- "yes"
          value.restrict.dotpot.num <- 40

          fluidRow(
            column(6,selectInput("restrict.dotpot","Restrict to top list",choices = c("no","yes"),selected = selected.restrict.dotpot)),
            column(6,numericInput("restrict.dotpot.num","Total genes to display:", value = value.restrict.dotpot.num)),

          )

        }

        else {
          fluidRow(
            column(6,selectInput("restrict.dotpot","Restrict to top list",choices = c("no","yes"),selected = "no")),
            column(6,numericInput("restrict.dotpot.num","Total genes to display:", value = 40)),
          )
        }

      }

      # df3.meta <- df3.meta[!grepl("RNA",df3.meta) & !grepl("BCR",df3.meta) & !grepl("TCR",df3.meta)& !grepl("_gene",df3.meta) & !grepl("allele",df3.meta) & !grepl("percent",df3.meta) & !grepl("cdr3",df3.meta)]



    })

    ## prioritising and automating the analysis ------

    output$Simple_workflow_step1 <- renderPrint({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "upload file")
      )

      TCR_Expanded_Df <- TCR_Expanded_fun(sc,(input$Samp_col),(input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      length.samp.ID <- length(unique(TCR_Expanded_Df$ID_Column))
      # print(length.samp.ID)
      mat <- mat_sum(sc,input$Samp_col,input$V_gene_sc)
      if(max(mat$TotalSamps)==1 && length.samp.ID ==1) {
        print("one individual and one sample")

        TCR_Expanded_Df <- TCR_Expanded_fun(sc,(input$Samp_col),(input$V_gene_sc))
        TCR_Expanded_Df$obs <- 1
        TCR_Expanded_Df2 <-  subset(TCR_Expanded_Df,TCR_Expanded_Df$percent>input$cut.off_percent_rep)
        observations <- sum(TCR_Expanded_Df2$obs)
        percentage <- sum(TCR_Expanded_Df2$percent)
        if(observations>0) {
          print("ImmunoDom")
          # observations <-  Frequency_expanded_df$obs[Frequency_expanded_df$percent>10 & Frequency_expanded_df$Frequency_expanded %in% "5. Gigantic (0.1 > X <= 0.5)"]
          # percent <-  Frequency_expanded_df$percent[Frequency_expanded_df$percent>10 & Frequency_expanded_df$Frequency_expanded %in% "5. Gigantic (0.1 > X <= 0.5)"]

          print(paste0("There are ",print(observations)," immuno dominant clonotype(s) that account for ",round(percentage,2),"% of the repertoire"))
          # print(paste0("This process will download (1) UMAP count and top plots, bar plot of immunodominat (>10% of repertoire)"))

          # cat("test")
        } else {
          print("Polyclonal")
        }

      } else if (max(mat$TotalSamps)>1 || length.samp.ID >1) {
        print("multiple individuals or samples")
      } else {

        print("other")
      }

    })

    ### UI outputs -----
    output$Module_case_statements <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "upload file")
      )

      TCR_Expanded_Df <- TCR_Expanded_fun(sc,(input$Samp_col),(input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      TCR_Expanded_Df2 <-  subset(TCR_Expanded_Df,TCR_Expanded_Df$percent>input$cut.off_percent_rep)
      observations <- sum(TCR_Expanded_Df2$obs)
      length.samp.ID <- length(unique(TCR_Expanded_Df$ID_Column))
      print(length.samp.ID)

      mat <- mat_sum(sc,input$Samp_col,input$V_gene_sc)

      if(max(mat$TotalSamps)==1 && length.samp.ID ==1) {
        print("one individual and one sample")
        if(observations>0) {
          fluidRow(
            column(12,actionButton("ImmDom_download_buttonOneOne","Download ImmunoDom (1 & 1) analysis"))
          )
        } else {
          fluidRow(
            # column(12,selectInput("AnalysisType","Preset parameters for", = c("ImmunoDom"))),
            column(12,actionButton("Poly_download_buttonOneOne","Download Polyclonal (1 & 1) analysis"))
          )
        }
      } else if (max(mat$TotalSamps)>1 || length.samp.ID >1) {

        BD_sum <- Top_clonotypes_multiCounts()

        if(dim(BD_sum)[1]>10) {
          BD_sum <- BD_sum[1:11,]
        }


        fluidRow(
          column(6,numericInput("cut.off_percent_repMulti","Priority cut-off",value = 1,step = 0.001, min = 0, max = 1)),
          column(12,actionButton("Multi_download_button","Download multi analysis"))
        )
      } else {
        # print("multiple individuals and multiple samples")
      }
    })

    output$Default_priority_cutoffAG <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )
      md <- sc@meta.data
      x = today()
      if (length(input.data_sc_clusTCR_AG())>0) {
        df7 <- AG_cluster()

        df7 <- df7[order(df7$priority,decreasing = F),]

        if (max(df7$Updated_order)>10) {
          prior <- subset(df7,df7$Updated_order > 11)
          numericInput("priority_cutoff","Priority cut-off (AG)",value = min(prior$priority),step = 0.01, min = 0,max = 1)
        } else {
          prior <- df7
          numericInput("priority_cutoff","Priority cut-off (AG)",value = 1,step = 0.01, min = 0,max = 1)
        }
      }

      else {
        # numericInput("priority_cutoff","Priority cut-off (AG)",value = 1,step = 0.01, min = 0)
      }
    })
    output$Default_priority_cutoffBD <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )
      md <- sc@meta.data
      x = today()
      if (length(input.data_sc_clusTCR_BD())>0) {
        df7 <- BD_cluster()

        df7 <- df7[order(df7$priority,decreasing = F),]
        if (max(df7$Updated_order)>10) {
          prior <- subset(df7,df7$Updated_order > 11)
          numericInput("priority_cutoffBD","Priority cut-off (BD)",value = min(prior$priority),step = 0.01, min = 0,max = 1)
        } else {
          prior <- df7
          numericInput("priority_cutoffBD","Priority cut-off (BD)",value = 0.25,step = 0.01, min = 0,max = 1)
        }
      } else {
      }

    })

    # one sample one indiv ImmunoDom -----
    observeEvent(input$ImmDom_download_buttonOneOne,{
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "upload file")
      )

      x = today()
      message("Downloading the Summary table...")
      top.name.clonotypes <- paste("Prioritisation/ImmunoDom/Expansion_summary_table_",x,".csv",sep="")
      write.csv(Top_clonotype_df2(),top.name.clonotypes, row.names = F)

      message("Downloading the top UMAP...")
      top.name.clonotypes.count_png <- paste("Prioritisation/ImmunoDom/Expansion_UMAP_top_",x,".png",sep="")
      png(top.name.clonotypes.count_png, width = input$width_png_TCR.UMAP_top,height = input$height_png_TCR.UMAP_top,res = input$resolution_PNG_TCR.UMAP_top)
      plot(Topclonotypes())
      dev.off()

      message("Downloading the count UMAP...")
      top.name.clonotypes.top_png <- paste("Prioritisation/ImmunoDom/Expansion_UMAP_count_",x,".png",sep="")
      png(top.name.clonotypes.top_png, width = input$width_png_TCR.UMAP,height = input$height_png_TCR.UMAP,res = input$resolution_PNG_TCR.UMAP)
      plot(UMAP.TCRclonalit2())
      dev.off()

      message("Downloading the Dom bar plot...")
      top_clonotype_bar_code_immdom()

      message("Downloading Dom stats files and dot plot...")
      top_clone_FindMaker_looped()

      message("Downloading AG cluster table")

      top.name.clonotypes <- paste("Prioritisation/ImmunoDom/Cluster_summary_table_AG",x,".csv",sep="")
      write.csv(clusTCR2_df(),top.name.clonotypes, row.names = F)

    })
    top_clonotype_bar_code_immdom <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload")
      )
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$V_gene_sc]

      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload")
      )

      req(input$Samp_col,input$V_gene_sc)
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$V_gene_sc]

      BD_sum <- Top_clonotype_df2()
      BD_sum <- subset(BD_sum,BD_sum$Total_count>1)

      TCR_Expanded_Df <- TCR_Expanded_fun(sc,(input$Samp_col),(input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      TCR_Expanded_Df2 <-  subset(TCR_Expanded_Df,TCR_Expanded_Df$percent>input$cut.off_percent_rep)
      observations <- sum(TCR_Expanded_Df2$obs)
      percentage <- sum(TCR_Expanded_Df2$percent)

      withProgress(message = 'Performing ImmDom Analysis', value = 0, {
        for (i in 1:observations) {

          incProgress(1/observations, detail = paste("Clone", i,"of",observations))
          message(BD_sum$cluster_name[i])
          name.clone <- BD_sum$cluster_name[i]
          top_BD_clonotype <- df3.meta[df3.meta$cluster_name %in% name.clone,]
          # print(top_BD_clonotype)

          dtop_clonotype_bar_code <- top_BD_clonotype

          # req(input$Graph_split_order)

          dtop_clonotype_bar_code$Selected_group <- dtop_clonotype_bar_code[,names(dtop_clonotype_bar_code) %in% input$Split_group_by_]
          num <- 1
          # num <- as.data.frame(num[complete.cases(num)==T,])
          as.data.frame(length(num))
          if (input$colourtype == "default") {
            colorblind_vector <- gg_fill_hue(num)
          } else if (input$colourtype == "hcl.colors") {
            colorblind_vector <- c(hcl.colors(num, palette = "viridis"))
          } else if (input$colourtype == "topo.colors") {
            colorblind_vector <- c(topo.colors(num))
          } else if (input$colourtype == "heat.colors") {
            colorblind_vector <- c(heat.colors(num))
          } else if (input$colourtype == "terrain.colors") {
            colorblind_vector <- c(terrain.colors(num))
          } else if (input$colourtype == "rainbow") {
            colorblind_vector <- c(rainbow((num)))
          } else if (input$colourtype == "random") {
            colorblind_vector <- distinctColorPalette(num)

          }  else {

          }


          colorblind_vector <- as.data.frame(colorblind_vector)
          names(colorblind_vector) <- "cols"

          dtop_clonotype_bar_code$Selected_chain2 <- dtop_clonotype_bar_code[,names(dtop_clonotype_bar_code) %in% input$V_gene_sc]
          dtop_clonotype_bar_code$Selected_chain3 <- gsub("_"," ",dtop_clonotype_bar_code$Selected_chain2)
          dtop_clonotype_bar_code$Selected_chain3 <- gsub("[.]"," ",dtop_clonotype_bar_code$Selected_chain3)

          dtop_clonotype_bar_code <- dtop_clonotype_bar_code[dtop_clonotype_bar_code$Selected_group %in% input$Graph_split_order,]
          dtop_clonotype_bar_code$Selected_group <- factor(dtop_clonotype_bar_code$Selected_group,levels = input$Graph_split_order)

          ggplot_plot <- ggplot(dtop_clonotype_bar_code, aes(x=Selected_group, fill=Selected_chain3,colour = Selected_chain3, label = Selected_chain3)) +
            geom_bar() +
            theme_bw()+
            scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector$cols, na.value=input$NA_col_analysis)+
            scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector$cols, na.value=input$NA_col_analysis)+
            # scale_alpha_manual(values = rep(1,length(unique(dtop_clonotype_bar_code$Selected_chain))), na.value=0.5)+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
            theme(
              axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
              axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
              axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=90),
              axis.title.x = element_blank(),
              legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.position = input$legend_position,
              legend.title = element_blank()
            ) +
            guides(color = "none", size = "none")

          x = today()
          top.name.clonotypes.top_png <- paste("Prioritisation/ImmunoDom/",i,"_top_clone_",gsub("[/]","",gsub("&","",name.clone)),"_",x,".png",sep="")

          num_width <- length(unique(dtop_clonotype_bar_code$Selected_group))

          png(top.name.clonotypes.top_png, width = (num_width*100+400),height = input$height_png_TCR.UMAP,res = input$resolution_PNG_TCR.UMAP)
          plot(ggplot_plot)
          dev.off()

        }
      })
    })
    top_clone_FindMaker_looped <- reactive({

      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload")
      )
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$V_gene_sc]
      sc@meta.data$Vgene <- sc@meta.data[,names(sc@meta.data) %in% input$V_gene_sc]

      BD_sum <- Top_clonotype_df2()
      BD_sum <- subset(BD_sum,BD_sum$Total_count>1)
      TCR_Expanded_Df <- TCR_Expanded_fun(sc,(input$Samp_col),(input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      TCR_Expanded_Df2 <-  subset(TCR_Expanded_Df,TCR_Expanded_Df$percent>input$cut.off_percent_rep)
      observations <- sum(TCR_Expanded_Df2$obs)
      percentage <- sum(TCR_Expanded_Df2$percent)

      withProgress(message = 'Performing FindMarker ImmDom Analysis', value = 0, {


        for (i in 1:observations) {
          incProgress(1/observations, detail = paste("Clone", i,"of",observations))
          name.clone <- BD_sum$cluster_name[i]

          message(paste0("Downloading Dom stats files and dot plot...",name.clone))

          sc@meta.data$Gene_select <- ifelse(sc@meta.data$Vgene %in% name.clone,name.clone,"other")
          sc@meta.data
          Idents(object = sc) <- sc@meta.data$Gene_select

          min.pct.expression<- input$min_point_ #standard setting: 0.25
          min.logfc<-  input$LogFC_ #0.25 is standard
          # p.val.cutoff <-  input$pval_top #(1/10^3) is standard, use (1/10^0) to ignore

          cluster.names <- unique(Idents(sc))[order(unique(Idents(sc)))]
          # print(paste0("calculating markers for cluster ",name.clone,". Total: ",length(cluster.names)," clusters"))
          markers.fm.list <- FindMarkers(sc, ident.1 = name.clone, min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
          markers.fm.list2 <- subset(markers.fm.list,markers.fm.list$p_val_adj < input$pval.ex.filter)

          x=today()
          clonotype.name.stats <- paste("Prioritisation/ImmunoDom/",i,"_",gsub("[/]","",gsub("&","",name.clone)),"_stats_table","_",today(), ".csv", sep = "")
          write.csv(markers.fm.list2,clonotype.name.stats,row.names = T)

          message(paste0("Saved csv",name.clone))
          list.names <- rownames(markers.fm.list2)

          if (length(rownames(markers.fm.list2))>40) {
            list.names <- list.names[1:40]
          }

          else {
            list.names <- rownames(markers.fm.list2)
          }

          size_legend = input$Bar_legend_size-2



          plotdotplot <- DotPlot(sc, features = list.names) +
            RotatedAxis() +
            theme(
              axis.title.y = element_blank(),
              axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
              axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size, angle = 90),
              axis.title.x = element_blank(),
              legend.title = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.text = element_text(colour="black", size=size_legend,family=input$font_type),
              legend.position = input$legend_position,
            ) +
            scale_colour_gradient2(low = input$low.dotplot, mid = input$middle.dotplot, high = input$high.dotplot) +
            scale_x_discrete(labels = label_wrap(20)) +
            scale_y_discrete(labels = label_wrap(20))


          file.name.clone <- paste("Prioritisation/ImmunoDom/",i,"_",gsub("[/]","",gsub("&","",name.clone)),"_dotplot_plot","_",today(), ".png", sep = "")

          ### download the dot plot -------
          png(file.name.clone, width = input$width_png_all_expression_dotplot_top, height = input$height_png_all_expression_dotplot_top,res = input$resolution_PNG_all_expression_dotplot_top)
          plot(plotdotplot)
          dev.off()


          ##### download the OverRep ------

          df <- sc@meta.data
          # require()

          geneSet <- read.csv(system.file("OverRep","GeneSets.csv",package = "STEGO.R"),header = T)

          background.genes.name <- as.data.frame(rownames(sc@assays$RNA$scale.data))
          names(background.genes.name) <- "V1"
          background.genes <- length(rownames(sc@assays$RNA$scale.data))


          geneSet$background.genes <- background.genes

          DEx.genes <- as.data.frame(rownames(markers.fm.list2))
          names(DEx.genes) <- "V1"

          total.sig <- length(DEx.genes$V1)
          geneSet$total.sig <- length(DEx.genes$V1)

          geneSet$background.geneset <- NA
          geneSet$background.geneset.name <- NA
          geneSet$in.geneset <- NA
          geneSet$in.geneset.name <- NA

          if(input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
            geneSet$GeneSet <- gsub("-",".",geneSet$GeneSet)
          }

          if(input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
            require(stringr)
            geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
          }

          for (j in 1:dim(geneSet)[1]) {
            # listed GeneSet
            message(paste("GeneSet: ", j))
            Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet,";")[j])
            names(Gene.set.testing) <- "V1"
            Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
            names(Gene.set.testing2) <- "V1"
            background.overlap <- merge(Gene.set.testing2,background.genes.name,by= "V1")
            geneSet$background.geneset[j] <- length(background.overlap$V1)
            geneSet$background.geneset.name[j] <- as.character(paste(unlist(background.overlap[1]), collapse=';'))
            # in sig gene list
            overlap <- merge(background.overlap,DEx.genes,by= "V1")

            geneSet$in.geneset[j] <- length(overlap$V1)
            geneSet$in.geneset.name[j] <- as.character(paste(unlist(overlap[1]), collapse=';'))

          }

          geneSet2 <- subset(geneSet,geneSet$in.geneset>0)

          for (k in 1:dim(geneSet2)[1]) {
            tota.gene.set <- geneSet2$background.geneset[k] # genes that are identified in background
            tota.gene.set
            in.geneset <-  geneSet2$in.geneset[k]# DEx in geneset

            background.genes
            not.in.total <- background.genes - tota.gene.set
            not.in.geneset.sig <- total.sig - in.geneset
            d <- data.frame( gene.in.interest=c( in.geneset, not.in.geneset.sig),gene.not.interest=c( tota.gene.set, not.in.total))
            row.names(d) <- c("In_category", "not_in_category")

            if (in.geneset>0) {
              geneSet2$p.val[k] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
              geneSet2$lowerCI[k] <-  unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
              geneSet2$upperCI[k] <-unlist(fisher.test(d)$conf.int)[2]
              geneSet2$OR[k] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1],3)
            }

            else {
              geneSet2$p.value[k] <- "-"
              geneSet2$lowerCI[k] <-  "-"
              geneSet2$upperCI[k] <- "-"
              geneSet2$OR[k] <- "-"
            }
          }

          geneSet2 <- geneSet2[order(geneSet2$p.val,decreasing = F),]
          geneSet2 <- subset(geneSet2,geneSet2$in.geneset>=input$in.geneset.cutoff_top)
          geneSet2 <- subset(geneSet2,geneSet2$p.val<=input$p.val_cutoff_top)
          geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
          geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")

          file.name.clone <- paste("Prioritisation/ImmunoDom/",i,"_",gsub("[/]","",gsub("&","",name.clone)),"_dotplot_plot","_",today(), ".png", sep = "")
          top.name.overrep <- paste("Prioritisation/ImmunoDom/",i,"_",gsub("[/]","",gsub("&","",name.clone)),"_OverRep","_",today(), ".csv", sep = "")
          write.csv(geneSet2,top.name.overrep, row.names = F)



        }

      })
    })

    # one sample one indiv Polyclonal -----
    observeEvent(input$Poly_download_buttonOneOne,{
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "upload file")
      )

      TCR_Expanded_Df <- TCR_Expanded_fun(sc,(input$Samp_col),(input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      TCR_Expanded_Df2 <-  subset(TCR_Expanded_Df,TCR_Expanded_Df$percent>input$cut.off_percent_rep)
      observations <- sum(TCR_Expanded_Df2$obs)

      x = today()
      message("Downloading contig summary table...")
      top.name.clonotypes <- paste("Prioritisation/PolyClonal/Expansion_summary_table_",x,".csv",sep="")
      write.csv(Top_clonotype_df2(),top.name.clonotypes, row.names = F)

      message("Downloading the count UMAP...")
      top.name.clonotypes.top_png <- paste("Prioritisation/PolyClonal/Expansion_UMAP_count_",x,".png",sep="")
      png(top.name.clonotypes.top_png, width = input$width_png_TCR.UMAP,height = input$height_png_TCR.UMAP,res = input$resolution_PNG_TCR.UMAP)
      plot(UMAP.TCRclonalit2())
      dev.off()

      message("Downloading the expansion UMAP...")
      top.name.clonotypes.top_png <- paste("Prioritisation/PolyClonal/Exp_UMAP_cutoff_count.",input$cut.off_expanded,".and.freq",input$cutoff.expanded,"_",x,".png",sep="")
      # png(top.name.clonotypes.top_png, width = input$width_png_TCR.UMAP,height = input$height_png_TCR.UMAP,res = input$resolution_PNG_TCR.UMAP)
      png(top.name.clonotypes.top_png, width = input$width_png_UMAP_Expanded,
          height = input$height_png_UMAP_Expanded ,
          res = input$resolution_PNG_UMAP_Expanded )
      plot(UMAP_Expanded_plot())
      dev.off()

      message("Downloading stats table...")
      Exp_stats_cutoff_count.name <- paste("Prioritisation/PolyClonal/Exp_stats_cutoff_count.",input$cut.off_expanded,".and.freq",input$cutoff.expanded,"_",x,".csv",sep="")
      write.csv(Vals_expanded.stats(),Exp_stats_cutoff_count.name, row.names = T)

      message("Downloading the expansion stats dotplot...")
      top.name.clonotypes.top_png <- paste("Prioritisation/PolyClonal/Exp_dotplot.",input$cut.off_expanded,".and.freq",input$cutoff.expanded,"_",x,".png",sep="")
      png(top.name.clonotypes.top_png, width = input$width_png_all_expression_dotplot_ex,
          height = input$height_png_all_expression_dotplot_ex,
          res = input$resolution_PNG_all_expression_dotplot_ex)
      plot(relative_expression_plot_ex())
      dev.off()

      message("Downloading overrep table...")
      Exp_stats_cutoff_count.name <- paste("Prioritisation/PolyClonal/Exp_OverRep_cutoff_count.",input$cut.off_expanded,".and.freq",input$cutoff.expanded,"_",x,".csv",sep="")
      write.csv(Over_rep_Exp(),Exp_stats_cutoff_count.name, row.names = F)


    })


    # multiple samples (either individuals or samples) ------
    Upset_plot_overlap_Multi <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "upload file")
      )

      df <- sc@meta.data
      df <- as.data.frame(df)
      unique.df <- unique(df[,names(df) %in% c(input$Samp_col,input$V_gene_sc) ])
      names(unique.df) <- c("group","chain")
      unique.df <- subset(unique.df,unique.df$chain != "NA")
      unique.df <- subset(unique.df,unique.df$group != "NA")
      unique.df$cloneCount <- 1
      mat <- acast(unique.df, chain~group, value.var="cloneCount")
      mat[is.na(mat)] <- 0
      Count_data <- as.data.frame(rowSums(mat))
      names(Count_data) <- "V1"
      unique.df <- (df[,names(df) %in% c(input$Samp_col,input$V_gene_sc) ])
      names(unique.df) <- c("group","chain")
      unique.df <- subset(unique.df,unique.df$chain != "NA")
      unique.df <- subset(unique.df,unique.df$group != "NA")
      unique.df$cloneCount <- 1
      mat <- acast(unique.df, chain~group, value.var="cloneCount")
      mat[is.na(mat)] <- 0
      sum_data <- as.data.frame(rowSums(mat))
      names(sum_data) <- "V1"
      mat <- as.data.frame(mat)
      mat$TotalSamps <-Count_data$V1
      mat$CloneTotal <-sum_data$V1
      mat <- mat[order(mat$CloneTotal, decreasing = T),]
      mat <- mat[order(mat$TotalSamps, decreasing = T),]
      mat
    })

    Upset_plot_multi <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "upload file")
      )
      df <- sc@meta.data
      df <- as.data.frame(df)
      unique.df <- unique(df[,names(df) %in% c(input$Samp_col,input$V_gene_sc) ])
      names(unique.df) <- c("group","chain")
      unique.df <- unique.df[unique.df$group %in% input$ID_Column_factor,]

      unique.df <- subset(unique.df,unique.df$chain != "NA")
      unique.df <- subset(unique.df,unique.df$group != "NA")
      unique.df$cloneCount <- 1
      mat <- acast(unique.df, chain~group, value.var="cloneCount")
      mat[is.na(mat)] <- 0
      mat <- as.data.frame(mat)
      df.x <- make_comb_mat(mat)
      list.names <- as.character(input$ID_Column_factor)

      ht = draw(UpSet(df.x,
                      pt_size = unit(5, "mm"),
                      lwd = 1,
                      row_names_gp =  gpar(fontfamily = input$font_type, fontsize = 12),
                      column_names_gp = gpar(fontfamily = input$font_type,fontsize =12),
                      top_annotation = upset_top_annotation(df.x,
                                                            add_numbers = T,
                                                            numbers_gp = gpar(fontfamily = input$font_type,fontsize =12),
                                                            annotation_name_gp = gpar(fontfamily = input$font_type, fontsize = 12),
                                                            gp = gpar(fill = "black"),
                      ),
                      right_annotation = upset_right_annotation(df.x,
                                                                add_numbers = T,
                                                                numbers_gp = gpar(fontfamily = input$font_type,fontsize =12),
                                                                annotation_name_gp = gpar(fontfamily = input$font_type,fontsize=12),
                                                                gp = gpar(fill = "black"),
                      ),
                      set_order  =  list.names

      ), padding = unit(c(20, 20, 20, 20), "mm"))
      ht
    })
    clonal_plot_multi <- reactive({
      df4 <- TCR_Expanded()
      df4
      names(df4)[names(df4) %in% input$Samp_col] <- "ID_Column"
      df4 <- df4[df4$ID_Column %in% input$ID_Column_factor,]
      df4$ID_Column <- as.character(df4$ID_Column)
      df4$ID_Column <- factor(df4$ID_Column,levels = input$ID_Column_factor)

      # df4 <- TCR_Expanded()
      df4 <- df4[order(df4[,names(df4) %in% input$Graph_type_bar]),]
      df4
      col.df <- as.data.frame(unique(df4[,names(df4) %in% input$Graph_type_bar]))
      names(col.df) <- "V1"

      colorblind_vector <-as.data.frame(unlist(colors_clonal_plot()))

      if (dim(colorblind_vector)[1]==0) {
        num <- length(col.df$V1)

        if (input$colourtype == "default") {
          colorblind_vector <- c(gg_fill_hue(num))
        } else if (input$colourtype == "hcl.colors") {
          colorblind_vector <- c(hcl.colors(num, palette = "viridis"))
        } else if (input$colourtype == "topo.colors") {
          colorblind_vector <- c(topo.colors(num))
        } else if (input$colourtype == "heat.colors") {
          colorblind_vector <- c(heat.colors(num))
        } else if (input$colourtype == "terrain.colors") {
          colorblind_vector <- c(terrain.colors(num))
        } else if (input$colourtype == "rainbow") {
          colorblind_vector <- c(rainbow(num))
        } else if (input$colourtype == "random") {
          colorblind_vector <- distinctColorPalette(num)

        }  else {

        }

      }
      col.df$col <- colorblind_vector
      col.df


      # ggplot(df4, aes(y = frequency, x = Sample_Name, fill=Clonality, label = Clonality)) +
      #   # ggplot(df4,aes(x=Sample_Name,y=frequency,fill=Clonality))+
      #   geom_bar(stat="identity")+
      #   theme_bw() +
      #   scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 20),values=alpha(heat.colors(5), 1))

      ggplot(df4,aes(x=ID_Column,y=frequency,fill=get(input$Graph_type_bar),colour= get(input$Graph_type_bar),label=get(input$Graph_type_bar)))+
        geom_bar(stat="identity")+
        theme_bw() +
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 50),values=alpha(col.df$col, 1),na.value = input$NA_col_analysis) +
        scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 50),values = alpha(col.df$col, 1)) +
        theme(
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=90),
          axis.title.x = element_blank(),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )

    })

    Top_clonotypes_multiCounts <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "upload file")
      )

      BD_sum <- Upset_plot_overlap_Multi()
      BD_sum <- subset(BD_sum,BD_sum$CloneTotal>2)
      BD_sum <- subset(BD_sum,BD_sum$TotalSamps>1)
      BD_sum$priority <- 1/(BD_sum$CloneTotal * BD_sum$TotalSamps)
      BD_sum$same <- ifelse(BD_sum$CloneTotal==BD_sum$TotalSamps,"NEx","Ex")
      BD_sum <- subset(BD_sum,BD_sum$same=="Ex")
      BD_sum <- BD_sum[,!names(BD_sum) %in% "same"]
      BD_sum$cluster_name <- rownames(BD_sum)
      BD_sum

    })


    Top_clonotypes_multiCounts_barplot <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload")
      )
      x = today()
      req(input$Samp_col,input$V_gene_sc,input$cut.off_percent_repMulti)
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$V_gene_sc]


      BD_sum <- Top_clonotypes_multiCounts()
      BD_sum$obs <- 1
      BD_sum <-  subset(BD_sum,BD_sum$priority<input$cut.off_percent_repMulti)

      top.name.clonotypes.top_png <- paste("Prioritisation/Multi/PublicLike/","Selected_clones_",x,".csv",sep="")
      write.csv(BD_sum,top.name.clonotypes.top_png)
      observations <- sum(BD_sum$obs)

      withProgress(message = 'Performing Multi Overlap Analysis (barplots)', value = 0, {
        for (i in 1:observations) {

          incProgress(1/observations, detail = paste("Clone", i,"of",observations))
          message(BD_sum$cluster_name[i])
          name.clone <- BD_sum$cluster_name[i]
          top_BD_clonotype <- df3.meta[df3.meta$cluster_name %in% name.clone,]
          # print(top_BD_clonotype)

          dtop_clonotype_bar_code <- top_BD_clonotype

          # req(input$Graph_split_order)

          dtop_clonotype_bar_code$Selected_group <- dtop_clonotype_bar_code[,names(dtop_clonotype_bar_code) %in% input$Split_group_by_]
          num <- 1
          # num <- as.data.frame(num[complete.cases(num)==T,])
          as.data.frame(length(num))
          if (input$colourtype == "default") {
            colorblind_vector <- gg_fill_hue(num)
          } else if (input$colourtype == "hcl.colors") {
            colorblind_vector <- c(hcl.colors(num, palette = "viridis"))
          } else if (input$colourtype == "topo.colors") {
            colorblind_vector <- c(topo.colors(num))
          } else if (input$colourtype == "heat.colors") {
            colorblind_vector <- c(heat.colors(num))
          } else if (input$colourtype == "terrain.colors") {
            colorblind_vector <- c(terrain.colors(num))
          } else if (input$colourtype == "rainbow") {
            colorblind_vector <- c(rainbow((num)))
          } else if (input$colourtype == "random") {
            colorblind_vector <- distinctColorPalette(num)

          }  else {

          }


          colorblind_vector <- as.data.frame(colorblind_vector)
          names(colorblind_vector) <- "cols"

          dtop_clonotype_bar_code$Selected_chain2 <- dtop_clonotype_bar_code[,names(dtop_clonotype_bar_code) %in% input$V_gene_sc]
          dtop_clonotype_bar_code$Selected_chain3 <- gsub("_"," ",dtop_clonotype_bar_code$Selected_chain2)
          dtop_clonotype_bar_code$Selected_chain3 <- gsub("[.]"," ",dtop_clonotype_bar_code$Selected_chain3)

          dtop_clonotype_bar_code <- dtop_clonotype_bar_code[dtop_clonotype_bar_code$Selected_group %in% input$Graph_split_order,]
          dtop_clonotype_bar_code$Selected_group <- factor(dtop_clonotype_bar_code$Selected_group,levels = input$Graph_split_order)

          ggplot_plot <- ggplot(dtop_clonotype_bar_code, aes(x=Selected_group, fill=Selected_chain3,colour = Selected_chain3, label = Selected_chain3)) +
            geom_bar() +
            theme_bw()+
            scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector$cols, na.value=input$NA_col_analysis)+
            scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector$cols, na.value=input$NA_col_analysis)+
            # scale_alpha_manual(values = rep(1,length(unique(dtop_clonotype_bar_code$Selected_chain))), na.value=0.5)+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
            theme(
              axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
              axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
              axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=90),
              axis.title.x = element_blank(),
              legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.position = input$legend_position,
              legend.title = element_blank()
            ) +
            guides(color = "none", size = "none")


          top.name.clonotypes.top_png <- paste("Prioritisation/Multi/PublicLike/",i,"_top_clone_",gsub("[/]","",gsub("&","",name.clone)),"_",x,".png",sep="")

          num_width <- length(unique(dtop_clonotype_bar_code$Selected_group))

          png(top.name.clonotypes.top_png, width = (num_width*100+500),height = input$height_png_TCR.UMAP,res = input$resolution_PNG_TCR.UMAP)
          plot(ggplot_plot)
          dev.off()

        }
      })

    })

    top_clone_FindMaker_looped_Multi <- reactive({

      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload")
      )
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$V_gene_sc]
      sc@meta.data$Vgene <- sc@meta.data[,names(sc@meta.data) %in% input$V_gene_sc]

      BD_sum <- Top_clonotypes_multiCounts()
      BD_sum <-  subset(BD_sum,BD_sum$priority<input$cut.off_percent_repMulti)
      BD_sum$obs <- 1
      observations <- sum(BD_sum$obs)


      withProgress(message = 'Performing Multi Analysis (FindMarkers, Dotplot, OverRep)', value = 0, {


        for (i in 1:observations) {
          incProgress(1/observations, detail = paste("Clone", i,"of",observations))
          name.clone <- BD_sum$cluster_name[i]

          message(paste0("Downloading Dom stats files and dot plot...",name.clone))

          sc@meta.data$Gene_select <- ifelse(sc@meta.data$Vgene %in% name.clone,name.clone,"other")
          sc@meta.data
          Idents(object = sc) <- sc@meta.data$Gene_select

          min.pct.expression<- input$min_point_ #standard setting: 0.25
          min.logfc<-  input$LogFC_ #0.25 is standard
          # p.val.cutoff <-  input$pval_top #(1/10^3) is standard, use (1/10^0) to ignore

          cluster.names <- unique(Idents(sc))[order(unique(Idents(sc)))]
          # print(paste0("calculating markers for cluster ",name.clone,". Total: ",length(cluster.names)," clusters"))
          markers.fm.list <- FindMarkers(sc, ident.1 = name.clone, min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
          markers.fm.list2 <- subset(markers.fm.list,markers.fm.list$p_val_adj < input$pval.ex.filter)

          x=today()
          clonotype.name.stats <- paste("Prioritisation/Multi/PublicLike/",i,"_",gsub("[/]","",gsub("&","",name.clone)),"_stats_table","_",today(), ".csv", sep = "")
          write.csv(markers.fm.list2,clonotype.name.stats,row.names = T)

          message(paste0("Saved csv ",name.clone))
          list.names <- rownames(markers.fm.list2)

          if (length(rownames(markers.fm.list2))>40) {
            list.names <- list.names[1:40]
          }

          else {
            list.names <- rownames(markers.fm.list2)
          }

          size_legend = input$Bar_legend_size-2


          plotdotplot <- DotPlot(sc, features = list.names) +
            RotatedAxis() +
            theme(
              axis.title.y = element_blank(),
              axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
              axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size, angle = 90),
              axis.title.x = element_blank(),
              legend.title = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.text = element_text(colour="black", size=size_legend,family=input$font_type),
              legend.position = input$legend_position,
            ) +
            scale_colour_gradient2(low = input$low.dotplot, mid = input$middle.dotplot, high = input$high.dotplot) +
            scale_x_discrete(labels = label_wrap(20)) +
            scale_y_discrete(labels = label_wrap(20))


          file.name.clone <- paste("Prioritisation/Multi/PublicLike/",i,"_",gsub("[/]","",gsub("&","",name.clone)),"_dotplot_plot","_",today(), ".png", sep = "")

          ### download the dot plot -------
          png(file.name.clone, width = input$width_png_all_expression_dotplot_top, height = input$height_png_all_expression_dotplot_top,res = input$resolution_PNG_all_expression_dotplot_top)
          plot(plotdotplot)
          dev.off()


          ##### download the OverRep ------
          geneSet <- read.csv(system.file("OverRep","GeneSets.csv",package = "STEGO.R"),header = T)

          background.genes.name <- as.data.frame(rownames(sc@assays$RNA$scale.data))
          names(background.genes.name) <- "V1"
          background.genes <- length(rownames(sc@assays$RNA$scale.data))


          geneSet$background.genes <- background.genes

          DEx.genes <- as.data.frame(rownames(markers.fm.list2))
          names(DEx.genes) <- "V1"

          total.sig <- length(DEx.genes$V1)
          geneSet$total.sig <- length(DEx.genes$V1)

          geneSet$background.geneset <- NA
          geneSet$background.geneset.name <- NA
          geneSet$in.geneset <- NA
          geneSet$in.geneset.name <- NA

          if(input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
            geneSet$GeneSet <- gsub("-",".",geneSet$GeneSet)
          }

          if(input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
            require(stringr)
            geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
          }
          message(paste(i, "Performing Over rep analysis"))
          for (j in 1:dim(geneSet)[1]) {
            # listed GeneSet

            Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet,";")[j])
            names(Gene.set.testing) <- "V1"
            Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
            names(Gene.set.testing2) <- "V1"
            background.overlap <- merge(Gene.set.testing2,background.genes.name,by= "V1")
            geneSet$background.geneset[j] <- length(background.overlap$V1)
            geneSet$background.geneset.name[j] <- as.character(paste(unlist(background.overlap[1]), collapse=';'))
            # in sig gene list
            overlap <- merge(background.overlap,DEx.genes,by= "V1")

            geneSet$in.geneset[j] <- length(overlap$V1)
            geneSet$in.geneset.name[j] <- as.character(paste(unlist(overlap[1]), collapse=';'))

          }

          geneSet2 <- subset(geneSet,geneSet$in.geneset>0)

          for (k in 1:dim(geneSet2)[1]) {
            tota.gene.set <- geneSet2$background.geneset[k] # genes that are identified in background
            tota.gene.set
            in.geneset <-  geneSet2$in.geneset[k]# DEx in geneset

            background.genes
            not.in.total <- background.genes - tota.gene.set
            not.in.geneset.sig <- total.sig - in.geneset
            d <- data.frame( gene.in.interest=c( in.geneset, not.in.geneset.sig),gene.not.interest=c( tota.gene.set, not.in.total))
            row.names(d) <- c("In_category", "not_in_category")

            if (in.geneset>0) {
              geneSet2$p.val[k] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
              geneSet2$lowerCI[k] <-  unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
              geneSet2$upperCI[k] <-unlist(fisher.test(d)$conf.int)[2]
              geneSet2$OR[k] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1],3)
            }

            else {
              geneSet2$p.value[k] <- "-"
              geneSet2$lowerCI[k] <-  "-"
              geneSet2$upperCI[k] <- "-"
              geneSet2$OR[k] <- "-"
            }
          }

          geneSet2 <- geneSet2[order(geneSet2$p.val,decreasing = F),]
          geneSet2 <- subset(geneSet2,geneSet2$in.geneset>=input$in.geneset.cutoff_top)
          geneSet2 <- subset(geneSet2,geneSet2$p.val<=input$p.val_cutoff_top)
          geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
          geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")

          top.name.overrep <- paste("Prioritisation/Multi/PublicLike/",i,"_",gsub("[/]","",gsub("&","",name.clone)),"_OverRep","_",today(), ".csv", sep = "")
          write.csv(geneSet2,top.name.overrep, row.names = F)

        }

      })
    })

    # download the top x per individual
    Top_clonotypes_Private <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "upload file")
      )
      BD_sum <- Upset_plot_overlap_Multi()
      BD_sum <- subset(BD_sum,BD_sum$CloneTotal>2)
      BD_sum <- subset(BD_sum,BD_sum$TotalSamps==1)
      BD_sum$cluster_name <- rownames(BD_sum)
      BD_sum

    })


    top_clone_FindMaker_looped_Private <- reactive({

      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload")
      )
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$V_gene_sc]
      sc@meta.data$Vgene <- sc@meta.data[,names(sc@meta.data) %in% input$V_gene_sc]

      BD_sum <- Top_clonotypes_Private()
      BD_sum$obs <- 1
      observations <- sum(BD_sum$obs)

      withProgress(message = 'Performing Multi Analysis (FindMarkers, Dotplot, OverRep)', value = 0, {

        for (i in 1:observations) {
          incProgress(1/observations, detail = paste("Clone", i,"of",observations))
          name.clone <- BD_sum$cluster_name[i]

          message(paste0("Downloading Dom stats files and dot plot...",name.clone))

          sc@meta.data$Gene_select <- ifelse(sc@meta.data$Vgene %in% name.clone,name.clone,"other")
          sc@meta.data
          Idents(object = sc) <- sc@meta.data$Gene_select

          min.pct.expression<- input$min_point_ #standard setting: 0.25
          min.logfc<-  input$LogFC_ #0.25 is standard
          # p.val.cutoff <-  input$pval_top #(1/10^3) is standard, use (1/10^0) to ignore

          cluster.names <- unique(Idents(sc))[order(unique(Idents(sc)))]
          # print(paste0("calculating markers for cluster ",name.clone,". Total: ",length(cluster.names)," clusters"))
          markers.fm.list <- FindMarkers(sc, ident.1 = name.clone, min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
          markers.fm.list2 <- subset(markers.fm.list,markers.fm.list$p_val_adj < input$pval.ex.filter)

          x=today()
          clonotype.name.stats <- paste("Prioritisation/Multi/Unique/",i,"_",gsub("[/]","",gsub("&","",name.clone)),"_stats_table","_",today(), ".csv", sep = "")
          write.csv(markers.fm.list2,clonotype.name.stats,row.names = T)

          message(paste0("Saved csv ",name.clone))
          list.names <- rownames(markers.fm.list2)

          if (length(rownames(markers.fm.list2))>40) {
            list.names <- list.names[1:40]
          }

          else {
            list.names <- rownames(markers.fm.list2)
          }

          size_legend = input$Bar_legend_size-2


          plotdotplot <- DotPlot(sc, features = list.names) +
            RotatedAxis() +
            theme(
              axis.title.y = element_blank(),
              axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
              axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size, angle = 90),
              axis.title.x = element_blank(),
              legend.title = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.text = element_text(colour="black", size=size_legend,family=input$font_type),
              legend.position = input$legend_position,
            ) +
            scale_colour_gradient2(low = input$low.dotplot, mid = input$middle.dotplot, high = input$high.dotplot) +
            scale_x_discrete(labels = label_wrap(20)) +
            scale_y_discrete(labels = label_wrap(20))


          file.name.clone <- paste("Prioritisation/Multi/Unique/",i,"_",gsub("[/]","",gsub("&","",name.clone)),"_dotplot_plot","_",today(), ".png", sep = "")

          ### download the dot plot -------
          png(file.name.clone, width = input$width_png_all_expression_dotplot_top, height = input$height_png_all_expression_dotplot_top,res = input$resolution_PNG_all_expression_dotplot_top)
          plot(plotdotplot)
          dev.off()


          ##### download the OverRep ------
          geneSet <- read.csv(system.file("OverRep","GeneSets.csv",package = "STEGO.R"),header = T)

          background.genes.name <- as.data.frame(rownames(sc@assays$RNA$scale.data))
          names(background.genes.name) <- "V1"
          background.genes <- length(rownames(sc@assays$RNA$scale.data))


          geneSet$background.genes <- background.genes

          DEx.genes <- as.data.frame(rownames(markers.fm.list2))
          names(DEx.genes) <- "V1"

          total.sig <- length(DEx.genes$V1)
          geneSet$total.sig <- length(DEx.genes$V1)

          geneSet$background.geneset <- NA
          geneSet$background.geneset.name <- NA
          geneSet$in.geneset <- NA
          geneSet$in.geneset.name <- NA

          if(input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
            geneSet$GeneSet <- gsub("-",".",geneSet$GeneSet)
          }

          if(input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
            require(stringr)
            geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
          }
          message(paste(i, "Performing Over rep analysis"))
          for (j in 1:dim(geneSet)[1]) {
            # listed GeneSet

            Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet,";")[j])
            names(Gene.set.testing) <- "V1"
            Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
            names(Gene.set.testing2) <- "V1"
            background.overlap <- merge(Gene.set.testing2,background.genes.name,by= "V1")
            geneSet$background.geneset[j] <- length(background.overlap$V1)
            geneSet$background.geneset.name[j] <- as.character(paste(unlist(background.overlap[1]), collapse=';'))
            # in sig gene list
            overlap <- merge(background.overlap,DEx.genes,by= "V1")

            geneSet$in.geneset[j] <- length(overlap$V1)
            geneSet$in.geneset.name[j] <- as.character(paste(unlist(overlap[1]), collapse=';'))

          }

          geneSet2 <- subset(geneSet,geneSet$in.geneset>0)

          for (k in 1:dim(geneSet2)[1]) {
            tota.gene.set <- geneSet2$background.geneset[k] # genes that are identified in background
            tota.gene.set
            in.geneset <-  geneSet2$in.geneset[k]# DEx in geneset

            background.genes
            not.in.total <- background.genes - tota.gene.set
            not.in.geneset.sig <- total.sig - in.geneset
            d <- data.frame( gene.in.interest=c( in.geneset, not.in.geneset.sig),gene.not.interest=c( tota.gene.set, not.in.total))
            row.names(d) <- c("In_category", "not_in_category")

            if (in.geneset>0) {
              geneSet2$p.val[k] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
              geneSet2$lowerCI[k] <-  unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
              geneSet2$upperCI[k] <-unlist(fisher.test(d)$conf.int)[2]
              geneSet2$OR[k] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1],3)
            }

            else {
              geneSet2$p.value[k] <- "-"
              geneSet2$lowerCI[k] <-  "-"
              geneSet2$upperCI[k] <- "-"
              geneSet2$OR[k] <- "-"
            }
          }

          geneSet2 <- geneSet2[order(geneSet2$p.val,decreasing = F),]
          geneSet2 <- subset(geneSet2,geneSet2$in.geneset>=input$in.geneset.cutoff_top)
          geneSet2 <- subset(geneSet2,geneSet2$p.val<=input$p.val_cutoff_top)
          geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
          geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")

          top.name.overrep <- paste("Prioritisation/Multi/Unique/",i,"_",gsub("[/]","",gsub("&","",name.clone)),"_OverRep","_",today(), ".csv", sep = "")
          write.csv(geneSet2,top.name.overrep, row.names = F)

        }

      })
    })

    observeEvent(input$Multi_download_button,{
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "upload file")
      )
      x = today()
      message("Downloading the Summary table...")
      top.name.clonotypes <- paste("Prioritisation/Multi/Overlap_Table_",x,".csv",sep="")
      write.csv(Upset_plot_overlap_Multi(),top.name.clonotypes, row.names = T)

      message("Downloading the Upset Plot...")
      df <- sc@meta.data
      df <- as.data.frame(df)
      unique.df <- unique(df[,names(df) %in% c(input$Samp_col,input$V_gene_sc) ])
      names(unique.df) <- c("group","chain")
      unique.df <- unique.df[unique.df$group %in% input$ID_Column_factor,]
      if (length(unique(unique.df$group))<31) {

        top.name.clonotypes.count_png <- paste("Prioritisation/Multi/Overlap_Upset_plot",x,".png",sep="")
        png(top.name.clonotypes.count_png, width = input$width_png_TCR.UMAP_top,height = input$height_png_TCR.UMAP_top,res = input$resolution_PNG_TCR.UMAP_top)
        plot(Upset_plot_multi())
        dev.off()
      } else {
        message("More than 31 groups, you will need to create the plot manually")
      }

      message("Downloading the stacked barplot Plot...")
      top.name.clonotypes.count_png <- paste("Prioritisation/Multi/Stacked_bar_plot",x,".png",sep="")

      df4 <- TCR_Expanded()
      df4
      names(df4)[names(df4) %in% input$Samp_col] <- "ID_Column"
      num_indiv <- length(unique(df4$ID_Column))

      png(top.name.clonotypes.count_png, width = (num_indiv*100+600),height = input$height_png_TCR.UMAP_top,res = input$resolution_PNG_TCR.UMAP_top)
      plot(clonal_plot_multi())
      dev.off()

      message("Downloading Multi bar plots...")
      Top_clonotypes_multiCounts_barplot()

      message("Downloading Multi Marker and OverRep analysis ")
      top_clone_FindMaker_looped_Multi()

      message("Downloading Private Marker and OverRep analysis summary table ...")
      top.name.clonotypes <- paste("Prioritisation/Multi/Unique_Table_",x,".csv",sep="")
      write.csv(Top_clonotypes_Private(),top.name.clonotypes, row.names = T)

      message("Downloading Private Marker and OverRep analysis ")
      top_clone_FindMaker_looped_Private()

    })

    output$Test_table_1 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload")
      )
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$V_gene_sc]
      sc@meta.data$Vgene <- sc@meta.data[,names(sc@meta.data) %in% input$V_gene_sc]

      BD_sum <- Top_clonotypes_Private()
      BD_sum$obs <- 1
      observations <- sum(BD_sum$obs)
      name.clone <- BD_sum$cluster_name[1]
      sc@meta.data$Gene_select <- ifelse(sc@meta.data$Vgene %in% name.clone,name.clone,"other")
      sc@meta.data


      # num_width <- length(unique(dtop_clonotype_bar_code$Selected_group))
    })

    # add in upset plot and table to show clonal expansion -> then make expansion overall with consideration to the sample. Also add in count, not just percentage.
    output$Number_of_clonotypes_to_ <- renderPrint({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload")
      )
      x = today()
      req(input$Samp_col,input$V_gene_sc,input$cut.off_percent_repMulti)
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$V_gene_sc]

      BD_sum <- Top_clonotypes_multiCounts()
      BD_sum$obs <- 1
      BD_sum <-  subset(BD_sum,BD_sum$priority<input$cut.off_percent_repMulti)
      observations <- sum(BD_sum$obs)



      BD_sum <- Top_clonotypes_Private()
      BD_sum$obs <- 1
      observations2 <- sum(BD_sum$obs)
      df1 <- paste("The analysis will be limited to the top",observations,"Public/MultiSample Clones")
      df2 <- paste("The analysis will be limited to the top",observations2,"Private Clones")
      df3 <- rbind(df1,df2)
      rownames(df3) <- c(1,2)
      df3
      # paste()
    })

    # Clustering priority  -----
    observeEvent(input$ClusterDownload_download_buttonOneOne,{
      x = today()
      if (length(AG_cluster())>0) {
        message("Downloading AG cluster table...")
        Exp_stats_cutoff_count.name <- paste("Prioritisation/Clustering/Cluster_summary_table_AG_",x,".csv",sep="")
        AG_cluster <- AG_cluster()

        write.csv(AG_cluster,Exp_stats_cutoff_count.name, row.names = F)
        message("Downloading AG cluster analysis...")
        req(input$priority_cutoff)
        ggPlotUMAPClusterAG()
      } else {

      }

    })
    observeEvent(input$ClusterDownload_download_buttonOneOne,{
      x = today()

      if (length(BD_cluster())>0) {
        message("Downloading BD cluster table...")
        Exp_stats_cutoff_count.name <- paste("Prioritisation/Clustering/Cluster_summary_table_BD_",x,".csv",sep="")
        BD_cluster <- BD_cluster()

        write.csv(BD_cluster,Exp_stats_cutoff_count.name, row.names = F)
        message("Downloading BD cluster analysis...")
        ggPlotUMAPClusterBD()
      } else {
      }

    })


    output$number_clusters_to_analyse_AG <- renderPrint({
      df1 <- AG_cluster()
      validate(
        need(nrow(df1)>0,
             "Upload ClusTCR file")
      )
      df1 <- subset(df1,df1$priority < input$priority_cutoff)

      paste("The analysis will be limited to the top", max(df1$Updated_order),"AG cluster(s)")
    })
    output$number_clusters_to_analyse_BD <- renderPrint({
      df1 <- BD_cluster()
      validate(
        need(nrow(df1)>0,
             "Upload ClusTCR file")
      )
      df1 <- subset(df1,df1$priority < input$priority_cutoffBD)

      paste("The analysis will be limited to the top", max(df1$Updated_order),"BD cluster(s)")
    })

    ggPlotUMAPClusterAG <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload Files")
      )
      md <- sc@meta.data
      clusterAG2 <- AG_cluster()
      validate(
        need(nrow(clusterAG2)>0,
             "Upload clusTCR AG table")
      )
      req(clusterAG2,input$Clusters_to_dis_PIE,input$Colour_By_this,input$priority_cutoff)
      clusterAG <- subset(clusterAG2,clusterAG2$priority < input$priority_cutoff)

      len.order <- length(unique(clusterAG$Updated_order))
      # clusterAG <- subset(cluster,cluster$priority<input$priority_cutoff)
      withProgress(message = 'Performing AG cluster analysis', value = 0, {

        for (i in 1:len.order) {
          incProgress(1/len.order, detail = paste("AG cluster", i,"of",len.order))
          ## ggplot UMAP -----
          cluster <- clusterAG
          names(cluster)[names(cluster) %in% input$Samp_col_cluster] <- "ID_Column"

          cluster <- cluster[cluster$Updated_order %in% i,]
          cluster$colour <- cluster[,names(cluster) %in% input$Colour_By_this]
          cluster$colour <- gsub("_"," ",cluster$colour)
          cluster$colour <- factor(cluster$colour, levels = unique(cluster$colour))
          cluster$colour <- gsub("NA",NA,cluster$colour)

          len.colour <- length(unique(cluster$colour))
          col.df <- as.data.frame(unique(cluster$colour))

          num <- (length(unlist(col.df)))

          if (input$colourtype == "default") {
            colorblind_vector <- c(gg_fill_hue(num))
          } else if (input$colourtype == "hcl.colors") {
            colorblind_vector <- c(hcl.colors(num, palette = "viridis"))
          } else if (input$colourtype == "topo.colors") {
            colorblind_vector <- c(topo.colors(num))
          } else if (input$colourtype == "heat.colors") {
            colorblind_vector <- c(heat.colors(num))
          } else if (input$colourtype == "terrain.colors") {
            colorblind_vector <- c(terrain.colors(num))
          } else if (input$colourtype == "rainbow") {
            colorblind_vector <- c(rainbow(num))
          } else if (input$colourtype == "random") {
            colorblind_vector <- distinctColorPalette(num)

          }  else {

          }


          col.df$col <- colorblind_vector

          figure <- ggplot(data=cluster,aes(x=UMAP_1,UMAP_2,colour=colour))+
            geom_point(size = input$size.dot.umap)+
            scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 20),values = col.df$col,na.value=input$NA_col_analysis) +
            theme_bw()+
            theme(
              legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.title = element_blank(),
              legend.position = input$legend_position,
              # strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
              axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
              axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
              axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=0),
              axis.title.x = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
            )

          message(paste(i," Downloading the count UMAP"))
          x = today()
          top.name.clonotypes.top_png <- paste("Prioritisation/Clustering/",i,"_AG_cluster_UMAP_",x,".png",sep="")
          png(top.name.clonotypes.top_png, width = input$width_png_TCR.UMAP,height = input$height_png_TCR.UMAP,res = input$resolution_PNG_TCR.UMAP)
          plot(figure)
          dev.off()
          ## Motif plot -----

          Network_df <- cluster[order(cluster$Updated_order),]
          Network_df <- Network_df %>% distinct(CDR3_Vgene, .keep_all = TRUE) # make Unique
          motifplot <- Motif_from_cluster_file(Network_df,Clust_selected = i,selected_cluster_column = "Updated_order")

          message(paste(i," Downloading motif plot"))
          top.name.clonotypes.top_png <- paste("Prioritisation/Clustering/",i,"_AG_motif_",x,".png",sep="")
          png(top.name.clonotypes.top_png, width = input$width_png_Motif_ClusTCR2_cluster,
              height = input$height_png_Motif_ClusTCR2_cluster,
              res = input$resolution_PNG_Motif_ClusTCR2_cluster)
          plot(motifplot)
          dev.off()

          # Stats table ------
          cluster <- clusterAG

          names(cluster)[names(cluster) %in% input$Samp_col_cluster] <- "ID_Column"
          cluster <- cluster[order(cluster$Updated_order),]

          rownames(cluster) <- cluster$Cell_Index

          checking <- cluster[,names(cluster) %in% c("Updated_order","Cell_Index")]
          md.checking <- merge(md,checking,by="Cell_Index",all.x=T)
          md.checking <- md.checking[order(md.checking$order),]
          rownames(md.checking) <- md.checking$Cell_Index

          md.checking$Clust_selected <- ifelse(md.checking$Updated_order == i,i,"NS")
          md.checking$Clust_selected[is.na(md.checking$Clust_selected)] <- "NS"
          md.checking <- md.checking[order(md.checking$order),]

          sc@meta.data <- md.checking
          Idents(object = sc) <- sc@meta.data$Clust_selected

          name.check.clust <- i
          min.pct.expression<- input$min_point_ #standard setting: 0.25
          min.logfc<-  input$LogFC_ #0.25 is standard

          markers.fm.list <- FindMarkers(sc, ident.1 = name.check.clust, min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
          markers.fm.list2 <- subset(markers.fm.list,markers.fm.list$p_val_adj < input$pval.ex.filter)
          message(paste(i," Downloading stats table"))
          Exp_stats_cutoff_count.name <- paste("Prioritisation/Clustering/",i,"_AG_cluster_statsTab_",x,".csv",sep="")

          write.csv(markers.fm.list2,Exp_stats_cutoff_count.name, row.names = T)

          # stats dotplot ----

          if (length(rownames(markers.fm.list2))<40) {
            list.names <- rownames(markers.fm.list2)
          } else {
            list.names <- rownames(markers.fm.list2)
            list.names <- list.names[1:40]
          }

          size_legend = input$Bar_legend_size-2

          dotplotClust <- DotPlot(sc, features = list.names) +
            RotatedAxis() +
            theme(
              axis.title.y = element_blank(),
              axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
              axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size, angle = 90),
              axis.title.x = element_blank(),
              legend.title = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.text = element_text(colour="black", size=size_legend,family=input$font_type),
              legend.position = input$legend_position,
            ) +
            scale_colour_gradient2(low = input$low.dotplot.clust, mid = input$middle.dotplot.clust, high = input$high.dotplot.clust)+
            scale_x_discrete(labels = label_wrap(20)) +
            scale_y_discrete(labels = label_wrap(20))

          top.name.clonotypes.top_png <- paste("Prioritisation/Clustering/",i,"_AG_dot.plot_",x,".png",sep="")
          png(top.name.clonotypes.top_png, width = input$width_png_all_expression_dotplot_clust,
              height = input$height_png_all_expression_dotplot_clust,
              res = input$resolution_PNG_all_expression_dotplot_clust)
          plot(dotplotClust)
          dev.off()

          # column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_all_expression_dotplot_clust','Download PDF')),
          # column(2,numericInput("width_png_all_expression_dotplot_clust","Width of PNG", value = 2400)),
          # column(2,numericInput("height_png_all_expression_dotplot_clust","Height of PNG", value = 700)),
          # column(2,numericInput("resolution_PNG_all_expression_dotplot_clust","Resolution of PNG", value = 144)),
          # column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_all_expression_dotplot_clust','Download PNG'))

          # stats OverRep analysis ----
          geneSet <- read.csv(system.file("OverRep","GeneSets.csv",package = "STEGO.R"),header = T)

          background.genes.name <- as.data.frame(rownames(sc@assays$RNA$scale.data))
          names(background.genes.name) <- "V1"
          background.genes <- length(rownames(sc@assays$RNA$scale.data))

          #
          geneSet$background.genes <- background.genes

          DEx.genes <- as.data.frame(rownames(markers.fm.list2))
          names(DEx.genes) <- "V1"
          total.sig <- length(DEx.genes$V1)
          geneSet$total.sig <- length(DEx.genes$V1)
          # geneSet
          geneSet$background.geneset <- NA
          geneSet$background.geneset.name <- NA
          geneSet$in.geneset <- NA
          geneSet$in.geneset.name <- NA

          if(input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
            geneSet$GeneSet <- gsub("-",".",geneSet$GeneSet)
          }

          if(input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
            require(stringr)
            geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
          }
          message(paste("Starting OverRep analysis of cluster ", i))
          for (j in 1:dim(geneSet)[1]) {
            # listed GeneSet

            Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet,";")[j])
            names(Gene.set.testing) <- "V1"
            Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
            names(Gene.set.testing2) <- "V1"
            background.overlap <- merge(Gene.set.testing2,background.genes.name,by= "V1")
            geneSet$background.geneset[j] <- length(background.overlap$V1)
            geneSet$background.geneset.name[j] <- as.character(paste(unlist(background.overlap[1]), collapse=';'))
            # in sig gene list
            overlap <- merge(background.overlap,DEx.genes,by= "V1")

            geneSet$in.geneset[j] <- length(overlap$V1)
            geneSet$in.geneset.name[j] <- as.character(paste(unlist(overlap[1]), collapse=';'))

          }
          geneSet2 <- subset(geneSet,geneSet$in.geneset>0)
          for (j in 1:dim(geneSet2)[1]) {
            tota.gene.set <- geneSet2$background.geneset[j] # genes that are identified in background
            in.geneset <-  geneSet2$in.geneset[j]# DEx in geneset
            not.in.total <- background.genes - tota.gene.set
            not.in.geneset.sig <- total.sig - in.geneset
            d <- data.frame( gene.in.interest=c( in.geneset, not.in.geneset.sig),gene.not.interest=c(tota.gene.set, not.in.total))
            row.names(d) <- c("In_category", "not_in_category")

            if (in.geneset>0) {
              geneSet2$p.val[j] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
              geneSet2$lowerCI[j] <-  unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
              geneSet2$upperCI[j] <-unlist(fisher.test(d)$conf.int)[2]
              geneSet2$OR[j] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1],3)
            } else {
              geneSet2$p.value[j] <- "-"
              geneSet2$lowerCI[j] <-  "-"
              geneSet2$upperCI[j] <- "-"
              geneSet2$OR[j] <- "-"
            }
          }
          geneSet2 <- geneSet2[order(geneSet2$p.val,decreasing = F),]
          geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
          geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
          message("Downloading the Summary table...")
          top.name.clonotypes <- paste("Prioritisation/Clustering/",i,"_AG_OverRep_",x,".csv",sep="")
          write.csv(geneSet2,top.name.clonotypes, row.names = F)

        }
      })

    })

    ggPlotUMAPClusterBD <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload Files")
      )
      md <- sc@meta.data
      clusterBD2 <- BD_cluster()
      validate(
        need(nrow(clusterBD2)>0,
             "Upload clusTCR BD table")
      )
      req(clusterBD2,input$Clusters_to_dis_PIE,input$Colour_By_this,input$priority_cutoffBD)
      clusterBD <- subset(clusterBD2,clusterBD2$priority < input$priority_cutoffBD)

      len.order <- length(unique(clusterBD$Updated_order))
      # clusterAG <- subset(cluster,cluster$priority<input$priority_cutoff)
      withProgress(message = 'Performing BD cluster analysis', value = 0, {

        for (i in 1:len.order) {
          incProgress(1/len.order, detail = paste("AG cluster", i,"of",len.order))
          cluster <- clusterBD
          names(cluster)[names(cluster) %in% input$Samp_col_cluster] <- "ID_Column"
          cluster <- cluster[cluster$Updated_order %in% i,]
          ## ggplot UMAP -----
          cluster$colour <- cluster[,names(cluster) %in% input$Colour_By_this]
          cluster$colour <- gsub("_"," ",cluster$colour)
          cluster$colour <- factor(cluster$colour, levels = unique(cluster$colour))
          cluster$colour <- gsub("NA",NA,cluster$colour)

          len.colour <- length(unique(cluster$colour))
          col.df <- as.data.frame(unique(cluster$colour))

          num <- (length(unlist(col.df)))

          if (input$colourtype == "default") {
            colorblind_vector <- c(gg_fill_hue(num))
          } else if (input$colourtype == "hcl.colors") {
            colorblind_vector <- c(hcl.colors(num, palette = "viridis"))
          } else if (input$colourtype == "topo.colors") {
            colorblind_vector <- c(topo.colors(num))
          } else if (input$colourtype == "heat.colors") {
            colorblind_vector <- c(heat.colors(num))
          } else if (input$colourtype == "terrain.colors") {
            colorblind_vector <- c(terrain.colors(num))
          } else if (input$colourtype == "rainbow") {
            colorblind_vector <- c(rainbow(num))
          } else if (input$colourtype == "random") {
            colorblind_vector <- distinctColorPalette(num)

          }  else {

          }


          col.df$col <- colorblind_vector

          figure <- ggplot(data=cluster,aes(x=UMAP_1,UMAP_2,colour=colour))+
            geom_point(size = input$size.dot.umap)+
            scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 20),values = col.df$col,na.value=input$NA_col_analysis) +
            theme_bw()+
            theme(
              legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.title = element_blank(),
              legend.position = input$legend_position,
              # strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
              axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
              axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
              axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=0),
              axis.title.x = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
            )

          message(paste(i," Downloading the count UMAP"))
          x = today()
          top.name.clonotypes.top_png <- paste("Prioritisation/Clustering/",i,"_BD_cluster_UMAP_",x,".png",sep="")
          png(top.name.clonotypes.top_png, width = input$width_png_TCR.UMAP,height = input$height_png_TCR.UMAP,res = input$resolution_PNG_TCR.UMAP)
          plot(figure)
          dev.off()
          ## Motif plot -----
          Network_df <- cluster[order(cluster$Updated_order),]
          Network_df <- Network_df %>% distinct(CDR3_Vgene, .keep_all = TRUE) # make Unique
          motifplot <- Motif_from_cluster_file(Network_df,Clust_selected = i,selected_cluster_column = "Updated_order")

          message(paste(i," Downloading motif plot"))
          top.name.clonotypes.top_png <- paste("Prioritisation/Clustering/",i,"_BD_motif_",x,".png",sep="")
          png(top.name.clonotypes.top_png, width = input$width_png_Motif_ClusTCR2_cluster,
              height = input$height_png_Motif_ClusTCR2_cluster,
              res = input$resolution_PNG_Motif_ClusTCR2_cluster)
          plot(motifplot)
          dev.off()

          # Stats table ------
          cluster <- clusterBD

          names(cluster)[names(cluster) %in% input$Samp_col_cluster] <- "ID_Column"
          cluster <- cluster[order(cluster$Updated_order),]

          rownames(cluster) <- cluster$Cell_Index

          checking <- cluster[,names(cluster) %in% c("Updated_order","Cell_Index")]
          md.checking <- merge(md,checking,by="Cell_Index",all.x=T)
          md.checking <- md.checking[order(md.checking$order),]
          rownames(md.checking) <- md.checking$Cell_Index

          md.checking$Clust_selected <- ifelse(md.checking$Updated_order == i,i,"NS")
          md.checking$Clust_selected[is.na(md.checking$Clust_selected)] <- "NS"
          md.checking <- md.checking[order(md.checking$order),]

          sc@meta.data <- md.checking
          Idents(object = sc) <- sc@meta.data$Clust_selected

          name.check.clust <- i
          min.pct.expression<- input$min_point_ #standard setting: 0.25
          min.logfc<-  input$LogFC_ #0.25 is standard

          markers.fm.list <- FindMarkers(sc, ident.1 = name.check.clust, min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
          markers.fm.list2 <- subset(markers.fm.list,markers.fm.list$p_val_adj < input$pval.ex.filter)
          message(paste(i," Downloading stats table"))
          Exp_stats_cutoff_count.name <- paste("Prioritisation/Clustering/",i,"_BD_cluster_statsTab_",x,".csv",sep="")
          write.csv(markers.fm.list2,Exp_stats_cutoff_count.name, row.names = T)

          # stats dotplot ----

          if (length(rownames(markers.fm.list2))<40) {
            list.names <- rownames(markers.fm.list2)
          } else {
            list.names <- rownames(markers.fm.list2)
            list.names <- list.names[1:40]
          }

          size_legend = input$Bar_legend_size-2

          dotplotClust <- DotPlot(sc, features = list.names) +
            RotatedAxis() +
            theme(
              axis.title.y = element_blank(),
              axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
              axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size, angle = 90),
              axis.title.x = element_blank(),
              legend.title = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.text = element_text(colour="black", size=size_legend,family=input$font_type),
              legend.position = input$legend_position,
            ) +
            scale_colour_gradient2(low = input$low.dotplot.clust, mid = input$middle.dotplot.clust, high = input$high.dotplot.clust)+
            scale_x_discrete(labels = label_wrap(20)) +
            scale_y_discrete(labels = label_wrap(20))

          top.name.clonotypes.top_png <- paste("Prioritisation/Clustering/",i,"_BD_dot.plot_",x,".png",sep="")
          png(top.name.clonotypes.top_png, width = input$width_png_all_expression_dotplot_clust,
              height = input$height_png_all_expression_dotplot_clust,
              res = input$resolution_PNG_all_expression_dotplot_clust)
          plot(dotplotClust)
          dev.off()

          # stats OverRep analysis ----

          if(dim(markers.fm.list2)[1]>0) {

            geneSet <- read.csv(system.file("OverRep","GeneSets.csv",package = "STEGO.R"),header = T)

            background.genes.name <- as.data.frame(rownames(sc@assays$RNA$scale.data))
            names(background.genes.name) <- "V1"
            background.genes <- length(rownames(sc@assays$RNA$scale.data))

            #
            geneSet$background.genes <- background.genes

            DEx.genes <- as.data.frame(rownames(markers.fm.list2))
            names(DEx.genes) <- "V1"
            total.sig <- length(DEx.genes$V1)
            geneSet$total.sig <- length(DEx.genes$V1)
            # geneSet
            geneSet$background.geneset <- NA
            geneSet$background.geneset.name <- NA
            geneSet$in.geneset <- NA
            geneSet$in.geneset.name <- NA

            if(input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
              geneSet$GeneSet <- gsub("-",".",geneSet$GeneSet)
            }

            if(input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
              require(stringr)
              geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
            }
            message(paste("Starting OverRep analysis of cluster ", i))
            for (j in 1:dim(geneSet)[1]) {
              # listed GeneSet

              Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet,";")[j])
              names(Gene.set.testing) <- "V1"
              Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
              names(Gene.set.testing2) <- "V1"
              background.overlap <- merge(Gene.set.testing2,background.genes.name,by= "V1")
              geneSet$background.geneset[j] <- length(background.overlap$V1)
              geneSet$background.geneset.name[j] <- as.character(paste(unlist(background.overlap[1]), collapse=';'))
              # in sig gene list
              overlap <- merge(background.overlap,DEx.genes,by= "V1")

              geneSet$in.geneset[j] <- length(overlap$V1)
              geneSet$in.geneset.name[j] <- as.character(paste(unlist(overlap[1]), collapse=';'))

            }

            geneSet
            geneSet2 <- subset(geneSet,geneSet$in.geneset>0)
            for (j in 1:dim(geneSet2)[1]) {
              tota.gene.set <- geneSet2$background.geneset[j] # genes that are identified in background
              in.geneset <-  geneSet2$in.geneset[j]# DEx in geneset
              not.in.total <- background.genes - tota.gene.set
              not.in.geneset.sig <- total.sig - in.geneset
              d <- data.frame( gene.in.interest=c( in.geneset, not.in.geneset.sig),gene.not.interest=c(tota.gene.set, not.in.total))
              row.names(d) <- c("In_category", "not_in_category")

              if (in.geneset>0) {
                geneSet2$p.val[j] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
                geneSet2$lowerCI[j] <-  unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
                geneSet2$upperCI[j] <-unlist(fisher.test(d)$conf.int)[2]
                geneSet2$OR[j] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1],3)
              } else {
                geneSet2$p.value[j] <- "-"
                geneSet2$lowerCI[j] <-  "-"
                geneSet2$upperCI[j] <- "-"
                geneSet2$OR[j] <- "-"
              }
            }
            geneSet2 <- geneSet2[order(geneSet2$p.val,decreasing = F),]
            geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
            geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
            message("Downloading the Summary table...")
            top.name.clonotypes <- paste("Prioritisation/Clustering/",i,"_BD_OverRep_",x,".csv",sep="")
            write.csv(geneSet2,top.name.clonotypes, row.names = F)

          }
        }
      })
    })

    output$Cluster_dowload_button_prior <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "upload file")
      )

      TCR_Expanded_Df <- TCR_Expanded_fun(sc,(input$Samp_col),(input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      TCR_Expanded_Df2 <-  subset(TCR_Expanded_Df,TCR_Expanded_Df$percent>input$cut.off_percent_rep)
      observations <- sum(TCR_Expanded_Df2$obs)
      length.samp.ID <- length(unique(TCR_Expanded_Df$ID_Column))

      fluidRow(
        # column(12,selectInput("AnalysisType","Preset parameters for", = c("ImmunoDom"))),
        column(12,actionButton("ClusterDownload_download_buttonOneOne","Download cluster analysis",style = "color: white; background-color:#4682b4"))
      )
    })

    ### clustering table -----
    output$PriorClustTB_Tab <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload Files")
      )
      df1 <- AG_cluster()
      df1 <- subset(df1,df1$priority < input$priority_cutoff)
      df1
    })


    ### test table -----
    output$colors.top_dt <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      sc <- UMAP_metadata_with_labs()
      validate(
        need(nrow(sc)>0,
             "Upload File")
      )
      md <- sc@meta.data
      x = today()

      if (length(input.data_sc_clusTCR_BD())>0) {
        clust <- input.data_sc_clusTCR_BD()

        if(input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") {
          names(md)[names(md) %in% "v_gene_BD"] <- "Selected_V_BD"
          names(md)[names(md) %in% "junction_aa_BD"] <- "AminoAcid_BD"
        } else {
          names(md)[names(md) %in% "v_gene_BD"] <- "Selected_V_BD"
          names(md)[names(md) %in% "cdr3_BD"] <- "AminoAcid_BD"
        }
      }
      md

      # names(cluster)[names(cluster) %in% input$Samp_col_cluster] <- "ID_Column"
      # cluster <- cluster[order(cluster$Updated_order),]
      #
      # rownames(cluster) <- cluster$Cell_Index
      #
      # checking <- cluster[,names(cluster) %in% c("Updated_order","Cell_Index")]
      # checking
      # md.checking <- merge(md,checking,by="Cell_Index",all.x=T)
      # md.checking$Clust_selected <- ifelse(md.checking$Updated_order == 1,1,"NS")
      # md.checking$Clust_selected[is.na(md.checking$Clust_selected)] <- "NS"
      # md.checking <- md.checking[order(md.checking$order),]
      #
      # sc@meta.data <- md.checking
      # Idents(object = sc) <- sc@meta.data$Clust_selected
      #
      #
      # md.checking <- md.checking[order(md.checking$order),]
      # rownames(md.checking) <- md.checking$Cell_Index
      #
      # #####
      # sc@meta.data <- md.checking
      # Idents(object = sc) <- sc@meta.data$Clust_selected
      #
      # name.check.clust <- 1
      # min.pct.expression<- input$min_point_ #standard setting: 0.25
      # min.logfc<-  input$LogFC_ #0.25 is standard
      #
      # markers.fm.list <- FindMarkers(sc, ident.1 = name.check.clust, min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
      # markers.fm.list2 <- subset(markers.fm.list,markers.fm.list$p_val_adj < input$pval.ex.filter)
      # markers.fm.list2
      #
      # geneSet <- read.csv(system.file("OverRep","GeneSets.csv",package = "STEGO.R"),header = T)
      # if(input$SeuratVersion == "Version 4") {
      #   background.genes.name <- as.data.frame(rownames(sc@assays$RNA$scale.data))
      #   names(background.genes.name) <- "V1"
      #   background.genes <- length(rownames(sc@assays$RNA$scale.data))
      # } else {
      #   background.genes.name <- as.data.frame(rownames(sc@assays$RNA$scale.data))
      #   names(background.genes.name) <- "V1"
      #   background.genes <- length(rownames(sc@assays$RNA$scale.data))
      # }
      #
      # #
      # geneSet$background.genes <- background.genes
      #
      # DEx.genes <- as.data.frame(rownames(markers.fm.list2))
      # names(DEx.genes) <- "V1"
      # total.sig <- length(DEx.genes$V1)
      # geneSet$total.sig <- length(DEx.genes$V1)
      # # geneSet
      # geneSet$background.geneset <- NA
      # geneSet$background.geneset.name <- NA
      # geneSet$in.geneset <- NA
      # geneSet$in.geneset.name <- NA
      #
      # if(input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
      #   geneSet$GeneSet <- gsub("-",".",geneSet$GeneSet)
      # }
      #
      # if(input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
      #   require(stringr)
      #   geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
      # }
      #
      # for (j in 1:dim(geneSet)[1]) {
      #   # listed GeneSet
      #   message(paste("GeneSet: ", j))
      #   Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet,";")[j])
      #   names(Gene.set.testing) <- "V1"
      #   Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
      #   names(Gene.set.testing2) <- "V1"
      #   background.overlap <- merge(Gene.set.testing2,background.genes.name,by= "V1")
      #   geneSet$background.geneset[j] <- length(background.overlap$V1)
      #   geneSet$background.geneset.name[j] <- as.character(paste(unlist(background.overlap[1]), collapse=';'))
      #   # in sig gene list
      #   overlap <- merge(background.overlap,DEx.genes,by= "V1")
      #
      #   geneSet$in.geneset[j] <- length(overlap$V1)
      #   geneSet$in.geneset.name[j] <- as.character(paste(unlist(overlap[1]), collapse=';'))
      #
      # }
      #
      # geneSet
      # geneSet2 <- subset(geneSet,geneSet$in.geneset>0)
      # for (j in 1:dim(geneSet2)[1]) {
      #   tota.gene.set <- geneSet2$background.geneset[j] # genes that are identified in background
      #   in.geneset <-  geneSet2$in.geneset[j]# DEx in geneset
      #   not.in.total <- background.genes - tota.gene.set
      #   not.in.geneset.sig <- total.sig - in.geneset
      #   d <- data.frame( gene.in.interest=c( in.geneset, not.in.geneset.sig),gene.not.interest=c(tota.gene.set, not.in.total))
      #   row.names(d) <- c("In_category", "not_in_category")
      #
      #   if (in.geneset>0) {
      #     geneSet2$p.val[j] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
      #     geneSet2$lowerCI[j] <-  unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
      #     geneSet2$upperCI[j] <-unlist(fisher.test(d)$conf.int)[2]
      #     geneSet2$OR[j] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1],3)
      #   } else {
      #     geneSet2$p.value[j] <- "-"
      #     geneSet2$lowerCI[j] <-  "-"
      #     geneSet2$upperCI[j] <- "-"
      #     geneSet2$OR[j] <- "-"
      #   }
      # }
      # geneSet2 <- geneSet2[order(geneSet2$p.val,decreasing = F),]
      # # geneSet2 <- subset(geneSet2,geneSet2$in.geneset>=input$in.geneset.cutoff_Clust)
      # # geneSet2 <- subset(geneSet2,geneSet2$p.val<=input$p.val_cutoff_Clust)
      # geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
      # geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
      # geneSet2

    })
    ### end -----
  }
  shinyApp(ui, server)

}
