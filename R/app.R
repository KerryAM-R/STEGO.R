#' STEGO application
#' @name STEGO
#' @export

runSTEGO <- function(){
  ### packages ------
  suppressWarnings(require("DescTools"))
  suppressMessages(require("bslib"))
  suppressMessages(require("circlize")) # colorRamp2
  suppressMessages(suppressWarnings(require("ClusTCR2")))
  suppressMessages(require("ComplexHeatmap"))
  suppressMessages(require("corrplot"))
  suppressMessages(require("cowplot"))
  suppressMessages(require("doParallel"))
  suppressMessages(require("dplyr"))
  suppressMessages(require("DT"))
  suppressMessages(require("extrafont"))
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
  suppressMessages(require("grid"))
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
  # suppressMessages(require("SeuratDisk"))
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
  # needs to be called before the colour picker due to a dependency issue
  suppressWarnings(suppressMessages(require("shinyjs")))

  # needs to be loaded after shinyjs to prevent the colourpicker issue...
  suppressMessages(require("colourpicker")) # select visual colour

  suppressMessages(require("chisq.posthoc.test"))
  message("All packages loaded")

  # font ------
  fonts <- fonttable()
  font <- as.data.frame(unique(fonts$FamilyName))
  names(font) <- "Fonts"

  #####

  if (dir.exists("custom_db/")) {

    suppressWarnings(source(system.file("scGATE","custom_df_scGATE.R",package = "STEGO.R")))

    suppressWarnings(source(system.file("scGATE","custom_annotation_models.R",package = "STEGO.R")))
  } else {
    message("To use the custom scGate models, please open STEGO.R in the Directory_for_project")
    message("Open the STEGO.Rproj, as this sets the working directory for your project to the current folder. ")

  }

  if (!dir.exists("Prioritization/")) {
    message("Creating folders for the Prioritization, as the program wasn't loaded from the STEGO.Rproj file")
    dir.create("Prioritization")
    dir.create("Prioritization/Multi")
    dir.create("Prioritization/Multi/PublicLike")
    dir.create("Prioritization/Multi/Unique")
    dir.create("Prioritization/Clustering")
    dir.create("Prioritization/Clustering/AG")
    dir.create("Prioritization/Clustering/BD")
    dir.create("Prioritization/EpitopePred")
    } else {
  }

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

  #### check if OLGA is installed ----
  require(reticulate)
  olgafunction_BD <- function(y) {
    olga <- system2('olga-compute_pgen', args=c("--humanTRB ",
                                                y),
                    wait=TRUE, stdout=TRUE)
    olga

  }

  olgafunction_AG <- function(y) {
    olga <- system2('olga-compute_pgen', args=c("--humanTRA ",
                                                y),
                    wait=TRUE, stdout=TRUE)
    olga
  }


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

  ###################

  # UI page -----
  ui <- fluidPage(

    # progress bar colouring, position --------
    tags$style(HTML("
  .dataTables_wrapper .dataTable td {
    white-space: nowrap;
  }
")),

# check box size and colouring ------

tags$style("input[type=checkbox] {
                    transform: scale(1.5);
           }

           "),
tags$head(
  tags$style(HTML('
      input[type="checkbox"]:before {
          border: 2px solid;
          color: #E9C2FF;
          background-color: white;
          content: "";
          height:  15px;
          left: 0;
          position: absolute;
          top: 0;
          width:  15px;
           vertical-align: top;
      }

      input[type="checkbox"]:checked:after {
          border: 2px solid;
          color: #6F00B0;
          background-color: #E9C2FF;
          content: "✓";
          font-size: smaller;
          vertical-align: middle;
          text-align: center;
          height:  15px;
          left: 0;
          position: absolute;
          top: 0;
          width:  15px;
          vertical-align: top;
      }
      label.checkbox {

    line-height: 50px;
    margin: 2px 0;
    display: block;
    height: 50px;} '
  )
  )
),


# change spinner position -----
tags$head(
  tags$style(HTML(
    ".centered-spinner {
        position: fixed;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
      }"
  ))
),

# change highlighted text in dropdown menue ------
tags$head(tags$style(
  HTML("
                  .selectize-dropdown-content .active {
                  background: #E9C2FF !important;
                  color: #6F00B0 !important;
                  }
                  .select-dropdown-content .active {
                  background: #E9C2FF !important;
                  color: #6F00B0 !important;
                  }
       .selectize-dropdown-content {
                  background: #6F00B0;
                  color: #E9C2FF;
                  }
                  .select-dropdown-content {
                  background: #6F00B0;
                  color: #E9C2FF;
                  }
       .selectize-dropdown-content .selected {
              background: #d8ffc2 !important; /* Red background for selected item */
              color: #41b000  !important; /* White text color for selected item */
            }


")
)),

# if I want to update the grey colour of the side bar panel
tags$head(tags$style(
  HTML('
         .well {

        }

        ')
)),

# change the button colours and hover ------
tags$head(
  tags$style(HTML("
                  .btn {
                    border: 2px solid #6F00B0;
                    color: #6F00B0;
                    background-color: #E9C2FF;
                    border-color: #6F00B0
                    }

                    "))),
tags$head(tags$style(HTML("
                  .btn:hover {
                      border: 2px solid #6F00B0;
                      color: #E9C2FF;
                      background-color: #6F00B0;
                  }
"))),

tags$head(tags$style(HTML("
                  .btn:active {
                        border: 2px solid #41b000 !important;
                        background: #d8ffc2 !important;
                        color: #41b000 !important;

                  }
"))),

tags$head(
  tags$style(HTML('
      .nav > li.disabled-1 > a {

        background: #d8ffc2 !important;
        color: #41b000 !important;
        pointer-events: none;
      }
    '))
),

tags$head(
  tags$style(HTML('
      .nav > li.disabled-2 > a {

        background: #d8ffc2 !important;
        color: #41b000 !important;
        pointer-events: none;
      }
    '))
),

tags$head(
  tags$style(HTML('
      .nav > li.disabled-3 > a {

        background: #d8ffc2 !important;
        color: #41b000 !important;
        pointer-events: none;
      }
    '))
),


tags$head(
  tags$style(HTML("
    .my-selectize .selectize-input {
      width: 900px !important;  /* Set the width as needed */
    }
  "))
),

# progress bar for to be purple and in the centre ------
tags$head(tags$style(HTML(".progress-bar {background-color: #6F00B0}"))),
tags$head(
  tags$style(
    HTML(
      ".shiny-notification {
              height: 200px;
              width: 800px;
              position:fixed;
              top: calc(50% - 50px);
              left: calc(50% - 400px);
              font-size: 150%;
              text-align: center;
              background-color: #E9C2FF;
             color: black
            }"
    )
  )
),
tags$style(HTML("
    .tabbable > .nav > li > a {background-color: #E9C2FF; color:#6F00B0}
    .tabbable > .nav > li.active > a {background-color: #6F00B0; color:#E9C2FF}
  ")),

tags$head(
  tags$style(type='text/css',"
               .navbar-nav {font-size: 16px}
               .navbar-nav > li > a {background-color: #E9C2FF; color:#6F00B0}

               ")),

tags$head(tags$style(HTML('
  .navbar-nav > li > .dropdown-menu {font-size: 16px}
  .navbar-nav > li > .dropdown-menu {background-color: #E9C2FF}
  .navbar-nav > li.active > .dropdown-menu {background-color: #E9C2FF; color:black}
  .dropdown-menu > li.hover > a {background-color:#6F00B0;}
  .dropdown-menu > li.active > a {background-color:#6F00B0;}
  .dropdown-menu > li > a:hover {
    color: white !important;
    background-color: #6F00B0 !important;
                          }
  '))),

tags$head(tags$style(HTML('
  ".navbar-default .navbar-nav>li>a {background-color: white}"
  '))),

# input Text boxes headers -----

tags$style(HTML("
    /* Styling for text input with class name-BD */
    .name-BD input[type='text'] {
    background-color: #E9C2FF; /* Change to your desired color */
    color: #6F00B0; /* Change 'red' to your desired color */
    border: 2px solid #6F00B0
    }
  ")),

tags$style(HTML("
    .name-header {
      color: #6F00B0;
      margin-bottom: -2px;
      font-size: 16px
    }
  ")),

tags$style(HTML("
    .name-header2 {
      color: #6F00B0;
      font-size: 16px
    }
  ")),

tags$style(HTML("
    .name-header3 {
      color: #6F00B0;
    font-size: 16px
    }
  ")),
# add in disable function ----



#####
navbarPage(
  title = "STEGO.R",


  #### step 1 ------


  navbarMenu(
    "STEP 1.",
    # 10x_Genomics ----
    tabPanel(
      "10x_Genomics",
      sidebarLayout(
        sidebarPanel(
          id = "tPanel4", style = "overflow-y:scroll; max-height: 800px; position:relative;", width = 3,
          conditionalPanel(
            condition = "input.panel_10x != 'auto_10x'",
            h4("File name"),
            fluidRow(
              column(6, div(class = "name-BD", textInput("group10x", "Add Group name", ""))),
              column(6, div(class = "name-BD", textInput("Indiv10x", "Add Individual name", "")))
            ),
            selectInput("Source_type_10x", "Input types", choices = c("Raw", ".h5")),
            conditionalPanel(
              condition = "input.Source_type_10x=='Raw'",
              fileInput("file_calls_10x", "Barcode file (.tsv.gz or .tsv)"),
              fileInput("file_features_10x", "Features file (.tsv.gz or .tsv)"),
              fileInput("file_matrix_10x", "Matrix file (.mtx.gz or .mtx)"),
            ),
            selectInput("csv_contig_file", "format of the contig file", choices = c("csv/csv.gz", "tsv")),
            fileInput("file_TCR_10x", "filtered contig annotations"),
            selectInput("BCR_TCR_10x", "Type of data", choices = c("TCR only", "BCR only"))),
          conditionalPanel(
            condition = "input.panel_10x == 'auto_10x'",


          )

        ),
        # 10x main panel -----
        mainPanel(
          width = 9,
          tabsetPanel(
            id = "panel_10x",
            tabPanel("Uploaded data",
                     value = 1,
                     conditionalPanel(
                       condition = "input.Source_type_10x=='Raw'",
                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                       div(DT::DTOutput("Barcode_10x_file_upload", height = "200px")),
                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                       div(DT::DTOutput("Features_10x_file_upload", height = "200px")),
                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                       div(DT::DTOutput("Matrix_10x_file_upload", height = "200px")),
                     ),
                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     div(DT::DTOutput("Contig_10x_file_upload", height = "200px")),

            ),
            tabPanel(
              "TCRex",
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              # fluidRow(
              downloadButton("downloaddf_TCRex_10x", "Download table"),
              div(DT::DTOutput("tb_TCRex_10x_df"),height = "200px"),
              # )

            ),
            tabPanel("Seurat QC",
                     value = 2,
                     tags$head(tags$style("#tb_10x_matrix2  {white-space: nowrap;  }")),
                     fluidRow(
                       column(3, downloadButton("downloadtb_10x_matrix2", "Download matrix"),style = "padding-top: 25px;"),
                       column(3, downloadButton("downloadtb_10x_metadata2", "Download metadata"),style = "padding-top: 25px;")
                     ),
                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     div(DT::DTOutput("tb_10x_matrix2")),
                     tags$head(tags$style("#sum_tb_10x1  {white-space: nowrap;  }")),

                     div(DT::DTOutput("sum_tb_10x1")),
                     div(DT::DTOutput("tb_10x_meta1")),

            ),
            tabPanel("ClusTCR",
                     value = 3,
                     tags$head(tags$style("#tb_10x_contigues1  {white-space: nowrap;  }")),
                     fluidRow(
                       column(3,selectInput("chain_clusTCR2_10x", "Select to download", choices = c("AG", "BD", "IgH", "IgLK")), ),
                       column(3, downloadButton("downloadtb_10x_contigues1", "Download clusTCR"),style = "padding-top: 25px;")
                     ),

                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     div(DT::DTOutput("tb_10x_contigues1")),

            ),
            tabPanel("TCR_Explore",
                     value = 4,
                     tags$head(tags$style("#dt_TCR_Explore_10x  {white-space: nowrap;  }")),
                     fluidRow(
                       column(3,downloadButton("downloaddt_TCR_Explore_10x", "Download TCR_Explore"),style = "padding-top: 25px;"),
                     ),
                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     div(DT::DTOutput("dt_TCR_Explore_10x")),
            ),
            tabPanel("Auto",value = "auto_10x",
                     h4("Set up the directories to have folders named sample_condition with the barcode, features, matrix and contig file"),
                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     fluidRow(
                       column(3,div(actionButton("checkFiles", "Check Files"), style = "padding-top: 25px;")),
                       column(3,div(actionButton("automateProcess", "Automate Process"), style = "padding-top: 25px;")),
                     ),

                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     verbatimTextOutput("outputText")
            )
          ),
        )
      )
    ),
    # 10x_Genomics end -----
    # BD Rhapsody  ------

    tabPanel(
      "BD rhapsody data",
      sidebarLayout(
        sidebarPanel(
          id = "tPanel4", style = "overflow-y:scroll; max-height: 1000px; position:relative;", width = 3,
          # UPLOAD the three files...
          # selectInput("dataset_BD", "Choose a dataset:", choices = c("test_data_BD", "own_data_BD")),
          div(class = "name-BD",
              h4("Add File Name", class = "name-header"),
              textInput("name_BD",h5("e.g., Samp1_Tumor", class = "name-header"), value = "")),
          fluidRow(
            column(6, numericInput("no_lines_skip_Tags", "If needed, skip first 7 lines of Sample Tag", value = 7, min = 0, max = 20, step = 7), ),
            column(6, fileInput("file_calls_BD", "Sample Tag Calls (.csv)",
                                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
            ), )
          ),
          h5("Upload the matrix files"),
          selectInput("Format_bd", "Format type", choices = c("cellXgene", "Barcode_features_matrix"), selected = "Barcode_features_matrix"),
          conditionalPanel(
            condition = "input.Format_bd=='Barcode_features_matrix'",
            fileInput("file_barcode_bd", "Barcode file (.tsv.gz or .tsv)"),
            fileInput("file_features_bd", "Features file (.tsv.gz or .tsv)"),
            fileInput("file_matrix_bd", "Matrix file"),
          ),
          conditionalPanel(
            condition = "input.Format_bd=='cellXgene'",
            fluidRow(
              column(6, numericInput("no_lines_skip_counts", "If needed, skip first 7 lines of Count Matrix", value = 7, min = 0, max = 10, step = 7), ),
              column(6, fileInput("file_counts_BD", "Counts (.csv)",
                                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
              ), ),
            ),
          ),
          h5("Upload the TCR Contig file"),
          selectInput("filtered_list", "Contig Format", choices = c("Paired", "Dominant", "Unfiltered"), selected = "Dominant"),
          conditionalPanel(
            condition = "input.filtered_list=='Paired'",
            fluidRow(
              column(6, numericInput("no_lines_skip_TCR", "If needed, skip first 7 lines of VDJ Contig file", value = 7, min = 0, max = 10, step = 7), ),
              column(6, fileInput("file_TCR_BD", "TCR file (.csv)",
                                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
              )),
            ),
          ),
          conditionalPanel(
            condition = "input.filtered_list=='Dominant' || input.filtered_list=='Unfiltered'",
            fileInput("file_TCR_bd2", "Contig AIRR file (tsv or tsv.gz)",
                      accept = c(".tsv", "tsv",".tsv.gz","tsv.gz")
            )
          ),
          ### filter out non-function TCR and un-paired TCR
          conditionalPanel(
            condition = "input.filtered_list=='Dominant' || input.filtered_list=='Unfiltered'",
            selectInput("locus_column", h5("Chain (e.g. locus)"), ""),
          ),
          conditionalPanel(
            condition = "input.filtered_list=='Paired'",
            fluidRow(
              column(6, selectInput("V_gene_AG_BDrhap", h6("Alpha/Gamma V gene"), ""), ),
              column(6, selectInput("Junction_AG_BDrhap", h6("Alpha/Gamma junction"), ""), ),
              column(6, selectInput("V_gene_BD_BDrhap", h6("Beta/Delta V gene"), "")),
              column(6, selectInput("Junction_BD_BDrhap", h6("Beta/Delta junction"), ""))
            ),
          ),

          fluidRow(
            column(6, checkboxInput("filtering_TCR", p("paired chains?", class = "name-header3"), value = FALSE, width = NULL)),
            column(6, checkboxInput("TCR_present",  p("TCR present?", class = "name-header3"), value = TRUE, width = NULL)),
            column(6, checkboxInput("BCR_present",  p("BCR present?", class = "name-header3"), value = FALSE, width = NULL)),
          ),
        ),
        # main panel ------
        mainPanel(
          width = 9,
          tabsetPanel(
            tabPanel(
              "Imported data",
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              div(DT::DTOutput("test.files")),
              conditionalPanel(
                condition = "input.Format_bd=='cellXgene'",
                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                div(DT::DTOutput("test.files3"))
              ),
              conditionalPanel(
                condition = "input.Format_bd=='Barcode_features_matrix'",
                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                div(DT::DTOutput("test.files.bd1")),
                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                div(DT::DTOutput("test.files.bd2")),
                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                div(DT::DTOutput("test.files.bd3")),
              ),
              conditionalPanel(
                condition = "input.filtered_list=='Paired'",
                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                div(DT::DTOutput("test.files2"))
              ),
              conditionalPanel(
                condition = "input.filtered_list=='Dominant' || input.filtered_list=='Unfiltered'",
                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                div(DT::DTOutput("test.files.bd4")),
              ),
            ),
            # tabPanel("Checking Merge",
            #
            #          div(DT::DTOutput("Check_table")),
            #
            #          ),
            tabPanel(
              "clusTCR2",
              tags$head(tags$style("#tb_clusTCR  {white-space: nowrap;  }")),
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              div(DT::DTOutput("tb_clusTCR")),
              selectInput("chain_clusTCR2_bd", "Select to download", choices = c("AG", "BD", "IgH", "IgLK")),
              downloadButton("downloaddf_clusTCR", "Download table")
            ),
            tabPanel(
              "TCRex",
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              div(DT::DTOutput("tb_TCRex_BDrap_df")),
              downloadButton("downloaddf_TCRex_BDrap", "Download table")
            ),
            tabPanel(
              "For Seurat",
              tags$head(tags$style("#tb_count_matrix  {white-space: nowrap;  }")),
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              div(DT::DTOutput("tb_count_matrix")),
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              div(DT::DTOutput("tb_metadata_sc")),
              fluidRow(
                column(3, downloadButton("downloadtb_count_matrix", "Download count table")),
                column(3),
                column(3, downloadButton("downloadtb_metadata_sc", "Download meta.data table")),
              ),
            ),
            tabPanel(
              "TCR_Explore",
              tags$head(tags$style("#tb_TCR_Explore  {white-space: nowrap;  }")),
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              div(DT::DTOutput("tb_TCR_Explore")),
              downloadButton("downloadtb_TCR_Explore", "Download table")
            ),
            tabPanel(
              "Multi-TCR",
              div(DT::DTOutput("tb_multiTCR")),
              downloadButton("downloadtb_multiTCR", "Download Multi-TCR table")
            ),
            tabPanel(
              "Create Sample Tags file",
              tags$head(tags$style("#tb_sample_tags_created  {white-space: nowrap;  }")),
              div(class = "name-BD",textInput("sample_tags_name", "Name of sample", value = "BD EA splenocyte")),
              div(DT::DTOutput("tb_sample_tags_created")),
              downloadButton("downloadtb_sample_tags", "Download Tags")
            ),
          )
        )
      ),
    ),
    # BD rhapsody end ----

    # array format -----
    tabPanel(
      "Array",
      sidebarLayout(
        sidebarPanel(
          id = "tPanel4", style = "overflow-y:scroll; max-height: 800px; position:relative;", width = 3,
          fileInput("file_calls_Array", "Matrix file (.txt.gz)",
                    accept = c("text/csv", "text/comma-separated-values,text/plain", "txt.gz")
          ),
          fileInput("file_contig_Array", "TCR file (.txt.gz)",
                    accept = c("text/csv", "text/comma-separated-values,text/plain", "txt.gz")
          ),
          checkboxInput("pairing_TCR_Array", "Paired only?"),
          div(class = "name-BD",textInput("sample_name_Array", "Add sample name", "Treatment_group")),
          div(class = "name-BD",textInput("name.array", "Add file name", "FileName"))
        ),
        mainPanel(
          width = 9,
          tabsetPanel(
            tabPanel(
              "Check files",
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              div(DT::DTOutput("test.files_array_Matrix")),
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              div(DT::DTOutput("test.files_array_contig")),
            ),
            tabPanel(
              "Filtering",
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              div(DT::DTOutput("test.files_array_contig_Filtered")),
            ),
            tabPanel(
              "clusTCR",
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              div(DT::DTOutput("test.files_ClusTCR2_array")),
              downloadButton("download_ClusTCR2_labs_array"),
            ),
            tabPanel("TCRex"),
            tabPanel(
              "Seurat",
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              div(DT::DTOutput("tb_array_matrix2")),
              div(DT::DTOutput("tb_Array_meta1")),
              fluidRow(
                column(3, downloadButton("downloadtb_array_matrix2", "Download matrix")),
                column(3, downloadButton("downloadtb_array_metadata2", "Download metadata"))
              ),
            ),
            tabPanel("TCR_Explore")
          )
        )
      )
    ),




    # Convert from .h5Seurat to .rds --------
    tabPanel(
      "Convert format",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          div(class = "name-BD",textInput("project_name5", h4("Name of Project", class = "name-header2"), value = "")),
          # textInput("", "Name of Project", value = ""),
          fileInput("file1_h5Seurat.file",
                    "Choose .h5Seurat files from directory",
                    multiple = TRUE,
                    accept = c(".h5Seurat", "h5Seurat")
          ),
          div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),

          downloadButton("downloaddf_SeruatObj_Convert_to_RDS", "Download converted .rds Seurat Obj")
        ),
        mainPanel(
          width = 9,
          tabsetPanel(
            id = "Converting_formatting",
            tabPanel("Converting",
                     value = "converting_PA",
                     h4("Convert .h5Seurat (V4 Seurat) to .rds"),
                     verbatimTextOutput("Convert_to_RDS_out")
            )
          )
        )
      )
    ),
    # End Nav menue data -----
  ), # NavbarMenu
  navbarMenu(
    "STEP 2.",
    # TCRex merge files ----
    tabPanel(
      "2a. TCRex merge",
      sidebarLayout(
        sidebarPanel(
          id = "tPanel4", style = "overflow-y:scroll; max-height: 800px; position:relative;", width = 3,
          # selectInput("dataset2", "Choose a dataset:", choices = c("test_data_clusTCR2","own_data_clusTCR2")),
          fileInput("file2_TCRexMerge", "Select files to merge",
                    multiple = TRUE,
          ),
          downloadButton("downloaddf_TCRexFiltered", "Download table")
        ),
        mainPanel(
          width = 9,
          tabsetPanel(
            id = "TCRex_tabs",
            tabPanel("Merge Multiple Files",
                     value = "merge",
                     div(DT::DTOutput("DEx_TCRexFiltered")),
                     # downloadButton('downloaddf_multiple_ClusTCR2','Download table')
            ),
          )
        )
      )
    ),

    # TCR clustering with ClusTCR2 -----
    tabPanel(
      "2b. ClusTCR2",
      sidebarLayout(
        sidebarPanel(
          id = "tPanel4", style = "overflow-y:scroll; max-height: 800px; position:relative;", width = 3,
          # selectInput("dataset2", "Choose a dataset:", choices = c("test_data_clusTCR2","own_data_clusTCR2")),

          selectInput("Clust_lab_tab_output", "Add prefix to file", choices = c("AG", "BD", "IgH", "IgKL")),
          conditionalPanel(
            condition = "input.clusTCR2_tabs=='merge'",
            fileInput("file1_ClusTCR2_multiple", "Select files to merge",
                      multiple = TRUE,
                      accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
            )
          ),
          conditionalPanel(
            condition = "input.clusTCR2_tabs != 'merge'",

            conditionalPanel(
              condition = "input.clusTCR2_tabs2 == 'processing1'",
              fileInput("file2_ClusTCR2", "Select file for Clustering",
                        accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
              ),

              fluidRow(
                column(6, selectInput("clusTCR2_names", label = h5("CDR3"), "")),
                column(6, selectInput("clusTCR2_Vgene", label = h5("V gene"), "")),
              ),
              fluidRow(
                column(6, checkboxInput("allele_ClusTCR2", "Remove allele *00?", value = T)),
                # column(6, numericInput("cores_ClusTCR2","Number of cores to parallel",value=1))
              ),
            ),
            conditionalPanel(
              condition = "input.clusTCR2_tabs2 != 'processing1'",
              h5("Click to run/update clustering"),
              div(actionButton("run_ClusTCR2", "Run Clustering"), style = "padding-top: 25px;")
            ),
          ),
          conditionalPanel(
            condition = "input.clusTCR2_tabs2 != 'processing1'",

            conditionalPanel(
              condition = "input.run_ClusTCR2 > 0",
              div(downloadButton("download_ClusTCR_labels", "Download Cluster table"),style = "padding-top: 25px;")
            ),
          ),
        ),
        mainPanel(
          width = 9,
          tabsetPanel(
            id = "clusTCR2_tabs",
            tabPanel("Merge Multiple Files",
                     value = "merge",

                     div(downloadButton("downloaddf_multiple_ClusTCR2", "Download table") ,style = "padding-top: 25px;"),

                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),

                     div(DT::DTOutput("DEx_multiple_ClusTCR2")),

            ),
            tabPanel("Clustering",
                     tabsetPanel(id = "clusTCR2_tabs2",
                                 tabPanel("Inputs",
                                          value = "processing1",
                                          div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                          div(DT::DTOutput("clust_dt2_table")),
                                 ),
                                 tabPanel(
                                   "Processing",
                                   div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                   verbatimTextOutput("ClusTCR2_Time"),
                                   verbatimTextOutput("verbatum_ClusTCR2")
                                 ),
                                 tabPanel(
                                   "Table for analysis",
                                   div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                   div(DT::DTOutput("ClusTCR2_lab")),
                                   p(" "),
                                 ),
                                 tabPanel(
                                   "Figures",
                                   # cluster number
                                   fluidRow(
                                     column(3, numericInput("selected_Cluster", "Selected cluster", value = 1)),
                                     column(3, numericInput("filter_connections", "Keep connections >", value = 1)),
                                     # Name (CDR3_V_gene_Cluster), cluster, CDR3, V_gene, Len (length of CDR3 sequence),CDR3_selected,Name_selected,cluster_selected, (_selected only prints names of the chosen cluster), None
                                     column(3, selectInput("lab_clust_by", "Label cluster by:", choices = c("Name", "cluster", "CDR3", "V_gene", "Len", "CDR3_selected", "Name_selected", "cluster_selected", "None"), selected = "cluster")),
                                     # column(3, selectInput("Clust_size_order", "Order of cluster", choices = c("cluster", "Original_cluster", "Clust_size_order"), selected = "Clust_size_order")),
                                     column(3, selectInput("colour_ClusTCR2", "Type of colouring", choices = c("color_all", "color_test"), selected = "color_test")),
                                     column(3, numericInput("text_size1", "Text size of selected cluster", value = 4)),
                                     column(3, numericInput("text_size2", "Text size of non-selected cluster", value = 2)),
                                   ),
                                   fluidRow(
                                     conditionalPanel(
                                       condition = "input.colour_ClusTCR2 == 'color_all'",
                                       column(3, selectInput("colour_ClusTCR2_types", "Colour panel", choices = c("default", "rainbow", "random", "heat.colors", "terrain.colors", "topo.colors", "hcl.colors"))),
                                     ),
                                     conditionalPanel(
                                       condition = "input.colour_ClusTCR2 == 'color_test'",
                                       fluidRow(
                                         column(3, colourInput("sel_colour_netplot", "Selected colour", "purple")),
                                         column(3, colourInput("nonsel_colour_netplot", "Non-selected colour", "grey80")),
                                       )
                                     )
                                   ),
                                   fluidRow(
                                     column(8, plotOutput("NP_ClusTCR", height = "600px")),
                                     column(4, plotOutput("MP_ClusTCR", height = "200px")),
                                   ),
                                   fluidRow(
                                     column(1, numericInput("width_Network_plot2", "Width of PDF", value = 10)),
                                     column(1, numericInput("height_Network_plot2", "Height of PDF", value = 8)),
                                     column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_Network_plot2", "Download PDF")),
                                     column(2, numericInput("width_png_Network_plot2", "Width of PNG", value = 1200)),
                                     column(2, numericInput("height_png_Network_plot2", "Height of PNG", value = 1000)),
                                     column(2, numericInput("resolution_PNG_Network_plot2", "Resolution of PNG", value = 144)),
                                     column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_Network_plot2", "Download PNG")),
                                   ),
                                   fluidRow(
                                     column(1, numericInput("width_Motif_plot2", "Width of PDF", value = 10)),
                                     column(1, numericInput("height_Motif_plot2", "Height of PDF", value = 3.5)),
                                     column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_Motif_plot2", "Download PDF")),
                                     column(2, numericInput("width_png_Motif_plot2", "Width of PNG", value = 1200)),
                                     column(2, numericInput("height_png_Motif_plot2", "Height of PNG", value = 600)),
                                     column(2, numericInput("resolution_PNG_Motif_plot2", "Resolution of PNG", value = 144)),
                                     column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_Motif_plot2", "Download PNG"))
                                   ),
                                 ),
                                 # tabPanel("Time",

                                 # )
                     )
            )
          )
        )
      )
    ),
  ),
  # end TCR clustering ------
  # STEP 3Quality control side bar panel -----
  navbarMenu(
    "STEP 3.",
    # STEP 3A ------
    tabPanel(
      "3a. Seurat QC",
      sidebarLayout(
        sidebarPanel(shinyjs::useShinyjs(),
                     id = "side-panel",
                     style = "overflow-y:scroll; max-height: 800px; position:relative;", width = 3,
                     div(class = "name-BD",textInput("project_name", h4("Sample/Project Name",class = "name-header2"), value = "")),
                     # selectInput("dataset_sc", "Choose a dataset:", choices = c("test_data_sc", "own_data_sc")),
                     # upload the file
                     fileInput("file_SC", "Load count matrix (csv for BDrhap, csv.gz for 10x, .h5 (10x)", ),
                     fileInput("file_SC_meta", "Upload file meta.data file (.csv.gz or .csv)", ),

                     # selectInput("species","Species",choices = c("human","mouse","other")),
                     selectInput("df_seruatobj_type", "Data type", choices = c("10x_Genomics (raw)", "10x_Genomics (.h5)","BD Rhapsody (Human Immune panel)","BD Rhapsody (Full panel)", "BD Rhapsody (Mouse)", "Array")),
                     selectInput("stored_in_expression", "Does the .h5 object has multiple part?", choices = c("no", "yes")),
                     uiOutput("feature_input"),
                     actionButton('run_violin', 'Filter', onclick = "$(tab).removeClass('disabled-1')"),

                     ######
                     # actionButton("run_violin", "Update Violin plot", onclick = "$(tab1).removeClass('disabled_after')"),
                     conditionalPanel(
                       condition = ("input.run_violin != 0"),
                       fluidRow(
                         column(6, numericInput("dimension_sc", "Max dimensions for clustering", value = 15)),
                         column(6, numericInput("resolution", "Resolution of clusters", value = 1)),
                       ),

                       div(actionButton("run_reduction", "Run clustering", onclick = "$(tab).removeClass('disabled-2')"), style = "padding-top: 25px;"),

                       conditionalPanel(
                         condition = ("input.run_reduction != 0"),
                         div(actionButton("run_metadata", "Impute metadata after clustering", onclick = "$(tab).removeClass('disabled-3')"), style = "padding-top: 25px;"),
                       ),
                     ),
                     conditionalPanel(
                       condition = ("input.run_metadata != 0"),
                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                       downloadButton("downloaddf_SeruatObj", "Download Seurat (.rds)"), style = "padding-top: 25px;"),

                     div(actionButton("clearDataBtn", "Clear Data"), # Add clear data button here
                     ),
                     tags$hr(),
                     # actionButton("reset_input", "Reset inputs")
        ),
        # QC main panel -----
        mainPanel(
          width = 9,
          tabsetPanel(
            # tabPanel("Instructions"
            #          ),
            tabPanel(
              "Uploaded Files",
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              div(DT::DTOutput("DEx_header_name_check.dt")),
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              div(DT::DTOutput("DEx_view.meta.dt")),
            ),
            ######
            tabPanel(
              "Violin and correlation",
              tabsetPanel(
                id = "Before_after",
                type = "tabs",
                tabPanel(
                  "Before",
                  div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                  plotOutput("before_plot_sc", height = "600px"),
                  fluidRow(
                    column(1, numericInput("width_before_plot_sc", "Width of PDF", value = 10)),
                    column(1, numericInput("height_before_plot_sc", "Height of PDF", value = 8)),
                    column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_before_plot_sc", "Download PDF")),
                    column(2, numericInput("width_png_before_plot_sc", "Width of PNG", value = 1200)),
                    column(2, numericInput("height_png_before_plot_sc", "Height of PNG", value = 1000)),
                    column(2, numericInput("resolution_PNG_before_plot_sc", "Resolution of PNG", value = 144)),
                    column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_before_plot_sc", "Download PNG")),
                  ),
                ),
                tabPanel(
                  "After",value = "after_filtering_violin",
                  p("hit 'update Violin plot' to check cut-offs"),
                  div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                  plotOutput("after_plot_sc", height = "600px"),
                  fluidRow(
                    column(1, numericInput("width_after_plot_sc", "Width of PDF", value = 10)),
                    column(1, numericInput("height_after_plot_sc", "Height of PDF", value = 8)),
                    column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_after_plot_sc", "Download PDF")),
                    column(2, numericInput("width_png_after_plot_sc", "Width of PNG", value = 1200)),
                    column(2, numericInput("height_png_after_plot_sc", "Height of PNG", value = 1000)),
                    column(2, numericInput("resolution_PNG_after_plot_sc", "Resolution of PNG", value = 144)),
                    column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_after_plot_sc", "Download PNG")),
                  ),
                  tags$script(
                    '
  var tab2 = $(\'a[data-value="after_filtering_violin"]\').parent();

  $(function(){
    $(tab2).addClass("disabled-1");

    $(tab2.parent()).on("click", "li.disabled-1", function(e) {
      e.preventDefault();
      return false;
    });
  });
  $(document).on("shiny:connected", function(e) {
    $("#run_violin").on("click", function() {
      $(tab2).removeClass("disabled-1");
    });
  });

  // Add this block to re-disable the tab upon clicking clearDataBtn
  $(document).on("shiny:connected", function(e) {
    $("#clearDataBtn").on("click", function() {
      $(tab2).addClass("disabled-1");
    });
  });
  '
                  ),
                )
              )
            ),
  # Elbow and heatmap  -----
  tabPanel(
    "Elbow Plot", value = "Elbow_plot",
    div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
    plotOutput("create_elbowPlot_sc", height = "600px"),
    fluidRow(
      column(1, numericInput("width_create_elbowPlot_sc", "Width of PDF", value = 10)),
      column(1, numericInput("height_create_elbowPlot_sc", "Height of PDF", value = 8)),
      column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_create_elbowPlot_sc", "Download PDF")),
      column(2, numericInput("width_png_create_elbowPlot_sc", "Width of PNG", value = 1200)),
      column(2, numericInput("height_png_create_elbowPlot_sc", "Height of PNG", value = 1000)),
      column(2, numericInput("resolution_PNG_create_elbowPlot_sc", "Resolution of PNG", value = 144)),
      column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_create_elbowPlot_sc", "Download PNG")),
    ),
    tags$script(
      '
    var tab4 = $(\'a[data-value="Elbow_plot"]\').parent();

    $(function(){
      $(tab4).addClass("disabled-1");

      $(tab4.parent()).on("click", "li.disabled-1", function(e) {
        e.preventDefault();
        return false;
      });
    });
    $(document).on("shiny:connected", function(e) {
      $("#run_violin").on("click", function() {
        $(tab4).removeClass("disabled-1");
      });
    });

      // Add this block to re-disable the tab upon clicking clearDataBtn
  $(document).on("shiny:connected", function(e) {
    $("#clearDataBtn").on("click", function() {
      $(tab4).addClass("disabled-1");
    });
  });


    '
    ),
  ),

  # Variable features -----
  tabPanel(
    "Top variable features", value = "Top_var_feat",
    div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
    plotOutput("plot_10_features_sc", height = "600px"),
    fluidRow(
      column(1, numericInput("width_plot_10_features_sc", "Width of PDF", value = 10)),
      column(1, numericInput("height_plot_10_features_sc", "Height of PDF", value = 8)),
      column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_plot_10_features_sc", "Download PDF")),
      column(2, numericInput("width_png_plot_10_features_sc", "Width of PNG", value = 1200)),
      column(2, numericInput("height_png_plot_10_features_sc", "Height of PNG", value = 1000)),
      column(2, numericInput("resolution_PNG_plot_10_features_sc", "Resolution of PNG", value = 144)),
      column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_plot_10_features_sc", "Download PNG")),
    ),

    tags$script(
      '
    var tab3 = $(\'a[data-value="Top_var_feat"]\').parent();

    $(function(){
      $(tab3).addClass("disabled-2");

      $(tab3.parent()).on("click", "li.disabled-2", function(e) {
        e.preventDefault();
        return false;
      });
    });
    $(document).on("shiny:connected", function(e) {
      $("#run_reduction").on("click", function() {
        $(tab3).removeClass("disabled-2");
      });
    });
          // Add this block to re-disable the tab upon clicking clearDataBtn
  $(document).on("shiny:connected", function(e) {
    $("#clearDataBtn").on("click", function() {
      $(tab3).addClass("disabled-2");
    });
  });
    '
    ),

  ),

  # UMAP  -----
  tabPanel(
    "UMAP", value = "sc_UMAP_QC",
    div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
    plotOutput("create_UMAP_sc", height = "600px"),
    tags$script(
      '
    var tab5 = $(\'a[data-value="sc_UMAP_QC"]\').parent();

    $(function(){
      $(tab5).addClass("disabled-2");

      $(tab5.parent()).on("click", "li.disabled-2", function(e) {
        e.preventDefault();
        return false;
      });
    });
    $(document).on("shiny:connected", function(e) {
      $("#run_reduction").on("click", function() {
        $(tab5).removeClass("disabled-2");
      });
    });
// Add this block to re-disable the tab upon clicking clearDataBtn
  $(document).on("shiny:connected", function(e) {
    $("#clearDataBtn").on("click", function() {
      $(tab5).addClass("disabled-2");
    });
  });
    '
    ),

  ), # Export a table with meta.data and expression.
tabPanel(
  "Add metadata", value = "Add_meta_data_to_sc",

  div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
  div(DT::DTOutput("DEx_table_meta.data")),

  tags$script(
    '
    var tab6 = $(\'a[data-value="Add_meta_data_to_sc"]\').parent();

    $(function(){
      $(tab6).addClass("disabled-2");

      $(tab6.parent()).on("click", "li.disabled-3", function(e) {
        e.preventDefault();
        return false;
      });
    });
    $(document).on("shiny:connected", function(e) {
      $("#run_reduction").on("click", function() {
        $(tab6).removeClass("disabled-2");
      });
    });
    // Add this block to re-disable the tab upon clicking clearDataBtn
  $(document).on("shiny:connected", function(e) {
    $("#clearDataBtn").on("click", function() {
      $(tab6).addClass("disabled-2");
    });
  });
    '
  ),
),
          ),
        ),
# end of QC -----
      ),
    ),
###################
# Merge multiple Seurat objects -----
tabPanel(
  "3b. Merge & Batch correction",
  sidebarLayout(
    # Sidebar with a slider input
    sidebarPanel(
      id = "tPanel5", style = "overflow-y:scroll; max-height: 800px; position:relative;", width = 3,
      div(class = "name-BD",textInput("project_name2", h4("Name of Project", class = "name-header2"), value = "")),
      conditionalPanel(
        condition = "input.Merging_and_batching == 'Merging_Harmony'",
        selectInput("sample.type.source_merging", "Species", choices = c("hs", "mm")),
        fileInput("file1_rds.file",
                  "Choose .rds files from merging",
                  multiple = TRUE,
                  accept = c("rds", ".rds")
        ),
        checkboxInput("include_additional_genes",p("Add additional genes?",class = "name-header2"),value = F),

        fileInput("file_user_genes",
                  "Choose .csv files that has required genes",
                  multiple = F,
                  accept = c("csv", ".csv","comma")
        ),

        downloadButton("downloaddf_SeruatObj_merged2", "Download Merged Seurat")
      ),
      conditionalPanel(
        condition = "input.Merging_and_batching != 'Merging_Harmony'",
        selectInput("Seruat_version_merge", "Seurat Version", choices = c("V4", "V5"), selected = "V5"),
        selectInput("sample.type.source", "Species", choices = ""),
        fileInput("file1_rds.Merged_data_for_harmony",
                  "Upload .rds file for Batch correction with Harmony",
                  multiple = F,
                  accept = c(".rds", "rds")
        ),
        downloadButton("downloaddf_SeruatObj_merged", "Download Batch corrected Seurat")
      )
    ),

    # Show a plot of the generated distribution
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "Merging_and_batching",
        tabPanel("Merge Files",
                 value = "Merging_Harmony",
                 div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                 tabsetPanel(
                   tabPanel(
                     "uploaded",
                     verbatimTextOutput("testing_mult"),
                   ),
                   tabPanel(
                     "merging",
                     actionButton("run_merging", "run merging"),
                     verbatimTextOutput("testing_mult2")
                   ),
                 ),
        ),

        tabPanel("Batch correction",
                 tabsetPanel(


                   tabPanel(
                     "Variable data",
                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     # verbatimTextOutput("Scaling_check_output"),
                     verbatimTextOutput("testing_mult3"),
                     actionButton("run_var", "Run VariableFeatures"),



                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     verbatimTextOutput("var_harmony_verbrose")
                   ),
                   tabPanel(
                     "Scale data",
                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     actionButton("run_scale", "Run Scale"),
                     div(DT::DTOutput("Tb_scaling_features_for_annotation")),
                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     verbatimTextOutput("scale_harmony_verbrose")
                   ),
                   tabPanel(
                     "PCA",
                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     actionButton("run_PCA", "Run PCA"),
                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     verbatimTextOutput("PCA_harmony_verbrose")
                   ),
                   tabPanel(
                     "harmony",
                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     actionButton("run_harmony", "Run Harmony"),
                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     verbatimTextOutput("harmony_verbrose"),
                   ),
                   tabPanel(
                     "Dimentional reduction",
                     fluidRow(
                       column(3, numericInput("dimension_Merged", "Max number of dimensions", value = 30)),
                       column(6, numericInput("res_merged", "Resolution of clusters", value = 0.5)),
                     ),
                     actionButton("run_reduction_harmony", "Run Dimentional reduction"),
                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     verbatimTextOutput("testing_mult4"),
                   ),
                   tabPanel(
                     "UMAP",
                     div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                     plotOutput("create_UMAP_merged", height = "600px"),
                     fluidRow(
                       column(1, numericInput("width_sc_merged", "Width of PDF", value = 10)),
                       column(1, numericInput("height_sc_merged", "Height of PDF", value = 8)),
                       column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_sc_merged", "Download Network PDF")),
                       column(2, numericInput("width_png_sc_merged", "Width of PNG", value = 1200)),
                       column(2, numericInput("height_png_sc_merged", "Height of PNG", value = 1000)),
                       column(2, numericInput("resolution_PNG_sc_merged", "Resolution of PNG", value = 144)),
                       column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_sc_merged", "Download Network PNG")),
                     ),
                   ),
                 )
        )
      ))
  )
),

###################

# Add annotations -----
tabPanel(
  "3c. Annotations",
  sidebarLayout(
    sidebarPanel(
      id = "tPanel4", style = "overflow-y:scroll; max-height: 800px; position:relative;", width = 3,
      # uiOutput("Detect_version"),
      div(class = "name-BD",textInput("project_name3", h4("Name of Project", class = "name-header2"), value = "")),
      selectInput("Data_types", "Source", choices = c("10x_HS", "BD_HS.Immune.Panel", "BD_HS.Full.Panel", "10x_MM", "BD_MM_Full.Panel", "BD_MM_Immune.Panel", "TCR-seq")),
      selectInput("sample.type.source.markers", "Species", choices = ""),
      fileInput("file1_rds.file2",
                "Choose merged or single .rds files from directory",
                multiple = TRUE,
                accept = c(".rds", "rds")
      ),

      selectInput("Require_custom_geneset", "Require custom genes?", choices = c("no", "yes")),
      uiOutput("scGate_cutoffs"),
      downloadButton("downloaddf_SeruatObj_annotated", "Download Annotated Seurat"),
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel(
          "Upload",
          div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
          verbatimTextOutput("testing_mult_anno")
        ),
        tabPanel(
          "scGATE",
          selectInput("reduction_anno", "Reduction to use", choices = c("calculate", "pca", "umap", "harmony"), selected = "harmony"),

          # custom annotations databases -----
          conditionalPanel(
            condition = "input.Require_custom_geneset == 'yes'",
            fluidRow(
              column(2, checkboxInput("GeneSet1_scGate", "GeneSet1", value = F)),
              column(2, textInput("geneset1_name", "Name", value = "GeneSet1")),
              column(2, checkboxInput("GeneSet2_scGate", "GeneSet2", value = F)),
              column(2, textInput("geneset2_name", "Name", value = "GeneSet2")),
              column(2, checkboxInput("GeneSet3_scGate", "GeneSet3", value = F)),
              column(2, textInput("geneset3_name", "Name", value = "GeneSet3")),
              column(2, checkboxInput("GeneSet4_scGate", "GeneSet4", value = F)),
              column(2, textInput("geneset4_name", "Name", value = "GeneSet4")),
              column(2, checkboxInput("GeneSet5_scGate", "GeneSet5", value = F)),
              column(2, textInput("geneset5_name", "Name", value = "GeneSet5")),
              column(2, checkboxInput("GeneSet6_scGate", "GeneSet6", value = F)),
              column(2, textInput("geneset6_name", "Name", value = "GeneSet6")),
              column(2, checkboxInput("GeneSet7_scGate", "GeneSet7", value = F)),
              column(2, textInput("geneset7_name", "Name", value = "GeneSet7")),
              column(2, checkboxInput("GeneSet8_scGate", "GeneSet8", value = F)),
              column(2, textInput("geneset8_name", "Name", value = "GeneSet8")),
              column(2, checkboxInput("GeneSet9_scGate", "GeneSet9", value = F)),
              column(2, textInput("geneset9_name", "Name", value = "GeneSet9")),
            )
          ),


          # human 10x annotations -----
          conditionalPanel(
            condition = "input.Data_types == '10x_HS' || input.Data_types == 'BD_HS.Full.Panel' || 'BD_HS.Immune.Panel'",
            div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
            fluidRow(
              div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
              column(3, checkboxInput("hs_function_scGATE", "Function (Human)", value = F)),
              column(3, checkboxInput("hs_generic_scGATE", "Generic (Human)", value = F)),
              column(3, checkboxInput("hs_exhausted_scGATE", "Exhausted (Human)", value = F)),
              column(3, checkboxInput("hs_senescence_scGATE", "Senescence (Human)", value = F)),
              column(3, checkboxInput("hs_cycling_scGATE", "Cycling (Human)", value = F)),
              column(3, checkboxInput("hs_TCRseq_scGATE", "TCR-seq (Human)", value = F)),
            )
          ),

          # BD rhapsody MM full panel ----
          conditionalPanel(
            "input.Data_types == 'BD_MM_Full.Panel' || input.Data_types =='10x_MM'",
            h5("Under development"),
            div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
            h6("Curated from Sharland lab"),
            fluidRow(
              column(3, checkboxInput("BDrhapsody_scGATE.MM.Tcell", "Major T cell popualtions", value = F)),
              column(3, checkboxInput("BDrhapsody_scGATE.MM.Memory", "Memory", value = F)),
              column(3, checkboxInput("BDrhapsody_scGATE.MM.signatures", "signatures", value = F)),
              column(3, checkboxInput("BDrhapsody_scGATE.MM.Innate.NK", "Innate & NK", value = F)),
            ),
            h6("From the human BD Rhapsody searches"),
            fluidRow(
              column(3, checkboxInput("BDrhapsody_scGATE.MM.TNF.IFNg", "TNF.IFNg", value = F)),
              column(3, checkboxInput("BDrhapsody_scGATE.MM.subtypes", "Subtypes", value = F)),
              column(3, checkboxInput("BDrhapsody_scGATE.MM.other", "other", value = F)),
            ),
          ),
          # BD rhapsody MM immune panel ----
          conditionalPanel(
            "input.Data_types == 'BD_MM_Immune.Panel'",
            h5("Under development")
          ),
          conditionalPanel(
            condition = "input.Require_custom_geneset == 'yes'",
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
          conditionalPanel(
            condition = "input.Data_types == '10x_HS' || input.Data_types == 'BD_HS.Full.Panel' || 'BD_HS.Immune.Panel'",
            div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
            verbatimTextOutput("scGATE_verbatum_function"),
            div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
            verbatimTextOutput("scGATE_verbatum_generic2"),
            div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
            verbatimTextOutput("scGATE_verbatum_exhausted"),
            div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
            verbatimTextOutput("scGATE_verbatum_senescence"),
            div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
            verbatimTextOutput("scGATE_verbatum_cycling"),
            div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
            verbatimTextOutput("scGATE_verbatum_TCRseq"),
          ),
          conditionalPanel(
            "input.Data_types == 'BD_MM_Full.Panel' || input.Data_types =='10x_MM'",
            verbatimTextOutput("scGATE_verbatum_BDrhapsody_MM.FP.Tcell"),
            verbatimTextOutput("scGATE_verbatum_BDrhapsody_MM.FP.Memory"),
            verbatimTextOutput("scGATE_verbatum_BDrhapsody_MM.FP.signatures"),
            verbatimTextOutput("scGATE_verbatum_BDrhapsody_MM.FP.Innate.NK"),
            verbatimTextOutput("scGATE_verbatum_BDrhapsody_MM.FP.TNF.IFNg"),
            verbatimTextOutput("scGATE_verbatum_BDrhapsody_MM.FP.subtypes"),
            verbatimTextOutput("scGATE_verbatum_BDrhapsody_MM.FP.other")
          ),
        ),
        tabPanel(
          "Check custom scGate annotations",

          fluidRow(
            column(4, plotOutput("create_custom_1", height = "600px")),
            column(4, plotOutput("create_custom_2", height = "600px"), ),
            column(4, plotOutput("create_custom_3", height = "600px")),
            column(4, plotOutput("create_custom_4", height = "600px"), ),
            column(4, plotOutput("create_custom_5", height = "600px")),
            column(4, plotOutput("create_custom_6", height = "600px"), ),
            column(4, plotOutput("create_custom_7", height = "600px")),
            column(4, plotOutput("create_custom_8", height = "600px"), ),
            column(4, plotOutput("create_custom_9", height = "600px")),
          )
          # plotOutput("create_custom_1", height = "600px"),
        ),
        # classification based on TCR_seq -----
        # display metadata -----
        tabPanel(
          "Feature plot",
          #   fluidRow(
          #     column(3, numericInput("min.ptc.sc", "minimum point", value = 0.25)),
          #     column(3, numericInput("logfc.ptc.sc", "Log fold change", value = 0.25)),
          #     column(3, selectInput("normalN", "Type of Differnetial expression",
          #       choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")
          #     )),
          #     column(3, style = "margin-top: 25px;", actionButton("run_differental.exp", "run differental expression"), )
          # ),
          #   fluidRow(
          #     column(4, selectInput("multiple_group_sc", "Include group comparison", choices = c("no", "yes"))),
          #     column(4, selectInput("meta_data_sc_clust", "Cluster by", choices = "")),
          #     column(4, selectInput("meta_data_sc_", "Add group", choices = "")),
          #
          # ),
          # tabsetPanel(
          # id = "Panel_DEX",
          # Cluster table -----
          # tabPanel("Checking files",
          #          div(DT::DTOutput("list_of_genes")),
          #
          #          # verbatimTextOutput("checking_files_markers_featurePlot_sc"),
          #          ),
          # tabPanel(
          p("   "),
          actionButton("run_string.data3", "View Feature plot"),
          fluidRow(column(12, selectInput("string.data3", "column names for summary", "", multiple = T, width = "1200px"))),
          fluidRow(
            column(2, checkboxInput("label_is_true_features", "Add plot lables", value = T)),
            column(2, selectInput("norm_expression_for_all", "Set Maximum", choices = c("no", "yes"))),
            column(2, numericInput("max_norm_FP", "Set maximum scale value", value = 10, step = 1, min = 1)),
            column(2, colourInput("lower_col_FP", "Min (Colour)", value = "grey90")),
            column(2, colourInput("upper_col_FP", "Max (colour)", value = "Darkblue"))
          ),
          plotOutput("markers_featurePlot_sc", height = "600px"),
          fluidRow(
            column(3, numericInput("width_markers_featurePlot_sc", "Width of PDF", value = 10)),
            column(3, numericInput("height_markers_featurePlot_sc", "Height of PDF", value = 8)),
            column(3),
            column(3, style = "margin-top: 25px;", downloadButton("downloadPlot_markers_featurePlot_sc", "Download PDF"))
          ),
          fluidRow(
            column(3, numericInput("width_png_markers_featurePlot_sc", "Width of PNG", value = 1200)),
            column(3, numericInput("height_png_markers_featurePlot_sc", "Height of PNG", value = 1000)),
            column(3, numericInput("resolution_PNG_markers_featurePlot_sc", "Resolution of PNG", value = 144)),
            column(3, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_markers_featurePlot_sc", "Download PNG"))
          )
          # ),
          # differential expression within clusters ----
          # tabPanel("Treatment differences within clusters",
          #   value = 55,
          #   actionButton("run_update_clust", "Update comparisons"),
          #   fluidRow(
          #     column(4, selectInput("unique.Idents1", "comaprison 1", choices = "")),
          #     column(4, selectInput("unique.Idents2", "comaprison 2", choices = "")),
          #   ),
          #   tabsetPanel(
          #     tabPanel(
          #       "Table",
          #       div(DT::DTOutput("DEx_table_comparison")),
          #       downloadButton("downloaddf_DEx_sc", "Download Table (.csv)"),
          #       downloadButton("downloaddf_DEx_sc_ggVolcanoR", "Download ggVolcanoR compatible table (.csv)")
          #     ),
          #     tabPanel(
          #       "Plot",
          #       plotOutput("volc_plot_cluster", height = "600px")
          #     )
          #   ),
          # ),
          # tabPanel("Cluster differences (All markers)",
          #   value = 5,
          #   div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
          #   div(DT::DTOutput("DEx_table_clusters")),
          #   downloadButton("downloaddf_DEx_table_clusters", "Download Table (.csv)")
          # ),
          # )
        ),

        # meta data table ------
        tabPanel(
          "Meta data table",
          fluidRow(
            # column(3,checkboxInput("add.kmeans","Add K-means classification", value = F)),
            # column(3,checkboxInput("add.scGATE","Add scGATE classifications", value = T))
          ),
          div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
          div(DT::DTOutput("DEx_table_TcellClass_scGATE")),
        )
      ),
    )
  )
),
# remove cells based on one factor -----
tabPanel(
  "3d. Remove/Edit Samps",
  sidebarLayout(
    sidebarPanel(
      id = "tPanelSamps", style = "max-height: 800px; position:relative;", width = 3,
      div(class = "name-BD",textInput("project_name4", h4("Name of Project", class = "name-header2"), value = "")),
      fileInput("file1_rds.fileSampsRemove",
                "Upload .rds file",
                multiple = F,
                accept = c(".rds", "rds")
      ),
      selectInput("Samp_col_SampToRemove", "Column name", choices = ""),
      downloadButton("downloaddf_SeruatObj_annotated_SampToKeep", "Download .rds"),
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel(
          "Remove Samps",
          h5("Before Samples are removed"),
          actionButton("run_remove_samps", "Remove samples"),
          verbatimTextOutput("Preliminary_samp_to_remove"),
          selectInput("DownVColumn", "Chose subset type:", choices = c("Meta_data", "Down_sampling")),
          numericInput("downsamp_limit", "Down sampling limit", value = 1000),
          selectInput("ID_Column_factor_SampToRemove", "Order of graph", choices = "", multiple = T, width = "1400px"),
          h5("After Samples are removed"),
          verbatimTextOutput("Filtered_samp_to_remove"),
        ),
        tabPanel("Update Meta data")
      )
    )
  )
),
# tabPanel("STEP3. Markdown"),
  ),

###################
# Analysis (UI side panel) ---------
tabPanel(
  "STEP 4. Analysis", id = "step4_analysis",

  # side bar layout ------
  sidebarLayout(
    sidebarPanel(
      ######
      id = "tPanel4", style = "overflow-y:scroll; max-height: 1000px; position:relative;", width = 3,
      conditionalPanel(
        condition = "input.check_up_files== 'up'",
        fileInput("file_SC_pro", "Upload seurat file",
                  accept = c("rds", ".rds", "rds")
        ),
        selectInput("add_additional_lables", "add additional labels", choices = c("no", "yes")),
        conditionalPanel(
          condition = "input.add_additional_lables== 'yes'",
          p("The .csv file first column should be label 'ID' and match the selected column"),
          selectInput("Samp_col2", "Sample column name", choices = ""),
          fileInput("file_Labels_to_add", "Upload other identifiers (.csv)",
                    accept = c(".csv", "csv")
          ),
        ),
        selectInput("Type_of_receptor", "Type of receptor", choices = c("TCR", "BCR"), selected = "TCR"),
        conditionalPanel(
          condition = "input.Type_of_receptor== 'TCR'",
          fileInput("file_cluster_file_AG", "Upload AG clusTCR2 file (.csv)",
                    accept = c(".csv", "csv")
          ),
          fileInput("file_cluster_file_BD", "Upload BD clusTCR2 file (.csv)",
                    accept = c(".csv", "csv")
          ),
        ),
        conditionalPanel(
          condition = "input.Type_of_receptor== 'BCR'",
          fileInput("file_cluster_file_IgH", "Upload IgH clusTCR2 file (.csv)",
                    accept = c(".csv", "csv")
          ),
          fileInput("file_cluster_file_IgKL", "Upload IgKL clusTCR2 file (.csv)",
                    accept = c(".csv", "csv")
          ),
        ),
        numericInput("skip_TCRex_up", "Skip # of lines for TCRex file", value = 7),
        fileInput("upload_TCRex_file", "Upload TCRex (.tsv)",
                  accept = c("tsv", ".tsv")
        ),
        selectInput("datasource", "Data source", choices = ""),
        selectInput("species_analysis", "Species", choices = ""),
      ),
      selectInput("V_gene_sc", "V gene with/without CDR3", choices = ""),

      conditionalPanel(
        condition = "input.check_up_files != 'up'",
        fluidRow(
          column(6,  selectInput("Samp_col", "Selected Individual", choices = "")),
          column(6, selectInput("Split_by_group", "Display by Selected Individual", choices = c("no", "yes"))),
          column(6, numericInput("wrap_row", h5("Number of plot rows"), value = 2)),
          conditionalPanel(
            condition = "input.Panel_TCRUMAP != 'TCR_and_GEX_tb'",
            column(6, selectInput("Split_group_by_", "Split graph by:", choices = ""))

          ),

        ),
        selectInput("colourtype", "Colouring Palettes", choices = c("default", "rainbow", "random", "heat.colors", "terrain.colors", "topo.colors", "hcl.colors", "one")),
      ),

      conditionalPanel(
        condition = "input.check_up_files != 'up'",
        conditionalPanel(
          condition = "input.Panel_TCRUMAP != 'Expanded'",
          conditionalPanel(
            condition = "input.check_up_files != 'up2'",
            uiOutput("classification_to_add")
          )
        )
      ),
      # Need to check the colouring by, may need to reduce to 1?
      conditionalPanel(
        condition = "input.check_up_files == 'up2'",
        uiOutput("classification_to_add_overview")
      ),

      # side bar layout expansion priority UI -------
      conditionalPanel(
        condition = "input.check_up_files == 'Prior' || input.check_up_files == 'TCR_and_GEX_tb' ",
        conditionalPanel(
          condition = "input.Panel_TCRUMAP == 'Expanded'",
          h4("Expanded cut-offs and colouring"),
          fluidRow(

            column(6, numericInput("cutoff.expanded", "Cut off greater than", value = 0.5, step = 0.01, min = 0, max = 0.99)),
            column(6, uiOutput("cut.off_expanded2")),
            column(6, uiOutput("classification_to_add2")),
          ),
          fluidRow(
            column(6, numericInput("cut.off_percent_rep", "Percent of Repertoire", value = 1, step = 1, min = 1, max = 100)),
            column(6, numericInput("size.dot.umap", "size of UMAP dot's", value = 2, step = 1, min = 1))
          )
        ),
      ),
      conditionalPanel(
        condition = "input.PriorTBMods == 'PriorRepertoireTB' || input.Panel_TCRUMAP == 'Expanded'",
        fluidRow(
          column(12, uiOutput("Expanded.dotplot.cutoffs")),
          column(12, selectizeInput("selected_Indiv_Ex_1", "Samp 1", choices = "", multiple = T)),
          column(12, selectizeInput("selected_Indiv_Ex_2", "Samp 2", choices = "", multiple = T)),
        )
      ),
      conditionalPanel(
        condition = "input.PriorTBMods == 'PriorClustTB' || input.Panel_TCRUMAP == 'ClusTCR2'",
        selectInput("Clusters_to_dis_PIE", "Clusters to display", choices = "", multiple = F)
      ),
      conditionalPanel(
        condition = "input.Panel_TCRUMAP=='ClusTCR2'",
        selectInput("chain_TCR", "Chain to display", choices = c("TRAG", "TRBD", "IgH", "IgKL")),
      ),
      #
      # Expanded stat cut-offs -----

      conditionalPanel(
        condition = "input.PriorTBMods == 'PriorRepertoireTB' || input.check_up_files == 'TCR_and_GEX_tb' ",
        fluidRow(
          column(4, numericInput("min_point_", "Min point cut off", value = 0.25)),
          column(4, numericInput("LogFC_", "Min LogFC cut off", value = 0.25)),
          column(4, numericInput("pval.ex.filter", "adj.p-val cut-off", value = 0.1)),
        )
      ),
      conditionalPanel(
        condition = "input.Panel_TCRUMAP == 'Marker'",
        column(12, selectInput("col_marker_scale", "Colour scale", choices = col_markers, selected = col_markers[1])),
      ),
      conditionalPanel(
        condition = "input.check_up_files != 'up' ",
        fluidRow(
          column(6, colourInput("one.colour.default", "One colour", "grey50")),
          column(6, colourInput("NA_col_analysis", "NA colour", "grey90"), )
        ),
        conditionalPanel(
          condition = "input.check_up_files == 'up2' || input.check_up_files == 'Prior'",

          fluidRow(column(12, selectInput("Graph_type_bar", "Type of graph", choices = c("Number_expanded", "Frequency_expanded", "Top_clonotypes")))),
          conditionalPanel(
            condition = "input.Graph_type_bar == 'Top_clonotypes'",
            fluidRow(
              column(12, numericInput("top_no_clonotypes", "Top clonotypes per group", value = 1, step = 1, min = 0, max = 20))
            )),
        ),

        h4("What individuals to include"),

        fluidRow(
          column(6, selectInput("by_indiv_pie_epi", "Display one individual?", choices = c("no", "yes"))),
          column(6, selectInput("selected_Indiv", "Display one individual", choices = ""), )
        ),
        h4("Plot parameters (all)"),
        fluidRow(
          column(6, numericInput("text_size", "Axis number size", value = 16)),
          column(6, numericInput("title.text.sizer2", "Axis text size", value = 30)),
          column(6, numericInput("Legend_size", "Legend text size", value = 12)),
          column(6, selectInput("legend_position", "Legend location", choices = c("top", "bottom", "left", "right", "none"), selected = "right")),
        ),
        fluidRow(
          column(6, numericInput("Strip_text_size", "Strip text size (e.g., grey bars)", value = 16)),
          column(6, numericInput("anno_text_size", "Annotation text size", value = 6)),
        ),
        selectInput("font_type", label = "Type of font", choices = font, selected = "Times New Roman"),
      ),
      ###### state of the analysis section as a csv file --------

      # Upload button
      # fileInput("uploadData", "Upload State")
    ),

    # add in clustering  (why did I add this comment?) -----
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "check_up_files",
        tabPanel("Uploaded data",
                 value = "up",
                 div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                 fluidRow(
                   column(12, div(DT::DTOutput("meta.data_check_upload"))),
                   column(12, div(DT::DTOutput("Sample_names_merging_sc"))),
                   # column(12,  div(DT::DTOutput("Tb_TCR_clonotypes.Umap"))),
                   column(12, div(DT::DTOutput("Tb_ClusTCR_test"))),
                   column(12, div(DT::DTOutput("Tb_tcrex_test")))
                 ),
        ),

        # UMAP -> TCR -----
        tabPanel("Overview",
                 value = "up2",
                 fluidRow(
                   div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                   column(12, selectInput("ID_Column_factor", "Order of graph (Selected Individual)", choices = "", multiple = T, width = "1200px")),
                 ),

                 tabsetPanel(
                   id = "QC_panel",
                   # T cell classification ------
                   tabPanel("GEX",
                            value = "GEX_panel",
                            tabsetPanel(
                              id = "Panel_class",
                              tabPanel("Percentage",
                                       value = 16,
                                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                       verbatimTextOutput("Percent_tab"),
                                       downloadButton("downloaddf_Percent_tab", "Download table")
                              ),
                              tabPanel("UMAP plot",
                                       value = 14,
                                       fluidRow(
                                         column(3, selectInput("show_selected", "Show all labels?", choices = c("All", "Selected_list"))),
                                         column(9, uiOutput("SiteNumInput", width = "900px")),
                                       ),
                                       fluidRow(
                                         column(2, numericInput("Filter_lower_UMAP1_marker_GEX", "UMAP_1 >", value = -20)),
                                         column(2, numericInput("Filter_lower_UMAP1_marker2_GEX", "UMAP_1 <", value = 20)),
                                         column(2, numericInput("Filter_lower_UMAP2_marker_GEX", "UMAP_2 >", value = -20)),
                                         column(2, numericInput("Filter_lower_UMAP2_marker2_GEX", "UMAP_2 <", value = 20)),
                                       ),
                                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                       fluidRow(
                                         column(
                                           3,
                                           wellPanel(
                                             id = "tPanel23", style = "overflow-y:scroll; max-height: 600px",
                                             uiOutput("cols_UMAP_all_classification")
                                           )
                                         ),
                                         column(9, plotOutput("UMAP_all_classification2", height = "600px"))
                                       ),
                                       fluidRow(
                                         column(1, numericInput("width_UMAP_all_classification", "Width of PDF", value = 10)),
                                         column(1, numericInput("height_UMAP_all_classification", "Height of PDF", value = 8)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_UMAP_all_classification", "Download PDF")),
                                         column(2, numericInput("width_png_UMAP_all_classification", "Width of PNG", value = 1200)),
                                         column(2, numericInput("height_png_UMAP_all_classification", "Height of PNG", value = 1000)),
                                         column(2, numericInput("resolution_PNG_UMAP_all_classification", "Resolution of PNG", value = 144)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_UMAP_all_classification", "Download PNG"))
                                       ),
                              ),
                              tabPanel("Pie chart",
                                       value = 15,
                                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                       fluidRow(
                                         # column(3),
                                         # column(9, h5("pie segments (Colour by:) & separation (Split graph by:)"))
                                       ),
                                       fluidRow(
                                         column(
                                           3,
                                           wellPanel(
                                             id = "tPanel23", style = "overflow-y:scroll; max-height: 600px",
                                             uiOutput("myPanel_pie")
                                           )
                                         ), # pie chart
                                         column(9, plotOutput("Classification_clonotype_pie", height = "600px"))
                                       ),
                                       # fluidRow(
                                       #   div(DT::DTOutput("table_pie")),
                                       # ),

                                       fluidRow(
                                         column(1, numericInput("width_Classification_clonotype_pie", "Width of PDF", value = 10)),
                                         column(1, numericInput("height_Classification_clonotype_pie", "Height of PDF", value = 8)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_Classification_clonotype_pie", "Download PDF")),
                                         column(2, numericInput("width_png_Classification_clonotype_pie", "Width of PNG", value = 1200)),
                                         column(2, numericInput("height_png_Classification_clonotype_pie", "Height of PNG", value = 1000)),
                                         column(2, numericInput("resolution_PNG_Classification_clonotype_pie", "Resolution of PNG", value = 144)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_Classification_clonotype_pie", "Download PNG"))
                                       ),
                              ),
                            ),
                   ),
                   # TCR overview on top of UMAP -----
                   tabPanel(
                     "TCR",
                     tabsetPanel(
                       id = "TCR",
                       tabPanel(
                         "Overlap",
                         tabsetPanel(
                           tabPanel(
                             "Summary Table",
                             selectInput("other_selected_summary_columns","Add other summary columns", multiple = T, "",width = "1200px"),
                             div(DT::DTOutput("Summary_TCR_tb")),
                             downloadButton("downloaddf_Summary_TCR_tb", "Download table")
                           ),

                           tabPanel(
                             "Upset_table",
                             div(DT::DTOutput("Upset_plot_overlap_Tb")),
                             downloadButton("downloaddf_Upset_plot_overlap_Tb", "Download table")
                           ),
                           tabPanel(
                             "Upset Plot",
                             plotOutput("Upset_plot_overlap"),
                             fluidRow(
                               column(1, numericInput("width_Upset_plot_overlap", "Width of PDF", value = 10)),
                               column(1, numericInput("height_Upset_plot_overlap", "Height of PDF", value = 8)),
                               column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_Upset_plot_overlap", "Download PDF")),
                               column(2, numericInput("width_png_Upset_plot_overlap", "Width of PNG", value = 1200)),
                               column(2, numericInput("height_png_Upset_plot_overlap", "Height of PNG", value = 1000)),
                               column(2, numericInput("resolution_PNG_Upset_plot_overlap", "Resolution of PNG", value = 144)),
                               column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_Upset_plot_overlap", "Download PNG"))
                             ),
                           ),
                         )
                       ),
                       ##### line graph ----
                       tabPanel("Line Graph",
                                p(" "),
                                fluidRow(column(2,actionButton("load_samp_name_list","Load samples"))),
                                fluidRow(
                                  column(2,checkboxInput("is_a_time_series",h4("Time series?"),value = T)),
                                  column(2,selectInput("separator_input", "Select Separators:",
                                                       choices = c("_" = "_.*", "-" = "-.*", "." = "\\..*", "|" = "\\|.*", "#" = "#*", "^" = "\\^*", "&" = "&.*"),
                                                       multiple = F),),
                                  column(2,selectInput("comparison_operator", "Choose comparison operator:",
                                                       choices = c("Equal to" = "==", "Greater than or equal to" = ">="),
                                                       selected = c("Greater than or equal to" = ">=")),),
                                  column(2,numericInput("cutoff_upset", "Enter cutoff value:", value = 2),),
                                  column(2,numericInput("max_number_lines_to","Maximum to display",value = 20)),
                                  column(2,numericInput("Total_count_Cutoff","Min count threshold",value = 1)),

                                ),

                                fluidRow(

                                  column(2,conditionalPanel(
                                    condition = "input.is_a_time_series",
                                    selectizeInput("Group_for_line_graph","Display multi-sample clones for: ",""))),
                                ),

                                tabsetPanel(
                                  tabPanel("Table",
                                           div(DT::DTOutput("Line_graph_table")),
                                  ),
                                  tabPanel("Line graph",
                                           fluidRow(
                                             column(2,sliderInput("number_of_conditions","Number of conditions",value = 2, min = 2 , max = 3)),
                                             column(2,selectInput("separator_input2", "Select Separators:",
                                                                  choices = c("_" = "_", "-" = "-", "." = "\\.", "|" = "\\|", "#" = "#", "^" = "\\^", "&" = "&"),
                                                                  selected = "_",
                                                                  multiple = T)),
                                             column(2,selectInput("display_all_samps_line","Display all?",choices = c("no","yes"))),

                                             column(3,
                                                    conditionalPanel(
                                                      condition = "input.number_of_conditions == 3",textInput("shape_legend_name","Shape legend name",value = ""))
                                             ),
                                           ),

                                           conditionalPanel(
                                             condition = "input.display_all_samps_line == 'no'",

                                             plotOutput("line_graph_output"),

                                           ),

                                           conditionalPanel(
                                             condition = "input.display_all_samps_line == 'yes'",

                                             plotOutput("line_graph_all_output"),

                                           ),



                                           fluidRow(
                                             column(1, numericInput("width_line_graph_output", "Width of PDF", value = 10)),
                                             column(1, numericInput("height_line_graph_output", "Height of PDF", value = 8)),
                                             column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_line_graph_output", "Download PDF")),
                                             column(2, numericInput("width_png_line_graph_output", "Width of PNG", value = 1200)),
                                             column(2, numericInput("height_png_line_graph_output", "Height of PNG", value = 1000)),
                                             column(2, numericInput("resolution_PNG_line_graph_output", "Resolution of PNG", value = 144)),
                                             column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_line_graph_output", "Download PNG"))
                                           ),
                                  )
                                ),

                       ),
                       #####
                       tabPanel("Clonal expansion plots",
                                value = 2,
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                fluidRow(
                                  conditionalPanel(
                                    condition = "input.Graph_type_bar=='Number_expanded' || input.Graph_type_bar=='Frequency_expanded'",
                                    fluidRow(
                                      column(
                                        3,
                                        wellPanel(
                                          id = "tPanel23", style = "overflow-y:scroll; max-height: 600px",
                                          uiOutput("myPanel_clonal_plot")
                                        )
                                      ),
                                      column(9, plotOutput("clonality.bar.graph", height = "600px"))
                                    )
                                  ),
                                  # selected clones to display ----
                                  conditionalPanel(
                                    condition = "input.Graph_type_bar=='Top_clonotypes'",
                                    fluidRow(
                                      column(
                                        3,
                                        wellPanel(
                                          id = "tPanel23", style = "overflow-y:scroll; max-height: 600px",
                                          uiOutput("myPanel_top_clonal_plot")
                                        )
                                      ),
                                      column(9, plotOutput("clonality.bar.graph2", height = "600px"))
                                    ),
                                  ),
                                ),
                                fluidRow(
                                  column(1, numericInput("width_clonality.bar.graph", "Width of PDF", value = 10)),
                                  column(1, numericInput("height_clonality.bar.graph", "Height of PDF", value = 8)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_clonaity.bar.graph", "Download PDF")),
                                  column(2, numericInput("width_png_clonality.bar.graph", "Width of PNG", value = 1200)),
                                  column(2, numericInput("height_png_clonality.bar.graph", "Height of PNG", value = 1000)),
                                  column(2, numericInput("resolution_PNG_clonality.bar.graph", "Resolution of PNG", value = 144)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_clonaity.bar.graph", "Download PNG"))
                                ),
                       ),
                       # UMAP clonality -> TCR -----
                       tabPanel("Conality (counts)",
                                value = 3,
                                tabsetPanel(
                                  tabPanel(
                                    "Test.table",
                                    div(DT::DTOutput("Tb_TCR_clonotypes.table"))
                                  ),
                                  tabPanel(
                                    "UMAP",
                                    div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                    fluidRow(
                                      column(3, selectInput("filter_umap_expand", "Filter plot", choices = c("no", "yes"))),
                                    ),
                                    conditionalPanel(
                                      condition = "input.filter_umap_expand == 'yes'",
                                      fluidRow(
                                        column(3, numericInput("UMAP_1x", "UMAP_1 <", value = 10)),
                                        column(3, numericInput("UMAP_1y", "UMAP_1 >", value = -10)),
                                        column(3, numericInput("UMAP_2x", "UMAP_2 <", value = 10)),
                                        column(3, numericInput("UMAP_2y", "UMAP_2 >", value = -10))
                                      ),
                                    ),
                                    conditionalPanel(
                                      condition = "input.Graph_type_bar=='Number_expanded' || input.Graph_type_bar=='Frequency_expanded'",
                                      fluidRow(
                                        column(
                                          3,
                                          wellPanel(
                                            id = "tPanel23", style = "overflow-y:scroll; max-height: 600px",
                                            uiOutput("cols_UMAP_clonal_plot")
                                          )
                                        ),
                                        column(9, plotOutput("clonality.TCR.UMAP", height = "600px"))
                                      ),
                                      fluidRow(
                                        column(1, numericInput("width_TCR.UMAP", "Width of PDF", value = 10)),
                                        column(1, numericInput("height_TCR.UMAP", "Height of PDF", value = 8)),
                                        column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_TCR.UMAP", "Download PDF")),
                                        column(2, numericInput("width_png_TCR.UMAP", "Width of PNG", value = 1600)),
                                        column(2, numericInput("height_png_TCR.UMAP", "Height of PNG", value = 1000)),
                                        column(2, numericInput("resolution_PNG_TCR.UMAP", "Resolution of PNG", value = 144)),
                                        column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_TCR.UMAP", "Download PNG"))
                                      ),
                                    ),
                                    conditionalPanel(
                                      condition = "input.Graph_type_bar=='Top_clonotypes'",
                                      #

                                      fluidRow(
                                        # column(4,selectInput("Split_by_group","Include group comparison",choices=c("no","yes"))),
                                        column(4, selectInput("display_all_samps", "Display all sample", choices = c("yes", "no"))),
                                      ),
                                      conditionalPanel(
                                        condition = "input.display_all_samps == 'no'",
                                        column(8, selectInput("ID_Column_metadata", "Select to display", choices = "", multiple = T, width = "800px"))
                                      ),
                                      div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                      fluidRow(
                                        column(
                                          3,
                                          wellPanel(
                                            id = "tPanel23", style = "overflow-y:scroll; max-height: 600px",
                                            uiOutput("cols_UMAP_Topclonotypes")
                                          )
                                        ),
                                        column(9, plotOutput("clonality.TCR.UMAP.top", height = "600px")),
                                        fluidRow(
                                          column(1, numericInput("width_TCR.UMAP_top", "Width of PDF", value = 10)),
                                          column(1, numericInput("height_TCR.UMAP_top", "Height of PDF", value = 8)),
                                          column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_TCR.UMAP_top", "Download PDF")),
                                          column(2, numericInput("width_png_TCR.UMAP_top", "Width of PNG", value = 1800)),
                                          column(2, numericInput("height_png_TCR.UMAP_top", "Height of PNG", value = 1000)),
                                          column(2, numericInput("resolution_PNG_TCR.UMAP_top", "Resolution of PNG", value = 144)),
                                          column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_TCR.UMAP_top", "Download PNG"))
                                        ),
                                      ),
                                    ),
                                  )
                                ),
                       ),
                     ),
                   ),



                   # tabPanel("Clonotype overlap per cluster (upset plot)",
                   #          # add in upset plot per cluster
                   #          )
                 )
        ),
        # end of differential expression -----
        # TCR and GEX analysis section-----
        tabPanel("TCR -> GEX",
                 value = "TCR_and_GEX_tb",
                 fluidRow(
                   column(
                     12,
                     # conditionalPanel(condition="input.Panel_TCRUMAP=='top_clone' || input.Panel_TCRUMAP=='Epitope'",
                     selectInput("Graph_split_order", "Order of split by:", choices = "", multiple = T, width = "1400px")
                   ),
                   column(3,colourInput("min_FC_col", "Zero colour", value = "white")),
                   column(3,colourInput("med_FC_col", "from one colour", value = "#E9C2FF")),
                   column(3,colourInput("max_FC_col", "Max colour", value = "#6F00B0")),
                   # ),
                 ),
                 # Classification to include ------
                 tabsetPanel(
                   id = "Panel_TCRUMAP",
                   # top clonotypes plot -----
                   tabPanel("Clone abundance",
                            value = "top_clone",
                            fluidRow(
                              column(
                                3,
                                conditionalPanel(
                                  condition = "input.Panel_TCRUMAP=='top_clone'",
                                  checkboxInput("limit_to_top_clones", "Limit to the top clones", value = T)
                                ),
                              ),
                              column(
                                3,
                                conditionalPanel(
                                  condition = "input.Panel_TCRUMAP=='top_clone' && input.limit_to_top_clones",
                                  numericInput("max_top_clone_limit", "1 to max value", value = 50)
                                ),
                              )
                            ),
                            fluidRow(
                              column(
                                12,
                                conditionalPanel(
                                  condition = "input.Panel_TCRUMAP=='top_clone'",
                                  selectInput("Selected_clonotype", "Select clonotype:", choices = "", width = "1400px")
                                ),
                              ),
                            ),
                            tabsetPanel(
                              tabPanel(
                                "Summary table",
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                div(DT::DTOutput("Top_clonotype_sum")),
                                downloadButton("download_Top_clonotype_sum", "Download table")
                              ),
                              tabPanel(
                                "Bar graph",
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                fluidRow(
                                  column(
                                    3,
                                    wellPanel(
                                      id = "tPanel23", style = "overflow-y:scroll; max-height: 600px",
                                      uiOutput("myPanel_Top_bar_clonotype")
                                    )
                                  ),
                                  column(9, plotOutput("top_clonotype", height = "600px"))
                                ),
                                fluidRow(
                                  column(1, numericInput("width_top_clonotype", "Width of PDF", value = 10)),
                                  column(1, numericInput("height_top_clonotype", "Height of PDF", value = 8)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_top_clonotype", "Download PDF")),
                                  column(2, numericInput("width_png_top_clonotype", "Width of PNG", value = 1200)),
                                  column(2, numericInput("height_png_top_clonotype", "Height of PNG", value = 1000)),
                                  column(2, numericInput("resolution_PNG_top_clonotype", "Resolution of PNG", value = 144)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_top_clonotype", "Download PNG"))
                                ),
                              ),
                              # tabPanel("Pie labs",
                              #          add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "200px",width = "200px", color = "#6F00B0"),
                              #          div(DT::DTOutput("Top_clonotype_Labs")),
                              #
                              # ),

                              tabPanel("heatmap",
                                       value = "TopHeat",
                                       p("Colour by: = x-axis; Split graph by: = y-axis for this heatmap"),
                                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                       plotOutput("heatmap_topclone_plot", height = "600px"),
                                       fluidRow(
                                         column(1, numericInput("width_heatmap_topclone_plot", "Width of PDF", value = 10)),
                                         column(1, numericInput("height_heatmap_topclone_plot", "Height of PDF", value = 8)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_heatmap_topclone_plot", "Download PDF")),
                                         column(2, numericInput("width_png_heatmap_topclone_plot", "Width of PNG", value = 1200)),
                                         column(2, numericInput("height_png_heatmap_topclone_plot", "Height of PNG", value = 1000)),
                                         column(2, numericInput("resolution_heatmap_topclone_plot", "Resolution of PNG", value = 144)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_heatmap_topclone_plot", "Download PNG"))
                                       ),
                              ),
                              #####.
                              tabPanel(
                                "Pie/UMAP chart",
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                fluidRow(
                                  column(3, selectInput("Plot_type_selected", "Plot", choices = c("pie", "UMAP"))),
                                  column(3, numericInput("size_selected_top", "Size of Point", value = 2)),
                                ),
                                fluidRow(
                                  column(
                                    3,
                                    wellPanel(
                                      id = "tPanel23", style = "overflow-y:scroll; max-height: 600px",
                                      uiOutput("myPanel_Top_pie_clonotype")
                                    )
                                  ),
                                  column(9, plotOutput("top_clonotype_pie", height = "600px")),
                                ),
                                fluidRow(
                                  column(1, numericInput("width_top_clonotype_pie", "Width of PDF", value = 10)),
                                  column(1, numericInput("height_top_clonotype_pie", "Height of PDF", value = 8)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_top_clonotype_pie", "Download PDF")),
                                  column(2, numericInput("width_png_top_clonotype_pie", "Width of PNG", value = 1200)),
                                  column(2, numericInput("height_png_top_clonotype_pie", "Height of PNG", value = 1000)),
                                  column(2, numericInput("resolution_PNG_top_clonotype_pie", "Resolution of PNG", value = 144)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_top_clonotype_pie", "Download PNG"))
                                ),
                              ),
                              # add in find marker for comparing population to other for top clonotype
                              # tabPanel("FindMarker"),

                              # tabPanel("Expression",



                              # tabsetPanel(
                              tabPanel(
                                "Ridge/Violin plots",
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                actionButton("run_string.data_Exp_top", "View Ridge plot"),
                                fluidRow(column(12, selectInput("string.data_Exp_top", "column names for summary", "", multiple = F, width = "1200px"))),
                                fluidRow(
                                  column(3, checkboxInput("restric_ex", "Restrict to above a threshold?", value = F)),
                                  column(3, numericInput("Gre_ex", "Expression above:", value = 0)),
                                  column(3, selectInput("plot_type_ridgvi", "Plot type", choices = c("Ridge (selected clonotype)", "Ridge (compare)", "Violin (selected clonotype)", "Violin (compare)"))),
                                ),
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                fluidRow(
                                  column(6, div(DT::DTOutput("Ridge_chart_alpha_gamma_stat"))),
                                  column(6, plotOutput("Ridge_chart_alpha_gamma_plot_out", height = "600px"))
                                ),
                                column(6, downloadButton("downloaddf_clusTCR_GEx", "Download stats")),
                                fluidRow(
                                  column(1, numericInput("width_Ridge_chart_alpha_gamma_plot_out", "Width of PDF", value = 10)),
                                  column(1, numericInput("height_Ridge_chart_alpha_gamma_plot_out", "Height of PDF", value = 8)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_Ridge_chart_alpha_gamma_plot_out", "Download PDF")),
                                  column(2, numericInput("width_png_Ridge_chart_alpha_gamma_plot_out", "Width of PNG", value = 1200)),
                                  column(2, numericInput("height_png_Ridge_chart_alpha_gamma_plot_out", "Height of PNG", value = 1000)),
                                  column(2, numericInput("resolution_PNG_Ridge_chart_alpha_gamma_plot_out", "Resolution of PNG", value = 144)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_Ridge_chart_alpha_gamma_plot_out", "Download PNG"))
                                ),
                              ),
                              tabPanel(
                                "Stats",
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                div(DT::DTOutput("Ridge_chart_alpha_gamma_stat_comp")),
                                downloadButton("downloaddf_FindMarker_Top", "Download stat (Right)")
                              ),
                              # dotplot top-----
                              tabPanel(
                                "Dotplot",
                                fluidRow(
                                  column(2, colourInput("low.dotplot", "Lower color:", "darkblue")),
                                  column(2, colourInput("middle.dotplot", "Middle color:", "white")),
                                  column(2, colourInput("high.dotplot", "High color:", "darkred")),
                                  column(2, checkboxInput("restict_no_points", "Restrict Label", value = F)),
                                  column(2, numericInput("pval.ex.top_genes", "Top genes to display", value = 40)),
                                ),
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                plotOutput("all_expression_dotplot_top", height = "400px"),
                                textInput("name_clonotype_selected", "Name of clone", "clone 1"),
                                fluidRow(
                                  column(1, numericInput("width_all_expression_dotplot_top", "Width of PDF", value = 20)),
                                  column(1, numericInput("height_all_expression_dotplot_top", "Height of PDF", value = 4)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_all_expression_dotplot_top", "Download PDF")),
                                  column(2, numericInput("width_png_all_expression_dotplot_top", "Width of PNG", value = 2400)),
                                  column(2, numericInput("height_png_all_expression_dotplot_top", "Height of PNG", value = 700)),
                                  column(2, numericInput("resolution_PNG_all_expression_dotplot_top", "Resolution of PNG", value = 144)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_all_expression_dotplot_top", "Download PNG"))
                                ),
                              ),
                              tabPanel(
                                "Over-representation",
                                fluidRow(
                                  column(3, numericInput("in.geneset.cutoff_top", "Min number of genes in GeneSet", value = 1, min = 0, step = 1, max = 60000)),
                                  column(3, numericInput("p.val_cutoff_top", "p-val cut-off", value = 0.05, min = 0, max = 1)),
                                  # column(3,numericInput("adjust_cutoff_top","BH cut-off",value = 1, min = 0, max = 1)),
                                ),
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                div(DT::DTOutput("Over_rep_Top_clones_Tab")),
                                downloadButton("downloadtb_over.rep.Top_Ex", "Download table")
                                #
                              ),
                              # ),
                              # ),
                            ),
                   ),

                   # expanded phenotype -----
                   tabPanel("Expanded",
                            value = "Expanded",
                            tabsetPanel(
                              id = "ExPan",
                              tabPanel(
                                "Table",
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                div(DT::DTOutput("Expansion_check")),
                              ),
                              tabPanel("ex.UMAP",
                                       value = "ExPan_UMAP",
                                       fluidRow(
                                         div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                         column(
                                           3,
                                           wellPanel(
                                             id = "tPanel23", style = "overflow-y:scroll; max-height: 600px",
                                             uiOutput("cols_UMAP_Expanded")
                                           )
                                         ),
                                         column(9, plotOutput("UMAP_Expanded", height = "600px"))
                                       ),
                                       fluidRow(
                                         column(1, numericInput("width_UMAP_Expanded", "Width of PDF", value = 10)),
                                         column(1, numericInput("height_UMAP_Expanded", "Height of PDF", value = 8)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_UMAP_Expanded", "Download PDF")),
                                         column(2, numericInput("width_png_UMAP_Expanded", "Width of PNG", value = 1200)),
                                         column(2, numericInput("height_png_UMAP_Expanded", "Height of PNG", value = 1000)),
                                         column(2, numericInput("resolution_PNG_UMAP_Expanded", "Resolution of PNG", value = 144)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_UMAP_Expanded", "Download PNG"))
                                       ),
                              ),
                              tabPanel("Stats",
                                       value = "ExPan_stat",
                                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                       div(DT::DTOutput("compare.stat_Ex")),
                                       downloadButton("downloadtb_compare.stat_Ex", "Download table")
                              ),
                              tabPanel("Dotplot",
                                       value = "ExPan_dot",
                                       fluidRow(
                                         column(2, colourInput("low.dotplot.ex", "Lower color:", "darkblue")),
                                         column(2, colourInput("middle.dotplot.ex", "Middle color:", "white")),
                                         column(2, colourInput("high.dotplot.ex", "High color:", "darkred")),
                                       ),
                                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                       plotOutput("relative_expression_dotplot_ex", height = "600px"),
                                       fluidRow(
                                         column(1, numericInput("width_all_expression_dotplot_ex", "Width of PDF", value = 20)),
                                         column(1, numericInput("height_all_expression_dotplot_ex", "Height of PDF", value = 4)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_all_expression_dotplot_ex", "Download PDF")),
                                         column(2, numericInput("width_png_all_expression_dotplot_ex", "Width of PNG", value = 2400)),
                                         column(2, numericInput("height_png_all_expression_dotplot_ex", "Height of PNG", value = 600)),
                                         column(2, numericInput("resolution_PNG_all_expression_dotplot_ex", "Resolution of PNG", value = 144)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_all_expression_dotplot_ex", "Download PNG"))
                                       ),
                              ),
                              tabPanel("Over-representation",
                                       value = "ExPan_OvRep",
                                       fluidRow(
                                         column(3, numericInput("in.geneset.cutoff_Exp", "Min number of genes in GeneSet", value = 1, min = 0, step = 1, max = 60000)),
                                         column(3, numericInput("p.val_cutoff_Exp", "p-val cut-off", value = 0.05, min = 0, max = 1)),
                                       ),
                                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                       div(DT::DTOutput("Over_rep_Exp_Tab")),
                                       downloadButton("downloadtb_over.rep_Exp", "Download table")
                              )
                            )
                   ),

                   # ClusTCR2 Analysis -----
                   tabPanel("ClusTCR2",
                            value = "ClusTCR2",
                            tabsetPanel(
                              tabPanel(
                                "Table.Clust",
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                div(DT::DTOutput("Tb_ClusTCR_selected")),
                                downloadButton("downloadtb_Tb_ClusTCR_selected", "Download table")
                              ),
                              tabPanel(
                                "UMAP",
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                fluidRow(
                                  column(
                                    3,
                                    wellPanel(
                                      id = "tPanel23", style = "overflow-y:scroll; max-height: 600px",
                                      uiOutput("myPanel_cols_clust_UMAP")
                                    )
                                  ),
                                  column(9, plotOutput("UMAP_ClusTCR2_plot", height = "600px"))
                                ),
                                fluidRow(
                                  column(1, numericInput("width_UMAP_ClusTCR2_plot", "Width of PDF", value = 10)),
                                  column(1, numericInput("height_UMAP_ClusTCR2_plot", "Height of PDF", value = 8)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_UMAP_ClusTCR2_plot", "Download PDF")),
                                  column(2, numericInput("width_png_UMAP_ClusTCR2_plot", "Width of PNG", value = 1200)),
                                  column(2, numericInput("height_png_UMAP_ClusTCR2_plot", "Height of PNG", value = 1000)),
                                  column(2, numericInput("resolution_PNG_UMAP_ClusTCR2_plot", "Resolution of PNG", value = 144)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_UMAP_ClusTCR2_plot", "Download PNG"))
                                ),
                              ),
                              tabPanel(
                                "motif",
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                plotOutput("Motif_ClusTCR2_cluster", height = "300px"),
                                verbatimTextOutput("print_unique_cases"),
                                fluidRow(
                                  column(1, numericInput("width_Motif_ClusTCR2_cluster", "Width of PDF", value = 10)),
                                  column(1, numericInput("height_Motif_ClusTCR2_cluster", "Height of PDF", value = 4)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_Motif_ClusTCR2_cluster", "Download PDF")),
                                  column(2, numericInput("width_png_Motif_ClusTCR2_cluster", "Width of PNG", value = 2400)),
                                  column(2, numericInput("height_png_Motif_ClusTCR2_cluster", "Height of PNG", value = 600)),
                                  column(2, numericInput("resolution_PNG_Motif_ClusTCR2_cluster", "Resolution of PNG", value = 144)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_Motif_ClusTCR2_cluster", "Download PNG"))
                                ),
                              ),
                              #####
                              tabPanel("heatmap",
                                       value = "ClustHeat",
                                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                       p("Colour by: = x-axis; Split graph by: = y-axis for this heatmap"),
                                       plotOutput("Pie_ClusTCR2_plot", height = "600px"),
                                       fluidRow(
                                         column(1, numericInput("width_Pie_ClusTCR2_plot", "Width of PDF", value = 10)),
                                         column(1, numericInput("height_Pie_ClusTCR2_plot", "Height of PDF", value = 8)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_Pie_ClusTCR2_plot", "Download PDF")),
                                         column(2, numericInput("width_png_Pie_ClusTCR2_plot", "Width of PNG", value = 1200)),
                                         column(2, numericInput("height_png_Pie_ClusTCR2_plot", "Height of PNG", value = 1000)),
                                         column(2, numericInput("resolution_PNG_Pie_ClusTCR2_plot", "Resolution of PNG", value = 144)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_Pie_ClusTCR2_plot", "Download PNG"))
                                       ),
                              ),

                              tabPanel("Stats",
                                       value = "ClusPan_stat",
                                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                       div(DT::DTOutput("compare.stat_Cluster_DT")),
                                       downloadButton("downloadtb_compare.stat_Cluster", "Download table")
                              ),

                              # dotplot Cluster ------
                              tabPanel("Dotplot",
                                       value = "ClusPan_dot",
                                       fluidRow(
                                         column(2, colourInput("low.dotplot.clust", "Lower color:", "darkblue")),
                                         column(2, colourInput("middle.dotplot.clust", "Middle color:", "white")),
                                         column(2, colourInput("high.dotplot.clust", "High color:", "darkred")),
                                         column(2, checkboxInput("restrict.dotpot.clust", "Restrict to top list", value = F)),
                                         column(2, numericInput("restrict.dotpot.num.clust", "Total genes to display:", value = 10))
                                       ),
                                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                       plotOutput("all_expression_dotplot_cluster", height = "400px"),
                                       fluidRow(
                                         column(1, numericInput("width_all_expression_dotplot_clust", "Width of PDF", value = 20)),
                                         column(1, numericInput("height_all_expression_dotplot_clust", "Height of PDF", value = 8)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_all_expression_dotplot_clust", "Download PDF")),
                                         column(2, numericInput("width_png_all_expression_dotplot_clust", "Width of PNG", value = 2400)),
                                         column(2, numericInput("height_png_all_expression_dotplot_clust", "Height of PNG", value = 700)),
                                         column(2, numericInput("resolution_PNG_all_expression_dotplot_clust", "Resolution of PNG", value = 144)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_all_expression_dotplot_clust", "Download PNG"))
                                       ),
                              ),
                              tabPanel("Over-representation",
                                       value = "ClusPan_OvRep",
                                       fluidRow(
                                         column(3, numericInput("in.geneset.cutoff_Clust", "Min number of genes in GeneSet", value = 1, min = 0, step = 1, max = 60000)),
                                         column(3, numericInput("p.val_cutoff_Clust", "p-val cut-off", value = 0.05, min = 0, max = 1)),
                                         # column(3,numericInput("adjust_cutoff_Clust","BH cut-off",value = 1, min = 0, max = 1)),
                                       ),
                                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                       div(DT::DTOutput("Over_rep_Cluster_Tab")),
                                       downloadButton("downloadtb_over.rep.cluster", "Download table")
                              )
                            )
                   ),
                   # epitope analysis -----
                   tabPanel("Epitope",
                            value = "Epitope",
                            conditionalPanel(
                              condition = "input.Panel_TCRUMAP=='Epitope'",
                              fluidRow(
                                column(3, uiOutput("classification_to_add_epitope")),
                                column(3, uiOutput("classification_to_add_epitope2")),
                                column(3, actionButton("Update_epi", "Add in Epitope list")),
                                column(3, selectInput("Epi_of_interest", "Epitope of interest", ""))
                              ),
                            ),
                            tabsetPanel(
                              id = "EpitipeTabs",
                              # tabPanel("Uploaded Epitope file",
                              #          add_busy_spinner(spin = "fading-circle",position = "top-right",margins = c(10,10),height = "200px",width = "200px", color = "#6F00B0"),
                              #          div(DT::DTOutput("MainTcell_Check")),
                              # ),
                              tabPanel(
                                "Summary Table",
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                div(DT::DTOutput("Pie_Epitope_dt")),
                                downloadButton("downloaddf_Pie_Epitope_dt", "Download table")
                              ),
                              tabPanel(
                                "Heatmap",
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                plotOutput("Heatmap_epi_plot", height = "600px"),
                                fluidRow(
                                  column(1, numericInput("width_Heatmap_epi_plot", "Width of PDF", value = 10)),
                                  column(1, numericInput("height_Heatmap_epi_plot", "Height of PDF", value = 8)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_Heatmap_epi_plot", "Download PDF")),
                                  column(2, numericInput("width_png_Heatmap_epi_plot", "Width of PNG", value = 1200)),
                                  column(2, numericInput("height_png_Heatmap_epi_plot", "Height of PNG", value = 1000)),
                                  column(2, numericInput("resolution_PNG_Heatmap_epi_plot", "Resolution of PNG", value = 144)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_Heatmap_epi_plot", "Download PNG"))
                                ),
                              ),
                              tabPanel(
                                "UMAP",
                                numericInput("value_size_epi_umap", "Size of epitope dots", value = 2),

                                # column(3,selectInput("epitope_umap_selected","Select",choices = c("beta","epitope","pathology"),selected = "pathology")),
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                fluidRow(
                                  column(
                                    3,
                                    wellPanel(
                                      id = "tPanel23", style = "overflow-y:scroll; max-height: 600px",
                                      uiOutput("myPanel_cols_epitope")
                                    )
                                  ),
                                  column(9, plotOutput("UMAP_Epitope_plot", height = "600px"))
                                ),
                                fluidRow(
                                  column(1, numericInput("width_UMAP_Epitope", "Width of PDF", value = 10)),
                                  column(1, numericInput("height_UMAP_Epitope", "Height of PDF", value = 8)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_UMAP_Epitope", "Download PDF")),
                                  column(2, numericInput("width_png_UMAP_Epitope", "Width of PNG", value = 1200)),
                                  column(2, numericInput("height_png_UMAP_Epitope", "Height of PNG", value = 1000)),
                                  column(2, numericInput("resolution_PNG_UMAP_Epitope", "Resolution of PNG", value = 144)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_UMAP_Epitope", "Download PNG"))
                                ),
                              ),
                              tabPanel(
                                "Pie (Expression)",
                                div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                fluidRow(
                                  column(
                                    3,
                                    wellPanel(
                                      id = "tPanel23", style = "overflow-y:scroll; max-height: 600px",
                                      uiOutput("myPanel_cols_epitope_pie")
                                    )
                                  ),
                                  column(9, plotOutput("Pie_Epitope_plot", height = "600px"))
                                ),
                                fluidRow(
                                  column(1, numericInput("width_Pie_Epitope", "Width of PDF", value = 10)),
                                  column(1, numericInput("height_Pie_Epitope", "Height of PDF", value = 8)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_Pie_Epitope", "Download PDF")),
                                  column(2, numericInput("width_png_Pie_Epitope", "Width of PNG", value = 1200)),
                                  column(2, numericInput("height_png_Pie_Epitope", "Height of PNG", value = 1000)),
                                  column(2, numericInput("resolution_PNG_Pie_Epitope", "Resolution of PNG", value = 144)),
                                  column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_Pie_Epitope", "Download PNG"))
                                ),
                              ),
                              tabPanel("Stats",
                                       value = "EpiPan_stat",
                                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                       div(DT::DTOutput("Epi_of_interest_DF")),
                                       div(DT::DTOutput("compare.stat_Epi_DT")),
                                       downloadButton("downloadtb_compare.stat_Epi", "Download table")
                              ),
                              tabPanel("Dotplot",
                                       value = "EpiPan_dot",
                                       fluidRow(
                                         column(2, colourInput("low.dotplot.epi", "Lower color:", "darkblue")),
                                         column(2, colourInput("middle.dotplot.epi", "Middle color:", "white")),
                                         column(2, colourInput("high.dotplot.epi", "High color:", "darkred")),
                                         column(2, checkboxInput("restrict.dotpot.epi", "Restrict to top list", value = F)),
                                         column(2, numericInput("restrict.dotpot.num.epi", "Total genes to display:", value = 10))
                                       ),
                                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                       plotOutput("all_expression_dotplot_epi", height = "400px"),
                                       fluidRow(
                                         column(1, numericInput("width_all_expression_dotplot_epi", "Width of PDF", value = 20)),
                                         column(1, numericInput("height_all_expression_dotplot_epi", "Height of PDF", value = 4)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_all_expression_dotplot_epi", "Download PDF")),
                                         column(2, numericInput("width_png_all_expression_dotplot_epi", "Width of PNG", value = 2400)),
                                         column(2, numericInput("height_png_all_expression_dotplot_epi", "Height of PNG", value = 700)),
                                         column(2, numericInput("resolution_PNG_all_expression_dotplot_epi", "Resolution of PNG", value = 144)),
                                         column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_all_expression_dotplot_epi", "Download PNG"))
                                       ),
                              ),
                              tabPanel("Over-representation",
                                       value = "EpiPan_OvRep",
                                       fluidRow(
                                         column(3, numericInput("in.geneset.cutoff_Epi", "Min number of genes in GeneSet", value = 1, min = 0, step = 1, max = 60000)),
                                         column(3, numericInput("p.val_cutoff_Epi", "p-val cut-off", value = 0.05, min = 0, max = 1)),
                                       ),
                                       div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                       div(DT::DTOutput("Over_rep_Epi_Tab")),
                                       downloadButton("downloadtb_over.rep.Epi", "Download table")
                                       #
                              )
                            ),
                   ),

                   # Overlap ------

                   # tabPanel("Clonotypes per cluster (Pie/bar plot)"),
                   # tabPanel("Upset plot")
                 ),
        ),
        ########
        # prioritization strategy ------

        tabPanel("Automation (TCR -> GEX)",
                 value = "Prior",
                 h5("under construction"),
                 tabsetPanel(
                   id = "PriorTBMods",
                   tabPanel("Analysis steps",
                            value = "ModstoConsid",
                            p("The prioritization section automates the analysis based on certain thresholds"),
                            p("There are three section"),
                            p("1. Clonotype"),
                            p("    - Single sample: immunodominant based on percentage threshold"),
                            p("2. Cluster"),
                            p("3. Epitope"),
                            p("Note: minimum threshold of expansion (n=3) per analysis. Cannot perform stats on n=1 or n=2"),
                   ),
                   # modules of priority Top clone ------
                   tabPanel("Clonotype",
                            value = "PriorRepertoireTB",
                            fluidRow(
                              column(3, checkboxInput("Download_public_overlapping", "Download Public-like", value = T)),
                              column(3, checkboxInput("Download_public_overlapping_bar", "Download Public-like Over-rep", value = F)),
                              column(3, checkboxInput("Download_private_overlapping", "Download Private clone analysis", value = F)),
                              column(3, checkboxInput("restrict_to_expanded", "Restrict to expanded", value = F))
                            ),
                            add_busy_spinner(spin = "fading-circle", position = "top-right",  height = "200px", width = "200px", color = "purple"),
                            verbatimTextOutput("Simple_workflow_step1"),
                            verbatimTextOutput("Number_of_clonotypes_to_"),
                            add_busy_spinner(spin = "fading-circle", position = "top-right",  height = "200px", width = "200px", color = "purple"),
                            div(DT::DTOutput("PriorClono_Tab")),
                            uiOutput("Module_case_statements"),
                            div(DT::DTOutput("private_clonotypes")),
                   ),
                   # modules of priority cluster ------
                   tabPanel("Cluster",
                            value = "PriorClustTB",
                            fluidRow(
                              column(3, uiOutput("Default_priority_cutoffAG")),
                              column(3, uiOutput("Default_priority_cutoffBD")),
                              # column(3,numericInput("CloneTotal_input","Select clones > ",value = 2, min = 2)),
                            ),
                            fluidRow(
                              column(3, numericInput("Sample_count_cluster", "AG Sample count >=", value = 2)),
                              column(3, numericInput("Total_cloneCount_cluster", "AG Clone count >=", value = 3, min = 3)),
                              column(3, numericInput("Sample_count_clusterBD", "BD Sample count >=", value = 3, min = 3)),
                              column(3, numericInput("Total_cloneCount_clusterBD", "BD Clone count >=", value = 2)),
                            ),
                            add_busy_spinner(spin = "fading-circle", position = "top-right",  height = "200px", width = "200px", color = "purple"),
                            verbatimTextOutput("number_clusters_to_analyse_AG"),
                            add_busy_spinner(spin = "fading-circle", position = "top-right",  height = "200px", width = "200px", color = "purple"),
                            verbatimTextOutput("number_clusters_to_analyse_BD"),
                            add_busy_spinner(spin = "fading-circle", position = "top-right",  height = "200px", width = "200px", color = "purple"),
                            uiOutput("Cluster_dowload_button_prior"),
                            # div(DT::DTOutput("colors.top_dt")),
                   ),
                   # modules of priority Epitope ------
                   tabPanel("Epitope/Annotation",
                            value = "PriorEpiTB",
                            add_busy_spinner(spin = "fading-circle", position = "top-right",  height = "200px", width = "200px", color = "purple"),
                            checkboxInput("epitope_uploaded", "Add in Epitope data", value = T),
                            fluidRow(
                              column(3, uiOutput("AddInEpiUI_1")),
                              column(3, uiOutput("AddInEpiUI_2")),
                            ),
                            actionButton("EpitopePrior_Download_Bt", "Automated Epitope analysis"),
                            div(DT::DTOutput("Test_table_1")),
                   ),

                   ### end priority
                 )
        ),

        # GEX -> TCR --------
        tabPanel(
          "GEX -> TCR",
          # h5("Under Development..."),
          tabsetPanel(
            ##### annotation GEX -> TCR --------
            tabPanel("Annotation",
                     value = "PanelAnno_GEXTCR",
                     h5("Under Active Development... Only Table is complete"),
                     fluidRow(
                       column(3, uiOutput("AddInAnnoUI_man_1")),
                       column(3, uiOutput("AddInAnnoUI_man_2")),
                       column(3, selectInput("Slected_annotation", "Select annotation", choices = ""))
                     ),
                     tabsetPanel(
                       tabPanel(
                         "Table",
                         div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                         div(DT::DTOutput("AnnoTable_perMarkers")),
                         downloadButton("downloaddf_AnnoTable_perMarkers", "Download table"),
                       ),
                       tabPanel("UMAP"),
                       tabPanel("TCR per Anno ID"),
                       tabPanel("Stats"),
                       tabPanel("Dot plot"),
                       tabPanel("Over representation")
                     )
            ),

            # marker specific TCR analysis --------
            tabPanel("Marker",
                     value = "Marker",
                     tabsetPanel(
                       id = "Marker_Panel",
                       tabPanel(
                         "Single marker",
                         actionButton("load_marker_genes_sig", "Load genes"),
                         fluidRow(
                           column(4, selectizeInput("Var_to_col_marker", "Marker col", "")),
                           column(2, numericInput("Filter_lower_UMAP1_marker", "UMAP_1 >", value = -20)),
                           column(2, numericInput("Filter_lower_UMAP1_marker2", "UMAP_1 <", value = 20)),
                           column(2, numericInput("Filter_lower_UMAP2_marker", "UMAP_2 >", value = -20)),
                           column(2, numericInput("Filter_lower_UMAP2_marker2", "UMAP_2 <", value = 20)),
                         ),
                         # conditionalPanel(condition="input.Marker_Panel =='Marker_Panel_plot_stats'",
                         fluidRow(
                           column(3, selectInput("Analysis_marker_stats_type", "Analysis type", choices = c("Population", "Expanded"))),
                           column(3, selectInput("Help_marer", "Definitions", choices = c("No", "Yes"))),
                           column(3, numericInput("cutoff_marker_gt", "Marker +ve cut-off (>)", value = 0, step = 0.1)),
                           column(3, numericInput("pos_expanded_cut_off", "clone count (>) for +ve", value = 2, min = 1)),
                         ),
                         conditionalPanel(
                           condition = "input.Help_marer == 'Yes'",
                           # p(strong("1. Total."), "Refers to the marker of interest compared to all negative cells"),
                           p(strong("1. Population."), "Refers to the the restricted UMAP selected population and compared the +ve to -ve marker population."),
                           p(strong("2. Expanded."), "Compares the clonal expaned (min 3) vs non-expanded population for all cells +ve for the marker in the restricted UMAP space."),
                         ),
                         tabsetPanel(
                           tabPanel("Table",
                                    value = "Marker_Panel_table",
                                    div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                    div(DT::DTOutput("marker_selected_tab")),
                           ),
                           tabPanel("UMAP plot",
                                    value = "Marker_Panel_plot_UMAP",
                                    h4("Select area of the plot to keep for the specific marker"),
                                    h6("Recommended to filter to broad populations based on UMAP e.g., CD4, CD8 or other"),
                                    div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                    numericInput("max_scale", "MAX scale", value = ""),
                                    plotOutput("marker_selected_UMAP_plot", height = "600px"),
                                    fluidRow(
                                      column(3, numericInput("width_marker_selected_UMAP_plot", "Width of PDF", value = 10)),
                                      column(3, numericInput("height_marker_selected_UMAP_plot", "Height of PDF", value = 8)),
                                      column(3),
                                      column(3, style = "margin-top: 25px;", downloadButton("downloadPlot_marker_selected_UMAP_plot", "Download PDF"))
                                    ),
                                    fluidRow(
                                      column(3, numericInput("width_png_marker_selected_UMAP_plot", "Width of PNG", value = 1200)),
                                      column(3, numericInput("height_png_marker_selected_UMAP_plot", "Height of PNG", value = 1000)),
                                      column(3, numericInput("resolution_PNG_marker_selected_UMAP_plot", "Resolution of PNG", value = 144)),
                                      column(3, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_marker_selected_UMAP_plot", "Download PNG"))
                                    )
                           ),
                           tabPanel("Violin/Ridge plot",
                                    value = "Marker_Panel_plot_VR",
                                    h4("Filter marker of interest based on threshold"),
                                    fluidRow(
                                      column(3, selectInput("select_plot_vio.ridge", "Plot type", choices = c("Violin", "Ridge")))
                                    ),
                                    div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                    plotOutput("marker_selected_VioRidge_plot", height = "600px")
                           ),
                           # tabPanel("UMAP"),
                           tabPanel("TCR/BCR mapped",
                                    value = "Marker_Panel_plot_TCR",
                                    h4("TCR and/or BCR seqeunces that are positive for that marker"),
                                    div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                    fluidRow(
                                      column(6, div(DT::DTOutput("TCR_marker_positive_count"))),
                                      column(6, div(DT::DTOutput("TCR_marker_neg_count")))
                                    ),
                                    div(DT::DTOutput("merged_marker_hist_table")),
                                    downloadButton("downloaddf_clonotype_distribution", "Download table"),
                                    div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                    plotOutput("marker_selected_histogram_plot", height = "600px")
                           ),
                           tabPanel("Stats",
                                    value = "MP_plot_stats",
                                    div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                                    div(DT::DTOutput("Compare.stat_marker")),
                                    downloadButton("downloaddf_Marker_stats", "Download table")

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
                       tabPanel(
                         "Dual marker analysis",
                         actionButton("load_marker_genes_dual", "Load genes"),
                         # marker 1 # cut-off 1
                         # marker 2 @ cut-off 2 (negative control)
                         fluidRow(
                           column(4, selectizeInput("Var_to_col_marker2", "X-axis Marker", "")),
                           column(4, selectizeInput("Var_to_col_marker3", "Y-axis Marker", ""))
                         ),
                         fluidRow(
                           column(2, numericInput("Filter_dual_UMAP1_marker", "UMAP_1 >", value = -20)),
                           column(2, numericInput("Filter_dual_UMAP1_marker2", "UMAP_1 <", value = 20)),
                           column(2, numericInput("Filter_dual_UMAP2_marker", "UMAP_2 >", value = -20)),
                           column(2, numericInput("Filter_dual_UMAP2_marker2", "UMAP_2 <", value = 20)),
                           # column(2),
                           column(2, numericInput("X_axis_dot_dual", "X-axis line", value = 0, step = 0.1)),
                           column(2, numericInput("Y_axis_dot_dual", "Y-axis line", value = 0, step = 0.1)),
                         ),
                         tabsetPanel(
                           tabPanel(
                             "table",
                             div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                             div(DT::DTOutput("meta_data_for_features_scale2_df")),
                           ),
                           tabPanel(
                             "Feature plots",
                             div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                             fluidRow(
                               column(3, numericInput("max_scale2", "MAX scale (left)", value = "")),
                               column(3, numericInput("max_scale3", "MAX scale (right)", value = "")),
                             ),
                             div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                             fluidRow(
                               column(6, plotOutput("marker_selected_UMAP_plot2", height = "600px")),
                               column(6, plotOutput("marker_selected_UMAP_plot3", height = "600px")),
                             )
                           ),
                           tabPanel(
                             "X by Y plot",
                             div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                             plotOutput("df_dotplot_marker_plot", height = "600px"),
                             fluidRow(
                               column(1, numericInput("width_df_dotplot_marker_plot", "Width of PDF", value = 10)),
                               column(1, numericInput("height_df_dotplot_marker_plot", "Height of PDF", value = 8)),
                               column(2, style = "margin-top: 25px;", downloadButton("downloadPlot_df_dotplot_marker_plot", "Download Network PDF")),
                               column(2, numericInput("width_png_df_dotplot_marker_plot", "Width of PNG", value = 1200)),
                               column(2, numericInput("height_png_df_dotplot_marker_plot", "Height of PNG", value = 1000)),
                               column(2, numericInput("resolution_PNG_df_dotplot_marker_plot", "Resolution of PNG", value = 144)),
                               column(2, style = "margin-top: 25px;", downloadButton("downloadPlotPNG_df_dotplot_marker_plot", "Download Network PNG")),
                             ),
                           ),
                           # tabPanel("Violin plot"),# will include splitting by T cell/B cell markers? Ig?
                           # tabPanel("UMAP"),
                           tabPanel(
                             "TCR table",
                             div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                             div(DT::DTOutput("dual_maker_TCR_Sum_DT")),
                             downloadButton("Dule_marker_TCRsummary_DT", "Download table")
                           ),
                           tabPanel(
                             "Stats",
                             div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
                             selectInput("quad_dualmarker", "Quadrant to compute", choices = c("Q1", "Q2", "Q3", "Q4")),
                             div(DT::DTOutput("dual_maker_TCR_statsTB")),
                             downloadButton("Dule_marker_statsTBDownload", "Download table")
                           ),
                         )
                       )
                     )
            ),

          ),

        ),

        ######
        # tabPanel("Parameters Table",
        #          div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
        #          div(DT::DTOutput("parameterTable")),
        #          # Download button
        #          downloadButton("downloadData", "Download State"),
        #
        #
        #
        #          ),
      )


      ######
    )
  )
),
# end of Integration -----

navbarMenu(
  "Post Analysis",

  # Extracting meta data -------
  tabPanel("Extract meta data",
           sidebarLayout(
             sidebarPanel(width = 3,
                          fileInput("file1_extract.metadata",
                                    "Upload .rds file",
                                    multiple = F,
                                    accept = c(".rds", "rds")
                          ),

                          downloadButton("download_extracting_md_from_sc_tb","download meta data")

             ),
             mainPanel(
               width = 9,
               div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
               div(DT::DTOutput("extracting_md_from_sc_tb")),
             )
           )
  ),

  # filtering CSV files ------
  tabPanel("Filtering",
           value = "",
           sidebarLayout(
             sidebarPanel(width = 3,
                          fileInput("file1_FilteringCluster",
                                    "Upload .csv file",
                                    multiple = F,
                                    accept = c(".csv", "csv","csv.gz",".gz")
                          ),

                          downloadButton("download_button_filtered_data_df", "Download Filtered Data"),
                          uiOutput("filter_checkboxes"),
             ),
             mainPanel(
               width = 9,
               # tabsetPanel(
               id = "Other_post_analysis",
               div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
               uiOutput("filter_inputs"),
               DTOutput("filtered_table")

               # )
             )
           )
           # p("Convert .h5Seurat to .rds")
  ),
  # summarising the data ------

  tabPanel("Summarizing",
           value = "",
           sidebarLayout(
             sidebarPanel(width = 3,
                          fileInput("file1_summarising_data",
                                    "Upload .csv file",
                                    multiple = F,
                                    accept = c(".csv", "csv","csv.gz",".gz")
                          ),
                          downloadButton("download_button_summarising_data", "Download summarised Data"),
             ),
             mainPanel(
               width = 9,
               # tabsetPanel(
               id = "summary_for_data_tab",
               div(id = "spinner-container",class = "centered-spinner",add_busy_spinner(spin = "fading-circle",height = "200px",width = "200px",color = "#6F00B0")),
               selectInput("columns_for_summary","Select columns for summarising",choices = "",multiple = T, width = "1200px"),
               DTOutput("chain_table_summary_tb"),
               # uiOutput("filter_inputs"),
               DTOutput("loaded_for_summary")

               # )
             )
           )
           # p("Convert .h5Seurat to .rds")
  ),



  # post analysis OLGA ------
  tabPanel("OLGA",
           fluidRow(

             add_busy_spinner(spin = "fading-circle", position = "top-right",  height = "200px", width = "200px", color = "#6F00B0"),
             column(3,actionButton("load_olga","Use OLGA if present")),
             conditionalPanel(
               condition = "input.load_olga",
               add_busy_spinner(spin = "fading-circle", position = "top-right",  height = "200px", width = "200px", color = "#6F00B0"),
               column(3,selectInput("Olga_installed","OLGA installed","")),

             )),

           conditionalPanel(
             condition = "input.Olga_installed == 'installed'",

             sidebarLayout(
               sidebarPanel(width = 3,

                            fileInput("file1_rds_OLGA",
                                      "Choose .rds files from merging",
                                      multiple = FALSE,
                                      accept = c("rds", ".rds")
                            ),
                            selectInput("chain_type_olga","Chain type: ",choices = c("TRB","TRA")),
                            downloadButton("downloaddf_pgen_dt", "Download Table (.csv)")
               ),

               mainPanel(
                 width = 9,
                 tabsetPanel(
                   id = "OLGA_analysis",
                   tabPanel("Loaded Table",
                            add_busy_spinner(spin = "fading-circle", position = "top-right",  height = "200px", width = "200px", color = "#6F00B0"),
                            div(DT::DTOutput("Pgen_Selected")),
                   ),
                   tabPanel("Pgen",
                            add_busy_spinner(spin = "fading-circle", position = "top-right",  height = "200px", width = "200px", color = "#6F00B0"),
                            div(DT::DTOutput("Pgen_BD")),
                   )
                 )
               )
             )
           ),
           conditionalPanel(
             condition = "input.Olga_installed == 'Unavailable' && input.load_olga",
             tabPanel("Instructions",

                      p("Require the user to install python and olga from the command line"),
                      p("Install OLGA via the following website"),
                      p("https://github.com/statbiophys/OLGA")

             ),
           )

  ),

  #####

),

navbarMenu("Info",
           tabPanel("Helpful resources"),
           tabPanel("Citing STEGO.R"),
           tabPanel("Session Information")



),
) # nav page
  )

  ########
  # server ------
  server <- function(input, output, session) {
    # convert ------
    Convert_to_RDS <- reactive({
      inFile_sc_pro2 <- input$file1_h5Seurat.file
      if (is.null(inFile_sc_pro2)) {
        return(NULL)
      } else {
        dataframe <- LoadH5Seurat(inFile_sc_pro2$datapath)
      }
    })

    output$Convert_to_RDS_out <- renderPrint({
      sc <- input$file1_h5Seurat.file
      validate(
        need(
          nrow(sc) > 0,
          "Upload files"
        )
      )
      df <- Convert_to_RDS()
      print(df)
    })

    output$downloaddf_SeruatObj_Convert_to_RDS <- downloadHandler(
      filename = function() {
        x <- today()
        if (input$project_name5 == "") {
          showModal(modalDialog(
            title = "Error",
            "Please Add File Name.",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {


          # paste(input$project_name,"_SC.obj_",x,".h5Seurat", sep = "")
          paste(input$project_name5, "_V4_Seurat", x, ".rds", sep = "")
        }
      },
      content = function(file) {
        SaveSeuratRds(Convert_to_RDS(), file)
        # SaveH5Seurat(vals_meta.sc$metadata_SCobj,file)
      }
    )

    # add UI ------
    # observeEvent(input$reset_input, {
    #   shinyjs::reset("side-panel")
    # })

    output$scGate_cutoffs <- renderUI({
      if (input$Data_types == "10x_HS" || input$Data_types == "BD_HS.Full.Panel" || input$Data_types == "10x_MM" || input$Data_types == "BD_MM_Full.Panel") {
        numericInput("threshold_scGate", "scGate threshold", value = 0.2)
      } else {
        numericInput("threshold_scGate", "scGate threshold", value = 0.5)
      }
      # c("10x_HS","BD_HS.Immune.Panel","BD_HS.Full.Panel","10x_MM","BD_MM_Full.Panel","BD_MM_Immune.Panel"
    })



    output$classification_to_add <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      df3.meta <- c(names(sc@meta.data))
      # df3.meta <- df3.meta[!grepl("RNA",df3.meta) & !grepl("BCR",df3.meta) & !grepl("TCR",df3.meta)& !grepl("_gene",df3.meta) & !grepl("allele",df3.meta) & !grepl("percent",df3.meta) & !grepl("cdr3",df3.meta)]

      # tags$style("#Colour_By_this {background-color:#e5e5e5;}")

      selectInput("Colour_By_this", "Colour by: ", choices = df3.meta, selected = "Sample_Name")
    })
    ## overview colouring
    output$classification_to_add_overview <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      df3.meta <- c(names(sc@meta.data))
      # df3.meta <- df3.meta[!grepl("RNA",df3.meta) & !grepl("BCR",df3.meta) & !grepl("TCR",df3.meta)& !grepl("_gene",df3.meta) & !grepl("allele",df3.meta) & !grepl("percent",df3.meta) & !grepl("cdr3",df3.meta)]
      # tags$style("#Colour_By_this_overview {background-color:#e5e5e5;}")
      selectInput("Colour_By_this_overview", "Colour by: ", choices = df3.meta, selected = "seurat_clusters")
    })
    output$classification_to_add_epitope <- renderUI({
      sc <- UMAP_metadata_with_labs()

      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )
      df3.meta <- sc@meta.data
      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta, epi, by = "CDR3_beta")
      selectInput("epitope_umap_selected", "Colour Pie by (hm = y-axis):", choices = names(df3.meta), selected = "beta")
    })

    output$classification_to_add_epitope2 <- renderUI({
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )
      df3.meta <- sc@meta.data
      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta, epi, by = "CDR3_beta")

      selectInput("epitope_umap_selected2", "Split Pie by (hm = x-axis):", choices = names(df3.meta), selected = "epitope")
    })

    # UI for Step 3a.

    # user interface parameters-----
    output$feature_input <- renderUI({
      if (input$df_seruatobj_type == "10x_Genomics (raw)") {
        fluidRow(
          column(6, numericInput("features.min", "minimum features (>)", value = 200)),
          column(6, numericInput("features.max", "Maximum features (<)", value = 6000)),
          column(6, numericInput("percent.mt", "Mitochondrial DNA cut-off (<)", value = 20)),
          column(6, numericInput("percent.rb", "Ribosomal RNA cut-off (>)", value = 5)),
        )
      } else if (input$df_seruatobj_type == "10x_Genomics (.h5)") {
        fluidRow(
          column(6, numericInput("features.min", "minimum features (>)", value = 200)),
          column(6, numericInput("features.max", "Maximum features (<)", value = 6000)),
          column(6, numericInput("percent.mt", "Mitochondrial DNA cut-off (<)", value = 20)),
          column(6, numericInput("percent.rb", "Ribosomal RNA cut-off (>)", value = 5)),
        )
      } else if (input$df_seruatobj_type == "BD Rhapsody (Mouse)") {
        fluidRow(
          column(6, numericInput("features.min", "minimum features (>)", value = 200)),
          column(6, numericInput("features.max", "Maximum features (<)", value = 6000)),
          column(6, numericInput("percent.mt", "Mitochondrial DNA cut-off (<)", value = 20)),
          column(6, numericInput("percent.rb", "Ribosomal RNA cut-off (>)", value = 0)),
        )
      } else if (input$df_seruatobj_type == "BD Rhapsody (Full panel)") {
        fluidRow(
          column(6, numericInput("features.min", "minimum features (>)", value = 200)),
          column(6, numericInput("features.max", "Maximum features (<)", value = 6000)),
          column(6, numericInput("percent.mt", "Mitochondrial DNA cut-off (<)", value = 20)),
          column(6, numericInput("percent.rb", "Ribosomal RNA cut-off (>)", value = 5)),
        )
      } else if (input$df_seruatobj_type == "BD Rhapsody (Human Immune panel)") {
        fluidRow(
          column(6, numericInput("features.min", "minimum features (>)", value = 45)),
          column(6, numericInput("features.max", "Maximum features (<)", value = 160)),
          column(6, numericInput("percent.mt", "Mitochondrial DNA cut-off (<)", value = 0)),
          column(6, numericInput("percent.rb", "Ribosomal RNA cut-off (>)", value = 0)),
        )
      } else {
        fluidRow(
          column(6, numericInput("features.min", "minimum features (>)", value = 200)),
          column(6, numericInput("features.max", "Maximum features (<)", value = 6000)),
          column(6, numericInput("percent.mt", "Mitochondrial DNA cut-off (<)", value = 20)),
          column(6, numericInput("percent.rb", "Ribosomal RNA cut-off (>)", value = 0)),
        )
      }
    })
    # three files required for BD data: Sample Tag calls, TCR file and count ----
    input.data.calls.bd <- reactive({
      inFile12 <- input$file_calls_BD
      if (is.null(inFile12)) {
        return(NULL)
      } else {
        dataframe <- read.csv(
          inFile12$datapath,
          skip = input$no_lines_skip_Tags, header = T
        )
      }
    })

    output$test.files <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      calls <- input.data.calls.bd()
      validate(
        need(
          nrow(calls) > 0,
          error_message_val1
        )
      )
      calls
    })

    input.data.TCR.BD <- reactive({
      inFileTCRBD <- input$file_TCR_BD
      if (is.null(inFileTCRBD)) {
        return(NULL)
      } else {
        dataframe <- read.csv(
          inFileTCRBD$datapath,
          skip = input$no_lines_skip_TCR, header = T
        )
      }
    })
    output$test.files2 <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      calls <- input.data.TCR.BD()
      validate(
        need(
          nrow(calls) > 0,
          error_message_val2
        )
      )
      calls
    })
    input.data.count.BD <- reactive({
      inFilecountBD <- input$file_counts_BD
      if (is.null(inFilecountBD)) {
        return(NULL)
      } else {
        dataframe <- read.csv(inFilecountBD$datapath, skip = input$no_lines_skip_counts, header = T)
      }
    })
    output$test.files3 <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      calls <- as.data.frame(input.data.count.BD())
      validate(
        need(
          nrow(calls) > 0,
          error_message_val3
        )
      )
      calls[1:6, 1:6]
    })

    # 10x format is now in BD rhapsody ----
    # barcode file -----
    input.data.barcode.bd <- reactive({
      inFile_bd_barcode <- input$file_barcode_bd
      if (is.null(inFile_bd_barcode)) {
        return(NULL)
      } else {
        dataframe <- read.table(
          inFile_bd_barcode$datapath
        )
      }
    })

    output$test.files.bd1 <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      calls <- input.data.barcode.bd()
      validate(
        need(
          nrow(calls) > 0,
          "Upload file"
        )
      )
      calls
    })

    # features file -----
    input.data.features.bd2 <- reactive({
      inFile_bd_features <- input$file_features_bd
      if (is.null(inFile_bd_features)) {
        return(NULL)
      } else {
        dataframe <- read.table(
          inFile_bd_features$datapath
        )
      }
    })

    output$test.files.bd2 <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      calls <- input.data.features.bd2()
      validate(
        need(
          nrow(calls) > 0,
          "upload file"
        )
      )
      calls
    })

    # Matrix file ----
    input.data.matrix.bd <- reactive({
      inFile_bd_matrix <- input$file_matrix_bd
      if (is.null(inFile_bd_matrix)) {
        return(NULL)
      } else {
        dataframe <- Matrix::readMM(inFile_bd_matrix$datapath)
      }
    })
    output$test.files.bd3 <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      calls <- as.data.frame(input.data.matrix.bd())
      validate(
        need(
          nrow(calls) > 0,
          "Upload Matrix"
        )
      )
      calls[1:10, 1:10]
    })



    # contig files ----
    input.data.TCR.bd2 <- reactive({
      inFile_bd2_TCR <- input$file_TCR_bd2
      if (is.null(inFile_bd2_TCR)) {
        return(NULL)
      } else {
        dataframe <- read.table(inFile_bd2_TCR$datapath, sep = "\t", header = T)
      }
    })
    output$test.files.bd4 <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      calls <- input.data.TCR.bd2()
      validate(
        need(
          nrow(calls) > 0,
          "Upload AIRR Contigs (Dominant)"
        )
      )
      calls
    })

    ## BD combining Cell x genes ----
    TCR_Filtering_Paired <- function() {
      TCR <- as.data.frame(input.data.TCR.BD())
      sample_tags <- input.data.calls.bd()

      validate(
        need(
          nrow(TCR) > 0 && nrow(sample_tags) > 0,
          "Upload files"
        )
      )

      calls_TCR <- merge(sample_tags, TCR, by = "Cell_Index")
      calls_TCR_count <- calls_TCR

      # filtering to paired TCR
      calls_TCR_count$TRAG <- ifelse(grepl("TR", calls_TCR_count$TCR_Alpha_Gamma_V_gene_Dominant), "TRAG", "Missing TRAG gene")
      calls_TCR_count$TRBD <- ifelse(grepl("TR", calls_TCR_count$TCR_Beta_Delta_V_gene_Dominant), "TRBD", "Missing TRBD gene")
      calls_TCR_count$TRAG_fun <- ifelse(grepl("[*]", calls_TCR_count$TCR_Alpha_Gamma_CDR3_Translation_Dominant), "Non-functional",
                                         ifelse(grepl("Missing", calls_TCR_count$TRAG), "Missing TCR",
                                                ifelse(grepl("TRGV10", calls_TCR_count$TCR_Alpha_Gamma_V_gene_Dominant), "pseudogene", "productive")
                                         )
      )

      calls_TCR_count$TRBD_fun <- ifelse(grepl("[*]", calls_TCR_count$TCR_Beta_Delta_CDR3_Translation_Dominant), "Non-functional",
                                         ifelse(grepl("Missing", calls_TCR_count$TRBD), "Missing TCR", "productive")
      )


      calls_TCR_count$paired_TCR <- ifelse(calls_TCR_count$TRAG == "TRAG" & calls_TCR_count$TRBD == "TRBD", "paired TCR",
                                           ifelse(calls_TCR_count$TRAG == "TRAG" | calls_TCR_count$TRBD == "TRBD", "Unpaired TCR", "No TCR")
      )
      calls_TCR_count$Productive_TCR <- ifelse(calls_TCR_count$TRBD_fun == "productive" & calls_TCR_count$TRAG_fun == "productive", "productive TCR",
                                               ifelse(calls_TCR_count$TRBD_fun == "productive" & calls_TCR_count$TRAG_fun == "Missing TCR", "productive TCR",
                                                      ifelse(calls_TCR_count$TRBD_fun == "Missing TCR" & calls_TCR_count$TRAG_fun == "productive", "productive TCR",
                                                             "unproductive"
                                                      )
                                               )
      )

      if (input$BCR_present == T) {
        calls_TCR_count$v_gene_IgL <- ifelse(grepl("IG", calls_TCR_count$BCR_Light_V_gene_Dominant), "v_gene_IgL", "Missing IgL gene")
        calls_TCR_count$v_gene_IgH <- ifelse(grepl("IG", calls_TCR_count$BCR_Heavy_V_gene_Dominant), "v_gene_IgH", "Missing IgL gene")
        calls_TCR_count$v_gene_IgL_fun <- ifelse(grepl("[*]", calls_TCR_count$BCR_Light_CDR3_Translation_Dominant), "Non-functional",
                                                 ifelse(grepl("Missing", calls_TCR_count$v_gene_IgL), "Missing BCR", "productive")
        )
        calls_TCR_count$v_gene_IgH_fun <- ifelse(grepl("[*]", calls_TCR_count$BCR_Heavy_CDR3_Translation_Dominant), "Non-functional",
                                                 ifelse(grepl("Missing", calls_TCR_count$v_gene_IgH), "Missing BCR", "productive")
        )
        #
        calls_TCR_count$paired_BCR <- ifelse(calls_TCR_count$v_gene_IgL == "v_gene_IgL" & calls_TCR_count$v_gene_IgH == "v_gene_IgH", "paired BCR",
                                             ifelse(calls_TCR_count$v_gene_IgL == "v_gene_IgL" | calls_TCR_count$v_gene_IgH == "v_gene_IgH", "Unpaired BCR", "No BCR")
        )


        calls_TCR_count$Productive_BCR <- ifelse(calls_TCR_count$v_gene_IgL_fun == "productive" & calls_TCR_count$v_gene_IgH_fun == "productive", "productive BCR",
                                                 ifelse(calls_TCR_count$v_gene_IgL_fun == "productive" & calls_TCR_count$v_gene_IgH_fun == "Missing BCR", "productive BCR",
                                                        ifelse(calls_TCR_count$v_gene_IgL_fun == "Missing BCR" & calls_TCR_count$v_gene_IgH_fun == "productive", "productive BCR",
                                                               "unproductive"
                                                        )
                                                 )
        )

        calls_TCR_count$BCR_TCR <- ifelse(calls_TCR_count$Productive_BCR == "productive BCR" & calls_TCR_count$Productive_TCR == "productive TCR", "Both TCR and BCR present", NA)
      } else {
        calls_TCR_count <- calls_TCR_count
      }

      # filter out non functional TCR
      if (input$BCR_present == T) {
        productive <- calls_TCR_count[calls_TCR_count$Productive_TCR %in% "productive TCR" | calls_TCR_count$Productive_BCR %in% "productive BCR", ]
      } else {
        productive <- calls_TCR_count[calls_TCR_count$Productive_TCR %in% "productive TCR", ]
      }
      # filter out un-paired TCR's
      if (input$filtering_TCR == T) {
        if (input$BCR_present == T) {
          paired <- productive[productive$paired_TCR %in% "paired TCR" | productive$paired_BCR %in% "paired BCR", ]
        } else {
          paired <- productive[productive$paired_TCR %in% "paired TCR", ]
        }
      } else {
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
          choices = names(input.data.TCR.BD()),
          selected = "chain"
        )
      } else {
        updateSelectInput(
          session,
          "locus_column",
          choices = names(input.data.TCR.bd2()),
          selected = "locus"
        )
      }
    })

    filtering_required_TCR <- reactive({
      contigs <- input.data.TCR.bd2()

      validate(
        need(
          nrow(contigs) > 0,
          "Upload files"
        )
      )

      TCR_unfiltered <- contigs
      TCR_unfiltered$seq_issue <- ifelse(grepl("[*]", TCR_unfiltered$sequence_alignment_aa), "stop-codon", "productive")
      TCR_unfiltered_prod <- subset(TCR_unfiltered, TCR_unfiltered$seq_issue == "productive")
      print(names(TCR_unfiltered_prod))
      if(length(TCR_unfiltered_prod$cdr3_length)>0) {
      } else {
        message("adding junction length")
        TCR_unfiltered_prod$cdr3_length <- nchar(TCR_unfiltered_prod$junction_aa)

      }

      TCR_unfiltered_prod <- subset(TCR_unfiltered_prod, TCR_unfiltered_prod$cdr3_length < 30)
      TCR_unfiltered_prod <- subset(TCR_unfiltered_prod, TCR_unfiltered_prod$cdr3_length > 5)
      df_unique <- TCR_unfiltered_prod[, names(TCR_unfiltered_prod) %in% c("locus", "junction_aa")]
      df_unique$Clonal_Expanded <- 1
      df_unique_sum <- ddply(df_unique, names(df_unique)[-c(3)], numcolwise(sum))
      TCR_unfiltered_prod <- merge(TCR_unfiltered_prod, df_unique_sum, by = c("locus", "junction_aa"))
      TCR_unfiltered_prod$juct_issue <- ifelse(grepl("CTCSAR", TCR_unfiltered_prod$junction_aa), "CTCSAR",
                                               ifelse(grepl(".TCSAR", TCR_unfiltered_prod$junction_aa), ".TCSAR",
                                                      ifelse(grepl("^C", TCR_unfiltered_prod$junction_aa) & grepl("F$", TCR_unfiltered_prod$junction_aa), "C___F",
                                                             ifelse(grepl("^C", TCR_unfiltered_prod$junction_aa), "C__",
                                                                    ifelse(TCR_unfiltered_prod$junction_aa == "", "blank",
                                                                           "other"
                                                                    )
                                                             )
                                                      )
                                               )
      )



      TCR_unfiltered_prod_junc <- subset(TCR_unfiltered_prod, TCR_unfiltered_prod$juct_issue != "blank")

      count.chain <- as.data.frame(table(TCR_unfiltered_prod_junc$cell_id))
      names(count.chain) <- c("cell_id", "seq_identified")

      TCR_unfiltered_prod_junc_seq <- merge(TCR_unfiltered_prod_junc, count.chain, by = "cell_id")
      sum_tab <- TCR_unfiltered_prod_junc_seq[, names(TCR_unfiltered_prod_junc_seq) %in% c("locus", "cell_id")]
      sum_tab$count <- 1
      sum_tab_2 <- as.data.frame(ddply(sum_tab, names(sum_tab)[1:2], numcolwise(sum)))
      count.chain <- as.data.frame(table(TCR_unfiltered_prod_junc_seq$cell_id, TCR_unfiltered_prod_junc_seq$locus))
      names(count.chain) <- c("cell_id", "chain", "seq_identified")
      TRA <- subset(count.chain, count.chain$chain == "TRA")
      names(TRA)[2:3] <- paste0(names(TRA), "_TRA")[2:3]
      TRB <- subset(count.chain, count.chain$chain == "TRB")
      names(TRB)[2:3] <- paste0(names(TRB), "_TRB")[2:3]
      TRG <- subset(count.chain, count.chain$chain == "TRG")
      names(TRG)[2:3] <- paste0(names(TRG), "_TRG")[2:3]
      TRD <- subset(count.chain, count.chain$chain == "TRD")
      names(TRD)[2:3] <- paste0(names(TRD), "_TRD")[2:3]

      TRAB <- merge(TRA, TRB, by = "cell_id")
      TRGD <- merge(TRG, TRD, by = "cell_id")
      TRAB_TRGD <- merge(TRAB, TRGD, by = "cell_id")
      names(TRAB_TRGD)
      head(TRAB_TRGD)
      TRAB_TRGD_2 <- TRAB_TRGD %>%
        mutate(Clonality = case_when(
          seq_identified_TRA > 0 & seq_identified_TRB > 0 & seq_identified_TRG > 0 & seq_identified_TRD > 0 ~ "AB GD", # all
          # three chains called
          seq_identified_TRA > 0 & seq_identified_TRB > 0 & seq_identified_TRG > 0 & seq_identified_TRD == 0 ~ "AB G",
          seq_identified_TRA > 0 & seq_identified_TRB > 0 & seq_identified_TRG == 0 & seq_identified_TRD > 0 ~ "AB D",
          seq_identified_TRA > 0 & seq_identified_TRB == 0 & seq_identified_TRG > 0 & seq_identified_TRD > 0 ~ "A GD",
          seq_identified_TRA == 0 & seq_identified_TRB > 0 & seq_identified_TRG > 0 & seq_identified_TRD > 0 ~ "B GD",
          # two chains called
          seq_identified_TRA > 0 & seq_identified_TRB > 0 & seq_identified_TRG == 0 & seq_identified_TRD == 0 ~ "AB",
          seq_identified_TRA == 0 & seq_identified_TRB == 0 & seq_identified_TRG > 0 & seq_identified_TRD > 0 ~ "GD",
          seq_identified_TRA == 0 & seq_identified_TRB > 0 & seq_identified_TRG > 0 & seq_identified_TRD == 0 ~ "B G",
          seq_identified_TRA == 0 & seq_identified_TRB > 0 & seq_identified_TRG == 0 & seq_identified_TRD > 0 ~ "B D",
          seq_identified_TRA > 0 & seq_identified_TRB == 0 & seq_identified_TRG == 0 & seq_identified_TRD > 0 ~ "A D",
          seq_identified_TRA > 0 & seq_identified_TRB == 0 & seq_identified_TRG > 0 & seq_identified_TRD == 0 ~ "A G",
          # one chains called
          seq_identified_TRA > 0 & seq_identified_TRB == 0 & seq_identified_TRG == 0 & seq_identified_TRD == 0 ~ "A",
          seq_identified_TRA == 0 & seq_identified_TRB > 0 & seq_identified_TRG == 0 & seq_identified_TRD == 0 ~ "B",
          seq_identified_TRA == 0 & seq_identified_TRB == 0 & seq_identified_TRG > 0 & seq_identified_TRD == 0 ~ "G",
          seq_identified_TRA == 0 & seq_identified_TRB == 0 & seq_identified_TRG == 0 & seq_identified_TRD > 0 ~ "D",
          TRUE ~ "Other"
        ))
      TRAB_TRGD_paired <- TRAB_TRGD_2
      TRAB_TRGD_paired <- TRAB_TRGD_paired[, !grepl("chain_", names(TRAB_TRGD_paired))]
      TRAB_TRGD_paired$sum <- rowSums(TRAB_TRGD_paired[2:5])
      names(TRAB_TRGD_paired)
      head(TRAB_TRGD_paired)
      TRAB_TRGD_paired$pairing_type <- ifelse(TRAB_TRGD_paired$Clonality == "AB" & TRAB_TRGD_paired$sum == 2, "standard AB",
                                              ifelse(TRAB_TRGD_paired$Clonality == "GD" & TRAB_TRGD_paired$sum == 2, "standard GD",
                                                     ifelse(TRAB_TRGD_paired$sum == 1, "Unpaired",
                                                            "non-standard"
                                                     )
                                              )
      )

      df <- merge(TCR_unfiltered_prod_junc_seq, TRAB_TRGD_paired, by = "cell_id")
      df_2 <- df
      Standard_AB <- subset(df_2, df_2$pairing_type == "standard AB")
      Standard_GD <- subset(df_2, df_2$pairing_type == "standard GD")
      non_standard <- subset(df_2, df_2$pairing_type == "non-standard")
      non_standard <- non_standard[order(non_standard$consensus_count, decreasing = T), ]
      non_standard$filter <- paste(non_standard$cell_id, non_standard$locus)
      non_standard2 <- non_standard[!duplicated(non_standard$filter), ]
      non_standard2 <- non_standard2[order(non_standard2$cell_id, decreasing = F), ]
      unique(non_standard2$Clonality)
      non_standard2$Standard_to_keep <- ifelse(non_standard2$Clonality == "AB G" & non_standard2$locus == "TRG", "remove_G",
                                               ifelse(non_standard2$Clonality == "AB D" & non_standard2$locus == "TRD", "remove_D",
                                                      ifelse(non_standard2$Clonality == "B GD" & non_standard2$locus == "TRB", "remove_B",
                                                             ifelse(non_standard2$Clonality == "A GD" & non_standard2$locus == "TRA", "remove_A",
                                                                    ifelse(non_standard2$Clonality == "AB GD" & c(non_standard2$locus == "TRG" | non_standard2$locus == "TRD"), "remove_GD", "keep")
                                                             )
                                                      )
                                               )
      )

      non_standard <- subset(non_standard2, non_standard2$Standard_to_keep == "keep")
      non_standard <- non_standard[, !names(non_standard) %in% c("filter", "Standard_to_keep")]
      dim(non_standard)


      Filtered_TCR_list <- rbind(Standard_AB, Standard_GD, non_standard)
      Filtered_TCR_list <- Filtered_TCR_list[order(Filtered_TCR_list$cell_id, decreasing = F), ]
      contigs <- Filtered_TCR_list
      contigs
    })

    tb_bd_meta.data_TCR <- function() {
      contigs <- input.data.TCR.bd2()
      sample_tags <- input.data.calls.bd()
      validate(
        need(
          nrow(contigs) > 0 && nrow(sample_tags) > 0,
          "Upload files"
        )
      )

      if (input$filtered_list == "Unfiltered") {
        contigs <- filtering_required_TCR()
      } else {
        contigs <- contigs
      }
      names(contigs)[names(contigs) %in% "cell_id"] <- "Cell_Index"

      contigs_merge <- merge(contigs, sample_tags, by = "Cell_Index")
      # remove non-functional sequences
      if (nrow(contigs_merge[-c(grep("[*]", contigs_merge$junction_aa)), ] > 0)) {
        contigs_merge <- contigs_merge[-c(grep("[*]", contigs_merge$junction_aa)), ]
      }

      # removed undefined sequences

      contigs_merge <- subset(contigs_merge, contigs_merge$Sample_Name != "Multiplet")
      contigs_merge <- subset(contigs_merge, contigs_merge$Sample_Name != "Undetermined")

      if (input$filtered_list == "Unfiltered") {
        contigs_lim <- contigs_merge[!names(contigs_merge) %in% c(
          "consensus_count", "sequence_id", "duplicate_count", "germline_alignment", "reads", "length", "cdr3", "rev_comp", "complete_vdj", names(contigs_merge[grep("fwr", names(contigs_merge))]), "cdr1", "cdr2", "productive", names(contigs_merge[grep("sequence", names(contigs_merge))]), names(contigs_merge[grep("cigar", names(contigs_merge))]), names(contigs_merge[grep("support", names(contigs_merge))]), "Sample_Tag", "sum", "dominant", "putative_cell", "juct_issue", "seq_issue", "c_call"
          # "sum","dominant","putative_cell","seq_issue","juct_issue"
        )]
      } else {
        contigs_lim <- contigs_merge[!names(contigs_merge) %in% c(
          "consensus_count", "sequence_id", "duplicate_count", "germline_alignment", "reads", "length", "cdr3", "rev_comp", "complete_vdj", names(contigs_merge[grep("fwr", names(contigs_merge))]), "cdr1", "cdr2", "productive", names(contigs_merge[grep("sequence", names(contigs_merge))]), names(contigs_merge[grep("cigar", names(contigs_merge))]), names(contigs_merge[grep("support", names(contigs_merge))]), "Sample_Tag"
          # "sum","dominant","putative_cell","seq_issue","juct_issue"
        )]
      }

      contigs_lim$v_gene <- gsub("[*]0.", "", contigs_lim$v_call)
      contigs_lim$j_gene <- gsub("[*]0.", "", contigs_lim$j_call)
      contigs_lim$d_gene <- gsub("[*]0.", "", contigs_lim$d_call)
      contigs_lim

      names(contigs_lim)[names(contigs_lim) %in% input$locus_column] <- "chain"

      contig_AG <- subset(contigs_lim, contigs_lim$chain == "TRA" | contigs_lim$chain == "TRG")

      if (input$filtered_list == "Unfiltered") {
        contig_AG <- contig_AG[!duplicated(contig_AG$Cell_Index), ]
        name.list <- names(contig_AG[c(
          names(contig_AG[grep("gene", names(contig_AG))]),
          names(contig_AG[grep("call", names(contig_AG))]),
          names(contig_AG[grep("cdr3", names(contig_AG))]),
          names(contig_AG[grep("cdr2", names(contig_AG))]),
          names(contig_AG[grep("cdr1", names(contig_AG))]),
          names(contig_AG[grep("junction", names(contig_AG))]),
          "chain",
          "cell_type_experimental",
          "Clonal_Expanded"
        )])
      } else {
        name.list <- names(contig_AG[c(
          names(contig_AG[grep("gene", names(contig_AG))]),
          names(contig_AG[grep("call", names(contig_AG))]),
          names(contig_AG[grep("cdr3", names(contig_AG))]),
          names(contig_AG[grep("cdr2", names(contig_AG))]),
          names(contig_AG[grep("cdr1", names(contig_AG))]),
          names(contig_AG[grep("junction", names(contig_AG))]),
          "chain",
          "cell_type_experimental"
        )])
      }

      contig_AG <- contig_AG %>%
        select(all_of(name.list), everything())
      names(contig_AG)[1:summary(name.list)[1]] <- paste(names(contig_AG[names(contig_AG) %in% name.list]), "_AG", sep = "")
      # contig_AG

      if (input$filtered_list == "Unfiltered") {
        if (length(contigs_lim$Clonality[contigs_lim$Clonality == "B D"]) > 0) {
          contig_D <- subset(contigs_lim, contigs_lim$chain == "TRD" & contigs_lim$Clonality == "B D")
          name.list <- names(contig_D[c(
            names(contig_D[grep("gene", names(contig_D))]),
            names(contig_D[grep("call", names(contig_D))]),
            names(contig_D[grep("cdr3", names(contig_D))]),
            names(contig_D[grep("cdr2", names(contig_D))]),
            names(contig_D[grep("cdr1", names(contig_D))]),
            names(contig_D[grep("junction", names(contig_D))]),
            "chain",
            "cell_type_experimental",
            "Clonal_Expanded"
          )])

          contig_D <- contig_D %>%
            select(all_of(name.list), everything())
          names(contig_D)[1:summary(name.list)[1]] <- paste(names(contig_D[names(contig_D) %in% name.list]), "_AG", sep = "")

          contig_AG <- rbind(contig_AG, contig_D)
        } else {
          contig_AG <- contig_AG
        }
      }

      contig_BD <- subset(contigs_lim, contigs_lim$chain == "TRB" | contigs_lim$chain == "TRD")

      if (input$filtered_list == "Unfiltered") {
        name.list <- names(contig_BD[c(
          names(contig_BD[grep("gene", names(contig_BD))]),
          names(contig_BD[grep("call", names(contig_BD))]),
          names(contig_BD[grep("cdr3", names(contig_BD))]),
          names(contig_BD[grep("cdr2", names(contig_BD))]),
          names(contig_BD[grep("cdr1", names(contig_BD))]),
          names(contig_BD[grep("junction", names(contig_BD))]),
          "chain",
          "cell_type_experimental",
          "Clonal_Expanded"
        )])
      } else {
        name.list <- names(contig_BD[c(
          names(contig_BD[grep("gene", names(contig_BD))]),
          names(contig_BD[grep("call", names(contig_BD))]),
          names(contig_BD[grep("cdr3", names(contig_BD))]),
          names(contig_BD[grep("cdr2", names(contig_BD))]),
          names(contig_BD[grep("cdr1", names(contig_BD))]),
          names(contig_BD[grep("junction", names(contig_BD))]),
          "chain",
          "cell_type_experimental"
        )])
      }

      contig_BD <- contig_BD %>%
        select(all_of(name.list), everything())
      names(contig_BD)[1:summary(name.list)[1]] <- paste(names(contig_BD[names(contig_BD) %in% name.list]), "_BD", sep = "")

      if (input$filtered_list == "Unfiltered") {
        if (length(contigs_lim$Clonality[contigs_lim$Clonality == "B D"]) > 0) {
          contig_B <- subset(contigs_lim, contigs_lim$chain == "TRB" & contigs_lim$Clonality == "B D")

          name.list <- names(contig_B[c(
            names(contig_B[grep("gene", names(contig_B))]),
            names(contig_B[grep("call", names(contig_B))]),
            names(contig_B[grep("cdr3", names(contig_B))]),
            names(contig_B[grep("cdr2", names(contig_B))]),
            names(contig_B[grep("cdr1", names(contig_B))]),
            names(contig_B[grep("junction", names(contig_B))]),
            "chain",
            "cell_type_experimental",
            "Clonal_Expanded"
          )])

          contig_B <- contig_B %>%
            select(all_of(name.list), everything())
          names(contig_B)[1:summary(name.list)[1]] <- paste(names(contig_B[names(contig_B) %in% name.list]), "_BD", sep = "")

          head(contig_B)
          contig_BD$other <- ifelse(contig_BD$Cell_Index %in% contig_B$Cell_Index & contig_BD$chain_BD == "TRD", "duplicate", "other")
          contig_BD <- subset(contig_BD, contig_BD$other == "other")
          contig_BD <- contig_BD[, !names(contig_BD) %in% "other"]
          contig_BD <- rbind(contig_BD, contig_B)
        }
      }

      if (input$filtered_list == "Unfiltered") {
        contig_paired <- merge(contig_AG, contig_BD, by = c("Cell_Index", "Sample_Name", names(contig_BD[grep("seq_", names(contig_BD))]), "Clonality", "pairing_type"), all = T)
      } else {
        contig_paired <- merge(contig_AG, contig_BD, by = c("Cell_Index", "Sample_Name"), all = T)
      }
      print("merged AG and BD")


      if (input$filtered_list == "Unfiltered") {
        contig_paired$pairing <- ifelse(contig_paired$Clonality == "AB", "abTCR Paired",
                                        ifelse(contig_paired$Clonality == "GD", "gdTCR Paired",
                                               ifelse(contig_paired$Clonality == "B D", "bdTCR Paired",
                                                      ifelse(contig_paired$Clonality == "B G", "bgTCR Paired",
                                                             ifelse(contig_paired$Clonality == "A D", "adTCR Paired",
                                                                    ifelse(contig_paired$Clonality == "B", "orphan B",
                                                                           ifelse(contig_paired$Clonality == "A", "orphan A",
                                                                                  ifelse(contig_paired$Clonality == "G", "orphan G",
                                                                                         ifelse(contig_paired$Clonality == "D", "orphan D", NA)
                                                                                  )
                                                                           )
                                                                    )
                                                             )
                                                      )
                                               )
                                        )
        )
        print(table(contig_paired$pairing))
        contig_paired$pairing[is.na(contig_paired$pairing)] <- "OTHER"
        print(table(contig_paired$pairing))
      }

      contig_paired_only <- contig_paired
      contig_paired_only[is.na(contig_paired_only)] <- "None"
      contig_paired_only[contig_paired_only == ""] <- "None"
      contig_paired_only$d_gene_BD[contig_paired_only$d_gene_BD == "None"] <- "_"
      # contig_paired_only$d_gene_AG[contig_paired_only$d_gene_AG == 'None'] <- "_"

      if (input$filtering_TCR == T) {
        contig_paired_only <- subset(contig_paired_only, contig_paired_only$junction_AG != "None")
        contig_paired_only <- subset(contig_paired_only, contig_paired_only$junction_BD != "None")
        print(dim(contig_paired_only))
      }

      contig_paired_only$vj_gene_AG <- paste(contig_paired_only$v_gene_AG, contig_paired_only$j_gene_AG, sep = ".")
      contig_paired_only$vj_gene_AG <- gsub("None.None", "", contig_paired_only$vj_gene_AG)

      contig_paired_only$vj_gene_BD <- paste(contig_paired_only$v_gene_BD, contig_paired_only$j_gene_BD, sep = ".")
      contig_paired_only$vj_gene_BD <- gsub(".NA.", ".", contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub("[.]None[.]", ".", contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub("None.None", "", contig_paired_only$vj_gene_BD)

      contig_paired_only$vdj_gene_BD <- paste(contig_paired_only$v_gene_BD, contig_paired_only$d_gene_BD, contig_paired_only$j_gene_BD, sep = ".")
      contig_paired_only$vdj_gene_BD <- gsub("[.]_[.]", ".", contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("None.None", ".", contig_paired_only$vdj_gene_BD)
      #
      contig_paired_only$vj_gene_cdr3_AG <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$junction_aa_AG, sep = "_")
      contig_paired_only$vj_gene_cdr3_AG <- gsub("_None", "", contig_paired_only$vj_gene_cdr3_AG)
      #
      contig_paired_only$vj_gene_cdr3_BD <- paste(contig_paired_only$vj_gene_BD, contig_paired_only$junction_aa_BD, sep = "_")
      contig_paired_only$vj_gene_cdr3_BD <- gsub("_None", "", contig_paired_only$vj_gene_cdr3_BD)
      #
      contig_paired_only$vdj_gene_cdr3_BD <- paste(contig_paired_only$vdj_gene_BD, contig_paired_only$junction_aa_BD, sep = "_")
      contig_paired_only$vdj_gene_cdr3_BD <- gsub("_None", "", contig_paired_only$vdj_gene_cdr3_BD)
      #
      contig_paired_only$vj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$vj_gene_BD, sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$vdj_gene_BD, sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- gsub("^ & ", "", contig_paired_only$vdj_gene_AG_BD)
      contig_paired_only$vdj_gene_AG_BD <- gsub(" & $", "", contig_paired_only$vdj_gene_AG_BD)
      #
      # #updating names to be consistant....
      contig_paired_only$vj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG, contig_paired_only$vj_gene_cdr3_BD, sep = " & ")
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub("^ & ", "", contig_paired_only$vj_gene_cdr3_AG_BD)
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub(" & $", "", contig_paired_only$vj_gene_cdr3_AG_BD)

      contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG, contig_paired_only$vdj_gene_cdr3_BD, sep = " & ")
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub("^ & ", "", contig_paired_only$vdj_gene_cdr3_AG_BD)
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub(" & $", "", contig_paired_only$vdj_gene_cdr3_AG_BD)
      # contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vdj_gene_cdr3_BD,sep = " & ")
      names(contig_paired_only)[names(contig_paired_only) %in% "cell_id"] <- "Cell_Index"
      contig_paired_only <- contig_paired_only[!duplicated(contig_paired_only$Cell_Index), ] # remove duplicates


      contig_paired_only <- merge(contig_paired_only, sample_tags, by = c("Cell_Index", "Sample_Name"), all = T)

      if (input$filtered_list == "Unfiltered") {
        contig_paired_only <- contig_paired_only %>%
          select(
            all_of(c(
              "Cell_Index", "Sample_Name", "Sample_Tag", "pairing", "pairing_type", "Clonality", "seq_identified", "seq_identified_TRA", "seq_identified_TRB", "seq_identified_TRG", "seq_identified_TRD", "cell_type_experimental_AG", "cell_type_experimental_BD", "chain_AG", "chain_BD",
              names(contig_paired_only[grep("call", names(contig_paired_only))])
            )),
            everything()
          )
      } else {
        contig_paired_only <- contig_paired_only %>%
          select(
            all_of(c(
              "Cell_Index", "Sample_Name", "Sample_Tag", "cell_type_experimental_AG", "cell_type_experimental_BD", "chain_AG", "chain_BD",
              names(contig_paired_only[grep("call", names(contig_paired_only))])
            )),
            everything()
          )
      }


      # names(contig_paired_only) <- gsub("_aa","",names(contig_paired_only) )
      contig_paired_only
      # merge(sample_tags,contig_paired_only,by.x="cell_id",by.y="Cell_Index",all=T)
    }

    filtering_required_BCR <- reactive({
      contigs <- input.data.TCR.bd2()

      validate(
        need(
          nrow(contigs) > 0,
          "Upload files"
        )
      )

      TCR_unfiltered <- contigs
      TCR_unfiltered$seq_issue <- ifelse(grepl("[*]", TCR_unfiltered$sequence_alignment_aa), "stop-codon", "productive")
      TCR_unfiltered_prod <- subset(TCR_unfiltered, TCR_unfiltered$seq_issue == "productive")
      print(names(TCR_unfiltered_prod))
      if(length(TCR_unfiltered_prod$cdr3_length)>0) {
      } else {
        message("adding junction length")
        TCR_unfiltered_prod$cdr3_length <- nchar(TCR_unfiltered_prod$junction_aa)
      }

      TCR_unfiltered_prod <- subset(TCR_unfiltered_prod, TCR_unfiltered_prod$cdr3_length < 30)
      TCR_unfiltered_prod <- subset(TCR_unfiltered_prod, TCR_unfiltered_prod$cdr3_length > 5)
      df_unique <- TCR_unfiltered_prod[, names(TCR_unfiltered_prod) %in% c("locus", "junction_aa")]
      df_unique$Clonal_Expanded <- 1
      df_unique_sum <- ddply(df_unique, names(df_unique)[-c(3)], numcolwise(sum))

      count.chain <- as.data.frame(table(TCR_unfiltered_prod_junc$cell_id))
      names(count.chain) <- c("cell_id", "seq_identified")

      TCR_unfiltered_prod_junc_seq <- merge(TCR_unfiltered_prod_junc, count.chain, by = "cell_id")
      sum_tab <- TCR_unfiltered_prod_junc_seq[, names(TCR_unfiltered_prod_junc_seq) %in% c("locus", "cell_id")]
      sum_tab$count <- 1
      sum_tab_2 <- as.data.frame(ddply(sum_tab, names(sum_tab)[1:2], numcolwise(sum)))
      count.chain <- as.data.frame(table(TCR_unfiltered_prod_junc_seq$cell_id, TCR_unfiltered_prod_junc_seq$locus))
      names(count.chain) <- c("cell_id", "chain", "seq_identified")

      IGH <- subset(count.chain, count.chain$chain == "IgH")
      names(IGH)[2:3] <- paste0(names(IGH), "_IGH")[2:3]
      IGL <- subset(count.chain, count.chain$chain == "IGL")
      names(IGL)[2:3] <- paste0(names(TRB), "_IGL")[2:3]
      IGK <- subset(count.chain, count.chain$chain == "IGK")
      names(IGK)[2:3] <- paste0(names(TRD), "_IGK")[2:3]

      IgHL <- merge(IGH, IGL, by = "cell_id")
      IgHK <- merge(IGH, IGK, by = "cell_id")
      IgHL_IgHK <- merge(IgHL, IgHK, by = "cell_id")
      # TRAB_TRGD_2 <- TRAB_TRGD %>%
      #   mutate(Clonality = case_when(
      #     seq_identified_TRA > 0 & seq_identified_TRB > 0 & seq_identified_TRG > 0 & seq_identified_TRD > 0 ~ "AB GD", # all
      #     # three chains called
      #     seq_identified_TRA > 0 & seq_identified_TRB > 0 & seq_identified_TRG > 0 & seq_identified_TRD == 0 ~ "AB G",
      #     seq_identified_TRA > 0 & seq_identified_TRB > 0 & seq_identified_TRG == 0 & seq_identified_TRD > 0 ~ "AB D",
      #     seq_identified_TRA > 0 & seq_identified_TRB == 0 & seq_identified_TRG > 0 & seq_identified_TRD > 0 ~ "A GD",
      #     seq_identified_TRA == 0 & seq_identified_TRB > 0 & seq_identified_TRG > 0 & seq_identified_TRD > 0 ~ "B GD",
      #     # two chains called
      #     seq_identified_TRA > 0 & seq_identified_TRB > 0 & seq_identified_TRG == 0 & seq_identified_TRD == 0 ~ "AB",
      #     seq_identified_TRA == 0 & seq_identified_TRB == 0 & seq_identified_TRG > 0 & seq_identified_TRD > 0 ~ "GD",
      #     seq_identified_TRA == 0 & seq_identified_TRB > 0 & seq_identified_TRG > 0 & seq_identified_TRD == 0 ~ "B G",
      #     seq_identified_TRA == 0 & seq_identified_TRB > 0 & seq_identified_TRG == 0 & seq_identified_TRD > 0 ~ "B D",
      #     seq_identified_TRA > 0 & seq_identified_TRB == 0 & seq_identified_TRG == 0 & seq_identified_TRD > 0 ~ "A D",
      #     seq_identified_TRA > 0 & seq_identified_TRB == 0 & seq_identified_TRG > 0 & seq_identified_TRD == 0 ~ "A G",
      #     # one chains called
      #     seq_identified_TRA > 0 & seq_identified_TRB == 0 & seq_identified_TRG == 0 & seq_identified_TRD == 0 ~ "A",
      #     seq_identified_TRA == 0 & seq_identified_TRB > 0 & seq_identified_TRG == 0 & seq_identified_TRD == 0 ~ "B",
      #     seq_identified_TRA == 0 & seq_identified_TRB == 0 & seq_identified_TRG > 0 & seq_identified_TRD == 0 ~ "G",
      #     seq_identified_TRA == 0 & seq_identified_TRB == 0 & seq_identified_TRG == 0 & seq_identified_TRD > 0 ~ "D",
      #     TRUE ~ "Other"
      #   ))
      # TRAB_TRGD_paired <- TRAB_TRGD_2
      # TRAB_TRGD_paired <- TRAB_TRGD_paired[, !grepl("chain_", names(TRAB_TRGD_paired))]
      # TRAB_TRGD_paired$sum <- rowSums(TRAB_TRGD_paired[2:5])
      # names(TRAB_TRGD_paired)
      # head(TRAB_TRGD_paired)
      # TRAB_TRGD_paired$pairing_type <- ifelse(TRAB_TRGD_paired$Clonality == "AB" & TRAB_TRGD_paired$sum == 2, "standard AB",
      #                                         ifelse(TRAB_TRGD_paired$Clonality == "GD" & TRAB_TRGD_paired$sum == 2, "standard GD",
      #                                                ifelse(TRAB_TRGD_paired$sum == 1, "Unpaired",
      #                                                       "non-standard"
      #                                                )
      #                                         )
      # )
      #
      # df <- merge(TCR_unfiltered_prod_junc_seq, TRAB_TRGD_paired, by = "cell_id")
      # df_2 <- df
      # Standard_AB <- subset(df_2, df_2$pairing_type == "standard AB")
      # Standard_GD <- subset(df_2, df_2$pairing_type == "standard GD")
      # non_standard <- subset(df_2, df_2$pairing_type == "non-standard")
      # non_standard <- non_standard[order(non_standard$consensus_count, decreasing = T), ]
      # non_standard$filter <- paste(non_standard$cell_id, non_standard$locus)
      # non_standard2 <- non_standard[!duplicated(non_standard$filter), ]
      # non_standard2 <- non_standard2[order(non_standard2$cell_id, decreasing = F), ]
      # unique(non_standard2$Clonality)
      # non_standard2$Standard_to_keep <- ifelse(non_standard2$Clonality == "AB G" & non_standard2$locus == "TRG", "remove_G",
      #                                          ifelse(non_standard2$Clonality == "AB D" & non_standard2$locus == "TRD", "remove_D",
      #                                                 ifelse(non_standard2$Clonality == "B GD" & non_standard2$locus == "TRB", "remove_B",
      #                                                        ifelse(non_standard2$Clonality == "A GD" & non_standard2$locus == "TRA", "remove_A",
      #                                                               ifelse(non_standard2$Clonality == "AB GD" & c(non_standard2$locus == "TRG" | non_standard2$locus == "TRD"), "remove_GD", "keep")
      #                                                        )
      #                                                 )
      #                                          )
      # )
      #
      # non_standard <- subset(non_standard2, non_standard2$Standard_to_keep == "keep")
      # non_standard <- non_standard[, !names(non_standard) %in% c("filter", "Standard_to_keep")]
      # dim(non_standard)
      #
      #
      # Filtered_TCR_list <- rbind(Standard_AB, Standard_GD, non_standard)
      # Filtered_TCR_list <- Filtered_TCR_list[order(Filtered_TCR_list$cell_id, decreasing = F), ]
      # contigs <- Filtered_TCR_list
      # contigs
    })

    tb_bd_meta.data_BCR <- function() {
      contigs <- input.data.TCR.bd2()
      sample_tags <- input.data.calls.bd()
      validate(
        need(
          nrow(contigs) > 0 && nrow(sample_tags) > 0,
          "Upload files"
        )
      )

      if (input$filtered_list == "Unfiltered") {
        contigs <- filtering_required_BCR()
      } else {
        contigs <- contigs
      }
      names(contigs)[names(contigs) %in% "cell_id"] <- "Cell_Index"
      contigs_merge <- merge(contigs, sample_tags, by = "Cell_Index")

      if (nrow(contigs_merge[-c(grep("[*]", contigs_merge$junction_aa)), ] > 0)) {
        contigs_merge <- contigs_merge[-c(grep("[*]", contigs_merge$junction_aa)), ]
      }

      # removed undefined sequences

      contigs_merge <- subset(contigs_merge, contigs_merge$Sample_Name != "Multiplet")
      contigs_merge <- subset(contigs_merge, contigs_merge$Sample_Name != "Undetermined")

      if (input$filtered_list == "Unfiltered") {
        contigs_lim <- contigs_merge[!names(contigs_merge) %in% c(
          "consensus_count", "sequence_id", "duplicate_count", "germline_alignment", "reads", "length", "cdr3", "rev_comp", "complete_vdj", names(contigs_merge[grep("fwr", names(contigs_merge))]), "cdr1", "cdr2", "productive", names(contigs_merge[grep("sequence", names(contigs_merge))]), names(contigs_merge[grep("cigar", names(contigs_merge))]), names(contigs_merge[grep("support", names(contigs_merge))]), "Sample_Tag", "sum", "dominant", "putative_cell", "juct_issue", "seq_issue", "c_call"
          # "sum","dominant","putative_cell","seq_issue","juct_issue"
        )]
      } else {
        contigs_lim <- contigs_merge[!names(contigs_merge) %in% c(
          "consensus_count", "sequence_id", "duplicate_count", "germline_alignment", "reads", "length", "cdr3", "rev_comp", "complete_vdj", names(contigs_merge[grep("fwr", names(contigs_merge))]), "cdr1", "cdr2", "productive", names(contigs_merge[grep("sequence", names(contigs_merge))]), names(contigs_merge[grep("cigar", names(contigs_merge))]), names(contigs_merge[grep("support", names(contigs_merge))]), "Sample_Tag"
          # "sum","dominant","putative_cell","seq_issue","juct_issue"
        )]
      }

      contigs_lim$v_gene <- gsub("[*]0.", "", contigs_lim$v_call)
      contigs_lim$j_gene <- gsub("[*]0.", "", contigs_lim$j_call)
      contigs_lim$d_gene <- gsub("[*]0.", "", contigs_lim$d_call)

      names(contigs_lim)[names(contigs_lim) %in% input$locus_column] <- "chain"

      contig_IgH <- subset(contigs_lim, contigs_lim$chain == "IGH")

      if (input$filtered_list == "Unfiltered") {
        contig_IgH <- contig_IgH[!duplicated(contig_IgH$Cell_Index), ]
        name.list <- names(contig_IgH[c(
          names(contig_IgH[grep("gene", names(contig_IgH))]),
          names(contig_IgH[grep("call", names(contig_IgH))]),
          names(contig_IgH[grep("cdr3", names(contig_IgH))]),
          names(contig_IgH[grep("cdr2", names(contig_IgH))]),
          names(contig_IgH[grep("cdr1", names(contig_IgH))]),
          names(contig_IgH[grep("junction", names(contig_IgH))]),
          "chain",
          "cell_type_experimental",
          "Clonal_Expanded"
        )])
      } else {
        name.list <- names(contig_IgH[c(
          names(contig_IgH[grep("gene", names(contig_IgH))]),
          names(contig_IgH[grep("call", names(contig_IgH))]),
          names(contig_IgH[grep("cdr3", names(contig_IgH))]),
          names(contig_IgH[grep("cdr2", names(contig_IgH))]),
          names(contig_IgH[grep("cdr1", names(contig_IgH))]),
          names(contig_IgH[grep("junction", names(contig_IgH))]),
          "chain",
          "cell_type_experimental"
        )])
      }

      contig_IgH <- contig_IgH %>%
        select(all_of(name.list), everything())
      names(contig_IgH)[1:summary(name.list)[1]] <- paste(names(contig_IgH[names(contig_IgH) %in% name.list]), "_IgH", sep = "")


      contig_IgLK <- subset(contigs_lim, contigs_lim$chain == "IGK" | contigs_lim$chain == "IGL")

      if (input$filtered_list == "Unfiltered") {
        name.list <- names(contig_IgLK[c(
          names(contig_IgLK[grep("gene", names(contig_IgLK))]),
          names(contig_IgLK[grep("call", names(contig_IgLK))]),
          names(contig_IgLK[grep("cdr3", names(contig_IgLK))]),
          names(contig_IgLK[grep("cdr2", names(contig_IgLK))]),
          names(contig_IgLK[grep("cdr1", names(contig_IgLK))]),
          names(contig_IgLK[grep("junction", names(contig_IgLK))]),
          "chain",
          "cell_type_experimental",
          "Clonal_Expanded"
        )])
      } else {
        name.list <- names(contig_IgLK[c(
          names(contig_IgLK[grep("gene", names(contig_IgLK))]),
          names(contig_IgLK[grep("call", names(contig_IgLK))]),
          names(contig_IgLK[grep("cdr3", names(contig_IgLK))]),
          names(contig_IgLK[grep("cdr2", names(contig_IgLK))]),
          names(contig_IgLK[grep("cdr1", names(contig_IgLK))]),
          names(contig_IgLK[grep("junction", names(contig_IgLK))]),
          "chain",
          "cell_type_experimental"
        )])
      }

      contig_IgLK <- contig_IgLK %>%
        select(all_of(name.list), everything())
      names(contig_IgLK)[1:summary(name.list)[1]] <- paste(names(contig_IgLK[names(contig_IgLK) %in% name.list]), "_IgLK", sep = "")

      if (input$filtered_list == "Unfiltered") {
        contig_paired <- merge(contig_IgH, contig_IgLK, by = c("Cell_Index", "Sample_Name"), all = T)
        # contig_paired <- merge(contig_AG, contig_BD, by = c("Cell_Index", "Sample_Name", names(contig_BD[grep("seq_", names(contig_BD))]), "Clonality", "pairing_type"), all = T)
      } else {
        contig_paired <- merge(contig_IgH, contig_IgLK, by = c("Cell_Index", "Sample_Name"), all = T)
      }
      print("merged IgH and IgLK")



      contig_paired_only <- contig_paired
      contig_paired_only[is.na(contig_paired_only)] <- "None"
      contig_paired_only[contig_paired_only == ""] <- "None"
      contig_paired_only$d_gene_IgH[contig_paired_only$d_gene_IgH == "None"] <- "_"
      print(names(contig_paired_only))
      if (input$filtering_TCR == T) {
        contig_paired_only <- subset(contig_paired_only, contig_paired_only$junction_IgH != "None")
        contig_paired_only <- subset(contig_paired_only, contig_paired_only$junction_IgLK != "None")
        print(dim(contig_paired_only))
      }


      # heavy chain
      contig_paired_only$vj_gene_IgH <- paste(contig_paired_only$v_gene_IgH, contig_paired_only$j_gene_IgH, sep = ".")
      if (length(gsub("None.None", "", contig_paired_only$vj_gene_IgLK))>0) {
        contig_paired_only$vj_gene_IgH <- gsub("None.None", "", contig_paired_only$vj_gene_IgH)
      }
      contig_paired_only$vj_gene_IgH <- paste(contig_paired_only$v_gene_IgH, contig_paired_only$j_gene_IgH, sep = ".")
      contig_paired_only$vj_gene_IgH <- gsub(".NA.", ".", contig_paired_only$vj_gene_IgH)
      contig_paired_only$vj_gene_IgH <- gsub("[.]None[.]", ".", contig_paired_only$vj_gene_IgH)
      contig_paired_only$vj_gene_IgH <- gsub("None.None", "", contig_paired_only$vj_gene_IgH)
      print(names(contig_paired_only))
      contig_paired_only$vj_gene_cdr3_IgH <- paste(contig_paired_only$vj_gene_IgH, contig_paired_only$junction_aa_IgH, sep = "_")
      contig_paired_only$vj_gene_cdr3_IgH <- gsub("_None", "", contig_paired_only$vj_gene_cdr3_IgH)
      print(names(contig_paired_only))
      contig_paired_only$vdj_gene_IgH <- paste(contig_paired_only$v_gene_IgH,contig_paired_only$d_gene_IgH, contig_paired_only$j_gene_IgH, sep = ".")
      contig_paired_only$vdj_gene_IgH <- gsub(".NA.", ".", contig_paired_only$vdj_gene_IgH)
      contig_paired_only$vdj_gene_IgH <- gsub("[.]None[.]", ".", contig_paired_only$vdj_gene_IgH)
      contig_paired_only$vdj_gene_IgH <- gsub("None.None", "", contig_paired_only$vdj_gene_IgH)
      print(names(contig_paired_only))
      contig_paired_only$vdj_gene_cdr3_IgH <- paste(contig_paired_only$vdj_gene_IgH, contig_paired_only$junction_aa_IgH, sep = "_")
      contig_paired_only$vdj_gene_cdr3_IgH <- gsub("_None", "", contig_paired_only$vdj_gene_cdr3_IgH)

      print(names(contig_paired_only))
      # light and kappa


      contig_paired_only$vj_gene_IgLK <- paste(contig_paired_only$v_gene_IgLK, contig_paired_only$j_gene_IgLK, sep = ".")
      if (length(gsub("None.None", "", contig_paired_only$vj_gene_IgLK))>0) {
        contig_paired_only$vj_gene_IgLK <- gsub("None.None", "", contig_paired_only$vj_gene_IgLK)
      }

      #
      contig_paired_only$vj_gene_IgLK <- paste(contig_paired_only$v_gene_IgLK, contig_paired_only$j_gene_IgLK, sep = ".")
      contig_paired_only$vj_gene_IgLK <- gsub(".NA.", ".", contig_paired_only$vj_gene_IgLK)
      contig_paired_only$vj_gene_IgLK <- gsub("[.]None[.]", ".", contig_paired_only$vj_gene_IgLK)
      contig_paired_only$vj_gene_IgLK <- gsub("None.None", "", contig_paired_only$vj_gene_IgLK)
      contig_paired_only$vj_gene_cdr3_IgLK <- paste(contig_paired_only$vj_gene_IgLK, contig_paired_only$junction_aa_IgLK, sep = "_")
      contig_paired_only$vj_gene_cdr3_IgLK <- gsub("_None", "", contig_paired_only$vj_gene_cdr3_IgLK)

      contig_paired_only$vj_gene_IgH_IgLK <- paste(contig_paired_only$vj_gene_IgH, contig_paired_only$vj_gene_IgLK, sep = " & ")
      contig_paired_only$vdj_gene_IgH_IgLK <- paste(contig_paired_only$vdj_gene_IgH, contig_paired_only$vj_gene_IgLK, sep = " & ")
      contig_paired_only$vdj_gene_IgH_IgLK <- gsub("^ & ", "", contig_paired_only$vdj_gene_IgH_IgLK)
      contig_paired_only$vdj_gene_IgH_IgLK <- gsub(" & $", "", contig_paired_only$vdj_gene_IgH_IgLK)

      #updating names to be consistant....
      contig_paired_only$vj_gene_cdr3_IgH_IgLK <- paste(contig_paired_only$vj_gene_cdr3_IgH, contig_paired_only$vj_gene_cdr3_IgLK, sep = " & ")
      contig_paired_only$vj_gene_cdr3_IgH_IgLK <- gsub("^ & ", "", contig_paired_only$vj_gene_cdr3_IgH_IgLK)
      contig_paired_only$vj_gene_cdr3_IgH_IgLK <- gsub(" & $", "", contig_paired_only$vj_gene_cdr3_IgH_IgLK)

      contig_paired_only$vdj_gene_cdr3_IgH_IgLK <- paste(contig_paired_only$vdj_gene_cdr3_IgH, contig_paired_only$vj_gene_cdr3_IgLK, sep = " & ")
      contig_paired_only$vdj_gene_cdr3_IgH_IgLK <- gsub("^ & ", "", contig_paired_only$vdj_gene_cdr3_IgH_IgLK)
      contig_paired_only$vdj_gene_cdr3_IgH_IgLK <- gsub(" & $", "", contig_paired_only$vdj_gene_cdr3_IgH_IgLK)

      names(contig_paired_only)[names(contig_paired_only) %in% "cell_id"] <- "Cell_Index"
      print(names(contig_paired_only))

      contig_paired_only <- contig_paired_only[!duplicated(contig_paired_only$Cell_Index), ] # remove duplicates

      # asdf
      contig_paired_only <- merge(contig_paired_only, sample_tags, by = c("Cell_Index", "Sample_Name"), all = T)

      if (input$filtered_list == "Unfiltered") {
        contig_paired_only <- contig_paired_only %>%
          select(
            all_of(c(
              "Cell_Index", "Sample_Name", "Sample_Tag", "pairing", "pairing_type", "Clonality", "seq_identified", "seq_identified_TRA", "seq_identified_TRB", "seq_identified_TRG", "seq_identified_TRD", "cell_type_experimental_IgH", "cell_type_experimental_IgLK", "chain_IgH", "chain_IgLK",
              names(contig_paired_only[grep("call", names(contig_paired_only))])
            )),
            everything()
          )
      } else {
        contig_paired_only <- contig_paired_only %>%
          select(
            all_of(c(
              "Cell_Index", "Sample_Name", "Sample_Tag", "cell_type_experimental_IgH", "cell_type_experimental_IgLK", "chain_IgH", "chain_IgLK",
              names(contig_paired_only[grep("call", names(contig_paired_only))])
            )),
            everything()
          )
      }

      contig_paired_only
      # merge(sample_tags,contig_paired_only,by.x="cell_id",by.y="Cell_Index",all=T)
    }

    BCR_TCR_md <- reactive({
      TCR <- tb_bd_meta.data_TCR()
      BCR <- tb_bd_meta.data_BCR()

      TCR_BCR <- merge(TCR,BCR,by=c("Cell_Index","Sample_Name","Sample_Tag"),all = T)
      TCR_BCR


    })

    ## create Sample_tags file =====

    samp.tags <- reactive({
      if (input$filtered_list == "Paired") {
        TCR <- as.data.frame(TCR_Filtering_Paired())
      } else {
        TCR <- input.data.TCR.bd2()
      }


      validate(
        need(
          nrow(TCR) > 0,
          error_message_val2
        )
      )
      if (input$filtered_list == "Paired") {
        Sample_Tags <- as.data.frame(TCR$Cell_Index)
      } else {
        Sample_Tags <- as.data.frame(unique(TCR$cell_id))
      }

      names(Sample_Tags) <- "Cell_Index"
      Sample_Tags$Sample_Tag <- "SampleTag01"
      #
      Sample_Tags$Sample_Name <- input$sample_tags_name

      Sample_Tags <- Sample_Tags[order(Sample_Tags$Cell_Index), ]
      # Sample_Tags
      Sample_Tags
    })

    output$tb_sample_tags_created <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      samp.tags()
    })

    output$downloadtb_sample_tags <- downloadHandler(
      filename = function() {
        if (input$name_BD == "") {
          showModal(modalDialog(
            title = "Error",
            "Please Add File Name.",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {

          x <- today()
          paste(input$name_BD, "_Sample_Tags_", x, ".csv", sep = "")
        }
      },
      content = function(file) {
        df <- as.data.frame(samp.tags())
        write.csv(df, file, row.names = F)
      }
    )

    ### testing table merge ------
    observe({
      updateSelectInput(
        session,
        "V_gene_AG_BDrhap",
        choices = names(input.data.TCR.BD()),
        selected = "TCR_Alpha_Gamma_V_gene_Dominant"
      )
    }) # junction sequence
    observe({
      updateSelectInput(
        session,
        "Junction_AG_BDrhap",
        choices = names(input.data.TCR.BD()),
        selected = "TCR_Alpha_Gamma_CDR3_Translation_Dominant"
      )
    }) # junction sequence
    observe({
      updateSelectInput(
        session,
        "V_gene_BD_BDrhap",
        choices = names(input.data.TCR.BD()),
        selected = "TCR_Beta_Delta_V_gene_Dominant"
      )
    }) # junction sequence
    observe({
      updateSelectInput(
        session,
        "Junction_BD_BDrhap",
        choices = names(input.data.TCR.BD()),
        selected = "TCR_Beta_Delta_CDR3_Translation_Dominant"
      )
    }) # junction sequence

    output$Check_table <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      df_nocouts <- TCR_Filtering_Paired()

      if (nrow(df_nocouts[df_nocouts$Total_VDJ_Read_Count != 0, ] > 0)) {
        df_nocouts2 <- df_nocouts[df_nocouts$Total_VDJ_Read_Count != 0, ]
      }
      df_nocouts2

      df_nocouts2_AG <- df_nocouts2[, names(df_nocouts2) %in% c("Sample_Name", input$V_gene_AG_BDrhap, input$Junction_AG_BDrhap)]
      names(df_nocouts2_AG) <- c("Sample_Name", "v_call", "junction_aa")
      df_nocouts2_AG

      df_nocouts2_BD <- df_nocouts2[, names(df_nocouts2) %in% c("Sample_Name", input$V_gene_BD_BDrhap, input$Junction_BD_BDrhap, "", "")]
      names(df_nocouts2_BD) <- c("Sample_Name", "v_call", "junction_aa")

      df_nocouts3 <- as.data.frame(rbind(df_nocouts2_AG, df_nocouts2_BD))
      df_nocouts3$v_call <- gsub("^$", NA, df_nocouts3$v_call)
      df_nocouts3 <- df_nocouts3[complete.cases(df_nocouts3$v_call), ]

      df_nocouts3$junction_aa <- gsub("^$", NA, df_nocouts3$junction_aa)
      df_nocouts3 <- df_nocouts3[complete.cases(df_nocouts3$junction_aa), ]

      if (nrow(df_nocouts3[-c(grep("[*]", df_nocouts3$junction_aa)), ] > 0)) {
        df_nocouts3 <- df_nocouts3[-c(grep("[*]", df_nocouts3$junction_aa)), ]
      }


      df_nocouts3 <- df_nocouts3[!duplicated(df_nocouts3$junction_aa), ]
      df_nocouts3
    })

    ## for the clusTCR ----
    df_clusTCR <- function() {
      df_nocouts <- TCR_Filtering_Paired()

      if (nrow(df_nocouts[df_nocouts$Total_VDJ_Read_Count != 0, ] > 0)) {
        df_nocouts2 <- df_nocouts[df_nocouts$Total_VDJ_Read_Count != 0, ]
      }

      #
      df_nocouts2_AG <- df_nocouts2[, names(df_nocouts2) %in% c("Sample_Name", input$V_gene_AG_BDrhap, input$Junction_AG_BDrhap)]
      names(df_nocouts2_AG) <- c("Sample_Name", "v_call", "junction_aa")
      df_nocouts2_AG

      df_nocouts2_BD <- df_nocouts2[, names(df_nocouts2) %in% c("Sample_Name", input$V_gene_BD_BDrhap, input$Junction_BD_BDrhap)]
      names(df_nocouts2_BD) <- c("Sample_Name", "v_call", "junction_aa")

      if (input$BCR_present == T) {
        df_nocouts2_IgL <- df_nocouts2[, names(df_nocouts2) %in% c("Sample_Name", "BCR_Light_CDR3_Translation_Dominant", "BCR_Light_V_gene_Dominant")]
        df_nocouts2_IgH <- df_nocouts2[, names(df_nocouts2) %in% c("Sample_Name", "BCR_Heavy_CDR3_Translation_Dominant", "BCR_Heavy_V_gene_Dominant")]
        names(df_nocouts2_IgL) <- c("Sample_Name", "v_call", "junction_aa")
        names(df_nocouts2_IgH) <- c("Sample_Name", "v_call", "junction_aa")
        df_nocouts3 <- as.data.frame(rbind(df_nocouts2_AG, df_nocouts2_BD, df_nocouts2_IgL, df_nocouts2_IgH))
      } else {
        df_nocouts3 <- as.data.frame(rbind(df_nocouts2_AG, df_nocouts2_BD))
      }

      df_nocouts3$v_call <- gsub("^$", NA, df_nocouts3$v_call)
      df_nocouts3 <- df_nocouts3[complete.cases(df_nocouts3$v_call), ]

      df_nocouts3$junction_aa <- gsub("^$", NA, df_nocouts3$junction_aa)
      df_nocouts3 <- df_nocouts3[complete.cases(df_nocouts3$junction_aa), ]

      if (nrow(df_nocouts3[-c(grep("[*]", df_nocouts3$junction_aa)), ] > 0)) {
        df_nocouts3 <- df_nocouts3[-c(grep("[*]", df_nocouts3$junction_aa)), ]
      }

      df_nocouts3 <- df_nocouts3[!duplicated(df_nocouts3$junction_aa), ]
      df_nocouts3
    }

    tb_bd_contigues_contig <- reactive({
      contigs <- input.data.TCR.bd2()
      sample_tags <- input.data.calls.bd()
      validate(
        need(
          nrow(contigs) > 0 && nrow(sample_tags) > 0,
          "Upload files"
        )
      )

      if (input$filtered_list == "Unfiltered") {
        df <- filtering_required_TCR()
        contigs <- contigs[contigs$cell_id %in% unique(df$cell_id), ]
      } else {
        contigs <- contigs
      }

      contigs_merge <- merge(contigs, sample_tags, by.x = "cell_id", by.y = "Cell_Index")

      contigs_merge <- subset(contigs_merge, contigs_merge$Sample_Name != "Multiplet")
      contigs_merge <- subset(contigs_merge, contigs_merge$Sample_Name != "Undetermined")

      contigs2 <- contigs[names(contigs) %in% c("v_call", "junction_aa")]
      contigs2$v_call <- gsub("[*]0.", "", contigs2$v_call)

      # contigs2$Sample_Name <- sample_tags$Sample_Name
      contigs2[is.na(contigs2)] <- "None"
      contigs2[contigs2 == ""] <- "None"

      if (nrow(contigs2[-c(grep("[*]", contigs2$junction_aa)), ] > 0)) {
        contigs2 <- contigs2[-c(grep("[*]", contigs2$junction_aa)), ]
      }

      contigs2 <- subset(contigs2, contigs2$junction_aa != "None")
      contigs2[!duplicated(contigs2[, c("v_call", "junction_aa")]), ]
    })



    df_clusTCR_Filtered <- reactive({
      if (input$filtered_list == "Paired") {
        df <- df_clusTCR()
      } else {
        df <- tb_bd_contigues_contig()
        df <- df[grep("^C", df$junction_aa), ]
      }

      if (input$chain_clusTCR2_bd == "AG") {
        (rbind(df[grep("TRAV", df$v_call), ], df[grep("TRGV", df$v_call), ]))
      } else if (input$chain_clusTCR2_bd == "IgH") {
        rbind(df[grep("IGH", df$v_call), ])
      } else if (input$chain_clusTCR2_bd == "IgLK") {
        (rbind(df[grep("IGL", df$v_call), ], df[grep("IGK", df$v_call), ]))
      } else {
        (rbind(df[grep("TRBV", df$v_call), ], df[grep("TRDV", df$v_call), ]))
      }
    })

    output$tb_clusTCR <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      if (input$filtered_list == "Paired") {
        df_clusTCR()
      } else {
        df_clusTCR_Filtered()
      }
    })

    output$downloaddf_clusTCR <- downloadHandler(
      filename = function() {
        if (input$name_BD == "") {
          showModal(modalDialog(
            title = "Error",
            "Please Add File Name.",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          x <- today()
          if (input$chain_clusTCR2_bd == "AG") {
            paste("AG_", input$name_BD, "", "_ClusTCR_", x, ".csv", sep = "")
          } else if (input$chain_clusTCR2_bd == "IgH") {
            paste("IgH_", input$name_BD, "_ClusTCR_", x, ".csv", sep = "")
          } else if (input$chain_clusTCR2_bd == "IgLK") {
            paste("IgLK_", input$name_BD, "_ClusTCR_", x, ".csv", sep = "")
          } else {
            paste("BD_", input$name_BD, "_ClusTCR_", x, ".csv", sep = "")
          }
        }


      },
      content = function(file) {
        df <- as.data.frame(df_clusTCR_Filtered())
        write.csv(df, file, row.names = F)
      }
    )

    ## count table for Seurat----
    tb_bd_matrix <- function() {
      barcode <- input.data.barcode.bd()
      features <- input.data.features.bd2()

      mmMat <- as.data.frame(input.data.matrix.bd())

      validate(
        need(
          nrow(barcode) > 0 & nrow(features) > 0 & nrow(mmMat) > 0,
          "Upload files"
        )
      )

      rownames(mmMat) <- make.unique(features$V2) # select which feature to label genes...
      names(mmMat) <- barcode$V1
      mmMat$Gene_Name <- rownames(mmMat)
      mmMat <- mmMat %>%
        select(all_of("Gene_Name"), everything())
      mmMat
    }

    df_count.matrix_bd <- reactive({
      mat <- as.data.frame(input.data.count.BD())

      validate(
        need(
          nrow(mat) > 0,
          "Upload file"
        )
      )
      mat <- as.data.frame(mat)
      head(mat)[1:6]
      rownames(mat) <- mat$Cell_Index

      mat <- subset(mat, select = -Cell_Index)
      head(mat)[1:6]

      mat <- as.data.frame(t(mat))
      head(mat)[1:6]

      rownames(mat) <- gsub("^X", "", rownames(mat))
      head(mat)[1:6]
      names(mat) <- gsub("^X", "", names(mat))
      mat
    })

    output$tb_count_matrix <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      if (input$Format_bd == "cellXgene") {
        (df_count.matrix_bd())[1:6, 1:6]
      } else {
        tb_bd_matrix()[1:6, 1:6]
      }
    })
    output$downloadtb_count_matrix <- downloadHandler(
      filename = function() {

        if (input$name_BD == "") {
          showModal(modalDialog(
            title = "Error",
            "Please Add File Name.",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          if (input$Format_bd == "cellXgene") {
            paste(input$name_BD, "_BD_Count_Matrix_", gsub("-", ".", Sys.Date()), ".csv", sep = "")
          } else {
            paste(input$name_BD, "_count-matrix_bd_", gsub("-", ".", Sys.Date()), ".csv.gz", sep = "")
          }
        }
      },
      content = function(file) {
        if (input$Format_bd == "cellXgene") {
          df <- as.data.frame(df_count.matrix_bd())
          # write.table(,file, row.names = T)
          write.csv(df, file, row.names = T)
        } else {
          df <- as.data.frame(tb_bd_matrix())
          # write.table(,file, row.names = T)
          write_csv(df, gzfile(file))
        }
      }
    )

    ### meta.data for seruat ----
    meta.data_for_Seuratobj <- function() {
      if (input$filtered_list == "Paired") {
        contig_paired_only <- TCR_Filtering_Paired()
        contig_paired_only
        contig_paired_only$v_gene_AG <- gsub("[*]0.", "", contig_paired_only$TCR_Alpha_Gamma_V_gene_Dominant)
        contig_paired_only$j_gene_AG <- gsub("[*]0.", "", contig_paired_only$TCR_Alpha_Gamma_J_gene_Dominant)
        contig_paired_only$cdr3_AG <- gsub("[*]0.", "", contig_paired_only[, names(contig_paired_only) %in% input$Junction_AG_BDrhap])
        contig_paired_only$v_gene_BD <- gsub("[*]0.", "", contig_paired_only$TCR_Beta_Delta_V_gene_Dominant)
        contig_paired_only$j_gene_BD <- gsub("[*]0.", "", contig_paired_only$TCR_Beta_Delta_J_gene_Dominant)
        contig_paired_only$d_gene_BD <- gsub("[*]0.", "", contig_paired_only$TCR_Beta_Delta_D_gene_Dominant)
        contig_paired_only$cdr3_BD <- gsub("[*]0.", "", contig_paired_only[, names(contig_paired_only) %in% input$Junction_BD_BDrhap])

        contig_paired_only$d_gene_BD <- gsub("^$", "NA", contig_paired_only$d_gene_BD)
        #
        contig_paired_only$vj_gene_AG <- paste(contig_paired_only$v_gene_AG, contig_paired_only$j_gene_AG, sep = ".")
        contig_paired_only$vj_gene_AG <- gsub("NA.NA", "", contig_paired_only$vj_gene_AG)
        #
        contig_paired_only$vj_gene_BD <- paste(contig_paired_only$v_gene_BD, contig_paired_only$j_gene_BD, sep = ".")
        contig_paired_only$vj_gene_BD <- gsub(".NA.", ".", contig_paired_only$vj_gene_BD)
        contig_paired_only$vj_gene_BD <- gsub(".None.", ".", contig_paired_only$vj_gene_BD)
        contig_paired_only$vj_gene_BD <- gsub("NA.NA", "", contig_paired_only$vj_gene_BD)
        #
        contig_paired_only$vdj_gene_BD <- paste(contig_paired_only$v_gene_BD, contig_paired_only$d_gene_BD, contig_paired_only$j_gene_BD, sep = ".")
        contig_paired_only$vdj_gene_BD <- gsub(".NA.", ".", contig_paired_only$vdj_gene_BD)
        contig_paired_only$vdj_gene_BD <- gsub("[.]NA[.]", ".", contig_paired_only$vdj_gene_BD)
        contig_paired_only$vdj_gene_BD <- gsub("NA.NA", "", contig_paired_only$vdj_gene_BD)
        #
        contig_paired_only$vj_gene_cdr3_AG <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$cdr3_AG, sep = "_")
        contig_paired_only$vj_gene_cdr3_AG <- gsub("_NA", "", contig_paired_only$vj_gene_cdr3_AG)
        #
        contig_paired_only$vj_gene_cdr3_BD <- paste(contig_paired_only$vj_gene_BD, contig_paired_only$cdr3_BD, sep = "_")
        contig_paired_only$vj_gene_cdr3_BD <- gsub("_NA", "", contig_paired_only$vj_gene_cdr3_BD)
        #
        contig_paired_only$vdj_gene_cdr3_BD <- paste(contig_paired_only$vdj_gene_BD, contig_paired_only$cdr3_BD, sep = "_")
        contig_paired_only$vdj_gene_cdr3_BD <- gsub("_NA", "", contig_paired_only$vdj_gene_cdr3_BD)
        #
        contig_paired_only$vj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$vj_gene_BD, sep = " & ")
        contig_paired_only$vdj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$vdj_gene_BD, sep = " & ")
        contig_paired_only$vdj_gene_AG_BD <- gsub("^ & ", "", contig_paired_only$vdj_gene_AG_BD)
        contig_paired_only$vdj_gene_AG_BD <- gsub(" & $", "", contig_paired_only$vdj_gene_AG_BD)
        #
        # #updating names to be consistant....
        contig_paired_only$vj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG, contig_paired_only$vj_gene_cdr3_BD, sep = " & ")
        contig_paired_only$vj_gene_cdr3_AG_BD <- gsub("^ & ", "", contig_paired_only$vj_gene_cdr3_AG_BD)
        contig_paired_only$vj_gene_cdr3_AG_BD <- gsub(" & $", "", contig_paired_only$vj_gene_cdr3_AG_BD)

        contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG, contig_paired_only$vdj_gene_cdr3_BD, sep = " & ")
        contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub("^ & ", "", contig_paired_only$vdj_gene_cdr3_AG_BD)
        contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub(" & $", "", contig_paired_only$vdj_gene_cdr3_AG_BD)


        contig_paired_only <- contig_paired_only[!names(contig_paired_only) %in% c("Total_VDJ_Read_Count", "Total_VDJ_Molecule_Count", "TCR_Alpha_Gamma_Read_Count", "TCR_Alpha_Gamma_Molecule_Count", "TCR_Beta_Delta_Read_Count", "TCR_Beta_Delta_Molecule_Count", "Sample_Tag", "TCR_Alpha_Gamma_C_gene_Dominant", "TCR_Beta_Delta_C_gene_Dominant", "TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant", "TCR_Beta_Delta_CDR3_Nucleotide_Dominant")]
        contig_paired_only
        if (input$BCR_present == T) {
          contig_paired_only$v_gene_IgL <- gsub("[*]0.", "", contig_paired_only$BCR_Light_V_gene_Dominant)
          contig_paired_only$j_gene_IgL <- gsub("[*]0.", "", contig_paired_only$BCR_Light_J_gene_Dominant)
          contig_paired_only$v_gene_IgH <- gsub("[*]0.", "", contig_paired_only$BCR_Heavy_V_gene_Dominant)
          contig_paired_only$d_gene_IgH <- gsub("[*]0.", "", contig_paired_only$BCR_Heavy_D_gene_Dominant)
          contig_paired_only$j_gene_IgH <- gsub("[*]0.", "", contig_paired_only$BCR_Heavy_J_gene_Dominant)

          contig_paired_only$v_gene_IgH_L <- paste(contig_paired_only$v_gene_IgH, contig_paired_only$v_gene_IgL, sep = " & ")

          contig_paired_only$v_gene_cdr3_IgH <- paste(contig_paired_only$v_gene_IgH, contig_paired_only$BCR_Heavy_CDR3_Translation_Dominant, sep = "_")
          contig_paired_only$v_gene_cdr3_IgL <- paste(contig_paired_only$v_gene_IgL, contig_paired_only$BCR_Light_CDR3_Translation_Dominant, sep = "_")
          contig_paired_only$v_gene_cdr3_IgH_L <- paste(contig_paired_only$v_gene_cdr3_IgH, contig_paired_only$v_gene_cdr3_IgL, sep = ".")

          names(contig_paired_only)[names(contig_paired_only) %in% c("BCR_Light_CDR3_Translation_Dominant")] <- "cdr3_IgL"
          names(contig_paired_only)[names(contig_paired_only) %in% c("BCR_Light_V_gene_Dominant")] <- "v_allele_IgL"
          names(contig_paired_only)[names(contig_paired_only) %in% c("BCR_Light_J_gene_Dominant")] <- "j_allele_IgL"
          names(contig_paired_only)[names(contig_paired_only) %in% c("BCR_Heavy_CDR3_Translation_Dominant")] <- "cdr3_IgH"
          names(contig_paired_only)[names(contig_paired_only) %in% c("BCR_Heavy_V_gene_Dominant")] <- "v_allele_IgH"
          names(contig_paired_only)[names(contig_paired_only) %in% c("BCR_Heavy_J_gene_Dominant")] <- "j_allele_IgH"
          names(contig_paired_only)[names(contig_paired_only) %in% c("BCR_Heavy_D_gene_Dominant")] <- "d_allele_IgH"
          contig_paired_only
        } else {
          contig_paired_only
        }
      } else {

        if(!input$BCR_present && input$TCR_present){
          tb_bd_meta.data_TCR()
        } else if (input$BCR_present && !input$TCR_present) {
          tb_bd_meta.data_BCR()
        } else if (input$BCR_present && input$TCR_present) {
          BCR <-as.data.frame(tb_bd_meta.data_BCR())

          TCR <- tb_bd_meta.data_TCR()

          merge(TCR,BCR,by = c("Cell_Index","Sample_Name","Sample_Tag"),all = T)


        } else {
          as.data.frame("Select TCR and/or BCR")

        }



      }
    }

    output$tb_metadata_sc <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      meta.data_for_Seuratobj()
    })

    output$downloadtb_metadata_sc <- downloadHandler(
      filename = function() {
        if (input$name_BD == "") {
          showModal(modalDialog(
            title = "Error",
            "Please Add File Name.",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          x <- today()
          paste(input$name_BD, "_Meta.data_", x, ".csv", sep = "")
        }
      },
      content = function(file) {
        df <- as.data.frame(meta.data_for_Seuratobj())
        write.csv(df, file, row.names = F)
      }
    )


    # for TCRex output -----
    TCRex_BDrap_df <- reactive({
      calls_TCR_paired.fun <- TCR_Filtering_Paired()


      if (nrow(calls_TCR_paired.fun[-c(grep("[*]", calls_TCR_paired.fun$TCR_Beta_Delta_CDR3_Translation_Dominant)), ] > 0)) {
        calls_TCR_paired.fun <- calls_TCR_paired.fun[-c(grep("[*]", calls_TCR_paired.fun$TCR_Beta_Delta_CDR3_Translation_Dominant)), ] # remove st
      }

      calls_TCR_paired.fun$TRBV_gene <- gsub("[*]0.", "", calls_TCR_paired.fun$TCR_Beta_Delta_V_gene_Dominant)
      calls_TCR_paired.fun$CDR3_beta <- paste("C", calls_TCR_paired.fun$TCR_Beta_Delta_CDR3_Translation_Dominant, "F", sep = "")
      calls_TCR_paired.fun$TRBJ_gene <- gsub("[*]0.", "", calls_TCR_paired.fun$TCR_Beta_Delta_J_gene_Dominant)
      calls_TCR_paired.fun$cloneCount <- 1
      calls_TCR_paired.fun2 <- calls_TCR_paired.fun[, names(calls_TCR_paired.fun) %in% c("TRBV_gene", "CDR3_beta", "TRBJ_gene", "cloneCount")]

      calls_TCR_paired.fun3 <- ddply(calls_TCR_paired.fun2, names(calls_TCR_paired.fun2)[-c(4)], numcolwise(sum))
      calls_TCR_paired.fun3 <- subset(calls_TCR_paired.fun3, calls_TCR_paired.fun3$CDR3_beta != "CF")

      calls_TCR_paired.fun3 <- calls_TCR_paired.fun3[!grepl("TRD", calls_TCR_paired.fun3$TRBV_gene), ]
      calls_TCR_paired.fun3
    })

    TCRex_Bd_df <- reactive({
      contigs <- input.data.TCR.bd2()
      sample_tags <- input.data.calls.bd()
      validate(
        need(
          nrow(contigs) > 0 && nrow(sample_tags) > 0,
          "Upload files"
        )
      )
      if (input$filtered_list == "Unfiltered") {
        df <- filtering_required_TCR()
        contigs <- contigs[contigs$cell_id %in% unique(df$cell_id), ]
      } else {
        contigs <- contigs
      }
      contigs_merge <- merge(contigs, sample_tags, by.x = "cell_id", by.y = "Cell_Index")

      contigs_merge <- subset(contigs_merge, contigs_merge$Sample_Name != "Multiplet")
      contigs_merge <- subset(contigs_merge, contigs_merge$Sample_Name != "Undetermined")

      contigs2 <- contigs_merge[, names(contigs_merge) %in% c("v_call", "j_call", "junction_aa")]
      contigs2[is.na(contigs2)] <- "None"
      contigs2[contigs2 == ""] <- "None"
      contigs2
      names(contigs2)[names(contigs2) %in% c("junction_aa")] <- "CDR3_beta"
      if (nrow(contigs2[-c(grep("[*]", contigs2$CDR3_beta)), ] > 0)) {
        contigs2 <- contigs2[-c(grep("[*]", contigs2$CDR3_beta)), ]
      }
      names(contigs2)[names(contigs2) %in% c("v_call")] <- "TRBV_gene"
      names(contigs2)[names(contigs2) %in% c("j_call")] <- "TRBJ_gene"
      contigs2 <- contigs2 %>%
        select(TRBV_gene, CDR3_beta, everything())
      contigs2 <- subset(contigs2, contigs2$CDR3_beta != "None")
      contigs2 <- contigs2
      contigs2$TRBV_gene <- gsub("[*]0.", "", contigs2$TRBV_gene)
      contigs2$TRBJ_gene <- gsub("[*]0.", "", contigs2$TRBJ_gene)
      contigs2$cloneCount <- 1
      calls_TCR_paired.fun3 <- ddply(contigs2, names(contigs2)[-c(4)], numcolwise(sum))
      calls_TCR_paired.fun3 <- calls_TCR_paired.fun3[grepl("TRB", calls_TCR_paired.fun3$TRBV_gene), ]
      calls_TCR_paired.fun3 <- calls_TCR_paired.fun3[,!names(calls_TCR_paired.fun3) %in% "cloneCount"]
      calls_TCR_paired.fun3
    })

    tb_TCRex_BDrap_df <- reactive({
      if (input$filtered_list == "Paired") {
        TCRex_BDrap_df()
      } else {
        TCRex_Bd_df()
      }
    })

    output$tb_TCRex_BDrap_df <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      tb_TCRex_BDrap_df()
    })


    output$downloaddf_TCRex_BDrap <- downloadHandler(
      filename = function() {
        if (input$name_BD == "") {
          showModal(modalDialog(
            title = "Error",
            "Please Add File Name.",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          x <- today()
          paste(input$name_BD, "_TCRex_", x, ".tsv", sep = "")
        }
      },
      content = function(file) {
        df <- as.data.frame(tb_TCRex_BDrap_df())
        write.table(df, file, row.names = F, sep = "\t", quote = F)
      }
    )

    ## TCR_Explore compatible -----
    df_TCR_Explore <- function() {
      calls_TCR_paired.fun <- meta.data_for_Seuratobj()
      calls_TCR_paired.fun$cloneCount <- 1 # added for TCR_Explore (total number of cells with that gene and sequence)

      calls_TCR_paired.fun <- calls_TCR_paired.fun[!names(calls_TCR_paired.fun) %in% c("Total_VDJ_Read_Count", "Total_VDJ_Molecule_Count", "TCR_Alpha_Gamma_Read_Count", "TCR_Alpha_Gamma_Molecule_Count", "TCR_Beta_Delta_Read_Count", "TCR_Beta_Delta_Molecule_Count", "Cell_Type_Experimental")]


      calls_TCR_paired.fun[order(calls_TCR_paired.fun$vdj_gene_AG_BD), ]

      calls_TCR_paired.fun$vdj_gene_AG_BD <- gsub("^$", NA, calls_TCR_paired.fun$vdj_gene_AG_BD)
      calls_TCR_paired.fun <- calls_TCR_paired.fun[complete.cases(calls_TCR_paired.fun$vdj_gene_AG_BD), ]

      calls_TCR_paired.fun <- calls_TCR_paired.fun %>%
        select(cloneCount, everything())

      calls_TCR_paired.fun
    }

    output$tb_TCR_Explore <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      df_TCR_Explore()
    })

    output$downloadtb_TCR_Explore <- downloadHandler(
      filename = function() {
        if (input$name_BD == "") {
          showModal(modalDialog(
            title = "Error",
            "Please Add File Name.",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          x <- today()
          paste(input$name_BD, "TCR_Explore.compatible_", x, ".csv", sep = "")
        }
      },
      content = function(file) {
        df <- as.data.frame(df_TCR_Explore())
        write.csv(df, file, row.names = F)
      }
    )


    # multi-TCR files -----

    multi_TCR_BDrhap <- reactive({
      contigs <- input.data.TCR.bd2()
      sample_tags <- input.data.calls.bd()
      validate(
        need(
          nrow(contigs) > 0 && nrow(sample_tags) > 0,
          "Upload files"
        )
      )

      if (input$filtered_list == "Unfiltered") {
        sample_tags$order <- 1:dim(sample_tags)[1]
        head(sample_tags)
        names(contigs)
        contigs_merge <- merge(contigs, sample_tags, by.x = "cell_id", by.y = "Cell_Index") # need this for BD rhap to remove multiplets/undetermined

        # remove non-functional sequences
        if (nrow(contigs_merge[-c(grep("[*]", contigs_merge$junction_aa)), ] > 0)) {
          contigs_merge <- contigs_merge[-c(grep("[*]", contigs_merge$junction_aa)), ]
        }
        names(contigs_merge)

        contigs_merge <- contigs_merge[contigs_merge[, names(contigs_merge) %in% "productive"] == "True", ]

        contigs_merge <- contigs_merge[grep("^C", contigs_merge$junction_aa), ]
        contigs_merge <- subset(contigs_merge, contigs_merge$Sample_Name != "Multiplet")
        contigs_merge <- subset(contigs_merge, contigs_merge$Sample_Name != "Undetermined")
        contigs_lim <- contigs_merge[!names(contigs_merge) %in% c("consensus_count", "sequence_id", "duplicate_count", "germline_alignment", "reads", "length", "cdr3", "rev_comp", "complete_vdj", names(contigs_merge[grep("fwr", names(contigs_merge))]), names(contigs_merge[grep("cdr1", names(contigs_merge))]), "junction", names(contigs_merge[grep("length", names(contigs_merge))]), "productive", names(contigs_merge[grep("cdr2", names(contigs_merge))]), "productive", names(contigs_merge[grep("sequence", names(contigs_merge))]), names(contigs_merge[grep("cigar", names(contigs_merge))]), names(contigs_merge[grep("support", names(contigs_merge))]), "Sample_Tag", "sum", "dominant", "putative_cell", "juct_issue", "seq_issue", "c_call", "cell_type_experimental")]
        contigs_lim$v_gene <- gsub("[*]0.", "", contigs_lim$v_call)
        contigs_lim$j_gene <- gsub("[*]0.", "", contigs_lim$j_call)
        contigs_lim$d_gene <- gsub("[*]0.", "", contigs_lim$d_call)
        contigs_lim[grep("call", names(contigs_lim))]

        names(contigs_lim) %in% names(contigs_lim[grep("call", names(contigs_lim))])

        contigs_lim <- contigs_lim[, !names(contigs_lim) %in% names(contigs_lim[grep("call", names(contigs_lim))])]
        head(contigs_lim)
        contigs_lim$cdr3_aa <- gsub("$^", "none", contigs_lim$cdr3_aa)
        contigs_lim <- subset(contigs_lim, contigs_lim$cdr3_aa != "none")

        # contigs_lim$v_gene_BD <- gsub("[*]0.","",calls_TCR_paired.fun$TCR_Beta_Delta_V_gene_Dominant)
        # names(contigs_lim)   <-  gsub("_call","_gene",names(contigs_lim))
        names(contigs_lim)

        names(contigs_lim)[names(contigs_lim) %in% "locus"] <- "chain" # input$locus_column

        contig_A <- subset(contigs_lim, contigs_lim$chain == "TRA")
        head(contig_A)
        name.list <- names(contig_A[c(
          names(contig_A[grep("gene", names(contig_A))]),
          names(contig_A[grep("call", names(contig_A))]),
          names(contig_A[grep("cdr3", names(contig_A))]),
          names(contig_A[grep("junction", names(contig_A))]),
          "chain"
        )])

        contig_A <- contig_A %>%
          select(all_of(name.list), everything())
        names(contig_A)[1:summary(name.list)[1]] <- paste(names(contig_A[names(contig_A) %in% name.list]), "_A", sep = "")

        names(contig_A) %in% c("cell_id", "Sample_Name")

        contig_G <- subset(contigs_lim, contigs_lim$chain == "TRG")
        head(contig_G)
        name.list <- names(contig_G[c(
          names(contig_G[grep("gene", names(contig_G))]),
          names(contig_G[grep("call", names(contig_G))]),
          names(contig_G[grep("cdr3", names(contig_G))]),
          names(contig_G[grep("junction", names(contig_G))]),
          "chain"
        )])
        head(contig_A)
        contig_G <- contig_G %>%
          select(all_of(name.list), everything())
        names(contig_G)[1:summary(name.list)[1]] <- paste(names(contig_G[names(contig_G) %in% name.list]), "_G", sep = "")

        names(contig_G) %in% c("cell_id", "Sample_Name")


        # contig_A

        contig_B <- subset(contigs_lim, contigs_lim$chain == "TRB")
        head(contig_B)
        name.list <- names(contig_B[c(
          names(contig_B[grep("gene", names(contig_B))]),
          names(contig_B[grep("call", names(contig_B))]),
          names(contig_B[grep("cdr3", names(contig_B))]),
          names(contig_B[grep("junction", names(contig_B))]),
          "chain"
        )])

        contig_B <- contig_B %>%
          select(all_of(name.list), everything())
        names(contig_B)[1:summary(name.list)[1]] <- paste(names(contig_B[names(contig_B) %in% name.list]), "_B", sep = "")

        names(contig_B) %in% c("cell_id", "Sample_Name")


        contig_D <- subset(contigs_lim, contigs_lim$chain == "TRD")
        head(contig_D)
        name.list <- names(contig_D[c(
          names(contig_D[grep("gene", names(contig_D))]),
          names(contig_D[grep("call", names(contig_D))]),
          names(contig_D[grep("cdr3", names(contig_D))]),
          names(contig_D[grep("junction", names(contig_D))]),
          "chain"
        )])

        contig_D <- contig_D %>%
          select(all_of(name.list), everything())
        names(contig_D)[1:summary(name.list)[1]] <- paste(names(contig_D[names(contig_D) %in% name.list]), "_D", sep = "")

        names(contig_D) %in% c("cell_id", "Sample_Name")

        contig_paired_AB <- merge(contig_A, contig_B, by = c("cell_id", "Sample_Name", "order"), all = T)
        contig_paired_GD <- merge(contig_G, contig_D, by = c("cell_id", "Sample_Name", "order"), all = T)
        contig_paired <- merge(contig_paired_AB, contig_paired_GD, by = c("cell_id", "Sample_Name", "order"), all = T)

        contig_paired <- contig_paired[!duplicated(contig_paired[names(contig_paired) %in% c("cell_id", "Sample_Name", "order", "junction_aa_A", "junction_aa_B", "junction_aa_G", "junction_aa_D")]), ]
        dim(contig_paired)

        contig_paired[is.na(contig_paired)] <- "-"
        contig_paired <- as.data.frame(contig_paired)

        head(contig_paired)

        contig_paired$pairing <- ifelse(contig_paired$chain_B == "TRB" & contig_paired$chain_A == "TRA" & contig_paired$chain_D == "TRD" & contig_paired$chain_G == "TRG", "AB GD",
                                        ifelse(contig_paired$chain_B == "TRB" & contig_paired$chain_A == "TRA" & contig_paired$chain_D != "TRD" & contig_paired$chain_G == "TRG", "AB G",
                                               ifelse(contig_paired$chain_B == "TRB" & contig_paired$chain_A == "TRA" & contig_paired$chain_D == "TRD" & contig_paired$chain_G != "TRG", "AB D",
                                                      ifelse(contig_paired$chain_B != "TRB" & contig_paired$chain_A == "TRA" & contig_paired$chain_D == "TRD" & contig_paired$chain_G == "TRG", "AG D",
                                                             ifelse(contig_paired$chain_B == "TRB" & contig_paired$chain_A != "TRA" & contig_paired$chain_D == "TRD" & contig_paired$chain_G == "TRG", "B GD",
                                                                    ifelse(contig_paired$chain_B == "TRB" & contig_paired$chain_A == "TRA" & contig_paired$chain_D != "TRD" & contig_paired$chain_G != "TRG", "AB",
                                                                           ifelse(contig_paired$chain_B != "TRB" & contig_paired$chain_A != "TRA" & contig_paired$chain_D == "TRD" & contig_paired$chain_G == "TRG", "GD",
                                                                                  ifelse(contig_paired$chain_B == "TRB" & contig_paired$chain_A != "TRA" & contig_paired$chain_D != "TRD" & contig_paired$chain_G == "TRG", "B G",
                                                                                         ifelse(contig_paired$chain_B != "TRB" & contig_paired$chain_A == "TRA" & contig_paired$chain_D == "TRD" & contig_paired$chain_G != "TRG", "A D",
                                                                                                "other"
                                                                                         )
                                                                                  )
                                                                           )
                                                                    )
                                                             )
                                                      )
                                               )
                                        )
        )





        unique(contig_paired$pairing)
        contig_paired

        contig_paired_only <- contig_paired
        head(contig_paired_only)
        # contig_paired_only <- subset(contig_paired_only,contig_paired_only$pairing!="single")

        AB_G <- contig_paired_only[contig_paired_only$pairing %in% c("AB GD", "AB G", "B GD", "B G"), ]

        grep("junction_aa", (names(AB_G)))

        junction_only <- contig_paired_only[, c(names(contig_paired_only)[grep("junction_aa", (names(contig_paired_only)))], names(contig_paired_only)[grep("v_gene", (names(contig_paired_only)))], names(contig_paired_only)[grep("chain_", (names(contig_paired_only)))]), ]
        names.list <- names(junction_only)
        junction_only$cloneCount <- 1

        grep("^C", junction_only$junction_aa_A)

        junction_only_sum <- ddply(junction_only, names.list, numcolwise(sum))

        junction_only_sum$pairing <- ifelse(junction_only_sum$chain_B == "TRB" & junction_only_sum$chain_A == "TRA" & junction_only_sum$chain_D == "TRD" & junction_only_sum$chain_G == "TRG", "AB GD",
                                            ifelse(junction_only_sum$chain_B == "TRB" & junction_only_sum$chain_A == "TRA" & junction_only_sum$chain_D != "TRD" & junction_only_sum$chain_G == "TRG", "AB G",
                                                   ifelse(junction_only_sum$chain_B == "TRB" & junction_only_sum$chain_A == "TRA" & junction_only_sum$chain_D == "TRD" & junction_only_sum$chain_G != "TRG", "AB D",
                                                          ifelse(junction_only_sum$chain_B != "TRB" & junction_only_sum$chain_A == "TRA" & junction_only_sum$chain_D == "TRD" & junction_only_sum$chain_G == "TRG", "AG D",
                                                                 ifelse(junction_only_sum$chain_B == "TRB" & junction_only_sum$chain_A != "TRA" & junction_only_sum$chain_D == "TRD" & junction_only_sum$chain_G == "TRG", "B GD",
                                                                        ifelse(junction_only_sum$chain_B == "TRB" & junction_only_sum$chain_A == "TRA" & junction_only_sum$chain_D != "TRD" & junction_only_sum$chain_G != "TRG", "AB",
                                                                               ifelse(junction_only_sum$chain_B != "TRB" & junction_only_sum$chain_A != "TRA" & junction_only_sum$chain_D == "TRD" & junction_only_sum$chain_G == "TRG", "GD",
                                                                                      ifelse(junction_only_sum$chain_B == "TRB" & junction_only_sum$chain_A != "TRA" & junction_only_sum$chain_D != "TRD" & junction_only_sum$chain_G == "TRG", "B G",
                                                                                             ifelse(junction_only_sum$chain_B == "TRB" & junction_only_sum$chain_A != "TRA" & junction_only_sum$chain_D == "TRD" & junction_only_sum$chain_G != "TRG", "B D",
                                                                                                    ifelse(junction_only_sum$chain_B != "TRB" & junction_only_sum$chain_A == "TRA" & junction_only_sum$chain_D == "TRD" & junction_only_sum$chain_G != "TRG", "A D",
                                                                                                           ifelse(junction_only_sum$chain_B != "TRB" & junction_only_sum$chain_A == "TRA" & junction_only_sum$chain_D != "TRD" & junction_only_sum$chain_G == "TRG", "A G",
                                                                                                                  ifelse(junction_only_sum$chain_B != "TRB" & junction_only_sum$chain_A == "TRA" & junction_only_sum$chain_D != "TRD" & junction_only_sum$chain_G != "TRG", "A",
                                                                                                                         ifelse(junction_only_sum$chain_B == "TRB" & junction_only_sum$chain_A != "TRA" & junction_only_sum$chain_D != "TRD" & junction_only_sum$chain_G != "TRG", "B",
                                                                                                                                ifelse(junction_only_sum$chain_B != "TRB" & junction_only_sum$chain_A != "TRA" & junction_only_sum$chain_D != "TRD" & junction_only_sum$chain_G == "TRG", "G",
                                                                                                                                       ifelse(junction_only_sum$chain_B != "TRB" & junction_only_sum$chain_A != "TRA" & junction_only_sum$chain_D == "TRD" & junction_only_sum$chain_G != "TRG", "D",
                                                                                                                                              "other"
                                                                                                                                       )
                                                                                                                                )
                                                                                                                         )
                                                                                                                  )
                                                                                                           )
                                                                                                    )
                                                                                             )
                                                                                      )
                                                                               )
                                                                        )
                                                                 )
                                                          )
                                                   )
                                            )
        )

        head(junction_only_sum)
        junction_only_sum_selected <- subset(junction_only_sum, junction_only_sum$junction_aa_A == "CVLSASSSFSKLVF")
        head(junction_only_sum_selected)
        merged_junction_contig <- merge(contig_paired_only, junction_only_sum, by = c("junction_aa_A", "junction_aa_B", "junction_aa_G", "junction_aa_D", "v_gene_A", "v_gene_B", "v_gene_G", "v_gene_D", "chain_A", "chain_B", "chain_G", "chain_D", "pairing")) # "v_gene_B","v_gene_D","v_gene_G"
        name.list_all <- c("cell_id", "Sample_Name", "pairing", "order")
        merged_junction_contig <- merged_junction_contig %>%
          select(all_of(name.list_all), everything())

        # merged_junction_contig <- merged_junction_contig[!names(merged_junction_contig) %in% c(names(merged_junction_contig[grep("v_gene",names(merged_junction_contig))]))]
        merged_junction_contig <- merged_junction_contig[order(merged_junction_contig$cell_id), ]
        merged_junction_contig
      } else {
        as.data.frame("Upload unfilted AIRR file")
      }
    })

    output$tb_multiTCR <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      multi_TCR_BDrhap()
    })

    output$downloadtb_multiTCR <- downloadHandler(
      filename = function() {
        if (input$name_BD == "") {
          showModal(modalDialog(
            title = "Error",
            "Please Add File Name.",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          x <- today()
          paste(input$name_BD, "_Multi_TCR_", x, ".csv", sep = "")
        }
      },
      content = function(file) {
        df <- as.data.frame(multi_TCR_BDrhap())
        write.csv(df, file, row.names = F)
      }
    )

    # 10x_Genomics data -----
    ## barcode file -----
    input.data.barcode.10x <- reactive({
      inFile_10x_barcode <- input$file_calls_10x
      if (is.null(inFile_10x_barcode)) {
        return(NULL)
      } else {
        dataframe <- read.table(
          inFile_10x_barcode$datapath
        )
      }
    })

    output$Barcode_10x_file_upload <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE, scrollY = T), {
      calls <- input.data.barcode.10x()
      validate(
        need(
          nrow(calls) > 0,
          "Upload Barcode"
        )
      )
      calls
    })
    ## features file -----
    input.data.features.10x <- reactive({
      inFile_10x_features <- input$file_features_10x
      if (is.null(inFile_10x_features)) {
        return(NULL)
      } else {
        dataframe <- read.table(
          inFile_10x_features$datapath
        )
      }
    })

    output$Features_10x_file_upload <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE, scrollY = T), {
      calls <- input.data.features.10x()
      validate(
        need(
          nrow(calls) > 0,
          "Upload Features/Gene File"
        )
      )
      calls
    })

    # Matrix file
    input.data.matrix.10x <- reactive({
      inFile_10x_matrix <- input$file_matrix_10x
      if (is.null(inFile_10x_matrix)) {
        return(NULL)
      } else {
        dataframe <- Matrix::readMM(inFile_10x_matrix$datapath)
      }
    })
    output$Matrix_10x_file_upload <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      calls <- as.data.frame(input.data.matrix.10x())
      validate(
        need(
          nrow(calls) > 0,
          "Upload Matrix file"
        )
      )
      calls[1:10, 1:10]
    })

    ## contig files ----
    input.data.TCR.10x <- reactive({
      inFile_10x_TCR <- input$file_TCR_10x
      if (is.null(inFile_10x_TCR)) {
        return(NULL)
      } else {
        if (input$csv_contig_file == "csv/csv.gz") {
          dataframe <- read.csv(inFile_10x_TCR$datapath, na.strings = c("", "NA"))
        } else {
          dataframe <- read.table(inFile_10x_TCR$datapath, sep = "\t", header = T)
        }
        dataframe
      }
    })
    output$Contig_10x_file_upload <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE, scrollY = T), {
      if (input$BCR_TCR_10x == "TCR only") {
        calls <- input.data.TCR.10x()
        validate(
          need(
            nrow(calls) > 0,
            "Upload TCR Files (Filtered contig file)"
          )
        )
        tb_10x_meta.data_TCR()
      } else {
        calls <- input.data.TCR.10x()
        validate(
          need(
            nrow(calls) > 0,
            "Upload BCR Files (Filtered contig file)"
          )
        )
        tb_10x_meta.data_BCR()
      }
    })
    ## meta.data for seurat ----
    tb_10x_for.tcr.exploreAIRR <- reactive({
      contigs <- input.data.TCR.10x()
      validate(
        need(
          nrow(contigs) > 0,
          "Upload TCR Files (Filtered contig file)"
        )
      )
      contigs
      names(contigs) <- gsub("call", "gene", names(contigs))

      names(contigs) <- gsub("junction_aa", "cdr3", names(contigs))
      contigs
      contigs_lim <- contigs[!names(contigs) %in% c("is_cell", "contig_id", "high_confidence", "raw_consensus_id", "exact_subclonotype_id", "umis", "reads", "length", "cdr3_nt", "germline_alignment", "cdr3_nt_id", "cdr3_nt_alignment", "rev_comp", names(contigs[grep("_end", names(contigs))]), names(contigs[grep("cigar", names(contigs))]), names(contigs[grep("length", names(contigs))]), names(contigs[grep("count", names(contigs))]), names(contigs[grep("sequence", names(contigs))]), names(contigs[grep("fwr", names(contigs))]), names(contigs[grep("cdr1", names(contigs))]), names(contigs[grep("cdr2", names(contigs))]))]


      contigs_lim$chain <- ifelse(grepl("TRA", contigs_lim$v_gene), "TRA",
                                  ifelse(grepl("TRB", contigs_lim$v_gene), "TRB",
                                         ifelse(grepl("TRG", contigs_lim$v_gene), "TRG",
                                                ifelse(grepl("TRD", contigs_lim$v_gene), "TRD",
                                                       ""
                                                )
                                         )
                                  )
      )

      contigs_lim
      names(contigs_lim) <- gsub("junction", "cdr3_nt", names(contigs_lim))

      contig_AG <- subset(contigs_lim, contigs_lim$chain == "TRA" | contigs_lim$chain == "TRG")
      head(contig_AG)
      name.list <- names(contig_AG[c(
        names(contig_AG[grep("gene", names(contig_AG))]),
        names(contig_AG[grep("cdr3", names(contig_AG))]),
        "chain"
      )])
      name.list

      contig_AG <- contig_AG %>%
        select(all_of(name.list), everything())
      names(contig_AG)[1:summary(name.list)[1]] <- paste(names(contig_AG[names(contig_AG) %in% name.list]), "_AG", sep = "")
      contig_AG

      contig_BD <- subset(contigs_lim, contigs_lim$chain == "TRB" | contigs_lim$chain == "TRD")
      name.list <- names(contig_BD[c(
        names(contig_BD[grep("gene", names(contig_BD))]),
        names(contig_BD[grep("cdr3", names(contig_BD))]),
        "chain"
      )])
      contig_BD <- contig_BD %>%
        select(all_of(name.list), everything())


      names(contig_BD)[1:summary(name.list)[1]] <- paste(names(contig_BD[names(contig_BD) %in% name.list]), "_BD", sep = "")
      contig_BD
      # contig_paired <- merge(contig_AG,contig_BD, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      # contig_paired <- merge(contig_AG,contig_BD, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      contig_paired <- merge(contig_AG, contig_BD, by = c("cell_id", "clone_id"), all = T)

      contig_paired$pairing <- ifelse(contig_paired$chain_BD == "TRB" & contig_paired$chain_AG == "TRA", "abTCR Paired",
                                      ifelse(contig_paired$chain_BD == "TRD" & contig_paired$chain_AG == "TRG", "gdTCR Paired", NA)
      )
      contig_paired
      contig_paired$pairing[is.na(contig_paired$pairing)] <- "unpaired"
      contig_paired <- contig_paired[!names(contig_paired) %in% c("d_gene_AG")]
      contig_paired_only <- contig_paired
      contig_paired_only <- subset(contig_paired_only, contig_paired_only$cdr3_BD != "None")
      contig_paired_only <- subset(contig_paired_only, contig_paired_only$cdr3_AG != "None")
      contig_paired_only$d_gene_BD <- gsub("^$", "NA", contig_paired_only$d_gene_BD)
      #
      contig_paired_only$vj_gene_AG <- paste(contig_paired_only$v_gene_AG, contig_paired_only$j_gene_AG, sep = ".")
      contig_paired_only$vj_gene_AG <- gsub("NA.NA", "", contig_paired_only$vj_gene_AG)
      #
      contig_paired_only$vj_gene_BD <- paste(contig_paired_only$v_gene_BD, contig_paired_only$j_gene_BD, sep = ".")
      contig_paired_only$vj_gene_BD <- gsub(".NA.", ".", contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub(".None.", ".", contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub("NA.NA", "", contig_paired_only$vj_gene_BD)
      #
      contig_paired_only$vdj_gene_BD <- paste(contig_paired_only$v_gene_BD, contig_paired_only$d_gene_BD, contig_paired_only$j_gene_BD, sep = ".")
      contig_paired_only$vdj_gene_BD <- gsub(".NA.", ".", contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("[.]NA[.]", ".", contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("NA.NA", "", contig_paired_only$vdj_gene_BD)
      #
      contig_paired_only$vj_gene_cdr3_AG <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$cdr3_AG, sep = "_")
      contig_paired_only$vj_gene_cdr3_AG <- gsub("_NA", "", contig_paired_only$vj_gene_cdr3_AG)
      #
      contig_paired_only$vj_gene_cdr3_BD <- paste(contig_paired_only$vj_gene_BD, contig_paired_only$cdr3_BD, sep = "_")
      contig_paired_only$vj_gene_cdr3_BD <- gsub("_NA", "", contig_paired_only$vj_gene_cdr3_BD)
      #
      contig_paired_only$vdj_gene_cdr3_BD <- paste(contig_paired_only$vdj_gene_BD, contig_paired_only$cdr3_BD, sep = "_")
      contig_paired_only$vdj_gene_cdr3_BD <- gsub("_NA", "", contig_paired_only$vdj_gene_cdr3_BD)
      #
      contig_paired_only$vj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$vj_gene_BD, sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$vdj_gene_BD, sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- gsub("^ & ", "", contig_paired_only$vdj_gene_AG_BD)
      contig_paired_only$vdj_gene_AG_BD <- gsub(" & $", "", contig_paired_only$vdj_gene_AG_BD)
      #
      # #updating names to be consistant....
      contig_paired_only$vj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG, contig_paired_only$vj_gene_cdr3_BD, sep = " & ")
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub("^ & ", "", contig_paired_only$vj_gene_cdr3_AG_BD)
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub(" & $", "", contig_paired_only$vj_gene_cdr3_AG_BD)

      contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG, contig_paired_only$vdj_gene_cdr3_BD, sep = " & ")
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub("^ & ", "", contig_paired_only$vdj_gene_cdr3_AG_BD)
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub(" & $", "", contig_paired_only$vdj_gene_cdr3_AG_BD)


      names(contig_paired_only)[names(contig_paired_only) %in% "cell_id"] <- "Cell_Index"

      contig_paired_only$Sample_Name <- paste(input$Indiv10x, "_", input$group10x, sep = "")

      contig_paired_only <- contig_paired_only %>%
        select(all_of(c("Cell_Index", "Sample_Name")), everything())
      contig_paired_only
    })
    tb_10x_meta.data_TCR_CSV <- function() {
      contigs <- input.data.TCR.10x()
      validate(
        need(
          nrow(contigs) > 0,
          "Upload file"
        )
      )
      contigs_lim <- contigs[!names(contigs) %in% c("is_cell", "contig_id", "high_confidence", "raw_consensus_id", "exact_subclonotype_id", "umis", "reads", "length", "cdr3_nt", names(contigs[grep("fwr", names(contigs))]), names(contigs[grep("cdr1", names(contigs))]), names(contigs[grep("cdr2", names(contigs))]))]
      contigs_lim
      contig_AG <- subset(contigs_lim, contigs_lim$chain == "TRA" | contigs_lim$chain == "TRG")
      name.list <- names(contig_AG[c(
        names(contig_AG[grep("gene", names(contig_AG))]),
        names(contig_AG[grep("cdr3", names(contig_AG))]),
        "chain"
      )])
      contig_AG <- contig_AG %>%
        select(all_of(name.list), everything())
      names(contig_AG)[1:summary(name.list)[1]] <- paste(names(contig_AG[names(contig_AG) %in% name.list]), "_AG", sep = "")
      contig_AG

      contig_BD <- subset(contigs_lim, contigs_lim$chain == "TRB" | contigs_lim$chain == "TRD")
      name.list <- names(contig_BD[c(
        names(contig_BD[grep("gene", names(contig_BD))]),
        names(contig_BD[grep("cdr3", names(contig_BD))]),
        "chain"
      )])
      contig_BD <- contig_BD %>%
        select(all_of(name.list), everything())


      names(contig_BD)[1:summary(name.list)[1]] <- paste(names(contig_BD[names(contig_BD) %in% name.list]), "_BD", sep = "")
      contig_BD
      # contig_paired <- merge(contig_AG,contig_BD, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      # contig_paired <- merge(contig_AG,contig_BD, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      contig_paired <- merge(contig_AG, contig_BD, by = c("barcode", "full_length", "productive", "raw_clonotype_id"), all = T)

      contig_paired$pairing <- ifelse(contig_paired$chain_BD == "TRB" & contig_paired$chain_AG == "TRA", "abTCR Paired",
                                      ifelse(contig_paired$chain_BD == "TRD" & contig_paired$chain_AG == "TRG", "gdTCR Paired", NA)
      )
      contig_paired
      contig_paired$pairing[is.na(contig_paired$pairing)] <- "unpaired"
      contig_paired <- contig_paired[!names(contig_paired) %in% c("d_gene_AG")]
      contig_paired_only <- contig_paired
      contig_paired_only <- subset(contig_paired_only, contig_paired_only$cdr3_BD != "None")
      contig_paired_only <- subset(contig_paired_only, contig_paired_only$cdr3_AG != "None")
      contig_paired_only$d_gene_BD <- gsub("^$", "NA", contig_paired_only$d_gene_BD)
      #
      contig_paired_only$vj_gene_AG <- paste(contig_paired_only$v_gene_AG, contig_paired_only$j_gene_AG, sep = ".")
      contig_paired_only$vj_gene_AG <- gsub("NA.NA", "", contig_paired_only$vj_gene_AG)
      #
      contig_paired_only$vj_gene_BD <- paste(contig_paired_only$v_gene_BD, contig_paired_only$j_gene_BD, sep = ".")
      contig_paired_only$vj_gene_BD <- gsub(".NA.", ".", contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub(".None.", ".", contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub("NA.NA", "", contig_paired_only$vj_gene_BD)
      #
      contig_paired_only$vdj_gene_BD <- paste(contig_paired_only$v_gene_BD, contig_paired_only$d_gene_BD, contig_paired_only$j_gene_BD, sep = ".")
      contig_paired_only$vdj_gene_BD <- gsub(".NA.", ".", contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("[.]NA[.]", ".", contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("NA.NA", "", contig_paired_only$vdj_gene_BD)
      #
      contig_paired_only$vj_gene_cdr3_AG <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$cdr3_AG, sep = "_")
      contig_paired_only$vj_gene_cdr3_AG <- gsub("_NA", "", contig_paired_only$vj_gene_cdr3_AG)
      #
      contig_paired_only$vj_gene_cdr3_BD <- paste(contig_paired_only$vj_gene_BD, contig_paired_only$cdr3_BD, sep = "_")
      contig_paired_only$vj_gene_cdr3_BD <- gsub("_NA", "", contig_paired_only$vj_gene_cdr3_BD)
      #
      contig_paired_only$vdj_gene_cdr3_BD <- paste(contig_paired_only$vdj_gene_BD, contig_paired_only$cdr3_BD, sep = "_")
      contig_paired_only$vdj_gene_cdr3_BD <- gsub("_NA", "", contig_paired_only$vdj_gene_cdr3_BD)
      #
      contig_paired_only$vj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$vj_gene_BD, sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$vdj_gene_BD, sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- gsub("^ & ", "", contig_paired_only$vdj_gene_AG_BD)
      contig_paired_only$vdj_gene_AG_BD <- gsub(" & $", "", contig_paired_only$vdj_gene_AG_BD)
      #
      # #updating names to be consistant....
      contig_paired_only$vj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG, contig_paired_only$vj_gene_cdr3_BD, sep = " & ")
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub("^ & ", "", contig_paired_only$vj_gene_cdr3_AG_BD)
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub(" & $", "", contig_paired_only$vj_gene_cdr3_AG_BD)

      contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG, contig_paired_only$vdj_gene_cdr3_BD, sep = " & ")
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub("^ & ", "", contig_paired_only$vdj_gene_cdr3_AG_BD)
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub(" & $", "", contig_paired_only$vdj_gene_cdr3_AG_BD)
      # contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vdj_gene_cdr3_BD,sep = " & ")
      names(contig_paired_only)[names(contig_paired_only) %in% "barcode"] <- "Cell_Index"
      contig_paired_only <- contig_paired_only[!duplicated(contig_paired_only$Cell_Index), ] # remove duplicates
      contig_paired_only$Sample_Name <- paste(input$Indiv10x, "_", input$group10x, sep = "")

      contig_paired_only <- contig_paired_only %>%
        select(all_of(c("Cell_Index", "Sample_Name")), everything())

      contig_paired_only
    }
    tb_10x_meta.data_TCR <- reactive({
      if (input$csv_contig_file == "csv/csv.gz") {
        dataframe <- tb_10x_meta.data_TCR_CSV()
      } else {
        dataframe <- tb_10x_for.tcr.exploreAIRR()
      }
      dataframe
    })


    tb_10x_meta.data_BCR <- function() {
      contigs <- input.data.TCR.10x()
      validate(
        need(
          nrow(contigs) > 0,
          "Upload file"
        )
      )
      contigs <- contigs[order(contigs$umis, decreasing = T), ]

      contigs_lim <- contigs[!names(contigs) %in% c("is_cell", "contig_id", "high_confidence", "raw_consensus_id", "exact_subclonotype_id", "reads", "length", "cdr3_nt", names(contigs[grep("fwr", names(contigs))]), names(contigs[grep("cdr1", names(contigs))]), names(contigs[grep("cdr2", names(contigs))]))]
      contigs_lim
      contig_LK <- subset(contigs_lim, contigs_lim$chain == "IGL" | contigs_lim$chain == "IGK")
      name.list <- names(contig_LK[c(
        names(contig_LK[grep("gene", names(contig_LK))]),
        names(contig_LK[grep("cdr3", names(contig_LK))]),
        "chain"
      )])
      name.list
      contig_LK <- contig_LK %>%
        select(all_of(name.list), everything())
      names(contig_LK)[1:summary(name.list)[1]] <- paste(names(contig_LK[names(contig_LK) %in% name.list]), "_IgL", sep = "")
      head(contig_LK)

      contig_H <- subset(contigs_lim, contigs_lim$chain == "IGH")

      name.list <- names(contig_H[c(
        names(contig_H[grep("gene", names(contig_H))]),
        names(contig_H[grep("cdr3", names(contig_H))]),
        "chain"
      )])
      contig_H <- contig_H %>%
        select(all_of(name.list), everything())


      names(contig_H)[1:summary(name.list)[1]] <- paste(names(contig_H[names(contig_H) %in% name.list]), "_IgH", sep = "")
      head(contig_H)
      # contig_paired <- merge(contig_LK,contig_H, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      # contig_paired <- merge(contig_LK,contig_H, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
      contig_paired <- merge(contig_LK, contig_H, by = c("barcode", "full_length", "productive", "raw_clonotype_id"), all = T)

      contig_paired$pairing <- ifelse(contig_paired$chain_IgH == "IGH" & contig_paired$chain_IgL == "IGK", "IGK Paired",
                                      ifelse(contig_paired$chain_IgH == "IGH" & contig_paired$chain_IgL == "IGL", "IGL Paired", NA)
      )

      contig_paired
      contig_paired$pairing[is.na(contig_paired$pairing)] <- "unpaired"
      contig_paired <- contig_paired[!names(contig_paired) %in% c("d_gene_IgL")]
      contig_paired_only <- contig_paired
      contig_paired_only <- subset(contig_paired_only, contig_paired_only$cdr3_IgH != "None")
      contig_paired_only <- subset(contig_paired_only, contig_paired_only$cdr3_IgL != "None")
      dim(contig_paired_only)

      contig_paired_only$d_gene_IgH <- sub("^$", "NA", contig_paired_only$d_gene_IgH)
      #
      contig_paired_only$vj_gene_LK <- paste(contig_paired_only$v_gene_IgL, contig_paired_only$j_gene_IgL, sep = ".")
      contig_paired_only$vj_gene_LK <- gsub("NA.NA", "", contig_paired_only$vj_gene_LK)
      #
      contig_paired_only$vj_gene_H <- paste(contig_paired_only$v_gene_IgH, contig_paired_only$j_gene_IgH, sep = ".")
      contig_paired_only$vj_gene_H <- gsub(".NA.", ".", contig_paired_only$vj_gene_H)
      contig_paired_only$vj_gene_H <- gsub(".None.", ".", contig_paired_only$vj_gene_H)
      contig_paired_only$vj_gene_H <- gsub("NA.NA", "", contig_paired_only$vj_gene_H)
      #
      contig_paired_only$vdj_gene_H <- paste(contig_paired_only$v_gene_IgH, contig_paired_only$d_gene_IgH, contig_paired_only$j_gene_IgH, sep = ".")
      contig_paired_only$vdj_gene_H <- gsub(".NA.", ".", contig_paired_only$vdj_gene_H)
      contig_paired_only$vdj_gene_H <- gsub(".None.", ".", contig_paired_only$vdj_gene_H)
      contig_paired_only$vdj_gene_H <- gsub("NA.NA", "", contig_paired_only$vdj_gene_H)
      #
      contig_paired_only$vj_gene_cdr3_LK <- paste(contig_paired_only$vj_gene_LK, contig_paired_only$cdr3_IgL, sep = "_")
      contig_paired_only$vj_gene_cdr3_LK <- gsub("_NA", "", contig_paired_only$vj_gene_cdr3_LK)
      #
      contig_paired_only$vj_gene_cdr3_H <- paste(contig_paired_only$vj_gene_H, contig_paired_only$cdr3_IgH, sep = "_")
      contig_paired_only$vj_gene_cdr3_H <- gsub("_NA", "", contig_paired_only$vj_gene_cdr3_H)
      #
      contig_paired_only$vdj_gene_cdr3_H <- paste(contig_paired_only$vdj_gene_H, contig_paired_only$cdr3_IgH, sep = "_")
      contig_paired_only$vdj_gene_cdr3_H <- gsub("_NA", "", contig_paired_only$vdj_gene_cdr3_H)
      #
      contig_paired_only$vj_gene_LK_H <- paste(contig_paired_only$vj_gene_LK, contig_paired_only$vj_gene_H, sep = " & ")
      contig_paired_only$vdj_gene_LK_H <- paste(contig_paired_only$vj_gene_LK, contig_paired_only$vdj_gene_H, sep = " & ")
      contig_paired_only$vdj_gene_LK_H <- gsub("^ & ", "", contig_paired_only$vdj_gene_LK_H)
      contig_paired_only$vdj_gene_LK_H <- gsub(" & $", "", contig_paired_only$vdj_gene_LK_H)
      #
      # #updating names to be consistant....
      contig_paired_only$vj_gene_cdr3_LK_H <- paste(contig_paired_only$vj_gene_cdr3_LK, contig_paired_only$vj_gene_cdr3_H, sep = " & ")
      contig_paired_only$vj_gene_cdr3_LK_H <- gsub("^ & ", "", contig_paired_only$vj_gene_cdr3_LK_H)
      contig_paired_only$vj_gene_cdr3_LK_H <- gsub(" & $", "", contig_paired_only$vj_gene_cdr3_LK_H)

      contig_paired_only$vdj_gene_cdr3_LK_H <- paste(contig_paired_only$vj_gene_cdr3_LK, contig_paired_only$vdj_gene_cdr3_H, sep = " & ")
      contig_paired_only$vdj_gene_cdr3_LK_H <- gsub("^ & ", "", contig_paired_only$vdj_gene_cdr3_LK_H)
      contig_paired_only$vdj_gene_cdr3_LK_H <- gsub(" & $", "", contig_paired_only$vdj_gene_cdr3_LK_H)
      # contig_paired_only$vdj_gene_cdr3_LK_H <- paste(contig_paired_only$vj_gene_cdr3_LK,contig_paired_only$vdj_gene_cdr3_H,sep = " & ")
      names(contig_paired_only)[names(contig_paired_only) %in% "barcode"] <- "Cell_Index"
      dup <- contig_paired_only[duplicated(contig_paired_only$Cell_Index), ]
      contig_paired_only <- contig_paired_only[order(contig_paired_only$Cell_Index, contig_paired_only$umis.x, contig_paired_only$umis.y, decreasing = T), ]

      contig_paired_only_dup <- contig_paired_only[!duplicated(contig_paired_only$Cell_Index), ] # remove duplicate barcodes.
      names(contig_paired_only_dup)
      contig_paired_only_dup <- contig_paired_only_dup[!names(contig_paired_only_dup) %in% c("umis.x", "umis.y")]
      contig_paired_only_dup$Sample_Name <- paste(input$Indiv10x, "_", input$group10x, sep = "")

      contig_paired_only_dup <- contig_paired_only_dup %>%
        select(all_of(c("Cell_Index", "Sample_Name")), everything())

      contig_paired_only_dup

      contig_paired_only_dup
    }

    output$tb_10x_meta1 <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 5, scrollX = TRUE), {
      if (input$BCR_TCR_10x == "TCR only") {
        tb_10x_meta.data_TCR()
      } else if (input$BCR_TCR_10x == "BCR only") {
        tb_10x_meta.data_BCR()
      } else {
        # TCR <- tb_10x_meta.data_TCR()
        # BCR <- tb_10x_meta.data_BCR()
        # contig_paired <- names(BCR)[!grepl("_IgH",names(BCR)) & !grepl("_IgL",names(BCR))]
        # contig_paired <- merge(TCR,BCR,by = merge.names, all=T)
        # contig_paired
      }
    })
    output$downloadtb_10x_metadata2 <- downloadHandler(
      filename = function() {
        if (input$group10x == "" && input$Indiv10x == "") {
          showModal(modalDialog(
            title = "Error",
            "Add Group and Individual Name",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          x <- today()
          paste( input$Indiv10x, "_", input$group10x, "_metadata_10x_", x, ".csv", sep = "")
        }

      },
      content = function(file) {
        if (input$BCR_TCR_10x == "TCR only") {
          df <- tb_10x_meta.data_TCR()
        } else if (input$BCR_TCR_10x == "BCR only") {
          df <- tb_10x_meta.data_BCR()
        } else {

        }
        # write.table(,file, row.names = T)
        write.csv(df, file, row.names = F)
      }
    )

    # meta-data summary ----
    sum_tb_10x <- function() {
      contigs <- input.data.TCR.10x()
      validate(
        need(
          nrow(contigs) > 0,
          "Upload file"
        )
      )

      if (input$BCR_TCR_10x == "TCR only") {
        contig_paired <- tb_10x_meta.data_TCR()
      } else if (input$BCR_TCR_10x == "BCR only") {
        contig_paired <- tb_10x_meta.data_BCR()
      } else {
        TCR <- tb_10x_meta.data_TCR()
        BCR <- tb_10x_meta.data_BCR()
        contig_paired <- names(BCR)[!grepl("_IgH", names(BCR)) & !grepl("_IgL", names(BCR))]
        contig_paired <- merge(TCR, BCR, by = merge.names, all = T)
      }

      # contig_paired <- as.data.frame(tb_10x_meta.data())

      count.df <- contig_paired[names(contig_paired) %in% c(names(contig_paired)[grepl("chain", names(contig_paired))], "pairing")]
      count.df$count <- 1
      # ddply(count.df,names(count.df) ,numcolwise(sum))
      df1 <- ddply(count.df, names(count.df)[c(-4)], numcolwise(sum))
      df1
    }
    output$sum_tb_10x1 <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      sum_tb_10x()
    })

    ## ClusTCR T cell receptor sequences  -----
    tb_10x_contigues_contig <- reactive({
      contigs <- input.data.TCR.10x()
      validate(
        need(
          nrow(contigs) > 0,
          "Upload TCR Files (Filtered contig file)"
        )
      )
      if (input$csv_contig_file == "tsv") {
        names(contigs) <- gsub("call", "gene", names(contigs))
        names(contigs) <- gsub("junction_aa", "cdr3", names(contigs))
      } else {

      }

      contigs2 <- contigs[names(contigs) %in% c("v_gene", "cdr3")]
      contigs2$Sample_Name <- paste(  input$Indiv10x, "_", input$group10x,sep = "")
      names(contigs2)[names(contigs2) %in% c("cdr3")] <- "junction_aa"
      names(contigs2)[names(contigs2) %in% c("v_gene")] <- "v_call"
      if (nrow(contigs2[-c(grep("[*]", contigs2$junction_aa)), ] > 0)) {
        contigs2 <- contigs2[-c(grep("[*]", contigs2$junction_aa)), ]
      }

      contigs2$junction_aa <- gsub("^$", "None", contigs2$junction_aa)

      if (nrow(contigs2[-c(grep("None", contigs2$junction_aa)), ] > 0)) {
        contigs2 <- contigs2[-c(grep("None", contigs2$junction_aa)), ]
      }

      df <- contigs2[!duplicated(contigs2[, c("v_call", "junction_aa")]), ]

      if (input$chain_clusTCR2_10x == "AG") {
        (rbind(df[grep("TRAV", df$v_call), ], df[grep("TRGV", df$v_call), ]))
      } else if (input$chain_clusTCR2_10x == "IgH") {
        rbind(df[grep("IGH", df$v_call), ])
      } else if (input$chain_clusTCR2_10x == "IgLK") {
        (rbind(df[grep("IGL", df$v_call), ], df[grep("IGK", df$v_call), ]))
      } else {
        (rbind(df[grep("TRBV", df$v_call), ], df[grep("TRDV", df$v_call), ]))
      }
    })

    output$tb_10x_contigues1 <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      tb_10x_contigues_contig()
    })

    output$downloadtb_10x_contigues1 <- downloadHandler(
      filename = function() {
        if (input$group10x == "" && input$Indiv10x == "") {
          showModal(modalDialog(
            title = "Error",
            "Add Group and Individual Name",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          x <- today()
          paste(input$chain_clusTCR2_10x, "_", input$Indiv10x, "_", input$group10x,"_clusTCR_10x_", x, ".csv", sep = "")
        }

      },
      content = function(file) {
        df <- as.data.frame(tb_10x_contigues_contig())
        write.csv(df, file, row.names = F)
      }
    )


    ## TCR explore 10x -----
    TCR_Explore_10x <- function() {
      if (input$BCR_TCR_10x == "TCR only") {
        if (input$csv_contig_file == "csv/csv.gz") {
          contig_paired_only <- tb_10x_meta.data_TCR()
        } else {
          contig_paired_only <- tb_10x_for.tcr.exploreAIRR()
        }
      } else if (input$BCR_TCR_10x == "BCR only") {
        contig_paired_only <- tb_10x_meta.data_BCR()
      } else {
        TCR <- tb_10x_meta.data_TCR()
        BCR <- tb_10x_meta.data_BCR()
        contig_paired_only <- names(BCR)[!grepl("_IgH", names(BCR)) & !grepl("_IgL", names(BCR))]
        contig_paired_only <- merge(TCR, BCR, by = merge.names, all = T)
      }

      contig_paired_only$cloneCount <- 1
      contig_paired_only$group <- input$group10x
      contig_paired_only$Indiv <- input$Indiv10x
      contig_paired_only <- contig_paired_only %>%
        select(cloneCount, group, Indiv, everything())
      contig_paired_only
    }

    output$dt_TCR_Explore_10x <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      TCR_Explore_10x()
    })

    output$downloaddt_TCR_Explore_10x <- downloadHandler(
      filename = function() {
        if (input$group10x == "" && input$Indiv10x == "") {
          showModal(modalDialog(
            title = "Error",
            "Add Group and Individual Name",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          x <- today()
          paste(input$Indiv10x, "_", input$group10x , "_TCR_Explore_10x_", x, ".csv", sep = "")
        }

      },
      content = function(file) {
        df <- as.data.frame(TCR_Explore_10x())
        write.csv(df, file, row.names = F)
      }
    )



    ## count table ----
    tb_10x_matrix <- function() {
      barcode <- input.data.barcode.10x()
      features <- input.data.features.10x()
      mmMat <- as.data.frame(input.data.matrix.10x())

      validate(
        need(
          nrow(barcode) > 0 & nrow(features) > 0 & nrow(mmMat) > 0,
          "Upload files"
        )
      )

      rownames(mmMat) <- make.unique(features$V2) # select which feature to label genes...
      names(mmMat) <- barcode$V1
      mmMat$Gene_Name <- rownames(mmMat)
      mmMat <- mmMat %>%
        select(all_of("Gene_Name"), everything())
      mmMat
    }

    tb_10x_matrix_h5 <- function() {
      mmMat <- as.data.frame(input.data.h5.10x())

      validate(
        need(
          nrow(mmMat) > 0,
          "Upload h5 files"
        )
      )
      mmMat
    }

    output$tb_10x_matrix2 <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      tb_10x_matrix()[1:6, 1:6]
    })

    output$downloadtb_10x_matrix2 <- downloadHandler(
      filename = function() {
        if (input$group10x == "" && input$Indiv10x == "") {
          showModal(modalDialog(
            title = "Error",
            "Add Group and Individual Name",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          x <- today()
          paste(input$Indiv10x, "_", input$group10x, "_count-matrix_10x_", x, ".csv.gz", sep = "")
        }

      },
      content = function(file) {
        df <- as.data.frame(tb_10x_matrix())
        # write.table(,file, row.names = T)
        write_csv(df, gzfile(file))
      }
    )

    # for TCRex output -----
    TCRex_10x_df <- function() {
      contigs <- input.data.TCR.10x()
      validate(
        need(
          nrow(contigs) > 0,
          "Upload TCR Files (Filtered contig file)"
        )
      )

      if (input$csv_contig_file == "tsv") {
        names(contigs) <- gsub("call", "gene", names(contigs))
        names(contigs) <- gsub("junction_aa", "cdr3", names(contigs))
      } else {

      }

      contigs2 <- contigs[, names(contigs) %in% c("v_gene", "j_gene", "cdr3")]

      names(contigs2)[names(contigs2) %in% c("cdr3")] <- "CDR3_beta"
      if (nrow(contigs2[-c(grep("[*]", contigs2$CDR3_beta)), ] > 0)) {
        contigs2 <- contigs2[-c(grep("[*]", contigs2$CDR3_beta)), ]
      }
      #

      names(contigs2)[names(contigs2) %in% c("v_gene")] <- "TRBV_gene"
      names(contigs2)[names(contigs2) %in% c("j_gene")] <- "TRBJ_gene"
      contigs2 <- contigs2 %>%
        select(TRBV_gene, CDR3_beta, everything())


      contigs2$CDR3_beta <- gsub("^$", "None", contigs2$CDR3_beta)

      if (nrow(contigs2[-c(grep("None", contigs2$CDR3_beta)), ] > 0)) {
        contigs2 <- subset(contigs2, contigs2$CDR3_beta != "None")
      }

      contigs2$TRBJ_gene <- gsub("^$", "None", contigs2$TRBJ_gene)

      if (nrow(contigs2[-c(grep("None", contigs2$TRBJ_gene)), ] > 0)) {
        contigs2 <- subset(contigs2, contigs2$TRBJ_gene != "None")
      }


      contigs2 <- contigs2
      contigs2$TRBV_gene <- gsub("[*]0.", "", contigs2$TRBV_gene)
      contigs2$TRBJ_gene <- gsub("[*]0.", "", contigs2$TRBJ_gene)
      contigs2$cloneCount <- 1
      calls_TCR_paired.fun3 <- ddply(contigs2, names(contigs2)[-c(4)], numcolwise(sum))
      calls_TCR_paired.fun3 <- calls_TCR_paired.fun3[grepl("TRB", calls_TCR_paired.fun3$TRBV_gene), ]
      calls_TCR_paired.fun3
    }

    output$tb_TCRex_10x_df <- DT::renderDT(filter = list(position = "top", clear = FALSE), escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      TCRex_10x_df()
    })

    output$downloaddf_TCRex_10x <- downloadHandler(
      filename = function() {
        if (input$group10x == "" && input$Indiv10x == "") {
          showModal(modalDialog(
            title = "Error",
            "Add Group and Individual Name",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          x <- today()
          paste(input$Indiv10x, "_", input$group10x, "_TCRex_", x, ".tsv", sep = "")
        }
      },
      content = function(file) {
        df <- as.data.frame(TCRex_10x_df())
        df <- df[, names(df) %in% c("TRBV_gene", "CDR3_beta", "TRBJ_gene")]
        write.table(df, file, row.names = F, sep = "\t", quote = F)
      }
    )

    # automated 10x process if the folder is of the correct structure

    checkFiles_dt <- reactive({
      req(input$checkFiles)

      if (input$checkFiles) {
        FN <- tempfile()
        zz <- file(FN, open = "wt")
        sink(zz, type = "output")
        sink(zz, type = "message")
        STEGO.R::preprocessing_10x(downloadTCRex = FALSE, downloadClusTCR = FALSE, downloadSeurat = FALSE, downloadTCR_Explore = FALSE)
        sink(type = "message")
        sink(type = "output")
        cat(readLines(FN), sep = "\n")
      } else {

      }
    })

    output$outputText <- renderPrint({
      checkFiles_dt()

    })

    observeEvent(input$automateProcess, {

      withProgress(message = "Pre processing underway", value = 0,
                   {
                     suppressWarnings(
                       preprocessing_10x(downloadTCRex = T, downloadClusTCR = T, downloadSeurat = T, downloadTCR_Explore = T)
                     )
                   }
      )
    })

    # Array ------
    input.data.calls.array <- reactive({
      inFile_arrayM <- input$file_calls_Array
      if (is.null(inFile_arrayM)) {
        return(NULL)
      } else {
        dataframe <- read.table(inFile_arrayM$datapath)
      }
    })
    output$test.files_array_Matrix <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      calls <- input.data.calls.array()
      validate(
        need(
          nrow(calls) > 0,
          "Upload Matrix file"
        )
      )
      head(calls)[1:6]
    })

    input.data.calls.array_contig <- reactive({
      inFile_arrayC <- input$file_contig_Array
      if (is.null(inFile_arrayC)) {
        return(NULL)
      } else {
        dataframe <- read.table(inFile_arrayC$datapath, sep = "\t", header = T, row.names = 1)
      }
    })
    output$test.files_array_contig <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      calls <- input.data.calls.array_contig()
      validate(
        need(
          nrow(calls) > 0,
          "Upload Contig/TCR file"
        )
      )
      calls
    })


    Filtering_TCR_Array <- reactive({
      df <- input.data.calls.array_contig()
      df$Cell_Index <- rownames(df)
      df_filtering <- subset(df, df$clonotype != "None")

      if (input$pairing_TCR_Array == T) {
        # "Single pair"     "None"            "Orphan beta"     "Extra alpha"     "Extra beta"      "Orphan alpha"    "Two full chains"
        df_filtering <- subset(df_filtering, df_filtering$chain_pairing != "Orphan beta")
        df_filtering <- subset(df_filtering, df_filtering$chain_pairing != "Orphan alpha")
      } else {
        df_filtering
      }
      df_filtering <- df_filtering[, !grepl("_2_", names(df_filtering))]
      df_filtering <- df_filtering[, !grepl("_junction_ins", names(df_filtering))]
      df_filtering <- df_filtering[, !grepl("_expr", names(df_filtering))]
      df_filtering <- df_filtering[, !grepl("clon", names(df_filtering))]
      df_filtering <- df_filtering[, !grepl("multi_chain", names(df_filtering))]
      df_filtering <- df_filtering[, !grepl("TRA_1_d_gene", names(df_filtering))]
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

      contig_paired_only$vj_gene_AG <- paste(contig_paired_only$v_gene_AG, contig_paired_only$j_gene_AG, sep = ".")
      contig_paired_only$vj_gene_AG <- gsub("NA.NA", "", contig_paired_only$vj_gene_AG)
      contig_paired_only$vj_gene_AG <- gsub("None.None", "", contig_paired_only$vj_gene_AG)
      #
      contig_paired_only$vj_gene_BD <- paste(contig_paired_only$v_gene_BD, contig_paired_only$j_gene_BD, sep = ".")
      contig_paired_only$vj_gene_BD <- gsub(".NA.", ".", contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub(".None.", ".", contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub("NA.NA", "", contig_paired_only$vj_gene_BD)
      #
      contig_paired_only$vdj_gene_BD <- paste(contig_paired_only$v_gene_BD, contig_paired_only$d_gene_BD, contig_paired_only$j_gene_BD, sep = ".")
      contig_paired_only$vdj_gene_BD <- gsub(".NA.", ".", contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub(".None.", ".", contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("NA.NA", "", contig_paired_only$vdj_gene_BD)
      #
      contig_paired_only$vj_gene_cdr3_AG <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$cdr3_AG, sep = "_")
      contig_paired_only$vj_gene_cdr3_AG <- gsub("_NA", "", contig_paired_only$vj_gene_cdr3_AG)
      contig_paired_only$vj_gene_cdr3_AG <- gsub("_None", "", contig_paired_only$vj_gene_cdr3_AG)
      #
      contig_paired_only$vj_gene_cdr3_BD <- paste(contig_paired_only$vj_gene_BD, contig_paired_only$cdr3_BD, sep = "_")
      contig_paired_only$vj_gene_cdr3_BD <- gsub("_NA", "", contig_paired_only$vj_gene_cdr3_BD)
      #
      contig_paired_only$vdj_gene_cdr3_BD <- paste(contig_paired_only$vdj_gene_BD, contig_paired_only$cdr3_BD, sep = "_")
      contig_paired_only$vdj_gene_cdr3_BD <- gsub("_NA", "", contig_paired_only$vdj_gene_cdr3_BD)
      #
      contig_paired_only$vj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$vj_gene_BD, sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$vdj_gene_BD, sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- gsub("^ & ", "", contig_paired_only$vdj_gene_AG_BD)
      contig_paired_only$vdj_gene_AG_BD <- gsub(" & $", "", contig_paired_only$vdj_gene_AG_BD)
      #
      # #updating names to be consistant....
      contig_paired_only$vj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG, contig_paired_only$vj_gene_cdr3_BD, sep = " & ")
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub("^& ", "", contig_paired_only$vj_gene_cdr3_AG_BD)
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub(" & $", "", contig_paired_only$vj_gene_cdr3_AG_BD)

      contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG, contig_paired_only$vdj_gene_cdr3_BD, sep = " & ")
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub("^ & ", "", contig_paired_only$vdj_gene_cdr3_AG_BD)
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub(" & $", "", contig_paired_only$vdj_gene_cdr3_AG_BD)
      # contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG,contig_paired_only$vdj_gene_cdr3_BD,sep = " & ")
      names(contig_paired_only)[names(contig_paired_only) %in% "barcode"] <- "Cell_Index"
      contig_paired_only$Sample_Name <- input$sample_name_Array

      contig_paired_only <- contig_paired_only %>%
        select(all_of(c("Cell_Index", "Sample_Name")), everything())

      contig_paired_only
    })


    output$test.files_array_contig_Filtered <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      calls <- Filtering_TCR_Array()
      validate(
        need(
          nrow(calls) > 0,
          "Upload Contig TCR file"
        )
      )
      calls
    })

    # clusTCR array -----
    ClusTCR2_array <- reactive({
      df <- Filtering_TCR_Array()
      df_nocouts2_AG <- df[, names(df) %in% c("Sample_Name", "cdr3_AG", "v_gene_AG")]
      df_nocouts2_BD <- df[, names(df) %in% c("Sample_Name", "cdr3_BD", "v_gene_BD")]
      names(df_nocouts2_AG) <- c("Sample_Name", "junction_aa", "v_call")
      names(df_nocouts2_BD) <- c("Sample_Name", "junction_aa", "v_call")



      df_nocouts3 <- as.data.frame(rbind(df_nocouts2_AG, df_nocouts2_BD))
      if (nrow(df_nocouts3[-c(grep("[*]", df_nocouts3$junction_aa)), ] > 0)) {
        df_nocouts3 <- df_nocouts3[-c(grep("[*]", df_nocouts3$junction_aa)), ]
      }
      df_nocouts3 <- subset(df_nocouts3, df_nocouts3$v_call != "None")

    })

    output$test.files_ClusTCR2_array <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      calls <- ClusTCR2_array()
      calls
    })
    output$download_ClusTCR2_labs_array <- downloadHandler(
      filename = function() {

        paste("ClusTCR2_output_Array", gsub("-", ".", Sys.Date()), ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(ClusTCR2_array())
        write.csv(df, file, row.names = F)
      }
    )
    # TCRex array  -----

    ## seurat count table ----
    tb_Array_matrix <- function() {
      mmMat <- as.data.frame(input.data.calls.array())

      validate(
        need(
          nrow(mmMat) > 0,
          "Upload files"
        )
      )
      mmMat
    }

    output$tb_array_matrix2 <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      tb_Array_matrix()[1:6, 1:6]
    })

    output$downloadtb_array_matrix2 <- downloadHandler(
      filename = function() {
        if (input$name.array == "") {
          showModal(modalDialog(
            title = "Error",
            "Add Group and Individual Name",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          x <- today()
          paste(input$name.array, "_count-matrix_array_", x, ".csv", sep = "")
        }


      },
      content = function(file) {
        df <- as.data.frame(t(tb_Array_matrix()))
        # write.table(,file, row.names = T)
        write.csv(df, file)
      }
    )

    ## meta.data for seurat ----
    tb_Array_meta.data_TCR <- function() {
      contigs <- Filtering_TCR_Array()
      validate(
        need(
          nrow(contigs) > 0,
          "Upload file"
        )
      )
      contigs$Cell_Index <- rownames(contigs)
      contigs$Cell_Index <- gsub("[:]", "-", contigs$Cell_Index)
      contigs$Cell_Index <- gsub("[.]", "-", contigs$Cell_Index)
      contigs$Sample_Name <- input$sample_name_Array

      contig_paired_only <- contigs %>%
        select(all_of(c("Cell_Index", "Sample_Name")), everything())

      contig_paired_only
    }

    output$tb_Array_meta1 <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 5, scrollX = TRUE), {
      tb_Array_meta.data_TCR()
    })
    output$downloadtb_array_metadata2 <- downloadHandler(
      filename = function() {
        paste(input$name.array, "_metadata_array_", gsub("-", ".", Sys.Date()), ".csv", sep = "")
      },
      content = function(file) {
        df <- tb_Array_meta.data_TCR()
        write.csv(df, file, row.names = F)
      }
    )

    # TCRex Merge  ------
    input.data_TCRexMerge <- reactive({
      inFile2_TCRexMerge <- input$file2_TCRexMerge
      validate(
        need(
          nrow(inFile2_TCRexMerge) > 0,
          "Upload mutliple TCRex files to merge"
        )
      )
      num <- dim(inFile2_TCRexMerge)[1]
      samples_list <- vector("list", length = num)
      samples_list
      for (i in 1:num) {
        sc <- read.table(input$file2_TCRexMerge[[i, "datapath"]], sep = "\t", header = T)
        samples_list[[i]] <- sc
      }

      samples_list
    })

    merged_TCRexFiltered <- reactive({
      myfiles <- input.data_TCRexMerge()
      df <- rbind(myfiles[[1]])
      for (i in 2:length(myfiles)) {
        df <- rbind(df, myfiles[[i]])
      }

      if (nrow(df[-c(grep("[*]", df$CDR3_beta)), ] > 0)) {
        df <- df[-c(grep("[*]", df$CDR3_beta)), ]
      }

      df$CDR3_beta <- gsub("^$", "None", df$CDR3_beta)

      if (nrow(df[-c(grep("None", df$CDR3_beta)), ] > 0)) {
        df <- df[-c(grep("None", df$CDR3_beta)), ]
      }

      df$TRBJ_gene <- gsub("^$", "None", df$TRBJ_gene)

      if (nrow(df[-c(grep("None", df$TRBJ_gene)), ] > 0)) {
        df <- subset(df, df$TRBJ_gene != "None")
      }

      df[!duplicated(df$CDR3_beta), ]
    })

    output$DEx_TCRexFiltered <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 5, scrollX = TRUE), {
      df <- merged_TCRexFiltered()

      df
    })


    output$downloaddf_TCRexFiltered <- downloadHandler(
      filename = function() {
        x <- today()
        paste("TCRex_merged_", x, ".tsv", sep = "")
      },
      content = function(file) {
        df <- merged_TCRexFiltered()
        write.table(df, file, row.names = F, sep = "\t")
      }
    )

    # observeEvent(input$downloaddf_TCRexFiltered,{
    #   write.table(merged_TCRexFiltered(),paste("TCRex_merged_",today(),".tsv", sep = ""), row.names = F,sep = "\t")
    # }
    #
    #             )



    # ClusTCR2 ------
    input.data_ClusTCR2 <- reactive({
      inFile2_ClusTCR2 <- input$file2_ClusTCR2
      if (is.null(inFile2_ClusTCR2)) {
        return(NULL)
      } else {
        dataframe <- read.csv(
          inFile2_ClusTCR2$datapath,
          header = TRUE
        )
      }
    })

    input.data_ClusTCR2_multiple <- reactive({
      inFile.seq <- input$file1_ClusTCR2_multiple
      num <- dim(inFile.seq)[1]
      samples_list <- vector("list", length = num)
      samples_list
      for (i in 1:num) {
        sc <- read.csv(input$file1_ClusTCR2_multiple[[i, "datapath"]])
        samples_list[[i]] <- sc
      }
      samples_list
    })

    merged_Clust_filtered <- reactive({
      inFile2_ClusTCR2 <- input$file1_ClusTCR2_multiple

      validate(
        need(
          nrow(inFile2_ClusTCR2) > 0,
          "Upload mutliple ClusTCR2 files to merge"
        )
      )

      myfiles <- input.data_ClusTCR2_multiple()
      df <- rbind(myfiles[[1]])
      for (i in 2:length(myfiles)) {
        df <- rbind(df, myfiles[[i]])
      }

      if (nrow(df[-c(grep("[*]", df$junction_aa)), ] > 0)) {
        df <- df[-c(grep("[*]", df$junction_aa)), ]
      }

      df$CDR3_beta <- gsub("^$", "None", df$junction_aa)

      if (nrow(df[-c(grep("None", df$CDR3_beta)), ] > 0)) {
        df <- df[-c(grep("None", df$CDR3_beta)), ]
      }

      df[!duplicated(df$junction_aa), ]
    })

    output$DEx_multiple_ClusTCR2 <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 5, scrollX = TRUE), {
      df <- merged_Clust_filtered()

      df
    })

    output$downloaddf_multiple_ClusTCR2 <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Clust_lab_tab_output, "_Multi_ClusTCR_", x, ".csv", sep = "")
      },
      content = function(file) {
        df <- merged_Clust_filtered()
        write.csv(df, file, row.names = F)
      }
    )

    ## clustering images ------

    observe({
      updateSelectInput(
        session,
        "clusTCR2_names",
        choices = names(input.data_ClusTCR2()),
        selected = "junction_aa"
      )
    }) # junction sequence
    observe({
      updateSelectInput(
        session,
        "clusTCR2_Vgene",
        choices = names(input.data_ClusTCR2()),
        selected = "v_call"
      )
    }) # V_gene

    vals_ClusTCR2_old <- reactiveValues(output_dt_old = NULL)
    clust_dt2 <- reactive({
      df1 <- input.data_ClusTCR2()
      validate(
        need(
          nrow(df1) > 0,
          "Upload ClusTCR file"
        )
      )
      df1
    })
    output$clust_dt2_table <- DT::renderDT({
      df1 <- input.data_ClusTCR2()
      validate(
        need(
          nrow(df1) > 0,
          "Upload ClusTCR file"
        )
      )
      df1
    })

    # ## run clustering ------

    vals_ClusTCR2 <- reactiveValues(output_dt2 = NULL)

    observeEvent(input$run_ClusTCR2, {
      df1 <- input.data_ClusTCR2()
      validate(
        need(
          nrow(df1) > 0,
          "Upload ClusTCR file"
        )
      )

      req(input$clusTCR2_names, input$clusTCR2_Vgene)

      clust_dt_DATA_5 <- df1[, names(df1) %in% c(input$clusTCR2_names, input$clusTCR2_Vgene)]
      ptm <- proc.time()

      df_cluster <- ClusTCR(clust_dt_DATA_5, allele = input$allele_ClusTCR2, v_gene = input$clusTCR2_Vgene)

      if (dim(df_cluster)[1] > 1) {
        cluster_lab <- mcl_cluster(df_cluster, expansion = 1, inflation = 1)
        end <- proc.time() - ptm
        cluster_lab[[3]] <- end
        vals_ClusTCR2$output_dt2 <- cluster_lab
        vals_ClusTCR2$output_dt2
      }
    })

    output$verbatum_ClusTCR2 <- renderPrint({
      df1 <- input.data_ClusTCR2()
      validate(
        need(
          nrow(df1) > 0,
          "Upload ClusTCR file"
        )
      )
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")

      # message("Clusters exist")
      if (is.null(vals_ClusTCR2$output_dt2)) {
        message("No clusters == 1 edit distance and therefore MCL not performed")
      } else {

      }

      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })



    ClusTCR2_lab_df <- reactive({
      df1 <- input.data_ClusTCR2()
      validate(
        need(
          nrow(df1) > 0,
          "Upload ClusTCR file"
        )
      )
      req(vals_ClusTCR2$output_dt2)
      output_dt <- vals_ClusTCR2$output_dt2
      df <- output_dt[[1]]

      df <- df[,names(df) %in% c("count","Clust_size_order","CDR3_Vgene")]
      df
    })


    output$ClusTCR2_lab <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      df1 <- input.data_ClusTCR2()
      validate(
        need(
          nrow(df1) > 0,
          "Upload ClusTCR file"
        )
      )
      ClusTCR2_lab_df()
    })

    output$download_ClusTCR_labels <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Clust_lab_tab_output, "_ClusTCR2_output_", x, ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(ClusTCR2_lab_df())
        write.csv(df, file, row.names = F)
      }
    )

    output$ClusTCR2_Time <- renderPrint({
      df1 <- input.data_ClusTCR2()
      validate(
        need(
          nrow(df1) > 0,
          "Upload ClusTCR file"
        )
      )
      req(vals_ClusTCR2$output_dt2)

      output_dt <- vals_ClusTCR2$output_dt2

      if (is.null(output_dt)) {
        return(NULL)
      }
      output_dt[[3]]
    })


    # create plots
    Network_plot_clusTCR2 <- reactive({
      df1 <- input.data_ClusTCR2()
      validate(
        need(
          nrow(df1) > 0,
          "Upload ClusTCR file"
        )
      )
      req(vals_ClusTCR2$output_dt2, input$filter_connections, input$selected_Cluster, input$lab_clust_by,  input$text_size1, input$text_size2)
      Network_df <- vals_ClusTCR2$output_dt2
      set.seed(123)
      # ?netplot
      netplot_ClusTCR2(Network_df,
                       filter_plot = input$filter_connections,
                       Clust_selected = input$selected_Cluster,
                       label = input$lab_clust_by,
                       Clust_column_name = "Clust_size_order",
                       colour = input$colour_ClusTCR2,
                       selected_text_size = input$text_size1,
                       non_selected_text_size = input$text_size2,
                       alpha_selected = 1, alpha_non_selected = 1,
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
        need(
          nrow(df1) > 0,
          "Upload ClusTCR file"
        )
      )
      Network_plot_clusTCR2()
    })

    motif_plot_clusTCR2 <- reactive({
      df1 <- input.data_ClusTCR2()
      validate(
        need(
          nrow(df1) > 0,
          "Upload ClusTCR file"
        )
      )
      req(vals_ClusTCR2$output_dt2)
      Network_df <- vals_ClusTCR2$output_dt2
      set.seed(123)
      motif_plot(Network_df, Clust_selected = input$selected_Cluster, Clust_column_name = input$Clust_size_order)
    })

    output$MP_ClusTCR <- renderPlot({
      df1 <- input.data_ClusTCR2()
      validate(
        need(
          nrow(df1) > 0,
          "Upload ClusTCR file"
        )
      )

      motif_plot_clusTCR2()
    })

    ## downloading plot -----
    output$downloadPlot_Network_plot2 <- downloadHandler(
      filename = function() {
        x <- today()
        paste("clusTCR2_Network_plot_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_Network_plot2, height = input$height_Network_plot2, onefile = FALSE) # open the pdf device
        plot(Network_plot_clusTCR2())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_Network_plot2 <- downloadHandler(
      filename = function() {
        x <- today()
        paste("clusTCR2_Network_plot_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_Network_plot2,
            height = input$height_png_Network_plot2,
            res = input$resolution_PNG_Network_plot2
        )
        plot(Network_plot_clusTCR2())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    # downlaod Motif
    output$downloadPlot_Motif_plot2 <- downloadHandler(
      filename = function() {
        x <- today()
        paste("clusTCR2_Motif_plot", " Cluster ", input$selected_Cluster, " ", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_Motif_plot2, height = input$height_Motif_plot2, onefile = FALSE) # open the pdf device
        plot(motif_plot_clusTCR2())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_Motif_plot2 <- downloadHandler(
      filename = function() {
        x <- today()
        paste("clusTCR2_Motif_plot", " Cluster ", input$selected_Cluster, " ", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Motif_plot2, height = input$height_png_Motif_plot2, res = input$resolution_PNG_Motif_plot2)
        plot(motif_plot_clusTCR2())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )


    # Seurat -----
    ## uploading the raw files ----

    # Define a reactive value to store the data
    data_sc <- reactiveVal(NULL)

    # Update the reactive value with the file data
    observeEvent(input$file_SC, {
      inFile_sc <- input$file_SC
      if (is.null(inFile_sc)) {
        data_sc(NULL)  # Clear the data if no file is uploaded
      } else {
        if (input$df_seruatobj_type == "10x_Genomics (raw)") {
          dataframe <- read.csv(inFile_sc$datapath)
        } else if (input$df_seruatobj_type == "10x_Genomics (.h5)") {
          dataframe <- suppressMessages(Read10X_h5(inFile_sc$datapath, use.names = TRUE, unique.features = TRUE))
        } else if (input$df_seruatobj_type == "BD Rhapsody (Mouse)") {
          dataframe <- read.csv(inFile_sc$datapath, row.names = 1)
        } else if (input$df_seruatobj_type == "BD Rhapsody (Full panel)") {
          dataframe <- read.csv(inFile_sc$datapath, row.names = 1)
        }else if (input$df_seruatobj_type == "BD Rhapsody (Human Immune panel)") {
          dataframe <- read.csv(inFile_sc$datapath, row.names = 1)
        } else {
          dataframe <- read.csv(inFile_sc$datapath, row.names = 1)
        }
        data_sc(dataframe)  # Update the reactive value with the new data
      }
    })

    # Use the reactive value to display data in the table
    output$DEx_header_name_check.dt <- DT::renderDT({
      sc <- data_sc()

      validate(
        need(
          !is.null(sc),
          "Upload data"
        )
      )

      sc
    })

    observeEvent(input$clearDataBtn, {
      # Clear the reactive value to remove the data
      data_sc(NULL)
      # Reset the text input field
      shinyjs::reset("project_name")
      # Use shinyjs to clear file input
      shinyjs::reset("file_SC")
    })

    observeEvent(input$clearDataBtn, {
      data_sc_meta(NULL)
      # Use shinyjs to clear file inputs
      shinyjs::reset("file_SC_meta")
    })


    output$DEx_header_name_check.dt <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 5, scrollX = TRUE), {
      df.test <- data_sc()

      validate(
        need(
          nrow(df.test) > 0,
          error_message_val_sc
        )
      )
      req(df.test)

      if (input$df_seruatobj_type == "10x_Genomics (raw)") {
        names(df.test) <- gsub("[.]1", "-1", names(df.test))
        rownames(df.test) <- make.unique(df.test$Gene_Name)
        df.test2 <- df.test[, !names(df.test) %in% c("Gene_Name")]
      } else if (input$df_seruatobj_type == "10x_Genomics (.h5)") {
        df.test2 <- as.data.frame("file possibly too big and will not be rendered")
      } else if (input$df_seruatobj_type == "BD Rhapsody (Mouse)") {
        names(df.test) <- gsub("X", "", names(df.test))
        df.test2 <- df.test[!rownames(df.test) %in% c("Cell_Index"), ]
      } else if (input$df_seruatobj_type == "BD Rhapsody (Full panel)") {
        names(df.test) <- gsub("X", "", names(df.test))
        df.test2 <- df.test[!rownames(df.test) %in% c("Cell_Index"), ]
      } else if (input$df_seruatobj_type == "BD Rhapsody (Human Immune panel)") {
        names(df.test) <- gsub("X", "", names(df.test))
        df.test2 <- df.test[!rownames(df.test) %in% c("Cell_Index"), ]
      } else {
        names(df.test) <- gsub("[.]", "-", names(df.test))
        rownames(df.test) <- gsub("[.]", "-", rownames(df.test))
        df.test2 <- df.test[!rownames(df.test) %in% c("Cell_Index"), ]
      }
      head(df.test2)[1:6]
    })

    ## reading in 10x and BD data ----
    df_seruatobj <- reactive({
      suppressWarnings({
        df.test <- data_sc()
        validate(
          need(
            length(df.test) > 0,
            error_message_val_sc
          )
        )

        validate(
          need(input$project_name != "", "Please enter a file name.")
        )

        if (input$df_seruatobj_type == "10x_Genomics (raw)") {
          names(df.test) <- gsub("[.]1", "-1", names(df.test))
          rownames(df.test) <- make.unique(df.test$Gene_Name)
          df.test2 <- df.test[, !names(df.test) %in% c("Gene_Name")]

          sc <- CreateSeuratObject(counts = df.test2, assay = "RNA", project = input$project_name)
          sc <- PercentageFeatureSet(sc, pattern = "^MT-", col.name = "mtDNA")
          sc <- PercentageFeatureSet(sc, pattern = "^RP[SL]", col.name = "rRNA")
          sc
        } else if (input$df_seruatobj_type == "10x_Genomics (.h5)") {
          # rownames(df.test$`Gene Expression`) <- gsub("GRCh38___","",rownames(df.test$`Gene Expression`))
          # sc <- CreateSeuratObject(counts = df.test[[1]], project = input$project_name)
          sc <- CreateSeuratObject(counts = df.test, project = input$project_name)
          sc <- PercentageFeatureSet(sc, pattern = "^MT-", col.name = "mtDNA")
          sc <- PercentageFeatureSet(sc, pattern = "^RP[SL]", col.name = "rRNA")
          sc
        } else if (input$df_seruatobj_type == "BD Rhapsody (Mouse)") {
          names(df.test) <- as.character(gsub("X", "", names(df.test)))


          sc <- CreateSeuratObject(counts = df.test, assay = "RNA", project = input$project_name)
          sc <- PercentageFeatureSet(sc, pattern = "^Mt", col.name = "mtDNA")
          sc <- PercentageFeatureSet(sc, pattern = "Rp[sl]", col.name = "rRNA")
          sc
        } else if (input$df_seruatobj_type == "BD Rhapsody (Full panel)") {
          names(df.test) <- as.character(gsub("X", "", names(df.test)))
          sc <- CreateSeuratObject(counts = df.test, assay = "RNA", project = input$project_name)
          sc <- PercentageFeatureSet(sc, pattern = "^MT-", col.name = "mtDNA")
          sc <- PercentageFeatureSet(sc, pattern = "^RP[SL]", col.name = "rRNA")
          sc
        } else if (input$df_seruatobj_type == "BD Rhapsody (Human Immune panel)") {
          names(df.test) <- as.character(gsub("X", "", names(df.test)))
          sc <- CreateSeuratObject(counts = df.test, assay = "RNA", project = input$project_name)
          sc <- PercentageFeatureSet(sc, pattern = "^MT-", col.name = "mtDNA")
          sc <- PercentageFeatureSet(sc, pattern = "^RP[SL]", col.name = "rRNA")
          sc
        } else {
          names(df.test) <- gsub("[.]", "-", names(df.test))
          rownames(df.test) <- gsub("[.]", "-", rownames(df.test))
          sc <- CreateSeuratObject(counts = df.test, assay = "RNA", project = input$project_name)
          sc <- PercentageFeatureSet(sc, pattern = "^MT-", col.name = "mtDNA")
          sc <- PercentageFeatureSet(sc, pattern = "^RP[SL]", col.name = "rRNA")
        }
      })
    })
    before_plot <- reactive({
      sc <- df_seruatobj()
      validate(
        need(
          nrow(sc) > 0,
          "Upload files"
        )
      )
      suppressWarnings({
        VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "mtDNA", "rRNA"), ncol = 2)
      })


    })

    output$before_plot_sc <- renderPlot({
      withProgress(
        message = "Figure is being generated...",
        detail = "",
        value = 0,
        {
          test_fun()
        }
      )
      before_plot()
    })


    output$downloadPlot_before_plot_sc <- downloadHandler(
      filename = function() {
        if (input$project_name == "") {
          showModal(modalDialog(
            title = "Error",
            "Add File Name",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          x <- today()
          paste(input$project_name, "_before_plot_sc_", x, ".pdf", sep = "")
        }
      },
      content = function(file) {
        pdf(file, width = input$width_before_plot_sc, height = input$height_before_plot_sc, onefile = FALSE) # open the pdf device
        plot(before_plot())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_before_plot_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name, "_before_plot_sc_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_before_plot_sc, height = input$height_png_before_plot_sc, res = input$resolution_PNG_before_plot_sc)
        plot(before_plot())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    ##### after filtering plot -------


    vals2 <- reactiveValues(after_violin_plot = NULL)
    observeEvent(input$run_violin, {
      suppressWarnings({
        sc <- df_seruatobj()
        req(sc)
        validate(
          need(input$project_name != "", "Please enter a file name.")
        )

        tryCatch({
          suppressWarnings({
            vals2$after_violin_plot <- subset(sc, subset = nFeature_RNA >= input$features.min & nFeature_RNA <= input$features.max & mtDNA <= input$percent.mt & rRNA >= input$percent.rb)
          })

        }, error = function(e) {
          # Handle the error by displaying a notification
          showNotification("No cells remain. Please check your input parameters.", type = "error")
        })
      })
    })

    output$after_plot_sc <- renderPlot({
      sc <- vals2$after_violin_plot
      validate(
        need(
          nrow(sc) > 0,
          "Run Filtering"
        )
      )
      VlnPlot(vals2$after_violin_plot, features = c("nFeature_RNA", "nCount_RNA", "mtDNA", "rRNA"), ncol = 2)
    })

    output$downloadPlot_after_plot_sc <- downloadHandler(
      filename = function() {
        if (input$project_name == "") {
          showModal(modalDialog(
            title = "Error",
            "Add File Name",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          x <- today()
          paste(input$project_name, "_after_plot_sc_", x, ".pdf", sep = "")
        }
      },
      content = function(file) {
        pdf(file, width = input$width_after_plot_sc, height = input$height_after_plot_sc, onefile = FALSE) # open the pdf device
        plot_after <- VlnPlot(vals2$after_violin_plot, features = c("nFeature_RNA", "nCount_RNA", "mtDNA", "rRNA"), ncol = 2)
        plot(plot_after)
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_after_plot_sc <- downloadHandler(
      filename = function() {
        x <- today()

        paste(input$project_name, "_after_plot_sc_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_after_plot_sc, height = input$height_png_after_plot_sc, res = input$resolution_PNG_after_plot_sc)
        plot_after <- VlnPlot(vals2$after_violin_plot, features = c("nFeature_RNA", "nCount_RNA", "mtDNA", "rRNA"), ncol = 2)
        plot(plot_after)
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    ### normalisationa and feature plot ------
    feature_serartobj <- reactive({
      sc <- vals2$after_violin_plot
      validate(
        need(
          nrow(sc) > 0,
          "Run Clustering"
        )
      )
      sc <- NormalizeData(sc)
      sc <- FindVariableFeatures(sc, selection.method = "vst")
    })

    plot_10_features <- reactive({
      sc <- feature_serartobj()
      validate(
        need(
          nrow(sc) > 0,
          "Run Clustering"
        )
      )
      # Identify the 10 most highly variable genes
      top10 <- head(VariableFeatures(sc), 10)
      # plot variable features with and without labels
      plot1 <- VariableFeaturePlot(sc)
      plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

      plot2
    })
    output$plot_10_features_sc <- renderPlot({
      plot_10_features()
    })

    ###### download variable features (QC) -----
    output$downloadPlot_plot_10_features_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name, "_10_features_sc_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_plot_10_features_sc, height = input$height_plot_10_features_sc, onefile = FALSE) # open the pdf device
        plot(plot_10_features())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_plot_10_features_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name, "_10_features_sc_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_plot_10_features_sc, height = input$height_png_plot_10_features_sc, res = input$resolution_PNG_plot_10_features_sc)
        plot(plot_10_features())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )


    ## PCA and choosing # of dimensions to reduce ----
    create_PCA <- reactive({
      sc <- feature_serartobj()
      validate(
        need(
          nrow(sc) > 0,
          "Run Clustering"
        )
      )
      all.genes <- rownames(sc)
      sc <- ScaleData(sc, features = all.genes)
      sc <- RunPCA(sc, features = VariableFeatures(object = sc))
      sc
    })

    ##### create elbow plot ----
    create_elbowPlot <- reactive({
      ElbowPlot(create_PCA())
    })
    output$create_elbowPlot_sc <- renderPlot({
      create_elbowPlot()
    })

    #######  download Elbow plot -----
    output$downloadPlot_create_elbowPlot_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name, "_elbowPlot_sc_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_create_elbowPlot_sc, height = input$height_create_elbowPlot_sc, onefile = FALSE) # open the pdf device
        plot(create_elbowPlot())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_create_elbowPlot_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name, "_elbowPlot_sc_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_create_elbowPlot_sc, height = input$height_png_create_elbowPlot_sc, res = input$resolution_PNG_create_elbowPlot_sc)
        plot(create_elbowPlot())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )


    ######  UMAP plot -----
    vals_clust <- reactiveValues(sc_clustering = NULL)

    observeEvent(input$run_reduction, {
      sc <- create_PCA()

      validate(
        need(
          nrow(sc) > 0,
          "Run Clustering"
        )
      )

      sc <- FindNeighbors(sc, dims = 1:input$dimension_sc)
      sc <- FindClusters(sc, resolution = input$resolution)
      sc <- RunUMAP(sc, dims = 1:input$dimension_sc)
      vals_clust$sc_clustering <- sc
      vals_clust$sc_clustering
    })

    create_UMAP <- reactive({
      sc <- vals_clust$sc_clustering

      validate(
        need(
          nrow(sc) > 0,
          "Run Clustering"
        )
      )

      DimPlot(sc, reduction = "umap")
    })
    output$create_UMAP_sc <- renderPlot({
      create_UMAP()
    })

    ## Differential expression with two conditions -----

    # Define a reactive value to store the data
    data_sc_meta <- reactiveVal(NULL)

    observeEvent(input$file_SC_meta, {
      inFile_sc_meta <- input$file_SC_meta
      req(inFile_sc_meta)
      dataframe <- read.csv(inFile_sc_meta$datapath)
      data_sc_meta(dataframe)
    })

    # Use the reactive value to display data in the table
    output$DEx_view.meta.dt <- DT::renderDT({
      sc <- data_sc_meta()

      validate(
        need(
          !is.null(sc),
          "Upload metadata"
        )
      )

      sc
    })



    vals_meta.sc <- reactiveValues(metadata_SCobj = NULL)

    observeEvent(input$run_metadata, {
      sc <- vals_clust$sc_clustering
      validate(
        need(
          nrow(sc) > 0,
          "Run clustering or add in metadata"
        )
      )

      meta.data.import <- data_sc_meta()
      validate(
        need(
          nrow(meta.data.import) > 0,
          "Impute Metadata"
        )
      )
      sc@meta.data$Cell_Index <- rownames(sc@meta.data)

      sc@meta.data$order <- 1:length(sc@meta.data$Cell_Index)
      scMeta.data <- sc@meta.data
      meta.data2 <- merge(scMeta.data, meta.data.import, by = "Cell_Index", all.x = T)
      sc@meta.data <- meta.data2
      rownames(sc@meta.data) <- sc@meta.data$Cell_Index
      sc@meta.data <- sc@meta.data[order((sc@meta.data$order)), ]
      vals_meta.sc$metadata_SCobj <- sc
      vals_meta.sc$metadata_SCobj
    })

    output$DEx_table_meta.data <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      calls <- vals_meta.sc$metadata_SCobj
      validate(
        need(
          nrow(calls) > 0,
          "Impute metadata"
        )
      )
      calls@meta.data$Cell_Index <- gsub("[.]", "-", calls@meta.data$Cell_Index)
      calls@meta.data
    })

    # save Seurat object -----
    output$downloaddf_SeruatObj <- downloadHandler(
      filename = function() {
        x <- today()
        # paste(input$project_name,"_SC.obj_",x,".h5Seurat", sep = "")
        paste(input$project_name, "_SC.obj_", x, ".rds", sep = "")
      },
      content = function(file) {
        SaveSeuratRds(vals_meta.sc$metadata_SCobj, file)
        # SaveH5Seurat(vals_meta.sc$metadata_SCobj,file)
      }
    )

    # merging multiple Seurat Obj -----

    data_user_genes <- reactive({
      inFile_user_genes <- input$file_user_genes
      if (is.null(inFile_user_genes)) {
        return(NULL)
      } else {
        dataframe <- read.csv(inFile_user_genes$datapath,header = T)
      }
    })

    getData <- reactive({
      inFile.seq <- input$file1_rds.file

      num <- dim(inFile.seq)[1]
      seurat_object_list <- vector("list", length = num)
      list.sc <- list()

      for (i in 1:num) {
        message("reading in file ", i)
        list.sc[[i]] <- LoadSeuratRds(input$file1_rds.file[[i, "datapath"]])

        message("Reducing file size for file ", i)
        if (input$sample.type.source_merging == "hs") {


          if (length(data_user_genes())>0 & input$include_additional_genes) {
            features.var.needed <- read.csv(system.file("Kmean", "human.variable.features.csv", package = "STEGO.R"))
            user_required_genes <- data_user_genes()
            names(user_required_genes) <- "V1"

            features.var.needed <- merge(features.var.needed,user_required_genes,by = "V1",all = T)
            print(dim(features.var.needed))
          } else {

            features.var.needed <- read.csv(system.file("Kmean", "human.variable.features.csv", package = "STEGO.R"))
          }

          list.sc[[i]] <- subset(list.sc[[i]], features = features.var.needed$V1)
        } else {
          features.var.needed <- read.csv(system.file("Kmean", "human.variable.features.csv", package = "STEGO.R"))
          features.var.needed$V1 <- str_to_title(features.var.needed$V1)
          list.sc[[i]] <- subset(list.sc[[i]], features = features.var.needed$V1)
        }
        list.sc[[i]]@project.name <- "SeuratProject"
        message("Updating cell Index ID for ", i)
        list.sc[[i]]@meta.data$Cell_Index_old <- list.sc[[i]]@meta.data$Cell_Index
        sl <- object.size(list.sc[[i]])
        message(i, " object is ", round(sl[1] / 1000^3, 1), " Gb in R env.")
      }


      list.sc
    })

    output$testing_mult <- renderPrint({
      sc <- input$file1_rds.file
      validate(
        need(
          nrow(sc) > 0,
          "Upload files"
        )
      )
      df <- getData()
      print(df)
    })

    merging_sc_ob <- reactiveValues(Val2 = NULL)

    observeEvent(input$run_merging, {
      sc <- input$file1_rds.file
      validate(
        need(
          nrow(sc) > 0,
          "Upload files"
        )
      )
      list_seurat <- getData()
      num <- length(list_seurat)
      validate(
        need(input$project_name2 != "", "Please enter a project name.")
      )

      if (num > 1) {
        merged_object <- reduce(list_seurat, function(x, y) {
          merge(x = x, y = y, merge.data = TRUE, project = input$project_name2)
        })
        sl <- object.size(merged_object)
        message("The merged object is ", round(sl[1] / 1000^3, 1), " Gb in R env.")
        getData
        merging_sc_ob$Val2 <- merged_object
      }
    })
    merging_sc <- reactive({
      merging_sc_ob$Val2
    })

    output$testing_mult2 <- renderPrint({
      merging_sc()
    })

    output$downloaddf_SeruatObj_merged2 <- downloadHandler(
      filename = function() {
        if (input$project_name2 == "") {
          showModal(modalDialog(
            title = "Error",
            "Add File Name",
            size = "s",
            easyClose = TRUE
          ))
          return(NULL)
        } else {
          x <- today()
          paste(input$project_name2, "_merged_", x, ".rds", sep = "")
        }


      },
      content = function(file) {
        SaveSeuratRds(merging_sc_ob$Val2, file)
      }
    )


    ### performing harmony
    Merged_data_for_harmony <- reactive({
      inFile_sc_Merged_data_for_harmony <- input$file1_rds.Merged_data_for_harmony
      if (is.null(inFile_sc_Merged_data_for_harmony)) {
        return(NULL)
      }

      LoadSeuratRds(inFile_sc_Merged_data_for_harmony$datapath)
    })
    output$testing_mult3 <- renderPrint({
      Merged_data_for_harmony()
    })

    # how to identify if mouse or human
    observe({
      sc <- Merged_data_for_harmony()
      validate(
        need(
          nrow(sc) > 0,
          "Run Variable"
        )
      )
      gene.names <- rownames(sc@assays$RNA$counts)

      gene.names_hs <- gene.names[str_detect(gene.names, "^[A-Z][A-Z/0-9][A-Z/0-9]") &
                                    !str_detect(gene.names, "^X[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^A[0-9][0-9][0-9]") & !str_detect(gene.names, "^AC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^B[0-9][0-9][0-9]") & !str_detect(gene.names, "^BC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^C[0-9][0-9][0-9]") & !str_detect(gene.names, "^BC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^D[0-9][0-9][0-9]") & !str_detect(gene.names, "^E[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^F[0-9][0-9][0-9]") & !str_detect(gene.names, "^G[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^H[0-9][0-9][0-9]") & !str_detect(gene.names, "^I[0-9][0-9][0-9]")]

      gene.names_mm <- gene.names[str_detect(gene.names, "^[A-Z][a-z/0-9][a-z/0-9]") &
                                    !str_detect(gene.names, "^X[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^A[0-9][0-9][0-9]") & !str_detect(gene.names, "^AC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^B[0-9][0-9][0-9]") & !str_detect(gene.names, "^BC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^C[0-9][0-9][0-9]") & !str_detect(gene.names, "^BC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^D[0-9][0-9][0-9]") & !str_detect(gene.names, "^E[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^F[0-9][0-9][0-9]") & !str_detect(gene.names, "^G[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^H[0-9][0-9][0-9]") & !str_detect(gene.names, "^I[0-9][0-9][0-9]")]

      if (length(gene.names_hs) > 0) {
        len_hs <- length(gene.names_hs)
      } else {
        len_hs <- 1
      }

      if (length(gene.names_mm) > 0) {
        len_mm <- length(gene.names_mm)
      } else {
        len_mm <- 1
      }

      ratio_hs.mm <- len_hs / len_mm

      if (ratio_hs.mm > 1) {
        updateSelectInput(
          session,
          "sample.type.source",
          choices = c("hs", "mm"),
          selected = "hs"
        )
      } else if (ratio_hs.mm < 1) {
        updateSelectInput(
          session,
          "sample.type.source",
          choices = c("hs", "mm"),
          selected = "mm"
        )
      } else {
        updateSelectInput(
          session,
          "sample.type.source",
          choices = c("hs", "mm"),
          selected = ""
        )
      }
    })


    # selectInput("","Species",choices = c("hs","mm")),

    Vals_norm <- reactiveValues(Norm1 = NULL)
    Vals_norm1 <- reactiveValues(Norm1 = NULL)
    Vals_norm2 <- reactiveValues(Norm1 = NULL)
    Vals_norm3 <- reactiveValues(Norm1 = NULL)
    Vals_norm4 <- reactiveValues(Norm1 = NULL)
    Vals_normk <- reactiveValues(anno = NULL)
    # Vals_normk <- reactiveValues(anno=NULL)

    output$Scaling_check_output <- renderPrint({
      sc <- Merged_data_for_harmony()
      validate(
        need(
          nrow(sc) > 0,
          "Run Scale"
        )
      )
      req(Vals_norm4$Norm1)
      kmeans <- read.csv(system.file("Kmean", "Kmeans.requires.annotation.csv", package = "STEGO.R"))
      data <- Merged_data_for_harmony()
      VersionControl <- data@version
      VersionControl
    })

    # find variable features
    observeEvent(input$run_var, {
      sc <- Merged_data_for_harmony()
      validate(
        need(
          nrow(sc) > 0,
          "Run Variable"
        )
      )

      if (input$Seruat_version_merge == "V4") {
        # var.genes <- as.data.frame(sc@assays$RNA@var.features)
      } else if (input$Seruat_version_merge == "V5") {
        sc <- JoinLayers(sc,  assay = "RNA")
      } else {

      }

      all.genes <- rownames(sc)

      if (length(all.genes) < 3001) {
        sc <- FindVariableFeatures(sc, selection.method = "vst")
      } else {
        sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)
      }
      Vals_norm4$Norm1 <- sc
    })

    output$var_harmony_verbrose <- renderPrint({
      sc <- Vals_norm4$Norm1
      validate(
        need(
          nrow(sc) > 0,
          "Run Variable"
        )
      )
      sc
    })
    # scale data

    observeEvent(input$run_scale, {
      sc <- Vals_norm4$Norm1
      validate(
        need(
          nrow(sc) > 0,
          "Run Scale"
        )
      )
      req(Vals_norm4$Norm1)
      kmeans <- read.csv(system.file("Kmean", "Kmeans.requires.annotation.csv", package = "STEGO.R"))

      if (input$Seruat_version_merge == "V4") {
        var.genes <- as.data.frame(sc@assays$RNA@var.features)
      } else if (input$Seruat_version_merge == "V5") {
        var.genes <- sc@assays$RNA@meta.data$var.features
        var.genes <- as.data.frame(var.genes[!is.na(var.genes)])
        var.genes
      } else {

      }

      names(var.genes) <- "V1"
      if (input$sample.type.source == "hs") {
        kmeans2 <- as.data.frame(kmeans$Human)
        names(kmeans2) <- "V1"
      } else {
        kmeans2 <- as.data.frame(kmeans$Mouse)
        names(kmeans2) <- "V1"
      }
      Vals_normk$anno <- unique(rbind(var.genes, kmeans2))
    })



    output$Tb_scaling_features_for_annotation <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 5, scrollX = TRUE), {
      sc <- Vals_norm4$Norm1
      validate(
        need(
          nrow(sc) > 0,
          "Run Scale"
        )
      )
      req(Vals_norm4$Norm1)
      kmeans <- read.csv(system.file("Kmean", "Kmeans.requires.annotation.csv", package = "STEGO.R"))

      if (input$Seruat_version_merge == "V4") {
        var.genes <- as.data.frame(sc@assays$RNA@var.features)
      } else if (input$Seruat_version_merge == "V5") {
        var.genes <- sc@assays$RNA@meta.data$var.features
        var.genes <- as.data.frame(var.genes[!is.na(var.genes)])
        var.genes
      } else {

      }

      names(var.genes) <- "V1"
      if (input$sample.type.source == "hs") {
        kmeans2 <- as.data.frame(kmeans$Human)
        names(kmeans2) <- "V1"
      } else {
        kmeans2 <- as.data.frame(kmeans$Mouse)
        names(kmeans2) <- "V1"
      }
      unique(rbind(var.genes, kmeans2))
    })


    observeEvent(input$run_scale, {
      sc <- Vals_norm4$Norm1
      validate(
        need(
          nrow(sc) > 0,
          "Run Scale"
        )
      )
      all.genes <- Vals_normk$anno
      sc <- ScaleData(sc, features = all.genes$V1)
      Vals_norm1$Norm1 <- sc
    })

    output$scale_harmony_verbrose <- renderPrint({
      sc <- Vals_norm1$Norm1

      validate(
        need(
          nrow(sc) > 0,
          "Run Scale"
        )
      )
      sc
    })

    observeEvent(input$run_PCA, {
      sc <- Vals_norm1$Norm1
      validate(
        need(
          nrow(sc) > 0,
          "Run PCA"
        )
      )
      sc <- RunPCA(sc)
      Vals_norm2$Norm1 <- sc
    })

    output$PCA_harmony_verbrose <- renderPrint({
      sc <- Vals_norm2$Norm1
      validate(
        need(
          nrow(sc) > 0,
          "Run PCA"
        )
      )
      head(sc@meta.data)
    })

    observeEvent(input$run_harmony, {
      sc <- Vals_norm2$Norm1
      validate(
        need(
          nrow(sc) > 0,
          "Run Harmony"
        )
      )
      sc <- RunHarmony(sc, "orig.ident", plot_convergence = TRUE)
      Vals_norm3$Norm1 <- sc
    })


    output$harmony_verbrose <- renderPrint({
      sc <- Vals_norm3$Norm1
      validate(
        need(
          nrow(sc) > 0,
          "Run Harmony"
        )
      )
      sc
    })

    observeEvent(input$run_reduction_harmony, {
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

    output$testing_mult4 <- renderPrint({
      df <- Vals_norm$Norm1
      validate(
        need(
          nrow(df) > 0,
          "Run reduction"
        )
      )
      df
    })
    output$create_UMAP_merged <- renderPlot({
      sc <- Vals_norm$Norm1
      validate(
        need(
          nrow(sc) > 0,
          "Run reduction"
        )
      )
      DimPlot(sc, reduction = "umap", group.by = "orig.ident", pt.size = 1)
    })

    # download Harmony merged ----
    output$downloadPlot_sc_merged <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name2, "_SC_Merged_UMAP", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_sc_merged, height = input$height_sc_merged, onefile = FALSE) # open the pdf device
        plot(DimPlot(Vals_norm$Norm1, reduction = "umap", group.by = "orig.ident", pt.size = 1))
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_sc_merged <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name2, "_SC_Merged_UMAP", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_sc_merged, height = input$height_png_sc_merged, res = input$resolution_PNG_sc_merged)
        plot(DimPlot(Vals_norm$Norm1, reduction = "umap", group.by = "orig.ident", pt.size = 1))
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    output$downloaddf_SeruatObj_merged <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name2, "_harmony_", x, ".rds", sep = "")
      },
      content = function(file) {
        SaveSeuratRds(Vals_norm$Norm1, file)
      }
    )

    # remove unwanted cells from scOBJ ------
    getData_SampRemove <- reactive({
      inFile_sc_SampRemove <- input$file1_rds.fileSampsRemove
      if (is.null(inFile_sc_SampRemove)) {
        return(NULL)
      }

      LoadSeuratRds(inFile_sc_SampRemove$datapath)
    })
    output$Preliminary_samp_to_remove <- renderPrint({
      inFile_sc_SampRemove <- input$file1_rds.fileSampsRemove
      if (is.null(inFile_sc_SampRemove)) {
        return(NULL)
      }

      LoadSeuratRds(inFile_sc_SampRemove$datapath)
    })
    observe({
      sc <- getData_SampRemove()
      validate(
        need(
          nrow(sc) > 0,
          "Upload .h5Seurat object"
        )
      )

      df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "Samp_col_SampToRemove",
        choices = names(df3.meta),
        selected = "Sample_Name"
      )
    })

    select_group_metadata_SampToRemove <- reactive({
      sc <- getData_SampRemove()

      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )
      df <- sc@meta.data
      df2 <- as.data.frame(unique(df[names(df) %in% input$Samp_col_SampToRemove]))
      df2 <- as.data.frame(df2)
      df2
    })

    observe({
      df2 <- select_group_metadata_SampToRemove()
      validate(
        need(
          nrow(df2) > 0,
          error_message_val1
        )
      )
      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- as.data.frame(df2[order(df2$V1), ])
      names(df2) <- "V1"
      df2
      # df2 <- subset(df2,df2$V1 != "NA")

      df3 <- subset(df2, df2$V1 != "NA")

      updateSelectInput(
        session,
        "ID_Column_factor_SampToRemove",
        choices = df2$V1,
        selected = df3$V1
      )
    })

    Samps_to_remove <- reactiveValues(Samp1 = NULL)

    observeEvent(input$run_remove_samps, {
      sc <- input$file1_rds.fileSampsRemove
      validate(
        need(
          nrow(sc) > 0,
          "Upload files"
        )
      )

      if (input$DownVColumn == "Meta_data") {
        sc <- getData_SampRemove()

        sc@meta.data$selected <- sc@meta.data[, names(sc@meta.data) %in% input$Samp_col_SampToRemove]
        sc@meta.data$keep <- ifelse(sc@meta.data$selected %in% c(input$ID_Column_factor_SampToRemove), "keep", "NS")
        sc <- subset(x = sc, subset = keep == "keep")
        sc@meta.data <- sc@meta.data[, !names(sc@meta.data) %in% c("selected", "keep")]
        Samps_to_remove$Samp1 <- sc
      } else {
        sc <- getData_SampRemove()
        sc@meta.data$orig.ident_old <- sc@meta.data$orig.ident
        Idents(sc) <- sc@meta.data[, names(sc@meta.data) %in% input$Samp_col_SampToRemove]

        sc2 <- subset(x = sc, downsample = input$downsamp_limit)
        Idents(sc2) <- sc2@meta.data$orig.ident_old
        Samps_to_remove$Samp1 <- sc2
      }
    })

    Filtered_samp_to_remove_process <- reactive({
      Samps_to_remove$Samp1
    })

    output$Filtered_samp_to_remove <- renderPrint({
      sc <- input$file1_rds.fileSampsRemove
      validate(
        need(
          nrow(sc) > 0,
          "Upload files"
        )
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
      filename = function() {
        x <- today()
        paste(input$project_name4, "_keep_", x, ".rds", sep = "")
      },
      content = function(file) {
        sc <- Filtered_samp_to_remove_process()
        SaveSeuratRds(sc, file)
      }
    )

    # Add cell annotations to merged Seurat object -----
    getData_2 <- reactive({
      inFile_sc_pro2 <- input$file1_rds.file2
      if (is.null(inFile_sc_pro2)) {
        return(NULL)
      } else {
        dataframe <- LoadSeuratRds(inFile_sc_pro2$datapath)
      }
    })

    output$testing_mult_anno <- renderPrint({
      sc <- input$file1_rds.file2
      validate(
        need(
          nrow(sc) > 0,
          "Upload files"
        )
      )
      df <- getData_2()

      print(df)
    })

    observe({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      df3.meta <- sc@meta.data

      updateSelectInput(
        session,
        "V_gene_Class_2",
        choices = names(df3.meta),
        selected = "vdj_gene_cdr3_AG_BD"
      )
    })

    Vals_norm2 <- reactiveValues(Norm2 = NULL)
    Annotation <- reactiveValues(LengthofAnno = NULL)
    ## add classification based on TCR-seq ----
    # how to identify if mouse or human
    observe({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Run Variable"
        )
      )
      gene.names <- rownames(sc@assays$RNA$counts)

      gene.names_hs <- gene.names[str_detect(gene.names, "^[A-Z][A-Z/0-9][A-Z/0-9]") &
                                    !str_detect(gene.names, "^X[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^A[0-9][0-9][0-9]") & !str_detect(gene.names, "^AC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^B[0-9][0-9][0-9]") & !str_detect(gene.names, "^BC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^C[0-9][0-9][0-9]") & !str_detect(gene.names, "^BC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^D[0-9][0-9][0-9]") & !str_detect(gene.names, "^E[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^F[0-9][0-9][0-9]") & !str_detect(gene.names, "^G[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^H[0-9][0-9][0-9]") & !str_detect(gene.names, "^I[0-9][0-9][0-9]")]

      gene.names_mm <- gene.names[str_detect(gene.names, "^[A-Z][a-z/0-9][a-z/0-9]") &
                                    !str_detect(gene.names, "^X[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^A[0-9][0-9][0-9]") & !str_detect(gene.names, "^AC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^B[0-9][0-9][0-9]") & !str_detect(gene.names, "^BC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^C[0-9][0-9][0-9]") & !str_detect(gene.names, "^BC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^D[0-9][0-9][0-9]") & !str_detect(gene.names, "^E[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^F[0-9][0-9][0-9]") & !str_detect(gene.names, "^G[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^H[0-9][0-9][0-9]") & !str_detect(gene.names, "^I[0-9][0-9][0-9]")]

      if (length(gene.names_hs) > 0) {
        len_hs <- length(gene.names_hs)
      } else {
        len_hs <- 1
      }

      if (length(gene.names_mm) > 0) {
        len_mm <- length(gene.names_mm)
      } else {
        len_mm <- 1
      }

      ratio_hs.mm <- len_hs / len_mm

      if (ratio_hs.mm > 1) {
        updateSelectInput(
          session,
          "sample.type.source.markers",
          choices = c("hs", "mm"),
          selected = "hs"
        )
      } else if (ratio_hs.mm < 1) {
        updateSelectInput(
          session,
          "sample.type.source.markers",
          choices = c("hs", "mm"),
          selected = "mm"
        )
      } else {
        updateSelectInput(
          session,
          "sample.type.source.markers",
          choices = c("hs", "mm"),
          selected = ""
        )
      }
    })

    TCR_seq_classification <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload files"
        )
      )

      if (input$sample.type.source.markers == "hs") {
        sc@meta.data$unconventional <- ifelse(sc@meta.data$vj_gene_AG == "TRAV1-2.TRAJ33", "MAIT",
                                              ifelse(sc@meta.data$vj_gene_AG == "TRAV1-2.TRAJ12", "MAIT",
                                                     ifelse(sc@meta.data$vj_gene_AG == "TRAV1-2.TRAJ23", "MAIT",
                                                            ifelse(sc@meta.data$vj_gene_AG == "TRAV10.TRAJ18", "iNKT",
                                                                   ifelse(sc@meta.data$v_gene_BD == "TRBV4-1" & sc@meta.data$v_gene_AG == "TRAV17", "CD1b-restricted",
                                                                          ifelse(sc@meta.data$v_gene_BD == "TRBV4-1" & sc@meta.data$v_gene_AG != "TRAV17", "CD1c-restricted (possible)",
                                                                                 # ifelse(sc@meta.data$chain_AG == 'TRG' & sc@meta.data$chain_BD == 'TRB',"gb T cell",
                                                                                 ifelse(sc@meta.data$chain_AG == "TRG" & sc@meta.data$chain_BD == "TRD", "gd T cell", NA)
                                                                          )
                                                                   )
                                                            )
                                                     )
                                              )
        )
      } else {

      }


      sc
    })

    output$TCR_seq_classification_df <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 5, scrollX = TRUE), {
      sc <- TCR_seq_classification()
      sc@meta.data
    })



    # Human annotations -------
    scGATE_anno_generic <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      req(input$threshold_scGate)

      len <- length(rownames(sc@assays$RNA$scale.data))

      if (input$hs_generic_scGATE) {
        scGate_models_DB <- custom_db_scGATE(system.file("scGATE", "human/generic", package = "STEGO.R"))
        models.list <- scGate_models_DB

        sc <- scGate(sc,
                     model = models.list,
                     pos.thr = input$threshold_scGate,
                     neg.thr = input$threshold_scGate,
                     nfeatures = len,
                     reduction = input$reduction_anno,
                     ncores = 8, min.cells = 1
        )

        sc@meta.data$generic <- sc@meta.data$scGate_multi
        sc@meta.data <- sc@meta.data[!grepl("_UCell", names(sc@meta.data))]
        sc@meta.data <- sc@meta.data[!grepl("is.pure_", names(sc@meta.data))]
        sc@meta.data <- sc@meta.data[!grepl("scGate_multi", names(sc@meta.data))]
        sc
      } else {
        sc
      }
      sc
    })

    output$scGATE_verbatum_generic2 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$hs_generic_scGATE) {
        scGATE_anno_generic()
      } else {
        print("Generic not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })

    scGATE_anno_exhausted <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      req(input$threshold_scGate)

      len <- length(rownames(sc@assays$RNA$scale.data))

      if (input$hs_exhausted_scGATE) {
        scGate_models_DB <- custom_db_scGATE(system.file("scGATE", "human/exhausted", package = "STEGO.R"))

        models.list <- scGate_models_DB

        sc <- scGate(sc,
                     model = models.list,
                     pos.thr = input$threshold_scGate,
                     neg.thr = input$threshold_scGate,
                     nfeatures = len,
                     reduction = input$reduction_anno,
                     ncores = 8, min.cells = 1
        )

        sc@meta.data$exhausted <- sc@meta.data$scGate_multi
        sc@meta.data <- sc@meta.data[!grepl("_UCell", names(sc@meta.data))]
        sc@meta.data <- sc@meta.data[!grepl("is.pure_", names(sc@meta.data))]
        sc@meta.data <- sc@meta.data[!grepl("scGate_multi", names(sc@meta.data))]
        sc
      } else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_exhausted <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$hs_exhausted_scGATE) {
        scGATE_anno_exhausted()
      } else {
        print("exhausted not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })

    scGATE_anno_senescence <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      req(input$threshold_scGate)

      len <- length(rownames(sc@assays$RNA$scale.data))

      if (input$hs_senescence_scGATE) {
        scGate_models_DB <- custom_db_scGATE(system.file("scGATE", "human/senescence", package = "STEGO.R"))
        models.list <- scGate_models_DB

        sc <- scGate(sc,
                     model = models.list,
                     pos.thr = input$threshold_scGate,
                     neg.thr = input$threshold_scGate,
                     nfeatures = len,
                     reduction = input$reduction_anno,
                     ncores = 8, min.cells = 1
        )

        sc@meta.data$senescence <- sc@meta.data$scGate_multi
        # sc@meta.data$stress <- ifelse(grepl("Exhau",sc@meta.data$Stress),"Exhausted",sc@meta.data$Stress)
        # sc@meta.data$stress <- ifelse(grepl("Multi",sc@meta.data$Stress),"Exhausted.Senescence",sc@meta.data$Stress)
        sc@meta.data <- sc@meta.data[!grepl("_UCell", names(sc@meta.data))]
        sc@meta.data <- sc@meta.data[!grepl("is.pure_", names(sc@meta.data))]
        sc@meta.data <- sc@meta.data[!grepl("scGate_multi", names(sc@meta.data))]
        sc
      } else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_senescence <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$hs_senescence_scGATE) {
        scGATE_anno_senescence()
      } else {
        print("senescence not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })

    scGATE_anno_function <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      req(input$threshold_scGate)

      len <- length(rownames(sc@assays$RNA$scale.data))

      if (input$hs_function_scGATE) {
        scGate_models_DB <- custom_db_scGATE(system.file("scGATE", "human/function", package = "STEGO.R"))

        models.list <- scGate_models_DB

        sc <- scGate(sc,
                     model = models.list,
                     pos.thr = input$threshold_scGate,
                     neg.thr = input$threshold_scGate,
                     nfeatures = len,
                     reduction = input$reduction_anno,
                     ncores = 8, min.cells = 1
        )

        sc@meta.data$Tcell.function <- sc@meta.data$scGate_multi
        sc@meta.data <- sc@meta.data[!grepl("_UCell", names(sc@meta.data))]
        sc@meta.data <- sc@meta.data[!grepl("is.pure_", names(sc@meta.data))]
        sc@meta.data <- sc@meta.data[!grepl("scGate_multi", names(sc@meta.data))]
        sc
      } else {
        sc
      }
      sc
    })

    output$scGATE_verbatum_function <- renderPrint({
      if (input$hs_function_scGATE) {
        sc <- scGATE_anno_function()
        table(sc@meta.data$Tcell.function)
      } else {
        print("Function not run")
      }
    })

    scGATE_anno_cycling <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      req(input$threshold_scGate)

      len <- length(rownames(sc@assays$RNA$scale.data))

      if (input$hs_cycling_scGATE) {
        scGate_models_DB <- custom_db_scGATE(system.file("scGATE", "human/cycling", package = "STEGO.R"))

        models.list <- scGate_models_DB
        print(models.list)
        sc <- scGate(sc,
                     model = models.list,
                     pos.thr = input$threshold_scGate,
                     neg.thr = input$threshold_scGate,
                     nfeatures = len,
                     reduction = input$reduction_anno,
                     ncores = 8, min.cells = 1
        )

        sc@meta.data$cycling <- sc@meta.data$scGate_multi
        sc@meta.data <- sc@meta.data[!grepl("_UCell", names(sc@meta.data))]
        sc@meta.data <- sc@meta.data[!grepl("is.pure_", names(sc@meta.data))]
        sc@meta.data <- sc@meta.data[!grepl("scGate_multi", names(sc@meta.data))]
        sc
      } else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_cycling <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$hs_cycling_scGATE) {
        scGATE_anno_cycling()
      } else {
        print("cycling not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })

    scGATE_anno_TCRseq <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload files"
        )
      )
      if (input$hs_TCRseq_scGATE) {
        if (input$sample.type.source.markers == "hs") {
          sc@meta.data$unconventional <- ifelse(sc@meta.data$vj_gene_AG == "TRAV1-2.TRAJ33", "MAIT",
                                                ifelse(sc@meta.data$vj_gene_AG == "TRAV1-2.TRAJ12", "MAIT",
                                                       ifelse(sc@meta.data$vj_gene_AG == "TRAV1-2.TRAJ23", "MAIT",
                                                              ifelse(sc@meta.data$vj_gene_AG == "TRAV1-2.TRAJ9", "CD1b-restricted(poss)",
                                                                     ifelse(sc@meta.data$vj_gene_AG == "TRAV10.TRAJ18", "iNKT",
                                                                            ifelse(sc@meta.data$v_gene_BD == "TRBV4-1" & sc@meta.data$v_gene_AG == "TRAV17", "CD1b-restricted(poss)",
                                                                                   ifelse(sc@meta.data$v_gene_BD == "TRBV4-1" & sc@meta.data$v_gene_AG != "TRAV17", "CD1c-restricted(poss)",
                                                                                          # ifelse(sc@meta.data$chain_AG == 'TRG' & sc@meta.data$chain_BD == 'TRB',"gb T cell",
                                                                                          ifelse(sc@meta.data$chain_AG == "TRG" & sc@meta.data$chain_BD == "TRD", "gd T cell",
                                                                                                 ifelse(sc@meta.data$chain_AG == "TRA" & sc@meta.data$chain_BD == "TRB", "ab T cell", "")
                                                                                          )
                                                                                   )
                                                                            )
                                                                     )
                                                              )
                                                       )
                                                )
          )
          sc
        } else {
        }
      }
      sc
    })
    output$scGATE_verbatum_TCRseq <- renderPrint({
      if (input$hs_TCRseq_scGATE) {
        sc <- scGATE_anno_TCRseq()
        table(sc@meta.data$unconventional)
      } else {
        print("TCR-seq not run")
      }
    })

    # BD rhapsody gates Mouse Full panel -----

    scGATE_anno_BD_MM.FP_T.cell <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )

      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len > 2000) {
        len <- 2000
      }

      if (input$BDrhapsody_scGATE.MM.Tcell == T) {
        models.list <- list()
        my_scGate_model <- gating_model(name = "Cd8b1", level = 1, signature = c("Cd8b1", "Cd8a"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "CD4neg", level = 2, signature = c("Cd4"), positive = F)
        models.list$CD8 <- my_scGate_model

        my_scGate_model <- gating_model(name = "Cd8b1", level = 1, signature = c("Cd8b1", "Cd8a"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "CD4neg", level = 2, signature = c("Cd4"), positive = T)
        models.list$CD4 <- my_scGate_model

        my_scGate_model <- gating_model(name = "Cd8b1", level = 1, signature = c("Cd8b1", "Cd8a"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "CD4neg", level = 2, signature = c("Cd4"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "CD3E", level = 3, signature = c("Cd3e"), positive = T)
        models.list$DN <- my_scGate_model
        obj <- scGate(sc, model = models.list, pos.thr = input$threshold_scGate, neg.thr = input$threshold_scGate, nfeatures = len, ncores = 8)

        sc@meta.data$T_cells <- obj@meta.data$scGate_multi
        sc
      } else {
        sc
      }
    })
    output$scGATE_verbatum_BDrhapsody_MM.FP.Tcell <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$BDrhapsody_scGATE.MM.Tcell == T) {
        scGATE_anno_BD_MM.FP_T.cell()
      } else {
        print("MM T cell not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })

    scGATE_anno_BD_MM.FP_Memory <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )

      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len > 2000) {
        len <- 2000
      }
      if (input$BDrhapsody_scGATE.MM.Memory == T) {
        models.list <- list()
        my_scGate_model <- gating_model(name = "Cd44", level = 1, signature = c("Cd44"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "Sell", level = 2, signature = c("Sell"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Ly6c2", level = 3, signature = c("Ly6c2"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "Itga4", level = 4, signature = c("Itga4"), positive = F)
        models.list$Naive <- my_scGate_model

        my_scGate_model <- gating_model(name = "Cd44", level = 1, signature = c("Cd44"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Sell", level = 2, signature = c("Sell"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Ly6c2", level = 3, signature = c("Ly6c2"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Itga4", level = 4, signature = c("Itga4"), positive = F)
        models.list$VM <- my_scGate_model

        my_scGate_model <- gating_model(name = "Cd44", level = 1, signature = c("Cd44"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Sell", level = 2, signature = c("Sell"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Ly6c2", level = 3, signature = c("Ly6c2"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Itga4", level = 4, signature = c("Itga4"), positive = T)
        models.list$CM <- my_scGate_model

        my_scGate_model <- gating_model(name = "Cd44", level = 1, signature = c("Cd44"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Sell", level = 2, signature = c("Sell"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "Itga4", level = 3, signature = c("Itga4"), positive = T)
        models.list$Eff.or.RM <- my_scGate_model

        obj <- scGate(sc, model = models.list, pos.thr = input$threshold_scGate, neg.thr = input$threshold_scGate, nfeatures = len, ncores = 8)

        sc@meta.data$Memory <- obj@meta.data$scGate_multi
      } else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_BDrhapsody_MM.FP.Memory <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$BDrhapsody_scGATE.MM.Memory == T) {
        scGATE_anno_BD_MM.FP_Memory()
      } else {
        print("MM memory not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })

    scGATE_anno_BD_MM.FP_signatures <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )

      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len > 2000) {
        len <- 2000
      }

      if (input$BDrhapsody_scGATE.MM.signatures == T) {
        models.list <- list()
        my_scGate_model <- gating_model(name = "Gzma", level = 1, signature = c("Gzma"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Gzmb", level = 2, signature = c("Gzmb"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Prf1", level = 3, signature = c("Prf1"), positive = T)
        models.list$Eff <- my_scGate_model

        my_scGate_model <- gating_model(name = "Sell", level = 1, signature = c("Sell"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Ccr7", level = 2, signature = c("Ccr7"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Cxcr6", level = 3, signature = c("Cxcr6"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "Itga1", level = 4, signature = c("Itga1"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "Prdm1", level = 5, signature = c("Prdm1"), positive = F)
        models.list$Recirculation <- my_scGate_model

        my_scGate_model <- gating_model(name = "Sell", level = 1, signature = c("Sell"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "Ccr7", level = 2, signature = c("Ccr7"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "Cxcr6", level = 3, signature = c("Cxcr6"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Itga1", level = 4, signature = c("Itga1"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Prdm1", level = 5, signature = c("Prdm1"), positive = T)
        models.list$Tiss_res <- my_scGate_model

        my_scGate_model <- gating_model(name = "Il2ra", level = 2, signature = c("Il2ra"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Foxp3", level = 2, signature = c("Foxp3"), positive = T)
        models.list$Treg <- my_scGate_model

        obj <- scGate(sc, model = models.list, pos.thr = input$threshold_scGate, neg.thr = input$threshold_scGate, nfeatures = len, ncores = 8)
        sc@meta.data$Signature <- obj@meta.data$scGate_multi
        sc@meta.data$Signature <- ifelse(sc@meta.data$Signature == "NA", NA, sc@meta.data$Signature)
      } else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_BDrhapsody_MM.FP.signatures <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$BDrhapsody_scGATE.MM.signatures == T) {
        scGATE_anno_BD_MM.FP_signatures()
      } else {
        print("MM signatures not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })

    scGATE_anno_BD_MM.FP_Innate.NK <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )

      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len > 2000) {
        len <- 2000
      }

      if (input$BDrhapsody_scGATE.MM.Innate.NK == T) {
        models.list <- list()
        my_scGate_model <- gating_model(name = "Gzma", level = 1, signature = c("Gzma"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Gzmb", level = 2, signature = c("Gzmb"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Prf1", level = 3, signature = c("Prf1"), positive = T)
        models.list$Eff <- my_scGate_model

        my_scGate_model <- gating_model(name = "Sell", level = 1, signature = c("Sell"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Ccr7", level = 2, signature = c("Ccr7"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Cxcr6", level = 3, signature = c("Cxcr6"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "Itga1", level = 4, signature = c("Itga1"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "Prdm1", level = 5, signature = c("Prdm1"), positive = F)
        models.list$Recirculation <- my_scGate_model

        my_scGate_model <- gating_model(name = "Sell", level = 1, signature = c("Sell"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "Ccr7", level = 2, signature = c("Ccr7"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "Cxcr6", level = 3, signature = c("Cxcr6"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Itga1", level = 4, signature = c("Itga1"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Prdm1", level = 5, signature = c("Prdm1"), positive = T)
        models.list$Tiss_res <- my_scGate_model

        my_scGate_model <- gating_model(name = "Il2ra", level = 2, signature = c("Il2ra"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Foxp3", level = 2, signature = c("Foxp3"), positive = T)
        models.list$Treg <- my_scGate_model

        obj <- scGate(sc, model = models.list, pos.thr = input$threshold_scGate, neg.thr = input$threshold_scGate, nfeatures = len, ncores = 8)
        sc@meta.data$Innate.NK <- obj@meta.data$scGate_multi
        sc@meta.data$Innate.NK <- ifelse(sc@meta.data$Innate.NK == "NA", NA, sc@meta.data$Innate.NK)
      } else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_BDrhapsody_MM.FP.Innate.NK <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$BDrhapsody_scGATE.MM.Innate.NK == T) {
        scGATE_anno_BD_MM.FP_Innate.NK()
      } else {
        print("MM Innate.NK not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })

    scGATE_anno_BD_MM.FP_TNF.IFNg <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )

      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len > 2000) {
        len <- 2000
      }

      if (input$BDrhapsody_scGATE.MM.TNF.IFNg == T) {
        models.list <- list()
        models.list

        my_scGate_model <- gating_model(name = "Ifng", level = 1, signature = c("Ifng"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Tnf", level = 2, signature = c("Tnf"), positive = T)
        models.list$IFNG.TNF <- my_scGate_model

        my_scGate_model <- gating_model(name = "Ifng", level = 1, signature = c("Ifng"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Tnf", level = 2, signature = c("Tnf"), positive = F)
        models.list$IFNG <- my_scGate_model

        my_scGate_model <- gating_model(name = "Ifng", level = 1, signature = c("Ifng"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "Tnf", level = 2, signature = c("Tnf"), positive = T)
        models.list$TNF <- my_scGate_model


        obj <- scGate(sc, model = models.list, pos.thr = input$threshold_scGate, neg.thr = input$threshold_scGate, nfeatures = len, ncores = 8)

        sc@meta.data$TNF.IFNg <- obj@meta.data$scGate_multi
        sc@meta.data$TNF.IFNg <- ifelse(sc@meta.data$TNF.IFNg == "NA", NA, sc@meta.data$TNF.IFNg)
      } else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_BDrhapsody_MM.FP.TNF.IFNg <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$BDrhapsody_scGATE.MM.TNF.IFNg == T) {
        scGATE_anno_BD_MM.FP_TNF.IFNg()
      } else {
        print("MM TNF.IFNg not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })

    scGATE_anno_BD_MM.FP_subtypes <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )

      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len > 2000) {
        len <- 2000
      }

      if (input$BDrhapsody_scGATE.MM.subtypes == T) {
        models.list <- list()
        my_scGate_model <- gating_model(name = "Ccr4", level = 1, signature = c("Ccr4"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Foxp3", level = 2, signature = c("Foxp3"), positive = F)
        models.list$Th2_like <- my_scGate_model

        my_scGate_model <- gating_model(name = "Cxcr3", level = 1, signature = c("Cxcr3"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Tbx21", level = 2, signature = c("Tbx21"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "CD8", level = 3, signature = c("Cd8a", "Cd8b1"), positive = F)
        models.list$Th1_like <- my_scGate_model
        my_scGate_model <- gating_model(name = "RORC", level = 1, signature = c("Rorc"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Foxp3", level = 2, signature = c("Foxp3"), positive = F)
        models.list$Th17_like <- my_scGate_model
        my_scGate_model <- gating_model(name = "Cxcr5", level = 1, signature = c("Cxcr5"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Foxp3", level = 2, signature = c("Foxp3"), positive = F)
        models.list$Tfh_like <- my_scGate_model
        my_scGate_model <- gating_model(name = "Foxp3", level = 1, signature = c("Foxp3"), positive = T)
        models.list$Treg_like <- my_scGate_model

        my_scGate_model <- gating_model(name = "Klrk1", level = 1, signature = c("Klrk1"), positive = T) # https://pubmed.ncbi.nlm.nih.gov/36705564/
        my_scGate_model <- gating_model(my_scGate_model, name = "Il7r", level = 2, signature = c("Il7r"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Cxcr3", level = 3, signature = c("Cxcr3"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "Tbx21", level = 3, signature = c("Tbx21"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "Cd8a", level = 4, signature = c("Cd8a"), positive = T)
        models.list$KILR_CD8_effector <- my_scGate_model

        my_scGate_model <- gating_model(name = "Gzma", level = 1, signature = c("Gzma"), positive = T) # https://pubmed.ncbi.nlm.nih.gov/36705564/
        my_scGate_model <- gating_model(my_scGate_model, name = "Gzmk", level = 2, signature = c("Gzmk"), positive = T)
        models.list$Effector_CD8 <- my_scGate_model

        my_scGate_model <- gating_model(name = "CD3e", level = 1, signature = c("Cd3e"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "Cd14", level = 2, signature = c("Cd14"), positive = T)
        models.list$Mono.Tcell.Complex <- my_scGate_model

        obj <- scGate(sc, model = models.list, pos.thr = input$threshold_scGate, neg.thr = input$threshold_scGate, nfeatures = len, ncores = 8)
        sc@meta.data$subtypes <- obj@meta.data$scGate_multi
        sc@meta.data$subtypes <- ifelse(sc@meta.data$subtypes == "NA", NA, sc@meta.data$subtypes)
      } else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_BDrhapsody_MM.FP.subtypes <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$BDrhapsody_scGATE.MM.subtypes == T) {
        scGATE_anno_BD_MM.FP_subtypes()
      } else {
        print("MM subtypes not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })

    scGATE_anno_BD_MM.FP_other <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )

      req(input$threshold_scGate)
      len <- length(rownames(sc@assays$RNA$scale.data))
      if (len > 2000) {
        len <- 2000
      }

      if (input$BDrhapsody_scGATE.MM.other == T) {
        models.list <- list()
        my_scGate_model <- gating_model(name = "Pdcd1", level = 1, signature = c("Pdcd1", "Tigit"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "B3GAT1", level = 2, signature = c("B3gat1", "Klrg1"), positive = T)
        models.list$exhausted_and_senescent <- my_scGate_model

        my_scGate_model <- gating_model(name = "PDCD1", level = 1, signature = c("Pdcd1", "Tigit"), positive = T)
        my_scGate_model <- gating_model(my_scGate_model, name = "B3GAT1", level = 2, signature = c("B3gat1", "Klrg1"), positive = F)
        models.list$exhausted <- my_scGate_model

        my_scGate_model <- gating_model(name = "PDCD1", level = 1, signature = c("Pdcd1", "Tigit"), positive = F)
        my_scGate_model <- gating_model(my_scGate_model, name = "B3GAT1", level = 2, signature = c("B3gat1", "Klrg1"), positive = T)
        models.list$senescent <- my_scGate_model
        my_scGate_model <- gating_model(name = "MKI67", level = 1, signature = c("Mki67", "Top2a"), positive = T)
        models.list$CellCycling <- my_scGate_model

        obj <- scGate(sc, model = models.list, pos.thr = input$threshold_scGate, neg.thr = input$threshold_scGate, nfeatures = len, ncores = 8)
        sc@meta.data$other <- obj@meta.data$scGate_multi
        sc@meta.data$other <- ifelse(sc@meta.data$other == "NA", NA, sc@meta.data$other)
      } else {
        sc
      }
      sc
    })
    output$scGATE_verbatum_BDrhapsody_MM.FP.other <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$BDrhapsody_scGATE.MM.other == T) {
        scGATE_anno_BD_MM.FP_other()
      } else {
        print("MM other not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })
    # BD rhapsody gates Mouse immune panel -----

    # verbatium outputs -----

    # scGATE user_Custom -------
    # geneset1
    scGate_anno_GeneSet1 <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      if (input$GeneSet1_scGate) {
        req(input$threshold_scGate)
        req(scGate_models_DB_geneset1)

        len <- length(rownames(sc@assays$RNA$scale.data))

        models.list <- scGate_models_DB_geneset1

        sc <- scGate(sc,
                     model = models.list,
                     pos.thr = input$threshold_scGate,
                     neg.thr = input$threshold_scGate,
                     nfeatures = len,
                     reduction = input$reduction_anno,
                     ncores = 8, min.cells = 1
        )

        sc@meta.data$geneSet1 <- sc@meta.data$scGate_multi

        sc
      } else {
        sc
      }
      sc
    })

    create_UMAP_custom_1 <- reactive({
      sc <- scGate_anno_GeneSet1()

      validate(
        need(
          nrow(sc) > 0,
          "Run Clustering"
        )
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet1", pt.size = 1)
    })
    output$create_custom_1 <- renderPlot({
      sc <- scGate_anno_GeneSet1()
      req(sc)
      create_UMAP_custom_1()
    })
    output$scGATE_verbatum_GeneSet1 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$GeneSet1_scGate) {
        scGate_anno_GeneSet1()
      } else {
        print("Geneset1 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })
    # geneset2
    scGate_anno_GeneSet2 <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      if (input$GeneSet2_scGate) {
        req(input$threshold_scGate)
        req(scGate_models_DB_geneset2)

        len <- length(rownames(sc@assays$RNA$scale.data))

        models.list <- scGate_models_DB_geneset2
        sc <- scGate(sc,
                     model = models.list,
                     pos.thr = input$threshold_scGate,
                     neg.thr = input$threshold_scGate,
                     nfeatures = len,
                     reduction = input$reduction_anno,
                     ncores = 8, min.cells = 1
        )
        sc@meta.data$geneSet2 <- sc@meta.data$scGate_multi

        sc
      } else {
        sc
      }
      sc
    })

    create_UMAP_custom_2 <- reactive({
      sc <- scGate_anno_GeneSet2()

      validate(
        need(
          nrow(sc) > 0,
          "Run Clustering"
        )
      )
      DimPlot(sc, reduction = "umap", group.by = "geneSet2", pt.size = 1)
    })
    output$create_custom_2 <- renderPlot({
      sc <- scGate_anno_GeneSet2()
      req(sc)
      create_UMAP_custom_2()
    })

    output$scGATE_verbatum_GeneSet2 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$GeneSet2_scGate) {
        scGate_anno_GeneSet2()
      } else {
        print("Geneset2 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })
    # geneset3
    scGate_anno_GeneSet3 <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      if (input$GeneSet3_scGate) {
        req(input$threshold_scGate)
        req(scGate_models_DB_geneset3)
        len <- length(rownames(sc@assays$RNA$scale.data))

        models.list <- scGate_models_DB_geneset3
        sc <- scGate(sc,
                     model = models.list,
                     pos.thr = input$threshold_scGate,
                     neg.thr = input$threshold_scGate,
                     nfeatures = len,
                     reduction = input$reduction_anno,
                     ncores = 8, min.cells = 1
        )

        sc@meta.data$geneSet3 <- sc@meta.data$scGate_multi

        sc
      } else {
        sc
      }
      sc
    })

    create_UMAP_custom_3 <- reactive({
      sc <- scGate_anno_GeneSet3()

      validate(
        need(
          nrow(sc) > 0,
          "Run Clustering"
        )
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet3", pt.size = 1)
    })
    output$create_custom_3 <- renderPlot({
      create_UMAP_custom_3()
    })

    output$scGATE_verbatum_GeneSet3 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$GeneSet3_scGate) {
        scGate_anno_GeneSet3()
      } else {
        print("Geneset3 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })
    # geneset4
    scGate_anno_GeneSet4 <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      if (input$GeneSet4_scGate) {
        req(input$threshold_scGate)
        len <- length(rownames(sc@assays$RNA$scale.data))
        models.list <- scGate_models_DB_geneset4
        sc <- scGate(sc,
                     model = models.list,
                     pos.thr = input$threshold_scGate,
                     neg.thr = input$threshold_scGate,
                     nfeatures = len,
                     reduction = input$reduction_anno,
                     ncores = 8, min.cells = 1
        )
        sc@meta.data$geneSet4 <- obj@meta.data$scGate_multi

        sc
      } else {
        sc
      }
      sc
    })

    create_UMAP_custom_4 <- reactive({
      sc <- scGate_anno_GeneSet4()

      validate(
        need(
          nrow(sc) > 0,
          "Run Clustering"
        )
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet4", pt.size = 1)
    })
    output$create_custom_4 <- renderPlot({
      create_UMAP_custom_4()
    })

    output$scGATE_verbatum_GeneSet4 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$GeneSet4_scGate) {
        scGate_anno_GeneSet4()
      } else {
        print("Geneset4 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })
    # geneset5
    scGate_anno_GeneSet5 <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      if (input$GeneSet5_scGate) {
        req(input$threshold_scGate)

        len <- length(rownames(sc@assays$RNA$scale.data))

        models.list <- scGate_models_DB_geneset5
        sc <- scGate(sc,
                     model = models.list,
                     pos.thr = input$threshold_scGate,
                     neg.thr = input$threshold_scGate,
                     nfeatures = len,
                     reduction = input$reduction_anno,
                     ncores = 8, min.cells = 1
        )
        sc@meta.data$geneSet5 <- obj@meta.data$scGate_multi

        sc
      } else {
        sc
      }
      sc
    })

    create_UMAP_custom_5 <- reactive({
      sc <- scGate_anno_GeneSet5()

      validate(
        need(
          nrow(sc) > 0,
          "Run Clustering"
        )
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet5", pt.size = 1)
    })
    output$create_custom_5 <- renderPlot({
      create_UMAP_custom_5()
    })

    output$scGATE_verbatum_GeneSet5 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$GeneSet5_scGate) {
        scGate_anno_GeneSet5()
      } else {
        print("Geneset5 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })
    # geneset6
    scGate_anno_GeneSet6 <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      if (input$GeneSet6_scGate) {
        req(input$threshold_scGate)
        len <- length(rownames(sc@assays$RNA$scale.data))

        models.list <- scGate_models_DB_geneset6
        sc <- scGate(sc,
                     model = models.list,
                     pos.thr = input$threshold_scGate,
                     neg.thr = input$threshold_scGate,
                     nfeatures = len,
                     reduction = input$reduction_anno,
                     ncores = 8, min.cells = 1
        )
        sc@meta.data$geneSet6 <- sc@meta.data$scGate_multi

        sc
      } else {
        sc
      }
      sc
    })

    create_UMAP_custom_6 <- reactive({
      sc <- scGate_anno_GeneSet6()

      validate(
        need(
          nrow(sc) > 0,
          "Run Clustering"
        )
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet6", pt.size = 1)
    })
    output$create_custom_6 <- renderPlot({
      create_UMAP_custom_6()
    })

    output$scGATE_verbatum_GeneSet6 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$GeneSet6_scGate) {
        scGate_anno_GeneSet6()
      } else {
        print("Geneset6 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })
    # geneset7
    scGate_anno_GeneSet7 <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      if (input$GeneSet7_scGate) {
        req(input$threshold_scGate)
        len <- length(rownames(sc@assays$RNA$scale.data))
        if (len > 2000) {
          len <- 2000
        }
        models.list <- scGate_models_DB_geneset7
        obj <- obj <- obj <- scGate(sc, model = models.list, pos.thr = input$threshold_scGate, neg.thr = input$threshold_scGate, nfeatures = len, ncores = 8)
        sc@meta.data$geneSet7 <- obj@meta.data$scGate_multi

        sc
      } else {
        sc
      }
      sc
    })

    create_UMAP_custom_7 <- reactive({
      sc <- scGate_anno_GeneSet7()

      validate(
        need(
          nrow(sc) > 0,
          "Run Clustering"
        )
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet7", pt.size = 1)
    })
    output$create_custom_7 <- renderPlot({
      create_UMAP_custom_7()
    })

    output$scGATE_verbatum_GeneSet7 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$GeneSet7_scGate) {
        scGate_anno_GeneSet7()
      } else {
        print("Geneset7 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })
    # geneset8
    scGate_anno_GeneSet8 <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      if (input$GeneSet8_scGate) {
        req(input$threshold_scGate)
        len <- length(rownames(sc@assays$RNA$scale.data))
        if (len > 2000) {
          len <- 2000
        }
        models.list <- scGate_models_DB_geneset8
        obj <- obj <- obj <- scGate(sc, model = models.list, pos.thr = input$threshold_scGate, neg.thr = input$threshold_scGate, nfeatures = len, ncores = 8)
        sc@meta.data$geneSet8 <- obj@meta.data$scGate_multi

        sc
      } else {
        sc
      }
      sc
    })

    create_UMAP_custom_8 <- reactive({
      sc <- scGate_anno_GeneSet8()

      validate(
        need(
          nrow(sc) > 0,
          "Run Clustering"
        )
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet8", pt.size = 1)
    })
    output$create_custom_8 <- renderPlot({
      create_UMAP_custom_8()
    })
    output$scGATE_verbatum_GeneSet8 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$GeneSet8_scGate) {
        scGate_anno_GeneSet8()
      } else {
        print("Geneset8 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })
    # geneset8
    scGate_anno_GeneSet9 <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      if (input$GeneSet9_scGate) {
        req(input$threshold_scGate)
        len <- length(rownames(sc@assays$RNA$scale.data))
        if (len > 2000) {
          len <- 2000
        }
        models.list <- scGate_models_DB_geneset9
        obj <- obj <- obj <- scGate(sc, model = models.list, pos.thr = input$threshold_scGate, neg.thr = input$threshold_scGate, nfeatures = len, ncores = 8)
        sc@meta.data$geneSet9 <- obj@meta.data$scGate_multi

        sc
      } else {
        sc
      }
      sc
    })

    create_UMAP_custom_9 <- reactive({
      sc <- scGate_anno_GeneSet9()

      validate(
        need(
          nrow(sc) > 0,
          "Run Clustering"
        )
      )

      DimPlot(sc, reduction = "umap", group.by = "geneSet9", pt.size = 1)
    })
    output$create_custom_9 <- renderPlot({
      create_UMAP_custom_9()
    })

    output$scGATE_verbatum_GeneSet9 <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz, type = "output")
      sink(zz, type = "message")
      if (input$GeneSet9_scGate) {
        scGate_anno_GeneSet9()
      } else {
        print("Geneset9 not annotated yet")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep = "\n")
    })


    # creating the final object ------
    scGATE_anno_HS <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )

      if (input$hs_function_scGATE) {
        obj <- scGATE_anno_function()
        sc@meta.data$Tcell.function <- obj@meta.data$Tcell.function
      }
      if (input$hs_generic_scGATE) {
        obj <- scGATE_anno_generic()
        sc@meta.data$generic <- obj@meta.data$generic
      }


      if (input$hs_exhausted_scGATE) {
        obj <- scGATE_anno_exhausted()
        sc@meta.data$exhausted <- obj@meta.data$exhausted
      }
      if (input$hs_senescence_scGATE) {
        obj <- scGATE_anno_senescence()
        sc@meta.data$senescence <- obj@meta.data$senescence
      }
      if (input$hs_cycling_scGATE) {
        obj <- scGATE_anno_cycling()
        sc@meta.data$cycling <- obj@meta.data$cycling
      }
      if (input$hs_TCRseq_scGATE) {
        obj <- scGATE_anno_TCRseq()
        sc@meta.data$unconventional <- obj@meta.data$unconventional
      }

      if (input$GeneSet1_scGate) {
        obj <- scGate_anno_GeneSet1()
        sc@meta.data$geneSet1 <- obj@meta.data$geneSet1
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet1"] <- input$geneset1_name
      }
      if (input$GeneSet2_scGate) {
        obj <- scGate_anno_GeneSet2()
        sc@meta.data$geneSet2 <- obj@meta.data$geneSet2
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet2"] <- input$geneset2_name
      }
      if (input$GeneSet3_scGate) {
        obj <- scGate_anno_GeneSet3()
        sc@meta.data$geneSet3 <- obj@meta.data$geneSet3
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet3"] <- input$geneset3_name
      }
      if (input$GeneSet4_scGate) {
        obj <- scGate_anno_GeneSet4()
        sc@meta.data$geneSet4 <- obj@meta.data$geneSet4
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet4"] <- input$geneset4_name
      }
      if (input$GeneSet5_scGate) {
        obj <- scGate_anno_GeneSet5()
        sc@meta.data$geneSet5 <- obj@meta.data$geneSet5
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet5"] <- input$geneset5_name
      }
      if (input$GeneSet6_scGate) {
        obj <- scGate_anno_GeneSet6()
        sc@meta.data$geneSet6 <- obj@meta.data$geneSet6
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet6"] <- input$geneset6_name
      }
      if (input$GeneSet7_scGate) {
        obj <- scGate_anno_GeneSet7()
        sc@meta.data$geneSet7 <- obj@meta.data$geneSet7
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet7"] <- input$geneset7_name
      }
      if (input$GeneSet8_scGate) {
        obj <- scGate_anno_GeneSet8()
        sc@meta.data$geneSet8 <- obj@meta.data$geneSet8
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet8"] <- input$geneset8_name
      }
      if (input$GeneSet9_scGate) {
        obj <- scGate_anno_GeneSet9()
        sc@meta.data$geneSet9 <- obj@meta.data$geneSet9
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet9"] <- input$geneset9_name
      } else {
        sc
      }

      sc
    })


    scGATE_anno_BD_MM.FP <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )
      if (input$BDrhapsody_scGATE.MM.Tcell) {
        obj <- scGATE_anno_BD_MM.FP_T.cell()
        sc@meta.data$T_cells <- obj@meta.data$T_cells
      }
      if (input$BDrhapsody_scGATE.MM.Memory) {
        obj <- scGATE_anno_BD_MM.FP_Memory()
        sc@meta.data$Memory <- obj@meta.data$Memory
      }
      if (input$BDrhapsody_scGATE.MM.signatures) {
        obj <- scGATE_anno_BD_MM.FP_signatures()
        sc@meta.data$signatures <- obj@meta.data$signatures
      }
      if (input$BDrhapsody_scGATE.MM.Innate.NK) {
        obj <- scGATE_anno_BD_MM.FP_Innate.NK()
        sc@meta.data$Innate.NK <- obj@meta.data$Innate.NK
      }
      if (input$BDrhapsody_scGATE.MM.TNF.IFNg) {
        obj <- scGATE_anno_BD_MM.FP_TNF.IFNg()
        sc@meta.data$TNF.IFNg <- obj@meta.data$TNF.IFNg
      }
      if (input$BDrhapsody_scGATE.MM.subtypes) {
        obj <- scGATE_anno_BD_MM.FP_subtypes()
        sc@meta.data$subtypes <- obj@meta.data$subtypes
      }
      if (input$BDrhapsody_scGATE.MM.other) {
        obj <- scGATE_anno_BD_MM.FP_other()
        sc@meta.data$other <- obj@meta.data$other
      }

      if (input$GeneSet1_scGate) {
        obj <- scGate_anno_GeneSet1()
        sc@meta.data$geneSet1 <- obj@meta.data$geneSet1
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet1"] <- input$geneset1_name
      }
      if (input$GeneSet2_scGate) {
        obj <- scGate_anno_GeneSet2()
        sc@meta.data$geneSet2 <- obj@meta.data$geneSet2
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet2"] <- input$geneset2_name
      }
      if (input$GeneSet3_scGate) {
        obj <- scGate_anno_GeneSet3()
        sc@meta.data$geneSet3 <- obj@meta.data$geneSet3
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet3"] <- input$geneset3_name
      }
      if (input$GeneSet4_scGate) {
        obj <- scGate_anno_GeneSet4()
        sc@meta.data$geneSet4 <- obj@meta.data$geneSet4
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet4"] <- input$geneset4_name
      }
      if (input$GeneSet5_scGate) {
        obj <- scGate_anno_GeneSet5()
        sc@meta.data$geneSet5 <- obj@meta.data$geneSet5
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet5"] <- input$geneset5_name
      }
      if (input$GeneSet6_scGate) {
        obj <- scGate_anno_GeneSet6()
        sc@meta.data$geneSet6 <- obj@meta.data$geneSet6
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet6"] <- input$geneset6_name
      }
      if (input$GeneSet7_scGate) {
        obj <- scGate_anno_GeneSet7()
        sc@meta.data$geneSet7 <- obj@meta.data$geneSet7
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet7"] <- input$geneset7_name
      }
      if (input$GeneSet8_scGate) {
        obj <- scGate_anno_GeneSet8()
        sc@meta.data$geneSet8 <- obj@meta.data$geneSet8
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet8"] <- input$geneset8_name
      }
      if (input$GeneSet9_scGate) {
        obj <- scGate_anno_GeneSet9()
        sc@meta.data$geneSet9 <- obj@meta.data$geneSet9
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet9"] <- input$geneset9_name
      } else {
        sc
      }

      sc
    })
    scGATE_anno_BD_MM.IP <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )

      if (input$GeneSet1_scGate) {
        obj <- scGate_anno_GeneSet1()
        sc@meta.data$geneSet1 <- obj@meta.data$geneSet1
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet1"] <- input$geneset1_name
      }
      if (input$GeneSet2_scGate) {
        obj <- scGate_anno_GeneSet2()
        sc@meta.data$geneSet2 <- obj@meta.data$geneSet2
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet2"] <- input$geneset2_name
      }
      if (input$GeneSet3_scGate) {
        obj <- scGate_anno_GeneSet3()
        sc@meta.data$geneSet3 <- obj@meta.data$geneSet3
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet3"] <- input$geneset3_name
      }
      if (input$GeneSet4_scGate) {
        obj <- scGate_anno_GeneSet4()
        sc@meta.data$geneSet4 <- obj@meta.data$geneSet4
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet4"] <- input$geneset4_name
      }
      if (input$GeneSet5_scGate) {
        obj <- scGate_anno_GeneSet5()
        sc@meta.data$geneSet5 <- obj@meta.data$geneSet5
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet5"] <- input$geneset5_name
      }
      if (input$GeneSet6_scGate) {
        obj <- scGate_anno_GeneSet6()
        sc@meta.data$geneSet6 <- obj@meta.data$geneSet6
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet6"] <- input$geneset6_name
      }
      if (input$GeneSet7_scGate) {
        obj <- scGate_anno_GeneSet7()
        sc@meta.data$geneSet7 <- obj@meta.data$geneSet7
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet7"] <- input$geneset7_name
      }
      if (input$GeneSet8_scGate) {
        obj <- scGate_anno_GeneSet8()
        sc@meta.data$geneSet8 <- obj@meta.data$geneSet8
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet8"] <- input$geneset8_name
      }
      if (input$GeneSet9_scGate) {
        obj <- scGate_anno_GeneSet9()
        sc@meta.data$geneSet9 <- obj@meta.data$geneSet9
        names(sc@meta.data)[names(sc@meta.data) %in% "geneSet9"] <- input$geneset9_name
      } else {
        sc
      }

      sc
    })

    # all.annotations added -----
    output$DEx_table_TcellClass_scGATE <- DT::renderDT(
      escape = FALSE,
      filter = list(position = "top", clear = FALSE),
      options = list(
        autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100),
        pageLength = 5, scrollX = TRUE
      ),
      {
        # c("10x_HS","BD_HS.Immune.Panel","BD_HS.Full.Panel","10x_MM","BD_MM_Full.Panel","BD_MM_Immune.Panel",)
        sc <- getData_2()
        validate(
          need(
            nrow(sc) > 0,
            "Upload file for annotation"
          )
        )
        if (input$Data_types == "10x_HS") {
          sc <- scGATE_anno_HS()
          as.data.frame(sc@meta.data)
        } else if (input$Data_types == "BD_HS.Immune.Panel") {
          sc <- scGATE_anno_HS()
          as.data.frame(sc@meta.data)
        } else if (input$Data_types == "BD_HS.Full.Panel") {
          sc <- scGATE_anno_HS()
          as.data.frame(sc@meta.data)
        } else if (input$Data_types == "10x_MM") {
          sc <- scGATE_anno_BD_MM.FP()
          as.data.frame(sc@meta.data)
        } else if (input$Data_types == "BD_MM_Full.Panel") {
          sc <- scGATE_anno_BD_MM.FP()
          as.data.frame(sc@meta.data)
        } else if (input$Data_types == "BD_MM_Immune.Panel") {
          sc <- scGATE_anno_BD_MM.IP()
          as.data.frame(sc@meta.data)
        } else if (input$Data_types == "TCR-seq") {
          sc <- TCR_seq_classification()
          as.data.frame(sc@meta.data)
        } else {
          df <- as.data.frame("Other panels are under development")
          names(df) <- "V1"
          df
        }
      }
    )

    output$downloaddf_SeruatObj_annotated <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$project_name3, "_annotated_", x, ".rds", sep = "")
      },
      content = function(file) {
        if (input$Data_types == "10x_HS") {
          sc <- scGATE_anno_HS()
          SaveSeuratRds(sc, file)
        } else if (input$Data_types == "BD_HS.Immune.Panel") {
          sc <- scGATE_anno_HS()
          SaveSeuratRds(sc, file)
        } else if (input$Data_types == "BD_HS.Full.Panel") {
          sc <- scGATE_anno_HS()
          SaveSeuratRds(sc, file)
        }
        # mouse panels
        else if (input$Data_types == "10x_MM") {
          sc <- scGATE_anno_BD_MM.FP()
          SaveSeuratRds(sc, file)
        } else if (input$Data_types == "BD_MM_Full.Panel") {
          sc <- scGATE_anno_BD_MM.FP()
          SaveSeuratRds(sc, file)
        } else if (input$Data_types == "BD_MM_Immune.Panel") {
          sc <- scGATE_anno_BD_MM.IP()
          SaveSeuratRds(sc, file)
        } else if (input$Data_types == "TCR-seq") {
          sc <- TCR_seq_classification()
          SaveSeuratRds(sc, file)
        } else {

        }
      }
    )

    ## differential expression -----
    observe({
      df.test <- getData_2()
      validate(
        need(
          nrow(df.test) > 0,
          error_message_val_sc
        )
      )
      updateSelectInput(
        session,
        "meta_data_sc_clust",
        choices = names(df.test@meta.data),
        selected = c("seurat_clusters")
      )
    })
    observe({
      df.test <- getData_2()
      validate(
        need(
          nrow(df.test) > 0,
          error_message_val_sc
        )
      )
      updateSelectInput(
        session,
        "meta_data_sc_",
        choices = names(df.test@meta.data),
        selected = c("orig.ident")
      )
    })
    df_sc_clust <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )

      if (input$multiple_group_sc == "yes") {
        Idents(sc) <- paste(
          sc@meta.data[[input$meta_data_sc_clust]],
          sc@meta.data[[input$meta_data_sc_]]
        )
      } else {
        Idents(sc) <- paste(sc@meta.data[[input$meta_data_sc_clust]])
      }

      sc
    })


    list_of_genes <- reactive({
      sc <- getData_2()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )

      kmeans <- read.csv(system.file("Kmean", "Kmeans.requires.annotation.csv", package = "STEGO.R"))
      var.genes <- as.data.frame(unique(rownames(sc@assays$RNA$scale.data)))
      names(var.genes) <- "V1"

      if (input$sample.type.source.markers == "hs") {
        kmeans2 <- as.data.frame(kmeans$Human)
        names(kmeans2) <- "V1"
      } else {
        kmeans2 <- as.data.frame(kmeans$Mouse)
        names(kmeans2) <- "V1"
      }
      kmeans2
      var.genes

      name.df <- unique(c(kmeans$V1, var.genes$V1))
      name.df <- as.data.frame(name.df)
      names(name.df) <- "Gene_Name"
      name.df <- as.data.frame(name.df[order(name.df$Gene_Name), ])
      names(name.df) <- "Gene_Name"

      name.df
    })




    output$list_of_genes_df <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      df.test <- getData_2()
      validate(
        need(
          nrow(df.test) > 0,
          error_message_val_sc
        )
      )
      list_of_genes()
    })


    output$meta_data_comaprison_check <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      calls <- sc@meta.data
      calls
    })

    vals_clust_markers <- reactiveValues(markers_for_table = NULL)
    observeEvent(input$run_differental.exp, {
      sc <- df_sc_clust()

      validate(
        need(
          nrow(sc) > 0,
          "Run DEx"
        )
      )


      sc.markers <- FindAllMarkers(df_sc_clust(),
                                   only.pos = TRUE,
                                   min.pct = input$min.ptc.sc,
                                   logfc.threshold = input$logfc.ptc.sc,
                                   test.use = input$normalN
      )
      vals_clust_markers$markers_for_table <- sc.markers
    })

    output$DEx_table_clusters <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      calls <- vals_clust_markers$markers_for_table
      calls
    })

    output$downloaddf_DEx_table_clusters <- downloadHandler(
      filename = function() {
        paste("DEx_all_", gsub("-", ".", Sys.Date()), ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(vals_clust_markers$markers_for_table)
        write.csv(df, file, row.names = T)
      }
    )

    #
    # markers_heatmap <- reactive({
    #   sc.markers %>%
    #     group_by(cluster) %>%
    #     top_n(n = 5, wt = avg_log2FC) -> top10
    #   DoHeatmap(sc, features = top10$gene) + NoLegend()
    # })
    #

    ## View featurePlot markers ----
    observeEvent(input$run_string.data3, {
      df.test <- getData_2()
      validate(
        need(
          nrow(df.test) > 0,
          error_message_val_sc
        )
      )
      name.df <- as.data.frame(list_of_genes())
      names(name.df) <- "Gene_Name"
      updateSelectInput(
        session,
        "string.data3",
        choices = name.df$Gene_Name,
        selected = c("GZMB", "CD4", "CD8A","SELL","CD27")
      )
    })

    #
    markers_featurePlot <- reactive({
      sc <- getData_2()
      Feture_plots <- list()
      feature_name <- c(input$string.data3)
      req(input$string.data3)
      x <- length(feature_name)

      if (input$norm_expression_for_all == "yes") {
        for (i in 1:x) {
          Feture_plots[[i]] <- FeaturePlot(sc, features = feature_name[i], raster = FALSE, label = input$label_is_true_features) +
            scale_color_gradientn(colours = c(input$lower_col_FP, input$upper_col_FP), limits = c(0, input$max_norm_FP)) +
            theme(plot.background = element_rect(fill = "white", color = NA)) +
            theme(plot.title = element_text(size = 10))
        }
      } else {
        for (i in 1:x) {
          Feture_plots[[i]] <- FeaturePlot(sc, features = feature_name[i], raster = FALSE, label = input$label_is_true_features) +
            scale_color_gradientn(colours = c(input$lower_col_FP, input$upper_col_FP)) +
            theme(plot.background = element_rect(fill = "white", color = NA)) +
            theme(plot.title = element_text(size = 10))
        }
      }

      n <- length(Feture_plots)
      nCol <- round(sqrt(n), 0)
      do.call("grid.arrange", c(Feture_plots, ncol = nCol))
    })
    output$markers_featurePlot_sc <- renderPlot({
      markers_featurePlot()
    })

    # download markers_featurePlot_sc
    output$downloadPlot_markers_featurePlot_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste("TCR_Explore_markers_featurePlot_sc_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_markers_featurePlot_sc, height = input$height_markers_featurePlot_sc, onefile = FALSE) # open the pdf device

        plot(markers_featurePlot())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_markers_featurePlot_sc <- downloadHandler(
      filename = function() {
        x <- today()
        paste("TCR_Explore_markers_featurePlot_sc_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_markers_featurePlot_sc, height = input$height_png_markers_featurePlot_sc, res = input$resolution_PNG_markers_featurePlot_sc)
        plot(markers_featurePlot())
        dev.off()
      }, contentType = "application/png"
    )



    # Analysis -----
    # observe({
    #
    #
    #   # if (input$colourtype == "Paired") {
    #     updateSelectInput(
    #       session,
    #       "colourtype",
    #       choices = names(input.data.TCR.BD()),
    #       selected = "chain"
    #     )
    #   # } else {
    #
    #   # }
    # })
    #

    observe({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val1
        )
      )

      meta.data <- sc@meta.data


      if ("junction_aa_BD" %in% names(sc@meta.data)) {
        updateSelectInput(
          session,
          "datasource",
          choices = c("10x_Genomics", "BD_Rhapsody_Paired", "BD_Rhapsody_AIRR"),
          selected = "BD_Rhapsody_AIRR"
        )
      } else if ("cdr3_BD" %in% names(sc@meta.data)) {
        sc@meta.data$orientation <- (ifelse(grepl("^C", sc@meta.data$cdr3_BD) & grepl("F$", sc@meta.data$cdr3_BD), "C___F", NA))
        if ("C___F" %in% sc@meta.data$orientation) {
          updateSelectInput(
            session,
            "datasource",
            choices = c("10x_Genomics", "BD_Rhapsody_Paired", "BD_Rhapsody_AIRR"),
            selected = "10x_Genomics"
          )
        } else {
          updateSelectInput(
            session,
            "datasource",
            choices = c("10x_Genomics", "BD_Rhapsody_Paired", "BD_Rhapsody_AIRR"),
            selected = "BD_Rhapsody_Paired"
          )
        }
      } else {
        updateSelectInput(
          session,
          "datasource",
          choices = c("10x_Genomics", "BD_Rhapsody_Paired", "BD_Rhapsody_AIRR"),
          selected = ""
        )
      }
    })
    ## uploading seruat obj ----
    data_sc_pro <- reactive({
      inFile_sc_pro <- input$file_SC_pro
      if (is.null(inFile_sc_pro)) {
        return(NULL)
      } else {
        dataframe <- LoadSeuratRds(inFile_sc_pro$datapath)
      }
    })
    data_sc_clusTCR_AG <- reactive({
      inFile_cluster_fileAG <- input$file_cluster_file_AG
      if (is.null(inFile_cluster_fileAG)) {
        return(NULL)
      } else {
        dataframe <- read.csv(inFile_cluster_fileAG$datapath)
      }
    })


    data_sc_clusTCR_BD <- reactive({
      inFile_cluster_fileBD <- input$file_cluster_file_BD
      if (is.null(inFile_cluster_fileBD)) {
        return(NULL)
      } else {
        dataframe <- read.csv(inFile_cluster_fileBD$datapath)
      }
    })


    data_sc_clusTCR_IgH <- reactive({
      inFile_cluster_fileIgH <- input$file_cluster_file_IgH
      if (is.null(inFile_cluster_fileIgH)) {
        return(NULL)
      } else {
        dataframe <- read.csv(inFile_cluster_fileIgH$datapath)
      }
    })

    data_sc_clusTCR_IgKL <- reactive({
      inFile_cluster_fileIgKL <- input$file_cluster_file_IgKL
      if (is.null(inFile_cluster_fileIgKL)) {
        return(NULL)
      } else {
        dataframe <- read.csv(inFile_cluster_fileIgH$datapath)
      }
    })

    ## Add in additional sample identifers - matched to Sample_Name or Orig.indent
    data_sc_Labels_to_add <- reactive({
      inFile_Labels_to_add <- input$file_Labels_to_add
      if (is.null(inFile_Labels_to_add)) {
        return(NULL)
      } else {
        dataframe <- read.csv(inFile_Labels_to_add$datapath)
      }
    })
    ## add additional sample labels to file. -----
    UMAP_metadata_with_labs <- reactive({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      req(input$Samp_col2, input$datasource, input$species_analysis, input$V_gene_sc)


      sc@meta.data$Cell_Index_old <- sc@meta.data$Cell_Index
      sc@meta.data$Cell_Index <- rownames(sc@meta.data)
      sc@meta.data$order <- 1:dim(sc@meta.data)[1]
      reduction <- (sc@reductions$umap)
      UMAP <- as.data.frame(reduction@cell.embeddings)
      names(UMAP)[1:2] <- toupper(names(UMAP)[1:2])
      UMAP$Cell_Index <- rownames(UMAP)
      meta.data <- as.data.frame(sc@meta.data)
      meta.data <- meta.data
      umap.meta <- merge(UMAP, meta.data, by = "Cell_Index", sort = F)

      if (input$add_additional_lables == "yes") {
        labs <- data_sc_Labels_to_add()
        validate(
          need(
            nrow(labs) > 0,
            "Upload lab file"
          )
        )
        names(umap.meta)[names(umap.meta) %in% input$Samp_col2] <- "ID"
        umap.meta <- merge(labs, umap.meta, by = "ID", sort = F)
        names(umap.meta)[names(umap.meta) %in% "ID"] <- input$Samp_col2
      }
      rownames(umap.meta) <- umap.meta$Cell_Index
      umap.meta <- umap.meta[order(umap.meta$order), ]
      sc@meta.data <- umap.meta
      sc
    })

    output$Sample_names_merging_sc <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      umap.meta <- sc@meta.data
      names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
      names(umap.meta)[names(umap.meta) %in% input$V_gene_sc] <- "v_gene_selected"
      if (nrow(umap.meta) > 0) {
        umap.meta
      } else {
        as.data.frame("Isssue with the file...")
      }
    })


    # checking issue with analysis ------
    output$meta.data_check_upload <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )

      umap.meta <- sc@meta.data
      if (length(input$Samp_col) > 0) {
        # message(paste("the ID column selected is: ",input$Samp_col))
        names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
      }

      if (length(input$V_gene_sc) > 0) {
        names(umap.meta)[names(umap.meta) %in% input$V_gene_sc] <- "v_gene_selected"
      }
      umap.meta
      req(TCR_Expanded())
      sc_merged <- merge(umap.meta, TCR_Expanded(), by = c("v_gene_selected", "ID_Column"), all.x = T)
      sc_merged

      if (nrow(sc_merged) > 0) {
        sc_merged
      } else {
        as.data.frame("Possible issue with Cell Index, check merging")
      }
    })
    # upload TCRex file ----
    data_sc_TCRex <- reactive({
      inupload_TCRex_file <- input$upload_TCRex_file
      if (is.null(inupload_TCRex_file)) {
        return(NULL)
      } else {
        dataframe <- read.table(inupload_TCRex_file$datapath, skip = input$skip_TCRex_up, header = T, sep = "\t")
      }
    })

    ## uploaded clusTCR table -----
    output$Tb_ClusTCR_test <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      AG_calls <- data_sc_clusTCR_AG()
      BD_calls <- data_sc_clusTCR_BD()
      validate(
        need(
          nrow(AG_calls) > 0 & nrow(BD_calls) > 0,
          "Upload clusTCR table, which is needed for TCR -> UMAP section"
        )
      )
      rbind(AG_calls, BD_calls)
    })

    output$Tb_tcrex_test <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      calls <- data_sc_TCRex()
      validate(
        need(
          nrow(calls) > 0,
          "Upload TCRex table"
        )
      )
      calls
    })

    ### observe  Samp_col V_gene_sc Samp_col2----
    observe({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      req(sc)

      df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "Samp_col",
        choices = names(df3.meta),
        selected = "Sample_Name"
      )
    })
    observe({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      req(sc)
      df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "V_gene_sc",
        choices = names(df3.meta),
        selected = "vdj_gene_cdr3_AG_BD"
      )
    })
    observe({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      req(sc)
      df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "Samp_col2",
        choices = names(df3.meta),
        selected = "Sample_Name"
      )
    })

    output$Tb_TCR_clonotypes.Umap <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 5, scrollX = TRUE), {
      calls <- select_group_metadata()
      validate(
        need(
          nrow(calls) > 0,
          error_message_val_UMAP
        )
      )
      calls
      UMAP.wt.clonality2 <- For_col_top()
      colorblind_vector <- as.data.frame(unlist(colors_UMAP_Topclonotypes()))

      if (dim(colorblind_vector)[1] == 0) {
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
        } else {

          colorblind_vector <- rep(input$one.colour.default,dim(num)[1])
        }
      }
      colorblind_vector <- as.data.frame(colorblind_vector)


      topclones_col <- as.data.frame(unique(UMAP.wt.clonality2$topclones))
      names(topclones_col) <- "topclones"
      topclones_col$col <- colorblind_vector
      topclones_col
    })


    ### summary table for publications -----

    observe({
      sc <- UMAP_metadata_with_labs()
      req(sc)
      md <- sc@meta.data
      # if (input$type.chain == 'ab') {
      updateSelectInput(
        session,
        "other_selected_summary_columns",
        choices=names(md),
        selected = c("v_gene_AG","j_gene_AG","junction_aa_AG","v_gene_BD","j_gene_BD","junction_aa_BD"))
      # }
    })

    Summary_TCR_table <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      req(sc,input$V_gene_sc,input$Samp_col,input$other_selected_summary_columns)

      md <- as.data.frame(sc@meta.data)

      md$cloneCount <- 1

      df <- as.data.frame(md)

      TCR <- df[!is.na(df[,names(df) %in% input$V_gene_sc]),]
      TCR$ID_Column <- TCR[,names(TCR) %in% input$Samp_col]

      select_cols <- c(input$V_gene_sc, input$other_selected_summary_columns, "ID_Column")
      print(select_cols)
      TCR_input <- TCR[, names(TCR) %in% c("cloneCount",select_cols)]
      TCR_input <- TCR_input %>%
        select("cloneCount","ID_Column", everything())
      # print( head(TCR_input))

      TCR_total <- as.data.frame(ddply(TCR_input,c(input$V_gene_sc),numcolwise(sum)))
      meta2.names <- names(TCR_input)

      total.condition <- as.data.frame(ddply(TCR_input, "ID_Column", numcolwise(sum)))
      # print(head(total.condition))
      emtpy <- matrix(nrow = dim(TCR_input)[1], ncol = dim(total.condition)[1])
      # print(dim(emtpy))

      for (i in 1:dim(TCR_input)[1]) {
        emtpy[i, ] <- ifelse(TCR_input$ID_Column[i] == total.condition$ID_Column[1:dim(total.condition)[1]],
                             total.condition[total.condition$ID_Column == total.condition$ID_Column[1:dim(total.condition)[1]], 2], F
        )
      }
      # as.data.frame(emtpy)
      # print(rowSums(emtpy))
      TCR_input$row_sum <- rowSums(emtpy)

      TCR_input$frequency <- TCR_input$cloneCount / TCR_input$row_sum
      TCR_input$percent <- TCR_input$frequency * 100

      df3 <- as.data.frame(ddply(TCR_input, c(select_cols,"row_sum"), numcolwise(sum)))
      df3 <- df3[order(df3$frequency, decreasing = T), ]
      df3 <- df3[order(df3$ID_Column, decreasing = T), ]
      df3$count_div_total <-paste("freq: ",df3$cloneCount,"/",df3$row_sum,sep = "")
      df3

    })

    output$Summary_TCR_tb <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 5, scrollX = TRUE), {
      Summary_TCR_table()
    })

    output$downloaddf_Summary_TCR_tb <- downloadHandler(
      filename = function() {
        paste("Summary_TCR_Table.csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(Summary_TCR_table())
        write.csv(df, file, row.names = F)
      }
    )



    ## TCR expansion ----
    TCR_Expanded <- reactive({
      sc <- UMAP_metadata_with_labs()

      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )

      req(input$V_gene_sc, input$Samp_col)

      df3.meta <- sc@meta.data
      df3.meta2 <- df3.meta[, names(df3.meta) %in% c(input$Samp_col, input$V_gene_sc)]
      df3.meta2$ID_Column <- df3.meta2[, names(df3.meta2) %in% input$Samp_col]
      df3.meta2$v_gene_selected <- df3.meta2[, names(df3.meta2) %in% input$V_gene_sc]
      # names(df3.meta2)[names(df3.meta2) %in% input$Samp_col] <- "ID_Column"
      # names(df3.meta2)[names(df3.meta2) %in% input$V_gene_sc] <- "v_gene_selected"
      # }
      df3.meta3 <- df3.meta2
      df3.meta3$v_gene_selected <- ifelse(df3.meta3$v_gene_selected == "_._", "Unknown", df3.meta3$v_gene_selected)
      df3.meta3$v_gene_selected <- ifelse(df3.meta3$v_gene_selected == "NA_NA & NA_NA", "Unknown", df3.meta3$v_gene_selected)
      df3.meta3$v_gene_selected <- ifelse(df3.meta3$v_gene_selected == "", "Unknown", df3.meta3$v_gene_selected)
      df3.meta3$v_gene_selected <- ifelse(df3.meta3$v_gene_selected == "NA", "Unknown", df3.meta3$v_gene_selected)
      df3.meta3$v_gene_selected[is.na(df3.meta3$v_gene_selected)] <- "Unknown"
      df3.meta3
      if (nrow(df3.meta3[-c(grep("Unknown", df3.meta3$v_gene_selected)), ] > 0)) {
        df3.meta3 <- df3.meta3[-c(grep("Unknown", df3.meta3$v_gene_selected)), ]
      }
      df3.meta3
      meta2.names <- names(df3.meta3)
      df3.meta3$samp.count <- 1
      total.condition <- as.data.frame(ddply(df3.meta3, "ID_Column", numcolwise(sum)))
      total.condition
      emtpy <- matrix(nrow = dim(df3.meta3)[1], ncol = dim(total.condition)[1])

      for (i in 1:dim(df3.meta3)[1]) {
        emtpy[i, ] <- ifelse(df3.meta3$ID_Column[i] == total.condition$ID_Column[1:dim(total.condition)[1]],
                             total.condition[total.condition$ID_Column == total.condition$ID_Column[1:dim(total.condition)[1]], 2], F
        )
      }
      as.data.frame(emtpy)
      #
      df3.meta3$frequency <- 1 / rowSums(emtpy)
      df3.meta3$percent <- 1 / rowSums(emtpy) * 100

      df3 <- as.data.frame(ddply(df3.meta3, meta2.names, numcolwise(sum)))
      df3 <- df3[order(df3$samp.count, decreasing = T), ]

      df4 <- df3 %>%
        mutate(Frequency_expanded = case_when(
          frequency <= 1e-4 ~ "1. Rare (0 > X < 1e-4)",
          frequency <= 0.001 ~ "2. Small (1e-4 > X <= 0.001)",
          frequency <= 0.01 ~ "3. Medium (0.001 > X <= 0.01)",
          frequency <= 0.10 ~ "4. Large (0.01 > X <= 0.1)",
          frequency <= 0.50 ~ "5. Gigantic (0.1 > X <= 0.5)",
          frequency <= 1 ~ "6. Hyperexpanded (0.5 > X <= 1)",
          TRUE ~ "Other"
        ))

      df4 <- df4 %>%
        mutate(Number_expanded = case_when(
          samp.count <= 1 ~ "1. Single (0 < X <= 1)",
          samp.count <= 5 ~ "2. Small (1 < X <= 5)",
          samp.count <= 20 ~ "3. Medium (5 < X <= 20)",
          samp.count <= 100 ~ "4. Large (20 < X <= 100)",
          samp.count <= 500 ~ "5. Hyperexpanded (100 < X <= 500)",
          TRUE ~ "6. Hyperexpanded (>500)"
        ))
      df4
    })
    ## Umap -----
    create_UMAP2 <- reactive({
      sc <- sc()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      umap.meta <- sc@meta.data
      umap.meta$ID_Column <- umap.meta[, names(umap.meta) %in% input$Samp_col]
      umap.meta$v_gene_selected <- umap.meta[, names(umap.meta) %in% input$V_gene_sc]
      umap.meta


      sc <- merge(umap.meta, TCR_Expanded(), by = c("v_gene_selected", "ID_Column"), all.x = T)

      ggplot(sc, aes(x = UMAP_1, UMAP_2, colour = seurat_clusters)) +
        geom_point() +
        scale_color_manual(values = rainbow(length(unique(sc$seurat_clusters))), na.value = input$NA_col_analysis) +
        theme_bw() +
        theme(
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 0),
          axis.title.x = element_text(colour = "black", angle = 0, vjust = .5, face = "plain", family = input$font_type, size = input$title.text.sizer2),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
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
        paste("TCR_Explore_UMAP_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_UMAP2, height = input$height_UMAP2, onefile = FALSE) # open the pdf device
        plot(create_UMAP2())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_UMAP2 <- downloadHandler(
      filename = function() {
        x <- today()
        paste("TCR_Explore_UMAP_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_UMAP2,
            height = input$height_png_UMAP2,
            res = input$resolution_PNG_UMAP2
        )
        plot(create_UMAP2())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    # clonotypes Plot -----
    cols_clonal_plot <- reactive({
      df4 <- TCR_Expanded()
      df4 <- df4[order(df4[, names(df4) %in% input$Graph_type_bar]), ]
      df_col <- unique(df4[, names(df4) %in% input$Graph_type_bar])

      num <- as.data.frame(unique(df4[, names(df4) %in% input$Graph_type_bar]))
      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep = "_"), paste(num[i, ]), col.gg[i])
        })
      } else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep = "_"), paste(num[i, ]), col.hcl[i])
        })
      } else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep = "_"), paste(num[i, ]), col.topo[i])
        })
      } else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep = "_"), paste(num[i, ]), heat_col[i])
        })
      } else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep = "_"), paste(num[i, ]), col.terrain[i])
        })
      } else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep = "_"), paste(num[i, ]), palette_rainbow[i])
        })
      } else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep = "_"), paste(num[i, ]), palette1[i])
        })
      } else if (input$colourtype == "one") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype", i, sep = "_"), paste(num[i, ]), input$one.colour.default)
        })
      } # one colour
    })
    output$myPanel_clonal_plot <- renderUI({
      cols_clonal_plot()
    })
    colors_clonal_plot <- reactive({
      df4 <- TCR_Expanded()
      df4 <- df4[order(df4[, names(df4) %in% input$Graph_type_bar]), ]
      df_col <- unique(df4[, names(df4) %in% input$Graph_type_bar])

      num <- as.data.frame(unique(df4[, names(df4) %in% input$Graph_type_bar]))
      # num <- unique(top_BD_cluster$Selected_function)

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.clonotype", i, sep = "_")]]
      })
    })
    clonal_plot <- reactive({
      df4 <- TCR_Expanded()
      df4$ID_Column <- df4[, names(df4) %in% input$Samp_col]
      message("Added ID")
      df4 <- df4[df4$ID_Column %in% input$ID_Column_factor, ]
      df4$ID_Column <- as.character(df4$ID_Column)
      df4$ID_Column <- factor(df4$ID_Column, levels = input$ID_Column_factor)
      df.col.1 <- unlist(colors_clonal_plot())
      req(df4, input$Samp_col)

      ggplot(df4, aes(x = ID_Column, y = frequency, fill = get(input$Graph_type_bar), colour = get(input$Graph_type_bar), label = get(input$Graph_type_bar))) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 50), values = alpha(df.col.1, 1), na.value = input$NA_col_analysis) +
        scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 50), values = alpha(df.col.1, 1)) +
        theme(
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
          axis.title.x = element_blank(),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )
    })

    ### top clonotypes -----
    vals <- reactiveValues(top10 = NULL)
    cols_top_clonal_plot <- reactive({
      df4 <- TCR_Expanded()
      top10 <- df4 %>%
        group_by(ID_Column) %>%
        top_n(n = input$top_no_clonotypes, wt = frequency)

      top10 <- top10[order(top10$frequency, decreasing = F), ]
      unique.top <- unique(top10$v_gene_selected)
      top10$v_gene_selected <- factor(top10$v_gene_selected,
                                      levels = unique.top
      )

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
          colourInput(paste("col.clonotype_top", i, sep = "_"), paste(num[i, ]), col.gg[i])
        })
      } else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep = "_"), paste(num[i, ]), col.hcl[i])
        })
      } else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep = "_"), paste(num[i, ]), col.topo[i])
        })
      } else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep = "_"), paste(num[i, ]), heat_col[i])
        })
      } else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep = "_"), paste(num[i, ]), col.terrain[i])
        })
      } else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep = "_"), paste(num[i, ]), palette_rainbow[i])
        })
      } else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep = "_"), paste(num[i, ]), palette1[i])
        })
      } else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.clonotype_top", i, sep = "_"), paste(num[i, ]), input$one.colour.default)
        })
      }
    })

    output$myPanel_top_clonal_plot <- renderUI({
      cols_top_clonal_plot()
    })
    colors_top_clonal_plot <- reactive({
      df4 <- TCR_Expanded()
      df4
      top10 <- df4 %>%
        group_by(ID_Column) %>%
        top_n(n = input$top_no_clonotypes, wt = frequency)

      top10 <- top10[order(top10$frequency, decreasing = F), ]
      unique.top <- unique(top10$v_gene_selected)
      top10$v_gene_selected <- factor(top10$v_gene_selected,
                                      levels = unique.top
      )

      num <- as.data.frame(unique(top10$v_gene_selected))

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.clonotype_top", i, sep = "_")]]
      })
    })

    top_clonal_plot <- reactive({
      df4 <- TCR_Expanded()
      df4
      top10 <- df4 %>%
        group_by(ID_Column) %>%
        top_n(n = input$top_no_clonotypes, wt = frequency)

      top10 <- top10[order(top10$frequency, decreasing = F), ]
      unique.top <- unique(top10$v_gene_selected)
      print(unique.top)
      top10$v_gene_selected <- factor(top10$v_gene_selected,
                                      levels = unique.top
      )

      df.col.1 <- unlist(colors_top_clonal_plot())
      unique.top$col <- as.list(df.col.1)

      top10_2 <- top10[,names(top10) %in% c("ID_Column","frequency","v_gene_selected")]

      ggplot(top10_2, aes(x = ID_Column, y = frequency, fill = v_gene_selected, label = v_gene_selected, colour = "black")) +
        geom_bar(stat = "identity") +
        theme_bw()

      scale_fill_manual(values = alpha(unique.top$col, 1), na.value = input$NA_col_analysis, labels = ~ stringr::str_wrap(.x, width = 20)) +
        scale_colour_manual(values = "black", guide = "none") +
        theme(
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
          axis.title.x = element_blank(),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )

      # ggplot(top10,aes(x=Sample_Name,y=frequency,fill=v_gene_cdr3_AB_GD, colour="Black"))+
      #   geom_bar()+
      #   scale_fill_manual(values=alpha(df.col.1, 1)) +
      #   scale_colour_manual(values = "black",guide = "none")
    })

    output$clonality.bar.graph <- renderPlot({
      if (input$Graph_type_bar == "Frequency_expanded") {
        clonal_plot()
      } else if (input$Graph_type_bar == "Number_expanded") {
        clonal_plot()
      } else {
        top_clonal_plot()
      }
    })

    output$clonality.bar.graph2 <- renderPlot({
      if (input$Graph_type_bar == "Frequency_expanded") {
        clonal_plot()
      } else if (input$Graph_type_bar == "Number_expanded") {
        clonal_plot()
      } else {
        top_clonal_plot()
      }
    })
    # Downloading the bar plot -------
    output$downloadPlot_clonaity.bar.graph <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$V_gene_sc, "_clonal_expanded_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_clonality.bar.graph, height = input$height_clonality.bar.graph, onefile = FALSE) # open the pdf device
        if (input$Graph_type_bar == "Frequency_expanded") {
          plot(clonal_plot())
        } else if (input$Graph_type_bar == "Number_expanded") {
          plot(clonal_plot())
        } else {
          plot(top_clonal_plot())
        }
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_clonaity.bar.graph <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$V_gene_sc, "_clonal_expanded_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_clonality.bar.graph,
            height = input$height_png_clonality.bar.graph,
            res = input$resolution_PNG_clonality.bar.graph
        )
        if (input$Graph_type_bar == "Clonality") {
          plot(clonal_plot())
        } else if (input$Graph_type_bar == "Number_expanded") {
          plot(clonal_plot())
        } else {
          plot(top_clonal_plot())
        }
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )


    # UMAP clonotype -> TCR -----
    UMAP.TCRclonalit <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      umap.meta <- as.data.frame(sc@meta.data)
      umap.meta
      umap.meta$ID_Column <- umap.meta[, names(umap.meta) %in% input$Samp_col]
      umap.meta$v_gene_selected <- umap.meta[, names(umap.meta) %in% input$V_gene_sc]
      TCR_Expanded <- as.data.frame(TCR_Expanded())

      UMAP.wt.clonality <- merge(umap.meta, TCR_Expanded, by = c("v_gene_selected", "ID_Column", input$Samp_col, input$V_gene_sc), all.x = T)
      UMAP.wt.clonality

      if (input$Graph_type_bar == "Number_expanded") {
        UMAP.wt.clonality$TYPE.clonality <- UMAP.wt.clonality$Number_expanded
      } else {
        UMAP.wt.clonality$TYPE.clonality <- paste(UMAP.wt.clonality$Frequency_expanded)
      }

      UMAP.wt.clonality[is.na(UMAP.wt.clonality)] <- "unknown"
      UMAP.wt.clonality <- subset(UMAP.wt.clonality, UMAP.wt.clonality$TYPE.clonality != "unknown")
      UMAP.wt.clonality
    })


    output$Tb_TCR_clonotypes.table <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      umap.meta <- sc@meta.data
      umap.meta$ID_Column <- umap.meta[, names(umap.meta) %in% input$Samp_col]
      umap.meta$v_gene_selected <- umap.meta[, names(umap.meta) %in% input$V_gene_sc]
      umap.meta
      FILE_MERGED <- merge(umap.meta, TCR_Expanded(), by = c("v_gene_selected", "ID_Column", input$Samp_col, input$V_gene_sc), all.x = T)
      FILE_MERGED$ID_Column[FILE_MERGED$v_gene_selected == "NA"] <- "unknown"
      subset(FILE_MERGED, FILE_MERGED$v_gene_selected != "unknown")
      UMAP.TCRclonalit()
    })



    cols_UMAP_clonal_plot2 <- reactive({
      UMAP.wt.clonality <- UMAP.TCRclonalit()
      UMAP.wt.clonality <- UMAP.wt.clonality[order(UMAP.wt.clonality$TYPE.clonality), ]
      UMAP.wt.clonality$TYPE.clonality <- factor(UMAP.wt.clonality$TYPE.clonality, levels = unique(UMAP.wt.clonality$TYPE.clonality))
      num <- as.data.frame(unique(UMAP.wt.clonality$TYPE.clonality))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep = "_"), paste(num[i, ]), col.gg[i])
        })
      } else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep = "_"), paste(num[i, ]), col.hcl[i])
        })
      } else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep = "_"), paste(num[i, ]), col.topo[i])
        })
      } else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep = "_"), paste(num[i, ]), heat_col[i])
        })
      } else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep = "_"), paste(num[i, ]), col.terrain[i])
        })
      } else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep = "_"), paste(num[i, ]), palette_rainbow[i])
        })
      } else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep = "_"), paste(num[i, ]), palette1[i])
        })
      } else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep = "_"), paste(num[i, ]), input$one.colour.default)
        })
      } # one colour
    })


    output$cols_UMAP_clonal_plot <- renderUI({
      cols_UMAP_clonal_plot2()
    })
    colors_UMAP_clonal_plot <- reactive({
      UMAP.wt.clonality <- UMAP.TCRclonalit()
      UMAP.wt.clonality <- UMAP.wt.clonality[order(UMAP.wt.clonality$TYPE.clonality), ]
      UMAP.wt.clonality$TYPE.clonality <- factor(UMAP.wt.clonality$TYPE.clonality, levels = unique(UMAP.wt.clonality$TYPE.clonality))
      # UMAP.wt.clonality <- UMAP.wt.clonality[order(UMAP.wt.clonality$TYPE.clonality),]
      num <- as.data.frame(unique(UMAP.wt.clonality$TYPE.clonality))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.UMAP_clonotype", i, sep = "_")]]
      })
    })

    observe({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )

      df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "RDS_V_gene_A",
        choices = names(df3.meta),
        selected = "Av_gene"
      )

      updateSelectInput(
        session,
        "RDS_V_gene_B",
        choices = names(df3.meta),
        selected = "Bv_gene"
      )

      updateSelectInput(
        session,
        "RDS_cdr3_A",
        choices = names(df3.meta),
        selected = "Acdr3"
      )
      updateSelectInput(
        session,
        "RDS_cdr3_B",
        choices = names(df3.meta),
        selected = "Bcdr3"
      )
    })



    UMAP.TCRclonalit2 <- reactive({
      UMAP.wt.clonality <- UMAP.TCRclonalit()
      UMAP.wt.clonality <- UMAP.wt.clonality[order(UMAP.wt.clonality$TYPE.clonality), ]

      colorblind_vector <- as.data.frame(unlist(colors_UMAP_clonal_plot()))

      if (dim(colorblind_vector)[1] == 0) {
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
        } else {

        }
      }
      colorblind_vector <- as.data.frame(colorblind_vector)
      names(colorblind_vector) <- "cols"
      UMAP.wt.clonality$ID_Column <- UMAP.wt.clonality[, names(UMAP.wt.clonality) %in% input$Samp_col]

      if (input$Split_by_group == "yes") {
        UMAP.wt.clonality <- UMAP.wt.clonality[UMAP.wt.clonality$ID_Column %in% input$ID_Column_factor, ]
      }

      UMAP.wt.clonality$ID_Column <- factor(UMAP.wt.clonality$ID_Column, levels = input$ID_Column_factor)

      if (input$filter_umap_expand == "yes") {
        UMAP.wt.clonality <- subset(UMAP.wt.clonality, UMAP.wt.clonality$UMAP_1 < input$UMAP_1x & UMAP.wt.clonality$UMAP_1 > input$UMAP_1y)
        UMAP.wt.clonality <- subset(UMAP.wt.clonality, UMAP.wt.clonality$UMAP_2 < input$UMAP_2x & UMAP.wt.clonality$UMAP_2 > input$UMAP_2y)
      }

      plot <- ggplot(UMAP.wt.clonality, aes(x = UMAP_1, UMAP_2, colour = TYPE.clonality, alpha = TYPE.clonality, label = TYPE.clonality)) +
        geom_point() +
        scale_color_manual(values = colorblind_vector$cols, na.value = input$NA_col_analysis, labels = ~ stringr::str_wrap(.x, width = 40)) +
        scale_alpha_manual(values = rep(1, length(unique(UMAP.wt.clonality$TYPE.clonality))), na.value = 0.5, labels = ~ stringr::str_wrap(.x, width = 40)) +
        theme_bw() +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 0),
          axis.title.x = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )


      if (input$Split_by_group == "no") {
        plot
      } else {
        plot + facet_wrap(~ID_Column, nrow = input$wrap_row)
      }
    })
    output$clonality.TCR.UMAP <- renderPlot({
      UMAP.TCRclonalit2()
    })

    output$downloadPlot_TCR.UMAP <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$V_gene_sc, "_TCR.UMAP_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_TCR.UMAP, height = input$height_TCR.UMAP, onefile = FALSE) # open the pdf device

        plot(UMAP.TCRclonalit2())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_TCR.UMAP <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$V_gene_sc, "_TCR.UMAP_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_TCR.UMAP, height = input$height_png_TCR.UMAP, res = input$resolution_PNG_TCR.UMAP)
        plot(UMAP.TCRclonalit2())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )
    # Freq ------
    observe({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "Samp_col_UMAP",
        choices = names(df3.meta),
        selected = "ID_Column"
      )
    }) # ID_Column

    For_col_top <- reactive({
      df4 <- TCR_Expanded()
      df4 <- subset(df4, df4$samp.count != 1)
      top10 <- df4 %>%
        group_by(ID_Column) %>%
        top_n(n = input$top_no_clonotypes, wt = frequency)

      top10 <- top10[order(top10$frequency, decreasing = F), ]
      unique.top <- unique(top10$v_gene_selected)
      top10$v_gene_selected <- factor(top10$v_gene_selected,
                                      levels = unique.top
      )

      top10.2 <- as.data.frame(top10$v_gene_selected)
      names(top10.2) <- "v_gene_selected"
      top10.2$ID_Column <- top10$ID_Column
      top10.2$topclones <- top10.2$v_gene_selected

      UMAP.wt.clonality <- UMAP.TCRclonalit()
      UMAP.wt.clonality2 <- merge(UMAP.wt.clonality, top10.2, by = c("v_gene_selected", "ID_Column"))
      UMAP.wt.clonality2 <- UMAP.wt.clonality2[order(UMAP.wt.clonality2$Cell_Index), ]
      UMAP.wt.clonality2 %>% mutate(across(where(is.factor), as.character)) -> UMAP.wt.clonality2
      UMAP.wt.clonality2$topclones[is.na(UMAP.wt.clonality2$topclones)] <- "not selected"
      UMAP.wt.clonality2$topclones <- factor(UMAP.wt.clonality2$topclones, levels = c("not selected", unique(as.character(top10.2$v_gene_selected))))
      UMAP.wt.clonality2 <- UMAP.wt.clonality2[order(UMAP.wt.clonality2$topclones), ]
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
          colourInput(paste("col.UMAP_clonotype_top", i, sep = "_"), paste(num[i, ]), col.gg[i])
        })
      } else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype_top", i, sep = "_"), paste(num[i, ]), col.hcl[i])
        })
      } else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype_top", i, sep = "_"), paste(num[i, ]), col.topo[i])
        })
      } else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype_top", i, sep = "_"), paste(num[i, ]), heat_col[i])
        })
      } else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype_top", i, sep = "_"), paste(num[i, ]), col.terrain[i])
        })
      } else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype_top", i, sep = "_"), paste(num[i, ]), palette_rainbow[i])
        })
      } else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype_top", i, sep = "_"), paste(num[i, ]), palette1[i])
        })
      } else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_clonotype", i, sep = "_"), paste(num[i, ]), input$one.colour.default)
        })
      } # one colour
    })

    output$cols_UMAP_Topclonotypes <- renderUI({
      cols_UMAP_Topclonotypes()
    })

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
        input[[paste("col.UMAP_clonotype_top", i, sep = "_")]]
      })
    })

    select_group_metadata <- reactive({
      sc <- UMAP_metadata_with_labs()

      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )
      df <- sc@meta.data
      df2 <- as.data.frame(unique(df[names(df) %in% input$Samp_col]))
      df2 <- as.data.frame(df2)
      df2
    })

    observe({
      df2 <- select_group_metadata()

      validate(
        need(
          nrow(df2) > 0,
          error_message_val1
        )
      )

      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- as.data.frame(df2[order(df2$V1), ])
      names(df2) <- "V1"
      df2

      updateSelectInput(
        session,
        "ID_Column_metadata",
        choices = df2$V1,
        selected = df2$V1[1]
      )
    })
    observe({
      df2 <- select_group_metadata()
      validate(
        need(
          nrow(df2) > 0,
          error_message_val1
        )
      )
      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- as.data.frame(df2[order(df2$V1), ])
      names(df2) <- "V1"
      df3 <- subset(df2, df2$V1 != "NA")

      updateSelectInput(
        session,
        "ID_Column_factor",
        choices = df2$V1,
        selected = df3$V1
      )
    })
    observe({
      df2 <- select_group_metadata()
      validate(
        need(
          nrow(df2) > 0,
          error_message_val1
        )
      )
      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- as.data.frame(df2[order(df2$V1), ])
      names(df2) <- "V1"
      df2 <- subset(df2, df2$V1 != "NA")

      updateSelectInput(
        session,
        "selected_Indiv",
        choices = df2$V1,
        selected = df2$V1[1]
      )
    })

    Topclonotypes <- reactive({
      UMAP.wt.clonality2 <- For_col_top()
      print(names(UMAP.wt.clonality2))
      UMAP.wt.clonality2 <- UMAP.wt.clonality2[UMAP.wt.clonality2$ID_Column %in% input$ID_Column_factor, ]
      UMAP.wt.clonality2$ID_Column <- factor(UMAP.wt.clonality2$ID_Column, levels = input$ID_Column_factor)

      colorblind_vector <- as.data.frame(unlist(colors_UMAP_Topclonotypes()))

      if (dim(colorblind_vector)[1] == 0) {
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
        } else {

        }
      }
      colorblind_vector <- as.data.frame(colorblind_vector)


      names(colorblind_vector) <- "col"
      topclones_col <- as.data.frame(unique(UMAP.wt.clonality2$topclones))
      names(topclones_col) <- "topclones"
      topclones_col$col <- colorblind_vector$col
      topclones_col



      if (input$display_all_samps == "yes" & input$Split_by_group == "no") {
        topclones_col
      } else if (input$display_all_samps == "yes" & input$Split_by_group == "yes") {
        topclones_col
      } else if (input$display_all_samps == "no" & input$Split_by_group == "yes") {
        UMAP.wt.clonality2 <- UMAP.wt.clonality2[UMAP.wt.clonality2$ID_Column %in% input$ID_Column_metadata, ]
        topclones_col <- topclones_col[topclones_col$topclones %in% unique(UMAP.wt.clonality2$topclones), ]
        topclones_col
      } else {
        UMAP.wt.clonality2 <- UMAP.wt.clonality2[UMAP.wt.clonality2$ID_Column %in% input$ID_Column_metadata, ]
        topclones_col <- topclones_col[topclones_col$topclones %in% unique(UMAP.wt.clonality2$topclones), ]
        topclones_col
      }

      UMAP.wt.clonality2$topclones2 <- gsub(" & ", " ", UMAP.wt.clonality2$topclones)

      # names(UMAP.wt.clonality2)[names(UMAP.wt.clonality2) %in% input$Samp_col] <- "ID_Column"
      UMAP.wt.clonality2 <- UMAP.wt.clonality2[UMAP.wt.clonality2$ID_Column %in% input$ID_Column_factor, ]
      UMAP.wt.clonality2$ID_Column <- factor(UMAP.wt.clonality2$ID_Column, levels = input$ID_Column_factor)



      plot <- ggplot(data = UMAP.wt.clonality2, aes(x = UMAP_1, UMAP_2, colour = topclones2, label = topclones2)) +
        geom_point() +
        scale_color_manual(values = topclones_col$col, na.value = input$NA_col_analysis, labels = ~ stringr::str_wrap(.x, width = 50)) +
        # scale_color_manual(values = colorblind_vector)+
        scale_alpha_manual(values = 1, na.value = 0.1, labels = ~ stringr::str_wrap(.x, width = 50)) +
        theme_bw() +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 0),
          axis.title.x = element_text(colour = "black", angle = 0, vjust = .5, face = "plain", family = input$font_type, size = input$title.text.sizer2),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )

      if (input$Split_by_group == "no") {
        plot
      } else {
        plot + facet_wrap(~ID_Column, nrow = input$wrap_row)
      }
    })

    output$clonality.TCR.UMAP.top <- renderPlot({
      Topclonotypes()
    })

    output$downloadPlot_TCR.UMAP_top <- downloadHandler(
      filename = function() {
        x <- today()
        paste("TCR_Explore_TCR.UMAP_top_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_TCR.UMAP_top, height = input$height_TCR.UMAP_top, onefile = FALSE) # open the pdf device

        plot(Topclonotypes())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_TCR.UMAP_top <- downloadHandler(
      filename = function() {
        x <- today()
        paste("TCR_Explore_TCR.UMAP_top_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_TCR.UMAP_top, height = input$height_png_TCR.UMAP_top, res = input$resolution_PNG_TCR.UMAP_top)
        plot(Topclonotypes())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )


    ### check Idents ----

    vals_meta.sc2 <- reactiveValues(AddIndents_SCobj = NULL)
    vals_meta.sc3 <- reactiveValues(AddIndents_SCobj2 = NULL)

    observeEvent(input$run_update_clust, {
      sc <- df_sc_clust()
      df.to.order <- as.data.frame(unique(Idents(df_sc_clust())))
      names(df.to.order) <- "V1"
      df.to.order2 <- as.data.frame(do.call(rbind, strsplit(as.character(df.to.order$V1), " ")))
      df.to.order2$uniqueID <- df.to.order$V1

      if (input$multiple_group_sc == "yes") {
        df.to.order2 <- df.to.order2 %>%
          arrange(V1, V2)
      } else {
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
        choices = unique.Idents$uniqueID,
        selected = unique.Idents$uniqueID[1]
      )
    }) # ident.1
    observe({
      unique.Idents <- vals_meta.sc3$AddIndents_SCobj2
      updateSelectInput(
        session,
        "unique.Idents2",
        choices = unique.Idents$uniqueID,
        selected = unique.Idents$uniqueID[2]
      )
    }) # ident.2

    DEx_sc <- reactive({
      b.interferon.response <- FindMarkers(df_sc_clust(),
                                           ident.1 = input$unique.Idents1,
                                           ident.2 = input$unique.Idents2,
                                           logfc.threshold = 0.25, # Give me ALL results
                                           min.pct = 0.25
      )
      b.interferon.response
    })

    output$DEx_table_comparison <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 5, scrollX = TRUE), {
      DEx_sc()
    })

    output$downloaddf_DEx_sc <- downloadHandler(
      filename = function() {
        paste(" Differential expresion of ", input$unique.Idents1, " vs. ", input$unique.Idents2, gsub("-", ".", Sys.Date()), ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(DEx_sc())
        write.csv(df, file, row.names = T)
      }
    )


    output$downloaddf_DEx_sc_ggVolcanoR <- downloadHandler(
      filename = function() {
        paste(" Differential expresion of ", input$unique.Idents1, " vs. ", input$unique.Idents2, gsub("-", ".", Sys.Date()), ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(DEx_sc())
        df1 <- df
        df1$ID <- rownames(df1)
        df1$Pvalue <- df1$p_val_adj
        df1$logFC <- df1$avg_log2FC
        df2 <- df1[, names(df1) %in% c("ID", "logFC", "Pvalue")]
        write.csv(df2, file, row.names = F)
      }
    )

    ggplot_DEx_sc <- reactive({
      dat <- DEx_sc()
      dat <- as.data.frame(dat)
      dat <- dat[order(dat$p_val_adj), ]
      dat$ID <- rownames(dat)
      neg <- -1 * 0.58
      pos <- 0.58

      maximum <- input$max

      mutateddf <- mutate(dat, sig = ifelse(dat$p_val_adj < 0.05, "Pvalue<0.05", "Not Sig"))
      sig <- subset(dat, dat$p_val_adj < 0.05 & abs(dat$avg_log2FC) > 0.58)
      top <- sig[(1:30), ]

      gene_list <- top$ID

      mutateddf.gene <- mutate(mutateddf, top = ifelse(mutateddf$ID %in% gene_list, "top", "other"))
      mutateddf.gene

      # no labels -----
      sub.mutateddf.gene <- mutate(mutateddf.gene,
                                   colour = ifelse(mutateddf.gene$p_val_adj < 0.05 & mutateddf.gene$avg_log2FC > pos, "sig_up",
                                                   ifelse(mutateddf.gene$p_val_adj < 0.05 & mutateddf.gene$avg_log2FC < neg, "sig_down", "NS")
                                   ),
                                   alpha = ifelse(mutateddf.gene$p_val_adj < 0.05 & mutateddf.gene$avg_log2FC > pos, 0.25,
                                                  ifelse(mutateddf.gene$p_val_adj < 0.05 & mutateddf.gene$avg_log2FC < neg, 0.25, 0.1)
                                   ),
                                   shape = ifelse(mutateddf.gene$p_val_adj < 0.05 & mutateddf.gene$avg_log2FC > pos, 19,
                                                  ifelse(mutateddf.gene$p_val_adj < 0.05 & mutateddf.gene$avg_log2FC < neg, 19, 1)
                                   ),
                                   size = ifelse(mutateddf.gene$p_val_adj < 0.05 & mutateddf.gene$avg_log2FC > pos, 3,
                                                 ifelse(mutateddf.gene$p_val_adj < 0.05 & mutateddf.gene$avg_log2FC < neg, 3, 1)
                                   )
      )
      # range of genes -----
      sub.mutateddf.gene2 <- mutate(mutateddf.gene,
                                    colour = ifelse(mutateddf.gene$ID %in% gene_list & mutateddf.gene$avg_log2FC > pos & mutateddf.gene$p_val_adj < 0.05, "Labelled_up",
                                                    ifelse(mutateddf.gene$ID %in% gene_list & mutateddf.gene$avg_log2FC < neg & mutateddf.gene$p_val_adj < 0.05, "Labelled_down", ifelse(mutateddf.gene$p_val_adj < 0.05 & mutateddf.gene$avg_log2FC > pos, "Significant-up",
                                                                                                                                                                                         ifelse(mutateddf.gene$p_val_adj < 0.05 & mutateddf.gene$avg_log2FC < neg, "Significant-down", "Non-significant")
                                                    ))
                                    ),
                                    size = ifelse(mutateddf.gene$ID %in% gene_list & mutateddf.gene$avg_log2FC > pos & mutateddf.gene$p_val_adj < 0.05, 3,
                                                  ifelse(mutateddf.gene$ID %in% gene_list & mutateddf.gene$avg_log2FC < neg & mutateddf.gene$p_val_adj < 0.05, 3, ifelse(mutateddf.gene$p_val_adj < 0.05 & mutateddf.gene$avg_log2FC > pos, 3,
                                                                                                                                                                         ifelse(mutateddf.gene$p_val_adj < 0.05 & mutateddf.gene$avg_log2FC < neg, 3, 1)
                                                  ))
                                    ),
                                    shape = ifelse(mutateddf.gene$ID %in% gene_list & mutateddf.gene$avg_log2FC > pos & mutateddf.gene$p_val_adj < 0.05, 19,
                                                   ifelse(mutateddf.gene$ID %in% gene_list & mutateddf.gene$avg_log2FC < neg & mutateddf.gene$p_val_adj < 0.05, 19, ifelse(mutateddf.gene$p_val_adj < 0.05 & mutateddf.gene$avg_log2FC > pos, 15,
                                                                                                                                                                           ifelse(mutateddf.gene$p_val_adj < 0.05 & mutateddf.gene$avg_log2FC < neg, 15, 1)
                                                   ))
                                    ),
      )

      y_lable1 <- bquote("-" ~ Log[10] ~ (.("p_val_adj")))
      y_lable1

      colour_class3 <- c("Significant-down", "Significant-up", "Labelled_down", "Labelled_up", "labelled-Non-significant", "Non-significant")
      colour.df3 <- as.data.frame(c("Significant-down", "Significant-up", "Labelled_down", "Labelled_up", "labelled-Non-significant", "Non-significant"))
      names(colour.df3) <- "label"
      # colour.df3$V1 <- c(input$down,input$up,input$col_lab1,input$col_lab2,input$col_lab3,input$NS)
      # colour.df3$shape <- c(input$shape2,input$shape1.1,input$shape1,input$shape1,input$shape1,input$shape3)
      # colour.df3$size <- c(input$size2,input$size1.1,input$size1,input$size1,input$size1,input$size3)
      # colour.df3$alpha <- c(input$alpha2,input$alpha2,input$alpha1,input$alpha1,input$alpha1,input$alpha3)
      # colour.class4 <- colour.df3[colour.df3$label %in% unique(sub.mutateddf.gene2$colour),]
      # sub.mutateddf.gene2$colour <- factor(sub.mutateddf.gene2$colour, levels = colour.class4$label)

      ggplot() +
        geom_point(aes(
          x = sub.mutateddf.gene2$avg_log2FC, y = -log10(sub.mutateddf.gene2$p_val_adj),
          col = sub.mutateddf.gene2$colour,
          # shape=sub.mutateddf.gene2$shape,
          # alpha=sub.mutateddf.gene2$colour,
          size = sub.mutateddf.gene2$size
        ), ) +
        # scale_color_manual(name="legend",values=colour.class4$V1, labels = colour.class4$label) +
        # scale_shape_manual(name="legend",values=colour.class4$shape, labels=colour.class4$label)+
        # scale_size_manual(name="legend",values=colour.class4$size, labels=colour.class4$label)+
        # scale_alpha_manual(name="legend",values=colour.class4$alpha, labels=colour.class4$label) +

        geom_text_repel(
          data = sub.mutateddf.gene2[sub.mutateddf.gene2$ID %in% gene_list, ],
          aes(
            x = sub.mutateddf.gene2$avg_log2FC[sub.mutateddf.gene2$ID %in% gene_list],
            y = -log10(sub.mutateddf.gene2$p_val_adj)[sub.mutateddf.gene2$ID %in% gene_list],
            label = sub.mutateddf.gene2$ID[sub.mutateddf.gene2$ID %in% gene_list]
          ),
          size = 6,
          family = "serif",
          segment.alpha = 0.5,
          show.legend = F, box.padding = unit(1, "lines"),
          max.overlaps = Inf
        ) +
        guides(shape = guide_legend(override.aes = list(size = 5))) +
        theme_bw(base_size = 18) +
        theme(
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")
        ) +
        geom_vline(xintercept = pos, linetype = "dashed", color = "darkorange") +
        geom_vline(xintercept = neg, linetype = "dashed", color = "darkorange") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkorange") +
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
      labs(
        y = y_lable1,
        x = expression(Log[2] ~ Fold ~ Change),
        title = ""
      )
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
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      meta.data <- sc@meta.data
      totals <- meta.data[, names(meta.data) %in% c(input$Samp_col, input$Colour_By_this_overview)]

      totals <- totals %>%
        select(all_of(input$Samp_col), everything())

      names(totals) <- c("groups", "Function")
      tab <- table(totals$Function, totals$groups)
      as.matrix(t(tab) / colSums(tab) * 100)
    })

    # output$Percent_tab <- renderPrint(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
    output$Percent_tab <- renderPrint({
      percent <- Percent_tab_df()

      round(percent,4)

    })
    output$downloaddf_Percent_tab <- downloadHandler(
      filename = function() {
        x <- today()
        paste("GEx_percent_", input$Colour_By_this_overview, "_", x, ".tsv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(Percent_tab_df())
        write.table(df, file, row.names = F, sep = "\t", quote = F)
      }
    )

    output$SiteNumInput <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )

      df_class <- sc@meta.data
      df_class$selected <- df_class[, names(df_class) %in% input$Colour_By_this_overview]
      selectInput("SiteNumInput", "Show on graph", choices = unique(df_class$selected), selected = NULL, multiple = TRUE)
    })

    # Umap classification plot -------
    cols_UMAP_all_classification <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )

      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this_overview]
      top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function), NA, top_BD_cluster$Selected_function)
      top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function, levels = unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep = "_"), paste(num[i, ]), col.gg[i])
        })
      } else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep = "_"), paste(num[i, ]), col.hcl[i])
        })
      } else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep = "_"), paste(num[i, ]), col.topo[i])
        })
      } else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep = "_"), paste(num[i, ]), heat_col[i])
        })
      } else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep = "_"), paste(num[i, ]), col.terrain[i])
        })
      } else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep = "_"), paste(num[i, ]), palette_rainbow[i])
        })
      } else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep = "_"), paste(num[i, ]), palette1[i])
        })
      } else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_all_classification", i, sep = "_"), paste(num[i, ]), input$one.colour.default)
        })
      } # one colour
    })

    output$cols_UMAP_all_classification <- renderUI({
      cols_UMAP_all_classification()
    })

    colors_UMAP_all_classification <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val1
        )
      )
      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this_overview]
      top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function), NA, top_BD_cluster$Selected_function)
      top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function, levels = unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.UMAP_all_classification", i, sep = "_")]]
      })
    })

    UMAP_all_classification <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val1
        )
      )
      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this_overview]
      top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function), NA, top_BD_cluster$Selected_function)
      # top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
      col.file <- as.data.frame(unique(top_BD_cluster$Selected_function))
      col.file <- as.data.frame(col.file[complete.cases(col.file) == T, ])
      names(col.file) <- "V1"
      if (input$show_selected == "Selected_list") {
        col.file$col <- unlist(colors_UMAP_all_classification())
        col.file <- col.file[col.file$V1 %in% input$SiteNumInput, ]
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function, levels = col.file$V1)
      } else {
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function, levels = unique(top_BD_cluster$Selected_function))

        col.file$col <- unlist(colors_UMAP_all_classification())
        col.file
      }

      names(top_BD_cluster)[names(top_BD_cluster) %in% input$Samp_col] <- "ID_Column"

      if (input$by_indiv_pie_epi == "yes") {
        top_BD_cluster <- top_BD_cluster[top_BD_cluster$ID_Column %in% input$selected_Indiv, ]
      }

      # top_BD_cluster$alpha_val <- ifelse(is.na(top_BD_cluster$Selected_function)==T,0.1,1)


      top_BD_cluster <- top_BD_cluster[top_BD_cluster$ID_Column %in% input$ID_Column_factor, ]
      top_BD_cluster$ID_Column <- factor(top_BD_cluster$ID_Column, levels = input$ID_Column_factor)
      md <- top_BD_cluster
      md <- subset(md, md$UMAP_1 > input$Filter_lower_UMAP1_marker_GEX)
      md <- subset(md, md$UMAP_1 < input$Filter_lower_UMAP1_marker2_GEX)

      md <- subset(md, md$UMAP_2 > input$Filter_lower_UMAP2_marker_GEX)
      md <- subset(md, md$UMAP_2 < input$Filter_lower_UMAP2_marker2_GEX)
      top_BD_cluster <- md

      df <- ggplot(top_BD_cluster, aes(x = UMAP_1, UMAP_2, colour = Selected_function, alpha = Selected_function, label = Selected_function)) +
        geom_point() +
        scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 20), values = col.file$col, na.value = input$NA_col_analysis) +
        scale_alpha_manual(labels = ~ stringr::str_wrap(.x, width = 20), values = rep(1, length(unique(top_BD_cluster$Selected_function))), na.value = 0.1) +
        theme_bw() +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 0),
          axis.title.x = element_text(colour = "black", angle = 0, vjust = .5, face = "plain", family = input$font_type, size = input$title.text.sizer2),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )

      if (input$Split_by_group == "no") {
        df <- df
      } else {
        df <- df + facet_wrap(~ID_Column, nrow = input$wrap_row)
      }
      df
    })

    output$UMAP_all_classification2 <- renderPlot({
      UMAP_all_classification()
    })

    output$downloadPlot_UMAP_all_classification <- downloadHandler(
      filename = function() {
        x <- today()
        paste("UMAP_GEx_", input$Colour_By_this_overview, "_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_UMAP_all_classification, height = input$height_UMAP_all_classification, onefile = FALSE) # open the pdf device
        plot(UMAP_all_classification())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_UMAP_all_classification <- downloadHandler(
      filename = function() {
        x <- today()
        paste("UMAP_GEx_", input$Colour_By_this_overview, "_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_UMAP_all_classification,
            height = input$height_png_UMAP_all_classification,
            res = input$resolution_PNG_UMAP_all_classification
        )
        plot(UMAP_all_classification())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    #### Pie chart -----

    observe({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val1
        )
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
        choices = names(meta.data),
        selected = "orig.ident"
      )
      # }
    })
    # pie colouring  ----
    cols_pie <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val1
        )
      )
      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this_overview]
      top_BD_cluster$Selected_function <- ifelse(top_BD_cluster$Selected_function == "NA", "-", top_BD_cluster$Selected_function)
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function), ]

      # top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep = "_"), paste(num[i, ]), col.gg[i])
        })
      } else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep = "_"), paste(num[i, ]), col.hcl[i])
        })
      } else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep = "_"), paste(num[i, ]), col.topo[i])
        })
      } else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep = "_"), paste(num[i, ]), heat_col[i])
        })
      } else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep = "_"), paste(num[i, ]), col.terrain[i])
        })
      } else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep = "_"), paste(num[i, ]), palette_rainbow[i])
        })
      } else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep = "_"), paste(num[i, ]), palette1[i])
        })
      } else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep = "_"), paste(num[i, ]), input$one.colour.default)
        })
      } # one colour
    })
    output$myPanel_pie <- renderUI({
      cols_pie()
    })
    colors_pie <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val1
        )
      )
      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this_overview]
      top_BD_cluster$Selected_function <- ifelse(top_BD_cluster$Selected_function == "NA", "-", top_BD_cluster$Selected_function)
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function), ]
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.pie", i, sep = "_")]]
      })
    })

    output$table_pie <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val1
        )
      )
      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this_overview]
      top_BD_cluster$Selected_function <- ifelse(top_BD_cluster$Selected_function == "NA", "-", top_BD_cluster$Selected_function)
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function), ]
      top_BD_cluster$Selected_group <- top_BD_cluster[, names(top_BD_cluster) %in% input$Split_group_by_]
      top_BD_cluster$Selected_group <- factor(top_BD_cluster$Selected_group, level = input$Graph_split_order)
      top_BD_cluster$Selected_indiv <- top_BD_cluster[, names(top_BD_cluster) %in% input$Samp_col]
      top_BD_cluster
      df.col <- unlist(colors_pie())

      if (input$by_indiv_pie_epi == "yes") {
        top_BD_cluster <- top_BD_cluster[top_BD_cluster$Selected_indiv %in% input$selected_Indiv, ]
      }

      top_BD_cluster
    })

    Pie_chart_Class <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val1
        )
      )
      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this_overview]
      top_BD_cluster$Selected_function <- ifelse(top_BD_cluster$Selected_function == "NA", "-", top_BD_cluster$Selected_function)
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function), ]
      top_BD_cluster$Selected_group <- top_BD_cluster[, names(top_BD_cluster) %in% input$Split_group_by_]

      df.col <- unlist(colors_pie())
      names(top_BD_cluster)[names(top_BD_cluster) %in% input$Samp_col] <- "ID_Column"

      if (input$by_indiv_pie_epi == "yes") {
        top_BD_cluster <- top_BD_cluster[top_BD_cluster$ID_Column %in% input$selected_Indiv, ]
      }

      df3.meta3 <- as.data.frame(table(top_BD_cluster$Selected_group, top_BD_cluster$Selected_function))
      total.condition <- as.data.frame(ddply(df3.meta3, "Var1", numcolwise(sum)))
      emtpy <- matrix(nrow = dim(df3.meta3)[1], ncol = dim(total.condition)[1])

      for (j in 1:dim(total.condition)[1]) {
        for (i in 1:dim(df3.meta3)[1]) {
          emtpy[i, j] <- ifelse(df3.meta3$Var1[i] == total.condition$Var1[j], total.condition$Freq[j], F)
        }
      }

      df3.meta3$n <- df3.meta3$Freq / rowSums(emtpy)

      ggplot(df3.meta3, aes(x = "", y = n, fill = Var2, group = Var1)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        theme_void(20) +
        facet_wrap(~Var1, nrow = input$wrap_row) +
        theme(
          legend.key.size = unit(1, "cm")
        ) +
        scale_fill_manual(values = df.col, na.value = input$NA_col_analysis) +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          legend.text = element_text(size = input$Legend_size, family = input$font_type),
          legend.title = element_blank()
        )

      # theme(strip.background =element_rect(fill=input$strip.colour.tree))+
      # theme(strip.text = element_text(colour = input$strip.text.colour.tree))
      # scale_pattern_fill_manual(values = class_col.t)
    })

    output$Classification_clonotype_pie <- renderPlot({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val1
        )
      )
      Pie_chart_Class()
    })
    output$downloadPlot_Classification_clonotype_pie <- downloadHandler(
      filename = function() {
        x <- today()
        paste("GEx_pie_", input$Colour_By_this_overview, "_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_Classification_clonotype_pie, height = input$height_Classification_clonotype_pie, onefile = FALSE) # open the pdf device
        plot(Pie_chart_Class())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_Classification_clonotype_pie <- downloadHandler(
      filename = function() {
        x <- today()
        paste("GEx_pie_", input$Colour_By_this_overview, "_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_Classification_clonotype_pie,
            height = input$height_png_Classification_clonotype_pie,
            res = input$resolution_PNG_Classification_clonotype_pie
        )
        plot(Pie_chart_Class())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    # chi-square of expression -----
    chi_squ <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      meta.data <- sc@meta.data
      totals <- meta.data[, names(meta.data) %in% c(input$Samp_col, input$Split_group_by_)]

      totals %>%
        select(all_of(input$Samp_col), everything())

      names(totals) <- c("samps", "split")
      totals$split <- as.character(totals$split)
      # totals <- totals[!totals$Function %in% c("NA"),]
      # totals <- totals[!totals$groups %in% c("NA"),]
      # totals <- totals[totals$groups %in% input$ID_Column_factor,]
      # totals$samps <- factor(totals$samps,levels = c(input$ID_Column_factor))
      totals
    })

    output$Chi_tab_before <- renderPrint({
      chi_squ()
      totals <- chi_squ()
      tb_totals <- table(totals$samps, totals$split)
      as.data.frame(tb_totals)
      chisq <- chisq.test(tb_totals)
      chisq
    })

    output$Post_hoc_chi <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      totals <- chi_squ()
      tb_totals <- table(totals$samps, totals$split)
      df <- tb_totals
      posthoc <- chisq.posthoc.test(df)
      resid <- subset(posthoc, posthoc$Value == "Residuals")
      pval <- subset(posthoc, posthoc$Value == "p values")
      as.data.frame(pval)
    })

    Chi_square_plot2 <- reactive({
      totals <- as.data.frame(chi_squ())
      totals <- totals[totals$samps %in% input$ID_Column_factor, ]
      tb_totals <- table(totals$samps, totals$split)
      chisq <- chisq.test(tb_totals)
      chisq
      res <- chisq$residuals
      res <- setNames(melt(res), c("y", "x", "residuals"))


      plot_chi <- ggplot(res, aes(x = x, y = y, size = residuals, fill = residuals))

      plot_chi <- plot_chi +
        geom_point(shape = 21, colour = "black") +
        # scale_size_area(max_size = 10) +
        scale_fill_gradient2(
          low = input$lower_col_chi,
          mid = input$mid_col_chi,
          high = input$high_col_chi,
          space = "Lab",
          na.value = "grey90",
          guide = "colourbar",
          aesthetics = "fill"
        ) +
        theme_bw() +
        theme(
          axis.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          axis.text.x = element_text(colour = "black", size = input$Legend_size, family = input$font_type, angle = 90),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          axis.line.x = element_blank(),
          panel.grid = element_blank(),
          legend.position = "right",
        ) +
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
        pdf(file, width = input$width_Chi_square_plot, height = input$height_Chi_square_plot, onefile = FALSE) # open the pdf device
        plot(Chi_square_plot2())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_Chi_square_plot <- downloadHandler(
      filename = function() {
        paste("_Chi_square_plot_", ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Chi_square_plot, height = input$height_png_Chi_square_plot, res = input$resolution_PNG_Chi_square_plot)
        plot(Chi_square_plot2())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )


    # Select top clonotype ----
    Top_clonotype_df2 <- reactive({
      sc <- UMAP_metadata_with_labs()
      req(input$selected_Indiv, input$Samp_col, input$V_gene_sc)
      validate(
        need(
          nrow(sc) > 0,
          error_message_val1
        )
      )
      dataframe_one <- sc@meta.data
      names(dataframe_one)[names(dataframe_one) %in% input$Samp_col] <- "ID_Column"

      if (input$by_indiv_pie_epi == "yes") {
        df3.meta <- dataframe_one[dataframe_one$ID_Column %in% input$selected_Indiv, ]
      } else {
        df3.meta <- dataframe_one
      }

      # summarising the clonality
      df3.meta$cluster_name <- df3.meta[, names(df3.meta) %in% input$V_gene_sc]
      df3.meta$cloneCount <- 1

      BD <- df3.meta[, names(df3.meta) %in% c("cluster_name", "cloneCount")]
      BD
      BD_sum <- ddply(BD, names(BD)[-c(2)], numcolwise(sum))
      as.data.frame(BD_sum)
      BD_sum <- BD_sum[!BD_sum$cluster_name == "_", ]
      BD_sum <- BD_sum[!BD_sum$cluster_name %in% "NA", ]
      names(BD_sum)[2] <- "Total_count"
      BD_sum$frequency <- BD_sum$Total_count / sum(BD_sum[, c("Total_count")], na.rm = T)

      # BD_sum$frequency <- BD_sum$Total_count/sum(BD_sum$Total_count)
      BD_sum <- BD_sum[order(BD_sum$Total_count, decreasing = T), ]

      BD_sum
    })

    output$Top_clonotype_sum <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      Top_clonotype_df2()
    })


    output$download_Top_clonotype_sum <- downloadHandler(
      filename = function() {
        x <- today()
        paste("clonotype_summary_table_", x, ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(Top_clonotype_df2())
        write.csv(df, file, row.names = F)
      }
    )

    # top clonotypes observe events -----
    observe({
      BD_sum <- Top_clonotype_df2()
      validate(
        need(
          nrow(BD_sum) > 0,
          error_message_val1
        )
      )

      # if (input$analysis_set.type)
      BD_sum <- subset(BD_sum, BD_sum$Total_count > 2)
      # if(input$restict_vGene) {
      #   BD_sum$Topclone <- BD_sum[,names(BD_sum) %in% input$v_gene_top_clone]
      #   BD_sum <- subset(BD_sum,BD_sum$Topclone == input$top_clone_selected)
      # }
      if (input$limit_to_top_clones) {
        BD_sum <- BD_sum[1:input$max_top_clone_limit, ]
      }

      updateSelectInput(
        session,
        "Selected_clonotype",
        choices = BD_sum$cluster_name,
        selected = BD_sum$cluster_name[1]
      )
    })





    select_split_by <- reactive({
      sc <- UMAP_metadata_with_labs()
      # df <- top_clonotype_bar_code()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
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
        need(
          nrow(df2) > 0,
          error_message_val1
        )
      )
      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- (df2[order(df2$V1), ])

      updateSelectInput(
        session,
        "Graph_split_order",
        choices = df2,
        selected = df2
      )
    })

    # top clonptype bar graph -----
    top_clonotype_bar_code <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload"
        )
      )
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[, names(df3.meta) %in% input$V_gene_sc]
      top_BD_clonotype <- df3.meta[df3.meta$cluster_name %in% input$Selected_clonotype, ]
      top_BD_clonotype
    })

    ##### bar graph of top clonotypes -----
    cols_Top_bar_clonotype <- reactive({
      dtop_clonotype_bar_code <- top_clonotype_bar_code()
      dtop_clonotype_bar_code$Selected_chain <- dtop_clonotype_bar_code[, names(dtop_clonotype_bar_code) %in% input$V_gene_sc]

      # top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(unique(dtop_clonotype_bar_code$Selected_chain))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")
      one <- rep(input$one.colour.default,dim(num)[1])

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_bar_clonotype", i, sep = "_"), paste(num[i, ]), col.gg[i])
        })
      } else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_bar_clonotype", i, sep = "_"), paste(num[i, ]), col.hcl[i])
        })
      } else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_bar_clonotype", i, sep = "_"), paste(num[i, ]), col.topo[i])
        })
      } else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_bar_clonotype", i, sep = "_"), paste(num[i, ]), heat_col[i])
        })
      } else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_bar_clonotype", i, sep = "_"), paste(num[i, ]), col.terrain[i])
        })
      } else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_bar_clonotype", i, sep = "_"), paste(num[i, ]), palette_rainbow[i])
        })
      } else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_bar_clonotype", i, sep = "_"), paste(num[i, ]), palette1[i])
        })
      } else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.pie", i, sep = "_"), paste(num[i, ]), one[i])
        })
      } # one colour
    })
    output$myPanel_Top_bar_clonotype <- renderUI({
      cols_Top_bar_clonotype()
    })
    colors_cols_Top_bar_clonotype <- reactive({
      dtop_clonotype_bar_code <- top_clonotype_bar_code()
      dtop_clonotype_bar_code$Selected_chain <- dtop_clonotype_bar_code[, names(dtop_clonotype_bar_code) %in% input$V_gene_sc]
      # top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(unique(dtop_clonotype_bar_code$Selected_chain))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.Top_bar_clonotype", i, sep = "_")]]
      })
    })

    ggplot_top_BD_clonotype_vs_SC <- reactive({
      dtop_clonotype_bar_code <- top_clonotype_bar_code()

      req(input$Graph_split_order)

      dtop_clonotype_bar_code$Selected_group <- dtop_clonotype_bar_code[, names(dtop_clonotype_bar_code) %in% input$Split_group_by_]
      colorblind_vector <- as.data.frame(unlist(colors_cols_Top_bar_clonotype()))

      if (dim(colorblind_vector)[1] == 0) {
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
        } else {

        }
      }
      colorblind_vector <- as.data.frame(colorblind_vector)
      names(colorblind_vector) <- "cols"

      dtop_clonotype_bar_code$Selected_chain2 <- dtop_clonotype_bar_code[, names(dtop_clonotype_bar_code) %in% input$V_gene_sc]
      dtop_clonotype_bar_code$Selected_chain3 <- gsub("_", " ", dtop_clonotype_bar_code$Selected_chain2)
      dtop_clonotype_bar_code$Selected_chain3 <- gsub("[.]", " ", dtop_clonotype_bar_code$Selected_chain3)

      dtop_clonotype_bar_code <- dtop_clonotype_bar_code[dtop_clonotype_bar_code$Selected_group %in% input$Graph_split_order, ]
      dtop_clonotype_bar_code$Selected_group <- factor(dtop_clonotype_bar_code$Selected_group, levels = input$Graph_split_order)

      ggplot(dtop_clonotype_bar_code, aes(x = Selected_group, fill = Selected_chain3, colour = Selected_chain3, label = Selected_chain3)) +
        geom_bar() +
        theme_bw() +
        scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector$cols, na.value = input$NA_col_analysis) +
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector$cols, na.value = input$NA_col_analysis) +
        # scale_alpha_manual(values = rep(1,length(unique(dtop_clonotype_bar_code$Selected_chain))), na.value=0.5)+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
          axis.title.x = element_blank(),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
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
        paste(input$Selected_clonotype, "top_clonotype_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_top_clonotype, height = input$height_top_clonotype, onefile = FALSE) # open the pdf device
        plot(ggplot_top_BD_clonotype_vs_SC())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_top_clonotype <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Selected_clonotype, "top_clonotype_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_top_clonotype,
            height = input$height_png_top_clonotype,
            res = input$resolution_PNG_top_clonotype
        )
        plot(ggplot_top_BD_clonotype_vs_SC())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    #### top clonotype pie -----
    cols_Top_pie_clonotype <- reactive({
      top_BD_cluster <- top_clonotype_bar_code()
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this]
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function), ]
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this]
      top_BD_cluster$Selected_group <- top_BD_cluster[, names(top_BD_cluster) %in% input$Split_group_by_]
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
          colourInput(paste("col.Top_pie_clonotype", i, sep = "_"), paste(num[i, ]), col.gg[i])
        })
      } else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep = "_"), paste(num[i, ]), col.hcl[i])
        })
      } else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep = "_"), paste(num[i, ]), col.topo[i])
        })
      } else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep = "_"), paste(num[i, ]), heat_col[i])
        })
      } else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep = "_"), paste(num[i, ]), col.terrain[i])
        })
      } else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep = "_"), paste(num[i, ]), palette_rainbow[i])
        })
      } else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep = "_"), paste(num[i, ]), palette1[i])
        })
      } else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.Top_pie_clonotype", i, sep = "_"), paste(num[i, ]), input$one.colour.default)
        })
      } # one colour
    })
    output$myPanel_Top_pie_clonotype <- renderUI({
      cols_Top_pie_clonotype()
    })
    colors_cols_Top_pie_clonotype <- reactive({
      top_BD_cluster <- top_clonotype_bar_code()
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this]
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function), ]
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this]
      top_BD_cluster$Selected_group <- top_BD_cluster[, names(top_BD_cluster) %in% input$Split_group_by_]
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[!num == "NA"])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.Top_pie_clonotype", i, sep = "_")]]
      })
    })

    pie_ag_bd_table <- reactive({
      top_BD_cluster <- top_clonotype_bar_code()
      top_BD_cluster
      # req(input$Colour_By_this,input$Split_group_by_,input$Graph_split_order)

      top_BD_cluster$Selected_function <- as.character((top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this]))
      top_BD_cluster$Selected_group <- as.character((top_BD_cluster[, names(top_BD_cluster) %in% input$Split_group_by_]))
      top_BD_cluster <- top_BD_cluster[top_BD_cluster$Selected_group %in% input$Graph_split_order, ]
      top_BD_cluster$Selected_group <- factor(top_BD_cluster$Selected_group, levels = input$Graph_split_order)

      df3.meta3 <- top_BD_cluster[, names(top_BD_cluster) %in% c("Selected_function", "Selected_group")]
      df3.meta3$count <- 1
      df3.meta3
      total.condition <- ddply(df3.meta3, "Selected_group", numcolwise(sum))
      total.condition
      #
      total.group.condition <- ddply(df3.meta3, c("Selected_group", "Selected_function"), numcolwise(sum))
      total.group.condition
      #
      nrow_count <- dim(total.group.condition)[1]
      ncol_count <- dim(total.condition)[1]
      #
      emtpy_mat <- matrix(nrow = nrow_count, ncol = ncol_count)

      for (i in 1:ncol_count) {
        emtpy_mat[, i] <- (ifelse(total.group.condition$Selected_group == total.condition$Selected_group[i], total.condition$count[i], F))
      }
      as.data.frame(emtpy_mat)
      #
      total.group.condition$n <- total.group.condition$count / rowSums(emtpy_mat)
      total.group.condition$Selected_function <- ifelse(total.group.condition$Selected_function == "NA", NA, total.group.condition$Selected_function)
      total.group.condition
    })

    output$Top_clonotype_Labs <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      pie_ag_bd_table()
    })

    Pie_chart_alpha_gamma <- reactive({
      df.col <- unlist(colors_cols_Top_pie_clonotype())

      total.group.condition <- pie_ag_bd_table()


      total.group.condition$Selected_function <- gsub("_", " ", total.group.condition$Selected_function)
      total.group.condition$Selected_function <- gsub("[.]", " ", total.group.condition$Selected_function)

      # ggplot(dtop_clonotype_bar_code, aes(x=Selected_group, fill=Selected_chain3,colour = Selected_chain3, label = Selected_chain3)) +
      #   geom_bar() +
      #   theme_bw()+
      #   scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector, na.value=input$NA_col_analysis)+
      #   scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector, na.value=input$NA_col_analysis)+

      ggplot(total.group.condition, aes(x = "", y = n, fill = as.character(Selected_function), group = as.character(Selected_group), label = Selected_function)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        theme_void(20) +
        facet_wrap(~Selected_group, nrow = input$wrap_row) +
        theme(
          legend.key.size = unit(1, "cm"),
          legend.title = element_blank()
        ) +
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = df.col, na.value = input$NA_col_analysis) +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.position = input$legend_position,
          legend.title = element_blank()
        )
    })

    output$top_clonotype_pie <- renderPlot({
      if (input$Plot_type_selected == "pie") {
        Pie_chart_alpha_gamma()
      } else {
        UMAP_chart_alpha_gamma()
      }
    })

    #### top  UMAP -----
    UMAP_chart_alpha_gamma <- reactive({
      top_BD_cluster <- top_clonotype_bar_code()
      req(input$Colour_By_this, input$Split_group_by_, input$Graph_split_order)
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this]
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function), ]

      top_BD_cluster$Selected_function <- as.character(top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this])

      top_BD_cluster$Selected_group <- as.character(top_BD_cluster[, names(top_BD_cluster) %in% input$Split_group_by_])
      #
      top_BD_cluster <- top_BD_cluster[top_BD_cluster$Selected_group %in% input$Graph_split_order, ]
      top_BD_cluster$Selected_group <- factor(top_BD_cluster$Selected_group, levels = input$Graph_split_order)
      top_BD_cluster$Selected_function <- ifelse(top_BD_cluster$Selected_function == "NA", NA, top_BD_cluster$Selected_function)

      df.col <- unlist(colors_cols_Top_pie_clonotype())

      plot_pie <- ggplot(data = top_BD_cluster, aes(x = UMAP_1, UMAP_2, colour = Selected_function)) +
        geom_point() +
        scale_color_manual(values = df.col, na.value = input$NA_col_analysis) +
        # scale_size_manual(values = rep(input$size_selected_top,dim(top_BD_cluster)[1]),na.value = 1) +
        theme_bw() +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        ) +
        scale_size(guide = "none")

      # if (input$Split_by_group2=="no") {
      #   plot_pie
      # }
      # else {

      plot_pie <- plot_pie + facet_wrap(~ top_BD_cluster$Selected_group, nrow = input$wrap_row)
      #

      plot_pie
    })


    #### top heatmap ------
    Selected_top_plot <- reactive({
      cluster <- top_clonotype_bar_code()
      validate(
        need(
          nrow(cluster) > 0,
          "Upload Seurat file"
        )
      )
      cluster$ID_Column2 <- cluster[, names(cluster) %in% input$Split_group_by_]
      Network_df <- cluster
      Network_df$Selected <- Network_df[, names(Network_df) %in% input$Colour_By_this]
      Network_df$cloneCount <- 1
      df.1 <- acast(Network_df, ID_Column2 ~ Selected, value.var = "cloneCount")
      df.1[is.na(df.1)] <- 0
      min.FC <- 0
      med.FC <- 1
      max.FC <- max(df.1)

      gp_labs <- input$text_size - 2

      ht <- Heatmap(df.1,
                    heatmap_legend_param = list(
                      title = "Count",
                      title_gp = gpar(
                        fontsize = input$Legend_size,
                        fontface = "bold", fontfamily = input$font_type
                      ),
                      labels_gp = gpar(fontsize = input$Legend_size, fontfamily = input$font_type)
                    ),
                    col = colorRamp2(c(min.FC,med.FC, max.FC), c(input$min_FC_col, input$med_FC_col, input$max_FC_col)),
                    row_names_gp = grid::gpar(fontsize = input$text_size, fontfamily = input$font_type),
                    column_names_gp = grid::gpar(fontsize = input$text_size, fontfamily = input$font_type),
      )

      draw(ht, padding = unit(c(5, 5, 5, 5), "mm"))
    })



    output$heatmap_topclone_plot <- renderPlot({
      # Pie_chart_ClusTCR2()
      Selected_top_plot()
    })


    output$downloadPlot_heatmap_topclone_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Selected_clonotype, "_heatmap_topclone_plot_", input$Split_group_by_, "_", input$Colour_By_this, "_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_heatmap_topclone_plot, height = input$height_heatmap_topclone_plot, onefile = FALSE) # open the pdf device
        plot(Selected_top_plot())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_heatmap_topclone_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Selected_clonotype, "_heatmap_topclone_plot_", input$Split_group_by_, "_", input$Colour_By_this, "_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_heatmap_topclone_plot,
            height = input$height_png_heatmap_topclone_plot,
            res = input$resolution_heatmap_topclone_plot
        )
        plot(Selected_top_plot())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    # plot pie and UMAP ----

    output$downloadPlot_top_clonotype_pie <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Selected_clonotype, "_top_clonotype_pie_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_top_clonotype_pie, height = input$height_top_clonotype_pie, onefile = FALSE) # open the pdf device
        if (input$Plot_type_selected == "pie") {
          plot(Pie_chart_alpha_gamma())
        } else {
          plot(UMAP_chart_alpha_gamma())
        }
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_top_clonotype_pie <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Selected_clonotype, "_top_clonotype_pie_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_top_clonotype_pie,
            height = input$height_png_top_clonotype_pie,
            res = input$resolution_PNG_top_clonotype_pie
        )

        if (input$Plot_type_selected == "pie") {
          plot(Pie_chart_alpha_gamma())
        } else {
          plot(UMAP_chart_alpha_gamma())
        }


        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )



    ### Add in the TCR clustering to the seurat object ----

    #### Ridge plot (expression per T cell) ----
    vals_Ridge_top <- reactiveValues(output_stats = NULL)

    compare.stat <- reactive({
      sc <- UMAP_metadata_with_labs()

      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )

      df <- sc@meta.data
      df
      unique.df <- (df[, names(df) %in% c(input$Samp_col, input$V_gene_sc)])
      names(unique.df) <- c("group", "chain")

      unique.df <- subset(unique.df, unique.df$chain != "NA")
      unique.df <- subset(unique.df, unique.df$group != "NA")
      unique.df$cloneCount <- 1
      unique.df
      df_unique_sum <- ddply(unique.df, names(unique.df)[-c(3)], numcolwise(sum))
      df_unique_sum <- df_unique_sum[order(df_unique_sum$cloneCount, decreasing = T), ]

      sc@meta.data$Vgene <- sc@meta.data[, names(sc@meta.data) %in% input$V_gene_sc]

      name.clone <- input$Selected_clonotype
      sc@meta.data$Gene_select <- ifelse(sc@meta.data$Vgene %in% name.clone, name.clone, "unselected")
      sc@meta.data
      unique(sc@meta.data$Gene_select)
      Idents(object = sc) <- sc@meta.data$Gene_select
      as.data.frame(Idents(object = sc))
      min.pct.expression <- input$min_point_ # standard setting: 0.25§
      min.logfc <- input$LogFC_ # 0.25 is standard
      p.val.cutoff <- input$pval_top # (1/10^3) is standard, use (1/10^0) to ignore

      cluster.names <- unique(Idents(sc))[order(unique(Idents(sc)))]
      as.data.frame(cluster.names)
      # # print(paste0("calculating markers for cluster ",name.clone,". Total: ",length(cluster.names)," clusters"))
      markers.fm.list <- FindMarkers(sc, ident.1 = name.clone, min.pct = min.pct.expression, logfc.threshold = min.logfc, only.pos = TRUE)
      markers.fm.list
      # markers.fm.list2 <- subset(markers.fm.list,markers.fm.list$p_val_adj < input$pval.ex.filter)
      # as.data.frame(markers.fm.list2)
    })
    output$Ridge_chart_alpha_gamma_stat_comp <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      sc <- data_sc_pro()

      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )
      compare.stat()
    })

    output$downloaddf_FindMarker_Top <- downloadHandler(
      filename = function() {
        paste("Stats_", gsub("-", ".", Sys.Date()), ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(compare.stat())
        write.csv(df, file, row.names = T)
      }
    )


    observeEvent(input$run_string.data_Exp_top, {
      # df <- compare.stat()

      sc <- data_sc_pro()

      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )

      # if(grepl(5, sc@version)) {
      #
      #
      # }

      df <- as.data.frame(rownames(sc@assays$RNA$scale.data))
      names(df) <- "V1"

      req(df)
      updateSelectInput(
        session,
        "string.data_Exp_top",
        choices = df$V1,
        selected = df$V1[1]
      )
    })


    output$test.table_ridge <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      sc <- data_sc_pro()

      validate(
        need(
          nrow(sc) > 0,
          "Upload file for annotation"
        )
      )

      df <- as.data.frame(sc@assays$RNA$scale.data)
      req(df)
      as.data.frame(df)[1:6]
    })

    # ridge plot ------
    Ridge_chart_alpha_gamma_df <- reactive({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )

      df <- as.data.frame(sc@assays$RNA$scale.data)
      req(df)

      req(input$string.data_Exp_top)

      MainTcell <- as.data.frame(t(df))
      meta.data <- as.data.frame(sc@meta.data)
      # names(MainTcell) <- ifelse(grepl("pAbO",names(MainTcell)),names(MainTcell),toupper(names(MainTcell)))
      MainTcell$Cell_Index <- rownames(MainTcell)
      # names(MainTcell) <- gsub("CD8B.","CD8B",names(MainTcell))

      gene_df <- MainTcell[, names(MainTcell) %in% c("Cell_Index", input$string.data_Exp_top)]

      top_BD_cluster <- top_clonotype_bar_code()
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Split_group_by_]
      top_BD_cluster <- top_BD_cluster[top_BD_cluster$Selected_function %in% input$Graph_split_order, ]
      top_BD_cluster$Selected_function <- as.character(top_BD_cluster[, names(top_BD_cluster) %in% input$Split_group_by_])
      top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function, levels = input$Graph_split_order)

      merge(top_BD_cluster, gene_df, by = "Cell_Index")
    })

    Ridge_chart_alpha_gamma_plot <- reactive({
      df <- Ridge_chart_alpha_gamma_df()
      df$expressed <- ifelse(df[, names(df) %in% input$string.data_Exp_top] > input$Gre_ex, "expressed", "Not expressed")

      if (input$restric_ex == T) {
        df2 <- subset(df, df$expressed == "expressed")
      } else {
        (
          df2 <- df
        )
      }

      df2$Selected_function <- factor(df2$Selected_function, levels = input$Graph_split_order)

      ggplot(df2, aes(x = get(input$string.data_Exp_top), y = Selected_function, fill = Selected_function)) +
        geom_density_ridges() +
        theme_ridges() +
        theme(
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          # legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none"
        ) +
        ggtitle(input$string.data_Exp_top)
    })

    Ridge_chart_alpha_gamma_stat_table <- reactive({
      df <- Ridge_chart_alpha_gamma_df()
      df$expressed <- ifelse(df[, names(df) %in% input$string.data_Exp_top] > input$Gre_ex, "expressed", "Not expressed")

      if (input$restric_ex == T) {
        df2 <- subset(df, df$expressed == "expressed")
      } else {
        (
          df2 <- df
        )
      }
      df2[is.na(df2)] <- "-"
      at <- TukeyHSD(aov(get(input$string.data_Exp_top) ~ get(input$Split_group_by_), data = df2))
      tab <- as.data.frame(at[1])
      names(tab) <- c("diff", "lwr", "upr", "p.adj")
      tab$stat <- ifelse(tab$p.adj < 0.0001, "****",
                         ifelse(tab$p.adj < 0.001, "***",
                                ifelse(tab$p.adj < 0.01, "**",
                                       ifelse(tab$p.adj < 0.05, "*", "NS")
                                )
                         )
      )
      tab
    })

    output$Ridge_chart_alpha_gamma_stat <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      Ridge_chart_alpha_gamma_stat_table()
    })

    output$downloaddf_clusTCR_GEx <- downloadHandler(
      filename = function() {
        paste("Stats_", gsub("-", ".", Sys.Date()), ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(Ridge_chart_alpha_gamma_stat_table())
        write.csv(df, file, row.names = T)
      }
    )
    ### Violin plots top clonotype -----
    Violin_plot_table <- reactive({
      sc <- data_sc_pro()
      need(
        nrow(sc) > 0,
        error_message_val_sc
      )

      df <- as.data.frame(sc@assays$RNA$scale.data)
      req(df)

      MainTcell <- as.data.frame(t(df))
      meta.data <- as.data.frame(sc@meta.data)
      # names(MainTcell) <- ifelse(grepl("pAbO",names(MainTcell)),names(MainTcell),toupper(names(MainTcell)))
      MainTcell$Cell_Index <- rownames(MainTcell)
      # names(MainTcell) <- gsub("CD8B.","CD8B",names(MainTcell))

      gene_df <- MainTcell[, names(MainTcell) %in% c("Cell_Index", input$string.data_Exp_top)]

      top_BD_cluster <- top_clonotype_bar_code()
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Split_group_by_]
      top_BD_cluster <- top_BD_cluster[top_BD_cluster$Selected_function %in% input$Graph_split_order, ]
      top_BD_cluster$Selected_function <- as.character(top_BD_cluster[, names(top_BD_cluster) %in% input$Split_group_by_])
      top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function, levels = input$Graph_split_order)

      merge(top_BD_cluster, gene_df, by = "Cell_Index")
    })


    Violin_chart_alpha_gamma_plot <- reactive({
      df <- Ridge_chart_alpha_gamma_df()
      # df <- df[df$Selected_group %in% input$Graph_split_order,]
      # top_BD_cluster$Selected_group <- factor(top_BD_cluster$Selected_group,levels = input$Graph_split_order)
      df$expressed <- ifelse(df[, names(df) %in% input$string.data_Exp_top] > input$Gre_ex, "expressed", "Not expressed")

      if (input$restric_ex == T) {
        df2 <- subset(df, df$expressed == "expressed")
      } else {
        (
          df2 <- df
        )
      }
      df2$Selected_function <- factor(df2$Selected_function, levels = input$Graph_split_order)
      ggplot(df2, aes(y = get(input$string.data_Exp_top), x = Selected_function, fill = Selected_function)) +
        geom_violin() +
        geom_jitter(height = 0, width = 0.1) +
        theme(legend.position = "none", ) +
        theme_bw() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
          axis.title.x = element_blank(),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_text(colour = "black", size = 20, family = input$font_type),
          legend.position = "none",
        ) +
        ggtitle(input$string.data_Exp_top)
    })

    # vs non-selected clonotypes -----

    Ridge_chart_alpha_gamma_df_all <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )
      df <- as.data.frame(sc@assays$RNA$scale.data)
      req(df)

      MainTcell <- as.data.frame(t(df))
      meta.data <- as.data.frame(sc@meta.data)
      MainTcell$Cell_Index <- rownames(MainTcell)
      gene_df <- MainTcell[, names(MainTcell) %in% c("Cell_Index", input$string.data_Exp_top)]
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[, names(df3.meta) %in% input$V_gene_sc]
      df3.meta$selected_top_clonotype <- ifelse(df3.meta$cluster_name == input$Selected_clonotype, input$Selected_clonotype, "other")
      top_BD_clonotype2 <- merge(df3.meta, gene_df, by = "Cell_Index")
      top_BD_clonotype2
    })
    Ridge_chart_alpha_gamma_plot_comp <- reactive({
      df <- Ridge_chart_alpha_gamma_df_all()
      df$expressed <- ifelse(df[, names(df) %in% input$string.data_Exp_top] > input$Gre_ex, "expressed", "Not expressed")

      if (input$restric_ex == T) {
        df2 <- subset(df, df$expressed == "expressed")
      } else {
        (
          df2 <- df
        )
      }

      ggplot(df2, aes(x = get(input$string.data_Exp_top), y = selected_top_clonotype, fill = selected_top_clonotype)) +
        geom_density_ridges() +
        theme_ridges() +
        theme(
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          # axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none"
        )
    })

    Violin_chart_alpha_gamma_plot_comp <- reactive({
      df <- Ridge_chart_alpha_gamma_df_all()
      df$expressed <- ifelse(df[, names(df) %in% input$string.data_Exp_top] > input$Gre_ex, "expressed", "Not expressed")

      if (input$restric_ex == T) {
        df2 <- subset(df, df$expressed == "expressed")
      } else {
        (
          df2 <- df
        )
      }

      ggplot(df2, aes(y = get(input$string.data_Exp_top), x = selected_top_clonotype, fill = selected_top_clonotype)) +
        geom_violin() +
        geom_jitter(height = 0, width = 0.1) +
        theme(legend.position = "none", ) +
        theme_bw() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )
    })
    output$Ridge_chart_alpha_gamma_plot_out <- renderPlot({
      if (input$plot_type_ridgvi == "Ridge (selected clonotype)") {
        Ridge_chart_alpha_gamma_plot()
      } else if (input$plot_type_ridgvi == "Ridge (compare)") {
        Ridge_chart_alpha_gamma_plot_comp()
      } else if (input$plot_type_ridgvi == "Violin (compare)") {
        Violin_chart_alpha_gamma_plot_comp()
      } else {
        Violin_chart_alpha_gamma_plot()
      }
    })

    output$downloadPlot_Ridge_chart_alpha_gamma_plot_out <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$string.data_Exp_top, "_", input$Selected_clonotype, "_", input$plot_type_ridgvi, "_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_Ridge_chart_alpha_gamma_plot_out, height = input$height_Ridge_chart_alpha_gamma_plot_out, onefile = FALSE) # open the pdf device
        if (input$plot_type_ridgvi == "Ridge (selected clonotype)") {
          df <- Ridge_chart_alpha_gamma_plot()
        } else if (input$plot_type_ridgvi == "Ridge (compare)") {
          df <- Ridge_chart_alpha_gamma_plot_comp()
        } else if (input$plot_type_ridgvi == "Violin (compare)") {
          df <- Violin_chart_alpha_gamma_plot_comp()
        } else {
          df <- Violin_chart_alpha_gamma_plot()
        }
        plot(df)
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_Ridge_chart_alpha_gamma_plot_out <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$string.data_Exp_top, "_", input$Selected_clonotype, "_", input$plot_type_ridgvi, "_", x, ".png", sep = "")
        # paste("_",input$epitope_umap_selected,"_",input$epitope_umap_selected2,"_Heatmap_", x, "", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_Ridge_chart_alpha_gamma_plot_out,
            height = input$height_png_Ridge_chart_alpha_gamma_plot_out,
            res = input$resolution_PNG_Ridge_chart_alpha_gamma_plot_out
        )
        if (input$plot_type_ridgvi == "Ridge (selected clonotype)") {
          df <- Ridge_chart_alpha_gamma_plot()
        } else if (input$plot_type_ridgvi == "Ridge (compare)") {
          df <- Ridge_chart_alpha_gamma_plot_comp()
        } else if (input$plot_type_ridgvi == "Violin (compare)") {
          df <- Violin_chart_alpha_gamma_plot_comp()
        } else {
          df <- Violin_chart_alpha_gamma_plot()
        }
        plot(df)
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    ### top dot plot -----
    all_expression_plot_top <- reactive({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )

      df <- sc@meta.data
      unique.df <- (df[, names(df) %in% c(input$Samp_col, input$V_gene_sc)])
      names(unique.df) <- c("group", "chain")

      # unique.df <- unique.df[unique.df$group %in% "LTR35-12",]

      unique.df <- subset(unique.df, unique.df$chain != "NA")
      unique.df <- subset(unique.df, unique.df$group != "NA")
      unique.df$cloneCount <- 1

      df_unique_sum <- ddply(unique.df, names(unique.df)[-c(3)], numcolwise(sum))
      df_unique_sum <- df_unique_sum[order(df_unique_sum$cloneCount, decreasing = T), ]

      sc@meta.data$Vgene <- sc@meta.data[, names(sc@meta.data) %in% input$V_gene_sc]

      name.clone <- input$Selected_clonotype

      sc@meta.data$Gene_select <- ifelse(sc@meta.data$Vgene %in% name.clone, input$name_clonotype_selected, "unselected")
      unique(sc@meta.data$Gene_select)
      Idents(object = sc) <- sc@meta.data$Gene_select
      if (input$restict_no_points == F) {
        list.names <- rownames(compare.stat())
      } else {
        list.names <- rownames(compare.stat())
        list.names <- list.names[1:input$pval.ex.top_genes]
      }

      size_legend <- input$Legend_size - 2

      DotPlot(sc, features = list.names) +
        RotatedAxis() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
          axis.title.x = element_blank(),
          legend.title = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.text = element_text(colour = "black", size = size_legend, family = input$font_type),
          legend.position = input$legend_position,
        ) +
        scale_colour_gradient2(low = input$low.dotplot, mid = input$middle.dotplot, high = input$high.dotplot) +
        scale_x_discrete(labels = label_wrap(20)) +
        scale_y_discrete(labels = label_wrap(20))
    })

    output$all_expression_dotplot_top <- renderPlot({
      all_expression_plot_top()
    })

    output$downloadPlot_all_expression_dotplot_top <- downloadHandler(
      filename = function() {
        paste(input$Selected_clonotype, "_dotplot", "_", today(), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_all_expression_dotplot_top, height = input$height_all_expression_dotplot_top, onefile = FALSE) # open the pdf device
        df <- all_expression_plot_top()
        plot(df)
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_all_expression_dotplot_top <- downloadHandler(
      filename = function() {
        paste(input$Selected_clonotype, "_dotplot", "_", today(), ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_all_expression_dotplot_top,
            height = input$height_png_all_expression_dotplot_top,
            res = input$resolution_PNG_all_expression_dotplot_top
        )
        df <- all_expression_plot_top()

        plot(df)
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    # Over representation analysis for top clonotypes -----
    observe({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          "Run Variable"
        )
      )
      gene.names <- rownames(sc@assays$RNA$counts)

      gene.names_hs <- gene.names[str_detect(gene.names, "^[A-Z][A-Z/0-9][A-Z/0-9]") &
                                    !str_detect(gene.names, "^X[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^A[0-9][0-9][0-9]") & !str_detect(gene.names, "^AC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^B[0-9][0-9][0-9]") & !str_detect(gene.names, "^BC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^C[0-9][0-9][0-9]") & !str_detect(gene.names, "^BC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^D[0-9][0-9][0-9]") & !str_detect(gene.names, "^E[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^F[0-9][0-9][0-9]") & !str_detect(gene.names, "^G[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^H[0-9][0-9][0-9]") & !str_detect(gene.names, "^I[0-9][0-9][0-9]")]

      gene.names_mm <- gene.names[str_detect(gene.names, "^[A-Z][a-z/0-9][a-z/0-9]") &
                                    !str_detect(gene.names, "^X[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^A[0-9][0-9][0-9]") & !str_detect(gene.names, "^AC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^B[0-9][0-9][0-9]") & !str_detect(gene.names, "^BC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^C[0-9][0-9][0-9]") & !str_detect(gene.names, "^BC[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^D[0-9][0-9][0-9]") & !str_detect(gene.names, "^E[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^F[0-9][0-9][0-9]") & !str_detect(gene.names, "^G[0-9][0-9][0-9]") &
                                    !str_detect(gene.names, "^H[0-9][0-9][0-9]") & !str_detect(gene.names, "^I[0-9][0-9][0-9]")]

      if (length(gene.names_hs) > 0) {
        len_hs <- length(gene.names_hs)
      } else {
        len_hs <- 1
      }

      if (length(gene.names_mm) > 0) {
        len_mm <- length(gene.names_mm)
      } else {
        len_mm <- 1
      }

      ratio_hs.mm <- len_hs / len_mm

      if (ratio_hs.mm > 1) {
        updateSelectInput(
          session,
          "species_analysis",
          choices = c("hs", "mm"),
          selected = "hs"
        )
      } else if (ratio_hs.mm < 1) {
        updateSelectInput(
          session,
          "species_analysis",
          choices = c("hs", "mm"),
          selected = "mm"
        )
      } else {
        updateSelectInput(
          session,
          "species_analysis",
          choices = c("hs", "mm"),
          selected = "hs"
        )
      }
    })



    Over_rep_Top_clones_old <- reactive({
      sc <- data_sc_pro()

      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )
      req(input$species_analysis, input$datasource)

      df <- sc@meta.data
      # require()

      geneSet <- read.csv(system.file("OverRep", "GeneSets.csv", package = "STEGO.R"), header = T)

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

      if (input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        geneSet$GeneSet <- gsub("-", ".", geneSet$GeneSet)
      }

      if (input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        require(stringr)
        geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
      }

      for (i in 1:dim(geneSet)[1]) {
        # listed GeneSet
        message(paste("GeneSet: ", i))
        Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet, ";")[i])
        names(Gene.set.testing) <- "V1"
        Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
        names(Gene.set.testing2) <- "V1"
        # message(paste(dim(Gene.set.testing2)[1],"GeneSet total"))
        # background genes
        background.overlap <- merge(Gene.set.testing2, background.genes.name, by = "V1")
        # message(paste(dim(background.overlap)[1],"in Background"))
        geneSet$background.geneset[i] <- length(background.overlap$V1)
        geneSet$background.geneset.name[i] <- as.character(paste(unlist(background.overlap[1]), collapse = ";"))
        # in sig gene list
        overlap <- merge(background.overlap, DEx.genes, by = "V1")
        # message(paste(dim(overlap)[1],"# Sig genes"))
        geneSet$in.geneset[i] <- length(overlap$V1)
        geneSet$in.geneset.name[i] <- as.character(paste(unlist(overlap[1]), collapse = ";"))
      }

      geneSet2 <- subset(geneSet, geneSet$in.geneset > 0)

      for (i in 1:dim(geneSet2)[1]) {
        tota.gene.set <- geneSet2$background.geneset[i] # genes that are identified in background
        tota.gene.set
        in.geneset <- geneSet2$in.geneset[i] # DEx in geneset

        background.genes
        not.in.total <- background.genes - tota.gene.set
        not.in.geneset.sig <- total.sig - in.geneset
        d <- data.frame(gene.in.interest = c(in.geneset, not.in.geneset.sig), gene.not.interest = c(tota.gene.set, not.in.total))
        row.names(d) <- c("In_category", "not_in_category")

        if (in.geneset > 0) {
          geneSet2$p.val[i] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
          geneSet2$lowerCI[i] <- unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
          geneSet2$upperCI[i] <- unlist(fisher.test(d)$conf.int)[2]
          geneSet2$OR[i] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1], 3)
        } else {
          geneSet2$p.value[i] <- "-"
          geneSet2$lowerCI[i] <- "-"
          geneSet2$upperCI[i] <- "-"
          geneSet2$OR[i] <- "-"
        }
        # message(print(d))
        # message(print(round(prop.table(d),3)))
      }

      geneSet2
    })


    Over_rep_Top_clones <- reactive({
      geneSet2 <- Over_rep_Top_clones_old()

      validate(
        need(
          nrow(geneSet2) > 0,
          error_message_val_sc
        )
      )
      geneSet2 <- geneSet2[order(geneSet2$p.val, decreasing = F), ]
      geneSet2 <- subset(geneSet2, geneSet2$in.geneset >= input$in.geneset.cutoff_top)
      geneSet2 <- subset(geneSet2, geneSet2$p.val <= input$p.val_cutoff_top)
      geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
      geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
      # name.list <- c("Geneset_ID","p.val","FDR","Bonferroni","OR","lowerCI","upperCI","in.geneset.name","in.geneset","background.geneset","total.sig","background.genes","background.geneset.name")
      # geneSet2 <- geneSet2 %>%
      #     select(all_of(name.list), everything())

      geneSet2
    })

    output$Over_rep_Top_clones_Tab <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1, 2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      Over_rep_Top_clones()
    })

    output$downloadtb_over.rep.Top_Ex <- downloadHandler(
      filename = function() {
        paste0(input$Selected_clonotype, "_", today(), "_over_rep.csv")
      },
      content = function(file) {
        df <- as.data.frame(Over_rep_Top_clones())
        write.csv(df, file, row.names = F)
      }
    )

    # Expanded TCR interrogation (regarless of TCR sequence ) -------

    Expansion_check_tb <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      req(input$Split_group_by_, input$Graph_split_order, input$V_gene_sc, input$cut.off_expanded)
      message(paste("Calculating expansion"))
      if (input$Split_group_by_ != "") {
        Expansion.meta.data <- sc@meta.data
        Expansion.meta.data$Samp <- Expansion.meta.data[, names(Expansion.meta.data) %in% input$Split_group_by_] # changing to the Selected Individual column

        if (length(input$Graph_split_order) > 0) {
          Expansion.meta.data <- Expansion.meta.data[Expansion.meta.data$Samp %in% c(input$Graph_split_order), ]

          if (dim(Expansion.meta.data)[1] > 0) {
            CDR3_Vgene_all <- as.data.frame(Expansion.meta.data[, names(Expansion.meta.data) %in% input$V_gene_sc])
            names(CDR3_Vgene_all) <- "V1"
            CDR3_Vgene_all[is.na(CDR3_Vgene_all)] <- "Unknown"
            CDR3_Vgene_all$count <- 1
            req(CDR3_Vgene_all)
            total.condition <- as.data.frame(ddply(CDR3_Vgene_all, "V1", numcolwise(sum)))
            req(total.condition)
            total.condition <- total.condition[total.condition$V1 != "Unknown", ]
            req(total.condition)
            total.condition
            CDR3_Vgene_all2 <- as.data.frame(Expansion.meta.data[, names(Expansion.meta.data) %in% c(input$V_gene_sc, "Samp")])
            names(CDR3_Vgene_all2)[names(CDR3_Vgene_all2) %in% input$V_gene_sc] <- "V1"
            CDR3_Vgene_all2[is.na(CDR3_Vgene_all2)] <- "Unknown"
            names(CDR3_Vgene_all2)[names(CDR3_Vgene_all2) %in% "Samp"] <- "V2"
            CDR3_Vgene_all2$count.group <- 1
            total.condition.group <- as.data.frame(ddply(CDR3_Vgene_all2, c("V1", "V2"), numcolwise(sum)))
            total.condition.group <- total.condition.group[total.condition.group$V1 != "Unknown", ]
            total.condition.group <- total.condition.group[total.condition.group$V2 != "Unknown", ]
            total.condition.group
            total.condition.group.counts <- merge(total.condition.group, total.condition, by = "V1")
            total.condition.group.counts <- total.condition.group.counts[order(total.condition.group.counts$count.group, decreasing = T), ]
            total.condition.group.counts$expand.singlets <- ifelse(total.condition.group.counts$count.group > input$cut.off_expanded, "Ex", "NEx")
            names(total.condition.group.counts)[names(total.condition.group.counts) %in% "V1"] <- input$V_gene_sc
            names(total.condition.group.counts)[names(total.condition.group.counts) %in% "V2"] <- "Samp"
            total.condition.group.counts
            umap.meta2 <- merge(Expansion.meta.data, total.condition.group.counts, by = c("Samp", input$V_gene_sc), sort = F, all.x = T)
            umap.meta2
            rownames(umap.meta2) <- umap.meta2$Cell_Index
            umap.meta2
            message(paste("Finished calculating expansion"))
            umap.meta2$Selected_Status <- umap.meta2$Samp
            umap.meta2$expansion.status <- paste0(umap.meta2$Selected_Status, "_", umap.meta2$expand.singlets)
            umap.meta2 <- umap.meta2[!is.na(umap.meta2$Samp), ]
            umap.meta2 <- umap.meta2[!is.na(umap.meta2$expand.singlets), ]
            umap.meta2
          }
        }
      }
    })

    output$Expansion_check <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )

      Expansion_check_tb()
    })

    select_group_metadata_ex <- reactive({
      df <- Expansion_check_tb()
      validate(
        need(
          nrow(df) > 0,
          error_message_val1
        )
      )
      req(df)
      df2 <- df$expansion.status
      df2 <- as.data.frame(df2)
      df2
    })
    observe({
      df2 <- select_group_metadata_ex()

      validate(
        need(
          nrow(df2) > 0,
          error_message_val1
        )
      )
      req(df2)

      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- as.data.frame(df2[order(df2$V1), ])
      names(df2) <- "V1"
      df2
      df2 <- subset(df2, df2$V1 != "NA")
      df3 <- unique(df2$V1)
      updateSelectizeInput(
        session,
        "selected_Indiv_Ex_1",
        choices = df3,
        selected = df3[1]
      )
    })
    observe({
      df2 <- select_group_metadata_ex()
      validate(
        need(
          nrow(df2) > 0,
          error_message_val1
        )
      )
      req(df2)
      df2 <- as.data.frame(df2)
      names(df2) <- "V1"
      df2 <- as.data.frame(df2[order(df2$V1), ])
      names(df2) <- "V1"
      df2
      df2 <- subset(df2, df2$V1 != "NA")
      df3 <- unique(df2$V1)
      updateSelectizeInput(
        session,
        "selected_Indiv_Ex_2",
        choices = df3,
        selected = df3[2]
      )
    })
    # #
    # # ### table check ------
    Vals_Expans_table <- reactiveValues(output_expTab = NULL)
    # # ### add in colouring specific to Expanded
    # #
    output$classification_to_add2 <- renderUI({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      req(Expansion_check_tb())
      df3.meta <- c(names(Expansion_check_tb()))

      if (length(df3.meta) > 0) {
        selectInput("Colour_By_this_Expanded", "Expanded Colour by: ", choices = df3.meta, selected = "expand.singlets")
      } else {
        selectInput("Colour_By_this_Expanded", "Expanded Colour by: ", choices = "expand.singlets", selected = "expand.singlets")
      }
    })
    #
    #
    # ## Expansion UMAP plot -----
    cols_UMAP_Expanded <- reactive({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      req(input$Colour_By_this_Expanded)

      top_BD_cluster <- Expansion_check_tb()
      req(top_BD_cluster)
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this_Expanded]

      if (input$Colour_By_this_Expanded == "expand.singlets") {
        top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function), NA, top_BD_cluster$Selected_function)
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function, levels = c("NEx", "Ex"))
        top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function, decreasing = F), ]
      } else {
        top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function), NA, top_BD_cluster$Selected_function)
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function, levels = unique(top_BD_cluster$Selected_function))
      }
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num) == T, ])
      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$Colour_By_this_Expanded == "expand.singlets") {
        col.gg[1] <- "grey90"
        col.gg[2] <- "#00BFC4"

        palette_rainbow[1] <- "grey90"
        heat_col[1] <- "grey90"
        col.terrain[1] <- "grey90"
        col.topo[1] <- "grey90"
        col.hcl[1] <- "grey90"
      } else {

      }


      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep = "_"), paste(num[i, ]), col.gg[i])
        })
      } else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep = "_"), paste(num[i, ]), col.hcl[i])
        })
      } else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep = "_"), paste(num[i, ]), col.topo[i])
        })
      } else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep = "_"), paste(num[i, ]), heat_col[i])
        })
      } else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep = "_"), paste(num[i, ]), col.terrain[i])
        })
      } else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep = "_"), paste(num[i, ]), palette_rainbow[i])
        })
      } else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep = "_"), paste(num[i, ]), palette1[i])
        })
      } else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.UMAP_Expanded", i, sep = "_"), paste(num[i, ]), input$one.colour.default)
        })
      } # one colour
    })
    output$cols_UMAP_Expanded <- renderUI({
      cols_UMAP_Expanded()
    })
    colors_Expanded <- reactive({
      sc <- data_sc_pro()

      validate(
        need(
          nrow(sc) > 0,
          error_message_val1
        )
      )
      top_BD_cluster <- Expansion_check_tb()
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this_Expanded]



      if (input$Colour_By_this_Expanded == "expand.singlets") {
        top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function), NA, top_BD_cluster$Selected_function)
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function, levels = c("NEx", "Ex"))
        top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function, decreasing = F), ]
      } else {
        top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function), NA, top_BD_cluster$Selected_function)
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function, levels = unique(top_BD_cluster$Selected_function))
      }

      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.UMAP_Expanded", i, sep = "_")]]
      })
    })
    UMAP_Expanded_plot <- reactive({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val1
        )
      )
      req(input$Colour_By_this_Expanded, input$wrap_row)
      top_BD_cluster <- Expansion_check_tb()
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this_Expanded]

      if (input$Colour_By_this_Expanded == "expand.singlets") {
        top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function), NA, top_BD_cluster$Selected_function)
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function, levels = c("NEx", "Ex"))
        top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function, decreasing = F), ]
      } else {
        top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function), NA, top_BD_cluster$Selected_function)
        top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function, levels = unique(top_BD_cluster$Selected_function))
      }

      col.file <- as.data.frame(unique(top_BD_cluster$Selected_function))
      col.file <- as.data.frame(col.file[complete.cases(col.file) == T, ])
      names(col.file) <- "V1"

      col.file$col <- unlist(colors_Expanded())
      colorblind_vector <- as.data.frame(unlist(colors_Expanded()))

      if (dim(colorblind_vector)[1] == 0) {
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
        } else {

        }
      }
      colorblind_vector <- as.data.frame(colorblind_vector)

      names(colorblind_vector) <- "cols"
      top_BD_cluster$ID_Column <- top_BD_cluster[, names(top_BD_cluster) %in% input$Split_group_by_]
      top_BD_cluster <- top_BD_cluster[top_BD_cluster$ID_Column %in% input$Graph_split_order, ] # subset to a specific individuals
      top_BD_cluster$ID_Column <- factor(top_BD_cluster$ID_Column, levels = input$Graph_split_order)

      df <- ggplot(top_BD_cluster, aes(x = UMAP_1, UMAP_2, colour = Selected_function, alpha = Selected_function, label = Selected_function))
      df <- df + geom_point()
      df <- df + scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 20), values = colorblind_vector$col, na.value = input$NA_col_analysis) +
        scale_alpha_manual(labels = ~ stringr::str_wrap(.x, width = 20), values = rep(1, length(unique(top_BD_cluster$Selected_function))), na.value = 0.1) +
        theme_bw() +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 0),
          axis.title.x = element_text(colour = "black", angle = 0, vjust = .5, face = "plain", family = input$font_type, size = input$title.text.sizer2),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )

      if (input$Split_by_group == "no") {
        df <- df
      } else {
        df <- df + facet_wrap(~ID_Column, nrow = input$wrap_row)
      }
      df
    })
    output$UMAP_Expanded <- renderPlot({
      UMAP_Expanded_plot()
    })
    #
    # # download UMAP expanded ------
    output$downloadPlot_UMAP_Expanded <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_Expanded_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_UMAP_Expanded, height = input$height_UMAP_Expanded, onefile = FALSE) # open the pdf device
        plot(UMAP_Expanded_plot())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_UMAP_Expanded <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_Expanded_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_UMAP_Expanded,
            height = input$height_png_UMAP_Expanded,
            res = input$resolution_PNG_UMAP_Expanded
        )
        plot(UMAP_Expanded_plot())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    # ## Expanded Violin plot ------
    #
    compare.stat_Expanded <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      req(input$Samp_col, Expansion_check_tb())
      md <- sc@meta.data
      ex.md <- Expansion_check_tb()
      ex.md2 <- ex.md[, names(ex.md) %in% c("Cell_Index", "expanded.singlets", "Selected_Status", "expansion.status")]
      md.ex <- merge(md, ex.md2, sort = F)
      rownames(md.ex) <- md.ex$Cell_Index
      sc@meta.data <- md.ex
      sc
    })
    # Vals_expanded.stats <- reactiveValues(output_ex1=NULL)

    Vals_expanded.stats <- reactive({
      sc <- data_sc_pro()

      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )

      sc <- compare.stat_Expanded()
      req(input$selected_Indiv_Ex_1, input$selected_Indiv_Ex_2, input$pval.ex.filter)
      # head(sc@meta.data)
      message("Updating Ident")
      Idents(object = sc) <- sc@meta.data$expansion.status

      min.pct.expression <- input$min_point_ # standard setting: 0.25
      min.logfc <- input$LogFC_ # 0.25 is standard
      message(paste0(" Calculating markers for cluster ", c(input$selected_Indiv_Ex_1), " vs ", c(input$selected_Indiv_Ex_2)))
      markers.fm.list <- FindMarkers(sc, ident.1 = input$selected_Indiv_Ex_1, ident.2 = c(input$selected_Indiv_Ex_2), min.pct = min.pct.expression, logfc.threshold = min.logfc, only.pos = TRUE)
      markers.fm.list2 <- subset(markers.fm.list, markers.fm.list$p_val_adj < input$pval.ex.filter)
      Vals_expanded.stats <- as.data.frame(markers.fm.list2)
      Vals_expanded.stats
    })

    output$compare.stat_Ex <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      sc <- data_sc_pro()

      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )
      as.data.frame(Vals_expanded.stats())
    })
    #
    output$downloadtb_compare.stat_Ex <- downloadHandler(
      filename = function() {
        paste("compare.stat_Expanded_", input$Split_group_by_, ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(Vals_expanded.stats())
        write.csv(df, file, row.names = T)
      }
    )
    ## dot plot Expanded ----
    relative_expression_plot_ex <- reactive({
      sc <- Vals_expanded.stats()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )

      if (input$restrict.dotpot == "yes") {
        list.names <- rownames(Vals_expanded.stats())[1:input$restrict.dotpot.num]
      } else {
        list.names <- rownames(Vals_expanded.stats())
      }

      sc <- compare.stat_Expanded()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )

      Idents(object = sc) <- sc@meta.data$expansion.status

      size_legend <- input$Legend_size - 2

      DotPlot(sc, features = list.names) +
        RotatedAxis() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
          axis.title.x = element_blank(),
          legend.title = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.text = element_text(colour = "black", size = size_legend, family = input$font_type),
          legend.position = input$legend_position,
        ) +
        scale_colour_gradient2(low = input$low.dotplot.ex, mid = input$middle.dotplot.ex, high = input$high.dotplot.ex) +
        scale_x_discrete(labels = label_wrap(20)) +
        scale_y_discrete(labels = label_wrap(20))
    })

    output$relative_expression_dotplot_ex <- renderPlot({
      relative_expression_plot_ex()
    })


    output$downloadPlot_all_expression_dotplot_ex <- downloadHandler(
      filename = function() {
        paste("Expanded_dotplot", input$Split_group_by_, "_", today(), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_all_expression_dotplot_ex, height = input$height_all_expression_dotplot_ex, onefile = FALSE) # open the pdf device
        df <- relative_expression_plot_ex()
        plot(df)
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_all_expression_dotplot_ex <- downloadHandler(
      filename = function() {
        paste("Expanded_dotplot", input$Split_group_by_, "_", today(), ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_all_expression_dotplot_ex,
            height = input$height_png_all_expression_dotplot_ex,
            res = input$resolution_PNG_all_expression_dotplot_ex
        )
        df <- relative_expression_plot_ex()

        plot(df)
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    # expression for one gene -------

    Violin_chart_Expanded <- reactive({
      df <- Ridge_chart_alpha_gamma_df()

      df2$Selected_function <- factor(df2$Selected_function, levels = input$Graph_split_order)
      ggplot(df2, aes(y = get(input$string.data_Exp_top), x = Selected_function, fill = Selected_function)) +
        geom_violin() +
        geom_jitter(height = 0, width = 0.1) +
        theme(legend.position = "none", ) +
        theme_bw() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
          axis.title.x = element_blank(),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_text(colour = "black", size = 20, family = input$font_type),
          legend.position = "none",
        ) +
        ggtitle(input$string.data_Exp_top)
    })




    # Over representation analysis for Expanded  -----

    Over_rep_Exp_old <- reactive({
      sc <- data_sc_pro()

      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )

      req(input$datasource, input$in.geneset.cutoff_Exp)

      geneSet <- read.csv(system.file("OverRep", "GeneSets.csv", package = "STEGO.R"), header = T)

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

      if (input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        geneSet$GeneSet <- gsub("-", ".", geneSet$GeneSet)
      }

      if (input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        require(stringr)
        geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
      }

      for (i in 1:dim(geneSet)[1]) {
        # listed GeneSet
        message(paste("GeneSet: ", i))
        Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet, ";")[i])
        names(Gene.set.testing) <- "V1"
        Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
        names(Gene.set.testing2) <- "V1"
        background.overlap <- merge(Gene.set.testing2, background.genes.name, by = "V1")
        # message(paste(dim(background.overlap)[1],"in Background"))
        geneSet$background.geneset[i] <- length(background.overlap$V1)
        geneSet$background.geneset.name[i] <- as.character(paste(unlist(background.overlap[1]), collapse = ";"))
        # in sig gene list
        overlap <- merge(background.overlap, DEx.genes, by = "V1")
        # message(paste(dim(overlap)[1],"# Sig genes"))
        geneSet$in.geneset[i] <- length(overlap$V1)
        geneSet$in.geneset.name[i] <- as.character(paste(unlist(overlap[1]), collapse = ";"))
      }

      geneSet2 <- subset(geneSet, geneSet$in.geneset > 0)

      for (i in 1:dim(geneSet2)[1]) {
        tota.gene.set <- geneSet2$background.geneset[i] # genes that are identified in background
        tota.gene.set
        in.geneset <- geneSet2$in.geneset[i] # DEx in geneset

        background.genes
        not.in.total <- background.genes - tota.gene.set
        not.in.geneset.sig <- total.sig - in.geneset
        d <- data.frame(gene.in.interest = c(in.geneset, not.in.geneset.sig), gene.not.interest = c(tota.gene.set, not.in.total))
        row.names(d) <- c("In_category", "not_in_category")

        if (in.geneset > 0) {
          geneSet2$p.val[i] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
          geneSet2$lowerCI[i] <- unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
          geneSet2$upperCI[i] <- unlist(fisher.test(d)$conf.int)[2]
          geneSet2$OR[i] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1], 3)
        } else {
          geneSet2$p.value[i] <- "-"
          geneSet2$lowerCI[i] <- "-"
          geneSet2$upperCI[i] <- "-"
          geneSet2$OR[i] <- "-"
        }
        # message(print(d))
        # message(print(round(prop.table(d),3)))
      }

      geneSet2
    })

    Over_rep_Exp <- reactive({
      geneSet2 <- Over_rep_Exp_old()

      validate(
        need(
          nrow(geneSet2) > 0,
          error_message_val_sc
        )
      )

      req(input$datasource, input$in.geneset.cutoff_Exp)

      geneSet2 <- geneSet2[order(geneSet2$p.val, decreasing = F), ]
      geneSet2 <- subset(geneSet2, geneSet2$in.geneset >= input$in.geneset.cutoff_Exp)

      # geneSet2 <- subset(geneSet2,geneSet2$pval.BH.adj<=input$adjust_cutoff_Exp)

      geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
      geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
      geneSet2 <- subset(geneSet2, geneSet2$p.val <= input$p.val_cutoff_Exp)
      # name.list <- c("Geneset_ID","p.val","FDR","Bonferroni","OR","lowerCI","upperCI","in.geneset.name","in.geneset","background.geneset","total.sig","background.genes","background.geneset.name")
      # geneSet2 <- geneSet2 %>%
      #   select(all_of(name.list), everything())
      geneSet2
    })
    output$Over_rep_Exp_Tab <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1, 2, 5, 10, 20, 50, 100), pageLength = 20, scrollX = TRUE), {
      Over_rep_Exp()
    })

    output$downloadtb_over.rep_Exp <- downloadHandler(
      filename = function() {
        paste0(input$Selected_clonotype, "_", today(), "_over_rep.csv")
      },
      content = function(file) {
        df <- as.data.frame(Over_rep_Exp())
        write.csv(df, file, row.names = F)
      }
    )
    #### Epitope upload -----
    df_tcrex <- reactive({
      epi <- data_sc_TCRex()
      epi[!(duplicated(epi$CDR3_beta)), ]
    })

    output$MainTcell_Check <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      df_tcrex()
    })
    #### Epitope heatmap -----
    heatmap_epitope <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )

      df3.meta <- sc@meta.data

      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta, epi, by = "CDR3_beta")
      df3.meta$selected <- df3.meta[, names(df3.meta) %in% input$epitope_umap_selected]
      df3.meta <- df3.meta[order(df3.meta$selected, decreasing = F), ]
      df3.meta$selected <- factor(df3.meta$selected, levels = unique(df3.meta$selected))
      df3.meta$Selected_group <- df3.meta[, names(df3.meta) %in% input$epitope_umap_selected2]
      df3.meta$cloneCount <- 1

      # df <- as.data.frame(ddply(df3.meta,(c("selected","Selected_group","cloneCount")),numcolwise(sum)))

      df.1 <- acast(df3.meta, selected ~ Selected_group, value.var = "cloneCount")
      df.1
      head(df.1)
      df.1[is.na(df.1)] <- 0
      dim(df.1)
      min.FC <- 0
      med.FC <- 1
      max.FC <- max(df.1)

      ht <- Heatmap(df.1,
                    heatmap_legend_param = list(
                      title = "Unique\nTCR's",
                      title_gp = gpar(
                        fontsize = input$Legend_size,
                        fontface = "bold", fontfamily = input$font_type
                      ),
                      labels_gp = gpar(fontsize = input$Legend_size, fontfamily = input$font_type)
                    ),
                    col = colorRamp2(c(min.FC,med.FC, max.FC), c(input$min_FC_col, input$med_FC_col, input$max_FC_col)),
                    row_names_gp = grid::gpar(fontsize = input$text_size, fontfamily = input$font_type),
                    column_names_gp = grid::gpar(fontsize = input$text_size, fontfamily = input$font_type),
      )

      # draw(ht, padding = unit(c(5, 5, 5, 5), "mm"))

      draw(ht, padding = unit(c(10, 10, 10, 10), "mm"))
    })

    output$Heatmap_epi_plot <- renderPlot({
      heatmap_epitope()
    })

    output$downloadPlot_Heatmap_epi_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_", input$epitope_umap_selected, "_", input$epitope_umap_selected2, "_Heatmap_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_Heatmap_epi_plot, height = input$height_Heatmap_epi_plot, onefile = FALSE) # open the pdf device
        plot(heatmap_epitope())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_Heatmap_epi_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_", input$epitope_umap_selected, "_", input$epitope_umap_selected2, "_Heatmap_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_Heatmap_epi_plot,
            height = input$height_png_Heatmap_epi_plot,
            res = input$resolution_PNG_Heatmap_epi_plot
        )
        plot(heatmap_epitope())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    #### epitope UMAP interrogation -----
    cols_epitope <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )
      df3.meta <- sc@meta.data

      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta, epi, by = "CDR3_beta", all.x = T)
      df3.meta$selected <- df3.meta[, names(df3.meta) %in% input$epitope_umap_selected]
      df3.meta <- df3.meta[order(df3.meta$selected, decreasing = F), ]
      df3.meta$selected <- factor(df3.meta$selected, levels = unique(df3.meta$selected))

      num <- as.data.frame(unique(df3.meta$selected))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep = "_"), paste(num[i, ]), col.gg[i])
        })
      } else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep = "_"), paste(num[i, ]), col.hcl[i])
        })
      } else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep = "_"), paste(num[i, ]), col.topo[i])
        })
      } else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep = "_"), paste(num[i, ]), heat_col[i])
        })
      } else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep = "_"), paste(num[i, ]), col.terrain[i])
        })
      } else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep = "_"), paste(num[i, ]), palette_rainbow[i])
        })
      } else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep = "_"), paste(num[i, ]), palette1[i])
        })
      } else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope", i, sep = "_"), paste(num[i, ]), input$one.colour.default)
        })
      } # one colour
    })

    output$myPanel_cols_epitope <- renderUI({
      cols_epitope()
    })
    colors_cols_epitope <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )
      df3.meta <- sc@meta.data

      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta, epi, by = "CDR3_beta", all.x = T)
      df3.meta$selected <- df3.meta[, names(df3.meta) %in% input$epitope_umap_selected]


      df3.meta <- df3.meta[order(df3.meta$selected, decreasing = F), ]
      df3.meta$selected <- factor(df3.meta$selected, levels = unique(df3.meta$selected))

      num <- as.data.frame(unique(df3.meta$selected))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.cols_epitope", i, sep = "_")]]
      })
    })

    UMAP_Epitope <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )
      df3.meta <- sc@meta.data

      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta, epi, by = "CDR3_beta")
      df3.meta$selected <- df3.meta[, names(df3.meta) %in% input$epitope_umap_selected]
      df3.meta <- df3.meta[order(df3.meta$selected, decreasing = F), ]
      df3.meta$selected <- factor(df3.meta$selected, levels = unique(df3.meta$selected))


      names(df3.meta)[names(df3.meta) %in% input$Split_group_by_] <- "ID_Column"
      # df3.meta$ID_Column<- df3.meta[df3.meta$ID_Column %in% input$ID_Column_factor,]
      df3.meta$ID_Column <- factor(df3.meta$ID_Column, levels = input$Graph_split_order)

      num <- as.data.frame(unique(df3.meta$selected))
      num <- as.data.frame(num[complete.cases(num) == T, ])
      palette_rainbow <- unlist(colors_cols_epitope())

      df <- ggplot(data = df3.meta, aes(x = UMAP_1, UMAP_2, colour = selected, size = selected, alpha = selected)) +
        geom_point() +
        scale_color_manual(na.value = input$NA_col_analysis, values = palette_rainbow) +
        scale_size_manual(na.value = 0.25, values = rep(input$value_size_epi_umap, dim(num)[1])) +
        scale_alpha_manual(na.value = 0.25, values = rep(1, dim(num)[1])) +
        theme(
          legend.text = element_text(colour = "black", size = 12, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        ) +
        theme_bw()
      if (input$Split_by_group == "no") {
        df <- df
      } else {
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
        paste("_", input$epitope_umap_selected, "_UMAP_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_UMAP_Epitope, height = input$height_UMAP_Epitope, onefile = FALSE) # open the pdf device
        plot(UMAP_Epitope())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_UMAP_Epitope <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_", input$epitope_umap_selected, "_UMAP_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_UMAP_Epitope,
            height = input$height_png_UMAP_Epitope,
            res = input$resolution_PNG_UMAP_Epitope
        )
        plot(UMAP_Epitope())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    # Epitope pie chart function -----

    cols_epitope_pie <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )
      df3.meta <- sc@meta.data
      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta, epi, by = "CDR3_beta")
      df3.meta$selected <- df3.meta[, names(df3.meta) %in% input$epitope_umap_selected]


      df3.meta <- df3.meta[order(df3.meta$selected, decreasing = F), ]
      df3.meta$selected <- factor(df3.meta$selected, levels = unique(df3.meta$selected))

      num <- as.data.frame(unique(df3.meta$selected))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep = "_"), paste(num[i, ]), col.gg[i])
        })
      } else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep = "_"), paste(num[i, ]), col.hcl[i])
        })
      } else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep = "_"), paste(num[i, ]), col.topo[i])
        })
      } else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep = "_"), paste(num[i, ]), heat_col[i])
        })
      } else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep = "_"), paste(num[i, ]), col.terrain[i])
        })
      } else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep = "_"), paste(num[i, ]), palette_rainbow[i])
        })
      } else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep = "_"), paste(num[i, ]), palette1[i])
        })
      } else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_epitope_pie", i, sep = "_"), paste(num[i, ]), input$one.colour.default)
        })
      } # one colour
    })

    output$myPanel_cols_epitope_pie <- renderUI({
      cols_epitope_pie()
    })

    colors_cols_epitope_pie <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )
      df3.meta <- sc@meta.data

      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta, epi, by = "CDR3_beta")
      df3.meta$selected <- df3.meta[, names(df3.meta) %in% input$epitope_umap_selected]

      df3.meta <- df3.meta[order(df3.meta$selected, decreasing = F), ]
      df3.meta$selected <- factor(df3.meta$selected, levels = unique(df3.meta$selected))

      num <- as.data.frame(unique(df3.meta$selected))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.cols_epitope_pie", i, sep = "_")]]
      })
    })

    Pie_chart_Epitope <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )
      df3.meta <- sc@meta.data
      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta, epi, by = "CDR3_beta")
      df3.meta$Selected_function <- df3.meta[, names(df3.meta) %in% input$epitope_umap_selected]
      df3.meta <- df3.meta[order(df3.meta$Selected_function, decreasing = F), ]
      df3.meta$Selected_function <- factor(df3.meta$Selected_function, levels = unique(df3.meta$Selected_function))

      top_BD_cluster <- df3.meta
      top_BD_cluster$Selected_group <- top_BD_cluster[, names(top_BD_cluster) %in% input$epitope_umap_selected2]
      top_BD_cluster$Selected_indiv <- top_BD_cluster[, names(top_BD_cluster) %in% input$Samp_col]

      if (input$by_indiv_pie_epi == "yes") {
        top_BD_cluster <- top_BD_cluster[top_BD_cluster$Selected_indiv %in% input$selected_Indiv, ]
      }

      df.col <- unlist(colors_cols_epitope_pie())

      df3.meta3 <- as.data.frame(table(top_BD_cluster$Selected_group, top_BD_cluster$Selected_function))
      total.condition <- as.data.frame(ddply(df3.meta3, "Var1", numcolwise(sum)))
      dim(total.condition)[1]
      dim(df3.meta3)[1]
      emtpy <- matrix(nrow = dim(df3.meta3)[1], ncol = dim(total.condition)[1])

      for (i in 1:dim(df3.meta3)[1]) {
        emtpy[i, ] <- ifelse(df3.meta3$Var1[i] == total.condition$Var1[1:dim(total.condition)[1]],
                             total.condition[total.condition$Var1 == total.condition$Var1[1:dim(total.condition)[1]], 2], F
        )
      }
      df3.meta3$n <- df3.meta3$Freq / rowSums(emtpy)



      ggplot(df3.meta3, aes(x = "", y = n, fill = Var2, group = Var1)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        theme_void(20) +
        facet_wrap(~Var1, nrow = input$wrap_row) +
        theme(
          legend.key.size = unit(1, "cm"),
          legend.title = element_blank()
        ) +
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 20), values = df.col, na.value = input$NA_col_analysis) +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.position = input$legend_position,
          legend.title = element_blank()
        )
    })

    output$Pie_Epitope_plot <- renderPlot({
      Pie_chart_Epitope()
    })

    Pie_Epitope_dt_process <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )
      df3.meta <- sc@meta.data

      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }

      req(input$epitope_umap_selected, input$epitope_umap_selected2, df3.meta$CDR3_beta)

      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta, epi, by = "CDR3_beta")

      df3.meta$Selected_function <- df3.meta[, names(df3.meta) %in% input$epitope_umap_selected]
      top_BD_cluster <- df3.meta
      top_BD_cluster$Selected_group <- top_BD_cluster[, names(top_BD_cluster) %in% input$epitope_umap_selected2]
      top_BD_cluster$Selected_indiv <- top_BD_cluster[, names(top_BD_cluster) %in% input$Samp_col]
      top_BD_cluster$cloneCount <- 1
      top_BD_cluster2 <- top_BD_cluster[, names(top_BD_cluster) %in% c("cloneCount", "Selected_function", "Selected_group", "Selected_indiv")]
      top_BD_cluster2 <- top_BD_cluster2 %>%
        select(cloneCount, everything())
      df2 <- as.data.frame(ddply(top_BD_cluster2, names(top_BD_cluster2)[-c(1)], numcolwise(sum)))
      df2$fraction <- df2$cloneCount / sum(df2$cloneCount)
      df2$Percent <- round(df2$cloneCount / sum(df2$cloneCount) * 100, 2)
      df2$Selected_function <- ifelse(grepl("NA", df2$Selected_function), "-", df2$Selected_function)
      df2[order(df2$Percent, decreasing = T), ]
    })

    output$Pie_Epitope_dt <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 5, scrollX = TRUE), {
      Pie_Epitope_dt_process()
    })




    # Epitope stat -----
    Epitope_of_interest <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )

      df3.meta <- sc@meta.data
      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }

      checking1 <- merge(df3.meta, epi, by = "CDR3_beta")

      checking2 <- checking1[, names(checking1) %in% c("epitope", "pathology")]
      names.check <- names(checking2)
      checking2$cloneCount <- 1
      sum.check2 <- ddply(checking2, names.check, numcolwise(sum))
      sum.check2_Morethan2 <- subset(sum.check2, sum.check2$cloneCount > 1)
      sum.check2_Morethan2 <- sum.check2_Morethan2[order(sum.check2_Morethan2$cloneCount, decreasing = T), ]
      sum.check2_Morethan2
    })

    output$Epi_of_interest_DF <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 5, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )
      as.data.frame(Epitope_of_interest())
    })

    observeEvent(input$Update_epi, {
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )

      sum.check2_Morethan2 <- as.data.frame(Epitope_of_interest())

      updateSelectInput(
        session,
        "Epi_of_interest",
        choices = sum.check2_Morethan2$epitope,
        selected = sum.check2_Morethan2$epitope[1]
      )
    })

    compare.stat_Epi <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )
      req(input$Epi_of_interest)
      df3.meta <- sc@meta.data
      as.data.frame(Epitope_of_interest())
      epi <- epi[epi$epitope %in% input$Epi_of_interest, ]
      epi <- epi[names(epi) %in% c("CDR3_beta", "epitope", "pathology")]
      epi <- unique(epi)

      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      df3.meta
      checking <- merge(df3.meta, epi, by = "CDR3_beta", all.x = T)
      rownames(checking) <- checking$Cell_Index
      name.check.epi <- unlist(unique(epi$epitope))
      checking$epi_selected <- ifelse(as.character(checking$epitope) == name.check.epi, name.check.epi, "not selected")
      rownames(checking) <- (checking$Cell_Index)
      checking <- checking[order(checking$order, decreasing = F), ]
      checking
      sc@meta.data <- checking
      Idents(object = sc) <- sc@meta.data$epi_selected

      min.pct.expression <- input$min_point_ # standard setting: 0.25
      min.logfc <- input$LogFC_ # 0.25 is standard

      cluster.names <- unique(Idents(sc))[order(unique(Idents(sc)))]

      markers.fm.list <- FindMarkers(sc, ident.1 = name.check.epi, min.pct = min.pct.expression, logfc.threshold = min.logfc, only.pos = TRUE)
      markers.fm.list2 <- subset(markers.fm.list, markers.fm.list$p_val_adj < input$pval.ex.filter)
      as.data.frame(markers.fm.list2)
    })

    output$compare.stat_Epi_DT <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )
      compare.stat_Epi()
    })

    # Epitope dotplot ----
    all_expression_plot_epi <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )
      req(input$Epi_of_interest)
      df3.meta <- sc@meta.data
      as.data.frame(Epitope_of_interest())
      epi <- epi[epi$epitope %in% input$Epi_of_interest, ]
      epi <- epi[names(epi) %in% c("CDR3_beta", "epitope", "pathology")]
      epi <- unique(epi)

      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      df3.meta
      checking <- merge(df3.meta, epi, by = "CDR3_beta", all.x = T)

      checking

      rownames(checking) <- checking$Cell_Index
      name.check.epi <- as.character(unlist(unique(epi$epitope)))
      checking$epitope <- as.character(checking$epitope)

      checking$epi_selected <- ifelse(as.character(checking$epitope) != name.check.epi, "NS", name.check.epi)
      checking$epi_selected[is.na(checking$epi_selected)] <- "NS"

      checking <- checking[order(checking$order, decreasing = F), ]
      checking
      sc@meta.data <- checking
      Idents(object = sc) <- sc@meta.data$epi_selected

      if (input$restrict.dotpot.epi == F) {
        list.names <- rownames(compare.stat_Epi())
      } else {
        list.names <- rownames(compare.stat_Epi())
        list.names <- list.names[1:input$restrict.dotpot.num.epi]
      }

      size_legend <- input$Legend_size - 2

      DotPlot(sc, features = list.names) +
        RotatedAxis() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
          axis.title.x = element_blank(),
          legend.title = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.text = element_text(colour = "black", size = size_legend, family = input$font_type),
          legend.position = input$legend_position,
        ) +
        scale_colour_gradient2(low = input$low.dotplot.epi, mid = input$middle.dotplot.epi, high = input$high.dotplot.epi) +
        scale_x_discrete(labels = label_wrap(20)) +
        scale_y_discrete(labels = label_wrap(20))
    })

    output$checking_epi_dot_issue <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1, 2, 5, 10, 20, 50, 100), pageLength = 20, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )
      req(input$Epi_of_interest)
      df3.meta <- sc@meta.data
      as.data.frame(Epitope_of_interest())
      epi <- epi[epi$epitope %in% input$Epi_of_interest, ]
      epi <- epi[names(epi) %in% c("CDR3_beta", "epitope", "pathology")]
      epi <- unique(epi)

      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      df3.meta
      checking <- merge(df3.meta, epi, by = "CDR3_beta", all.x = T)

      checking

      rownames(checking) <- checking$Cell_Index
      name.check.epi <- as.character(unlist(unique(epi$epitope)))
      checking$epitope <- as.character(checking$epitope)

      checking$epi_selected <- ifelse(as.character(checking$epitope) != name.check.epi, "NS", name.check.epi)
      checking$epi_selected[is.na(checking$epi_selected)] <- "NS"

      checking <- checking[order(checking$order, decreasing = F), ]
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
        paste(input$Epi_of_interest, "_dotplot", "_", today(), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_all_expression_dotplot_epi, height = input$height_all_expression_dotplot_epi, onefile = FALSE) # open the pdf device
        df <- all_expression_plot_epi()
        plot(df)
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_all_expression_dotplot_epi <- downloadHandler(
      filename = function() {
        paste(input$Epi_of_interest, "_dotplot", "_", today(), ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_all_expression_dotplot_epi,
            height = input$height_png_all_expression_dotplot_epi,
            res = input$resolution_PNG_all_expression_dotplot_epi
        )
        df <- all_expression_plot_ex()

        plot(df)
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )
    # Over representation analysis for Epitope  -----
    Over_rep_Epi_old <- reactive({
      sc <- UMAP_metadata_with_labs()
      epi <- data_sc_TCRex()
      validate(
        need(
          nrow(sc) > 0 & nrow(epi) > 0,
          "Upload Files"
        )
      )
      req(input$Epi_of_interest)
      df3.meta <- sc@meta.data
      as.data.frame(Epitope_of_interest())
      epi <- epi[epi$epitope %in% input$Epi_of_interest, ]
      epi <- epi[names(epi) %in% c("CDR3_beta", "epitope", "pathology")]
      epi <- unique(epi)

      if (input$datasource == "BD_Rhapsody_Paired") {
        df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
      } else if (input$datasource == "BD_Rhapsody_AIRR") {
        df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
      } else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      df3.meta
      checking <- merge(df3.meta, epi, by = "CDR3_beta", all.x = T)

      checking

      rownames(checking) <- checking$Cell_Index
      name.check.epi <- as.character(unlist(unique(epi$epitope)))
      checking$epitope <- as.character(checking$epitope)

      checking$epi_selected <- ifelse(as.character(checking$epitope) != name.check.epi, "NS", name.check.epi)
      checking$epi_selected[is.na(checking$epi_selected)] <- "NS"

      checking <- checking[order(checking$order, decreasing = F), ]
      sc@meta.data <- checking
      Idents(object = sc) <- sc@meta.data$epi_selected

      geneSet <- read.csv(system.file("OverRep", "GeneSets.csv", package = "STEGO.R"), header = T)

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

      if (input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        geneSet$GeneSet <- gsub("-", ".", geneSet$GeneSet)
      }

      if (input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        require(stringr)
        geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
      }

      for (i in 1:dim(geneSet)[1]) {
        # listed GeneSet
        message(paste("GeneSet: ", i))
        Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet, ";")[i])
        names(Gene.set.testing) <- "V1"
        Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
        names(Gene.set.testing2) <- "V1"
        # message(paste(dim(Gene.set.testing2)[1],"GeneSet total"))
        # background genes
        background.overlap <- merge(Gene.set.testing2, background.genes.name, by = "V1")
        # message(paste(dim(background.overlap)[1],"in Background"))
        geneSet$background.geneset[i] <- length(background.overlap$V1)
        geneSet$background.geneset.name[i] <- as.character(paste(unlist(background.overlap[1]), collapse = ";"))
        # in sig gene list
        overlap <- merge(background.overlap, DEx.genes, by = "V1")
        # message(paste(dim(overlap)[1],"# Sig genes"))
        geneSet$in.geneset[i] <- length(overlap$V1)
        geneSet$in.geneset.name[i] <- as.character(paste(unlist(overlap[1]), collapse = ";"))
      }

      geneSet2 <- subset(geneSet, geneSet$in.geneset > 0)

      for (i in 1:dim(geneSet2)[1]) {
        tota.gene.set <- geneSet2$background.geneset[i] # genes that are identified in background
        tota.gene.set
        in.geneset <- geneSet2$in.geneset[i] # DEx in geneset

        background.genes
        not.in.total <- background.genes - tota.gene.set
        not.in.geneset.sig <- total.sig - in.geneset
        d <- data.frame(gene.in.interest = c(in.geneset, not.in.geneset.sig), gene.not.interest = c(tota.gene.set, not.in.total))
        row.names(d) <- c("In_category", "not_in_category")

        if (in.geneset > 0) {
          geneSet2$p.val[i] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
          geneSet2$lowerCI[i] <- unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
          geneSet2$upperCI[i] <- unlist(fisher.test(d)$conf.int)[2]
          geneSet2$OR[i] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1], 3)
        } else {
          geneSet2$p.value[i] <- "-"
          geneSet2$lowerCI[i] <- "-"
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
        need(
          nrow(geneSet2) > 0,
          "Upload Files"
        )
      )
      geneSet2 <- geneSet2[order(geneSet2$p.val, decreasing = F), ]
      geneSet2 <- subset(geneSet2, geneSet2$in.geneset >= input$in.geneset.cutoff_Epi)
      geneSet2 <- subset(geneSet2, geneSet2$p.val <= input$p.val_cutoff_Epi)

      geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
      geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
      name.list <- c("Geneset_ID", "p.val", "FDR", "Bonferroni", "OR", "lowerCI", "upperCI", "in.geneset.name", "in.geneset", "background.geneset", "total.sig", "background.genes", "background.geneset.name")
      geneSet2 <- geneSet2 %>%
        select(all_of(name.list), everything())
    })

    output$Over_rep_Epi_Tab <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1, 2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      Over_rep_Epi()
    })

    output$downloadtb_over.rep.Epi <- downloadHandler(
      filename = function() {
        paste0("Epitope", "_", today(), "_over_rep.csv")
      },
      content = function(file) {
        df <- as.data.frame(Over_rep_Epi())
        write.csv(df, file, row.names = F)
      }
    )


    #### download Epitope files ------
    output$downloadPlot_Pie_Epitope <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$epitope_Pie_Epitope, "_Pie_epi_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_Pie_Epitope, height = input$height_Pie_Epitope, onefile = FALSE) # open the pdf device
        plot(Pie_chart_Epitope())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_Pie_Epitope <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$epitope_umap_selected, "_Pie_epi_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_Pie_Epitope,
            height = input$height_png_Pie_Epitope,
            res = input$resolution_PNG_Pie_Epitope
        )
        plot(Pie_chart_Epitope())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )


    output$downloadPlot_UMAP_Epitope <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$epitope_umap_selected, "_UMAP_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_UMAP_Epitope, height = input$height_UMAP_Epitope, onefile = FALSE) # open the pdf device
        plot(UMAP_Epitope())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_UMAP_Epitope <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$epitope_umap_selected, "_UMAP_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_UMAP_Epitope,
            height = input$height_png_UMAP_Epitope,
            res = input$resolution_PNG_UMAP_Epitope
        )
        plot(UMAP_Epitope())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    output$downloaddf_Pie_Epitope_dt <- downloadHandler(
      filename = function() {
        paste("Summary_table_", gsub("-", ".", Sys.Date()), ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(Pie_Epitope_dt_process())
        write.csv(df, file, row.names = F)
      }
    )



    #### clusTCR2 seting up the dataframe -----
    AG_cluster <- reactive({
      x <- today()
      clust <- data_sc_clusTCR_AG()
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0 & nrow(clust) > 0,
          "Upload File"
        )
      )
      req(sc)
      md <- sc@meta.data

      req(clust, md)
      if (input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") {
        names(md)[names(md) %in% "v_gene_AG"] <- "Selected_V_AG"
        names(md)[names(md) %in% "junction_aa_AG"] <- "AminoAcid_AG"
      } else {
        names(md)[names(md) %in% "v_gene_AG"] <- "Selected_V_AG"
        names(md)[names(md) %in% "cdr3_AG"] <- "AminoAcid_AG"
      }

      req(md$Selected_V_AG, md$AminoAcid_AG, input$Samp_col)

      md$CDR3_Vgene <- paste(md$AminoAcid_AG, md$Selected_V_AG, sep = "_")
      df <- merge(md, clust, by = "CDR3_Vgene")

      # Total clone Count -----
      df2 <- as.data.frame(df$Clust_size_order)
      names(df2) <- "Clust_size_order"
      df2$Total_cloneCount <- 1
      df3 <- as.data.frame(ddply(df2, "Clust_size_order", numcolwise(sum)))
      df3 <- subset(df3, df3$Total_cloneCount > 2)
      df4 <- merge(df3, df, by = "Clust_size_order")
      df4
      # updated cluster count
      clusterCount <- df4[, names(df4) %in% c("Clust_size_order", "CDR3_Vgene")]
      clusterCount <- clusterCount[!duplicated(clusterCount), ]
      clusterCount$ClusterCount <- 1
      clusterCount2 <- as.data.frame(ddply(clusterCount, c("Clust_size_order"), numcolwise(sum)))
      clusterCount2 <- subset(clusterCount2, clusterCount2$ClusterCount > 1)
      df4_clusterCount2 <- merge(df4, clusterCount2, by = c("Clust_size_order"))
      df4_clusterCount2
      # Sample Count
      SampCount <- df4_clusterCount2[, names(df4_clusterCount2) %in% c("Clust_size_order", input$Samp_col)]
      SampCount <- SampCount[!duplicated(SampCount), ]
      SampCount$Sample_count <- 1
      df9 <- as.data.frame(ddply(SampCount, c("Clust_size_order"), numcolwise(sum)))
      df10 <- merge(df9, df4_clusterCount2, by = "Clust_size_order")

      # Calculating the priority
      df10$priority <- 1 / (df10$Total_cloneCount * df10$ClusterCount * df10$Sample_count)
      df10 <- df10[order(df10$priority, decreasing = F), ]
      df10
      df10$priority[is.na(df10$priority)] <- 0
      df10 <- subset(df10, df10$priority > 0)
      df10
      # updated order
      updatedOrder <- df10[, names(df10) %in% c("priority", "Clust_size_order")]
      updatedOrder <- updatedOrder[!duplicated(updatedOrder), ]
      updatedOrder <- updatedOrder[order(updatedOrder$priority, decreasing = F), ]
      updatedOrder$Updated_order <- 1:dim(updatedOrder)[1]
      updatedOrder

      # final data frame
      df7 <- merge(updatedOrder, df10, by = c("Clust_size_order", "priority"))
      df7 <- df7[order(df7$priority, decreasing = F), ]


      clusterAG <- df7 %>%
        select(all_of(c(input$Samp_col, "Sample_count", "Total_cloneCount", "ClusterCount", "priority", "Updated_order")), everything())
    })
    BD_cluster <- reactive({
      sc <- UMAP_metadata_with_labs()
      clust <- data_sc_clusTCR_BD()

      validate(
        need(
          nrow(sc) > 0 & nrow(clust) > 0,
          "Upload File"
        )
      )
      md <- sc@meta.data
      x <- today()

      if (input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") {
        names(md)[names(md) %in% "v_gene_BD"] <- "Selected_V_BD"
        names(md)[names(md) %in% "junction_aa_BD"] <- "AminoAcid_BD"
      } else {
        names(md)[names(md) %in% "v_gene_BD"] <- "Selected_V_BD"
        names(md)[names(md) %in% "cdr3_BD"] <- "AminoAcid_BD"
      }

      req(md$Selected_V_BD, md$AminoAcid_BD, input$Samp_col)
      md$CDR3_Vgene <- paste(md$AminoAcid_BD, md$Selected_V_BD, sep = "_")
      df <- merge(md, clust, by = "CDR3_Vgene")

      # Total clone Count -----
      df2 <- as.data.frame(df$Clust_size_order)
      names(df2) <- "Clust_size_order"
      df2$Total_cloneCount <- 1
      df3 <- as.data.frame(ddply(df2, "Clust_size_order", numcolwise(sum)))
      df3 <- subset(df3, df3$Total_cloneCount > 2)
      df4 <- merge(df3, df, by = "Clust_size_order")
      df4
      # updated cluster count
      clusterCount <- df4[, names(df4) %in% c("Clust_size_order", "CDR3_Vgene")]
      clusterCount <- clusterCount[!duplicated(clusterCount), ]
      clusterCount$ClusterCount <- 1
      clusterCount2 <- as.data.frame(ddply(clusterCount, c("Clust_size_order"), numcolwise(sum)))
      clusterCount2 <- subset(clusterCount2, clusterCount2$ClusterCount > 1)
      df4_clusterCount2 <- merge(df4, clusterCount2, by = c("Clust_size_order"))
      df4_clusterCount2
      # Sample Count
      SampCount <- df4_clusterCount2[, names(df4_clusterCount2) %in% c("Clust_size_order", input$Samp_col)]
      SampCount <- SampCount[!duplicated(SampCount), ]
      SampCount$Sample_count <- 1
      df9 <- as.data.frame(ddply(SampCount, c("Clust_size_order"), numcolwise(sum)))
      df10 <- merge(df9, df4_clusterCount2, by = "Clust_size_order")

      # Calculating the priority
      df10$priority <- 1 / (df10$Total_cloneCount * df10$ClusterCount * df10$Sample_count)
      df10 <- df10[order(df10$priority, decreasing = F), ]
      df10$priority[is.na(df10$priority)] <- 0
      df10 <- subset(df10, df10$priority > 0)
      # updated order
      updatedOrder <- df10[, names(df10) %in% c("priority", "Clust_size_order")]
      updatedOrder <- updatedOrder[!duplicated(updatedOrder), ]
      updatedOrder <- updatedOrder[order(updatedOrder$priority, decreasing = F), ]
      updatedOrder$Updated_order <- 1:dim(updatedOrder)[1]
      updatedOrder

      # final data frame
      df7 <- merge(updatedOrder, df10, by = c("Clust_size_order", "priority"))
      df7 <- df7[order(df7$priority, decreasing = F), ]

      clusterBD <- df7 %>%
        select(all_of(c(input$Samp_col, "Sample_count", "Total_cloneCount", "ClusterCount", "priority", "Updated_order")), everything())
      clusterBD
    })

    clusTCR2_df <- reactive({
      if (input$chain_TCR == "TRAG") {
        if (length(data_sc_clusTCR_AG()) > 0) {
          req(AG_cluster())
          AG_cluster()
        }
      } else if (input$chain_TCR == "TRBD") {
        if (length(data_sc_clusTCR_BD()) > 0) {
          req(BD_cluster())
          BD_cluster()
        }
      } else {

      }
    })

    output$Tb_ClusTCR_selected <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 5, scrollX = TRUE), {
      cluster <- clusTCR2_df()
      validate(
        need(
          nrow(cluster) > 0,
          "Upload clusTCR table, which is needed for TCR -> UMAP section"
        )
      )
      cluster
    })

    output$downloadtb_Tb_ClusTCR_selected <- downloadHandler(
      filename = function() {
        paste0(input$Clusters_to_dis_PIE, "_Cluster_", today(), "_over_rep.csv")
      },
      content = function(file) {
        df <- as.data.frame(clusTCR2_df())
        write.csv(df, file, row.names = F)
      }
    )

    # cluster to display ------
    observe({
      cluster <- clusTCR2_df()

      validate(
        need(
          nrow(cluster) > 0,
          "upload clustering"
        )
      )
      cluster$ID_Column <- cluster[, names(cluster) %in% input$Samp_col]
      req(input$Samp_col)
      # cluster <- cluster[cluster$Clust_size_order %in% input$lower_cluster:input$upper_cluster,]
      updateSelectInput(
        session,
        "Clusters_to_dis_PIE",
        choices = unique(cluster$Updated_order),
        selected = unique(cluster$Updated_order)[1]
      )
    }) # cluster to display


    # cluster UMAP (1 to display) ClusTCR -----
    cols_clust_UMAP <- reactive({
      cluster <- clusTCR2_df()
      validate(
        need(
          nrow(cluster) > 0,
          "Upload clusTCR table"
        )
      )
      req(cluster, input$Colour_By_this, input$Samp_col)
      cluster$ID_Column <- cluster[, names(cluster) %in% input$Samp_col]
      cluster <- cluster[order(cluster$ID_Column), ]

      cluster$colour <- cluster[, names(cluster) %in% input$Colour_By_this]
      cluster$colour <- gsub("_", " ", cluster$colour)
      cluster$colour <- factor(cluster$colour, levels = unique(cluster$colour))
      cluster$colour <- gsub("NA", NA, cluster$colour)

      num <- as.data.frame(unique(cluster$colour))
      num <- as.data.frame(num[complete.cases(num) == T, ])
      num

      col.gg <- gg_fill_hue(dim(num)[1])
      palette_rainbow <- rainbow(dim(num)[1])
      heat_col <- heat.colors(dim(num)[1])
      col.terrain <- terrain.colors(dim(num)[1])
      col.topo <- topo.colors(dim(num)[1])
      col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")

      if (input$colourtype == "default") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep = "_"), paste(num[i, ]), col.gg[i])
        })
      } else if (input$colourtype == "hcl.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep = "_"), paste(num[i, ]), col.hcl[i])
        })
      } else if (input$colourtype == "topo.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep = "_"), paste(num[i, ]), col.topo[i])
        })
      } else if (input$colourtype == "heat.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep = "_"), paste(num[i, ]), heat_col[i])
        })
      } else if (input$colourtype == "terrain.colors") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep = "_"), paste(num[i, ]), col.terrain[i])
        })
      } else if (input$colourtype == "rainbow") {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep = "_"), paste(num[i, ]), palette_rainbow[i])
        })
      } else if (input$colourtype == "random") {
        palette1 <- distinctColorPalette(dim(num)[1])
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep = "_"), paste(num[i, ]), palette1[i])
        })
      } else {
        lapply(1:dim(num)[1], function(i) {
          colourInput(paste("col.cols_clust_UMAP", i, sep = "_"), paste(num[i, ]), input$one.colour.default)
        })
      } # one colour
    })
    output$myPanel_cols_clust_UMAP <- renderUI({
      cols_clust_UMAP()
    })
    colors_cols_cols_clust_UMAP <- reactive({
      cluster <- clusTCR2_df()
      validate(
        need(
          nrow(cluster) > 0,
          "Upload clusTCR table"
        )
      )
      req(cluster, input$Colour_By_this)
      cluster$ID_Column <- cluster[, names(cluster) %in% input$Samp_col]
      cluster <- cluster[order(cluster$ID_Column), ]
      # cluster <- cluster[cluster$Updated_order %in% input$Clusters_to_dis_PIE,]
      cluster$colour <- cluster[, names(cluster) %in% input$Colour_By_this]
      cluster$colour <- gsub("_", " ", cluster$colour)
      cluster$colour <- factor(cluster$colour, levels = unique(cluster$colour))
      cluster$colour <- gsub("NA", NA, cluster$colour)

      num <- as.data.frame(unique(cluster$colour))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.cols_clust_UMAP", i, sep = "_")]]
      })
    })



    UMAP_ClusTCR2 <- reactive({
      cluster <- clusTCR2_df()
      validate(
        need(
          nrow(cluster) > 0,
          "Upload clusTCR table, which is needed for TCR -> UMAP section"
        )
      )
      req(cluster, input$Clusters_to_dis_PIE, input$Colour_By_this)
      cluster$ID_Column <- cluster[, names(cluster) %in% input$Samp_col]
      cluster <- cluster[order(cluster$ID_Column), ]
      # cluster <- cluster[cluster$Updated_order %in% input$Clusters_to_dis_PIE,]
      cluster$colour <- cluster[, names(cluster) %in% input$Colour_By_this]
      cluster$colour <- gsub("_", " ", cluster$colour)
      num <- as.data.frame(unique(cluster$colour))
      num <- as.data.frame(num[complete.cases(num) == T, ])
      names(num) <- "ID"
      num$col <- unlist(colors_cols_cols_clust_UMAP())


      cluster$colour <- gsub("NA", NA, cluster$colour)
      cluster <- cluster[cluster$Updated_order %in% input$Clusters_to_dis_PIE, ]

      num <- num[num$ID %in% unique(cluster$colour), ]
      print(num)

      cluster$colour <- factor(cluster$colour, levels = num$ID)


      figure <- ggplot(data = cluster, aes(x = UMAP_1, UMAP_2, colour = colour)) +
        geom_point(size = input$size.dot.umap) +
        scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 20), values = num$col, na.value = input$NA_col_analysis) +
        theme_bw() +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
          # strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 0),
          axis.title.x = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
        )

      if (input$Split_by_group == "no") {
        figure
      } else {
        figure + facet_wrap(~ID_Column, nrow = input$wrap_row)
      }
    })

    output$UMAP_ClusTCR2_plot <- renderPlot({
      UMAP_ClusTCR2()
    })

    output$downloadPlot_UMAP_ClusTCR2_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_UMAP_ClusTCR2_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_UMAP_ClusTCR2_plot, height = input$height_UMAP_ClusTCR2_plot, onefile = FALSE) # open the pdf device
        plot(UMAP_ClusTCR2())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_UMAP_ClusTCR2_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("_UMAP_ClusTCR2_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_UMAP_ClusTCR2_plot,
            height = input$height_png_UMAP_ClusTCR2_plot,
            res = input$resolution_PNG_UMAP_ClusTCR2_plot
        )
        plot(UMAP_ClusTCR2())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    ##### Clustering motif plot ----

    motif_plot_sc <- reactive({
      cluster <- clusTCR2_df()
      validate(
        need(
          nrow(cluster) > 0,
          "Upload clusTCR table"
        )
      )
      cluster$ID_Column <- cluster[, names(cluster) %in% input$Samp_col]
      Network_df <- cluster[order(cluster$Updated_order), ]
      # Network_df2 <-
      Network_df %>% distinct(CDR3_Vgene, .keep_all = TRUE) # make Unique
      Network_df <- Network_df %>% distinct(CDR3_Vgene, .keep_all = TRUE) # make Unique
      Motif_from_cluster_file(Network_df, Clust_selected = input$Clusters_to_dis_PIE, selected_cluster_column = "Updated_order")
      # ?Motif_from_cluster_file
    })
    output$Motif_ClusTCR2_cluster <- renderPlot({
      motif_plot_sc()
    })

    # render which cases were contributing to the cluster
    output$print_unique_cases <- renderPrint({
      cluster <- clusTCR2_df()
      validate(
        need(
          nrow(cluster) > 0,
          "Upload clusTCR table, which is needed for TCR -> UMAP section"
        )
      )
      cluster$ID_Column <- cluster[, names(cluster) %in% input$Samp_col]

      df <- cluster[cluster$Updated_order %in% input$Clusters_to_dis_PIE, ]


      if (input$chain_TCR == "TRAG") {
        cat("This Motif if from: ", noquote(unique(df$ID_Column)), "; the TCR is:", noquote(unique(df$Selected_V_AG)), noquote(unique(df$j_gene_AG)))
      } else if (input$chain_TCR == "TRBD") {
        cat("This Motif if from:", noquote(unique(df$ID_Column)), "; the TCR is:", noquote(unique(df$Selected_V_BD)), noquote(unique(df$d_gene_BD)), noquote(unique(df$j_gene_BD)))
      } else { # BCR repertoire
      }
    })
    # download the motif plot -----
    output$downloadPlot_Motif_ClusTCR2_cluster <- downloadHandler(
      filename = function() {
        x <- today()
        paste("Motif_ClusTCR2_cluster_plot_", input$chain_TCR, "_", input$Clusters_to_dis_PIE, "_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_Motif_ClusTCR2_cluster, height = input$height_Motif_ClusTCR2_cluster, onefile = FALSE) # open the pdf device
        plot(motif_plot_sc())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_Motif_ClusTCR2_cluster <- downloadHandler(
      filename = function() {
        x <- today()
        paste("Motif_ClusTCR2_cluster_plot_", input$chain_TCR, "_", input$Clusters_to_dis_PIE, "_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_Motif_ClusTCR2_cluster,
            height = input$height_png_Motif_ClusTCR2_cluster,
            res = input$resolution_PNG_Motif_ClusTCR2_cluster
        )
        plot(motif_plot_sc())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    ## heatmap for each cluster ----
    Selected_cluster_plot <- reactive({
      cluster <- clusTCR2_df()
      validate(
        need(
          nrow(cluster) > 0,
          "Upload clusTCR table, which is needed for TCR -> UMAP section"
        )
      )
      cluster$ID_Column2 <- cluster[, names(cluster) %in% input$Split_group_by_]
      Network_df <- cluster[order(cluster$Updated_order), ]
      Network_df <- Network_df[Network_df$Updated_order %in% input$Clusters_to_dis_PIE, ]
      Network_df$Selected <- Network_df[, names(Network_df) %in% input$Colour_By_this]
      Network_df$cloneCount <- 1
      df.1 <- acast(Network_df, ID_Column2 ~ Selected, value.var = "cloneCount")
      df.1[is.na(df.1)] <- 0
      min.FC <- min(df.1)
      med.FC <- max(df.1) / 2
      max.FC <- max(df.1)

      ht <- Heatmap(df.1,
                    heatmap_legend_param = list(
                      title = "Count",
                      title_gp = gpar(
                        fontsize = input$Legend_size,
                        fontface = "bold", fontfamily = input$font_type
                      ),
                      labels_gp = gpar(fontsize = input$Legend_size, fontfamily = input$font_type)
                    ),
                    col = colorRamp2(c(min.FC,med.FC, max.FC), c(input$min_FC_col, input$med_FC_col, input$max_FC_col)),
                    row_names_gp = grid::gpar(fontsize = input$text_size, fontfamily = input$font_type),
                    column_names_gp = grid::gpar(fontsize = input$text_size, fontfamily = input$font_type),
      )

      draw(ht, padding = unit(c(10, 10, 10, 10), "mm"))
    })


    #### cluster pie chart function -----
    # cols_clusTCR2_pie <- reactive({
    #   cluster <- clusTCR2_df()
    #   validate(
    #     need(nrow(cluster)>0,
    #          "Upload clusTCR table, which is needed for TCR -> UMAP section")
    #   )
    #   cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]
    #   Network_df <- cluster[order(cluster$Updated_order),]
    #   Network_df <- Network_df[Network_df$Updated_order  %in% input$Clusters_to_dis_PIE,]
    #   Network_df$Selected <- Network_df[,names(Network_df) %in% input$Colour_By_this]
    #   Network_df$Selected <- gsub("NA",NA,Network_df$Selected)
    #
    #   num <- as.data.frame(unique(Network_df$Selected))
    #   num <- as.data.frame(num[complete.cases(num)==T,])
    #
    #
    #   col.gg <- gg_fill_hue(dim(num)[1])
    #   palette_rainbow <- rainbow(dim(num)[1])
    #   heat_col <- heat.colors(dim(num)[1])
    #   col.terrain <- terrain.colors(dim(num)[1])
    #   col.topo <- topo.colors(dim(num)[1])
    #   col.hcl <- hcl.colors(dim(num)[1], palette = "viridis")
    #
    #   if (input$colourtype == "default") {
    #     lapply(1:dim(num)[1], function(i) {
    #       colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), col.gg[i])
    #     })
    #   }
    #   else if (input$colourtype == "hcl.colors") {
    #     lapply(1:dim(num)[1], function(i) {
    #       colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), col.hcl[i])
    #     })
    #   }
    #   else if (input$colourtype == "topo.colors") {
    #     lapply(1:dim(num)[1], function(i) {
    #       colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), col.topo[i])
    #     })
    #   }
    #   else if (input$colourtype == "heat.colors") {
    #     lapply(1:dim(num)[1], function(i) {
    #       colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), heat_col[i])
    #     })
    #   }
    #   else if (input$colourtype == "terrain.colors") {
    #     lapply(1:dim(num)[1], function(i) {
    #       colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), col.terrain[i])
    #     })
    #   }
    #
    #   else if (input$colourtype == "rainbow") {
    #     lapply(1:dim(num)[1], function(i) {
    #       colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), palette_rainbow[i])
    #     }) }
    #   else if (input$colourtype == "random") {
    #     palette1 <- distinctColorPalette(dim(num)[1])
    #     lapply(1:dim(num)[1], function(i) {
    #       colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), palette1[i])
    #     })
    #
    #   }
    #   else {
    #     lapply(1:dim(num)[1], function(i) {
    #       colourInput(paste("col.cols_clusTCR2_pie", i, sep="_"), paste(num[i,]), input$one.colour.default)
    #     })
    #
    #
    #   } # one colour
    #
    # })
    #
    # output$myPanel_cols_clusTCR2_pie <- renderUI({cols_clusTCR2_pie()})
    #
    # colors_cols_clusTCR2_pie <- reactive({
    #   cluster <- clusTCR2_df()
    #   validate(
    #     need(nrow(cluster)>0,
    #          "Upload clusTCR table, which is needed for TCR -> UMAP section")
    #   )
    #   cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]
    #   Network_df <- cluster[order(cluster$Updated_order),]
    #   Network_df <- Network_df[Network_df$Updated_order  %in% input$Clusters_to_dis_PIE,]
    #   Network_df$Selected <- Network_df[,names(Network_df) %in% input$Colour_By_this]
    #   Network_df$Selected <- gsub("NA",NA,Network_df$Selected)
    #
    #   num <- as.data.frame(unique(Network_df$Selected))
    #   num <- as.data.frame(num[complete.cases(num)==T,])
    #
    #
    #   lapply(1:dim(num)[1], function(i) {
    #     input[[paste("col.cols_clusTCR2_pie", i, sep="_")]]
    #   })
    # })
    #
    # Pie_chart_ClusTCR2 <- reactive({
    #   cluster <- clusTCR2_df()
    #   validate(
    #     need(nrow(cluster)>0,
    #          "Upload clusTCR table, which is needed for TCR -> UMAP section")
    #   )
    #   cluster$ID_Column <- cluster[,names(cluster) %in% input$Samp_col]
    #   Network_df <- cluster[order(cluster$Updated_order),]
    #   Network_df <- Network_df[Network_df$Updated_order  %in% input$Clusters_to_dis_PIE,]
    #   Network_df$Selected <- Network_df[,names(Network_df) %in% input$Colour_By_this]
    #   Network_df$Selected <- gsub("NA",NA,Network_df$Selected)
    #
    #   df3.meta3 <-  as.data.frame(table(Network_df$ID_Column,Network_df$Selected))
    #   total.condition <- as.data.frame(ddply(df3.meta3,"Var1",numcolwise(sum)))
    #   dim(total.condition)[1]
    #   dim(df3.meta3)[1]
    #   emtpy <- matrix(nrow =dim(df3.meta3)[1],ncol=dim(total.condition)[1])
    #
    #   for (i in 1:dim(df3.meta3)[1]) {
    #
    #     emtpy[i,] <- ifelse(df3.meta3$Var1[i]==total.condition$Var1[1:dim(total.condition)[1]],
    #         total.condition[total.condition$Var1==total.condition$Var1[1:dim(total.condition)[1]],2],F)
    #   }
    #   df3.meta3$n <- df3.meta3$Freq/rowSums(emtpy)
    #   # df3.meta3
    #   df.col <- unlist(colors_cols_clusTCR2_pie())
    #
    #   ggplot(df3.meta3,aes(x="", y=n, fill=Var2, group = Var1)) +
    #     geom_bar(stat="identity", width=1) +
    #     coord_polar("y", start=0)  +
    #     theme_void(20) +
    #     facet_wrap(~Var1, nrow = input$wrap_row) +
    #     theme(
    #       legend.key.size = unit(1, 'cm'),
    #       legend.title = element_blank()) +
    #     scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 20),values = df.col, na.value = input$NA_col_analysis) +
    #     theme(strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
    #           legend.text = element_text(colour="black", size=input$Legend_size,family=input$font_type),
    #           legend.position = input$legend_position,
    #           legend.title = element_blank()
    #     )
    #
    # })
    ######
    output$Pie_ClusTCR2_plot <- renderPlot({
      # Pie_chart_ClusTCR2()
      Selected_cluster_plot()
    })

    output$downloadPlot_Pie_ClusTCR2_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("Pie_ClusTCR2_plot_", input$chain_TCR, "_", input$Clusters_to_dis_PIE, "_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_Pie_ClusTCR2_plot, height = input$height_Pie_ClusTCR2_plot, onefile = FALSE) # open the pdf device
        plot(Selected_cluster_plot())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_Pie_ClusTCR2_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("Pie_ClusTCR2_plot_", input$chain_TCR, "_", input$Clusters_to_dis_PIE, "_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_Pie_ClusTCR2_plot,
            height = input$height_png_Pie_ClusTCR2_plot,
            res = input$resolution_PNG_Pie_ClusTCR2_plot
        )
        plot(Selected_cluster_plot())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )


    # Cluster stats   -----

    compare.stat_Cluster <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )

      # req(input$V_call_clust_sc,input$junction_clust_sc)
      md <- sc@meta.data

      cluster <- clusTCR2_df()

      cluster$ID_Column <- cluster[, names(cluster) %in% input$Samp_col]
      cluster <- cluster[order(cluster$Updated_order), ]

      rownames(cluster) <- cluster$Cell_Index

      checking <- cluster[, names(cluster) %in% c("Updated_order", "Cell_Index")]
      md.checking <- merge(md, checking, by = "Cell_Index", all.x = T)
      md.checking <- md.checking[order(md.checking$order), ]
      rownames(md.checking) <- md.checking$Cell_Index

      md.checking$Clust_selected <- ifelse(md.checking$Updated_order == input$Clusters_to_dis_PIE, input$Clusters_to_dis_PIE, "NS")
      md.checking$Clust_selected[is.na(md.checking$Clust_selected)] <- "NS"
      md.checking <- md.checking[order(md.checking$order), ]
      md.checking
      sc@meta.data <- md.checking
      Idents(object = sc) <- sc@meta.data$Clust_selected

      name.check.clust <- input$Clusters_to_dis_PIE
      min.pct.expression <- input$min_point_ # standard setting: 0.25
      min.logfc <- input$LogFC_ # 0.25 is standard

      markers.fm.list <- FindMarkers(sc, ident.1 = name.check.clust, min.pct = min.pct.expression, logfc.threshold = min.logfc, only.pos = TRUE)
      markers.fm.list2 <- subset(markers.fm.list, markers.fm.list$p_val_adj < input$pval.ex.filter)
      as.data.frame(markers.fm.list2)
    })

    output$compare.stat_Cluster_DT <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )

      # req(input$V_call_clust_sc,input$junction_clust_sc)
      compare.stat_Cluster()
    })

    output$downloaddf_FindMarker_Cluster <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Clusters_to_dis_PIE, "_Cluster_Stats_", x, ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(compare.stat_Cluster())
        write.csv(df, file, row.names = T)
      }
    )

    # Cluster dot plot -----

    output$Cluster_of_interest_DF <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )

      # req(input$V_call_clust_sc,input$junction_clust_sc)
      md <- sc@meta.data

      cluster <- clusTCR2_df()

      cluster$ID_Column <- cluster[, names(cluster) %in% input$Samp_col]
      cluster <- cluster[order(cluster$Updated_order), ]

      rownames(cluster) <- cluster$Cell_Index

      checking <- cluster[, names(cluster) %in% c("Updated_order", "Cell_Index")]
      md.checking <- merge(md, checking, by = "Cell_Index", all.x = T)
      md.checking <- md.checking[order(md.checking$order), ]
      rownames(md.checking) <- md.checking$Cell_Index

      md.checking$Clust_selected <- ifelse(md.checking$Updated_order == input$Clusters_to_dis_PIE, input$Clusters_to_dis_PIE, "NS")
      md.checking$Clust_selected[is.na(md.checking$Clust_selected)] <- "NS"
      md.checking <- md.checking[order(md.checking$order), ]
      md.checking
    })

    all_expression_plot_cluster <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )

      # req(input$V_call_clust_sc,input$junction_clust_sc)
      md <- sc@meta.data

      cluster <- clusTCR2_df()

      cluster$ID_Column <- cluster[, names(cluster) %in% input$Samp_col]
      cluster <- cluster[order(cluster$Updated_order), ]

      rownames(cluster) <- cluster$Cell_Index

      checking <- cluster[, names(cluster) %in% c("Updated_order", "Cell_Index")]
      md.checking <- merge(md, checking, by = "Cell_Index", all.x = T)
      md.checking <- md.checking[order(md.checking$order), ]
      rownames(md.checking) <- md.checking$Cell_Index

      md.checking$Clust_selected <- ifelse(md.checking$Updated_order == input$Clusters_to_dis_PIE, input$Clusters_to_dis_PIE, "NS")
      md.checking$Clust_selected[is.na(md.checking$Clust_selected)] <- "NS"
      md.checking <- md.checking[order(md.checking$order), ]
      md.checking

      sc@meta.data <- md.checking
      Idents(object = sc) <- sc@meta.data$Clust_selected

      if (input$restrict.dotpot.clust == F) {
        list.names <- rownames(compare.stat_Cluster())
      } else {
        list.names <- rownames(compare.stat_Cluster())
        list.names <- list.names[1:input$restrict.dotpot.num.clust]
      }

      size_legend <- input$Legend_size - 2

      DotPlot(sc, features = list.names) +
        RotatedAxis() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
          axis.title.x = element_blank(),
          legend.title = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.text = element_text(colour = "black", size = size_legend, family = input$font_type),
          legend.position = input$legend_position,
        ) +
        scale_colour_gradient2(low = input$low.dotplot.clust, mid = input$middle.dotplot.clust, high = input$high.dotplot.clust) +
        scale_x_discrete(labels = label_wrap(20)) +
        scale_y_discrete(labels = label_wrap(20))
    })

    output$all_expression_dotplot_cluster <- renderPlot({
      all_expression_plot_cluster()
    })

    output$downloadPlot_all_expression_dotplot_clust <- downloadHandler(
      filename = function() {
        paste(input$Clusters_to_dis_PIE, "_cluster_dotplot", "_", today(), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file,
            width = input$width_all_expression_dotplot_clust,
            height = input$height_all_expression_dotplot_clust, onefile = FALSE
        ) # open the pdf device
        df <- all_expression_plot_cluster()
        plot(df)
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_all_expression_dotplot_clust <- downloadHandler(
      filename = function() {
        paste(input$Clusters_to_dis_PIE, "_cluster_dotplot", "_", today(), ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_all_expression_dotplot_clust,
            height = input$height_png_all_expression_dotplot_clust,
            res = input$resolution_PNG_all_expression_dotplot_clust
        )
        df <- all_expression_plot_cluster()

        plot(df)
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )




    # Over representation analysis for Cluster  -----
    Over_rep_cluster_old <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )

      # req(input$V_call_clust_sc,input$junction_clust_sc)
      md <- sc@meta.data

      cluster <- clusTCR2_df()

      cluster$ID_Column <- cluster[, names(cluster) %in% input$Samp_col]
      cluster <- cluster[order(cluster$Updated_order), ]

      rownames(cluster) <- cluster$Cell_Index

      checking <- cluster[, names(cluster) %in% c("Updated_order", "Cell_Index")]
      md.checking <- merge(md, checking, by = "Cell_Index", all.x = T)
      md.checking <- md.checking[order(md.checking$order), ]
      rownames(md.checking) <- md.checking$Cell_Index

      md.checking$Clust_selected <- ifelse(md.checking$Updated_order == input$Clusters_to_dis_PIE, input$Clusters_to_dis_PIE, "NS")
      md.checking$Clust_selected[is.na(md.checking$Clust_selected)] <- "NS"
      md.checking <- md.checking[order(md.checking$order), ]
      md.checking

      sc@meta.data <- md.checking
      Idents(object = sc) <- sc@meta.data$Clust_selected
      sc@meta.data <- checking
      # require()

      geneSet <- read.csv(system.file("OverRep", "GeneSets.csv", package = "STEGO.R"), header = T)

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

      if (input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        geneSet$GeneSet <- gsub("-", ".", geneSet$GeneSet)
      }

      if (input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        require(stringr)
        geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
      }

      for (i in 1:dim(geneSet)[1]) {
        # listed GeneSet
        message(paste("GeneSet: ", i))
        Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet, ";")[i])
        names(Gene.set.testing) <- "V1"
        Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
        names(Gene.set.testing2) <- "V1"
        # message(paste(dim(Gene.set.testing2)[1],"GeneSet total"))
        # background genes
        background.overlap <- merge(Gene.set.testing2, background.genes.name, by = "V1")
        # message(paste(dim(background.overlap)[1],"in Background"))
        geneSet$background.geneset[i] <- length(background.overlap$V1)
        geneSet$background.geneset.name[i] <- as.character(paste(unlist(background.overlap[1]), collapse = ";"))
        # in sig gene list
        overlap <- merge(background.overlap, DEx.genes, by = "V1")
        # message(paste(dim(overlap)[1],"# Sig genes"))
        geneSet$in.geneset[i] <- length(overlap$V1)
        geneSet$in.geneset.name[i] <- as.character(paste(unlist(overlap[1]), collapse = ";"))
      }

      geneSet2 <- subset(geneSet, geneSet$in.geneset > 0)

      for (i in 1:dim(geneSet2)[1]) {
        tota.gene.set <- geneSet2$background.geneset[i] # genes that are identified in background
        tota.gene.set
        in.geneset <- geneSet2$in.geneset[i] # DEx in geneset

        background.genes
        not.in.total <- background.genes - tota.gene.set
        not.in.geneset.sig <- total.sig - in.geneset
        d <- data.frame(gene.in.interest = c(in.geneset, not.in.geneset.sig), gene.not.interest = c(tota.gene.set, not.in.total))
        row.names(d) <- c("In_category", "not_in_category")

        if (in.geneset > 0) {
          geneSet2$p.val[i] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
          geneSet2$lowerCI[i] <- unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
          geneSet2$upperCI[i] <- unlist(fisher.test(d)$conf.int)[2]
          geneSet2$OR[i] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1], 3)
        } else {
          geneSet2$p.value[i] <- "-"
          geneSet2$lowerCI[i] <- "-"
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
        need(
          nrow(geneSet2) > 0,
          "Upload Files"
        )
      )

      geneSet2 <- geneSet2[order(geneSet2$p.val, decreasing = F), ]
      geneSet2 <- subset(geneSet2, geneSet2$in.geneset >= input$in.geneset.cutoff_Clust)
      geneSet2 <- subset(geneSet2, geneSet2$p.val <= input$p.val_cutoff_Clust)
      geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
      geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
      geneSet2
    })

    output$Over_rep_Cluster_Tab <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1, 2, 5, 10, 20, 50, 100), pageLength = 20, scrollX = TRUE), {
      Over_rep_cluster()
    })

    output$downloadtb_over.rep.cluster <- downloadHandler(
      filename = function() {
        paste0(input$Clusters_to_dis_PIE, "_Cluster_", today(), "_over_rep.csv")
      },
      content = function(file) {
        df <- as.data.frame(Over_rep_cluster())
        write.csv(df, file, row.names = F)
      }
    )


    ##### Overlap -----
    #### upset plot -----

    Upset_plot <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )
      df <- sc@meta.data
      df <- as.data.frame(df)
      unique.df <- unique(df[, names(df) %in% c(input$Samp_col, input$V_gene_sc)])
      names(unique.df) <- c("group", "chain")
      unique.df <- unique.df[unique.df$group %in% input$ID_Column_factor, ]

      unique.df <- subset(unique.df, unique.df$chain != "NA")
      unique.df <- subset(unique.df, unique.df$group != "NA")
      unique.df$cloneCount <- 1
      mat <- acast(unique.df, chain ~ group, value.var = "cloneCount")
      mat[is.na(mat)] <- 0
      mat <- as.data.frame(mat)
      df.x <- make_comb_mat(mat)
      list.names <- as.character(input$ID_Column_factor)

      ht <- draw(UpSet(df.x,
                       pt_size = unit(5, "mm"),
                       lwd = 1,
                       row_names_gp = gpar(fontfamily = input$font_type, fontsize = 12),
                       column_names_gp = gpar(fontfamily = input$font_type, fontsize = 12),
                       top_annotation = upset_top_annotation(df.x,
                                                             add_numbers = T,
                                                             numbers_gp = gpar(fontfamily = input$font_type, fontsize = 12),
                                                             annotation_name_gp = gpar(fontfamily = input$font_type, fontsize = 12),
                                                             gp = gpar(fill = "black"),
                       ),
                       right_annotation = upset_right_annotation(df.x,
                                                                 add_numbers = T,
                                                                 numbers_gp = gpar(fontfamily = input$font_type, fontsize = 12),
                                                                 annotation_name_gp = gpar(fontfamily = input$font_type, fontsize = 12),
                                                                 gp = gpar(fill = "black"),
                       ),
                       set_order = list.names
      ), padding = unit(c(20, 20, 20, 20), "mm"))
      ht
    })
    output$Upset_plot_overlap <- renderPlot({
      Upset_plot()
    })

    output$downloadPlot_Upset_plot_overlap <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Samp_col, "_", input$V_gene_sc, "_overlap_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_Upset_plot_overlap, height = input$height_Upset_plot_overlap, onefile = FALSE) # open the pdf device
        plot(Upset_plot())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_Upset_plot_overlap <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Samp_col, "_", input$V_gene_sc, "_overlap_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Upset_plot_overlap, height = input$height_png_Upset_plot_overlap, res = input$resolution_PNG_Upset_plot_overlap)
        plot(Upset_plot())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    # upset plot table -----
    Upset_plot_overlap <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )

      df <- sc@meta.data
      df <- as.data.frame(df)
      unique.df <- unique(df[, names(df) %in% c(input$Samp_col, input$V_gene_sc)])
      names(unique.df) <- c("group", "chain")
      unique.df <- subset(unique.df, unique.df$chain != "NA")
      unique.df <- subset(unique.df, unique.df$group != "NA")
      unique.df$cloneCount <- 1
      mat <- acast(unique.df, chain ~ group, value.var = "cloneCount")
      mat[is.na(mat)] <- 0
      Count_data <- as.data.frame(rowSums(mat))
      names(Count_data) <- "V1"
      unique.df <- (df[, names(df) %in% c(input$Samp_col, input$V_gene_sc)])
      names(unique.df) <- c("group", "chain")
      unique.df <- subset(unique.df, unique.df$chain != "NA")
      unique.df <- subset(unique.df, unique.df$group != "NA")
      unique.df$cloneCount <- 1
      mat <- acast(unique.df, chain ~ group, value.var = "cloneCount")
      mat[is.na(mat)] <- 0
      sum_data <- as.data.frame(rowSums(mat))
      names(sum_data) <- "V1"
      mat <- as.data.frame(mat)
      mat$TotalSamps <- Count_data$V1
      mat$CloneTotal <- sum_data$V1
      mat
    })

    output$Upset_plot_overlap_Tb <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 5, scrollX = TRUE), {
      Upset_plot_overlap()
    })

    output$downloaddf_Upset_plot_overlap_Tb <- downloadHandler(
      filename = function() {
        paste(input$name_BD, "_Overlap_", gsub("-", ".", Sys.Date()), ".tsv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(Upset_plot_overlap())
        write.csv(df, file)
      }
    )


    # line graph from upset plot  ----

    check_sep <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )

      clones <- Upset_plot_overlap()
      req(clones)

      if (input$comparison_operator == "==") {
        original_data <- subset(clones, TotalSamps == input$cutoff_upset)
      } else if (input$comparison_operator == ">=") {
        original_data <- subset(clones, TotalSamps >= input$cutoff_upset)
      }

      # Get the group names from the column names
      group_names <- unique(gsub(input$separator_input,"", colnames(original_data)))
      group_names <- group_names[!(group_names %in% c("background", "TotalSamps", "CloneTotal"))]
      print(group_names)
      group_names

    })

    select_top_five <- reactive({

      clones <- Upset_plot_overlap()
      req(clones)

      if (input$comparison_operator == "==") {
        original_data <- subset(clones, TotalSamps == input$cutoff_upset)
      } else if (input$comparison_operator == ">=") {
        original_data <- subset(clones, TotalSamps >= input$cutoff_upset)
      }
      print(dim(original_data))
      # Get the group names from the column names
      group_names <- unique(gsub(input$separator_input,"", colnames(original_data)))  # Assuming first two columns are not part of the groups

      len <- length(unlist(group_names)) - 2

      # Create an empty list to store top 5 data for each group
      top_5_data_list <- list()

      if(input$is_a_time_series) {
        for (i in 1:len) {
          # Subset the data for the current group
          group_data <- original_data[,names(original_data) %in% c(colnames(original_data)[grepl(group_names[i], colnames(original_data))])]

          # Check if group_data is empty or not a data frame
          if (!is.data.frame(group_data) || nrow(group_data) == 0) {
            cat("Group data for", group_names[i], "is empty or not a data frame. Skipping.\n")
            next  # Skip to the next iteration if group_data is empty or not a data frame
          }

          len_group_data <- dim(group_data)[2]
          print(len_group_data)
          # Calculate row sums
          group_data$CloneTotal <- rowSums(group_data)
          group_data <- subset(group_data,group_data$CloneTotal >= input$Total_count_Cutoff )

          if (dim(group_data)[2]==0) {
            cat(group_names[i], "has not clones. Skipping.\n")
            next  # Skip to the next iteration if group_data is empty or not a data frame
          }

          # Order by CloneTotal in descending order
          group_data <- group_data[order(-group_data$CloneTotal), ]

          top_5_group <- group_data %>%
            slice_max(CloneTotal, n = input$max_number_lines_to)

          # Remove the CloneTotal column
          top_5_group <- top_5_group[, !grepl("^CloneTotal", colnames(top_5_group))]
          print(top_5_group)
          # Store the top 5 data for the current group in the list object
          top_5_data_list[[group_names[i]]] <- top_5_group

        }

      } else {
        print(names(original_data))
        # print(check_sep())
        group_data <- original_data[,!names(original_data) %in% c("background", "TotalSamps", "CloneTotal")]
        print(head(group_data))
        group_data$CloneTotal <- rowSums(group_data)
        print(head(group_data))
        group_data <- subset(group_data,group_data$CloneTotal >= input$Total_count_Cutoff)
        group_data <- group_data[order(-group_data$CloneTotal), ]
        top_5_group <- group_data %>%
          slice_max(CloneTotal, n = input$max_number_lines_to)
        top_5_group <- top_5_group[, !grepl("^CloneTotal", colnames(top_5_group))]
        print(top_5_group)
        top_5_data_list <- top_5_group
        print(top_5_data_list)
      }

      top_5_data_list

    })

    observeEvent(input$load_samp_name_list,{

      # Get the group names from the column names
      group_names <- check_sep()
      req(group_names)

      if (length(group_names)>0) {
        updateSelectInput(
          session,
          "Group_for_line_graph",
          choices = group_names

        )
      } else {
        updateSelectInput(
          session,
          "Group_for_line_graph",
          choices = ""

        )
      }

    })

    output$Line_graph_table <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1, 2, 5, 10, 20, 50, 100), pageLength = 20, scrollX = TRUE), {

      list.df <- select_top_five()

      if(input$is_a_time_series) {
        df <- as.data.frame(list.df[[input$Group_for_line_graph]])
        df } else {

          select_top_five()
        }
    })

    Line_graph_for_tracing <- reactive({

      top_5_data_list <- select_top_five()
      print(top_5_data_list)
      if(input$is_a_time_series) {

        print(top_5_data_list)
        req(top_5_data_list)
        # Remove empty data frames from top_5_data_list
        top_5_data_list <- top_5_data_list[sapply(top_5_data_list, function(x) nrow(x) > 0)]

      } else {
        top_5_data_list <- select_top_five()

      }
      # Find the maximum count value across all datasets
      max_count <- max(sapply(top_5_data_list, function(df) max(df)))
      print(max_count)
      # Round the maximum count value to the nearest 5 and then add 5
      max_count <- ceiling(max_count / 5) * 5

      if(input$is_a_time_series) {
        # Create an empty list to store plots
        plot_list <- list()

        # Iterate over each year in top_5_data_list
        for (year in names(top_5_data_list)) {
          # Transpose the data frame and convert to data.frame
          top_5_data <- top_5_data_list[[year]]

          top_5_transposed <- as.data.frame(t(top_5_data), stringsAsFactors = FALSE)

          # Add Year column
          top_5_transposed$ID <- year

          # Convert row names into a regular column
          top_5_transposed$Sample_ID <- rownames(top_5_transposed)
          # Split Sample_ID into separate columns for Group and Time

          # as.data.frame(do.call(rbind, strsplit(as.character(top_5_transposed$Sample_ID), "[.]")))

          sample_parts <- strsplit(top_5_transposed$Sample_ID, input$separator_input2)
          print(sample_parts)
          #   # Loop over each part and assign them to new columns

          for (i in 1:length(names(top_5_data_list))) {  # Assuming there are 3 parts after splitting
            # Create new column names V1, V2, V3
            new_col_name <- paste0("V", i)

            # Extract the ith part from each split
            top_5_transposed[[new_col_name]] <- sapply(sample_parts, function(x) ifelse(length(x) >= i, x[i], NA))
          }

          print(head(top_5_transposed))
          all_names <- names(top_5_transposed)
          print(all_names)
          v_gene_names <- all_names[grep("TRAV|TRBV|TRGV|TRDV", all_names)]
          print(v_gene_names)

          # Reshape the data into long format
          data_long <- pivot_longer(top_5_transposed,
                                    cols = v_gene_names,   # Exclude the Sample_ID, Group, and Time columns
                                    names_to = "VDJ",  # New column name for time points
                                    values_to = "Count")     # New column name for values
          # data_long
          print(data_long)
          # Replace underscores with spaces in the VDJ variable
          data_long$VDJ <- gsub("_", " ", data_long$VDJ)
          data_long$VDJ <- gsub(" & ", " ", data_long$VDJ)
          #   # Wrap the Time_Point variable based on spaces
          data_long$VDJ <- str_wrap(data_long$VDJ, width = 10)  # Adjust width as needed
          #
          #   # Determine unique levels of VDJ
          unique_vdj <- as.data.frame(unique(data_long$VDJ))
          names(unique_vdj) <- "unique_vdj"

          set.seed(200)
          unique_vdj$shape <- sample(1:25, nrow(unique_vdj))
          data_long

          # Create the plot

          if(input$number_of_conditions == 2) {

            req(length(data_long$V2>0))

            p <- ggplot(data_long, aes(x = V2, y = Count, color = VDJ, shape = VDJ)) +
              geom_point(size = 7) +  # Increased point size
              geom_line(aes(group = paste(VDJ, V1)), linewidth = 1.25) +
              scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +
              scale_shape_manual(values = unique_vdj$shape) +  # Use default shapes
              labs(x = "", y = "", title = "", color = year, shape = year) +
              theme_minimal() +
              theme(legend.title = element_text(face = "bold", size = 16, family = input$font_type),
                    legend.text = element_text(size = input$Legend_size, family = input$font_type),
                    axis.text = element_text(size = 16, family = input$font_type),
                    axis.title = element_blank(),
                    plot.title = element_blank()
              ) +
              guides(color = guide_legend(
                title.theme = element_text(margin = margin(b = 0)),  # Increase margin between title and items
                label.theme = element_text(margin = margin(t = 15)),
                override.aes = list(size = 7)# Increase margin between items
              )) +
              ylim(0, max_count)  # Set y-axis limits

            # Store the plot in the list
            plot_list[[year]] <- p

          } else if (input$number_of_conditions == 3) {

            req(length(data_long$V2>0), length(data_long$V3>0))

            p <- ggplot(data_long, aes(x = V2, y = Count, color = VDJ, shape = V3)) +
              geom_point(size = 7) +  # Increased point size
              geom_line(aes(group = paste(VDJ, V1)), linewidth = 1.25) +
              scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +
              scale_shape_manual(values = unique_vdj$shape) +  # Use default shapes
              labs(x = "", y = "", title = "", color = year, shape = input$shape_legend_name) +
              theme_minimal() +
              theme(legend.title = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
                    legend.text = element_text(size = input$Legend_size, family = input$font_type),
                    axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
                    axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 0),
                    axis.title = element_blank(),
                    plot.title = element_blank()
              ) +
              guides(colour = guide_legend(order = 1,
                                           title.theme = element_text(margin = margin(b = 0)),  # Increase margin between title and items
                                           label.theme = element_text(margin = margin(t = 15)),

              ),
              shape = guide_legend(order = 1,
                                   title.theme = element_text(margin = margin(b = 0)),  # Increase margin between title and items
                                   label.theme = element_text(margin = margin(t = 15)),
              )) +
              ylim(0, max_count)  # Set y-axis limits

            # Store the plot in the list
            plot_list[[year]] <- p
          }
        }
      } else {
        top_5_transposed <- as.data.frame(t(top_5_data_list), stringsAsFactors = FALSE)
        print(top_5_transposed)

        # Convert row names into a regular column
        top_5_transposed$Sample_ID <- rownames(top_5_transposed)
        # Split Sample_ID into separate columns for Group and Time

        print(head(top_5_transposed))
        all_names <- names(top_5_transposed)
        print(all_names)
        # Find the maximum count value across all datasets
        max_count <- max(sapply(top_5_data_list, function(df) max(df)))
        print(max_count)
        # Round the maximum count value to the nearest 5 and then add 5
        max_count <- ceiling(max_count / 5) * 5
        top_5_transposed <- as.data.frame(t(top_5_data_list), stringsAsFactors = FALSE)
        # Convert row names into a regular column
        top_5_transposed$Sample_ID <- rownames(top_5_transposed)
        top_5_transposed
        # Split Sample_ID into separate columns for Group and Time
        all_names <- names(top_5_transposed)
        print(all_names)
        v_gene_names <- all_names[grep("TRAV|TRBV|TRGV|TRDV", all_names)]
        v_gene_names

        # Reshape the data into long format
        data_long <- pivot_longer(top_5_transposed,
                                  cols = v_gene_names,   # Exclude the Sample_ID, Group, and Time columns
                                  names_to = "VDJ",  # New column name for time points
                                  values_to = "Count")     # New column name for values
        #   # data_long

        #   # Replace underscores with spaces in the VDJ variable
        data_long$VDJ <- gsub("_", " ", data_long$VDJ)
        data_long$VDJ <- gsub(" & ", " ", data_long$VDJ)
        #   #   # Wrap the Time_Point variable based on spaces
        data_long$VDJ <- str_wrap(data_long$VDJ, width = 10)  # Adjust width as needed

        #   #   # Determine unique levels of VDJ
        unique_vdj <- as.data.frame(unique(data_long$VDJ))
        names(unique_vdj) <- "unique_vdj"
        unique_vdj
        set.seed(200)
        unique_vdj$shape <- sample(1:25, nrow(unique_vdj))
        print(data_long)
        #
        #
        plot_list <- ggplot(data_long, aes(x = Sample_ID, y = Count, color = VDJ, shape = VDJ)) +
          geom_point(size = 7) +  # Increased point size
          geom_line(aes(group = paste(VDJ)), linewidth = 1.25) +
          scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +
          scale_shape_manual(values = unique_vdj$shape) +  # Use default shapes
          labs(x = "TCR counts", y = "", title = "", color = "VDJ", shape = "VDJ") +
          theme_minimal() +
          theme(
            legend.title = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
            legend.text = element_text(size = input$Legend_size, family = input$font_type),
            axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
            axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 0),
            axis.title.y  = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
            axis.title.x  = element_blank()
          ) +
          guides(color = guide_legend(
            title.theme = element_text(margin = margin(b = 0)),  # Increase margin between title and items
            label.theme = element_text(margin = margin(t = 15)),
            override.aes = list(size = 7)# Increase margin between items
          )) +
          ylim(0, max_count)  # Set y-axis limits
        # }

      }
      plot_list
    })


    output$line_graph_output <- renderPlot({
      plot_list <- Line_graph_for_tracing()
      req(plot_list)


      if(input$is_a_time_series) {
        print(names(plot_list))
        plot_list[[input$Group_for_line_graph]]

      } else {
        print(plot_list)

        # p <- plot_list[[1]]
        plot_list
      }
    })



    output$line_graph_all_output <- renderPlot({
      require(cowplot)
      require(grid)
      require(gridExtra)
      plot_list <- Line_graph_for_tracing()
      req(plot_list)
      combined_plots <- plot_grid(plotlist = plot_list, nrow = input$wrap_row)
      y.grob <- textGrob("TCR counts",
                         gp=gpar(fontface = "bold", col="black", family = input$font_type, fontsize=30), rot=90)

      x.grob <- textGrob("Time",
                         gp=gpar(fontface = "bold", col="black", family = input$font_type, fontsize=30))

      grid.arrange(arrangeGrob(combined_plots, left = y.grob, bottom = x.grob))
    })


    # overlap table with UMAP and expression -----
    overlap_table <- reactive({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )
      reduction <- (sc@reductions$umap)
      UMAP <- as.data.frame(reduction@cell.embeddings)
      names(UMAP)[1:2] <- toupper(names(UMAP)[1:2])
      UMAP$Cell_Index <- rownames(UMAP)
      meta.data <- as.data.frame(sc@meta.data)
      umap.meta <- merge(UMAP, meta.data, by = "Cell_Index")

      umap.meta$chain <- umap.meta[,names(umap.meta) %in% input$V_gene_sc]
      umap.meta$ID_Column <- umap.meta[,names(umap.meta) %in% input$Samp_col]
      umap.meta$Selected_function <- umap.meta[,names(umap.meta) %in% input$Colour_By_this]
      umap.meta$Selected_group <- umap.meta[,names(umap.meta) %in% input$Split_group_by_]

      # names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
      # names(umap.meta)[names(umap.meta) %in% input$Colour_By_this] <- "Selected_function"
      # names(umap.meta)[names(umap.meta) %in% input$Split_group_by_] <- "Selected_group"
      # sc <- merge(umap.meta,TCR_Expanded(),by=c("v_gene_selected","ID_Column"),all.x=T)

      Upset_plot_overlap <- Upset_plot_overlap()
      Upset_plot_overlap_top <- subset(Upset_plot_overlap, Upset_plot_overlap$sum > 1)
      Upset_plot_overlap_top$chain <- rownames(Upset_plot_overlap_top)
      umap.meta.overlap <- merge(Upset_plot_overlap_top, umap.meta, by = c("chain",map.meta$ID_Column))
      print(names(umap.meta.overlap))
      umap.meta.overlap
    })

    # overlap UMAP plot
    create_UMAP_overlap <- reactive({
      umap.meta.overlap <- overlap_table()
      validate(
        need(
          nrow(umap.meta.overlap) > 0,
          error_message_val_UMAP
        )
      )

      ggplot(umap.meta.overlap, aes(x = UMAP_1, UMAP_2, colour = Selected_function)) +
        geom_point() +
        scale_color_manual(values = rainbow(length(unique(umap.meta.overlap$Selected_function))), na.value = input$NA_col_analysis) +
        theme_bw() +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.title.x = element_text(colour = "black", angle = 0, vjust = .5, face = "plain", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 0),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
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

      df3.meta3 <- as.data.frame(table(top_BD_cluster$Selected_group, top_BD_cluster$Selected_function))
      total.condition <- as.data.frame(ddply(df3.meta3, "Var1", numcolwise(sum)))
      emtpy <- matrix(nrow = dim(df3.meta3)[1], ncol = dim(total.condition)[1])
      for (i in 1:dim(df3.meta3)[1]) {
        emtpy[i, ] <- ifelse(df3.meta3$Var1[i] == total.condition$Var1[1:dim(total.condition)[1]],
                             total.condition[total.condition$Var1 == total.condition$Var1[1:dim(total.condition)[1]], 2], F
        )
      }

      df3.meta3$n <- df3.meta3$Freq / rowSums(emtpy)


      ggplot(df3.meta3, aes(x = "", y = n, fill = Var2, group = Var1)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        theme_void(20) +
        facet_wrap(~Var1, nrow = input$wrap_row) +
        theme(
          legend.key.size = unit(1, "cm")
        ) +
        # scale_fill_manual(values = df.col, na.value = input$NA_col_analysis) +
        theme(
          strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
          legend.text = element_text(size = input$Legend_size, family = input$font_type),
          legend.title = element_blank()
        )
    })

    output$Classification_Overlap_pie <- renderPlot({
      Overlap_Pie_chart_Class()
    })

    # Over representation analysis for Cluster  -----

    Over_rep_Overlap <- reactive({
      sc <- data_sc_pro()

      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )

      df <- sc@meta.data
      # require()

      geneSet <- read.csv(system.file("OverRep", "GeneSets.csv", package = "STEGO.R"), header = T)

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

      if (input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        geneSet$GeneSet <- gsub("-", ".", geneSet$GeneSet)
      }

      if (input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        require(stringr)
        geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
      }

      for (i in 1:dim(geneSet)[1]) {
        # listed GeneSet
        message(paste("GeneSet: ", i))
        Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet, ";")[i])
        names(Gene.set.testing) <- "V1"
        Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
        names(Gene.set.testing2) <- "V1"
        # message(paste(dim(Gene.set.testing2)[1],"GeneSet total"))
        # background genes
        background.overlap <- merge(Gene.set.testing2, background.genes.name, by = "V1")
        # message(paste(dim(background.overlap)[1],"in Background"))
        geneSet$background.geneset[i] <- length(background.overlap$V1)
        geneSet$background.geneset.name[i] <- as.character(paste(unlist(background.overlap[1]), collapse = ";"))
        # in sig gene list
        overlap <- merge(background.overlap, DEx.genes, by = "V1")
        # message(paste(dim(overlap)[1],"# Sig genes"))
        geneSet$in.geneset[i] <- length(overlap$V1)
        geneSet$in.geneset.name[i] <- as.character(paste(unlist(overlap[1]), collapse = ";"))
      }

      geneSet2 <- subset(geneSet, geneSet$in.geneset > 0)

      for (i in 1:dim(geneSet2)[1]) {
        tota.gene.set <- geneSet2$background.geneset[i] # genes that are identified in background
        tota.gene.set
        in.geneset <- geneSet2$in.geneset[i] # DEx in geneset

        background.genes
        not.in.total <- background.genes - tota.gene.set
        not.in.geneset.sig <- total.sig - in.geneset
        d <- data.frame(gene.in.interest = c(in.geneset, not.in.geneset.sig), gene.not.interest = c(tota.gene.set, not.in.total))
        row.names(d) <- c("In_category", "not_in_category")

        if (in.geneset > 0) {
          geneSet2$p.val[i] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
          geneSet2$lowerCI[i] <- unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
          geneSet2$upperCI[i] <- unlist(fisher.test(d)$conf.int)[2]
          geneSet2$OR[i] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1], 3)
        } else {
          geneSet2$p.value[i] <- "-"
          geneSet2$lowerCI[i] <- "-"
          geneSet2$upperCI[i] <- "-"
          geneSet2$OR[i] <- "-"
        }
        # message(print(d))
        # message(print(round(prop.table(d),3)))
      }
      geneSet2 <- geneSet2[order(geneSet2$p.val, decreasing = F), ]
      geneSet2 <- geneSet2[order(geneSet2$p.val, decreasing = F), ]
      geneSet2 <- subset(geneSet2, geneSet2$in.geneset >= input$in.geneset.cutoff_Exp)
      geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
      geneSet2$pval.bonferroni.adj <- p.adjust(geneSet2$p.val, method = "bonferroni")
      geneSet2 <- subset(geneSet2, geneSet2$p.val <= input$p.val_cutoff_Exp)
      # geneSet2 <- subset(geneSet2,geneSet2$pval.BH.adj<=input$adjust_cutoff_top)

      geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
      geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
      name.list <- c("Geneset_ID", "p.val", "FDR", "Bonferroni", "OR", "lowerCI", "upperCI", "in.geneset.name", "in.geneset", "background.geneset", "total.sig", "background.genes", "background.geneset.name")
      geneSet2 <- geneSet2 %>%
        select(all_of(name.list), everything())
    })

    output$Over_rep_overlap_Tab <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1, 2, 5, 10, 20, 50, 100), pageLength = 20, scrollX = TRUE), {
      Over_rep_Overlap()
    })

    output$downloadtb_over.rep.overlap <- downloadHandler(
      filename = function() {
        paste0("Overlap", "_", today(), "_over_rep.csv")
      },
      content = function(file) {
        df <- as.data.frame(Over_rep_Overlap())
        write.csv(df, file, row.names = F)
      }
    )

    # marker analysis -----

    ### single marker feature plot -----

    MainTcell_counts <- reactive({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )
      md <- sc@meta.data

      df1 <- as.data.frame(sc@assays$RNA$counts[rownames(sc@assays$RNA$counts) %in% c(input$Var_to_col_marker, input$Var_to_col_marker2, input$Var_to_col_marker3, "CD4", "CD8A", "Cd4", "Cd8", "CD3E", "Cd3e"), ])
      req(df1)
      names(df1) <- colnames(sc@assays$RNA$counts)
      df1 <- as.data.frame(t(df1))
      df1$Cell_Index <- rownames(df1)
      df1
    })

    MainTcell_counts_names <- reactive({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )
      req(sc)

      df <- rownames(sc@assays$RNA$scale.data)[rowSums(sc@assays$RNA$scale.data) != 0]

      df <- as.data.frame(df)
      names(df) <- "V1"
      df
    })

    observeEvent(input$load_marker_genes_sig, {
      sc <- (MainTcell_counts_names())
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
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
        need(
          nrow(MainTcell) > 0,
          "Upload File"
        )
      )

      MainTcell$Cell_Index <- rownames(MainTcell)
      MainTcell
    })

    meta_data_for_features <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )

      md <- sc@meta.data
      md <- merge(md, selected_gene(), by = "Cell_Index")
      md$log2.count <- log2(as.numeric(md[, names(md) %in% input$Var_to_col_marker]))
      md$log2.count[is.na(md$log2.count)] <- "-Inf"
      md
      md <- subset(md, md$UMAP_1 > input$Filter_lower_UMAP1_marker)
      md <- subset(md, md$UMAP_1 < input$Filter_lower_UMAP1_marker2)

      md <- subset(md, md$UMAP_2 > input$Filter_lower_UMAP2_marker)
      md <- subset(md, md$UMAP_2 < input$Filter_lower_UMAP2_marker2)

      subset(md, md$v_gene_AG != "NA")
    })

    # subset based on UMAP location, subset bashed on norm value for the marker of interest threshold as well. -----

    MainTcell_scale <- reactive({
      sc <- data_sc_pro()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )
      md <- sc@meta.data
      df1 <- as.data.frame(sc@assays$RNA$scale.data[rownames(sc@assays$RNA$scale.data) %in% c(input$Var_to_col_marker, input$Var_to_col_marker2, input$Var_to_col_marker3, "CD4", "CD8A", "Cd4", "Cd8", "CD3E", "Cd3e", "JCHAIN", "Jchain", "CD19", "Cd19", "MS4A1", "Ms4a1"), ])

      names(df1) <- colnames(sc@assays$RNA$scale.data)
      df1 <- as.data.frame(t(df1))
      df1$Cell_Index <- rownames(df1)
      df1
    })

    selected_scale <- reactive({
      MainTcell <- MainTcell_scale()
      validate(
        need(
          nrow(MainTcell) > 0,
          "Upload File"
        )
      )

      MainTcell$Cell_Index <- rownames(MainTcell)
      MainTcell
    })

    meta_data_for_features_scale <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )

      md <- sc@meta.data
      md <- merge(md, selected_scale(), by = "Cell_Index")
      md$scale <- as.numeric(md[, names(md) %in% input$Var_to_col_marker])

      md <- subset(md, md$UMAP_1 > input$Filter_lower_UMAP1_marker)
      md <- subset(md, md$UMAP_1 < input$Filter_lower_UMAP1_marker2)

      md <- subset(md, md$UMAP_2 > input$Filter_lower_UMAP2_marker)
      md <- subset(md, md$UMAP_2 < input$Filter_lower_UMAP2_marker2)
      subset(md, md$v_gene_AG != "NA")
    })


    # umap scale -------
    marker_selected_UMAP <- reactive({
      umap.meta <- meta_data_for_features_scale()
      validate(
        need(
          nrow(umap.meta) > 0,
          "Upload File"
        )
      )

      umap.meta <- meta_data_for_features_scale()
      umap.meta_pos <- subset(umap.meta, umap.meta$scale != "-Inf")

      umap.meta$scale <- ifelse(umap.meta$scale == "-Inf", NA, umap.meta$scale)
      umap.meta$scale <- as.numeric(umap.meta$scale)
      umap.meta_pos$scale <- as.numeric(umap.meta_pos$scale)

      ggplot(umap.meta, aes(x = UMAP_1, y = UMAP_2, color = scale)) +
        geom_point(size = 1) +
        geom_point(data = umap.meta_pos, aes(x = UMAP_1, y = UMAP_2, color = scale), size = 1) +
        # scale_color_continuous(type = "viridis",na.value="grey75") +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
        ) +
        scale_color_distiller(direction = 1, palette = input$col_marker_scale, na.value = "grey85", limits = c(min(0), max(input$max_scale))) +
        theme(
          plot.title = element_text(colour = "black", family = input$font_type, size = 24, face = "bold", vjust = .5),
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.title.x = element_text(colour = "black", angle = 0, vjust = .5, face = "plain", family = input$font_type, size = input$title.text.sizer2),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        ) +
        labs(title = input$Var_to_col_marker)
    })

    output$marker_selected_UMAP_plot <- renderPlot({
      umap.meta <- MainTcell_counts_names()
      validate(
        need(
          nrow(umap.meta) > 0,
          "Upload File"
        )
      )
      marker_selected_UMAP()
    })


    output$downloadPlot_marker_selected_UMAP_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Var_to_col_marker, "_marker_selected_UMAP_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_marker_selected_UMAP_plot, height = input$height_marker_selected_UMAP_plot, onefile = FALSE) # open the pdf device

        plot(marker_selected_UMAP())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_marker_selected_UMAP_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Var_to_col_marker, "_marker_selected_UMAP_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_marker_selected_UMAP_plot, height = input$height_png_marker_selected_UMAP_plot, res = input$resolution_PNG_marker_selected_UMAP_plot)
        plot(marker_selected_UMAP())
        dev.off()
      }, contentType = "application/png"
    )

    # scale violin plot ----
    marker_selected_Violine.ridge <- reactive({
      # umap.meta <- meta_data_for_features_scale()
      umap.meta <- meta_data_for_features()
      validate(
        need(
          nrow(umap.meta) > 0,
          "Upload File"
        )
      )

      umap.meta <- meta_data_for_features_scale()

      if (input$select_plot_vio.ridge == "Violin") {
        ggplot(umap.meta, aes(y = scale, x = Sample_Name)) +
          geom_violin() +
          geom_jitter() +
          theme(
            legend.position = "none",
          ) +
          theme_bw() +
          geom_hline(yintercept = input$cutoff_marker_gt, color = "orange", size = 3)
        # geom_vline(xintercept=neg, linetype="dashed", color = "darkorange")
      } else {
        ggplot(umap.meta, aes(x = scale, y = Sample_Name)) +
          geom_density_ridges() +
          theme_ridges() +
          geom_vline(xintercept = input$cutoff_marker_gt, color = "orange", size = 3)
      }
    })


    output$marker_selected_VioRidge_plot <- renderPlot({
      umap.meta <- MainTcell_counts_names()
      validate(
        need(
          nrow(umap.meta) > 0,
          "Upload File"
        )
      )
      marker_selected_Violine.ridge()
    })


    filtered_positive_marker_TCRsum <- reactive({
      # umap.meta <- meta_data_for_features_scale()
      umap.meta <- meta_data_for_features()
      validate(
        need(
          nrow(umap.meta) > 0,
          "Upload File"
        )
      )

      umap.meta$log2.count <- ifelse(umap.meta$log2.count == "-Inf", -2, umap.meta$log2.count)
      umap.meta$log2.count <- as.numeric(umap.meta$log2.count)

      umap.meta_marker_pos <- subset(umap.meta, umap.meta$log2.count > input$cutoff_marker_gt)
      umap.meta_marker_pos$cloneCount <- 1
      umap.meta_marker_pos_TCR <- as.data.frame(umap.meta_marker_pos[, names(umap.meta_marker_pos) %in% c(input$V_gene_sc, "cloneCount")])
      file.names <- input$V_gene_sc

      df_unique_sum <- ddply(umap.meta_marker_pos_TCR, file.names, numcolwise(sum))
      df_unique_sum$type <- "Pos"
      df_unique_sum[order(df_unique_sum$cloneCount, decreasing = T), ]
    })

    filtered_negative_marker_TCRsum <- reactive({
      umap.meta <- meta_data_for_features()
      validate(
        need(
          nrow(umap.meta) > 0,
          "Upload File"
        )
      )

      umap.meta$log2.count <- ifelse(umap.meta$log2.count == "-Inf", -2, umap.meta$log2.count)
      umap.meta$log2.count <- as.numeric(umap.meta$log2.count)

      umap.meta_marker_neg <- subset(umap.meta, umap.meta$log2.count < input$cutoff_marker_gt)
      umap.meta_marker_neg$cloneCount <- 1
      umap.meta_marker_neg_TCR <- as.data.frame(umap.meta_marker_neg[, names(umap.meta_marker_neg) %in% c(input$V_gene_sc, "cloneCount")])
      file.names <- input$V_gene_sc

      df_unique_sum <- ddply(umap.meta_marker_neg_TCR, file.names, numcolwise(sum))
      df_unique_sum$type <- "Neg"
      df_unique_sum[order(df_unique_sum$cloneCount, decreasing = T), ]
    })

    output$TCR_marker_positive_count <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1, 2, 5, 10, 20, 50, 100), pageLength = 1, scrollX = TRUE), {
      umap.meta <- filtered_positive_marker_TCRsum()
      umap.meta
    })

    output$TCR_marker_neg_count <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1, 2, 5, 10, 20, 50, 100), pageLength = 1, scrollX = TRUE), {
      umap.meta <- filtered_negative_marker_TCRsum()
      umap.meta
    })

    # histogram of counts marker and TCR -----

    merged_marker_hist <- reactive({
      a <- filtered_negative_marker_TCRsum()
      a$expand <- ifelse(a$cloneCount >= 2, "Ex", "NEx")
      names(a)[2:4] <- paste(names(a)[2:4], "_neg", sep = "")

      b <- filtered_positive_marker_TCRsum()
      b$expand <- ifelse(b$cloneCount >= 2, "Ex", "NEx")
      names(b)[2:4] <- paste(names(b)[2:4], "_pos", sep = "")

      by_ab <- names(a)[1]
      ab <- merge(b, a, by = by_ab, all = T)



      ab[is.na(ab)] <- 0

      ab$expand_pos <- ifelse(ab$expand_pos != 0, ab$expand_pos, "Abs")
      ab$expand_neg <- ifelse(ab$expand_neg != 0, ab$expand_neg, "Abs")
      ab <- ab[, !names(ab) %in% c("type_neg", "type_pos")]
      ab$diff <- ab$cloneCount_pos - ab$cloneCount_neg
      ab$ratio <- ab$cloneCount_pos / ab$cloneCount_neg
      ab[order(ab$diff, decreasing = T), ]
    })

    output$merged_marker_hist_table <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1, 2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      umap.meta <- merged_marker_hist()
      umap.meta
    })

    marker_histogram <- reactive({
      a <- filtered_negative_marker_TCRsum()
      b <- filtered_positive_marker_TCRsum()
      ab <- rbind(a, b)

      # ab <- subset(ab,ab$cloneCount>2)

      ab$log10_ab <- log10(ab$cloneCount)
      ggplot(ab, aes(x = cloneCount, color = type)) +
        geom_histogram(binwidth = 0.1, fill = "white", alpha = 0.5, position = "identity") +
        theme_classic() +
        # geom_density(alpha=.2) +
        theme(legend.position = "top") +
        facet_wrap(~type)
    })

    output$marker_selected_histogram_plot <- renderPlot({
      umap.meta <- MainTcell_counts_names()
      validate(
        need(
          nrow(umap.meta) > 0,
          "Upload File"
        )
      )
      marker_histogram()
    })

    stats_marker_one <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )

      md <- sc@meta.data
      md <- merge(md, selected_scale(), by = "Cell_Index", all.x = T, sort = F)
      rownames(md) <- md$Cell_Index
      sc@meta.data <- md
      sc@meta.data$scale <- (as.numeric(md[, names(md) %in% input$Var_to_col_marker]))

      sc <- subset(sc, UMAP_1 > input$Filter_lower_UMAP1_marker)
      sc <- subset(sc, UMAP_1 < input$Filter_lower_UMAP1_marker2)

      sc <- subset(sc, UMAP_2 > input$Filter_lower_UMAP2_marker)
      sc <- subset(sc, UMAP_2 < input$Filter_lower_UMAP2_marker2)
      # sc <- subset(sc,v_gene_AG != "NA")

      sc@meta.data$pos_neg <- ifelse(sc@meta.data$scale > input$cutoff_marker_gt, "pos", "neg")
      sc
    })

    output$marker_selected_tab <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1, 2, 5, 10, 20, 50, 100), pageLength = 20, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )

      md <- sc@meta.data
      md <- merge(md, selected_scale(), by = "Cell_Index", all.x = T, sort = F)
      rownames(md) <- md$Cell_Index
      sc@meta.data <- md
      sc@meta.data$scale <- (as.numeric(md[, names(md) %in% input$Var_to_col_marker]))

      sc <- subset(sc, UMAP_1 > input$Filter_lower_UMAP1_marker)
      sc <- subset(sc, UMAP_1 < input$Filter_lower_UMAP1_marker2)

      sc <- subset(sc, UMAP_2 > input$Filter_lower_UMAP2_marker)
      sc <- subset(sc, UMAP_2 < input$Filter_lower_UMAP2_marker2)
      sc2 <- suppressWarnings(subset(sc, v_gene_AG != "NA"))

      sc2@meta.data$pos_neg <- ifelse(sc2@meta.data$scale > input$cutoff_marker_gt, "pos", "neg")
      sc2@meta.data
    })


    compare.stat_marker <- reactive({
      sc <- stats_marker_one()

      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )

      Idents(object = sc) <- sc@meta.data$pos_neg
      min.pct.expression <- input$min_point_ # standard setting: 0.25
      min.logfc <- input$LogFC_ # 0.25 is standard

      markers.fm.list <- FindMarkers(sc, ident.1 = "pos", min.pct = min.pct.expression, logfc.threshold = min.logfc, only.pos = TRUE)
      as.data.frame(markers.fm.list)
    })

    output$Compare.stat_marker <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      sc <- data_sc_pro()

      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )
      compare.stat_marker()
    })

    # download tables ----

    output$downloaddf_clonotype_distribution <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Var_to_col_marker, "_clonotypes_", input$pos_expanded_cut_off, "_stats", x, ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(merged_marker_hist())
        write.csv(df, file, row.names = T)
      }
    )

    output$downloaddf_Marker_stats <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Var_to_col_marker, "_Exp_cuttoff", input$pos_expanded_cut_off, "_stats", x, ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(compare.stat_marker())
        write.csv(df, file, row.names = T)
      }
    )

    # display T cells assoicated with that marker. ----
    output$AddInAnnoUI_man_1 <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload Files"
        )
      )
      req(sc)
      df3.meta <- sc@meta.data
      selectInput("AddInAnnoUI_man_1", "Selected Function", choices = names(df3.meta), selected = "T_cells")
    })

    output$AddInAnnoUI_man_2 <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload Files"
        )
      )
      req(sc)
      df3.meta <- sc@meta.data
      selectInput("AddInAnnoUI_man_2", "Selected Group", choices = names(df3.meta), selected = "Sample_Name")
    })

    selected_anno_col <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload Files"
        )
      )
      req(sc, input$AddInAnnoUI_man_1)
      md <- sc@meta.data
      md.selected <- unique(md[, names(md) %in% input$AddInAnnoUI_man_1])
      md.selected
    })

    observe({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload Files"
        )
      )
      req(sc, input$AddInAnnoUI_man_1)
      df3.meta <- selected_anno_col()

      updateSelectizeInput(
        session,
        "Slected_annotation",
        choices = df3.meta,
        selected = unlist(df3.meta)[1]
      )
    })

    Common_Annotation_code_man <- reactive({
      sc <- UMAP_metadata_with_labs()

      validate(
        need(
          nrow(sc) > 0,
          "Upload Files"
        )
      )
      req(sc, input$AddInAnnoUI_man_1, input$AddInAnnoUI_man_2, input$Samp_col, input$V_gene_sc)
      top_BD_cluster <- sc@meta.data
      top_BD_cluster$Selected_function <- top_BD_cluster[, names(top_BD_cluster) %in% input$AddInAnnoUI_man_1]
      top_BD_cluster$Selected_group <- top_BD_cluster[, names(top_BD_cluster) %in% input$AddInAnnoUI_man_2]
      top_BD_cluster$Selected_TCR <- top_BD_cluster[, names(top_BD_cluster) %in% input$V_gene_sc]
      top_BD_cluster$Selected_indiv <- top_BD_cluster[, names(top_BD_cluster) %in% input$Samp_col]
      top_BD_cluster$colouring_col <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this]
      top_BD_cluster$ID_Column <- top_BD_cluster[, names(top_BD_cluster) %in% input$Split_group_by_]

      top_BD_cluster$cloneCount <- 1
      rownames(top_BD_cluster) <- top_BD_cluster$Cell_Index
      top_BD_cluster
    })


    Prior_Anno_dt_process_man <- reactive({
      top_BD_cluster <- Common_Annotation_code_man()
      top_BD_cluster2 <- top_BD_cluster[, names(top_BD_cluster) %in% c("cloneCount", "Selected_function", "Selected_group", "Selected_indiv", "Selected_TCR")]
      top_BD_cluster2 <- top_BD_cluster2 %>%
        select(cloneCount, everything())
      df2 <- as.data.frame(ddply(top_BD_cluster2, names(top_BD_cluster2)[-c(1)], numcolwise(sum)))
      df2 <- df2[!is.na(df2$Selected_TCR), ]

      # df2 <- df2[complete.cases(df2)==T,]
      df2$fraction <- df2$cloneCount / sum(df2$cloneCount)
      df2$Percent <- round(df2$cloneCount / sum(df2$cloneCount) * 100, 2)
      df2$Selected_function <- ifelse(grepl("NA", df2$Selected_function), "-", df2$Selected_function)
      df2[order(df2$Percent, decreasing = T), ]
    })



    output$AnnoTable_perMarkers <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload Files"
        )
      )
      md <- Prior_Anno_dt_process_man()
      req(md)
      md
    })

    output$downloaddf_AnnoTable_perMarkers <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$Var_to_col_marker, "_clonotypes_", input$pos_expanded_cut_off, "_stats", x, ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(merged_marker_hist())
        write.csv(df, file, row.names = T)
      }
    )
    # Dual marker analysis -----

    observeEvent(input$load_marker_genes_dual, {
      sc <- (MainTcell_counts_names())
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )
      req(sc)

      updateSelectizeInput(
        session,
        "Var_to_col_marker2",
        choices = sc$V1,
        select = "CD8B"
      )
    })

    observeEvent(input$load_marker_genes_dual, {
      sc <- (MainTcell_counts_names())
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
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
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )

      md <- sc@meta.data
      md <- merge(md, selected_scale(), by = "Cell_Index", all.x = T, sort = F)
      rownames(md) <- md$Cell_Index
      md$markerX <- as.numeric(md[, names(md) %in% input$Var_to_col_marker2])
      md$markerX <- ifelse(md$markerX == "-Inf", NA, md$markerX)
      md$markerY <- as.numeric(md[, names(md) %in% input$Var_to_col_marker3])
      md$markerY <- ifelse(md$markerY == "-Inf", NA, md$markerY)
      md
      sc@meta.data <- md

      sc <- subset(sc, UMAP_1 > input$Filter_lower_UMAP1_marker)
      sc <- subset(sc, UMAP_1 < input$Filter_lower_UMAP1_marker2)

      sc <- subset(sc, UMAP_2 > input$Filter_lower_UMAP2_marker)
      sc <- subset(sc, UMAP_2 < input$Filter_lower_UMAP2_marker2)


      sc@meta.data$Quad <- ifelse(sc@meta.data$markerX < input$X_axis_dot_dual & sc@meta.data$markerY > input$Y_axis_dot_dual, "Q1",
                                  ifelse(sc@meta.data$markerX > input$X_axis_dot_dual & sc@meta.data$markerY > input$Y_axis_dot_dual, "Q2",
                                         ifelse(sc@meta.data$markerX < input$X_axis_dot_dual & sc@meta.data$markerY < input$Y_axis_dot_dual, "Q3",
                                                ifelse(sc@meta.data$markerX > input$X_axis_dot_dual & sc@meta.data$markerY < input$Y_axis_dot_dual, "Q4", "other")
                                         )
                                  )
      )

      sc2 <- suppressWarnings(subset(sc, v_gene_AG != "NA"))
      sc2
    })

    output$meta_data_for_features_scale2_df <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1, 2, 5, 10, 20, 50, 100), pageLength = 1, scrollX = TRUE), {
      sc <- meta_data_for_features_scale2()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )

      md <- sc@meta.data
      md
    })

    # umap dual scale -------
    marker_selected_UMAP_scale2 <- reactive({
      sc2 <- meta_data_for_features_scale2()
      validate(
        need(
          nrow(sc2) > 0,
          "Upload File"
        )
      )

      umap.meta <- sc2@meta.data
      umap.meta_pos <- subset(umap.meta, umap.meta$markerX != "-Inf")

      umap.meta$markerX <- ifelse(umap.meta$markerX == "-Inf", NA, umap.meta$markerX)
      umap.meta$markerX <- as.numeric(umap.meta$markerX)
      umap.meta_pos$markerX <- as.numeric(umap.meta_pos$markerX)

      ggplot(umap.meta, aes(x = UMAP_1, y = UMAP_2, color = markerX)) +
        geom_point(size = 1) +
        geom_point(data = umap.meta_pos, aes(x = UMAP_1, y = UMAP_2, color = markerX), size = 1) +
        # scale_color_continuous(type = "viridis",na.value="grey75") +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
        ) +
        scale_color_distiller(direction = 1, palette = input$col_marker_scale, na.value = "grey85", limits = c(min(0), max(input$max_scale2))) +
        theme(
          plot.title = element_text(colour = "black", family = input$font_type, size = 24, face = "bold", vjust = .5),
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.title.x = element_text(colour = "black", angle = 0, vjust = .5, face = "plain", family = input$font_type, size = input$title.text.sizer2),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        ) +
        labs(title = input$Var_to_col_marker2)
    })

    marker_selected_UMAP_scale3 <- reactive({
      sc2 <- meta_data_for_features_scale2()
      validate(
        need(
          nrow(sc2) > 0,
          "Upload File"
        )
      )

      umap.meta <- sc2@meta.data
      umap.meta_pos <- subset(umap.meta, umap.meta$markerY != "-Inf")

      umap.meta$markerY <- ifelse(umap.meta$markerY == "-Inf", NA, umap.meta$markerY)
      umap.meta$markerY <- as.numeric(umap.meta$markerY)
      umap.meta_pos$markerY <- as.numeric(umap.meta_pos$markerY)

      ggplot(umap.meta, aes(x = UMAP_1, y = UMAP_2, color = markerY)) +
        geom_point(size = 1) +
        geom_point(data = umap.meta_pos, aes(x = UMAP_1, y = UMAP_2, color = markerY), size = 1) +
        # scale_color_continuous(type = "viridis",na.value="grey75") +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
        ) +
        scale_color_distiller(direction = 1, palette = input$col_marker_scale, na.value = "grey85", limits = c(min(0), max(input$max_scale3))) +
        theme(
          plot.title = element_text(colour = "black", family = input$font_type, size = 24, face = "bold", vjust = .5),
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.title.x = element_text(colour = "black", angle = 0, vjust = .5, face = "plain", family = input$font_type, size = input$title.text.sizer2),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        ) +
        labs(title = input$Var_to_col_marker3)
    })

    output$marker_selected_UMAP_plot2 <- renderPlot({
      umap.meta <- MainTcell_counts_names()
      validate(
        need(
          nrow(umap.meta) > 0,
          "Upload File"
        )
      )
      marker_selected_UMAP_scale2()
    })
    output$marker_selected_UMAP_plot3 <- renderPlot({
      umap.meta <- MainTcell_counts_names()
      validate(
        need(
          nrow(umap.meta) > 0,
          "Upload File"
        )
      )
      marker_selected_UMAP_scale3()
    })

    # umap dual scale -------
    df_dotplot_marker <- reactive({
      sc2 <- meta_data_for_features_scale2()
      validate(
        need(
          nrow(sc2) > 0,
          "Upload File"
        )
      )

      umap.meta <- sc2@meta.data

      selected.col <- umap.meta[, names(umap.meta) %in% input$Colour_By_this]

      umap.meta$selected <- umap.meta[, names(umap.meta) %in% input$Colour_By_this]
      umap.meta$selected <- as.character(umap.meta$selected)
      umap.meta$selected <- gsub("NA", NA, umap.meta$selected)

      len.colour <- length(unique(selected.col))
      colour_markersby <- rainbow(len.colour)

      umap.meta$markerX <- umap.meta[, names(umap.meta) %in% input$Var_to_col_marker2]
      umap.meta$markerY <- umap.meta[, names(umap.meta) %in% input$Var_to_col_marker3]

      Q1 <- umap.meta[umap.meta$markerX < input$X_axis_dot_dual & umap.meta$markerY > input$Y_axis_dot_dual, ]
      Q2 <- umap.meta[umap.meta$markerX > input$X_axis_dot_dual & umap.meta$markerY > input$Y_axis_dot_dual, ]
      Q3 <- umap.meta[umap.meta$markerX < input$X_axis_dot_dual & umap.meta$markerY < input$Y_axis_dot_dual, ]
      Q4 <- umap.meta[umap.meta$markerX > input$X_axis_dot_dual & umap.meta$markerY < input$Y_axis_dot_dual, ]

      Q1.per <- dim(Q1)[1] / dim(umap.meta)[1] * 100
      Q2.per <- dim(Q2)[1] / dim(umap.meta)[1] * 100
      Q3.per <- dim(Q3)[1] / dim(umap.meta)[1] * 100
      Q4.per <- dim(Q4)[1] / dim(umap.meta)[1] * 100

      annotations <- data.frame(
        xpos = c(-Inf, -Inf, Inf, Inf),
        ypos = c(-Inf, Inf, -Inf, Inf),
        annotateText = c(
          paste0("Q3: ", round(Q3.per, 2), "%"), paste0("Q1: ", round(Q1.per, 2), "%"),
          paste0("Q4: ", round(Q4.per, 2), "%"), paste0("Q2: ", round(Q2.per, 2), "%")
        ),
        hjustvar = c(-0.1, -0.1, 1.1, 1.2),
        vjustvar = c(-0.4, 1.4, -0.4, 1.4)
      ) #<- adjust

      plot <- ggplot() +
        geom_point(data = umap.meta, aes(x = get(input$Var_to_col_marker2), y = get(input$Var_to_col_marker3), colour = selected)) +
        theme_bw() +
        geom_hline(yintercept = input$Y_axis_dot_dual) +
        geom_vline(xintercept = input$X_axis_dot_dual) +
        xlab(input$Var_to_col_marker2) +
        ylab(input$Var_to_col_marker3) +
        scale_colour_manual(values = colour_markersby, na.value = "grey85") +
        geom_text(data = annotations, aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar, label = annotateText), colour = "black", family = input$font_type, size = input$anno_text_size) +
        theme(
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 0),
          axis.title.x = element_text(colour = "black", angle = 0, vjust = .5, face = "plain", family = input$font_type, size = input$title.text.sizer2),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
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
        need(
          nrow(umap.meta) > 0,
          "Upload File"
        )
      )
      df_dotplot_marker()
    })

    output$downloadPlot_df_dotplot_marker_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("Dual_marker_plot_", input$Var_to_col_marker2, "_", input$Var_to_col_marker3, "_", x, ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = input$width_df_dotplot_marker_plot, height = input$height_df_dotplot_marker_plot, onefile = FALSE) # open the pdf device
        plot(df_dotplot_marker())
        dev.off()
      }, contentType = "application/pdf"
    )

    output$downloadPlotPNG_df_dotplot_marker_plot <- downloadHandler(
      filename = function() {
        x <- today()
        paste("Dual_marker_plot_", input$Var_to_col_marker2, "_", input$Var_to_col_marker3, "_", x, ".png", sep = "")
      },
      content = function(file) {
        png(file,
            width = input$width_png_df_dotplot_marker_plot,
            height = input$height_png_df_dotplot_marker_plot,
            res = input$resolution_PNG_df_dotplot_marker_plot
        )
        plot(df_dotplot_marker())
        dev.off()
      }, contentType = "application/png" # MIME type of the image
    )

    # summary of TCR per quad -----
    dual_maker_TCR_Sum <- reactive({
      sc2 <- meta_data_for_features_scale2()
      validate(
        need(
          nrow(sc2) > 0,
          "Upload File"
        )
      )
      req(sc2)
      umap.meta <- sc2@meta.data

      blank <- as.data.frame(names(umap.meta))
      names(blank) <- "V1"
      rownames(blank) <- blank$V1
      blank$blank1 <- "blank"
      blank$blank2 <- "blank"
      blank$blank3 <- "blank"
      blank$blank4 <- "blank"
      blank <- blank[, -c(1)]
      blank <- as.data.frame(t(blank))
      head(blank)

      list1 <- c(input$Y_axis_dot_dual - 1, input$Y_axis_dot_dual - 1, input$Y_axis_dot_dual + 1, input$Y_axis_dot_dual + 1)
      blank$markerY <- list1

      list2 <- c(input$X_axis_dot_dual - 1, input$X_axis_dot_dual + 1, input$X_axis_dot_dual - 1, input$X_axis_dot_dual + 1)
      blank$markerX <- list2
      blank

      umap.meta <- rbind(umap.meta, blank)
      umap.meta

      Q1 <- umap.meta[umap.meta$markerX < input$X_axis_dot_dual & umap.meta$markerY > input$Y_axis_dot_dual, ]
      Q1$Q1 <- 1
      Q2 <- umap.meta[umap.meta$markerX > input$X_axis_dot_dual & umap.meta$markerY > input$Y_axis_dot_dual, ]
      Q2$Q2 <- 1
      Q3 <- umap.meta[umap.meta$markerX < input$X_axis_dot_dual & umap.meta$markerY < input$Y_axis_dot_dual, ]
      Q3$Q3 <- 1
      Q4 <- umap.meta[umap.meta$markerX > input$X_axis_dot_dual & umap.meta$markerY < input$Y_axis_dot_dual, ]
      Q4$Q4 <- 1

      # as.data.frame(is.data.frame(subset(umap.meta, markerX>input$X_axis_dot_dual & markerY> input$Y_axis_dot_dual)))


      Q1_quad <- Q1[, names(Q1) %in% c(input$V_gene_sc, "Q1")]
      Q1_quad <- ddply(Q1_quad, input$V_gene_sc, numcolwise(sum))

      Q2_quad <- Q2[, names(Q2) %in% c(input$V_gene_sc, "Q2")]
      Q2_quad <- ddply(Q2_quad, input$V_gene_sc, numcolwise(sum))

      Q3_quad <- Q3[, names(Q3) %in% c(input$V_gene_sc, "Q3")]
      Q3_quad <- ddply(Q3_quad, input$V_gene_sc, numcolwise(sum))

      Q4_quad <- Q4[, names(Q4) %in% c(input$V_gene_sc, "Q4")]
      Q4_quad <- ddply(Q4_quad, input$V_gene_sc, numcolwise(sum))

      Q1_Q2 <- merge(Q1_quad, Q2_quad, by = input$V_gene_sc, all = T)
      Q3_Q4 <- merge(Q3_quad, Q4_quad, by = input$V_gene_sc, all = T)
      umeta.TCR <- merge(Q1_Q2, Q3_Q4, by = input$V_gene_sc, all = T)

      umeta.TCR[is.na(umeta.TCR)] <- 0
      umeta.TCR_needed.names <- names(umeta.TCR)
      umeta.TCR$total <- rowSums(umeta.TCR[2:5])
      rownames(umeta.TCR) <- umeta.TCR[, c(1)]
      umeta.TCR <- umeta.TCR[!(row.names(umeta.TCR) %in% c("blank")), ]
      umeta.TCR <- umeta.TCR[, -c(1)]
      names(umeta.TCR) <- c(
        paste0(input$Var_to_col_marker3, "pos", input$Var_to_col_marker2, "neg"),
        paste0(input$Var_to_col_marker3, "pos", input$Var_to_col_marker2, "pos"),
        paste0(input$Var_to_col_marker3, "neg", input$Var_to_col_marker2, "neg"),
        paste0(input$Var_to_col_marker3, "neg", input$Var_to_col_marker2, "pos"),
        "total"
      )

      umeta.TCR
    })


    output$dual_maker_TCR_Sum_DT <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1, 2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      umeta.TCR <- dual_maker_TCR_Sum()
      umeta.TCR
    })

    output$Dule_marker_TCRsummary_DT <- downloadHandler(
      filename = function() {
        x <- today()
        paste("Dule_marker_TCRsummary_", input$Var_to_col_marker2, "_", input$Var_to_col_marker3, "_", x, ".csv", sep = "")
      },
      content = function(file) {
        df <- dual_maker_TCR_Sum()
        write.csv(df, file)
      }
    )

    #### STATs for selected qudrant ------

    stats_dual_markers_tb1 <- reactive({
      sc <- meta_data_for_features_scale2()

      validate(
        need(
          nrow(sc) > 0,
          error_message_val_sc
        )
      )

      Idents(object = sc) <- sc@meta.data$Quad
      min.pct.expression <- input$min_point_ # standard setting: 0.25
      min.logfc <- input$LogFC_ # 0.25 is standard

      markers.fm.list <- FindMarkers(sc, ident.1 = input$quad_dualmarker, min.pct = min.pct.expression, logfc.threshold = min.logfc, only.pos = TRUE)
      as.data.frame(markers.fm.list)
    })


    output$dual_maker_TCR_statsTB <- DT::renderDT(escape = FALSE, filter = list(position = "top", clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(1, 2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      umeta.TCR <- stats_dual_markers_tb1()
      umeta.TCR
    })

    output$Dule_marker_statsTBDownload <- downloadHandler(
      filename = function() {
        x <- today()
        paste(input$quad_dualmarker, "Dule_marker_stats_from_", input$Var_to_col_marker2, "_", input$Var_to_col_marker3, "_", x, ".csv", sep = "")
      },
      content = function(file) {
        df <- stats_dual_markers_tb1()
        write.csv(df, file)
      }
    )


    # Add ridge plot for the distribution...

    #

    # updating the UI for prioritization -------
    output$cut.off_expanded2 <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val_UMAP
        )
      )

      req(input$Samp_col, input$V_gene_sc)
      TCR_Expanded_Df <- TCR_Expanded_fun(sc, (input$Samp_col), (input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1

      total.condition <- ddply(TCR_Expanded_Df, c("samp.count"), numcolwise(sum))
      total.condition$count_total <- total.condition$samp.count * total.condition$obs
      total.condition$frequency <- total.condition$count_total / sum(total.condition$count_total)
      total.condition <- total.condition %>% mutate(cum_freq = cumsum(frequency))
      total.condition <- subset(total.condition, total.condition$cum_freq > input$cutoff.expanded)
      num.cut.off <- total.condition$samp.count[1]

      # df3.meta <- df3.meta[!grepl("RNA",df3.meta) & !grepl("BCR",df3.meta) & !grepl("TCR",df3.meta)& !grepl("_gene",df3.meta) & !grepl("allele",df3.meta) & !grepl("percent",df3.meta) & !grepl("cdr3",df3.meta)]
      numericInput("cut.off_expanded", "Cut off greater than", value = num.cut.off, step = 1, min = 1)
    })
    output$Expanded.dotplot.cutoffs <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )

      TCR_Expanded_Df <- TCR_Expanded_fun(sc, (input$Samp_col), (input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      TCR_Expanded_Df2 <- subset(TCR_Expanded_Df, TCR_Expanded_Df$percent > input$cut.off_percent_rep)
      observations <- sum(TCR_Expanded_Df2$obs)
      percentage <- sum(TCR_Expanded_Df2$percent)
      length.samp.ID <- length(unique(TCR_Expanded_Df$ID_Column))
      # print(length.samp.ID)
      mat <- mat_sum(sc, input$Samp_col, input$V_gene_sc)

      if (observations > 0) {
        fluidRow(
          column(6, selectInput("restrict.dotpot", "Restrict to top list", choices = c("no", "yes"), selected = "no")),
          column(6, numericInput("restrict.dotpot.num", "Total genes to display:", value = 40)),
        )
      } else {
        sc <- UMAP_metadata_with_labs()
        validate(
          need(
            nrow(sc) > 0,
            error_message_val_UMAP
          )
        )

        markers.fm.list2 <- Vals_expanded.stats()
        req(markers.fm.list2)
        if (length(rownames(markers.fm.list2)) > 40) {
          selected.restrict.dotpot <- "yes"
          value.restrict.dotpot.num <- 40

          fluidRow(
            column(6, selectInput("restrict.dotpot", "Restrict to top list", choices = c("no", "yes"), selected = selected.restrict.dotpot)),
            column(6, numericInput("restrict.dotpot.num", "Total genes to display:", value = value.restrict.dotpot.num)),
          )
        } else {
          fluidRow(
            column(6, selectInput("restrict.dotpot", "Restrict to top list", choices = c("no", "yes"), selected = "no")),
            column(6, numericInput("restrict.dotpot.num", "Total genes to display:", value = 40)),
          )
        }
      }

      # df3.meta <- df3.meta[!grepl("RNA",df3.meta) & !grepl("BCR",df3.meta) & !grepl("TCR",df3.meta)& !grepl("_gene",df3.meta) & !grepl("allele",df3.meta) & !grepl("percent",df3.meta) & !grepl("cdr3",df3.meta)]
    })

    ## prioritising and automating the analysis ------

    output$Simple_workflow_step1 <- renderPrint({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )

      TCR_Expanded_Df <- TCR_Expanded_fun(sc, (input$Samp_col), (input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      length.samp.ID <- length(unique(TCR_Expanded_Df$ID_Column))
      # print(length.samp.ID)
      mat <- mat_sum(sc, input$Samp_col, input$V_gene_sc)
      if (max(mat$TotalSamps) == 1 && length.samp.ID == 1) {
        print("one individual and one sample")

        TCR_Expanded_Df <- TCR_Expanded_fun(sc, (input$Samp_col), (input$V_gene_sc))
        TCR_Expanded_Df$obs <- 1
        TCR_Expanded_Df2 <- subset(TCR_Expanded_Df, TCR_Expanded_Df$percent > input$cut.off_percent_rep)
        observations <- sum(TCR_Expanded_Df2$obs)
        percentage <- sum(TCR_Expanded_Df2$percent)
        if (observations > 0) {
          print("ImmunoDom")
          # observations <-  Frequency_expanded_df$obs[Frequency_expanded_df$percent>10 & Frequency_expanded_df$Frequency_expanded %in% "5. Gigantic (0.1 > X <= 0.5)"]
          # percent <-  Frequency_expanded_df$percent[Frequency_expanded_df$percent>10 & Frequency_expanded_df$Frequency_expanded %in% "5. Gigantic (0.1 > X <= 0.5)"]

          print(paste0("There are ", print(observations), " immuno dominant clonotype(s) that account for ", round(percentage, 2), "% of the repertoire"))
          # print(paste0("This process will download (1) UMAP count and top plots, bar plot of immunodominat (>10% of repertoire)"))

          # cat("test")
        } else {
          print("Polyclonal")
        }
      } else if (max(mat$TotalSamps) > 1 || length.samp.ID > 1) {
        print("multiple individuals or samples")
      } else {
        print("other")
      }
    })

    ### UI outputs -----
    output$Module_case_statements <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )

      TCR_Expanded_Df <- TCR_Expanded_fun(sc, (input$Samp_col), (input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      TCR_Expanded_Df2 <- subset(TCR_Expanded_Df, TCR_Expanded_Df$percent > input$cut.off_percent_rep)
      observations <- sum(TCR_Expanded_Df2$obs)
      length.samp.ID <- length(unique(TCR_Expanded_Df$ID_Column))
      print(length.samp.ID)

      mat <- mat_sum(sc, input$Samp_col, input$V_gene_sc)

      if (max(mat$TotalSamps) == 1 && length.samp.ID == 1) {
        print("one individual and one sample")
        if (observations > 0) {
          fluidRow(
            column(12, actionButton("ImmDom_download_buttonOneOne", "Download ImmunoDom (1 & 1) analysis"))
          )
        } else {
          fluidRow(
            # column(12,selectInput("AnalysisType","Preset parameters for", = c("ImmunoDom"))),
            column(12, actionButton("Poly_download_buttonOneOne", "Download Polyclonal (1 & 1) analysis"))
          )
        }
      } else if (max(mat$TotalSamps) > 1 || length.samp.ID > 1) {
        # BD_sum <- Top_clonotypes_multiCounts()

        fluidRow(
          column(3, numericInput("cut.off_percent_repMulti", "Priority cut-off", value = 1, step = 0.001, min = 0, max = 1)),
          column(3, numericInput("CloneTotal_input_top2", "Select clones > (Private)", value = 2, min = 2)),
          column(12, actionButton("Multi_download_button", "Download multi analysis"))
        )
      } else {
        # print("multiple individuals and multiple samples")
      }
    })

    output$Default_priority_cutoffAG <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )
      md <- sc@meta.data
      x <- today()
      if (length(data_sc_clusTCR_AG()) > 0) {
        df7 <- AG_cluster()

        df7 <- df7[order(df7$priority, decreasing = F), ]

        if (max(df7$Updated_order) > 10) {
          prior <- subset(df7, df7$Updated_order > 11)
          numericInput("priority_cutoff", "Priority cut-off (AG)", value = min(prior$priority), step = 0.01, min = 0, max = 1)
        } else {
          prior <- df7
          numericInput("priority_cutoff", "Priority cut-off (AG)", value = 1, step = 0.01, min = 0, max = 1)
        }
      } else {
        # numericInput("priority_cutoff","Priority cut-off (AG)",value = 1,step = 0.01, min = 0)
      }
    })
    output$Default_priority_cutoffBD <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload File"
        )
      )
      md <- sc@meta.data
      x <- today()
      if (length(data_sc_clusTCR_BD()) > 0) {
        df7 <- BD_cluster()

        df7 <- df7[order(df7$priority, decreasing = F), ]
        if (max(df7$Updated_order) > 10) {
          prior <- subset(df7, df7$Updated_order > 11)
          numericInput("priority_cutoffBD", "Priority cut-off (BD)", value = min(prior$priority), step = 0.01, min = 0, max = 1)
        } else {
          prior <- df7
          numericInput("priority_cutoffBD", "Priority cut-off (BD)", value = 1, step = 0.01, min = 0, max = 1)
        }
      } else {
      }
    })

    # one sample one indiv ImmunoDom -----
    observeEvent(input$ImmDom_download_buttonOneOne, {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )

      x <- today()
      message("Downloading the Summary table...")
      top.name.clonotypes <- paste("prioritization/ImmunoDom/Expansion_summary_table_", x, ".csv", sep = "")
      write.csv(Top_clonotype_df2(), top.name.clonotypes, row.names = F)

      message("Downloading the top UMAP...")
      top.name.clonotypes.count_png <- paste("prioritization/ImmunoDom/Expansion_UMAP_top_", x, ".png", sep = "")
      png(top.name.clonotypes.count_png, width = input$width_png_TCR.UMAP_top, height = input$height_png_TCR.UMAP_top, res = input$resolution_PNG_TCR.UMAP_top)
      plot(Topclonotypes())
      dev.off()

      message("Downloading the count UMAP...")
      top.name.clonotypes.top_png <- paste("prioritization/ImmunoDom/Expansion_UMAP_count_", x, ".png", sep = "")
      png(top.name.clonotypes.top_png, width = input$width_png_TCR.UMAP, height = input$height_png_TCR.UMAP, res = input$resolution_PNG_TCR.UMAP)
      plot(UMAP.TCRclonalit2())
      dev.off()

      message("Downloading the Dom bar plot...")
      top_clonotype_bar_code_immdom()

      message("Downloading Dom stats files and dot plot...")
      top_clone_FindMaker_looped()

      message("Downloading AG cluster table")

      top.name.clonotypes <- paste("prioritization/ImmunoDom/Cluster_summary_table_AG", x, ".csv", sep = "")
      write.csv(clusTCR2_df(), top.name.clonotypes, row.names = F)
    })
    top_clonotype_bar_code_immdom <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload"
        )
      )
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[, names(df3.meta) %in% input$V_gene_sc]

      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload"
        )
      )

      req(input$Samp_col, input$V_gene_sc)
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[, names(df3.meta) %in% input$V_gene_sc]

      BD_sum <- Top_clonotype_df2()
      BD_sum <- subset(BD_sum, BD_sum$Total_count > 1)

      TCR_Expanded_Df <- TCR_Expanded_fun(sc, (input$Samp_col), (input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      TCR_Expanded_Df2 <- subset(TCR_Expanded_Df, TCR_Expanded_Df$percent > input$cut.off_percent_rep)
      observations <- sum(TCR_Expanded_Df2$obs)
      percentage <- sum(TCR_Expanded_Df2$percent)

      withProgress(message = "Performing ImmDom Analysis", value = 0, {
        for (i in 1:observations) {
          incProgress(1 / observations, detail = paste("Clone", i, "of", observations))
          message(BD_sum$cluster_name[i])
          name.clone <- BD_sum$cluster_name[i]
          top_BD_clonotype <- df3.meta[df3.meta$cluster_name %in% name.clone, ]
          # print(top_BD_clonotype)

          dtop_clonotype_bar_code <- top_BD_clonotype

          # req(input$Graph_split_order)

          dtop_clonotype_bar_code$Selected_group <- dtop_clonotype_bar_code[, names(dtop_clonotype_bar_code) %in% input$Split_group_by_]
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
          } else {

          }


          colorblind_vector <- as.data.frame(colorblind_vector)
          names(colorblind_vector) <- "cols"

          dtop_clonotype_bar_code$Selected_chain2 <- dtop_clonotype_bar_code[, names(dtop_clonotype_bar_code) %in% input$V_gene_sc]
          dtop_clonotype_bar_code$Selected_chain3 <- gsub("_", " ", dtop_clonotype_bar_code$Selected_chain2)
          dtop_clonotype_bar_code$Selected_chain3 <- gsub("[.]", " ", dtop_clonotype_bar_code$Selected_chain3)

          dtop_clonotype_bar_code <- dtop_clonotype_bar_code[dtop_clonotype_bar_code$Selected_group %in% input$Graph_split_order, ]
          dtop_clonotype_bar_code$Selected_group <- factor(dtop_clonotype_bar_code$Selected_group, levels = input$Graph_split_order)

          ggplot_plot <- ggplot(dtop_clonotype_bar_code, aes(x = Selected_group, fill = Selected_chain3, colour = Selected_chain3, label = Selected_chain3)) +
            geom_bar() +
            theme_bw() +
            scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector$cols, na.value = input$NA_col_analysis) +
            scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector$cols, na.value = input$NA_col_analysis) +
            # scale_alpha_manual(values = rep(1,length(unique(dtop_clonotype_bar_code$Selected_chain))), na.value=0.5)+
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(
              axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
              axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
              axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
              axis.title.x = element_blank(),
              legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
              legend.position = input$legend_position,
              legend.title = element_blank()
            ) +
            guides(color = "none", size = "none")

          x <- today()
          top.name.clonotypes.top_png <- paste("prioritization/ImmunoDom/", i, "_top_clone_", gsub("[/]", "", gsub("&", "", name.clone)), "_", x, ".png", sep = "")

          num_width <- length(unique(dtop_clonotype_bar_code$Selected_group))

          png(top.name.clonotypes.top_png, width = (num_width * 100 + 400), height = input$height_png_TCR.UMAP, res = input$resolution_PNG_TCR.UMAP)
          plot(ggplot_plot)
          dev.off()
        }
      })
    })
    top_clone_FindMaker_looped <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload"
        )
      )
      x <- today()
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[, names(df3.meta) %in% input$V_gene_sc]
      sc@meta.data$Vgene <- sc@meta.data[, names(sc@meta.data) %in% input$V_gene_sc]

      BD_sum <- Top_clonotype_df2()
      BD_sum <- subset(BD_sum, BD_sum$Total_count > 1)
      TCR_Expanded_Df <- TCR_Expanded_fun(sc, (input$Samp_col), (input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      TCR_Expanded_Df2 <- subset(TCR_Expanded_Df, TCR_Expanded_Df$percent > input$cut.off_percent_rep)
      observations <- sum(TCR_Expanded_Df2$obs)
      percentage <- sum(TCR_Expanded_Df2$percent)

      withProgress(message = "Performing FindMarker ImmDom Analysis", value = 0, {
        for (i in 1:observations) {
          incProgress(1 / observations, detail = paste("Clone", i, "of", observations))
          name.clone <- BD_sum$cluster_name[i]

          message(paste0("Downloading Dom stats files and dot plot...", name.clone))

          sc@meta.data$Gene_select <- ifelse(sc@meta.data$Vgene %in% name.clone, name.clone, "other")
          sc@meta.data
          Idents(object = sc) <- sc@meta.data$Gene_select

          min.pct.expression <- input$min_point_ # standard setting: 0.25
          min.logfc <- input$LogFC_ # 0.25 is standard
          # p.val.cutoff <-  input$pval_top #(1/10^3) is standard, use (1/10^0) to ignore

          cluster.names <- unique(Idents(sc))[order(unique(Idents(sc)))]
          # print(paste0("calculating markers for cluster ",name.clone,". Total: ",length(cluster.names)," clusters"))
          markers.fm.list <- FindMarkers(sc, ident.1 = name.clone, min.pct = min.pct.expression, logfc.threshold = min.logfc, only.pos = TRUE)
          markers.fm.list2 <- subset(markers.fm.list, markers.fm.list$p_val_adj < input$pval.ex.filter)

          if (length(markers.fm.list2$p_val_adj) > 0) {
            ### download the dot plot -------

            clonotype.name.stats <- paste("prioritization/ImmunoDom/", i, "_", gsub("[/]", "", gsub("&", "", name.clone)), "_stats_table", "_", today(), ".csv", sep = "")
            write.csv(markers.fm.list2, clonotype.name.stats, row.names = T)

            message(paste0("Saved csv", name.clone))
            list.names <- rownames(markers.fm.list2)

            if (length(rownames(markers.fm.list2)) > 40) {
              list.names <- list.names[1:40]
            } else {
              list.names <- rownames(markers.fm.list2)
            }

            size_legend <- input$Legend_size - 2

            plotdotplot <- DotPlot(sc, features = list.names) +
              RotatedAxis() +
              theme(
                axis.title.y = element_blank(),
                axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
                axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
                axis.title.x = element_blank(),
                legend.title = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
                legend.text = element_text(colour = "black", size = size_legend, family = input$font_type),
                legend.position = input$legend_position,
              ) +
              scale_colour_gradient2(low = input$low.dotplot, mid = input$middle.dotplot, high = input$high.dotplot) +
              scale_x_discrete(labels = label_wrap(20)) +
              scale_y_discrete(labels = label_wrap(20))


            file.name.clone <- paste("prioritization/ImmunoDom/", i, "_", gsub("[/]", "", gsub("&", "", name.clone)), "_dotplot_plot", "_", today(), ".png", sep = "")


            png(file.name.clone, width = input$width_png_all_expression_dotplot_top, height = input$height_png_all_expression_dotplot_top, res = input$resolution_PNG_all_expression_dotplot_top)
            plot(plotdotplot)
            dev.off()


            ##### download the OverRep ------

            df <- sc@meta.data
            # require()

            geneSet <- read.csv(system.file("OverRep", "GeneSets.csv", package = "STEGO.R"), header = T)

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

            if (input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
              geneSet$GeneSet <- gsub("-", ".", geneSet$GeneSet)
            }

            if (input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
              require(stringr)
              geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
            }

            for (j in 1:dim(geneSet)[1]) {
              # listed GeneSet
              message(paste("GeneSet: ", j))
              Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet, ";")[j])
              names(Gene.set.testing) <- "V1"
              Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
              names(Gene.set.testing2) <- "V1"
              background.overlap <- merge(Gene.set.testing2, background.genes.name, by = "V1")
              geneSet$background.geneset[j] <- length(background.overlap$V1)
              geneSet$background.geneset.name[j] <- as.character(paste(unlist(background.overlap[1]), collapse = ";"))
              # in sig gene list
              overlap <- merge(background.overlap, DEx.genes, by = "V1")

              geneSet$in.geneset[j] <- length(overlap$V1)
              geneSet$in.geneset.name[j] <- as.character(paste(unlist(overlap[1]), collapse = ";"))
            }

            geneSet2 <- subset(geneSet, geneSet$in.geneset > 0)
            message(paste(i, "has", length(geneSet2$in.geneset), "GeneSets in ImmunoDom"))

            if (length(geneSet2$in.geneset) > 0) {
              for (k in 1:dim(geneSet2)[1]) {
                tota.gene.set <- geneSet2$background.geneset[k] # genes that are identified in background
                tota.gene.set
                in.geneset <- geneSet2$in.geneset[k] # DEx in geneset

                background.genes
                not.in.total <- background.genes - tota.gene.set
                not.in.geneset.sig <- total.sig - in.geneset
                d <- data.frame(gene.in.interest = c(in.geneset, not.in.geneset.sig), gene.not.interest = c(tota.gene.set, not.in.total))
                row.names(d) <- c("In_category", "not_in_category")

                if (in.geneset > 0) {
                  geneSet2$p.val[k] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
                  geneSet2$lowerCI[k] <- unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
                  geneSet2$upperCI[k] <- unlist(fisher.test(d)$conf.int)[2]
                  geneSet2$OR[k] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1], 3)
                } else {
                  geneSet2$p.value[k] <- "-"
                  geneSet2$lowerCI[k] <- "-"
                  geneSet2$upperCI[k] <- "-"
                  geneSet2$OR[k] <- "-"
                }
              }

              geneSet2 <- geneSet2[order(geneSet2$p.val, decreasing = F), ]
              geneSet2 <- subset(geneSet2, geneSet2$in.geneset >= input$in.geneset.cutoff_top)
              geneSet2 <- subset(geneSet2, geneSet2$p.val <= input$p.val_cutoff_top)
              geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
              geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")

              file.name.clone <- paste("prioritization/ImmunoDom/", i, "_", gsub("[/]", "", gsub("&", "", name.clone)), "_dotplot_plot", "_", today(), ".png", sep = "")
              top.name.overrep <- paste("prioritization/ImmunoDom/", i, "_", gsub("[/]", "", gsub("&", "", name.clone)), "_OverRep", "_", today(), ".csv", sep = "")
              write.csv(geneSet2, top.name.overrep, row.names = F)
            }
          }
        }
      })
    })

    # one sample one indiv Polyclonal -----
    observeEvent(input$Poly_download_buttonOneOne, {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )

      TCR_Expanded_Df <- TCR_Expanded_fun(sc, (input$Samp_col), (input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      TCR_Expanded_Df2 <- subset(TCR_Expanded_Df, TCR_Expanded_Df$percent > input$cut.off_percent_rep)
      observations <- sum(TCR_Expanded_Df2$obs)

      x <- today()
      message("Downloading contig summary table...")
      top.name.clonotypes <- paste("prioritization/PolyClonal/Expansion_summary_table_", x, ".csv", sep = "")
      write.csv(Top_clonotype_df2(), top.name.clonotypes, row.names = F)

      message("Downloading the count UMAP...")
      top.name.clonotypes.top_png <- paste("prioritization/PolyClonal/Expansion_UMAP_count_", x, ".png", sep = "")
      png(top.name.clonotypes.top_png, width = input$width_png_TCR.UMAP, height = input$height_png_TCR.UMAP, res = input$resolution_PNG_TCR.UMAP)
      plot(UMAP.TCRclonalit2())
      dev.off()

      message("Downloading the expansion UMAP...")
      top.name.clonotypes.top_png <- paste("prioritization/PolyClonal/Exp_UMAP_cutoff_count.", input$cut.off_expanded, ".and.freq", input$cutoff.expanded, "_", x, ".png", sep = "")
      # png(top.name.clonotypes.top_png, width = input$width_png_TCR.UMAP,height = input$height_png_TCR.UMAP,res = input$resolution_PNG_TCR.UMAP)
      png(top.name.clonotypes.top_png,
          width = input$width_png_UMAP_Expanded,
          height = input$height_png_UMAP_Expanded,
          res = input$resolution_PNG_UMAP_Expanded
      )
      plot(UMAP_Expanded_plot())
      dev.off()

      message("Downloading stats table...")
      Exp_stats_cutoff_count.name <- paste("prioritization/PolyClonal/Exp_stats_cutoff_count.", input$cut.off_expanded, ".and.freq", input$cutoff.expanded, "_", x, ".csv", sep = "")
      write.csv(Vals_expanded.stats(), Exp_stats_cutoff_count.name, row.names = T)

      message("Downloading the expansion stats dotplot...")
      top.name.clonotypes.top_png <- paste("prioritization/PolyClonal/Exp_dotplot.", input$cut.off_expanded, ".and.freq", input$cutoff.expanded, "_", x, ".png", sep = "")
      png(top.name.clonotypes.top_png,
          width = input$width_png_all_expression_dotplot_ex,
          height = input$height_png_all_expression_dotplot_ex,
          res = input$resolution_PNG_all_expression_dotplot_ex
      )
      plot(relative_expression_plot_ex())
      dev.off()

      message("Downloading overrep table...")
      Exp_stats_cutoff_count.name <- paste("prioritization/PolyClonal/Exp_OverRep_cutoff_count.", input$cut.off_expanded, ".and.freq", input$cutoff.expanded, "_", x, ".csv", sep = "")
      write.csv(Over_rep_Exp(), Exp_stats_cutoff_count.name, row.names = F)
    })


    # multiple samples (either individuals or samples) ------
    Upset_plot_overlap_Multi <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )

      df <- sc@meta.data
      df <- as.data.frame(df)
      unique.df <- unique(df[, names(df) %in% c(input$Samp_col, input$V_gene_sc)])
      unique.df <- unique.df %>%
        select(input$Samp_col, everything())

      # count data of unique
      names(unique.df) <- c("group", "chain")
      unique.df <- subset(unique.df, unique.df$chain != "NA")
      unique.df <- subset(unique.df, unique.df$group != "NA")
      unique.df$cloneCount <- 1
      mat <- acast(unique.df, chain ~ group, value.var = "cloneCount")
      mat[is.na(mat)] <- 0
      Count_data <- as.data.frame(rowSums(mat))
      names(Count_data) <- "V1"

      # Total count
      unique.df <- (df[, names(df) %in% c(input$Samp_col, input$V_gene_sc)])
      unique.df <- unique.df %>%
        select(input$Samp_col, everything())

      names(unique.df) <- c("group", "chain")
      unique.df <- subset(unique.df, unique.df$chain != "NA")
      unique.df <- subset(unique.df, unique.df$group != "NA")
      unique.df$cloneCount <- 1
      mat <- acast(unique.df, chain ~ group, value.var = "cloneCount")
      mat[is.na(mat)] <- 0
      sum_data <- as.data.frame(rowSums(mat))
      names(sum_data) <- "V1"
      mat <- as.data.frame(mat)
      mat$TotalSamps <- Count_data$V1
      mat$CloneTotal <- sum_data$V1
      mat <- mat[order(mat$CloneTotal, decreasing = T), ]
      mat <- mat[order(mat$TotalSamps, decreasing = T), ]
      mat
    })

    Upset_plot_multi <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )
      df <- sc@meta.data
      df <- as.data.frame(df)
      unique.df <- unique(df[, names(df) %in% c(input$Samp_col, input$V_gene_sc)])
      unique.df <- unique.df %>%
        select(input$Samp_col, everything())

      names(unique.df) <- c("group", "chain")
      unique.df <- unique.df[unique.df$group %in% input$ID_Column_factor, ]

      unique.df <- subset(unique.df, unique.df$chain != "NA")
      unique.df <- subset(unique.df, unique.df$group != "NA")
      unique.df$cloneCount <- 1
      mat <- acast(unique.df, chain ~ group, value.var = "cloneCount")
      mat[is.na(mat)] <- 0
      mat <- as.data.frame(mat)
      df.x <- make_comb_mat(mat)
      list.names <- as.character(input$ID_Column_factor)

      ht <- draw(UpSet(df.x,
                       pt_size = unit(5, "mm"),
                       lwd = 1,
                       row_names_gp = gpar(fontfamily = input$font_type, fontsize = 12),
                       column_names_gp = gpar(fontfamily = input$font_type, fontsize = 12),
                       top_annotation = upset_top_annotation(df.x,
                                                             add_numbers = T,
                                                             numbers_gp = gpar(fontfamily = input$font_type, fontsize = 12),
                                                             annotation_name_gp = gpar(fontfamily = input$font_type, fontsize = 12),
                                                             gp = gpar(fill = "black"),
                       ),
                       right_annotation = upset_right_annotation(df.x,
                                                                 add_numbers = T,
                                                                 numbers_gp = gpar(fontfamily = input$font_type, fontsize = 12),
                                                                 annotation_name_gp = gpar(fontfamily = input$font_type, fontsize = 12),
                                                                 gp = gpar(fill = "black"),
                       ),
                       set_order = list.names
      ), padding = unit(c(20, 20, 20, 20), "mm"))
      ht
    })

    clonal_plot_multi <- reactive({
      df4 <- TCR_Expanded()
      df4 <- df4[df4$ID_Column %in% input$ID_Column_factor, ]
      df4$ID_Column <- as.character(df4$ID_Column)
      df4$ID_Column <- factor(df4$ID_Column, levels = input$ID_Column_factor)

      # df4 <- TCR_Expanded()
      df4 <- df4[order(df4[, names(df4) %in% input$Graph_type_bar]), ]
      df4
      col.df <- as.data.frame(unique(df4[, names(df4) %in% input$Graph_type_bar]))
      names(col.df) <- "V1"

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
      } else {

      }

      col.df$col <- colorblind_vector

      ggplot(df4, aes(x = ID_Column, y = frequency, fill = get(input$Graph_type_bar), colour = get(input$Graph_type_bar), label = get(input$Graph_type_bar))) +
        geom_bar(stat = "identity") +
        theme_bw() +
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 50), values = alpha(col.df$col, 1), na.value = input$NA_col_analysis) +
        scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 50), values = alpha(col.df$col, 1)) +
        theme(
          axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
          axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
          axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
          axis.title.x = element_blank(),
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )
    })

    Top_clonotypes_multiCounts <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )

      BD_sum <- Upset_plot_overlap_Multi()
      BD_sum <- subset(BD_sum, BD_sum$CloneTotal > 1)
      BD_sum <- subset(BD_sum, BD_sum$TotalSamps > 1)
      BD_sum$priority <- 1 / (BD_sum$CloneTotal * BD_sum$TotalSamps)
      if (input$restrict_to_expanded) {
        BD_sum$same <- ifelse(BD_sum$CloneTotal == BD_sum$TotalSamps, "NEx", "Ex")
        BD_sum <- subset(BD_sum, BD_sum$same == "Ex")
        BD_sum <- BD_sum[, !names(BD_sum) %in% "same"]
      }

      BD_sum$cluster_name <- rownames(BD_sum)
      BD_sum
    })

    Top_clonotypes_multiCounts_barplot <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload"
        )
      )
      x <- today()
      req(input$Samp_col, input$V_gene_sc, input$cut.off_percent_repMulti)
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[, names(df3.meta) %in% input$V_gene_sc]

      BD_sum <- Top_clonotypes_multiCounts()

      if (dim(BD_sum)[1] > 0) {
        BD_sum$obs <- 1
        BD_sum <- subset(BD_sum, BD_sum$priority < input$cut.off_percent_repMulti)

        top.name.clonotypes.top_png <- paste("prioritization/Multi/PublicLike/", "Selected_clones_", x, ".csv", sep = "")
        write.csv(BD_sum, top.name.clonotypes.top_png)
        observations <- sum(BD_sum$obs)

        withProgress(message = "Performing Multi Overlap Analysis (barplots)", value = 0, {
          for (i in 1:observations) {
            incProgress(1 / observations, detail = paste("Clone", i, "of", observations))
            message(BD_sum$cluster_name[i])
            name.clone <- BD_sum$cluster_name[i]
            top_BD_clonotype <- df3.meta[df3.meta$cluster_name %in% name.clone, ]
            # print(top_BD_clonotype)

            dtop_clonotype_bar_code <- top_BD_clonotype

            # req(input$Graph_split_order)

            dtop_clonotype_bar_code$Selected_group <- dtop_clonotype_bar_code[, names(dtop_clonotype_bar_code) %in% input$Split_group_by_]
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
            } else {

            }


            colorblind_vector <- as.data.frame(colorblind_vector)
            names(colorblind_vector) <- "cols"

            dtop_clonotype_bar_code$Selected_chain2 <- dtop_clonotype_bar_code[, names(dtop_clonotype_bar_code) %in% input$V_gene_sc]
            dtop_clonotype_bar_code$Selected_chain3 <- gsub("_", " ", dtop_clonotype_bar_code$Selected_chain2)
            dtop_clonotype_bar_code$Selected_chain3 <- gsub("[.]", " ", dtop_clonotype_bar_code$Selected_chain3)

            dtop_clonotype_bar_code <- dtop_clonotype_bar_code[dtop_clonotype_bar_code$Selected_group %in% input$Graph_split_order, ]
            dtop_clonotype_bar_code$Selected_group <- factor(dtop_clonotype_bar_code$Selected_group, levels = input$Graph_split_order)

            ggplot_plot <- ggplot(dtop_clonotype_bar_code, aes(x = Selected_group, fill = Selected_chain3, colour = Selected_chain3, label = Selected_chain3)) +
              geom_bar() +
              theme_bw() +
              scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector$cols, na.value = input$NA_col_analysis) +
              scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector$cols, na.value = input$NA_col_analysis) +
              # scale_alpha_manual(values = rep(1,length(unique(dtop_clonotype_bar_code$Selected_chain))), na.value=0.5)+
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
              theme(
                axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
                axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
                axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
                axis.title.x = element_blank(),
                legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
                legend.position = input$legend_position,
                legend.title = element_blank()
              ) +
              guides(color = "none", size = "none")


            top.name.clonotypes.top_png <- paste("prioritization/Multi/PublicLike/", i, "_top_clone_", gsub("[/]", "", gsub("&", "", name.clone)), "_", x, ".png", sep = "")

            num_width <- length(unique(dtop_clonotype_bar_code$Selected_group))

            png(top.name.clonotypes.top_png, width = (num_width * 100 + 500), height = input$height_png_TCR.UMAP, res = input$resolution_PNG_TCR.UMAP)
            plot(ggplot_plot)
            dev.off()
          }
        })
      } else {
        withProgress(message = "No multi-sample bar plots")
      }
    })

    top_clone_FindMaker_looped_Multi <- reactive({
      sc <- UMAP_metadata_with_labs()
      x <- today()

      validate(
        need(
          nrow(sc) > 0,
          "Upload"
        )
      )
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[, names(df3.meta) %in% input$V_gene_sc]
      sc@meta.data$Vgene <- sc@meta.data[, names(sc@meta.data) %in% input$V_gene_sc]
      BD_sum <- Top_clonotypes_multiCounts()
      BD_sum <- subset(BD_sum, BD_sum$CloneTotal > 2)

      if (dim(BD_sum)[1] > 0) {
        BD_sum <- subset(BD_sum, BD_sum$priority < input$cut.off_percent_repMulti)
        BD_sum$obs <- 1
        observations <- sum(BD_sum$obs)

        withProgress(message = "Performing Multi Analysis (FindMarkers, Dotplot, OverRep)", value = 0, {
          for (i in 1:observations) {
            incProgress(1 / observations, detail = paste("Clone", i, "of", observations))
            name.clone <- BD_sum$cluster_name[i]

            message(paste0("Downloading Dom stats files and dot plot...", name.clone))

            sc@meta.data$Gene_select <- ifelse(sc@meta.data$Vgene %in% name.clone, name.clone, "other")
            sc@meta.data
            Idents(object = sc) <- sc@meta.data$Gene_select

            min.pct.expression <- input$min_point_ # standard setting: 0.25
            min.logfc <- input$LogFC_ # 0.25 is standard
            # p.val.cutoff <-  input$pval_top #(1/10^3) is standard, use (1/10^0) to ignore

            cluster.names <- unique(Idents(sc))[order(unique(Idents(sc)))]
            # print(paste0("calculating markers for cluster ",name.clone,". Total: ",length(cluster.names)," clusters"))
            markers.fm.list <- FindMarkers(sc, ident.1 = name.clone, min.pct = min.pct.expression, logfc.threshold = min.logfc, only.pos = TRUE)
            markers.fm.list2 <- subset(markers.fm.list, markers.fm.list$p_val_adj < input$pval.ex.filter)
            message(paste(length(markers.fm.list2$p_val_adj), "total markers for cluster", i))
            if (length(markers.fm.list2$p_val_adj) > 0) {
              clonotype.name.stats <- paste("prioritization/Multi/PublicLike/", i, "_", gsub("[/]", "", gsub("&", "", name.clone)), "_stats_table", "_", today(), ".csv", sep = "")
              write.csv(markers.fm.list2, clonotype.name.stats, row.names = T)

              message(paste0("Saved csv ", name.clone))
              list.names <- rownames(markers.fm.list2)

              if (length(rownames(markers.fm.list2)) > 40) {
                list.names <- list.names[1:40]
              } else {
                list.names <- rownames(markers.fm.list2)
              }

              size_legend <- input$Legend_size - 2


              plotdotplot <- DotPlot(sc, features = list.names) +
                RotatedAxis() +
                theme(
                  axis.title.y = element_blank(),
                  axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
                  axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
                  axis.title.x = element_blank(),
                  legend.title = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
                  legend.text = element_text(colour = "black", size = size_legend, family = input$font_type),
                  legend.position = input$legend_position,
                ) +
                scale_colour_gradient2(low = input$low.dotplot, mid = input$middle.dotplot, high = input$high.dotplot) +
                scale_x_discrete(labels = label_wrap(20)) +
                scale_y_discrete(labels = label_wrap(20))


              file.name.clone <- paste("prioritization/Multi/PublicLike/", i, "_", gsub("[/]", "", gsub("&", "", name.clone)), "_dotplot_plot", "_", today(), ".png", sep = "")

              ### download the dot plot -------
              png(file.name.clone, width = input$width_png_all_expression_dotplot_top, height = input$height_png_all_expression_dotplot_top, res = input$resolution_PNG_all_expression_dotplot_top)
              plot(plotdotplot)
              dev.off()


              ##### download the OverRep ------
              geneSet <- read.csv(system.file("OverRep", "GeneSets.csv", package = "STEGO.R"), header = T)

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

              if (input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
                geneSet$GeneSet <- gsub("-", ".", geneSet$GeneSet)
              }

              if (input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
                require(stringr)
                geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
              }
              message(paste(i, "Performing Over rep analysis"))
              for (j in 1:dim(geneSet)[1]) {
                # listed GeneSet

                Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet, ";")[j])
                names(Gene.set.testing) <- "V1"
                Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
                names(Gene.set.testing2) <- "V1"
                background.overlap <- merge(Gene.set.testing2, background.genes.name, by = "V1")
                geneSet$background.geneset[j] <- length(background.overlap$V1)
                geneSet$background.geneset.name[j] <- as.character(paste(unlist(background.overlap[1]), collapse = ";"))
                # in sig gene list
                overlap <- merge(background.overlap, DEx.genes, by = "V1")

                geneSet$in.geneset[j] <- length(overlap$V1)
                geneSet$in.geneset.name[j] <- as.character(paste(unlist(overlap[1]), collapse = ";"))
              }

              geneSet2 <- subset(geneSet, geneSet$in.geneset > 0)
              message(paste(i, "has", length(geneSet2$in.geneset), "GeneSets in Multi"))

              if (length(geneSet2$in.geneset) > 0) {
                for (k in 1:dim(geneSet2)[1]) {
                  tota.gene.set <- geneSet2$background.geneset[k] # genes that are identified in background
                  tota.gene.set
                  in.geneset <- geneSet2$in.geneset[k] # DEx in geneset

                  background.genes
                  not.in.total <- background.genes - tota.gene.set
                  not.in.geneset.sig <- total.sig - in.geneset
                  d <- data.frame(gene.in.interest = c(in.geneset, not.in.geneset.sig), gene.not.interest = c(tota.gene.set, not.in.total))
                  row.names(d) <- c("In_category", "not_in_category")

                  if (in.geneset > 0) {
                    geneSet2$p.val[k] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
                    geneSet2$lowerCI[k] <- unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
                    geneSet2$upperCI[k] <- unlist(fisher.test(d)$conf.int)[2]
                    geneSet2$OR[k] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1], 3)
                  } else {
                    geneSet2$p.value[k] <- "-"
                    geneSet2$lowerCI[k] <- "-"
                    geneSet2$upperCI[k] <- "-"
                    geneSet2$OR[k] <- "-"
                  }
                }

                geneSet2 <- geneSet2[order(geneSet2$p.val, decreasing = F), ]
                geneSet2 <- subset(geneSet2, geneSet2$in.geneset >= input$in.geneset.cutoff_top)
                geneSet2 <- subset(geneSet2, geneSet2$p.val <= input$p.val_cutoff_top)
                geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
                geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")

                top.name.overrep <- paste("prioritization/Multi/PublicLike/", i, "_", gsub("[/]", "", gsub("&", "", name.clone)), "_OverRep", "_", today(), ".csv", sep = "")
                write.csv(geneSet2, top.name.overrep, row.names = F)
              }
            }
          }
        })
      } else {
        withProgress(message("No Markers"))
      }
    })

    # download the top x per individual
    Top_clonotypes_Private <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )
      BD_sum <- Upset_plot_overlap_Multi()
      BD_sum <- subset(BD_sum, BD_sum$CloneTotal > input$CloneTotal_input_top2)
      BD_sum <- subset(BD_sum, BD_sum$TotalSamps == 1)
      BD_sum$cluster_name <- rownames(BD_sum)
      BD_sum
    })
    top_clone_FindMaker_looped_Private <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload"
        )
      )
      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[, names(df3.meta) %in% input$V_gene_sc]
      sc@meta.data$Vgene <- sc@meta.data[, names(sc@meta.data) %in% input$V_gene_sc]

      BD_sum <- Top_clonotypes_Private()
      BD_sum$obs <- 1
      observations <- sum(BD_sum$obs)

      withProgress(message = "Performing Multi Analysis (FindMarkers, Dotplot, OverRep)", value = 0, {
        for (i in 1:observations) {
          incProgress(1 / observations, detail = paste("Clone", i, "of", observations))
          name.clone <- BD_sum$cluster_name[i]

          message(paste0("Downloading Dom stats files and dot plot...", name.clone))
          ###### Stats table ------
          sc@meta.data$Gene_select <- ifelse(sc@meta.data$Vgene %in% name.clone, name.clone, "other")
          sc@meta.data
          Idents(object = sc) <- sc@meta.data$Gene_select

          min.pct.expression <- input$min_point_ # standard setting: 0.25
          min.logfc <- input$LogFC_ # 0.25 is standard
          # p.val.cutoff <-  input$pval_top #(1/10^3) is standard, use (1/10^0) to ignore

          cluster.names <- unique(Idents(sc))[order(unique(Idents(sc)))]
          # print(paste0("calculating markers for cluster ",name.clone,". Total: ",length(cluster.names)," clusters"))
          markers.fm.list <- FindMarkers(sc, ident.1 = name.clone, min.pct = min.pct.expression, logfc.threshold = min.logfc, only.pos = TRUE)
          markers.fm.list2 <- subset(markers.fm.list, markers.fm.list$p_val_adj < input$pval.ex.filter)
          x <- today()

          message(paste(length(markers.fm.list2$p_val_adj), "total markers for cluster", i))

          if (length(markers.fm.list2$p_val_adj) > 0) {
            clonotype.name.stats <- paste("prioritization/Multi/Unique/", i, "_", gsub("[/]", "", gsub("&", "", name.clone)), "_stats_table", "_", today(), ".csv", sep = "")
            write.csv(markers.fm.list2, clonotype.name.stats, row.names = T)

            message(paste0("Saved csv ", name.clone))
            list.names <- rownames(markers.fm.list2)

            if (length(rownames(markers.fm.list2)) > 40) {
              list.names <- list.names[1:40]
            } else {
              list.names <- rownames(markers.fm.list2)
            }

            size_legend <- input$Legend_size - 2


            plotdotplot <- DotPlot(sc, features = list.names) +
              RotatedAxis() +
              theme(
                axis.title.y = element_blank(),
                axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
                axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
                axis.title.x = element_blank(),
                legend.title = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
                legend.text = element_text(colour = "black", size = size_legend, family = input$font_type),
                legend.position = input$legend_position,
              ) +
              scale_colour_gradient2(low = input$low.dotplot, mid = input$middle.dotplot, high = input$high.dotplot) +
              scale_x_discrete(labels = label_wrap(20)) +
              scale_y_discrete(labels = label_wrap(20))


            file.name.clone <- paste("prioritization/Multi/Unique/", i, "_", gsub("[/]", "", gsub("&", "", name.clone)), "_dotplot_plot", "_", today(), ".png", sep = "")

            ### download the dot plot -------
            png(file.name.clone, width = input$width_png_all_expression_dotplot_top, height = input$height_png_all_expression_dotplot_top, res = input$resolution_PNG_all_expression_dotplot_top)
            plot(plotdotplot)
            dev.off()


            ##### download the OverRep ------
            geneSet <- read.csv(system.file("OverRep", "GeneSets.csv", package = "STEGO.R"), header = T)

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

            if (input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
              geneSet$GeneSet <- gsub("-", ".", geneSet$GeneSet)
            }

            if (input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
              require(stringr)
              geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
            }
            message(paste(i, "Performing Over rep analysis"))
            for (j in 1:dim(geneSet)[1]) {
              # listed GeneSet

              Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet, ";")[j])
              names(Gene.set.testing) <- "V1"
              Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
              names(Gene.set.testing2) <- "V1"
              background.overlap <- merge(Gene.set.testing2, background.genes.name, by = "V1")
              geneSet$background.geneset[j] <- length(background.overlap$V1)
              geneSet$background.geneset.name[j] <- as.character(paste(unlist(background.overlap[1]), collapse = ";"))
              # in sig gene list
              overlap <- merge(background.overlap, DEx.genes, by = "V1")

              geneSet$in.geneset[j] <- length(overlap$V1)
              geneSet$in.geneset.name[j] <- as.character(paste(unlist(overlap[1]), collapse = ";"))
            }

            geneSet2 <- subset(geneSet, geneSet$in.geneset > 0)

            message(paste(i, "has", length(geneSet2$in.geneset), "GeneSets in Unique (private)"))

            if (length(geneSet2$in.geneset) > 0) {
              for (k in 1:dim(geneSet2)[1]) {
                tota.gene.set <- geneSet2$background.geneset[k] # genes that are identified in background
                tota.gene.set
                in.geneset <- geneSet2$in.geneset[k] # DEx in geneset

                background.genes
                not.in.total <- background.genes - tota.gene.set
                not.in.geneset.sig <- total.sig - in.geneset
                d <- data.frame(gene.in.interest = c(in.geneset, not.in.geneset.sig), gene.not.interest = c(tota.gene.set, not.in.total))
                row.names(d) <- c("In_category", "not_in_category")

                if (in.geneset > 0) {
                  geneSet2$p.val[k] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
                  geneSet2$lowerCI[k] <- unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
                  geneSet2$upperCI[k] <- unlist(fisher.test(d)$conf.int)[2]
                  geneSet2$OR[k] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1], 3)
                } else {
                  geneSet2$p.value[k] <- "-"
                  geneSet2$lowerCI[k] <- "-"
                  geneSet2$upperCI[k] <- "-"
                  geneSet2$OR[k] <- "-"
                }
              }

              geneSet2 <- geneSet2[order(geneSet2$p.val, decreasing = F), ]
              geneSet2 <- subset(geneSet2, geneSet2$in.geneset >= input$in.geneset.cutoff_top)
              geneSet2 <- subset(geneSet2, geneSet2$p.val <= input$p.val_cutoff_top)
              geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
              geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")

              top.name.overrep <- paste("prioritization/Multi/Unique/", i, "_", gsub("[/]", "", gsub("&", "", name.clone)), "_OverRep", "_", today(), ".csv", sep = "")
              write.csv(geneSet2, top.name.overrep, row.names = F)
            }
          }
        }
      })
    })

    observeEvent(input$Multi_download_button, {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )
      x <- today()
      message("Downloading the Summary table...")
      top.name.clonotypes <- paste("prioritization/Multi/Overlap_Table_", x, ".csv", sep = "")
      write.csv(Upset_plot_overlap_Multi(), top.name.clonotypes, row.names = T)

      message("Downloading the Upset Plot...")
      df <- sc@meta.data
      df <- as.data.frame(df)
      unique.df <- unique(df[, names(df) %in% c(input$Samp_col, input$V_gene_sc)])
      unique.df <- unique.df %>%
        select(input$Samp_col, everything())

      names(unique.df) <- c("group", "chain")
      unique.df <- unique.df[unique.df$group %in% input$ID_Column_factor, ]
      if (length(unique(unique.df$group)) < 31) {
        top.name.clonotypes.count_png <- paste("prioritization/Multi/Overlap_Upset_plot", x, ".png", sep = "")
        png(top.name.clonotypes.count_png, width = input$width_png_TCR.UMAP_top, height = input$height_png_TCR.UMAP_top, res = input$resolution_PNG_TCR.UMAP_top)
        plot(Upset_plot_multi())
        dev.off()
      } else {
        message("More than 31 groups, you will need to create the plot manually")
      }

      message("Downloading the stacked barplot Plot...")
      top.name.clonotypes.count_png <- paste("prioritization/Multi/Stacked_bar_plot", x, ".png", sep = "")

      df4 <- TCR_Expanded()
      num_indiv <- length(unique(df4$ID_Column))
      png(top.name.clonotypes.count_png, width = (num_indiv * 100 + 600), height = input$height_png_TCR.UMAP_top, res = input$resolution_PNG_TCR.UMAP_top)
      plot(clonal_plot_multi())
      dev.off()

      if (input$Download_public_overlapping) {
        message("Downloading Multi bar plots...")

        Top_clonotypes_multiCounts_barplot()
        message("Downloading Multi Marker and OverRep analysis ")
        if (input$Download_public_overlapping_bar) {
          top_clone_FindMaker_looped_Multi()
        }
      }

      if (input$Download_private_overlapping) {
        message("Downloading Private Marker and OverRep analysis summary table ...")
        top.name.clonotypes <- paste("prioritization/Multi/Unique_Table_", x, ".csv", sep = "")
        write.csv(Top_clonotypes_Private(), top.name.clonotypes, row.names = T)

        message("Downloading Private Marker and OverRep analysis ")
        top_clone_FindMaker_looped_Private()
      }

      message("Finished")
    })

    # add in upset plot and table to show clonal expansion -> then make expansion overall with consideration to the sample. Also add in count, not just percentage.
    output$Number_of_clonotypes_to_ <- renderPrint({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload"
        )
      )
      x <- today()
      req(input$Samp_col, input$V_gene_sc, input$cut.off_percent_repMulti)

      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[, names(df3.meta) %in% input$V_gene_sc]

      BD_sum <- Top_clonotypes_multiCounts()
      dim(BD_sum)
      if (dim(BD_sum)[1] > 1) {
        BD_sum$obs <- 1
        BD_sum <- subset(BD_sum, BD_sum$priority < input$cut.off_percent_repMulti)
        observations <- sum(BD_sum$obs)
        df1 <- paste("The analysis will be limited to the top", observations, "Public/MultiSample Clones")
      } else {
        df1 <- paste("There are no public/multisample clonotypes")
      }

      BD_sum <- Top_clonotypes_Private()
      BD_sum$obs <- 1
      observations2 <- sum(BD_sum$obs)

      df2 <- paste("The analysis will be limited to the top", observations2, "Private Clones")
      df3 <- rbind(df1, df2)
      rownames(df3) <- c(1, 2)
      df3
    })

    output$private_clonotypes <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      BD_sum <- Top_clonotypes_multiCounts()
      req(BD_sum)
      BD_sum
    })

    output$PriorClono_Tab <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload"
        )
      )

      x <- today()
      req(input$Samp_col, input$V_gene_sc, input$cut.off_percent_repMulti)

      df3.meta <- sc@meta.data
      df3.meta$cluster_name <- df3.meta[, names(df3.meta) %in% input$V_gene_sc]

      BD_sum <- Top_clonotypes_multiCounts()
      req(BD_sum)

      if (dim(BD_sum)[1] > 0) {
        BD_sum$obs <- 1
        BD_sum2 <- BD_sum[, names(BD_sum) %in% c("priority", "obs")]
        total.condition <- as.data.frame(ddply(BD_sum2, "priority", numcolwise(sum)))
        total.condition <- total.condition %>% mutate(cum_freq = cumsum(obs))
        total.condition
      } else {
        as.data.frame("There are no public/multisample clonotypes")
      }
    })

    # Clustering priority  -----
    observeEvent(input$ClusterDownload_download_buttonOneOne, {
      x <- today()
      if (length(AG_cluster()) > 0) {
        message("Downloading AG cluster table...")
        Exp_stats_cutoff_count.name <- paste("prioritization/Clustering/Cluster_summary_table_AG_", x, ".csv", sep = "")
        AG_cluster <- AG_cluster()

        write.csv(AG_cluster, Exp_stats_cutoff_count.name, row.names = F)
        message("Downloading AG cluster analysis...")
        req(input$priority_cutoff)
        ggPlotUMAPClusterAG()
      } else {

      }
    })
    observeEvent(input$ClusterDownload_download_buttonOneOne, {
      x <- today()

      if (length(BD_cluster()) > 0) {
        message("Downloading BD cluster table...")
        Exp_stats_cutoff_count.name <- paste("prioritization/Clustering/Cluster_summary_table_BD_", x, ".csv", sep = "")
        BD_cluster <- BD_cluster()

        write.csv(BD_cluster, Exp_stats_cutoff_count.name, row.names = F)
        message("Downloading BD cluster analysis...")
        ggPlotUMAPClusterBD()
      } else {
      }
    })

    output$number_clusters_to_analyse_AG <- renderPrint({
      df1 <- AG_cluster()
      validate(
        need(
          nrow(df1) > 0,
          "Upload ClusTCR file"
        )
      )
      df1 <- subset(df1, df1$Sample_count >= input$Sample_count_cluster)
      df1 <- subset(df1, df1$Total_cloneCount >= input$Total_cloneCount_cluster)
      df1 <- subset(df1, df1$priority < input$priority_cutoff)
      len <- length(unique(df1$Updated_order))
      df2 <- paste("The analysis will be limited to the top", len, "AG cluster(s)")
      df3 <- paste("The cut-off used min clone clount:", min(df1$Total_cloneCount), "The cut-off used min sample clount:", min(df1$Sample_count))
      df4 <- rbind(df2, df3)
      rownames(df4) <- c(1, 2)
      df4
    })

    output$number_clusters_to_analyse_BD <- renderPrint({
      df1 <- BD_cluster()
      validate(
        need(
          nrow(df1) > 0,
          "Upload ClusTCR file"
        )
      )
      df1 <- subset(df1, df1$Sample_count >= input$Sample_count_clusterBD)
      df1 <- subset(df1, df1$Total_cloneCount >= input$Total_cloneCount_clusterBD)
      df1 <- subset(df1, df1$priority < input$priority_cutoffBD)
      len <- length(unique(df1$Updated_order))
      df2 <- paste("The analysis will be limited to the top", len, "BD cluster(s)")
      df3 <- paste("The cut-off used min clone clount:", min(df1$Total_cloneCount), "The cut-off used min sample clount:", min(df1$Sample_count))
      df4 <- rbind(df2, df3)
      rownames(df4) <- c(1, 2)
      df4
    })

    ggPlotUMAPClusterAG <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload Files"
        )
      )
      md <- sc@meta.data
      clusterAG2 <- AG_cluster()
      validate(
        need(
          nrow(clusterAG2) > 0,
          "Upload clusTCR AG table"
        )
      )
      req(clusterAG2, input$Clusters_to_dis_PIE, input$Colour_By_this, input$priority_cutoff)
      clusterAG <- subset(clusterAG2, clusterAG2$priority < input$priority_cutoff)
      clusterAG <- subset(clusterAG, clusterAG$Sample_count >= input$Sample_count_cluster)
      clusterAG <- subset(clusterAG, clusterAG$Total_cloneCount >= input$Total_cloneCount_cluster)

      len.order <- (unique(clusterAG$Updated_order))
      # clusterAG <- subset(cluster,cluster$priority<input$priority_cutoff)
      withProgress(message = "Performing AG cluster analysis", value = 0, {
        for (i in len.order) {
          incProgress(1 / length(len.order), detail = paste("AG cluster", i, "of", max(len.order)))

          ## ggplot UMAP -----
          cluster <- clusterAG
          cluster$ID_Column <- cluster[, names(cluster) %in% input$Samp_col]
          # names(cluster)[names(cluster) %in% _cluster] <- "ID_Column"

          cluster <- cluster[cluster$Updated_order %in% i, ]
          Vgene <- as.data.frame(strsplit(cluster$CDR3_Vgene, "_"))[2, 1]
          Vgene <- gsub("/", "", Vgene)
          sample_size <- unique(cluster$Sample_count)
          dirName <- paste0("prioritization/Clustering/AG/", i, "_", Vgene, "_", sample_size, "/")
          # dirName <- paste0("prioritization/Clustering/AG/")

          if (file.exists(dirName)) {
            cluster$colour <- cluster[, names(cluster) %in% input$Colour_By_this]
            cluster$colour <- gsub("_", " ", cluster$colour)
            cluster$colour <- factor(cluster$colour, levels = unique(cluster$colour))
            cluster$colour <- gsub("NA", NA, cluster$colour)

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
            } else {

            }


            col.df$col <- colorblind_vector

            figure <- ggplot(data = cluster, aes(x = UMAP_1, UMAP_2, colour = colour)) +
              geom_point(size = input$size.dot.umap) +
              scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 20), values = col.df$col, na.value = input$NA_col_analysis) +
              theme_bw() +
              theme(
                legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
                legend.title = element_blank(),
                legend.position = input$legend_position,
                # strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
                axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
                axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
                axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 0),
                axis.title.x = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
              )

            message(paste(i, " Downloading the count UMAP"))
            x <- today()
            top.name.clonotypes.top_png <- paste(dirName, i, "_", Vgene, "_AG_cluster_UMAP_", x, ".png", sep = "")
            png(top.name.clonotypes.top_png, width = input$width_png_TCR.UMAP, height = input$height_png_TCR.UMAP, res = input$resolution_PNG_TCR.UMAP)
            plot(figure)
            dev.off()
          } else {
            dir.create(dirName)
            cluster$colour <- cluster[, names(cluster) %in% input$Colour_By_this]
            cluster$colour <- gsub("_", " ", cluster$colour)
            cluster$colour <- factor(cluster$colour, levels = unique(cluster$colour))
            cluster$colour <- gsub("NA", NA, cluster$colour)

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
            } else {

            }


            col.df$col <- colorblind_vector

            figure <- ggplot(data = cluster, aes(x = UMAP_1, UMAP_2, colour = colour)) +
              geom_point(size = input$size.dot.umap) +
              scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 20), values = col.df$col, na.value = input$NA_col_analysis) +
              theme_bw() +
              theme(
                legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
                legend.title = element_blank(),
                legend.position = input$legend_position,
                # strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
                axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
                axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
                axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 0),
                axis.title.x = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
              )

            message(paste(i, " Downloading the count UMAP"))
            x <- today()
            top.name.clonotypes.top_png <- paste(dirName, i, "_", Vgene, "_AG_cluster_UMAP_", x, ".png", sep = "")
            png(top.name.clonotypes.top_png, width = input$width_png_TCR.UMAP, height = input$height_png_TCR.UMAP, res = input$resolution_PNG_TCR.UMAP)
            plot(figure)
            dev.off()
            ## Motif plot -----

            Network_df <- cluster[order(cluster$Updated_order), ]
            Network_df <- Network_df %>% distinct(CDR3_Vgene, .keep_all = TRUE) # make Unique
            motifplot <- Motif_from_cluster_file(Network_df, Clust_selected = i, selected_cluster_column = "Updated_order")

            message(paste(i, " Downloading motif plot"))
            top.name.clonotypes.top_png <- paste(dirName, i, "_", Vgene, "_AG_motif_", x, ".png", sep = "")
            png(top.name.clonotypes.top_png,
                width = input$width_png_Motif_ClusTCR2_cluster,
                height = input$height_png_Motif_ClusTCR2_cluster,
                res = input$resolution_PNG_Motif_ClusTCR2_cluster
            )
            plot(motifplot)
            dev.off()

            # Stats table ------
            cluster <- clusterAG

            names(cluster)[names(cluster) %in% input$Samp_col] <- "ID_Column"
            cluster <- cluster[order(cluster$Updated_order), ]

            rownames(cluster) <- cluster$Cell_Index

            checking <- cluster[, names(cluster) %in% c("Updated_order", "Cell_Index")]
            md.checking <- merge(md, checking, by = "Cell_Index", all.x = T)
            md.checking <- md.checking[order(md.checking$order), ]
            rownames(md.checking) <- md.checking$Cell_Index

            md.checking$Clust_selected <- ifelse(md.checking$Updated_order == i, i, "NS")
            md.checking$Clust_selected[is.na(md.checking$Clust_selected)] <- "NS"
            md.checking <- md.checking[order(md.checking$order), ]

            sc@meta.data <- md.checking
            Idents(object = sc) <- sc@meta.data$Clust_selected

            name.check.clust <- i
            min.pct.expression <- input$min_point_ # standard setting: 0.25
            min.logfc <- input$LogFC_ # 0.25 is standard

            markers.fm.list <- FindMarkers(sc, ident.1 = name.check.clust, min.pct = min.pct.expression, logfc.threshold = min.logfc, only.pos = TRUE)
            markers.fm.list2 <- subset(markers.fm.list, markers.fm.list$p_val_adj < input$pval.ex.filter)

            message(paste(length(markers.fm.list2$p_val_adj), "total markers for AG cluster", i))

            if (length(markers.fm.list2$p_val_adj) > 0) {
              message(paste(i, " Downloading stats table"))
              Exp_stats_cutoff_count.name <- paste(dirName, i, "_", Vgene, "_AG_cluster_statsTab_", x, ".csv", sep = "")

              write.csv(markers.fm.list2, Exp_stats_cutoff_count.name, row.names = T)

              # stats dotplot ----

              if (length(rownames(markers.fm.list2)) < 40) {
                list.names <- rownames(markers.fm.list2)
              } else {
                list.names <- rownames(markers.fm.list2)
                list.names <- list.names[1:40]
              }

              size_legend <- input$Legend_size - 2

              dotplotClust <- DotPlot(sc, features = list.names) +
                RotatedAxis() +
                theme(
                  axis.title.y = element_blank(),
                  axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
                  axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
                  axis.title.x = element_blank(),
                  legend.title = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
                  legend.text = element_text(colour = "black", size = size_legend, family = input$font_type),
                  legend.position = input$legend_position,
                ) +
                scale_colour_gradient2(low = input$low.dotplot.clust, mid = input$middle.dotplot.clust, high = input$high.dotplot.clust) +
                scale_x_discrete(labels = label_wrap(20)) +
                scale_y_discrete(labels = label_wrap(20))

              top.name.clonotypes.top_png <- paste(dirName, i, "_", Vgene, "_AG_dot.plot_", x, ".png", sep = "")
              png(top.name.clonotypes.top_png,
                  width = input$width_png_all_expression_dotplot_clust,
                  height = input$height_png_all_expression_dotplot_clust,
                  res = input$resolution_PNG_all_expression_dotplot_clust
              )
              plot(dotplotClust)
              dev.off()

              # column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_all_expression_dotplot_clust','Download PDF')),
              # column(2,numericInput("width_png_all_expression_dotplot_clust","Width of PNG", value = 2400)),
              # column(2,numericInput("height_png_all_expression_dotplot_clust","Height of PNG", value = 700)),
              # column(2,numericInput("resolution_PNG_all_expression_dotplot_clust","Resolution of PNG", value = 144)),
              # column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_all_expression_dotplot_clust','Download PNG'))

              # stats OverRep analysis ----
              geneSet <- read.csv(system.file("OverRep", "GeneSets.csv", package = "STEGO.R"), header = T)

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

              if (input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
                geneSet$GeneSet <- gsub("-", ".", geneSet$GeneSet)
              }

              if (input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
                require(stringr)
                geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
              }
              message(paste(i, "Starting OverRep analysis of Epitope "))
              for (j in 1:dim(geneSet)[1]) {
                # listed GeneSet

                Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet, ";")[j])
                names(Gene.set.testing) <- "V1"
                Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
                names(Gene.set.testing2) <- "V1"
                background.overlap <- merge(Gene.set.testing2, background.genes.name, by = "V1")
                geneSet$background.geneset[j] <- length(background.overlap$V1)
                geneSet$background.geneset.name[j] <- as.character(paste(unlist(background.overlap[1]), collapse = ";"))
                # in sig gene list
                overlap <- merge(background.overlap, DEx.genes, by = "V1")

                geneSet$in.geneset[j] <- length(overlap$V1)
                geneSet$in.geneset.name[j] <- as.character(paste(unlist(overlap[1]), collapse = ";"))
              }
              geneSet2 <- subset(geneSet, geneSet$in.geneset > 0)
              message(paste(i, "has", length(geneSet2$in.geneset), "GeneSets in Cluster"))

              if (length(geneSet2$in.geneset) > 0) {
                for (j in 1:dim(geneSet2)[1]) {
                  tota.gene.set <- geneSet2$background.geneset[j] # genes that are identified in background
                  in.geneset <- geneSet2$in.geneset[j] # DEx in geneset
                  not.in.total <- background.genes - tota.gene.set
                  not.in.geneset.sig <- total.sig - in.geneset
                  d <- data.frame(gene.in.interest = c(in.geneset, not.in.geneset.sig), gene.not.interest = c(tota.gene.set, not.in.total))
                  row.names(d) <- c("In_category", "not_in_category")

                  if (in.geneset > 0) {
                    geneSet2$p.val[j] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
                    geneSet2$lowerCI[j] <- unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
                    geneSet2$upperCI[j] <- unlist(fisher.test(d)$conf.int)[2]
                    geneSet2$OR[j] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1], 3)
                  } else {
                    geneSet2$p.value[j] <- "-"
                    geneSet2$lowerCI[j] <- "-"
                    geneSet2$upperCI[j] <- "-"
                    geneSet2$OR[j] <- "-"
                  }
                }
                geneSet2 <- geneSet2[order(geneSet2$p.val, decreasing = F), ]
                geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
                geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
                message("Downloading the Summary table...")
                top.name.clonotypes <- paste(dirName, i, "_", Vgene, "_AG_OverRep_", x, ".csv", sep = "")
                write.csv(geneSet2, top.name.clonotypes, row.names = F)
              }
            }
          }
        }
      })
    })
    ggPlotUMAPClusterBD <- reactive({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload Files"
        )
      )
      md <- sc@meta.data
      clusterBD2 <- BD_cluster()
      validate(
        need(
          nrow(clusterBD2) > 0,
          "Upload clusTCR BD table"
        )
      )
      req(clusterBD2, input$Clusters_to_dis_PIE, input$Colour_By_this, input$priority_cutoffBD)
      clusterBD <- subset(clusterBD2, clusterBD2$priority < input$priority_cutoffBD)
      clusterBD <- subset(clusterBD, clusterBD$Sample_count >= input$Sample_count_clusterBD)
      clusterBD <- subset(clusterBD, clusterBD$Total_cloneCount >= input$Total_cloneCount_clusterBD)

      len.order <- (unique(clusterBD$Updated_order))

      withProgress(message = "Performing BD cluster analysis", value = 0, {
        for (i in len.order) {
          # incProgress(1/len.order, detail = paste("BD cluster", i,"of",len.order))
          incProgress(1 / length(len.order), detail = paste("BD cluster", i, "of", max(len.order)))
          cluster <- clusterBD
          cluster$ID_Column <- cluster[, names(cluster) %in% input$Samp_col]
          cluster <- cluster[cluster$Updated_order %in% i, ]

          Vgene <- as.data.frame(strsplit(cluster$CDR3_Vgene, "_"))[2, 1]
          sample_size <- unique(cluster$Sample_count)
          dirName <- paste0("prioritization/Clustering/BD/", i, "_", Vgene, "_", sample_size, "/")
          if (file.exists(dirName)) {
            ## ggplot UMAP -----
            cluster$colour <- cluster[, names(cluster) %in% input$Colour_By_this]
            cluster$colour <- gsub("_", " ", cluster$colour)
            cluster$colour <- factor(cluster$colour, levels = unique(cluster$colour))
            cluster$colour <- gsub("NA", NA, cluster$colour)

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
            } else {

            }


            col.df$col <- colorblind_vector

            figure <- ggplot(data = cluster, aes(x = UMAP_1, UMAP_2, colour = colour)) +
              geom_point(size = input$size.dot.umap) +
              scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 20), values = col.df$col, na.value = input$NA_col_analysis) +
              theme_bw() +
              theme(
                legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
                legend.title = element_blank(),
                legend.position = input$legend_position,
                # strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
                axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
                axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
                axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 0),
                axis.title.x = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
              )

            message(paste(i, " Downloading the count UMAP"))
            x <- today()
            top.name.clonotypes.top_png <- paste(dirName, i, "_", Vgene, "_BD_cluster_UMAP_", x, ".png", sep = "")
            png(top.name.clonotypes.top_png, width = input$width_png_TCR.UMAP, height = input$height_png_TCR.UMAP, res = input$resolution_PNG_TCR.UMAP)
            plot(figure)
            dev.off()
            ## Motif plot -----
            Network_df <- cluster[order(cluster$Updated_order), ]
            Network_df <- Network_df %>% distinct(CDR3_Vgene, .keep_all = TRUE) # make Unique
            motifplot <- Motif_from_cluster_file(Network_df, Clust_selected = i, selected_cluster_column = "Updated_order")

            message(paste(i, " Downloading motif plot"))
            top.name.clonotypes.top_png <- paste(dirName, i, "_", Vgene, "_BD_motif_", x, ".png", sep = "")
            png(top.name.clonotypes.top_png,
                width = input$width_png_Motif_ClusTCR2_cluster,
                height = input$height_png_Motif_ClusTCR2_cluster,
                res = input$resolution_PNG_Motif_ClusTCR2_cluster
            )
            plot(motifplot)
            dev.off()
          } else {
            # if(file.exists(dirName)) break;
            dir.create(dirName)

            ## ggplot UMAP -----
            cluster$colour <- cluster[, names(cluster) %in% input$Colour_By_this]
            cluster$colour <- gsub("_", " ", cluster$colour)
            cluster$colour <- factor(cluster$colour, levels = unique(cluster$colour))
            cluster$colour <- gsub("NA", NA, cluster$colour)

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
            } else {

            }


            col.df$col <- colorblind_vector

            figure <- ggplot(data = cluster, aes(x = UMAP_1, UMAP_2, colour = colour)) +
              geom_point(size = input$size.dot.umap) +
              scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 20), values = col.df$col, na.value = input$NA_col_analysis) +
              theme_bw() +
              theme(
                legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
                legend.title = element_blank(),
                legend.position = input$legend_position,
                # strip.text = element_text(size = input$Strip_text_size, family = input$font_type),
                axis.title.y = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
                axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
                axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 0),
                axis.title.x = element_text(colour = "black", family = input$font_type, size = input$title.text.sizer2),
              )

            message(paste(i, " Downloading the count UMAP"))
            x <- today()
            top.name.clonotypes.top_png <- paste(dirName, i, "_", Vgene, "_BD_cluster_UMAP_", x, ".png", sep = "")
            png(top.name.clonotypes.top_png, width = input$width_png_TCR.UMAP, height = input$height_png_TCR.UMAP, res = input$resolution_PNG_TCR.UMAP)
            plot(figure)
            dev.off()
            ## Motif plot -----
            Network_df <- cluster[order(cluster$Updated_order), ]
            Network_df <- Network_df %>% distinct(CDR3_Vgene, .keep_all = TRUE) # make Unique
            motifplot <- Motif_from_cluster_file(Network_df, Clust_selected = i, selected_cluster_column = "Updated_order")

            message(paste(i, " Downloading motif plot"))
            top.name.clonotypes.top_png <- paste(dirName, i, "_", Vgene, "_BD_motif_", x, ".png", sep = "")
            png(top.name.clonotypes.top_png,
                width = input$width_png_Motif_ClusTCR2_cluster,
                height = input$height_png_Motif_ClusTCR2_cluster,
                res = input$resolution_PNG_Motif_ClusTCR2_cluster
            )
            plot(motifplot)
            dev.off()

            # Stats table ------
            cluster <- clusterBD

            names(cluster)[names(cluster) %in% input$Samp_col] <- "ID_Column"
            cluster <- cluster[order(cluster$Updated_order), ]

            rownames(cluster) <- cluster$Cell_Index

            checking <- cluster[, names(cluster) %in% c("Updated_order", "Cell_Index")]
            md.checking <- merge(md, checking, by = "Cell_Index", all.x = T)
            md.checking <- md.checking[order(md.checking$order), ]
            rownames(md.checking) <- md.checking$Cell_Index

            md.checking$Clust_selected <- ifelse(md.checking$Updated_order == i, i, "NS")
            md.checking$Clust_selected[is.na(md.checking$Clust_selected)] <- "NS"
            md.checking <- md.checking[order(md.checking$order), ]

            sc@meta.data <- md.checking
            Idents(object = sc) <- sc@meta.data$Clust_selected

            name.check.clust <- i
            min.pct.expression <- input$min_point_ # standard setting: 0.25
            min.logfc <- input$LogFC_ # 0.25 is standard

            markers.fm.list <- FindMarkers(sc, ident.1 = name.check.clust, min.pct = min.pct.expression, logfc.threshold = min.logfc, only.pos = TRUE)
            markers.fm.list2 <- subset(markers.fm.list, markers.fm.list$p_val_adj < input$pval.ex.filter)

            message(paste(i, "BD cluster has", length(markers.fm.list2$p_val_adj), "total markers"))

            if (length(markers.fm.list2$p_val_adj) > 0) {
              message(paste(i, " Downloading stats table"))
              Exp_stats_cutoff_count.name <- paste(dirName, i, "_", Vgene, "_BD_cluster_statsTab_", x, ".csv", sep = "")
              write.csv(markers.fm.list2, Exp_stats_cutoff_count.name, row.names = T)

              # stats dotplot ----

              if (length(rownames(markers.fm.list2)) < 40) {
                list.names <- rownames(markers.fm.list2)
              } else {
                list.names <- rownames(markers.fm.list2)
                list.names <- list.names[1:40]
              }

              size_legend <- input$Legend_size - 2

              dotplotClust <- DotPlot(sc, features = list.names) +
                RotatedAxis() +
                theme(
                  axis.title.y = element_blank(),
                  axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
                  axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
                  axis.title.x = element_blank(),
                  legend.title = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
                  legend.text = element_text(colour = "black", size = size_legend, family = input$font_type),
                  legend.position = input$legend_position,
                ) +
                scale_colour_gradient2(low = input$low.dotplot.clust, mid = input$middle.dotplot.clust, high = input$high.dotplot.clust) +
                scale_x_discrete(labels = label_wrap(20)) +
                scale_y_discrete(labels = label_wrap(20))

              top.name.clonotypes.top_png <- paste(dirName, i, "_", Vgene, "_BD_dot.plot_", x, ".png", sep = "")
              png(top.name.clonotypes.top_png,
                  width = input$width_png_all_expression_dotplot_clust,
                  height = input$height_png_all_expression_dotplot_clust,
                  res = input$resolution_PNG_all_expression_dotplot_clust
              )
              plot(dotplotClust)
              dev.off()

              # stats OverRep analysis ----

              if (dim(markers.fm.list2)[1] > 0) {
                geneSet <- read.csv(system.file("OverRep", "GeneSets.csv", package = "STEGO.R"), header = T)

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

                if (input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
                  geneSet$GeneSet <- gsub("-", ".", geneSet$GeneSet)
                }

                if (input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
                  require(stringr)
                  geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
                }
                message(paste("Starting OverRep analysis of cluster ", i))
                for (j in 1:dim(geneSet)[1]) {
                  # listed GeneSet

                  Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet, ";")[j])
                  names(Gene.set.testing) <- "V1"
                  Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
                  names(Gene.set.testing2) <- "V1"
                  background.overlap <- merge(Gene.set.testing2, background.genes.name, by = "V1")
                  geneSet$background.geneset[j] <- length(background.overlap$V1)
                  geneSet$background.geneset.name[j] <- as.character(paste(unlist(background.overlap[1]), collapse = ";"))
                  # in sig gene list
                  overlap <- merge(background.overlap, DEx.genes, by = "V1")

                  geneSet$in.geneset[j] <- length(overlap$V1)
                  geneSet$in.geneset.name[j] <- as.character(paste(unlist(overlap[1]), collapse = ";"))
                }

                geneSet2 <- subset(geneSet, geneSet$in.geneset > 0)
                message(paste(i, "has", length(geneSet2$in.geneset), "GeneSets in Epitope"))
                if (length(geneSet2$in.geneset) > 0) {
                  for (j in 1:dim(geneSet2)[1]) {
                    tota.gene.set <- geneSet2$background.geneset[j] # genes that are identified in background
                    in.geneset <- geneSet2$in.geneset[j] # DEx in geneset
                    not.in.total <- background.genes - tota.gene.set
                    not.in.geneset.sig <- total.sig - in.geneset
                    d <- data.frame(gene.in.interest = c(in.geneset, not.in.geneset.sig), gene.not.interest = c(tota.gene.set, not.in.total))
                    row.names(d) <- c("In_category", "not_in_category")

                    if (in.geneset > 0) {
                      geneSet2$p.val[j] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
                      geneSet2$lowerCI[j] <- unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
                      geneSet2$upperCI[j] <- unlist(fisher.test(d)$conf.int)[2]
                      geneSet2$OR[j] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1], 3)
                    } else {
                      geneSet2$p.value[j] <- "-"
                      geneSet2$lowerCI[j] <- "-"
                      geneSet2$upperCI[j] <- "-"
                      geneSet2$OR[j] <- "-"
                    }
                  }
                  geneSet2 <- geneSet2[order(geneSet2$p.val, decreasing = F), ]
                  geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
                  geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
                  message("Downloading over-rep table...")
                  top.name.clonotypes <- paste(dirName, i, "_", Vgene, "_BD_OverRep_", x, ".csv", sep = "")
                  write.csv(geneSet2, top.name.clonotypes, row.names = F)
                }
              }
            }
          }
        }
      })
    })
    output$Cluster_dowload_button_prior <- renderUI({
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )
      TCR_Expanded_Df <- TCR_Expanded_fun(sc, (input$Samp_col), (input$V_gene_sc))
      TCR_Expanded_Df$obs <- 1
      TCR_Expanded_Df2 <- subset(TCR_Expanded_Df, TCR_Expanded_Df$percent > input$cut.off_percent_rep)
      observations <- sum(TCR_Expanded_Df2$obs)
      length.samp.ID <- length(unique(TCR_Expanded_Df$ID_Column))

      fluidRow(
        # column(12,selectInput("AnalysisType","Preset parameters for", = c("ImmunoDom"))),
        column(12, actionButton("ClusterDownload_download_buttonOneOne", "Download cluster analysis", style = "color: white; background-color:#4682b4"))
      )
    })

    ### clustering table -----
    output$PriorClustTB_Tab <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload Files"
        )
      )
      df1 <- AG_cluster()
      df1 <- subset(df1, df1$priority < input$priority_cutoff)
      df1
    })


    # epitope prioritization ------
    observeEvent(input$EpitopePrior_Download_Bt, {
      x <- today()
      message(paste("Downloading sumamry table of", input$Prior_AddInEpiUI_1, "vs", input$Prior_AddInEpiUI_2))
      Exp_stats_cutoff_count.name <- paste("prioritization/EpitopePred/Epitope_summary_table_", input$Prior_AddInEpiUI_1, "vs", input$Prior_AddInEpiUI_2, "_", x, ".csv", sep = "")
      write.csv(Prior_Epitope_dt_process(), Exp_stats_cutoff_count.name, row.names = F)

      message(paste("Downloading heatmap", input$Prior_AddInEpiUI_1, " vs ", input$Prior_AddInEpiUI_2))

      name.epitope <- paste("prioritization/EpitopePred/", input$Prior_AddInEpiUI_1, "_", input$Prior_AddInEpiUI_2, "_Heatmap_", x, ".png", sep = "")

      Common_Epitope_code <- Common_Epitope_code()
      total_X <- length(unique(Common_Epitope_code$Selected_function))
      total_Y <- length(unique(Common_Epitope_code$Selected_group))

      png(name.epitope, width = 20 * total_Y + 400, height = 20 * total_X + 400, res = 144)
      plot(Prior_heatmap_epitope())
      dev.off()

      message(paste("Downloading UMAP", input$Prior_AddInEpiUI_1, " vs ", input$Prior_AddInEpiUI_2))
      name.epitope.UMAP <- paste("prioritization/EpitopePred/", "Epitope_UMAP_coloredby.", input$Prior_AddInEpiUI_1, "_", x, ".png", sep = "")

      len.colour <- length(unique(Common_Epitope_code$ID_Column))
      nCol <- round(sqrt(len.colour), 0)

      png(name.epitope.UMAP, width = 400 * nCol + 400, height = 400 * nCol, res = 144)
      plot(Prior_UMAP_Epitope())
      dev.off()

      message(paste("Starting Stats download..."))
      Prior_Stats_Epitope()
    })

    # epitope code for automated analysis ------
    df_tcrex2 <- reactive({
      epi <- data_sc_TCRex()
      req(epi)
      epi[!(duplicated(epi$CDR3_beta)), ]
    })
    # UI Prior Epi
    output$AddInEpiUI_1 <- renderUI({
      sc <- UMAP_metadata_with_labs()

      if (input$epitope_uploaded) {
        epi <- data_sc_TCRex()
        validate(
          need(
            nrow(sc) > 0 & nrow(epi) > 0,
            "Upload Files"
          )
        )
        df3.meta <- sc@meta.data
        if (input$datasource == "BD_Rhapsody_Paired") {
          df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
        } else if (input$datasource == "BD_Rhapsody_AIRR") {
          df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
        } else {
          df3.meta$CDR3_beta <- df3.meta$cdr3_BD
        }
        req(df3.meta$CDR3_beta)
        epi$beta <- epi$CDR3_beta

        df3.meta <- merge(df3.meta, epi, by = "CDR3_beta")
        selectInput("Prior_AddInEpiUI_1", "Colour Pie by (hm = y-axis):", choices = names(df3.meta), selected = "pathology")
      } else {
        df3.meta <- sc@meta.data
        selectInput("Prior_AddInEpiUI_1", "Colour Pie by (hm = y-axis):", choices = names(df3.meta), selected = "Sample_Name")
      }
    })
    output$AddInEpiUI_2 <- renderUI({
      sc <- UMAP_metadata_with_labs()

      if (input$epitope_uploaded) {
        epi <- data_sc_TCRex()
        validate(
          need(
            nrow(sc) > 0 & nrow(epi) > 0,
            "Upload Files"
          )
        )
        df3.meta <- sc@meta.data

        if (input$datasource == "BD_Rhapsody_Paired") {
          df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
        } else if (input$datasource == "BD_Rhapsody_AIRR") {
          df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
        } else {
          df3.meta$CDR3_beta <- df3.meta$cdr3_BD
        }
        req(df3.meta$CDR3_beta)
        epi$beta <- epi$CDR3_beta
        df3.meta <- merge(df3.meta, epi, by = "CDR3_beta")
        selectInput("Prior_AddInEpiUI_2", "Split Pie by (hm = x-axis):", choices = names(df3.meta), selected = "epitope")
      } else {
        df3.meta <- sc@meta.data
        selectInput("Prior_AddInEpiUI_2", "Split Pie by (hm = x-axis):", choices = names(df3.meta), selected = "Tcellfunction")
      }
    })

    Common_Epitope_code <- reactive({
      sc <- UMAP_metadata_with_labs()

      if (input$epitope_uploaded) {
        epi <- df_tcrex2()
        validate(
          need(
            nrow(sc) > 0 & nrow(epi) > 0,
            "Upload Files"
          )
        )
        df3.meta <- sc@meta.data

        if (input$datasource == "BD_Rhapsody_Paired") {
          df3.meta$CDR3_beta <- paste("C", df3.meta$cdr3_BD, "F", sep = "")
        } else if (input$datasource == "BD_Rhapsody_AIRR") {
          df3.meta$CDR3_beta <- df3.meta$junction_aa_BD
        } else {
          df3.meta$CDR3_beta <- df3.meta$cdr3_BD
        }

        req(df3.meta$CDR3_beta)

        epi$beta <- epi$CDR3_beta
        df3.meta <- merge(df3.meta, epi, by = "CDR3_beta", all.x = T)
      } else {
        df3.meta <- sc@meta.data
      }

      df3.meta$Selected_function <- df3.meta[, names(df3.meta) %in% input$Prior_AddInEpiUI_1]

      top_BD_cluster <- df3.meta
      top_BD_cluster$Selected_group <- top_BD_cluster[, names(top_BD_cluster) %in% input$Prior_AddInEpiUI_2]
      top_BD_cluster$Selected_indiv <- top_BD_cluster[, names(top_BD_cluster) %in% input$Samp_col]

      top_BD_cluster$Selected_Function_group <- paste(top_BD_cluster$Selected_function, top_BD_cluster$Selected_group, sep = "_") # merged both
      top_BD_cluster$colouring_col <- top_BD_cluster[, names(top_BD_cluster) %in% input$Colour_By_this]
      top_BD_cluster$ID_Column <- top_BD_cluster[, names(top_BD_cluster) %in% input$Split_group_by_]

      top_BD_cluster$cloneCount <- 1
      rownames(top_BD_cluster) <- top_BD_cluster$Cell_Index
      top_BD_cluster
    })

    Prior_Epitope_dt_process <- reactive({
      top_BD_cluster <- Common_Epitope_code()
      top_BD_cluster2 <- top_BD_cluster[, names(top_BD_cluster) %in% c("cloneCount", "Selected_function", "Selected_group", "Selected_indiv", "Selected_Function_group")]
      top_BD_cluster2 <- top_BD_cluster2 %>%
        select(cloneCount, everything())
      df2 <- as.data.frame(ddply(top_BD_cluster2, names(top_BD_cluster2)[-c(1)], numcolwise(sum)))
      df2 <- df2[complete.cases(df2) == T, ]

      df2 <- subset(df2, df2$Selected_function != "BackgroundData")
      df2 <- subset(df2, df2$Selected_function != "Background")
      df2 <- subset(df2, df2$Selected_group != "BackgroundData")
      df2 <- subset(df2, df2$Selected_group != "Background")

      df2$fraction <- df2$cloneCount / sum(df2$cloneCount)
      df2$Percent <- df2$cloneCount / sum(df2$cloneCount) * 100
      df2$Selected_function <- ifelse(grepl("NA", df2$Selected_function), "-", df2$Selected_function)
      df2[order(df2$Percent, decreasing = T), ]
    })

    Prior_heatmap_epitope <- reactive({
      df3.meta <- Common_Epitope_code()
      req(df3.meta)
      # df3.meta$Selected_function <- df3.meta[,names(df3.meta) %in% input$Prior_AddInEpiUI_1]
      # top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% input$Prior_AddInEpiUI_2]
      df3.meta2 <- df3.meta[!is.na(df3.meta$Selected_function), ]
      df3.meta2 <- subset(df3.meta2, df3.meta2$Selected_function != "BackgroundData")
      df3.meta2 <- subset(df3.meta2, df3.meta2$Selected_function != "Background")
      df3.meta3 <- df3.meta2[!is.na(df3.meta2$Selected_group), ]
      df3.meta3 <- subset(df3.meta3, df3.meta3$Selected_group != "BackgroundData")
      df3.meta3 <- subset(df3.meta3, df3.meta3$Selected_group != "Background")


      df.1 <- acast(df3.meta3, Selected_function ~ Selected_group, value.var = "cloneCount")
      df.1
      head(df.1)
      df.1[is.na(df.1)] <- 0
      dim(df.1)
      min.FC <- min(df.1)
      med.FC <- max(df.1) / 2
      max.FC <- max(df.1)

      ht <- Heatmap(df.1,
                    heatmap_legend_param = list(
                      title = "Unique\nTCR's",
                      title_gp = gpar(
                        fontsize = input$Legend_size,
                        fontface = "bold", fontfamily = input$font_type
                      ),
                      labels_gp = gpar(fontsize = input$Legend_size, fontfamily = input$font_type)
                    ),
                    col = colorRamp2(c(min.FC,med.FC, max.FC), c(input$min_FC_col, input$med_FC_col, input$max_FC_col)),
                    row_names_gp = grid::gpar(fontsize = input$text_size, fontfamily = input$font_type),
                    column_names_gp = grid::gpar(fontsize = input$text_size, fontfamily = input$font_type),
      )

      draw(ht, padding = unit(c(10, 10, 10, 10), "mm"))
    })

    Prior_UMAP_Epitope <- reactive({
      df3.meta <- Common_Epitope_code()
      req(df3.meta)
      df3.meta$selected <- df3.meta$Selected_Function_group
      df3.meta$selected <- ifelse(grepl("Background", df3.meta$selected), "Background", df3.meta$selected)
      # df3.meta <- subset(df3.meta,df3.meta$selected != "BackgroundData")
      df3.meta <- subset(df3.meta, df3.meta$selected != "Background")


      df3.meta <- df3.meta[order(df3.meta$selected, decreasing = F), ]
      df3.meta$selected <- factor(df3.meta$selected, levels = unique(df3.meta$selected))
      df3.meta$ID_Column <- factor(df3.meta$ID_Column, levels = unique(df3.meta$ID_Column))


      num <- as.data.frame(unique(df3.meta$selected))
      num <- as.data.frame(num[complete.cases(num) == T, ])

      col.df <- as.data.frame(unique(df3.meta$selected))

      num <- (length(unlist(col.df)))
      num <- num[!is.na(num)]

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
      } else {

      }


      col.df$col <- colorblind_vector



      df <- ggplot(data = df3.meta, aes(x = UMAP_1, UMAP_2, colour = selected, size = selected, alpha = selected)) +
        geom_point() +
        scale_color_manual(na.value = input$NA_col_analysis, values = col.df$col) +
        scale_size_manual(na.value = 0.25, values = rep(2, num)) +
        scale_alpha_manual(na.value = 0.25, values = rep(1, num)) +
        theme(
          legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        ) +
        theme_bw()

      len.colour <- length(unique(df3.meta$ID_Column))
      nCol <- round(sqrt(len.colour), 0)

      df <- df + facet_wrap(~ID_Column, ncol = nCol)

      df
    })

    Prior_Stats_Epitope <- reactive({
      x <- today()
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload Files"
        )
      )
      md <- Common_Epitope_code()
      req(md)
      top_BD_cluster2 <- md[, names(md) %in% c("cloneCount", "Selected_function", "Selected_group", "Selected_Function_group")]
      top_BD_cluster2 <- top_BD_cluster2 %>% select(cloneCount, everything())
      df2 <- as.data.frame(ddply(top_BD_cluster2, c("Selected_function", "Selected_group", "Selected_Function_group"), numcolwise(sum)))
      df2$selected <- df2$Selected_Function_group
      df2$selected <- ifelse(grepl("Background", df2$selected), "Background", df2$selected)

      df2 <- df2[complete.cases(df2) == T, ]
      df2 <- subset(df2, df2$selected != "Background")

      df2$fraction <- df2$cloneCount / sum(df2$cloneCount)
      df2$Percent <- df2$cloneCount / sum(df2$cloneCount) * 100
      df2$Selected_function <- ifelse(grepl("NA", df2$Selected_function), "-", df2$Selected_function)
      df2 <- subset(df2, df2$cloneCount > 2)
      df2 <- df2[order(df2$Percent, decreasing = T), ]

      Selected_group <- (df2$Selected_Function_group)

      len.order <- length(Selected_group)
      if (length(df2$Percent) > 0) {
        withProgress(message = "Peforming stats on selected Epitope", value = 0, {
          for (i in 1:len.order) {
            ##### UMAP ------
            message(Selected_group[i])
            incProgress(1 / len.order, detail = paste("Epitope", i, "of", len.order))

            Epi.selected <- subset(md, md$Selected_Function_group == df2$Selected_Function_group[i])

            # dirName <- paste0("prioritization/Clustering/BD/",i,"_",Vgene,"_",sample_size,"/")

            selected.epitope.name <- df2$Selected_group[i]
            Epi.selected.function.name <- df2$Selected_function[i]

            top.name.clonotypes.top_csv <- paste("prioritization/EpitopePred/", i, "_", selected.epitope.name, "_", Epi.selected.function.name, "_Epitope_summary_table_", x, ".csv", sep = "")

            write.csv(Epi.selected, top.name.clonotypes.top_csv)

            num <- as.data.frame(unique(Epi.selected$colouring_col))
            num <- as.data.frame(num[complete.cases(num) == T, ])

            col.df <- num

            num <- (length(unlist(col.df)))
            num <- num[!is.na(num)]

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
            } else {

            }


            col.df$col <- colorblind_vector

            df <- ggplot(data = Epi.selected, aes(x = UMAP_1, UMAP_2, colour = colouring_col, size = colouring_col, alpha = colouring_col)) +
              geom_point() +
              scale_color_manual(na.value = input$NA_col_analysis, values = col.df$col) +
              scale_size_manual(na.value = 0.25, values = rep(2, num)) +
              scale_alpha_manual(na.value = 0.25, values = rep(1, num)) +
              theme(
                legend.text = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
                legend.title = element_blank(),
                legend.position = input$legend_position,
              ) +
              theme_bw()

            message(paste(i, " Downloading the count UMAP"))
            x <- today()
            top.name.clonotypes.top_png <- paste("prioritization/EpitopePred/", i, "_", selected.epitope.name, "_", Epi.selected.function.name, "_Epitope_UMAP_", x, ".png", sep = "")
            png(top.name.clonotypes.top_png, width = input$width_png_TCR.UMAP, height = input$height_png_TCR.UMAP, res = input$resolution_PNG_TCR.UMAP)
            plot(df)
            dev.off()

            # Stats table ------
            md.checking <- Common_Epitope_code()
            md.checking <- md.checking[order(md.checking$order), ]
            rownames(md.checking) <- md.checking$Cell_Index


            md.checking$selected <- ifelse(md.checking$Selected_Function_group == df2$Selected_Function_group[i], df2$Selected_Function_group[i], "NS")
            md.checking$selected[is.na(md.checking$selected)] <- "NS"
            md.checking <- md.checking[order(md.checking$order), ]
            md.checking
            sc@meta.data <- md.checking
            Idents(object = sc) <- sc@meta.data$selected

            name.check.epi <- df2$Selected_Function_group[i]
            min.pct.expression <- input$min_point_ # standard setting: 0.25
            min.logfc <- input$LogFC_ # 0.25 is standard

            markers.fm.list <- FindMarkers(sc, ident.1 = name.check.epi, min.pct = min.pct.expression, logfc.threshold = min.logfc, only.pos = TRUE)
            markers.fm.list2 <- subset(markers.fm.list, markers.fm.list$p_val_adj < input$pval.ex.filter)

            message(paste(i, selected.epitope.name, "has", length(markers.fm.list2$p_val_adj), "total markers"))

            if (length(markers.fm.list2$p_val_adj) > 0) {
              message(paste(i, " Downloading stats table"))
              Exp_stats_cutoff_count.name <- paste("prioritization/EpitopePred/", i, "_", selected.epitope.name, "_", Epi.selected.function.name, "_Epitope_statsTab_", x, ".csv", sep = "")
              write.csv(markers.fm.list2, Exp_stats_cutoff_count.name, row.names = T)

              #   # stats dotplot ----
              if (length(rownames(markers.fm.list2)) < 40) {
                list.names <- rownames(markers.fm.list2)
              } else {
                list.names <- rownames(markers.fm.list2)
                list.names <- list.names[1:40]
              }

              size_legend <- input$Legend_size - 2

              dotplotClust <- DotPlot(sc, features = list.names) +
                RotatedAxis() +
                theme(
                  axis.title.y = element_blank(),
                  axis.text.y = element_text(colour = "black", family = input$font_type, size = input$text_size),
                  axis.text.x = element_text(colour = "black", family = input$font_type, size = input$text_size, angle = 90),
                  axis.title.x = element_blank(),
                  legend.title = element_text(colour = "black", size = input$Legend_size, family = input$font_type),
                  legend.text = element_text(colour = "black", size = size_legend, family = input$font_type),
                  legend.position = input$legend_position,
                ) +
                scale_colour_gradient2(low = input$low.dotplot.clust, mid = input$middle.dotplot.clust, high = input$high.dotplot.clust) +
                scale_x_discrete(labels = label_wrap(20)) +
                scale_y_discrete(labels = label_wrap(20))

              top.name.clonotypes.top_png <- paste("prioritization/EpitopePred/", i, "_", selected.epitope.name, "_", Epi.selected.function.name, "_Epitope_dot.plot_", x, ".png", sep = "")
              png(top.name.clonotypes.top_png,
                  width = input$width_png_all_expression_dotplot_clust,
                  height = input$height_png_all_expression_dotplot_clust,
                  res = input$resolution_PNG_all_expression_dotplot_clust
              )
              plot(dotplotClust)
              dev.off()
              # stats OverRep analysis ----
              if (dim(markers.fm.list2)[1] > 0) {
                geneSet <- read.csv(system.file("OverRep", "GeneSets.csv", package = "STEGO.R"), header = T)

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

                if (input$datasource == "BD_Rhapsody_Paired" || input$datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
                  geneSet$GeneSet <- gsub("-", ".", geneSet$GeneSet)
                }

                if (input$species_analysis == "mm") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
                  require(stringr)
                  geneSet$GeneSet <- str_to_title(geneSet$GeneSet)
                }
                message(paste("Starting OverRep analysis of cluster ", i))
                for (j in 1:dim(geneSet)[1]) {
                  # listed GeneSet

                  Gene.set.testing <- as.data.frame(strsplit(geneSet$GeneSet, ";")[j])
                  names(Gene.set.testing) <- "V1"
                  Gene.set.testing2 <- as.data.frame(unique(Gene.set.testing$V1))
                  names(Gene.set.testing2) <- "V1"
                  background.overlap <- merge(Gene.set.testing2, background.genes.name, by = "V1")
                  geneSet$background.geneset[j] <- length(background.overlap$V1)
                  geneSet$background.geneset.name[j] <- as.character(paste(unlist(background.overlap[1]), collapse = ";"))
                  # in sig gene list
                  overlap <- merge(background.overlap, DEx.genes, by = "V1")

                  geneSet$in.geneset[j] <- length(overlap$V1)
                  geneSet$in.geneset.name[j] <- as.character(paste(unlist(overlap[1]), collapse = ";"))
                }

                geneSet2 <- subset(geneSet, geneSet$in.geneset > 0)
                message(paste("there are", length(geneSet2$in.geneset), "GeneSets in Epitope", i))

                if (length(geneSet2$in.geneset) > 0) {
                  for (j in 1:dim(geneSet2)[1]) {
                    tota.gene.set <- geneSet2$background.geneset[j] # genes that are identified in background
                    in.geneset <- geneSet2$in.geneset[j] # DEx in geneset
                    not.in.total <- background.genes - tota.gene.set
                    not.in.geneset.sig <- total.sig - in.geneset
                    d <- data.frame(gene.in.interest = c(in.geneset, not.in.geneset.sig), gene.not.interest = c(tota.gene.set, not.in.total))
                    row.names(d) <- c("In_category", "not_in_category")

                    if (in.geneset > 0) {
                      geneSet2$p.val[j] <- unlist(fisher.test(d, alternative = "greater")$p.value)[1]
                      geneSet2$lowerCI[j] <- unlist(fisher.test(d, alternative = "greater")$conf.int)[1]
                      geneSet2$upperCI[j] <- unlist(fisher.test(d)$conf.int)[2]
                      geneSet2$OR[j] <- round(unlist(fisher.test(d, alternative = "greater")$estimate)[1], 3)
                    } else {
                      geneSet2$p.value[j] <- "-"
                      geneSet2$lowerCI[j] <- "-"
                      geneSet2$upperCI[j] <- "-"
                      geneSet2$OR[j] <- "-"
                    }
                  }
                  geneSet2 <- geneSet2[order(geneSet2$p.val, decreasing = F), ]
                  geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
                  geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")
                  message("Downloading the Summary table...")
                  # top.name.clonotypes.top_png <- paste("prioritization/EpitopePred/",i,"_",selected.epitope.name,"_dot.plot_",x,".png",sep="")
                  top.name.clonotypes <- paste("prioritization/EpitopePred/", i, "_", selected.epitope.name, "_", Epi.selected.function.name, "_Epitope_OverRep_", x, ".csv", sep = "")
                  write.csv(geneSet2, top.name.clonotypes, row.names = F)
                }
              }
            } # loop to find
            ######
          }
        })
      } else {
        message("No epitopes to analyse")
      }
    })

    output$Test_table_1 <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload Files"
        )
      )
      md <- Common_Epitope_code()
      req(md)
      top_BD_cluster2 <- md[, names(md) %in% c("cloneCount", "Selected_function", "Selected_group", "Selected_Function_group")]
      top_BD_cluster2 <- top_BD_cluster2 %>% select(cloneCount, everything())
      df2 <- as.data.frame(ddply(top_BD_cluster2, c("Selected_function", "Selected_group", "Selected_Function_group"), numcolwise(sum)))
      df2 <- df2[complete.cases(df2) == T, ]
      df2$Selected_Function_group <- ifelse(grepl("Background", df2$Selected_Function_group), "Background", df2$Selected_Function_group)
      df2 <- subset(df2, df2$Selected_Function_group != "Background")
      df2$fraction <- df2$cloneCount / sum(df2$cloneCount)
      df2$Percent <- df2$cloneCount / sum(df2$cloneCount) * 100
      df2$Selected_function <- ifelse(grepl("NA", df2$Selected_function), "-", df2$Selected_function)
      df2 <- subset(df2, df2$cloneCount > 2)
      df2 <- df2[order(df2$Percent, decreasing = T), ]
      df2
    })
    # epitope prioritization ------
    observeEvent(input$EpitopePrior_Download_Bt, {
      x <- today()
      message(paste("Downloading sumamry table of", input$Prior_AddInEpiUI_1, "vs", input$Prior_AddInEpiUI_2))
      Exp_stats_cutoff_count.name <- paste("prioritization/EpitopePred/Epitope_summary_table_", input$Prior_AddInEpiUI_1, "vs", input$Prior_AddInEpiUI_2, "_", x, ".csv", sep = "")
      write.csv(Prior_Epitope_dt_process(), Exp_stats_cutoff_count.name, row.names = F)

      message(paste("Downloading heatmap", input$Prior_AddInEpiUI_1, " vs ", input$Prior_AddInEpiUI_2))

      name.epitope <- paste("prioritization/EpitopePred/", input$Prior_AddInEpiUI_1, "_", input$Prior_AddInEpiUI_2, "_Heatmap_", x, ".png", sep = "")

      Common_Epitope_code <- Common_Epitope_code()
      total_X <- length(unique(Common_Epitope_code$Selected_function))
      total_Y <- length(unique(Common_Epitope_code$Selected_group))

      png(name.epitope, width = 20 * total_Y + 600, height = 20 * total_X + 400, res = 144)
      plot(Prior_heatmap_epitope())
      dev.off()

      message(paste("Downloading UMAP", input$Prior_AddInEpiUI_1, " vs ", input$Prior_AddInEpiUI_2))
      name.epitope.UMAP <- paste("prioritization/EpitopePred/", "Epitope_UMAP_coloredby.", input$Prior_AddInEpiUI_1, "_", x, ".png", sep = "")

      len.colour <- length(unique(Common_Epitope_code$ID_Column))
      nCol <- round(sqrt(len.colour), 0)

      png(name.epitope.UMAP, width = 400 * nCol + 400, height = 400 * nCol, res = 144)
      plot(Prior_UMAP_Epitope())
      dev.off()

      message(paste("Starting Stats download..."))
      Prior_Stats_Epitope()
    })
    ### test table -----
    output$colors.top_dt <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 10, scrollX = TRUE), {
      sc <- UMAP_metadata_with_labs()
      validate(
        need(
          nrow(sc) > 0,
          "Upload Files"
        )
      )
      md <- sc@meta.data
      clusterBD2 <- BD_cluster()
      validate(
        need(
          nrow(clusterBD2) > 0,
          "Upload clusTCR BD table"
        )
      )
      clusterBD2
      cluster <- clusterBD2

      names(cluster)[names(cluster) %in% input$Samp_col] <- "ID_Column"
      cluster <- cluster[order(cluster$Updated_order), ]
      rownames(cluster) <- cluster$Cell_Index

      checking <- cluster[, names(cluster) %in% c("Updated_order", "Cell_Index")]
      checking
      md.checking <- merge(md, checking, by = "Cell_Index", all.x = T)
      md.checking <- md.checking[order(md.checking$order), ]
      rownames(md.checking) <- md.checking$Cell_Index

      md.checking$Clust_selected <- ifelse(md.checking$Updated_order == 10, 10, "NS")
      md.checking$Clust_selected[is.na(md.checking$Clust_selected)] <- "NS"
      md.checking <- md.checking[order(md.checking$order), ]
      md.checking
      sc@meta.data <- md.checking
      Idents(object = sc) <- sc@meta.data$Clust_selected

      name.check.clust <- 10
      min.pct.expression <- input$min_point_ # standard setting: 0.25
      min.logfc <- input$LogFC_ # 0.25 is standard

      markers.fm.list <- FindMarkers(sc, ident.1 = name.check.clust, min.pct = min.pct.expression, logfc.threshold = min.logfc, only.pos = TRUE)
      markers.fm.list2 <- subset(markers.fm.list, markers.fm.list$p_val_adj < input$pval.ex.filter)
      as.data.frame(length(markers.fm.list2$p_val_adj))
    })
    ### end -----
    ### pgen -----

    observeEvent(input$load_olga,{
      require(reticulate)

      if(py_module_available("olga")) {
        message("OLGA is installed")

        updateSelectInput(
          session,
          "Olga_installed",
          choices = "installed"
        )

      } else {
        message("Please Install OLGA, if you want to do the post analysis section")
        updateSelectInput(
          session,
          "Olga_installed",
          choices = "Unavailable"
        )


      }
      hide(id = "load_olga")
      hide(id = "Olga_installed")

    })


    OLGA_data <- reactive({
      inFile_sc_OLGA <- input$file1_rds_OLGA
      if (is.null(inFile_sc_OLGA)) {
        return(NULL)
      }

      LoadSeuratRds(inFile_sc_OLGA$datapath)
    })

    observe({
      sc <- OLGA_data()
      validate(
        need(
          nrow(sc) > 0,
          error_message_val1
        )
      )
      meta.data <- sc@meta.data
      updateSelectInput(
        session,
        "Samp_col_OLGA",
        choices = names(meta.data),
        selected = "Sample_Name"
      )
    })


    output$Pgen_Selected <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      sc <- OLGA_data()
      validate(
        need(
          nrow(sc) > 0,
          "Upload .RDS file"
        )
      )
      md <- sc@meta.data
      req(input$chain_type_olga)
      md2 <- md[,names(md) %in% c("junction_aa_BD","cdr3_BD","junction_aa_AG","cdr3_AG","chain_AG","chain_BD")]
      print(names(md2))
      print("cdr3_BD"  %in% names(md2))
      if("cdr3_BD"  %in% names(md2)) {
        if(input$chain_type_olga == "TRB") {
          md2 <- subset(md2,md2$chain_BD == "TRB")
          md2 <- md2[,names(md2) %in% c("chain_BD","cdr3_BD")]
        } else {
          md2 <- subset(md2,md2$chain_AG == "TRA")
          md2 <- md2[,names(md2) %in% c("chain_AG","cdr3_AG")]
        }

      } else {
        if(input$chain_type_olga == "TRB") {
          md2 <- subset(md2,md2$chain_BD == "TRB")
          md2 <- md2[,names(md2) %in% c("chain_BD","junction_aa_BD")]
        } else {
          md2 <- subset(md2,md2$chain_AG == "TRA")
          md2 <- md2[,names(md2) %in% c("chain_AG","junction_aa_AG")]
        }
      }



      md2 <- md2[complete.cases(md2),]
      as.data.frame(md2)
      # df <-  md2[,names(md2) %in% c(input$cdr3_aa_olga)]

    })


    pgen_dt <- reactive({
      sc <- OLGA_data()
      validate(
        need(
          nrow(sc) > 0,
          "upload file"
        )
      )

      md <- sc@meta.data
      req(input$chain_type_olga)

      md2 <- md[,names(md) %in% c("junction_aa_BD","cdr3_BD","junction_aa_AG","cdr3_AG","chain_AG","chain_BD")]
      if("cdr3_BD"  %in% names(md2)) {

        if(input$chain_type_olga == "TRB") {
          md2 <- subset(md2,md2$chain_BD == "TRB")
          md2 <- md2[,names(md2) %in% c("chain_BD","cdr3_BD")]
          df <- unique(md2$cdr3_BD)
          pgen_dat <- as.data.frame(matrix(ncol = 2, nrow = length(df)))

        } else {
          md2 <- subset(md2,md2$chain_AG == "TRA")
          md2 <- md2[,names(md2) %in% c("chain_AG","cdr3_AG")]
          df <- unique(md2$cdr3_AG)
          pgen_dat <- as.data.frame(matrix(ncol = 2, nrow = length(df)))
        }

      }

      else {
        if(input$chain_type_olga == "TRB") {
          md2 <- subset(md2,md2$chain_BD == "TRB")
          md2 <- md2[,names(md2) %in% c("chain_BD","junction_aa_BD")]
          df <- unique(md2$junction_aa_BD)
          pgen_dat <- as.data.frame(matrix(ncol = 2, nrow = length(df)))

        } else {
          md2 <- subset(md2,md2$chain_AG == "TRA")
          md2 <- md2[,names(md2) %in% c("chain_AG","junction_aa_AG")]
          df <- unique(md2$junction_aa_AG)
          pgen_dat <- as.data.frame(matrix(ncol = 2, nrow = length(df)))

        }
      }

      colnames(pgen_dat) = c("cdr3_aa","pgen")

      for (i in 1:length(df)) {

        if(input$chain_type_olga == "TRB") {
          dfoutput <- foreach(dfrow=iter(df[i], by='row'), .combine=rbind) %do%
            olgafunction_BD(dfrow)
        } else {
          dfoutput <- foreach(dfrow=iter(df[i], by='row'), .combine=rbind) %do%
            olgafunction_AG(dfrow)
        }
        text <- dfoutput[2]
        split_substring <- str_split(text, " ")[[1]]

        # Remove the colon character from the first element
        split_substring[7] <- str_replace_all(split_substring[7], ":", "")
        pgen_dat[i,1] <- split_substring[7]
        pgen_dat[i,2] <- as.numeric(split_substring[8])
        message(paste("Compeleted",i,"of",length(df)))
        # print(pgen_dat[i,])
      }

      pgen_dat

    })

    output$Pgen_BD <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 20, scrollX = TRUE), {
      sc <- OLGA_data()
      validate(
        need(
          nrow(sc) > 0,
          "Upload .RDS file"
        )
      )
      pgen_dt()

    })

    output$downloaddf_pgen_dt <- downloadHandler(
      filename = function() {
        paste(input$chain_type_olga,"_pgen", ".csv", sep = "")
      },
      content = function(file) {
        df <- as.data.frame(pgen_dt())
        write.csv(df, file, row.names = F)
      }
    )


    # extract meta-data from sc object -----

    getData_extract_md <- reactive({
      inFile_extract.metadata <- input$file1_extract.metadata
      if (is.null(inFile_extract.metadata)) {
        return(NULL)
      }
      LoadSeuratRds(inFile_extract.metadata$datapath)
    })

    extracting_md_from_sc <- reactive({
      sc <- getData_extract_md()
      req(sc)

      md <- sc@meta.data
      md

    })

    output$extracting_md_from_sc_tb <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 20, scrollX = TRUE), {
      calls <- extracting_md_from_sc()
      req(calls)
      calls
    })
    # Add a download button
    output$download_extracting_md_from_sc_tb <- downloadHandler(
      filename = function() {
        paste("Meta_data_",today(), ".csv", sep = "")
      },
      content = function(file) {
        filtered_data_df <- extracting_md_from_sc()
        write.csv(filtered_data_df, file, row.names = F)
      }
    )


    ### filtering clustering table in post analysis ------
    getData_FilteringCluster <- reactive({
      inFile_sc_FilteringCluster <- input$file1_FilteringCluster
      if (is.null(inFile_sc_FilteringCluster)) {
        return(NULL)
      }

      read.csv(inFile_sc_FilteringCluster$datapath,header = T)
    })

    output$filter_checkboxes <- renderUI({

      data_filtered <- getData_FilteringCluster()
      req(data_filtered)

      selectizeInput(
        inputId = "selected_filters",
        label = "Select Filters to Include:",
        choices = names(data_filtered),
        selected = NULL,
        multiple = TRUE,
        options = list(placeholder = "Select filters")
      )
    })

    output$filter_inputs <- renderUI({
      # Create a list to store selectizeInput and checkbox widgets
      filter_inputs <- list()
      data_filtered <- getData_FilteringCluster()
      req(data_filtered)

      # Loop through the selected filters
      selected_filters <- input$selected_filters


      for (i in seq_along(selected_filters)) {
        selected_filter <- selected_filters[i]
        print(paste0("filter_ui_", selected_filter))

        filter_inputs[[selected_filter]] <- tagList(
          column(
            width = 2,
            checkboxInput(
              inputId = paste0("all_none_", selected_filter),
              label = paste0("All/None ", selected_filter),
              value = FALSE
            )
          ),


          column(
            width = 10,  # Set width to 10 units (adjust as needed)
            uiOutput(paste0("filter_ui_", selected_filter))  # Ensure unique IDs
          )
        )
      }

      # Wrap the list of selectizeInput and checkbox widgets in a tagList
      tagList(filter_inputs)
    })

    # Update the server logic to dynamically generate selectizeInput based on checkbox value
    observe({
      data_filtered <- getData_FilteringCluster()
      req(data_filtered)

      selected_filters <- input$selected_filters

      # Use lapply to generate all output elements at once
      lapply(seq_along(selected_filters), function(i) {
        selected_filter <- selected_filters[i]
        print(selected_filter)
        output[[paste0("filter_ui_", selected_filter)]] <- renderUI({
          checkbox_value <- input[[paste0("all_none_", selected_filter)]]
          selectizeInput(
            inputId = paste0("filter_", selected_filter),  # Ensure unique IDs
            label = selected_filter,
            choices = unique(data_filtered[[selected_filter]]),
            selected = if (checkbox_value) unique(data_filtered[[selected_filter]]) else NULL,
            multiple = TRUE,
            width = "900px",  # Set the width here
            options = list(
              placeholder = "Select values"
            )
          )
        })
      })
    })

    # Update the filtered_data reactive to handle the checkbox value
    filtered_data <- reactive({
      data_filtered <- getData_FilteringCluster()
      req(data_filtered)

      selected_filters <- input$selected_filters
      # print(selected_filters)
      # If no filters selected, return the original data
      if (length(selected_filters) == 0) {
        return(data_filtered)
      }

      # Apply filters based on user selection
      filtered_data <- data_filtered
      for (i in seq_along(selected_filters)) {
        filter <- selected_filters[i]
        checkbox_value <- input[[paste0("all_none_", filter)]]
        selected_values <- input[[paste0("filter_", filter)]]
        # print(selected_values)
        if (!checkbox_value && length(selected_values) > 0) {
          filtered_data <- filtered_data[filtered_data[[filter]] %in% selected_values, ]
        } else if (checkbox_value && length(selected_values) > 0) {

          filtered_data <- filtered_data[filtered_data[[filter]] %in% selected_values, ]
        }
      }

      filtered_data
    })



    # Render the filtered table
    output$filtered_table <- renderDT({
      filtered_data()
    }, options = list(
      scrollX = TRUE,
      columnDefs = list(list(targets = "_all", className = "dt-no-wrap"))
    ))

    # Add a download button
    output$download_button_filtered_data_df <- downloadHandler(
      filename = function() {
        selected_filters <- paste(input$selected_filters, collapse = "_")
        paste("filtered_data_", selected_filters, "_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        filtered_data_df <- filtered_data()
        write.csv(filtered_data_df, file, row.names = FALSE)
      }
    )



    # output$download_filters_button <- downloadHandler(
    #   filename = function() {
    #     paste("selected_filters", Sys.Date(), ".txt", sep = "_")
    #   },
    #   content = function(file) {
    #     selected_filters <- input$selected_filters
    #     writeLines(selected_filters, file)
    #   }
    # )
    # summarizing a csv file ------
    getData_summarising_data <- reactive({
      inFile_summarising_data <- input$file1_summarising_data
      if (is.null(inFile_summarising_data)) {
        return(NULL)
      }
      read.csv(inFile_summarising_data$datapath)
    })

    observe({
      dat <- getData_summarising_data()
      updateSelectInput(
        session,
        "columns_for_summary",
        choices = names(dat),
        selected = "Sample_Name"
      )
    })

    chain_table_summary <- reactive({
      df <- getData_summarising_data()
      req(df)

      df <- as.data.frame(df)
      df$cloneCount <- 1
      df2 <- df[,c("cloneCount",input$columns_for_summary)]
      df2 <- as.data.frame(df2)
      df3 <- as.data.frame(ddply(df2,input$columns_for_summary,numcolwise(sum)))
      df3
    })

    output$loaded_for_summary <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 2, scrollX = TRUE), {
      calls <- getData_summarising_data()
      req(calls)
      calls
    })

    output$chain_table_summary_tb <- DT::renderDT(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2, 5, 10, 20, 50, 100), pageLength = 20, scrollX = TRUE), {

      calls <- chain_table_summary()
      req(calls)
      calls
    })

    # Add a download button
    output$download_button_summarising_data <- downloadHandler(
      filename = function() {
        list_of_ <- paste(input$columns_for_summary, collapse = " ")
        # print(list_of_)
        paste("summary_",c(list_of_),"_",today(), ".csv", sep = "")
      },
      content = function(file) {
        df <- chain_table_summary()
        write.csv(df, file, row.names = F)
      }
    )
    # extract meta data from filtered object -----

    # extracting parameters and replacing defaults ------

    # Function to get all input values within a specific tab panel
    # Function to get all input values
    getAllInputValues <- function(session) {
      inputs <- list()
      for (name in names(session$input)) {
        inputs[[name]] <- isolate(input[[name]])
      }
      inputs <- Filter(function(x) !is.null(x), inputs)
      return(inputs)
    }

    # Data to be downloaded
    data_step4_params <- reactive({
      df <- data.frame(Parameter = names(getAllInputValues(session)),
                       Value = unlist(getAllInputValues(session)))
      return(df)
    })

    # output$parameterTable <- renderDT({
    #   df <- data.frame(Parameter = names(getAllInputValues(session)),
    #                    Value = unlist(getAllInputValues(session)))
    #   return(df)
    # })

    # Download handler
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("parameters_state", Sys.Date(), ".csv", sep = "_")
      },
      content = function(file) {
        write.csv(data_step4_params(), file, row.names = FALSE)
      }
    )

    # Upload handler
    observeEvent(input$uploadData, {
      req(input$uploadData)
      inFile <- input$uploadData
      if (!is.null(inFile)) {
        uploaded_data <- read.csv(inFile$datapath, stringsAsFactors = FALSE)
        for (i in 1:nrow(uploaded_data)) {
          param_name <- uploaded_data$Parameter[i]
          param_value <- uploaded_data$Value[i]
          # Check input type and update accordingly
          if (is.numeric(input[[param_name]])) {
            param_value <- as.numeric(param_value)
          } else if (is.character(input[[param_name]])) {
            # Nothing to do for character inputs
          } else if (is.logical(input[[param_name]])) {
            param_value <- as.logical(param_value)
          }
          updateInput(session, param_name, value = param_value)
        }
      }
    })

  }
  ########
  # run the app in browser -----
  shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
}
