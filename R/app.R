#' Run STEGO application.
#' @name STEGO
#' @importFrom stats chisq.test
#' @import BiocParallel
#' @import ClusTCR2
#' @export

runSTEGO <- function(){
  source(system.file("Global","required_functions.R",package = "STEGO.R"))
  suppressWarnings(source(system.file("scGATE","custom_df_scGATE.R",package = "STEGO.R")))
  # UI page -----
  ui <- fluidPage(
    theme=bs_theme(version = 5, bootswatch = "default"),
    navbarPage(title = "STEGO_R",
               theme=bs_theme(version = 5, bootswatch = "default"),
               navbarMenu("Quality control",
                          ## 10x Genomics ----
                          tabPanel("10x genomics",
                                   sidebarLayout(
                                     sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                                  # selectInput("dataset_10x", "Choose a dataset:", choices = c("test_data_10x", "own_data_10x")),
                                                  fileInput('file_calls_10x', 'Barcode file (.tsv.gz or .tsv)',
                                                            accept=c('.tsv','.tsv.gz')),
                                                  fileInput('file_features_10x', 'Features file (.tsv.gz or .tsv)',
                                                            accept=c('.tsv','.tsv.gz')),
                                                  fileInput('file_matrix_10x', 'Matrix file (.mtx.gz or .mtx)',
                                                            accept=c('.mtx.gz','.mtx')),
                                                  fileInput('file_TCR_10x', 'filtered contig annotations (.csv)',
                                                            accept=c('.csv')),
                                                  textInput("name.10x","Name added to files",value = ""),
                                                  textInput("sample_name_10x","Add sample name","Treatment_group"),
                                                  h6("TCR_explore name"),
                                                  fluidRow(
                                                    column(6,textInput("group10x","Add sample name","Group")),
                                                    column(6,textInput("Indiv10x","Add sample name","Indiv"))
                                                  ),
                                                  selectInput("BCR_TCR_10x","Type of data",choices = c("TCR only","BCR only","both")),

                                     ),
                                     ### 10x main panel -----
                                     mainPanel(
                                       tabsetPanel(id = "panel_10x",
                                                   tabPanel("Uploaded data",value = 1,
                                                            div(DT::dataTableOutput("test.files.10x1")),
                                                            div(DT::dataTableOutput("test.files.10x2")),
                                                            # div(DT::dataTableOutput("test.files.10x3")),
                                                            div(DT::dataTableOutput("test.files.10x4")),
                                                   ),
                                                   tabPanel("TCRex",
                                                            add_busy_spinner(spin = "fading-circle"),
                                                            div(DT::dataTableOutput("tb_TCRex_10x_df")),
                                                            downloadButton('downloaddf_TCRex_10x','Download table')

                                                   ),
                                                   tabPanel("For Seurat QC",value = 2,
                                                            tags$head(tags$style("#tb_10x_matrix2  {white-space: nowrap;  }")),
                                                            add_busy_spinner(spin = "fading-circle"),
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
                          # 10x genomics end -----


                          ## BD Rhapsody  ------
                          tabPanel("BD rhapsody data",
                                   sidebarLayout(
                                     sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                                  # UPLOAD the three files...
                                                  # selectInput("dataset_BD", "Choose a dataset:", choices = c("test_data_BD", "own_data_BD")),
                                                  fileInput('file_calls_BD', 'Sample Tag Calls (.csv)',
                                                            accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                                  numericInput("no_lines_skip_Tags","Information present? skip first 7 lines",value = 0,min=0,max=10,step=7),
                                                  fileInput('file_TCR_BD', 'TCR file (.csv)',
                                                            accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                                  numericInput("no_lines_skip_TCR","Information present? skip first 7 lines",value = 0,min=0,max=10,step=7),
                                                  fileInput('file_counts_BD', 'Counts (.csv)',
                                                            accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                                  numericInput("no_lines_skip_counts","Information present? skip first 7 lines",value = 0,min=0,max=10,step=7),
                                                  ### filter out non-function TCR and un-paired TCR
                                                  checkboxInput("filter_zero_expression", "remove cells with no mRNA expression?", value = FALSE, width = NULL),
                                                  checkboxInput("filtering_TCR", "Keep functional paired chains?", value = FALSE, width = NULL),
                                                  checkboxInput("filter_non_function", "Remove non-functional TCR/BCR?", value = FALSE, width = NULL),
                                                  checkboxInput("BCR_present", "BCR present?", value = FALSE, width = NULL),

                                                  textInput("name.BD","Name added to files",value = ""),
                                     ),
                                     mainPanel(
                                       shiny::tabsetPanel(
                                         tabPanel("Imported data",
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("test.files")),
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("test.files2")),
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("test.files3")),

                                         ),
                                         tabPanel("Filtering",
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("Filtering_BD")),

                                         ),
                                         tabPanel("clusTCR2",
                                                  tags$head(tags$style("#tb_clusTCR  {white-space: nowrap;  }")),
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("tb_clusTCR")),
                                                  downloadButton('downloaddf_clusTCR','Download table')
                                         ),

                                         tabPanel("TCRex",
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("tb_TCRex_BDrap_df")),
                                                  downloadButton('downloaddf_TCRex_BDrap','Download table')

                                         ),

                                         tabPanel("For Seurat",
                                                  tags$head(tags$style("#tb_count_matrix  {white-space: nowrap;  }")),
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("tb_count_matrix")),
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("tb_clusTCR_sum")),
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("tb_metadata_sc")),
                                                  fluidRow(
                                                    column(3,downloadButton('downloadtb_count_matrix','Download count table')),
                                                    column(3),
                                                    column(3,downloadButton('downloadtb_metadata_sc','Download meta.data table')),
                                                  ),
                                         ),
                                         tabPanel("TCR_Explore",
                                                  tags$head(tags$style("#tb_TCR_Explore  {white-space: nowrap;  }")),
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("tb_TCR_Explore")),
                                                  downloadButton('downloadtb_TCR_Explore','Download table')
                                         ),
                                         tabPanel("Create Sample Tags file",
                                                  tags$head(tags$style("#tb_sample_tags_created  {white-space: nowrap;  }")),
                                                  textInput("sample_tags_name","Name of sample",value = "BD EA splenocyte"),
                                                  div(DT::dataTableOutput("tb_sample_tags_created")),
                                                  downloadButton('downloadtb_sample_tags','Download Tags')
                                         ),
                                       )
                                     )

                                   )
                          ),
                          # BD rhapsody end ----
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
                                         tabPanel("Upload files",
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("test.files_array_Matrix")),
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("test.files_array_contig")),

                                         ),
                                         tabPanel("Filtering",
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("test.files_array_contig_Filtered")),
                                         ),
                                         tabPanel("clusTCR",
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("test.files_ClusTCR2_array")),
                                                  downloadButton("download_ClusTCR2_labs_array"),
                                         ),
                                         tabPanel("TCRex"),
                                         tabPanel("Seurat",
                                                  add_busy_spinner(spin = "fading-circle"),
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
               ### TCR clustering with ClusTCR2 -----
               tabPanel("ClusTCR2",
                        sidebarLayout(
                          sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                       # selectInput("dataset2", "Choose a dataset:", choices = c("test_data_clusTCR2","own_data_clusTCR2")),
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
                                                          column(6, numericInput("cores_ClusTCR2","Number of cores to parallel",value=4))
                                                        )
                                       ),
                          ),
                          mainPanel(
                            tabsetPanel(id = "clusTCR2_tabs",
                                        tabPanel("Merge Multiple Files",value = "merge",
                                                 div(DT::dataTableOutput("DEx_multiple_ClusTCR2")),
                                                 downloadButton('downloaddf_multiple_ClusTCR2','Download table')
                                        ),
                                        tabPanel("Uploaded file",value = "processing1",
                                                 div(DT::dataTableOutput("clust_dt2")),
                                        ),
                                        tabPanel("Outputs",value = "processing2",
                                                 tabsetPanel(
                                                   tabPanel('Cluster Labels',
                                                            add_busy_spinner(spin = "fading-circle"),
                                                            div(DT::dataTableOutput("ClusTCR2_lab")),
                                                            downloadButton('download_ClusTCR_labels','Download Cluster table'),
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
                                                   tabPanel("Time",
                                                            add_busy_spinner(spin = "fading-circle"),
                                                            verbatimTextOutput('ClusTCR2_Time'),
                                                            # div(DT::dataTableOutput("")),
                                                   )
                                                 ))

                            )
                          )
                        )
               ),
               ### end TCR clustering ------
               tabPanel("Seurat QC",
                        sidebarLayout(
                          sidebarPanel(id = "tPanel2",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                       # selectInput("dataset_sc", "Choose a dataset:", choices = c("test_data_sc", "own_data_sc")),
                                       # upload the file
                                       fileInput('file_SC', 'Load csv (BDrhap) or csv.gz (10x)',
                                                 accept=c('text/csv','.csv','.csv.gz','gz','csv.gz','h5','.h5'),
                                       ),
                                       textInput("project_name","Name of sample",value = ""),
                                       # selectInput("species","Species",choices = c("human","mouse","other")),
                                       selectInput("df_seruatobj_type","Data type", choices = c("10x Genomics","BD Rhapsody","Array")),
                                       uiOutput("feature_input"),
                                       actionButton("run","Update Violin plot"),
                                       fluidRow(
                                         column(6,numericInput("dimension_heatmap.min","View heatmap dimensions.min", value = 1)),
                                         column(6,numericInput("dimension_heatmap.max","View heatmap dimensions.max", value = 10)),
                                         column(6,numericInput("numberofcells","Number of cells to use for heatmap", value = 500))
                                       ),

                                       actionButton("run_reduction","Run clustering"),
                                       fluidRow(
                                         column(6,numericInput("dimension_sc","Max dimensions for clustering", value = 15)),
                                         column(6,numericInput("resolution","Resolution of clusters", value = 1)),
                                       ),

                                       fileInput('file_SC_meta', 'Upload file meta.data file (.csv.gz or .csv)',
                                                 accept=c('.csv','.csv.gz')),
                                       actionButton("run_metadata","Impute metadata after clustering"),
                                       # selectInput("save_type","Type of output",choices = c(".h5Seurat",".rds")),
                                       add_busy_spinner(spin = "fading-circle"),
                                       downloadButton('downloaddf_SeruatObj','Download Seurat')
                                       # column(4,numericInput("percent.mt","Mictochondrial DNA cut-off", value = 20)),
                          ),
                          mainPanel(
                            tabsetPanel(
                              # QC panel -----
                              tabPanel("QC plots",
                                       tabsetPanel(
                                         tabPanel("Header check",
                                                  div(DT::dataTableOutput("DEx_header_name_check.dt")),
                                         ),
                                         tabPanel("Violin and correlation",
                                                  tabsetPanel(
                                                    tabPanel("Before",
                                                             add_busy_spinner(spin = "fading-circle"),
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
                                                             add_busy_spinner(spin = "fading-circle"),
                                                             plotOutput("after_plot_sc", height = "600px"),
                                                             fluidRow(
                                                               column(3,numericInput("width_Before.after_plot_sc", "Width of PDF", value=10)),
                                                               column(3,numericInput("height_Before.after_plot_sc", "Height of PDF", value=8)),
                                                               column(3),
                                                               column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_Before.after_plot_sc','Download PDF'))),
                                                    )
                                                  ),
                                         ),
                                         #### Variable features -----
                                         tabPanel("Top variable features",
                                                  add_busy_spinner(spin = "fading-circle"),
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
                                         #### Elbow and heatmap  -----
                                         tabPanel("Elbow Plot",
                                                  plotOutput("create_elbowPlot_sc", height = "600px"),
                                                  fluidRow(
                                                    column(1,numericInput("width_create_elbowPlot_sc", "Width of PDF", value=10)),
                                                    column(1,numericInput("height_create_elbowPlot_sc", "Height of PDF", value=8)),
                                                    column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_create_elbowPlot_scc','Download PDF')),
                                                    column(2,numericInput("width_png_create_elbowPlot_sc","Width of PNG", value = 1200)),
                                                    column(2,numericInput("height_png_create_elbowPlot_sc","Height of PNG", value = 1000)),
                                                    column(2,numericInput("resolution_PNG_create_elbowPlot_sc","Resolution of PNG", value = 144)),
                                                    column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_create_elbowPlot_sc','Download PNG')),
                                                  ),

                                         ),
                                         tabPanel("DimHeatmap",
                                                  add_busy_spinner(spin = "fading-circle"),

                                                  plotOutput("create_PCA_heatmap_sc", height = "1000px")
                                         ),
                                         # tabPanel("Resolution plot"),
                                         #### UMAP  -----
                                         tabPanel("UMAP",
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  plotOutput("create_UMAP_sc", height = "600px")
                                         ), # Export a table with meta.data and expression.
                                         tabPanel("Add metadata",
                                                  div(DT::dataTableOutput("DEx_view.meta.dt")),
                                                  div(DT::dataTableOutput("DEx_table_meta.data")),
                                                  # selectInput("","comaprison 1")
                                         ),


                                       ),


                              ),
                              # end of QC -----
                            ),
                          )
                        )
               ),
               # Merge multiple Seurat objects -----
               tabPanel("Merge SC",
                        sidebarLayout(
                          # Sidebar with a slider input
                          sidebarPanel(id = "tPanel5",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                       fileInput("file1_h5Seurat.file",
                                                 "Choose .h5Seurat files from directory",
                                                 multiple = TRUE,
                                                 accept=c("h5Seurat",".h5Seurat")),
                                       textInput("project_name2","Name of Project",value = ""),
                                       downloadButton('downloaddf_SeruatObj_merged','Download Merged Seurat')
                          ),

                          # Show a plot of the generated distribution
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Upload files",
                                       add_busy_spinner(spin = "fading-circle"),
                                       verbatimTextOutput("testing_mult"),
                                       # verbatimTextOutput("testing_mult2")
                              ),
                              tabPanel("Variable data",
                                       add_busy_spinner(spin = "fading-circle"),
                                       actionButton("run_var","Run VariableFeatures"),
                                       add_busy_spinner(spin = "fading-circle"),
                                       verbatimTextOutput("var_harmony_verbrose")
                              ),
                              tabPanel("Scale data",
                                       add_busy_spinner(spin = "fading-circle"),
                                       actionButton("run_scale","Run Scale"),
                                       add_busy_spinner(spin = "fading-circle"),
                                       verbatimTextOutput("scale_harmony_verbrose")
                              ),
                              tabPanel("PCA",
                                       add_busy_spinner(spin = "fading-circle"),
                                       actionButton("run_PCA","Run PCA"),
                                       add_busy_spinner(spin = "fading-circle"),
                                       verbatimTextOutput("PCA_harmony_verbrose")
                              ),
                              tabPanel("harmony",
                                       add_busy_spinner(spin = "fading-circle"),
                                       actionButton("run_harmony","Run Harmony"),
                                       add_busy_spinner(spin = "fading-circle"),
                                       verbatimTextOutput("harmony_verbrose"),
                              ),
                              tabPanel("Dimentional reduction",
                                       fluidRow(
                                         column(3,numericInput("dimension_Merged","Max number of dimensions", value = 30)),
                                         column(6,numericInput("res_merged","Resolution of clusters", value = 0.5)),
                                       ),
                                       actionButton("run_reduction_harmony","Run Dimentional reduction"),
                                       add_busy_spinner(spin = "fading-circle"),
                                       verbatimTextOutput("testing_mult3")

                              ),
                              tabPanel("UMAP",
                                       add_busy_spinner(spin = "fading-circle"),
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
               # Add annotations -----
               tabPanel("Annotations",
                        sidebarLayout(
                          sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                       fileInput("file1_h5Seurat.file2",
                                                 "Choose merged or single .h5Seurat files from directory",
                                                 multiple = TRUE,
                                                 accept=c('.h5Seurat','h5Seurat')),
                                       textInput("project_name3","Name of Project",value = ""),
                                       downloadButton('downloaddf_SeruatObj_annotated','Download Annotated Seurat')
                          ),

                          mainPanel(
                            tabsetPanel(
                              tabPanel("Upload",
                                       add_busy_spinner(spin = "fading-circle"),
                                       verbatimTextOutput("testing_mult_anno")
                              ),
                              tabPanel("scGATE",
                                       tabsetPanel(
                                         tabPanel("Annotations to add",
                                                  selectInput("GenericID_scGATE","Generic To include",choices = c("Bcell","CD4T","CD8T","CD8TIL" ,"Erythrocyte" ,"Megakaryocyte" , "MoMacDC","Myeloid","NK","PanBcell","panDC","PlasmaCell","Tcell","Tcell.alphabeta"), selected = c("Bcell","MoMacDC","NK","CD8T","CD4T","PlasmaCell"), multiple = T),
                                                  fluidRow(column(3,checkboxInput("generic_scGATE","Generic (Human)", value = F)),
                                                           column(3,checkboxInput("CD4_scGATE","CD4 T cell (Human)", value = F)),
                                                           column(3,checkboxInput("CD8_scGATE","CD8 T cell (Human)", value = F)),
                                                           column(3,checkboxInput("cellTypist_higher_scGATE","CellTypist (higher) T cell (Human)", value = F)),
                                                           column(3,checkboxInput("cellTypist_lowerBcells_scGATE","CellTypist (lower) B cell (Human)", value = F)),
                                                           column(3,checkboxInput("cellTypist_lowerCD4T_scGATE","CellTypist (lower) CD4 T cell (Human)", value = F)),
                                                           column(3,checkboxInput("cellTypist_lowerCD8T_scGATE","CellTypist (lower) CD8 T cell (Human)", value = F)),
                                                           column(3,checkboxInput("cellTypist_lowerOtherT_scGATE","CellTypist (lower) Other T cell (Human)", value = F)),
                                                           column(3,checkboxInput("Ex.Sen_scGATE","Exhausted/Senescence (Human)", value = F)),
                                                           column(3,checkboxInput("COVID_scGATE","COVID markers (Human)", value = F)),
                                                           column(3,checkboxInput("Activation_scGATE","Activation markers (Human)", value = F)),
                                                           column(3,checkboxInput("IFNgTNFa_scGATE","IFNg and TNFa (Human)", value = F)),
                                                           column(3,checkboxInput("GNLY.PFR1.GZMB_scGATE","GNLY.PFR1.GZMB markers (Human)", value = F)),
                                                           column(3,checkboxInput("Interlukin_scGATE","Interlukin markers (Human)", value = F)),
                                                  ),
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  fluidRow(
                                                    column(6,verbatimTextOutput("scGATE_verbatum_Generic")),
                                                    column(6,verbatimTextOutput("scGATE_verbatum_CD4")),
                                                    column(6,verbatimTextOutput("scGATE_verbatum_CD8")),
                                                    column(6,verbatimTextOutput("scGATE_verbatum_cellTypist_higher")),
                                                    column(6,verbatimTextOutput("scGATE_verbatum_cellTypist_lower_Bcells")),
                                                    column(6,verbatimTextOutput("scGATE_verbatum_cellTypist_lower_CD4T")),
                                                    column(6,verbatimTextOutput("scGATE_verbatum_cellTypist_lower_CD8T")),
                                                    column(6,verbatimTextOutput("scGATE_verbatum_cellTypist_lower_OtherT")),
                                                    column(6,verbatimTextOutput("scGATE_verbatum_Ex.Sen")),
                                                    column(6,verbatimTextOutput("scGATE_verbatum_COVID")),
                                                    column(6,verbatimTextOutput("scGATE_verbatum_Activation")),
                                                    column(6,verbatimTextOutput("scGATE_verbatum_IFNgTNFa")),
                                                    column(6,verbatimTextOutput("scGATE_verbatum_GNLY.PFR1.GZMB")),
                                                    column(6,verbatimTextOutput("scGATE_verbatum_Interlukin")),
                                                  )
                                         ),
                                         tabPanel("Table",
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("DEx_table_TcellClass_scGATE")),
                                         )
                                       )
                              ),
                              tabPanel("K-means clustering",
                                       tabsetPanel(
                                         tabPanel("Checking annotation",
                                                  selectInput("V_gene_Class_2","V gene with/without CDR3",choices = ""),
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("DEx_table_TcellClass_3")),
                                         ),
                                         tabPanel("Classification to include",
                                                  p("Curated list"),
                                                  fluidRow(column(3,checkboxInput("add.classification_T_cell_2","General Adaptive cell markers", value = T)),
                                                           column(3,checkboxInput("add.classification_T_cell_Function_2","Function", value = T)),
                                                           column(3,checkboxInput("add.classification_T_cell_Function_CD4_2","Function (CD4)", value = T)),
                                                           column(3,checkboxInput("add.classification_T_cell_Function_CD8_2","Function (CD8)", value = T)),
                                                           column(3,checkboxInput("add.classification_T_cell_Function_CD4_CD8pos_2","Function (CD4+CD8+)", value = T)),
                                                           column(3,checkboxInput("add.classification_T_cell_Function_CD4_CD8neg_2","Function (CD4-CD8-)", value = T)),
                                                           column(3,checkboxInput("add.classification_B_cell_Function_2","Function B cell", value = T)),

                                                           column(3,checkboxInput("add.classification_T_cell_Memory_2","Memory", value = T)),
                                                           column(3,checkboxInput("add.classification_T_cell_Activation_2","Activation status", value = T)),),
                                                  h5("Cell Typist based list"),
                                                  fluidRow(
                                                    column(3,checkboxInput("add.classification_CellTypist_list_overview_2","CellTypist Overview", value = T)),
                                                    column(3,checkboxInput("add.classification_CellTypist_list_2","CellTypist T cell classification", value = T)),
                                                    column(3,checkboxInput("add.classification_CellTypist_cycling_2","CellTypist Cell cycle", value = T)),
                                                    column(3,checkboxInput("add.classification_CellTypist_list_overview_3","CellTypist & curated list", value = T)),
                                                    # column(3,)
                                                  ),

                                         ),

                                         # tabPanel("Data2Talk",
                                         #          p("Upload",tags$a(href="https://talk2data.bioturing.com/predict", "Data2Talk prediction")),
                                         #
                                         #
                                         #
                                         #          ),

                                         tabPanel("Table",
                                                  add_busy_spinner(spin = "fading-circle"),
                                                  div(DT::dataTableOutput("DEx_table_TcellClass_2")),

                                         ),


                                       ),
                              ),
                              tabPanel("Meta data table",
                                       fluidRow(
                                         column(3,checkboxInput("add.kmeans","Add K-means classification", value = T)),
                                         column(3,checkboxInput("add.scGATE","Add scGATE classifications", value = T))
                                       ),
                                       add_busy_spinner(spin = "fading-circle"),
                                       div(DT::dataTableOutput("DEx_table_TcellClass_scGATE.kmeans")),
                              )
                            )


                          ),

                        )

               ),


               ## Analysis (UI side panel)---------
               tabPanel("Analysis",
                        sidebarLayout(
                          sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,

                                       selectInput("STEGO_R_pro","QC processed",choices = c("STEGO_R (.h5Seurat)")), #,"Seurat (.rds)"
                                       textInput("name.file_clust","Name added to files",value = ""),
                                       conditionalPanel(condition="input.check_up_files== 'up'",


                                                        fileInput('file_SC_pro', 'Upload seurat file',
                                                                  accept=c("h5Seurat",'.h5Seurat','rds')),
                                                        fileInput('file_cluster_file', 'Upload clustering file from clusTCR2 (.csv)',
                                                                  accept=c('csv')),

                                                        numericInput("skip_TCRex_up","Skip # of lines for TCRex file",value = 7),
                                                        fileInput('upload_TCRex_file', 'Upload TCRex (.tsv)',
                                                                  accept=c('tsv','.tsv')),

                                       ),
                                       conditionalPanel(condition="input.STEGO_R_pro == 'STEGO_R (.h5Seurat)'",
                                                        selectInput("datasource", "Data source",choices=c("10x Genomics","BD Rhapsody")),
                                       ),
                                       fluidRow(column(6,selectInput("Samp_col","Sample column name",choices = "")),
                                                column(6,selectInput("Split_by_group","Include group comparison",choices=c("no","yes"))),
                                                column(6,conditionalPanel(condition="input.STEGO_R_pro == 'STEGO_R (.h5Seurat)'",
                                                                          selectInput("V_gene_sc","V gene with/without CDR3",choices = "")))
                                       ),

                                       fluidRow(
                                         column(4,selectInput("colourtype","Type of colouring",choices = c("default","rainbow","random","heat.colors","terrain.colors","topo.colors","hcl.colors","one colour"))),
                                         column(4,colourInput("one.colour.default","One colour","grey")),
                                         column(4,colourInput("NA_col_analysis","NA colour","grey75"),)

                                       ),

                                       conditionalPanel(condition="input.QC_panel==2 || input.QC_panel==3",
                                                        fluidRow(column(12,selectInput("Graph_type_bar","Type of graph",choices = c("Number_expanded","Clonality","Top_clonotypes")))),
                                                        fluidRow(
                                                          # column(4,selectInput("Graph_type_bar","Type of graph",choices = c("Number_expanded","Clonality","Top_clonotypes")),),

                                                          # column(4,numericInput("number_colour_req","Number of colours needed",value = 5))
                                                          conditionalPanel(condition="input.Graph_type_bar =='Top_clonotypes'",
                                                                           column(6,numericInput("top_no_clonotypes","Top clonotypes per group",value = 1,step = 1, min = 0, max = 20)),
                                                          ),
                                                        ),
                                       ),

                                       conditionalPanel(condition="input.STEGO_R_pro=='Seurat (.rds)'",
                                                        h4("Add in chain information"),
                                                        fluidRow(
                                                          column(6,selectInput("RDS_V_gene_A","A/G gene column",choices = "")),
                                                          column(6,selectInput("RDS_V_gene_B","B/D gene column",choices = "")),
                                                          column(6,selectInput("RDS_cdr3_A","A/G cdr3 column",choices = "")),
                                                          column(6,selectInput("RDS_cdr3_B","B/D cdr3 column",choices = "")),
                                                        ),
                                       ),

                                       column(12,selectInput("ID_Column_factor","Order of graph",choices = "", multiple = T, width = "1200px")),


                                       selectInput("font_type",label = h4("Type of font"),choices = font,selected = "serif"),

                                       # conditionalPanel(condition="input.QC_panel==2 || input.QC_panel==3",
                                       fluidRow(
                                         column(6,numericInput("text_size","Size of #",value=16)),
                                         column(6,numericInput("title.text.sizer2","Axis text size",value=30)),

                                         column(6,numericInput("Bar_legend_size","Legend text size",value=16)),
                                         column(4, selectInput("legend_position","Legend location",choices = c("top","bottom","left","right","none"),selected = "right")),
                                       ),
                                       # ),
                                       p(""),
                                       # selectInput("dataset_cluster_file", "Choose a dataset:", choices = c("test_data_clusTCR", "own_data_clusTCR")),
                                       # numericInput("num", label = h3("Numeric input"), value = 1),
                                       # materialSwitch(inputId = "mode", label = icon("moon"),
                                       #                right=TRUE,status = "success")
                          ),

                          # add in clustering  -----
                          mainPanel(
                            tabsetPanel(id = "check_up_files",
                                        tabPanel("Check files uploaded",value = 'up',
                                                 add_busy_spinner(spin = "fading-circle"),
                                                 fluidRow(column(12,  div(DT::dataTableOutput("Tb_TCR_clonotypes.Umap"))),
                                                          column(12,   div(DT::dataTableOutput("Tb_ClusTCR_test"))),

                                                 ),
                                                 fluidRow(

                                                   column(12,   div(DT::dataTableOutput("Tb_tcrex_test")))
                                                 ),
                                        ),
                                        ### UMAP -> TCR -----
                                        tabPanel("Overview of TCR",
                                                 tabsetPanel(id = "QC_panel",
                                                             tabPanel("Check UMAP", value = 1,
                                                                      add_busy_spinner(spin = "fading-circle"),
                                                                      plotOutput("create_UMAP_sc2", height = "600px"),

                                                                      fluidRow(
                                                                        column(3,numericInput("width_UMAP2", "Width of PDF", value=10)),
                                                                        column(3,numericInput("height_UMAP2", "Height of PDF", value=8)),
                                                                        column(3),
                                                                        column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_UMAP2','Download PDF'))),
                                                                      fluidRow(
                                                                        column(3,numericInput("width_png_UMAP2","Width of PNG", value = 1200)),
                                                                        column(3,numericInput("height_png_UMAP2","Height of PNG", value = 1000)),
                                                                        column(3,numericInput("resolution_PNG_UMAP2","Resolution of PNG", value = 144)),
                                                                        column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_UMAP2','Download PNG'))

                                                                      ),
                                                             ),
                                                             tabPanel("Clonal expansion plots",value = 2,
                                                                      add_busy_spinner(spin = "fading-circle"),


                                                                      fluidRow(
                                                                        # uiOutput('myPanel_top_clonal_plot'),
                                                                        conditionalPanel(condition="input.Graph_type_bar=='Number_expanded' || input.Graph_type_bar=='Clonality'",
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
                                                                      # plotOutput("", height = "600px"),

                                                                      fluidRow(
                                                                        column(3,numericInput("width_clonality.bar.graph", "Width of PDF", value=10)),
                                                                        column(3,numericInput("height_clonality.bar.graph", "Height of PDF", value=8)),
                                                                        column(3),
                                                                        column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_clonaity.bar.graph','Download PDF'))),

                                                                      fluidRow(
                                                                        column(3,numericInput("width_png_clonality.bar.graph","Width of PNG", value = 2400)),
                                                                        column(3,numericInput("height_png_clonality.bar.graph","Height of PNG", value = 1400)),
                                                                        column(3,numericInput("resolution_PNG_clonality.bar.graph","Resolution of PNG", value = 144)),
                                                                        column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_clonaity.bar.graph','Download PNG'))
                                                                      )
                                                             ),
                                                             #### UMAP clonality -> TCR -----
                                                             tabPanel("UMAP with clonality (counts)", value = 3,
                                                                      add_busy_spinner(spin = "fading-circle"),
                                                                      # column(4,selectInput("Split_by_group","Include group comparison",choices=c("no","yes"))),
                                                                      conditionalPanel(condition="input.Graph_type_bar=='Number_expanded' || input.Graph_type_bar=='Clonality'",
                                                                                       fluidRow(
                                                                                         column(3,
                                                                                                wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                          uiOutput('cols_UMAP_clonal_plot'))),
                                                                                         column(9, plotOutput("clonality.TCR.UMAP",height="600px"))
                                                                                       ),

                                                                                       fluidRow(
                                                                                         column(3,numericInput("width_TCR.UMAP", "Width of PDF", value=10)),
                                                                                         column(3,numericInput("height_TCR.UMAP", "Height of PDF", value=5)),
                                                                                         column(3),
                                                                                         column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_TCR.UMAP','Download PDF'))),

                                                                                       fluidRow(
                                                                                         column(3,numericInput("width_png_TCR.UMAP","Width of PNG", value = 2400)),
                                                                                         column(3,numericInput("height_png_TCR.UMAP","Height of PNG", value = 1200)),
                                                                                         column(3,numericInput("resolution_PNG_TCR.UMAP","Resolution of PNG", value = 144)),
                                                                                         column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_TCR.UMAP','Download PNG'))
                                                                                       ),

                                                                      ),

                                                                      conditionalPanel(condition="input.Graph_type_bar=='Top_clonotypes'",
                                                                                       # column(12,   div(DT::dataTableOutput("Tb_For_colouring_check"))),

                                                                                       fluidRow(
                                                                                         # column(4,selectInput("Split_by_group","Include group comparison",choices=c("no","yes"))),
                                                                                         column(4,selectInput("display_all_samps","Display all sample",choices=c("yes","no"))),




                                                                                       ),

                                                                                       column(12,selectInput("ID_Column_metadata","Select to display",choices = "", multiple = T, width = "1200px")),

                                                                                       add_busy_spinner(spin = "fading-circle"),
                                                                                       fluidRow(
                                                                                         column(3,
                                                                                                wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                          uiOutput('cols_UMAP_Topclonotypes'))),
                                                                                         column(9, plotOutput("clonality.TCR.UMAP.top",height="600px")),

                                                                                         fluidRow(
                                                                                           column(3,numericInput("width_TCR.UMAP_top", "Width of PDF", value=10)),
                                                                                           column(3,numericInput("height_TCR.UMAP_top", "Height of PDF", value=8)),
                                                                                           column(3),
                                                                                           column(3,style = "margin-top: 25px;",downloadButton('downloadPlot_TCR.UMAP_top','Download PDF'))),

                                                                                         fluidRow(
                                                                                           column(3,numericInput("width_png_TCR.UMAP_top","Width of PNG", value = 1200)),
                                                                                           column(3,numericInput("height_png_TCR.UMAP_top","Height of PNG", value = 1400)),
                                                                                           column(3,numericInput("resolution_PNG_TCR.UMAP_top","Resolution of PNG", value = 144)),
                                                                                           column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_TCR.UMAP_top','Download PNG'))
                                                                                         ),

                                                                                       ),
                                                                      ),

                                                             ),
                                                             # tabPanel("Clonotype overlap per cluster (upset plot)",
                                                             #          # add in upset plot per cluster
                                                             #          )
                                                 )
                                        ),
                                        ### UMAP with TCR expression end -----
                                        tabPanel("Differential expression",

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
                                                             #### Cluster table -----
                                                             tabPanel("Cluster differences (Feature plot)",
                                                                      actionButton("run_string.data3","View Feature plot"),
                                                                      fluidRow(column(12, selectInput("string.data3","column names for summary","",multiple = T, width = "1200px") )),
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
                                                             #### differential expression within clusters ----
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
                                                                      add_busy_spinner(spin = "fading-circle"),
                                                                      div(DT::dataTableOutput("DEx_table_clusters")),


                                                                      downloadButton('downloaddf_DEx_table_clusters','Download Table (.csv)')



                                                             ),
                                                             #### Cluster table plot -----

                                                             # tabPanel("Cluster differences (Heatmap/expression dot plot) "),
                                                             # heatmap (Z-score)
                                                             # clusters, mean expression and percentage with expression (Top x genes per cluster?)


                                                             #### check meta data table ----


                                                 )

                                        ),
                                        ### end of differential expression -----
                                        ### TCR and GEX analysis section-----
                                        tabPanel("TCR and GEX",

                                                 fluidRow(
                                                   column(4,uiOutput("classification_to_add")),
                                                   column(4, selectInput("clust_group","Select cluster:",choices ="")),
                                                   column(4, selectInput("Gene_V_top","Select V with/without CDR3",choices =""))),
                                                 conditionalPanel(condition="input.Panel_TCRUMAP=='top_clone'",
                                                                  fluidRow(
                                                                    column(4,actionButton("run_top","Update selected clonotype")),
                                                                    column(4, selectInput("Selected_clonotype","Select clonotype:",choices ="")),
                                                                    column(4, selectInput("Selected_chain","Select chain:",choices ="")),
                                                                  )),
                                                 conditionalPanel(condition="input.Panel_TCRUMAP=='ClusTCR2'",
                                                                  fluidRow(
                                                                    column(3,selectInput("chain_TCR","Chains included",choices = c("TCR","BCR","both")),),
                                                                    # column(3,selectInput("V_gene_clust","V gene",choices = c("v_gene_AG","v_gene_BD","v_gene_IgL","v_gene_IgH")),),
                                                                    # column(3, selectInput("cdr3_clust","CDR3",choices = c("cdr3_AG","cdr3_BD","cdr3_IgL","cdr3_IgH")),)
                                                                  ),
                                                 ),

                                                 conditionalPanel(condition="input.Panel_TCRUMAP=='Epitope'",
                                                                  fluidRow(
                                                                    column(4,uiOutput("classification_to_add_epitope")),
                                                                    column(4,uiOutput("classification_to_add_epitope2")),
                                                                  ),
                                                 ),




                                                 ##### Classification to include ------
                                                 tabsetPanel(id = "Panel_TCRUMAP",
                                                             # T cell classification ------
                                                             tabPanel("GEX",
                                                                      tabsetPanel(id = "Panel_class",
                                                                                  tabPanel("Percentage",value = 16,
                                                                                           add_busy_spinner(spin = "fading-circle"),
                                                                                           div(DT::dataTableOutput("Percent_tab")),
                                                                                           downloadButton('downloaddf_Percent_tab','Download table')
                                                                                  ),
                                                                                  tabPanel("UMAP plot",value = 14,
                                                                                           fluidRow(column(3,selectInput("show_selected","Show all labels?",choices=c("All","Selected_list"))),
                                                                                                    column(9,uiOutput("SiteNumInput",width = "900px")),
                                                                                           ),
                                                                                           add_busy_spinner(spin = "fading-circle"),
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
                                                                                           add_busy_spinner(spin = "fading-circle"),
                                                                                           fluidRow(
                                                                                             column(3,numericInput("strip_size","Size of header",value = 10)),
                                                                                             column(9,p("Fill = Select function type; Group = Select cluster"))
                                                                                           ),
                                                                                           fluidRow(column(3,
                                                                                                           wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                                                     uiOutput('myPanel_pie'))), # pie chart
                                                                                                    column(9, plotOutput("Classification_clonotype_pie",height="600px"))),

                                                                                           fluidRow(
                                                                                             column(1,numericInput("width_Classification_clonotype_pie", "Width of PDF", value=10)),
                                                                                             column(1,numericInput("height_Classification_clonotype_pie", "Height of PDF", value=8)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_Classification_clonotype_pie','Download PDF')),
                                                                                             column(2,numericInput("width_png_Classification_clonotype_pie","Width of PNG", value = 1200)),
                                                                                             column(2,numericInput("height_png_Classification_clonotype_pie","Height of PNG", value = 1000)),
                                                                                             column(2,numericInput("resolution_PNG_Classification_clonotype_pie","Resolution of PNG", value = 144)),
                                                                                             column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_Classification_clonotype_pie','Download PNG'))),

                                                                                  ),
                                                                                  tabPanel("Longitudinal (χ2) ",
                                                                                           h5("balloonplot of χ2"),
                                                                                           p("The plots include the pearsons residuals for each cell (aka standardized residuals), for the proportion of cell contribution. The Contrib = contribution in percentage (%)"),p("See:",tags$a(href="http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r", " for Information")),


                                                                                           selectInput("type_res","Type of comparison",choices = c("Residuals","Contrib")),
                                                                                           conditionalPanel(condition="input.type_res=='Residuals'",
                                                                                                            fluidRow(
                                                                                                              column(3,colourInput("lower_col_chi","Lower colour:","purple")),
                                                                                                              column(3,colourInput("mid_col_chi","Mid colour:","White")),
                                                                                                              column(3,colourInput("high_col_chi","Higher colour:","Gold")),
                                                                                                            )),

                                                                                           conditionalPanel(condition="input.type_res=='Contrib'",
                                                                                                            fluidRow(
                                                                                                              column(3,colourInput("lower_col_chi2","Lower colour:","purple")),
                                                                                                              # column(3,colourInput("mid_col_chi","Mid colour:","White")),
                                                                                                              column(3,colourInput("high_col_chi2","Higher colour:","Gold")),
                                                                                                            )),

                                                                                           tabsetPanel(
                                                                                             tabPanel("Table",
                                                                                                      div(DT::dataTableOutput("Chi_tab_before")),
                                                                                             ),
                                                                                             tabPanel("balloonplot",
                                                                                                      plotOutput("Chi_square_plot",height="600px"),
                                                                                                      fluidRow(
                                                                                                        column(1,numericInput("width_Chi_square_plot", "Width of PDF", value=10)),
                                                                                                        column(1,numericInput("height_Chi_square_plot", "Height of PDF", value=8)),
                                                                                                        column(2,style = "margin-top: 25px;",downloadButton('downloadPlot_Chi_square_plot','Download PDF')),
                                                                                                        column(2,numericInput("width_png_Chi_square_plot","Width of PNG", value = 1200)),
                                                                                                        column(2,numericInput("height_png_Chi_square_plot","Height of PNG", value = 1000)),
                                                                                                        column(2,numericInput("resolution_PNG_Chi_square_plot","Resolution of PNG", value = 144)),
                                                                                                        column(2,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_Chi_square_plot','Download PNG'))),
                                                                                             )
                                                                                           ),
                                                                                  )

                                                                      ),
                                                             ),
                                                             # top clonotypes plot -----
                                                             tabPanel("Top clonotypes", value = "top_clone",
                                                                      # fluidRow(column(3,selectInput("View_AG_or_BD_top_clone","Select top chain:",choices = c("Alpha or Gamma","Beta or Delta"))),
                                                                      #
                                                                      #          ),
                                                                      tabsetPanel(
                                                                        tabPanel("Table",
                                                                                 add_busy_spinner(spin = "fading-circle"),
                                                                                 div(DT::dataTableOutput("Top_clonotype_df")),

                                                                        ),
                                                                        tabPanel("Bar graph",
                                                                                 add_busy_spinner(spin = "fading-circle"),

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
                                                                        tabPanel("Pie/UMAP chart",

                                                                                 add_busy_spinner(spin = "fading-circle"),
                                                                                 fluidRow(
                                                                                   column(3,selectInput("Plot_type_selected","Plot",choices = c("pie","UMAP"))),
                                                                                   column(3,numericInput("strip_size_pie","Size of header",value = 10)),
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
                                                                        tabPanel("Expression",
                                                                                 actionButton("run_string.data_Exp_top","View Ridge plot"),
                                                                                 fluidRow(column(12, selectInput("string.data_Exp_top","column names for summary","",multiple = F, width = "1200px") )),
                                                                                 tabsetPanel(
                                                                                   tabPanel("Plots",
                                                                                            div(DT::dataTableOutput("Violin_chart_alpha_gamma")),
                                                                                            fluidRow(
                                                                                              column(3, checkboxInput("restric_ex","Restrict to above a threshold?", value = F )),
                                                                                              column(3, numericInput("Gre_ex","Expression above:", value = 0 )),
                                                                                              column(3, selectInput("plot_type_ridgvi","Plot type", choices = c("Ridge (selected clonotype)","Ridge (compare)","Violin (selected clonotype)", "Violin (compare)"))),
                                                                                            ),

                                                                                            add_busy_spinner(spin = "fading-circle"),
                                                                                            plotOutput("Ridge_chart_alpha_gamma_plot_out",height="600px"),
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
                                                                                            add_busy_spinner(spin = "fading-circle"),
                                                                                            h5("Left = relative to all; Right based on function"),
                                                                                            fluidRow(
                                                                                              column(6, div(DT::dataTableOutput("Ridge_chart_alpha_gamma_stat"))),
                                                                                              column(6, div(DT::dataTableOutput("Ridge_chart_alpha_gamma_stat_comp")))
                                                                                            ),

                                                                                            downloadButton('downloaddf_clusTCR_GEx','Download stat (left)')
                                                                                   ),
                                                                                 ),



                                                                        ),
                                                                      ),
                                                             ),
                                                             # clusTCR specific analysis
                                                             # tabPanel("clusTCR"),
                                                             # epitope analysis -----
                                                             tabPanel("Epitope", value = "Epitope",
                                                                      tabsetPanel(

                                                                        tabPanel("Uploaded Epitope file",
                                                                                 add_busy_spinner(spin = "fading-circle"),
                                                                                 div(DT::dataTableOutput("MainTcell_Check")),
                                                                        ),

                                                                        tabPanel("Summary Table",
                                                                                 add_busy_spinner(spin = "fading-circle"),
                                                                                 div(DT::dataTableOutput("Pie_Epitope_dt")),
                                                                                 downloadButton('downloaddf_Pie_Epitope_dt','Download table')


                                                                        ),

                                                                        tabPanel("Heatmap",
                                                                                 fluidRow(
                                                                                   column(3,selectInput("epitope_hm","y-axis",choices = c("CDR3_beta","epitope","pathology"),selected = "epitope")),
                                                                                   column(3,selectInput("pathology_hm","x-axis",choices = c("CDR3_beta","epitope","pathology"),selected = "pathology")),
                                                                                 ),

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
                                                                        )
                                                                      ),

                                                             ),
                                                             # ClusTCR2 Analysis -----
                                                             tabPanel("ClusTCR2",value = "ClusTCR2",
                                                                      tabsetPanel(
                                                                        tabPanel("Table",

                                                                                 add_busy_spinner(spin = "fading-circle"),
                                                                                 div(DT::dataTableOutput("Tb_ClusTCR_selected")),

                                                                        ),
                                                                        tabPanel("UMAP",
                                                                                 fluidRow(
                                                                                   column(3,selectInput("ClusTCR_display","Colour by:",choices = c("all","Selected"))),
                                                                                   column(3,conditionalPanel(condition="input.ClusTCR_display=='Selected'",selectInput("Clusters_to_dis","Clusters to display",
                                                                                                                                                                       choices = "",multiple = T))),

                                                                                   column(3,selectInput("Clusters_to_dis_motif","Clusters to display (motif)",
                                                                                                        choices = "",multiple = F))
                                                                                 ),
                                                                                 fluidRow(
                                                                                   # column(3,
                                                                                   #        wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                   #                  uiOutput('myPanel_cols_epitope'))),
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
                                                                        )
                                                                      )

                                                             ),
                                                             # Overlap ------
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
                                                                        tabPanel("UMAP")
                                                                      )
                                                             ),
                                                             # tabPanel("Clonotypes per cluster (Pie/bar plot)"),
                                                             # tabPanel("Upset plot")
                                                 ),
                                        ),
                                        ### Longitudinal  -----
                                        #
                                        # tabPanel("Longitudinal data",
                                        #          tabsetPanel(
                                        #            tabPanel("Gene expression",
                                        #                     tabsetPanel(
                                        #                       tabPanel("Violine Plot"),
                                        #                       tabPanel("Ridge Plot"),
                                        #                       tabPanel("Stat")
                                        #                     )),
                                        #
                                        #            tabPanel("plot_trajectories"),
                                        #            # Module Expression
                                        #            # possibely represented as a heatmap? Downloadable table for pathway analysis
                                        #            tabPanel("Series of treemaps/chord diagrams?"),
                                        #            tabPanel("Correlation plot?"),
                                        #            tabPanel("Shannon_entropy vs time"), # changes in diveristy
                                        #            tabPanel("Module")
                                        #          )
                                        # ),
                            )
                          )
                        )
               ),
               # end of Integration -----
               # tabPanel("Epitope")
    ) # nav page
  )

  ########
  # server ------
  server <- function(input, output,session) {
    # add UI ------
    output$classification_to_add <- renderUI({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )

      df3.meta <- c(names(sc@meta.data))
      df3.meta <- df3.meta[!grepl("RNA",df3.meta) & !grepl("BCR",df3.meta) & !grepl("TCR",df3.meta)& !grepl("_gene",df3.meta) & !grepl("allele",df3.meta) & !grepl("percent",df3.meta) & !grepl("cdr3",df3.meta)]
      selectInput("clust_names_top","Select function type:",choices = df3.meta,selected="seurat_clusters")

    })


    output$classification_to_add_epitope <- renderUI({
      df3.meta <-  Add.UMAP.reduction()
      # summarising the clonality
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(df3.meta)>0 & nrow(epi)>0,
             "Upload Files")
      )
      if(input$datasource == "BD rhapsody") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

      }
      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta,epi,by="CDR3_beta")

      selectInput("epitope_umap_selected","Colour Pie by (hm = y-axis):",choices = names(df3.meta),selected="beta")

    })


    output$classification_to_add_epitope2 <- renderUI({
      df3.meta <-  Add.UMAP.reduction()
      # summarising the clonality
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(df3.meta)>0 & nrow(epi)>0,
             "Upload Files")
      )
      if(input$datasource == "BD rhapsody") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

      }
      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta,epi,by="CDR3_beta")

      selectInput("epitope_umap_selected2","Split Pie by (hm = x-axis):",choices = names(df3.meta),selected="beta")

    })

    # column(3,selectInput("epitope_umap_selected","Select",choices = c("beta","epitope","pathology"),selected = "pathology")),

    # user interface parameters-----
    output$feature_input <- renderUI({
      if (input$df_seruatobj_type=="10x Genomics") {
        fluidRow(
          column(6,numericInput("features.min","minimum features (<)", value = 200)),
          column(6,numericInput("features.max","Maximum features (<)", value = 6000)),
          column(6,numericInput("percent.mt","Mictochondrial DNA cut-off (<)", value = 20)),
          column(6,numericInput("percent.rb","Ribosomal RNA cut-off (>)", value = 5)),
        )
      }



      else if (input$df_seruatobj_type=="BD Rhapsody") {
        fluidRow(
          column(6,numericInput("features.min","minimum features (<)", value = 45)),
          column(6,numericInput("features.max","Maximum features (<)", value = 160)),
          column(6,numericInput("percent.mt","Mictochondrial DNA cut-off (<)", value = 20)),
          column(6,numericInput("percent.rb","Ribosomal RNA cut-off (>)", value = 0)),
        )
      }

      else {
        fluidRow(
          column(6,numericInput("features.min","minimum features (<)", value = 200)),
          column(6,numericInput("features.max","Maximum features (<)", value = 6000)),
          column(6,numericInput("percent.mt","Mictochondrial DNA cut-off (<)", value = 20)),
          column(6,numericInput("percent.rb","Ribosomal RNA cut-off (>)", value = 0)),
        )
      }

    })
    # human BD rhapsody data -----
    ## three files required for BD data: Sample Tag calls, TCR file and count ----
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

    ## combining all data ----
    df <- function() {
      calls <- input.data.calls.bd();
      TCR <- as.data.frame(input.data.TCR.BD())
      counts <- as.data.frame(input.data.count.BD())

      validate(
        need(nrow(counts)>0 & nrow(TCR)>0 & nrow(calls)>0,
             error_message_val4)
      )

      calls_TCR <- merge(calls,TCR, by ="Cell_Index")
      calls_TCR_count <- merge(calls_TCR,counts, by ="Cell_Index")
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
      if (input$filter_non_function ==T && input$BCR_present ==T) {
        productive <- calls_TCR_count[calls_TCR_count$Productive_TCR %in% "productive TCR" | calls_TCR_count$Productive_BCR %in% "productive BCR",] }
      else if (input$filter_non_function ==T && input$BCR_present ==F) {

        productive <- calls_TCR_count[calls_TCR_count$Productive_TCR %in% "productive TCR",]
      }
      else if (input$filter_non_function ==F && input$BCR_present ==T) {
        productive <- calls_TCR_count
      }
      else {
        productive <- calls_TCR_count
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
    output$Filtering_BD <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      df()


    })

    ## create Sample_tags file =====

    samp.tags <- reactive({
      TCR <- as.data.frame(input.data.TCR.BD())
      validate(
        need(nrow(TCR)>0,
             error_message_val2)
      )
      Sample_Tags <- as.data.frame(TCR$Cell_Index)
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
        paste(input$name.BD,"_Sample_Tags_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(samp.tags())
        write.csv(df,file, row.names = F)
      } )

    ## for the clusTCR ----
    df_clusTCR <- function () {
      calls_TCR_paired.fun <- df()
      counts <- as.data.frame(input.data.count.BD())
      df_nocouts <- calls_TCR_paired.fun[!names(calls_TCR_paired.fun) %in% c(names(counts),"cloneCount") ]
      df_nocouts2 <- df_nocouts[df_nocouts$Total_VDJ_Read_Count != 0, ]
      df_nocouts2_AG <- df_nocouts2[,names(df_nocouts2) %in% c("Sample_Name","TCR_Alpha_Gamma_CDR3_Translation_Dominant","TCR_Alpha_Gamma_V_gene_Dominant")]
      df_nocouts2_BD <- df_nocouts2[,names(df_nocouts2) %in% c("Sample_Name","TCR_Beta_Delta_CDR3_Translation_Dominant","TCR_Beta_Delta_V_gene_Dominant")]
      names(df_nocouts2_AG) <- c("Sample_Name","v_call","junction_aa")
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

      df_nocouts3 <- df_nocouts3[-c(grep("^$",df_nocouts3$junction_aa)),] # remove blank columns
      df_nocouts3 <- df_nocouts3[-c(grep("[*]",df_nocouts3$junction_aa)),] # remove st

      df_nocouts3

    }
    output$tb_clusTCR <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      df_clusTCR()
    })
    output$downloaddf_clusTCR <- downloadHandler(
      filename = function(){
        paste(input$name.BD," ClusTCR ",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(df_clusTCR())
        write.csv(df,file, row.names = T)
      } )

    ### summary of missing TCR -----
    sum_tb_BD <- function () {
      calls_TCR_paired.fun <- df()
      counts <- as.data.frame(input.data.count.BD())
      df_nocouts <- calls_TCR_paired.fun[!names(calls_TCR_paired.fun) %in% c(names(counts),"cloneCount") ]
      contig_paired <- df_nocouts[df_nocouts$Total_VDJ_Read_Count != 0, ]
      contig_paired$TRAG <-ifelse(grepl("TR",contig_paired$TCR_Alpha_Gamma_V_gene_Dominant),"TRAG","Missing TRAG gene")
      contig_paired$TRBD <-ifelse(grepl("TR",contig_paired$TCR_Beta_Delta_V_gene_Dominant),"TRBD","Missing TRBD gene")

      contig_paired$TRAG_fun <-ifelse(grepl("[*]",contig_paired$TCR_Alpha_Gamma_CDR3_Translation_Dominant),"Non-functional",
                                      ifelse(grepl("Missing",contig_paired$TRAG),"Missing TCR",

                                             ifelse(grepl("TRGV10",contig_paired$TCR_Alpha_Gamma_V_gene_Dominant),"pseudogene","productive")))

      contig_paired$TRBD_fun <- ifelse(grepl("[*]",contig_paired$TCR_Beta_Delta_CDR3_Translation_Dominant),"Non-functional",
                                       ifelse(grepl("Missing",contig_paired$TRBD),"Missing TCR","productive"))

      if (input$BCR_present ==T) {

        contig_paired$v_gene_IgL <- ifelse(grepl("IG",contig_paired$BCR_Light_V_gene_Dominant),"v_gene_IgL","Missing IgL gene")

        contig_paired$v_gene_IgH <- ifelse(grepl("IG",contig_paired$BCR_Heavy_V_gene_Dominant),"v_gene_IgL","Missing IgL gene")

        contig_paired$v_gene_IgL_fun <- ifelse(grepl("[*]",contig_paired$BCR_Light_CDR3_Translation_Dominant),"Non-functional",
                                               ifelse(grepl("Missing",contig_paired$v_gene_IgL),"Missing BCR","productive"))
        contig_paired$v_gene_IgH_fun <- ifelse(grepl("[*]",contig_paired$BCR_Heavy_CDR3_Translation_Dominant),"Non-functional",
                                               ifelse(grepl("Missing",contig_paired$v_gene_IgH),"Missing BCR","productive"))

        df.names.test <- c('TRAG',"TRAG_fun","TRBD","TRBD_fun","v_gene_IgL","v_gene_IgL_fun","v_gene_IgH","v_gene_IgH_fun")
        df3 <- contig_paired[names(contig_paired) %in% df.names.test]
        df3$cloneCount <- 1
        df1 <- ddply(df3,names(df3)[-c(9)] ,numcolwise(sum))
        df1


      }
      else {
        df.names.test <- c('TRAG',"TRAG_fun","TRBD","TRBD_fun","cloneCount")
        df3 <- contig_paired[names(contig_paired) %in% df.names.test]
        df3$cloneCount <- 1
        df1 <- ddply(df3,names(df3)[-c(5)] ,numcolwise(sum))
        df1
      }
    }

    sum_tb_BD_2 <- function () {
      calls_TCR_paired.fun <- df()
      counts <- as.data.frame(input.data.count.BD())
      df_nocouts <- calls_TCR_paired.fun[!names(calls_TCR_paired.fun) %in% c(names(counts),"cloneCount") ]
      contig_paired <- df_nocouts[df_nocouts$Total_VDJ_Read_Count != 0, ]
      contig_paired$TRAG <-ifelse(grepl("TR",contig_paired$TCR_Alpha_Gamma_V_gene_Dominant),"TRAG","Missing TRAG gene")
      contig_paired$TRBD <-ifelse(grepl("TR",contig_paired$TCR_Beta_Delta_V_gene_Dominant),"TRBD","Missing TRBD gene")

      contig_paired$TRAG_fun <-ifelse(grepl("[*]",contig_paired$TCR_Alpha_Gamma_CDR3_Translation_Dominant),"Non-functional",
                                      ifelse(grepl("Missing",contig_paired$TRAG),"Missing TCR",

                                             ifelse(grepl("TRGV10",contig_paired$TCR_Alpha_Gamma_V_gene_Dominant),"pseudogene","productive")))

      contig_paired$TRBD_fun <- ifelse(grepl("[*]",contig_paired$TCR_Beta_Delta_CDR3_Translation_Dominant),"Non-functional",
                                       ifelse(grepl("Missing",contig_paired$TRBD),"Missing TCR","productive"))

      if (input$BCR_present ==T) {

        contig_paired$v_gene_IgL <- ifelse(grepl("IG",contig_paired$BCR_Light_V_gene_Dominant),"v_gene_IgL","Missing IgL gene")

        contig_paired$v_gene_IgH <- ifelse(grepl("IG",contig_paired$BCR_Heavy_V_gene_Dominant),"v_gene_IgL","Missing IgL gene")

        contig_paired$v_gene_IgL_fun <- ifelse(grepl("[*]",contig_paired$BCR_Light_CDR3_Translation_Dominant),"Non-functional",
                                               ifelse(grepl("Missing",contig_paired$v_gene_IgL),"Missing BCR","productive"))
        contig_paired$v_gene_IgH_fun <- ifelse(grepl("[*]",contig_paired$BCR_Heavy_CDR3_Translation_Dominant),"Non-functional",
                                               ifelse(grepl("Missing",contig_paired$v_gene_IgH),"Missing BCR","productive"))

        contig_paired
      }
      else {
        contig_paired
      }

    }

    output$tb_clusTCR_sum <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      sum_tb_BD()
    })


    ## Differential expression files ----
    df_count.matrix <- function () {
      calls_TCR_paired.fun <- df()
      counts <- as.data.frame(input.data.count.BD())

      calls_TCR_paired.fun$v_gene_AG <- gsub("[*]0.","",calls_TCR_paired.fun$TCR_Alpha_Gamma_V_gene_Dominant)
      calls_TCR_paired.fun$v_gene_BD <- gsub("[*]0.","",calls_TCR_paired.fun$TCR_Beta_Delta_V_gene_Dominant)
      calls_TCR_paired.fun$cdr3_AG <- paste(calls_TCR_paired.fun$v_gene_AG,calls_TCR_paired.fun$TCR_Alpha_Gamma_CDR3_Translation_Dominant,sep = "_")
      calls_TCR_paired.fun$cdr3_BD <- paste(calls_TCR_paired.fun$v_gene_BD,calls_TCR_paired.fun$TCR_Beta_Delta_CDR3_Translation_Dominant,sep = "_")
      calls_TCR_paired.fun$v_gene_cdr3_AB_GD <- paste(calls_TCR_paired.fun$cdr3_AG,calls_TCR_paired.fun$cdr3_BD,sep = ".")

      rownames(calls_TCR_paired.fun) <- paste(calls_TCR_paired.fun$Sample_Name,"*",calls_TCR_paired.fun$cdr3_AG.cdr3_BD," ",calls_TCR_paired.fun$Cell_Index,sep="")
      df_nocouts <- calls_TCR_paired.fun[names(calls_TCR_paired.fun) %in% c(names(counts))]
      rownames(df_nocouts) <- df_nocouts$Cell_Index
      df_nocouts
    }

    output$tb_count_matrix <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      as.data.frame(t(df_count.matrix()))[1:6,1:6]
    })
    output$downloadtb_count_matrix <- downloadHandler(
      filename = function(){
        paste(input$name.BD," BD_Count_Matrix_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(t(df_count.matrix()))
        # write.table(,file, row.names = T)
        write.csv(df, file, row.names = T)

      })

    ### meta.data for seruat ----

    meta.data_for_Seuratobj <- function () {

      calls_TCR_paired.fun <- df()
      counts <- input.data.count.BD()

      calls_TCR_paired.fun$v_gene_AG <- gsub("[*]0.","",calls_TCR_paired.fun$TCR_Alpha_Gamma_V_gene_Dominant)
      calls_TCR_paired.fun$j_gene_AG <- gsub("[*]0.","",calls_TCR_paired.fun$TCR_Alpha_Gamma_J_gene_Dominant)
      calls_TCR_paired.fun$v_gene_BD <- gsub("[*]0.","",calls_TCR_paired.fun$TCR_Beta_Delta_V_gene_Dominant)
      calls_TCR_paired.fun$v_gene_AG_BD <- paste(calls_TCR_paired.fun$v_gene_AG,calls_TCR_paired.fun$v_gene_BD,sep = " & ")
      calls_TCR_paired.fun$v_gene_cdr3_AG <- paste(calls_TCR_paired.fun$v_gene_AG,calls_TCR_paired.fun$TCR_Alpha_Gamma_CDR3_Translation_Dominant,sep = "_")
      calls_TCR_paired.fun$v_gene_cdr3_BD <- paste(calls_TCR_paired.fun$v_gene_BD,calls_TCR_paired.fun$TCR_Beta_Delta_CDR3_Translation_Dominant,sep = "_")
      calls_TCR_paired.fun$v_gene_cdr3_AB_GD <- paste(calls_TCR_paired.fun$v_gene_cdr3_AG,calls_TCR_paired.fun$v_gene_cdr3_BD,sep = ".")

      calls_TCR_paired.fun <- calls_TCR_paired.fun[!names(calls_TCR_paired.fun) %in% c(names(counts[-c(1)]),"Total_VDJ_Read_Count","Total_VDJ_Molecule_Count","TCR_Alpha_Gamma_Read_Count","TCR_Alpha_Gamma_Molecule_Count","TCR_Beta_Delta_Read_Count","TCR_Beta_Delta_Molecule_Count","Sample_Tag","TCR_Alpha_Gamma_C_gene_Dominant","TCR_Beta_Delta_C_gene_Dominant","TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant","TCR_Beta_Delta_CDR3_Nucleotide_Dominant") ]

      if (input$BCR_present ==T) {

        calls_TCR_paired.fun$v_gene_IgL <- gsub("[*]0.","",calls_TCR_paired.fun$BCR_Light_V_gene_Dominant)
        calls_TCR_paired.fun$j_gene_IgL <- gsub("[*]0.","",calls_TCR_paired.fun$BCR_Light_J_gene_Dominant)
        calls_TCR_paired.fun$v_gene_IgH <- gsub("[*]0.","",calls_TCR_paired.fun$BCR_Heavy_V_gene_Dominant)
        calls_TCR_paired.fun$d_gene_IgH <- gsub("[*]0.","",calls_TCR_paired.fun$BCR_Heavy_D_gene_Dominant)
        calls_TCR_paired.fun$j_gene_IgH <- gsub("[*]0.","",calls_TCR_paired.fun$BCR_Heavy_J_gene_Dominant)

        calls_TCR_paired.fun$v_gene_IgH_L <- paste(calls_TCR_paired.fun$v_gene_IgH,calls_TCR_paired.fun$v_gene_IgL,sep = " & ")

        calls_TCR_paired.fun$v_gene_cdr3_IgH <- paste(calls_TCR_paired.fun$v_gene_IgH,calls_TCR_paired.fun$BCR_Heavy_CDR3_Translation_Dominant,sep = "_")
        calls_TCR_paired.fun$v_gene_cdr3_IgL <- paste(calls_TCR_paired.fun$v_gene_IgL,calls_TCR_paired.fun$BCR_Light_CDR3_Translation_Dominant,sep = "_")
        calls_TCR_paired.fun$v_gene_cdr3_IgH_L <- paste(calls_TCR_paired.fun$v_gene_cdr3_IgH,calls_TCR_paired.fun$v_gene_cdr3_IgL,sep = ".")

        names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("BCR_Light_CDR3_Translation_Dominant")] <- "cdr3_IgL"
        names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("BCR_Light_V_gene_Dominant")] <- "v_allele_IgL"
        names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("BCR_Light_J_gene_Dominant")] <- "j_allele_IgL"
        names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("BCR_Heavy_CDR3_Translation_Dominant")] <- "cdr3_IgH"
        names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("BCR_Heavy_V_gene_Dominant")] <- "v_allele_IgH"
        names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("BCR_Heavy_J_gene_Dominant")] <- "j_allele_IgH"
        names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("BCR_Heavy_D_gene_Dominant")] <- "d_allele_IgH"
        calls_TCR_paired.fun

      }

      else {
        calls_TCR_paired.fun
      }

      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TCR_Alpha_Gamma_CDR3_Translation_Dominant")] <- "cdr3_AG"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TCR_Alpha_Gamma_J_gene_Dominant")] <- "j_gene_AG"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TCR_Alpha_Gamma_V_gene_Dominant")] <- "v_allele_AG"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TCR_Beta_Delta_CDR3_Translation_Dominant")] <- "cdr3_BD"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TCR_Beta_Delta_V_gene_Dominant")] <- "v_allele_BD"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TCR_Beta_Delta_J_gene_Dominant")] <- "j_gene_BD"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TCR_Beta_Delta_D_gene_Dominant")] <- "d_gene_BD"
      calls_TCR_paired.fun

      # calls_TCR_paired.fun %>%
      #   select(Cell_Type_Experimental, everything())
    }

    output$tb_metadata_sc <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      meta.data_for_Seuratobj()
    })

    output$downloadtb_metadata_sc <- downloadHandler(
      filename = function(){
        paste(input$name.BD,"_Meta.data_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(meta.data_for_Seuratobj())
        write.csv(df,file, row.names = F)
      } )


    # for TCRex output -----
    TCRex_BDrap_df <- function () {

      calls_TCR_paired.fun <- df()
      calls_TCR_paired.fun <- calls_TCR_paired.fun[-c(grep("[*]",calls_TCR_paired.fun$TCR_Beta_Delta_CDR3_Translation_Dominant)),] # remove st
      calls_TCR_paired.fun$TRBV_gene <- gsub("[*]0.","",calls_TCR_paired.fun$TCR_Beta_Delta_V_gene_Dominant)
      calls_TCR_paired.fun$CDR3_beta <- paste("C",calls_TCR_paired.fun$TCR_Beta_Delta_CDR3_Translation_Dominant,"F",sep="")
      calls_TCR_paired.fun$TRBJ_gene<- gsub("[*]0.","",calls_TCR_paired.fun$TCR_Beta_Delta_J_gene_Dominant)
      calls_TCR_paired.fun$cloneCount <- 1
      calls_TCR_paired.fun2 <- calls_TCR_paired.fun[,names(calls_TCR_paired.fun) %in% c("TRBV_gene","CDR3_beta","TRBJ_gene","cloneCount")]

      calls_TCR_paired.fun3 <- ddply(calls_TCR_paired.fun2,names(calls_TCR_paired.fun2)[-c(4)] ,numcolwise(sum))
      calls_TCR_paired.fun3 <- subset(calls_TCR_paired.fun3,calls_TCR_paired.fun3$CDR3_beta != "CF")

      calls_TCR_paired.fun3 <- calls_TCR_paired.fun3[!grepl("TRD", calls_TCR_paired.fun3$TRBV_gene),]
      calls_TCR_paired.fun3


    }

    output$tb_TCRex_BDrap_df <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      TCRex_BDrap_df()
    })

    output$downloaddf_TCRex_BDrap <- downloadHandler(
      filename = function(){
        paste(input$name.BD,"_TCRex_",gsub("-", ".", Sys.Date()),".tsv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(TCRex_BDrap_df())
        df <- df[,names(df) %in% c("TRBV_gene","CDR3_beta","TRBJ_gene")]
        write.table(df,file, row.names = F,sep="\t", quote = F)
      } )

    # count matrix download ----
    output$downloadtb_count_matrix <- downloadHandler(
      filename = function(){
        paste(input$name.BD,"_BD_Count_Matrix_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(t(df_count.matrix()))
        # write.table(,file, row.names = T)
        write.csv(df, file, row.names = T)

      })
    ## TCR_Explore compatible -----
    df_TCR_Explore <- function () {
      calls_TCR_paired.fun <- df()
      counts <- input.data.count.BD()
      calls_TCR_paired.fun$v_gene_AG <- gsub("[*]0.","",calls_TCR_paired.fun$TCR_Alpha_Gamma_V_gene_Dominant)
      calls_TCR_paired.fun$v_gene_BD <- gsub("[*]0.","",calls_TCR_paired.fun$TCR_Beta_Delta_V_gene_Dominant)
      calls_TCR_paired.fun$v_gene_AG_BD <- paste(calls_TCR_paired.fun$v_gene_AG,calls_TCR_paired.fun$v_gene_BD,sep = " & ")

      calls_TCR_paired.fun$v_gene_cdr3_AG <- paste(calls_TCR_paired.fun$v_gene_AG,calls_TCR_paired.fun$TCR_Alpha_Gamma_CDR3_Translation_Dominant,sep = "_")
      calls_TCR_paired.fun$v_gene_cdr3_BD <- paste(calls_TCR_paired.fun$v_gene_BD,calls_TCR_paired.fun$TCR_Beta_Delta_CDR3_Translation_Dominant,sep = "_")
      calls_TCR_paired.fun$v_gene_cdr3_AG.cdr3_BD <- paste(calls_TCR_paired.fun$v_gene_cdr3_AG,calls_TCR_paired.fun$v_gene_cdr3_BD,sep = ".")
      calls_TCR_paired.fun$cloneCount <- 1 # added for TCR_Explore (total number of cells with that gene and sequence)

      calls_TCR_paired.fun <- calls_TCR_paired.fun[!names(calls_TCR_paired.fun) %in% c(names(counts[-c(1)]),"Total_VDJ_Read_Count","Total_VDJ_Molecule_Count","TCR_Alpha_Gamma_Read_Count","TCR_Alpha_Gamma_Molecule_Count","TCR_Beta_Delta_Read_Count","TCR_Beta_Delta_Molecule_Count","Cell_Type_Experimental") ]

      if (input$BCR_present ==T) {
        calls_TCR_paired.fun$v_gene_IgL <- gsub("[*]0.","",calls_TCR_paired.fun$BCR_Light_V_gene_Dominant)
        calls_TCR_paired.fun$j_gene_IgL <- gsub("[*]0.","",calls_TCR_paired.fun$BCR_Light_J_gene_Dominant)
        calls_TCR_paired.fun$v_gene_IgH <- gsub("[*]0.","",calls_TCR_paired.fun$BCR_Heavy_V_gene_Dominant)
        calls_TCR_paired.fun$d_gene_IgH <- gsub("[*]0.","",calls_TCR_paired.fun$BCR_Heavy_D_gene_Dominant)
        calls_TCR_paired.fun$j_gene_IgH <- gsub("[*]0.","",calls_TCR_paired.fun$BCR_Heavy_J_gene_Dominant)

        calls_TCR_paired.fun$v_gene_IgH_L <- paste(calls_TCR_paired.fun$v_gene_IgH,calls_TCR_paired.fun$v_gene_IgL,sep = " & ")

        calls_TCR_paired.fun$v_gene_cdr3_IgH <- paste(calls_TCR_paired.fun$v_gene_IgH,calls_TCR_paired.fun$BCR_Heavy_CDR3_Translation_Dominant,sep = "_")
        calls_TCR_paired.fun$v_gene_cdr3_IgL <- paste(calls_TCR_paired.fun$v_gene_IgL,calls_TCR_paired.fun$BCR_Light_CDR3_Translation_Dominant,sep = "_")
        calls_TCR_paired.fun$v_gene_cdr3_IgH_L <- paste(calls_TCR_paired.fun$v_gene_cdr3_IgH,calls_TCR_paired.fun$v_gene_cdr3_IgL,sep = ".")

        names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("BCR_Light_CDR3_Translation_Dominant")] <- "cdr3_IgL"
        names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("BCR_Light_V_gene_Dominant")] <- "v_allele_IgL"
        names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("BCR_Light_J_gene_Dominant")] <- "j_allele_IgL"
        names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("BCR_Heavy_CDR3_Translation_Dominant")] <- "cdr3_IgH"
        names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("BCR_Heavy_V_gene_Dominant")] <- "v_allele_IgH"
        names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("BCR_Heavy_J_gene_Dominant")] <- "j_allele_IgH"
        names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("BCR_Heavy_D_gene_Dominant")] <- "d_allele_IgH"
        calls_TCR_paired.fun
      }

      else {
        calls_TCR_paired.fun
      }

      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TCR_Alpha_Gamma_CDR3_Translation_Dominant")] <- "cdr3_AG"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TCR_Alpha_Gamma_J_gene_Dominant")] <- "j_gene_AG"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TCR_Alpha_Gamma_V_gene_Dominant")] <- "v_allele_AG"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TCR_Beta_Delta_CDR3_Translation_Dominant")] <- "cdr3_BD"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TCR_Beta_Delta_V_gene_Dominant")] <- "v_allele_BD"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TCR_Beta_Delta_J_gene_Dominant")] <- "j_gene_BD"
      names(calls_TCR_paired.fun)[names(calls_TCR_paired.fun) %in% c("TCR_Beta_Delta_D_gene_Dominant")] <- "d_gene_BD"
      calls_TCR_paired.fun


      calls_TCR_paired.fun[order(calls_TCR_paired.fun$v_allele_AG),]
      calls_TCR_paired.fun %>%
        select(cloneCount, everything())
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


    # 10x genomics data -----
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
        dataframe <- read.csv(inFile_10x_TCR$datapath)}
    })
    output$test.files.10x4 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      calls <- input.data.TCR.10x()
      validate(
        need(nrow(calls)>0,
             error_message_val_10x_features)
      )
      calls
    })

    ## meta.data for seurat ----
    tb_10x_meta.data_TCR <- function () {
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
      contig_paired_only$d_gene_BD <- sub("^$","NA", contig_paired_only$d_gene_BD)
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
      contig_paired_only$vdj_gene_BD <- gsub(".None.",".",contig_paired_only$vdj_gene_BD)
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
      contig_paired_only$Sample_Name <- input$sample_name_10x

      contig_paired_only <- contig_paired_only %>%
        select(all_of(c("Cell_Index","Sample_Name")), everything())

      contig_paired_only
    }
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
      contig_paired_only$vdj_gene_cdr3_H <- paste(contig_paired_only$vdj_gene_H,contig_paired_only$cdr3_H,sep = "_")
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
      contig_paired_only_dup$Sample_Name <- input$sample_name_10x

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
        TCR <- tb_10x_meta.data_TCR()
        BCR <- tb_10x_meta.data_BCR()
        contig_paired <- names(BCR)[!grepl("_IgH",names(BCR)) & !grepl("_IgL",names(BCR))]
        contig_paired <- merge(TCR,BCR,by = merge.names, all=T)
        contig_paired
      }

    })
    output$downloadtb_10x_metadata2 <- downloadHandler(
      filename = function(){
        paste(input$name.10x,"_metadata_10x_",gsub("-", ".", Sys.Date()),".csv", sep = "")
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
        write.csv(df, file)

      } )

    # meta-data summary
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
      contigs2 <- contigs[names(contigs) %in% c("v_gene","cdr3")]
      contigs2$Sample_Name <- input$sample_name_10x
      names(contigs2)[names(contigs2) %in% c("cdr3")] <- "junction_aa"
      names(contigs2)[names(contigs2) %in% c("v_gene")] <- "v_call"
      contigs2 <- contigs2[-c(grep("[*]",contigs2$junction_aa)),]
      contigs2 <- subset(contigs2,contigs2$junction_aa!= "None")
      contigs2[!duplicated(contigs2[,c('v_call','junction_aa')]),]
    })

    output$tb_10x_contigues1 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      tb_10x_contigues_contig()
    })

    output$downloadtb_10x_contigues1 <- downloadHandler(
      filename = function(){
        paste(input$name.10x,"_clusTCR_10x_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(tb_10x_contigues_contig())
        write.csv(df,file, row.names = F)
      } )


    ## TCR explore 10x -----
    TCR_Explore_10x <- function () {
      if (input$BCR_TCR_10x=="TCR only") {
        contig_paired_only <-  tb_10x_meta.data_TCR()

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
        paste(input$name.10x,"_TCR_Explore_10x_",gsub("-", ".", Sys.Date()),".csv", sep = "")
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

    output$tb_10x_matrix2 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 2, scrollX = TRUE),{
      tb_10x_matrix()[1:6,1:6]
    })

    output$downloadtb_10x_matrix2 <- downloadHandler(
      filename = function(){
        paste(input$name.10x,"_count-matrix_10x_",gsub("-", ".", Sys.Date()),".csv.gz", sep = "")
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
      contigs2 <- contigs[,names(contigs) %in% c("v_gene","j_gene","cdr3")]

      names(contigs2)[names(contigs2) %in% c("cdr3")] <- "CDR3_beta"
      contigs2 <- contigs2[-c(grep("[*]",contigs2$CDR3_beta)),]

      names(contigs2)[names(contigs2) %in% c("v_gene")] <- "TRBV_gene"
      names(contigs2)[names(contigs2) %in% c("j_gene")] <- "TRBJ_gene"
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
    }

    output$tb_TCRex_10x_df <- DT::renderDataTable(filter = list(position = 'top', clear = FALSE), escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      TCRex_10x_df()
    })

    output$downloaddf_TCRex_10x <- downloadHandler(
      filename = function(){
        paste(input$name.10x," TCRex ",gsub("-", ".", Sys.Date()),".tsv", sep = "")
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

      df_nocouts3 <- df_nocouts3[-c(grep("[*]",df_nocouts3$junction_aa)),] # remove st
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
        paste(input$name.array,"_count-matrix_array_",gsub("-", ".", Sys.Date()),".csv", sep = "")
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
      myfiles <- input.data_ClusTCR2_multiple()
      df <- rbind(myfiles[[1]])
      for (i in 2:length(myfiles)) {
        df <- rbind(df,myfiles[[i]])
      }
      # df <- df[-c(grep("[*]",df$junction_aa)),]
      df
      df[!duplicated(df$junction_aa,df$v_call),]
    })

    output$DEx_multiple_ClusTCR2 <-  DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      df <- merged_Clust_filtered()
      # validate(
      #   need(nrow(df)>0,
      #        "Upload mutliple CLusTCR2 files to merge")
      # )
      df
    })

    output$downloaddf_multiple_ClusTCR2 <- downloadHandler(
      filename = function(){
        paste("Multi_ClusTCR",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- merged_Clust_filtered()
        write.csv(df,file, row.names = F)
      })

    # rendering the table

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
    # observe({
    #   updateSelectInput(
    #     session,
    #     "clusTCR2_Group_column",
    #     choices=names(input.data()),
    #     selected = "group")
    #
    # }) # group

    clust_dt2 <- reactive({
      df1 <- input.data_ClusTCR2()
      validate(
        need(nrow(df1)>0,
             "Upload ClusTCR file")
      )

      df1
      df1 <- df1 %>%
        select(input$clusTCR2_names,input$clusTCR2_Vgene,everything())

      df2 <- df1
      names(df2)[1:2] <- c("V1","V2")
      df2$V1 <- ifelse(grepl("^$",df2$V1),"Missing TCR",df2$V1)
      df2$V2 <- ifelse(grepl("^$",df2$V2),"Missing TCR",df2$V2)
      head(df2)
      df3 <- subset(df2, df2$V1 != "Missing TCR" | df2$V2 != "Missing TCR")
      head(df3)
      names(df3)[1:2] <- c(input$clusTCR2_names,input$clusTCR2_Vgene)
      df3
    })
    output$clust_dt2 <- DT::renderDataTable({
      clust_dt2()
    })

    ## run clustering
    vals_ClusTCR2 <- reactiveValues(output_dt2=NULL)
    observeEvent(input$run_ClusTCR2,{
      clust_dt2 <- clust_dt2()
      validate(
        need(nrow(clust_dt2)>0,
             "Run ClusTCR2")
      )
      ptm <- proc.time()

      # clust_dt2 <- clust_dt2()

      df_cluster <- ClusTCR(clust_dt2, allele =input$allele_ClusTCR2, v_gene = input$clusTCR2_Vgene, cores_selected = input$cores_ClusTCR2)
      cluster_lab <- mcl_cluster(df_cluster,expansion = 1,inflation = 1)
      end <- proc.time() - ptm
      cluster_lab[[3]] <- end
      vals_ClusTCR2$output_dt2 <- cluster_lab
    })

    ClusTCR2_lab_df <- reactive({
      output_dt <- vals_ClusTCR2$output_dt2
      output_dt[[1]]

    })

    output$ClusTCR2_lab <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      ClusTCR2_lab_df()
    })

    output$download_ClusTCR_labels <- downloadHandler(
      filename = function(){
        paste("ClusTCR2_output",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(ClusTCR2_lab_df())
        write.csv(df,file, row.names = F)
      } )

    output$ClusTCR2_Time <- renderPrint({
      output_dt <- vals_ClusTCR2$output_dt2

      if (is.null(output_dt)){return(NULL)}
      output_dt[[3]]
    })


    # create plots
    Network_plot_clusTCR2 <- reactive({
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
                       non_selected_text_col = "black")

    })
    output$NP_ClusTCR <- renderPlot({
      Network_plot_clusTCR2()
    })

    motif_plot_clusTCR2 <- reactive({
      Network_df <- vals_ClusTCR2$output_dt2
      set.seed(123)
      motif_plot(Network_df,Clust_selected = input$selected_Cluster,Clust_column_name = input$Clust_size_order)

    })

    output$MP_ClusTCR <- renderPlot({
      motif_plot_clusTCR2()
    })

    ## downloading plot -----
    output$downloadPlot_Network_plot2 <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste("clusTCR2_Network_plot_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Network_plot2,height=input$height_Network_plot2, onefile = FALSE) # open the pdf device
        plot(Network_plot_clusTCR2())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Network_plot2 <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste("clusTCR2_Network_plot_", gsub("/", "-", x), ".png", sep = "")
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
        x <- gsub(":", ".", Sys.time())
        paste("clusTCR2_Motif_plot"," Cluster ",input$selected_Cluster," ",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Motif_plot2,height=input$height_Motif_plot2, onefile = FALSE) # open the pdf device
        plot(motif_plot_clusTCR2())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Motif_plot2 <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste("clusTCR2_Motif_plot"," Cluster ",input$selected_Cluster," ", gsub("/", "-", x), ".png", sep = "")
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
        if (input$df_seruatobj_type =="10x Genomics") {
          dataframe = read.csv(inFile_sc$datapath)
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
      if (input$df_seruatobj_type =="10x Genomics") {
        names(df.test) <- gsub("[.]1","-1",names(df.test) )
        rownames(df.test) <- make.unique(df.test$Gene_Name)
        df.test2 <- df.test[,!names(df.test) %in% c("Gene_Name")]
      }

      else if (input$df_seruatobj_type=="BD Rhapsody") {
        names(df.test) <- gsub("X","",names(df.test))
        df.test2 <- df.test[!rownames(df.test) %in% c("Cell_Index"),]
      }

      else {
        names(df.test) <- gsub("[.]","-",names(df.test))
        df.test2 <- df.test[!rownames(df.test) %in% c("Cell_Index"),]

      }
      head(df.test2)[1:6]
    })

    ## reading in 10x and BD data ----
    df_seruatobj <- reactive({
      df.test <- input.data_sc()
      validate(
        need(nrow(df.test)>0,
             error_message_val_sc)
      )

      if (input$df_seruatobj_type =="10x Genomics") {
        names(df.test) <- gsub("[.]1","-1",names(df.test) )
        rownames(df.test) <- make.unique(df.test$Gene_Name)
        df.test2 <- df.test[,!names(df.test) %in% c("Gene_Name")]
        sc <- CreateSeuratObject(counts = df.test2, project = input$project_name)
        sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
        sc[["percent.rb"]] <- PercentageFeatureSet(sc, pattern = "^RP[SL]")
        sc
      }
      else if (input$df_seruatobj_type=="BD Rhapsody") {
        names(df.test) <- gsub("X","",names(df.test))
        df.test2 <- df.test[!rownames(df.test) %in% c("Cell_Index"),]
        sc <- CreateSeuratObject(counts = df.test2, project = input$project_name)
        sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
        sc[["percent.rb"]] <- PercentageFeatureSet(sc, pattern = "^RP[SL]")
        sc
      }

      else {
        names(df.test) <- gsub("[.]","-",names(df.test))
        sc <- CreateSeuratObject(counts = df.test, project = input$project_name)
        sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
        sc[["percent.rb"]] <- PercentageFeatureSet(sc, pattern = "^RP[SL]")
        sc
      }
    })
    before_plot <- reactive({
      sc <- df_seruatobj()
      validate(
        need(nrow(sc)>0,
             "Upload files")
      )
      VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 2)
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
        x <- gsub(":", ".", Sys.time())
        paste(input$project_name2,"_before_plot_sc_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_before_plot_sc,height=input$height_before_plot_sc, onefile = FALSE) # open the pdf device
        plot(before_plot())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_before_plot_sc <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$project_name2,"_before_plot_sc_", gsub("/", "-", x), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_before_plot_sc, height = input$height_png_before_plot_sc, res = input$resolution_PNG_before_plot_sc)
        plot(before_plot())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    ##### after filtering plot


    vals2 <- reactiveValues(after_violin_plot=NULL)
    observeEvent(input$run,{
      sc <- df_seruatobj()
      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )
      vals2$after_violin_plot <- subset(sc, subset = nFeature_RNA >= input$features.min & nFeature_RNA <= input$features.max & percent.mt <= input$percent.mt & percent.rb >= input$percent.rb)
    })

    output$after_plot_sc <- renderPlot({
      sc <- vals2$after_violin_plot
      validate(
        need(nrow(sc)>0,
             "Run Filtering")
      )
      VlnPlot(vals2$after_violin_plot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 2)
    })

    output$downloadPlot_Before.after_plot_sc <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste("inv.simpson.index.",gsub("/", "-", x), ".pdf", sep = "")
      }, content = function(file) {
        plot1 <- before_plot()
        # plot2 <-
        pdf(file, width=input$width_Before.after_plot_sc,height=input$height_Before.after_plot_sc, onefile = FALSE) # open the pdf device
        print(VlnPlot(vals2$after_violin_plot, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 2))
        # group.diversity1()
        dev.off()},
      contentType = "application/pdf" )

    ### normalisationa and feature plot ------
    feature_serartobj <- reactive({
      sc <- vals2$after_violin_plot
      validate(
        need(nrow(sc)>0,
             "Run Clustering")
      )
      sc <- NormalizeData(sc)
      sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
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
        x <- gsub(":", ".", Sys.time())
        paste(input$project_name,"_10_features_sc",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_plot_10_features_sc,height=input$height_plot_10_features_sc, onefile = FALSE) # open the pdf device
        plot(plot_10_features())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_sc_merged <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$project_name,"_SC_Merged_UMAP", gsub("/", "-", x), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_plot_10_features_sc, height = input$height_png_plot_10_features_sc, res = input$resolution_PNG_plot_10_features_sc)
        plot( plot_10_features())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )


    ## PCA and chosing # of dimentions to reduce ----
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
        x <- gsub(":", ".", Sys.time())
        paste(input$project_name,"_elbowPlot_sc",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_create_elbowPlot_sc,height=input$height_plot_create_elbowPlot_sc, onefile = FALSE) # open the pdf device
        plot(plot_10_features())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_create_elbowPlot_sc <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$project_name,"_elbowPlot_sc", gsub("/", "-", x), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_create_elbowPlot_sc, height = input$height_png_create_elbowPlot_sc, res = input$resolution_PNG_create_elbowPlot_sc)
        plot( plot_10_features())
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
      dataframe = read.csv(inFile_sc_meta$datapath)
    })

    output$DEx_view.meta.dt <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
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
             "Add Metadata")
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

    output$DEx_table_meta.data <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
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
        paste(input$project_name,"_SC.obj_",gsub("-", ".", Sys.Date()),".h5Seurat", sep = "")
      },
      content = function(file){
        SaveH5Seurat(vals_meta.sc$metadata_SCobj,file)
      })

    # merging multiple Seurat Obj -----
    getData <- reactive({
      inFile.seq <- input$file1_h5Seurat.file
      num <- dim(inFile.seq)[1]
      samples_list <- vector("list", length = num)
      samples_list
      for (i in 1:num) {
        sc <- LoadH5Seurat(input$file1_h5Seurat.file[[i, 'datapath']])
        samples_list[[i]] <- sc
      }
      samples_list
    })
    output$testing_mult <- renderPrint({
      sc <- input$file1_h5Seurat.file
      validate(
        need(nrow(sc)>0,
             "Upload files")
      )
      df <- getData()
      print(df)

    })
    merging_sc <- reactive({
      sc <- input$file1_h5Seurat.file
      validate(
        need(nrow(sc)>0,
             "Upload files")
      )
      samples_list <- getData()
      num <- length(samples_list)
      pbmc.normalized <- (samples_list[[1]])

      for (i in 2:num ) {
        pbmc.normalized <- merge(pbmc.normalized, y = samples_list[[i]], add.cell.ids = c("", paste("df",i,sep = "")), project = input$project_name2,merge.data = TRUE)
        pbmc.normalized@meta.data
      }
      pbmc.normalized
    })
    output$testing_mult2 <- renderPrint({
      pbmc.normalized <- merging_sc()
      table(pbmc.normalized$orig.ident)

    })

    Vals_norm <- reactiveValues(Norm1=NULL)
    Vals_norm1 <- reactiveValues(Norm1=NULL)
    Vals_norm2 <- reactiveValues(Norm1=NULL)
    Vals_norm3 <- reactiveValues(Norm1=NULL)
    Vals_norm4 <- reactiveValues(Norm1=NULL)

    observeEvent(input$run_var,{
      sc <- merging_sc()
      validate(
        need(nrow(sc)>0,
             "Run Variable")
      )
      sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
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

    observeEvent(input$run_scale,{
      sc <- Vals_norm4$Norm1
      validate(
        need(nrow(sc)>0,
             "Run Scale")
      )

      kmeans <- read.csv("../inst/Kmeans.requires.annotation.csv")

      var.genes <- as.data.frame(sc@assays$RNA@var.features)
      names(var.genes) <- "V1"
      all.genes <- rbind(kmeans,all.genes)
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
      sc <- RunPCA(sc, features = VariableFeatures(object = sc))
      Vals_norm2$Norm1 <- sc
    })

    output$PCA_harmony_verbrose <- renderPrint({
      sc <- Vals_norm2$Norm1
      validate(
        need(nrow(sc)>0,
             "Run Harmony")
      )
      sc
    })

    observeEvent(input$run_harmony,{
      sc <- Vals_norm2$Norm1
      validate(
        need(nrow(sc)>0,
             "Run Harmony")
      )
      sc <- sc %>%
        RunHarmony("orig.ident", plot_convergence = TRUE)
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

    # download Harmony merged ----
    output$downloadPlot_sc_merged <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$project_name2,"_SC_Merged_UMAP",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_sc_merged,height=input$height_sc_merged, onefile = FALSE) # open the pdf device
        plot(DimPlot(Vals_norm$Norm1, reduction = "umap", group.by = "orig.ident", pt.size = 1))
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_sc_merged <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$project_name2,"_SC_Merged_UMAP", gsub("/", "-", x), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_sc_merged, height = input$height_png_sc_merged, res = input$resolution_PNG_sc_merged)
        plot(DimPlot(Vals_norm$Norm1, reduction = "umap", group.by = "orig.ident", pt.size = 1))
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    output$downloaddf_SeruatObj_merged <- downloadHandler(
      filename = function(){
        paste(input$project_name2,"_merged_",gsub("-", ".", Sys.Date()),".h5Seurat", sep = "")
      },
      content = function(file){
        as.h5Seurat(Vals_norm$Norm1,file)
      } )


    # Add cell annotations to merged Seurat object -----
    getData_2 <- reactive({
      inFile_sc_pro2 <- input$file1_h5Seurat.file2
      if (is.null(inFile_sc_pro2)) return(NULL)
      else {
        dataframe = LoadH5Seurat(inFile_sc_pro2$datapath)
      }

    })

    output$testing_mult_anno <- renderPrint({
      sc <- input$file1_h5Seurat.file2
      validate(
        need(nrow(sc)>0,
             "Upload files")
      )
      df <- getData_2()

      print(df)

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

    test.data_anno <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Update K-means")
      )
      df= as.data.frame(sc[["RNA"]]@scale.data)

      MainTcell <- as.data.frame(t(df))
      MainTcell$Cell_Index <- rownames(MainTcell)
      MainTcell <- MainTcell %>%
        select(Cell_Index, everything())
      head(MainTcell)[1:6]
      md <- sc@meta.data
      md$Motif_gene <- md[,names(md) %in% input$V_gene_Class_2]

      md$Chain2 <- ifelse(grepl("_._", md$Motif_gene),0,
                          ifelse(md$Motif_gene=="_",0,
                                 ifelse(md$Motif_gene=="",0,
                                        ifelse(md$Motif_gene==" & ",0,
                                               ifelse(grepl("TRBV",md$Motif_gene) |
                                                        grepl("TRAV",md$Motif_gene),1,
                                                      ifelse(grepl("TRGV",md$Motif_gene) |
                                                               grepl("TRDV",md$Motif_gene),-1,0))))))


      md2 <- md[names(md) %in% c("Cell_Index","Chain2")]
      head(md2)
      MainTcell <- merge(md2,MainTcell,by="Cell_Index",all.y=T)
      rownames(MainTcell) <- MainTcell$Cell_Index
      MainTcell
    })

    output$DEx_table_TcellClass_3 <-  DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      sc <- test.data_anno()
      head(sc)[1:6]
    })
    # add T cell classification
    add.classification_T_cell_Df_2 <- reactive({
      MainTcell <- test.data_anno()

      if (input$add.classification_T_cell_2==T) {

        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("PTPRC","CD4","CD8A","CD8B","CD3D","CD3E","CD3G","MS4A1","JCHAIN","CLEC10A","FCGR3A","S100A9","TRAC","TRDC","TRGV9", "TRDV2","TRDV1")]
        head(MainTcell_test)
        df <- MainTcell_test
        nms <- c("CD4","CD8A","CD8B","CD3D","CD3E","CD3G","MS4A1","JCHAIN","CLEC10A","FCGR3A","S100A9","TRDC","TRAC","TRGV9", "TRDV2","TRDV1","PTPRC")   # Vector of columns you want in this data
        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]

        centre_size <- round(dim(MainTcell_test)[1]/20,0)

        if(centre_size>500) {
          centre_size=500
        }

        else {

          centre_size <- round(dim(MainTcell_test)[1]/20,0)
        }

        set.seed(123)
        kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 1)

        df_centres <- as.data.frame(kmeans10$centers)




        df_centres <- df_centres %>%
          mutate(classify.Adaptive.cell = case_when(
            TRAC>0 & TRDC <0 & c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 >0 & CD8B > 0 & CD8A > 0 ~ "CD4+CD8ab+ ab T cells",
            TRAC<0 & TRDC >0 & c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 >0 & CD8B > 0 & CD8A > 0 ~ "CD4+CD8ab+ gd T cells",
            c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 >0 & CD8B < 0 & CD8A > 0 ~ "CD4+CD8a+",
            c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 >0 & CD8B < 0 & CD8A < 0 & TRAC> 0 ~ "CD4+ ab T cell",
            c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B > 0 & CD8A > 0 & TRDV1>0 ~ "CD8ab+ Vd1 T cells",
            c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B > 0 & CD8A > 0 ~ "CD8ab+",
            c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B < 0 & CD8A > 0 & TRDC>0 & TRDV2>0 & TRGV9>0 ~ "CD8a+ Vd2g9 T cell",
            c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B < 0 & CD8A > 0 & TRDC>0 & TRDV1>0 ~ "CD8a+ Vd1 T cell",
            c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B < 0 & CD8A > 0 & TRDC>0 & TRDV2>0  ~ "CD8a+ Vd2 T cell",
            c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B < 0 & CD8A > 0 & TRDC>0 ~ "CD8a+ T cell",
            c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B < 0 & CD8A > 0 ~ "CD8a+",
            c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B < 0 & CD8A < 0 & TRDC>0 & TRAC<0 ~ "CD4-CD8- gd T cell",
            c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B < 0 & CD8A < 0 & TRAC>0 & TRDC<0 ~ "CD4-CD8- ab T cell",
            c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B < 0 & CD8A < 0   ~ "CD4-CD8- T cell",
            # CD3D < 0 & CD3E<0 & CD3G<0 & CD4 <0 & CD8B < 0 & CD8A < 0~ "Not a T cell",
            CD8A>0 & TRDC>0 ~ "CD8a+ gd T cell",
            TRDC>0 ~ "gd T cell",
            c(CD3D > 0 | CD3E>0 | CD3G>0) & TRAC>0 ~ "ab T cell",
            CD3D > 0 | CD3E>0 | CD3G>0 ~ "T cell",
            MS4A1>0 ~ "B cell",
            JCHAIN>0 ~ "Plasma cell",
            CLEC10A>0 ~ "DC",
            FCGR3A>0 ~ "NK-like",
            S100A9>0 ~ "Monocytes",
            CD4 >0 & CD8B < 0 & CD8A < 0 & TRAC> 0 ~ "CD4+ ab T cell",
            CD4 >0 & CD8B < 0 & CD8A < 0 ~ "CD4+ (TCR downreg)",
            CD4 <0 & CD8B > 0 & CD8A > 0 ~ "CD8ab+ (TCR downreg)",
            TRAC>0  ~ "ab T cell",
            PTPRC>0 ~ "T/B cell",
            TRUE ~ NA_character_))

        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","classify.Adaptive.cell")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","classify.Adaptive.cell")]
        Cluster_Df_class2 <- Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
        Cluster_Df_class2

      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2
      }

    })
    add.classification_T.B_cell_Function_Df_2 <- reactive({
      MainTcell <- test.data_anno()

      if (input$add.classification_T_cell_Function_2==T) {

        rownames(MainTcell) <- MainTcell$Cell_Index
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                                                             "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM","MS4A1","IGHA2","IGHA1","JCHAIN","IGHG1","IGHG2","IGHG3","IGHG4","IGHE","IGHD")]



        df <- MainTcell_test
        # Vector of columns you want in this data
        nms <- c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                 "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM","MS4A1","IGHA2","IGHA1","JCHAIN","IGHG1","IGHG2","IGHG3","IGHG4","IGHE","IGHD")

        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]

        d_frame <- na.omit(MainTcell_test)
        centre_size <- round(dim(MainTcell_test)[1]/10,0)
        if(centre_size>500) {
          centre_size=500
        }
        else {
          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }

        set.seed(123)
        kmeans10 <- kmeans(d_frame, centers = centre_size, nstart = 1)
        df_centres <- as.data.frame(kmeans10$centers)

        df_centres <- df_centres %>%
          mutate(classify.T.cell_Function = case_when(
            CD4 > 0 & c(CD8A<0 | CD8B<0) & FOXP3 >0 & IL2RA >0 ~ "CD4+ Treg FOXP3+CD25+",
            CD4 < 0 & c(CD8A>0 | CD8B>0) & TGFB1>0 & FOXP3 >0 ~ "CD8+ Tregs FOXP3+TGFB1+",
            CD4 >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & IL2>0 & IFNG>0 ~ "IL2+IFNg+ (Th1)",
            CD4 >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & TNF >0 & IFNG>0 ~ "IFNg+TNF+ (Th1)",
            CD4 >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & TNF >0 & IFNG<0 ~ "IFNg-TNF+ (Th1)",
            CD4 >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & TNF <0 & IFNG>0 ~ "IFNg+TNF- (Th1)",
            CD4 >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & IL2>0 ~ "Th1 IL2+",
            CD4 >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & IFNG>0 ~ "Th1 IFNg+",
            # CD4 >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) ~ "Th1 (other)",
            CD4 >0 & IL4>0 ~ "CD4+IL4+ (Th2)",
            CD4 >0 & c(IL9>0) ~ "CD4+IL9+ (Th9)",
            CD4 >0 & c(IL17A>0 | IL17F>0) ~ "CD4+IL17+ (Th17)",
            CD4 >0 & c(IL22>0) ~ "IL22+ (Th22)",
            CD4 >0 & c(IL21>0) ~ "IL21+ (Tfh)",
            CD8A> 0 & CCR6>0 & KLRB1>0 ~ "CD8a+CCR6+ KLRB1+ (MAIT)",
            CD4<0 & c(CD8A>0 | CD8B>0) & GZMB < 0 & GNLY>0 & PRF1<0 ~ "GNLY+CD8+ Cytotoxic",
            CD4<0 & c(CD8A>0 | CD8B>0) & GZMB > 0 & GNLY<0 & PRF1>0 ~ "GZMB+PRF1+CD8+ Cytotoxic",
            CD4<0 & c(CD8A>0 | CD8B>0) & c(GZMB > 0 | GNLY>0 | PRF1 >0) ~ "CD8+ Cytotoxic",
            CD4>0 & CD8A<0 & CD8B<0 & c(GZMB > 0 | GNLY>0 | PRF1 >0) ~ "CD4+ Cytotoxic",
            CD4<0 & CD8A <0 & CD8B <0 & c(GZMB > 0 | GNLY>0 | PRF1 >0) ~ "CD4-CD8- Cytotoxic",
            CD4>0 & c(CD8A >0 | CD8B >0) & c(GZMB > 0 | GNLY>0 | PRF1 >0) ~ "CD4+CD8+ Cytotoxic",
            CD4>0 & TNF>0 ~ "TNF+CD4+",
            CD4>0 & IFNG>0 ~ "IFNg+CD4+",
            # B cell markers
            IGHM>0 & MS4A1>0 ~ "IgM+ B cell",
            IGHM<0 & MS4A1>0~ "IgM- B cell",
            c(IGHA2>0 | IGHA1>0) & JCHAIN>0 ~ "IgA+ Plasma cell",
            c(IGHG1>0 | IGHG2>0| IGHG3>0| IGHG4>0) & JCHAIN>0 ~ "IgG+ Plasma cell",
            IGHE>0 & JCHAIN>0 ~ "IgE+ Plasma cell",
            IGHD>0 & MS4A1>0 ~ "IgD+ B cell",
            TRUE ~ NA_character_))
        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","classify.T.cell_Function")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        Cluster_Df
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","classify.T.cell_Function")]
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]

      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }
    })
    add.classification_T_cell_Function_CD4_df_2 <- reactive({
      MainTcell <- test.data_anno()

      if (input$add.classification_T_cell_Function_CD4_2==T) {
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                                                             "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM","MS4A1","IGHA2","IGHA1","JCHAIN","IGHG1","IGHG2","IGHG3","IGHG4","IGHE","IGHD")]
        df <- MainTcell_test
        # Vector of columns you want in this data
        nms <- c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                 "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM","MS4A1","IGHA2","IGHA1","JCHAIN","IGHG1","IGHG2","IGHG3","IGHG4","IGHE","IGHD")
        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]
        MainTcell_test[is.na(MainTcell_test)] <- 0
        MainTcell_test
        centre_size <- round(dim(MainTcell_test)[1]/10,0)

        if(centre_size>500) {
          centre_size=500
        }

        else {

          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }

        set.seed(123)
        kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 1)
        kmeans10$centers
        df_centres <- as.data.frame(kmeans10$centers)
        df_centres <- df_centres %>%
          mutate(classify.T.cell_Function_CD4 = case_when(
            CD4 > 0 & c(CD8A<0 | CD8B<0) & FOXP3 >0 & IL2RA >0 ~ "CD4+ Treg FOXP3+CD25+",
            CD4 >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & IL2>0 & IFNG>0 ~ "IL2+IFNg+ (Th1)",
            CD4 >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & TNF >0 & IFNG>0 ~ "IFNg+TNF+ (Th1)",
            CD4 >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & TNF >0 & IFNG<0 ~ "IFNg-TNF+ (Th1)",
            CD4 >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & TNF <0 & IFNG>0 ~ "IFNg+TNF- (Th1)",
            CD4 >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & IL2>0 ~ "Th1 IL2+",
            CD4 >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & IFNG>0 ~ "Th1 IFNg+",
            CD4 >0 & IL4>0 ~ "CD4+IL4+ (Th2)",
            CD4 >0 & c(IL9>0) ~ "CD4+IL9+ (Th9)",
            CD4 >0 & c(IL17A>0 | IL17F>0) ~ "CD4+IL17+ (Th17)",
            CD4 >0 & c(IL22>0) ~ "CD4+IL22+ (Th22)",
            CD4 >0 & c(IL21>0) ~ "CD4+IL21+ (Tfh)",
            CD8A> 0 & CCR6>0 & KLRB1>0 ~ "CD8a+CCR6+ KLRB1+ (MAIT)",
            CD4>0 & GZMB < 0 & GNLY>0 & PRF1<0 ~ "GNLY+CD4+ Cytotoxic",
            CD4>0 & GZMB > 0 & GNLY<0 & PRF1>0 ~ "GZMB+PRF1+CD4+ Cytotoxic",
            # CD4>0 & c(CD8A >0 | CD8B >0) & c(GZMB > 0 | GNLY>0 | PRF1 >0) ~ "CD4+CD8+ Cytotoxic",
            CD4>0 & c(GZMB > 0 | GNLY>0 | PRF1 >0) ~ "CD4+ Cytotoxic",
            # CD4>0 & MKI67 >0 & TOP2A > 0  ~ "Cycling CD4+ T cells",
            CD4>0 & TNF>0 ~ "TNF+CD4+",
            CD4>0 & IFNG>0 ~ "IFNg+CD4+",
            CD4 >0 & CD8A<0 & CD8B<0 ~ "CD4+",
            TRUE ~ NA_character_))

        df_centres$clust_name <- rownames(df_centres)
        (df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","classify.T.cell_Function_CD4")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","classify.T.cell_Function_CD4")]
        Cluster_Df_class2 <- Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
        Cluster_Df_class2
      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }


    })
    add.classification_T_cell_Function_CD8_df_2 <- reactive({

      MainTcell <- test.data_anno()

      if (input$add.classification_T_cell_Function_CD8_2==T) {
        rownames(MainTcell) <- MainTcell$Cell_Index
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                                                             "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM","MS4A1","IGHA2","IGHA1","JCHAIN","IGHG1","IGHG2","IGHG3","IGHG4","IGHE","IGHD")]
        df <- MainTcell_test
        # Vector of columns you want in this data
        nms <- c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                 "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM","MS4A1","IGHA2","IGHA1","JCHAIN","IGHG1","IGHG2","IGHG3","IGHG4","IGHE","IGHD")
        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]

        MainTcell_test[is.na(MainTcell_test)] <- 0
        centre_size <- round(dim(MainTcell_test)[1]/10,0)

        if(centre_size>500) {
          centre_size=500
        }

        else {

          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }
        set.seed(123)
        kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 1)
        df_centres <- as.data.frame(kmeans10$centers)

        df_centres <- df_centres %>%
          mutate(classify.T.cell_Function_CD8 = case_when(
            CD4 < 0 & c(CD8A>0 | CD8B>0) & TGFB1>0 & FOXP3 >0 ~ "CD8+ Tregs FOXP3+TGFB1+",
            CD8A> 0 & CCR6>0 & KLRB1>0 ~ "CD8a+CCR6+KLRB1+ (MAIT)",
            CD4<0 & c(CD8A>0 | CD8B>0) & GZMB < 0 & GNLY>0 & PRF1<0 ~ "GNLY+CD8+ Cytotoxic",
            CD4<0 & c(CD8A>0 | CD8B>0) & GZMB > 0 & GNLY<0 & PRF1>0 ~ "GZMB+PRF1+CD8+ Cytotoxic",
            CD4<0 & c(CD8A>0 | CD8B>0) & c(GZMB > 0 | GNLY>0 | PRF1 >0) ~ "CD8+ Cytotoxic",
            CD4<0 &c(CD8A >0 | CD8B >0) & IL6>0 ~ "CD8+IL6+",
            CD4<0 &c(CD8A >0 | CD8B >0) & IFNG>0 ~ "CD8+IFNg+",
            CD4<0 &c(CD8A >0 | CD8B >0) & TNF>0 ~ "CD8+TNF+",
            CD4<0 & c(CD8A >0 | CD8B >0) & IL2>0 ~ "CD8+IL2+",
            CD4<0 &c(CD8A >0 | CD8B >0)  ~ "CD8+ T cells",
            TRUE ~ NA_character_))

        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","classify.T.cell_Function_CD8")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        Cluster_Df
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","classify.T.cell_Function_CD8")]
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }


    })
    add.classification_T_cell_Function_CD4_CD8pos_df_2 <- reactive({
      MainTcell <- test.data_anno()

      if (input$add.classification_T_cell_Function_CD4_CD8pos_2==T) {
        rownames(MainTcell) <- MainTcell$Cell_Index
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                                                             "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM","MS4A1","IGHA2","IGHA1","JCHAIN","IGHG1","IGHG2","IGHG3","IGHG4","IGHE","IGHD")]
        df <- MainTcell_test
        # Vector of columns you want in this data
        nms <- c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                 "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM","MS4A1","IGHA2","IGHA1","JCHAIN","IGHG1","IGHG2","IGHG3","IGHG4","IGHE","IGHD")
        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]

        MainTcell_test[is.na(MainTcell_test)] <- 0
        centre_size <- round(dim(MainTcell_test)[1]/10,0)

        if(centre_size>500) {
          centre_size=500
        }

        else {

          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }
        set.seed(123)
        kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 1)
        df_centres <- as.data.frame(kmeans10$centers)

        df_centres <- df_centres %>%
          mutate(classify.T.cell_Function_CD4_CD8pos = case_when(
            CD4>0 & c(CD8A >0 | CD8B >0) & c(GZMB > 0 | GNLY>0 | PRF1 >0) ~ "CD4+CD8+ Cytotoxic",
            CD4>0 & c(CD8A >0 | CD8B >0) & TNF>0 ~ "TNF+CD4+CD8+",
            CD4>0 & c(CD8A >0 | CD8B >0) & IFNG>0 ~ "IFNg+CD4+CD8+",
            CD4>0 & c(CD8A >0 | CD8B >0) ~ "CD4+CD8+",
            TRUE ~ NA_character_))

        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","classify.T.cell_Function_CD4_CD8pos")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        Cluster_Df
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","classify.T.cell_Function_CD4_CD8pos")]
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }


    })
    add.classification_T_cell_Function_CD4_CD8neg_df_2 <- reactive({
      MainTcell <- test.data_anno()
      if (input$add.classification_T_cell_Function_CD4_CD8neg_2==T) {
        rownames(MainTcell) <- MainTcell$Cell_Index
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F",
                                                             "IL22","IL21","CCR6","GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","MKI67","TOP2A")]
        df <- MainTcell_test
        # Vector of columns you want in this data
        nms <- c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6","GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","MKI67","TOP2A")
        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]

        MainTcell_test[is.na(MainTcell_test)] <- 0
        centre_size <- round(dim(MainTcell_test)[1]/10,0)

        if(centre_size>500) {
          centre_size=500
        }

        else {

          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }
        set.seed(123)
        kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 1)
        df_centres <- as.data.frame(kmeans10$centers)

        df_centres <- df_centres %>%
          mutate(classify.T.cell_Function_CD4_CD8neg = case_when(
            c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4<0 & CD8A <0 & CD8B <0 & c(GZMB > 0 | GNLY>0 | PRF1 >0) ~ "CD4-CD8- Cytotoxic",
            c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B < 0 & CD8A < 0 ~ "CD4-CD8-",
            TRUE ~ NA_character_))

        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","classify.T.cell_Function_CD4_CD8neg")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        Cluster_Df
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","classify.T.cell_Function_CD4_CD8neg")]

        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }


    })
    add.classification_B_cell_Function_Df_2 <- reactive({

      MainTcell <- test.data_anno()

      if (input$add.classification_B_cell_Function_2==T) {

        rownames(MainTcell) <- MainTcell$Cell_Index
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                                                             "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM","MS4A1","IGHA2","IGHA1","JCHAIN","IGHG1","IGHG2","IGHG3","IGHG4","IGHE","IGHD")]
        df <- MainTcell_test
        # Vector of columns you want in this data
        nms <- c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                 "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM","MS4A1","IGHA2","IGHA1","JCHAIN","IGHG1","IGHG2","IGHG3","IGHG4","IGHE","IGHD")

        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]

        d_frame <- na.omit(MainTcell_test)
        centre_size <- round(dim(MainTcell_test)[1]/10,0)

        if(centre_size>500) {
          centre_size=500
        }

        else {

          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }
        set.seed(123)
        kmeans10 <- kmeans(d_frame, centers = centre_size, nstart = 1)
        df_centres <- as.data.frame(kmeans10$centers)

        df_centres <- df_centres %>%
          mutate(classify.B.cell_Function = case_when(

            # B cell markers
            IGHM>0 & MS4A1>0 ~ "IgM+ B cell",
            IGHM<0 & MS4A1>0~ "IgM- B cell",
            c(IGHA2>0 | IGHA1>0) & JCHAIN>0 ~ "IgA+ Plasma cell",
            c(IGHG1>0 | IGHG2>0| IGHG3>0| IGHG4>0) & JCHAIN>0 ~ "IgG+ Plasma cell",
            IGHE>0 & JCHAIN>0 ~ "IgE+ Plasma cell",
            IGHD>0 & MS4A1>0 ~ "IgD+ B cell",

            TRUE ~ NA_character_))


        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","classify.B.cell_Function")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        Cluster_Df
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","classify.B.cell_Function")]
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]

      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }


    })
    add.classification_T_cell_Memory_df_2 <- reactive({
      MainTcell <- test.data_anno()
      if (input$add.classification_T_cell_Memory_2==T) {
        rownames(MainTcell) <- MainTcell$Cell_Index
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD8A","CD8B","CD4","CD27","CCR7","SELL", "CXCR3","IL7R","IL2RA","CX3CR1","CD2","CD5","CD7","PDCD1","CD69","KLRG1",
                                                             "CD44","ITGAE","CCR8","KLRC1","CTLA4","ICOS","JCHAIN","MS4A1")]

        df <- MainTcell_test
        nms <- c("CD8A","CD8B","CD4","CD27","CCR7","SELL", "CXCR3","IL7R","IL2RA","CX3CR1","CD2","CD5","CD7","PDCD1","CD69","KLRG1",
                 "CD44","ITGAE","CCR8","KLRC1","CTLA4","ICOS","JCHAIN","MS4A1")  # Vector of columns you want in this data
        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]
        head(MainTcell_test)
        MainTcell_test[is.na(MainTcell_test)] <- 0
        centre_size <- round(dim(MainTcell_test)[1]/10,0)
        if(centre_size>500) {
          centre_size=500
        }
        else {
          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }
        set.seed(123)
        kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 10)
        kmeans10$centers
        df_centres <- as.data.frame(kmeans10$centers)

        df_centres <- df_centres %>%
          mutate(classify.T.cell_Memory = case_when(
            CD27>0 & SELL>0 & IL7R>0 & CCR7>0 ~ "Naive (CD27+CD62L+CD127+CCR7+)",
            IL2RA>0 & IL7R<0 & CCR7<0  ~"Effector CD25+CD127-CCR7-", #CD25+CD127-
            IL2RA>0 & CD69>0 & KLRG1>0 ~ "Effector CD25+CD69+KLRG1+",
            IL2RA<0 & IL7R>0 & CCR7<0~ "Effector Memory CD25-CD127+CCR7-", #CD25 IL2RA)-CD127 (IL7R)+
            IL2RA>0 & IL7R>0  ~ "Central Memory CD25+CD127+", #CD25 IL2RA)-CD127 (IL7R)+
            SELL>0 & CD44>0 & CCR7>0 ~ "Central memory CD62L+CD44+CCR7+",
            PDCD1>0 & IL2RA>0 & c(CD8A>0 | CD8B>0)~ "CD25+PD−1+CD8+",
            IL2RA>0 & KLRG1 >0 & IL7R>0 & c(CD8A>0 | CD8B>0) & CD4<0 ~ "Other CD8+ memory",
            CX3CR1>0 ~ "Peripheral memory",
            CD69>0 & ITGAE>0 & CCR8>0 ~ "Skin-Tissue resident", #
            CD69>0 & ITGAE>0 & KLRC1>0 & CTLA4>0 & ICOS>0 ~ "Lung-Tissue resident", #
            CD69>0 & ITGAE>0  ~ "Other resident", #
            CD2>0 & CD5>0 & CD7>0 ~ 'Developing T cells (CD2+CD5+CD7+)',
            CD27>0 & JCHAIN>0 ~ "Memory Plasma cell",
            CD27>0 & MS4A1>0 ~ "Memory B cell",
            TRUE ~ NA_character_))

        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","classify.T.cell_Memory")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","classify.T.cell_Memory")]
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }


    })
    add.classification_T_cell_Activation_df_2 <- reactive({
      MainTcell <- test.data_anno()
      if (input$add.classification_T_cell_Activation_2==T) {
        rownames(MainTcell) <- MainTcell$Cell_Index
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD8A","CD8B","CD4","IL2RA","CD27","TNFRSF8","DPP4","CD40LG","TNFRSF4","ICAM1","PDCD1","B3GAT1","IL6","JCHAIN","MS4A1")]

        df <- MainTcell_test
        nms <- c("CD8A","CD8B","CD4","IL2RA","CD27","TNFRSF8","DPP4","CD40LG","TNFRSF4","ICAM1","PDCD1","B3GAT1","IL6","JCHAIN","MS4A1")   # Vector of columns you want in this data
        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]
        head(MainTcell_test)
        MainTcell_test[is.na(MainTcell_test)] <- 0
        centre_size <- round(dim(MainTcell_test)[1]/10,0)

        if(centre_size>500) {
          centre_size=500
        }

        else {

          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }
        set.seed(123)
        kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 10)
        kmeans10$centers
        df_centres <- as.data.frame(kmeans10$centers)

        df_centres <- df_centres %>%
          mutate(classify.T.B.cell_Activation = case_when(
            # https://pubmed.ncbi.nlm.nih.gov/22120733/
            # CD69>0 ~ "CD69+",
            IL2RA>0 & CD27>0 ~ "CD25+CD27+",
            IL2RA>0 & TNFRSF8 >0 ~ "CD25+TNFRSF8+",
            IL2RA>0 & DPP4 >0 ~ "CD25+DPP4+",
            IL2RA>0 & CD40LG >0 ~ "CD25+CD40LG+",
            IL2RA>0 & TNFRSF4 >0 ~ "CD25+TNFRSF4+",
            IL2RA>0 & ICAM1> 0 ~  "CD25+CD54+",
            PDCD1>0 ~ "Exhaused (PD1+)", # PDCD1= PD1 AND B3GAT1=CD57
            B3GAT1>0 ~ "Senescence (CD57+)", # PDCD1= PD1 AND B3GAT1=CD57
            JCHAIN>0 & IL6>0 ~ "IL6+ Plasma cell",
            MS4A1>0 & IL6>0 ~ "IL6+ B cell",
            TRUE ~ NA_character_))

        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","classify.T.B.cell_Activation")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","classify.T.B.cell_Activation")]
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }


    })
    add.classification_CellTypist_list_lower_2 <- reactive({
      MainTcell <- test.data_anno()
      if (input$add.classification_CellTypist_list_overview_2==T) {
        rownames(MainTcell) <- MainTcell$Cell_Index
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("MKI67","TOP2A","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD","MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","S100A9","EBI3","CCL19","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","FCGR3A","S100A13","TLE1","AREG","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","C1QA","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1","Chain2")]
        # -----
        df <- MainTcell_test
        # Vector of columns you want in this data
        nms <- c("MKI67","TOP2A","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD","MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","S100A9","EBI3","CCL19","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","FCGR3A","S100A13","TLE1","AREG","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","C1QA","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1","Chain2")

        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]
        MainTcell_test[is.na(MainTcell_test)] <- 0
        centre_size <- round(dim(MainTcell_test)[1]/10,0)

        if(centre_size>500) {
          centre_size=500
        }

        else {

          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }
        set.seed(123)
        kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 1)
        df_centres <- as.data.frame(kmeans10$centers)

        df_centres <- df_centres %>%
          mutate(CellTypist_list_lower = case_when(
            #### B cells -----
            FCRL2>0 & ITGAX>0 & TBX21>0     ~ "Age-associated B cells", # B cells
            CD79A>0 & MS4A1>0 & CD19>0      ~ "B cell", # B cell
            CXCR5>0 & TNFRSF13B>0 & CD22>0  ~ "Follicular B cells",
            POU2AF1>0 & CD40>0 & SUGCT>0    ~ "Germinal center B cells",
            CR2>0 & CD27>0 & MS4A1>0        ~ "Memory B cells",
            IGHM>0 & IGHD>0 & TCL1A>0       ~ "Naive B cells",
            MKI67>0 & SUGCT>0 & AICDA>0     ~ "Proliferative germinal center B cells",
            CD24>0 & MYO1C>0 & MS4A1>0      ~ "Transitional B cells",
            MME>0 & CD24>0 & MKI67>0        ~ "Large pre-B cells",
            IL7R>0 & ZCCHC7>0 & RAG1>0      ~ "Pre-pro-B cells",
            MME>0 &  DNTT>0 & GLL1>0        ~ "Pro-B cells",
            MME>0 & CD24>0 & IGLL5>0        ~ "Small pre-B cells",
            JCHAIN>0 & MZB1>0 & XBP1>0      ~ "Plasma cell",
            # MKI67 >0 & TOP2A > 0 & Chain2==1 ~ "Cycling ab T cells",
            # MKI67 >0 & TOP2A > 0 & Chain2== -1 ~ "Cycling gd T cells",
            #### T cells ----
            ACY3>0 & CD34>0 & SPINK2>0 ~ "ETP",
            FXYD2>0 & HES1>0 & CD99>0 ~ "DN thymocytes",
            CD1A>0 & CD8A>0 & SMPD3>0 ~ "DP thymocytes",
            ZNF683 >0 & GNG4 >0 & PDCD1 >0 ~ "CD8aa T cells",
            TOX2 >0 & SATB1 >0 & CCR9 >0 ~ "CD8a/b(entry)",
            PDCD1 >0 & ICOS >0 & CXCR5>0  ~"Follicular helper T cells",
            KLRB1 >0 & SLC4A10 >0 & `TRAV1-2`>0  ~ "KLRB1+ MAIT cells (TRAV1-2+)",
            GZMK >0 & CD4>0 & IL10>0 ~ "Memory CD4+ cytotoxic T cells",
            NKG7>0 & GNLY>0 & CD8A>0 ~ "NK CD8+ T cells",
            CTLA4>0 & IL2RA >0 & FOXP3>0 ~ "Regulatory T cells",
            MIR155HG>0 & BIRC3>0 & SMS>0 ~ "T(agonist)",
            CD8A>0 & CCR7>0 & SELL>0 ~ "Tcm/Naive cytotoxic T cells",
            CD4>0 & CCR7>0 & SELL>0 ~ "Tcm/Naive helper T cells",
            KLRB1>0 & AQP3>0 & ITGB1>0 ~ "Tem/Effector helper T cells",
            PDCD1>0 & CD4>0 & CTLA4>0 ~ "Tem/Effector helper T cells PD1+",
            CX3CR1>0 & GZMB>0 & GNLY>0 ~ "Tem/Temra cytotoxic T cells",
            GZMK>0 & CD8A>0 & CCL5>0 ~ "Tem/Trm cytotoxic CD8+ T cells",
            CD27>0 & CCR7>0 & IKZF2>0 ~ "Treg(diff)",
            ITGA1>0 & ITGAE>0 & CXCR6>0 ~"Trm cytotoxic T cells",
            CCL5>0 & CXCR3>0 & TBX21>0  ~ "Type 1 helper T cells",
            IL7R>0 & CCR6>0 & ZBTB16>0 ~ "Type 17 helper T cells",
            TRDC>0 & TRGC1>0 & CCL5>0 ~ "Gamma-delta T cells",

            #### DC -----
            CD1C>0 & FCER1A>0 & CLEC10A>0 ~ "DC",
            BATF3>0 & CADM1>0 & CLEC9A>0 ~ "DC1",
            CLEC10A>0 & FCER1A>0 & CD1C>0 ~ "DC2",
            CLEC10A>0 & FCER1A>0 & S100A9>0 ~ "DC3",
            EBI3>0 & CCR7>0 & CCL19>0 ~ "Migratory DC",
            CLEC10A>0 & KLF4>0 & AXL>0 ~ " Transitional DC",
            IRF8>0 & CLEC10A>0 & PCLAF>0 ~ "DC precursor",
            #### ILC cells ----
            GNLY>0 & FCGR3A>0 & NKG7>0 ~ "CD16+ NK cells",
            GNLY>0 & CD160>0 & NKG7>0 ~ "CD16- NK cells",
            S100A13>0 & TLE1>0 & AREG>0 ~ "ILC",
            CXCR3>0 & CD3D>0 & IKZF3>0 ~ "ILC1",
            GATA3>0 & KLRG1>0 & HPGDS>0 ~ "ILC2",
            IL4I1>0 & RORC>0 & KIT>0 ~ "ILC3",
            GNLY>0 & XCL2>0 & NKG7>0 ~ "NK cells",
            CCL5>0 & GZMK>0 & FCGR3A>0 ~ "Traditional NK cells",
            #### Monocytes -----
            TYROBP>0 & C1QC>0 & HMOX1>0 ~ "Mono-mac",
            LYZ>0 & VCAN>0 & S100A9>0 ~ "Monocyte precursor",
            S100A9>0 & CD14>0 & S100A12>0 ~ "Classical monocytes",
            S100A9>0 & LYZ>0 & FCN1>0 ~ "Monocytes",
            FCGR3A>0 & C1QA>0 & CX3CR1>0 ~ "Non-classical monocytes",


            TRUE ~ NA_character_))

        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","CellTypist_list_lower")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        Cluster_Df
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","CellTypist_list_lower")]
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }


    })
    add.classification_CellTypist_list_higher_2 <- reactive({
      MainTcell <- test.data_anno()
      if (input$add.classification_CellTypist_list_overview_2==T) {
        rownames(MainTcell) <- MainTcell$Cell_Index
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("MKI67","TOP2A","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD","MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","S100A9","EBI3","CCL19","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","FCGR3A","S100A13","TLE1","AREG","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","C1QA","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1","Chain2")]
        # -----
        df <- MainTcell_test
        # Vector of columns you want in this data
        nms <- c("MKI67","TOP2A","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD","MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","S100A9","EBI3","CCL19","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","FCGR3A","S100A13","TLE1","AREG","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","C1QA","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1","Chain2")
        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]

        MainTcell_test[is.na(MainTcell_test)] <- 0
        centre_size <- round(dim(MainTcell_test)[1]/10,0)
        if(centre_size>500) {
          centre_size=500
        }
        else {
          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }
        set.seed(123)
        kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 1)
        df_centres <- as.data.frame(kmeans10$centers)

        df_centres <- df_centres %>%
          mutate(CellTypist_list_higher= case_when(
            #### B cells -----
            FCRL2>0 & ITGAX>0 & TBX21>0     ~ "B cell", # B cell
            CD79A>0 & MS4A1>0 & CD19>0      ~ "B cell", # B cell
            CXCR5>0 & TNFRSF13B>0 & CD22>0  ~ "B cell", # B cell
            POU2AF1>0 & CD40>0 & SUGCT>0    ~ "B cell", # B cell
            CR2>0 & CD27>0 & MS4A1>0        ~ "B cell", # B cell
            IGHM>0 & IGHD>0 & TCL1A>0       ~ "B cell", # B cell
            MKI67>0 & SUGCT>0 & AICDA>0     ~ "B cell", # B cell
            CD24>0 & MYO1C>0 & MS4A1>0      ~ "B cell", # B cell
            MME>0 & CD24>0 & MKI67>0        ~ "B cell", # B cell
            IL7R>0 & ZCCHC7>0 & RAG1>0      ~ "B cell", # B cell
            MME>0 &  DNTT>0 & GLL1>0        ~ "B cell", # B cell
            MME>0 & CD24>0 & IGLL5>0        ~ "B cell", # B cell
            JCHAIN>0 & MZB1>0 & XBP1>0      ~ "Plasma cell",
            # MKI67 >0 & TOP2A > 0 & Chain2==1 ~ "Cycling ab T cells",
            # MKI67 >0 & TOP2A > 0 & Chain2== -1 ~ "Cycling gd T cells",
            #### T cells ----
            ACY3>0 & CD34>0 & SPINK2>0 ~ "T cells",
            FXYD2>0 & HES1>0 & CD99>0 ~ "T cells",
            CD1A>0 & CD8A>0 & SMPD3>0 ~ "T cells",
            ZNF683 >0 & GNG4 >0 & PDCD1 >0 ~ "T cells",
            TOX2 >0 & SATB1 >0 & CCR9 >0 ~ "T cells",
            PDCD1 >0 & ICOS >0 & CXCR5>0  ~"T cells",
            KLRB1 >0 & SLC4A10 >0 & `TRAV1-2`>0  ~ "T cells",
            GZMK >0 & CD4>0 & IL10>0 ~ "T cells",
            NKG7>0 & GNLY>0 & CD8A>0 ~ "T cells",
            CTLA4>0 & IL2RA >0 & FOXP3>0 ~ "T cells",
            MIR155HG>0 & BIRC3>0 & SMS>0 ~ "T cells",
            CD8A>0 & CCR7>0 & SELL>0 ~ "T cells",
            CD4>0 & CCR7>0 & SELL>0 ~ "T cells",
            KLRB1>0 & AQP3>0 & ITGB1>0 ~ "T cells",
            PDCD1>0 & CD4>0 & CTLA4>0 ~ "T cells",
            CX3CR1>0 & GZMB>0 & GNLY>0 ~ "T cells",
            GZMK>0 & CD8A>0 & CCL5>0 ~ "T cells",
            CD27>0 & CCR7>0 & IKZF2>0 ~ "T cells",
            ITGA1>0 & ITGAE>0 & CXCR6>0 ~"T cells",
            CCL5>0 & CXCR3>0 & TBX21>0  ~ "T cells",
            IL7R>0 & CCR6>0 & ZBTB16>0 ~ "T cells",
            TRDC>0 & TRGC1>0 & CCL5>0 ~ "T cells",
            #### DC -----
            CD1C>0 & FCER1A>0 & CLEC10A>0 ~ "DC",
            BATF3>0 & CADM1>0 & CLEC9A>0 ~ "DC",
            CLEC10A>0 & FCER1A>0 & CD1C>0 ~ "DC",
            CLEC10A>0 & FCER1A>0 & S100A9>0 ~ "DC",
            EBI3>0 & CCR7>0 & CCL19>0 ~"DC",
            CLEC10A>0 & KLF4>0 & AXL>0 ~ "DC",
            IRF8>0 & CLEC10A>0 & PCLAF>0 ~ "DC",
            #### ILC cells ----
            GNLY>0 & FCGR3A>0 & NKG7>0 ~ "ILC",
            GNLY>0 & CD160>0 & NKG7>0 ~ "ILC",
            S100A13>0 & TLE1>0 & AREG>0 ~ "ILC",
            CXCR3>0 & CD3D>0 & IKZF3>0 ~ "ILC",
            GATA3>0 & KLRG1>0 & HPGDS>0 ~"ILC",
            IL4I1>0 & RORC>0 & KIT>0 ~ "ILC",
            GNLY>0 & XCL2>0 & NKG7>0 ~ "ILC",
            CCL5>0 & GZMK>0 & FCGR3A>0 ~ "ILC",
            #### Monocytes -----
            TYROBP>0 & C1QC>0 & HMOX1>0 ~ "Monocytes",
            LYZ>0 & VCAN>0 & S100A9>0 ~ "Monocytes",
            S100A9>0 & CD14>0 & S100A12>0 ~ "Monocytes",
            S100A9>0 & LYZ>0 & FCN1>0 ~ "Monocytes",
            FCGR3A>0 & C1QA>0 & CX3CR1>0 ~ "Monocytes",


            TRUE ~ NA_character_))

        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","CellTypist_list_higher")]
        Cluster_Df <- as.data.frame(kmeans10$cluster)
        Cluster_Df
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","CellTypist_list_higher")]
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }
    })
    add.classification_CellTypist_list_lower_Tcell <- reactive({
      MainTcell <- test.data_anno()
      if (input$add.classification_CellTypist_list_2==T) {
        rownames(MainTcell) <- MainTcell$Cell_Index
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("MKI67","TOP2A","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD","MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","S100A9","EBI3","CCL19","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","FCGR3A","S100A13","TLE1","AREG","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","C1QA","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1","Chain2")]
        # -----
        df <- MainTcell_test
        # Vector of columns you want in this data
        nms <- c("MKI67","TOP2A","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD","MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","S100A9","EBI3","CCL19","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","FCGR3A","S100A13","TLE1","AREG","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","C1QA","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1","Chain2")

        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]
        MainTcell_test[is.na(MainTcell_test)] <- 0
        centre_size <- round(dim(MainTcell_test)[1]/10,0)

        if(centre_size>500) {
          centre_size=500
        }

        else {

          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }
        set.seed(123)
        kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 1)
        df_centres <- as.data.frame(kmeans10$centers)

        df_centres <- df_centres %>%
          mutate(CellTypist_list_lower_T_Cell = case_when(
            #### B cells -----
            FCRL2>0 & ITGAX>0 & TBX21>0     ~ "", # B cells
            CD79A>0 & MS4A1>0 & CD19>0      ~ "", # B cells
            CXCR5>0 & TNFRSF13B>0 & CD22>0  ~ "", # B cells
            POU2AF1>0 & CD40>0 & SUGCT>0    ~ "", # B cells
            CR2>0 & CD27>0 & MS4A1>0        ~ "", # B cells
            IGHM>0 & IGHD>0 & TCL1A>0       ~ "", # B cells
            MKI67>0 & SUGCT>0 & AICDA>0     ~ "", # B cells
            CD24>0 & MYO1C>0 & MS4A1>0      ~ "", # B cells
            MME>0 & CD24>0 & MKI67>0        ~ "", # B cells
            IL7R>0 & ZCCHC7>0 & RAG1>0      ~ "", # B cells
            MME>0 &  DNTT>0 & GLL1>0        ~ "", # B cells
            MME>0 & CD24>0 & IGLL5>0        ~ "", # B cells
            JCHAIN>0 & MZB1>0 & XBP1>0      ~ "", # B cells

            #### T cells ----
            ACY3>0 & CD34>0 & SPINK2>0 ~ "ETP",
            FXYD2>0 & HES1>0 & CD99>0 ~ "DN thymocytes",
            CD1A>0 & CD8A>0 & SMPD3>0 ~ "DP thymocytes",
            ZNF683 >0 & GNG4 >0 & PDCD1 >0 ~ "CD8aa T cells",
            TOX2 >0 & SATB1 >0 & CCR9 >0 ~ "CD8a/b(entry)",
            PDCD1 >0 & ICOS >0 & CXCR5>0  ~"Follicular helper T cells",
            KLRB1 >0 & SLC4A10 >0 & `TRAV1-2`>0  ~ "KLRB1+ MAIT cells (TRAV1-2+)",
            GZMK >0 & CD4>0 & IL10>0 ~ "Memory CD4+ cytotoxic T cells",
            NKG7>0 & GNLY>0 & CD8A>0 ~ "NK CD8+ T cells",
            CTLA4>0 & IL2RA >0 & FOXP3>0 ~ "Regulatory T cells",
            MIR155HG>0 & BIRC3>0 & SMS>0 ~ "T(agonist)",
            CD8A>0 & CCR7>0 & SELL>0 ~ "Tcm/Naive cytotoxic T cells",
            CD4>0 & CCR7>0 & SELL>0 ~ "Tcm/Naive helper T cells",
            KLRB1>0 & AQP3>0 & ITGB1>0 ~ "Tem/Effector helper T cells",
            PDCD1>0 & CD4>0 & CTLA4>0 ~ "Tem/Effector helper T cells PD1+",
            CX3CR1>0 & GZMB>0 & GNLY>0 ~ "Tem/Temra cytotoxic T cells",
            GZMK>0 & CD8A>0 & CCL5>0 ~ "Tem/Trm cytotoxic CD8+ T cells",
            CD27>0 & CCR7>0 & IKZF2>0 ~ "Treg(diff)",
            ITGA1>0 & ITGAE>0 & CXCR6>0 ~"Trm cytotoxic T cells",
            CCL5>0 & CXCR3>0 & TBX21>0  ~ "Type 1 helper T cells",
            IL7R>0 & CCR6>0 & ZBTB16>0 ~ "Type 17 helper T cells",
            TRDC>0 & TRGC1>0 & CCL5>0 ~ "Gamma-delta T cells",

            #### DC -----
            CD1C>0 & FCER1A>0 & CLEC10A>0 ~ "", # dc
            BATF3>0 & CADM1>0 & CLEC9A>0 ~ "", # dc
            CLEC10A>0 & FCER1A>0 & CD1C>0 ~ "", # dc
            CLEC10A>0 & FCER1A>0 & S100A9>0 ~ "", # dc
            EBI3>0 & CCR7>0 & CCL19>0 ~ "", # dc
            CLEC10A>0 & KLF4>0 & AXL>0 ~ "", # dc
            IRF8>0 & CLEC10A>0 & PCLAF>0 ~ "", # dc
            #### ILC cells ----
            GNLY>0 & FCGR3A>0 & NKG7>0 ~ "", # dc
            GNLY>0 & CD160>0 & NKG7>0 ~ "", # dc
            S100A13>0 & TLE1>0 & AREG>0 ~ "", # dc
            CXCR3>0 & CD3D>0 & IKZF3>0 ~ "", # dc
            GATA3>0 & KLRG1>0 & HPGDS>0 ~ "", # dc
            IL4I1>0 & RORC>0 & KIT>0 ~ "", # dc
            GNLY>0 & XCL2>0 & NKG7>0 ~ "", # dc
            CCL5>0 & GZMK>0 & FCGR3A>0 ~ "", # dc
            #### Monocytes -----
            TYROBP>0 & C1QC>0 & HMOX1>0 ~ "", # dc
            LYZ>0 & VCAN>0 & S100A9>0 ~ "", # dc
            S100A9>0 & CD14>0 & S100A12>0 ~ "", # dc
            S100A9>0 & LYZ>0 & FCN1>0 ~ "", # dc
            FCGR3A>0 & C1QA>0 & CX3CR1>0 ~ "", # dc


            TRUE ~ NA_character_))

        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","CellTypist_list_lower_T_Cell")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        Cluster_Df
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","CellTypist_list_lower_T_Cell")]
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }


    })
    add.classification_CellTypist_cycling_df_2 <- reactive({
      MainTcell <- test.data_anno()
      if (input$add.classification_CellTypist_cycling_2==T) {
        rownames(MainTcell) <- MainTcell$Cell_Index
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("MKI67","TOP2A","CD3D","CD3E","CD3G","CLEC10A","GNLY","S100A9","TRAC","TRDC")]

        df <- MainTcell_test
        # Vector of columns you want in this data
        nms <- c("MKI67","TOP2A","CD3D","CD3E","CD3G","CLEC10A","GNLY","S100A9","TRAC","TRDC")
        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]

        MainTcell_test[is.na(MainTcell_test)] <- 0
        centre_size <- round(dim(MainTcell_test)[1]/10,0)

        if(centre_size>500) {
          centre_size=500
        }

        else {

          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }
        set.seed(123)

        kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 1)
        df_centres <- as.data.frame(kmeans10$centers)

        df_centres <- df_centres %>%
          mutate(CellTypist_list_Cell_cycling = case_when(
            MKI67 >0 & TOP2A > 0 & TRAC==1 ~ "Cycling ab T cells",
            MKI67 >0 & TOP2A > 0 & TRDC== -1 ~ "Cycling gd T cells",
            MKI67 >0 & TOP2A > 0 & c(CD3D>0 | CD3E>0 | CD3G>0)  ~ "Cycling T cell",
            MKI67 >0 & TOP2A>0 & CLEC10A>0 ~ "Cycling DC",
            MKI67 >0 & TOP2A>0 & GNLY>0 ~ "Cycling NK cells",
            MKI67 >0 & TOP2A>0 & S100A9>0 ~ "Cycling Monocytes",
            TRUE ~ NA_character_))

        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","CellTypist_list_Cell_cycling")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        Cluster_Df
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","CellTypist_list_Cell_cycling")]
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }


    })
    add.classification_CellTypist_list_lower_curated_2 <- reactive({
      MainTcell <- test.data_anno()
      if (input$add.classification_CellTypist_list_overview_3==T) {
        rownames(MainTcell) <- MainTcell$Cell_Index
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("MKI67","TOP2A","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD","MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","S100A9","EBI3","CCL19","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","FCGR3A","S100A13","TLE1","AREG","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","C1QA","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1","Chain2")]
        # -----
        df <- MainTcell_test
        # Vector of columns you want in this data
        nms <- c("MKI67","TOP2A","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD","MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","S100A9","EBI3","CCL19","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","FCGR3A","S100A13","TLE1","AREG","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","C1QA","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1","Chain2")

        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]
        MainTcell_test[is.na(MainTcell_test)] <- 0
        centre_size <- round(dim(MainTcell_test)[1]/10,0)
        if(centre_size>500) {
          centre_size=500
        }

        else {

          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }
        set.seed(123)
        kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 1)
        df_centres <- as.data.frame(kmeans10$centers)

        df_centres <- df_centres %>%
          mutate(CellTypist_list_lower_curated = case_when(
            #### T cells ----
            ACY3>0 & CD34>0 & SPINK2>0 ~ "ETP",
            FXYD2>0 & HES1>0 & CD99>0 ~ "DN thymocytes",
            CD1A>0 & CD8A>0 & SMPD3>0 ~ "DP thymocytes",
            ZNF683 >0 & GNG4 >0 & PDCD1 >0 ~ "CD8aa T cells",
            TOX2 >0 & SATB1 >0 & CCR9 >0 ~ "CD8a/b(entry)",
            PDCD1 >0 & ICOS >0 & CXCR5>0  ~"Follicular helper T cells",
            KLRB1 >0 & SLC4A10 >0 & `TRAV1-2`>0  ~ "KLRB1+ MAIT cells (TRAV1-2+)",
            KLRB1 >0 & `TRAV1-2`>0 & CD8A>0 ~ "KLRB1+ MAIT cells (TRAV1-2+)",
            GZMK >0 & CD4>0 & IL10>0 ~ "Memory CD4+ cytotoxic T cells",
            NKG7>0 & GNLY>0 & CD8A>0 ~ "NK CD8+ T cells",
            CTLA4>0 & IL2RA >0 & FOXP3>0 ~ "Regulatory T cells",
            MIR155HG>0 & BIRC3>0 & SMS>0 ~ "T(agonist)",
            CD8A>0 & CCR7>0 & SELL>0 ~ "Tcm/Naive cytotoxic T cells",
            CD4>0 & CCR7>0 & SELL>0 ~ "Tcm/Naive helper T cells",
            KLRB1>0 & AQP3>0 & ITGB1>0 ~ "Tem/Effector helper T cells",
            PDCD1>0 & CD4>0 & CTLA4>0 ~ "Tem/Effector helper T cells PD1+",
            CX3CR1>0 & GZMB>0 & GNLY>0 ~ "Tem/Temra cytotoxic T cells",
            GZMK>0 & CD8A>0 & CCL5>0 ~ "Tem/Trm cytotoxic CD8+ T cells",
            CD27>0 & CCR7>0 & IKZF2>0 ~ "Treg(diff)",
            ITGA1>0 & ITGAE>0 & CXCR6>0 ~"Trm cytotoxic T cells",
            CCL5>0 & CXCR3>0 & TBX21>0  ~ "Th1 T cells",
            IL7R>0 & CCR6>0 & ZBTB16>0 ~ "Th17 T cells",
            # my annotation
            TRDC>0 & TRGC1>0 & CCL5>0 ~ "Gamma-delta T cells (CCL5+)",
            TRDC>0 & TRGC1>0 & TOX2>0 ~ "Gamma-delta T cells (TOX2+)",
            CD4>0 & CD8A>0 & Chain2>0 ~ "DP thymocytes",
            CD4>0 & FOXP3>0 & CD27>0 ~ "Regulatory T cells",
            CD4>0 &	NKG7>0 &	GNLY>0 ~ "CD4+NKGL7+GNLY+ T cell",
            CD4>0 & GNLY>0  ~ "CD4+GNLY+ T cell",
            CD4>0 & CCL5>0 ~ "Th1 T cells",
            CD4>0 & IL10>0 ~ "IL10+CD4+ T cell",
            CD4>0 & MKI67>0 &	TOP2A>0 ~ "Cycling CD4+ T cell",
            CD4>0 & c(CD3E>0 | CD3D>0 | CD3G>0) ~ "CD4+ T cell",
            CD4>0 & Chain2>0 ~ "CD4+ T cell",
            CD8A>0 & GNLY>0  ~ "CD8+GNLY+ T cell",
            TRDC>0 & TRGC1>0 ~ "Gamma-delta T cells",
            CD4<0 & CD8A<0 & GNLY>0 & Chain2>0 ~ "DN GNLY+ T cell",
            CD4<0 & CD8A<0 & IL2RA>0 & Chain2>0 ~ "DN Tregs cell", #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3230797/
            CD4<0 & CD8A<0 & Chain2>0 ~ "DN T cells",


            #### B cells -----
            FCRL2>0 & ITGAX>0 & TBX21>0     ~ "Age-associated B cells", # B cells
            CD79A>0 & MS4A1>0 & CD19>0      ~ "B cell", # B cell
            CXCR5>0 & TNFRSF13B>0 & CD22>0  ~ "Follicular B cells",
            POU2AF1>0 & CD40>0 & SUGCT>0    ~ "Germinal center B cells",
            CR2>0 & CD27>0 & MS4A1>0        ~ "Memory B cells",
            IGHM>0 & IGHD>0 & TCL1A>0       ~ "Naive B cells",
            MKI67>0 & SUGCT>0 & AICDA>0     ~ "Proliferative germinal center B cells",
            CD24>0 & MYO1C>0 & MS4A1>0      ~ "Transitional B cells",
            MME>0 & CD24>0 & MKI67>0        ~ "Large pre-B cells",
            IL7R>0 & ZCCHC7>0 & RAG1>0      ~ "Pre-pro-B cells",
            MME>0 &  DNTT>0 & GLL1>0        ~ "Pro-B cells",
            MME>0 & CD24>0 & IGLL5>0        ~ "Small pre-B cells",
            JCHAIN>0 & MZB1>0 & XBP1>0      ~ "Plasma cell",
            #### DC -----
            CD1C>0 & FCER1A>0 & CLEC10A>0 ~ "DC",
            BATF3>0 & CADM1>0 & CLEC9A>0 ~ "DC1",
            CLEC10A>0 & FCER1A>0 & CD1C>0 ~ "DC2",
            CLEC10A>0 & FCER1A>0 & S100A9>0 ~ "DC3",
            EBI3>0 & CCR7>0 & CCL19>0 ~ "Migratory DC",
            CLEC10A>0 & KLF4>0 & AXL>0 ~ " Transitional DC",
            IRF8>0 & CLEC10A>0 & PCLAF>0 ~ "DC precursor",
            #### ILC cells ----
            GNLY>0 & FCGR3A>0 & NKG7>0 ~ "CD16+ NK cells",
            GNLY>0 & CD160>0 & NKG7>0 ~ "CD16- NK cells",
            S100A13>0 & TLE1>0 & AREG>0 ~ "ILC",
            CXCR3>0 & CD3D>0 & IKZF3>0 ~ "ILC1",
            GATA3>0 & KLRG1>0 & HPGDS>0 ~ "ILC2",
            IL4I1>0 & RORC>0 & KIT>0 ~ "ILC3",
            GNLY>0 & XCL2>0 & NKG7>0 ~ "NK cells",
            CCL5>0 & GZMK>0 & FCGR3A>0 ~ "Traditional NK cells",
            #### Monocytes -----
            TYROBP>0 & C1QC>0 & HMOX1>0 ~ "Mono-mac",
            LYZ>0 & VCAN>0 & S100A9>0 ~ "Monocyte precursor",
            S100A9>0 & CD14>0 & S100A12>0 ~ "Classical monocytes",
            S100A9>0 & LYZ>0 & FCN1>0 ~ "Monocytes",
            FCGR3A>0 & C1QA>0 & CX3CR1>0 ~ "Non-classical monocytes",


            TRUE ~ NA_character_))

        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","CellTypist_list_lower_curated")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        Cluster_Df
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","CellTypist_list_lower_curated")]
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }


    })
    add.classification_CellTypist_list_CD4_lower_curated_2 <- reactive({
      MainTcell <- test.data_anno()
      if (input$add.classification_CellTypist_list_overview_3==T) {
        rownames(MainTcell) <- MainTcell$Cell_Index
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("MKI67","TOP2A","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD","MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","S100A9","EBI3","CCL19","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","FCGR3A","S100A13","TLE1","AREG","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","C1QA","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1","Chain2")]
        # -----
        df <- MainTcell_test
        # Vector of columns you want in this data
        nms <- c("MKI67","TOP2A","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD","MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","S100A9","EBI3","CCL19","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","FCGR3A","S100A13","TLE1","AREG","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","C1QA","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1","Chain2")

        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]
        MainTcell_test[is.na(MainTcell_test)] <- 0
        centre_size <- round(dim(MainTcell_test)[1]/10,0)
        if(centre_size>500) {
          centre_size=500
        }
        else {

          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }
        set.seed(123)
        kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 1)
        df_centres <- as.data.frame(kmeans10$centers)

        df_centres <- df_centres %>%
          mutate(CellTypist_list_CD4_lower_curated = case_when(
            #### T cells ----
            ACY3>0 & CD34>0 & SPINK2>0 ~ "Non-CD4+ T cell",
            FXYD2>0 & HES1>0 & CD99>0 ~ "Non-CD4+ T cell",
            CD1A>0 & CD8A>0 & SMPD3>0 ~ "Non-CD4+ T cell",
            ZNF683 >0 & GNG4 >0 & PDCD1 >0 ~ "Non-CD4+ T cell",
            TOX2 >0 & SATB1 >0 & CCR9 >0 ~ "Non-CD4+ T cell",
            PDCD1 >0 & ICOS >0 & CXCR5>0  ~"Non-CD4+ T cell",
            KLRB1 >0 & SLC4A10 >0 & `TRAV1-2`>0  ~ "Non-CD4+ T cell",
            KLRB1 >0 & `TRAV1-2`>0 & CD8A>0 ~ "Non-CD4+ T cell",
            GZMK >0 & CD4>0 & IL10>0 ~ "Memory CD4+ cytotoxic T cells",
            NKG7>0 & GNLY>0 & CD8A>0 ~ "Non-CD4+ T cell",
            CTLA4>0 & IL2RA >0 & FOXP3>0 ~ "Regulatory T cells",
            MIR155HG>0 & BIRC3>0 & SMS>0 ~ "Non-CD4+ T cell",
            CD8A>0 & CCR7>0 & SELL>0 ~ "Non-CD4+ T cell",
            CD4>0 & CCR7>0 & SELL>0 ~ "Tcm/Naive helper T cells",
            KLRB1>0 & AQP3>0 & ITGB1>0 ~ "Tem/Effector helper T cells",
            PDCD1>0 & CD4>0 & CTLA4>0 ~ "Tem/Effector helper T cells PD1+",
            CX3CR1>0 & GZMB>0 & GNLY>0 ~ "Tem/Temra cytotoxic T cells",
            GZMK>0 & CD8A>0 & CCL5>0 ~ "Non-CD4+ T cell",
            CD27>0 & CCR7>0 & IKZF2>0 ~ "Treg(diff)",
            ITGA1>0 & ITGAE>0 & CXCR6>0 ~"Non-CD4+ T cell",
            CCL5>0 & CXCR3>0 & TBX21>0  ~ "Th1 T cells",
            IL7R>0 & CCR6>0 & ZBTB16>0 ~ "Th17 T cells",
            # my annotation
            TRDC>0 & TRGC1>0 & CCL5>0 ~ "Non-CD4+ T cell",
            TRDC>0 & TRGC1>0 & TOX2>0 ~ "Non-CD4+ T cell",
            CD4>0 & CD8A>0 & Chain2>0 ~ "Non-CD4+ T cell",
            CD4>0 & FOXP3>0 & CD27>0 ~ "Regulatory T cells",
            CD4>0 &	NKG7>0 &	GNLY>0 ~ "CD4+NKGL7+GNLY+ T cell",
            CD4>0 & GNLY>0  ~ "CD4+GNLY+ T cell",
            CD4>0 & CCL5>0 ~ "Th1 T cells",
            CD4>0 & IL10>0 ~ "IL10+CD4+ T cell",
            CD4>0 & MKI67>0 &	TOP2A>0 ~ "Cycling CD4+ T cell",
            CD4>0 & c(CD3E>0 | CD3D>0 | CD3G>0) ~ "CD4+ T cell",
            CD4>0 & Chain2>0 ~ "CD4+ T cell",
            CD8A>0 & GNLY>0  ~ "Non-CD4+ T cell",
            TRDC>0 & TRGC1>0 ~ "Non-CD4+ T cell",
            CD4<0 & CD8A<0 & GNLY>0 & Chain2>0 ~ "Non-CD4+ T cell",
            CD4<0 & CD8A<0 & IL2RA>0 & Chain2>0 ~ "Non-CD4+ T cell",
            CD4<0 & CD8A<0 & Chain2>0 ~ "Non-CD4+ T cell",


            #### B cells -----
            FCRL2>0 & ITGAX>0 & TBX21>0     ~ "Non-CD4+ T cell",
            CD79A>0 & MS4A1>0 & CD19>0      ~ "Non-CD4+ T cell",
            CXCR5>0 & TNFRSF13B>0 & CD22>0  ~ "Non-CD4+ T cell",
            POU2AF1>0 & CD40>0 & SUGCT>0    ~ "Non-CD4+ T cell",
            CR2>0 & CD27>0 & MS4A1>0        ~ "Non-CD4+ T cell",
            IGHM>0 & IGHD>0 & TCL1A>0       ~ "Non-CD4+ T cell",
            MKI67>0 & SUGCT>0 & AICDA>0     ~ "Non-CD4+ T cell",
            CD24>0 & MYO1C>0 & MS4A1>0      ~ "Non-CD4+ T cell",
            MME>0 & CD24>0 & MKI67>0        ~ "Non-CD4+ T cell",
            IL7R>0 & ZCCHC7>0 & RAG1>0      ~ "Non-CD4+ T cell",
            MME>0 &  DNTT>0 & GLL1>0        ~ "Non-CD4+ T cell",
            MME>0 & CD24>0 & IGLL5>0        ~ "Non-CD4+ T cell",
            JCHAIN>0 & MZB1>0 & XBP1>0      ~ "Non-CD4+ T cell",
            #### DC -----
            CD1C>0 & FCER1A>0 & CLEC10A>0 ~"Non-CD4+ T cell",
            BATF3>0 & CADM1>0 & CLEC9A>0 ~ "Non-CD4+ T cell",
            CLEC10A>0 & FCER1A>0 & CD1C>0 ~ "Non-CD4+ T cell",
            CLEC10A>0 & FCER1A>0 & S100A9>0 ~ "Non-CD4+ T cell",
            EBI3>0 & CCR7>0 & CCL19>0 ~ "Non-CD4+ T cell",
            CLEC10A>0 & KLF4>0 & AXL>0 ~ "Non-CD4+ T cell",
            IRF8>0 & CLEC10A>0 & PCLAF>0 ~ "Non-CD4+ T cell",
            #### ILC cells ----
            GNLY>0 & FCGR3A>0 & NKG7>0 ~ "Non-CD4+ T cell",
            GNLY>0 & CD160>0 & NKG7>0 ~ "Non-CD4+ T cell",
            S100A13>0 & TLE1>0 & AREG>0 ~ "Non-CD4+ T cell",
            CXCR3>0 & CD3D>0 & IKZF3>0 ~ "Non-CD4+ T cell",
            GATA3>0 & KLRG1>0 & HPGDS>0 ~ "Non-CD4+ T cell",
            IL4I1>0 & RORC>0 & KIT>0 ~ "Non-CD4+ T cell",
            GNLY>0 & XCL2>0 & NKG7>0 ~ "Non-CD4+ T cell",
            CCL5>0 & GZMK>0 & FCGR3A>0 ~ "Non-CD4+ T cell",
            #### Monocytes -----
            TYROBP>0 & C1QC>0 & HMOX1>0 ~ "Non-CD4+ T cell",
            LYZ>0 & VCAN>0 & S100A9>0 ~ "Non-CD4+ T cell",
            S100A9>0 & CD14>0 & S100A12>0 ~ "Non-CD4+ T cell",
            S100A9>0 & LYZ>0 & FCN1>0 ~ "Non-CD4+ T cell",
            FCGR3A>0 & C1QA>0 & CX3CR1>0 ~ "Non-CD4+ T cell",


            TRUE ~ NA_character_))

        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","CellTypist_list_CD4_lower_curated")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        Cluster_Df
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","CellTypist_list_CD4_lower_curated")]
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }


    })
    add.classification_CellTypist_list_curated_T_cell_2 <- reactive({
      MainTcell <- test.data_anno()
      if (input$add.classification_CellTypist_list_overview_3==T) {
        rownames(MainTcell) <- MainTcell$Cell_Index
        MainTcell_test <- MainTcell[,names(MainTcell) %in% c("MKI67","TOP2A","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD","MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","S100A9","EBI3","CCL19","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","FCGR3A","S100A13","TLE1","AREG","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","C1QA","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1","Chain2")]
        # -----
        df <- MainTcell_test
        # Vector of columns you want in this data
        nms <- c("MKI67","TOP2A","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD","MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","S100A9","EBI3","CCL19","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","FCGR3A","S100A13","TLE1","AREG","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","C1QA","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1","Chain2")

        Missing <- setdiff(nms, names(df))  # Find names of missing columns
        df[Missing] <- 0                    # Add them, filled with '0's
        MainTcell_test <- df[nms]
        MainTcell_test[is.na(MainTcell_test)] <- 0
        centre_size <- round(dim(MainTcell_test)[1]/10,0)
        if(centre_size>500) {
          centre_size=500
        }

        else {

          centre_size <- round(dim(MainTcell_test)[1]/10,0)
        }
        set.seed(123)
        kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 1)
        df_centres <- as.data.frame(kmeans10$centers)

        df_centres <- df_centres %>%
          mutate(CellTypist_list_curated_T_cell = case_when(
            #### T cells ----
            ACY3>0 & CD34>0 & SPINK2>0 ~ "ETP",
            FXYD2>0 & HES1>0 & CD99>0 ~ "DN thymocytes",
            CD1A>0 & CD8A>0 & SMPD3>0 ~ "DP thymocytes",
            ZNF683 >0 & GNG4 >0 & PDCD1 >0 ~ "CD8aa T cells",
            TOX2 >0 & SATB1 >0 & CCR9 >0 ~ "CD8a/b(entry)",
            PDCD1 >0 & ICOS >0 & CXCR5>0  ~"Follicular helper T cells",
            KLRB1 >0 & SLC4A10 >0 & `TRAV1-2`>0  ~ "KLRB1+ MAIT cells (TRAV1-2+)",
            KLRB1 >0 & `TRAV1-2`>0 & CD8A>0 ~ "KLRB1+ MAIT cells (TRAV1-2+)",
            GZMK >0 & CD4>0 & IL10>0 ~ "Memory CD4+ cytotoxic T cells",
            NKG7>0 & GNLY>0 & CD8A>0 ~ "NK CD8+ T cells",
            CTLA4>0 & IL2RA >0 & FOXP3>0 ~ "Regulatory T cells",
            MIR155HG>0 & BIRC3>0 & SMS>0 ~ "T(agonist)",
            CD8A>0 & CCR7>0 & SELL>0 ~ "Tcm/Naive cytotoxic T cells",
            CD4>0 & CCR7>0 & SELL>0 ~ "Tcm/Naive helper T cells",
            KLRB1>0 & AQP3>0 & ITGB1>0 ~ "Tem/Effector helper T cells",
            PDCD1>0 & CD4>0 & CTLA4>0 ~ "Tem/Effector helper T cells PD1+",
            CX3CR1>0 & GZMB>0 & GNLY>0 ~ "Tem/Temra cytotoxic T cells",
            GZMK>0 & CD8A>0 & CCL5>0 ~ "Tem/Trm cytotoxic CD8+ T cells",
            CD27>0 & CCR7>0 & IKZF2>0 ~ "Treg(diff)",
            ITGA1>0 & ITGAE>0 & CXCR6>0 ~"Trm cytotoxic T cells",
            CCL5>0 & CXCR3>0 & TBX21>0  ~ "Th1 T cells",
            IL7R>0 & CCR6>0 & ZBTB16>0 ~ "Th17 T cells",
            # my annotation
            TRDC>0 & TRGC1>0 & CCL5>0 ~ "Gamma-delta T cells (CCL5+)",
            TRDC>0 & TRGC1>0 & TOX2>0 ~ "Gamma-delta T cells (TOX2+)",
            CD4>0 & CD8A>0 & Chain2>0 ~ "DP thymocytes",
            CD4>0 & FOXP3>0 & CD27>0 ~ "Regulatory T cells",
            CD4>0 &	NKG7>0 &	GNLY>0 ~ "CD4+NKGL7+GNLY+ T cell",
            CD4>0 & GNLY>0  ~ "CD4+GNLY+ T cell",
            CD4>0 & CCL5>0 ~ "Th1 T cells",
            CD4>0 & IL10>0 ~ "IL10+CD4+ T cell",
            CD4>0 & MKI67>0 &	TOP2A>0 ~ "Cycling CD4+ T cell",
            CD4>0 & c(CD3E>0 | CD3D>0 | CD3G>0) ~ "CD4+ T cell",
            CD4>0 & Chain2>0 ~ "CD4+ T cell",
            CD8A>0 & GNLY>0  ~ "CD8+GNLY+ T cell",
            TRDC>0 & TRGC1>0 ~ "Gamma-delta T cells",
            CD4<0 & CD8A<0 & GNLY>0 & Chain2>0 ~ "DN GNLY+ T cell",
            CD4<0 & CD8A<0 & IL2RA>0 & Chain2>0 ~ "DN Tregs cell", #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3230797/
            CD4<0 & CD8A<0 & Chain2>0 ~ "DN T cells",

            TRUE ~ NA_character_))

        df_centres$clust_name <- rownames(df_centres)
        head(df_centres)
        df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","CellTypist_list_curated_T_cell")]
        Cluster_Df <- as.data.frame( kmeans10$cluster)
        Cluster_Df
        names(Cluster_Df) <- "clust_name"
        Cluster_Df$Cell_Index <- rownames(Cluster_Df)
        Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
        Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","CellTypist_list_curated_T_cell")]
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }
      else {
        Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
        names(Cluster_Df_class2) <- "Cell_Index"
        Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      }


    })
    Vals_norm3 <- reactiveValues(Norm3=NULL)
    anno_final_kmeans <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Run K-means")
      )

      annotation <- cbind(add.classification_T_cell_Df_2(),
                          add.classification_T.B_cell_Function_Df_2(),
                          add.classification_T_cell_Function_CD4_df_2(),
                          add.classification_T_cell_Function_CD8_df_2(),
                          add.classification_T_cell_Function_CD4_CD8pos_df_2(),
                          add.classification_T_cell_Function_CD4_CD8neg_df_2(),
                          add.classification_B_cell_Function_Df_2(),
                          add.classification_T_cell_Activation_df_2(),
                          add.classification_T_cell_Memory_df_2(),
                          add.classification_CellTypist_list_higher_2(),
                          add.classification_CellTypist_list_lower_2(),
                          add.classification_CellTypist_list_lower_Tcell(),
                          add.classification_CellTypist_cycling_df_2(),
                          add.classification_CellTypist_list_lower_curated_2(),
                          add.classification_CellTypist_list_curated_T_cell_2(),
                          add.classification_CellTypist_list_CD4_lower_curated_2()
      )

      rownames(annotation) <- annotation$Cell_Index
      annotation <- annotation[,!grepl("Cell_Index",names(annotation))]
      annotation <- annotation[,!grepl("add.classi",names(annotation))]
      annotation$Cell_Index <- rownames(annotation)
      annotation
    })

    output$DEx_table_TcellClass_2 <-  DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "K-means processing")
      )
      anno_final_kmeans()

    })

    # scGATE annotations -------
    scGATE_anno_generic <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )

      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      if (input$generic_scGATE==T) {
        scGate_models_DB <- custom_db_scGATE(system.file("scGATE","human/generic",package = "STEGO.R"))
        models.list <- scGate_models_DB[c(input$GenericID_scGATE)]
        obj <- scGate(sc, model = models.list)
        df$generic <- obj@meta.data$scGate_multi
      }
      else {
        df
      }
      df
    })
    scGATE_anno_CD4 <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      if (input$CD4_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/CD4_TIL",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list)
        df$CD4_TIL <- obj@meta.data$scGate_multi
      }
      else {
        df
      }
      df
    })
    scGATE_anno_CD8 <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      if (input$CD8_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/CD8_TIL",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list)
        df$CD8_TIL <- obj@meta.data$scGate_multi
      }
      else {
        df
      }
      df
    })
    scGATE_anno_cellTypist_higher <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      if (input$cellTypist_higher_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/CellTypist_higher",package = "STEGO.R")))
        models.list <- scGate_models_DB[c("Bcell","CD4Tcell","CD8Tcell","DC","DNTcell", "MonoMac","NK","Plasma")]
        obj <- scGate(sc, model = models.list)
        df$CellTypist_higher <- obj@meta.data$scGate_multi
      }
      else {
        df
      }
      df
    })
    scGATE_anno_cellTypist_lowerBcells <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      if (input$cellTypist_lowerBcells_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/CellTypist_lower_B_cells",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list)
        df$CellTypist_lower_B_cells <- obj@meta.data$scGate_multi
      }
      else {
        df
      }
      df
    })
    scGATE_anno_cellTypist_lowerCD4T <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      if (input$cellTypist_lowerCD4T_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/CellTypist_lower_CD4",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list)
        df$CellTypist_lower_CD4 <- obj@meta.data$scGate_multi
      }
      else {
        df
      }
      df
    })
    scGATE_anno_cellTypist_lowerCD8T <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      if (input$cellTypist_lowerCD8T_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/CellTypist_lower_CD8",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list)
        df$CellTypist_lower_CD8 <- obj@meta.data$scGate_multi
      }
      else {
        df
      }
      df
    })
    scGATE_anno_cellTypist_lowerOtherT <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      if (input$cellTypist_lowerOtherT_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/CellTypist_lower_Other",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list)
        df$CellTypist_lower_Other <- obj@meta.data$scGate_multi
      }
      else {
        df
      }
      df
    })
    scGATE_anno_Ex.Sen <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      if (input$Ex.Sen_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/Exhausted_Senescence",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list)
        df$Exhausted_Senescence <- obj@meta.data$scGate_multi
      }
      else {
        df
      }
      df
    })
    scGATE_anno_COVID <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      if (input$COVID_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/COVID",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list)
        df$COVID <- obj@meta.data$scGate_multi
      }
      else {
        df
      }
      df
    })
    scGATE_anno_Activation <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      if (input$Activation_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/Activation",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list)
        df$Activation <- obj@meta.data$scGate_multi
      }
      else {
        df
      }
      df
    })
    scGATE_anno_IFNgTNFa <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      if (input$IFNgTNFa_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/IFNgTNFa",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list)
        df$IFNgTNFa <- obj@meta.data$scGate_multi
      }
      else {
        df
      }
      df
    })
    scGATE_anno_GNLY.PFR1.GZMB <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      if (input$GNLY.PFR1.GZMB_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/GNLY.PFR1.GZMB",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list)
        df$GNLY.PFR1.GZMB <- obj@meta.data$scGate_multi
      }
      else {
        df
      }
      df
    })
    scGATE_anno_Interlukin <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
      names(df) <- "Cell_Index"
      if (input$Interlukin_scGATE==T) {
        scGate_models_DB <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/Interlukins",package = "STEGO.R")))
        models.list <- scGate_models_DB
        obj <- scGate(sc, model = models.list)
        df$Interlukins <- obj@meta.data$scGate_multi
      }
      else {
        df
      }
      df
    })

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
    output$scGATE_verbatum_cellTypist_higher <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$cellTypist_higher_scGATE==T) {
        scGATE_anno_cellTypist_higher()}
      else{
        print("cellTypist_higher not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    output$scGATE_verbatum_cellTypist_lower_Bcells <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$cellTypist_lowerBcells_scGATE==T) {
        scGATE_anno_cellTypist_lowerBcells()}
      else{
        print("cellTypist lower Bcells not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    output$scGATE_verbatum_cellTypist_lower_CD4T <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$cellTypist_lowerCD4T_scGATE==T) {
        scGATE_anno_cellTypist_lowerCD4T()}
      else{
        print("cellTypist lower CD4 T cells not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    output$scGATE_verbatum_cellTypist_lower_CD8T <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$cellTypist_lowerCD8T_scGATE==T) {
        scGATE_anno_cellTypist_lowerCD8T()}
      else{
        print("cellTypist lower CD8 T cells not run")
      }
      sink(type = "message")
      sink(type = "output")
      cat(readLines(FN), sep="\n")
    })
    output$scGATE_verbatum_cellTypist_lower_OtherT <- renderPrint({
      FN <- tempfile()
      zz <- file(FN, open = "wt")
      sink(zz ,type = "output")
      sink(zz, type = "message")
      if (input$cellTypist_lowerOtherT_scGATE==T) {
        scGATE_anno_cellTypist_lowerOtherT()}
      else{
        print("cellTypist lower Other T cells not run")
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
      if (input$Ex.Sen_scGATE==T) {
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

    scGATE_anno <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      df <- cbind(scGATE_anno_generic(),
                  scGATE_anno_CD4(),
                  scGATE_anno_CD8(),
                  scGATE_anno_cellTypist_lowerBcells(),
                  scGATE_anno_cellTypist_lowerCD4T(),
                  scGATE_anno_cellTypist_lowerCD8T(),
                  scGATE_anno_cellTypist_lowerOtherT(),
                  scGATE_anno_Ex.Sen(),
                  scGATE_anno_COVID(),
                  scGATE_anno_Activation(),
                  scGATE_anno_IFNgTNFa(),
                  scGATE_anno_GNLY.PFR1.GZMB(),
                  scGATE_anno_Interlukin()
      )

      rownames(df) <- df$Cell_Index
      df <- df[,!grepl("Cell_Index",names(df))]
      df$Cell_Index <- rownames(df)
      as.data.frame(df)
    })

    output$DEx_table_TcellClass_scGATE <-  DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )
      scGATE_anno()
    })
    # all.annotations added -----
    final_anno <- reactive({
      sc <- getData_2()
      validate(
        need(nrow(sc)>0,
             "Upload file for annotation")
      )


      if (input$add.scGATE==T && input$add.kmeans==T) {
        annotation <- merge(scGATE_anno(),anno_final_kmeans(),by = "Cell_Index", sort=F,all=T)
      }

      else if (input$add.scGATE==T && input$add.kmeans==F) {
        annotation <- cbind(scGATE_anno())
        annotation
      }
      else if (input$add.scGATE==F && input$add.kmeans==T) {
        annotation <- cbind(anno_final_kmeans())
        annotation
      }
      else {
        annotation <- as.data.frame(sc@meta.data[,names(sc@meta.data) %in% c("Cell_Index")])
        names(annotation) <- "Cell_Index"
        annotation
      }
      rownames(annotation) <- annotation$Cell_Index
      as.data.frame(annotation)
      md <- sc@meta.data
      md.anno <- merge(md,annotation,by = "Cell_Index", sort=F,all=T)
      sc@meta.data <- md.anno
      sc
    })

    output$DEx_table_TcellClass_scGATE.kmeans <-  DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),
                                                                      options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100),
                                                                                     pageLength = 5, scrollX = TRUE),{
                                                                                       sc <- getData_2()
                                                                                       validate(
                                                                                         need(nrow(sc)>0,
                                                                                              "K-means processing")
                                                                                       )
                                                                                       # sc <- final_anno()
                                                                                       # sc@meta.data
                                                                                       sc <- final_anno()
                                                                                       sc@meta.data
                                                                                     })

    output$downloaddf_SeruatObj_annotated <- downloadHandler(
      filename = function(){
        paste(input$project_name3,"_annotated_",gsub("-", ".", Sys.Date()),".h5Seurat", sep = "")
      },
      content = function(file){
        as.h5Seurat(final_anno(),file)
      } )

    #
    # Analysis -----
    ## uploading seruat obj ----
    input.data_sc_pro <- reactive({
      inFile_sc_pro <- input$file_SC_pro
      if (is.null(inFile_sc_pro)) return(NULL)
      else {

        if (input$STEGO_R_pro=="STEGO_R (.h5Seurat)") {
          dataframe = LoadH5Seurat(inFile_sc_pro$datapath)
        }

        else {
          sc <- readRDS(inFile_sc_pro$datapath)
        }

      }

    })

    # input.data_sc_clusTCR <- reactive({switch(input$dataset_sc_pro,"test_data_sc_pro" = test.data_sc_clusTCR(),"own_data_sc_pro" = own.data_sc_clusTCR())})
    # # input.data_sc_clusTCR <- reactive({switch(input$dataset_cluster_file,"test_data_clusTCR" = test.data_sc_clusTCR(),"own_data_clusTCR" = own.data_sc_clusTCR())})
    # test.data_sc_clusTCR <- reactive({
    #   dataframe = read.csv(system.file("extdata","BDrhap/clusTCR/ClusTCR2.csv",package = "STEGO.R"))
    #
    #   # dataframe = read.csv("../Test.Cluster/Cluster_table  2022.12.15.csv")
    # })
    input.data_sc_clusTCR <- reactive({
      inFile_cluster_file <- input$file_cluster_file
      if (is.null(inFile_cluster_file)) return(NULL)
      else {
        dataframe = read.csv(inFile_cluster_file$datapath)
      }

    })

    # upload TCRex file ----
    # input.data_sc_TCRex <- reactive({switch(input$dataset_sc_pro,"test_data_sc_pro" = test.data_sc_TCRex(),"own_data_sc_pro" = own.data_sc_TCRex())})
    # # input.data_sc_clusTCR <- reactive({switch(input$dataset_cluster_file,"test_data_clusTCR" = test.data_sc_clusTCR(),"own_data_clusTCR" = own.data_sc_clusTCR())})
    # test.data_sc_TCRex <- reactive({
    #
    #   dataframe = read.table(system.file("extdata","BDrhap/TCRex/tcrex_nsjhx29ivo.tsv",package = "STEGO.R"),skip=7,header = T,sep="\t")
    #
    #   # read.table("../Public data/Bd Rhapsody/TCRex/tcrex_nsjhx29ivo.tsv",skip = 7,header = T,sep="\t")
    # })

    input.data_sc_TCRex <- reactive({
      inupload_TCRex_file <- input$upload_TCRex_file
      if (is.null(inupload_TCRex_file)) return(NULL)
      else {
        dataframe = read.table(inupload_TCRex_file$datapath,skip = input$skip_TCRex_up,header = T,sep="\t")
      }

    })

    ## uploaded clusTCR table -----
    output$Tb_ClusTCR_test <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{


      calls <- input.data_sc_clusTCR()
      validate(
        need(nrow(calls)>0,
             "Upload clusTCR table, which is needed for TCR -> UMAP section")
      )

      calls
    })



    output$Tb_tcrex_test <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      calls <- input.data_sc_TCRex()
      validate(
        need(nrow(calls)>0,
             "Upload TCRex table")
      )

      calls
    })

    ### observe ----
    observe({
      sc <- input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )

      df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "Samp_col",
        choices=names(df3.meta),
        selected = "orig.ident")
    })
    #
    observe({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )

      df3.meta <- sc@meta.data


      if (input$datasource == "BD rhapsody" & input$STEGO_R_pro== "STEGO_R (.h5Seurat)") {
        updateSelectInput(
          session,
          "V_gene_sc",
          choices=names(df3.meta),
          selected = "v_gene_cdr3_AB_GD")
      }

      else if (input$datasource == "10x Genomics" & input$STEGO_R_pro== "STEGO_R (.h5Seurat)") {
        updateSelectInput(
          session,
          "V_gene_sc",
          choices=names(df3.meta),
          selected = "vdj_gene_cdr3_AG_BD")
      }

      else if (input$STEGO_R_pro== "Seurat (.rds)") {
        updateSelectInput(
          session,
          "V_gene_sc",
          choices=names(df3.meta),
          selected = "Av_gene")
      }

      else {
        updateSelectInput(
          session,
          "V_gene_sc",
          choices=names(df3.meta),
          selected = "vdj_gene_cdr3_AB_GD")
      }


    })


    output$Tb_TCR_clonotypes.Umap <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      calls <- select_group_metadata()
      validate(
        need(nrow(calls)>0,
             error_message_val_UMAP)
      )

      UMAP.wt.clonality2 <- For_col_top()
      colorblind_vector <-as.data.frame(unlist(colors_UMAP_Topclonotypes()))

      topclones_col <- as.data.frame(unique(UMAP.wt.clonality2$topclones))
      names(topclones_col) <- "topclones"
      # topclones_col$col <- colorblind_vector()
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
    TCR_Expansion <- reactive({
      sc <- input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )

      df3.meta <- as.data.frame(sc@meta.data)

      if (input$STEGO_R_pro== "Seurat (.rds)") {
        df3.meta$Cell_Index <- rownames(df3.meta)
        df3.meta$v_gene_AG <- df3.meta[,names(df3.meta) %in% input$RDS_V_gene_A]
        df3.meta$v_gene_BD <- df3.meta[,names(df3.meta) %in% input$RDS_V_gene_B]
        df3.meta$cdr3_AG <- df3.meta[,names(df3.meta) %in% input$RDS_cdr3_A]
        df3.meta$cdr3_BD <- df3.meta[,names(df3.meta) %in% input$RDS_cdr3_B]
        df3.meta$v_gene_cdr3_BD <- paste(df3.meta$v_gene_BD,df3.meta$cdr3_BD,sep="_")
        df3.meta$v_gene_cdr3_AG<- paste(df3.meta$v_gene_AG,df3.meta$cdr3_AG,sep="_")
        df3.meta$v_gene_cdr3_AG_BD <- paste(df3.meta$v_gene_cdr3_AG,df3.meta$v_gene_cdr3_BD,sep=" & ")

        df3.meta2 <- df3.meta[names(df3.meta) %in% c(input$Samp_col,"v_gene_cdr3_AG_BD")]
        names(df3.meta2)[names(df3.meta2) %in% input$Samp_col] <- "ID_Column"
        names(df3.meta2)[names(df3.meta2) %in% "v_gene_cdr3_AG_BD"] <- "v_gene_selected"
        df3.meta2
      }
      else {
        df3.meta2 <- df3.meta[names(df3.meta) %in% c(input$Samp_col,input$V_gene_sc)]
        names(df3.meta2)[names(df3.meta2) %in% input$Samp_col] <- "ID_Column"
        names(df3.meta2)[names(df3.meta2) %in% input$V_gene_sc] <- "v_gene_selected"
      }
      df3.meta2
      df3.meta3 <- subset(df3.meta2,df3.meta2$v_gene_selected!="_._")
      df3.meta3 <- subset(df3.meta3,df3.meta3$v_gene_selected!="NA_NA & NA_NA")
      df3.meta3$v_gene_selected <- ifelse(df3.meta3$v_gene_selected=="",NA,df3.meta3$v_gene_selected)
      df3.meta3$v_gene_selected <- ifelse(df3.meta3$v_gene_selected=="NA",NA,df3.meta3$v_gene_selected)
      df3.meta3 <- df3.meta3[complete.cases(df3.meta3),]
      df3.meta3
      meta2.names <- names(df3.meta3)
      df3.meta3$samp.count <- 1
      total.condition <- as.data.frame(ddply(df3.meta3,"ID_Column",numcolwise(sum)))
      emtpy <- matrix(nrow =dim(df3.meta3)[1],ncol=dim(total.condition)[1])

      for (i in 1:dim(df3.meta3)[1]) {

        emtpy[i,] <- ifelse(df3.meta3$ID_Column[i]==total.condition$ID_Column[1:dim(total.condition)[1]],
                            total.condition[total.condition$ID_Column==total.condition$ID_Column[1:dim(total.condition)[1]],2],F)
      }
      as.data.frame(emtpy)

      df3.meta3$frequency <- 1/rowSums(emtpy)
      df3.meta3$percent <- 1/rowSums(emtpy)*100

      df3 <- as.data.frame(ddply(df3.meta3,meta2.names,numcolwise(sum)))
      df3 <- df3[order(df3$samp.count,decreasing = T),]

      df4 <- df3 %>%
        mutate(Clonality = case_when(
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

    ## Umap -----
    create_UMAP2 <- reactive({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      reduction <- (sc@reductions$umap)
      UMAP <- as.data.frame(reduction@cell.embeddings)
      UMAP$Cell_Index <- rownames(UMAP)
      meta.data <- as.data.frame(sc@meta.data)
      if(input$STEGO_R_pro == "STEGO_R (.h5Seurat)") {
        meta.data <- meta.data
        umap.meta <- merge(UMAP,meta.data,by="Cell_Index")
        names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
        names(umap.meta)[names(umap.meta) %in% input$V_gene_sc] <- "v_gene_selected"
      }
      else {
        meta.data$Cell_Index <- rownames(meta.data)
        meta.data$v_gene_AG <- meta.data[,names(meta.data) %in% input$RDS_V_gene_A]
        meta.data$v_gene_BD <- meta.data[,names(meta.data) %in% input$RDS_V_gene_B]
        meta.data$cdr3_AG <- meta.data[,names(meta.data) %in% input$RDS_cdr3_A]
        meta.data$cdr3_BD <- meta.data[,names(meta.data) %in% input$RDS_cdr3_B]
        meta.data$v_gene_cdr3_BD <- paste(meta.data$v_gene_BD,meta.data$cdr3_BD,sep="_")
        meta.data$v_gene_cdr3_AG<- paste(meta.data$v_gene_AG,meta.data$cdr3_AG,sep="_")
        meta.data$v_gene_cdr3_AG_BD <- paste(meta.data$v_gene_cdr3_AG,meta.data$v_gene_cdr3_BD,sep=" & ")
        umap.meta <- merge(UMAP,meta.data,by="Cell_Index")
        names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
        names(umap.meta)[names(umap.meta) %in% "v_gene_cdr3_AG_BD"] <- "v_gene_selected"
      }

      sc <- merge(umap.meta,TCR_Expansion(),by=c("v_gene_selected","ID_Column"),all.x=T)
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
        x <- gsub(":", ".", Sys.time())
        paste("TCR_Explore_UMAP_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_UMAP2,height=input$height_UMAP2, onefile = FALSE) # open the pdf device
        plot(create_UMAP2())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_UMAP2 <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste("TCR_Explore_UMAP_", gsub("/", "-", x), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_UMAP2,
            height = input$height_png_UMAP2,
            res = input$resolution_PNG_UMAP2)
        plot(create_UMAP2())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    output$Tb_TCR_clonotypes.table <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      UMAP.wt.clonality <- UMAP.TCRclonalit()

      UMAP.wt.clonality <- subset(UMAP.wt.clonality)
      num <- as.data.frame(unique(UMAP.wt.clonality$TYPE.clonality))
      num <- as.data.frame(num[complete.cases(num)==T,])
      num
    })

    # clonotypes Plot -----

    cols_clonal_plot <- reactive({
      df4 <- TCR_Expansion()
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
      df4 <- TCR_Expansion()
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

      df4 <- TCR_Expansion()
      names(df4)[names(df4) %in% input$Samp_col] <- "ID_Column"
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
      df4 <- TCR_Expansion()
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
      df4 <- TCR_Expansion()
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
      df4 <- TCR_Expansion()
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
      if (input$Graph_type_bar=="Clonality") {
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
      if (input$Graph_type_bar=="Clonality") {
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
        x <- gsub(":", ".", Sys.time())
        paste("TCR_Explore_motif_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_clonality.bar.graph,height=input$height_clonality.bar.graph, onefile = FALSE) # open the pdf device
        if (input$Graph_type_bar=="Clonality") {
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
        x <- gsub(":", ".", Sys.time())
        paste("TCR_Explore_motif_", gsub("/", "-", x), ".png", sep = "")
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
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      reduction <- (sc@reductions$umap)
      UMAP <- as.data.frame(reduction@cell.embeddings)
      UMAP$Cell_Index <- rownames(UMAP)
      meta.data <- as.data.frame(sc@meta.data)
      if(input$STEGO_R_pro == "STEGO_R (.h5Seurat)") {
        meta.data <- meta.data
        umap.meta <- merge(UMAP,meta.data,by="Cell_Index")
        head(umap.meta)

        names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
        names(umap.meta)[names(umap.meta) %in% input$V_gene_sc] <- "v_gene_selected"
        UMAP.wt.clonality <- merge(umap.meta,TCR_Expansion(),by=c("v_gene_selected","ID_Column"),all.x=T)
        UMAP.wt.clonality$Chain <- ifelse(grepl("_._", UMAP.wt.clonality$v_gene_selected),NA,
                                          ifelse(UMAP.wt.clonality$v_gene_selected=="_",NA,
                                                 ifelse(UMAP.wt.clonality$v_gene_selected=="",NA,
                                                        ifelse(UMAP.wt.clonality$v_gene_selected==" & ",NA,
                                                               ifelse(grepl("TRBV",UMAP.wt.clonality$v_gene_selected) |
                                                                        grepl("TRAV",UMAP.wt.clonality$v_gene_selected),"abTCR",
                                                                      ifelse(grepl("TRGV",UMAP.wt.clonality$v_gene_selected) |
                                                                               grepl("TRDV",UMAP.wt.clonality$v_gene_selected),"gdTCR","Other"))))))

      }
      else {
        meta.data$Cell_Index <- rownames(meta.data)
        meta.data$v_gene_AG <- meta.data[,names(meta.data) %in% input$RDS_V_gene_A]
        meta.data$v_gene_BD <- meta.data[,names(meta.data) %in% input$RDS_V_gene_B]
        meta.data$cdr3_AG <- meta.data[,names(meta.data) %in% input$RDS_cdr3_A]
        meta.data$cdr3_BD <- meta.data[,names(meta.data) %in% input$RDS_cdr3_B]
        meta.data$v_gene_cdr3_BD <- paste(meta.data$v_gene_BD,meta.data$cdr3_BD,sep="_")
        meta.data$v_gene_cdr3_AG<- paste(meta.data$v_gene_AG,meta.data$cdr3_AG,sep="_")
        meta.data$v_gene_cdr3_AG_BD <- paste(meta.data$v_gene_cdr3_AG,meta.data$v_gene_cdr3_BD,sep=" & ")
        umap.meta <- merge(UMAP,meta.data,by="Cell_Index")
        head(umap.meta)

        names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
        names(umap.meta)[names(umap.meta) %in% "v_gene_cdr3_AG_BD"] <- "v_gene_selected"

        UMAP.wt.clonality <- merge(umap.meta,TCR_Expansion(),by=c("v_gene_selected","ID_Column"),all.x=T)
        UMAP.wt.clonality$Chain <- ifelse(grepl("_._", UMAP.wt.clonality$v_gene_selected),NA,
                                          ifelse(UMAP.wt.clonality$v_gene_selected=="_",NA,
                                                 ifelse(UMAP.wt.clonality$v_gene_selected=="",NA,
                                                        ifelse(UMAP.wt.clonality$v_gene_selected==" & ",NA,
                                                               ifelse(grepl("TRBV",UMAP.wt.clonality$v_gene_selected) |
                                                                        grepl("TRAV",UMAP.wt.clonality$v_gene_selected),"abTCR",
                                                                      ifelse(grepl("TRGV",UMAP.wt.clonality$v_gene_selected) |
                                                                               grepl("TRDV",UMAP.wt.clonality$v_gene_selected),"gdTCR","Other"))))))

      }
      # umap.meta <- merge(UMAP,meta.data,by="Cell_Index")
      # head(umap.meta)
      #
      # names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
      # names(umap.meta)[names(umap.meta) %in% input$V_gene_sc] <- "v_gene_selected"




      if (input$Graph_type_bar== "Number_expanded") {
        UMAP.wt.clonality$TYPE.clonality <- paste(UMAP.wt.clonality$Chain,UMAP.wt.clonality$Number_expanded)
      }
      else {
        UMAP.wt.clonality$TYPE.clonality <- paste(UMAP.wt.clonality$Chain,UMAP.wt.clonality$Clonality)

      }

      UMAP.wt.clonality$TYPE.clonality<- ifelse(grepl("NA", UMAP.wt.clonality$TYPE.clonality),NA,UMAP.wt.clonality$TYPE.clonality)

      UMAP.wt.clonality

    })

    UMAP.TCRclonalit2 <- reactive({
      UMAP.wt.clonality <- UMAP.TCRclonalit()
      UMAP.wt.clonality$TYPE.clonality <- factor(UMAP.wt.clonality$TYPE.clonality,levels = unique(UMAP.wt.clonality$TYPE.clonality))
      # UMAP.wt.clonality$TYPE.clonality <- factor(UMAP.wt.clonality$TYPE.clonality,levels =unique(UMAP.wt.clonality$TYPE.clonality))
      # list <- unique(UMAP.wt.clonality$TYPE.clonality)[-is.na(unique(UMAP.wt.clonality$TYPE.clonality))]
      # UMAP.wt.clonality$TYPE.clonality <- factor(UMAP.wt.clonality$TYPE.clonality,levels = unique(UMAP.wt.clonality$TYPE.clonality), exclude = NA, ordered = T)

      colorblind_vector <- unlist(colors_UMAP_clonal_plot())
      names(UMAP.wt.clonality)[names(UMAP.wt.clonality) %in% input$Samp_col] <- "ID_Column"
      UMAP.wt.clonality$ID_Column <- factor(UMAP.wt.clonality$ID_Column,levels = input$ID_Column_factor)

      plot <- ggplot(UMAP.wt.clonality,aes(x=UMAP_1,UMAP_2,colour=TYPE.clonality,alpha = TYPE.clonality,label=TYPE.clonality))+
        geom_point()+
        scale_color_manual(values = colorblind_vector, na.value=input$NA_col_analysis,labels = ~ stringr::str_wrap(.x, width = 40))+
        scale_alpha_manual(values = rep(1,length(unique(UMAP.wt.clonality$TYPE.clonality))), na.value=0.5,labels = ~ stringr::str_wrap(.x, width = 40))+
        theme_bw() +
        theme(
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
        # UMAP.wt.clonality2$ID.names <- UMAP.wt.clonality2[,names(UMAP.wt.clonality2) %in% input$meta_data_sc_]
        # UMAP.wt.clonality2$ID.names <- UMAP.wt.clonality2[,names(UMAP.wt.clonality2) %in% input$Samp_col_UMAP]

        plot+ facet_wrap(~ID_Column)
      }

    })
    output$clonality.TCR.UMAP <- renderPlot({
      UMAP.TCRclonalit2()
    })

    output$downloadPlot_TCR.UMAP <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste("TCR.UMAP_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_TCR.UMAP,height=input$height_TCR.UMAP, onefile = FALSE) # open the pdf device

        plot(UMAP.TCRclonalit2())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_TCR.UMAP <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste("TCR.UMAP_", gsub("/", "-", x), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_TCR.UMAP,height = input$height_png_TCR.UMAP,res = input$resolution_PNG_TCR.UMAP)
        plot(UMAP.TCRclonalit2())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )
    # clonality ------
    observe({
      sc <- input.data_sc_pro()
      # sc$ID_Column <-"ID_Column"
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
      df4 <- TCR_Expansion()
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
      UMAP.wt.clonality2 <- merge(UMAP.wt.clonality,top10.2,by=c("v_gene_selected","ID_Column"),all.x=T)
      UMAP.wt.clonality2 <- UMAP.wt.clonality2[order(UMAP.wt.clonality2$Cell_Index),]
      UMAP.wt.clonality2 %>% mutate(across(where(is.factor), as.character)) -> UMAP.wt.clonality2
      UMAP.wt.clonality2$topclones[is.na(UMAP.wt.clonality2$topclones)]<- "not selected"
      UMAP.wt.clonality2$topclones <- factor(UMAP.wt.clonality2$topclones,levels = c("not selected",unique(as.character(top10.2$v_gene_selected))))
      # UMAP.wt.clonality2 <- UMAP.wt.clonality2[order(UMAP.wt.clonality2$Cell_Index),]
      UMAP.wt.clonality2 <- UMAP.wt.clonality2[order(UMAP.wt.clonality2$topclones),]
      UMAP.wt.clonality2

    })

    cols_UMAP_Topclonotypes <- reactive({
      UMAP.wt.clonality2 <- For_col_top()

      num <- as.data.frame(unique(UMAP.wt.clonality2$topclones))

      col.gg <- c("grey",gg_fill_hue(dim(num)[1]))
      palette_rainbow <- c("grey",rainbow(dim(num)[1]))
      heat_col <- c("grey",heat.colors(dim(num)[1]))
      col.terrain <- c("grey",terrain.colors(dim(num)[1]))
      col.topo <- c("grey",topo.colors(dim(num)[1]))
      col.hcl <- c("grey",hcl.colors(dim(num)[1], palette = "viridis"))


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

      col.gg <- c("grey",gg_fill_hue(dim(num)[1]))
      palette_rainbow <- c("grey",rainbow(dim(num)[1]))
      heat_col <- c("grey",heat.colors(dim(num)[1]))
      col.terrain <- c("grey",terrain.colors(dim(num)[1]))
      col.topo <- c("grey",topo.colors(dim(num)[1]))
      col.hcl <- c("grey",hcl.colors(dim(num)[1], palette = "viridis"))

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.UMAP_clonotype_top", i, sep="_")]]
      })
    })

    select_group_metadata <- function () {
      sc <- input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      df <- as.data.frame(sc@meta.data)
      df2 <- as.data.frame(unique(df[names(df) %in% input$Samp_col]))
      df2 <- as.data.frame(df2)
      df2

    }
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
        choices=df2
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
        "ID_Column_factor",
        choices=df2,
        selected = df2
      )


    })

    Topclonotypes <- reactive({
      UMAP.wt.clonality2 <- For_col_top()
      UMAP.wt.clonality2$ID_Column <- factor(UMAP.wt.clonality2$ID_Column, levels = input$ID_Column_factor)

      colorblind_vector <-as.data.frame(unlist(colors_UMAP_Topclonotypes()))
      names(colorblind_vector) <- "col"
      topclones_col <- as.data.frame(unique(UMAP.wt.clonality2$topclones))
      names(topclones_col) <- "topclones"
      topclones_col$col <- colorblind_vector$col
      topclones_col


      if (input$display_all_samps=="yes" & input$Split_by_group=="no") {
        topclones_col
      }

      else if (input$display_all_samps=="yes" & input$Split_by_group=="yes") {
        topclones_col
      }

      else if (input$display_all_samps=="no" & input$Split_by_group=="yes") {
        UMAP.wt.clonality2 <- UMAP.wt.clonality2[UMAP.wt.clonality2$ID_Column %in% input$ID_Column_metadata,]
        topclones_col <- topclones_col[topclones_col$topclones %in% unique(UMAP.wt.clonality2$topclones),]
        topclones_col
      }

      else {
        UMAP.wt.clonality2 <- UMAP.wt.clonality2[UMAP.wt.clonality2$ID_Column %in% input$ID_Column_metadata,]
        topclones_col <- topclones_col[topclones_col$topclones %in% unique(UMAP.wt.clonality2$topclones),]
        topclones_col
      }

      UMAP.wt.clonality2$topclones2 <- gsub(" & "," ",UMAP.wt.clonality2$topclones)

      names(UMAP.wt.clonality2)[names(UMAP.wt.clonality2) %in% input$Samp_col] <- "ID_Column"
      UMAP.wt.clonality2$ID_Column <- factor(UMAP.wt.clonality2$ID_Column,levels = input$ID_Column_factor)

      plot <- ggplot(data=UMAP.wt.clonality2,aes(x=UMAP_1,UMAP_2,colour=topclones2,label =topclones2 ))+
        geom_point()+
        scale_color_manual(values = topclones_col$col, na.value=input$NA_col_analysis,labels = ~ stringr::str_wrap(.x, width = 50))+
        # scale_color_manual(values = colorblind_vector)+
        scale_alpha_manual(values = 1, na.value = 0.1,labels = ~ stringr::str_wrap(.x, width = 50))+
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

      if (input$Split_by_group=="no") {
        plot
      }
      else {
        # UMAP.wt.clonality2$ID.names <- UMAP.wt.clonality2[,names(UMAP.wt.clonality2) %in% input$meta_data_sc_]
        # UMAP.wt.clonality2$ID.names <- UMAP.wt.clonality2[,names(UMAP.wt.clonality2) %in% input$Samp_col_UMAP]

        plot+ facet_wrap(~ID_Column)
      }
    })

    output$clonality.TCR.UMAP.top <- renderPlot({
      Topclonotypes()
    })

    output$downloadPlot_TCR.UMAP_top <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste("TCR_Explore_TCR.UMAP_top_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_TCR.UMAP_top,height=input$height_TCR.UMAP_top, onefile = FALSE) # open the pdf device

        plot(Topclonotypes())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_TCR.UMAP_top <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste("TCR_Explore_TCR.UMAP_top_", gsub("/", "-", x), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_TCR.UMAP_top,height = input$height_png_TCR.UMAP_top,res = input$resolution_PNG_TCR.UMAP_top)
        plot(Topclonotypes())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    ## differential expression -----
    observe({
      df.test <- input.data_sc_pro()
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
      df.test <- input.data_sc_pro()
      validate(
        need(nrow(df.test)>0,
             error_message_val_sc)
      )
      updateSelectInput(
        session,
        "meta_data_sc_",
        choices=names(df.test@meta.data),
        selected = c("Sample_Name"))
    })
    df_sc_clust <- reactive({
      sc <- input.data_sc_pro()
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
        paste(input$name.file_clust,"DEx_all_",gsub("-", ".", Sys.Date()),".csv", sep = "")
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
      df.test <- input.data_sc_pro()
      validate(
        need(nrow(df.test)>0,
             error_message_val_sc)
      )
      name.df <- as.data.frame(unique(df.test@assays$RNA@var.features))

      names(name.df) <- "Gene_Name"
      name.df <- as.data.frame(name.df[order(name.df$Gene_Name),])
      names(name.df) <- "Gene_Name"
      updateSelectInput(
        session,
        "string.data3",
        choices=name.df$Gene_Name,
        selected = c("GZMB"))

    })
    #
    markers_featurePlot <- reactive({
      FeaturePlot(df_sc_clust(), features = c(input$string.data3), min.cutoff = 'q10', label = T)
    })
    output$markers_featurePlot_sc <- renderPlot({
      markers_featurePlot()
    })

    # download markers_featurePlot_sc
    output$downloadPlot_markers_featurePlot_sc <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"TCR_Explore_markers_featurePlot_sc_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_markers_featurePlot_sc,height=input$height_markers_featurePlot_sc, onefile = FALSE) # open the pdf device

        plot(markers_featurePlot())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_markers_featurePlot_sc <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"TCR_Explore_markers_featurePlot_sc_", gsub("/", "-", x), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_markers_featurePlot_sc,height = input$height_png_markers_featurePlot_sc,res = input$resolution_PNG_markers_featurePlot_sc)
        plot(markers_featurePlot())
        dev.off()},   contentType = "application/png")


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
        paste(input$name.file_clust," Differential expresion of ",input$unique.Idents1, " vs. ", input$unique.Idents2,gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(DEx_sc())
        write.csv(df,file, row.names = T)
      } )


    output$downloaddf_DEx_sc_ggVolcanoR <- downloadHandler(
      filename = function(){
        paste(input$name.file_clust," Differential expresion of ",input$unique.Idents1, " vs. ", input$unique.Idents2,gsub("-", ".", Sys.Date()),".csv", sep = "")
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
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )

      umap.meta.classification <-  sc@meta.data
      umap.meta.classification$Motif_gene <- umap.meta.classification[,names(umap.meta.classification) %in% input$V_gene_sc]
      umap.meta.classification$Chain <- ifelse(grepl("TRBV",umap.meta.classification$Motif_gene) |
                                                 grepl("TRAV",umap.meta.classification$Motif_gene),"abTCR",
                                               ifelse(grepl("TRGV",umap.meta.classification$Motif_gene) |
                                                        grepl("TRDV",umap.meta.classification$Motif_gene),"gdTCR","unknown"))

      umap.meta.classification$umap.meta.classification.chain <- paste(umap.meta.classification$Chain,umap.meta.classification$classify.T.cell)
      #
      umap.meta.classification$umap.meta.classification.chain <- ifelse(grepl("unknown", umap.meta.classification$umap.meta.classification.chain),"unknown",umap.meta.classification$umap.meta.classification.chain)
      umap.meta.classification$umap.meta.classification.chain <- ifelse(grepl(" NA", umap.meta.classification$umap.meta.classification.chain),"",umap.meta.classification$umap.meta.classification.chain)

      umap.meta.classification$Selected_function <- umap.meta.classification[,names(umap.meta.classification) %in% input$clust_names_top]
      umap.meta.classification$cloneCount <- 1
      top_BD_cluster2 <-  umap.meta.classification[,names(umap.meta.classification) %in% c("cloneCount","Selected_function")]
      top_BD_cluster2 <- top_BD_cluster2 %>%
        select(cloneCount, everything())
      df2 <- as.data.frame(ddply(top_BD_cluster2,names(top_BD_cluster2)[-c(1)],numcolwise(sum)))
      df2$fraction <- df2$cloneCount/sum(df2$cloneCount)
      df2$Percent <- round(df2$cloneCount/sum(df2$cloneCount)*100,2)
      df2$Selected_function <- ifelse(grepl("NA", df2$Selected_function),NA,df2$Selected_function)

      df2[order(df2$Percent,decreasing = T),]
    })

    output$Percent_tab <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{

      Percent_tab_df()
    })
    output$downloaddf_Percent_tab <- downloadHandler(
      filename = function(){
        paste(input$name.file_clust,"_Percent_",input$clust_names_top,"_",gsub("-", ".", Sys.Date()),".tsv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(Percent_tab_df())
        write.table(df,file, row.names = F,sep="\t", quote = F)
      } )

    Add.UMAP.reduction <- reactive ({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )
      reduction <- (sc@reductions$umap)
      UMAP <- as.data.frame(reduction@cell.embeddings)
      UMAP$Cell_Index <- rownames(UMAP)
      meta.data <- as.data.frame(sc@meta.data)
      if(input$STEGO_R_pro == "STEGO_R (.h5Seurat)") {
        meta.data <- meta.data
      }
      else {
        meta.data$Cell_Index <- rownames(meta.data)
      }

      meta.data$Motif_gene <- meta.data[,names(meta.data) %in% input$V_gene_sc]

      meta.data$Chain <- ifelse(grepl("TRBV",meta.data$Motif_gene) |
                                  grepl("TRAV",meta.data$Motif_gene),"abTCR",
                                ifelse(grepl("TRGV",meta.data$Motif_gene) |
                                         grepl("TRDV",meta.data$Motif_gene),"gdTCR","unknown"))

      meta.data$umap.meta.classification.chain <- paste(meta.data$Chain,meta.data$classify.T.cell)
      #
      meta.data$umap.meta.classification.chain <- ifelse(grepl("unknown", meta.data$umap.meta.classification.chain),"unknown",meta.data$umap.meta.classification.chain)
      meta.data$umap.meta.classification.chain <- ifelse(grepl(" NA", meta.data$umap.meta.classification.chain),"",meta.data$umap.meta.classification.chain)

      meta.data <- meta.data[order(as.numeric(meta.data$Cell_Index),decreasing = F),]

      umap.meta <- merge(UMAP,meta.data,by="Cell_Index")

      umap.meta


    }) # Add.UMAP.reduction

    output$SiteNumInput <- renderUI({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP))

      df_class <- sc@meta.data
      df_class$selected <- df_class[,names(df_class) %in% input$clust_names_top]
      selectInput("SiteNumInput", "Show on graph", choices = unique(df_class$selected), selected = NULL, multiple = TRUE)
    })

    # Umap classification plot -------
    cols_UMAP_all_classification <- reactive({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )

      top_BD_cluster <-  Add.UMAP.reduction()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
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
      top_BD_cluster <-  Add.UMAP.reduction()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
      top_BD_cluster$Selected_function <- ifelse(grepl("NA", top_BD_cluster$Selected_function),NA,top_BD_cluster$Selected_function)
      top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num)==T,])

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.UMAP_all_classification", i, sep="_")]]
      })
    })

    UMAP_all_classification <- reactive({
      top_BD_cluster <-  Add.UMAP.reduction()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
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
      # top_BD_cluster$alpha_val <- ifelse(is.na(top_BD_cluster$Selected_function)==T,0.1,1)

      names(top_BD_cluster)[names(top_BD_cluster) %in% input$Samp_col] <- "ID_Column"
      top_BD_cluster$ID_Column <- factor(top_BD_cluster$ID_Column,levels = input$ID_Column_factor)

      df <- ggplot(top_BD_cluster,aes(x=UMAP_1,UMAP_2,colour=Selected_function,alpha = Selected_function,label = Selected_function))+
        geom_point()+
        scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 20),values = col.file$col, na.value=input$NA_col_analysis)+
        scale_alpha_manual(labels = ~ stringr::str_wrap(.x, width = 20),values = rep(1,length(unique(top_BD_cluster$Selected_function))), na.value=0.1) +
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

      if (input$Split_by_group=="no") {
        df <- df
      }
      else {
        df <- df + facet_wrap(~ID_Column)
      }
      df
    })

    output$UMAP_all_classification2 <- renderPlot({
      UMAP_all_classification()
    })

    output$downloadPlot_UMAP_all_classification  <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"top_clonotype_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_UMAP_all_classification ,height=input$height_UMAP_all_classification , onefile = FALSE) # open the pdf device
        plot(UMAP_all_classification())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_UMAP_all_classification  <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"top_clonotype_", gsub("/", "-", x), ".png", sep = "")
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
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      meta.data <- as.data.frame(sc@meta.data)
      updateSelectInput(
        session,
        "clust_group",
        choices=names(meta.data),
        selected = "seurat_clusters")
    })
    # pie colouring  ----
    cols_pie <- reactive({
      top_BD_cluster <-  Add.UMAP.reduction()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
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
      top_BD_cluster <-  Add.UMAP.reduction()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num)==T,])

      lapply(1:dim(num)[1], function(i) {
        input[[paste("col.pie", i, sep="_")]]
      })
    })

    Pie_chart_Class <- reactive({
      top_BD_cluster <-  Add.UMAP.reduction()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
      top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_group]
      df.col <- unlist(colors_pie())

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
        facet_wrap(~Var1) +
        theme(
          legend.key.size = unit(1, 'cm'))+
        scale_fill_manual(values = df.col, na.value = input$NA_col_analysis) +
        theme(strip.text = element_text(size = input$strip_size, family = input$font_type),
              legend.text = element_text(size = input$Bar_legend_size, family = input$font_type),
              legend.title = element_blank()

        )
      # theme(strip.background =element_rect(fill=input$strip.colour.tree))+
      # theme(strip.text = element_text(colour = input$strip.text.colour.tree))
      # scale_pattern_fill_manual(values = class_col.t)
    })

    output$Classification_clonotype_pie <- renderPlot({
      Pie_chart_Class()
    })
    output$downloadPlot_Classification_clonotype_pie <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"_Classification_clonotype_pie_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Classification_clonotype_pie,height=input$height_Classification_clonotype_pie, onefile = FALSE) # open the pdf device
        plot(Pie_chart_Class())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Classification_clonotype_pie <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"_Classification_clonotype_pie_", gsub("/", "-", x), ".png", sep = "")
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
      meta.data <-  Add.UMAP.reduction()
      totals <- meta.data[,names(meta.data) %in% c(input$Samp_col,input$clust_names_top)]
      names(totals) <- c("groups","Function")
      totals <- totals[!totals$Function %in% c("NA"),]
      totals <- totals[!totals$groups %in% c("NA"),]
      totals$groups <- factor(totals$groups,levels = input$ID_Column_factor)
      totals
    })

    output$Chi_tab_before <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      chi_squ()
    })
    # Chi_square <- reactiveValues(plot_chi=NULL)

    chi.seq.residuals <- reactive({

    })



    Chi_square_plot2 <- reactive({
      totals <- chi_squ()
      tb_totals <- table(totals$groups,totals$Function)
      chisq <- chisq.test(tb_totals)
      chisq
      if (input$type_res=="Residuals") {
        res <- chisq$residuals
        res <- setNames(melt(res), c('x', 'y', 'residuals'))
        res <- res[res$x %in% input$ID_Column_factor,]

        plot_chi <- ggplot(res, aes(x = x, y = y,size = residuals, fill = residuals)) +
          geom_point(shape = 21, colour = "black")+
          scale_size_area(max_size = 10) +
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

      }
      else {
        contrib <- 100*chisq$residuals^2/chisq$statistic
        contrib <- setNames(melt(contrib), c('x', 'y', 'Percentage'))
        contrib <- contrib[contrib$x %in% input$ID_Column_factor,]
        plot_chi <- ggplot(contrib, aes(x = x, y = y,size = Percentage, fill = Percentage)) +
          geom_point(shape = 21, colour = "black")+
          scale_size_area(max_size = 10) +
          scale_fill_gradient2(
            # low = input$lower_col_chi2,
            mid = input$lower_col_chi2,
            high = input$high_col_chi2,
            space = "Lab",
            na.value = "grey90",
            guide = "colourbar",
            aesthetics = "fill" )+
          theme_bw()+
          theme(
            axis.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
            axis.text.x = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type,angle = 90),
            legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
            # legend.title = element_blank(),
            axis.line.x = element_blank(),
            panel.grid = element_blank()
            # legend.position = input$legend_position,
          )+
          guides(size = "none")
      }
      plot_chi
    })

    output$Chi_square_plot <- renderPlot({
      Chi_square_plot2()
    })


    output$downloadPlot_Chi_square_plot <- downloadHandler(
      filename = function() {
        paste(input$name.file_clust,"_Chi_square_plot", ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Chi_square_plot,height=input$height_Chi_square_plot, onefile = FALSE) # open the pdf device
        plot(Chi_square_plot2())
        dev.off()
      }, contentType = "application/pdf")

    output$downloadPlotPNG_Chi_square_plot <- downloadHandler(
      filename = function() {
        paste(input$name.file_clust,"_Chi_square_plot_",".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Chi_square_plot,height = input$height_png_Chi_square_plot,res = input$resolution_PNG_Chi_square_plot)
        plot(Chi_square_plot2())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )


    # Select top clonotype ----
    observe({
      sc <- input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val1)
      )
      meta.data <- as.data.frame(sc@meta.data)
      # meta.data <- meta.data[,names(meta.data) %in% c(grep("v",names(meta.data)), grep("cdr3",names(meta.data)))]

      if (input$datasource == "BD rhapsody") {
        updateSelectInput(
          session,
          "Gene_V_top",
          choices=names(meta.data),
          selected = "v_gene_cdr3_AB_GD")
      }

      else {
        updateSelectInput(
          session,
          "Gene_V_top",
          choices=names(meta.data),
          selected = "vdj_gene_cdr3_AG_BD")
      }

      # updateSelectInput(
      #  session,
      #  "",
      #  choices=names(meta.data),
      #  selected = "v_gene_cdr3_AG")

    })
    Top_clonotype_df2 <- reactive({
      df3.meta <- Add.UMAP.reduction()
      # summarising the clonality
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$Gene_V_top]
      df3.meta$cloneCount <- 1

      BD <- df3.meta[,names(df3.meta) %in% c("cluster_name","cloneCount")]
      BD
      BD_sum <- ddply(BD,names(BD)[-c(2)] ,numcolwise(sum))
      as.data.frame(BD_sum)
      BD_sum <- BD_sum[!BD_sum$cluster_name=="_",]
      names(BD_sum)[2] <- "Total_count"

      BD_sum$frequency <- BD_sum$Total_count/sum(BD_sum$Total_count)
      BD_sum<- BD_sum[order(BD_sum$frequency,decreasing = T),]
      BD_sum2 <- BD_sum[!BD_sum$cluster_name %in% "NA",]
      BD_sum2
    })
    observeEvent(input$run_top,{
      BD_sum <- Top_clonotype_df2()
      validate(
        need(nrow(BD_sum)>0,
             error_message_val1)
      )


      updateSelectInput(
        session,
        "Selected_clonotype",
        choices=BD_sum$cluster_name[1:50],
        selected = BD_sum$cluster_name[1])
    })
    observeEvent(input$run_top,{
      BD_sum <- Add.UMAP.reduction()
      validate(
        need(nrow(BD_sum)>0,
             error_message_val1)
      )

      updateSelectInput(
        session,
        "Selected_chain",
        choices=names(BD_sum),
        selected = "Motif_gene")
    })
    output$Top_clonotype_df<- DT::renderDataTable({
      Top_clonotype_df2()

      df3.meta <- Add.UMAP.reduction()
      # summarising the clonality
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$Gene_V_top]
      top_BD_clonotype <- df3.meta[df3.meta$cluster_name %in% input$Selected_clonotype,]
      top_BD_clonotype
      datatable(top_BD_clonotype, extensions = "Buttons", options = list(searching = TRUE,
                                                                         ordering = TRUE,
                                                                         buttons = c('copy','csv', 'excel'),
                                                                         dom = 'Bfrtip',
                                                                         pageLength=10,
                                                                         lengthMenu=c(2,5,10,20,50,100),
                                                                         scrollX = TRUE
      ))
    }, server = FALSE)

    # top clonptype bar graph -----
    cols_Top_bar_clonotype <- reactive({
      dtop_clonotype_bar_code <- top_clonotype_bar_code()
      dtop_clonotype_bar_code$Selected_chain <- dtop_clonotype_bar_code[,names(dtop_clonotype_bar_code) %in% input$Selected_chain]
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
      dtop_clonotype_bar_code$Selected_chain <- dtop_clonotype_bar_code[,names(dtop_clonotype_bar_code) %in% input$Selected_chain]
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
    top_clonotype_bar_code <- reactive({
      df3.meta <- Add.UMAP.reduction()
      # summarising the clonality
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$Gene_V_top]

      top_BD_clonotype <- df3.meta[df3.meta$cluster_name %in% input$Selected_clonotype,]
      top_BD_clonotype

    })
    ggplot_top_BD_clonotype_vs_SC <- reactive({
      dtop_clonotype_bar_code <- top_clonotype_bar_code()
      dtop_clonotype_bar_code$Selected_group <- dtop_clonotype_bar_code[,names(dtop_clonotype_bar_code) %in% input$clust_group]
      colorblind_vector <- unlist(colors_cols_Top_bar_clonotype())

      dtop_clonotype_bar_code$Selected_chain2 <- dtop_clonotype_bar_code[,names(dtop_clonotype_bar_code) %in% input$Selected_chain]

      dtop_clonotype_bar_code$Selected_chain3 <- gsub("_"," ",dtop_clonotype_bar_code$Selected_chain2)
      dtop_clonotype_bar_code$Selected_chain3 <- gsub("[.]"," ",dtop_clonotype_bar_code$Selected_chain3)

      ggplot(dtop_clonotype_bar_code, aes(x=Selected_group, fill=Selected_chain3,colour = Selected_chain3, label = Selected_chain3)) +
        geom_bar() +
        theme_bw()+
        scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector, na.value=input$NA_col_analysis)+
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector, na.value=input$NA_col_analysis)+
        # scale_alpha_manual(values = rep(1,length(unique(dtop_clonotype_bar_code$Selected_chain))), na.value=0.5)+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        theme(
          axis.title.y = element_text(colour="black",family=input$font_type,size = input$title.text.sizer2),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=0),
          axis.title.x = element_text(colour="black",angle=0,vjust=.5,face="plain",family=input$font_type,size = input$title.text.sizer2),
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
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"top_clonotype_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_top_clonotype,height=input$height_top_clonotype, onefile = FALSE) # open the pdf device
        plot(ggplot_top_BD_clonotype_vs_SC())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_top_clonotype <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"top_clonotype_", gsub("/", "-", x), ".png", sep = "")
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
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
      top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_group]
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
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
      top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_group]
      # top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(unique(top_BD_cluster$Selected_function))
      num <- as.data.frame(num[complete.cases(num)==T,])

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
    Pie_chart_alpha_gamma <- reactive({
      top_BD_cluster <-  top_clonotype_bar_code()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
      top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_group]
      df.col <- unlist(colors_cols_Top_pie_clonotype())

      df3.meta3 <-  as.data.frame(table(top_BD_cluster$Selected_group,top_BD_cluster$Selected_function))
      total.condition <- as.data.frame(ddply(df3.meta3,"Var1",numcolwise(sum)))
      dim(total.condition)[1]
      dim(df3.meta3)[1]
      emtpy <- matrix(nrow =dim(df3.meta3)[1],ncol=dim(total.condition)[1])
      length(empty)

      for (i in 1:dim(df3.meta3)[1]) {

        emtpy[i,] <- ifelse(df3.meta3$Var1[i]==total.condition$Var1[1:dim(total.condition)[1]],
                            total.condition[total.condition$Var1==total.condition$Var1[1:dim(total.condition)[1]],2],F)
      }
      as.data.frame(emtpy)

      df3.meta3$n <- df3.meta3$Freq/rowSums(emtpy)

      ggplot(df3.meta3,aes(x="", y=n, fill=Var2, group = Var1)) +
        geom_bar(stat="identity", width=1)+
        coord_polar("y", start=0)  +
        theme_void(20) +
        facet_wrap(~Var1) +
        theme(
          legend.key.size = unit(1, 'cm'),
          legend.title = element_blank()) +
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 20),values = df.col, na.value = input$NA_col_analysis) +
        theme(strip.text = element_text(size = input$strip_size_pie, family = input$font_type),
              legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.position = input$legend_position,
              legend.title = element_blank()
        )
      # theme(strip.background =element_rect(fill=input$strip.colour.tree))+
      # theme(strip.text = element_text(colour = input$strip.text.colour.tree))
      # scale_pattern_fill_manual(values = class_col.t)

    })

    output$top_clonotype_pie <- renderPlot({
      if (input$Plot_type_selected=="pie") {
        Pie_chart_alpha_gamma()
      }
      else {

        if (input$Split_by_group=="no") {
          UMAP_chart_alpha_gamma()
        }
        else {

          UMAP_chart_alpha_gamma() + facet_wrap(~get(input$Samp_col))
        }

      }
    })

    output$downloadPlot_top_clonotype_pie <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"top_clonotype_pie_",gsub("/", "-", x), ".pdf", sep = "")
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
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"top_clonotype_", gsub("/", "-", x), ".png", sep = "")
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

    #### top  UMAP -----
    UMAP_chart_alpha_gamma <- reactive({

      df3.meta <- Add.UMAP.reduction()

      top_BD_cluster <-  top_clonotype_bar_code()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
      top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
      top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_group]
      df.col <- unlist(colors_cols_Top_pie_clonotype())
      df3.meta$Selected_function <- NA

      names(top_BD_cluster)[names(top_BD_cluster) %in% input$Samp_col] <- "ID_Column"
      top_BD_cluster$ID_Column <- factor(top_BD_cluster$ID_Column,levels = input$ID_Column_factor)

      ggplot()+
        geom_point(data = df3.meta,aes(x=UMAP_1,UMAP_2,colour=Selected_function))+
        geom_point(data = top_BD_cluster,aes(x=UMAP_1,UMAP_2,colour=Selected_function))+
        scale_color_manual(values = df.col, na.value=input$NA_col_analysis) +
        # scale_size_manual(values = rep(input$size_selected_top,dim(top_BD_cluster)[1]),na.value = 1) +
        theme_bw()+
        theme(
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
        )+
        scale_size(guide = 'none')

    })

    ### Add in the TCR clustering to the seurat object ----

    #### Ridge plot (expression per T cell) ----
    vals_Ridge_top <- reactiveValues(output_stats=NULL)

    compare.stat <- reactive({
      sc <-input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )
      # df = as.data.frame(sc[["RNA"]]@scale.data)
      df = as.data.frame(sc[["RNA"]]@counts)
      md <- sc@meta.data
      MainTcell <- as.data.frame(t(df))
      MainTcell$Cell_Index <- rownames(MainTcell)
      # names(MainTcell) <- gsub("CD8B.","CD8B",names(MainTcell))

      MainTcell$Cell_Index == md$Cell_Index
      unique(md$v_allele_AG)

      md$Gene_select <- md[,names(md) %in% input$Gene_V_top]

      MainTcell$selected <- ifelse(md$Gene_select == input$Selected_clonotype,"Selected","other")

      MainTcell_selected <-  subset(MainTcell,MainTcell$selected=="Selected")
      MainTcell_other<-  subset(MainTcell,MainTcell$selected=="other")
      num <- dim(MainTcell)[2]-2
      mat_emp <- matrix(nrow=num,ncol = 6)

      for (i in 1:num) {
        mat_emp[i,1] <- round(t.test(MainTcell_other[,i], MainTcell_selected[,i])$p.value, digits = 3)
        mat_emp[i,2] <- round(t.test(MainTcell_other[,i], MainTcell_selected[,i])$conf.int[1],2)
        mat_emp[i,3] <- round(t.test(MainTcell_other[,i], MainTcell_selected[,i])$conf.int[2],2)
        mat_emp[i,4] <- round(t.test(MainTcell_other[,i], MainTcell_selected[,i])$estimate[1],2)
        mat_emp[i,5] <- round(t.test(MainTcell_other[,i], MainTcell_selected[,i])$estimate[2],2)
        mat_emp[i,6] <- wilcox.test(MainTcell_other[,i], MainTcell_selected[,i])$p.value
      }

      mat_emp <- as.data.frame(mat_emp)
      rownames(mat_emp) <- names(MainTcell)[1:num]
      names(mat_emp) <- c("Pval","lowCI","upperCI","mean of other","mean of selected","Wilcoxon")

      mat_emp_sub <- subset(mat_emp,mat_emp$Pval!="NaN")

      mat_emp_sub$bonferroni <- p.adjust(mat_emp_sub$Pval, method = "bonferroni", n = length(mat_emp_sub$Pval))
      mat_emp_sub$BH_t.test <- p.adjust(mat_emp_sub$Pval, method = "BH", n = length(mat_emp_sub$Pval))
      mat_emp_sub$BH_Wil <- p.adjust(mat_emp_sub$Wilcoxon, method = "BH", n = length(mat_emp_sub$Pval))
      mat_emp_sub[order(mat_emp_sub$BH_Wil,decreasing = F),]


    })
    output$Ridge_chart_alpha_gamma_stat_comp <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength =10, scrollX = TRUE),{

      sc <-input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )
      compare.stat()

    })

    observeEvent(input$run_string.data_Exp_top,{
      df <- compare.stat()
      validate(
        need(nrow(df)>0,
             error_message_val_sc)
      )
      df <- subset(df,df$BH_Wil<0.05)

      updateSelectInput(
        session,
        "string.data_Exp_top",
        choices=rownames(df),
        selected = rownames(df)[1])

    })

    Ridge_chart_alpha_gamma_df <- reactive({
      sc <-input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )
      df= as.data.frame(sc[["RNA"]]@scale.data)

      MainTcell <- as.data.frame(t(df))
      meta.data <- as.data.frame(sc@meta.data)
      # names(MainTcell) <- ifelse(grepl("pAbO",names(MainTcell)),names(MainTcell),toupper(names(MainTcell)))
      MainTcell$Cell_Index <- rownames(MainTcell)
      # names(MainTcell) <- gsub("CD8B.","CD8B",names(MainTcell))

      gene_df <- MainTcell[,names(MainTcell) %in% c("Cell_Index",input$string.data_Exp_top)]

      top_BD_cluster <-  top_clonotype_bar_code()
      top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_group]
      top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
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

      ggplot(df2, aes(x = get(input$string.data_Exp_top), y = get(input$clust_group), fill = get(input$clust_group))) +
        geom_density_ridges() +
        theme_ridges() +
        theme(legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              # legend.title = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = input$legend_position)+
        ggtitle(input$string.data_Exp_top)

    })



    output$Ridge_chart_alpha_gamma_stat <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength =10, scrollX = TRUE),{

      df <- Ridge_chart_alpha_gamma_df()
      df$expressed <- ifelse(df[,names(df) %in% input$string.data_Exp_top] >input$Gre_ex,"expressed","Not expressed")

      if(input$restric_ex == T) {
        df2 <- subset(df,df$expressed=="expressed")
      }
      else(
        df2 <- df
      )
      df2[is.na(df2)] <- "Unknown"
      at <- TukeyHSD(aov(get(input$string.data_Exp_top)~ get(input$clust_group),data = df2))
      tab <- as.data.frame(at[1])
      names(tab) <- c("diff" ,"lwr","upr","p.adj")
      tab$stat <- ifelse(tab$p.adj<0.0001,"****",
                         ifelse(tab$p.adj<0.001,"***",
                                ifelse(tab$p.adj<0.01,"**",
                                       ifelse(tab$p.adj<0.05,"*","NS"
                                       ))))
      tab
    })

    output$downloaddf_clusTCR_GEx <- downloadHandler(
      filename = function(){
        paste(input$name.file_clust,"Stats_",gsub("-", ".", Sys.Date()),".csv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(Ridge_chart_alpha_gamma_stat_comp())
        write.csv(df,file, row.names = T)
      } )




    ### Violin plots -----

    Violin_chart_alpha_gamma_plot <- reactive({
      df <- Ridge_chart_alpha_gamma_df()
      df$expressed <- ifelse(df[,names(df) %in% input$string.data_Exp_top] >input$Gre_ex,"expressed","Not expressed")

      if(input$restric_ex == T) {
        df2 <- subset(df,df$expressed=="expressed")
      }
      else(
        df2 <- df
      )

      ggplot(df2, aes(y = get(input$string.data_Exp_top), x = get(input$clust_group), fill = get(input$clust_group))) +
        geom_violin() +
        geom_jitter(height = 0, width = 0.1)+
        theme(legend.position = "none",

        )+
        theme_bw() +
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_text(colour="black",family=input$font_type,size = input$text_size),
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=0),
          axis.title.x = element_blank(),
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_text(colour="black",size=20,family=input$font_type),
          legend.position = "none",
        )+
        ggtitle(input$string.data_Exp_top)

    })

    # vs non-selected clonotypes -----

    Ridge_chart_alpha_gamma_df_all <- reactive({
      sc <-input.data_sc_pro()
      validate(
        need(nrow(sc)>0,
             error_message_val_sc)
      )
      df= as.data.frame(sc[["RNA"]]@scale.data)
      MainTcell <- as.data.frame(t(df))
      meta.data <- as.data.frame(sc@meta.data)
      MainTcell$Cell_Index <- rownames(MainTcell)
      gene_df <- MainTcell[,names(MainTcell) %in% c("Cell_Index",input$string.data_Exp_top)]
      df3.meta <- Add.UMAP.reduction()
      df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% input$Gene_V_top]
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
              # axis.text.x = element_blank(),
              # axis.text.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = input$legend_position)

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
          axis.text.x = element_text(colour="black",family=input$font_type,size = input$text_size,angle=0),
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
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"_",input$plot_type_ridgvi,"_",gsub("/", "-", x), ".pdf", sep = "")
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
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"_",input$epitope_hm,"_",input$pathology_hm,"_Heatmap_", gsub("/", "-", x), ".png", sep = "")
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


    #### Epitope upload -----

    df_tcrex <- reactive({
      tcrex <- input.data_sc_TCRex()
    })

    output$MainTcell_Check <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{

      df_tcrex()

    })

    output$downloaddf_TCRex_meta <- downloadHandler(
      filename = function(){
        paste(input$name.file_clust,"_","TCRex_metadata_",gsub("-", ".", Sys.Date()),".tsv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(df_tcrex())
        write.csv(df,file, row.names = F,quote = F)
      } )

    #### Epitope heatmap -----
    heatmap_epitope <- reactive({
      df3.meta <-  Add.UMAP.reduction()
      # summarising the clonality
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(df3.meta)>0 & nrow(epi)>0,
             "Upload Files")
      )
      if(input$datasource == "BD rhapsody") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

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
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"_",input$epitope_hm,"_",input$pathology_hm,"_Heatmap_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Heatmap_epi_plot ,height=input$height_Heatmap_epi_plot , onefile = FALSE) # open the pdf device
        plot(heatmap_epitope())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Heatmap_epi_plot  <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"_",input$epitope_hm,"_",input$pathology_hm,"_Heatmap_", gsub("/", "-", x), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Heatmap_epi_plot,
            height = input$height_png_Heatmap_epi_plot,
            res = input$resolution_PNG_Heatmap_epi_plot)
        plot(heatmap_epitope())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    #### epitope interrogation -----
    cols_epitope <- reactive({
      df3.meta <-  Add.UMAP.reduction()
      # summarising the clonality
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(df3.meta)>0 & nrow(epi)>0,
             "Upload Files")
      )
      if(input$datasource == "BD rhapsody") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

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
      df3.meta <-  Add.UMAP.reduction()
      # summarising the clonality
      epi <- input.data_sc_TCRex()

      validate(
        need(nrow(df3.meta)>0 & nrow(epi)>0,
             "Upload Files")
      )

      if(input$datasource == "BD rhapsody") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

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
      df3.meta <-  Add.UMAP.reduction()
      # summarising the clonality
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(df3.meta)>0 & nrow(epi)>0,
             "Upload Files")
      )
      if(input$datasource == "BD rhapsody") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")
      }
      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta,epi,by="CDR3_beta",all.x=T)
      df3.meta$selected <- df3.meta[,names(df3.meta) %in% input$epitope_umap_selected]
      df3.meta <- df3.meta[order(df3.meta$selected,decreasing = F),]
      df3.meta$selected <- factor(df3.meta$selected,levels = unique(df3.meta$selected))


      names(df3.meta)[names(df3.meta) %in% input$Samp_col] <- "ID_Column"
      df3.meta$ID_Column <- factor(df3.meta$ID_Column,levels = input$ID_Column_factor)

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
        df <- df + facet_wrap(~ID_Column)
      }
      df
    })
    output$UMAP_Epitope_plot <- renderPlot({
      UMAP_Epitope()
    })

    output$downloadPlot_UMAP_Epitope <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"_",input$epitope_umap_selected,"_UMAP_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_UMAP_Epitope,height=input$height_UMAP_Epitope, onefile = FALSE) # open the pdf device
        plot(UMAP_Epitope())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_UMAP_Epitope <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"_",input$epitope_umap_selected,"_UMAP_", gsub("/", "-", x), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_UMAP_Epitope,
            height = input$height_png_UMAP_Epitope,
            res = input$resolution_PNG_UMAP_Epitope)
        plot(UMAP_Epitope())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    # pie chart function -----
    cols_epitope_pie <- reactive({
      df3.meta <-  Add.UMAP.reduction()
      # summarising the clonality
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(df3.meta)>0 & nrow(epi)>0,
             "Upload Files")
      )
      if(input$datasource == "BD rhapsody") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

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
      df3.meta <-  Add.UMAP.reduction()
      # summarising the clonality
      epi <- input.data_sc_TCRex()

      validate(
        need(nrow(df3.meta)>0 & nrow(epi)>0,
             "Upload Files")
      )

      if(input$datasource == "BD rhapsody") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")

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
      df3.meta <-  Add.UMAP.reduction()
      # summarising the clonality
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(df3.meta)>0 & nrow(epi)>0,
             "Upload Files")
      )
      if(input$datasource == "BD rhapsody") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")
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
        facet_wrap(~Var1) +
        theme(
          legend.key.size = unit(1, 'cm'),
          legend.title = element_blank()) +
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 20),values = df.col, na.value = input$NA_col_analysis) +
        theme(strip.text = element_text(size = input$strip_size_pie, family = input$font_type),
              legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
              legend.position = input$legend_position,
              legend.title = element_blank()
        )

    })

    output$Pie_Epitope_plot <- renderPlot({
      Pie_chart_Epitope()
    })

    Pie_Epitope_dt_process <- reactive({
      df3.meta <-  Add.UMAP.reduction()
      epi <- input.data_sc_TCRex()
      validate(
        need(nrow(df3.meta)>0 & nrow(epi)>0,
             "Upload Files")
      )
      if(input$datasource == "BD rhapsody") {
        df3.meta$CDR3_beta <- paste("C",df3.meta$cdr3_BD,"F",sep="")
      }
      else {
        df3.meta$CDR3_beta <- df3.meta$cdr3_BD
      }
      epi$beta <- epi$CDR3_beta
      df3.meta <- merge(df3.meta,epi,by="CDR3_beta")
      df3.meta$Selected_function <- df3.meta[,names(df3.meta) %in% input$epitope_umap_selected]
      top_BD_cluster <-  df3.meta
      top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% input$epitope_umap_selected2]
      top_BD_cluster$cloneCount <- 1
      top_BD_cluster2 <-  top_BD_cluster[,names(top_BD_cluster) %in% c("cloneCount","Selected_function","Selected_group")]
      top_BD_cluster2 <- top_BD_cluster2 %>%
        select(cloneCount, everything())
      df2 <- as.data.frame(ddply(top_BD_cluster2,names(top_BD_cluster2)[-c(1)],numcolwise(sum)))
      df2$fraction <- df2$cloneCount/sum(df2$cloneCount)
      df2$Percent <- round(df2$cloneCount/sum(df2$cloneCount)*100,2)
      df2$Selected_function <- ifelse(grepl("NA", df2$Selected_function),NA,df2$Selected_function)
      df2[order(df2$Percent,decreasing = T),]
    })

    output$Pie_Epitope_dt <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{

      Pie_Epitope_dt_process()

    })

    #### download Epitope files ------
    output$downloadPlot_Pie_Epitope <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$epitope_Pie_Epitope,"_Pie_epi_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Pie_Epitope,height=input$height_Pie_Epitope, onefile = FALSE) # open the pdf device
        plot(Pie_chart_Epitope())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Pie_Epitope <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$epitope_umap_selected,"_Pie_epi_", gsub("/", "-", x), ".png", sep = "")
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
        x <- gsub(":", ".", Sys.time())
        paste(input$epitope_umap_selected,"_UMAP_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_UMAP_Epitope,height=input$height_UMAP_Epitope, onefile = FALSE) # open the pdf device
        plot(UMAP_Epitope())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_UMAP_Epitope <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$epitope_umap_selected,"_UMAP_", gsub("/", "-", x), ".png", sep = "")
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

    #### clusTCR2 figure -----
    clusTCR2_df <- reactive({
      clust <- input.data_sc_clusTCR()
      md <- Add.UMAP.reduction()
      validate(
        need(nrow(clust)>0 & nrow(md)>0,
             "Upload clusTCR table, which is needed for TCR -> UMAP section")
      )
      if (input$chain_TCR=="TCR") {
        md$CDR3_Vgene <- paste(md$cdr3_AG,md$v_gene_AG,sep="_")
        split <- as.data.frame(do.call(rbind, strsplit(as.character(md$CDR3_Vgene), "-")))[1]
        md$CDR3_Vgene <- split$V1
        df <- merge(md,clust,by = "CDR3_Vgene")

        md$CDR3_Vgene <- paste(md$cdr3_BD,md$v_gene_BD,sep="_")
        split <- as.data.frame(do.call(rbind, strsplit(as.character(md$CDR3_Vgene), "-")))[1]
        md$CDR3_Vgene <- split$V1
        df2 <- merge(md,clust,by = "CDR3_Vgene")
        cluster <- rbind(df,df2)
      }
      else if (input$chain_TCR=="BCR") {


        md$CDR3_Vgene <- paste(md$cdr3_IgL,md$v_gene_IgL,sep="_")
        split <- as.data.frame(do.call(rbind, strsplit(as.character(md$CDR3_Vgene), "-")))[1]
        md$CDR3_Vgene <- split$V1
        df3 <- merge(md,clust,by = "CDR3_Vgene")

        md$CDR3_Vgene <- paste(md$cdr3_IgH,md$v_gene_IgH,sep="_")
        split <- as.data.frame(do.call(rbind, strsplit(as.character(md$CDR3_Vgene), "-")))[1]
        md$CDR3_Vgene <- split$V1
        df4 <- merge(md,clust,by = "CDR3_Vgene")
        cluster <- rbind(df3,df4)
      }
      else {
        md$CDR3_Vgene <- paste(md$cdr3_AG,md$v_gene_AG,sep="_")
        split <- as.data.frame(do.call(rbind, strsplit(as.character(md$CDR3_Vgene), "-")))[1]
        md$CDR3_Vgene <- split$V1
        df <- merge(md,clust,by = "CDR3_Vgene")

        md$CDR3_Vgene <- paste(md$cdr3_BD,md$v_gene_BD,sep="_")
        split <- as.data.frame(do.call(rbind, strsplit(as.character(md$CDR3_Vgene), "-")))[1]
        md$CDR3_Vgene <- split$V1
        df2 <- merge(md,clust,by = "CDR3_Vgene")

        md$CDR3_Vgene <- paste(md$cdr3_IgL,md$v_gene_IgL,sep="_")
        split <- as.data.frame(do.call(rbind, strsplit(as.character(md$CDR3_Vgene), "-")))[1]
        md$CDR3_Vgene <- split$V1
        df3 <- merge(md,clust,by = "CDR3_Vgene")

        md$CDR3_Vgene <- paste(md$cdr3_IgH,md$v_gene_IgH,sep="_")
        split <- as.data.frame(do.call(rbind, strsplit(as.character(md$CDR3_Vgene), "-")))[1]
        md$CDR3_Vgene <- split$V1
        df4 <- merge(md,clust,by = "CDR3_Vgene")

        cluster <- rbind(df,df2,df3,df4)
      }
      cluster <- cluster[order(cluster$Clust_size_order),]
      cluster
    })





    output$Tb_ClusTCR_selected <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      clusTCR2_df()
    })

    # umap ClusTCR -----
    observe({
      clust <- input.data_sc_clusTCR()
      validate(
        need(nrow(clust)>0,
             error_message_val_sc)
      )
      clust <- clust[order(clust$Clust_size_order),]
      validate(
        need(nrow(clust)>0,
             "Upload clusTCR table, which is needed for TCR -> UMAP section")
      )
      updateSelectInput(
        session,
        "Clusters_to_dis",
        choices=unique(clust$Clust_size_order),
        selected = c(1,2)
      )
    }) # junction sequence

    output$Tb_ClusTCR_col <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
      cluster <- clusTCR2_df()
      md <- Add.UMAP.reduction()
      validate(
        need(nrow(cluster)>0 & nrow(md)>0,
             "Upload clusTCR table, which is needed for TCR -> UMAP section")
      )

      cluster <- cluster[order(cluster$Clust_size_order),]
      cluster$Clust_size_order <- factor(cluster$Clust_size_order, levels = unique(cluster$Clust_size_order))
      # cluster <- cluster[cluster$Clust_size_order %in% input$Clusters_to_dis,]
      col.df <- as.data.frame(unique(cluster$Clust_size_order))
      names(col.df) <- "V1"
      col.df$palette_rainbow <- rainbow(length(unique(cluster$Clust_size_order)))
      col.df2 <- col.df
      col.df2

    })


    UMAP_ClusTCR2 <- reactive({
      cluster <- clusTCR2_df()
      md <- Add.UMAP.reduction()
      if (input$ClusTCR_display == "all" ) {
        cluster <- cluster[order(cluster$Clust_size_order),]
        cluster$Clust_size_order <- factor(cluster$Clust_size_order, levels = unique(cluster$Clust_size_order))
        # cluster <- cluster[cluster$Clust_size_order %in% input$Clusters_to_dis,]
        col.df <- as.data.frame(unique(cluster$Clust_size_order))
        names(col.df) <- "V1"
        col.df$palette_rainbow <- rainbow(length(unique(cluster$Clust_size_order)))
        col.df2 <- col.df
        md$Clust_size_order <- NA
        md$Clust_size_order <- factor(md$Clust_size_order, levels = unique(cluster$Clust_size_order))

        figure <- ggplot()+
          geom_point(data=md,aes(x=UMAP_1,UMAP_2,color = Clust_size_order))+
          geom_point(data=cluster,aes(x=UMAP_1,UMAP_2,colour=Clust_size_order))+
          scale_color_manual(na.value="grey", values = c(col.df2$palette_rainbow),breaks = c(unique(col.df2$V1)))+
          # scale_size_manual(na.value=0.25,values = rep(3,dim(num)[1]))+
          theme_bw()+
          # labs(color=NULL,size = 12)+
          theme(
            legend.text = element_text(colour="black", size=24,family=input$font_type),
            legend.title = element_blank(),
            legend.position = "right",
          )

      }

      else {

        cluster <- cluster[order(cluster$Clust_size_order),]
        cluster$Clust_size_order <- factor(cluster$Clust_size_order, levels = unique(cluster$Clust_size_order))

        col.df <- as.data.frame(unique(cluster$Clust_size_order))
        names(col.df) <- "V1"
        col.df$palette_rainbow <- rainbow(length(unique(cluster$Clust_size_order)))
        col.df2 <- col.df[col.df$V1 %in% input$Clusters_to_dis,]

        md$Clust_size_order <- NA
        md$Clust_size_order <- factor(md$Clust_size_order, levels = unique(cluster$Clust_size_order))
        cluster <- cluster[cluster$Clust_size_order %in% input$Clusters_to_dis,]

        figure <- ggplot()+
          geom_point(data=md,aes(x=UMAP_1,UMAP_2,color = Clust_size_order))+
          geom_point(data=cluster,aes(x=UMAP_1,UMAP_2,colour=Clust_size_order))+
          scale_color_manual(na.value="grey", values = c(col.df2$palette_rainbow),breaks = c(unique(col.df2$V1)))+
          # scale_size_manual(na.value=0.25,values = rep(3,dim(num)[1]))+
          theme_bw()+
          # labs(color=NULL,size = 12)+
          theme(
            legend.text = element_text(colour="black", size=24,family=input$font_type),
            legend.title = element_blank(),
            legend.position = "right",
          )
      }

      figure

    })

    output$UMAP_ClusTCR2_plot <- renderPlot({
      UMAP_ClusTCR2()
    })
    output$downloadPlot_UMAP_ClusTCR2_plot <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste("_UMAP_ClusTCR2_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_UMAP_ClusTCR2_plot,height=input$height_UMAP_ClusTCR2_plot, onefile = FALSE) # open the pdf device
        plot(UMAP_ClusTCR2())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_UMAP_ClusTCR2_plot <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste("_UMAP_ClusTCR2_", gsub("/", "-", x), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_UMAP_ClusTCR2_plot,
            height = input$height_png_UMAP_ClusTCR2_plot,
            res = input$resolution_PNG_UMAP_ClusTCR2_plot)
        plot(UMAP_ClusTCR2())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    ##### motif plot ----
    observe({
      clust <- input.data_sc_clusTCR()
      validate(
        need(nrow(clust)>0,
             "upload clustering")
      )

      clust <- clust[order(clust$Clust_size_order),]
      updateSelectInput(
        session,
        "Clusters_to_dis_motif",
        choices=unique(clust$Clust_size_order),
        selected = 5
      )
    }) # cluster to display
    motif_plot_sc <- reactive({
      Network_df <- input.data_sc_clusTCR()
      validate(
        need(nrow(Network_df)>0,
             "upload ClusTCR2 file")
      )
      Motif_from_cluster_file(Network_df,Clust_selected = input$Clusters_to_dis_motif)
    })

    output$Motif_ClusTCR2_cluster <- renderPlot({
      motif_plot_sc()
    })

    # Ridge ClusTCR -----
    Ridge_clusTCR2 <- reactive({
      cluster <- clusTCR2_df()
      md <- Add.UMAP.reduction()
      if (input$ClusTCR_display == "all" ) {
        cluster <- cluster[order(cluster$Clust_size_order),]
        cluster$Clust_size_order <- factor(cluster$Clust_size_order, levels = unique(cluster$Clust_size_order))
        # cluster <- cluster[cluster$Clust_size_order %in% input$Clusters_to_dis,]
        col.df <- as.data.frame(unique(cluster$Clust_size_order))
        names(col.df) <- "V1"
        col.df$palette_rainbow <- rainbow(length(unique(cluster$Clust_size_order)))
        col.df2 <- col.df
        md$Clust_size_order <- NA
        md$Clust_size_order <- factor(md$Clust_size_order, levels = unique(cluster$Clust_size_order))

        figure <- ggplot()+
          geom_point(data=md,aes(x=UMAP_1,UMAP_2,color = Clust_size_order))+
          geom_point(data=cluster,aes(x=UMAP_1,UMAP_2,colour=Clust_size_order))+
          scale_color_manual(na.value="grey", values = c(col.df2$palette_rainbow),breaks = c(unique(col.df2$V1)))+
          # scale_size_manual(na.value=0.25,values = rep(3,dim(num)[1]))+
          theme_bw()+
          # labs(color=NULL,size = 12)+
          theme(
            legend.text = element_text(colour="black", size=24,family=input$font_type),
            legend.title = element_blank(),
            legend.position = "right",
          )

      }
      else {

        cluster <- cluster[order(cluster$Clust_size_order),]
        cluster$Clust_size_order <- factor(cluster$Clust_size_order, levels = unique(cluster$Clust_size_order))

        col.df <- as.data.frame(unique(cluster$Clust_size_order))
        names(col.df) <- "V1"
        col.df$palette_rainbow <- rainbow(length(unique(cluster$Clust_size_order)))
        col.df2 <- col.df[col.df$V1 %in% input$Clusters_to_dis,]

        md$Clust_size_order <- NA
        md$Clust_size_order <- factor(md$Clust_size_order, levels = unique(cluster$Clust_size_order))
        cluster <- cluster[cluster$Clust_size_order %in% input$Clusters_to_dis,]

        figure <- ggplot()+
          geom_point(data=md,aes(x=UMAP_1,UMAP_2,color = Clust_size_order))+
          geom_point(data=cluster,aes(x=UMAP_1,UMAP_2,colour=Clust_size_order))+
          scale_color_manual(na.value="grey", values = c(col.df2$palette_rainbow),breaks = c(unique(col.df2$V1)))+
          # scale_size_manual(na.value=0.25,values = rep(3,dim(num)[1]))+
          theme_bw()+
          # labs(color=NULL,size = 12)+
          theme(
            legend.text = element_text(colour="black", size=24,family=input$font_type),
            legend.title = element_blank(),
            legend.position = "right",
          )
      }

      figure

    })

    output$ridge_ClusTCR2_plot <- renderPlot({
      Pie_ClusTCR2()
    })

    ###

    ##### Overlap -----
    #### upset plot -----

    Upset_plot <- reactive({
      df <-  Add.UMAP.reduction()
      validate(
        need(nrow(df)>0,
             "upload file")
      )

      df <-  Add.UMAP.reduction()
      df <- as.data.frame(df)
      unique.df <- unique(df[,names(df) %in% c(input$Samp_col,input$V_gene_sc) ])
      names(unique.df) <- c("group","chain")
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
                      column_names_gp = gpar(fontfamily = input$font_type),
                      top_annotation = upset_top_annotation(df.x,
                                                            add_numbers = T,
                                                            annotation_name_gp = gpar(fontfamily = input$font_type),
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
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"_Overlap_",gsub("/", "-", x), ".pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width=input$width_Upset_plot_overlap,height=input$height_Upset_plot_overlap, onefile = FALSE) # open the pdf device
        plot(Upset_plot())
        dev.off()}, contentType = "application/pdf" )

    output$downloadPlotPNG_Upset_plot_overlap <- downloadHandler(
      filename = function() {
        x <- gsub(":", ".", Sys.time())
        paste(input$name.file_clust,"_Overlap_", gsub("/", "-", x), ".png", sep = "")
      },
      content = function(file) {
        png(file, width = input$width_png_Upset_plot_overlap,height = input$height_png_Upset_plot_overlap,res = input$resolution_PNG_Upset_plot_overlap)
        plot(Upset_plot())
        dev.off()},   contentType = "application/png" # MIME type of the image
    )

    # upset plot table -----

    Upset_plot_overlap <- reactive({
      df <-  Add.UMAP.reduction()
      validate(
        need(nrow(df)>0,
             "upload file")
      )

      df <-  Add.UMAP.reduction()
      df <- as.data.frame(df)
      unique.df <- unique(df[,names(df) %in% c(input$Samp_col,input$V_gene_sc) ])
      names(unique.df) <- c("group","chain")
      unique.df <- subset(unique.df,unique.df$chain != "NA")
      unique.df <- subset(unique.df,unique.df$group != "NA")
      unique.df$cloneCount <- 1
      unique.df
      mat <- acast(unique.df, chain~group, value.var="cloneCount")
      mat[is.na(mat)] <- 0
      sum_data <- as.data.frame(rowSums(mat))
      names(sum_data) <- "V1"

      mat <- as.data.frame(mat)

      mat$sum <-sum_data$V1
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


    # overlap UMAP plot

    create_UMAP_overlap <- reactive({
      sc <- input.data_sc_pro()

      validate(
        need(nrow(sc)>0,
             error_message_val_UMAP)
      )

      reduction <- (sc@reductions$umap)
      UMAP <- as.data.frame(reduction@cell.embeddings)
      UMAP$Cell_Index <- rownames(UMAP)
      meta.data <- as.data.frame(sc@meta.data)
      if(input$STEGO_R_pro == "STEGO_R (.h5Seurat)") {
        meta.data <- meta.data
        umap.meta <- merge(UMAP,meta.data,by="Cell_Index")
        names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
        names(umap.meta)[names(umap.meta) %in% input$V_gene_sc] <- "v_gene_selected"

      }
      else {
        meta.data$Cell_Index <- rownames(meta.data)
        meta.data$v_gene_AG <- meta.data[,names(meta.data) %in% input$RDS_V_gene_A]
        meta.data$v_gene_BD <- meta.data[,names(meta.data) %in% input$RDS_V_gene_B]
        meta.data$cdr3_AG <- meta.data[,names(meta.data) %in% input$RDS_cdr3_A]
        meta.data$cdr3_BD <- meta.data[,names(meta.data) %in% input$RDS_cdr3_B]
        meta.data$v_gene_cdr3_BD <- paste(meta.data$v_gene_BD,meta.data$cdr3_BD,sep="_")
        meta.data$v_gene_cdr3_AG<- paste(meta.data$v_gene_AG,meta.data$cdr3_AG,sep="_")
        meta.data$v_gene_cdr3_AG_BD <- paste(meta.data$v_gene_cdr3_AG,meta.data$v_gene_cdr3_BD,sep=" & ")

        umap.meta <- merge(UMAP,meta.data,by="Cell_Index")
        names(umap.meta)[names(umap.meta) %in% input$Samp_col] <- "ID_Column"
        names(umap.meta)[names(umap.meta) %in% "v_gene_cdr3_AG_BD"] <- "v_gene_selected"

      }


      sc <- merge(umap.meta,TCR_Expansion(),by=c("v_gene_selected","ID_Column"),all.x=T)



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
    output$create_UMAP_overlap <- renderPlot({
      create_UMAP_overlap()
    })

    ### end -----
  }
  shinyApp(ui, server)

}
