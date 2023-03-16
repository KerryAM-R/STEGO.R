#' Run STEGO application.
#' @name STEGO
#' @export runSTEGO



runSTEGO <- function(...)  {

  options(shiny.maxRequestSize = 20000*1024^2)
# ?numericInput
# UI page -----
ui <- fluidPage(
  theme=bs_theme(version = 5, bootswatch = "default"),
  navbarPage(title = "STEGO_R",
             theme=bs_theme(version = 5, bootswatch = "default"),
             navbarMenu("Quality control",
## BD Rhapsody  ------
                        tabPanel("BD rhapsody data",
                                 sidebarLayout(
                                   sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                     # UPLOAD the three files...
                                     selectInput("dataset_BD", "Choose a dataset:", choices = c("test_data_BD", "own_data_BD")),
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
                                       tabPanel("clusTCR",
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
## 10x Genomics ----
                        tabPanel("10x genomics",
                                 sidebarLayout(
                                   sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                     selectInput("dataset_10x", "Choose a dataset:", choices = c("test_data_10x", "own_data_10x")),
                                     fileInput('file_calls_10x', 'Barcode file (.tsv.gz or .tsv)',
                                               accept=c('.tsv','.tsv.gz')),
                                     fileInput('file_features_10x', 'Features file (.tsv.gz or .tsv)',
                                               accept=c('.tsv','.tsv.gz')),
                                     fileInput('file_matrix_10x', 'Matrix file (.mtx.gz or .mtx)',
                                               accept=c('.mtx.gz','.mtx')),
                                     fileInput('file_TCR_10x', 'filtered contig annotations (.csv)',
                                               accept=c('.csv')),
                                     conditionalPanel(condition="input.panel_10x==2 | input.panel_10x==3",
                                     textInput("sample_name_10x","Add sample name","Treatment_group")
                                                      ),

                                     conditionalPanel(condition="input.panel_10x==4",
                                                      fluidRow(
                                                        column(6,textInput("group10x","Add sample name","Group")),
                                                        column(6,textInput("Indiv10x","Add sample name","Indiv"))
                                                              )
                                     ),

                                     conditionalPanel(condition="input.panel_10x==2",
                                                      fluidRow(column(6,downloadButton('downloadtb_10x_matrix2','Download matrix')),
                                                               column(6,downloadButton('downloadtb_10x_metadata2','Download metadata'))
                                                               )
                                                      ),
                                     conditionalPanel(condition="input.panel_10x==3",
                                                      downloadButton('downloadtb_10x_contigues1','Download clusTCR')
                                     ),
                                     conditionalPanel(condition="input.panel_10x==4",
                                                      downloadButton('downloaddt_TCR_Explore_10x','Download TCR_Explore')
                                     ),
                                     selectInput("BCR_TCR_10x","Type of data",choices = c("TCR only","BCR only","both")),

                                     textInput("name.10x","Name added to files",value = ""),

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

                                                ),
                                       tabPanel("ClusTCR",value = 3,
                                                tags$head(tags$style("#tb_10x_contigues1  {white-space: nowrap;  }")),
                                                div(DT::dataTableOutput("tb_10x_contigues1")),
                                                ),
                                       tabPanel("TCR_Explore",value = 4,
                                                tags$head(tags$style("#dt_TCR_Explore_10x  {white-space: nowrap;  }")),
                                                div(DT::dataTableOutput("dt_TCR_Explore_10x")),

                                                ),
                                     ),
                                   )
                                 )
                        ),
# 10x genomics end -----
             ), # NavbarMenu
### TCR clustering with ClusTCR2 -----
              tabPanel("ClusTCR2",
                sidebarLayout(

                    sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,
                                 selectInput("dataset2", "Choose a dataset:", choices = c("test_data_clusTCR2","own_data_clusTCR2")),
                                 # selectInput("dataset_clusTCR2", "Choose a dataset:", choices = c("test_data_clusTCR2","own_data_clusTCR2")),
                                 fileInput('file2_ClusTCR2', 'Select file for single samples',
                                           accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
                                 #
                                 # fluidRow(
                                 #   column(6,radioButtons('sep_clusTCR2', 'Separator', c( Tab='\t', Comma=','), ',')),
                                 #   column(6,radioButtons('quote_clusTCR2', 'Quote', c(None='', 'Double Quote'='"', 'Single Quote'="'"), '"'))
                                 # ),
                                 actionButton("run_ClusTCR2","Update Clustering"),
                                 fluidRow(
                                   column(6,selectInput( "clusTCR2_names",label = h5("CDR3"),"")),
                                   column(6,selectInput( "clusTCR2_Vgene",label = h5("V gene"),"")),
                                   # column(6,selectInput( "clusTCR2_group_column",label = h5("Group"),"")),
                                   # checkboxInput("include_V_gene", "Include V gene for clustering?", value = T, width = NULL),
                                 ),
                                 fluidRow(
                                   column(6,checkboxInput("allele_ClusTCR2","Remove allele *00?", value = T)),
                                   column(6, numericInput("cores_ClusTCR2","Number of cores to parallel",value=4))
                                 )

                ),
                mainPanel(
                  tabsetPanel(
                    tabPanel("Uploaded file",
                             div(DT::dataTableOutput("clust_dt2")),
                             ),
                    tabPanel("Outputs",
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
                          selectInput("dataset_sc", "Choose a dataset:", choices = c("test_data_sc", "own_data_sc")),
                          # upload the file
                          fileInput('file_SC', 'Unprocessed filtered file (.csv.gz or .csv)',
                                    accept=c('.csv','.csv.gz')),
                          # selectInput("species","Species",choices = c("human","mouse","other")),
                          selectInput("df_seruatobj_type","Data type", choices = c("10x","BD rhapsody"),selected = "BD rhapsody"),

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
                          selectInput("save_type","Type of output",choices = c(".h5",".rds")),
                           downloadButton('downloaddf_SeruatObj','Download Seurat')
                            # column(4,numericInput("percent.mt","Mictochondrial DNA cut-off", value = 20)),
                        ),
                        mainPanel(
                          tabsetPanel(
# QC panel -----
                            tabPanel("QC plots",
                                     tabsetPanel(
                                       tabPanel("Violin and correlation",
                                               tabsetPanel(
                                                 tabPanel("Before",
                                                          add_busy_spinner(spin = "fading-circle"),
                                                          plotOutput("before_plot_sc", height = "600px"),
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

                                               # # fluidRow(
                                               # #   column(3,numericInput("width_png_UMAP2","Width of PNG", value = 1200)),
                                               # #   column(3,numericInput("height_png_UMAP2","Height of PNG", value = 1000)),
                                               # #   column(3,numericInput("resolution_PNG_UMAP2","Resolution of PNG", value = 144)),
                                               # #   column(3,style = "margin-top: 25px;",downloadButton('downloadPlotPNG_UMAP2','Download PNG'))
                                               #
                                               # ),

                                       ),
#### Variable features -----
                                       tabPanel("Top variable features",
                                                add_busy_spinner(spin = "fading-circle"),
                                                plotOutput("plot_10_features_sc", height = "600px")
                                       ),
#### Elbow and heatmap  -----
                                       tabPanel("Elbow Plot",
                                                plotOutput("create_elbowPlot_sc", height = "600px")

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
                                    tabPanel("Cell annotations",
                                             tabsetPanel(
                                               tabPanel("CellTypist",
                                                        fluidRow(

                                                        column(3,checkboxInput("cellTypist_add","Add in CellTypist classification (human)", value = F)),
                                                        column(9, selectInput("cellTypistModels_selected","Cell typist Models",
                                                                              choices = cellTypistModels, selected = "Immune_All_Low.pkl",multiple = T))
                                                        ),
                                                        add_busy_spinner(spin = "fading-circle"),
                                                        div(DT::dataTableOutput("DEx_table_cellTypist")),

                                                        ),

                                               tabPanel("K-means clustering",
                                                        tabsetPanel(
                                                          tabPanel("Classification to include",
                                                                  selectInput("V_gene_Class","V gene with/without CDR3",choices = ""),
                                                                  p("Lituature Curated list"),
                                                                  fluidRow(column(3,checkboxInput("add.classification_T_cell","General T cell markers", value = T)),
                                                                           column(3,checkboxInput("add.classification_T_cell_Function","Function", value = T)),
                                                                           column(3,checkboxInput("add.classification_T_cell_Function_CD4","Function (CD4)", value = T)),
                                                                           column(3,checkboxInput("add.classification_T_cell_Function_CD8","Function (CD8)", value = T)),
                                                                           column(3,checkboxInput("add.classification_T_cell_Function_CD4_CD8pos","Function (CD4+CD8+)", value = T)),
                                                                           column(3,checkboxInput("add.classification_T_cell_Function_CD4_CD8neg","Function (CD4-CD8-)", value = T)),
                                                                           column(3,checkboxInput("add.classification_B_cell_Function","Function B cell", value = T)),
                                                                           column(3,checkboxInput("add.classification_T_cell_Memory","Memory", value = T)),
                                                                           column(3,checkboxInput("add.classification_T_cell_Activation","Activation status", value = T)),),
                                                                  h5("Cell Typist based list"),
                                                                  fluidRow(

                                                                    column(3,checkboxInput("add.classification_CellTypist_list_overview","CellTypist Overview", value = T)),

                                                                    column(3,checkboxInput("add.classification_CellTypist_list","CellTypist based classification", value = T)),
                                                                    column(3,checkboxInput("add.classification_CellTypist_list_CD8","CellTypist CD8 classification", value = T)),
                                                                    column(3,checkboxInput("add.classification_CellTypist_list_CD4","CellTypist CD4 classification", value = T)),
                                                                    column(3,checkboxInput("add.classification_CellTypist_cycling","CellTypist Cell cycle", value = T)),

                                                                  ),
                                                                 add_busy_spinner(spin = "fading-circle"),
                                                                 div(DT::dataTableOutput("DEx_table_TcellClass")),


                                                                  ),
                                                            ),
                                                        ),
                                               tabPanel("ProjecTILs") # need to work out how to add this.
                                             )
),
# end of QC -----
                          ),
                        )
                      )
             ),
## Analysis (UI side panel)---------
             tabPanel("Analysis",
                      sidebarLayout(
                        sidebarPanel(id = "tPanel4",style = "overflow-y:scroll; max-height: 800px; position:relative;", width=3,

                        selectInput("STEGO_R_pro","QC processed",choices = c("STEGO_R (.h5Seurat)","Seurat (.rds)")),
                        selectInput("dataset_sc_pro", "Choose a dataset:", choices = c("test_data_sc_pro", "own_data_sc_pro"), selected = "test_data_sc_pro"),

                         conditionalPanel(condition="input.check_up_files== 'up'",


                          fileInput('file_SC_pro', 'Upload seurat file',
                                    accept=c('.h5Seurat','.rds')),
                          fileInput('file_cluster_file', 'Upload clustering file from clusTCR (.csv)',
                                    accept=c('.csv')),
                          fileInput('upload_TCRex_file', 'Upload TCRex file from clusTCR (.tsv)',
                                    accept=c('.tsv')),
                          textInput("name.file_clust","Name added to files",value = ""),
                                     ),
                        conditionalPanel(condition="input.STEGO_R_pro == 'STEGO_R (.h5Seurat)'",
                          selectInput("datasource", "Data source",choices=c("BD rhapsody","10x Genomics")),
                        ),
                          fluidRow(column(6,selectInput("Samp_col","Sample column name",choices = "")),
                          conditionalPanel(condition="input.STEGO_R_pro == 'STEGO_R (.h5Seurat)'",
                                   column(6,selectInput("V_gene_sc","V gene with/without CDR3",choices = "")))
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




                          selectInput("font_type",label = h4("Type of font"),choices = font,selected = "serif"),

                          # conditionalPanel(condition="input.QC_panel==2 || input.QC_panel==3",
                                           fluidRow(
                                             column(6,numericInput("text_size","Size of #",value=16)),
                                             column(6,numericInput("title.text.sizer2","Axis text size",value=30)),

                                             column(6,numericInput("Bar_legend_size","Legend text size",value=8)),
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
                                                column(4,selectInput("Split_by_group","Include group comparison",choices=c("no","yes"))),
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
### TCR -> UMAP TCR structure -> back to expression data -----


                                  tabPanel("TCR interrogation",

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

##### Classification to include ------
                                     tabsetPanel(id = "Panel_TCRUMAP",

# T cell classification ------
                                       tabPanel("Overview",
                                                tabsetPanel(id = "Panel_class",
                                                            tabPanel("Meta data",value = 16,
                                                                     add_busy_spinner(spin = "fading-circle"),
                                                                     div(DT::dataTableOutput("UMAP_tb_download")),
                                                            ),
                                                 tabPanel("Percentage",value = 16,


                                                          add_busy_spinner(spin = "fading-circle"),
                                                                     div(DT::dataTableOutput("Percent_tab")),
                                                            ),


                                                  # tabPanel("Table",value = 13,
                                                  #          add_busy_spinner(spin = "fading-circle"),
                                                  #          div(DT::dataTableOutput("test.files_classification")),),
                                                  tabPanel("UMAP plot",value = 14,
                                                           fluidRow(column(3,selectInput("show_selected","Show all labels?",choices=c("All","Selected_list"))),
                                                                    column(9,uiOutput("SiteNumInput")),
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


                                                             ),
                                                           fluidRow(column(3,
                                                                           wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                                     uiOutput('myPanel_pie'))), # pie chart
                                                                    column(9, plotOutput("Classification_clonotype_pie",height="600px"))),

                                                           # plotOutput("",height="600px")
                                                           ),


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

                                                           ),
                                                           fluidRow(
                                                             column(3,
                                                                    wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                                                              uiOutput('myPanel_Top_pie_clonotype'))),
                                                             column(9, plotOutput("top_clonotype_pie",height="600px"))

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
                                                           div(DT::dataTableOutput("Violin_chart_alpha_gamma")),
                                                           fluidRow(
                                                             column(3, checkboxInput("restric_ex","Restrict to above a threshold?", value = F )),
                                                             column(3, numericInput("Gre_ex","Expression above:", value = 0 )),
                                                             column(3, selectInput("plot_type_ridgvi","Plot type", choices = c("Ridge (selected clonotype)","Ridge (compare)","Violin (selected clonotype)", "Violin (compare)"))),
                                                           ),
                                                           actionButton("run_string.data_Exp_top","View Ridge plot"),
                                                           fluidRow(column(12, selectInput("string.data_Exp_top","column names for summary","",multiple = F, width = "1200px") )),
                                                           add_busy_spinner(spin = "fading-circle"),
                                                           plotOutput("Ridge_chart_alpha_gamma_plot_out",height="600px"),
                                                           conditionalPanel(condition="input.plot_type_ridgvi=='Ridge (selected clonotype)' | input.plot_type_ridgvi=='Violin (selected clonotype)'",
                                                                            div(DT::dataTableOutput("Ridge_chart_alpha_gamma_stat"))),
                                                           conditionalPanel(condition="input.plot_type_ridgvi=='Ridge (compare)' | input.plot_type_ridgvi=='Violin (compare)'",
                                                                            div(DT::dataTableOutput("Ridge_chart_alpha_gamma_stat_comp"))),
                                                  ),
                                                    ),
                                                  ),
# clusTCR specific analysis
                                        # tabPanel("clusTCR"),
# epitope analysis -----
tabPanel("Epitope",
         tabsetPanel(
           tabPanel("Uploaded Epitope file",
                    add_busy_spinner(spin = "fading-circle"),


                    div(DT::dataTableOutput("MainTcell_Check")),
           ),

           tabPanel("Heatmap",
                    fluidRow(
                      column(3,selectInput("epitope_hm","y-axis",choices = c("CDR3_beta","epitope","pathology"),selected = "epitope")),
                      column(3,selectInput("pathology_hm","x-axis",choices = c("CDR3_beta","epitope","pathology"),selected = "pathology")),

                             ),

                    plotOutput("Heatmap_epi_plot",height="600px"),

                    ),
           tabPanel("UMAP",
                    column(3,selectInput("epitope_umap_selected","Select",choices = c("beta","epitope","pathology"),selected = "pathology")),
                    fluidRow(
                      column(3,
                             wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                                       uiOutput('myPanel_cols_epitope'))),
                      column(9, plotOutput("UMAP_Epitope_plot",height="600px"))
                    ),
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
                    # div(DT::dataTableOutput("Tb_ClusTCR_col")),


                    fluidRow(
                      # column(3,
                      #        wellPanel(id = "tPanel23",style = "overflow-y:scroll; max-height: 600px",
                      #                  uiOutput('myPanel_cols_epitope'))),
                      column(9, plotOutput("UMAP_ClusTCR2_plot",height="600px"))
                    ),

           ),
           tabPanel("motif",
                    plotOutput("Motif_ClusTCR2_cluster",height="600px")
                    )



         )

)
                                       # tabPanel("Clonotypes per cluster (Pie/bar plot)"),
                                       # tabPanel("Upset plot")
                                      )
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
    df3.meta <- df3.meta[!grepl("RNA",df3.meta) & !grepl("BCR",df3.meta) &!grepl("orig.ident",df3.meta) & !grepl("TCR",df3.meta)& !grepl("gene",df3.meta) & !grepl("allele",df3.meta) & !grepl("percent",df3.meta) & !grepl("cdr3",df3.meta)]
    selectInput("clust_names_top","Select function type:",choices = df3.meta,selected="classify.T.cell")

  })


  # user interface parameters-----
  output$feature_input <- renderUI({
    if (input$df_seruatobj_type=="10x") {
      fluidRow(
        column(6,numericInput("features.min","minimum features (<)", value = 200)),
        column(6,numericInput("features.max","Maximum features (<)", value = 6000)),
        column(6,numericInput("percent.mt","Mictochondrial DNA cut-off (<)", value = 20)),
        column(6,numericInput("percent.rb","Ribosomal RNA cut-off (>)", value = 0)),
      )
    }

    else {
      fluidRow(
        column(6,numericInput("features.min","minimum features (<)", value = 45)),
        column(6,numericInput("features.max","Maximum features (<)", value = 160)),
        column(6,numericInput("percent.mt","Mictochondrial DNA cut-off (<)", value = 20)),
        column(6,numericInput("percent.rb","Ribosomal RNA cut-off (>)", value = 0)),
      )
    }

  })
# human BD rhapsody data -----
  ## three files required for BD data: Sample Tag calls, TCR file and count ----
    input.data.calls.bd <- reactive({switch(input$dataset_BD,"test_data_BD" = test.data.calls.bd(), "own_data_BD" = own.data.calls.bd())})
    test.data.calls.bd <- reactive({
        dataframe = read.csv(system.file("extdata","BDrhap/demo/Sample_Tags 2023.02.27.csv",package = "STEGO.R"))

       # read.csv(system.file("../Public data/Bd Rhapsody/BD-Demo-VDJ-2/001493_BD-Demo-VDJ/Sample_Tags 2023.02.27.csv"))
      })
    own.data.calls.bd <- reactive({
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

    input.data.TCR.BD <- reactive({switch(input$dataset_BD,"test_data_BD" = test.data.TCR.bd(), "own_data_BD" = own.data.TCR.bd())})
    test.data.TCR.bd <- reactive({
      # dataframe = read.csv("../Nic BD data/gd-T-cell-CBZ_VDJ_perCell.csv")
        dataframe = read.csv(system.file("extdata","BDrhap/demo/RhapVDJDemo-BCR_VDJ_perCell.csv",package = "STEGO.R"),skip = 7)
      # dataframe = read.csv("../Public data/Bd Rhapsody/BD-Demo-VDJ-2/001493_BD-Demo-VDJ/",skip = 7)
      })
    own.data.TCR.bd <- reactive({
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

    input.data.count.BD <- reactive({switch(input$dataset_BD,"test_data_BD" = test.data.count.bd(), "own_data_BD" = own.data.count.bd())})
    test.data.count.bd <- reactive({
      # dataframe = read.csv("../Nic BD data/Combined_gd-T-cell-CBZ_DBEC_MolsPerCell.csv")
        dataframe = read.csv(system.file("extdata","BDrhap/demo/RhapVDJDemo-BCR_DBEC_MolsPerCell.csv",package = "STEGO.R"),skip = 7)
        # dataframe = read.csv("../Public data/Bd Rhapsody/BD-Demo-VDJ-2/001493_BD-Demo-VDJ/RhapVDJDemo-BCR_DBEC_MolsPerCell.csv",skip = 7)
    })
    own.data.count.bd <- reactive({
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

      paired <- paired[!names(paired) %in% c("TRAG", "TRBD","TRAG_fun" , "TRBD_fun"   ),]
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
        paste(input$name.BD," Sample_Tags ",gsub("-", ".", Sys.Date()),".csv", sep = "")
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


    # output$downloadtb_tb_clusTCR_sum <- downloadHandler(
    #   filename = function(){
    #     paste(input$name.BD," BD__",gsub("-", ".", Sys.Date()),".csv", sep = "")
    #   },
    #   content = function(file){
    #     df <- as.data.frame(sum_tb_BD_2())
    #     # write.table(,file, row.names = T)
    #     write.csv(df, file, row.names = T)
    #
    #   })

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
        paste(input$name.BD," Meta.data ",gsub("-", ".", Sys.Date()),".csv", sep = "")
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
        paste(input$name.BD," TCRex ",gsub("-", ".", Sys.Date()),".tsv", sep = "")
      },
      content = function(file){
        df <- as.data.frame(TCRex_BDrap_df())
        df <- df[,names(df) %in% c("TRBV_gene","CDR3_beta","TRBJ_gene")]
        write.table(df,file, row.names = F,sep="\t", quote = F)
      } )

    # count matrix download ----
    output$downloadtb_count_matrix <- downloadHandler(
      filename = function(){
        paste(input$name.BD," BD_Count_Matrix_",gsub("-", ".", Sys.Date()),".csv", sep = "")
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
  input.data.barcode.10x <- reactive({switch(input$dataset_10x,"test_data_10x" = test.data.barcode.10x(), "own_data_10x" = own.data.barcode.10x())})
  test.data.barcode.10x <- reactive({

    dataframe = read.table(system.file("extdata","10x/RAW/GSM4455933_K409_LN_GEX_barcodes.tsv.gz",package = "STEGO.R"))
    # dataframe = read.table("../RHM003/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz")
  })
  own.data.barcode.10x <- reactive({
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
  input.data.features.10x <- reactive({switch(input$dataset_10x,"test_data_10x" = test.data.features.10x(), "own_data_10x" = own.data.features.10x())})
  test.data.features.10x <- reactive({
    dataframe = read.table(system.file("extdata","10x/RAW/GSM4455933_K409_LN_GEX_genes.tsv.gz",package = "STEGO.R"))
    # dataframe = read.table("../RHM003/count/sample_filtered_feature_bc_matrix/features.tsv.gz")
  })
  own.data.features.10x <- reactive({
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
  input.data.matrix.10x <- reactive({switch(input$dataset_10x,"test_data_10x" = test.data.matrix.10x(), "own_data_10x" = own.data.matrix.10x())})
  test.data.matrix.10x <- reactive({
    dataframe = Matrix::readMM(system.file("extdata","10x/RAW/GSM4455933_K409_LN_GEX_matrix.mtx.gz",package = "STEGO.R"))
  })
  own.data.matrix.10x <- reactive({
    inFile_10x_matrix <- input$file_matrix_10x
    if (is.null(inFile_10x_matrix)) return(NULL)

    else {
      dataframe <- Matrix::readMM(
        inFile_10x_matrix$datapath)}
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
  input.data.TCR.10x <- reactive({switch(input$dataset_10x,"test_data_10x" = test.data.TCR.10x(), "own_data_10x" = own.data.TCR.10x())})
  test.data.TCR.10x <- reactive({
    dataframe = read.csv(system.file("extdata","10x/RAW/GSM4455934_K409_LN_VDJ_filtered_contig_annotations.csv.gz",package = "STEGO.R"))
  })
  own.data.TCR.10x <- reactive({
    inFile_10x_TCR <- input$file_TCR_10x
    if (is.null(inFile_10x_TCR)) return(NULL)

    else {
      dataframe <- read.csv(
        inFile_10x_TCR$datapath)}
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
    names(contig_LK)[1:summary(name.list)[1]] <-paste(names(contig_LK[names(contig_LK) %in% name.list]),"_LK",sep="")
    contig_LK

    contig_H <- subset(contigs_lim,contigs_lim$chain=="IGH")

    name.list <- names(contig_H[c(names(contig_H[grep("gene",names(contig_H))]),
                                  names(contig_H[grep("cdr3",names(contig_H))]),
                                  "chain")])
    contig_H <- contig_H %>%
      select(all_of(name.list), everything())


    names(contig_H)[1:summary(name.list)[1]] <-paste(names(contig_H[names(contig_H) %in% name.list]),"_H",sep="")
    contig_H
    # contig_paired <- merge(contig_LK,contig_H, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
    # contig_paired <- merge(contig_LK,contig_H, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
    contig_paired <- merge(contig_LK,contig_H, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)

    contig_paired$pairing <- ifelse(contig_paired$chain_H=="IGH" & contig_paired$chain_LK=="IGK","IGK Paired",
                                    ifelse(contig_paired$chain_H=="IGH" & contig_paired$chain_LK=="IGL","IGL Paired",NA
                                    ))

    contig_paired
    contig_paired$pairing[is.na(contig_paired$pairing)] <- "unpaired"
    contig_paired <- contig_paired[!names(contig_paired) %in% c("d_gene_LK")]

    dim(contig_paired)

    contig_paired_only <- contig_paired
    contig_paired_only <- subset(contig_paired_only,contig_paired_only$cdr3_H!="None")

    contig_paired_only <- subset(contig_paired_only,contig_paired_only$cdr3_LK!="None")
    dim(contig_paired_only)

    contig_paired_only$d_gene_H <- sub("^$","NA", contig_paired_only$d_gene_H)
    #
    contig_paired_only$vj_gene_LK <- paste(contig_paired_only$v_gene_LK,contig_paired_only$j_gene_LK,sep = ".")
    contig_paired_only$vj_gene_LK <- gsub("NA.NA","",contig_paired_only$vj_gene_LK)
    #
    contig_paired_only$vj_gene_H <- paste(contig_paired_only$v_gene_H,contig_paired_only$j_gene_H,sep = ".")
    contig_paired_only$vj_gene_H <- gsub(".NA.",".",contig_paired_only$vj_gene_H)
    contig_paired_only$vj_gene_H <- gsub(".None.",".",contig_paired_only$vj_gene_H)
    contig_paired_only$vj_gene_H <- gsub("NA.NA","",contig_paired_only$vj_gene_H)
    #
    contig_paired_only$vdj_gene_H <- paste(contig_paired_only$v_gene_H,contig_paired_only$d_gene_H,contig_paired_only$j_gene_H,sep = ".")
    contig_paired_only$vdj_gene_H <- gsub(".NA.",".",contig_paired_only$vdj_gene_H)
    contig_paired_only$vdj_gene_H <- gsub(".None.",".",contig_paired_only$vdj_gene_H)
    contig_paired_only$vdj_gene_H <- gsub("NA.NA","",contig_paired_only$vdj_gene_H)
    #
    contig_paired_only$vj_gene_cdr3_LK <- paste(contig_paired_only$vj_gene_LK,contig_paired_only$cdr3_LK,sep = "_")
    contig_paired_only$vj_gene_cdr3_LK <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_LK)
    #
    contig_paired_only$vj_gene_cdr3_H <- paste(contig_paired_only$vj_gene_H,contig_paired_only$cdr3_H,sep = "_")
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
    write.csv(contig_paired_only,"contig_paired_only.csv")
    dup <- contig_paired_only[duplicated(contig_paired_only$Cell_Index),]
    contig_paired_only <- contig_paired_only[order(contig_paired_only$Cell_Index, contig_paired_only$umis.x,contig_paired_only$umis.y,decreasing = T),]

    contig_paired_only_dup <- contig_paired_only[!duplicated(contig_paired_only$Cell_Index),] # remove duplicate barcodes.
    names(contig_paired_only_dup)
    contig_paired_only_dup <- contig_paired_only_dup[!names(contig_paired_only_dup) %in% c("umis.x","umis.y")]
    contig_paired_only_dup$Sample_Name <- input$sample_name_10x

    contig_paired_only_dup <- contig_paired_only_dup %>%
      select(all_of(c("Cell_Index","Sample_Name")), everything())

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

   }

  })
  output$downloadtb_10x_metadata2 <- downloadHandler(
    filename = function(){
      paste(input$name.10x," metadata_10x_",gsub("-", ".", Sys.Date()),".csv", sep = "")
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
      write_csv(df, file)

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
      merge.names <- names(contig_paired_only_dup)[!grepl("_H",names(contig_paired_only_dup)) & !grepl("_LK",names(contig_paired_only_dup))]
      TCR <- tb_10x_meta.data_TCR()
      BCR <- tb_10x_meta.data_BCR()
      merge(TCR,BCR,by = merge.names, all=T)

    }

    # contig_paired <- as.data.frame(tb_10x_meta.data())

    count.df <- contig_paired_only_dup[names(contig_paired_only_dup) %in% c(names(contig_paired_only_dup)[grepl("chain",names(contig_paired_only_dup))],"pairing")]
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
        contigs2 <- subset(contigs2,contigs2$junction_aa!= "None")
        contigs2[!duplicated(contigs2[,c('v_call','junction_aa')]),]
  })

  output$tb_10x_contigues1 <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
      tb_10x_contigues_contig()
    })

  output$downloadtb_10x_contigues1 <- downloadHandler(
    filename = function(){
      paste(input$name.10x," clusTCR_10x_",gsub("-", ".", Sys.Date()),".csv", sep = "")
    },
    content = function(file){
      df <- as.data.frame(tb_10x_contigues_contig())
      write.csv(df,file, row.names = F)
    } )


  ## TCR explore 10x -----
  TCR_Explore_10x <- function () {
    contig_paired_only <- tb_10x_meta.data()
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
      paste(input$name.10x,"TCR_Explore_10x  ",gsub("-", ".", Sys.Date()),".csv", sep = "")
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
      paste(input$name.10x," count-matrix_10x_",gsub("-", ".", Sys.Date()),".csv.gz", sep = "")
    },
    content = function(file){
      df <- as.data.frame(tb_10x_matrix())
      # write.table(,file, row.names = T)
      write_csv(df, file)

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


# ClusTCR2 ------
  input.data_ClusTCR2 <- reactive({switch(input$dataset2,"test_data_clusTCR2" = test.data_clusTCR2(),"own_data_clusTCR2" = own.data_clusTCR2())})
  test.data_clusTCR2 <- reactive({
      # system.file("extdata","inst/extdata/clusTCR/cdr3.csv",package = "STEGO.R")
      dataframe = read.csv(system.file("extdata","clusTCR/cdr3.csv",package = "STEGO.R"))
    # dataframe = read.csv("../cdr3.csv")
  })
  own.data_clusTCR2  <- reactive({
    inFile2_ClusTCR2 <- input$file2_ClusTCR2
    if (is.null(inFile2_ClusTCR2)) return(NULL)
    else {
      dataframe <- read.csv(
        inFile2_ClusTCR2$datapath,
        header=TRUE,
        # sep=input$sep_ClusTCR2,
        # quote=input$quote_ClusTCR2
        )}
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
    # registerDoParallel(cores=8)
    ptm <- proc.time()
    clust_dt2 <- clust_dt2()

    df_cluster <- ClusTCR(clust_dt2, allele =input$allele_ClusTCR2, v_gene = input$clusTCR2_Vgene, cores_selected = cores_ClusTCR2)
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
      paste("Cluster_table  ",gsub("-", ".", Sys.Date()),".csv", sep = "")
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
   netplot(Network_df,
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
  input.data_sc <- reactive({switch(input$dataset_sc,"test_data_sc" = test.data_sc(),"own_data_sc" = own.data_sc())})
  test.data_sc <- reactive({
    if (input$df_seruatobj_type =="10x") {
        # dataframe = read.csv(system.file("extdata","clusTCR/cdr3.csv",package = "STEGO.R"))
      dataframe = read.csv(system.file("extdata","10x/Mahuron_Melanoma_2020/LN/K409_LN.count-matrix_10x_2023.03.08.csv.gz",package = "STEGO.R"))
      # head(dataframe)[1:6,1:6]
      # dataframe = read.csv("../filename.csv.gz")
    }
    else {
      dataframe = read.csv(system.file("extdata","BDrhap/Seurat/BD_Count_Matrix_2023.03.07.csv",package = "STEGO.R"),row.names = 1)
      # dataframe = read.csv("../Public data/Bd Rhapsody/QC_output/BD_Count_Matrix_2023.02.27.csv",row.names = 1)
    }
  })
  own.data_sc <- reactive({
    inFile_sc <- input$file_SC
    if (is.null(inFile_sc)) return(NULL)
    else {
      if (input$df_seruatobj_type =="10x") {
        dataframe = read.csv(inFile_sc$datapath)
      }
      else {
        dataframe = read.csv(inFile_sc$datapath,row.names = 1)
      }

    }

  })
    ## reading in 10x and BD data ----
  df_seruatobj <- reactive({
    df.test <- input.data_sc()
    validate(
      need(nrow(df.test)>0,
           error_message_val_sc)
    )

    if (input$df_seruatobj_type =="10x") {
      names(df.test) <- gsub("[.]1","-1",names(df.test) )
      rownames(df.test) <- make.unique(df.test$Gene_Name)
      df.test2 <- df.test[,!names(df.test) %in% c("Gene_Name","Gene_ID")]
      sc <- CreateSeuratObject(counts = df.test2, project = "rhm003")
      sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
      sc[["percent.rb"]] <- PercentageFeatureSet(sc, pattern = "^RP[SL]")
      sc
    }
    else {
      names(df.test) <- gsub("X","",names(df.test))
      df.test2 <- df.test[!rownames(df.test) %in% c("Cell_Index"),]
      sc <- CreateSeuratObject(counts = df.test2, project = "Test.data")
      sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
      sc[["percent.rb"]] <- PercentageFeatureSet(sc, pattern = "^RP[SL]")
      sc
    }



  })
  before_plot <- reactive({
    sc <- df_seruatobj()
    VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 2)

  })

  output$before_plot_sc <- renderPlot({
    withProgress(message = 'Figure is being generated...',
                 detail = '', value = 0, {
                   test_fun()
                 })
    before_plot()
  })

  vals2 <- reactiveValues(after_violin_plot=NULL)
  observeEvent(input$run,{
      sc <- df_seruatobj()
      vals2$after_violin_plot <- subset(sc, subset = nFeature_RNA >= input$features.min & nFeature_RNA <= input$features.max & percent.mt <= input$percent.mt & percent.rb >= input$percent.rb)
      vals2$after_violin_plot

    })

  output$after_plot_sc <- renderPlot({

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

  # output$downloadPlotPNG_simpson.inv <- downloadHandler(
  #   filename = function() {
  #     x <- gsub(":", ".", Sys.time())
  #     paste("inv.simpson.index.", gsub("/", "-", x), ".png", sep = "")
  #   },
  #   content = function(file) {
  #
  #     png(file, width = input$width_png_simpson.inv, height = input$height_png_simpson.inv, res = input$resolution_PNG_simpson.inv)
  #     grid.arrange(print(group.diversity1()),print(group.diversity2()),ncol=2)
  #     # group.diversity1()
  #     dev.off()}, contentType = "application/png" # MIME type of the image
  # )

  ### normalisationa and feature plot ------
  feature_serartobj <- function (){
    sc <- vals2$after_violin_plot
    sc <- NormalizeData(sc)
    sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
  }

  plot_10_features <- reactive({
    sc <- feature_serartobj()
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



    ## PCA and chosing # of dimentions to reduce ----
  create_PCA <- reactive({
    sc <- feature_serartobj()
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

  create_elbowPlot <- reactive({
    ElbowPlot(create_PCA())
  })
  output$create_elbowPlot_sc<- renderPlot({
    create_elbowPlot()
  })

  vals_clust <- reactiveValues(sc_clustering=NULL)

  observeEvent(input$run_reduction,{
      sc <- create_PCA()
      sc<- FindNeighbors(sc, dims = 1:input$dimension_sc)
      sc <- FindClusters(sc, resolution = input$resolution)
      sc <- RunUMAP(sc, dims = 1:input$dimension_sc)
      vals_clust$sc_clustering <- sc
      vals_clust$sc_clustering
  })

  create_UMAP <- reactive({
  DimPlot(vals_clust$sc_clustering, reduction = "umap")
  })
  output$create_UMAP_sc<- renderPlot({
    create_UMAP()
  })

    ## Differential expression with two conditions -----
  input.data_sc_meta <- reactive({switch(input$dataset_sc,"test_data_sc" = test.data_sc_meta(),"own_data_sc" = own.data_sc_meta())})
  test.data_sc_meta <- reactive({
    if (input$df_seruatobj_type =="10x") {
        # dataframe = read.csv(system.file("extdata","clusTCR/cdr3.csv",package = "STEGO.R"))
      dataframe = read.csv(system.file("extdata","10x/Mahuron_Melanoma_2020/LN/K409_LN metadata_10x_2023.03.08.csv",package = "STEGO.R"))
      # dataframe = read.csv("../Test.metadata/metadata_10x_2022.12.14.csv")
    }
    else {
      # dataframe = read.csv(system.file("extdata","BDrhap/Seurat/BD_Count_Matrix_2023.03.07.csv",package = "STEGO.R"),row.names = 1)
      dataframe = read.csv(system.file("extdata","BDrhap/Seurat/Meta.data 2023.03.07.csv",package = "STEGO.R"))
      # dataframe = read.csv("../Public data/Bd Rhapsody/QC_output/Meta.data 2023.02.27.csv")
    }
  })

  own.data_sc_meta <- reactive({
    inFile_sc_meta <- input$file_SC_meta
    if (is.null(inFile_sc_meta)) return(NULL)
    else {
      if (input$df_seruatobj_type =="10x") {
        dataframe = read.csv(inFile_sc_meta$datapath)
      }
      else {
        dataframe = read.csv(inFile_sc_meta$datapath)
      }

    }

  })

  output$DEx_view.meta.dt <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    input.data_sc_meta()
  })

  vals_meta.sc <- reactiveValues(metadata_SCobj=NULL)

  observeEvent(input$run_metadata,{
      sc <- vals_clust$sc_clustering
      meta.data.import <- input.data_sc_meta()
      validate(
        need(nrow(meta.data.import)>0,
             error_message_val_sc)
      )
      sc@meta.data$Cell_Index <- rownames(sc@meta.data)
      scMeta.data <- sc@meta.data
      meta.data2 <- merge(scMeta.data,meta.data.import,by="Cell_Index",all.x=T)
      sc@meta.data <- meta.data2
      rownames(sc@meta.data) <- sc@meta.data$Cell_Index
      sc@meta.data <- sc@meta.data[order(as.numeric(sc@meta.data$Cell_Index)),]
      vals_meta.sc$metadata_SCobj <- sc
      vals_meta.sc$metadata_SCobj
  })

  output$DEx_table_meta.data <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    calls <- vals_meta.sc$metadata_SCobj
    calls@meta.data
  })
# Add cell annotations to Seurat object -----

  # add.CellTypist <- reactive({
  #   sc <- vals_meta.sc$metadata_SCobj
  #
  #   if(input$cellTypist_add==T){
  #     sc <- RIRA::RunCellTypist(sc, runCelltypistUpdate = F, modelName = input$cellTypistModels_selected)
  #     sc@meta.data$CellTypist <- sc@meta.data$majority_voting
  #
  #   }
  #   else(
  #     sc <- sc
  #   )
  #   sc
  # })
# cell Typist
  output$DEx_table_cellTypist <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
    calls <- add.CellTypist()
    calls@meta.data
  })

  #
  observe({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           error_message_val_UMAP)
    )
    df3.meta <- sc@meta.data
      updateSelectInput(
        session,
        "V_gene_Class",
        choices=names(df3.meta),
        selected = "v_gene_cdr3_AB_GD")
  })

  add.classification <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )

    df= as.data.frame(sc[["RNA"]]@scale.data)
    MainTcell <- as.data.frame(t(df))
    MainTcell
    # names(MainTcell) <- ifelse(grepl("pAbO",names(MainTcell)),names(MainTcell),toupper(names(MainTcell)))
    MainTcell$Cell_Index <- rownames(MainTcell)
    # names(MainTcell) <- gsub("CD8B.","CD8B",names(MainTcell))
    md <- sc@meta.data
    md$Motif_gene <- md[,names(md) %in% input$V_gene_Class]

    md$Chain2 <- ifelse(grepl("_._", md$Motif_gene),0,
                        ifelse(md$Motif_gene=="_",0,
                               ifelse(md$Motif_gene=="",0,
                                      ifelse(md$Motif_gene==" & ",0,
                                             ifelse(grepl("TRBV",md$Motif_gene) |
                                                      grepl("TRAV",md$Motif_gene),1,
                                                    ifelse(grepl("TRGV",md$Motif_gene) |
                                                             grepl("TRDV",md$Motif_gene),-1,0))))))

    md2 <- md[names(md) %in% c("Cell_Index","Chain2")]
    MainTcell <- merge(md2,MainTcell,by="Cell_Index",all.y=T)
    rownames(MainTcell) <- MainTcell$Cell_Index
    MainTcell

  })
  add.classification_T_cell_Df <- reactive({

    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )


    MainTcell <- add.classification()
    if (input$add.classification_T_cell==T) {

      MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD4","CD8A","CD8B","CD3D","CD3E","CD3G","MS4A1")]
      head(MainTcell_test)
      df <- MainTcell_test
      nms <- c("CD4","CD8A","CD8B","CD3D","CD3E","CD3G","MS4A1")   # Vector of columns you want in this data
      Missing <- setdiff(nms, names(df))  # Find names of missing columns
      df[Missing] <- 0                    # Add them, filled with '0's
      MainTcell_test <- df[nms]

      MainTcell_test[is.na(MainTcell_test)] <- 0
      set.seed(123)
      kmeans10 <- kmeans(MainTcell_test, centers = 20, nstart = 10)
      kmeans10$centers
      df_centres <- as.data.frame(kmeans10$centers)
      df_centres <- df_centres %>%
        mutate(classify.T.cell = case_when(
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 >0 & CD8B > 0 & CD8A > 0 ~ "CD4+CD8ab+",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 >0 & CD8B > 0 & CD8A < 0 ~ "CD4+CD8b+",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 >0 & CD8B < 0 & CD8A > 0 ~ "CD4+CD8a+",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 >0 & CD8B < 0 & CD8A < 0 ~ "CD4+",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B > 0 & CD8A > 0 ~ "CD8ab+",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B < 0 & CD8A > 0 ~ "CD8a+",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B > 0 & CD8A < 0 ~ "CD8b+",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4 <0 & CD8B < 0 & CD8A < 0 ~ "CD4-CD8-",
          # CD3D < 0 & CD3E<0 & CD3G<0 & CD4 <0 & CD8B < 0 & CD8A < 0~ "Not a T cell",
          CD3D > 0 | CD3E>0 | CD3G>0 ~ "T cell",
          MS4A1>0 ~ "B cell",
          JCHAIN>0 ~ "Plasma cell",

          TRUE ~ NA_character_))

      df_centres$clust_name <- rownames(df_centres)
      head(df_centres)
      df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","classify.T.cell")]
      Cluster_Df <- as.data.frame( kmeans10$cluster)
      names(Cluster_Df) <- "clust_name"
      Cluster_Df$Cell_Index <- rownames(Cluster_Df)
      Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
      Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","classify.T.cell")]
      Cluster_Df_class2 <- Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
      Cluster_Df_class2

    }
    else {
      Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
      names(Cluster_Df_class2) <- "Cell_Index"
      Cluster_Df_class2
    }

  })
  add.classification_T.B_cell_Function_Df <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )
    MainTcell <- add.classification()

    if (input$add.classification_T_cell_Function==T) {

      rownames(MainTcell) <- MainTcell$Cell_Index
      MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                                                           "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM")]
      df <- MainTcell_test
      # Vector of columns you want in this data
      nms <- c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
               "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM")
      Missing <- setdiff(nms, names(df))  # Find names of missing columns
      df[Missing] <- 0                    # Add them, filled with '0's
      MainTcell_test <- df[nms]

      d_frame <- na.omit(MainTcell_test)
      set.seed(123)
      kmeans10 <- kmeans(d_frame, centers = 500, nstart = 1)
      df_centres <- as.data.frame(kmeans10$centers)

      df_centres <- df_centres %>%
        mutate(classify.Adaptive.cell_Function = case_when(
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
          IGHM>0 & MS4A1>0 ~ "IgM+ B cell",
          IGHM<0 & MS4A1>0~ "IgM- B cell",
          c(IGHA2>0 | IGHA1>0) & JCHAIN>0 ~ "IgA+ Plasma cell",
          c(IGHG1>0 | IGHG2>0| IGHG3>0| IGHG4>0) & JCHAIN>0 ~ "IgG+ Plasma cell",
          IGHE>0 & JCHAIN>0 ~ "IgG+ Plasma cell",
          IGHD>0 & MS4A1>0 ~ "IgD+ B cell",
          TRUE ~ NA_character_))

      df_centres$clust_name <- rownames(df_centres)
      head(df_centres)
      df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","classify.Adaptive.cell_Function")]
      Cluster_Df <- as.data.frame( kmeans10$cluster)
      Cluster_Df
      names(Cluster_Df) <- "clust_name"
      Cluster_Df$Cell_Index <- rownames(Cluster_Df)
      Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
      Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","classify.Adaptive.cell_Function")]
      Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]

    }
    else {
      Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
      names(Cluster_Df_class2) <- "Cell_Index"
      Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
    }


  })
  add.classification_T_cell_Function_CD4_df <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )
    MainTcell <- as.data.frame(add.classification())

    if (input$add.classification_T_cell_Function_CD4==T) {
      MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                                                           "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM")]
      df <- MainTcell_test
      # Vector of columns you want in this data
      nms <- c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
               "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM")
      Missing <- setdiff(nms, names(df))  # Find names of missing columns
      df[Missing] <- 0                    # Add them, filled with '0's
      MainTcell_test <- df[nms]
      MainTcell_test[is.na(MainTcell_test)] <- 0
      MainTcell_test
      set.seed(123)
      kmeans10 <- kmeans(MainTcell_test, centers = 50, nstart = 1)
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
          CD4>0 & MKI67 >0 & TOP2A > 0  ~ "Cycling CD4+ T cells",
          CD4>0 & TNF>0 ~ "TNF+CD4+",
          CD4>0 & IFNG>0 ~ "IFNg+CD4+",
          CD4 >0 & CD8A<0 & CD8B<0  ~ "CD4+",
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
  add.classification_T_cell_Function_CD8_df <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )
    MainTcell <- add.classification()

    if (input$add.classification_T_cell_Function_CD8==T) {
      rownames(MainTcell) <- MainTcell$Cell_Index
      MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                                                           "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM")]
      df <- MainTcell_test
      # Vector of columns you want in this data
      nms <- c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
               "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM")
      Missing <- setdiff(nms, names(df))  # Find names of missing columns
      df[Missing] <- 0                    # Add them, filled with '0's
      MainTcell_test <- df[nms]

      MainTcell_test[is.na(MainTcell_test)] <- 0
      # kmeans10 <-pamk(d_frame, krange= 500:1000, ns=1000)
      set.seed(123)
      kmeans10 <- kmeans(MainTcell_test, centers = 500, nstart = 1)
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
  add.classification_T_cell_Function_CD4_CD8pos_df <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )
    MainTcell <- add.classification()

    if (input$add.classification_T_cell_Function_CD4_CD8pos==T) {
      rownames(MainTcell) <- MainTcell$Cell_Index
      MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                                                           "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM")]
      df <- MainTcell_test
      # Vector of columns you want in this data
      nms <- c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
               "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM")
      Missing <- setdiff(nms, names(df))  # Find names of missing columns
      df[Missing] <- 0                    # Add them, filled with '0's
      MainTcell_test <- df[nms]

      MainTcell_test[is.na(MainTcell_test)] <- 0
      set.seed(123)
      kmeans10 <- kmeans(MainTcell_test, centers = 500, nstart = 1)
      df_centres <- as.data.frame(kmeans10$centers)

      df_centres <- df_centres %>%
        mutate(classify.T.cell_Function_CD4_CD8pos = case_when(
          CD4>0 &c(CD8A >0 | CD8B >0) & MKI67 >0 & TOP2A > 0  ~ "Cycling CD4+CD8+ T cells",
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
  add.classification_T_cell_Function_CD4_CD8neg_df <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )
    MainTcell <- add.classification()
    if (input$add.classification_T_cell_Function_CD4_CD8neg==T) {
      rownames(MainTcell) <- MainTcell$Cell_Index
      MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                                                           "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM")]
      df <- MainTcell_test
      # Vector of columns you want in this data
      nms <- c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
               "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM")
      Missing <- setdiff(nms, names(df))  # Find names of missing columns
      df[Missing] <- 0                    # Add them, filled with '0's
      MainTcell_test <- df[nms]

      MainTcell_test[is.na(MainTcell_test)] <- 0
      set.seed(123)
      kmeans10 <- kmeans(MainTcell_test, centers = 500, nstart = 1)
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


  add.classification_B_cell_Function_Df <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )
    MainTcell <- add.classification()

    if (input$add.classification_B_cell_Function==T) {

      rownames(MainTcell) <- MainTcell$Cell_Index
      MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
                                                           "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM")]
      df <- MainTcell_test
      # Vector of columns you want in this data
      nms <- c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21","CCR6",
               "GZMB","PRF1","GNLY","IL6","KLRB1","TGFB1","IGHM")
      Missing <- setdiff(nms, names(df))  # Find names of missing columns
      df[Missing] <- 0                    # Add them, filled with '0's
      MainTcell_test <- df[nms]

      d_frame <- na.omit(MainTcell_test)
      set.seed(123)
      kmeans10 <- kmeans(d_frame, centers = 500, nstart = 1)
      df_centres <- as.data.frame(kmeans10$centers)

      df_centres <- df_centres %>%
        mutate(classify.B.cell_Function = case_when(
          IGHM>0 & MS4A1>0 ~ "IgM+ B cell",
          IGHM<0 & MS4A1>0~ "IgM- B cell",
          c(IGHA2>0 | IGHA1>0) & JCHAIN>0 ~ "IgA+ Plasma cell",
          c(IGHG1>0 | IGHG2>0| IGHG3>0| IGHG4>0) & JCHAIN>0 ~ "IgG+ Plasma cell",
          IGHE>0 & JCHAIN>0 ~ "IgG+ Plasma cell",
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

  add.classification_T_cell_Memory_df <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )
    MainTcell <- add.classification()
    if (input$add.classification_T_cell_Memory==T) {
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
      set.seed(123)
      kmeans10 <- kmeans(MainTcell_test, centers = 1000, nstart = 10)
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
  add.classification_T_cell_Activation_df <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )
    MainTcell <- add.classification()
    if (input$add.classification_T_cell_Activation==T) {
      rownames(MainTcell) <- MainTcell$Cell_Index
      MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD8A","CD8B","CD4","IL2RA","CD27","TNFRSF8","DPP4","CD40LG","TNFRSF4","ICAM1","PDCD1","B3GAT1")]

      df <- MainTcell_test
      nms <- c("CD8A","CD8B","CD4","IL2RA","CD27","TNFRSF8","DPP4","CD40LG","TNFRSF4","ICAM1","PDCD1","B3GAT1")   # Vector of columns you want in this data
      Missing <- setdiff(nms, names(df))  # Find names of missing columns
      df[Missing] <- 0                    # Add them, filled with '0's
      MainTcell_test <- df[nms]
      head(MainTcell_test)
      MainTcell_test[is.na(MainTcell_test)] <- 0
      set.seed(123)
      kmeans10 <- kmeans(MainTcell_test, centers = 500, nstart = 10)
      kmeans10$centers
      df_centres <- as.data.frame(kmeans10$centers)

      df_centres <- df_centres %>%
        mutate(classify.T.cell_Activation = case_when(
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
          TRUE ~ NA_character_))

      df_centres$clust_name <- rownames(df_centres)
      head(df_centres)
      df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","classify.T.cell_Activation")]
      Cluster_Df <- as.data.frame( kmeans10$cluster)
      names(Cluster_Df) <- "clust_name"
      Cluster_Df$Cell_Index <- rownames(Cluster_Df)
      Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
      Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","classify.T.cell_Activation")]
      Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
    }
    else {
      Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
      names(Cluster_Df_class2) <- "Cell_Index"
      Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
    }


  })
  add.classification_CellTypist_list_lower <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )
    MainTcell <- add.classification()
    if (input$add.classification_CellTypist_list_overview==T) {
      rownames(MainTcell) <- MainTcell$Cell_Index
      MainTcell_test <- MainTcell[,names(MainTcell) %in% c("MKI67","TOP2A","Chain2","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD", "IKZF2", "MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","TBX21","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","CD27","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","CLEC10A","S100A9","EBI3","CCR7","CCL19","CLEC10A","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","GNLY","FCGR3A","NKG7","S100A13","TLE1","AREG","CXCR3","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","GZMK","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","FCGR3A","C1QA","CX3CR1","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1")]
      # -----
      df <- MainTcell_test
      # Vector of columns you want in this data
      nms <- c("MKI67","TOP2A","Chain2","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD", "IKZF2", "MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","TBX21","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","CD27","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","CLEC10A","S100A9","EBI3","CCR7","CCL19","CLEC10A","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","GNLY","FCGR3A","NKG7","S100A13","TLE1","AREG","CXCR3","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","GZMK","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","FCGR3A","C1QA","CX3CR1","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1")

      Missing <- setdiff(nms, names(df))  # Find names of missing columns
      df[Missing] <- 0                    # Add them, filled with '0's
      MainTcell_test <- df[nms]
      MainTcell_test[is.na(MainTcell_test)] <- 0
      set.seed(123)
      kmeans10 <- kmeans(MainTcell_test, centers = 200, nstart = 1)
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
          JCHAIN>0 & MZB1>0 & XBP1>0      ~ "Plasma cell"

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
  add.classification_CellTypist_list_higher <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )
    MainTcell <- add.classification()
    if (input$add.classification_CellTypist_list_overview==T) {
      rownames(MainTcell) <- MainTcell$Cell_Index
      MainTcell_test <- MainTcell[,names(MainTcell) %in% c("MKI67","TOP2A","Chain2","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD", "IKZF2", "MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","TBX21","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","CD27","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","CLEC10A","S100A9","EBI3","CCR7","CCL19","CLEC10A","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","GNLY","FCGR3A","NKG7","S100A13","TLE1","AREG","CXCR3","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","GZMK","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","FCGR3A","C1QA","CX3CR1","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1")]
      # -----
      df <- MainTcell_test
      # Vector of columns you want in this data
      nms <- c("MKI67","TOP2A","Chain2","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD", "IKZF2", "MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","TBX21","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","CD27","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","CLEC10A","S100A9","EBI3","CCR7","CCL19","CLEC10A","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","GNLY","FCGR3A","NKG7","S100A13","TLE1","AREG","CXCR3","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","GZMK","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","FCGR3A","C1QA","CX3CR1","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1")
      Missing <- setdiff(nms, names(df))  # Find names of missing columns
      df[Missing] <- 0                    # Add them, filled with '0's
      MainTcell_test <- df[nms]

      MainTcell_test[is.na(MainTcell_test)] <- 0
      set.seed(123)
      kmeans10 <- kmeans(MainTcell_test, centers = 200, nstart = 1)
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

          # MKI67 >0 & TOP2A > 0 & Chain2==1 ~ "Cycling ab T cells",
          # MKI67 >0 & TOP2A > 0 & Chain2== -1 ~ "Cycling gd T cells",
          #### T cells ----
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

  add.classification_CellTypist_list_df <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )
    MainTcell <- add.classification()
    if (input$add.classification_CellTypist_list==T) {
      rownames(MainTcell) <- MainTcell$Cell_Index
      MainTcell_test <- MainTcell[,names(MainTcell) %in% c("MKI67","TOP2A","Chain2","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD", "IKZF2", "MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G")]

      df <- MainTcell_test
      # Vector of columns you want in this data
      nms <- c("MKI67","TOP2A","Chain2","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD", "IKZF2", "MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G")
      Missing <- setdiff(nms, names(df))  # Find names of missing columns
      df[Missing] <- 0                    # Add them, filled with '0's
      MainTcell_test <- df[nms]

      MainTcell_test[is.na(MainTcell_test)] <- 0
      # kmeans10 <-pamk(d_frame, krange= 500:1000, ns=1000)
      kmeans10 <- kmeans(MainTcell_test, centers = 1000, nstart = 1)
      df_centres <- as.data.frame(kmeans10$centers)

      df_centres <- df_centres %>%
        mutate(CellTypist_list = case_when(
          # MKI67 >0 & TOP2A > 0 & Chain2==1 ~ "Cycling ab T cells",
          # MKI67 >0 & TOP2A > 0 & Chain2== -1 ~ "Cycling gd T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & PDCD1 >0 & ZNF683 >0 & CD8A >0 ~ "CD8aa T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & PDCD1 >0 & ICOS >0 & CXCR5>0  ~"Follicular helper T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & GZMK >0 & CD4>0 & IL10>0 ~ "Memory GZMK+IL10+CD4+ cytotoxic T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & GZMK >0 & CD4>0 & IL10>0 ~ "Memory GZMK+IL10-CD4+ cytotoxic T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & NKG7>0 & GNLY>0 & CD8A>0 & Chain2==1 ~ "NK CD8+ ab T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & NKG7>0 & GNLY>0 & CD8A>0 & Chain2== -1 ~ "NK CD8+ gd T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & NKG7>0 & GNLY>0 & CD4>0 & Chain2==1 ~ "NK CD4+ ab T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & NKG7>0 & GNLY>0 & CD8A>0 & Chain2== -1 ~ "NK CD8- gd T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & CD4>0 &CTLA4>0 & IL2RA >0 & FOXP3>0 ~ "Regulatory T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &KLRB1 >0 & SLC4A10 >0 & `TRAV1-2`>0  ~ "KLRB1+SLC4A10+ MAIT cells (TRAV1-2+)",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &KLRB1 >0 & SLC4A10 >0 & `TRAV1-2`>0  ~ "KLRB1+ MAIT cells (TRAV1-2+)",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &KLRB1 >0 & SLC4A10 >0 & `TRAV1-2`==0 ~ "MAIT cells (TRAV1-2neg)",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &CD8A>0 & CCR7>0 & SELL>0 ~ "Tcm/Naive cytotoxic T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &CD4>0 & CCR7>0 & SELL>0 ~ "Tcm/Naive helper T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &KLRB1>0 & AQP3>0 & ITGB1>0 ~ "Tem/Effector helper T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &PDCD1>0 & CD4>0 & CTLA4>0 ~ "Tem/Effector helper T cells PD1+",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &CX3CR1>0 & GZMB>0 & GNLY>0 ~ "Tem/Temra cytotoxic T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &GZMK>0 & CD8A>0 & CCL5>0 ~ "Tem/Trm cytotoxic CD8+ T cells",
          (CD3D>0 | CD3E >0 |CD3G>0 ) &GZMK>0 & CD4>0 & CCL5>0 ~ "Tem/Trm cytotoxic CD4+ T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &CD27>0 & CCR7>0 & IKZF2>0 ~ "Treg(diff)",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &ITGA1>0 & ITGAE>0 & CXCR6>0 ~"Trm cytotoxic T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &CCL5>0 & CXCR3>0 & TBX21>0 & CD4>0 ~ "Type 1 helper CD4+ ab T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &CCL5>0 & CXCR3>0 & TBX21>0 & Chain2== -1 ~ "Type 1 helper gd T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &IL7R>0 & CCR6>0 & ZBTB16>0 ~ "Type 17 helper T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &ITGAD>0 & Chain2== -1 & IKZF2>0 ~ "CRTAM+ gamma-delta T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &MIR155HG>0 & BIRC3>0 & SMS>0 ~ "T(agonist)",

          TRUE ~ NA_character_))

      df_centres$clust_name <- rownames(df_centres)
      head(df_centres)
      df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","CellTypist_list")]
      Cluster_Df <- as.data.frame( kmeans10$cluster)
      Cluster_Df
      names(Cluster_Df) <- "clust_name"
      Cluster_Df$Cell_Index <- rownames(Cluster_Df)
      Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
      Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","CellTypist_list")]
      Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
    }
    else {
      Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
      names(Cluster_Df_class2) <- "Cell_Index"
      Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
    }


  })
  add.classification_CellTypist_list_CD8_df <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )
    MainTcell <- add.classification()
    if (input$add.classification_CellTypist_list_CD8==T) {
      rownames(MainTcell) <- MainTcell$Cell_Index
      MainTcell_test <- MainTcell[,names(MainTcell) %in% c("MKI67","TOP2A","Chain2","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD", "IKZF2", "MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G")]

      df <- MainTcell_test
      # Vector of columns you want in this data
      nms <- c("MKI67","TOP2A","Chain2","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD", "IKZF2", "MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G")
      Missing <- setdiff(nms, names(df))  # Find names of missing columns
      df[Missing] <- 0                    # Add them, filled with '0's
      MainTcell_test <- df[nms]

      MainTcell_test[is.na(MainTcell_test)] <- 0
      set.seed(123)
      kmeans10 <- kmeans(MainTcell_test, centers = 1000, nstart = 1)
      df_centres <- as.data.frame(kmeans10$centers)

      df_centres <- df_centres %>%
        mutate(CellTypist_list_CD8 = case_when(
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &PDCD1 >0 & ZNF683 >0 & CD8A >0 ~ "CD8aa T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &NKG7>0 & GNLY>0 & CD8A>0 & Chain2==1 ~ "NK CD8+ ab T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & NKG7>0 & GNLY>0 & CD8A>0 & Chain2== -1 ~ "NK CD8+ gd T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &CTLA4>0 & IL2RA >0 & FOXP3>0 & CD8A>0 ~ "Regulatory CD8+ T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &CD8A>0 & CCR7>0 & SELL>0 ~ "Tcm/Naive cytotoxic T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &KLRB1>0 & AQP3>0 & ITGB1>0 ~ "Tem/Effector helper T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & CX3CR1>0 & GZMB>0 & GNLY>0 ~ "Tem/Temra cytotoxic T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &GZMK>0 & CD8A>0 & CCL5>0 ~ "Tem/Trm cytotoxic CD8+ T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &ITGA1>0 & ITGAE>0 & CXCR6>0 ~"Trm cytotoxic T cells",
          TRUE ~ NA_character_))

      df_centres$clust_name <- rownames(df_centres)
      head(df_centres)
      df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","CellTypist_list_CD8")]
      Cluster_Df <- as.data.frame( kmeans10$cluster)
      Cluster_Df
      names(Cluster_Df) <- "clust_name"
      Cluster_Df$Cell_Index <- rownames(Cluster_Df)
      Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
      Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","CellTypist_list_CD8")]
      Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
    }
    else {
      Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
      names(Cluster_Df_class2) <- "Cell_Index"
      Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
    }


  })
  add.classification_CellTypist_list_CD4_df <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )
    MainTcell <- add.classification()
    if (input$add.classification_CellTypist_list_CD4==T) {
      rownames(MainTcell) <- MainTcell$Cell_Index
      MainTcell_test <- MainTcell[,names(MainTcell) %in% c("MKI67","TOP2A","Chain2","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD", "IKZF2", "MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G")]

      df <- MainTcell_test
      # Vector of columns you want in this data
      nms <- c("MKI67","TOP2A","Chain2","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD", "IKZF2", "MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G")
      Missing <- setdiff(nms, names(df))  # Find names of missing columns
      df[Missing] <- 0                    # Add them, filled with '0's
      MainTcell_test <- df[nms]

      MainTcell_test[is.na(MainTcell_test)] <- 0
      set.seed(123)
      kmeans10 <- kmeans(MainTcell_test, centers = 500, nstart = 1)
      df_centres <- as.data.frame(kmeans10$centers)

      df_centres <- df_centres %>%
        mutate(CellTypist_list_CD4 = case_when(
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &PDCD1 >0 & ICOS >0 & CXCR5>0  ~"Follicular helper T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) & GZMK >0 & CD4>0 & IL10>0 ~ "Memory GZMK+IL10+CD4+ cytotoxic T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &GZMK >0 & CD4>0 & IL10>0 ~ "Memory GZMK+IL10-CD4+ cytotoxic T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &NKG7>0 & GNLY>0 & CD4>0 & Chain2==1 ~ "NK CD4+ ab T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &CTLA4>0 & IL2RA >0 & FOXP3>0 & CD4>0 ~ "Regulatory T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &KLRB1 >0 & SLC4A10 >0 & `TRAV1-2`>0  ~ "MAIT cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &CD4>0 & CCR7>0 & SELL>0 ~ "Tcm/Naive helper T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &KLRB1>0 & AQP3>0 & ITGB1>0 ~ "Tem/Effector helper T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &PDCD1>0 & CD4>0 & CTLA4>0 ~ "Tem/Effector helper T cells PD1+",
          (CD3D>0 | CD3E >0 |CD3G>0 ) &GZMK>0 & CD4>0 & CCL5>0 ~ "Tem/Trm cytotoxic CD4+ T cells",
          # c(CD3D>0 | CD3E >0 |CD3G>0 ) & CX3CR1>0 & GZMB>0 & GNLY>0 ~ "Tem/Temra cytotoxic T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &CD27>0 & CCR7>0 & IKZF2>0 ~ "Treg(diff)",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &ITGA1>0 & ITGAE>0 & CXCR6>0 ~"Trm cytotoxic T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &CCL5>0 & CXCR3>0 & TBX21>0 & CD4>0 ~ "Type 1 helper CD4+ ab T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &IL7R>0 & CCR6>0 & ZBTB16>0 ~ "Type 17 helper T cells",
          c(CD3D>0 | CD3E >0 |CD3G>0 ) &MIR155HG>0 & BIRC3>0 & SMS>0 ~ "T(agonist)",

          TRUE ~ NA_character_))

      df_centres$clust_name <- rownames(df_centres)
      head(df_centres)
      df_centres2 <- df_centres[names(df_centres) %in% c("clust_name","CellTypist_list_CD4")]
      Cluster_Df <- as.data.frame( kmeans10$cluster)
      Cluster_Df
      names(Cluster_Df) <- "clust_name"
      Cluster_Df$Cell_Index <- rownames(Cluster_Df)
      Cluster_Df_class <- merge(df_centres2,Cluster_Df,by="clust_name")
      Cluster_Df_class2 <- Cluster_Df_class[,names(Cluster_Df_class) %in% c("Cell_Index","CellTypist_list_CD4")]
      Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
    }
    else {
      Cluster_Df_class2 <- as.data.frame(MainTcell[,names(MainTcell) %in% c("Cell_Index")])
      names(Cluster_Df_class2) <- "Cell_Index"
      Cluster_Df_class2[order(Cluster_Df_class2$Cell_Index),]
    }


  })
  add.classification_CellTypist_cycling_df <- reactive({
    sc <- vals_meta.sc$metadata_SCobj
    validate(
      need(nrow(sc)>0,
           "Imput metadata required")
    )
    MainTcell <- add.classification()
    if (input$add.classification_CellTypist_cycling==T) {
      rownames(MainTcell) <- MainTcell$Cell_Index
      MainTcell_test <- MainTcell[,names(MainTcell) %in% c("MKI67","TOP2A","Chain2","CD3D","CD3E","CD3G","CLEC10A","GNLY","S100A9")]

      df <- MainTcell_test
      # Vector of columns you want in this data
      nms <- c("MKI67","TOP2A","Chain2","CD3D","CD3E","CD3G","CLEC10A","GNLY","S100A9")
      Missing <- setdiff(nms, names(df))  # Find names of missing columns
      df[Missing] <- 0                    # Add them, filled with '0's
      MainTcell_test <- df[nms]

      MainTcell_test[is.na(MainTcell_test)] <- 0
      set.seed(123)
      kmeans10 <- kmeans(MainTcell_test, centers = 100, nstart = 1)
      df_centres <- as.data.frame(kmeans10$centers)

      df_centres <- df_centres %>%
        mutate(CellTypist_list_Cell_cycling = case_when(
          MKI67 >0 & TOP2A > 0 & Chain2==1 ~ "Cycling ab T cells",
          MKI67 >0 & TOP2A > 0 & Chain2== -1 ~ "Cycling gd T cells",
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
    ### kmean class ----
    # # SELL = CD62L; CCR7; IL7R == CD127 ;IL2RA== CD25 ; ITGA1 == CD49a (not in dataset)
    # # https://www.thermofisher.com/be/en/home/life-science/cell-analysis/cell-analysis-learning-center/immunology-at-work/cytotoxic-t-cell-overview.html
  ## download the classification -----
  classification_SC <- reactive({
    annotation <- cbind(add.classification_T_cell_Df(),
                        add.classification_T.B_cell_Function_Df(),
                        add.classification_T_cell_Function_CD4_df(),
                        add.classification_T_cell_Function_CD8_df(),
                        add.classification_T_cell_Function_CD4_CD8pos_df(),
                        add.classification_T_cell_Function_CD4_CD8neg_df(),

                        add.classification_T_cell_Activation_df(),
                        add.classification_T_cell_Memory_df(),

                        add.classification_CellTypist_list_higher(),
                        add.classification_CellTypist_list_lower(),
                        add.classification_CellTypist_list_df(),
                        add.classification_CellTypist_list_CD8_df(),
                        add.classification_CellTypist_list_CD4_df(),
                        add.classification_CellTypist_cycling_df()
        )

    rownames(annotation) <- annotation$Cell_Index
    annotation <- annotation[,!grepl("Cell_Index",names(annotation))]
    annotation$Cell_Index <- rownames(annotation)

    sc <- add.CellTypist()

    md <- sc@meta.data

    md.anno <-  merge(md,annotation,by = "Cell_Index",sort=F,all.x=T)
    rownames(md.anno) <- md.anno$Cell_Index
    sc@meta.data <- md.anno
    sc
  })

  ## output ----

  output$DEx_table_TcellClass<- DT::renderDataTable( {
   sc <- classification_SC()
   datatable(sc@meta.data, extensions = "Buttons", options = list(searching = TRUE,
                                                                       ordering = TRUE,
                                                                       buttons = c('copy','csv', 'excel'),
                                                                       dom = 'Bfrtip',
                                                                       pageLength=10,
                                                                       lengthMenu=c(2,5,10,20,50,100),
                                                                       scrollX = TRUE
    ))
  }, server = FALSE)

  # output$DEx_table_TcellClass <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE),  options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 5, scrollX = TRUE),{
  #
  #
  #
  #   # calls <- calls[, !duplicated(names(calls))]
  #   # calls <- as.data.frame(calls)
  #
  # })

  # save Seurat object -----


  output$downloaddf_SeruatObj <- downloadHandler(
    filename = function(){
      if (input$save_type ==".h5") {

        paste("Seurat Obj ",gsub("-", ".", Sys.Date()),".h5Seurat", sep = "")
        }

      else {
        paste("Seurat Obj ",gsub("-", ".", Sys.Date()),".rds", sep = "")
    }
    },
    content = function(file){
      if (input$save_type ==".h5") {
        as.h5Seurat(classification_SC(),file)
      }

      else {
        saveRDS(classification_SC(),file)
      }

    } )

# Analysis -----
  ## uploading seruat obj ----
  input.data_sc_pro <- reactive({switch(input$dataset_sc_pro,"test_data_sc_pro" = test.data_sc_pro(),"own_data_sc_pro" = own.data_sc_pro())})
  test.data_sc_pro <- reactive({
    if (input$STEGO_R_pro=="STEGO_R (.h5Seurat)") {
        # dataframe = read.csv(system.file("extdata","clusTCR/cdr3.csv",package = "STEGO.R"))

      dataframe = LoadH5Seurat(system.file("extdata","BDrhap/Analysis/Seurat Obj 2023.02.28.h5Seurat",package = "STEGO.R"))
    }

    else {
      dataframe = readRDS("../KULSMC35CRC16highQualCD8.rds")
    }


  })

  own.data_sc_pro <- reactive({
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

  input.data_sc_clusTCR <- reactive({switch(input$dataset_sc_pro,"test_data_sc_pro" = test.data_sc_clusTCR(),"own_data_sc_pro" = own.data_sc_clusTCR())})
  # input.data_sc_clusTCR <- reactive({switch(input$dataset_cluster_file,"test_data_clusTCR" = test.data_sc_clusTCR(),"own_data_clusTCR" = own.data_sc_clusTCR())})
  test.data_sc_clusTCR <- reactive({
      dataframe = read.csv(system.file("extdata","BDrhap/clusTCR/ClusTCR2.csv",package = "STEGO.R"))

    # dataframe = read.csv("../Test.Cluster/Cluster_table  2022.12.15.csv")
  })
  own.data_sc_clusTCR <- reactive({
    inFile_cluster_file <- input$file_cluster_file
    if (is.null(inFile_cluster_file)) return(NULL)
    else {
      dataframe = read.csv(inFile_cluster_file$datapath)
    }

  })

  # upload TCRex file
  input.data_sc_TCRex <- reactive({switch(input$dataset_sc_pro,"test_data_sc_pro" = test.data_sc_TCRex(),"own_data_sc_pro" = own.data_sc_TCRex())})
  # input.data_sc_clusTCR <- reactive({switch(input$dataset_cluster_file,"test_data_clusTCR" = test.data_sc_clusTCR(),"own_data_clusTCR" = own.data_sc_clusTCR())})
  test.data_sc_TCRex <- reactive({

    dataframe = read.table(system.file("extdata","BDrhap/TCRex/tcrex_nsjhx29ivo.tsv",package = "STEGO.R"),skip=7,header = T,sep="\t")

    # read.table("../Public data/Bd Rhapsody/TCRex/tcrex_nsjhx29ivo.tsv",skip = 7,header = T,sep="\t")
  })

  own.data_sc_TCRex <- reactive({
    inupload_TCRex_file <- input$upload_TCRex_file
    if (is.null(inupload_TCRex_file)) return(NULL)
    else {
      dataframe = read.table(inupload_TCRex_file$datapath,skip = 7,header = T,sep="\t")
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

    if (input$STEGO_R_pro== "Seurat (.rds)") {
      updateSelectInput(
        session,
        "Samp_col",
        choices=names(df3.meta),
        selected = "orig.ident")
    }


    else {
      updateSelectInput(
        session,
        "Samp_col",
        choices=names(df3.meta),
        selected = "Sample_Name")
    }

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


    ### top clonotypes -----
  vals <- reactiveValues(top10=NULL)
  cols_top_clonal_plot <- reactive({
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


    plot <- ggplot(UMAP.wt.clonality,aes(x=UMAP_1,UMAP_2,colour=TYPE.clonality,alpha = TYPE.clonality,label=TYPE.clonality))+
      geom_point()+
      scale_color_manual(values = colorblind_vector, na.value=input$NA_col_analysis,labels = ~ stringr::str_wrap(.x, width = 20))+
      scale_alpha_manual(values = rep(1,length(unique(UMAP.wt.clonality$TYPE.clonality))), na.value=0.5,labels = ~ stringr::str_wrap(.x, width = 20))+
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

      plot+ facet_wrap(~UMAP.wt.clonality$ID_Column)
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
  df4
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

  Topclonotypes <- reactive({
    UMAP.wt.clonality2 <- For_col_top()
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


    plot <- ggplot(data=UMAP.wt.clonality2,aes(x=UMAP_1,UMAP_2,colour=topclones2,label =topclones2 ))+
      geom_point()+
      scale_color_manual(values = topclones_col$col, na.value=input$NA_col_analysis,labels = ~ stringr::str_wrap(.x, width = 20))+
      # scale_color_manual(values = colorblind_vector)+
      scale_alpha_manual(values = 1, na.value = 0.1,labels = ~ stringr::str_wrap(.x, width = 20))+
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

      plot+ facet_wrap(~UMAP.wt.clonality2$ID_Column)
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

  ## differental expression -----
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
      paste(input$name.file_clust," ClusTCR ",gsub("-", ".", Sys.Date()),".csv", sep = "")
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
      paste("TCR_Explore_markers_featurePlot_sc_",gsub("/", "-", x), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width=input$width_markers_featurePlot_sc,height=input$height_markers_featurePlot_sc, onefile = FALSE) # open the pdf device

      plot(markers_featurePlot())
      dev.off()}, contentType = "application/pdf" )

  output$downloadPlotPNG_markers_featurePlot_sc <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("TCR_Explore_markers_featurePlot_sc_", gsub("/", "-", x), ".png", sep = "")
    },
    content = function(file) {
      png(file, width = input$width_png_markers_featurePlot_sc,height = input$height_png_markers_featurePlot_sc,res = input$resolution_PNG_markers_featurePlot_sc)
      plot(markers_featurePlot())
      dev.off()},   contentType = "application/png")


    ### check Idents ----

  vals_meta.sc2 <- reactiveValues(AddIndents_SCobj=NULL)
  vals_meta.sc3 <- reactiveValues(AddIndents_SCobj2=NULL)

  observeEvent(input$run_update_clust,{
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

  output$UMAP_tb_download<- DT::renderDataTable( {
    sc <- input.data_sc_pro()
    validate(
      need(nrow(sc)>0,
           error_message_val_UMAP)
    )



    datatable(UMAP.TCRclonalit(), extensions = "Buttons", options = list(searching = TRUE,
                                                                       ordering = TRUE,
                                                                       buttons = c('copy','csv', 'excel'),
                                                                       dom = 'Bfrtip',
                                                                       pageLength=10,
                                                                       lengthMenu=c(2,5,10,20,50,100),
                                                                       scrollX = TRUE
    ))
  }, server = FALSE)


  output$Percent_tab <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
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

  output$test.files_classification <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{
    Add.UMAP.reduction()
  })


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
           error_message_val_UMAP)
    )

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

    ggplot(top_BD_cluster,aes(x=UMAP_1,UMAP_2,colour=Selected_function,alpha = Selected_function,label = Selected_function))+
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

  })

  output$UMAP_all_classification2 <- renderPlot({
    UMAP_all_classification()
  })

  output$downloadPlot_UMAP_all_classification  <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("top_clonotype_",gsub("/", "-", x), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width=input$width_UMAP_all_classification ,height=input$height_UMAP_all_classification , onefile = FALSE) # open the pdf device
      plot(UMAP_all_classification())
      dev.off()}, contentType = "application/pdf" )

  output$downloadPlotPNG_UMAP_all_classification  <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("top_clonotype_", gsub("/", "-", x), ".png", sep = "")
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
    # top_BD_cluster <-  Add.UMAP.reduction()
    # top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
    # top_BD_cluster$Selected_function[is.na(top_BD_cluster$Selected_function )] <- "unknown"
    # top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
    # top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function,levels = unique(top_BD_cluster$Selected_function))
    # num <- as.data.frame(unique(top_BD_cluster$Selected_function))
    # num
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

    top_BD_cluster %>%
      count(Selected_group,Selected_function) %>%
      group_by(Selected_group) %>%
      mutate(n = n/sum(n)) %>%
      ggplot(aes(x="", y=n, fill=Selected_function, group = Selected_group)) +
      geom_bar(stat="identity", width=1)+
      coord_polar("y", start=0)  +
      theme_void(20) +
      facet_wrap(~Selected_group) +
      theme(
        legend.key.size = unit(1, 'cm'))+
      # geom_col_pattern(aes(pattern=Selected_function,fill=Selected_function,pattern_fill = Selected_function),
      #                  # colour                   = 'black',
      #                  # pattern_aspect_ratio = 0.25,
      #                  pattern_density          = 0.5,
      #                  pattern_key_scale_factor = 0.75)
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

  # Select top clonotype ----

  observe({
    sc <- input.data_sc_pro()
    validate(
      need(nrow(sc)>0,
           error_message_val1)
    )
    meta.data <- as.data.frame(sc@meta.data)
    # meta.data <- meta.data[,names(meta.data) %in% c(grep("v",names(meta.data)), grep("cdr3",names(meta.data)))]
    updateSelectInput(
      session,
      "Gene_V_top",
      choices=names(meta.data),
      selected = "v_gene_cdr3_AG")

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
    BD_sum

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
      choices=BD_sum$cluster_name,
      selected = "TRGV9_ALWEVKELGKKIKV")
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

  output$Top_clonotype_df<- DT::renderDataTable( {
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
      paste("top_clonotype_",gsub("/", "-", x), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width=input$width_top_clonotype,height=input$height_top_clonotype, onefile = FALSE) # open the pdf device
      plot(ggplot_top_BD_clonotype_vs_SC())
      dev.off()}, contentType = "application/pdf" )

  output$downloadPlotPNG_top_clonotype <- downloadHandler(
    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("top_clonotype_", gsub("/", "-", x), ".png", sep = "")
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

    top_BD_cluster %>%
      count(Selected_group,Selected_function) %>%
      group_by(Selected_group) %>%
      mutate(n = n/sum(n)) %>%
      ggplot(aes(x="", y=n, fill=Selected_function, group = Selected_group)) +
      geom_bar(stat="identity", width=1)+
      coord_polar("y", start=0)  +
      theme_void(20) +
      facet_wrap(~Selected_group) +
      theme(
        legend.key.size = unit(1, 'cm'))+
      # geom_col_pattern(aes(pattern=Selected_function,fill=Selected_function,pattern_fill = Selected_function),
      #                  # colour                   = 'black',
      #                  # pattern_aspect_ratio = 0.25,
      #                  pattern_density          = 0.5,
      #                  pattern_key_scale_factor = 0.75)
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
  UMAP_chart_alpha_gamma()
}
  })

  output$downloadPlot_top_clonotype_pie <- downloadHandler(



    filename = function() {
      x <- gsub(":", ".", Sys.time())
      paste("top_clonotype_pie_",gsub("/", "-", x), ".pdf", sep = "")
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
      paste("top_clonotype_", gsub("/", "-", x), ".png", sep = "")
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
    top_BD_cluster <-  top_clonotype_bar_code()
    top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
    top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
    top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_names_top]
    top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% input$clust_group]
    df.col <- unlist(colors_cols_Top_pie_clonotype())

      ggplot(top_BD_cluster,aes(x=UMAP_1,UMAP_2,colour=Selected_function))+
        geom_point(aes(size = 5))+
        scale_color_manual(values = df.col, na.value=input$NA_col_analysis) +
        theme(
          legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
          legend.title = element_blank(),
          legend.position = input$legend_position,
             )+
        theme_bw()

  })

  ### Add in the TCR clustering to the seurat object ----

  #### Ridge plot (expression per T cell) ----

  observeEvent(input$run_string.data_Exp_top,{
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
      "string.data_Exp_top",
      choices=name.df$Gene_Name,
      selected = c("GZMB"))

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
            legend.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = input$legend_position)

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


output$Ridge_chart_alpha_gamma_stat_comp <- DT::renderDataTable(escape = FALSE, options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength =10, scrollX = TRUE),{

  sc <-input.data_sc_pro()

  validate(
    need(nrow(sc)>0,
         error_message_val_sc)
  )
  df = as.data.frame(sc[["RNA"]]@scale.data)
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
    mat_emp[i,1] <- t.test(MainTcell_other[,i], MainTcell_selected[,i])$p.value
    mat_emp[i,2] <- t.test(MainTcell_other[,i], MainTcell_selected[,i])$conf.int[1]
    mat_emp[i,3] <- t.test(MainTcell_other[,i], MainTcell_selected[,i])$conf.int[2]
    mat_emp[i,4] <- t.test(MainTcell_other[,i], MainTcell_selected[,i])$estimate[1]
    mat_emp[i,5] <- t.test(MainTcell_other[,i], MainTcell_selected[,i])$estimate[2]
    mat_emp[i,6] <- wilcox.test(MainTcell_other[,i], MainTcell_selected[,i])$p.value
  }

  mat_emp <- as.data.frame(mat_emp)
  rownames(mat_emp) <- names(MainTcell)[1:num]
  names(mat_emp) <- c("Pval","lowCI","upperCI","mean of other","mean of selected","Wilcoxon")

  mat_emp_sub <- subset(mat_emp,mat_emp$Pval!="NaN")

  mat_emp_sub$bonferroni <- p.adjust(mat_emp_sub$Pval, method = "bonferroni", n = length(mat_emp_sub$Pval))
  mat_emp_sub$BH_t.test <- p.adjust(mat_emp_sub$Pval, method = "BH", n = length(mat_emp_sub$Pval))
  mat_emp_sub$BH_Wil <- p.adjust(mat_emp_sub$Wilcoxon, method = "BH", n = length(mat_emp_sub$Pval))
  mat_emp_sub

  })

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
        legend.title = element_blank(),
        legend.position = input$legend_position,
      )

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

#### Epitope upload -----

  df_tcrex <- reactive({
    tcrex <- input.data_sc_TCRex()
  })

  output$MainTcell_Check <- DT::renderDataTable(escape = FALSE, filter = list(position = 'top', clear = FALSE), options = list(autoWidth = FALSE, lengthMenu = c(2,5,10,20,50,100), pageLength = 10, scrollX = TRUE),{

    df_tcrex()

  })

  output$downloaddf_TCRex_meta <- downloadHandler(
    filename = function(){
      paste(" TCRex_metadata ",gsub("-", ".", Sys.Date()),".tsv", sep = "")
    },
    content = function(file){
      df <- as.data.frame(df_tcrex())
      write.csv(df,file, row.names = F,quote = F)
    } )

  #### Epitope heatmap -----
  heatmap_epitope <- reactive({
    df <- input.data_sc_TCRex()
    dat <-df
    dat <- as.data.frame(dat)
    dat$cloneCount <- 1
    df <- as.data.frame(ddply(dat,(c("CDR3_beta","epitope","pathology","cloneCount")),numcolwise(sum)))

    df.1 <- acast(df, get(input$epitope_hm)~get(input$pathology_hm), value.var="cloneCount")
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
      top_BD_cluster$CDR3_beta <- top_BD_cluster$cdr3_BD
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

      top_BD_cluster$CDR3_beta <- top_BD_cluster$cdr3_BD

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
      top_BD_cluster$CDR3_beta <- top_BD_cluster$cdr3_BD
    }
    epi$beta <- epi$CDR3_beta
    df3.meta <- merge(df3.meta,epi,by="CDR3_beta",all.x=T)
    df3.meta$selected <- df3.meta[,names(df3.meta) %in% input$epitope_umap_selected]
    df3.meta <- df3.meta[order(df3.meta$selected,decreasing = F),]
    df3.meta$selected <- factor(df3.meta$selected,levels = unique(df3.meta$selected))
    num <- as.data.frame(unique(df3.meta$selected))
    num <- as.data.frame(num[complete.cases(num)==T,])
    palette_rainbow <- unlist(colors_cols_epitope())

    ggplot(data=df3.meta,aes(x=UMAP_1,UMAP_2,colour=selected,size= selected))+
      geom_point()+
      scale_color_manual(na.value=input$NA_col_analysis, values = palette_rainbow)+
      scale_size_manual(na.value=0.25,values = rep(3,dim(num)[1]))+
      theme(
        legend.text = element_text(colour="black", size=12,family=input$font_type),
        legend.title = element_blank(),
        legend.position = input$legend_position,
      )+
      theme_bw()


    # ggplot(top_BD_cluster,aes(x=UMAP_1,UMAP_2,colour=Selected_function))+
    #   geom_point()+
    #   scale_color_manual(na.value=, values = df.col)+
    #   theme(
    #     legend.text = element_text(colour="black", size=input$Bar_legend_size,family=input$font_type),
    #     legend.title = element_blank(),
    #     legend.position = ,
    #   )+
    #   theme_bw()
  })
  output$UMAP_Epitope_plot <- renderPlot({
    UMAP_Epitope()
  })


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
          legend.text = element_text(colour="black", size=24,family="sans"),
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
          legend.text = element_text(colour="black", size=24,family="sans"),
          legend.title = element_blank(),
          legend.position = "right",
        )
    }

    figure

  })

  output$UMAP_ClusTCR2_plot <- renderPlot({
    UMAP_ClusTCR2()
  })
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
          legend.text = element_text(colour="black", size=24,family="sans"),
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
          legend.text = element_text(colour="black", size=24,family="sans"),
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

### end -----
}
shinyApp(ui, server)
}
