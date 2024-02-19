
options(shiny.maxRequestSize = 100000*1024^2)

suppressMessages(require("bslib"))
suppressMessages(require("circlize")) # colorRamp2
# suppressMessages(suppressWarnings(require("ClusTCR2")))
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
suppressMessages(require("magrittr"))
# suppressMessages(require("motifStack")) # function
suppressMessages(require("network"))
suppressMessages(require("plyr"))
suppressMessages(require("purrr")) # merging Seurat objects
suppressMessages(require("randomcoloR"))
suppressMessages(require("RColorBrewer"))
suppressMessages(require("Rcpp"))
suppressMessages(require("readr"))
suppressMessages(require("reshape2")) # acast function
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


utils::globalVariables(c("CDR3_beta","Cell_Index","ClusTCR","Clust_size_order","CreateSeuratObject","DimHeatmap","DimPlot","ElbowPlot","FeaturePlot","FindAllMarkers","FindClusters","FindMarkers","FindNeighbors","FindVariableFeatures","Heatmap","ID_Column","Idents","Indiv","LabelPoints","LoadH5Seurat","Motif_from_cluster_file","NormalizeData","PercentageFeatureSet","Pie_ClusTCR2","Ridge_chart_alpha_gamma_stat_comp","RunHarmony","RunPCA","RunUMAP","ScaleData","Selected_chain3","Selected_function","Selected_group","TRBV_gene","TYPE.clonality","TukeyHSD","UMAP_1","UMAP_2","UpSet","V1","V2","VariableFeaturePlot","VariableFeatures","VlnPlot","acast","across","actionButton","add_busy_spinner","aes","all_of","alpha","aov","arrange","as.h5Seurat","bs_theme","case_when","checkboxInput","cloneCount","colorRamp2","colourInput","column","complete.cases","conditionalPanel","coord_polar","cores_ClusTCR2","count","datatable","ddply","dev.off","distinctColorPalette","div","downloadButton","downloadHandler","draw","element_blank","element_line","element_text","everything","facet_wrap","fileInput","fluidPage","fluidRow","frequency","geom_bar","geom_density_ridges","geom_hline","geom_jitter","geom_point","geom_text_repel","geom_violin","geom_vline","get_scGateDB","ggplot","ggtitle","gpar","group","group_by","guide_legend","guides","h4","h5","hcl","hcl.colors","head","heat.colors","incProgress","kmeans","labs","mainPanel","make_comb_mat","mcl_cluster","merge.names","motif_plot","mutate","n","nFeature_RNA","na.omit","navbarMenu","navbarPage","need","netplot","numcolwise","numericInput","observe","observeEvent","p","p.adjust","pdf","percent.mt","percent.rb","plotOutput","png","rainbow","reactive","reactiveValues","read.csv","read.table","renderPlot","renderPrint","renderUI","runif","scGate","scale_alpha_manual","scale_color_manual","scale_colour_manual","scale_fill_manual","scale_size","scale_size_manual","select","selectInput","selected","selected_top_clonotype","seurat_clusters","shinyApp","sidebarLayout","sidebarPanel","t.test","tabPanel","tabsetPanel","tags","terrain.colors","textInput","theme","theme_bw","theme_ridges","theme_void","top_n","topclones2","topo.colors","uiOutput","unit","updateSelectInput","upset_right_annotation","upset_top_annotation","v_gene_selected","validate","verbatimTextOutput","wellPanel","where","wilcox.test","withProgress","write.csv","write.table","write_csv","Idents<-","Read10X_h5","corrplot","input.data.barcode.10x","input.data.features.10x","input.data.matrix.10x","read_h5","renderText","shinyDirButton","shinyDirChoose","write_h5","x","Var1","Var2","custom_db_scGATE","error_message_val1","error_message_val2","error_message_val3","error_message_val4","error_message_val_10x_barcode","error_message_val_10x_features","error_message_val_UMAP","error_message_val_sc","font","gg_fill_hue","h6","test_fun","SaveH5Seurat","empty","netplot_ClusTCR2","Percentage","melt","residuals","scale_fill_gradient2","scale_size_area","setNames","y"))
