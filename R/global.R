#' @export
#'
# compatibility with Immunarch ----
options(shiny.maxRequestSize = 20000*1024^2)
suppressMessages(require("bslib"))
suppressMessages(require("circlize")) # colorRamp2
suppressMessages(suppressWarnings(require("ClusTCR2")))
suppressMessages(require("colourpicker")) # select visual colour
suppressMessages(require("ComplexHeatmap"))
suppressMessages(library("corrplot"))
suppressMessages(require("doParallel"))
suppressMessages(require("dplyr"))
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
suppressMessages(require("motifStack")) # function
suppressMessages(require("network"))
suppressMessages(require("plyr"))
suppressMessages(require("purrr"))
suppressMessages(require("randomcoloR"))
suppressMessages(require("RColorBrewer"))
suppressMessages(require("Rcpp"))
suppressMessages(require("readr"))
suppressMessages(require("reshape2")) # acast function
suppressMessages(require("scGate"))
suppressMessages(require("Seurat"))
suppressMessages(require("SeuratDisk"))
suppressMessages(require("SeuratData"))
suppressMessages(require("SeuratObject"))
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
suppressMessages(require("VLF")) ## aa.count.function

font_add_google("Gochi Hand", "gochi")
font_add_google("Schoolbell", "bell")
font_add_google("Press Start 2P", "Game")

showtext_auto()

# font ------
font <- as.data.frame(font_families())
font
names(font) <- "Fonts"

# set file upload limit

# reticulate::use_condaenv("/Users/kerrymullan/miniconda", required = TRUE) # changes to the location that the python conda packages are installed.

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
} # colouring function

# default colouring
gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

utils::globalVariables(c("%>%","CDR3_beta","Cell_Index","ClusTCR","Clust_size_order","CreateSeuratObject","DimHeatmap","DimPlot","ElbowPlot","FeaturePlot","FindAllMarkers","FindClusters","FindMarkers","FindNeighbors","FindVariableFeatures","Heatmap","ID_Column","Idents","Indiv","LabelPoints","LoadH5Seurat","Motif_from_cluster_file","NormalizeData","PercentageFeatureSet","Pie_ClusTCR2","Ridge_chart_alpha_gamma_stat_comp","RunHarmony","RunPCA","RunUMAP","ScaleData","Selected_chain3","Selected_function","Selected_group","TRBV_gene","TYPE.clonality","TukeyHSD","UMAP_1","UMAP_2","UpSet","V1","V2","VariableFeaturePlot","VariableFeatures","VlnPlot","acast","across","actionButton","add_busy_spinner","aes","all_of","alpha","aov","arrange","as.h5Seurat","bs_theme","case_when","checkboxInput","cloneCount","colorRamp2","colourInput","column","complete.cases","conditionalPanel","coord_polar","cores_ClusTCR2","count","datatable","ddply","dev.off","distinctColorPalette","div","downloadButton","downloadHandler","draw","element_blank","element_line","element_text","everything","facet_wrap","fileInput","fluidPage","fluidRow","frequency","geom_bar","geom_density_ridges","geom_hline","geom_jitter","geom_point","geom_text_repel","geom_violin","geom_vline","get_scGateDB","ggplot","ggtitle","gpar","group","group_by","guide_legend","guides","h4","h5","hcl","hcl.colors","head","heat.colors","incProgress","kmeans","labs","mainPanel","make_comb_mat","mcl_cluster","merge.names","motif_plot","mutate","n","nFeature_RNA","na.omit","navbarMenu","navbarPage","need","netplot","numcolwise","numericInput","observe","observeEvent","p","p.adjust","pdf","percent.mt","percent.rb","plotOutput","png","rainbow","reactive","reactiveValues","read.csv","read.table","renderPlot","renderPrint","renderUI","runif","scGate","scale_alpha_manual","scale_color_manual","scale_colour_manual","scale_fill_manual","scale_size","scale_size_manual","select","selectInput","selected","selected_top_clonotype","seurat_clusters","shinyApp","sidebarLayout","sidebarPanel","t.test","tabPanel","tabsetPanel","tags","terrain.colors","textInput","theme","theme_bw","theme_ridges","theme_void","top_n","topclones2","topo.colors","uiOutput","unit","updateSelectInput","upset_right_annotation","upset_top_annotation","v_gene_selected","validate","verbatimTextOutput","wellPanel","where","wilcox.test","withProgress","write.csv","write.table","write_csv","Idents<-"))
