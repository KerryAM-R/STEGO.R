#' @export
#'
# compatibility with Immunarch ----
suppressMessages(require("shiny"))
suppressMessages(require("shinyBS"))
suppressMessages(require("shinyWidgets"))
suppressMessages(require("Seurat"))
suppressMessages(require("plyr"))
suppressMessages(require("dplyr"))
suppressMessages(require("bslib"))
suppressMessages(require("Matrix"))
suppressMessages(require("DT"))
suppressMessages(require("shinybusy"))
suppressMessages(require("tidyverse"))
suppressMessages(require("SeuratDisk"))
suppressMessages(suppressWarnings(require("SeuratData")))
suppressMessages(require("igraph"))
suppressMessages(require("linkcomm")) # create the network graphs.
suppressMessages(require("ggpattern"))
# suppressMessages(require("plotly"))
suppressMessages(require("ggrepel"))
suppressMessages(require("showtext"))
suppressMessages(require("reshape2")) # acast function
suppressMessages(require("GGally"))
suppressMessages(require("ggnet"))
suppressMessages(require("network"))
suppressMessages(require("VLF")) ## aa.count.function
suppressMessages(require("motifStack")) # function
suppressMessages(require("ggseqlogo")) # create figure
suppressMessages(require("RIRA")) # cell typist
suppressMessages(require("colourpicker")) # select visual colour
suppressMessages(require("RColorBrewer"))
suppressMessages(require("randomcoloR"))
suppressMessages(require("ggridges")) # ridges plot
suppressMessages(require("fpc")) #
suppressMessages(require("ComplexHeatmap"))
suppressMessages(require("circlize")) # colorRamp2
suppressMessages(suppressWarnings(require("ClusTCR2")))
suppressMessages(require("doParallel"))
suppressMessages(require("RColorBrewer"))
suppressMessages(require("stringr"))
suppressMessages(require(gridExtra))

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
