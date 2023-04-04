#' Packages needed for STEGO.R
#' @param Load_required_packages load packages
#' @export

Load_required_packages <- function() {
  # options(shiny.maxRequestSize = 20000*1024^2)
  options(shiny.maxRequestSize = 100000*1024^2)
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
  suppressMessages(require("SeuratObject"))
  suppressMessages(require("shiny"))
  suppressMessages(require("shinyBS"))
  suppressMessages(require("shinybusy"))
  suppressMessages(require("shinyWidgets"))
  suppressMessages(require("shinyFiles"))
  suppressMessages(require("showtext"))
  suppressMessages(require("showtextdb"))
  suppressMessages(require("stringr"))
  suppressMessages(require("sysfonts"))
  suppressMessages(require("tibble"))
  suppressMessages(require("tidyr"))
  suppressMessages(require("VLF"))

}
