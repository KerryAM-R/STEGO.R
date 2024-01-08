#' Packages needed for STEGO.R
#' @param Load_required_packages load packages
#' @export

Load_required_packages <- function() {
  # options(shiny.maxRequestSize = 20000*1024^2)
  options(shiny.maxRequestSize = 100000*1024^2)

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
  suppressMessages(require("scales"))
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
}

#' Packages needed for STEGO.R
#' @param gg_fill_hue coloring function
#' @export

gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Packages needed for STEGO.R
#' @param TCR_Expanded_fun Feature set
#' @param sc seurat object description
#' @param Samp_col Name of the sample in the meta data
#' @param V_gene_sc
#' @export

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


#' Packages needed for STEGO.R
#' @param mat_sum Feature set
#' @param sc seurat object description
#' @param Samp_col Name of the sample in the meta data
#' @param V_gene_sc V gene
#' @export

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
