#' Generate Default Fill Colors for ggplot2 Using Hues
#' @name gg_fill_hue
#' @description
#' This function generates a sequence of default fill colors suitable for use in ggplot2 plots by varying hues.
#'
#' @param n Number of colors to generate.
#'
#' @return A vector of colors in HCL space.
#'
#' @export
#'
gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Update Seurat Object with Additional Labels and Embeddings
#' @name update_sc
#' @description
#' This function adds additional labels and embeddings to a Seurat object if needed.
#' It merges the new IDs with the existing metadata based on a specified column name.
#' This process is essential for incorporating additional information into the analysis.
#'
#' @importFrom plyr ddply mutate
#' @import Seurat
#'
#' @param sc Merged and annotated Seurat object.
#' @param add_additional_lables Boolean indicating whether to add additional labels. Default is TRUE.
#' @param updateID File with new IDs to be merged with metadata.
#' @param colName Column name for merging metadata with updated IDs. Default is "Sample_Name".
#' @return Updated Seurat object with additional labels and embeddings.
#'
#' @export

update_sc <- function (sc = sc,
                       add_additional_lables = TRUE,
                       updateID = updateID,
                       colName = "Sample_Name") {
  sc <- sc
  sc@meta.data$order <- 1:dim(sc@meta.data)[1]
  reduction <- (sc@reductions$umap)
  UMAP <- as.data.frame(reduction@cell.embeddings)
  names(UMAP)[1:2] <- toupper(names(UMAP)[1:2])
  UMAP$Cell_Index <- rownames(UMAP)
  meta.data <- as.data.frame(sc@meta.data)
  umap.meta <- merge(UMAP,meta.data,by="Cell_Index",sort = FALSE)

  if(add_additional_lables) {
    labs <- updateID
    names(umap.meta)[names(umap.meta) %in% colName] <- "ID"
    umap.meta <- merge(labs,umap.meta,by="ID", sort = FALSE)
    names(umap.meta)[names(umap.meta) %in% "ID"] <- colName
  }

  rownames(umap.meta) <- umap.meta$Cell_Index
  umap.meta <- umap.meta[order(umap.meta$order),]
  sc@meta.data <- umap.meta
  sc

}



#' Summarize TCR-seq Data Based on scRepertoire Categories
#' @name TCR_Expanded
#' @description
#' This function summarizes TCR-seq data based on scRepertoire categories, providing counts, frequencies, and percentages of expanded clones.
#'
#' @param sc Merged and annotated Seurat object.
#' @param Samp_col Sample column name. Default is "Sample_Name".
#' @param V_gene_sc Column name containing VDJ gene and CDR3 sequence. Default is "vdj_gene_cdr3_AG_BD".
#' @param font_type Font type for plot. Default is "Times New Roman".
#'
#' @return A data frame summarizing TCR-seq data based on scRepertoire categories.
#'
#' @import plyr
#' @import extrafont
#' @return Adding the expanded values to the meta-data and returns the seurat object
#' @export

TCR_Expanded <- function (sc = sc,
                          Samp_col = "Sample_Name",
                          V_gene_sc = "vdj_gene_cdr3_AG_BD",
                          font_type = "Times New Roman"
) {

  sc <- sc
  df3.meta <- sc@meta.data
  df3.meta2 <- df3.meta[, names(df3.meta) %in% c(Samp_col, V_gene_sc)]
  df3.meta2$ID_Column <- df3.meta2[, names(df3.meta2) %in% Samp_col]
  df3.meta2$v_gene_selected <- df3.meta2[, names(df3.meta2) %in% V_gene_sc]

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
      frequency <= 1e-5 ~ "1. Very rare (0 > X \u2264 1e-5)",
      frequency <= 1e-4 ~ "2. Rare (1e-5 > X \u2264 1e-4)",
      frequency <= 0.001 ~ "3. Small (1e-4 > X \u2264 0.001)",
      frequency <= 0.01 ~ "4. Medium (0.001 > X \u2264 0.01)",
      frequency <= 0.10 ~ "5. Large (0.01 > X \u2264 0.1)",
      frequency <= 0.50 ~ "6. Hyperexpanded (0.1 > X \u2264 0.5)",
      frequency <= 1 ~ "7. Hyperexpanded (0.5 > X \u2264 1)",
      TRUE ~ "Other"
    ))

  df4 <- df4 %>%
    mutate(Number_expanded = case_when(
      samp.count <= 1 ~ "1. Single (X = 1)",
      samp.count <= 5 ~ "2. Small (1 < X \u2264 5)",
      samp.count <= 20 ~ "3. Medium (5 < X \u2264 20)",
      samp.count <= 100 ~ "4. Large (20 < X \u2264 100)",
      samp.count <= 500 ~ "5. Hyperexpanded (100 < X \u2264 500)",
      samp.count > 500 ~ "6. Hyperexpanded (X > 500)",
      TRUE ~ "6. Hyperexpanded (>500)"
    ))
  df4
}
