


gg_fill_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Merging RDS Seurat objects
#' @name Updating_IDs
#' @description
#' This function will add in the additional labels required for the analysis if needed. this process is also required to add in the reductions to the meta.data
#'
#' @param sc Add the merged and annotated file
#' @param update_sc If you have added in the additional labels if required
#' @param updateID file with the new id's based on '[colName]
#' @param colName for the merging of the meta data with the updated ID, and will default to 'Sample_Name'
#' @export

update_sc <- function (sc = sc,
                       add_additional_lables = T,
                       updateID = updateID,
                       colName = "Sample_Name") {
  sc <- sc
  sc@meta.data$order <- 1:dim(sc@meta.data)[1]
  reduction <- (sc@reductions$umap)
  UMAP <- as.data.frame(reduction@cell.embeddings)
  names(UMAP)[1:2] <- toupper(names(UMAP)[1:2])
  UMAP$Cell_Index <- rownames(UMAP)
  meta.data <- as.data.frame(sc@meta.data)
  umap.meta <- merge(UMAP,meta.data,by="Cell_Index",sort = F)

  if(add_additional_lables) {
    labs <- updateID
    names(umap.meta)[names(umap.meta) %in% colName] <- "ID"
    umap.meta <- merge(labs,umap.meta,by="ID", sort = F)
    names(umap.meta)[names(umap.meta) %in% "ID"] <- colName
  }

  rownames(umap.meta) <- umap.meta$Cell_Index
  umap.meta <- umap.meta[order(umap.meta$order),]
  sc@meta.data <- umap.meta
  sc

}

#' TCR-seq summary table for both group and total expansion
#' @name Summary_TCR_table
#' @description
#' This function will add in the additional labels required for the analysis if needed. this process is also required to add in the reductions to the meta.data
#'
#' @param sc Add the merged and annotated file
#' @param Samp_col Sample column name with the default being "Sample_Name"
#' @param V_gene_sc The vdj and CDR3 sequence "vdj_gene_cdr3_AG_BD"
#' @param save_file save file to the Clonotypes folder
#' @export
#'

Summary_TCR_table <- function (sc = sc,
                               Samp_col = "Sample_Name",
                               V_gene_sc = "vdj_gene_cdr3_AG_BD",
                               save_file = T
) {
  library(magrittr)
  library(dplyr)
  sc <- sc
  df <- sc@meta.data
  df <- as.data.frame(df)
  unique.df <- unique(df[,names(df) %in% c( Samp_col, V_gene_sc) ])
  unique.df <- unique.df %>%
    select( Samp_col, everything())

  # count data of unique ------
  names(unique.df) <- c("group","chain")
  unique.df <- subset(unique.df,unique.df$chain != "NA")
  unique.df <- subset(unique.df,unique.df$group != "NA")
  unique.df$cloneCount <- 1
  mat <- acast(unique.df, chain~group, value.var="cloneCount")
  mat[is.na(mat)] <- 0
  Count_data <- as.data.frame(rowSums(mat))
  names(Count_data) <- "V1"

  # Total count ------
  unique.df <- (df[,names(df) %in% c(Samp_col, V_gene_sc) ])
  unique.df <- unique.df %>%
    select(all_of(Samp_col), everything())

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
  # print(as.data.frame(mat))
  if(save_file) {
    newfolder <- "Clonotypes"
    newpath <- file.path(newfolder)

    if (dir.exists(newpath) != T) {
      dir.create(file.path(dirname(path), newfolder))
    }

    write.csv(mat,"Clonotypes/TCR_summary_file.csv")

  }
  mat
}

#' Column factor ID function
#' @name ID_Column_factor_function
#' @description
#' This section is to add the factor order for the graphs
#'
#' @param sc Add the merged and annotated file
#' @param Samp_col Sample column name with the default being "Sample_Name"
#' @export
#'

ID_Column_factor_function <- function (sc = sc, Samp_col = "Sample_Name") {
  sc <- sc
  df <- sc@meta.data
  df2 <- as.data.frame(unique(df[names(df) %in% Samp_col]))
  df2 <- as.data.frame(df2)
  names(df2) <- "V1"
  df2 <- as.data.frame(df2[order(df2$V1), ])
  names(df2) <- "V1"
  df3 <- subset(df2, df2$V1 != "NA")
  df3

}

#' Creating the summary table in the form of an UPSET plot
#' @name Upset_plot_multi
#' @description
#' This section is to add the factor order for the graphs
#'
#' @param sc Add the merged and annotated file
#' @param Samp_col Sample column name with the default being "Sample_Name"
#' @param V_gene_sc column name vdj_gene_cdr3_AG_BD
#' @param font_type from the extrafont package; default: Times New Roman
#' @param save_plots save file to director or you can store as an object
#' @param width_px Width of plot 2000
#' @param height_px Height of plot 1200
#' @param resolution_px Width of plot 144
#' @export
#'

Upset_plot_multi <- function (sc = sc,
                              Samp_col = "Sample_Name",
                              V_gene_sc = "vdj_gene_cdr3_AG_BD",
                              font_type = "Times New Roman",
                              save_plots = T,
                              width_px = 2000,
                              height_px = 1200,
                              resolution_px = 144
) {
  require(extrafont)
  require(ComplexHeatmap)
  require(dplyr)
  require(reshape2)

  newfolder <- "Clonotypes"
  newpath <- file.path(newfolder)

  if (dir.exists(newpath) != T) {
    dir.create(file.path(dirname(path), newfolder))
  }

  Summary_TCR_table(sc,save_file = T)

  sc <- sc
  df <- sc@meta.data

  df <- as.data.frame(df)
  unique.df <- unique(df[,names(df) %in% c(Samp_col, V_gene_sc) ])
  unique.df <- unique.df %>%
    select(Samp_col, everything())
  names(unique.df) <- c("group","chain")
  ID_factor <- ID_Column_factor_function(sc)
  unique.df <- unique.df[unique.df$group %in%  ID_factor$V1,]
  unique.df <- subset(unique.df,unique.df$chain != "NA")
  unique.df <- subset(unique.df,unique.df$group != "NA")
  unique.df$cloneCount <- 1
  mat <- acast(unique.df, chain~group, value.var="cloneCount")
  mat[is.na(mat)] <- 0
  mat <- as.data.frame(mat)
  df.x <- make_comb_mat(mat)
  list.names <- as.character(ID_factor$V1)
  ht = draw(UpSet(df.x,
                  pt_size = unit(5, "mm"),
                  lwd = 1,
                  row_names_gp =  gpar(fontfamily = font_type, fontsize = 12),
                  column_names_gp = gpar(fontfamily = font_type,fontsize =12),
                  top_annotation = upset_top_annotation(df.x,
                                                        add_numbers = T,
                                                        numbers_gp = gpar(fontfamily = font_type,fontsize =12),
                                                        annotation_name_gp = gpar(fontfamily = font_type, fontsize = 12),
                                                        gp = gpar(fill = "black"),
                  ),
                  right_annotation = upset_right_annotation(df.x,
                                                            add_numbers = T,
                                                            numbers_gp = gpar(fontfamily = font_type,fontsize =12),
                                                            annotation_name_gp = gpar(fontfamily = font_type,fontsize=12),
                                                            gp = gpar(fill = "black"),
                  ),
                  set_order  =  list.names

  ), padding = unit(c(20, 20, 20, 20), "mm"))

  if(save_plots) {
    top.name.clonotypes.count_png <- paste("Clonotypes/Overlap_Upset_plot.png", sep = "")
    png(top.name.clonotypes.count_png, width = width_px, height = height_px, res = resolution_px)
    plot(ht)
    dev.off()

  }

  ht

}

#' Function for summarising the TCR-seq based on scRepertiore categories
#' @name TCR_Expanded
#' @description
#' This section is to add the factor order for the graphs
#'
#' @param sc Add the merged and annotated file
#' @param Samp_col Sample column name with the default being "Sample_Name"
#' @param V_gene_sc column name vdj_gene_cdr3_AG_BD
#' @param font_type from the extrafont package; default: Times New Roman
#' @import plyr
#' @import extrafont
#' @export
#'

TCR_Expanded <- function (sc = sc,
                          Samp_col = "Sample_Name",
                          V_gene_sc = "vdj_gene_cdr3_AG_BD",
                          font_type = "Times New Roman"
) {
  require(plyr)
  require(extrafont)
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
}

#' Frequency distribution or count
#' @name clonal_plot_multi
#' @description
#' This section is to add the factor order for the graphs
#'
#' @param sc Add the merged and annotated file
#' @param Samp_col Sample column name with the default being "Sample_Name"
#' @param V_gene_sc column name vdj_gene_cdr3_AG_BD
#' @param colourtype colour types include: default, hcl.colors, topo.colors, heat.colors, terrain.colors, rainbow, random
#' @param NA_col_analysis If NA's are present it will make them coloured as grey90
#' @param font_type from the extrafont package; default: Times New Roman
#' @param title_size Size of the graph plot in 20
#' @param text_size Size of the text 12
#' @param Legend_size Size of the text 12
#' @param height_px Height of plot 1200
#' @param resolution_px Width of plot 144
#' @param save_file save file to director or you can store as an object
#' @import ggplot2
#' @export
#'

clonal_plot_multi <- function (sc = sc,
                               Samp_col = "Sample_Name",
                               V_gene_sc = "vdj_gene_cdr3_AG_BD",
                               colourtype = "default",
                               Graph_type_bar = "Frequency_expanded",
                               NA_col_analysis = "grey90",
                               font_type = "Times New Roman",
                               legend_position = "right",
                               title_size = 20,
                               text_size = 12,
                               Legend_size = 12,
                               height_px = 1200,
                               resolution_px = 144,
                               save_file = F
) {

  require(ggplot2)
  sc <- sc

  # this is already loaded and
  df4 <- TCR_Expanded(sc, Samp_col = Samp_col,
                      V_gene_sc = V_gene_sc,
                      font_type = font_type)

  ID_Column_factor <- ID_Column_factor_function(sc)
  head(ID_Column_factor)
  df4 <- df4[df4$ID_Column %in% ID_Column_factor$V1,]
  df4$ID_Column <- as.character(df4$ID_Column)
  df4$ID_Column <- factor(df4$ID_Column,levels = ID_Column_factor$V1)
  df4$graphType <- df4[,names(df4) %in% Graph_type_bar]
  df4 <- df4[order(df4$graphType),]
  col.df <- as.data.frame(unique(df4$graphType))
  print(head(col.df))
  names(col.df) <- "V1"
  num <- length(col.df$V1)

  if (colourtype == "default") {
    colorblind_vector <- c(gg_fill_hue(num))
  } else if (colourtype == "hcl.colors") {
    colorblind_vector <- c(hcl.colors(num, palette = "viridis"))
  } else if (colourtype == "topo.colors") {
    colorblind_vector <- c(topo.colors(num))
  } else if (colourtype == "heat.colors") {
    colorblind_vector <- c(heat.colors(num))
  } else if (colourtype == "terrain.colors") {
    colorblind_vector <- c(terrain.colors(num))
  } else if (colourtype == "rainbow") {
    colorblind_vector <- c(rainbow(num))
  } else if (colourtype == "random") {
    colorblind_vector <- distinctColorPalette(num)

  }  else {

  }

  col.df$col <- colorblind_vector
  print(head(col.df))
  print(head(df4))

  pl_bar <-  ggplot(df4,aes(x=ID_Column,y=frequency,fill = graphType,label = graphType))+
    geom_bar(stat="identity")+
    theme_bw() +
    scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 50),values=alpha(col.df$col, 1),na.value = NA_col_analysis) +
    scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 50),values = alpha(col.df$col, 1)) +
    theme(
      axis.title.y = element_text(colour="black",family=font_type,size = title_size),
      axis.text.y = element_text(colour="black",family=font_type,size = text_size),
      axis.text.x = element_text(colour="black",family=font_type,size = text_size,angle=90),
      axis.title.x = element_blank(),
      legend.text = element_text(colour="black", size=Legend_size,family=font_type),
      legend.title = element_blank(),
      legend.position = legend_position,
    )

  if(save_file) {
    num_indiv <- length(unique(df4$ID_Column))
    top.name.clonotypes.count_png <- paste("Clonotypes/",Graph_type_bar,"_barplot.png",sep = "")
    png(top.name.clonotypes.count_png, width = (num_indiv * 100 + 600), height = height_px, res = resolution_px)
    plot(pl_bar)
    dev.off()
  }
  pl_bar
}

#' Selected clones function for the automated process
#' @name selected_clonotypes
#' @description
#' This section is to add the factor order for the graphs
#'
#' @param sc Add the merged and annotated file
#' @param Samp_col Sample column name with the default being "Sample_Name"
#' @param V_gene_sc column name vdj_gene_cdr3_AG_BD
#' @param colourtype colour types include: default, hcl.colors, topo.colors, heat.colors, terrain.colors, rainbow, random
#' @param NA_col_analysis If NA's are present it will make them coloured as grey90
#' @param font_type from the extrafont package; default: Times New Roman
#' @param title_size Size of the graph plot in 20
#' @param text_size Size of the text 12
#' @param Legend_size Size of the text 12
#' @param height_px Height of plot 1200
#' @param resolution_px Width of plot 144
#' @param save_file save file to director or you can store as an object
#' @import ggplot2
#' @export
#'


selected_clonotypes <- function (sc = sc, restrict_to_expanded = F) {
  sc <- sc

  BD_sum <- as.data.frame(Summary_TCR_table(sc))
  BD_sum <- subset(BD_sum,BD_sum$CloneTotal>1)
  BD_sum <- subset(BD_sum,BD_sum$TotalSamps>1)
  BD_sum$TotalSamps

  BD_sum$priority <- 1/(BD_sum$TotalSamps * BD_sum$CloneTotal)

  if(restrict_to_expanded) {
    BD_sum$same <- ifelse(BD_sum$CloneTotal == BD_sum$TotalSamps,"NEx","Ex")
    BD_sum <- subset(BD_sum,BD_sum$same == "Ex")
    print(head(BD_sum))
    BD_sum <- BD_sum[,names(BD_sum) != "same"]
  }

  BD_sum$cluster_name <- rownames(BD_sum)
  BD_sum

}

#' Automated function for the PublicLike clonotypes
#' @name Clonotypes_PublicLike
#' @description
#' This section is to add the factor order for the graphs
#' @param sc seurat file
#' @param restrict_to_expanded Restrict the analysis to include only expanded. removed 1 count in every individual.
#' @param Samp_col Sample column name with the default being "Sample_Name"
#' @param V_gene_sc column name vdj_gene_cdr3_AG_BD
#' @param species_analysis Different types of specifices of either human (hs) or mouse (mm)
#' @param colourtype colour types include: default, hcl.colors, topo.colors, heat.colors, terrain.colors, rainbow, random
#' @param Split_group_by_ split the bar graph by this variable. Defaulted to Tcellfunction
#' @param NA_col_analysis If NA's are present it will make them coloured as grey90
#' @param font_type from the extrafont package; default: Times New Roman
#' @param LogFC_ logFC variable for the findmarkers function
#' @param min_point_ Minimum point value
#' @param pval.ex.filter Find Marker P-value cut-off
#' @param low.dotplot lower colour of the dot plot
#' @param middle.dotplot middle colour of the dot plot
#' @param high.dotplot high colour of the dot plot
#' @param gene_set_min_count minumim number of genes for each gene set
#' @param p.val_cutoff_top p-value cut-off for the over-representation nalysis
#' @param title_size Size of the graph plot in 20
#' @param text_size Size of the text 12
#' @param Legend_size Size of the text 12
#' @param legend_position location of the legend: right, left, top, bottom or none
#' @param height_px Height of plot 1200
#' @param resolution_px Width of plot 144
#' @param save_file save file to director or you can store as an object
#' @import ggplot2
#' @export
#'

Clonotypes_PublicLike <- function (sc = sc,
                                   restrict_to_expanded = T,
                                   Samp_col = "Sample_Name",
                                   V_gene_sc = "vdj_gene_cdr3_AG_BD",
                                   species_analysis = "hs",
                                   colourtype = "default",
                                   Split_group_by_ = "Tcellfunction",
                                   cutoff_priority = 1,
                                   font_type = "Times New Roman",
                                   NA_col_analysis = "grey90",
                                   min_point_ = 0.25,
                                   LogFC_ = 0.25,
                                   pval.ex.filter = 0.05,
                                   low.dotplot = "darkblue", middle.dotplot = "white", high.dotplot = "darkred",
                                   gene_set_min_count = 1,
                                   p.val_cutoff_top = 1,
                                   title_size = 20,
                                   text_size = 12,
                                   Legend_size = 12,
                                   legend_position = "right",
                                   height_px = 1200,
                                   resolution_px = 144
) {
  require(Seurat)
  require(scales)
  sc <- sc

  newfolder = "Clonotypes/Multi/Private"
  newpath <- file.path(newfolder)

  if (dir.exists(newpath) != T) {
    newfolder <- "Clonotypes/Multi"
    newpath <- file.path(newfolder)

    newfolder <- "Clonotypes/Multi/PublicLike"
    dir.create(file.path(newfolder))

    newfolder <- "Clonotypes/Multi/Private"
    dir.create(file.path(newfolder))
  }

  df3.meta <- sc@meta.data
  df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% V_gene_sc]
  sc@meta.data$Vgene <- sc@meta.data[,names(sc@meta.data) %in% V_gene_sc]

  BD_sum <- selected_clonotypes(sc,restrict_to_expanded =  restrict_to_expanded)

  print(dim(BD_sum))
  if (dim(BD_sum)[1]>0) {
    BD_sum$obs <- 1
    BD_sum <-  subset(BD_sum,BD_sum$priority<cutoff_priority)

    top.name.clonotypes.top_png <- paste(newpath,"/Selected_clones.csv",sep="")
    write.csv(BD_sum,top.name.clonotypes.top_png)
    observations <- sum(BD_sum$obs)

    for (i in 1:observations) {
      message(1/observations, detail = paste("Clone", i,"of",observations))

      ##### bar graph -------
      message(BD_sum$cluster_name[i])
      name.clone <- BD_sum$cluster_name[i]
      top_BD_clonotype <- df3.meta[df3.meta$cluster_name %in% name.clone,]
      # print(top_BD_clonotype)

      dtop_clonotype_bar_code <- top_BD_clonotype

      dtop_clonotype_bar_code$Selected_group <- dtop_clonotype_bar_code[,names(dtop_clonotype_bar_code) %in% Split_group_by_]
      num <- 1
      # num <- as.data.frame(num[complete.cases(num)==T,])
      as.data.frame(length(num))
      if (colourtype == "default") {
        colorblind_vector <- gg_fill_hue(num)
      } else if (colourtype == "hcl.colors") {
        colorblind_vector <- c(hcl.colors(num, palette = "viridis"))
      } else if (colourtype == "topo.colors") {
        colorblind_vector <- c(topo.colors(num))
      } else if (colourtype == "heat.colors") {
        colorblind_vector <- c(heat.colors(num))
      } else if (colourtype == "terrain.colors") {
        colorblind_vector <- c(terrain.colors(num))
      } else if (colourtype == "rainbow") {
        colorblind_vector <- c(rainbow((num)))
      } else if (colourtype == "random") {
        colorblind_vector <- distinctColorPalette(num)

      }  else {

      }


      colorblind_vector <- as.data.frame(colorblind_vector)
      names(colorblind_vector) <- "cols"

      dtop_clonotype_bar_code$Selected_chain2 <- dtop_clonotype_bar_code[,names(dtop_clonotype_bar_code) %in% V_gene_sc]
      dtop_clonotype_bar_code$Selected_chain3 <- gsub("_"," ",dtop_clonotype_bar_code$Selected_chain2)
      dtop_clonotype_bar_code$Selected_chain3 <- gsub("[.]"," ",dtop_clonotype_bar_code$Selected_chain3)

      ggplot_plot <- ggplot(dtop_clonotype_bar_code, aes(x=Selected_group, fill=Selected_chain3,colour = Selected_chain3, label = Selected_chain3)) +
        geom_bar() +
        theme_bw()+
        scale_color_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector$cols, na.value = NA_col_analysis)+
        scale_fill_manual(labels = ~ stringr::str_wrap(.x, width = 10), values = colorblind_vector$cols, na.value = NA_col_analysis)+
        # scale_alpha_manual(values = rep(1,length(unique(dtop_clonotype_bar_code$Selected_chain))), na.value=0.5)+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        theme(
          axis.title.y = element_text(colour="black",family=font_type,size = title_size),
          axis.text.y = element_text(colour="black",family=font_type,size = text_size),
          axis.text.x = element_text(colour="black",family=font_type,size = text_size,angle=90),
          axis.title.x = element_blank(),
          legend.text = element_text(colour="black", size=Legend_size,family=font_type),
          legend.position = legend_position,
          legend.title = element_blank()
        ) +
        guides(color = "none", size = "none")



      newfolder = "Clonotypes/Multi/Private/"
      top.name.clonotypes.top_png <- paste(newfolder,i,"_top_clone_",gsub("[/]","",gsub("&","",name.clone)),".png",sep="")

      num_width <- length(unique(dtop_clonotype_bar_code$Selected_group))

      png(top.name.clonotypes.top_png, width = (num_width*100+500),height = height_px,res = resolution_px)
      plot(ggplot_plot)
      dev.off()

      ###### Find marker loop ------
      name.clone <- BD_sum$cluster_name[i]

      message(paste0("Downloading Dom stats files and dot plot...",name.clone))

      sc@meta.data$Gene_select <- ifelse(sc@meta.data$Vgene %in% name.clone,name.clone,"other")
      sc@meta.data
      Idents(object = sc) <- sc@meta.data$Gene_select

      min.pct.expression<- min_point_
      min.logfc<-  LogFC_

      cluster.names <- unique(Idents(sc))[order(unique(Idents(sc)))]
      markers.fm.list <- FindMarkers(sc, ident.1 = name.clone, min.pct = min.pct.expression,  logfc.threshold = min.logfc, only.pos=TRUE)
      markers.fm.list2 <- subset(markers.fm.list,markers.fm.list$p_val_adj < pval.ex.filter)
      message(paste(length(markers.fm.list2$p_val_adj),"total markers for cluster",i))

      ###### dotplot ------
      if(length(markers.fm.list2$p_val_adj)>0) {

        newfolder = "Clonotypes/Multi/Private/"
        clonotype.name.stats <- paste(newfolder,i,"_",gsub("[/]","",gsub("&","",name.clone)),"_stats_table.csv", sep = "")
        write.csv(markers.fm.list2,clonotype.name.stats,row.names = T)

        groupLen <- length(unique(sc@meta.data$Gene_select))

        message(paste0("Saved csv ",name.clone))
        list.names <- rownames(markers.fm.list2)

        if (length(rownames(markers.fm.list2))>40) {
          list.names <- list.names[1:40]
        }

        else {
          list.names <- rownames(markers.fm.list2)
        }

        size_legend = Legend_size-2

        plotdotplot <- DotPlot(sc, features = list.names) +
          RotatedAxis() +
          theme(
            axis.title.y = element_blank(),
            axis.text.y = element_text(colour="black",family=font_type,size = text_size),
            axis.text.x = element_text(colour="black",family=font_type,size = text_size, angle = 90),
            axis.title.x = element_blank(),
            legend.title = element_text(colour="black", size=Legend_size,family=font_type),
            legend.text = element_text(colour="black", size=size_legend,family=font_type),
            legend.position = legend_position,
          ) +
          scale_colour_gradient2(low = low.dotplot, mid = middle.dotplot, high = high.dotplot) +
          scale_x_discrete(labels = label_wrap(20)) +
          scale_y_discrete(labels = label_wrap(20))

        newfolder = "Clonotypes/Multi/Private/"
        file.name.clone <- paste(newfolder,i,"_",gsub("[/]","",gsub("&","",name.clone)),"_dotplot_plot.png", sep = "")

        ### download the dot plot -------
        png(file.name.clone, width = 2400, height = 600,res = 144)
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

        # if(datasource == "BD_Rhapsody_Paired" || datasource == "BD_Rhapsody_AIRR") { # selectInput("datasource", "Data source",choices=c("10x_Genomics","BD_Rhapsody_Paired","BD_Rhapsody_AIRR")),
        #   geneSet$GeneSet <- gsub("-",".",geneSet$GeneSet)
        # }

        if(species_analysis == "mm") {
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
        message(paste(i,"has", length(geneSet2$in.geneset),"GeneSets in Multi"))

        if(length(geneSet2$in.geneset)>0) {
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
          geneSet2 <- subset(geneSet2,geneSet2$in.geneset >= gene_set_min_count)
          geneSet2 <- subset(geneSet2,geneSet2$p.val <= p.val_cutoff_top)
          geneSet2$FDR <- p.adjust(geneSet2$p.val, method = "fdr")
          geneSet2$Bonferroni <- p.adjust(geneSet2$p.val, method = "bonferroni")

          newfolder = "Clonotypes/Multi/Private/"
          top.name.overrep <- paste(newfolder,i,"_",gsub("[/]","",gsub("&","",name.clone)),"_OverRep.csv", sep = "")
          write.csv(geneSet2,top.name.overrep, row.names = F)

        }
      }
    }
  } else {

  }
}
