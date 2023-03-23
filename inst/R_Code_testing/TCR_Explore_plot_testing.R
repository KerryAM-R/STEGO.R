sc <- LoadH5Seurat("~/Documents/Post-doc/PostDoc_Antwerp/R_files/1_Files_scTCRseq/Antwerp/data_mmr_dextramer/M6_2023.03.20.h5Seurat")

md <- sc@meta.data

write.csv(table(md$Sample_Name,md$orig.ident),"TCR_coverate.csv")

require(scGate)

reduction <- (sc@reductions$umap)
UMAP <- as.data.frame(reduction@cell.embeddings)
UMAP$Cell_Index <- rownames(UMAP)
meta.data <- as.data.frame(sc@meta.data)
meta.data$Cell_Index <- rownames(meta.data)

names(meta.data)




FeaturePlot(sc, features = c("CD4","PTPRC","CD3E","CD8A"), min.cutoff = 'q10', label = T)
FeaturePlot(sc, features = c("CD4"), min.cutoff = 'q1', label = T,split.by = "orig.ident",ncol = 4)


meta.data$Motif_gene <- meta.data[,names(meta.data) %in% "vdj_gene_cdr3_AG_BD"]

meta.data$Chain <- ifelse(grepl("TRBV",meta.data$Motif_gene) |
                            grepl("TRAV",meta.data$Motif_gene),"abTCR",
                          ifelse(grepl("TRGV",meta.data$Motif_gene) |
                                   grepl("TRDV",meta.data$Motif_gene),"gdTCR","unknown"))

meta.data$umap.meta.classification.chain <- paste(meta.data$Chain,meta.data$classify.T.cell)
#
meta.data$umap.meta.classification.chain <- ifelse(grepl("unknown", meta.data$umap.meta.classification.chain),"unknown",meta.data$umap.meta.classification.chain)
meta.data$umap.meta.classification.chain <- ifelse(grepl(" NA", meta.data$umap.meta.classification.chain),"",meta.data$umap.meta.classification.chain)

head(meta.data)

meta.data <- meta.data[order(meta.data$Cell_Index,decreasing = F),]

umap.meta <- merge(UMAP,meta.data,by="Cell_Index")

umap.meta


df3.meta <- umap.meta
df3.meta$Selected_group <- df3.meta[,names(df3.meta) %in% "orig.ident"]
unique(df3.meta$Selected_group)



df3.meta$Selected_group <- factor( df3.meta$Selected_group,levels = c("M6_d0" ,   "M6_d3"  , "M6_d7" , "M6_d21" ,"M6_d150" ,"M6_d365" ))



md <- df3.meta[,names(df3.meta) %in% c("orig.ident","vdj_gene_cdr3_AG_BD") ]

# upset.parameters <- function () {
#   df <- input.data2();
#
#   validate(
#     need(nrow(df)>0,
#          error_message_val1)
#   )


  df <- as.data.frame(md)
  head(df)
  unique.df <- df[,names(df) %in% c(input$Samp_col,input$V_gene_Class)]

  unique.df <- unique(df[c("orig.ident","vdj_gene_cdr3_AG_BD")])
  unique.df
  names(unique.df) <- c("group","chain")
  unique.df <- unique.df[unique.df$chain!="NA",]


  unique.df$group <- factor( unique.df$group,levels = c("M6_d0" ,   "M6_d3"  , "M6_d7" , "M6_d21" ,"M6_d150" ,"M6_d365" ))
  unique.df$cloneCount <- 1
  mat <- acast(unique.df, chain~group, value.var="cloneCount")
  mat[is.na(mat)] <- 0
  mat <- as.data.frame(mat)
  # a <- as.data.frame(unique(names(mat)))
  #   a$V1 <- distinctColorPalette(dim(a)[1])
 head( mat)

  df.x <- make_comb_mat(mat)

  list.names <- as.character(unique(  unique.df$group))
  list.names

  if (input$upset_anno == "Colour by degree") {

    your_list <- c(input$upset_colours_list1)
    your_list_df <- as.data.frame((unlist(strsplit(your_list, ','))))
    names(your_list_df) <- "ID"

    ht = draw(UpSet(df.x,
                    pt_size = unit(input$upset.point.size, "mm"),
                    lwd = input$upset.lwd,
                    row_names_gp =  gpar(fontfamily = input$font_type,fontsize = input$upset.text.size),#changes font size of "set size" labels
                    column_names_gp = gpar(fontfamily = input$font_type),
                    comb_col = c(your_list_df$ID)[comb_degree(df.x)],
                    top_annotation = upset_top_annotation(df.x,
                                                          numbers_gp = gpar(fontfamily = input$font_type,fontsize = input$font.size.anno.upset),
                                                          annotation_name_gp = gpar(fontfamily = input$font_type,fontsize=input$font.size.anno.upset),
                                                          gp = gpar(fill = input$top_annotation_colour),

                    ),
                    right_annotation = upset_right_annotation(df.x,
                                                              add_numbers = T,
                                                              numbers_gp = gpar(fontfamily = input$font_type,fontsize = input$font.size.anno.upset),
                                                              annotation_name_gp = gpar(fontfamily = input$font_type,fontsize=input$font.size.anno.upset),
                                                              gp = gpar(fill = input$right_annotation_colour),
                    ),


                    set_order  = c(input$order.of.group)
    ), padding = unit(c(20, 20, 20, 20), "mm"))


  }
  ?upset_top_annotation

    ht = draw(UpSet(df.x,
                    pt_size = unit(5, "mm"),
                    lwd = 1,
                    row_names_gp =  gpar(fontfamily = "sans", fontsize = 12),
                    column_names_gp = gpar(fontfamily = "sans"),
                    top_annotation = upset_top_annotation(df.x,
                                                          add_numbers = T,
                                                          annotation_name_gp = gpar(fontfamily = "sans"),
                                                          gp = gpar(fill = "black"),
                    ),


                    right_annotation = upset_right_annotation(df.x,
                                                              add_numbers = T,
                                                              numbers_gp = gpar(fontfamily = "sans",fontsize =12),
                                                              annotation_name_gp = gpar(fontfamily = "sans",fontsize=12),
                                                              gp = gpar(fill = "black"),
                    ),


                    set_order  =  list.names

    ), padding = unit(c(20, 20, 20, 20), "mm"))

    }










TCR_Explore <- read.csv("~/Documents/Post-doc/PostDoc_Antwerp/R_files/1_Files_scTCRseq/Antwerp/data_mmr_dextramer/TCR_Explore_overlap_3time.csv")
# CD4 <- read.csv("~/Desktop/TCR_CD4.csv")
# TCR_Explore <- subset(TCR_Explore,TCR_Explore$Type=="CD4")
head(TCR_Explore)
unique(TCR_Explore$Row.Labels)
TCR_Explore_2  <- TCR_Explore
# TCR_Explore_2 <- TCR_Explore[grepl("clone  4",TCR_Explore$clonotype_name),]
TCR_Explore_2
TCR_Explore_overlap <- merge(df3.meta,TCR_Explore_2,by.x="Motif_gene",by.y="Row.Labels")
# TCR_Explore_overlap <- merge(df3.meta,CD4,by.x="Motif_gene",by.y="vdj_gene_cdr3_AG_BD")
TCR_Explore_overlap

TCR_Explore_overlap$Samps <- paste(TCR_Explore_overlap$clonotype_name," (Time points = ",TCR_Explore_overlap$Time.points,"; Total = ",TCR_Explore_overlap$Sum,")",sep = "")
unique(TCR_Explore_overlap$Samps)

top_BD_cluster <-  TCR_Explore_overlap
top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% "Samps"]
# top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% "vdj_gene_cdr3_AG_BD"]
top_BD_cluster$Selected_function

top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
num= unique(top_BD_cluster$Selected_function)
num


top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% "orig.ident"]

top_BD_cluster$Selected_group <- factor( top_BD_cluster$Selected_group,levels = c("M6_d0" ,   "M6_d3"  , "M6_d7" , "M6_d21" ,"M6_d150" ,"M6_d365" ))
unique(top_BD_cluster$Selected_function)
num_leng <- length(unique(top_BD_cluster$Selected_function))
# ?facet_grid
# df.col <- unlist(colors_cols_Top_pie_clonotype())


df3.meta$Selected_function <- NA
df3.meta$Selected_function <- factor(df3.meta$Selected_function, levels = unique(TCR_Explore_overlap$Samps))
df3.meta$Selected_function

top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function, levels = unique(TCR_Explore_overlap$Samps))
top_BD_cluster$Selected_function
df3.meta$Selected_function
df3.meta_lim <- subset(df3.meta,df3.meta$UMAP_1> -5)
df3.meta_lim
palette_rainbow <- rainbow(length(num))

set.seed(1000)
palette_rainbow <- distinctColorPalette(length(num))



ggplot()+
  geom_point(data = df3.meta_lim,aes(x=UMAP_1,UMAP_2,colour=Selected_function,size=Selected_function,shape=Selected_function))+
  geom_point(data = top_BD_cluster,aes(x=UMAP_1,UMAP_2,colour=Selected_function,size=Selected_function,shape=Selected_function))+
  scale_color_manual(values = palette_rainbow, na.value="grey90") +
  scale_size_manual(values = rep(4,num_leng),na.value = 0.55) +
  scale_shape_manual(values = c(19,17,15,16,22,37,10,14),na.value = 1) +
  theme_bw()+
  theme(
    legend.text = element_text(colour="black", size=12,family="sans"),
    legend.title = element_blank(),
    legend.position = "right",
  )+
  # scale_size() +
  facet_wrap(~Selected_group,nrow = 3)

names(top_BD_cluster)
top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% "CellTypist_list_CD4_lower_curated"]
unique(top_BD_cluster$Selected_function )



top_BD_cluster$Selected_function <- gsub("Tem/Temra cytotoxic T cells","Non-CD4+ T cell",top_BD_cluster$Selected_function)


unique(top_BD_cluster$Selected_function)

top_BD_cluster$Selected_function <- factor(top_BD_cluster$Selected_function, levels = unique(top_BD_cluster$CellTypist_list_CD4_lower_curated))


num= unique(top_BD_cluster$Selected_function)
palette_rainbow <- terrain.colors(length(num))
top_BD_cluster %>%
  count(clonotype_name,Selected_function) %>%
  group_by(clonotype_name) %>%
  mutate(n = n/sum(n)) %>%
  ggplot(aes(x="", y=n, fill=Selected_function, group = clonotype_name)) +
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)  +
  theme_void(20) +
  facet_wrap(~clonotype_name,nrow=4) +
  theme(
    legend.key.size = unit(1, 'cm'))+
  # geom_col_pattern(aes(pattern=Selected_function,fill=Selected_function,pattern_fill = Selected_function),
  #                  # colour                   = 'black',
  #                  # pattern_aspect_ratio = 0.25,
  #                  pattern_density          = 0.5,
  #                  pattern_key_scale_factor = 0.75)
  scale_fill_manual(values = palette_rainbow, na.value = "grey90") +

theme(strip.text = element_text(size = 12, family = "sans"),
        legend.text = element_text(size = 12, family = "sans"),
        legend.title = element_blank()

  )

