sc <- LoadH5Seurat("~/Documents/Post-doc/PostDoc_Antwerp/R_files/1_Files_scTCRseq/Antwerp/data_mmr_dextramer/M62023.03.14.h5Seurat")

reduction <- (sc@reductions$umap)
UMAP <- as.data.frame(reduction@cell.embeddings)
UMAP$Cell_Index <- rownames(UMAP)
meta.data <- as.data.frame(sc@meta.data)
meta.data$Cell_Index <- rownames(meta.data)

names(meta.data)

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

TCR_Explore <- read.csv("~/Documents/Post-doc/PostDoc_Antwerp/R_files/1_Files_scTCRseq/Antwerp/data_mmr_dextramer/TCR_Explore_overlap.csv")
head(TCR_Explore)
# TCR_Explore <- TCR_Explore[1:2,]

TCR_Explore <- TCR_Explore[grepl("TRAV26-2.TRAJ54",TCR_Explore$Row.Labels),]


TCR_Explore
TCR_Explore_overlap <- merge(df3.meta,TCR_Explore,by.x="Motif_gene",by.y="Row.Labels")
TCR_Explore_overlap$Motif_gene <- paste(TCR_Explore_overlap$Motif_gene," (Samples = ",TCR_Explore_overlap$Time.points,"; Total = ",TCR_Explore_overlap$Sum,")",sep = "")


top_BD_cluster <-  TCR_Explore_overlap

names(top_BD_cluster)
top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% "CellTypist_list_lower"]
top_BD_cluster$Selected_function

top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
num= unique(top_BD_cluster$Selected_function)
num
palette_rainbow <- rainbow(length(num))

top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% "orig.ident"]

top_BD_cluster$Selected_group <- factor( top_BD_cluster$Selected_group,levels = c("M6_d0" ,   "M6_d3"  , "M6_d7" , "M6_d21" ,"M6_d150" ,"M6_d365" ))
num_leng <- length(unique(top_BD_cluster$Selected_function))
# ?facet_grid
# df.col <- unlist(colors_cols_Top_pie_clonotype())
df3.meta$Selected_function <- NA
set.seed(123)
ggplot()+
  geom_point(data = df3.meta,aes(x=UMAP_1,UMAP_2,colour=Selected_function,size = Selected_function))+
  geom_point(data = top_BD_cluster,aes(x=UMAP_1,UMAP_2,colour=Selected_function,size = Selected_function))+
  scale_color_manual(values = palette_rainbow, na.value="grey80") +
  scale_size_manual(values = rep(2,num_leng),na.value = 0.55,guide = 'none') +
  theme_bw()+
  theme(
    legend.text = element_text(colour="black", size=12,family="sans"),
    legend.title = element_blank(),
    legend.position = "right",
  )+
  # scale_size() +
  facet_wrap(~Selected_group)
