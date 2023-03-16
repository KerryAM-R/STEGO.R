source("R/global.R")
sc <- LoadH5Seurat("~/Documents/Post-doc/PostDoc_Antwerp/R_files/1_Files_scTCRseq/Antwerp/data_mmr_dextramer/M62023.03.14.h5Seurat")
sc <- LoadH5Seurat("~/Documents/Post-doc/PostDoc_Antwerp/R_files/1_Files_scTCRseq/Antwerp/data_mmr_dextramer/M62023.03.15.h5Seurat")

FeaturePlot(sc, features = "FCGR3A", min.cutoff = 'q10', label = T)


df= as.data.frame(sc[["RNA"]]@scale.data)
head(df)[1:6]
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


names(df3.meta)

df3.meta$cluster_name <- df3.meta[,names(df3.meta) %in% "cdr3_BD"]

top_BD_clonotype <- df3.meta[df3.meta$cluster_name %in% "CASSLKDGRQPQHF",]
top_BD_clonotype


top_BD_cluster <-  top_BD_clonotype
top_BD_cluster$Selected_function <- top_BD_cluster[,names(top_BD_cluster) %in% "classify.Adaptive.cell"]
top_BD_cluster <- top_BD_cluster[order(top_BD_cluster$Selected_function),]
num= unique(top_BD_cluster$Selected_function)
num
palette_rainbow <- rainbow(length(num))

top_BD_cluster$Selected_group <- top_BD_cluster[,names(top_BD_cluster) %in% "orig.ident"]
# df.col <- unlist(colors_cols_Top_pie_clonotype())
df3.meta$Selected_function <- NA
ggplot()+
  geom_point(data = df3.meta,aes(x=UMAP_1,UMAP_2,colour=Selected_function))+
  geom_point(data = top_BD_cluster,aes(x=UMAP_1,UMAP_2,colour=Selected_function))+
  scale_color_manual(values = palette_rainbow, na.value="grey90") +
  # scale_size_manual(values = rep(input$size_selected_top,dim(top_BD_cluster)[1]),na.value = 1) +
  theme_bw()+
  theme(
    legend.text = element_text(colour="black", size=16,family="sans"),
    legend.title = element_blank(),
    legend.position = "right",
  )+
  scale_size(guide = 'none') +
  facet_grid(~orig.ident)
