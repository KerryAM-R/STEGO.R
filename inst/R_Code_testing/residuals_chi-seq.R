source("R/global.R")
sc <- LoadH5Seurat("~/Documents/Post-doc/PostDoc_Antwerp/R_files/1_Files_scTCRseq/Antwerp/data_mmr_dextramer/M6_2023.03.17.h5Seurat")

reduction <- (sc@reductions$umap)
UMAP <- as.data.frame(reduction@cell.embeddings)
UMAP$Cell_Index <- rownames(UMAP)
meta.data <- as.data.frame(sc@meta.data)
meta.data$Cell_Index <- rownames(meta.data)
meta.data$Motif_gene <- meta.data[,names(meta.data) %in% "vdj_gene_cdr3_AG_BD"]
meta.data$cloneCount <- 1
BD <- meta.data[,names(meta.data) %in% c("Motif_gene","cloneCount")]
BD_sum <- ddply(BD,names(BD)[-c(2)] ,numcolwise(sum))
BD_sum2 <- BD_sum[!BD_sum$Motif_gene %in% "NA",]
head(BD_sum2)
names(BD_sum)[2] <- "Total_count"

BD_sum$frequency <- BD_sum$Total_count/sum(BD_sum$Total_count)
UMAP_md2<- BD_sum[order(BD_sum$frequency,decreasing = T),]
# BD_sum <- subset(BD_sum,BD_sum$cluster_name!="NA")
BD_sum

# chisq
unique(meta.data$orig.ident)

names(meta.data)
totals <- meta.data[,names(meta.data) %in% c("orig.ident","classify.Adaptive.cell")]
totals <- meta.data[,names(meta.data) %in% c("orig.ident","CellTypist_list_lower")]
totals <- meta.data[,names(meta.data) %in% c("orig.ident","CellTypist_list_curated_T_cell")]

# totals <- meta.data[,names(meta.data) %in% c("orig.ident","CellTypist_list_higher")]
# totals <- totals[!totals$orig.ident %in% c( "M6_d7"  , "M6_d21" ),]

names(totals) <- c("groups","Function")
names(totals)
unique(totals$Function)
totals <- totals[!totals$Function %in% c( "NA","Pre-pro-B cells","B cell","Monocyte precursor","DC","ILC1","Transitional B cells",
                                          "Naive B cells","DC1" ,"Mono-mac"),]

totals <- totals[!totals$Function %in% c( "NA","DN T cells" ,"NK CD8+ T cells","Tcm/Naive cytotoxic T cells","DN Tregs cell","DN GNLY+ T cell",
                                          "Tem/Trm cytotoxic CD8+ T cells","Gamma-delta T cells","CD8+GNLY+ T cell" ,"DP thymocytes" ,
                                          "DN thymocytes" ,"Gamma-delta T cells (CCL5+)","KLRB1+ MAIT cells (TRAV1-2+)","Tem/Temra cytotoxic T cells"),]
unique(totals$Function)



# totals <- totals[!totals$groups %in% c( "M6_d7" ),]

totals$groups <- factor(totals$groups,levels = c("M6_d0", "M6_d3" , "M6_d7" ,"M6_d21" ,"M6_d150", "M6_d365"))
totals$groups <- factor(totals$groups,levels = c("M6_d0", "M6_d3" , "M6_d7" ,"M6_d21" ,"M6_d150", "M6_d365"))
tb_totals <- table(totals$groups,totals$Function)
tb_totals
chisq <- chisq.test(tb_totals)
chisq
chisq$observed
chisq$expected
round(chisq$residuals, 3) # pearson residual is the standardised distances
library(corrplot)
# ?corrplot

corrplot(chisq$residuals, is.cor = FALSE,hclust.method ="ward.D2",tl.col = "black")
corrplot(chisq$residuals, is.cor = FALSE, hclust.method ="ward.D2",tl.col = "black",col = topo.colors(10)) # COL2('PRGn')

contrib <- 100*chisq$residuals^2/chisq$statistic
round(contrib, 3)

library(RColorBrewer)

corrplot(contrib, is.cor = FALSE,tl.col = "black",col = topo.colors(100) )
chisq$p.value
chisq$statistic

