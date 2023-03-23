sc <- LoadH5Seurat("~/Documents/Post-doc/PostDoc_Antwerp/R_files/1_Files_scTCRseq/Antwerp/data_mmr_dextramer/M6_2023.03.16.h5Seurat")
df= as.data.frame(sc[["RNA"]]@scale.data)
MainTcell <- as.data.frame(t(df))
# names(MainTcell) <- ifelse(grepl("pAbO",names(MainTcell)),names(MainTcell),toupper(names(MainTcell)))
MainTcell$Cell_Index <- rownames(MainTcell)

# names(MainTcell) <- gsub("CD8B.","CD8B",names(MainTcell))
md <- sc@meta.data
names(md)
md$Motif_gene <- md[,names(md) %in% "vj_gene_cdr3_AG_BD"]
unique(md$Motif_gene)

ifelse(grepl("TRBV",md$Motif_gene),1,0)

md$Chain2
md$Chain2 <- ifelse(grepl("_._", md$Motif_gene),0,
                    ifelse(md$Motif_gene=="_",0,
                           ifelse(md$Motif_gene=="",0,
                                  ifelse(md$Motif_gene==" & ",0,
                                         ifelse(grepl("TRBV",md$Motif_gene) |
                                                  grepl("TRAV",md$Motif_gene),1,
                                                ifelse(grepl("TRGV",md$Motif_gene) |
                                                         grepl("TRDV",md$Motif_gene),-1,0))))))


md$CD4_pos <- ifelse(md$seurat_clusters==0,2,
                     ifelse( md$seurat_clusters==1,2,
                             ifelse(md$seurat_clusters==2,2,-1)))

md$CD4_pos
md2 <- md[names(md) %in% c("Cell_Index","Chain2","CD4_pos")]
MainTcell <- merge(md2,MainTcell,by="Cell_Index",all.y=T)



rownames(MainTcell) <- MainTcell$Cell_Index

FeaturePlot(sc, features = "CD8A", min.cutoff = 'q10', label = T)
FeaturePlot(sc, features = "CD4", min.cutoff = 'q10', label = T)

MainTcell_test <- MainTcell[,names(MainTcell) %in% c("MKI67","TOP2A","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD","MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","S100A9","EBI3","CCL19","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","FCGR3A","S100A13","TLE1","AREG","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","C1QA","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1","Chain2","HLA-DRA", "HLA-DRB5","HLA-DRB1","CD4_pos")]
# -----
df <- MainTcell_test
# Vector of columns you want in this data
nms <- c("MKI67","TOP2A","PDCD1","ZNF683","CD8A","ICOS","CXCR5","GZMK","CD4","IL10","NKG7","GNLY","CTLA4","IL2RA","FOXP3","KLRB1","SLC4A10","TRAV1-2","CCR7","SELL","AQP3","ITGB1","CX3CR1","GZMB","CCL5","CD27","IKZF2","ITGA1","ITGAE","CXCR6","CXCR3","TBX21","IL7R","CCR6","ZBTB16","ITGAD","MIR155HG","BIRC3","SMS","CD3D","CD3E","CD3G","FCRL2","ITGAX","CD79A","MS4A1","CD19","TNFRSF13B","CD22","POU2AF1","SUGCT","CR2","IGHM","IGHD","TCL1A","AICDA","CD24","MYO1C","MME","ZCCHC7","RAG1","DNTT","GLL1","IGLL5","CD40","CD1C","FCER1A","CLEC10A","BATF3","CADM1","CLEC9A","S100A9","EBI3","CCL19","KLF4","AXL","IRF8","PCLAF","FXYD2","HES1","CD99","CD1A","SMPD3","ACY3","SPINK2","CD34","FCGR3A","S100A13","TLE1","AREG","IKZF3","GATA3","KLRG1","HPGDS","IL4I1","RORC","KIT","XCL2","TYROBP","C1QC","HMOX1","LYZ","VCAN","CD14","S100A12","FCN1","C1QA","CD160","TOX2","SATB1","CCR9","GNG4","TRDC","TRGC1","JCHAIN","MZB1","XBP1","Chain2","HLA-DRA", "HLA-DRB5","HLA-DRB1")

Missing <- setdiff(nms, names(df))  # Find names of missing columns
df[Missing] <- 0                    # Add them, filled with '0's
MainTcell_test <- df[nms]
MainTcell_test[is.na(MainTcell_test)] <- 0
centre_size <- round(dim(MainTcell_test)[1]/10,0)
set.seed(123)
kmeans10 <- kmeans(MainTcell_test, centers = centre_size, nstart = 1)
df_centres <- as.data.frame(kmeans10$centers)

df_centres <- df_centres %>%
  mutate(CellTypist_list_lower = case_when(

    #### T cells ----

    ZNF683 >0 & GNG4 >0 & PDCD1 >0 ~ "CD8aa T cells",
    TOX2 >0 & SATB1 >0 & CCR9 >0 ~ "CD8a/b(entry)",
    PDCD1 >0 & ICOS >0 & CXCR5>0  ~"Follicular helper T cells",
    KLRB1 >0 & SLC4A10 >0 & `TRAV1-2`>0  ~ "KLRB1+ MAIT cells (TRAV1-2+)",
    KLRB1 >0 & `TRAV1-2`>0 & CD8A>0 ~ "KLRB1+ MAIT cells (TRAV1-2+)",
    GZMK >0 & c(CD4>0 | CD4_pos>0) & IL10>0 ~ "Memory CD4+ cytotoxic T cells",
    NKG7>0 & GNLY>0 & CD8A>0 ~ "NK CD8+ T cells",
    CTLA4>0 & IL2RA >0 & FOXP3>0 ~ "Regulatory T cells",
    MIR155HG>0 & BIRC3>0 & SMS>0 ~ "T(agonist)",
    CD8A>0 & CCR7>0 & SELL>0 ~ "Tcm/Naive cytotoxic T cells",
    c(CD4>0 | CD4_pos>0) & CCR7>0 & SELL>0 ~ "Tcm/Naive helper T cells",
    KLRB1>0 & AQP3>0 & ITGB1>0 ~ "Tem/Effector helper T cells",
    PDCD1>0 & c(CD4>0 | CD4_pos>0) & CTLA4>0 ~ "Tem/Effector helper T cells PD1+",
    CX3CR1>0 & GZMB>0 & GNLY>0 ~ "Tem/Temra cytotoxic T cells",
    GZMK>0 & CD8A>0 & CCL5>0 ~ "Tem/Trm cytotoxic CD8+ T cells",
    CD27>0 & CCR7>0 & IKZF2>0 ~ "Treg(diff)",
    ITGA1>0 & ITGAE>0 & CXCR6>0 ~"Trm cytotoxic T cells",
    CCL5>0 & CXCR3>0 & TBX21>0  ~ "Th1 T cells",
    IL7R>0 & CCR6>0 & ZBTB16>0 ~ "Th17 T cells",
    # my annotation
    TRDC>0 & TRGC1>0 & CCL5>0 ~ "Gamma-delta T cells (CCL5+)",
    TRDC>0 & TRGC1>0 & TOX2>0 ~ "Gamma-delta T cells (TOX2+)",
    c(CD4>0 | CD4_pos>0) & CD8A>0 & Chain2>0 ~ "DP thymocytes",
    c(CD4>0 | CD4_pos>0) & FOXP3>0 & CD27>0 ~ "Regulatory T cells",
    c(CD4>0 | CD4_pos>0) &	NKG7>0 &	GNLY>0 ~ "CD4+NKGL7+GNLY+ T cell",
    c(CD4>0 | CD4_pos>0) & GNLY>0  ~ "CD4+GNLY+ T cell",
    c(CD4>0 | CD4_pos>0) & CCL5>0 ~ "Th1 T cells",
    c(CD4>0 | CD4_pos>0) & IL10>0 ~ "IL10+CD4+ T cell",
    c(CD4>0 | CD4_pos>0) & MKI67>0 &	TOP2A>0 ~ "Cycling CD4+ T cell",
    c(CD4>0 | CD4_pos>0) & c(CD3E>0 | CD3D>0 | CD3G>0) ~ "CD4+ T cell",
    c(CD4>0 | CD4_pos>0) & Chain2>0 ~ "CD4+ T cell",
    CD8A>0 & GNLY>0  ~ "CD8+GNLY+ T cell",
    TRDC>0 & TRGC1>0 ~ "Gamma-delta T cells",
    CD4<0 & CD8A<0 & GNLY>0 & Chain2>0 ~ "DN GNLY+ T cell",
    CD4<0 & CD8A<0 & IL2RA>0 & Chain2>0 ~ "DN Tregs cell", #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3230797/
    CD4<0 & CD8A<0 & Chain2>0 ~ "DN T cells",
    ACY3>0 & CD34>0 & SPINK2>0 ~ "ETP",
    FXYD2>0 & HES1>0 & CD99>0 ~ "DN thymocytes",
    CD1A>0 & CD8A>0 & SMPD3>0 ~ "DP thymocytes",

    ### end -----
    TRUE ~ NA_character_))
df_centres$clust_name <- rownames(df_centres)
head(df_centres)
write.csv(df_centres,"df_centres.csv")
