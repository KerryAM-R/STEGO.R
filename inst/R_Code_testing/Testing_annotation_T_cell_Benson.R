source("R/global.R")
sc <- LoadH5Seurat("~/Documents/Post-doc/PostDoc_Antwerp/R_files/1_Files_scTCRseq/Antwerp/data_mmr_dextramer/M6_2023.03.17.h5Seurat")
df= as.data.frame(sc[["RNA"]]@scale.data)
df= as.data.frame(sc[["RNA"]]@counts)
MainTcell <- as.data.frame(t(df))
# names(MainTcell) <- ifelse(grepl("pAbO",names(MainTcell)),names(MainTcell),toupper(names(MainTcell)))
MainTcell$Cell_Index <- rownames(MainTcell)



# names(MainTcell) <- gsub("CD8B.","CD8B",names(MainTcell))
md <- sc@meta.data
names(md)
md$Motif_gene <- md[,names(md) %in% "vj_gene_cdr3_AG_BD"]
names(md)
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
"seurat_clusters"
md2 <- md[names(md) %in% c("Cell_Index","Chain2","vdj_gene_cdr3_AG_BD","CD4_pos")]
MainTcell <- merge(md2,MainTcell,by="Cell_Index",all.y=T)
rownames(MainTcell) <- MainTcell$Cell_Index

MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD4","CD8A","CD8B","CD3E","CD3D","CD3G","vdj_gene_cdr3_AG_BD")]

paired <- read.csv("~/Downloads/paired_chain_CDR32023.03.17 (1).csv")
MainTcell_test_paired <- merge(MainTcell_test,paired,by = "vdj_gene_cdr3_AG_BD")

write.csv(MainTcell_test_paired,"MainTcell.tcr.count.csv")

names(MainTcell)[grepl("HLA",names(MainTcell))]

# FeaturePlot(sc, features = "HLA-DR", min.cutoff = 'q10', label = T)

MainTcell_test <- MainTcell[,names(MainTcell) %in% c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21",
                                                     "CCR6","GZMB","PRF1","GNLY","KLRB1","TGFB1","IGHM","MS4A1","IGHA2","IGHA1","JCHAIN","IGHG1","IGHG2","IGHG3","IGHG4","IGHE",
                                                     "IGHD","CD3E","CD3D","CD3G","CD19","IL6","CD27","IL10","CD4_pos")]
df <- MainTcell_test
# Vector of columns you want in this data
nms <- c("CD4","CD8A","CD8B","FOXP3","IL2RA","IL12A","IL18","CXCR3","IL2","IFNG","TNF","IL4","IL9","IL17A","IL17F","IL22","IL21",
         "CCR6","GZMB","PRF1","GNLY","KLRB1","TGFB1","IGHM","MS4A1","IGHA2","IGHA1","JCHAIN","IGHG1","IGHG2","IGHG3","IGHG4","IGHE",
         "IGHD","CD3E","CD3D","CD3G","CD19","IL6","CD27","IL10")

Missing <- setdiff(nms, names(df))  # Find names of missing columns
df[Missing] <- 0                    # Add them, filled with '0's
MainTcell_test <- df[nms]

d_frame <- na.omit(MainTcell_test)
centre_size <- round(dim(MainTcell_test)[1]/20,0)
set.seed(123)
kmeans10 <- kmeans(d_frame, centers = centre_size, nstart = 1)
df_centres <- as.data.frame(kmeans10$centers)

df_centres <- df_centres %>%
  mutate(classify.T.cell_Function = case_when(
    c(CD4 > 0 |CD4_pos>0) & c(CD8A<0 | CD8B<0) & FOXP3 >0 & IL2RA >0 ~ "CD4+ Treg FOXP3+CD25+",
    c(CD4 > 0 |CD4_pos>0) & c(CD8A>0 | CD8B>0) & TGFB1>0 & FOXP3 >0 ~ "CD8+ Tregs FOXP3+TGFB1+",
    c(CD4 > 0 |CD4_pos>0) & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & IL2>0 & IFNG>0 ~ "IL2+IFNg+ (Th1)",
    c(CD4 > 0 |CD4_pos>0) & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & TNF >0 & IFNG>0 ~ "IFNg+TNF+ (Th1)",
    c(CD4 > 0 |CD4_pos>0) >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & TNF >0 & IFNG<0 ~ "IFNg-TNF+ (Th1)",
    c(CD4 > 0 |CD4_pos>0) >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & TNF <0 & IFNG>0 ~ "IFNg+TNF- (Th1)",
    c(CD4 > 0 |CD4_pos>0) >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & IL2>0 ~ "Th1 IL2+",
    c(CD4 > 0 |CD4_pos>0) >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) & IFNG>0 ~ "CD4+ IFNg+",

    # CD4 >0 & CD8A<0 & CD8B<0 & c(IL12A >0 | IL18 > 0 | CXCR3>0) ~ "Th1 (other)",
    c(CD4 > 0 |CD4_pos>0) & IL4>0 ~ "CD4+IL4+ (Th2)",
    c(CD4 > 0 |CD4_pos>0) & c(IL9>0) ~ "CD4+IL9+ (Th9)",
    c(CD4 > 0 |CD4_pos>0) & c(IL17A>0 | IL17F>0) ~ "CD4+IL17+ (Th17)",
    c(CD4 > 0 |CD4_pos>0) & c(IL22>0) ~ "IL22+ (Th22)",
    c(CD3D>0 | CD3E >0 |CD3G>0 ) & c(IL22>0) ~ "IL22+",
    CD4 >0 & c(IL21>0) ~ "IL21+ (Tfh)",
    CD8A> 0 & CCR6>0 & KLRB1>0 ~ "CD8a+CCR6+ KLRB1+ (MAIT)",
    CD4<0 & c(CD8A>0 | CD8B>0) & GZMB < 0 & GNLY>0 & PRF1<0 ~ "GNLY+CD8+ Cytotoxic",
    CD4<0 & c(CD8A>0 | CD8B>0) & GZMB > 0 & GNLY<0 & PRF1>0 ~ "GZMB+PRF1+CD8+ Cytotoxic",
    CD4<0 & c(CD8A>0 | CD8B>0) & TNF > 0   & IFNG > 0~ "IFNg+TNF+CD8+ Cytotoxic",
    CD4<0 & c(CD8A>0 | CD8B>0) & TNF < 0   & IFNG > 0~ "IFNg+TNF-CD8+ Cytotoxic",
    CD4<0 & c(CD8A>0 | CD8B>0) & TNF > 0   & IFNG < 0~ "IFNg-TNF+CD8+ Cytotoxic",
    CD4<0 & c(CD8A>0 | CD8B>0) & c(GZMB > 0 | GNLY>0 | PRF1 >0) ~ "CD8+ Cytotoxic",
    CD4>0 & CD8A<0 & CD8B<0 & c(GZMB > 0 | GNLY>0 | PRF1 >0) ~ "CD4+ Cytotoxic",

    CD4<0 & CD8A <0 & CD8B <0 & GZMB > 0 & GNLY>0 & PRF1 >0 ~ "GZMB+PRF1+GNLY+CD4-CD8- Cytotoxic",
    CD4<0 & CD8A <0 & CD8B <0 & GZMB < 0 & GNLY>0 & PRF1 <0 ~ "GZMB-PRF1-GNLY+CD4-CD8- Cytotoxic",
    CD4<0 & CD8A <0 & CD8B <0 & GZMB > 0 & GNLY<0 & PRF1 >0 ~ "GZMB+PRF1+GNLY-CD4-CD8- Cytotoxic",

    c(CD3D>0 | CD3E >0 |CD3G>0 )&CD4<0 & CD8A <0 & CD8B <0 & IFNG > 0 & TNF>0 ~ "IFNg+TNF+ CD4-CD8- Cytotoxic",
    c(CD3D>0 | CD3E >0 |CD3G>0 )&CD4<0 & CD8A <0 & CD8B <0 & IFNG < 0 & TNF>0 ~ "IFNg-TNF+ CD4-CD8- Cytotoxic",
    c(CD3D>0 | CD3E >0 |CD3G>0 )&CD4<0 & CD8A <0 & CD8B <0 & IFNG > 0 & TNF<0 ~ "IFNg+TNF- CD4-CD8- Cytotoxic",

    c(CD3D>0 | CD3E >0 |CD3G>0 )&CD4<0 & CD8A <0 & CD8B <0 & c(GZMB > 0 | GNLY>0 | PRF1 >0) ~ "CD4-CD8- Cytotoxic",

    CD4>0 & c(CD8A >0 | CD8B >0) & c(GZMB > 0 | GNLY>0 | PRF1 >0) ~ "CD4+CD8+ Cytotoxic",
    CD4>0 & TNF>0 & IFNG>0 ~ "IFNg+TNF+CD4+",
    CD4>0 & TNF>0 & IFNG<0 ~ "IFNg-TNF+CD4+",
    CD4>0 & TNF<0 & IFNG>0 ~ "IFNg+TNF-CD4+",
    # B cell markers
    # CD19>0 & IGHM>0 & MS4A1>0 & IL6<0~ "IgM+ B cell", # ,"SDC1" = CD138 (plasma cell)
    TRUE ~ NA_character_))
df_centres$clust_name <- rownames(df_centres)
head(df_centres)
