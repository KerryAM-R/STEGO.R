require(STEGO.R)
# large file merging and annotating
# as the program window may disappear, we recommend merging and annotating using the following commands

###### merging seurat object ------
merging_multi_SeuratRDS(set_directory = "~/Desktop/STEGO_copy.R/", merge_RDS = F, pattern_RDS = "Anno.rds$") # Check that you are in the correct working directory with your RDS files

sc_merge <- merging_multi_SeuratRDS(set_directory = "~/Desktop/STEGO_copy.R/", merge_RDS = T, pattern_RDS = "Anno.rds$") # once that is check, switch merge_RDS to TRUE or T

unique(sc_merge@meta.data$Tcellfunction)
unique(sc_merge@meta.data$cycling)
unique(sc_merge@meta.data$generic)
unique(sc_merge@meta.data$stress)

sc_merge@meta.data <- sc_merge@meta.data[!grepl("_UCell",names(sc_merge@meta.data))]
sc_merge@meta.data <- sc_merge@meta.data[!grepl("is.pure_",names(sc_merge@meta.data))]
sc_merge@meta.data <- sc_merge@meta.data[!grepl("scGate_multi",names(sc_merge@meta.data))]
saveRDS(sc_merge,"~/Desktop/STEGO_copy.R/sc_merge.rds") # save the merged file - it will not have the scaled data and PCA stored any more due to the merging process. This is in case R crashes and it needs to be read in.

md <- sc_merge@meta.data
sc_merge@meta.data$Tcellfunction <- (ifelse(is.na(sc_merge@meta.data$Tcellfunction),"unknown",md$Tcellfunction))
Idents(sc) <- sc_merge@meta.data$Tcellfunction
sc_merge@meta.data$orig.ident <- sc_merge@meta.data$Tcellfunction

scDown <- subset(x = sc, downsample = 500)
scDown

# DimPlot(sc,group.by = "seurat_clusters", split.by = "Source", raster=FALSE, ncol = 4)
DimPlot(scDown,group.by = "Tcellfunction", split.by = "Tcellfunction", raster=FALSE, ncol = 9)


markerIdentified <- FindMarkers(sc,ident.1 = "unknown", only.pos = T,min.pct = 0.25)
head(markerIdentified)


md1 <- subset(md,md$Tcellfunction=="Tcell.Mono.Complex")
unique(md1$Tcellfunction)
## perform the harmony batch correction ------
# sc <- readRDS("sc_merge.rds")
sc <- harmony_batch_correction_1_variableFeatures(file = sc_merge)
sc <- harmony_batch_correction_2_Scaling(file = sc)
sc <- harmony_batch_correction_3_PC(file = sc)
sc <- harmony_batch_correction_4_Harmony(file = sc)
sc@meta.data <- sc@meta.data[!grepl("_UCell",names(sc@meta.data))]
sc@meta.data <- sc@meta.data[!grepl("is.pure_",names(sc@meta.data))]

saveRDS(sc,"~/Desktop/STEGO_copy.R/sc_harmony.rds")

sc
 #### annotating Seurat object -----
sc <- readRDS()
sc <- scGate_annotating(sc)
sc@meta.data$label <- "ONE"
#### Down sampling if required ----
list_scDown <- list()
Idents(sc) <- sc$seurat_clusters

scDownOne <- subset(x = sc, downsample = 5000)
if (dim(scDownOne@meta.data)[1]>120000) {
  message("Adjust down sampling to below 120K cells")
} else {
  list_scDown[[1]] <- scDownOne
  x <- colnames(scDownOne)
  sc$keep <- ifelse(sc$Cell_Index %in% setdiff(colnames(sc),x) ,"yes","no")
  table(sc$keep)
  scOther <- subset(sc, subset = keep=="yes")



  if( dim(scOther@meta.data)[1]>120000) {
    scDownTwo <- subset(x = scOther, downsample = 5000)
    list_scDown[[2]] <- scDownTwo

  }


}




if (dim(scOther@meta.data)[1]>120000) {
  scDownTwo <- subset(x = scDownOne, downsample = 22000)
  if (dim(scDownTwo@meta.data)[1]>120000) {
    message("Adjust down sampling to below 120K cells")
  }
}

### annotating Downsampled objects -----
sc <- harmony_batch_correction_1_variableFeatures(file = sc_merge)
sc <- harmony_batch_correction_2_Scaling(file = sc)
