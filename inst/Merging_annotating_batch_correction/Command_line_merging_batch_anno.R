require(STEGO.R)
# large file merging and annotating
# as the program window may disappear, we recommend merging and annotating using the following commands

###### merging seurat object ------
merging_multi_SeuratRDS(set_directory = "2_scObj/", merge_RDS = F, pattern_RDS = ".rds$") # Check that you are in the correct working directory with your RDS files
sc_merge <- merging_multi_SeuratRDS(set_directory = "2_scObj/", merge_RDS = T, pattern_RDS = ".rds$") # once that is check, switch merge_RDS to TRUE or T

saveRDS(sc_merge,"2_scObj/sc_merge.rds") # save the merged file - it will not have the scaled data and PCA stored any more due to the merging process. This is in case R crashes and it needs to be read in.

## perform the harmony batch correction ------
# sc <- readRDS("sc_merge.rds")

sc <- harmony_batch_correction_1_variableFeatures(file = sc_merge)
sc <- harmony_batch_correction_2_Scaling(file = sc)
sc <- harmony_batch_correction_3_PC(file = sc)
sc <- harmony_batch_correction_4_Harmony(file = sc)

saveRDS(sc,"2_scObj/sc_harmony.rds")


 #### annotating Seurat object -----
sc <- scGate_annotating(file = sc,
                     TcellFunction = T, # switch to T if wa
                     generic = T,
                     stress = T,
                     cycling = T,
                     threshold = 0.2, # change to 0.2 if you use the focused immuen panel from BD rhapsody
                     reductionType = "harmony")

saveRDS(sc,"3_analysis/sc_anno.rds")


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




  } else {
    list_scDown[[2]] <- scOther

  }


}




if (dim(scOther@meta.data)[1]>120000) {
  scDownTwo <- subset(x = scDownOne, downsample = 22000)
  if (dim(scDownTwo@meta.data)[1]>120000) {
    message("Adjust down sampling to below 120K cells")
  }
}

for (i in )

sc <- scGate_annotating(file = sc,
                        TcellFunction = T, # switch to T if wa
                        generic = T,
                        stress = T,
                        cycling = T,
                        threshold = 0.2, # change to 0.2 if you use the focused immuen panel from BD rhapsody
                        reductionType = "harmony")


### annotating Downsampled objects -----
sc <- harmony_batch_correction_1_variableFeatures(file = sc_merge)
sc <- harmony_batch_correction_2_Scaling(file = sc)
