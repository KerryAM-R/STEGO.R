require("harmony")
source("R/global.R")

fils1 <- list.files(path = "inst/extdata/SeuratObj", pattern = "h5Seurat", full.names = TRUE, ignore.case = TRUE)
samples_list <- vector("list", length = length(fils1))
samples_list
for (i in seq_along(fils1)) {
  sc <- LoadH5Seurat(fils1[i])
  samples_list[[i]] <- sc
}
length(fils1)

pbmc.normalized <- (samples_list[[1]])

for (i in 2:length(fils1) ) {
  pbmc.normalized <- merge(pbmc.normalized, y = samples_list[[i]], add.cell.ids = c("", paste("df",i,sep = "")), project = "PBMC12K",merge.data = TRUE)

}

table(pbmc.normalized$orig.ident)

# sc <- NormalizeData(pbmc.normalized)
sc <- FindVariableFeatures(pbmc.normalized, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes)
sc <- RunPCA(sc, features = VariableFeatures(object = sc))
sc <- sc %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)

p1 <- DimPlot(object = sc, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- DimPlot(object = sc, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p1|p2


sc <- sc %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.5) %>%
  identity()


DimPlot(sc, reduction = "umap", group.by = "orig.ident", pt.size = .1)
DimPlot(sc, reduction = "umap",group.by = "Cell_Type_Experimental")
DimPlot(sc, reduction = "umap",group.by = "classify.T.cell_Memory")
DimPlot(sc, reduction = "umap",group.by = "classify.T.cell_Memory")


