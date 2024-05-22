#' Automated preprocessing
#' @name automated_sc_filtering
#' @description
#' This function is to automate step 3a. We highly recommend that you check your base parameters in the interface and then apply them to all of the data. The default parameters are set for human 10x
#'
#' @param folder This is the location of the matrix and md files from the automated pre-processing
#' @param dataset_type Dataset types. Please use "10x" or "BD_rap". There is a slight difference in how the barcodes are labeled with 10x having strings while BD has numbers
#' @import Seurat
#' @import ggplot2
#' @param species either select hs or mm. This is due to the gene case with all caps for hs and proper case for mm.
#' @param features.min minimum number of features
#' @param features.max Maximum features
#' @param percent.mt Mitochondria DNA cut-off (<)
#' @param percent.rb Ribosomal RNA cut-off (>)
#' @param dimension_sc Set number of dimensions to use for the dimensional reduction.
#' @param resolution_sc Set the Seurat unsupervised clustering; however, this wont be used for annotation purposes in this pipeline.
#' @param limit_to_TCR_GEx Reduce the file to limit to the TCR seq and GEx only (recommended for large data sets)
#' @param save_plots Save the plots
#' @param output_dir Location of the process .rds Seurat object filtered files
#' @return Seurat object per file.
#' @export

automated_sc_filtering <- function(folder = "1_SeuratQC",
                                   dataset_type = "10x",
                                   species = "hs",
                                   features.min = 200,
                                   features.max = 6000, percent.mt = 20,
                                   percent.rb = 5,
                                   dimension_sc = 15,
                                   resolution_sc= 1,
                                   limit_to_TCR_GEx = F,
                                   save_plots = T,
                                   output_dir = "3_SCobj/3a/",
) {

  main_directory <- "1_SeuratQC"
  # main_directory <- main_directory
  main_folders <- list.files(paste(main_directory),full.names = T)
  samp_names <- main_folders[grepl("_count-matrix_10x.csv.gz", main_folders)]
  samp_names <- gsub("_count-matrix_10x.csv.gz","",samp_names)
  main_folders2 <- list.files(paste(main_directory),full.names = F)

  samp_names2 <- main_folders2[grepl("_count-matrix_10x.csv.gz", main_folders2)]
  samp_names2 <- gsub("_count-matrix_10x.csv.gz","",samp_names2)
  num <- length(samp_names2)
  for (i in 1:num) {
    message(paste("processing",i,"of",num))
    project_name <- samp_names[i]
    project_name
    project_name2 <- samp_names2[i]
    print(project_name2)
    sc_processed <- paste0(output_dir,project_name2,"_md_added",".rds")

    if(file.exists(sc_processed)) {

      message(sc_processed," already exists")
      next

    } else {
      message("Started processing",project_name2)
      if (dataset_type == "10x") {
        message("reading in matrix file")
        file.exists(paste0(samp_names[i],"_count-matrix_10x.csv.gz"))
        mat <- read.csv(paste0(samp_names[i],"_count-matrix_10x.csv.gz"))

        head(mat)[1:6]

        message("reading in meta data file")
        md <- read.csv(paste0(samp_names[i],"_metadata_10x.csv"))

        md <-  md[!duplicated(md$Cell_Index),]
        head(md)
      } else if (dataset_type == "BD_rap") {
        message("reading in matrix file")
        mat <- read.csv(paste0(samp_names[i],"_count-matrix_10x.csv.gz"), row.names = 1)
        message("reading in meta data file")
        md <- read.csv(paste0(samp_names[i],"_metadata_10x.csv"))
        # randomly removing duplication ----
        md <-  md[!duplicated(md$Cell_Index),]
      } else {
        message("Please select 10x or BD_rap")
      }

      # reading in the matrix data ---
      if (dataset_type == "10x" && species == "hs") {
        names(mat) <- gsub("[.]1", "-1", names(mat))
        rownames(mat) <- make.unique(mat$Gene_Name)
        mat2 <- mat[, !names(mat) %in% c("Gene_Name")]

        sc <- suppressWarnings(CreateSeuratObject(counts = mat2, assay = "RNA", project = project_name2))
        sc <- PercentageFeatureSet(sc, pattern = "^MT-", col.name = "mtDNA")
        sc <- PercentageFeatureSet(sc, pattern = "^RP[SL]", col.name = "rRNA")
        head(sc@meta.data)

      } else if (dataset_type == "10x" && species == "mm") {
        names(mat) <- gsub("[.]1", "-1", names(mat))
        rownames(mat) <- make.unique(mat$Gene_Name)
        mat2 <- mat[, !names(mat) %in% c("Gene_Name")]
        sc <- CreateSeuratObject(counts = mat, assay = "RNA", project = project_name2)
        rm(mat)
        rm(mat2)

        sc <- suppressMessages(PercentageFeatureSet(sc, pattern = "^Mt", col.name = "mtDNA"))
        sc <- suppressMessages(PercentageFeatureSet(sc, pattern = "Rp[sl]", col.name = "rRNA"))

      } else if (dataset_type == "BD_rap" && species == "mm") {
        names(mat) <- as.character(gsub("X", "", names(mat)))
        sc <- CreateSeuratObject(counts = mat, assay = "RNA", project = project_name2)
        sc <- PercentageFeatureSet(sc, pattern = "^Mt", col.name = "mtDNA")
        sc <- PercentageFeatureSet(sc, pattern = "Rp[sl]", col.name = "rRNA")

      } else if (dataset_type == "BD_rap" && species == "hs") {
        names(mat) <- as.character(gsub("X", "", names(mat)))

        sc <- suppressMessages(CreateSeuratObject(counts = mat, assay = "RNA", project = project_name2))
        rm(mat)
        sc <- suppressMessages(PercentageFeatureSet(sc, pattern = "^MT-", col.name = "mtDNA"))
        sc <- suppressMessages(PercentageFeatureSet(sc, pattern = "^RP[SL]", col.name = "rRNA"))

      } else {
        message("Please choose either 10x or BD_rap as well as define hs or mm")
      }
      print(sc)
      print(head(sc@meta.data))

      if(save_plots) {
        dir.create(paste0("Figures.Tables/QC_Figures/",project_name2))
        suppressWarnings({
          sc@meta.data$mtDNA[is.na(sc@meta.data$mtDNA)] <- 0
          sc@meta.data$rRNA[is.na(sc@meta.data$rRNA)] <- 0
            plot <- suppressWarnings(suppressMessages(VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "mtDNA", "rRNA"), ncol = 2)))
          file_name_before <- paste0("Figures.Tables/QC_Figures/",project_name2,"/",project_name2,"_1_before_filtering.png")
          png(file_name_before, width = 1000, height = 1200, res = 144)
          plot(plot)
          dev.off()
        })
      }

      sc <- suppressWarnings(suppressMessages(subset(sc, subset = nFeature_RNA >= features.min & nFeature_RNA <= features.max & mtDNA <= percent.mt & rRNA >= percent.rb)))

      if(save_plots) {
        suppressWarnings({
          if(sum(sc@meta.data$mtDNA)==0 || sum(sc@meta.data$rRNA)==0) {
            plot <- suppressWarnings(suppressMessages(VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)))
          } else {
            plot <- suppressWarnings(suppressMessages(VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "mtDNA", "rRNA"), ncol = 2)))
          }
          file_name_after <- paste0("Figures.Tables/QC_Figures/",project_name2,"/",project_name2,"_2_after_filtering.png")
          png(file_name_after, width = 1000, height = 1200, res = 144)
          plot(plot)
          dev.off()
        })
      }
      sc <- suppressMessages(NormalizeData(sc, verbose = F))

      sc <- suppressMessages(FindVariableFeatures(sc, selection.method = "vst", verbose = F))

      if(save_plots) {
        suppressWarnings({
          top10 <- head(VariableFeatures(sc, verbose = F), 10)
          # plot variable features with and without labels
          plot1 <- VariableFeaturePlot(sc)
          plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

          file_name_top10 <- paste0("Figures.Tables/QC_Figures/",project_name2,"/",project_name2,"_3_Top_features.png")
          png(file_name_top10, width = 1000, height = 1200, res = 144)
          plot(plot2)
          dev.off()
        })
      }

      all.genes <- rownames(sc)
      message("Scaling")
      sc <- suppressWarnings(suppressMessages(ScaleData(sc, features = all.genes, verbose = F)))
      message("PCA")
      sc <- suppressMessages(RunPCA(sc, features = VariableFeatures(object = sc),verbose = F))
      message("Finding Neighbours")
      sc <- suppressMessages(FindNeighbors(sc, dims = 1:dimension_sc, verbose = F))
      message("Finding clusters")
      sc <- suppressMessages(FindClusters(sc, resolution = resolution_sc, verbose = F))
      message("make UMAP")
      sc <- RunUMAP(sc, dims = 1:dimension_sc, verbose = F)

      if(save_plots) {
        suppressWarnings({
          plot <- ElbowPlot(sc)

          file_name_top10 <- paste0("Figures.Tables/QC_Figures/",project_name2,"/",project_name2,"_4_Elbow_plot.png")
          png(file_name_top10, width = 1000, height = 1200, res = 144)
          plot(plot)
          dev.off()
        })
      }

      if(save_plots) {
        suppressWarnings({
          plot <- DimPlot(sc, reduction = "umap")
          file_name_top10 <- paste0("Figures.Tables/QC_Figures/",project_name2,"/",project_name2,"_5_UMAP_plot.png")
          png(file_name_top10, width = 1000, height = 1200, res = 144)
          plot(plot)
          dev.off()
        })
      }

      message("adding in meta data")
      sc@meta.data$Cell_Index <- rownames(sc@meta.data)

      sc@meta.data$order <- 1:length(sc@meta.data$Cell_Index)
      scMeta.data <- sc@meta.data
      meta.data2 <- merge(scMeta.data, md, by = "Cell_Index", all.x = T)
      sc@meta.data <- meta.data2
      rownames(sc@meta.data) <- sc@meta.data$Cell_Index
      sc@meta.data <- sc@meta.data[order((sc@meta.data$order)), ]

      if (limit_to_TCR_GEx) {
        print(sc)
        sc <- subset(sc, subset = chain_AG == "TRA" | chain_AG == "TRG" | chain_AG == "TRD")
        print(sc)
      }

      sc_processed <- paste0(output_dir,project_name2,"_md_added",".rds")

      saveRDS(sc,sc_processed)
      message("Saved ", sc_processed)
      rm(sc)
    }
  }
}
