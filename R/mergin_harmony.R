#' Merging RDS Seurat objects
#' @name merging_multi_SeuratRDS
#' @description
#' This function is to aid im merging multiple Seurat object, which will then need to undergo harmony merging
#'
#' @param set_directory set the path to the directory with the. Default = current working directory. The program will reset the directory to the original to not impact the function of STEGO.R project file
#' @param merge_RDS Set to TRUE once you have check the directory with the .rds file
#' @param pattern_RDS uses the list.files function to identify
#' @param species Species: hs or mm, as the humans use upper case and the mouse gene use proper case.
#' @export

merging_multi_SeuratRDS <- function (set_directory = ".", merge_RDS = F, pattern_RDS = ".rds$",species = "hs") {
  x <- getwd()
  require(purrr)
  require(Seurat)
  setwd(set_directory)
  message(paste(getwd(),"is the current work directory"))
  temp = list.files(pattern=pattern_RDS)
  message(paste(temp,""))
  len.temp <- length(temp)
  message(paste("There are",len.temp,"rds files in the directory for merging"))
  if (merge_RDS) {
    list.sc <- list()
    features.var.needed <- read.csv(system.file("Kmean","human.variable.features.csv",package = "STEGO.R"))
    dim(features.var.needed)
    for (i in 1:len.temp) {
      model.name <- strsplit(temp,".rds")[[i]][1]
      message("reading in file ", model.name)
      list.sc[[model.name]] <- readRDS(temp[i])
      message("Reducing file size ", model.name)

      if (species == "hs") {
        features.var.needed <- read.csv(system.file("Kmean","human.variable.features.csv",package = "STEGO.R"))
      } else {
        features.var.needed <- read.csv(system.file("Kmean","human.variable.features.csv",package = "STEGO.R"))
        features.var.needed$V1 <- str_to_title(features.var.needed$V1)
      }
      list.sc[[i]] <- subset(list.sc[[i]],features=features.var.needed$V1)
      message("adding project name")
      list.sc[[i]]@project.name <- model.name
      message("Updating cell Index ID")
      list.sc[[i]]@meta.data$Cell_Index_old <- list.sc[[i]]@meta.data$Cell_Index
      print(list.sc[[i]])
      sl <- object.size(list.sc[[i]])
      message("stored object is ",round(sl[1]/1000^3,1)," Gb")
    }
  sl <- object.size(list.sc)
  message("stored object is ",round(sl[1]/1000^3,1)," Gb")
  merged_object <- reduce(list.sc, function(x, y) {
    merge(x = x, y = y, merge.data = TRUE)
  })
  merged_object@meta.data$Cell_Index_old <- merged_object@meta.data$Cell_Index
  merged_object@meta.data$Cell_Index <- rownames(merged_object@meta.data)
  sl <- object.size(merged_object)
  message("The merged object is ",round(sl[1]/1000^3,1)," Gb")

    setwd(x)
    merged_object
  } else {
    setwd(x)
    message("If these are the files you are looking for, set merge_RDS to T")
  }
}

#' Find Variable features
#' @name harmony_batch_correction_1_variableFeatures
#' @description
#' This function is to aid im merging multiple Seurat object, which will then need to undergo harmony merging
#'
#' @param file merged seurat object
#' @param feature_total total number of Fetures to find. If there are fewer then the selected total, it will use that length instead
#' @export
#'

harmony_batch_correction_1_variableFeatures <- function(file = sc, feature_total_limit = 3000) {
  require(Seurat)
  sc <- file
  all.genes <- rownames(sc)
  feature_total <- length(all.genes)
  num_featureTotal = feature_total_limit
  if(feature_total<num_featureTotal) {
    paste("There are less than",num_featureTotal,"Therefore, setting total Features to",feature_total)
    num_featureTotal <- feature_total
  } else {
    paste("There are more than",feature_total,"Therefore, setting total Features to",num_featureTotal)
  }
  sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = num_featureTotal)
  sc
  }

#' Create the scaling data
#' @name harmony_batch_correction_2_scaling
#' @description
#' Scales the Seurat object and ensure that the correct genes are present for annotating purposes
#'
#' @param file merged seurat object
#' @param Seruat_version Seurat version
#' @param species select either human 'hs' or mouse 'mm'
#' @export
#'

harmony_batch_correction_2_Scaling <- function(file = sc,  Seruat_version = "V4",species = "hs") {
  sc <- file
  kmeans <- read.csv(system.file("Kmean","Kmeans.requires.annotation.csv",package = "STEGO.R"))
  if(Seruat_version=="V4") {
    var.genes <- as.data.frame(sc@assays$RNA@var.features)
    print(head(var.genes))
  } else {
    var.genes <- sc@assays$RNA@meta.data$var.features
    var.genes <- as.data.frame(var.genes[!is.na(var.genes)])
    print(head(var.genes))
  }

  names(var.genes) <- "V1"
  if (species == 'hs') {
    kmeans2 <- as.data.frame(kmeans$Human)
    names(kmeans2) <- "V1"
  } else {
    kmeans2 <- as.data.frame(kmeans$Mouse)
    names(kmeans2) <- "V1"
  }
  Gene_to_scale <- unique(rbind(var.genes,kmeans2))
  names(Gene_to_scale) <- "V1"
  print(head(Gene_to_scale))
  print(dim(Gene_to_scale))

  sc <- ScaleData(sc, features = Gene_to_scale$V1)
  sc
}

#' Create the finding the PCA, identifical to the Seurat pipeline
#' @name harmony_batch_correction_3_PC
#' @description
#' This function is to aid im merging multiple Seurat object, which will then need to undergo harmony merging
#'
#' @param file merged seurat object
#' @export
#'

harmony_batch_correction_3_PC <- function(file = sc) {
  sc <- RunPCA(file)
  sc
}

#' Performign the Haromy integration
#' @name harmony_batch_correction_4_Harmony
#' @description
#' This function is to aid im merging multiple Seurat object, which will then need to undergo harmony merging
#'
#' @param file merged seurat object
#' @param selected_column_for_reduction select the column to use for the harmony e.g., orig.ident, Sample_name,
#' @param Maximum_PC_to_use Maximum number principle components (PC) to use. Default30
#' @param resolution_of_clusters Resolutions of the clusters - default 0.5
#' @export
#'

harmony_batch_correction_4_Harmony <- function(file = sc,selected_column_for_reduction = "orig.ident", Maximum_PC_to_use = 30, resolution_of_clusters = 0.5) {
  require("harmony")
  sc <- file
  sc <- sc %>%
    RunHarmony(selected_column_for_reduction, plot_convergence = TRUE)

  sc <- sc %>%
    RunUMAP(reduction = "harmony", dims = 1:Maximum_PC_to_use) %>%
    FindNeighbors(reduction = "harmony", dims = 1:Maximum_PC_to_use) %>%
    FindClusters(resolution = resolution_of_clusters) %>%
    identity()
  sc
}

#' Create the finding the PCA, identifical to the Seurat pipeline
#' @name Down_sampling_Seurat_OBJ
#' @description
#' This function is to aid im merging multiple Seurat object, which will then need to undergo harmony merging
#'
#' @param file merged seurat object
#' @param selected_column_for_reduction select the column to use for the harmony e.g., orig.ident, Sample_name,
#' @param Maximum_PC_to_use Maximum number principle components (PC) to use. Default30
#' @param resolution_of_clusters Resolutions of the clusters - default 0.5
#' @export
#'

Down_sampling_Seurat_OBJ <- function(file = file) {
  sc = file

}

