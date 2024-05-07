#' Merging RDS Seurat objects
#' @name merging_multi_SeuratRDS
#' @description
#' This function is to aid im merging multiple Seurat object, which will then need to undergo harmony merging after this step. We have included a loop that merges two Seurat objects at a time to ensure time efficiency  (exponentially gets smaller from the previous sequential merging process)
#'
#' @param set_directory set the path to the directory with the. Default = current working directory. The program will reset the directory to the original to not impact the function of STEGO.R project file. If you stop the script, you will need to manually reset the working directory with setwd("../"), to go back one directory.
#' @param merge_RDS Set to TRUE once you have check the directory with the .rds file
#' @param pattern_RDS uses the list.files function to identify the .rds objects for merging in the 2_scObj
#' @param species Species: hs or mm, as the humans (hs) use upper case and the mouse (mm) gene use proper case.
#' @import Seurat
#' @importFrom purrr reduce
#' @export

merging_multi_SeuratRDS <- function(set_directory = "2_scObj/", merge_RDS = FALSE, pattern_RDS = ".rds$", species = "hs") {

  x <- getwd()
  setwd(set_directory)
  message(paste(getwd(), "is the current work directory"))
  temp <- list.files(pattern = pattern_RDS)
  message(paste(temp, ""))
  len.temp <- length(temp)
  message(paste("There are", len.temp, "rds files in the directory for merging"))

  if (merge_RDS) {
    list.sc <- list()
    features.var.needed <- read.csv(system.file("Kmean", "human.variable.features.csv", package = "STEGO.R"))

    for (i in 1:len.temp) {
      model.name <- strsplit(temp[i], ".rds")[[1]]
      message("Reading in file ", i, " of ", len.temp, ": ", model.name)
      list.sc[[model.name]] <- readRDS(temp[i])
      message("Reducing file size ", model.name)

      if (species == "hs") {
        features.var.needed <- read.csv(system.file("Kmean", "human.variable.features.csv", package = "STEGO.R"))
      } else {
        features.var.needed <- read.csv(system.file("Kmean", "human.variable.features.csv", package = "STEGO.R"))
        features.var.needed$V1 <- tools::toTitleCase(features.var.needed$V1)
      }

      list.sc[[i]] <- subset(list.sc[[i]], features = features.var.needed$V1)
      message("Adding project name")
      list.sc[[i]]@project.name <- model.name
      message("Updating cell Index ID")
      list.sc[[i]]@meta.data$Cell_Index_old <- list.sc[[i]]@meta.data$Cell_Index
      print(list.sc[[i]])
      sl <- object.size(list.sc[[i]])
      message("Stored object is ", round(sl[1]/1000^3, 1), " Gb")
    }

    sl <- object.size(list.sc)
    message("stored object is ", round(sl[1]/1000^3, 1), " Gb")

    # Calculate the number of loops required
    num_objects <- length(list.sc)
    num_loops <- ceiling(log2(num_objects))

    # Print the estimated number of loops required
    message("Estimated number of loops required: ", num_loops)

    # Initialize loop count
    loop_count <- 0

    # Loop to merge pairs of Seurat objects until only one is left
    while (length(list.sc) > 1) {
      loop_count <- loop_count + 1
      message("Loop ",loop_count,"of estimated",num_loops)
      temp_list <- list()

      for (i in seq(1, length(list.sc), by = 2)) {


        if (i < length(list.sc)) {
          # Merge two Seurat objects
          merged_names <- paste(names(list.sc)[i], names(list.sc)[i + 1], sep = " and ")
          new_name <- paste0("Merged_", i)
          message("Merging ", merged_names, " to create ", new_name)
          merged_object <- suppressWarnings(merge(x = list.sc[[i]], y = list.sc[[i + 1]], merge.data = TRUE))
          temp_list[[new_name]] <- merged_object

          # Calculate size of merged object
          merged_size <- object.size(merged_object)
          message("Size of ", new_name, ": ", round(merged_size[1] / 1000^3, 1), " Gb")
        } else {
          # If there's an odd number of objects, keep the last one as is
          temp_list[[names(list.sc)[i]]] <- list.sc[[i]]
        }
      }

      list.sc <- temp_list
    }
    # Calculate total number of files merged
    total_merged <- length(temp) - 1

    merged_object <- list.sc[[1]]
    merged_object@meta.data$Cell_Index_old <- merged_object@meta.data$Cell_Index
    merged_object@meta.data$Cell_Index <- rownames(merged_object@meta.data)
    sl <- object.size(merged_object)
    message(paste(total_merged, "files were merged into 1. The merged object is ", round(sl[1]/1000^3, 2), "Gb"))

    #reset the working directory to the project directory
    setwd(x)
    return(merged_object)
  } else {
    setwd(x)
    message("If these are the files you are looking for, set merge_RDS to TRUE")
    # return(NULL)
  }
}

#' Find Variable features
#' @name harmony_batch_correction_1_variableFeatures
#' @description
#' This function is to aid im merging multiple Seurat object, which will then need to undergo harmony merging
#'
#' @param file merged seurat object
#' @param feature_total_limit limit number of features to create the scaled data for the harmony batch correction.
#' @import Seurat
#' @export
#'

harmony_batch_correction_1_variableFeatures <- function(file = sc, feature_total_limit = 3000) {
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
#' @import Seurat
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

#' Create the finding the PCA, identical to the Seurat pipeline
#' @name harmony_batch_correction_3_PC
#' @description
#' This function is to aid im merging multiple Seurat object, which will then need to undergo harmony merging
#'
#' @param file merged Seurat object
#' @import Seurat
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
#' @import Seurat
#' @import harmony
#' @param file merged seurat object
#' @param selected_column_for_reduction select the column to use for the harmony e.g., orig.ident, Sample_name,
#' @param Maximum_PC_to_use Maximum number principle components (PC) to use. Default30
#' @param resolution_of_clusters Resolutions of the clusters - default 0.5
#' @export
#'

harmony_batch_correction_4_Harmony <- function(file = sc,selected_column_for_reduction = "orig.ident", Maximum_PC_to_use = 30, resolution_of_clusters = 0.5) {

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

Down_sampling_Seurat_OBJ <- function(file = file, selected_column_for_reduction= "Sample_Name",Maximum_PC_to_use = 15,resolution_of_clusters=1) {
  sc = file
  selected_column_for_reduction
  Maximum_PC_to_use
  print("under development, not yet in use")


}

