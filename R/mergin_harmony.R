#' Merging the ClusTCR2 files
#' @name merging_ClusTCR2
#' @description
#' This function is to aid im merging multiple Seurat object, which will then need to undergo harmony merging after this step. We have included a loop that merges two Seurat objects at a time to ensure time efficiency  (exponentially gets smaller from the previous sequential merging process)
#'
#' @param directory location of the multiple clusTCR2 files from step 1. default is "1_ClusTCR"
#' @param chain Either select AG or BD. The AG are in the same column so that the A will pair with B and G with a D.
#' @param output_dir location of where the file will be saved. default is "2_ClusTCR2"
#' @return A merged file for clustering step based in ClusTCR2
#' @export

merging_ClusTCR2 <- function (directory = "1_ClusTCR",chain = "AG", output_dir = "2_ClusTCR2") {
  main_directory <- directory
  main_folders <- list.files(paste(main_directory), full.names = T)
  main_folders <- main_folders[grepl(chain,main_folders)]

  chain2 <- paste0(chain,"_")

  main_folders2 <- list.files(paste(main_directory), full.names = F)

  samp_names2 <- main_folders2[grepl(chain2,
                                     main_folders2)]
  samp_names2 <- gsub("_clusTCR_10x.csv", "", samp_names2)

  num <- length(samp_names2)

  myfiles <- list()

  for (i in 1:num) {
    clust <- read.csv(main_folders[i])
    head(clust)

    myfiles[[i]] <- clust
  }

  df <- rbind(myfiles[[1]])
  for (i in 2:length(myfiles)) {
    df <- rbind(df, myfiles[[i]])
  }

  if (nrow(df[-c(grep("[*]", df$junction_aa)), ] > 0)) {
    df <- df[-c(grep("[*]", df$junction_aa)), ]
  }

  df$CDR3_beta <- gsub("^$", "None", df$junction_aa)

  if (nrow(df[-c(grep("None", df$CDR3_beta)), ] > 0)) {
    df <- df[-c(grep("None", df$CDR3_beta)), ]
  }

  df <- df[!duplicated(df$junction_aa), ]
  print(head(df))

  name_file <- paste0(output_dir,"/",chain, "_Multi_ClusTCR2.csv", sep = "")
  write.csv(df,name_file,row.names = F)
  message("Saved", name_file)

}

#' Merging RDS Seurat objects
#' @name merging_multi_SeuratRDS
#' @description
#' This function is to aid im merging multiple Seurat object, which will then need to undergo harmony merging after this step. We have included a loop that merges two Seurat objects at a time to ensure time efficiency  (exponentially gets smaller from the previous sequential merging process)
#'
#' @param seurat_files set the path to the directory with the. Default = current working directory. The program will reset the directory to the original to not impact the function of STEGO.R project file. If you stop the script, you will need to manually reset the working directory with setwd("../"), to go back one directory.
#' @param merge_RDS Set to TRUE once you have check the directory with the .rds file
#' @param pattern_RDS uses the list.files function to identify the .rds objects for merging in the 2_scObj
#' @param species Species: hs or mm, as the humans (hs) use upper case and the mouse (mm) gene use proper case.
#' @param reduce_size Reduce the size of the matrix to the most variable TCR based on the 12 dataset.
#' @param own_features If you want to reduce the file size to your own features change to TRUE
#' @param own_features_df Add in a data frame of the list of your features of one column only.
#' @param Seurat_Version Current version of Seurat V5 requires joining layers in the merging process.
#' @param max_features set the maximum number of features to keep per file.
#' @import Seurat
#' @importFrom purrr reduce
#' @return A merged Seurat RDS object.
#' @export

merging_multi_SeuratRDS <- function(seurat_files = "3_SCobj/3a",
                                    merge_RDS = FALSE,
                                    pattern_RDS = ".rds$",
                                    species = "hs",
                                    reduce_size = TRUE,
                                    own_features = FALSE,
                                    Seurat_Version = "V5",
                                    max_features = 10000,
                                    own_features_df = own_features_df
                                    )  {

  main_directory <- seurat_files
  main_folders <- list.files(paste(main_directory), full.names = T)

  temp <- list.files(main_directory,pattern = pattern_RDS, full.names = T)
  temp2 <- list.files(main_directory,pattern = pattern_RDS, full.names = F)
  len.temp <- length(temp)

  message(paste("There are", len.temp, "rds files in the directory for merging"))

  if (merge_RDS) {
    list.sc <- list()
    features.var.needed <- read.csv(system.file("Kmean", "human.variable.features.csv", package = "STEGO.R"))

    for (i in 1:len.temp) {
      model.name <- strsplit(temp2[i], "_md_added.rds")[[1]]
      message("Reading in file ", i, " of ", len.temp, ": ", model.name)
      list.sc[[model.name]] <- readRDS(temp[i])
      if (reduce_size) {
      message("Reducing file size ", model.name)

      if (species == "hs") {
        list.sc[[model.name]] <- FindVariableFeatures(list.sc[[model.name]], selection.method = "vst", nfeatures = max_features)
        var.feat <- VariableFeatures(list.sc[[model.name]])
        var.feat <- as.data.frame(var.feat)
        names(var.feat) <- "V1"

        anno.features <- read.csv(system.file("Kmean", "Kmeans.requires.annotation.csv", package = "STEGO.R"))
        names(anno.features) <- "V1"
        anno.features2 <- merge(var.feat,anno.features,by = "V1",all.x = T)

        if (own_features) {
          ownfeatures <- as.data.frame(own_features_df)
          names(ownfeatures) <- "V1"

          features.var.needed <- merge(anno.features2,ownfeatures, by = "V1", all.x = T)
        }

      } else {
        features.var.needed <- read.csv(system.file("Kmean", "human.variable.features.csv", package = "STEGO.R"))
        features.var.needed$V1 <- tools::toTitleCase(features.var.needed$V1)

        if (own_features) {
          features.var.needed <- own_features_df
          features.var.needed$V1 <- tools::toTitleCase(features.var.needed$V1)
        }
      }
      message("removing previous scaled data")
      list.sc[[i]][["RNA"]]$scale.data <- NULL
      message("selecting the variable features")
      list.sc[[i]] <- subset(list.sc[[i]], features = features.var.needed$V1)
      }

      message("Adding project name")
      list.sc[[i]]@project.name <- model.name
      message("Updating cell Index ID")
      list.sc[[i]]@meta.data$Cell_Index_old <- list.sc[[i]]@meta.data$Cell_Index
      print(list.sc[[i]])
      sl <- object.size(list.sc[[i]])
      message("Stored object is ", round(sl[1]/1000^3, 3), " Gb")
    }

    print(length(list.sc))
    sl <- object.size(list.sc)
    message("stored object is ", round(sl[1]/1000^3, 3), " Gb")

    # remove potential graph issue.
    for (i in 1:length(list.sc)) {
      list.sc[[i]]@graphs <- list()
      list.sc[[i]]@misc <- list()


    }

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

          if(Seurat_Version == "V5") {
            message("Joining layers of ",new_name)
            merged_object <- JoinLayers(merged_object,  assay = "RNA")
          }

          temp_list[[new_name]] <- merged_object
          # Calculate size of merged object
          merged_size <- object.size(merged_object)
          message("Size of ", new_name, ": ", round(merged_size[1] / 1000^3, 1), " Gb")
        } else {
          # If there's an odd number of objects, keep the last one as is
          temp_list[[names(list.sc)[i]]] <- list.sc[[i]]
        }
      }
      print(temp_list)
      list.sc <- temp_list
    }
    # Calculate total number of files merged
    total_merged <- length(temp)

    merged_object <- list.sc[[1]]
    merged_object@meta.data$Cell_Index_old <- merged_object@meta.data$Cell_Index
    merged_object@meta.data$Cell_Index <- rownames(merged_object@meta.data)
    sl <- object.size(merged_object)
    message(paste(total_merged, "files were merged into 1. The merged object is ", round(sl[1]/1000^3, 2), "Gb"))

    #reset the working directory to the project directory
    return(merged_object)
  } else {
    message("If these are the files you are looking for, set merge_RDS to TRUE")
    # return(NULL)
  }
}

#' Find Variable features
#' @name harmony_batch_correction_1_variableFeatures
#' @description
#' This function is to aid im merging multiple Seurat object, which will then need to undergo harmony merging
#'
#' @param file merged Seurat object
#' @param feature_total_limit limit number of features to create the scaled data for the harmony batch correction. Default = 3000 features.
#' @import Seurat
#' @return Seurat object with the variable features.
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
#' This function scales a Seurat object and ensures that the correct genes are present for annotating purposes.
#' It is intended to be used sequentially following `harmony_batch_correction_1_variableFeatures`.
#' @import Seurat
#' @param file merged seurat object
#' @param Seruat_version Seurat version
#' @param species select either human 'hs' or mouse 'mm'
#' @return A Seurat object with scaled gene expression data.
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
#' This function is to aid in merging multiple Seurat object, which will then need to undergo harmony merging
#' @param file merged Seurat object
#' @import Seurat
#' @return A Seurat object containing the results of PCA.
#' @export
#'

harmony_batch_correction_3_PC <- function(file = sc) {
  sc <- RunPCA(file)
  sc
}

#' Performign the Haromy integration
#' @name harmony_batch_correction_4_Harmony
#' @description
#' This function facilitates the integration of multiple Seurat objects using the Harmony algorithm for batch correction.
#'
#' @param file Merged Seurat object.
#' @param selected_column_for_reduction The column used for reducing batch effects, e.g., 'orig.ident', 'Sample_name'.
#' @param Maximum_PC_to_use Maximum number of principal components (PCs) to use. Default is 30.
#' @param resolution_of_clusters Resolution parameter for clustering. Default is 0.5.
#'
#' @details
#' Harmony is a method for removing batch effects from single-cell RNA-seq data by aligning subpopulations of cells across different batches.
#' This function applies Harmony to the merged Seurat object provided, using the specified column for batch identification.
#'
#' @import Seurat
#' @import harmony
#'
#' @return A Seurat object after Harmony integration, UMAP dimensionality reduction, neighbor finding, and clustering.
#'
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
