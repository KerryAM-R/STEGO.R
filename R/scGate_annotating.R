#' scGate_annotating
#' @name scGate_annotating
#' @description
#' This function is to aid in annotating T cells expreiments with or without TCR-seq. There are six available T cell models including the T cell functions, generic annotations, exhausted, senescence, Th1_cytokines and cellular cycling (cell division) and a TCR-seq model.
#'
#' @param TcellFunction T or F; set to T if you want to include the current T cell model
#' @param file Seurat object file. This requires the file to have the scaled data available for annotation purposes.
#' @param generic generic annotations to identify T cells vs other immune cells.
#' @param exhausted T cell based stress models of exhaustion
#' @param senescence T cell based stress models of Senescence
#' @param Th1_cytokines T cell based stress models of IFNG and TNF
#' @param cycling uses TOP2A and MKI67
#' @param TCRseq Uses the meta data TCR-seq to call MAIT, iNKT, potential CD1 restricted, gd T cells and ab T cells
#' @param threshold Set the scGate threshold; default is 0.2 for full models and 0.5 for focused models.
#' @param reductionType Chose the time of dimensional reduction to use; default = harmony
#' @param chunk_size This is the total number of cells to perform the annotation model on. This is to prevent the out of memory issue when annotating. The cells are randomly selected. The same random seed is selected to ensure that random selection is consistent if this has to be repeated. The default is 50,000 cells (we recommend a maximum of 100,000 per loop). However, fewer cells may be required if your RAM is <32 Gb, to ensure each loop is completed.
#' @export

scGate_annotating <- function (file = file, TcellFunction = FALSE, generic = FALSE, exhausted = FALSE,
                               senescence = FALSE, cycling = FALSE, Th1_cytokines = FALSE, TCRseq = FALSE,
                               threshold = 0.2, reductionType = "harmony", chunk_size = 5000)
{
  set.seed(123) # Set a specific seed value, such as 123
  source(system.file("scGATE", "custom_df_scGATE.R", package = "STEGO.R"))
  require(Seurat)
  require(scGate)
  sc <- file
  len.obj <- dim(sc@meta.data)[1]
  threshold_scGate <- threshold

  # Define function to apply scGate to a chunk of data
  apply_scGate_to_chunk <- function(sc_chunk, models_list, threshold_scGate, reductionType) {
    sc_chunk <- scGate(sc_chunk, model = models_list,
                       pos.thr = threshold_scGate,
                       neg.thr = threshold_scGate,
                       nfeatures = ncol(sc_chunk),
                       reduction = reductionType,
                       min.cells = 1)
    sc_chunk@meta.data <- sc_chunk@meta.data[!grepl("_UCell", names(sc_chunk@meta.data))]
    sc_chunk@meta.data <- sc_chunk@meta.data[!grepl("is.pure_", names(sc_chunk@meta.data))]

    return(sc_chunk)
  }

  # Check if any function is selected
  if (!TcellFunction & !generic & !exhausted & !senescence & !cycling & !Th1_cytokines & !TCRseq) {
    message("No function selected. Please set at least one of the function flags to TRUE.")
  }

  # Store UMAP and scaling
  umap <- Embeddings(sc, reduction = "umap")
  harmony <- Embeddings(sc, reduction = "harmony")

  # Initialize a list to store results
  merged_sc_list <- list()

  # Calculate total number of cells and number of chunks
  total_cells <- ncol(sc)
  num_chunks <- ceiling(total_cells / chunk_size)

  # Print total number of cells and number of chunks
  cat("Total number of cells:", total_cells, "\n")
  cat("Number of chunks:", num_chunks, "\n")

  # Loop over chunks
  for (i in 1:num_chunks) {
    # Determine chunk size for this iteration
    if (i == num_chunks) {
      # Last chunk, use the remaining cells
      sc_chunk = sc
      print(sc_chunk)
    } else {
      current_chunk_size <- chunk_size
      sampled_idx <- sample(1:total_cells, current_chunk_size, replace = FALSE)
      # Extract chunk of data
      sc_chunk <- subset(sc, cells = sampled_idx)
      print(sc_chunk)

      # Remove sampled cells from the original object
      sc <- subset(sc, cells = -sampled_idx)
      print(sc)
    }

    # Print number of cells in the chunk
    cat("Chunk", i, "- Number of cells:", ncol(sc_chunk), "\n")

    # Apply scGate based on selected options
    if (TcellFunction) {
      models_list <- custom_db_scGATE(system.file("scGATE", "human/function", package = "STEGO.R"))
      sc_chunk <- apply_scGate_to_chunk(sc_chunk, models_list, threshold_scGate, reductionType)
      sc_chunk@meta.data$Tcellfunction <- sc_chunk@meta.data$scGate_multi
      print(table(sc_chunk@meta.data$Tcellfunction))

    }
    if (generic) {
      models_list <- custom_db_scGATE(system.file("scGATE", "human/generic", package = "STEGO.R"))
      sc_chunk <- apply_scGate_to_chunk(sc_chunk, models_list, threshold_scGate, reductionType)
      sc_chunk@meta.data$generic <- sc_chunk@meta.data$scGate_multi
      print(table(sc_chunk@meta.data$generic))
    }
    if (exhausted) {
      models_list <- custom_db_scGATE(system.file("scGATE", "human/exhausted", package = "STEGO.R"))
      sc_chunk <- apply_scGate_to_chunk(sc_chunk, models_list, threshold_scGate, reductionType)
      sc_chunk@meta.data$exhausted <- sc_chunk@meta.data$scGate_multi
      print(table(sc_chunk@meta.data$exhausted))
    }
    if (senescence) {
      models_list <- custom_db_scGATE(system.file("scGATE", "human/senescence", package = "STEGO.R"))
      sc_chunk <- apply_scGate_to_chunk(sc_chunk, models_list, threshold_scGate, reductionType)
      sc_chunk@meta.data$senescence <- sc_chunk@meta.data$scGate_multi
      print(table(sc_chunk@meta.data$senescence))
    }
    if (Th1_cytokines) {
      models_list <- custom_db_scGATE(system.file("scGATE", "human/Th1_cytokines", package = "STEGO.R"))
      sc_chunk <- apply_scGate_to_chunk(sc_chunk, models_list, threshold_scGate, reductionType)
      sc_chunk@meta.data$Th1_cytokines <- sc_chunk@meta.data$scGate_multi
      print(table(sc_chunk@meta.data$Th1_cytokines))
    }
    if (cycling) {
      models_list <- suppressWarnings(custom_db_scGATE(system.file("scGATE", "human/cycling", package = "STEGO.R")))
      sc_chunk <- apply_scGate_to_chunk(sc_chunk, models_list, threshold_scGate, reductionType)
      sc_chunk@meta.data$cycling <- sc_chunk@meta.data$scGate_multi
      print(table(sc_chunk@meta.data$cycling))
    }
    if (TCRseq) {
      sc_chunk@meta.data$TCRseq <- ifelse(sc_chunk@meta.data$vj_gene_AG == "TRAV1-2.TRAJ33", "MAIT",
                                                  ifelse(sc_chunk@meta.data$vj_gene_AG == "TRAV1-2.TRAJ12", "MAIT",
                                                         ifelse(sc_chunk@meta.data$vj_gene_AG == "TRAV1-2.TRAJ23", "MAIT",
                                                                ifelse(sc_chunk@meta.data$vj_gene_AG == "TRAV1-2.TRAJ9", "CD1b-restricted(poss)",
                                                                       ifelse(sc_chunk@meta.data$vj_gene_AG == "TRAV10.TRAJ18", "iNKT",
                                                                              ifelse(sc_chunk@meta.data$v_gene_BD == "TRBV4-1" & sc_chunk@meta.data$v_gene_AG == "TRAV17", "CD1b-restricted(poss)",
                                                                                     ifelse(sc_chunk@meta.data$v_gene_BD == "TRBV4-1" & sc_chunk@meta.data$v_gene_AG != "TRAV17", "CD1c-restricted(poss)",
                                                                                            ifelse(sc_chunk@meta.data$chain_AG == "TRG" & sc_chunk@meta.data$chain_BD == "TRD", "gd T cell",
                                                                                                   ifelse(sc_chunk@meta.data$chain_AG == "TRA" & sc_chunk@meta.data$chain_BD == "TRB", "ab T cell", "")))))))))

      print(table(sc_chunk@meta.data$TCRseq))

    }
    sc_chunk@meta.data <- sc_chunk@meta.data[!grepl("scGate_multi", names(sc_chunk@meta.data))]



    # Store result in list
    merged_sc_list[[length(merged_sc_list) + 1]] <- list(sc_chunk = sc_chunk, umap = umap, harmony = harmony)
  }

  # Loop to merge pairs of Seurat objects with arrays/matrices until only one is left

  # Loop to merge pairs of Seurat objects until only one is left
  while (length(merged_sc_list) > 1) {
    temp_list <- list()

    for (i in seq(1, length(merged_sc_list), by = 2)) {


      if (i < length(merged_sc_list)) {
        # Merge two Seurat objects
        merged_names <- paste(names(merged_sc_list)[i], names(merged_sc_list)[i + 1], sep = " and ")
        new_name <- paste0("Merged_", i)
        message("Merging ", merged_names, " to create ", new_name)
        merged_object <- merge(x = merged_sc_list[[i]]$sc_chunk, y = merged_sc_list[[i + 1]]$sc_chunk, merge.data = TRUE)
        # scalling <- suppressWarnings(merge(x = merged_sc_list[[i]]$scaling, y = merged_sc_list[[i + 1]]$scaling, merge.data = TRUE))
        temp_list[[new_name]] <- merged_object

        # Calculate size of merged object
        merged_size <- object.size(merged_object)
        message("Size of ", new_name, ": ", round(merged_size[1] / 1000^3, 1), " Gb")
      } else {
        # If there's an odd number of objects, keep the last one as is
        temp_list[[names(merged_sc_list)[i]]] <- merged_sc_list[[i]]
      }
    }

    merged_sc_list <- temp_list
  }

  # Extract the merged Seurat object
  merged_sc <- merged_sc_list[[1]]

  join_sc <- JoinLayers(merged_sc)
  join_sc@reductions$umap <- CreateDimReducObject(embeddings = umap, key = 'UMAP_', assay = 'RNA')
  join_sc@reductions$harmony <- CreateDimReducObject(embeddings = harmony, key = 'harmony_', assay = 'RNA')
  # join_sc_list <- list(sc = join_sc, umap = umap, harmony = harmony)
  return(join_sc)
}
