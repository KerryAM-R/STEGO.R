#' scGate Annotating Function
#'
#' @description
#' This function aids in annotating T cell experiments with or without TCR-seq data. It includes six available T cell models: T cell functions, generic annotations, immune checkpoint, senescence, Th1 cytokines, cellular cycling (cell division), and a TCR-seq model.
#'
#' @param TcellFunction Logical; set to TRUE if you want to include the current T cell model.
#' @param file Seurat object file. This requires the file to have the scaled data available for annotation purposes.
#' @param Threshold_test Logical; testing the scGate threshold for identifying the sub-populations. Recommended to use with BD Rhapsody immune panel.
#' @param signature_for_testing Character vector; signature for testing purposes. Default is c("CD8A", "CD8B").
#' @param threshold Numeric; set the scGate threshold. Default is 0.2 for full models and 0.55 for focused models from BD Rhapsody.
#' @param immune_checkpoint Logical; include T cell-based stress models of exhaustion.
#' @param senescence Logical; include T cell-based stress models of senescence.
#' @param Th1_cytokines Logical; include T cell-based stress models of IFNG and TNF.
#' @param cycling Logical; include models using TOP2A and MKI67.
#' @param TCRseq Logical; uses the metadata from TCR-seq to call MAIT, iNKT, potential CD1-restricted, gd T cells, and ab T cells.
#' @param reductionType Character; type of dimensional reduction to use. Default is "harmony".
#' @param Version Character; Seurat version. Either "V4", "V5", or "python". If converted from anndata to Seurat, use the python version.
#' @param chunk_size Integer; total number of cells to perform the annotation model on. Default is 50,000 cells.
#' @param output_dir Character; directory to store the outputs of this process. If running multiple times, it's best to change the name.
#' @import Seurat
#' @import scGate
#' @export
#'
#' @return Seurat object with annotations
#'

scGate_annotating <- function (file = file,
                               Threshold_test = FALSE,
                               signature_for_testing = c("CD8A","CD8B"),
                               threshold = 0.2,
                               TcellFunction = FALSE,
                               immune_checkpoint = FALSE,
                               senescence = FALSE, cycling = FALSE, Th1_cytokines = FALSE, TCRseq = FALSE,
                                reductionType = "harmony", chunk_size = 50000, output_dir = "output",
                               Version = c("V5","V4","python"))
{
  set.seed(123) # Set a specific seed value, such as 123
  source(system.file("scGATE", "custom_df_scGATE.R", package = "STEGO.R"))

  sc <- file
  len.obj <- dim(sc@meta.data)[1]
  threshold_scGate <- threshold
  output_dir = output_dir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Function to save table to a text file
  save_table <- function(table_data, output_file) {
    write.table(table_data, file = output_file, sep = "\t", row.names = TRUE, col.names = TRUE)
  }

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

  if  (Threshold_test) {
    message("Testing scGate threshold")
  } else if (!TcellFunction & !immune_checkpoint & !senescence & !cycling & !Th1_cytokines & !TCRseq) {
    message("No function selected. Please set at least one of the function flags to TRUE.")
  }


  if(Threshold_test) {
    my_scGate_model <- gating_model(name = "test", signature = signature_for_testing)

    sc <- scGate(sc, model = my_scGate_model,
                 nfeatures = ncol(sc),
                 neg.thr = threshold_scGate,
                 pos.thr = threshold_scGate,
                 reduction = reductionType, min.cells = 1)

    DimPlot(sc) | FeaturePlot(sc, features = signature_for_testing)
  } else {


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
      total_cells <- ncol(sc)
      # Determine chunk size for this iteration
      if (i == num_chunks) {
        # Last chunk, use the remaining cells
        sc_chunk = sc
        print(sc_chunk)
      } else {
        current_chunk_size <- chunk_size
        sampled_idx <- sample(1:total_cells, current_chunk_size)
        # Extract chunk of data
        sc_chunk <- subset(sc, cells = sampled_idx)
        print(sc_chunk)

        # Remove sampled cells from the original object
        sc <- subset(sc, cells = -sampled_idx)
        print(ncol(sc))
      }

      # Print number of cells in the chunk
      cat("Chunk", i,"of", num_chunks, "- Number of cells:", ncol(sc_chunk), "\n")

      # Apply scGate based on selected options
      if (TcellFunction) {
        models_list <- custom_db_scGATE(system.file("scGATE", "human/function", package = "STEGO.R"))
        sc_chunk <- apply_scGate_to_chunk(sc_chunk, models_list, threshold_scGate, reductionType)
        sc_chunk@meta.data$Tcellfunction <- sc_chunk@meta.data$scGate_multi

        tcell_function_table <- table(sc_chunk@meta.data$Tcellfunction)
        # print(tcell_function_table)
        save_table(tcell_function_table, paste0(output_dir,"/", "Tcellfunction_table_chunk_", i, ".txt"))


      }
      # if (generic) {
      #   models_list <- custom_db_scGATE(system.file("scGATE", "human/generic", package = "STEGO.R"))
      #   sc_chunk <- apply_scGate_to_chunk(sc_chunk, models_list, threshold_scGate, reductionType)
      #   sc_chunk@meta.data$generic <- sc_chunk@meta.data$scGate_multi
      #   # Print and save meta data
      #   tcell_function_table <- table(sc_chunk@meta.data$generic)
      #   # print(tcell_function_table)
      #   save_table(tcell_function_table, paste0(output_dir,"/", "generic_table_chunk_", i, ".txt"))
      # }
      if (immune_checkpoint) {
        models_list <- custom_db_scGATE(system.file("scGATE", "human/immune_checkpoint", package = "STEGO.R"))
        sc_chunk <- apply_scGate_to_chunk(sc_chunk, models_list, threshold_scGate, reductionType)
        sc_chunk@meta.data$immune_checkpoint <- sc_chunk@meta.data$scGate_multi

        tcell_function_table <- table(sc_chunk@meta.data$immune_checkpoint)
        # print(tcell_function_table)
        save_table(tcell_function_table, paste0(output_dir,"/", "exhausted_table_chunk_", i, ".txt"))

      }
      if (senescence) {
        models_list <- custom_db_scGATE(system.file("scGATE", "human/senescence", package = "STEGO.R"))
        sc_chunk <- apply_scGate_to_chunk(sc_chunk, models_list, threshold_scGate, reductionType)
        sc_chunk@meta.data$senescence <- sc_chunk@meta.data$scGate_multi
        # Print and save meta data
        tcell_function_table <- table(sc_chunk@meta.data$senescence)
        # print(tcell_function_table)
        save_table(tcell_function_table, paste0(output_dir,"/", "senescence_table_chunk_", i, ".txt"))
      }
      if (Th1_cytokines) {
        models_list <- custom_db_scGATE(system.file("scGATE", "human/Th1_cytokines", package = "STEGO.R"))
        sc_chunk <- apply_scGate_to_chunk(sc_chunk, models_list, threshold_scGate, reductionType)
        sc_chunk@meta.data$Th1_cytokines <- sc_chunk@meta.data$scGate_multi

        tcell_function_table <- table(sc_chunk@meta.data$Th1_cytokines)
        # print(tcell_function_table)
        save_table(tcell_function_table, paste0(output_dir,"/", "Th1_cytokines_table_chunk_", i, ".txt"))
      }
      if (cycling) {
        models_list <- suppressWarnings(custom_db_scGATE(system.file("scGATE", "human/cycling", package = "STEGO.R")))
        sc_chunk <- apply_scGate_to_chunk(sc_chunk, models_list, threshold_scGate, reductionType)
        sc_chunk@meta.data$cycling <- sc_chunk@meta.data$scGate_multi

        tcell_function_table <- table(sc_chunk@meta.data$cycling)
        # print(tcell_function_table)
        save_table(tcell_function_table, paste0(output_dir,"/", "cycling_table_chunk_", i, ".txt"))
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

        tcell_function_table <- table(sc_chunk@meta.data$TCRseq)
        # print(tcell_function_table)
        save_table(tcell_function_table, paste0(output_dir,"/", "TCRseq_table_chunk_", i, ".txt"))

      }
      sc_chunk@meta.data <- sc_chunk@meta.data[!grepl("scGate_multi", names(sc_chunk@meta.data))]



      # Store result in list
      merged_sc_list[[length(merged_sc_list) + 1]] <- sc_chunk
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
          merged_object <- merge(x = merged_sc_list[[i]], y = merged_sc_list[[i + 1]], merge.data = TRUE)
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
    if (num_chunks == 1) {
      join_sc <- merged_sc
    } else {

      if (Version == "V4" |Version == "Python") {
        join_sc <- merged_sc

        join_sc@reductions$umap <- CreateDimReducObject(embeddings = umap_reordered, key = 'UMAP_', assay = 'RNA')
        join_sc@reductions$harmony <- CreateDimReducObject(embeddings = harmony_reordered, key = 'harmony_', assay = 'RNA')
        join_sc@assays$RNA@scale.data  <- sc@assays$RNA@scale.data
        join_sc@assays$RNA@var.features <- sc@assays$RNA@var.features

      } else {
        join_sc <- JoinLayers(merged_sc)

        barcode_order <- rownames(join_sc@meta.data)
        umap_reordered <- umap[match(barcode_order, rownames(umap)), ]
        harmony_reordered <- umap[match(barcode_order, rownames(harmony)),]

        join_sc@reductions$umap <- CreateDimReducObject(embeddings = umap_reordered, key = 'UMAP_', assay = 'RNA')
        join_sc@reductions$harmony <- CreateDimReducObject(embeddings = harmony_reordered, key = 'harmony_', assay = 'RNA')
        join_sc@assays$RNA@layers$scale.data  <- sc@assays$RNA@layers$scale.data
        join_sc@assays$RNA@layers$var.features <- sc@assays$RNA@layers$var.features

        Idents(join_sc) <- join_sc@meta.data$seurat_clusters
        # join_sc_list <- list(sc = join_sc, umap = umap, harmony = harmony)


      }
    }

    return(join_sc)
  }
}
