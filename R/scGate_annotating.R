#' log_parameters_scGate
#'
#' @description
#' This function logs the parameters within the scGate annotation function. It should add the method text to the log file.
#'
#' @param log_file name of the log file called "method_file.txt". It will save in the project directory.
#' @param model_to_log annotation models used
#' @param params_to_log Other scGate parameters that were changed from the default.
#'
#' @export
#'
#' @return saves the log of the scGate prameters for the methods section.
#'

log_parameters_scGate <- function(log_file, model_to_log = NULL, params_to_log = NULL) {
  # Open the log_file in append mode
  con <- file(log_file, open = "a")

  # Get the versions of Seurat and STEGO.R
  scGate_version <- packageVersion("scGate")
  stego_version <- packageVersion("STEGO.R")

  # Start with a newline
  writeLines("\n", con)

  if (!is.null(params_to_log)) {
    non_default_params <- sapply(names(params_to_log), function(name) paste(name, "(", params_to_log[[name]], ")", sep = ""))
    non_default_params_str <- paste(non_default_params, collapse = ", ")
  } else {
    non_default_params_str <- ""
  }

  if (!is.null(model_to_log)) {
    models_used <- sapply(names(model_to_log), function(name) paste(name))
    if (length(models_used) > 1) {
      models_used_str <- paste(paste(models_used[-length(models_used)], collapse = ", "), "and", models_used[length(models_used)])
    } else {
      models_used_str <- models_used
    }
  } else {
    models_used_str <- ""
  }

  message <- sprintf("The merged object underwent annotating using the semi-supervised process with scGate (version %s)[ref]. The annotation models created in STEGO.R included %s. Additionally, the following scGate parameters were changed from the default values: %s. from the detu", scGate_version, models_used_str, non_default_params_str)

  if (grepl("chunk_size", message)) {
    message <- gsub("chunk_size", "number of cells per chunk ", message)
  }

  if (length(unique(names(params_to_log))) == 1) {
    message <- gsub("parameters", "parameter", message)
  }

  print(message)

  # Write the message to the log file
  writeLines(message, con)

  # Close the file connection
  close(con)
}

#' scGate Annotating Function
#'
#' @description
#' This function aids in annotating T cell experiments with or without TCR-seq data. It includes six available T cell models: T cell functions, generic annotations, immune checkpoint, senescence, Th1 cytokines, cellular cycling (cell division), and a TCR-seq model. Currently this function is for human only.
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

scGate_annotating <- function (
                               file = file,
                               Threshold_test = FALSE,
                               signature_for_testing = c("CD8A","CD8B"),
                               threshold = 0.25,
                               TcellFunction = FALSE,
                               SimpleFunction = FALSE,
                               immune_checkpoint = FALSE,
                               senescence = FALSE,
                               cycling = FALSE,
                               Th1_cytokines = FALSE,
                               TCRseq = FALSE,
                               reductionType = "harmony",
                               chunk_size = 50000,
                               output_dir = "output",
                               Version = c("V5","V4","python")
                               )
{
  set.seed(123) # Set a specific seed value, such as 123
  source(system.file("scGATE", "custom_df_scGATE.R", package = "STEGO.R"))

  default_models <- list(
    TcellFunction = F,
    SimpleFunction = F,
    immune_checkpoint = F,
    senescence = F,
    cycling = FALSE,
    Th1_cytokines = F,
    TCRseq = F
  )
  models <- list(
    TcellFunction = TcellFunction,
    SimpleFunction = SimpleFunction,
    immune_checkpoint = immune_checkpoint,
    senescence = senescence,
    cycling = cycling,
    Th1_cytokines = Th1_cytokines,
    TCRseq = TCRseq
  )
print(models)

  models_to_log <- list()
  for (param_name in names(model)) {
    if (!identical(model[[param_name]], default_models[[param_name]])) {
      models_to_log[[param_name]] <- model[[param_name]]
    }
  }
  models_to_log

  default_params <- list(
    threshold = 0.25,
    reductionType = "harmony",
    chunk_size = 50000
  )

  current_params <- list(
    threshold = threshold,
    reductionType = reductionType,
    chunk_size = chunk_size
  )

  params_to_log <- list()
  for (param_name in names(current_params)) {
      if (!identical(current_params[[param_name]], default_params[[param_name]])) {
        params_to_log[[param_name]] <- current_params[[param_name]]
      }
    }
  params_to_log

  log_file <- "log_file.txt"
  params_to_log

  if (length(params_to_log) > 0 && length(models_to_log)>0) {
    log_parameters_scGate(log_file, models_to_log,params_to_log)
  } else {
    log_parameters_scGate(log_file)
  }




#### adding in the models -----
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

        sc_chunk@meta.data$majorTCF <-ifelse(grepl("CD4",sc_chunk@meta.data$Tcellfunction),"CD4",
                                             ifelse(grepl("CD8ab",sc_chunk@meta.data$Tcellfunction),"CD8ab",
                                                    ifelse(grepl("CD8aa",sc_chunk@meta.data$Tcellfunction),"CD8ab",
                                                           ifelse(grepl("DN",sc_chunk@meta.data$Tcellfunction),"DN","other"
                                                           ))))

        tcell_function_table <- table(sc_chunk@meta.data$Tcellfunction)
        save_table(tcell_function_table, paste0(output_dir,"/", "Tcellfunction_table_chunk_", i, ".txt"))


      }
      if (SimpleFunction) {
        models_list <- custom_db_scGATE(system.file("scGATE", "human/simple_functions", package = "STEGO.R"))
        sc_chunk <- apply_scGate_to_chunk(sc_chunk, models_list, threshold_scGate, reductionType)
        sc_chunk@meta.data$TSimpleFunction <- sc_chunk@meta.data$scGate_multi
        sc_chunk@meta.data$majorSTF <-ifelse(grepl("CD8[.]CD4",sc_chunk@meta.data$TSimpleFunction),"DP",
                                             ifelse(grepl("CD4",sc_chunk@meta.data$TSimpleFunction),"CD4",
                                                    ifelse(grepl("Th",sc_chunk@meta.data$TSimpleFunction),"CD4",
                                                           ifelse(grepl("Treg",sc_chunk@meta.data$TSimpleFunction),"CD4",
                                                                  ifelse(grepl("CD8",sc_chunk@meta.data$TSimpleFunction),"CD8",
                                                                         ifelse(grepl("CD8",sc_chunk@meta.data$TSimpleFunction),"CD8",
                                                                                ifelse(grepl("DN",sc_chunk@meta.data$TSimpleFunction),"DN","other"
                                                                                )))))))

        tcell_function_table <- table(sc_chunk@meta.data$majorSTF)
        save_table(tcell_function_table, paste0(output_dir,"/", "TSimpleFunction_table_chunk_", i, ".txt"))
      }
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

    # saveRDS(merged_sc_list,"Annotation_chunks_obj.rds")
    # saveRDS(umap,"umap.rds")
    # saveRDS(harmony,"harmony.rds")
    print(names(merged_sc_list))
    merged_sc_list2 <- merged_sc_list
    names(merged_sc_list2) <-  paste0(rep("List_obj_name",num_chunks),1:num_chunks)

    # Loop to merge pairs of Seurat objects until only one is left
    # Run the loop until only one object remains in the list
    while (length(merged_sc_list2) > 1) {
      temp_list <- list()
      counter <- 1  # Reset the counter at the beginning of each iteration
      for (i in seq(1, length(merged_sc_list2), by = 2)) {
        if (i < length(merged_sc_list2)) {
          merged_names <- paste(names(merged_sc_list2)[i],
                                names(merged_sc_list2)[i + 1], sep = " and ")
          new_name <- paste0("Merged_", counter)  # Generate sequential name
          message("Merging ", merged_names, " to create ", new_name)
          merged_object <- merge(x = merged_sc_list2[[i]],
                                 y = merged_sc_list2[[i + 1]], merge.data = TRUE)

          if(Version == "V5") {
            message("Joining layers of ",new_name)
            merged_object <- JoinLayers(merged_object,  assay = "RNA")
          }

          temp_list[[new_name]] <- merged_object
          merged_size <- object.size(merged_object)
          message("Size of ", new_name, ": ", round(merged_size[1]/1000^3, 1), " Gb")
          counter <- counter + 1  # Increment the counter
        } else {
          # If there's an odd number of elements, keep the last one
          temp_list[[names(merged_sc_list2)[i]]] <- merged_sc_list2[[i]]
        }
      }
      # Replace the original list with the merged list
      merged_sc_list2 <- temp_list
    }

    # Extract the merged Seurat object
    merged_sc <- merged_sc_list2[[1]]

    message("Adding back in the UMAP and Harmony data")
    if (num_chunks == 1) {
      join_sc <- merged_sc
    } else {
      if (Version == "V4" |Version == "Python") {
        join_sc <- merged_sc

        join_sc@reductions$umap <- CreateDimReducObject(embeddings = umap_reordered, key = 'UMAP_', assay = 'RNA')
        join_sc@reductions$harmony <- CreateDimReducObject(embeddings = harmony_reordered, key = 'harmony_', assay = 'RNA')
        join_sc@assays$RNA@scale.data  <- file@assays$RNA@scale.data
        join_sc@assays$RNA@var.features <- file@assays$RNA@var.features

      } else {
        message("Joining layers")
        join_sc <- JoinLayers(merged_sc)
        message("reordering and adding back in the umap & harmony data")
        barcode_order <- rownames(join_sc@meta.data)
        umap_reordered <- umap[match(barcode_order, rownames(umap)), ]
        harmony_reordered <- umap[match(barcode_order, rownames(harmony)),]

        join_sc@reductions$umap <- CreateDimReducObject(embeddings = umap_reordered, key = 'UMAP_', assay = 'RNA')
        join_sc@reductions$harmony <- CreateDimReducObject(embeddings = harmony_reordered, key = 'harmony_', assay = 'RNA')

        message("Adding in seurat_clusters to the Idents column")
        join_sc@meta.data$seurat_clusters <- factor(join_sc@meta.data$seurat_clusters,levels = 0:100)
        print(join_sc)
        Idents(join_sc) <- join_sc@meta.data$seurat_clusters
        # join_sc_list <- list(sc = join_sc, umap = umap, harmony = harmony)
      }
    }
    join_sc
    return(join_sc)
  }
}
