#' scGate_annotating
#' @name scGate_annotating
#' @description
#' This function is to aid in annotating T cells expreiments with or without TCR-seq. There are four avaialbe T cell models including the T cell functions, generic annotations, T cell based stress and cellular cycling
#'
#' @param TcellFunction T or F; set to T if you want to include the current T cell model
#' @param file Seurat object file. This requires the file to have the scaled data available for annotation purposes.
#' @param generic generic annotations to identify T cells vs other immune cells.
#' @param stress T cell based stress models of exhaustion and Senescence
#' @param cycling uses TOP2A and MKI67
#' @param threshold Set the scGate threshold; default is 0.2 for full models and 0.5 for focused models.
#' @param reductionType Chose the time of dimentional reduction to use; default = harmony
#' @export

scGate_annotating <- function(file = file, TcellFunction = F,generic = F, stress = F, cycling = F, threshold = 0.2, reductionType = "harmony") {
  source(system.file("scGATE","custom_df_scGATE.R",package = "STEGO.R"))
  STEGO.R::Load_required_packages()
  sc <- file

  threshold_scGate = threshold
  len <- length(rownames(sc@assays$RNA$scale.data))
  if(TcellFunction) {
    models.list <- custom_db_scGATE(system.file("scGATE","human/function",package = "STEGO.R"))
    sc <- scGate(sc, model = models.list,pos.thr = threshold_scGate,neg.thr = threshold_scGate,nfeatures = len, ncores = 4, reduction = reductionType, min.cells = 1)
    sc@meta.data$Tcellfunction <- sc@meta.data$scGate_multi
    sc@meta.data <- sc@meta.data[!grepl("_UCell",names(sc@meta.data))]
    sc@meta.data <- sc@meta.data[!grepl("is.pure_",names(sc@meta.data))]
    sc@meta.data <- sc@meta.data[!grepl("scGate_multi",names(sc@meta.data))]
  }

  if (generic) {
    models.list <- custom_db_scGATE(system.file("scGATE","human/generic",package = "STEGO.R"))
    sc <- scGate(sc, model = models.list,pos.thr = threshold_scGate,neg.thr = threshold_scGate,nfeatures = len, ncores = 4, reduction = reductionType, min.cells = 1)
    sc@meta.data$generic <- sc@meta.data$scGate_multi
    sc@meta.data <- sc@meta.data[!grepl("_UCell",names(sc@meta.data))]
    sc@meta.data <- sc@meta.data[!grepl("is.pure_",names(sc@meta.data))]
    sc@meta.data <- sc@meta.data[!grepl("scGate_multi",names(sc@meta.data))]
  }

  if (stress)  {
    models.list <- custom_db_scGATE(system.file("scGATE","human/stress",package = "STEGO.R"))
    sc <- scGate(sc, model = models.list,pos.thr = threshold_scGate,neg.thr = threshold_scGate,nfeatures = len, reduction = reductionType, ncores = 4, min.cells = 1)
    sc@meta.data$stress <- sc@meta.data$scGate_multi
    sc@meta.data$stress <- ifelse(grepl("Exhau",sc@meta.data$stress),"Exhausted",sc@meta.data$stress)
    sc@meta.data$stress <- ifelse(grepl("Multi",sc@meta.data$stress),"Exhausted.Senescence",sc@meta.data$stress)
    sc@meta.data <- sc@meta.data[!grepl("_UCell",names(sc@meta.data))]
    sc@meta.data <- sc@meta.data[!grepl("is.pure_",names(sc@meta.data))]
    sc@meta.data <- sc@meta.data[!grepl("scGate_multi",names(sc@meta.data))]
  }

  if (cycling) {
    models.list <- suppressWarnings(custom_db_scGATE(system.file("scGATE","human/cycling",package = "STEGO.R")))
    models.list
    sc <- scGate(sc, model = models.list,pos.thr = threshold_scGate,neg.thr = threshold_scGate,nfeatures = len, ncores = 4, reduction = reductionType, min.cells = 1)
    sc@meta.data$cycling <- sc@meta.data$scGate_multi
    sc@meta.data <- sc@meta.data[!grepl("_UCell",names(sc@meta.data))]
    sc@meta.data <- sc@meta.data[!grepl("is.pure_",names(sc@meta.data))]
    sc@meta.data <- sc@meta.data[!grepl("scGate_multi",names(sc@meta.data))]
  }
  sc

}
