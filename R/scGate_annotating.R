#' scGate_annotating
#' @name scGate_annotating
#' @description
#' This function is to aid in annotating T cells expreiments with or without TCR-seq. There are four avaialbe T cell models including the T cell functions, generic annotations, T cell based stress and cellular cycling
#'
#' @param TcellFunction T or F; set to T if you want to include the current T cell model
#' @param file Seurat object file. This requires the file to have the scaled data available for annotation purposes.
#' @param generic generic annotations to identify T cells vs other immune cells.
#' @param exhausted T cell based stress models of exhaustion
#' @param senescence T cell based stress models of Senescence
#' @param Th1_cytokines T cell based stress models of IFNG and TNF
#' @param cycling uses TOP2A and MKI67
#' @param threshold Set the scGate threshold; default is 0.2 for full models and 0.5 for focused models.
#' @param reductionType Chose the time of dimentional reduction to use; default = harmony
#' @export

scGate_annotating <- function(file = file, TcellFunction = F,generic = F, exhausted = F, senescence =F, cycling = F, Th1_cytokines = F, threshold = 0.2, reductionType = "harmony") {
  source(system.file("scGATE","custom_df_scGATE.R",package = "STEGO.R"))
  require(Seurat)
  require(scGate)

  sc <- file
  len.obj <- dim(sc@meta.data)[1]

  # if (len.obj>120000) {
    # message(paste0("Cannot annotate larger >120K cells. Down sampling required"))
    sc
  # } else {

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

  if (exhausted)  {
    models.list <- custom_db_scGATE(system.file("scGATE","human/exhausted",package = "STEGO.R"))
    sc <- scGate(sc, model = models.list,pos.thr = threshold_scGate,neg.thr = threshold_scGate,nfeatures = len, reduction = reductionType, ncores = 4, min.cells = 1)
    sc@meta.data$exhausted <- sc@meta.data$scGate_multi
    sc@meta.data <- sc@meta.data[!grepl("_UCell",names(sc@meta.data))]
    sc@meta.data <- sc@meta.data[!grepl("is.pure_",names(sc@meta.data))]
    sc@meta.data <- sc@meta.data[!grepl("scGate_multi",names(sc@meta.data))]
  }

  if (senescence)  {
    models.list <- custom_db_scGATE(system.file("scGATE","human/senescence",package = "STEGO.R"))
    sc <- scGate(sc, model = models.list,pos.thr = threshold_scGate,neg.thr = threshold_scGate,nfeatures = len, reduction = reductionType, ncores = 4, min.cells = 1)
    sc@meta.data$senescence <- sc@meta.data$scGate_multi
    sc@meta.data <- sc@meta.data[!grepl("_UCell",names(sc@meta.data))]
    sc@meta.data <- sc@meta.data[!grepl("is.pure_",names(sc@meta.data))]
    sc@meta.data <- sc@meta.data[!grepl("scGate_multi",names(sc@meta.data))]
  }

  if (Th1_cytokines)  {
    models.list <- custom_db_scGATE(system.file("scGATE","human/Th1_cytokines",package = "STEGO.R"))
    sc <- scGate(sc, model = models.list,pos.thr = threshold_scGate,neg.thr = threshold_scGate,nfeatures = len, reduction = reductionType, ncores = 4, min.cells = 1)
    sc@meta.data$Th1_cytokines <- sc@meta.data$scGate_multi
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
  # }
}
