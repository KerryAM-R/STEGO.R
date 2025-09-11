#' @name custom_database_scGATE
#' @param repo_path.v path to one custom scGATE structure.
#' @export

custom_db_scGATE <- function(repo_path.v) {
  allfiles <- list.files(repo_path.v, recursive = TRUE,full.names = T)
  modelfiles <- grep("scGate_Model.tsv", allfiles, value = TRUE)
  master.table = read.table(paste0(repo_path.v,"/","master_table.tsv"),sep ="\t",header =T)
  df.models.toimpute <- list()
  files.to.impute <- list.files(file.path(repo_path.v),"_scGate_Model.tsv")

  for(f in 1:length(files.to.impute)){
    model.name <- strsplit(files.to.impute,"_scGate_Model.tsv")[[f]][1]
    model.name
    df.models.toimpute[[model.name]] <- suppressWarnings(read.table(paste0(repo_path.v,"/",files.to.impute[f]),sep ="\t",header =T))
  }

  for(f in 1:length(files.to.impute)){
    model.name <- strsplit(files.to.impute,"_scGate_Model.tsv")[[f]][1]
    message(model.name)
    df.models.toimpute[[model.name]] <- merge(df.models.toimpute[[model.name]][1:3],master.table,by=c("name"),sort = F)
  }
  df.models.toimpute
}
