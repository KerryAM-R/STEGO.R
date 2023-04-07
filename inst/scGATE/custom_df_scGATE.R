#' @name custom_database_scGATE
#' @param repo_path.v path to one custom scGATE structure.
#' @export

custom_db_scGATE <- function(repo_path.v) {
allfiles <- list.files(repo_path.v, recursive = TRUE)
modelfiles <- grep("scGate_Model.tsv", allfiles, value = TRUE)
uniq_dirs <- sort(unique(dirname(modelfiles)))
model_db <- list()
sub <- strsplit(uniq_dirs, split = "/")[[1]]
model_path <- file.path(repo_path.v, uniq_dirs)

master.table = read.table(paste0(model_path,"/","master_table.tsv"),sep ="\t",header =T)
df.models.toimpute <- list()
files.to.impute <- list.files(file.path(model_path),"_scGate_Model.tsv")
    for(f in 1:length(files.to.impute)){
      model.name <- strsplit(files.to.impute,"_scGate_Model.tsv")[[f]][1]
      model.name
      read.table(file.path(model_path,files.to.impute[[f]][1]),sep ="\t",header =T)
      df.models.toimpute[[model.name]] <- read.table(paste0(model_path,"/",files.to.impute[f]),sep ="\t",header =T)
    }
    for(f in 1:length(files.to.impute)){
      model.name <- strsplit(files.to.impute,"_scGate_Model.tsv")[[f]][1]
      model.name
      read.table(file.path(model_path,files.to.impute[[f]][1]),sep ="\t",header =T)

      df.models.toimpute[[model.name]] <- merge(read.table(paste0(model_path,"/",files.to.impute[f]),sep ="\t",header =T)[1:3],master.table,by=c("name"))

    }
    df.models.toimpute
}
