#' Automated preprocessing
#' @name preprocessing
#' @description
#' This function performs the pre-processing equivalent to STEP 1. This requires the user to name their folders indiv_group. Please use the _ to separate the values, this particularly important for the TCR_Explore as it will separate out the indiv from the group column.
#' Each folder should contain one of each that have the following words: barcode, features, matrix and contig
#'
#' @param downloadSeurat T or F; set to T if you want to include the current T cell model
#' @param downloadTCRex Seurat object file. This requires the file to have the scaled data available for annotation purposes.
#' @param downloadClusTCR generic annotations to identify T cells vs other immune cells.
#' @param downloadTCR_Explore T cell based stress models of exhaustion
#' @param csv_contig_file T cell based stress models of Senescence
#' @param main_directory This is the location of the raw data from 10x outputs
#' @export

preprocessing <- function (downloadTCRex=F, downloadClusTCR = F, downloadTCR_Explore = F, downloadSeurat = F, csv_contig_file = "csv", main_directory = "0_RAW_files/") {
  require(magrittr)
  require(plyr)
  require(dplyr)
  require(lubridate)
  main_directory = "0_RAW_files/"
  main_directory <- main_directory
  main_folders <- list.files(paste(main_directory))
  main_folders
  # Subdirectory within the main directory
  num <- length(main_folders)


  # Get list of files in the subdirectory
  for (i in 2:num) {

    sub_directory <- main_folders[i]
    print(sub_directory)
    files <- list.files(paste(main_directory, sub_directory, sep = ""),full.names = T)
    print(files)

    if (length(files[grepl("barcode",files)]) == 1 && length(files[grepl("contig",files)]) == 1) {
      message("files are present")
    } else {
      message("reformat directory to have Indiv_group name one of each file barcode, features, matrix and contig")
      next
    }
    # loading the files -----

    barcode <- read.table(files[grepl("barcode",files)])
    dim(barcode)
    print("barcode found")
    features <- read.table(files[grepl("feature",files)])
    print("features found")

    if (length(features)>0) {
      mat <- Matrix::readMM(files[grepl("matrix",files)])
      print("maxtrix found")
    }

    TCR <- read.csv(files[grepl("contig",files)])


    if (length(TCR)>0) {
      message("Can create the meta.data, TCRex, Clustering files and TCR_Explore")
    } else {
      message("missing TCR contig file, please restructure the folder.")
      next
    }
    #### TCRex -----
    if (downloadTCRex) {
      contigs <- TCR
      if (csv_contig_file == "tsv") {
        names(contigs) <- gsub("call", "gene", names(contigs))
        names(contigs) <- gsub("junction_aa", "cdr3", names(contigs))
      } else {

      }

      contigs2 <- contigs[, names(contigs) %in% c("v_gene", "j_gene", "cdr3")]

      names(contigs2)[names(contigs2) %in% c("cdr3")] <- "CDR3_beta"
      if (nrow(contigs2[-c(grep("[*]", contigs2$CDR3_beta)), ] > 0)) {
        contigs2 <- contigs2[-c(grep("[*]", contigs2$CDR3_beta)), ]
      }
      names(contigs2)[names(contigs2) %in% c("v_gene")] <- "TRBV_gene"
      names(contigs2)[names(contigs2) %in% c("j_gene")] <- "TRBJ_gene"
      contigs2 <- contigs2 %>%
        select(TRBV_gene, CDR3_beta, everything())


      contigs2$CDR3_beta <- gsub("^$", "None", contigs2$CDR3_beta)

      if (nrow(contigs2[-c(grep("None", contigs2$CDR3_beta)), ] > 0)) {
        contigs2 <- subset(contigs2, contigs2$CDR3_beta != "None")
      }

      contigs2$TRBJ_gene <- gsub("^$", "None", contigs2$TRBJ_gene)

      if (nrow(contigs2[-c(grep("None", contigs2$TRBJ_gene)), ] > 0)) {
        contigs2 <- subset(contigs2, contigs2$TRBJ_gene != "None")
      }

      contigs2 <- contigs2
      contigs2$TRBV_gene <- gsub("[*]0.", "", contigs2$TRBV_gene)
      contigs2$TRBJ_gene <- gsub("[*]0.", "", contigs2$TRBJ_gene)
      contigs2$cloneCount <- 1
      calls_TCR_paired.fun3 <- ddply(contigs2, names(contigs2)[-c(4)], numcolwise(sum))
      calls_TCR_paired.fun3 <- calls_TCR_paired.fun3[grepl("TRB", calls_TCR_paired.fun3$TRBV_gene), ]
      file_name_tcrex <- paste0("1_TCRex/",sub_directory,"_TCRex_",today(),".tsv")
      # file_name_tcrex
      write.table(calls_TCR_paired.fun3, file_name_tcrex, row.names = F, sep = "\t", quote = F)
      message(paste0("Downloaded TCRex for ",sub_directory))
    } else {
      message("Did not download TCRex, switch to TRUE to run.")
    }

    if (downloadClusTCR) {
      contigs <- TCR
      if (csv_contig_file == "tsv") {
        names(contigs) <- gsub("call", "gene", names(contigs))
        names(contigs) <- gsub("junction_aa", "cdr3", names(contigs))
      } else {

      }

      contigs2 <- contigs[names(contigs) %in% c("v_gene", "cdr3")]
      contigs2$Sample_Name <- sub_directory

      names(contigs2)[names(contigs2) %in% c("cdr3")] <- "junction_aa"
      names(contigs2)[names(contigs2) %in% c("v_gene")] <- "v_call"

      if (nrow(contigs2[-c(grep("[*]", contigs2$junction_aa)), ] > 0)) {
        contigs2 <- contigs2[-c(grep("[*]", contigs2$junction_aa)), ]
      }

      contigs2$junction_aa <- gsub("^$", "None", contigs2$junction_aa)

      if (nrow(contigs2[-c(grep("None", contigs2$junction_aa)), ] > 0)) {
        contigs2 <- contigs2[-c(grep("None", contigs2$junction_aa)), ]
      }

      df <- contigs2[!duplicated(contigs2[, c("v_call", "junction_aa")]), ]
      # if (input$chain_clusTCR2_10x == "AG") {
      TRAG <-  (rbind(df[grep("TRAV", df$v_call), ], df[grep("TRGV", df$v_call), ]))
      if(nrow(TRAG)>0) {
        file_name_clust <- paste0("1_ClusTCR/","AG_",sub_directory,"_clusTCR_10x_",today(),".csv")
        # file_name_tcrex
        write.csv(TRAG, file_name_clust, row.names = F)
      }
      # } else if (input$chain_clusTCR2_10x == "IgH") {
      IGH <-   rbind(df[grep("IGH", df$v_call), ])
      if(nrow(IGH)>0) {
        file_name_clust <- paste0("1_ClusTCR/","IgH_",sub_directory,"_clusTCR_10x_",today(),".csv")
        # file_name_tcrex
        write.csv(IGH, file_name_clust, row.names = F)
      }
      # } else if (input$chain_clusTCR2_10x == "IgLK") {
      IGLK <- (rbind(df[grep("IGL", df$v_call), ], df[grep("IGK", df$v_call), ]))
      if(nrow(IGLK)>0) {
        file_name_clust <- paste0("1_ClusTCR/","IgLK_",sub_directory,"_clusTCR_10x_",today(),".csv")
        # file_name_tcrex
        write.csv(IGLK, file_name_clust, row.names = F)
      }
      # } else {
      TRBD <-  (rbind(df[grep("TRBV", df$v_call), ], df[grep("TRDV", df$v_call), ]))
      if(nrow(TRBD)>0) {
        file_name_clust <- paste0("1_ClusTCR/","BD_",sub_directory,"_clusTCR_10x_",today(),".csv")
        # file_name_tcrex
        write.csv(TRBD, file_name_clust, row.names = F)
      }

    } else {
      message("if the files are present, change downloadClusTCR to T")
    }




    #### Seurat -----
    if (downloadSeurat) {

      ##### meta data ------
      contigs <- TCR
      names(contigs)[names(contigs) %in% "barcode"] <- "Cell_Index"
      contigs_lim <- contigs[!names(contigs) %in% c("is_cell", "contig_id", "high_confidence", "raw_consensus_id", "exact_subclonotype_id", "umis", "reads", "length", "germline_alignment", "cdr3_nt_id", "cdr3_nt_alignment", "rev_comp", names(contigs[grep("_end", names(contigs))]), names(contigs[grep("cigar", names(contigs))]), names(contigs[grep("length", names(contigs))]), names(contigs[grep("count", names(contigs))]), names(contigs[grep("sequence", names(contigs))]), names(contigs[grep("fwr", names(contigs))]), names(contigs[grep("_nt", names(contigs))]), "productive","raw_clonotype_id"
      )
      ]

      contigs_lim$chain <- ifelse(grepl("TRA", contigs_lim$v_gene), "TRA",
                                  ifelse(grepl("TRB", contigs_lim$v_gene), "TRB",
                                         ifelse(grepl("TRG", contigs_lim$v_gene), "TRG",
                                                ifelse(grepl("TRD", contigs_lim$v_gene), "TRD",
                                                       ""
                                                )
                                         )
                                  )
      )

      contig_AG <- subset(contigs_lim, contigs_lim$chain == "TRA" | contigs_lim$chain == "TRG")
      head(contig_AG)
      name.list <- names(contig_AG[c(
        names(contig_AG[grep("gene", names(contig_AG))]),
        names(contig_AG[grep("cdr1", names(contig_AG))]),
        names(contig_AG[grep("cdr2", names(contig_AG))]),
        names(contig_AG[grep("cdr3", names(contig_AG))]),
        "chain"
      )])
      contig_AG <- contig_AG %>%
        select(all_of(name.list), everything())
      names(contig_AG)[1:summary(name.list)[1]] <- paste(names(contig_AG[names(contig_AG) %in% name.list]), "_AG", sep = "")
      head(contig_AG)

      contig_BD <- subset(contigs_lim, contigs_lim$chain == "TRB" | contigs_lim$chain == "TRD")
      name.list <- names(contig_BD[c(
        names(contig_BD[grep("gene", names(contig_BD))]),
        names(contig_BD[grep("cdr1", names(contig_BD))]),
        names(contig_BD[grep("cdr2", names(contig_BD))]),
        names(contig_BD[grep("cdr3", names(contig_BD))]),
        "chain"
      )])
      contig_BD <- contig_BD %>%
        select(all_of(name.list), everything())

      names(contig_BD)[1:summary(name.list)[1]] <- paste(names(contig_BD[names(contig_BD) %in% name.list]), "_BD", sep = "")

      contig_paired <- merge(contig_AG, contig_BD, by = c("Cell_Index"), all = T)
      head(contig_paired)

      contig_paired$chain_BD == "TRB"

      contig_paired$pairing <- ifelse(contig_paired$chain_BD == "TRB" & contig_paired$chain_AG == "TRA", "abTCR Paired",
                                      ifelse(contig_paired$chain_BD == "TRD" & contig_paired$chain_AG == "TRG", "gdTCR Paired", "unpaired")
      )

      contig_paired$pairing[is.na(contig_paired$pairing)] <- "unpaired"
      table(contig_paired$pairing)

      contig_paired <- contig_paired[!names(contig_paired) %in% c("d_gene_AG")]
      contig_paired_only <- contig_paired
      contig_paired_only <- subset(contig_paired_only, contig_paired_only$cdr3_BD != "None")
      contig_paired_only <- subset(contig_paired_only, contig_paired_only$cdr3_AG != "None")
      contig_paired_only$d_gene_BD <- gsub("^$", "NA", contig_paired_only$d_gene_BD)
      #
      contig_paired_only$vj_gene_AG <- paste(contig_paired_only$v_gene_AG, contig_paired_only$j_gene_AG, sep = ".")
      contig_paired_only$vj_gene_AG <- gsub("NA.NA", "", contig_paired_only$vj_gene_AG)
      #
      contig_paired_only$vj_gene_BD <- paste(contig_paired_only$v_gene_BD, contig_paired_only$j_gene_BD, sep = ".")
      contig_paired_only$vj_gene_BD <- gsub(".NA.", ".", contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub(".None.", ".", contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub("NA.NA", "", contig_paired_only$vj_gene_BD)
      #
      contig_paired_only$vdj_gene_BD <- paste(contig_paired_only$v_gene_BD, contig_paired_only$d_gene_BD, contig_paired_only$j_gene_BD, sep = ".")
      contig_paired_only$vdj_gene_BD <- gsub(".NA.", ".", contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("[.]NA[.]", ".", contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("NA.NA", "", contig_paired_only$vdj_gene_BD)
      #
      contig_paired_only$vj_gene_cdr3_AG <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$cdr3_AG, sep = "_")
      contig_paired_only$vj_gene_cdr3_AG <- gsub("_NA", "", contig_paired_only$vj_gene_cdr3_AG)
      #
      contig_paired_only$vj_gene_cdr3_BD <- paste(contig_paired_only$vj_gene_BD, contig_paired_only$cdr3_BD, sep = "_")
      contig_paired_only$vj_gene_cdr3_BD <- gsub("_NA", "", contig_paired_only$vj_gene_cdr3_BD)
      #
      contig_paired_only$vdj_gene_cdr3_BD <- paste(contig_paired_only$vdj_gene_BD, contig_paired_only$cdr3_BD, sep = "_")
      contig_paired_only$vdj_gene_cdr3_BD <- gsub("_NA", "", contig_paired_only$vdj_gene_cdr3_BD)
      #
      contig_paired_only$vj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$vj_gene_BD, sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$vdj_gene_BD, sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- gsub("^ & ", "", contig_paired_only$vdj_gene_AG_BD)
      contig_paired_only$vdj_gene_AG_BD <- gsub(" & $", "", contig_paired_only$vdj_gene_AG_BD)
      #
      # #updating names to be consistant....
      contig_paired_only$vj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG, contig_paired_only$vj_gene_cdr3_BD, sep = " & ")
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub("^ & ", "", contig_paired_only$vj_gene_cdr3_AG_BD)
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub(" & $", "", contig_paired_only$vj_gene_cdr3_AG_BD)

      contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG, contig_paired_only$vdj_gene_cdr3_BD, sep = " & ")
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub("^ & ", "", contig_paired_only$vdj_gene_cdr3_AG_BD)
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub(" & $", "", contig_paired_only$vdj_gene_cdr3_AG_BD)
      dim(contig_paired_only)

      contig_paired_only_dup <- contig_paired_only[!duplicated(contig_paired_only$Cell_Index), ]
      # names(contig_paired_only)[names(contig_paired_only) %in% "cell_id"] <- "Cell_Index"

      contig_paired_only$Sample_Name <- sub_directory

      contig_paired_only <- contig_paired_only %>%
        select(all_of(c("Cell_Index", "Sample_Name")), everything())
      head(contig_paired_only)

      # paste(input$group10x, "_", input$Indiv10x, "", gsub("-", ".", Sys.Date()), ".csv", sep = "")
      file_name_md <- paste0("1_SeuratQC/",sub_directory,"_metadata_10x_",today(),".csv")
      write.csv(contig_paired_only, file_name_md, row.names = F)

      #### making the matrix -----
      mat <- as.data.frame(mat)
      rownames(mat) <- make.unique(features$V2) # select which feature to label genes...

      names(mat) <- barcode$V1

      mat$Gene_Name <- rownames(mat)
      mat <- mat %>%
        select(all_of("Gene_Name"), everything())

      file_name_mat <- paste0("1_SeuratQC/",sub_directory,"_count-matrix_10x_",today(),".csv.gz")
      write_csv(mat, gzfile(file_name_mat))

    } else {
      message("if the files are present, change downloadSeurat to T")
    }
    if (downloadTCR_Explore) {
      ##### meta data ------
      contigs <- TCR
      names(contigs)[names(contigs) %in% "barcode"] <- "Cell_Index"
      contigs_lim <- contigs[!names(contigs) %in% c("is_cell", "contig_id", "high_confidence", "raw_consensus_id", "exact_subclonotype_id", "umis", "reads", "length", "germline_alignment", "cdr3_nt_id", "cdr3_nt_alignment", "rev_comp", names(contigs[grep("_end", names(contigs))]), names(contigs[grep("cigar", names(contigs))]), names(contigs[grep("length", names(contigs))]), names(contigs[grep("count", names(contigs))]), names(contigs[grep("sequence", names(contigs))]), names(contigs[grep("fwr", names(contigs))]), names(contigs[grep("_nt", names(contigs))]), "productive","raw_clonotype_id"
      )
      ]

      contigs_lim$chain <- ifelse(grepl("TRA", contigs_lim$v_gene), "TRA",
                                  ifelse(grepl("TRB", contigs_lim$v_gene), "TRB",
                                         ifelse(grepl("TRG", contigs_lim$v_gene), "TRG",
                                                ifelse(grepl("TRD", contigs_lim$v_gene), "TRD",
                                                       ""
                                                )
                                         )
                                  )
      )

      contig_AG <- subset(contigs_lim, contigs_lim$chain == "TRA" | contigs_lim$chain == "TRG")
      head(contig_AG)
      name.list <- names(contig_AG[c(
        names(contig_AG[grep("gene", names(contig_AG))]),
        names(contig_AG[grep("cdr1", names(contig_AG))]),
        names(contig_AG[grep("cdr2", names(contig_AG))]),
        names(contig_AG[grep("cdr3", names(contig_AG))]),
        "chain"
      )])
      contig_AG <- contig_AG %>%
        select(all_of(name.list), everything())
      names(contig_AG)[1:summary(name.list)[1]] <- paste(names(contig_AG[names(contig_AG) %in% name.list]), "_AG", sep = "")
      head(contig_AG)

      contig_BD <- subset(contigs_lim, contigs_lim$chain == "TRB" | contigs_lim$chain == "TRD")
      name.list <- names(contig_BD[c(
        names(contig_BD[grep("gene", names(contig_BD))]),
        names(contig_BD[grep("cdr1", names(contig_BD))]),
        names(contig_BD[grep("cdr2", names(contig_BD))]),
        names(contig_BD[grep("cdr3", names(contig_BD))]),
        "chain"
      )])
      contig_BD <- contig_BD %>%
        select(all_of(name.list), everything())

      names(contig_BD)[1:summary(name.list)[1]] <- paste(names(contig_BD[names(contig_BD) %in% name.list]), "_BD", sep = "")

      contig_paired <- merge(contig_AG, contig_BD, by = c("Cell_Index"), all = T)

      contig_paired$chain_BD == "TRB"

      contig_paired$pairing <- ifelse(contig_paired$chain_BD == "TRB" & contig_paired$chain_AG == "TRA", "abTCR Paired",
                                      ifelse(contig_paired$chain_BD == "TRD" & contig_paired$chain_AG == "TRG", "gdTCR Paired", "unpaired")
      )

      contig_paired$pairing[is.na(contig_paired$pairing)] <- "unpaired"

      contig_paired <- contig_paired[!names(contig_paired) %in% c("d_gene_AG")]
      contig_paired_only <- contig_paired
      contig_paired_only <- subset(contig_paired_only, contig_paired_only$cdr3_BD != "None")
      contig_paired_only <- subset(contig_paired_only, contig_paired_only$cdr3_AG != "None")
      contig_paired_only$d_gene_BD <- gsub("^$", "NA", contig_paired_only$d_gene_BD)
      #
      contig_paired_only$vj_gene_AG <- paste(contig_paired_only$v_gene_AG, contig_paired_only$j_gene_AG, sep = ".")
      contig_paired_only$vj_gene_AG <- gsub("NA.NA", "", contig_paired_only$vj_gene_AG)
      #
      contig_paired_only$vj_gene_BD <- paste(contig_paired_only$v_gene_BD, contig_paired_only$j_gene_BD, sep = ".")
      contig_paired_only$vj_gene_BD <- gsub(".NA.", ".", contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub(".None.", ".", contig_paired_only$vj_gene_BD)
      contig_paired_only$vj_gene_BD <- gsub("NA.NA", "", contig_paired_only$vj_gene_BD)
      #
      contig_paired_only$vdj_gene_BD <- paste(contig_paired_only$v_gene_BD, contig_paired_only$d_gene_BD, contig_paired_only$j_gene_BD, sep = ".")
      contig_paired_only$vdj_gene_BD <- gsub(".NA.", ".", contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("[.]NA[.]", ".", contig_paired_only$vdj_gene_BD)
      contig_paired_only$vdj_gene_BD <- gsub("NA.NA", "", contig_paired_only$vdj_gene_BD)
      #
      contig_paired_only$vj_gene_cdr3_AG <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$cdr3_AG, sep = "_")
      contig_paired_only$vj_gene_cdr3_AG <- gsub("_NA", "", contig_paired_only$vj_gene_cdr3_AG)
      #
      contig_paired_only$vj_gene_cdr3_BD <- paste(contig_paired_only$vj_gene_BD, contig_paired_only$cdr3_BD, sep = "_")
      contig_paired_only$vj_gene_cdr3_BD <- gsub("_NA", "", contig_paired_only$vj_gene_cdr3_BD)
      #
      contig_paired_only$vdj_gene_cdr3_BD <- paste(contig_paired_only$vdj_gene_BD, contig_paired_only$cdr3_BD, sep = "_")
      contig_paired_only$vdj_gene_cdr3_BD <- gsub("_NA", "", contig_paired_only$vdj_gene_cdr3_BD)
      #
      contig_paired_only$vj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$vj_gene_BD, sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- paste(contig_paired_only$vj_gene_AG, contig_paired_only$vdj_gene_BD, sep = " & ")
      contig_paired_only$vdj_gene_AG_BD <- gsub("^ & ", "", contig_paired_only$vdj_gene_AG_BD)
      contig_paired_only$vdj_gene_AG_BD <- gsub(" & $", "", contig_paired_only$vdj_gene_AG_BD)
      #
      # #updating names to be consistant....
      contig_paired_only$vj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG, contig_paired_only$vj_gene_cdr3_BD, sep = " & ")
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub("^ & ", "", contig_paired_only$vj_gene_cdr3_AG_BD)
      contig_paired_only$vj_gene_cdr3_AG_BD <- gsub(" & $", "", contig_paired_only$vj_gene_cdr3_AG_BD)

      contig_paired_only$vdj_gene_cdr3_AG_BD <- paste(contig_paired_only$vj_gene_cdr3_AG, contig_paired_only$vdj_gene_cdr3_BD, sep = " & ")
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub("^ & ", "", contig_paired_only$vdj_gene_cdr3_AG_BD)
      contig_paired_only$vdj_gene_cdr3_AG_BD <- gsub(" & $", "", contig_paired_only$vdj_gene_cdr3_AG_BD)
      dim(contig_paired_only)

      contig_paired_only_dup <- contig_paired_only[!duplicated(contig_paired_only$Cell_Index), ]
      # names(contig_paired_only)[names(contig_paired_only) %in% "cell_id"] <- "Cell_Index"

      contig_paired_only$indiv.group <- sub_directory

      contig_paired_only <- contig_paired_only %>%
        select(all_of(c("Cell_Index", "indiv.group")), everything())
      head(contig_paired_only)

      contig_paired_only$cloneCount <- 1
      contig_paired_only$group <-  strsplit(sub_directory,split = "_")[[1]][2]
      contig_paired_only$indiv <-  strsplit(sub_directory,split = "_")[[1]][1]

      contig_paired_only <- contig_paired_only %>%
        select(cloneCount, indiv.group, group, indiv, everything())
      contig_paired_only

      file_name_tcrexplore <- paste0("1_TCR_Explore/",sub_directory,"_TCR_Explore_10x_",today(),".csv")
      write.csv(contig_paired_only, file_name_tcrexplore, row.names = F)


      message(paste0("Downloaded TCR_Explore for ",sub_directory))

    } else {
      message("if the files are present, change downloadTCR_Explore to T")
    }


  }
}









