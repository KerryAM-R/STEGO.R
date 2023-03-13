barcode = read.table("~/Documents/Post-doc/PostDoc_Antwerp/R_files/1_Files_scTCRseq/CZI/RawBcellfiles/GSE171703_RAW/GSM5231088_R125_barcodes.tsv.gz")
feature = read.table("~/Documents/Post-doc/PostDoc_Antwerp/R_files/1_Files_scTCRseq/CZI/RawBcellfiles/GSE171703_RAW/GSM5231088_R125_feature_barcode.csv.gz")
matrix = Matrix::readMM("~/Documents/Post-doc/PostDoc_Antwerp/R_files/1_Files_scTCRseq/CZI/RawBcellfiles/GSE171703_RAW/GSM5231088_R125_matrix.mtx.gz")
contigs = read.csv("~/Documents/Post-doc/PostDoc_Antwerp/R_files/1_Files_scTCRseq/CZI/RawBcellfiles/GSE171703_RAW/GSM5231088_R125_filtered_contig_annotations.csv.gz")

processed <- read.csv("inst/extdata/10x/B_cells/R125_BCR_count-matrix_10x_2023.03.13.csv.gz")

dim(processed)[2]

dim(contigs)
contigs <- contigs[order(contigs$umis,decreasing = T),]

contigs_lim <- contigs[!names(contigs) %in% c("is_cell","contig_id","high_confidence","raw_consensus_id","exact_subclonotype_id","reads","length","cdr3_nt",names(contigs[grep("fwr",names(contigs))]),names(contigs[grep("cdr1",names(contigs))]),names(contigs[grep("cdr2",names(contigs))])
)]
contigs_lim
contig_LK <- subset(contigs_lim,contigs_lim$chain=="IGL" | contigs_lim$chain=="IGK")
name.list <- names(contig_LK[c(names(contig_LK[grep("gene",names(contig_LK))]),
                               names(contig_LK[grep("cdr3",names(contig_LK))]),
                               "chain")])
name.list
contig_LK <- contig_LK %>%
  select(all_of(name.list), everything())
names(contig_LK)[1:summary(name.list)[1]] <-paste(names(contig_LK[names(contig_LK) %in% name.list]),"_LK",sep="")
contig_LK

contig_H <- subset(contigs_lim,contigs_lim$chain=="IGH")

name.list <- names(contig_H[c(names(contig_H[grep("gene",names(contig_H))]),
                               names(contig_H[grep("cdr3",names(contig_H))]),
                               "chain")])
contig_H <- contig_H %>%
  select(all_of(name.list), everything())


names(contig_H)[1:summary(name.list)[1]] <-paste(names(contig_H[names(contig_H) %in% name.list]),"_H",sep="")
contig_H
# contig_paired <- merge(contig_LK,contig_H, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
# contig_paired <- merge(contig_LK,contig_H, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)
contig_paired <- merge(contig_LK,contig_H, by=c("barcode", "full_length" ,"productive" ,"raw_clonotype_id"),all = T)

contig_paired$pairing <- ifelse(contig_paired$chain_H=="IGH" & contig_paired$chain_LK=="IGK","IGK Paired",
                                ifelse(contig_paired$chain_H=="IGH" & contig_paired$chain_LK=="IGL","IGL Paired",NA
                                ))

contig_paired
contig_paired$pairing[is.na(contig_paired$pairing)] <- "unpaired"
contig_paired <- contig_paired[!names(contig_paired) %in% c("d_gene_LK")]

dim(contig_paired)

contig_paired_only <- contig_paired
contig_paired_only <- subset(contig_paired_only,contig_paired_only$cdr3_H!="None")

contig_paired_only <- subset(contig_paired_only,contig_paired_only$cdr3_LK!="None")
dim(contig_paired_only)

contig_paired_only$d_gene_H <- sub("^$","NA", contig_paired_only$d_gene_H)
#
contig_paired_only$vj_gene_LK <- paste(contig_paired_only$v_gene_LK,contig_paired_only$j_gene_LK,sep = ".")
contig_paired_only$vj_gene_LK <- gsub("NA.NA","",contig_paired_only$vj_gene_LK)
#
contig_paired_only$vj_gene_H <- paste(contig_paired_only$v_gene_H,contig_paired_only$j_gene_H,sep = ".")
contig_paired_only$vj_gene_H <- gsub(".NA.",".",contig_paired_only$vj_gene_H)
contig_paired_only$vj_gene_H <- gsub(".None.",".",contig_paired_only$vj_gene_H)
contig_paired_only$vj_gene_H <- gsub("NA.NA","",contig_paired_only$vj_gene_H)
#
contig_paired_only$vdj_gene_H <- paste(contig_paired_only$v_gene_H,contig_paired_only$d_gene_H,contig_paired_only$j_gene_H,sep = ".")
contig_paired_only$vdj_gene_H <- gsub(".NA.",".",contig_paired_only$vdj_gene_H)
contig_paired_only$vdj_gene_H <- gsub(".None.",".",contig_paired_only$vdj_gene_H)
contig_paired_only$vdj_gene_H <- gsub("NA.NA","",contig_paired_only$vdj_gene_H)
#
contig_paired_only$vj_gene_cdr3_LK <- paste(contig_paired_only$vj_gene_LK,contig_paired_only$cdr3_LK,sep = "_")
contig_paired_only$vj_gene_cdr3_LK <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_LK)
#
contig_paired_only$vj_gene_cdr3_H <- paste(contig_paired_only$vj_gene_H,contig_paired_only$cdr3_H,sep = "_")
contig_paired_only$vj_gene_cdr3_H <- gsub("_NA","",contig_paired_only$vj_gene_cdr3_H)
#
contig_paired_only$vdj_gene_cdr3_H <- paste(contig_paired_only$vdj_gene_H,contig_paired_only$cdr3_H,sep = "_")
contig_paired_only$vdj_gene_cdr3_H <- gsub("_NA","",contig_paired_only$vdj_gene_cdr3_H)
#
contig_paired_only$vj_gene_LK_H <- paste(contig_paired_only$vj_gene_LK,contig_paired_only$vj_gene_H,sep = " & ")
contig_paired_only$vdj_gene_LK_H <- paste(contig_paired_only$vj_gene_LK,contig_paired_only$vdj_gene_H,sep = " & ")
contig_paired_only$vdj_gene_LK_H <- gsub("^ & ","",contig_paired_only$vdj_gene_LK_H)
contig_paired_only$vdj_gene_LK_H <- gsub(" & $","",contig_paired_only$vdj_gene_LK_H)
#
# #updating names to be consistant....
contig_paired_only$vj_gene_cdr3_LK_H <- paste(contig_paired_only$vj_gene_cdr3_LK,contig_paired_only$vj_gene_cdr3_H,sep = " & ")
contig_paired_only$vj_gene_cdr3_LK_H <- gsub("^ & ","",contig_paired_only$vj_gene_cdr3_LK_H)
contig_paired_only$vj_gene_cdr3_LK_H <- gsub(" & $","",contig_paired_only$vj_gene_cdr3_LK_H)

contig_paired_only$vdj_gene_cdr3_LK_H <- paste(contig_paired_only$vj_gene_cdr3_LK,contig_paired_only$vdj_gene_cdr3_H,sep = " & ")
contig_paired_only$vdj_gene_cdr3_LK_H <- gsub("^ & ","",contig_paired_only$vdj_gene_cdr3_LK_H)
contig_paired_only$vdj_gene_cdr3_LK_H <- gsub(" & $","",contig_paired_only$vdj_gene_cdr3_LK_H)
# contig_paired_only$vdj_gene_cdr3_LK_H <- paste(contig_paired_only$vj_gene_cdr3_LK,contig_paired_only$vdj_gene_cdr3_H,sep = " & ")
names(contig_paired_only)[names(contig_paired_only) %in% "barcode"] <- "Cell_Index"
write.csv(contig_paired_only,"contig_paired_only.csv")
dup <- contig_paired_only[duplicated(contig_paired_only$Cell_Index),]
contig_paired_only <- contig_paired_only[order(contig_paired_only$Cell_Index, contig_paired_only$umis.x,contig_paired_only$umis.y,decreasing = T),]

contig_paired_only_dup <- contig_paired_only[!duplicated(contig_paired_only$Cell_Index),] # remove duplicate barcodes.
names(contig_paired_only_dup)
contig_paired_only_dup <- contig_paired_only_dup[!names(contig_paired_only_dup) %in% c("umis.x","umis.y")]

contig_paired_only_dup$Sample_Name <- "Sample_name"

contig_paired_only_dup <- contig_paired_only_dup %>%
  select(all_of(c("Cell_Index","Sample_Name")), everything())

merge.names <- names(contig_paired_only_dup)[!grepl("_H",names(contig_paired_only_dup)) & !grepl("_LK",names(contig_paired_only_dup))]


contig_paired_only_dup[names(contig_paired_only_dup) %in% c(names(contig_paired_only_dup)[grepl("chain",names(contig_paired_only_dup))],"pairing")]

