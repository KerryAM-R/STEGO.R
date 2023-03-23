
temp = list.files(path = "~/Documents/GitHub/E-MTAB-9357/", pattern="*.txt.gz",full.names = T)
temp.name <-list.files(path = "~/Documents/GitHub/E-MTAB-9357/", pattern="*.txt.gz",full.names = F)
temp.name

myfiles = lapply(temp, read.delim)

for (i in 1:length(temp.name)) {
  myfiles[[i]]$Derived.Array.Data.File.3 <- temp.name[i]
}

View(myfiles[[186]])

df <- rbind(myfiles[[1]],myfiles[[2]])

for (i in 3:length(temp.name)) {
  df <- rbind(df,myfiles[[i]])
}

dim(df)
df_filtering <- subset(df,df$clonotype!="None")

dim(df_filtering)

df_filtering_beta <- subset(df_filtering,df_filtering$TRB_1_cdr3!="None")
dim(df_filtering_beta)
ID <- read.table("~/Documents/GitHub/E-MTAB-9357/ID.txt",sep= "\t",header = T)

ID_CD8 <- ID[,names(ID) %in% c("Source.Name","Characteristics.disease.","Derived.Array.Data.File.3")]
dim(ID_CD8)

ID_df_filtering_beta_CD8 <- merge(df_filtering_beta,ID_CD8,by="Derived.Array.Data.File.3")


ID_df_filtering_beta_CD8_2 <- ID_df_filtering_beta_CD8[,names(ID_df_filtering_beta_CD8) %in% c( "TRB_1_v_gene","TRB_1_cdr3","TRB_1_j_gene")]
head(ID_df_filtering_beta_CD8_2)
ID_df_filtering_beta_CD8_2 <- ID_df_filtering_beta_CD8_2 %>%
  select(TRB_1_v_gene,TRB_1_cdr3,TRB_1_j_gene, everything())
ID_df_filtering_beta_CD8_2$clone <- 1
names(ID_df_filtering_beta_CD8_2) <- c("TRBV_gene","CDR3_beta","TRBJ_gene","cloneCount")

calls_TCR_paired.fun3 <- ddply(ID_df_filtering_beta_CD8_2,names(ID_df_filtering_beta_CD8_2)[-c(4)] ,numcolwise(sum))
dim(calls_TCR_paired.fun3)
# write.table(calls_TCR_paired.fun3,"CD8_total.tsv", row.names = F,sep="\t", quote = F)


# ID_CD4 <- ID[,names(ID) %in% c("Source.Name","Characteristics.disease.","Derived.Array.Data.File.2")]
#
# ID_df_filtering_beta_CD4 <- merge(df_filtering_beta,ID_CD4,by.x="Derived.Array.Data.File.3",by.y="Derived.Array.Data.File.2")
#
# head(ID_df_filtering_beta_CD4)
#
# ID_df_filtering_beta_CD4_2 <- ID_df_filtering_beta_CD4[,names(ID_df_filtering_beta_CD4) %in% c("TRB_1_v_gene","TRB_1_cdr3","TRB_1_j_gene")]
# ID_df_filtering_beta_CD4_2 <- ID_df_filtering_beta_CD4_2 %>%
#   select(TRB_1_v_gene,TRB_1_cdr3, everything())
# ID_df_filtering_beta_CD4_2$clone <- 1
# names(ID_df_filtering_beta_CD4_2) <- c("TRBV_gene","CDR3_beta","TRBJ_gene","cloneCount")
#
# cm.colors(12)


ID_df_filtering_beta_CD8_2 <- ID_df_filtering_beta_CD8[,names(ID_df_filtering_beta_CD8) %in% c( "TRB_1_v_gene","TRB_1_cdr3","TRB_1_j_gene","Source.Name",
                                                                                                "Characteristics.disease.")]
head(ID_df_filtering_beta_CD8_2)
ID_df_filtering_beta_CD8_2 <- ID_df_filtering_beta_CD8_2 %>%
  select(TRB_1_v_gene,TRB_1_cdr3,TRB_1_j_gene, everything())
ID_df_filtering_beta_CD8_2$clone <- 1
names(ID_df_filtering_beta_CD8_2) <- c("TRBV_gene","CDR3_beta","TRBJ_gene","ID","Disease","cloneCount")

calls_TCR_paired.fun3 <- ddply(ID_df_filtering_beta_CD8_2,names(ID_df_filtering_beta_CD8_2)[-c(6)] ,numcolwise(sum))
dim(calls_TCR_paired.fun3)
# write.table(calls_TCR_paired.fun3,"CD8_Indiviudal.tsv", row.names = F,sep="\t", quote = F)

tcrex <- read.table("~/Downloads/tcrex_k847bxviei.tsv", header = T)
tcrex$TRBV_gene <- gsub("-0","-",tcrex$TRBV_gene)
tcrex$TRBV_gene <- gsub("V0","V",tcrex$TRBV_gene)

tcrex$TRBJ_gene <- gsub("-0","-",tcrex$TRBJ_gene)
tcrex$TRBJ_gene <- gsub("J0","J",tcrex$TRBJ_gene)
dim(tcrex)

tcrex_CD8 <- merge(calls_TCR_paired.fun3,tcrex,by=c("TRBV_gene","CDR3_beta","TRBJ_gene"),all.y=T)
dim(tcrex_CD8)
# write.csv(tcrex_CD8,"Individual.tcrex.cd8.csv")
