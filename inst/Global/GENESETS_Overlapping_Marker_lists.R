df <- readxl::read_xlsx("~/Documents/Post-doc/PostDoc_Antwerp/STEGO/sup. STEGO Summary of prioritisation process.xlsx",sheet = 15)

STEGO.R::Load_required_packages()
names(df)

GSE114724 <- as.data.frame(df$GSE114724)
names(GSE114724) <- "ID"
GSE114724$GSE114724 <- 1
GSE114724 <- GSE114724[complete.cases(GSE114724),]
rownames(GSE114724) <- GSE114724$ID

GSE139555 <- as.data.frame(df$GSE139555)
names(GSE139555) <- "ID"
GSE139555$GSE139555 <- 1
GSE139555 <- GSE139555[complete.cases(GSE139555),]
  rownames(GSE139555) <- GSE139555$ID

GSE144469 <- as.data.frame(df$GSE144469)
names(GSE144469) <- "ID"
GSE144469$GSE144469 <- 1
GSE144469 <- GSE144469[complete.cases(GSE144469),]
rownames(GSE144469) <- GSE144469$ID

GSE145370 <- as.data.frame(df$GSE145370)
names(GSE145370) <- "ID"
GSE145370$GSE145370 <- 1
GSE145370 <- GSE145370[complete.cases(GSE145370),]
rownames(GSE145370) <- GSE145370$ID

GSE148190 <- as.data.frame(df$GSE148190)
names(GSE148190) <- "ID"
GSE148190$GSE148190 <- 1
GSE148190 <- GSE148190[complete.cases(GSE148190),]
rownames(GSE148190) <- GSE148190$ID

GSE161192 <- as.data.frame(df$GSE161192)
names(GSE161192) <- "ID"
GSE161192$GSE161192 <- 1
GSE161192 <- GSE161192[complete.cases(GSE161192),]
rownames(GSE161192) <- GSE161192$ID

  GSE165499 <- as.data.frame(df$GSE165499)

names(GSE165499) <- "ID"
GSE165499$GSE165499 <- 1
GSE165499 <- GSE165499[complete.cases(GSE165499),]
rownames(GSE165499) <- GSE165499$ID

GSE168859 <- as.data.frame(df$GSE168859)
names(GSE168859) <- "ID"
GSE168859$GSE168859 <- 1
GSE168859 <- GSE168859[complete.cases(GSE168859),]
rownames(GSE168859) <- GSE168859$ID

GSE176201 <- as.data.frame(df$GSE176201)
names(GSE176201) <- "ID"
GSE176201$GSE176201 <- 1
GSE176201 <- GSE176201[complete.cases(GSE176201),]
rownames(GSE176201) <- GSE176201$ID

GSE180268 <- as.data.frame(df$GSE180268)
names(GSE180268) <- "ID"
GSE180268$GSE180268 <- 1
GSE180268 <- GSE180268[complete.cases(GSE180268),]
rownames(GSE180268) <- GSE180268$ID

GSE184703 <- as.data.frame(df$GSE184703)
names(GSE184703) <- "ID"
GSE184703$GSE184703 <- 1
GSE184703 <- GSE184703[complete.cases(GSE184703),]
rownames(GSE184703) <- GSE184703$ID

GSE185659 <- as.data.frame(df$GSE185659)
names(GSE185659) <- "ID"
GSE185659$GSE185659 <- 1
GSE185659 <- GSE185659[complete.cases(GSE185659),]
rownames(GSE185659) <- GSE185659$ID

# [1] "GSE114724" "GSE139555" "GSE144469" "GSE145370" "GSE148190" "GSE161192" "GSE165499" "GSE168859" "GSE176201" "GSE180268" "GSE184703" "GSE185659"

all <- merge(GSE114724,GSE139555,by = "ID", all = T)
all <- merge(all,GSE144469,by = "ID", all = T)
all <- merge(all,GSE145370,by = "ID", all = T)
all <- merge(all,GSE148190,by = "ID", all = T)
all <- merge(all,GSE161192,by = "ID", all = T)
all <- merge(all,GSE165499,by = "ID", all = T)
all <- merge(all,GSE168859,by = "ID", all = T)
all <- merge(all,GSE176201,by = "ID", all = T)
all <- merge(all,GSE180268,by = "ID", all = T)
all <- merge(all,GSE184703,by = "ID", all = T)
all <- merge(all,GSE185659,by = "ID", all = T)

rownames(all) <- gsub(38961,"SEPTIN6",all$ID)
all[is.na(all)] <- 0
all <- all[,!names(all) %in% "ID"]


mat <- as.data.frame(all)
df.x <- make_comb_mat(mat)
list.names <- as.character(rownames(all))

ht = draw(UpSet(df.x,
                pt_size = unit(5, "mm"),
                lwd = 1,
                # row_names_gp =  gpar(fontfamily = input$font_type, fontsize = 12),
                # column_names_gp = gpar(fontfamily = input$font_type),
                top_annotation = upset_top_annotation(df.x,
                                                      add_numbers = T,
                                                      # annotation_name_gp = gpar(fontfamily = input$font_type),
                                                      gp = gpar(fill = "black"),
                ),
                right_annotation = upset_right_annotation(df.x,
                                                          add_numbers = T,
                                                          # numbers_gp = gpar(fontfamily = input$font_type,fontsize =12),
                                                          # annotation_name_gp = gpar(fontfamily = input$font_type,fontsize=12),
                                                          gp = gpar(fill = "black"),
                ),
                # set_order  =  list.names

), padding = unit(c(20, 20, 20, 20), "mm"))
ht
3/6
7/11

4/8
all$sum <- rowSums(all)
list_row <- (rownames(subset(all,all$sum>=5)))
length(list_row)
as.character(paste(unlist(list_row), collapse=';'))

