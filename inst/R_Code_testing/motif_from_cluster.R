Network_df <- read.csv(system.file("extdata","BDrhap/clusTCR/ClusTCR2.csv",package = "STEGO.R"))

require(ClusTCR2)

motifplot_2 <- function(ClusTCR, Clust_selected=NULL) {

  motif_DF <- ClusTCR
  z1 <- as.data.frame(t(as.data.frame(strsplit(motif_DF$CDR3_Vgene,"_"))))
  motif_DF$motif <- z1$V1

  motif_DF$len1 <- nchar(motif_DF[,names(motif_DF) %in% "motif"])

  if (length(which(motif_DF$Clust_size_order==c(Clust_selected)))>0) {


    motif_DF_interest <- motif_DF[motif_DF$Clust_size_order %in% c(Clust_selected),]

    # subset(motif_DF,motif_DF$cluster==Clust_selected)
    motif <- as.data.frame(t(as.data.frame(strsplit(motif_DF_interest$motif, ""))))
    len <- nchar(unique(motif_DF_interest$motif))[1]

    motif_count <- aa.count.function(cbind(x=1,y=2,motif),len)

    motif_count1_aa<-pcm2pfm(motif_count)


    ggseqlogo(motif_count, seq_type='aa', method='p') +
      ylab('bits')+
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      theme(
        axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
        axis.text.y = element_text(colour="black",size=20,angle=0,hjust=1,vjust=0,face="plain",family="serif"),
        axis.title.x=element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain",family="serif"),
        axis.title.y = element_text(colour="black",size=20,angle=90,hjust=.5,vjust=.5,face="plain",family="serif"),
        legend.title  =element_blank(),
        legend.position = "right",
        legend.text = element_text(colour="black", size=12,family="serif")
      )

  }

  else {
    paste0("Cluster, ", Clust_selected,", does not exist")

  }
}

Motif_from_cluster_file(Network_df,Clust_selected = c(3))

motifplot_2(Network_df,Clust_selected = c(3))
