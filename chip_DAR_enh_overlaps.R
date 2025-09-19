











WT1_JNBERBOD_enh <- GenomicRanges::intersect(BERJNBOD_H3K4me1_H3K27ac, WT1_GR)
WT1_JNBERBOD_enh_anno <- annotatePeak(WT1_JNBERBOD_enh, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db")
WT1_JNBERBOD_enh_df <- as.data.frame(WT1_JNBERBOD_enh_anno)



#DAR
setwd("D:/1_ATAC Files/ATAC 2022 Files")
JN_DEpeaksAnno <- annotatePeak("JN_NP_DEpeaks.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0


JN_DEpeaksAnno<-as.data.frame(JN_DEpeaksAnno) #as data frame to manipulate easier
JN_DEpeaksAnno <- JN_DEpeaksAnno %>%  #rename columns! 
  rename("FC" = "V5") #I SWEAR TO GOD THIS IS DIFFERENT EVERY TIME ON HOW TO CHANGE THE NAME

JN_DEpeaks_DAR <- JN_DEpeaksAnno %>% filter(FC > 1 | FC < -1)  #filter by abs value; works great, nothing between -.99 to .99!  



JN_DEpeaks_DAR_GR <- JN_DEpeaks_DAR[,c(1,2,3,5)]
JN_DEpeaks_DAR_GR <- makeGRangesFromDataFrame(JN_DEpeaks_DAR_GR)


#WT1-DARs
WT1_DARs <- GenomicRanges::intersect(WT1_GR, JN_DEpeaks_DAR_GR)
WT1_DARs_df <- as.data.frame(WT1_DARs)

#WT1-JN_enhancers
WT1_JNH3K4me1_H3K27ac <- GenomicRanges::intersect(WT1_GR, JN_H3K4me1_H3K27ac)
WT1_JNH3K4me1_H3K27ac_df <- as.data.frame(WT1_JNH3K4me1_H3K27ac)

#DARs-enhancers
JN_DAR_H3K4me1_H3K27ac <- GenomicRanges::intersect(JN_DEpeaks_DAR_GR, JN_H3K4me1_H3K27ac)
JN_DAR_H3K4me1_H3K27ac_df <- as.data.frame(JN_DAR_H3K4me1_H3K27ac)

#WT1-DARs-enhancers
JN_WT1_DARs_enhancers <- GenomicRanges::intersect(WT1_DARs, JN_H3K4me1_H3K27ac)
JN_WT1_DARs_enhancers_df <- as.data.frame(JN_WT1_DARs_enhancers)





png("JN_ChIP_DARs.png",
    width = 8, height = 8, # Dimensions in inches
    units = "in",
    res = 300)

draw.pairwise.venn(
  area1 = 3265, # Total in binding
  area2 = 3478, # Total in DARs
  cross.area = 779, # overlap
  category = c("EWSR1::WT1", "JN-DSRCT DARs"),
  fill = c("blue", "orange"),
  euler.d = FALSE, # Disable euler diagrams for accurate representation
  scaled = FALSE, # Disable scaling for accurate representation
  print.mode = "raw", # Print raw numbers instead of percentages
  cat.pos = c(-50, 50), # Adjust category label positions
  cat.dist = c(0.17, 0.17), # Adjust category label distances
  margin = 0.1,
  cex = 1.5,      # Increase the size of the numbers inside the diagram
  cat.cex = 1.8 # Adjust margin around the diagram
)

dev.off()






png("JN_ChIP_DARs_enhancers2.png", 
    width = 8, height = 8, # Dimensions in inches
    units = "in", 
    res = 300)
draw.triple.venn(
  area1 = 3265, # Total in ChIP
  area2 = 3478, # Total in DARs
  area3 = 57842, # Total in enhancers
  n12 = 779,  # ChIP and DARs overlap
  n23 = 1968,  # DARs and enhancers overlap
  n13 = 2272,  # ChIP and enhancers overlap
  n123 = 581, # All three overlap
  category = c("ChIP", "DARs", "Enhancers"),
  fill = c("blue", "green", "red"),
  euler.d = TRUE, # Disable euler diagrams for accurate representation
  scaled = FALSE, # Disable scaling for accurate representation
  print.mode = "raw", # Print raw numbers instead of percentages
  cat.pos = c(330, 30, 180), # Adjust category label positions
  cat.dist = c(0.07, 0.07, 0.05), # Adjust category label distances
  margin = 0.1,
  cex = 1.5,     # Increase the size of the numbers inside the diagram
  cat.cex = 1.8 # Adjust margin around the diagram
)

dev.off() 



##4.13
#venns for chip and promoters/enhancers 


png("JN_ChIP_H3K4me3.png",
    width = 8, height = 8, # Dimensions in inches
    units = "in",
    res = 300)

draw.pairwise.venn(
  area1 = 3265, # Total in binding
  area2 = 20211, # Total in promoters
  cross.area = 1069, # overlap
  category = c("EWSR1::WT1", "Promoters"),
  fill = c("blue", "darkgreen"),
  euler.d = FALSE, # Disable euler diagrams for accurate representation
  scaled = FALSE, # Disable scaling for accurate representation
  print.mode = "raw", # Print raw numbers instead of percentages
  cat.pos = c(-50, 50), # Adjust category label positions
  cat.dist = c(0.17, 0.17), # Adjust category label distances
  margin = 0.1,
  cex = 1.5,      # Increase the size of the numbers inside the diagram
  cat.cex = 1.8 # Adjust margin around the diagram
)

dev.off()




png("JN_ChIP_H3K4me1_H3K27ac.png",
    width = 8, height = 8, # Dimensions in inches
    units = "in",
    res = 300)

draw.pairwise.venn(
  area1 = 3265, # Total in binding
  area2 = 57843, # Total in enhancers
  cross.area = 2272, # overlap
  category = c("EWSR1::WT1", "Promoters"),
  fill = c("blue", "lightgreen"),
  euler.d = FALSE, # Disable euler diagrams for accurate representation
  scaled = FALSE, # Disable scaling for accurate representation
  print.mode = "raw", # Print raw numbers instead of percentages
  cat.pos = c(-50, 50), # Adjust category label positions
  cat.dist = c(0.17, 0.17), # Adjust category label distances
  margin = 0.1,
  cex = 1.5,      # Increase the size of the numbers inside the diagram
  cat.cex = 1.8 # Adjust margin around the diagram
)

dev.off()

