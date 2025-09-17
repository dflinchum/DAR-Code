library(rtracklayer)
library(Sierra)
library(ChIPpeakAnno)
library(tidyverse)
library(ChIPseeker)
library(openxlsx)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene




setwd("D:/3_Cut&TAG_Files/hg38Aligned/4_macs2Output/1_peakFiles/BER")
BER_H3K4me1 <- import.bed("BER_H3K4me1_hg38.bed")
BER_H3K4me3<-import.bed("BER_H3K4me3_hg38.bed")
BER_H3K27ac<-import.bed("BER_H3K27ac_hg38.bed")
BER_H3K27me3<-import.bed("BER_H3K27me3_hg38.bed")
BER_enh <- import.bed("BER_H3K4me1_H3K27ac_hg38.bed")

setwd("D:/3_Cut&TAG_Files/hg38Aligned/4_macs2Output/1_peakFiles/JN")
JN_H3K4me1<-import.bed("JN_H3K4me1_hg38.bed")
JN_H3K4me3<-import.bed("JN_H3K4me3_hg38.bed")
JN_H3K27ac<-import.bed("JN_H3K27ac_hg38.bed")
JN_H3K27me3<-import.bed("JN_H3K27me3_hg38.bed")
JN_enh <- import.bed("JN_H3K4me1_H3K27ac_hg38.bed")

setwd("D:/3_Cut&TAG_Files/hg38Aligned/4_macs2Output/1_peakFiles/BOD")
BOD_H3K4me1<-import.bed("BOD_H3K4me1_hg38.bed")
BOD_H3K4me3<-import.bed("BOD_H3K4me3_hg38.bed")
BOD_H3K27ac<-import.bed("BOD_H3K27ac_hg38.bed")
BOD_H3K27me3<-import.bed("BOD_H3K27me3_hg38.bed")
BOD_enh <- import.bed("BOD_H3K4me1_H3K27ac_hg38.bed")


BER_H3K4me1_anno<- annotatePeak(BER_H3K4me1, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
BER_H3K4me3_anno <- annotatePeak(BER_H3K4me3, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
BER_H3K27ac_anno<- annotatePeak(BER_H3K27ac, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
BER_H3K27me3_anno<- annotatePeak(BER_H3K27me3, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
BER_enh_anno<- annotatePeak(BER_enh, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 

JN_H3K4me1_anno<- annotatePeak(JN_H3K4me1, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
JN_H3K4me3_anno <- annotatePeak(JN_H3K4me3, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
JN_H3K27ac_anno<- annotatePeak(JN_H3K27ac, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
JN_H3K27me3_anno<- annotatePeak(JN_H3K27me3, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
JN_enh_anno<- annotatePeak(JN_enh, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 

BOD_H3K4me1_anno<- annotatePeak(BOD_H3K4me1, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
BOD_H3K4me3_anno <- annotatePeak(BOD_H3K4me3, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
BOD_H3K27ac_anno<- annotatePeak(BOD_H3K27ac, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
BOD_H3K27me3_anno<- annotatePeak(BOD_H3K27me3, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
BOD_enh_anno<- annotatePeak(BOD_enh, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 


BER_H3K4me1_anno <- as.data.frame(BER_H3K4me1_anno)
BER_H3K4me3_anno <- as.data.frame(BER_H3K4me3_anno)
BER_H3K27ac_anno <- as.data.frame(BER_H3K27ac_anno)
BER_H3K27me3_anno <- as.data.frame(BER_H3K27me3_anno)
BER_enh_anno <- as.data.frame(BER_enh_anno)

JN_H3K4me1_anno <- as.data.frame(JN_H3K4me1_anno)
JN_H3K4me3_anno <- as.data.frame(JN_H3K4me3_anno)
JN_H3K27ac_anno <- as.data.frame(JN_H3K27ac_anno)
JN_H3K27me3_anno <- as.data.frame(JN_H3K27me3_anno)
JN_enh_anno <- as.data.frame(JN_enh_anno)

BOD_H3K4me1_anno <- as.data.frame(BOD_H3K4me1_anno)
BOD_H3K4me3_anno <- as.data.frame(BOD_H3K4me3_anno)
BOD_H3K27ac_anno <- as.data.frame(BOD_H3K27ac_anno)
BOD_H3K27me3_anno <- as.data.frame(BOD_H3K27me3_anno)
BOD_enh_anno <- as.data.frame(BOD_enh_anno)




##ber tss playing around
BER_H3K4me3_filter <- BER_H3K4me3_anno %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))

BER_H3K4me3_counts <- BER_H3K4me3_filter %>%  ##WORKED ITS PERFECT; 
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

BER_H3K4me3_TSS<- BER_H3K4me3_filter %>% filter(distanceToTSS > 0 & distanceToTSS < 2500) ##have to use "&" or "AND", not the "|"; this | is an "OR" function!!! 
BER_H3K4me3_TSS_small <- BER_H3K4me3_TSS %>% filter(distanceToTSS < 201)

BER_H3K4me3_TSS_small_distinct<- BER_H3K4me3_TSS_small %>% distinct(SYMBOL, .keep_all = TRUE) #removes any genes that have multiples


#jn tss playing around 
JN_H3K4me3_filter <- JN_H3K4me3_anno %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))

JN_H3K4me3_counts <- JN_H3K4me3_filter %>%  ##WORKED ITS PERFECT; 
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

JN_H3K4me3_TSS <- JN_H3K4me3_filter %>% filter(distanceToTSS > 0 & distanceToTSS < 2500)
JN_H3K4me3_TSS_small <- JN_H3K4me3_TSS %>% filter(distanceToTSS < 201)
JN_H3K4me3_TSS_small_distinct<- JN_H3K4me3_TSS_small %>% distinct(SYMBOL, .keep_all = TRUE) #removes any genes that have multiples 

#bod tss playing around
BOD_H3K4me3_filter <- BOD_H3K4me3_anno %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))

BOD_H3K4me3_counts <- BOD_H3K4me3_filter %>%  ##WORKED ITS PERFECT; 
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

BOD_H3K4me3_TSS <- BOD_H3K4me3_filter %>% filter(distanceToTSS > 0 & distanceToTSS < 2500) 
BOD_H3K4me3_TSS_small <- BOD_H3K4me3_TSS %>% filter(distanceToTSS < 201)
BOD_H3K4me3_TSS_small_distinct<- BOD_H3K4me3_TSS_small %>% distinct(SYMBOL, .keep_all = TRUE) #removes any genes that have multiples 







##extracting gene via column and then intersecting; works fine but maybe not best
BERtest <- BER_H3K4me3_TSS_small_distinct %>% dplyr::select(SYMBOL)#such a crude way to do this lol
JNtest <- JN_H3K4me3_TSS_small_distinct %>% dplyr::select(SYMBOL)
test2 <- intersect(BERtest, JNtest)
BODtest <- BOD_H3K4me3_TSS_small_distinct %>% dplyr::select(SYMBOL)
test3 <- intersect(BODtest, test2)


##try merging and then stripping everything away we don't care about
BERmerge <- merge(BER_H3K4me3_TSS_small_distinct, JN_H3K4me3_TSS_small_distinct, by="SYMBOL")
JNBER<- BERmerge %>% dplyr::select(SYMBOL, GENENAME.x, ENSEMBL.x, distanceToTSS.x)
write.xlsx(JNBER, "JNBER_H3K4me3_Promoters.xlsx")

JNBERBOD <- merge(BERmerge, BOD_H3K4me3_TSS_small_distinct, by="SYMBOL")
JNBERBODfinal <- JNBERBOD %>% dplyr::select(SYMBOL, GENENAME.x, ENSEMBL.x, distanceToTSS.x)
write.xlsx(JNBERBODfinal, "JNBERBOD_H3K4me3_Promoters.xlsx")
