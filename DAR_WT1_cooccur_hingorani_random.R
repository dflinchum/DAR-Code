library(readxl)
library(rtracklayer)
library(tidygraph)
library(ChIPseeker)
library(clusterProfiler)
library(readxl)
library(GenomicRanges)
library(ChIPpeakAnno)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


#bershwt1
berminus_co <- read_xlsx("BERshWT1minusDox_coOccurMotifs.xlsx", sheet = 2)
berplus_co <- read_xlsx("BERshWT1plusDox_coOccurMotifs.xlsx", sheet = 2)
#JNshwt1
jnminus_co <- read_xlsx("JNshWT1minusDox_coOccurMotifs.xlsx", sheet = 2)
jnplus_co <- read_xlsx("JNshWT1plusDox_coOccurMotifs.xlsx", sheet = 2)

setwd("D:/1_ATAC Files/ATAC 2020 Files/9 Co_occurMotifs")
##reading in WT1 coOccur motifs
##these are WT1 motif locations/ranges in 2020 JNBERBOD ATAC consensus peaks
WT1<-read_excel("2020ATAC_WT1_coOccur_motif_peaks.xlsx")

WT1<-data.frame(WT1)
WT1<-makeGRangesFromDataFrame(WT1)
WT1_motifs_df <- as.data.frame(WT1)
WT1_anno <- annotatePeak(WT1, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db")

###### want to cross WT1 with 



#JN
#need to annotate for the gene 
setwd("D:/5_R_generated_Data/JN/3. 2022 ATAC")
JN_DEpeaksAnno <- annotatePeak("JN_NP_DEpeaks.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0


JN_DEpeaksAnno<-as.data.frame(JN_DEpeaksAnno) #as data frame to manipulate easier; should tibble? 
JN_DEpeaksAnno <- JN_DEpeaksAnno %>%  #rename columns! 
  rename("FC" = "V5")
##JN_DEpeaksAnno files is the 3625 JN DAR Peaks (not genes!). want to filter by FC HERE before merging any gene stuff
JN_DEpeaks_DAR <- JN_DEpeaksAnno %>% filter(FC > 1 | FC < -1)  #filter by abs value; works great!  

#JN_DEpeaksAnno_DAR_summary <- JN_DEpeaks_DAR %>%    #remove duplicates by summing up FC and averaging out per gene
 # group_by(SYMBOL) %>%                          #doing this way basically adds the meanFC column to the end 
 # mutate(mean_FC = mean(FC)) %>%
 # ungroup()

JN_DEpeaksAnno_DAR_distinct<- JN_DEpeaks_DAR %>% distinct(SYMBOL, .keep_all = TRUE) #removes any genes that have multiples from ATAC anno
##OLD_____JN_DEpeaks_filtered <- JN_DEpeaksAnno_distinct %>% filter(mean_FC > 1 | mean_FC < -1)  #filter by abs value; works great!  

JN_ATAC_down<- JN_DEpeaksAnno_DAR_distinct %>% filter(FC > 1) #atac
JN_ATAC_up<- JN_DEpeaksAnno_DAR_distinct %>% filter(FC< -1) #atac



#BER
#need to annotate for the gene 
setwd("D:/5_R_generated_Data/BER/3. 2022 ATAC")
BER_DEpeaksAnno <- annotatePeak("BER_NP_DEpeaks.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0


BER_DEpeaksAnno<-as.data.frame(BER_DEpeaksAnno) #as data frame to manipulate easier
BER_DEpeaksAnno <- BER_DEpeaksAnno %>%  #rename columns! 
  rename("FC" = "V5")
##BER_DEpeaksAnno files is the 14072 BER DAR Peaks (not genes!). want to filter by FC HERE before merging any gene stuff
BER_DEpeaks_DAR <- BER_DEpeaksAnno %>% filter(FC > 1 | FC < -1)  #filter by abs value; works great!  

#BER_DEpeaksAnno_DAR_summary <- BER_DEpeaks_DAR %>%  #remove duplicates by summing up FC and averaging out per gene
 # group_by(SYMBOL) %>%   #doing this way basically adds the meanFC column to the end 
 # mutate(mean_FC = mean(V5)) %>%
#  ungroup()

BER_DEpeaksAnno_DAR_distinct<- BER_DEpeaks_DAR %>% distinct(SYMBOL, .keep_all = TRUE) #removes any genes that have multiples from ATAC anno
#BER_DEpeaks_filtered <- BER_DEpeaksAnno_distinct %>% filter(mean_FC > 1 | mean_FC < -1)  #filter by abs value; works great!  

BER_ATAC_down<- BER_DEpeaksAnno_DAR_distinct %>% filter(FC > 1) #atac
BER_ATAC_up<- BER_DEpeaksAnno_DAR_distinct %>% filter(FC < -1) #atac
############

JN_ATAC_down_simp <- JN_ATAC_down %>% dplyr::select("seqnames", "start", "end", "V4")
JN_ATAC_up_simp <- JN_ATAC_up %>% dplyr::select("seqnames", "start", "end", "V4")

BER_ATAC_down_simp <- BER_ATAC_down %>% dplyr::select("seqnames", "start", "end", "V4")
BER_ATAC_up_simp <- BER_ATAC_up %>% dplyr::select("seqnames", "start", "end", "V4")

JN_ATAC_down_simp <- makeGRangesFromDataFrame(JN_ATAC_down_simp)
JN_ATAC_up_simp <- makeGRangesFromDataFrame(JN_ATAC_up_simp)
BER_ATAC_down_simp <- makeGRangesFromDataFrame(BER_ATAC_down_simp)
BER_ATAC_up_simp <- makeGRangesFromDataFrame(BER_ATAC_up_simp)

JNBER_atac_down <- GenomicRanges::intersect(JN_ATAC_down_simp, BER_ATAC_down_simp)
WT1_JNBER_atac_down <- GenomicRanges::intersect(WT1, JNBER_atac_down)
WT1_JNBER_atac_down_Anno <- annotatePeak(WT1_JNBER_atac_down, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0
WT1_JNBER_atac_down_Anno <- as.data.frame(WT1_JNBER_atac_down_Anno)


JNBER_atac_up <- GenomicRanges::intersect(JN_ATAC_up_simp, BER_ATAC_up_simp)
WT1_JNBER_atac_up <- GenomicRanges::intersect(WT1, JNBER_atac_up)
WT1_JNBER_atac_up_anno <- annotatePeak(WT1_JNBER_atac_up, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0
WT1_JNBER_atac_up_anno <- as.data.frame(WT1_JNBER_atac_up_anno)

#####
setwd("D:/3_Cut&TAG_Files/hg38Aligned/4_macs2Output/1_peakFiles/BER")
BER_H3K4me1<-import.bed("BER_H3K4me1_hg38.bed")
BER_H3K4me3<-import.bed("BER_H3K4me3_hg38.bed")
BER_H3K27ac<-import.bed("BER_H3K27ac_hg38.bed")
BER_H3K27me3<-import.bed("BER_H3K27me3_hg38.bed")

setwd("D:/3_Cut&TAG_Files/hg38Aligned/4_macs2Output/1_peakFiles/JN")
JN_H3K4me1<-import.bed("JN_H3K4me1_hg38.bed")
JN_H3K4me3<-import.bed("JN_H3K4me3_hg38.bed")
JN_H3K27ac<-import.bed("JN_H3K27ac_hg38.bed")
JN_H3K27me3<-import.bed("JN_H3K27me3_hg38.bed")

setwd("D:/3_Cut&TAG_Files/hg38Aligned/4_macs2Output/1_peakFiles/BOD")
BOD_H3K4me1<-import.bed("BOD_H3K4me1_hg38.bed")
BOD_H3K4me3<-import.bed("BOD_H3K4me3_hg38.bed")
BOD_H3K27ac<-import.bed("BOD_H3K27ac_hg38.bed")
BOD_H3K27me3<-import.bed("BOD_H3K27me3_hg38.bed")


JNBER_H3K4me1 <- GenomicRanges::intersect(JN_H3K4me1,BER_H3K4me1)
JNBER_H3K4me3 <- GenomicRanges::intersect(JN_H3K4me3, BER_H3K4me3)
JNBER_H3K27ac <- GenomicRanges::intersect(JN_H3K27ac, BER_H3K27ac)
JNBER_H3K27me3 <- GenomicRanges::intersect(JN_H3K27me3, BER_H3K27me3)

JNBERBOD_H3K4me1 <- GenomicRanges::intersect(JNBER_H3K4me1, BOD_H3K4me1)
JNBERBOD_H3K4me3 <- GenomicRanges::intersect(JNBER_H3K4me3, BOD_H3K4me3)
JNBERBOD_H3K27ac <- GenomicRanges::intersect(JNBER_H3K27ac, BOD_H3K27ac)
JNBERBOD_H3K27me3 <- GenomicRanges::intersect(JNBER_H3K27me3, BOD_H3K27me3)

H3K27ac_H3K4me1 <- GenomicRanges::intersect(JNBERBOD_H3K27ac, JNBERBOD_H3K4me1)
H3K27ac_H3K4me1_df <- as.data.frame(H3K27ac_H3K4me1)

promoters <- JNBERBOD_H3K4me3
promoters_df <- as.data.frame(promoters)
JN_promoters <-JN_H3K4me3
#####
dsrct_enh <- H3K27ac_H3K4me1
dsrct_enh_anno <- annotatePeak(dsrct_enh, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db")
dsrct_enh_anno_df <- as.data.frame(dsrct_enh_anno)
dsrct_enh_df <-as.data.frame(dsrct_enh)

WT1_dsrct_enh <- GenomicRanges::intersect(dsrct_enh, WT1)
WT1_dsrct_enh_anno <- annotatePeak(WT1_dsrct_enh, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0
WT1_dsrct_enh_anno <- as.data.frame(WT1_dsrct_enh_anno)

setwd("D:/6_Data Files")
hingorani <- read_xlsx("Hingorani, 2020 WT1 ChIP-Seq Data.xlsx")
hingorani <-hingorani %>% rename("peaks" = "PeakID..cmd.annotatePeaks.pl.EWS.WT1_rep1_testp.7_macspeaks.txt.hg19.")
hingorani_GRanges <- hingorani %>% dplyr::select("Chr", "Start", "End", "peaks")
hingorani_GRanges <- makeGRangesFromDataFrame(hingorani_GRanges)

WT1motifs_hingorani <- GenomicRanges::intersect(WT1, hingorani_GRanges)
WT1motifs_hingorani_df <- as.data.frame(WT1motifs_hingorani)
WT1_motifs_hingorani_anno<- annotatePeak(WT1motifs_hingorani, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db")
WT1_motifs_hingorani_anno <- as.data.frame(WT1_motifs_hingorani_anno)


dsrct_enh_hingorani <- GenomicRanges::intersect(dsrct_enh, hingorani_GRanges)
dsrct_enh_hingorani_df <- as.data.frame(dsrct_enh_hingorani)
dsrct_promoters_hingorani <- GenomicRanges::intersect(promoters, hingorani_GRanges)
dsrct_promoters_hingorani_df <-as.data.frame(dsrct_promoters_hingorani)
dsrct_promoters_hingorani_anno <- annotatePeak(dsrct_promoters_hingorani, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db")


dsrct_JN_promoters_hingorani <- GenomicRanges::intersect(JN_promoters, hingorani_GRanges)
dsrct_JN_promoters_hingorani_df <-as.data.frame(dsrct_JN_promoters_hingorani)
dsrct_JN_promoters_hingorani_anno <- annotatePeak(dsrct_JN_promoters_hingorani, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db")



