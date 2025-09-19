library(openxlsx)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(ChIPpeakAnno)
library(ChIPseeker)
library(dplyr)

setwd("D:/6_Data Files/FANTOM")
fantom <- import.bed("Fantom_hg38_updated.bed") #enhancer peaks

#BER
BER_minus_WT1_peaks <- read.xlsx("BER_MinusDox_WT1_motif_peaks.xlsx") #load in only peaks that have WT1 motifs for each condition
BER_plus_WT1_peaks <- read.xlsx("BER_PlusDox_WT1_motif_peaks.xlsx")

BER_minus_WT1_peaks <- BER_minus_WT1_peaks %>% dplyr::select(Chr, Start, End, Strand) #editing format
BER_plus_WT1_peaks <- BER_plus_WT1_peaks %>% select(Chr, Start, End, Strand)

BER_minus_WT1_peaks <- makeGRangesFromDataFrame(BER_minus_WT1_peaks) #GRanges for intersection
BER_plus_WT1_peaks <- makeGRangesFromDataFrame(BER_plus_WT1_peaks)

BER_DEsamples <- read.csv("BER_DSRCTshWT1_Differential_Peak_Binding_NP.csv")
BER_DEpeaks <- BER_DEsamples[BER_DEsamples$FDR<0.05,c(2,3,4,13,10)]
BER_negDoxpeaks <- BER_DEsamples[BER_DEsamples$FDR<0.05 & BER_DEsamples$Fold > 0,c(2,3,4,13,10)]
BER_posDoxpeaks <- BER_DEsamples[BER_DEsamples$FDR<0.05 & BER_DEsamples$Fold < 0,c(2,3,4,13,10)]

BER_negDoxpeaks <- makeGRangesFromDataFrame(BER_negDoxpeaks)
BER_posDoxpeaks <- makeGRangesFromDataFrame(BER_posDoxpeaks)

BER_minus_WT1_peaks <-GenomicRanges::intersect(BER_posDoxpeaks, BER_minus_WT1_peaks, ignore.strand=TRUE) #have to make sure the FDR is appropriate/actually have DE peaks 
BER_plus_WT1_peaks <- GenomicRanges::intersect(BER_negDoxpeaks, BER_plus_WT1_peaks, ignore.strand=TRUE) #so load in DE peaks from chipseeker and intersect



fantom_BERminus_WT1 <- GenomicRanges::intersect(fantom, BER_minus_WT1_peaks, ignore.strand=TRUE) #Fantom enhancers intersected with each condition of WT1 peaks
fantom_BERplus_WT1 <- GenomicRanges::intersect(fantom, BER_plus_WT1_peaks, ignore.strand=TRUE)
fantom_BERminus_WT1_df <- as.data.frame(fantom_BERminus_WT1)
fantom_BERplus_WT1_df <- as.data.frame(fantom_BERplus_WT1)



export.bed(fantom_BERminus_WT1, "BER_WT1_minus_fantom.bed") #wrote this out and used it in bedtools for finding the unique regions
export.bed(fantom_BERplus_WT1, "BER_WT1_plus_fantom.bed")


BER_minusonly_WT1_enhancers <- import.bed("non_overlapping_in_BERminus.bed") #importing the unique regions from each after bedtools usage
BER_plusonly_WT1_enhancers <- import.bed("non_overlapping_in_BERplus.bed")

BER_minusonly_WT1_enhancers <- annotatePeak(BER_minusonly_WT1_enhancers, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
BER_plusonly_WT1_enhancers <- annotatePeak(BER_plusonly_WT1_enhancers, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db")

BER_minusonly_WT1_enhancers <- as.data.frame(BER_minusonly_WT1_enhancers)
BER_plusonly_WT1_enhancers <- as.data.frame(BER_plusonly_WT1_enhancers)



setwd("D:/3_Cut&TAG_Files/hg38Aligned/4_macs2Output/1_peakFiles/BER")
BER_H3K4me1<-import.bed("BER_H3K4me1_hg38.bed")
BER_H3K4me3<-import.bed("BER_H3K4me3_hg38.bed")
BER_H3K27ac<-import.bed("BER_H3K27ac_hg38.bed")
BER_H3K27me3<-import.bed("BER_H3K27me3_hg38.bed")



BER_H3K27ac_H3K4me1 <- GenomicRanges::intersect(BER_H3K27ac, BER_H3K4me1, ignore.strand=TRUE) #making "BER enhancers" from CUTTag stuff
BER_minusonly_WT1_enhancers_TAG <- GenomicRanges::intersect(BER_minusonly_WT1_enhancers, BER_H3K27ac_H3K4me1, ignore.strand=TRUE) #crossing CUTTag enh with our Fantom/WT1 peaks enhancers
BER_plusonly_WT1_enhancers_TAG <- GenomicRanges::intersect(BER_plusonly_WT1_enhancers, BER_H3K27ac_H3K4me1, ignore.strand=TRUE)

BER_minusonly_WT1_enhancers_TAG<- annotatePeak(BER_minusonly_WT1_enhancers_TAG, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
BER_plusonly_WT1_enhancers_TAG<- annotatePeak(BER_plusonly_WT1_enhancers_TAG, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db")

BER_minusonly_WT1_enhancers_TAG <- as.data.frame(BER_minusonly_WT1_enhancers_TAG)
BER_plusonly_WT1_enhancers_TAG <- as.data.frame(BER_plusonly_WT1_enhancers_TAG)


#JN
#JN
setwd("D:/1_ATAC Files/ATAC 2022 Files")
JN_minus_WT1_peaks <- read.xlsx("JN_MinusDox_WT1_motif_peaks.xlsx") #load in only peaks that have WT1 motifs for each condition
JN_plus_WT1_peaks <- read.xlsx("JN_PlusDox_WT1_motif_peaks.xlsx")

JN_minus_WT1_peaks <- JN_minus_WT1_peaks %>% select(Chr, Start, End, Strand) #editing format
JN_plus_WT1_peaks <- JN_plus_WT1_peaks %>% select(Chr, Start, End, Strand)

JN_minus_WT1_peaks <- makeGRangesFromDataFrame(JN_minus_WT1_peaks) #GRanges for intersection
JN_plus_WT1_peaks <- makeGRangesFromDataFrame(JN_plus_WT1_peaks)


#from Chipseeker
JN_DEsamples <- read.csv("JN_DSRCTshWT1_Differential_Peak_Binding_NP.csv")
#Create sets with differentially expressed peaks
JN_allpeaks <- JN_DEsamples[,c(2,3,4,13,10)]
JN_DEpeaks <- JN_DEsamples[JN_DEsamples$FDR<0.05,c(2,3,4,13,10)]
JN_negDoxpeaks <- JN_DEsamples[JN_DEsamples$FDR<0.05 & JN_DEsamples$Fold > 0,c(2,3,4,13,10)]
JN_posDoxpeaks <- JN_DEsamples[JN_DEsamples$FDR<0.05 & JN_DEsamples$Fold < 0,c(2,3,4,13,10)]
###

JN_negDoxpeaks <- makeGRangesFromDataFrame(JN_negDoxpeaks)
JN_posDoxpeaks <- makeGRangesFromDataFrame(JN_posDoxpeaks)
JN_minus_WT1_peaks <-GenomicRanges::intersect(JN_posDoxpeaks, JN_minus_WT1_peaks, ignore.strand=TRUE) #have to make sure the FDR is appropriate/actually have DE peaks 
JN_plus_WT1_peaks <- GenomicRanges::intersect(JN_negDoxpeaks, JN_plus_WT1_peaks, ignore.strand=TRUE) #so load in DE peaks from chipseeker and intersect



fantom_JNminus_WT1 <- GenomicRanges::intersect(fantom, JN_minus_WT1_peaks, ignore.strand=TRUE) #Fantom enhancers intersected with each condition of WT1 peaks
fantom_JNplus_WT1 <- GenomicRanges::intersect(fantom, JN_plus_WT1_peaks, ignore.strand=TRUE)
fantom_JNminus_WT1_df <- as.data.frame(fantom_JNminus_WT1)
fantom_JNplus_WT1_df <- as.data.frame(fantom_JNplus_WT1)



export.bed(fantom_JNminus_WT1, "JN_WT1_minus_fantom.bed") #wrote this out and used it in bedtools for finding the unique regions
export.bed(fantom_JNplus_WT1, "JN_WT1_plus_fantom.bed")


JN_minusonly_WT1_enhancers <- import.bed("non_overlapping_in_JNminus.bed") #importing the unique regions from each after bedtools usage
JN_plusonly_WT1_enhancers <- import.bed("non_overlapping_in_JNplus.bed") #CHECK TO MAKE SURE THE OUTPUT FROM BEDTOOLS AND THIS ARE THE SAME; RENAME THE BER ONES

JN_minusonly_WT1_enhancers <- annotatePeak(JN_minusonly_WT1_enhancers, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
JN_plusonly_WT1_enhancers <- annotatePeak(JN_plusonly_WT1_enhancers, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db")

JN_minusonly_WT1_enhancers <- as.data.frame(JN_minusonly_WT1_enhancers)
JN_plusonly_WT1_enhancers <- as.data.frame(JN_plusonly_WT1_enhancers)


setwd("D:/3_Cut&TAG_Files/hg38Aligned/4_macs2Output/1_peakFiles/JN")
JN_H3K4me1<-import.bed("JN_H3K4me1_hg38.bed")
JN_H3K4me3<-import.bed("JN_H3K4me3_hg38.bed")
JN_H3K27ac<-import.bed("JN_H3K27ac_hg38.bed")
JN_H3K27me3<-import.bed("JN_H3K27me3_hg38.bed")


JN_H3K27ac_H3K4me1 <- GenomicRanges::intersect(JN_H3K27ac, JN_H3K4me1, ignore.strand=TRUE) #making "JN enhancers" from CUTTag stuff
JN_minusonly_WT1_enhancers_TAG <- GenomicRanges::intersect(JN_minusonly_WT1_enhancers, JN_H3K27ac_H3K4me1, ignore.strand=TRUE) #crossing CUTTag enh with our Fantom/WT1 peaks enhancers
JN_plusonly_WT1_enhancers_TAG <- GenomicRanges::intersect(JN_plusonly_WT1_enhancers, JN_H3K27ac_H3K4me1, ignore.strand=TRUE)

JN_minusonly_WT1_enhancers_TAG<- annotatePeak(JN_minusonly_WT1_enhancers_TAG, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") 
JN_plusonly_WT1_enhancers_TAG<- annotatePeak(JN_plusonly_WT1_enhancers_TAG, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db")

JN_minusonly_WT1_enhancers_TAG <- as.data.frame(JN_minusonly_WT1_enhancers_TAG)
JN_plusonly_WT1_enhancers_TAG <- as.data.frame(JN_plusonly_WT1_enhancers_TAG)


common_in_both <- merge(JN_minusonly_WT1_enhancers_TAG, BER_minusonly_WT1_enhancers_TAG, by = "SYMBOL")


#checking out fantom and cut tag enhancer 
fantom_BER_enh <- GenomicRanges::intersect(fantom, BER_H3K27ac_H3K4me1, ignore.strand=TRUE)
fantom_JN_enh <- GenomicRanges::intersect(fantom, JN_H3K27ac_H3K4me1, ignore.strand=TRUE)
