






#load DSRCT consensus peaks
setwd("D:/1_ATAC Files/ATAC 2022 Files")
JN_DEsamples <- read.csv("JN_DSRCTshWT1_Differential_Peak_Binding_NP.csv")

#Create sets with differentially expressed peaks
JN_upreg <- JN_DEsamples[JN_DEsamples$FDR<0.05 & JN_DEsamples$Fold > 1 & JN_DEsamples$p.value <0.05, c(2,3,4,13,6)]
JN_downreg <- JN_DEsamples[JN_DEsamples$FDR<0.05 & JN_DEsamples$Fold < -1 & JN_DEsamples$p.value <0.05, c(2,3,4,13,6)]

JN_upreg <- makeGRangesFromDataFrame(JN_upreg, keep.extra.columns = TRUE)
JN_downreg <- makeGRangesFromDataFrame(JN_downreg, keep.extra.columns = TRUE)

export.bed(JN_upreg, "JN_upreg_HOMER.bed")
export.bed(JN_downreg, "JN_downreg_HOMER.bed")

#BER
BER_DEsamples <- read.csv("BER_DSRCTshWT1_Differential_Peak_Binding_NP.csv")

#Create sets with differentially expressed peaks
BER_upreg <- BER_DEsamples[BER_DEsamples$FDR<0.05 & BER_DEsamples$Fold > 1 & BER_DEsamples$p.value <0.05, c(2,3,4,13,6)]
BER_downreg <- BER_DEsamples[BER_DEsamples$FDR<0.05 & BER_DEsamples$Fold < -1 & BER_DEsamples$p.value <0.05, c(2,3,4,13,6)]



BER_upreg <- makeGRangesFromDataFrame(BER_upreg, keep.extra.columns = TRUE)
BER_downreg <- makeGRangesFromDataFrame(BER_downreg, keep.extra.columns = TRUE)

export.bed(BER_upreg, "BER_upreg_HOMER.bed")
export.bed(BER_downreg, "BER_downreg_HOMER.bed")


