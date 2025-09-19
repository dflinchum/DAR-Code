library(readxl)
library(openxlsx)
library(rtracklayer)
library(ChIPseeker)
library(ChIPpeakAnno)
library(VennDiagram)
library(dplyr)
library(clusterProfiler)
library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


#need JN and BER FC>1, padj or FDR <.05 in both ATAC and RNA

#RNA
setwd("D:/2_RNA Seq_Files/Excel Sheets")
fold_change <- read.xlsx("Regulation_JN_BER_SK2_BOD_KTS_ALL_GENES.xlsx")

BER_fold_change <- fold_change[,c(1,6,7)]
BER_fold_change <- BER_fold_change %>% dplyr::rename(Gene = gene_name) %>%
  na.omit()

JN_fold_change <- fold_change[,c(1,3,4)]
JN_fold_change <- JN_fold_change %>% dplyr::rename(Gene = gene_name) %>%
  na.omit()

JNBER_foldchange <- merge(JN_fold_change, BER_fold_change, by="Gene")

JNBER_upregulated <- JNBER_foldchange %>% filter(JN_shWT1_padj <0.05 & BER_shWT1_padj <0.05)
JNBER_upregulated <- JNBER_upregulated %>% filter(JNshWT1_log2FC > 1 & BERshWT1_log2FC > 1)

JNBER_downregulated <- JNBER_foldchange %>% filter(JN_shWT1_padj <0.05 & BER_shWT1_padj <0.05)
JNBER_downregulated <- JNBER_downregulated %>% filter(JNshWT1_log2FC < -1 & BERshWT1_log2FC < -1)


#JN
#need to annotate for the gene 
setwd("D:/1_ATAC Files/ATAC 2022 Files")
JN_DEpeaksAnno <- annotatePeak("JN_NP_DEpeaks.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0

JN_DEpeaksAnno<-as.data.frame(JN_DEpeaksAnno) #as data frame to manipulate easier
JN_DEpeaksAnno <- JN_DEpeaksAnno %>%  #rename columns! 
  rename("FC" = "V5") #I SWEAR TO GOD THIS IS DIFFERENT EVERY TIME ON HOW TO CHANGE THE NAME

##JN_DEpeaksAnno files is the 14072 JN DAR Peaks (not genes!). want to filter by FC HERE before merging any gene stuff
JN_DEpeaks_DAR <- JN_DEpeaksAnno %>% filter(FC > 1 | FC < -1)  #filter by abs value; works great, nothing between -.99 to .99!  


JN_DAR_closing<- JN_DEpeaks_DAR %>% filter(FC < -1) #atac
JN_DAR_opening<- JN_DEpeaks_DAR %>% filter(FC > 1) #atac


JN_DAR_closing_distinct <- JN_DAR_closing %>%
  group_by(SYMBOL) %>%
  slice_min(distanceToTSS, n = 1, with_ties = FALSE) %>%
  ungroup()

JN_DAR_opening_distinct <- JN_DAR_opening %>%
  group_by(SYMBOL) %>%
  slice_min(distanceToTSS, n = 1, with_ties = FALSE) %>%
  ungroup()    


#BER
#need to annotate for the gene 
setwd("D:/1_ATAC Files/ATAC 2022 Files")
BER_DEpeaksAnno <- annotatePeak("BER_NP_DEpeaks.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0

BER_DEpeaksAnno<-as.data.frame(BER_DEpeaksAnno) #as data frame to manipulate easier
BER_DEpeaksAnno <- BER_DEpeaksAnno %>%  #rename columns! 
  rename("FC" = "V5") #I SWEAR TO GOD THIS IS DIFFERENT EVERY TIME ON HOW TO CHANGE THE NAME

##BER_DEpeaksAnno files is the 14072 BER DAR Peaks (not genes!). want to filter by FC HERE before merging any gene stuff
BER_DEpeaks_DAR <- BER_DEpeaksAnno %>% filter(FC > 1 | FC < -1)  #filter by abs value; works great, nothing between -.99 to .99!  


BER_DAR_closing<- BER_DEpeaks_DAR %>% filter(FC < -1) #atac
BER_DAR_opening<- BER_DEpeaks_DAR %>% filter(FC > 1) #atac

BER_DAR_closing_distinct <- BER_DAR_closing %>% 
  group_by(SYMBOL) %>%
  slice_min(distanceToTSS, n= 1, with_ties = FALSE) %>%
  ungroup()
            
BER_DAR_opening_distinct <- BER_DAR_opening %>% 
  group_by(SYMBOL) %>%
  slice_min(distanceToTSS, n= 1, with_ties = FALSE) %>%
  ungroup()            
            


##JN-BER
JNBER_closing_distinct <- merge(JN_DAR_closing_distinct, BER_DAR_closing_distinct, by="SYMBOL")

JNBER_opening_distinct <- merge(JN_DAR_opening_distinct, BER_DAR_opening_distinct, by="SYMBOL")

#merge with RNA fold change
JNBER_closing_distinct_RNA_downreg <- merge(JNBER_closing_distinct, JNBER_downregulated, by.x="SYMBOL", by.y="Gene")
JNBER_opening_distinct_RNA_upreg <- merge(JNBER_opening_distinct, JNBER_upregulated, by.x="SYMBOL", by.y="Gene")

JNBER_closing_distinct_RNA_upreg <- merge(JNBER_closing_distinct, JNBER_upregulated, by.x="SYMBOL", by.y="Gene")
JNBER_opening_distinct_RNA_downreg <- merge(JNBER_opening_distinct, JNBER_downregulated, by.x="SYMBOL", by.y="Gene")

JNBER_closing_distinct_RNA_downreg_test <- JNBER_closing_distinct_RNA_downreg %>% distinct(SYMBOL, .keep_all = TRUE)
write.xlsx(JNBER_closing_distinct_RNA_downreg_test, "JNBER_closing_downregulated.xlsx")
JNBER_opening_distinct_RNA_upreg_test <- JNBER_opening_distinct_RNA_upreg %>% distinct(SYMBOL, .keep_all = TRUE)
write.xlsx(JNBER_opening_distinct_RNA_upreg_test, "JNBER_opening_upregulated.xlsx")
JNBER_closing_distinct_RNA_upreg_test <- JNBER_closing_distinct_RNA_upreg %>% distinct(SYMBOL, .keep_all = TRUE)
JNBER_opening_distinct_RNA_downreg_test <- JNBER_opening_distinct_RNA_downreg %>% distinct(SYMBOL, .keep_all = TRUE)


#graph
data <- data.frame(
  Regulation = rep(c("Repressed", "Induced"), each = 2),
  Accessibility = rep(c("Opening", "Closing"), 2),
  Genes = c(73, 24, 1, 76)
)

data$Regulation <- factor(data$Regulation, levels = c("Induced", "Repressed")) # Change the order here

ggplot(data, aes(x = Regulation, y = Genes, fill = Accessibility)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Genes), vjust = -0.5, position = position_dodge(width = 0.9), size = 10, fontface = "bold") +
  scale_fill_manual(values = c("darkred", "lightblue")) +
  labs(
    title = NULL,
    x = "",
    y = "Genes"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = c(0.5, 0.9),
    legend.title = element_blank(),
    panel.grid = element_blank(), # Add this line to remove gridlines
    legend.text = element_text(size = 30), # Increase legend font size
    axis.text = element_text(size = 30, face = "bold"),   # Increase axis text font size
    axis.title = element_text(size = 28, face="bold")  # Increase axis title font size
    )






