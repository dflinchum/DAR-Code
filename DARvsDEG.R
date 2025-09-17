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

#test hi
##2022 ATAC Data; Differentially accessible region stuff

#JN
#need to annotate for the gene 
setwd("D:/1_ATAC Files/ATAC 2022 Files")
JN_DEpeaksAnno <- annotatePeak("JN_NP_DEpeaks.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0

plotDistToTSS(JN_DEpeaksAnno,
              title="JN DAR: Distribution of transcription factor-binding loci relative to TSS") #have to do this here because needs to be AnnoObject thing


JN_DEpeaksAnno<-as.data.frame(JN_DEpeaksAnno) #as data frame to manipulate easier
JN_DEpeaksAnno <- JN_DEpeaksAnno %>%  #rename columns! 
  rename("FC" = "V5") #I SWEAR TO GOD THIS IS DIFFERENT EVERY TIME ON HOW TO CHANGE THE NAME

##JN_DEpeaksAnno files is the 14072 JN DAR Peaks (not genes!). want to filter by FC HERE before merging any gene stuff
JN_DEpeaks_DAR <- JN_DEpeaksAnno %>% filter(FC > 1 | FC < -1)  #filter by abs value; works great, nothing between -.99 to .99!  

#JN_DEpeaks_DAR <- JN_DEpeaks_DAR %>% filter(distanceToTSS > -10000 & distanceToTSS < 10000)  #filter by abs value; works great!  


JN_DEpeaksAnno_DAR_summary <- JN_DEpeaks_DAR %>%  #remove duplicates by summing up FC and averaging out per gene
  group_by(SYMBOL) %>%   #doing this way basically adds the meanFC column to the end 
  mutate(mean_FC = mean(FC)) %>%
  ungroup()

JN_DEpeaksAnno_DAR_distinct<- JN_DEpeaksAnno_DAR_summary %>% distinct(SYMBOL, .keep_all = TRUE) #removes any genes that have multiples from ATAC anno
JN_DEpeaksAnno_DAR_distinct_omit <- na.omit(JN_DEpeaksAnno_DAR_distinct)

#JN_DEpeaks_filtered <- JN_DEpeaksAnno_distinct %>% filter(mean_FC > 0 | mean_FC < 0)  #filter by abs value; works great!  

JN_DAR_downreg<- JN_DEpeaks_DAR %>% filter(FC < -1) #atac
JN_DAR_upreg<- JN_DEpeaks_DAR %>% filter(FC > 1) #atac

JN_DAR_downreg_omitNA <- na.omit(JN_DAR_downreg)
JN_DAR_upreg_omitNA <- na.omit(JN_DAR_upreg)

JN_DAR_downreg_surv <- merge(survival, JN_DAR_downreg, by.x="Gene", by.y="SYMBOL")
JN_DAR_downreg_surv <- JN_DAR_downreg_surv[,c(1,2,9,18)]

JN_DAR_upreg_surv <- merge(survival, JN_DAR_upreg, by.x="Gene", by.y="SYMBOL")
JN_DAR_upreg_surv<- JN_DAR_upreg_surv[,c(1,2,9,18)]



JN_DAR_gene_list <- list("FGF17, FBXO21, ZDHHC22, OTOGL") #these are open DARs w/sig survival and upregulated by EWSR1-WT1 


#BER
#need to annotate for the gene 
setwd("D:/1_ATAC Files/ATAC 2022 Files")
BER_DEpeaksAnno <- annotatePeak("BER_NP_DEpeaks.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0

plotDistToTSS(BER_DEpeaksAnno,
              title="BER DAR: Distribution of transcription factor-binding loci relative to TSS") #have to do this here because needs to be AnnoObject thing


BER_DEpeaksAnno<-as.data.frame(BER_DEpeaksAnno) #as data frame to manipulate easier
BER_DEpeaksAnno <- BER_DEpeaksAnno %>%  #rename columns! 
  rename("FC" = "V5")

##BER_DEpeaksAnno files is the 14072 BER DAR Peaks (not genes!). want to filter by FC HERE before merging any gene stuff
BER_DEpeaks_DAR <- BER_DEpeaksAnno %>% filter(FC > 1 | FC < -1)  #filter by abs value; works great, nothing between -.99 to .99!  

#BER_DEpeaks_DAR <- BER_DEpeaks_DAR %>% filter(distanceToTSS > -10000 & distanceToTSS < 10000)  #filter by abs value; works great!  


BER_DEpeaksAnno_DAR_summary <- BER_DEpeaks_DAR %>%  #remove duplicates by summing up FC and averaging out per gene
  group_by(SYMBOL) %>%   #doing this way basically adds the meanFC column to the end 
  mutate(mean_FC = mean(FC)) %>%
  ungroup()

BER_DEpeaksAnno_DAR_distinct<- BER_DEpeaksAnno_DAR_summary %>% distinct(SYMBOL, .keep_all = TRUE) #removes any genes that have multiples from ATAC anno
BER_DEpeaksAnno_DAR_distinct_omit <- na.omit(BER_DEpeaksAnno_DAR_distinct)

#BER_DEpeaks_filtered <- BER_DEpeaksAnno_distinct %>% filter(mean_FC > 0 | mean_FC < 0)  #filter by abs value; works great!  

BER_DAR_downreg<- BER_DEpeaks_DAR %>% filter(FC < -1) #atac
BER_DAR_upreg<- BER_DEpeaks_DAR %>% filter(FC > 1) #atac

BER_DAR_downreg_omitNA <- na.omit(BER_DAR_downreg)
BER_DAR_upreg_omitNA <-na.omit(BER_DAR_upreg)

BER_DAR_downreg_surv <- merge(survival, BER_DAR_downreg, by.x="Gene", by.y="SYMBOL")
BER_DAR_downreg_surv <- BER_DAR_downreg_surv[,c(1,2,9,18)]

BER_DAR_upreg_surv <- merge(survival, BER_DAR_upreg, by.x="Gene", by.y="SYMBOL")
BER_DAR_upreg_surv<- BER_DAR_upreg_surv[,c(1,2,9,18)]

BER_DAR_gene_list <- list("SAMD11, LPL, ARC, BDKRB2, GFRA2, SYT2, ZMAT3, FGF17, ITPR1, RAI2, OTOGL, JAK1, SRGAP3") #these are open DARs w/sig survival and upregulated by EWSR1-WT1 

###DOOOOOOOOOOO THIS ABOVE AND SEE WHERE WE LOSE THE FOXA2??????
##its p value? is shit.....

##get DAR upreg/downreg groups
JNBER_upreg <- merge(JN_DAR_upreg_omitNA, BER_DAR_upreg_omitNA, by="SYMBOL")
JNBER_upreg_surv <- merge(survival, JNBER_upreg, by.x="Gene", by.y="SYMBOL")

JNBER_downreg <- merge(JN_DAR_downreg_omitNA, BER_DAR_downreg_omitNA, by="SYMBOL")
JNBER_downreg_surv <- merge(survival, JNBER_downreg, by.x="Gene", by.y="SYMBOL")

#comparing DAR in -/+ dox conditions 
#opening
venn.diagram(x=list(JN_DAR_upreg_omitNA$SYMBOL, BER_DAR_upreg_omitNA$SYMBOL), 
             category.names=c("BER-DSRCT", "JN-DSRCT"), 
             filename='JN_BER_DAR_opening_Venn.png', 
             output=FALSE, 
             scaled=FALSE,
             cex=c(1.5),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = 1.5,
             margin = .1,
             #lwd = 2,
             fill = c("blue", "green")
)

#closing
venn.diagram(x=list(JN_DAR_downreg_omitNA$SYMBOL, BER_DAR_downreg_omitNA$SYMBOL), 
             category.names=c("BER-DSRCT", "JN-DSRCT"), 
             filename='JN_BER_DAR_closing_Venn.png', 
             output=FALSE, 
             scaled=FALSE,
             cex=c(1.5),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = 1.5,
             margin = .1,
             #lwd = 2,
             fill = c("blue", "green")
)

##this one works
##maybe sort the DARs by closest TSS? and see how that turns out 
venn.diagram(x=list(JN_DAR_downreg_omitNA$SYMBOL, BER_DAR_downreg_omitNA$SYMBOL, JN_DAR_upreg_omitNA$SYMBOL, BER_DAR_upreg_omitNA$SYMBOL), 
             category.names=c("JN DAR Downregulated", "BER DAR Downregulated", "JN DAR Upregulated", "BER DAR Upregulated"), 
             filename='JN_BER_DAR_up&down_Venn.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.27, 0.27, 0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green", "red", "yellow")
)
##
JNBER_atac_up <- merge(JN_ATAC_up, BER_ATAC_up, by="SYMBOL")

JNBER_atac_down <- merge(JN_ATAC_down, BER_ATAC_down, by="SYMBOL")

###making graphs of up/down etc
JN_ATAC_down <- JN_ATAC_down %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))
JN_ATAC_down_count <- JN_ATAC_down %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

JN_ATAC_down_count <- JN_ATAC_down_count %>% 
  rename("Annotation" = "annotation",
         "Peaks" = "count")


data2 <- JN_ATAC_down_count                                                 # Replicate original data
data2$Annotation <- factor(data2$Annotation,                                    # Factor levels in decreasing order
                           levels = data2$Annotation[order(data2$Peaks, decreasing = TRUE)])

test2 <- ggplot(data=data2, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen") +
  theme(text = element_text(size = 25),plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"))

test2 + ggtitle("JN atac down Differentially Accessible Region Distribution")


JN_ATAC_up <- JN_ATAC_up %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))
JN_ATAC_up_count <- JN_ATAC_up %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

JN_ATAC_up_count <- JN_ATAC_up_count %>% 
  rename("Annotation" = "annotation",
         "Peaks" = "count")


data3 <- JN_ATAC_up_count                                                 # Replicate original data
data3$Annotation <- factor(data3$Annotation,                                    # Factor levels in decreasing order
                           levels = data3$Annotation[order(data3$Peaks, decreasing = TRUE)])

test3 <- ggplot(data=data3, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen") +
  theme(text = element_text(size = 25),plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"))

test3 + ggtitle("JN atac up Differentially Accessible region Distribution")

BER_ATAC_down <- BER_ATAC_down %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))
BER_ATAC_down_count <- BER_ATAC_down %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

BER_ATAC_down_count <- BER_ATAC_down_count %>% 
  rename("Annotation" = "annotation",
         "Peaks" = "count")


data4 <- BER_ATAC_down_count                                                 # Replicate original data
data4$Annotation <- factor(data4$Annotation,                                    # Factor levels in decreasing order
                           levels = data4$Annotation[order(data4$Peaks, decreasing = TRUE)])

test4 <- ggplot(data=data4, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen") +
  theme(text = element_text(size = 25),plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"))

test4 + ggtitle("BER atac down Differentially Accessible Peak Distribution")


BER_ATAC_up <- BER_ATAC_up %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))
BER_ATAC_up_count <- BER_ATAC_up %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

BER_ATAC_up_count <- BER_ATAC_up_count %>% 
  rename("Annotation" = "annotation",
         "Peaks" = "count")


data5 <- BER_ATAC_up_count                                                 # Replicate original data
data5$Annotation <- factor(data5$Annotation,                                    # Factor levels in decreasing order
                           levels = data5$Annotation[order(data5$Peaks, decreasing = TRUE)])

test5 <- ggplot(data=data5, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen") +
  theme(text = element_text(size = 25),plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"))

test5 + ggtitle("BER atac up Differentially Accessible Peak Distribution")


##RNA-Seq DATA
#JN
setwd("D:/2_RNA Seq_Files/Excel Sheets/DEG")
JN_shWT1<-read.xlsx("JNshWT1_rsem_DGE_results.xlsx")
JN_shWT1_filtered<- JN_shWT1 %>% filter(log2FoldChange > 1 | log2FoldChange < -1, padj < 0.05)
JN_shWT1_filtered<-na.omit(JN_shWT1_filtered) #RNA
JN_up<-JN_shWT1_filtered %>% filter(log2FoldChange< -1, padj < 0.05) #RNA
JN_down<-JN_shWT1_filtered %>% filter(log2FoldChange>1, padj < 0.05) #RNA


#BER
BER_shWT1<-read.xlsx("Lee_BER_shWT1_4Days_rsem_DGE_results.xlsx")
BER_shWT1_filtered<- BER_shWT1 %>% filter(log2FoldChange > 1 | log2FoldChange < -1, padj < 0.05)
BER_shWT1_filtered<-na.omit(BER_shWT1_filtered) #RNA
BER_up<-BER_shWT1_filtered %>% filter(log2FoldChange< -1, padj < 0.05) #RNA  
BER_down<-BER_shWT1_filtered %>% filter(log2FoldChange>1, padj < 0.05) #RNA




##Venns

#this works   
setwd("D:/5_R_generated_Data")
venn.diagram(x=list(BER_DEpeaks_distinct$SYMBOL, BER_shWT1_filtered$gene_name), 
             category.names=c("BER DAR", "BER DEG"), 
             filename='BER_DAR_DEGnew.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green")
)

venn.diagram(x=list(BER_ATAC_down$SYMBOL, BER_down$gene_name), 
             category.names=c("BER Down ATAC", "BER Downregulated RNA"), 
             filename='BER_ATAC_RNA_DOWNnew.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green")
)

venn.diagram(x=list(BER_ATAC_up$SYMBOL, BER_up$gene_name), 
             category.names=c("BER UP ATAC", "BER Upregulated RNA"), 
             filename='BER_ATAC_RNA_UPnew.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green")
)



venn.diagram(x=list(BER_ATAC_down$SYMBOL, BER_up$gene_name), 
             category.names=c("BER ATAC Down", "BER upregulated RNA"), 
             filename='BER_ATACdown_RNAup.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green")
)

venn.diagram(x=list(BER_ATAC_up$SYMBOL, BER_down$gene_name), 
             category.names=c("BER ATAC Up", "BER downregulated RNA"), 
             filename='BER_ATACup_RNAdown.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green")
)

#JN  
venn.diagram(x=list(JN_DEpeaks_distinct$SYMBOL, JN_shWT1_filtered$gene_name), 
             category.names=c("JN DAR", "JN DEG"), 
             filename='JN_DAR_DEGnew.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green")
)

venn.diagram(x=list(JN_ATAC_down$SYMBOL, JN_down$gene_name), 
             category.names=c("JN Down ATAC", "JN Downregulated RNA"), 
             filename='JN_ATAC_RNA_DOWNnew.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green")
)

venn.diagram(x=list(JN_ATAC_up$SYMBOL, JN_up$gene_name), 
             category.names=c("JN UP ATAC", "JN Upregulated RNA"), 
             filename='JN_ATAC_RNA_UPnew.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green")
)




venn.diagram(x=list(JN_ATAC_up$SYMBOL, JN_down$gene_name), 
             category.names=c("JN ATAC Up", "JN downregulated RNA"), 
             filename='JN_ATACup_RNAdown.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green")
)

venn.diagram(x=list(JN_ATAC_down$SYMBOL, JN_up$gene_name), 
             category.names=c("JN ATAC Down", "JN upregulated RNA"), 
             filename='JN_ATACdown_RNAup.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green")
)


##JN and BER DAR vs DAR, and DEG vs DEG; or use something from Justin and the paper for DEG stuff; but can use these conditions here for similar approach

JN_DEpeaks_DAR_noNA <- na.omit(JN_DEpeaks_DAR)
BER_DEpeaks_DAR_noNA <- na.omit(BER_DEpeaks_DAR)


JN_DEpeaksAnno_DAR_distinct_noNA<- na.omit(JN_DEpeaksAnno_DAR_distinct)
BER_DEpeaksAnno_DAR_distinct_noNA<- na.omit(BER_DEpeaksAnno_DAR_distinct)

venn.diagram(x=list(JN_DEpeaksAnno_DAR_distinct_noNA$SYMBOL, BER_DEpeaksAnno_DAR_distinct_noNA$SYMBOL), 
             category.names=c("BER", "JN"), 
             filename='JN_BER_DAR_compare.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green")
)



venn.diagram(x=list(JN_shWT1_filtered$gene_name, BER_shWT1_filtered$gene_name), 
             category.names=c("BER", "JN"), 
             filename='JN_BER_DEG_compare.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green")
)



JN_BER_DAR<- merge(JN_DEpeaksAnno_DAR_distinct_noNA, BER_DEpeaksAnno_DAR_distinct_noNA, by.x="SYMBOL", by.y="SYMBOL" )
JN_BER_DEG<- merge(JN_shWT1_filtered, BER_shWT1_filtered, by.x="gene_name", by.y="gene_name")

JN_BER_DAR_DEG_list <- merge(JN_BER_DAR, JN_BER_DEG, by.x="SYMBOL", by.y="gene_name")
## write out/export the list and whichever columns we want for 
venn.diagram(x=list(JN_BER_DAR$SYMBOL, JN_BER_DEG$gene_name), 
             category.names=c("DAR", "DEG"), 
             filename='JN_BER_DAR_DEG_compare.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green")
)

write.xlsx(JN_BER_DAR_DEG_list, "JN_BER_DAR_DEG_newlist.xlsx")


##Merging 2022 DAR ATAC with UP/DOWN (DEG) shWT1 RNA data
#should see the same basic numbers as in venns 
JN_DAR_DEG<-merge(JN_DEpeaks_filtered, JN_shWT1_filtered, by.x="SYMBOL", by.y="gene_name")

JN_ATAC_RNA_down<-merge(JN_ATAC_down, JN_down, by.x="SYMBOL", by.y="gene_name")
JN_ATAC_RNA_up<-merge(JN_ATAC_up, JN_up, by.x="SYMBOL", by.y="gene_name")

JN_ATACup_RNAdown<-merge(JN_ATAC_up, JN_down, by.x="SYMBOL", by.y="gene_name")
JN_ATACdown_RNAup<-merge(JN_ATAC_down, JN_up, by.x="SYMBOL", by.y="gene_name")

BER_DAR_DEG<-merge(BER_DEpeaksAnno_DAR_distinct, BER_shWT1_filtered, by.x="SYMBOL", by.y="gene_name")

BER_ATAC_RNA_down<-merge(BER_ATAC_down, BER_down, by.x="SYMBOL", by.y="gene_name")
BER_ATAC_RNA_up<-merge(BER_ATAC_up, BER_up, by.x="SYMBOL", by.y="gene_name")

BER_ATACdown_RNAup<-merge(BER_ATAC_down, BER_up, by.x="SYMBOL", by.y="gene_name")
BER_ATACup_RNAdown<-merge(BER_ATAC_up, BER_down, by.x="SYMBOL", by.y="gene_name")

##WRITE THESE OUT (from above)
#JN
write.xlsx(JN_DAR_DEG, "JN_DAR_DEG.xlsx")
write.xlsx(JN_ATAC_RNA_down, "JN_ATAC_RNA_down.xlsx")
write.xlsx(JN_ATAC_RNA_up, "JN_ATAC_RNA_Up.xlsx")
write.xlsx(JN_ATACdown_RNAup, "JN_ATACdown_RNAup.xlsx")
write.xlsx(JN_ATACup_RNAdown, "JN_ATACup_RNAdown.xlsx")


write.xlsx(BER_DAR_DEG, "BER_DAR_DEG.xlsx")
write.xlsx(BER_ATAC_RNA_down, "BER_ATAC_RNA_down.xlsx")
write.xlsx(BER_ATAC_RNA_up, "BER_ATAC_RNA_Up.xlsx")
write.xlsx(BER_ATACdown_RNAup, "BER_ATACdown_RNAup.xlsx")
write.xlsx(BER_ATACup_RNAdown, "BER_ATACup_RNAdown.xlsx")



#JN and BER DAR DEG Overlap from gene names 

JN_BER_DAR_DEG<-merge(JN_DAR_DEG, BER_DAR_DEG, by="SYMBOL")

write.xlsx(JN_BER_DAR_DEG, "JN_BER_DAR_DEG.xlsx") 


venn.diagram(x=list(JN_DAR_DEG$SYMBOL, BER_DAR_DEG$SYMBOL), 
             category.names=c("BER DAR_DEG", "JN DAR_DEG"), 
             filename='JN_BER_DARDEGtest.png', 
             output=FALSE, 
             scaled=FALSE,
             #cex=c(1,1,1,1),
             #cat.pos = c(-35, -10, 80, 10),
             cat.dist = c(0.15, 0.15),
             cat.cex = .75,
             margin = .2,
             #lwd = 2,
             fill = c("blue", "green")
)




##GSEA DAR DEG

#JN DAR DEG GSEA
#playing around, got the graphs to produce but p value is still weird, needs more work 

setwd("D:/5_R_generated_Data/JN/4. 2020 ATAC x Cut&Tag/Melody hg38")

JN_DAR_DEG_FC=JN_DAR_DEG$FC
names(JN_DAR_DEG_FC)<-as.character(JN_DAR_DEG$geneId)
JN_DAR_DEG_FC <- sort(JN_DAR_DEG_FC, decreasing = TRUE)

JN_DAR_DEG_gsea <- gseGO(geneList = JN_DAR_DEG_FC, 
                     ont = 'BP', 
                     OrgDb = org.Hs.eg.db, 
                     keyType = "ENTREZID",
                     #nPerm = 1000, 
                     minGSSize = 20, 
                     pvalueCutoff = 0.05,
                     verbose = FALSE)
dotplot(JN_DAR_DEG_gsea, showCategory=15) + ggtitle("JN DAR DEG")
JN_DAR_DEG_gsea_results <- JN_DAR_DEG_gsea@result
write.xlsx(JJN_DAR_DEG_gsea_results, "JN_DAR_DEG_gsea_results.xlsx", quote=F)

JN_DAR_DEG_gseDO<-gseDO(JN_DAR_DEG_FC,
                     pvalueCutoff = 1)
dotplot(JN_DAR_DEG_gseDO, showCategory=15) + ggtitle("JN: ATAC/EnhancerMarks/shWT1 RNA_gseDO")
JN_DAR_DEG_gseDO_result <- JN_DAR_DEG_gseDO@result
write.xlsx(JN_gseaDO_enh_result, "JN_gseaDO_enh_result.xlsx", quote=F)


#JN all peaks gsea
jn_allpeaks<-import.bed("JN_NP_allpeaks.bed")
jn_allpeaks_anno<- annotatePeak(jn_allpeaks, TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db")
jn_allpeaks_anno<- as.data.frame(jn_allpeaks_anno)
jn_allpeaks_anno <- jn_allpeaks_anno %>%  #rename columns! 
  rename("FC" = "score")
jn_allpeaks_anno_FC=jn_allpeaks_anno$FC
names(jn_allpeaks_anno_FC)<-as.character(jn_allpeaks_anno$geneId)
jn_allpeaks_anno_FC <- sort(jn_allpeaks_anno_FC, decreasing = TRUE)

jn_allpeaks_gsea <- gseGO(geneList = jn_allpeaks_anno_FC, 
                         ont = 'BP', 
                         OrgDb = org.Hs.eg.db, 
                         keyType = "ENTREZID",
                         #nPerm = 1000, 
                         minGSSize = 20, 
                         pvalueCutoff = 0.05,
                         verbose = FALSE)
dotplot(jn_allpeaks_gsea, showCategory=15) + ggtitle("JN all peaks")





###JUNK CODE; TESTING AND PLAYING WITH DIFFERENT CONDITIONS
##need to clean up eventually 
##not needed at the moment
aggtest <- JN_DE_H3K27me3 %>% group_by(geneid) %>% 
  summarise(mean_foldchange=mean(foldchange),
            .groups = 'drop')






JN_up <- JNshWT1[(JNshWT1$padj<0.05 & JNshWT1$log2FoldChange<0),"gene_name"]
test3 <- JN_BER_BOD_promoter_merge[(JN_BER_BOD_promoter_merge$padj.x<0.05), "SYMBOL"]
head(test2)
test3<-as.data.frame(test3)
test3<-na.omit(test3)

#try to get jn up in non vector omg
dane16<-my_df %>% filter(d == "A", a > 0.5) 

#worked completely normally :)
otter<-JN_shWT1 %>% filter(log2FoldChange<0, padj > 0.5) 

####END JUNK



dsrct_targets <- read.table("dsrct_enhancer_targets.txt", header = TRUE, sep = ",", dec = ".")

dsrct_targets_DARDEG <- merge(JN_BER_DAR_DEG_list, dsrct_targets, by.x="SYMBOL", by.y="Gene.Symbol")
