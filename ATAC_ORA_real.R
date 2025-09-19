library(tidygraph)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)
library(readxl)
library(VennDiagram)
library(openxlsx)


#load DSRCT consensus peaks
setwd("D:/1_ATAC Files/ATAC 2022 Files")
JN_DEsamples <- read.csv("JN_DSRCTshWT1_Differential_Peak_Binding_NP.csv")

#Create sets with differentially expressed peaks
JN_upreg <- JN_DEsamples[JN_DEsamples$FDR<0.05 & JN_DEsamples$Fold > 1 & JN_DEsamples$p.value <0.05, c(2,3,4,13,10,11,12)]
JN_downreg <- JN_DEsamples[JN_DEsamples$FDR<0.05 & JN_DEsamples$Fold < -1 & JN_DEsamples$p.value <0.05, c(2,3,4,13,10,11,12)]
JN_all_genes <- JN_DEsamples[,c(2,3,4,13,10,11,12)] #for background

write.table(JN_upreg, "JN_upreg.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
JN_upreg <- annotatePeak("JN_upreg.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0
JN_upreg <- as.data.frame(JN_upreg)

write.table(JN_downreg, "JN_downreg.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
JN_downreg <- annotatePeak("JN_downreg.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0
JN_downreg <- as.data.frame(JN_downreg)

write.table(JN_all_genes, "JN_all_genes.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
JN_all_genes <- annotatePeak("JN_all_genes.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0
JN_all_genes <- as.data.frame(JN_all_genes)

#JN RNA
setwd("D:/2_RNA Seq_Files/Excel Sheets")
fold_change <- read.xlsx("Regulation_JN_BER_SK2_BOD_KTS_ALL_GENES.xlsx")
JN_fold_change <- fold_change[,c(1,3)]
JN_fold_change <- JN_fold_change %>% dplyr::rename(Gene = gene_name, log2_fold_change =JNshWT1_log2FC) %>%
  na.omit()

#unique genes
JN_upreg_unique <- JN_upreg %>% 
  distinct(SYMBOL, .keep_all = TRUE)

JN_upreg_unique <- merge(JN_upreg_unique, JN_fold_change, by.x="SYMBOL", by.y="Gene")
JN_upreg_unique <- JN_upreg_unique %>% filter(log2_fold_change >1)

JN_downreg_unique <- JN_downreg %>% 
  distinct(SYMBOL, .keep_all = TRUE)

JN_downreg_unique <- merge(JN_downreg_unique, JN_fold_change, by.x="SYMBOL", by.y="Gene")
JN_downreg_unique <- JN_downreg_unique %>% filter(log2_fold_change < -1)


JN_all_genes_unique <- JN_all_genes %>% 
  distinct(SYMBOL, .keep_all = TRUE)

##trying to get genes only in one or the other
##seemed to work pretty well, still have decent # of genes
JN_upreg_unique_diff <- setdiff(JN_upreg_unique$SYMBOL, JN_downreg_unique$SYMBOL)


JN_upreg_unique_diff_only <- JN_upreg_unique[JN_upreg_unique$SYMBOL %in% JN_upreg_unique_diff, ]


JN_downreg_unique_diff <- setdiff(JN_downreg_unique$SYMBOL, JN_upreg_unique$SYMBOL)

JN_downreg_unique_diff_only <- JN_downreg_unique[JN_downreg_unique$SYMBOL %in% JN_downreg_unique_diff, ]


#pulling IDs/gene names only for analysis
JN_upreg_unique_geneIDs <- as.character(JN_upreg_unique_diff_only$SYMBOL)
JN_downreg_unique_geneIDs <- as.character(JN_downreg_unique_diff_only$SYMBOL)
JN_all_genes_unique_geneIDs <- as.character(JN_all_genes_unique$SYMBOL)


#ORA on upreg genes
JN_upreg_unique_geneIDs_ego <- enrichGO(gene = JN_upreg_unique_geneIDs, #downregulated genes (i know name is weird but just notation)
                                       universe = JN_all_genes_unique_geneIDs, #all genes, similar notation stuff
                                       keyType = "SYMBOL",
                                       OrgDb = org.Hs.eg.db, 
                                       ont = "BP", #try MF and BP
                                       pAdjustMethod = "BH", 
                                       qvalueCutoff = 0.05, 
                                       readable = TRUE)

JN_upreg_unique_geneIDs_ego_summary_BP <- data.frame(JN_upreg_unique_geneIDs_ego) #just converting to df
dotplot(JN_upreg_unique_geneIDs_ego, showCategory=30, title="JN Upregulated ATAC after KD (BP)")



#ORA on downreg genes
JN_downreg_unique_geneIDs_ego <- enrichGO(gene = JN_downreg_unique_geneIDs, #downregulated genes (i know name is weird but just notation)
                                        universe = JN_all_genes_unique_geneIDs, #all genes, similar notation stuff
                                        keyType = "SYMBOL",
                                        OrgDb = org.Hs.eg.db, 
                                        ont = "BP", #try MF and BP
                                        pAdjustMethod = "BH", 
                                        qvalueCutoff = 0.05, 
                                        readable = TRUE)

JN_downreg_unique_geneIDs_ego_summary_BP <- data.frame(JN_downreg_unique_geneIDs_ego) #just converting to df
dotplot(JN_downreg_unique_geneIDs_ego, showCategory=30, title="JN downregulated ATAC after KD (BP)")



############################
#BER
setwd("D:/1_ATAC Files/ATAC 2022 Files")
#load DSRCT consensus peaks
BER_DEsamples <- read.csv("BER_DSRCTshWT1_Differential_Peak_Binding_NP.csv")

#Create sets with differentially expressed peaks
BER_upreg <- BER_DEsamples[BER_DEsamples$FDR<0.05 & BER_DEsamples$Fold > 1 & BER_DEsamples$p.value <0.05, c(2,3,4,13,10,11,12)]
BER_downreg <- BER_DEsamples[BER_DEsamples$FDR<0.05 & BER_DEsamples$Fold < -1 & BER_DEsamples$p.value <0.05, c(2,3,4,13,10,11,12)]
BER_all_genes <- BER_DEsamples[,c(2,3,4,13,10,11,12)] #for background

write.table(BER_upreg, "BER_upreg.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
BER_upreg <- annotatePeak("BER_upreg.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0
BER_upreg <- as.data.frame(BER_upreg)

write.table(BER_downreg, "BER_downreg.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
BER_downreg <- annotatePeak("BER_downreg.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0
BER_downreg <- as.data.frame(BER_downreg)

write.table(BER_all_genes, "BER_all_genes.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
BER_all_genes <- annotatePeak("BER_all_genes.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0
BER_all_genes <- as.data.frame(BER_all_genes)

#BER RNA
setwd("D:/2_RNA Seq_Files/Excel Sheets")
fold_change <- read.xlsx("Regulation_JN_BER_SK2_BOD_KTS_ALL_GENES.xlsx")
BER_fold_change <- fold_change[,c(1,6)]
BER_fold_change <- BER_fold_change %>% dplyr::rename(Gene = gene_name, log2_fold_change =BERshWT1_log2FC) %>%
  na.omit()


#unique genes
BER_upreg_unique <- BER_upreg %>% 
  distinct(SYMBOL, .keep_all = TRUE)

BER_upreg_unique <- merge(BER_upreg_unique, BER_fold_change, by.x="SYMBOL", by.y="Gene")
BER_upreg_unique <- BER_upreg_unique %>% filter(log2_fold_change >1)

BER_downreg_unique <- BER_downreg %>% 
  distinct(SYMBOL, .keep_all = TRUE)

BER_downreg_unique <- merge(BER_downreg_unique, BER_fold_change, by.x="SYMBOL", by.y="Gene")
BER_downreg_unique <- BER_downreg_unique %>% filter(log2_fold_change < -1)


BER_all_genes_unique <- BER_all_genes %>% 
  distinct(SYMBOL, .keep_all = TRUE)


##getting unique stuff not in both groups
BER_upreg_unique_diff <- setdiff(BER_upreg_unique$SYMBOL, BER_downreg_unique$SYMBOL)

BER_upreg_unique_diff_only <- BER_upreg_unique[BER_upreg_unique$SYMBOL %in% BER_upreg_unique_diff, ]


BER_downreg_unique_diff <- setdiff(BER_downreg_unique$SYMBOL, BER_upreg_unique$SYMBOL)

BER_downreg_unique_diff_only <- BER_downreg_unique[BER_downreg_unique$SYMBOL %in% BER_downreg_unique_diff, ]

#pulling IDs/gene names only for analysis
BER_upreg_unique_geneIDs <- as.character(BER_upreg_unique_diff_only$SYMBOL)
BER_downreg_unique_geneIDs <- as.character(BER_downreg_unique_diff_only$SYMBOL)
BER_all_genes_unique_geneIDs <- as.character(BER_all_genes_unique$SYMBOL)





#ORA on upreg genes
BER_upreg_unique_geneIDs_ego <- enrichGO(gene = BER_upreg_unique_geneIDs, #downregulated genes (i know name is weird but just notation)
                                        universe = BER_all_genes_unique_geneIDs, #all genes, similar notation stuff
                                        keyType = "SYMBOL",
                                        OrgDb = org.Hs.eg.db, 
                                        ont = "BP", #try MF and BP
                                        pAdjustMethod = "BH", 
                                        qvalueCutoff = 0.05, 
                                        readable = TRUE)

BER_upreg_unique_geneIDs_ego_summary_BP <- data.frame(BER_upreg_unique_geneIDs_ego) #just converting to df
dotplot(BER_upreg_unique_geneIDs_ego, showCategory=30, title="BER Upregulated ATAC after KD (BP)")



#ORA on downreg genes
BER_downreg_unique_geneIDs_ego <- enrichGO(gene = BER_downreg_unique_geneIDs, #downregulated genes (i know name is weird but just notation)
                                          universe = BER_all_genes_unique_geneIDs, #all genes, similar notation stuff
                                          keyType = "SYMBOL",
                                          OrgDb = org.Hs.eg.db, 
                                          ont = "BP", #try MF and BP
                                          pAdjustMethod = "BH", 
                                          qvalueCutoff = 0.05, 
                                          readable = TRUE)

BER_downreg_unique_geneIDs_ego_summary_BP <- data.frame(BER_downreg_unique_geneIDs_ego) #just converting to df
dotplot(BER_downreg_unique_geneIDs_ego, showCategory=30, title="BER downregulated ATAC after KD (BP)")



#
JNBER_summary <- merge(JN_upreg_unique_geneIDs_ego_summary_BP, BER_upreg_unique_geneIDs_ego_summary_BP, by="Description")


##looking for consensus stuff
JNBER_upreg <- merge(JN_upreg, BER_upreg, by="SYMBOL")
JNBER_downreg <- merge(JN_downreg, BER_downreg, by="SYMBOL")

JNBER_upreg <- JNBER_upreg %>% 
  distinct(SYMBOL, .keep_all = TRUE)

JNBER_downreg <- JNBER_downreg %>% 
  distinct(SYMBOL, .keep_all = TRUE)



##overlap of ATAC up and down
BER_upreg_unique_diff_only <- na.omit(BER_upreg_unique_diff_only)
BER_downreg_unique_diff_only <- na.omit(BER_downreg_unique_diff_only)
JN_upreg_unique_diff_only <- na.omit(JN_upreg_unique_diff_only)
JN_downreg_unique_diff_only <- na.omit(JN_downreg_unique_diff_only)


my_colors <- brewer.pal(4, "Set1")

BER_upreg_unique_naomit <- na.omit(BER_upreg_unique)
BER_downreg_unique_naomit <- na.omit(BER_downreg_unique)
JN_upreg_unique_naomit <- na.omit(JN_upreg_unique)
JN_downreg_unique_naomit <- na.omit(JN_downreg_unique)

venn.diagram(
  x = list(BER_upreg_unique_naomit$SYMBOL, BER_downreg_unique_naomit$SYMBOL, JN_upreg_unique_naomit$SYMBOL, JN_downreg_unique_naomit$SYMBOL),
  category.names = c("BER-DSRCT Increased", "BER-DSRCT Decreased", "JN-DSRCT Increased", "JN-DSRCT Decreased"),
  filename = "JN_BER_accessibility_overlap2.png",
  output = FALSE ,
  imagetype="png",
  main = "JN-BER Accessibility",
  main.cex = 1.5,
  cex = 1.5,
  cat.cex = 1,
  cat.pos = c(-10, 10, 0, 0),
  fill = my_colors
)


