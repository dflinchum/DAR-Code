library(RColorBrewer)
library(tidygraph)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(readxl)
library(ggplot2)
library(dplyr)




#BER
#need to annotate for the gene 
setwd("D:/1_ATAC Files/ATAC 2022 Files")
BER_DEpeaksAnno <- annotatePeak("BER_NP_DEpeaks.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0

BER_DEpeaksAnno<-as.data.frame(BER_DEpeaksAnno) #as data frame to manipulate easier
BER_DEpeaksAnno <- BER_DEpeaksAnno %>%  #rename columns! 
  rename("FC" = "V5")

#BER_DEpeaks_DAR <- BER_DEpeaksAnno %>% filter(FC > 1 | FC < -1)

BER_DEpeaksAnno <- BER_DEpeaksAnno %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))


BER_stable <- BER_DEpeaksAnno %>% filter(FC < 1 & FC > -1)
BER_closing <- BER_DEpeaksAnno %>% filter(FC < -1)
BER_opening <- BER_DEpeaksAnno %>% filter(FC > 1)


BER_stable_count <- BER_stable %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

BER_closing_count <- BER_closing %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

BER_opening_count <- BER_opening %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

BER_opening_closing <- merge(BER_closing_count, BER_opening_count, by="annotation")
BER_open_close_stable <- merge(BER_stable_count, BER_opening_closing, by="annotation")

BER_open_close_stable <- BER_open_close_stable[BER_open_close_stable$annotation != "5' UTR", ] #remove 5' and downstream rows
BER_open_close_stable <- BER_open_close_stable[BER_open_close_stable$annotation != "Downstream (<=300bp)", ] #remove 5' and downstream rows


BER_merged_anno_counts <- read_xlsx("BER_NegPosStable_Anno.xlsx", sheet = "Sheet2")


my_colors <- brewer.pal(5, "Set1")


ggplot(data = BER_merged_anno_counts, aes(x = Status, y = Count, fill = Annotation)) +
  geom_col(position = "fill") +
  labs(title = "BER", 
       x = "Peak Status", 
       y = "Percentage of Peaks") +
  scale_x_discrete(limit = unique(BER_merged_anno_counts$Status)) +
  theme_bw() +
  #scale_fill_manual(values = my_colors) +
  theme(
    # Legend title text size
    legend.title = element_text(size = 28),       # Adjust size as needed
    
    # Legend text size (labels)
    legend.text = element_text(size = 26),       # Adjust size as needed
    
    # X-axis label text size (e.g., "Status")
    axis.title.x = element_text(size = 28),       # Adjust size as needed
    
    # Y-axis label text size (e.g., "Percentage of Peaks")
    axis.title.y = element_text(size = 28),       # Adjust size as needed
    
    # X-axis tick text size (e.g., the Status categories)
    axis.text.x = element_text(size = 26),        # Adjust size as needed
    
    # Y-axis tick text size (e.g., the percentage values)
    axis.text.y = element_text(size = 26)         # Adjust size as needed
  )


#JN
#need to annotate for the gene 
setwd("D:/1_ATAC Files/ATAC 2022 Files")
JN_DEpeaksAnno <- annotatePeak("JN_NP_DEpeaks.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0

JN_DEpeaksAnno<-as.data.frame(JN_DEpeaksAnno) #as data frame to manipulate easier
JN_DEpeaksAnno <- JN_DEpeaksAnno %>%  #rename columns! 
  rename("FC" = "V5")

#JN_DEpeaks_DAR <- JN_DEpeaksAnno %>% filter(FC > 1 | FC < -1)

JN_DEpeaksAnno <- JN_DEpeaksAnno %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))


JN_stable <- JN_DEpeaksAnno %>% filter(FC < 1 & FC > -1)
JN_closing <- JN_DEpeaksAnno %>% filter(FC < -1)
JN_opening <- JN_DEpeaksAnno %>% filter(FC > 1)


JN_stable_count <- JN_stable %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

JN_closing_count <- JN_closing %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

JN_opening_count <- JN_opening %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

JN_opening_closing <- merge(JN_closing_count, JN_opening_count, by="annotation")
JN_open_close_stable <- merge(JN_stable_count, JN_opening_closing, by="annotation")

JN_open_close_stable <- JN_open_close_stable[JN_open_close_stable$annotation != "5' UTR", ] #remove 5' and downstream rows
JN_open_close_stable <- JN_open_close_stable[JN_open_close_stable$annotation != "Downstream (<=300bp)", ] #remove 5' and downstream rows


JN_merged_anno_counts <- read_xlsx("JN_NegPosStable_Anno.xlsx", sheet = "Sheet2")


my_colors <- brewer.pal(5, "Set1")


ggplot(data = JN_merged_anno_counts, aes(x = Status, y = Count, fill = Annotation)) +
  geom_col(position = "fill") +
  labs(title = "JN", 
       x = "Peak Status", 
       y = "Percentage of Peaks") +
  scale_x_discrete(limit = unique(JN_merged_anno_counts$Status)) +
  theme_bw() +
  #scale_fill_manual(values = my_colors) +
  theme(
    # Legend title text size
    legend.title = element_text(size = 28),       # Adjust size as needed
    
    # Legend text size (labels)
    legend.text = element_text(size = 26),       # Adjust size as needed
    
    # X-axis label text size (e.g., "Status")
    axis.title.x = element_text(size = 28),       # Adjust size as needed
    
    # Y-axis label text size (e.g., "Percentage of Peaks")
    axis.title.y = element_text(size = 28),       # Adjust size as needed
    
    # X-axis tick text size (e.g., the Status categories)
    axis.text.x = element_text(size = 26),        # Adjust size as needed
    
    # Y-axis tick text size (e.g., the percentage values)
    axis.text.y = element_text(size = 26)         # Adjust size as needed
  )
