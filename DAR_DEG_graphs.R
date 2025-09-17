library(tidyverse)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
install.packages("ggthemes") # Install 
library(ggthemes)


setwd("D:/5_R_generated_Data/JN/3. 2022 ATAC")
#JN All Peaks  Bar Graphs
JN_allpeaks <- annotatePeak("JN_NP_allpeaks.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0
JN_allpeaks<-as.data.frame(JN_allpeaks)
allpeaks<- JN_allpeaks %>% dplyr::select(annotation)
allpeaks<-as_tibble(allpeaks)
allpeaks <- allpeaks %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))
allpeaks_count <- allpeaks %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

allpeaks_count <- allpeaks_count %>% 
  rename("Annotation" = "annotation",
         "Peaks" = "count")


data1 <- allpeaks_count                                                 # Replicate original data
data1$Annotation <- factor(data1$Annotation,                                    # Factor levels in decreasing order
                           levels = data1$Annotation[order(data1$Peaks, decreasing = TRUE)])

test1 <- ggplot(data=data1, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen") +
  theme(text = element_text(size = 25),
        panel.background = element_rect(fill = "white"))

test1 + ggtitle("JN +/-dox All Peaks Distribution")




#JN Differentially Accessible Peaks Bar Graphs
JN_DE_peaks <- annotatePeak("JN_NP_DEpeaks.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0
JN_DE_peaks<-as.data.frame(JN_DE_peaks)
DEpeaks<- JN_DE_peaks %>% dplyr::select(annotation)
DEpeaks<-as_tibble(DEpeaks)
DEpeaks <- DEpeaks %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))
DEpeaks_count <- DEpeaks %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

DEpeaks_count <- DEpeaks_count %>% 
  rename("Annotation" = "annotation",
         "Peaks" = "count")


data2 <- DEpeaks_count                                                 # Replicate original data
data2$Annotation <- factor(data2$Annotation,                                    # Factor levels in decreasing order
                           levels = data2$Annotation[order(data2$Peaks, decreasing = TRUE)])

test2 <- ggplot(data=data2, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen") +
  theme(text = element_text(size = 25),plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"))

test2 + ggtitle("JN shWT1 Differentially Accessible Peak Distribution")
        




##DO BER BAR Graphs
setwd("D:/5_R_generated_Data/BER/3. 2022 ATAC")
#BER All Peaks 
BER_allpeaks <- annotatePeak("BER_NP_allpeaks.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0
BER_allpeaks<-as.data.frame(BER_allpeaks)
BER_allpeaks<- BER_allpeaks %>% dplyr::select(annotation)
BER_allpeaks<-as_tibble(BER_allpeaks)
BER_allpeaks <- BER_allpeaks %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))
BER_allpeaks_count <- BER_allpeaks %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

BER_allpeaks_count <- BER_allpeaks_count %>% 
  rename("annotation" = "Annotation",
         "count" = "Peaks")


BER_data1 <- BER_allpeaks_count                                                 # Replicate original data
BER_data1$Annotation <- factor(BER_data1$Annotation,                                    # Factor levels in decreasing order
                           levels = BER_data1$Annotation[order(BER_data1$Peaks, decreasing = TRUE)])

BER_test1 <- ggplot(data=BER_data1, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen") +
  theme(text = element_text(size = 25),
        panel.background = element_rect(fill = "white"))

BER_test1 + ggtitle("BER +/-dox All Peaks Distribution")



#BER Differentially Accessible Peaks Bar Graphs
BER_DE_peaks <- annotatePeak("BER_NP_DEpeaks.bed", TxDb=txdb, tssRegion=c(-3000,3000), annoDb="org.Hs.eg.db") #these have FDR <.05 and FC>/< 0
BER_DE_peaks<-as.data.frame(BER_DE_peaks)
BER_DEpeaks<- BER_DE_peaks %>% dplyr::select(annotation)
BER_DEpeaks<-as_tibble(BER_DEpeaks)
BER_DEpeaks <- BER_DEpeaks %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))
BER_DEpeaks_count <- BER_DEpeaks %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

BER_DEpeaks_count <- BER_DEpeaks_count %>% 
  rename("Annotation" = "annotation",
         "Peaks" = "count")


BER_data2 <- BER_DEpeaks_count                                                 # Replicate original data
BER_data2$Annotation <- factor(BER_data2$Annotation,                                    # Factor levels in decreasing order
                           levels = BER_data2$Annotation[order(BER_data2$Peaks, decreasing = TRUE)])

BER_test2 <- ggplot(data=BER_data2, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen") +
  theme(text = element_text(size = 25), plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"))

BER_test2 + ggtitle("BER shWT1 Differentially Accessible Peak Distribution")






#for log2 fold change graph with open/less open
#try making two different columns with the open vs closed
#and have FC be on the y axis 

JN_FCpeaks<- JN_DE_peaks %>% dplyr::select(V5)
JN_FCpeaks <- JN_FCpeaks %>% 
  rename("V5" = "Log2 Fold Change")






##ALLLLLLL FROM PRACTICE/WORKING STUFF OUT
allpeaks<-as_tibble(JN_allpeaks)
allpeaks <- allpeaks %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    TRUE ~ annotation
  ))

allpeaks_count <- allpeaks %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

#example
data1 <- data                                                 # Replicate original data
data1$x <- factor(data1$x,                                    # Change ordering manually
                  levels = c("B", "D", "E", "C", "A"))
#real
data1 <- allpeaks_count                                                 # Replicate original data
data1$Annotation <- factor(data1$Annotation,                            # Change ordering manually
                           levels = c("Distal Intergenic", "Intron", "Exon", "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "3' UTR", "5' UTR", "Downstream (<=300bp)"))

ggplot(data=data1, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen")




data2 <- data1                                                 # Replicate original data
data2$Annotation <- factor(data2$Annotation,                                    # Factor levels in increasing order
                           levels = data2$Annotation[order(data2$Peaks)])
ggplot(data=data2, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen")


data3 <- data1                                                 # Replicate original data
data3$Annotation <- factor(data3$Annotation,                                    # Factor levels in decreasing order
                           levels = data3$Annotation[order(data3$Peaks, decreasing = TRUE)])
ggplot(data=data3, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen")