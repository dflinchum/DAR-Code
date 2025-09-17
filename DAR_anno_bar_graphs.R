


JN_DEpeaks_DAR_all_up <- JN_DEpeaksAnno %>% filter(FC < 0)  #filter by abs value; works great, nothing between -.99 to .99!  
JN_DEpeaks_DAR_all_down <- JN_DEpeaksAnno %>% filter(FC > 0)  #filter by abs value; works great, nothing between -.99 to .99!


BER_DEpeaks_DAR_all_up <- BER_DEpeaksAnno %>% filter(FC < 0)  
BER_DEpeaks_DAR_all_down <- BER_DEpeaksAnno %>% filter(FC > 0)  


##this is basically all "differentially accessible regions" or DE peaks because there's no fold change filter here
#broken into "up" and "down" expression 
#jn up
JN_DEpeaks_DAR_all_up <- JN_DEpeaks_DAR_all_up %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))
JN_DEpeaks_DAR_all_up_count <- JN_DEpeaks_DAR_all_up %>%  
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

JN_DEpeaks_DAR_all_up_count <- JN_DEpeaks_DAR_all_up_count %>% 
  rename("Annotation" = "annotation",
         "Peaks" = "count")


otter1 <- JN_DEpeaks_DAR_all_up_count                                                 # Replicate original data
otter1$Annotation <- factor(otter1$Annotation,                                    # Factor levels in decreasing order
                           levels = otter1$Annotation[order(otter1$Peaks, decreasing = TRUE)])

bear1 <- ggplot(data=otter1, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen") +
  theme(text = element_text(size = 25),plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"))

bear1 + ggtitle("JN atac all up")


#jn down
JN_DEpeaks_DAR_all_down <- JN_DEpeaks_DAR_all_down %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))
JN_DEpeaks_DAR_all_down_count <- JN_DEpeaks_DAR_all_down %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

JN_DEpeaks_DAR_all_down_count <- JN_DEpeaks_DAR_all_down_count %>% 
  rename("Annotation" = "annotation",
         "Peaks" = "count")


otter2 <- JN_DEpeaks_DAR_all_down_count                                                 # Replicate original data
otter2$Annotation <- factor(otter2$Annotation,                                    # Factor levels in decreasing order
                            levels = otter2$Annotation[order(otter2$Peaks, decreasing = TRUE)])

bear2 <- ggplot(data=otter2, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen") +
  theme(text = element_text(size = 25),plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"))

bear2 + ggtitle("JN atac all down")


#ber up
BER_DEpeaks_DAR_all_up <- BER_DEpeaks_DAR_all_up %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))
BER_DEpeaks_DAR_all_up_count <- BER_DEpeaks_DAR_all_up %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

BER_DEpeaks_DAR_all_up_count <- BER_DEpeaks_DAR_all_up_count %>% 
  rename("Annotation" = "annotation",
         "Peaks" = "count")


otter3 <- BER_DEpeaks_DAR_all_up_count                                                 # Replicate original data
otter3$Annotation <- factor(otter3$Annotation,                                    # Factor levels in decreasing order
                            levels = otter3$Annotation[order(otter3$Peaks, decreasing = TRUE)])

bear3 <- ggplot(data=otter3, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen") +
  theme(text = element_text(size = 25),plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"))

bear3 + ggtitle("BER atac all up")

#ber down
BER_DEpeaks_DAR_all_down <- BER_DEpeaks_DAR_all_down %>%  #THIS WORKED
  mutate(annotation = case_when(
    str_detect(annotation, "Exon") ~ "Exon",
    str_detect(annotation, "Intron")  ~ "Intron",
    str_detect(annotation, "Promoter") ~ "Promoter",
    TRUE ~ annotation
  ))
BER_DEpeaks_DAR_all_down_count <- BER_DEpeaks_DAR_all_down %>%  ##WORKED ITS PERFECT; RENAME PROMOTER IF YOU WANT
  group_by(annotation) %>%
  summarize(
    count = n(),
  )

BER_DEpeaks_DAR_all_down_count <- BER_DEpeaks_DAR_all_down_count %>% 
  rename("Annotation" = "annotation",
         "Peaks" = "count")


otter4 <- BER_DEpeaks_DAR_all_down_count                                                 # Replicate original data
otter4$Annotation <- factor(otter4$Annotation,                                    # Factor levels in decreasing order
                            levels = otter4$Annotation[order(otter4$Peaks, decreasing = TRUE)])

bear4 <- ggplot(data=otter4, aes(x=Annotation, y=Peaks)) +
  geom_bar(stat="identity", fill = "darkgreen") +
  theme(text = element_text(size = 25),plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white"))

bear4 + ggtitle("BER atac all down")