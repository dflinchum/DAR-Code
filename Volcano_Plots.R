library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(openxlsx)

df <- read.xlsx("BER_diff_binding_4volcano.xlsx")
#colnames(df)[1] = "ID" ##renaming column 1 to ID

##this worked; created a basic scatterplot with black dots; good start 
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(0.8), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(0.8), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))
ggplot(data = df, aes(x = Fold, y = -log10(FDR))) +
  geom_point() 

##adds column with "no" in all the cells down
df$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$Fold > 0.6 & df$p.value < 0.05] <- "DOWN"


# if log2Foldchange < -0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$Fold < -0.6 & df$p.value < 0.05] <- "UP"

## shows all the different color names to choose for the graph
colors()

ggplot(data = df, aes(x = Fold, y = -log10(p.value), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "red2"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Downregulated", "Not significant", "Upregulated")) # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
##example below!
df$delabel <- ifelse(df$gene_symbol %in% head(df[order(df$padj), "gene_symbol"], 30), df$gene_symbol, NA)

##messing around with my stuff!
df$delabel <- ifelse(df$gene_symbol %in% head(df[order(df$FDR), "gene_symbol"], 20), df$gene_symbol, NA)

##trying to name individual ones
df$delabel <- ifelse(df$gene_symbol == "FOXA2", df$gene_symbol, NA)


##was using this for a while
ggplot(data = df, aes(x = Fold, y = -log10(p.value), col = diffexpressed, 
                      label = ifelse(abs(Fold) > 1, delabel, ""))) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "red2"), # to set the colours of our variable<br /><br /><br />
                     labels = c("EWSR1-WT1 'Off'", "Not significant", "EWSR1-WT1 'On'"))+ # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  #coord_cartesian(ylim = c(0, 35), xlim = c(-4, 4)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Expression:', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  ggtitle('BER Differentially Accessible Peaks') +
  geom_text_repel(max.overlaps = Inf)



##for plot below, mostly playing with different sig variables or names
ggplot(data = df, aes(x = Fold, y = -log10(FDR), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "black", "red2"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Downregulated", "Not significant", "Upregulated"))+ # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 75), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Expression:', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR")) +
  ggtitle('BER Differentially Accessible Peaks') +
  geom_text_repel(max.overlaps = Inf)



##
##JN
##


JN <- read.xlsx("JN_diff_bind_4volcano.xlsx")
#colnames(JN)[1] = "ID" ##renaming column 1 to ID

##this worked; created a basic scatterplot with black dots; good start 
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(0.8), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(0.8), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))
ggplot(data = JN, aes(x = Fold, y = -log10(FDR))) +
  geom_point() 

##adds column with "no" in all the cells down
JN$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "DOWN"
JN$diffexpressed[JN$Fold > 0.6 & JN$p.value < 0.05] <- "DOWN"


# if log2Foldchange < -0.6 and pvalue < 0.05, set as "UP"
JN$diffexpressed[JN$Fold < -0.6 & JN$p.value < 0.05] <- "UP"

ggplot(data = JN, aes(x = Fold, y = -log10(p.value), col = diffexpressed)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "red2"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Downregulated", "Not significant", "Upregulated")) # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
##example below!
JN$delabel <- ifelse(JN$gene_symbol %in% head(JN[order(JN$padj), "gene_symbol"], 30), JN$gene_symbol, NA)

##messing around with my stuff!
JN$delabel <- ifelse(JN$gene_symbol %in% head(JN[order(JN$FDR), "gene_symbol"], 30), JN$gene_symbol, NA)


ggplot(data = JN, aes(x = Fold, y = -log10(p.value), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "red2"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Downregulated", "Not significant", "Upregulated"))+ # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  #coord_cartesian(ylim = c(0, 250), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Expression:', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  ggtitle('JN Differentially Accessible Peaks') +
  geom_text_repel(max.overlaps = Inf)


##with the one below, I was trying to only label FC>1 so that the gray wasn't labeled. Seems to work!
ggplot(data = JN, aes(x = Fold, y = -log10(p.value), col = diffexpressed, 
                      label = ifelse(abs(Fold) > 1, delabel, ""))) + # Label only if |Fold| > 1
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "grey", "red2"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  labs(color = 'Expression:', x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
  ggtitle('JN Differentially Accessible Peaks') +
  geom_text_repel(max.overlaps = Inf)  






ggplot(data = JN, aes(x = Fold, y = -log10(FDR), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "gray", "red2"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))+ # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  #coord_cartesian(ylim = c(0, 30), xlim = c(-7, 7)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Expression:', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR")) +
  ggtitle('JN Differentially Accessible Peaks') +
  geom_text_repel(max.overlaps = Inf)

#10.31.23 note
##FDR vs pvalue doesn't seem super big difference, play around with both and pick what looks better
##and finalize the axis lengths or whatever

##may try swapping +/- stuff just for aesthetic purposes; put the upregulated on the right and downregulated on left

##11.2.23 note
#plots look good, happy with how they look and labeling (easy to edit now too)
#only thing to play around with is the pvalue/padj/FDR (can probably just use FDR)
#AND the fold change. Current setup has the summarized, mean FC from the differential stuff. Play around with highest sig FC instead of mean or whateverr
# or just use multiple and say hey that's how ATAC works? 


##split into "open/close" or +/-?? and then sum up?
##or do it by within a certain distance of TSS and sum up CLOSE ones? this seems like the most effective way
##or do all?? It works with all/multiple; only problem is all of them show up on graph labeled
##there's probably a way to get so that only ONE label shows up

# Biostatsquid theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))
