BER_downregulated_peaks <- BER_DEpeaks %>% filter(Fold > 1)
BER_downregulated_peaks <- BER_downregulated_peaks %>% dplyr::select(seqnames, start,end, Peak_Number)
BER_downregulated_peaks <- BER_downregulated_peaks %>% dplyr::rename("chr"=seqnames)

export.bed(BER_downregulated_peaks, "BER_downregulated_peaks_forHomer.bed")




##making annotation change graphs

data <- data.frame(
  Annotation = c("Intron", "Distal Intergenic", "Promoter", "Exon", "5' UTR", "3' UTR", "Downstream"),
  BER_Percent_Change = c(15.86016, 11.05671, -28.0928, 0.67844, 0.178106, 0.308945, 0.010436),
  JN_Percent_Change = c(-0.44154, -0.75181, -0.04979, 1.327506, 0, -0.065, -0.01937)
)
library(tidyr)
data_long <- data %>%
  pivot_longer(cols = c(BER_Percent_Change, JN_Percent_Change), names_to = "Cell_Line", values_to = "Percent_Change")

library(ggplot2)



ggplot(data_long, aes(x = Annotation, y = Percent_Change, fill = Cell_Line)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Percent Change in Peak Annotation by Cell Line",
       x = "Annotation",
       y = "Percent Change") +
  theme_minimal() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", alpha=0.5) +
  scale_fill_manual(values = c("BER_Percent_Change" = "skyblue", "JN_Percent_Change" = "lightcoral"))  # Corrected values


ggplot(data_long, aes(x = Annotation, y = Percent_Change, fill = Cell_Line)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Percent Change in Peak Annotation by Cell Line",
       x = "Annotation",
       y = "Percent Change") +
  theme_minimal() +
  geom_hline(yintercept=0, linetype="dashed", color = "black", alpha=0.5) +
  scale_fill_manual(values = c("BER_Percent_Change" = "skyblue", "JN_Percent_Change" = "lightcoral"),
                    labels = c("BER_Percent_Change" = "BER-DSRCT", "JN_Percent_Change" = "JN-DSRCT"))  # Modified labels
