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






JN_DAR<- JN_DEpeaks_filtered
BER_DAR<- BER_DEpeaks_filtered

JN_DAR_FC=JN_DAR$FC
names(JN_DAR_FC)<-as.character(JN_DAR$geneId)
JN_DAR_FC <- sort(JN_DAR_FC, decreasing = TRUE)

JN_DAR_gsea <- gseGO(geneList = JN_DAR_FC, 
                         ont = 'BP', 
                         OrgDb = org.Hs.eg.db, 
                         keyType = "ENTREZID",
                         #nPerm = 1000, 
                         minGSSize = 20, 
                         pvalueCutoff = 1,
                         verbose = FALSE)
dotplot(JN_DAR_gsea, showCategory=15) + ggtitle("JN DAR")
JN_DAR_gsea_results <- JN_DAR_gsea@result
write.xlsx(JN_DAR_gsea_results, "JN_DAR_gsea_results.xlsx", quote=F)

JN_DAR_DEG_gseDO<-gseDO(JN_DAR_DEG_FC,
                        pvalueCutoff = 1)
dotplot(JN_DAR_DEG_gseDO, showCategory=15) + ggtitle("JN: ATAC/EnhancerMarks/shWT1 RNA_gseDO")
JN_DAR_DEG_gseDO_result <- JN_DAR_DEG_gseDO@result
write.xlsx(JN_gseaDO_enh_result, "JN_gseaDO_enh_result.xlsx", quote=F)




BER_DAR_FC=BER_DAR$FC
names(BER_DAR_FC)<-as.character(BER_DAR$geneId)
BER_DAR_FC <- sort(BER_DAR_FC, decreasing = TRUE)

BER_DAR_gsea <- gseGO(geneList = BER_DAR_FC, 
                     ont = 'BP', 
                     OrgDb = org.Hs.eg.db, 
                     keyType = "ENTREZID",
                     #nPerm = 1000, 
                     minGSSize = 20, 
                     pvalueCutoff = .05,
                     verbose = FALSE)
dotplot(BER_DAR_gsea, showCategory=15) + ggtitle("BER DAR")
BER_DAR_gsea_results <- BER_DAR_gsea@result
write.xlsx(BER_DAR_gsea_results, "BER_DAR_gsea_results.xlsx", quote=F)

JN_DAR_DEG_gseDO<-gseDO(JN_DAR_DEG_FC,
                        pvalueCutoff = 1)
dotplot(JN_DAR_DEG_gseDO, showCategory=15) + ggtitle("JN: ATAC/EnhancerMarks/shWT1 RNA_gseDO")
JN_DAR_DEG_gseDO_result <- JN_DAR_DEG_gseDO@result
write.xlsx(JN_gseaDO_enh_result, "JN_gseaDO_enh_result.xlsx", quote=F)




