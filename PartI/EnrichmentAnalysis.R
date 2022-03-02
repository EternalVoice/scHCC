EnrichmentAnalysis <- function(DEmarkers,enrich.method,filePath,fileName){
  if(!requireNamespace("clusterProfiler")) BiocManager::install("clusterProfiler")
  if(!requireNamespace("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db")
  if(!requireNamespace("patchwork")) install.packages("patchwork")
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(patchwork)
  if(enrich.method == "GO"){
    # GO
    ego_ALL <- enrichGO(gene = rownames(DEmarkers), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    ego_CC <- enrichGO(gene = rownames(DEmarkers), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    ego_MF <- enrichGO(gene = rownames(DEmarkers), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    ego_BP <- enrichGO(gene = rownames(DEmarkers), OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
    ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
    ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)
    ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
    p_BP <- barplot(ego_BP, showCategory = 10) + ggtitle("Biological Process")
    p_CC <- barplot(ego_CC, showCategory = 10) + ggtitle("Cellular Component")
    p_MF <- barplot(ego_MF, showCategory = 10) + ggtitle("Molecular Function")
    pltc <- p_BP/p_CC/p_MF
    ggsave(paste(filePath,paste0(fileName,".GO.pdf"),sep="/"),plot=pltc, device = "pdf", width=12,height=10)
    save(ego_ALL,ego_BP,ego_MF,ego_CC,file = paste(filePath,paste0(fileName,".RData"),sep="/"))
    ego_results <- list(ego_ALL = ego_ALL,ego_BP = ego_BP, ego_CC = ego_CC, ego_MF = ego_MF)
    return(ego_results)
  } else if(enrich.method == "KEGG"){
    # kegg
    genelist <- bitr(rownames(DEmarkers),fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    genelist <- pull(genelist, ENTREZID)
    ekegg <- enrichKEGG(gene = genelist, organism = "hsa")
    p <- dotplot(ekegg, showCategory = 20)
    ggsave(paste(filePath,paste0(fileName,".KEGG.pdf"),sep="/"),plot = p, device = "pdf", width = 12,height = 10)
    save(ekegg,file = paste(filePath,paste0(fileName,".RData"), sep = "/"))
    return(ekegg)
  }
}