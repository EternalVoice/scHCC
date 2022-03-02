# >>>> 对于up/down的基因，结合细胞信息分析其最可能的细胞来源(iCAF, mCAF, stromal cell, 上皮细胞，肿瘤细胞，其他细胞等)

rm(list = ls())

library(Seurat)
library(ggplot2)
library(patchwork)

tumor <- readRDS("rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR.rds")

# liver stromal cell marker:VCAM1
p1 <- VlnPlot(tumor,features = "VCAM1",group.by = "CellType",pt.size = 0) + NoLegend()

current_ids <- c("B cell","CAF","Malignant cell","T cell","TAM","TEC","unclassified")
new_ids <- c("B cell","CAF","Malignant cell","T cell","TAM","TEC","Stromal cell")
tumor$CellType <- plyr::mapvalues(tumor$CellType, from = current_ids, to = new_ids)
current_ids <- c("B cell","CAF","CD4+ T cell","CD8+ T cell","Malignant cell","TAM","TEC","unclassified")
new_ids <- c("B cell","CAF","CD4+ T cell","CD8+ T cell","Malignant cell","TAM","TEC","Stromal cell")
tumor$CellType2 <- plyr::mapvalues(tumor$CellType2, from = current_ids, to = new_ids)

saveRDS(tumor,file = "rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR_2.rds")


##>>>>>>>>>>>>>>>>>>> UP 组 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

## CD3 UP
CD3.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",sheet = 1)
colnames(CD3.up)[1] <- 'Genes'
CD3.up.genes <- CD3.up$Genes

pal <- c("#0066A5","#46A040","#00AF99","#875691","#98D9E9","#F6313E","#FFA300")
p1 <- DimPlot(tumor, reduction = "tsne", group.by = "CellType", cols = pal)
for(gene in CD3.up.genes){
  p2 <- FeaturePlot(tumor, reduction = "tsne", features = gene, pt.size = 0.9, cols = c("#F3F3F4","#BF0032"))
  plt <- p1 / p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/UP/CD3T/FeatureScatter/",gene,".pdf"), plt, device = "pdf", width = 6, height = 8)
}

## CD4 UP
CD4.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD4.high_vs_low.DEmarkers.xlsx",sheet = 1)
colnames(CD4.up)[1] <- 'Genes'
CD4.up.genes <- CD4.up$Genes

pal <- c("#0066A5","#46A040","#00AF99","#875691","#98D9E9","#F6313E","#FFA300")
p1 <- DimPlot(tumor, reduction = "tsne", group.by = "CellType", cols = pal)
for(gene in CD4.up.genes){
  p2 <- FeaturePlot(tumor, reduction = "tsne", features = gene, pt.size = 0.9, cols = c("#F3F3F4","#BF0032"))
  plt <- p1 / p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/UP/CD4T/FeatureScatter/",gene,".pdf"), plt, device = "pdf", width = 6, height = 8)
}

## CD8 UP
CD8.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD8.high_vs_low.DEmarkers.xlsx",sheet = 1)
colnames(CD8.up)[1] <- 'Genes'
CD8.up.genes <- CD8.up$Genes

pal <- c("#0066A5","#46A040","#00AF99","#875691","#98D9E9","#F6313E","#FFA300")
p1 <- DimPlot(tumor, reduction = "tsne", group.by = "CellType", cols = pal)
for(gene in CD8.up.genes){
  p2 <- FeaturePlot(tumor, reduction = "tsne", features = gene, pt.size = 0.9, cols = c("#F3F3F4","#BF0032"))
  plt <- p1 / p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/UP/CD8T/FeatureScatter/",gene,".pdf"), plt, device = "pdf", width = 6, height = 8)
}


##>>>>>>>>>>>>>>>>>> DOWN组 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

## CD3 Down
CD3.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",sheet = 2)
colnames(CD3.down)[1] <- 'Genes'
CD3.down.genes <- CD3.down$Genes

pal <- c("#0066A5","#46A040","#00AF99","#875691","#98D9E9","#F6313E","#FFA300")
p1 <- DimPlot(tumor, reduction = "tsne", group.by = "CellType", cols = pal)
for(item in CD3.down.genes){
  p2 <- FeaturePlot(tumor,reduction = "tsne",features = item,pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  plt <- p1 / p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/Down/CD3T/FeatureScatter/",item,".pdf"),plt,device = "pdf",width = 6,height = 8)
}

## CD4 Down
CD4.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",sheet = 2)
colnames(CD4.down)[1] <- 'Genes'
CD4.down.genes <- CD4.down$Genes

pal <- c("#0066A5","#46A040","#00AF99","#875691","#98D9E9","#F6313E","#FFA300")
p1 <- DimPlot(tumor, reduction = "tsne", group.by = "CellType", cols = pal)
for(item in CD4.down.genes){
  p2 <- FeaturePlot(tumor,reduction = "tsne",features = item,pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  plt <- p1 / p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/Down/CD4T/FeatureScatter/",item,".pdf"),plt,device = "pdf",width = 6,height = 8)
}

## CD8 Down
CD8.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD8.high_vs_low.DEmarkers.xlsx",sheet = 2)
colnames(CD8.down)[1] <- 'Genes'
CD8.down.genes <- CD8.down$Genes

pal <- c("#0066A5","#46A040","#00AF99","#875691","#98D9E9","#F6313E","#FFA300")
p1 <- DimPlot(tumor, reduction = "tsne", group.by = "CellType", cols = pal)
for(item in CD8.down.genes){
  p2 <- FeaturePlot(tumor,reduction = "tsne",features = item,pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  plt <- p1 / p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/Down/CD8T/FeatureScatter/",item,".pdf"),plt,device = "pdf",width = 6,height = 8)
}
