library(Seurat)
library(ggplot2)
library(patchwork)

scRNAlist_sub_with_PVTT_MLN <- readRDS("rds/with_PVTT_MLN/1.scRNAlist_sub_with_PVTT_MLN.rds")
scRNA_merge_with_PVTT_MLN <- merge(scRNAlist_sub_with_PVTT_MLN[[1]], y = c(scRNAlist_sub_with_PVTT_MLN[[2]], scRNAlist_sub_with_PVTT_MLN[[3]], scRNAlist_sub_with_PVTT_MLN[[4]]))
saveRDS(scRNA_merge_with_PVTT_MLN, file = "rds/with_PVTT_MLN/2.scRNA_merge_with_PVTT_MLN.rds")
# preprocessing
scRNA_merge_with_PVTT_MLN <- readRDS("rds/with_PVTT_MLN/2.scRNA_merge_with_PVTT_MLN.rds")
# MT ratio
scRNA_merge_with_PVTT_MLN[["percent.mt"]] <- PercentageFeatureSet(scRNA_merge_with_PVTT_MLN, pattern = "^MT-")
# ribo ratio
scRNA_merge_with_PVTT_MLN[["percent.ribo"]] <- PercentageFeatureSet(scRNA_merge_with_PVTT_MLN, pattern = "^RP[SL]")
# red cell
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_match <- match(HB.genes, rownames(scRNA_merge_with_PVTT_MLN))
HB.genes <- rownames(scRNA_merge_with_PVTT_MLN)[match(HB.genes, rownames(scRNA_merge_with_PVTT_MLN))]
HB.genes <- HB.genes[!is.na(HB.genes)]
scRNA_merge_with_PVTT_MLN[["percent.hb"]] <- PercentageFeatureSet(scRNA_merge_with_PVTT_MLN, features = HB.genes)

col.num <- length(levels(Idents(scRNA_merge_with_PVTT_MLN)))
p1 <- VlnPlot(scRNA_merge_with_PVTT_MLN, features = "nFeature_RNA",pt.size = 0.0,ncol = 1,cols = rainbow(col.num)) + NoLegend()
p2 <- VlnPlot(scRNA_merge_with_PVTT_MLN, features = "nCount_RNA",pt.size = 0.0,ncol = 1,cols = rainbow(col.num)) + NoLegend()
p3 <- VlnPlot(scRNA_merge_with_PVTT_MLN, features = "percent.mt",pt.size = 0.0,ncol = 1,cols = rainbow(col.num)) + NoLegend()
ggsave("result/with_PVTT_MLN/1.scRNA_merge_vlnplot_before_QC_geneCnt.pdf",p1,device = "pdf",width = 8,height = 6)
ggsave("result/with_PVTT_MLN/1.scRNA_merge_vlnplot_before_QC_UMI.pdf",p2,device = "pdf",width = 8,height = 6)
ggsave("result/with_PVTT_MLN/1.scRNA_merge_vlnplot_before_QC_Pct.MT.pdf",p3,device = "pdf",width = 8,height = 6)

scRNA_merge_with_PVTT_MLN <- NormalizeData(scRNA_merge_with_PVTT_MLN) %>%
  FindVariableFeatures(selection.method="vst") %>% ScaleData()

scRNA_merge_with_PVTT_MLN <- RunPCA(scRNA_merge_with_PVTT_MLN,features = VariableFeatures(scRNA_merge_with_PVTT_MLN))
p1 <- DimPlot(scRNA_merge_with_PVTT_MLN,reduction = "pca",group.by = "orig.ident")
ggsave("result/with_PVTT_MLN/2.scRNA_merge_beforeQC_pca.pdf",p1,device = "pdf",width = 8,height = 6)
p2 <- ElbowPlot(scRNA_merge_with_PVTT_MLN,ndims = 30, reduction = "pca")
ggsave("result/with_PVTT_MLN/3.scRNA_merge_beforeQC_elbow.pdf",p2,device = "pdf",width = 6,height = 4)
pc.num = 1:22

# 聚类
scRNA_merge_with_PVTT_MLN <- FindNeighbors(scRNA_merge_with_PVTT_MLN, dims = pc.num) %>% FindClusters(resolution = 0.2)
scRNA_merge_with_PVTT_MLN <- RunUMAP(scRNA_merge_with_PVTT_MLN, reduction = "pca",dims = pc.num)
p1 <- DimPlot(scRNA_merge_with_PVTT_MLN, reduction = "umap")
p2 <- DimPlot(scRNA_merge_with_PVTT_MLN, reduction = "umap", group.by = "orig.ident") + theme(plot.title = element_blank())
pc <- p1 + p2
ggsave("result/with_PVTT_MLN/4.scRNA_merge_beforeQC_umap.pdf",pc,device = "pdf",width = 14,height = 4.5)

