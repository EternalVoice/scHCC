rm(list = ls());gc()
library(Seurat)
library(ggplot2)
library(dplyr)

scRNAlist_without_PVTT_MLN <- readRDS("rds/without_PVTT_MLN/1.scRNAlist_sub_without_PVTT_MLN.rds")

for (i in 1:length(scRNAlist_without_PVTT_MLN)){
  scRNAlist_without_PVTT_MLN[[i]] <- NormalizeData(scRNAlist_without_PVTT_MLN[[i]])
  scRNAlist_without_PVTT_MLN[[i]] <- FindVariableFeatures(scRNAlist_without_PVTT_MLN[[i]], selection.method = "vst")
}

# 整合
anchors <- FindIntegrationAnchors(object.list = scRNAlist_without_PVTT_MLN)
saveRDS(anchors,file = "rds/without_PVTT_MLN/anchors_without_PVTT_MLN.rds")
scRNA_integrated_without_PVTT_MLN <- IntegrateData(anchorset = anchors)
Idents(scRNA_integrated_without_PVTT_MLN) <- scRNA_integrated_without_PVTT_MLN$orig.ident
saveRDS(scRNA_integrated_without_PVTT_MLN,file = "rds/without_PVTT_MLN/total/1.scRNA_integrated_without_PVTT_MLN.rds")
col.num <- length(levels(as.factor(scRNA_integrated_without_PVTT_MLN$orig.ident)))
# 降维聚类
scRNA_integrated_without_PVTT_MLN <- ScaleData(scRNA_integrated_without_PVTT_MLN)
scRNA_integrated_without_PVTT_MLN <- RunPCA(scRNA_integrated_without_PVTT_MLN, features = VariableFeatures(scRNA_integrated_without_PVTT_MLN))
p1 <- DimPlot(scRNA_integrated_without_PVTT_MLN, reduction = "pca", group.by = "orig.ident",cols = rainbow(col.num))
ggsave("result/without_PVTT_MLN/total/integrated/1.scRNA_integrated_beforeQC_pca.pdf",p1,device="pdf",width=8,height=6)
p2 <- ElbowPlot(scRNA_integrated_without_PVTT_MLN,ndims=30,reduction="pca")
ggsave("result/without_PVTT_MLN/total/integrated/2.scRNA_integrated_beforeQC_elbow.pdf",p2,device="pdf",width=6,height=4)
pc.num = 1:24

# 聚类
scRNA_integrated_without_PVTT_MLN <- FindNeighbors(scRNA_integrated_without_PVTT_MLN, dims = pc.num) %>% FindClusters(resolution = 0.2)
scRNA_integrated_without_PVTT_MLN <- RunUMAP(scRNA_integrated_without_PVTT_MLN, reduction = "pca",dims = pc.num) %>% RunTSNE(dims=pc.num)

col.num2 <- length(levels(as.factor(scRNA_integrated_without_PVTT_MLN$seurat_clusters)))
p1 <- DimPlot(scRNA_integrated_without_PVTT_MLN, reduction = "tsne",cols = rainbow(col.num2))
p2 <- DimPlot(scRNA_integrated_without_PVTT_MLN, reduction = "tsne", group.by = "orig.ident",cols = rainbow(col.num)) + theme(plot.title = element_blank())
pc <- p1 + p2
ggsave("result/without_PVTT_MLN/total/integrated/3.scRNA_integrated_beforeQC_tsne.pdf",pc,device = "pdf",width = 14,height = 4.5)

# QC
DefaultAssay(scRNA_integrated_without_PVTT_MLN) <- "RNA"

# percent mt ribo
scRNA_integrated_without_PVTT_MLN[['percent.mt']] <- PercentageFeatureSet(scRNA_integrated_without_PVTT_MLN, pattern="^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA_integrated_without_PVTT_MLN@assays$RNA)) 
HB.genes <- rownames(scRNA_integrated_without_PVTT_MLN@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)]
scRNA_integrated_without_PVTT_MLN[['percent.hb']] <- PercentageFeatureSet(scRNA_integrated_without_PVTT_MLN, features=HB.genes)
scRNA_integrated_without_PVTT_MLN$project <- "Integrated"

violin <- VlnPlot(scRNA_integrated_without_PVTT_MLN, group.by = "project",
                  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                  pt.size = 0.00001, ncol = 3) + NoLegend()
ggsave("result/without_PVTT_MLN/total/integrated/4.scRNA_integrated_vlnplot_beforeQC.pdf",violin,device = "pdf",width = 8,height = 6)

# QC
minGene=200
maxGene=8000
pctMT=25

scRNA_integrated_filt_without_PVTT_MLN <- subset(scRNA_integrated_without_PVTT_MLN, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
violin <- VlnPlot(scRNA_integrated_filt_without_PVTT_MLN, group.by = "project",
                  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                  pt.size = 0.00001, ncol = 3) + NoLegend()
ggsave("result/without_PVTT_MLN/total/integrated/4.scRNA_integrated_vlnplot_afterQC.pdf",violin,device = "pdf",width = 8,height = 6)

p1 <- DimPlot(scRNA_integrated_filt_without_PVTT_MLN, reduction = "tsne",cols = rainbow(col.num2))
p2 <- DimPlot(scRNA_integrated_filt_without_PVTT_MLN, reduction = "tsne", group.by = "orig.ident",cols = rainbow(col.num)) + theme(plot.title = element_blank())
pc <- p1 + p2
ggsave("result/without_PVTT_MLN/total/integrated/3.scRNA_integrated_afterQC_tsne.pdf",pc,device = "pdf",width = 14,height = 4.5)
saveRDS(scRNA_integrated_filt_without_PVTT_MLN,file = "rds/without_PVTT_MLN/total/2.scRNA_integrated_DimRed_filt_without_PVTT_MLN.rds")

Idents(scRNA_integrated_filt_without_PVTT_MLN) <- scRNA_integrated_filt_without_PVTT_MLN$Sample
DefaultAssay(scRNA_integrated_filt_without_PVTT_MLN) <- "RNA"
table(scRNA_integrated_filt_without_PVTT_MLN$seurat_clusters,scRNA_integrated_filt_without_PVTT_MLN$Type)

current_cluster <- c(0:23)
new_cluster <- c("T cell","TAM","Malignant cell","TEC","Malignant cell","Malignant cell",
                 "T cell","Malignant cell","Malignant cell","TAM","CAF","Malignant cell",
                 "Malignant cell","B cell","B cell","Malignant cell","TAM","T cell","unclassified",
                 "unclassified","TEC","CAF","unclassified","TAM")
scRNA_integrated_filt_without_PVTT_MLN$CellType <- plyr::mapvalues(scRNA_integrated_filt_without_PVTT_MLN$seurat_clusters,
                                                                   from = current_cluster,to = new_cluster)
saveRDS(scRNA_integrated_filt_without_PVTT_MLN,file = "rds/without_PVTT_MLN/total/3.scRNA_integrated_filt_without_PVTT_MLN_celltypeAssigned.rds")

p1 <- DimPlot(scRNA_integrated_filt_without_PVTT_MLN,reduction = "tsne",group.by = "seurat_clusters",
              cols = rainbow(length(levels(as.factor(scRNA_integrated_filt_without_PVTT_MLN$seurat_clusters)))))
p2 <- DimPlot(scRNA_integrated_filt_without_PVTT_MLN,reduction = "tsne",group.by = "Sample",
              cols = rainbow(length(levels(as.factor(scRNA_integrated_filt_without_PVTT_MLN$Sample)))))
p3 <- DimPlot(scRNA_integrated_filt_without_PVTT_MLN,reduction = "tsne",group.by = "CellType",
              cols = rainbow(length(levels(as.factor(scRNA_integrated_filt_without_PVTT_MLN$CellType)))))
pc <- p1 + p2 + p3
ggsave("result/without_PVTT_MLN/total/integrated/5.scRNA_integrated_filt_tsne.pdf",pc,device = "pdf",width = 16,height = 4)
