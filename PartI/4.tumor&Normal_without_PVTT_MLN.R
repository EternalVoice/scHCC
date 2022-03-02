rm(list = ls());gc()
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

main_tiss <- readRDS("rds/without_PVTT_MLN/total/1.scRNA_integrated_without_PVTT_MLN.rds")

main_tiss$Tissue <- ifelse(
  main_tiss$Sample=="HCC03N" | main_tiss$Sample=="HCC04N" |
    main_tiss$Sample=="HCC05N" | main_tiss$Sample=="HCC06N" |
    main_tiss$Sample=="HCC07N" | main_tiss$Sample=="HCC08N" |
    main_tiss$Sample=="HCC09N" | main_tiss$Sample=="HCC10N",
  "normal", "tumor"
)

DefaultAssay(main_tiss) <- "RNA"
# 降维聚类
main_tiss <- FindVariableFeatures(main_tiss, selection.method = "vst", nfeatures = 2000)
DefaultAssay(main_tiss) <- "integrated"
main_tiss <- ScaleData(main_tiss, features = rownames(main_tiss))
main_tiss <- RunPCA(main_tiss, features = VariableFeatures(main_tiss))
p1 <- DimPlot(main_tiss, reduction = "pca", group.by = "orig.ident", cols = rainbow(length(levels(as.factor(main_tiss$orig.ident)))))
ggsave("result/without_PVTT_MLN/Main/1.main_tiss_beforeQC_pca.pdf", p1, device = "pdf", width = 8, height = 6)
p2 <- ElbowPlot(main_tiss, ndims = 30, reduction = "pca")
ggsave("result/without_PVTT_MLN/Main/2.main_tiss_beforeQC_elbow.pdf", p2, device = "pdf", width = 6, height = 4)

pc.num = 1:24

main_tiss <- FindNeighbors(main_tiss, dims = pc.num) %>% FindClusters(resolution = 0.2) %>% RunUMAP(reduction="pca",dims=pc.num) %>% RunTSNE(dims=pc.num)

p3 <- DimPlot(main_tiss, reduction = "tsne", cols = rainbow(length(levels(as.factor(main_tiss$seurat_clusters)))))
p4 <- DimPlot(main_tiss, reduction = "tsne", group.by = "orig.ident", cols = rainbow(length(levels(as.factor(main_tiss$orig.ident))))) + theme(plot.title = element_blank())
pc <- p3 + p4
ggsave("result/without_PVTT_MLN/Main/3.main_tiss_beforeQC_tsne.pdf", pc, device = "pdf", width = 14, height = 4.5)

DefaultAssay(main_tiss) <- "RNA"
main_tiss[["percent.mt"]] <- PercentageFeatureSet(main_tiss, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(main_tiss@assays$RNA))
HB.genes <- rownames(main_tiss@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
main_tiss[["percent.hb"]] <- PercentageFeatureSet(main_tiss, features = HB.genes)
main_tiss$project <- "Main"

violin <- VlnPlot(main_tiss, group.by = "project", features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0.00001, ncol=3) + NoLegend()
ggsave("result/without_PVTT_MLN/Main/4.main_tiss_violin_beforeQC.pdf", violin, device = "pdf", width = 8,height = 6)

# QC criteria
minGene = 200
maxGene = 8000
pctMT = 25

main_tiss <- subset(main_tiss, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)

violin <- VlnPlot(main_tiss, group.by = "project", features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0.00001,ncol=3) + NoLegend()
ggsave("result/without_PVTT_MLN/4.main_tiss_violin_afterQC.pdf", violin, device = "pdf", width = 8, height = 6)

p5 <- DimPlot(main_tiss, reduction = "tsne", cols = rainbow(length(levels(as.factor(main_tiss$seurat_clusters)))))
p6 <- DimPlot(main_tiss, reduction = "tsne", group.by = "orig.ident", cols = rainbow(length(levels(as.factor(main_tiss$orig.ident))))) + theme(plot.title=element_blank())
pc <- p5 + p6
ggsave("result/without_PVTT_MLN/Main/3.main_tiss_afterQC_tsne.pdf", pc, device = "pdf", width = 14, height = 4.5)

# cell type assignment
Idents(main_tiss) <- main_tiss$Sample
DefaultAssay(main_tiss) <- "RNA"
current_cluster <- c(0:23)
new_cluster <- c("TAM","T cell","Malignant cell","Malignant cell", "Malignant cell",
                 "Malignant cell","TEC","Malignant cell","CAF","Malignant cell",
                 "Malignant cell","B cell","T cell","Malignant cell","TAM",
                 "T cell","unclassified","B cell","TEC","CAF","unclassified","unclassified")
main_tiss$CellType <- plyr::mapvalues(main_tiss$seurat_clusters, from = current_cluster, to = new_cluster)
saveRDS(main_tiss, file = "rds/without_PVTT_MLN/Main/main_tiss_celltypeAssigned.rds")

p1 <- DimPlot(main_tiss,reduction = "tsne",group.by = "seurat_clusters",cols = rainbow(length(levels(as.factor(main_tiss$seurat_clusters)))))
p2 <- DimPlot(main_tiss,reduction = "tsne",group.by = "Sample",cols = rainbow(length(levels(as.factor(main_tiss$Sample)))))
p3 <- DimPlot(main_tiss,reduction = "tsne",group.by = "CellType",cols = rainbow(length(levels(as.factor(main_tiss$CellType)))))
pc <- p1 + p2 + p3
ggsave("result/without_PVTT_MLN/Main/5.tumor_tsne.pdf",pc,device = "pdf",width = 16,height = 4)


# 删除免疫细胞
