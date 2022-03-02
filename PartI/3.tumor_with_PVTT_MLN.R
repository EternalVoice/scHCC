rm(list = ls());gc()
library(Seurat)
library(ggplot2)
library(dplyr)

scRNA_integrated_with_PVTT_MLN <- readRDS("rds/with_PVTT_MLN/total/1.scRNA_integrated_with_PVTT_MLN.rds")

scRNA_integrated_with_PVTT_MLN$Tissue <- ifelse(
  scRNA_integrated_with_PVTT_MLN$Sample=="HCC03N" | scRNA_integrated_with_PVTT_MLN$Sample=="HCC04N" |
    scRNA_integrated_with_PVTT_MLN$Sample=="HCC05N" | scRNA_integrated_with_PVTT_MLN$Sample=="HCC06N" |
    scRNA_integrated_with_PVTT_MLN$Sample=="HCC07N" | scRNA_integrated_with_PVTT_MLN$Sample=="HCC08N" |
    scRNA_integrated_with_PVTT_MLN$Sample=="HCC09N" | scRNA_integrated_with_PVTT_MLN$Sample=="HCC10N",
  "normal", "tumor"
)

tumor_cells <- subset(scRNA_integrated_with_PVTT_MLN@meta.data, Tissue == "tumor")
tumor <- subset(scRNA_integrated_with_PVTT_MLN, cells = rownames(tumor_cells))

Idents(tumor) <- tumor$orig.ident
# 降维聚类
DefaultAssay(tumor) <- "RNA"
tumor <- FindVariableFeatures(tumor,selection.method = "vst",nfeatures = 2000)
DefaultAssay(tumor) <- "integrated"
tumor <- ScaleData(tumor,features = rownames(tumor))
tumor <- RunPCA(tumor, features = VariableFeatures(tumor))
p1 <- DimPlot(tumor, reduction = "pca", group.by = "orig.ident",cols=rainbow(length(levels(as.factor(tumor$orig.ident)))))
ggsave("result/with_PVTT_MLN/tumor/1.tumor_beforeQC_pca.pdf",p1,device="pdf",width=8,height=6)
p2 <- ElbowPlot(tumor,ndims=30,reduction="pca")
ggsave("result/with_PVTT_MLN/tumor/2.tumor_beforeQC_elbow.pdf",p2,device="pdf",width=6,height=4)

pc.num = 1:24
# 聚类
tumor <- FindNeighbors(tumor, dims = pc.num) %>% FindClusters(resolution = 0.2)
tumor <- RunUMAP(tumor,reduction="pca",dims=pc.num) %>% RunTSNE(dims=pc.num)

p1 <- DimPlot(tumor,reduction="tsne",cols=rainbow(length(levels(as.factor(tumor$seurat_clusters)))))
p2 <- DimPlot(tumor,reduction="tsne",group.by="orig.ident",cols=rainbow(length(levels(as.factor(tumor$orig.ident))))) + theme(plot.title = element_blank())
pc <- p1 + p2
ggsave("result/with_PVTT_MLN/tumor/3.tumor_beforeQC_tsne.pdf",pc,device = "pdf",width = 14,height = 4.5)

DefaultAssay(tumor) <- "RNA"
# percent mt ribo
tumor[['percent.mt']] <- PercentageFeatureSet(tumor, pattern="^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(tumor@assays$RNA)) 
HB.genes <- rownames(tumor@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)]
tumor[['percent.hb']] <- PercentageFeatureSet(tumor, features=HB.genes)
tumor$project <- "tumor"

violin <- VlnPlot(tumor, group.by = "project",
                  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                  pt.size = 0.00001, ncol = 3) + NoLegend()

ggsave("result/with_PVTT_MLN/tumor/4.tumor_vlnplot_beforeQC.pdf",violin,device = "pdf",width = 8,height = 6)

# QC
minGene=200
maxGene=8000
pctMT=25

tumor <- subset(tumor, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)

violin <- VlnPlot(tumor, group.by = "project",
                  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                  pt.size = 0.00001, ncol = 3) + NoLegend()

ggsave("result/with_PVTT_MLN/tumor/4.tumor_vlnplot_afterQC.pdf",violin,device = "pdf",width = 8,height = 6)

p1 <- DimPlot(tumor,reduction="tsne",cols=rainbow(length(levels(as.factor(tumor$seurat_clusters)))))
p2 <- DimPlot(tumor,reduction="tsne",group.by="orig.ident",cols=rainbow(length(levels(as.factor(tumor$orig.ident))))) + theme(plot.title=element_blank())
pc <- p1 + p2
ggsave("result/with_PVTT_MLN/tumor/3.tumor_afterQC_tsne.pdf",pc,device = "pdf",width = 14,height = 4.5)
saveRDS(tumor,file = "rds/with_PVTT_MLN/tumor/1.tumor_DimRed_filt_with_PVTT_MLN.rds")

Idents(tumor) <- tumor$Sample
DefaultAssay(tumor) <- "RNA"
table(tumor$seurat_clusters,tumor$Type)
# cell type assignment (preset)
current_cluster <- c(0:21)
new_cluster <- c("TAM","T cell","Malignant cell","Malignant cell", "Malignant cell",
                 "Malignant cell","Malignant cell","Malignant cell","TEC","CAF","Malignant cell",
                 "T cell","Malignant cell","B cell","TAM","TAM","T cell","B cell","TAM",
                 "Malignant cell","CAF","unclassified")
tumor$CellType <- plyr::mapvalues(tumor$seurat_clusters,from = current_cluster,to = new_cluster)
saveRDS(tumor,file = "rds/with_PVTT_MLN/tumor/2.tumor_filt_with_PVTT_MLN_celltypeAssigned.rds")

p1 <- DimPlot(tumor,reduction = "tsne",group.by = "seurat_clusters",cols = rainbow(length(levels(as.factor(tumor$seurat_clusters)))))
p2 <- DimPlot(tumor,reduction = "tsne",group.by = "Sample",cols = rainbow(length(levels(as.factor(tumor$Sample)))))
p3 <- DimPlot(tumor,reduction = "tsne",group.by = "CellType",cols = rainbow(length(levels(as.factor(tumor$CellType)))))
pc <- p1 + p2 + p3
ggsave("result/with_PVTT_MLN/tumor/5.tumor_tsne.pdf",pc,device = "pdf",width = 16,height = 4)

# CD3+: CD3D, CD3E, CD3G
# CD45+: PTPRC
# CD8+: CD8A, CD8B
# CD11b-: ITGAM
# B220: PTPRC
# CD11C-: ITGAX
# NK1.1-: KLRB1
# gdT-: CD27
markers <- c("CD3D","CD3E","CD3G","PTPRC","CD8A","CD8B","ITGAM","ITGAX","KLRB1","CD27")
Tcell <- DimPlot(tumor,group.by = "CellType",reduction = "tsne",cells.highlight = rownames(subset(tumor@meta.data,CellType=="T cell"))) + ggtitle("T cell") +NoLegend()
TAM <- DimPlot(tumor,group.by = "CellType",reduction = "tsne",cells.highlight = rownames(subset(tumor@meta.data,CellType=="TAM"))) + ggtitle("TAM") +NoLegend()
TEC <- DimPlot(tumor,group.by = "CellType",reduction = "tsne",cells.highlight = rownames(subset(tumor@meta.data,CellType=="TEC"))) + ggtitle("TEC") +NoLegend()
CAF <- DimPlot(tumor,group.by = "CellType",reduction = "tsne",cells.highlight = rownames(subset(tumor@meta.data,CellType=="CAF"))) + ggtitle("CAF") +NoLegend()
Bcell <- DimPlot(tumor,group.by = "CellType",reduction = "tsne",cells.highlight = rownames(subset(tumor@meta.data,CellType=="B cell"))) + ggtitle("B cell") +NoLegend()
MC <- DimPlot(tumor,group.by = "CellType",reduction = "tsne",cells.highlight = rownames(subset(tumor@meta.data,CellType=="Malignant cell"))) + ggtitle("Malignant cell") +NoLegend()
unclassified <- DimPlot(tumor,group.by = "CellType",reduction = "tsne",cells.highlight = rownames(subset(tumor@meta.data,CellType=="unclassified"))) + ggtitle("unclassified") +NoLegend()

T_sub_cells <- subset(tumor@meta.data, CellType == "T cell")
T_sub <- subset(tumor, cells = rownames(T_sub_cells))
