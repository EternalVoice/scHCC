#>>>>>>> 1.raw data preprocessing <<<<<<<<<<<<<<<<<<

library(ggplot2)
library(Seurat)
# 1. pre-processing

# input raw dataset 
total <- readRDS("../new/rds/without_PVTT_MLN/total/3.scRNA_integrated_filt_without_PVTT_MLN_celltypeAssigned.rds")

# p <- VlnPlot(tumor,features = "PTPRC",group.by = "CellType",pt.size = 0) + NoLegend()
# ggsave("0.PTPRC.pdf",p,width = 6,height = 3)

# Filtered unwanted patient samples --- H75
tumorFilt <- subset(total, cells = rownames(subset(total@meta.data, Sample != "H75")))
# p2 <- VlnPlot(tumorFilt, group.by = "Sample", features = "nFeature_RNA",pt.size = 0.001) + NoLegend()
# ggsave("../20210825/1.nFeature_RNA.pdf",p2,device = "pdf",width = 8,height = 4)

tumorFilt$Tissue <- ifelse(
  tumorFilt$Sample=="HCC03N" | tumorFilt$Sample=="HCC04N" | tumorFilt$Sample=="HCC05N" | tumorFilt$Sample=="HCC06N" |
    tumorFilt$Sample=="HCC07N" | tumorFilt$Sample=="HCC08N" | tumorFilt$Sample=="HCC09N" | tumorFilt$Sample=="HCC10N",
  "Normal", "Tumor"
)
tumorFilt@meta.data[,c(6:22)] <- NULL

saveRDS(tumorFilt, file = "rds/1.orignial_input.rds")


# CNP0000650
CNP0000650_matrix <- as.data.frame(data.table::fread("../new/input/new/CNP0000650/HCC_log_tpm_expression_matrix.txt.gz",header = T,sep = "\t",stringsAsFactors = F))

CNP0000650_matrix <- CNP0000650_matrix[!duplicated(CNP0000650_matrix[,1]),]
rownames(CNP0000650_matrix) <- CNP0000650_matrix[,1]
CNP0000650_matrix <- CNP0000650_matrix[,-1]

metadata <- read.table("../new/input/new/CNP0000650/HCC_cell_metadata.txt",header = T,sep = "\t",skip = 1)
rownames(metadata) <- metadata$type;
metadata$sample_attribute <- NULL
metadata$group <- stringr::str_split(metadata$group, "_",simplify = T)[,2]
metadata$type <- stringr::str_split(metadata$type, "_", simplify = T)[,1]

colnames(metadata)[1] <- "Sample";colnames(metadata)[3] <- "Tissue"
metadata$experiment_attribute <- NULL
metadata$Tissue <- gsub("Adjacent liver","Normal",metadata$Tissue)
save(CNP0000650_matrix,metadata,file = "rds/CNP0000650.Rdata")

CNP0000650 <- CreateSeuratObject(counts = CNP0000650_matrix, meta.data = metadata)

CNP0000650_barcode <- rownames(subset(CNP0000650@meta.data, group.2 != "Relapsed"))
CNP0000650 <- subset(CNP0000650, cells = rownames(subset(CNP0000650@meta.data, group.2 != "Relapsed")))
CNP0000650$group.2 <- as.character(CNP0000650$group.2)
saveRDS(CNP0000650,file = "rds/CNP0000650.filt.rds")

# Data Integration

rm(list = ls());gc()
library(Seurat)
library(ggplot2)
library(dplyr)

ori <- readRDS("rds/1.original_input.rds")
CNP0000650 <- readRDS("rds/CNP0000650.filt.rds")

ori <- NormalizeData(ori) %>% FindVariableFeatures(selection.method = "vst")
CNP0000650 <- NormalizeData(CNP0000650) %>% FindVariableFeatures(selection.method = "vst")

# 整合
anchors <- FindIntegrationAnchors(object.list = list(ori,CNP0000650))

saveRDS(anchors,file = "rds/anchors.rds")
scRNA_integrated <- IntegrateData(anchorset = anchors)

scRNA_integrated$group <- NULL
scRNA_integrated$group.2 <- NULL

saveRDS(scRNA_integrated,file = "rds/scRNA_integrated.rds")

scRNA_integrated <- readRDS("rds/scRNA_integrated.rds")

col.num <- length(levels(as.factor(scRNA_integrated$orig.ident)))


# 降维聚类
scRNA_integrated <- ScaleData(scRNA_integrated)
scRNA_integrated <- RunPCA(scRNA_integrated, features = VariableFeatures(scRNA_integrated))
p1 <- DimPlot(scRNA_integrated, reduction = "pca", group.by = "orig.ident",cols = rainbow(col.num))
ggsave("results/1.scRNA_integrated_beforeQC_pca.pdf",p1,device="pdf",width=8,height=6)
p2 <- ElbowPlot(scRNA_integrated,ndims=30,reduction="pca")
ggsave("results/2.scRNA_integrated_beforeQC_elbow.pdf",p2,device="pdf",width=6,height=4)

pc.num = 1:28

# 聚类
scRNA_integrated <- FindNeighbors(scRNA_integrated, dims = pc.num) %>% FindClusters(resolution = 0.5)
scRNA_integrated <- RunUMAP(scRNA_integrated, reduction = "pca",dims = pc.num) %>% RunTSNE(dims=pc.num)

col.num2 <- length(levels(as.factor(scRNA_integrated$seurat_clusters)))
p1 <- DimPlot(scRNA_integrated, reduction = "tsne",cols = rainbow(col.num2))
p2 <- DimPlot(scRNA_integrated, reduction = "tsne", group.by = "orig.ident",cols = rainbow(col.num)) + theme(plot.title = element_blank())
pc <- p1 + p2
ggsave("results/3.scRNA_integrated_beforeQC_tsne.pdf",pc,device = "pdf",width = 14,height = 4.5)

# QC
DefaultAssay(scRNA_integrated) <- "RNA"

# percent mt ribo
scRNA_integrated[['percent.mt']] <- PercentageFeatureSet(scRNA_integrated, pattern="^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA_integrated@assays$RNA)) 
HB.genes <- rownames(scRNA_integrated@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)]
scRNA_integrated[['percent.hb']] <- PercentageFeatureSet(scRNA_integrated, features=HB.genes)
scRNA_integrated$project <- "Integrated"

violin <- VlnPlot(scRNA_integrated, group.by = "project",
                  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                  pt.size = 0.00001, ncol = 3) + NoLegend()
ggsave("results/4.scRNA_integrated_vlnplot_beforeQC.pdf",violin,device = "pdf",width = 8,height = 6)

# QC
minGene=200
maxGene=10000
pctMT=25

scRNA_integrated <- subset(scRNA_integrated, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)
violin <- VlnPlot(scRNA_integrated, group.by = "project",
                  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                  pt.size = 0.00001, ncol = 3) + NoLegend()
ggsave("results/4.scRNA_integrated_vlnplot_afterQC.pdf",violin,device = "pdf",width = 8,height = 6)

p1 <- DimPlot(scRNA_integrated, reduction = "tsne",cols = rainbow(col.num2))
p2 <- DimPlot(scRNA_integrated, reduction = "tsne", group.by = "orig.ident",cols = rainbow(col.num)) + theme(plot.title = element_blank())
pc <- p1 + p2
ggsave("results/3.scRNA_integrated_afterQC_tsne.pdf",pc,device = "pdf",width = 14,height = 4.5)
saveRDS(scRNA_integrated,file = "rds/2.scRNA_integrated_DimRed_filt.rds")


#>>>>>>>>> TAM <<<<<<<<<<<<<<<<<<<

# TAM: CD45+CD3e-CD19-CD56-CD16-CD206+CD163+HLA-DR+CD11C+CD14+
# TAM: PTPRC,CD3E,CD19,NCAM1,FCGR3A,MRC1,CD163,HLA-DRA,ITGAX,CD14

markers <- c("PTPRC","CD3E","CD19","NCAM1","FCGR3A","MRC1","CD163","HLA-DRA","ITGAX","CD14")

for (marker in markers) {
  p <- VlnPlot(scRNA_integrated,features = marker, group.by = "seurat_clusters",pt.size = 0) + 
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank())
  ggsave(paste0("results/5.",marker,".pdf"),p,width = 10,height = 2)
}

p <- VlnPlot(scRNA_integrated,features = "MRC1", group.by = "seurat_clusters") + 
  NoLegend() +
  theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank())
ggsave(paste0("results/5.","MRC1",".pdf"),p,width = 10,height = 2)

p_Tsne <- DimPlot(scRNA_integrated,reduction = "tsne",label = T,label.size = 5) +NoLegend() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 1.5,colour = "black"))
ggsave("results/6.tsne.pdf",p_Tsne,height = 6,width = 6)

# TAM
clusters <- c(1,10,12,16,24,30,31)
cells <- NULL
for(cluster in clusters){
  cell <- rownames(subset(scRNA_integrated@meta.data,seurat_clusters == cluster))
  cells <- c(cells,cell)
}

TAM <- subset(scRNA_integrated,cells = cells)
TAM$integrated_snn_res.0.5 <- NULL;TAM$percent.mt <- NULL; TAM$percent.hb <- NULL
DefaultAssay(TAM) <- "integrated"
TAM <- ScaleData(TAM)
TAM <- RunPCA(TAM, features = VariableFeatures(TAM))
p2 <- ElbowPlot(TAM,ndims=30,reduction="pca")
pc.num = 1:24

# 聚类
TAM <- FindNeighbors(TAM, dims = pc.num) %>% FindClusters(resolution = 0.2)
TAM <- RunUMAP(TAM, reduction = "pca",dims = pc.num) %>% RunTSNE(dims=pc.num)

p1 <- DimPlot(TAM,reduction = "tsne",label = TRUE,label.size = 5) + NoLegend()
p2 <- DimPlot(TAM,reduction = "tsne",group.by = "Tissue") + 
  theme(plot.title = element_blank(),
        legend.position = "top")
p <- p1 + p2
ggsave("results/7.TAM.pdf",p,width = 8,height = 4)

saveRDS(TAM,file = "rds/TAM.rds")


#>>>>>>>>>>>>>>> MDSC <<<<<<<<<<<<<<<<<<<

# CD45+ CD3- B220- NK1.1- gdTCR- CD11b+ CD14- CD15+ CD33+ CD66b+


#>>>>>>>>>>>>>>> G-MDSC <<<<<<<<<<<<<<<<<<<

# CD14-, CD15+, CD66b+

# PTPRC+, CD3-, KLRB1-, ITGAM+, CD14-, FUT4+, CD33+, CEACAM8+
markers <- c("PTPRC","CD3D","CD3E","CD3G","KLRB1","ITGAM","CD14","FUT4","CD33","CEACAM8")

tmp1 <- c("CD3D","CD3E","CD3G","KLRB1","CD14")
for (marker in tmp1) {
  p <- VlnPlot(scRNA_integrated,features = marker, group.by = "seurat_clusters",pt.size = 0) + 
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank())
  ggsave(paste0("results/8.",marker,".pdf"),p,width = 10,height = 2)
}

tmp2 <- c("PTPRC","ITGAM","FUT4","CD33","CEACAM8")
for(marker in tmp2){
  p <- VlnPlot(scRNA_integrated,features = marker, group.by = "seurat_clusters") + 
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank())
  ggsave(paste0("results/8.",marker,".pdf"),p,width = 10,height = 2)
}

clusters <- c(3,4,7,9,13,14,17,18,26,36,37)

cells <- NULL
for(cluster in clusters){
  cell <- rownames(subset(scRNA_integrated@meta.data,seurat_clusters == cluster))
  cells <- c(cells,cell)
}

GMDSC <- subset(scRNA_integrated,cells = cells)
GMDSC$integrated_snn_res.0.5 <- NULL;GMDSC$percent.mt <- NULL; GMDSC$percent.hb <- NULL
DefaultAssay(GMDSC) <- "integrated"
GMDSC <- ScaleData(GMDSC)
GMDSC <- RunPCA(GMDSC, features = VariableFeatures(GMDSC))
p2 <- ElbowPlot(GMDSC,ndims=30,reduction="pca")
pc.num = 1:28

# 聚类
GMDSC <- FindNeighbors(GMDSC, dims = pc.num) %>% FindClusters(resolution = 0.2)
GMDSC <- RunUMAP(GMDSC, reduction = "pca",dims = pc.num) %>% RunTSNE(dims=pc.num)

p1 <- DimPlot(GMDSC,reduction = "tsne",label = TRUE,label.size = 5) + NoLegend()
p2 <- DimPlot(GMDSC,reduction = "tsne",group.by = "Tissue") + 
  theme(plot.title = element_blank(),
        legend.position = "top")
p <- p1 + p2
ggsave("results/GMDSC/8.GMDSC.pdf",p,width = 8,height = 4)

saveRDS(GMDSC,file = "rds/GMDSC.rds")


#>>>>>>>>>>>>>>> M-MDSC <<<<<<<<<<<<<<<<<<

# CD14+, CD15-, CD66b-
# PTPRC+, CD3-, KLRB1-, ITGAM+, CD14+, FUT4-, CD33+, CEACAM8-
markers <- c("PTPRC","CD3D","CD3E","CD3G","KLRB1","ITGAM","CD14","FUT4","CD33","CEACAM8")

tmp1 <- c("CD3D","CD3E","CD3G","KLRB1","FUT4","CEACAM8")
for (marker in tmp1) {
  p <- VlnPlot(scRNA_integrated,features = marker, group.by = "seurat_clusters",pt.size = 0) + 
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank())
  ggsave(paste0("results/MMDSC/9.",marker,".pdf"),p,width = 10,height = 2)
}

tmp2 <- c("PTPRC","ITGAM","CD14","CD33")
for(marker in tmp2){
  p <- VlnPlot(scRNA_integrated,features = marker, group.by = "seurat_clusters") + 
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank())
  ggsave(paste0("results/MMDSC/9.",marker,".pdf"),p,width = 10,height = 2)
}

clusters <- c(4,6,8,16,19,21,23,24,27,30,31,39)

cells <- NULL
for(cluster in clusters){
  cell <- rownames(subset(scRNA_integrated@meta.data,seurat_clusters == cluster))
  cells <- c(cells,cell)
}

MMDSC <- subset(scRNA_integrated,cells = cells)
MMDSC$integrated_snn_res.0.5 <- NULL;MMDSC$percent.mt <- NULL; MMDSC$percent.hb <- NULL
DefaultAssay(MMDSC) <- "integrated"
MMDSC <- ScaleData(MMDSC)
MMDSC <- RunPCA(MMDSC, features = VariableFeatures(MMDSC))
p2 <- ElbowPlot(MMDSC,ndims=30,reduction="pca")
pc.num = 1:25

# 聚类
MMDSC <- FindNeighbors(MMDSC, dims = pc.num) %>% FindClusters(resolution = 0.2)
MMDSC <- RunUMAP(MMDSC, reduction = "pca",dims = pc.num) %>% RunTSNE(dims=pc.num)

p1 <- DimPlot(MMDSC,reduction = "tsne",label = TRUE,label.size = 5) + NoLegend()
p2 <- DimPlot(MMDSC,reduction = "tsne",group.by = "Tissue") + 
  theme(plot.title = element_blank(),
        legend.position = "top")
p <- p1 + p2
ggsave("results/MMDSC/9.MMDSC.pdf",p,width = 8,height = 4)

saveRDS(MMDSC,file = "rds/MMDSC.rds")


#>>>>>>>>> TEX exhausted CD8+T cells <<<<<<<<<<<<<<<<<<

# CD8_Tex: CD45+CD8+CD3+CD223+ CD279+ CD366+ (TIM-3) Tox+
# CD8_Tex: PTPRC,CD8A,CD8B,CD3D,CD3E,CD3G,LAG3,PDCD1,HAVCR2,TOX

markers <- c("PTPRC","CD8A","CD8B","CD3D","CD3E","CD3G","LAG3","PDCD1","HAVCR2","TOX")

for(marker in markers){
  p <- VlnPlot(scRNA_integrated,features = marker, group.by = "seurat_clusters") + 
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank())
  ggsave(paste0("results/CD8_Tex/10.",marker,".pdf"),p,width = 10,height = 2)
}

clusters <- c(0,1,2,4,5,10,11,12,15,16,20,22,29,32,37,38)

cells <- NULL
for(cluster in clusters){
  cell <- rownames(subset(scRNA_integrated@meta.data,seurat_clusters == cluster))
  cells <- c(cells,cell)
}

CD8_Tex <- subset(scRNA_integrated,cells = cells)
CD8_Tex$integrated_snn_res.0.5 <- NULL;CD8_Tex$percent.mt <- NULL; CD8_Tex$percent.hb <- NULL
DefaultAssay(CD8_Tex) <- "integrated"
CD8_Tex <- ScaleData(CD8_Tex)
CD8_Tex <- RunPCA(CD8_Tex, features = VariableFeatures(CD8_Tex))
p2 <- ElbowPlot(CD8_Tex,ndims=30,reduction="pca")
pc.num = 1:26

# 聚类
CD8_Tex <- FindNeighbors(CD8_Tex, dims = pc.num) %>% FindClusters(resolution = 0.2)
CD8_Tex <- RunUMAP(CD8_Tex, reduction = "pca",dims = pc.num) %>% RunTSNE(dims=pc.num)

p1 <- DimPlot(CD8_Tex,reduction = "tsne",label = TRUE,label.size = 5) + NoLegend()
p2 <- DimPlot(CD8_Tex,reduction = "tsne",group.by = "Tissue") + 
  theme(plot.title = element_blank(),
        legend.position = "top")
p <- p1 + p2
ggsave("results/CD8_Tex/10.CD8_Tex.pdf",p,width = 8,height = 4)

saveRDS(CD8_Tex,file = "rds/CD8_Tex.rds")



#>>>>>>>>> TEX exhausted CD4+T cells <<<<<<<<<<<<<<<<<<

# CD4_Tex: CD45+CD4+CD3+CD223+ CD279+ CD366+ (TIM-3) Tox+
# CD4_Tex: PTPRC,CD4,CD3D,CD3E,CD3G,LAG3,PDCD1,HAVCR2,TOX
markers <- c("PTPRC","CD4","CD3D","CD3E","CD3G","LAG3","PDCD1","HAVCR2","TOX")

for(marker in markers){
  p <- VlnPlot(scRNA_integrated,features = marker, group.by = "seurat_clusters") + 
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank())
  ggsave(paste0("results/CD4_Tex/11.",marker,".pdf"),p,width = 10,height = 2)
}

clusters <- c(5,7,8,9,13,14,17,18,20,21,22,30,31,34,39,40)

cells <- NULL
for(cluster in clusters){
  cell <- rownames(subset(scRNA_integrated@meta.data,seurat_clusters == cluster))
  cells <- c(cells,cell)
}

CD4_Tex <- subset(scRNA_integrated,cells = cells)
CD4_Tex$integrated_snn_res.0.5 <- NULL;CD4_Tex$percent.mt <- NULL; CD4_Tex$percent.hb <- NULL
DefaultAssay(CD4_Tex) <- "integrated"
CD4_Tex <- ScaleData(CD4_Tex)
CD4_Tex <- RunPCA(CD4_Tex, features = VariableFeatures(CD4_Tex))
p2 <- ElbowPlot(CD4_Tex,ndims=30,reduction="pca")
pc.num = 1:26

# 聚类
CD4_Tex <- FindNeighbors(CD4_Tex, dims = pc.num) %>% FindClusters(resolution = 0.2)
CD4_Tex <- RunUMAP(CD4_Tex, reduction = "pca",dims = pc.num) %>% RunTSNE(dims=pc.num)

p1 <- DimPlot(CD4_Tex,reduction = "tsne",label = TRUE,label.size = 5) + NoLegend()
p2 <- DimPlot(CD4_Tex,reduction = "tsne",group.by = "Tissue") + 
  theme(plot.title = element_blank(),
        legend.position = "top")
p <- p1 + p2
ggsave("results/CD4_Tex/11.CD4_Tex.pdf",p,width = 8,height = 4)

saveRDS(CD4_Tex,file = "rds/CD4_Tex.rds")


#>>>>>>>>> CTL <<<<<<<<<<<<<<<<<<

# CTL: CD45+CD3+ CD8+ CD107a+ (LAMP-1) GranzymeB+ IFN-γ+
# CTL: PTPRC,CD3D,CD3E,CD3G,CD8A,CD8B,LAMP1,GZMB,IFNG

markers <- c("PTPRC","CD3D","CD3E","CD3G","CD8A","CD8B","LAMP1","GZMB","IFNG")

for(marker in markers){
  p <- VlnPlot(scRNA_integrated,features = marker, group.by = "seurat_clusters") + 
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank())
  ggsave(paste0("results/CTL/12.",marker,".pdf"),p,width = 10,height = 2)
}

clusters <- c(0,1,2,4,5,7,8,9,10,11,12,13,14,15,16,19,20,21,22,29,30,32,37,38)

cells <- NULL
for(cluster in clusters){
  cell <- rownames(subset(scRNA_integrated@meta.data,seurat_clusters == cluster))
  cells <- c(cells,cell)
}

CTL <- subset(scRNA_integrated,cells = cells)
CTL$integrated_snn_res.0.5 <- NULL;CTL$percent.mt <- NULL; CTL$percent.hb <- NULL
DefaultAssay(CTL) <- "integrated"
CTL <- ScaleData(CTL)
CTL <- RunPCA(CTL, features = VariableFeatures(CTL))
p2 <- ElbowPlot(CTL,ndims=30,reduction="pca")
pc.num = 1:28

# 聚类
CTL <- FindNeighbors(CTL, dims = pc.num) %>% FindClusters(resolution = 0.2)
CTL <- RunUMAP(CTL, reduction = "pca",dims = pc.num) %>% RunTSNE(dims=pc.num)

p1 <- DimPlot(CTL,reduction = "tsne",label = TRUE,label.size = 5) + NoLegend()
p2 <- DimPlot(CTL,reduction = "tsne",group.by = "Tissue") + 
  theme(plot.title = element_blank(),
        legend.position = "top")
p <- p1 + p2
ggsave("results/CTL/12.CTL.pdf",p,width = 8,height = 4)

saveRDS(CTL,file = "rds/CTL.rds")

#>>>>>>>>> Th1 <<<<<<<<<<<<<<<<<<

# Th1: CD45+CD3+ CD4+Tbet+ IFN-γ+CXCR3+
# Th1: PTPRC,CD3D,CD3E,CD3G,CD4,TBX21,IFNG,CXCR3

markers <- c("PTPRC","CD3D","CD3E","CD3G","CD4","TBX21","IFNG","CXCR3")

for(marker in markers){
  p <- VlnPlot(scRNA_integrated,features = marker, group.by = "seurat_clusters") + 
    NoLegend() +
    theme(plot.title = element_text(hjust = 0.5),axis.title.x = element_blank())
  ggsave(paste0("results/Th1/13.",marker,".pdf"),p,width = 10,height = 2)
}

markers <- c(0,1,2,3,4,5,12,13,14,15,16,20,21,22,29,30,32,37,38,39,40)

cells <- NULL
for(cluster in clusters){
  cell <- rownames(subset(scRNA_integrated@meta.data,seurat_clusters == cluster))
  cells <- c(cells,cell)
}

Th1 <- subset(scRNA_integrated,cells = cells)
Th1$integrated_snn_res.0.5 <- NULL;Th1$percent.mt <- NULL; Th1$percent.hb <- NULL
DefaultAssay(Th1) <- "integrated"
Th1 <- ScaleData(Th1)
Th1 <- RunPCA(Th1, features = VariableFeatures(Th1))
p2 <- ElbowPlot(Th1,ndims=30,reduction="pca")
pc.num = 1:26

# 聚类
Th1 <- FindNeighbors(Th1, dims = pc.num) %>% FindClusters(resolution = 0.2)
Th1 <- RunUMAP(Th1, reduction = "pca",dims = pc.num) %>% RunTSNE(dims=pc.num)

p1 <- DimPlot(Th1,reduction = "tsne",label = TRUE,label.size = 5) + NoLegend()
p2 <- DimPlot(Th1,reduction = "tsne",group.by = "Tissue") + 
  theme(plot.title = element_blank(),legend.position = "top")
p <- p1 + p2
ggsave("results/Th1/13.Th1.pdf",p,width = 8,height = 4)

saveRDS(Th1,file = "rds/Th1.rds")


##>>>>>>>>>>>> DE genes and Correlation Analysis <<<<<<<<<<<<<<<<<<<<


## >>>>>>>> Step 1. Find DE markers <<<<<<<<<<<<<<<<

DeGenePipeline <- function(object,cell) {
  DefaultAssay(object) <- "RNA"
  markers <- FindMarkers(object, ident.1 = "Tumor", ident.2 = "Normal", group.by = "Tissue", only.pos = TRUE)
  markers <- markers[order(markers$avg_logFC, decreasing = TRUE),]
  p1 <- DimPlot(object, group.by = "Tissue") + theme(plot.title = element_blank(), legend.position = "top")
  xlsx::write.xlsx(markers, file = paste0("results/",cell,"/markers.xlsx"))
  
  for(i in seq_along(rownames(markers))) {
    p2 <- FeaturePlot(object, features = rownames(markers)[i])
    p <- p1 + p2
    ggsave(paste0("results/",cell,"/DEmarkers/",i,"_",rownames(markers)[i],".pdf"),p,width = 8,height = 3.6)
  }
  
  highlight <- DimPlot(object,group.by = "Tissue",pt.size = 0.8,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
                       cells.highlight = rownames(subset(object@meta.data,Tissue=="Tumor"))) + NoLegend() +
    theme(panel.border = element_rect(colour = "black", size = 1), plot.title = element_blank())
  ggsave(paste0("results/",cell,"/tumor_highlight.pdf"),highlight,width = 6,height = 6)
}

#>>> TAM 
TAM <- readRDS("rds/TAM.rds")
DeGenePipeline(object = TAM,cell = "TAM")
rm(TAM);gc()

#>>> G-MDSC
GMDSC <- readRDS("rds/GMDSC.rds")
DeGenePipeline(object = GMDSC,cell = "GMDSC")
rm(GMDSC);gc()

#>>> M-MDSC
MMDSC <- readRDS("rds/MMDSC.rds")
DeGenePipeline(object = MMDSC,cell = "MMDSC")
rm(MMDSC);gc()

#>>> CD8_Tex
CD8_Tex <- readRDS("rds/CD8_Tex.rds")
DeGenePipeline(object = CD8_Tex,cell = "CD8_Tex")
rm(CD8_Tex);gc()

#>>> CD4_Tex
CD4_Tex <- readRDS("rds/CD4_Tex.rds")
DeGenePipeline(object = CD4_Tex,cell = "CD4_Tex")
rm(CD4_Tex);gc()

#>>> CTL
CTL <- readRDS("rds/CTL.rds")
DeGenePipeline(object = CTL,cell = "CTL")
rm(CTL);gc()

#>>> Th1
Th1 <- readRDS("rds/Th1.rds")
DeGenePipeline(object = Th1,cell = "Th1")
rm(Th1);gc()

## >>>>>>>> Step 2. Correlation Analysis <<<<<<<<<<<<<<<<

#>>> TAM 

## 分析TAM中高表达的基因，并分析它们的表达在数据集以及在TCGA中与CD163、CD206、CD8、IFNg的相关性。

rm(list = ls());gc()
scRNA <- readRDS("rds/2.scRNA_integrated_DimRed_filt.rds")

markers <- xlsx::read.xlsx("results/TAM/markers.xlsx",sheetIndex = 1)

gene <- c("CD163","MRC1","CD8A","CD8B","IFNG",as.character(markers[,1]))

# singlecell
expr <- GetAssayData(scRNA,assay = "RNA",slot = "data")[gene,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(scRNA), to = scRNA$Sample)
colnames(expr.T) <- c(gene,"Sample")

AvgExpPerSample <- function(expMtx){
  patient <- NULL
  for(i in 1:(ncol(expMtx)-1)){
    patient[i] <- data.frame(tapply(expMtx[,i],expMtx$Sample,sum))
  }
  names(patient) <- colnames(expMtx)[-ncol(expMtx)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expMtx[,1],expMtx$Sample,sum))
  Patient_sample <- as.data.frame(table(expMtx$Sample),stringsAsFactors = F)
  patient$Sample <- plyr::mapvalues(rownames(patient),from = Patient_sample$Var1,to = Patient_sample$Freq)
  patient <- as.data.frame(patient,stringsAsFactors = F)
  patient$Sample <- as.numeric(patient$Sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$Sample[j]
    }
  }
  patient$Sample <- NULL
  return(patient)
}

patient <- AvgExpPerSample(expr.T)
colnames(patient) <- gene
xlsx::write.xlsx(patient,file = "results/TAM/avgExpPerPatient.xlsx")

colnames(patient) <- gsub("-",".",colnames(patient))
pwd = "results/TAM/corr/singlecell/"
for(gene in colnames(patient)[1:5]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

for(i in 1:5){
  for(j in 6:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[i], y = colnames(patient)[j])) +
      geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,colnames(patient)[i],"/",colnames(patient)[j],"_vs_",colnames(patient)[i],".pdf"),p2,width = 6,height = 6)
  }
}

# TCGA
load("../new/refdb/TCGA_LIHC.rnaseq.RData")

genesTBL <- xlsx::read.xlsx("results/TAM/avgExpPerPatient.xlsx",sheetIndex = 1)

genes <- gsub("\\.","-",colnames(genesTBL)[-1])
genes <- genes[-grep("-",genes)]

setdiff(genes, colnames(LIHC.rnaseq))

genes <- genes[-grep("MTRNR2L12",genes)];genes <- genes[-grep("MTRNR2L8",genes)]
genes[grep("SELENOP",genes)] <- "SEPP1";genes[grep("PKM",genes)] <- "PKM2"
genes <- genes[-grep("MTRNR2L1",genes)];genes[grep("RACK1",genes)] <- "GNB2L1"
genes[grep("ERO1A",genes)] <- "ERO1L";genes <- genes[-grep("SNHG25",genes)]
genes[grep("SELENOK",genes)] <- "SELK";genes[grep("COPS9",genes)] <- "MYEOV2"
genes[grep("SELENOW",genes)] <- "SEPW1";genes[grep("LAMTOR2",genes)] <- "ROBLD3"
genes[grep("DNPH1",genes)] <- "C6orf108";genes[grep("SELENOS",genes)] <- "SELS"
genes[grep("MYDGF",genes)] <- "IL27";genes <- genes[-grep("HIST2H2AA4",genes)]
genes <- genes[-grep("SMIM4",genes)];genes[grep("MPC2",genes)] <- "BRP44"
genes[grep("SELENOH",genes)] <- "C11orf31";genes[grep("COX20",genes)] <- "FAM36A"
genes[grep("WDR45B",genes)] <- "WDR45L";genes[grep("MRNIP",genes)] <- "C5orf45"
genes[grep("GLMP",genes)] <- "C1orf85";genes[grep("CYSTM1",genes)] <- "C5orf32"
genes[grep("SELENOT",genes)] <- "SELT";genes[grep("SELENOF",genes)] <- "SEP15"
genes[grep("HYPK",genes)] <- "C15orf63";genes[grep("CCDC167",genes)] <- "C6orf129"
genes[grep("SLC50A1",genes)] <- "RAG1AP1";genes[grep("UQCC2",genes)] <- "C6orf125"

expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)
pwd = "results/TAM/corr/TCGA/"

for(gene in colnames(expr.norm)[1:5]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

colnames(expr.norm) <- gsub("-",".",colnames(expr.norm))

for(i in 1:5){
  for(j in 6:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x = colnames(expr.norm)[i], y = colnames(expr.norm)[j])) +
      geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,colnames(expr.norm)[i],"/",colnames(expr.norm)[j],"_vs_",colnames(expr.norm)[i],".pdf"),p2,width = 6,height = 6)
  }
}


#>>> G-MDSC 

## 分析G-MDSC中高表达的基因，并分析它们的表达在数据集以及在TCGA中与CD33,CD66b,CD11B,CD8,IFNg的相关性。

markers <- xlsx::read.xlsx("results/GMDSC/markers.xlsx",sheetIndex = 1)
gene <- c("CD33","CEACAM8","ITGAM","CD8A","CD8B","IFNG",as.character(markers[,1]))

# singlecell
expr <- GetAssayData(scRNA,assay = "RNA",slot = "data")[gene,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(scRNA), to = scRNA$Sample)
colnames(expr.T) <- c(gene,"Sample")

AvgExpPerSample <- function(expMtx){
  patient <- NULL
  for(i in 1:(ncol(expMtx)-1)){
    patient[i] <- data.frame(tapply(expMtx[,i],expMtx$Sample,sum))
  }
  names(patient) <- colnames(expMtx)[-ncol(expMtx)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expMtx[,1],expMtx$Sample,sum))
  Patient_sample <- as.data.frame(table(expMtx$Sample),stringsAsFactors = F)
  patient$Sample <- plyr::mapvalues(rownames(patient),from = Patient_sample$Var1,to = Patient_sample$Freq)
  patient <- as.data.frame(patient,stringsAsFactors = F)
  patient$Sample <- as.numeric(patient$Sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$Sample[j]
    }
  }
  patient$Sample <- NULL
  return(patient)
}

patient <- AvgExpPerSample(expr.T)
colnames(patient) <- gene
xlsx::write.xlsx(patient,file = "results/GMDSC/avgExpPerPatient.xlsx")

colnames(patient) <- gsub("-",".",colnames(patient))

pwd = "results/GMDSC/corr/singlecell/"
for(gene in colnames(patient)[1:6]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

for(i in 1:6){
  for(j in 7:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[i], y = colnames(patient)[j])) +
      geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,colnames(patient)[i],"/",colnames(patient)[j],"_vs_",colnames(patient)[i],".pdf"),p2,width = 6,height = 6)
  }
}

# TCGA
load("../new/refdb/TCGA_LIHC.rnaseq.RData")

genesTBL <- xlsx::read.xlsx("results/GMDSC/avgExpPerPatient.xlsx",sheetIndex = 1)

genes <- gsub("\\.","-",colnames(genesTBL)[-1])
genes <- genes[-grep("-",genes)]

setdiff(genes, colnames(LIHC.rnaseq))

genes <- genes[-grep("SNHG25",genes)];genes <- genes[-grep("MTRNR2L8",genes)]
genes[grep("MPC2",genes)] <- "BRP44";genes <- genes[-grep("MTRNR2L12",genes)]
genes <- genes[-grep("GAGE12H",genes)];genes[grep("PKM",genes)] <- "PKM2"
genes[grep("UQCC2",genes)] <- "C6orf125";genes <- genes[-grep("MTRNR2L10",genes)]
genes[grep("COA3",genes)] <- "CCDC56";genes[grep("COX20",genes)] <- "FAM36A"
genes[grep("SLC50A1",genes)] <- "RAG1AP1";genes <- genes[-grep("SMIM4",genes)]
genes[grep("SLIRP",genes)] <- "C14orf156";genes[grep("LAMTOR5",genes)] <- "HBXIP"
genes[grep("FUOM",genes)] <- "C10orf125";genes <- genes[-grep("LINC00261",genes)]
genes[grep("MZT2B",genes)] <- "FAM128B";genes[grep("COA4",genes)] <- "CHCHD8"
genes[grep("DNPH1",genes)] <- "C6orf108";genes <- genes[-grep("MTRNR2L1",genes)]
genes[grep("MSRB1",genes)] <- "SEPX1";genes <- genes[-grep("MT-ND1",genes)]
genes[grep("NDUFAF6",genes)] <- "C8orf38";genes[grep("CERS2",genes)] <- "LASS2"
genes[grep("COA6",genes)] <- "C1orf31";genes[grep("MZT2A",genes)] <- "FAM128A"
genes[grep("MPLKIP",genes)] <- "C7orf11";genes[grep("C7orf73",genes)] <- "LASS2"
genes[grep("MALSU1",genes)] <- "C7orf30";genes[grep("AK4",genes)] <- "AK3"
genes[grep("UQCC3",genes)] <- "C11orf83";genes[grep("CCDC167",genes)] <- "C6orf129"
genes[grep("NCBP2-AS2",genes)] <- "AK3";genes[grep("GLMP",genes)] <- "C1orf85"
genes[grep("COX14",genes)] <- "C12orf62";genes <- genes[-grep("MIR4458HG",genes)]
genes <- genes[-grep("SMLR1",genes)];genes[grep("SMCO4",genes)] <- "C11orf75"
genes[grep("WDR45B",genes)] <- "WDR45L";genes[grep("MSMO1",genes)] <- "SC4MOL"
genes <- genes[-grep("LINC01419",genes)];genes[grep("LAMTOR2",genes)] <- "ROBLD3"
genes[grep("NELFE",genes)] <- "RDBP";genes[grep("ECI2",genes)] <- "PECI"
genes[grep("ERO1A",genes)] <- "ERO1L";genes <- genes[-grep("SNHG19",genes)]
genes[grep("MRPL57",genes)] <- "MRP63";genes[grep("SWI5",genes)] <- "C9orf119"
genes[grep("HACD3",genes)] <- "PTPLAD1";genes <- genes[-grep("C8orf82",genes)]

expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)
pwd = "results/GMDSC/corr/TCGA/"

for(gene in colnames(expr.norm)[1:6]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

colnames(expr.norm) <- gsub("-",".",colnames(expr.norm))

for(i in 1:6){
  for(j in 7:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x = colnames(expr.norm)[i], y = colnames(expr.norm)[j])) +
      geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,colnames(expr.norm)[i],"/",colnames(expr.norm)[j],"_vs_",colnames(expr.norm)[i],".pdf"),p2,width = 6,height = 6)
  }
}


#>>> M-MDSC 

## ## 分析M-MDSC中高表达的基因，并分析它们的表达在数据集以及在TCGA中与CD33,CD11B,CD8的相关性。

markers <- xlsx::read.xlsx("results/MMDSC/markers.xlsx",sheetIndex = 1)
gene <- c("CD33","ITGAM","CD8A","CD8B",as.character(markers[,1]))
# singlecell
expr <- GetAssayData(scRNA,assay = "RNA",slot = "data")[gene,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(scRNA), to = scRNA$Sample)
colnames(expr.T) <- c(gene,"Sample")

AvgExpPerSample <- function(expMtx){
  patient <- NULL
  for(i in 1:(ncol(expMtx)-1)){
    patient[i] <- data.frame(tapply(expMtx[,i],expMtx$Sample,sum))
  }
  names(patient) <- colnames(expMtx)[-ncol(expMtx)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expMtx[,1],expMtx$Sample,sum))
  Patient_sample <- as.data.frame(table(expMtx$Sample),stringsAsFactors = F)
  patient$Sample <- plyr::mapvalues(rownames(patient),from = Patient_sample$Var1,to = Patient_sample$Freq)
  patient <- as.data.frame(patient,stringsAsFactors = F)
  patient$Sample <- as.numeric(patient$Sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$Sample[j]
    }
  }
  patient$Sample <- NULL
  return(patient)
}

patient <- AvgExpPerSample(expr.T)
colnames(patient) <- gene
xlsx::write.xlsx(patient,file = "results/MMDSC/avgExpPerPatient.xlsx")

colnames(patient) <- gsub("-",".",colnames(patient))
pwd = "results/MMDSC/corr/singlecell/"
for(gene in colnames(patient)[1:4]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

for(i in 1:4){
  for(j in 5:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[i], y = colnames(patient)[j])) +
      geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,colnames(patient)[i],"/",colnames(patient)[j],"_vs_",colnames(patient)[i],".pdf"),p2,width = 6,height = 6)
  }
}

# TCGA
load("../new/refdb/TCGA_LIHC.rnaseq.RData")

genesTBL <- xlsx::read.xlsx("results/MMDSC/avgExpPerPatient.xlsx",sheetIndex = 1)

genes <- gsub("\\.","-",colnames(genesTBL)[-1])
genes <- genes[-grep("-",genes)]

setdiff(genes, colnames(LIHC.rnaseq))

genes <- genes[-grep("IGLC3",genes)]; genes <- genes[-grep("IGLC2",genes)]
genes <- genes[-grep("IGKC",genes)]; genes <- genes[-grep("MTRNR2L8",genes)]
genes <- genes[-grep("MTRNR2L12",genes)]; genes[grep("SELENOP",genes)] <- "SEPP1"
genes[grep("MPC2",genes)] <- "BRP44"; genes <- genes[-grep("IGHM",genes)]
genes[grep("DNPH1",genes)] <- "C6orf108"; genes <- genes[-grep("SNHG25",genes)]
genes <- genes[-grep("IGHG1",genes)]; genes[grep("ECI2",genes)] <- "PECI"
genes[grep("COA3",genes)] <- "CCDC56"; genes <- genes[-grep("ANGPTL8",genes)]
genes <- genes[-grep("LINC01485",genes)];genes[grep("FUOM",genes)] <- "C10orf125"
genes[grep("AK4",genes)] <- "AK3"; genes[grep("CEBPZOS",genes)] <- "CEBPZ"
genes[grep("MSMO1",genes)] <- "SC4MOL"; genes[grep("UQCC2",genes)] <- "C6orf125"
genes[grep("EMC3",genes)] <- "TMEM111"; genes[grep("NAPRT",genes)] <- "NAPRT1"
genes[grep("CCDC167",genes)] <- "C6orf129"; genes[grep("MIEN1",genes)] <- "C17orf37"
genes[grep("LAMTOR2",genes)] <- "ROBLD3"; genes <- genes[-grep("MIR4458HG",genes)]
genes[grep("APMAP",genes)] <- "C20orf3"; genes[grep("METTL23",genes)] <- "C17orf95"
genes[grep("COPRS",genes)] <- "C17orf79"; genes <- genes[-grep("MTRNR2L1",genes)]
genes[grep("HACD3",genes)] <- "PTPLAD1"; genes <- genes[-grep("SMLR1",genes)]
genes[grep("PKM",genes)] <- "PKM2"; genes[grep("COA4",genes)] <- "CHCHD8"
genes[grep("COX20",genes)] <- "FAM36A"; genes[grep("MPLKIP",genes)] <- "C7orf11"
genes[grep("MRPL57",genes)] <- "MRP63"


expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)

pwd = "results/MMDSC/corr/TCGA/"

for(gene in colnames(expr.norm)[1:4]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

colnames(expr.norm) <- gsub("-",".",colnames(expr.norm))

for(i in 1:4){
  for(j in 5:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x = colnames(expr.norm)[i], y = colnames(expr.norm)[j])) +
      geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,colnames(expr.norm)[i],"/",colnames(expr.norm)[j],"_vs_",colnames(expr.norm)[i],".pdf"),p2,width = 6,height = 6)
  }
}

#>>> Tex CD8+T 

# 分析Tex CD8+T中高表达的基因，并分析它们的表达在我们的数据集以及在TCGA数据集中与
# CD8,T-bet,IFNg,GranzymeB,IL-2,PD-1,LAG3,TIM3,TOX的相关性

markers <- xlsx::read.xlsx("results/CD8_Tex/markers.xlsx",sheetIndex = 1)
gene <- c("CD8A","CD8B","TBX21","IFNG","GZMB","IL2","PDCD1","LAG3","HAVCR2","TOX",as.character(markers[,1]))
# singlecell
expr <- GetAssayData(scRNA,assay = "RNA",slot = "data")[gene,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(scRNA), to = scRNA$Sample)
colnames(expr.T) <- c(gene,"Sample")

AvgExpPerSample <- function(expMtx){
  patient <- NULL
  for(i in 1:(ncol(expMtx)-1)){
    patient[i] <- data.frame(tapply(expMtx[,i],expMtx$Sample,sum))
  }
  names(patient) <- colnames(expMtx)[-ncol(expMtx)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expMtx[,1],expMtx$Sample,sum))
  Patient_sample <- as.data.frame(table(expMtx$Sample),stringsAsFactors = F)
  patient$Sample <- plyr::mapvalues(rownames(patient),from = Patient_sample$Var1,to = Patient_sample$Freq)
  patient <- as.data.frame(patient,stringsAsFactors = F)
  patient$Sample <- as.numeric(patient$Sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$Sample[j]
    }
  }
  patient$Sample <- NULL
  return(patient)
}

patient <- AvgExpPerSample(expr.T)
colnames(patient) <- gene
xlsx::write.xlsx(patient,file = "results/CD8_Tex/avgExpPerPatient.xlsx")

colnames(patient) <- gsub("-",".",colnames(patient))
pwd = "results/CD8_Tex/corr/singlecell/"
for(gene in colnames(patient)[1:10]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

for(i in 1:10){
  for(j in 11:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[i], y = colnames(patient)[j])) +
      geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,colnames(patient)[i],"/",colnames(patient)[j],"_vs_",colnames(patient)[i],".pdf"),p2,width = 6,height = 6)
  }
}

# TCGA
load("../new/refdb/TCGA_LIHC.rnaseq.RData")

genesTBL <- xlsx::read.xlsx("results/CD8_Tex/avgExpPerPatient.xlsx",sheetIndex = 1)

genes <- gsub("\\.","-",colnames(genesTBL)[-1])
genes <- genes[-grep("-",genes)]

setdiff(genes, colnames(LIHC.rnaseq))

genes <- genes[-grep("IGLL5",genes)];genes <- genes[-grep("IGHG4",genes)]
genes <- genes[-grep("IGHG1",genes)];genes <- genes[-grep("MTRNR2L8",genes)]
genes <- genes[-grep("MTRNR2L12",genes)];genes <- genes[-grep("IGLC2",genes)]

genes[grep("MPC2",genes)] <- "BRP44";genes <- genes[-grep("IGLC3",genes)]
genes[grep("PKM",genes)] <- "PKM2";genes[grep("RACK1",genes)] <- "GNB2L1"
genes[grep("SELENOK",genes)] <- "SELK";genes <- genes[-grep("MTRNR2L1",genes)]
genes[grep("DNPH1",genes)] <- "C6orf108";genes <- genes[-grep("IGKC",genes)]
genes[grep("COPS9",genes)] <- "MYEOV2";genes <- genes[-grep("MZB1",genes)]
genes[grep("LAMTOR2",genes)] <- "ROBLD3";genes[grep("SELENOW",genes)] <- "SEPW1"
genes[grep("SELENOP",genes)] <- "SEPP1";genes[grep("PLPP1",genes)] <- "PPAP2A"
genes[grep("SELENOT",genes)] <- "SELT";genes[grep("CXCL8",genes)] <- "IL8"
genes[grep("SELENOH",genes)] <- "C11orf31";genes[grep("ERO1A",genes)] <- "ERO1L"
genes[grep("CNOT9",genes)] <- "RQCD1";genes[grep("SLC50A1",genes)] <- "RAG1AP1"
genes[grep("SELENOS",genes)] <- "SELS";genes[grep("SELENOF",genes)] <- "SEP15"
genes[grep("LRRC75A",genes)] <- "C17orf76";genes[grep("COA3",genes)] <- "CCDC56"
genes <- genes[-grep("HIST2H2AA4",genes)]; genes[grep("MINOS1",genes)] <- "C1orf151"
genes[grep("EMC4",genes)] <- "TMEM85";genes[grep("LAMTOR5",genes)] <- "HBXIP"
genes[grep("UQCC2",genes)] <- "C6orf125";genes[grep("COX20",genes)] <- "FAM36A"
genes <- genes[-grep("TMEM35B",genes)];genes[grep("SLIRP",genes)] <- "C14orf156"
genes[grep("MRNIP",genes)] <- "C5orf45";genes[grep("JCHAIN",genes)] <- "IGJ"
genes <- genes[-grep("BOLA2B",genes)];genes[grep("MYDGF",genes)] <- "IL27"
genes[grep("HYPK",genes)] <- "C15orf63";genes[grep("LBHD1",genes)] <- "C11orf48"
genes[grep("NSRP1",genes)] <- "CCDC55";genes <- genes[-grep("C11orf98",genes)]
genes <- genes[-grep("U2AF1L5",genes)]; genes[grep("SGO2",genes)] <- "SGOL2"
genes[grep("AKIP1",genes)] <- "C11orf17"; genes <- genes[-grep("NPIPB5",genes)]
genes[grep("COA6",genes)] <- "C1orf31";genes[grep("UTP11",genes)] <- "UTP11L"
genes[grep("CBWD7",genes)] <- "CBWD6";genes[grep("BBIP1",genes)] <- "NCRNA00081"
genes <- genes[-grep("NPIPB4",genes)]; genes[grep("FMC1",genes)] <- "C7orf55"
genes <- genes[-grep("LINC00998",genes)];genes[grep("TMEM258",genes)] <- "C11orf10"
genes[grep("CYSTM1",genes)] <- "C5orf32"

expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)
pwd = "results/CD8_Tex/corr/TCGA/"

for(gene in colnames(expr.norm)[1:10]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

colnames(expr.norm) <- gsub("-",".",colnames(expr.norm))

for(i in 1:10){
  for(j in 11:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x = colnames(expr.norm)[i], y = colnames(expr.norm)[j])) +
      geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,colnames(expr.norm)[i],"/",colnames(expr.norm)[j],"_vs_",colnames(expr.norm)[i],".pdf"),p2,width = 6,height = 6)
  }
}


#>>> Tex CD4+T

# 分析Tex CD4+T中高表达的基因，并分析它们的表达在我们的数据集以及在TCGA数据集中与
# CD4,T-bet,IFNg,GranzymeB,IL-2,PD-1,LAG3,TIM3,TOX的相关性

markers <- xlsx::read.xlsx("results/CD4_Tex/markers.xlsx",sheetIndex = 1)
gene <- c("CD4","TBX21","IFNG","GZMB","IL2","PDCD1","LAG3","HAVCR2","TOX",as.character(markers[,1]))
# singlecell
expr <- GetAssayData(scRNA,assay = "RNA",slot = "data")[gene,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(scRNA), to = scRNA$Sample)
colnames(expr.T) <- c(gene,"Sample")

AvgExpPerSample <- function(expMtx){
  patient <- NULL
  for(i in 1:(ncol(expMtx)-1)){
    patient[i] <- data.frame(tapply(expMtx[,i],expMtx$Sample,sum))
  }
  names(patient) <- colnames(expMtx)[-ncol(expMtx)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expMtx[,1],expMtx$Sample,sum))
  Patient_sample <- as.data.frame(table(expMtx$Sample),stringsAsFactors = F)
  patient$Sample <- plyr::mapvalues(rownames(patient),from = Patient_sample$Var1,to = Patient_sample$Freq)
  patient <- as.data.frame(patient,stringsAsFactors = F)
  patient$Sample <- as.numeric(patient$Sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$Sample[j]
    }
  }
  patient$Sample <- NULL
  return(patient)
}

patient <- AvgExpPerSample(expr.T)
colnames(patient) <- gene
xlsx::write.xlsx(patient,file = "results/CD4_Tex/avgExpPerPatient.xlsx")

colnames(patient) <- gsub("-",".",colnames(patient))
pwd = "results/CD4_Tex/corr/singlecell/"
for(gene in colnames(patient)[1:9]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

for(i in 1:9){
  for(j in 10:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[i], y = colnames(patient)[j])) +
      geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,colnames(patient)[i],"/",colnames(patient)[j],"_vs_",colnames(patient)[i],".pdf"),p2,width = 6,height = 6)
  }
}

# TCGA
load("../new/refdb/TCGA_LIHC.rnaseq.RData")

genesTBL <- xlsx::read.xlsx("results/CD4_Tex/avgExpPerPatient.xlsx",sheetIndex = 1)

genes <- gsub("\\.","-",colnames(genesTBL)[-1])
genes <- genes[-grep("-",genes)]

setdiff(genes, colnames(LIHC.rnaseq))

genes[grep("MPC2",genes)] <- "BRP44"
genes[grep("DNPH1",genes)] <- "C6orf108"
genes <- genes[-grep("MTRNR2L8",genes)]
genes <- genes[-grep("SNHG25",genes)]
genes <- genes[-grep("MTRNR2L12",genes)]
genes[grep("COA3",genes)] <- "CCDC56"
genes[grep("SLC50A1",genes)] <- "RAG1AP1"
genes[grep("MSMO1",genes)] <- "SC4MOL"
genes[grep("UQCC2",genes)] <- "C6orf125"
genes[grep("MPLKIP",genes)] <- "C7orf11"
genes[grep("SLIRP",genes)] <- "C14orf156"
genes[grep("FUOM",genes)] <- "C10orf125"
genes[grep("LAMTOR2",genes)] <- "ROBLD3"
genes[grep("AK4",genes)] <- "AK3"
genes[grep("PKM",genes)] <- "PKM2"
genes[grep("FAM213A",genes)] <- "C10orf58"
genes[grep("CYSTM1",genes)] <- "C5orf32"
genes[grep("CERS2",genes)] <- "LASS2"
genes[grep("ECI2",genes)] <- "PECI"
genes[grep("RACK1",genes)] <- "GNB2L1"
genes[grep("COA4",genes)] <- "CHCHD8"
genes <- genes[-grep("SF3B6",genes)]
genes[grep("MIEN1",genes)] <- "C17orf37"
genes[grep("SELENOP",genes)] <- "SEPP1"
genes[grep("COA6",genes)] <- "C1orf31"
genes[grep("HACD3",genes)] <- "PTPLAD1"
genes[grep("LAMTOR4",genes)] <- "C7orf59"
genes <- genes[-grep("MIR4458HG",genes)]
genes[grep("MRPL57",genes)] <- "MRP63"
genes[grep("ERO1A",genes)] <- "ERO1L"
genes <- genes[-grep("C7orf73",genes)]
genes[grep("NAPRT",genes)] <- "NAPRT1"
genes[grep("VMP1",genes)] <- "TMEM49"
genes[grep("UQCC3",genes)] <- "C11orf83"
genes[grep("GLMP",genes)] <- "C1orf85"
genes <- genes[-grep("LINC00493",genes)]
genes[grep("COPRS",genes)] <- "C17orf79"
genes <- genes[-grep("CHP1",genes)]
genes[grep("ATRAID",genes)] <- "C2orf28"
genes[grep("SELENOK",genes)] <- "SELK"
genes[grep("COX20",genes)] <- "FAM36A"
genes[grep("MZT2B",genes)] <- "FAM128B"
genes[grep("CXCL8",genes)] <- "IL8"
genes <- genes[-grep("MTRNR2L1",genes)]
genes[grep("LAMTOR5",genes)] <- "HBXIP"

genes[grep("CTSL",genes)] <- "CTSL1"
genes[grep("ADGRG6",genes)] <- "GPR126"
genes <- genes[-grep("LINC00116",genes)]


expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)
pwd = "results/CD4_Tex/corr/TCGA/"

for(gene in colnames(expr.norm)[1:9]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

colnames(expr.norm) <- gsub("-",".",colnames(expr.norm))

for(i in 1:9){
  for(j in 10:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x = colnames(expr.norm)[i], y = colnames(expr.norm)[j])) +
      geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,colnames(expr.norm)[i],"/",colnames(expr.norm)[j],"_vs_",colnames(expr.norm)[i],".pdf"),p2,width = 6,height = 6)
  }
}


#>>> CTL

# 分析CTL中高表达的基因，并分析它们的表达在我们的数据集以及在TCGA数据集中与
# CD45 CD3 CD8 CD107a (LAMP-1) GranzymeB IFN-γ的相关性

markers <- xlsx::read.xlsx("results/CTL/markers.xlsx",sheetIndex = 1)
gene <- c("PTPRC","CD3D","CD3E","CD3G","CD8A","CD8B","LAMP1","GZMB","IFNG",as.character(markers[,1]))
# singlecell
expr <- GetAssayData(scRNA,assay = "RNA",slot = "data")[gene,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(scRNA), to = scRNA$Sample)
colnames(expr.T) <- c(gene,"Sample")

AvgExpPerSample <- function(expMtx){
  patient <- NULL
  for(i in 1:(ncol(expMtx)-1)){
    patient[i] <- data.frame(tapply(expMtx[,i],expMtx$Sample,sum))
  }
  names(patient) <- colnames(expMtx)[-ncol(expMtx)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expMtx[,1],expMtx$Sample,sum))
  Patient_sample <- as.data.frame(table(expMtx$Sample),stringsAsFactors = F)
  patient$Sample <- plyr::mapvalues(rownames(patient),from = Patient_sample$Var1,to = Patient_sample$Freq)
  patient <- as.data.frame(patient,stringsAsFactors = F)
  patient$Sample <- as.numeric(patient$Sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$Sample[j]
    }
  }
  patient$Sample <- NULL
  return(patient)
}

patient <- AvgExpPerSample(expr.T)
colnames(patient) <- gene
xlsx::write.xlsx(patient,file = "results/CTL/avgExpPerPatient.xlsx")

colnames(patient) <- gsub("-",".",colnames(patient))
pwd = "results/CTL/corr/singlecell/"
for(gene in colnames(patient)[1:9]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

for(i in 1:9){
  for(j in 10:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[i], y = colnames(patient)[j])) +
      geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,colnames(patient)[i],"/",colnames(patient)[j],"_vs_",colnames(patient)[i],".pdf"),p2,width = 6,height = 6)
  }
}

# TCGA
load("../new/refdb/TCGA_LIHC.rnaseq.RData")

genesTBL <- xlsx::read.xlsx("results/CTL/avgExpPerPatient.xlsx",sheetIndex = 1)

genes <- gsub("\\.","-",colnames(genesTBL)[-1])
genes <- genes[-grep("-",genes)]

setdiff(genes, colnames(LIHC.rnaseq))

genes <- genes[-grep("IGHG4",genes)];genes <- genes[-grep("MTRNR2L8",genes)]
genes <- genes[-grep("IGHG1",genes)];genes[grep("MPC2",genes)] <- "BRP44"
genes <- genes[-grep("MTRNR2L12",genes)];genes <- genes[-grep("IGLC2",genes)]
genes <- genes[-grep("IGLL5",genes)];genes <- genes[-grep("IGLC3",genes)]
genes[grep("DNPH1",genes)] <- "C6orf108";genes[grep("PKM",genes)] <- "PKM2"
genes <- genes[-grep("IGKC",genes)];genes[grep("COA3",genes)] <- "CCDC56"
genes <- genes[-grep("MZB1",genes)];genes[grep("RACK1",genes)] <- "GNB2L1"
genes[grep("SELENOP",genes)] <- "SEPP1";genes[grep("SLC50A1",genes)] <- "RAG1AP1"
genes[grep("SELENOK",genes)] <- "SELK";genes[grep("UQCC2",genes)] <- "C6orf125"
genes[grep("COPS9",genes)] <- "MYEOV2";genes[grep("CYSTM1",genes)] <- "C5orf32"
genes[grep("MYDGF",genes)] <- "IL27";genes <- genes[-grep("SNHG25",genes)]
genes[grep("ERO1A",genes)] <- "ERO1L";genes[grep("SLIRP",genes)] <- "C14orf156"
genes[grep("PLPP1",genes)] <- "PPAP2A";genes[grep("LAMTOR5",genes)] <- "HBXIP"
genes[grep("SELENOW",genes)] <- "SEPW1"
genes[grep("FAM213A",genes)] <- "C10orf58";genes[grep("SELENOH",genes)] <- "C11orf31"
genes[grep("MSMO1",genes)] <- "SC4MOL";genes[grep("CERS2",genes)] <- "LASS2"
genes[grep("SELENOS",genes)] <- "SELS";genes[grep("AK4",genes)] <- "AK3"
genes[grep("MIEN1",genes)] <- "C17orf37";genes[grep("SELENOT",genes)] <- "SELT"
genes[grep("TMEM258",genes)] <- "C11orf10";genes[grep("MPLKIP",genes)] <- "C7orf11"
genes[grep("MZT2B",genes)] <- "FAM128B";genes[grep("FUOM",genes)] <- "C10orf125"
genes[grep("COA6",genes)] <- "C1orf31";genes[grep("EMC4",genes)] <- "TMEM85"
genes[grep("CCDC167",genes)] <- "C6orf129";genes[grep("CNOT9",genes)] <- "RQCD1"
genes[grep("COX20",genes)] <- "FAM36A";genes[grep("ECI2",genes)] <- "PECI"
genes[grep("MINOS1",genes)] <- "C1orf151";genes[grep("COA4",genes)] <- "CHCHD8"
genes[grep("SELENOF",genes)] <- "SEP15";genes[grep("COPRS",genes)] <- "C17orf79"


expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)
pwd = "results/CTL/corr/TCGA/"

for(gene in colnames(expr.norm)[1:9]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

colnames(expr.norm) <- gsub("-",".",colnames(expr.norm))

for(i in 1:9){
  for(j in 10:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x = colnames(expr.norm)[i], y = colnames(expr.norm)[j])) +
      geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,colnames(expr.norm)[i],"/",colnames(expr.norm)[j],"_vs_",colnames(expr.norm)[i],".pdf"),p2,width = 6,height = 6)
  }
}


#>>> Th1

# 分析Th1中高表达的基因，并分析它们的表达在我们的数据集以及在TCGA数据集中与
# CD45 CD3 CD4 Tbet IFN-γ CXCR3的相关性

markers <- xlsx::read.xlsx("results/Th1/markers.xlsx",sheetIndex = 1)
gene <- c("PTPRC","CD3D","CD3E","CD3G","CD4","TBX21","IFNG","CXCR3",as.character(markers[,1]))
# singlecell
expr <- GetAssayData(scRNA,assay = "RNA",slot = "data")[gene,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(scRNA), to = scRNA$Sample)
colnames(expr.T) <- c(gene,"Sample")

AvgExpPerSample <- function(expMtx){
  patient <- NULL
  for(i in 1:(ncol(expMtx)-1)){
    patient[i] <- data.frame(tapply(expMtx[,i],expMtx$Sample,sum))
  }
  names(patient) <- colnames(expMtx)[-ncol(expMtx)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expMtx[,1],expMtx$Sample,sum))
  Patient_sample <- as.data.frame(table(expMtx$Sample),stringsAsFactors = F)
  patient$Sample <- plyr::mapvalues(rownames(patient),from = Patient_sample$Var1,to = Patient_sample$Freq)
  patient <- as.data.frame(patient,stringsAsFactors = F)
  patient$Sample <- as.numeric(patient$Sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$Sample[j]
    }
  }
  patient$Sample <- NULL
  return(patient)
}

patient <- AvgExpPerSample(expr.T)
colnames(patient) <- gene
xlsx::write.xlsx(patient,file = "results/Th1/avgExpPerPatient.xlsx")

colnames(patient) <- gsub("-",".",colnames(patient))
pwd = "results/Th1/corr/singlecell/"
for(gene in colnames(patient)[1:8]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

for(i in 1:8){
  for(j in 9:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[i], y = colnames(patient)[j])) +
      geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,colnames(patient)[i],"/",colnames(patient)[j],"_vs_",colnames(patient)[i],".pdf"),p2,width = 6,height = 6)
  }
}

# TCGA
load("../new/refdb/TCGA_LIHC.rnaseq.RData")

genesTBL <- xlsx::read.xlsx("results/Th1/avgExpPerPatient.xlsx",sheetIndex = 1)

genes <- gsub("\\.","-",colnames(genesTBL)[-1])
genes <- genes[-grep("-",genes)]

setdiff(genes, colnames(LIHC.rnaseq))

genes <- genes[-grep("IGHG4",genes)];genes <- genes[-grep("MTRNR2L8",genes)]
genes <- genes[-grep("IGHG1",genes)];genes[grep("MPC2",genes)] <- "BRP44"
genes <- genes[-grep("MTRNR2L12",genes)];genes <- genes[-grep("IGLC2",genes)]
genes <- genes[-grep("IGLL5",genes)];genes <- genes[-grep("IGLC3",genes)]
genes[grep("DNPH1",genes)] <- "C6orf108";genes[grep("PKM",genes)] <- "PKM2"
genes <- genes[-grep("IGKC",genes)];genes[grep("COA3",genes)] <- "CCDC56"
genes <- genes[-grep("MZB1",genes)];genes[grep("RACK1",genes)] <- "GNB2L1"
genes[grep("SELENOP",genes)] <- "SEPP1";genes[grep("SLC50A1",genes)] <- "RAG1AP1"
genes[grep("SELENOK",genes)] <- "SELK";genes[grep("UQCC2",genes)] <- "C6orf125"
genes[grep("COPS9",genes)] <- "MYEOV2";genes[grep("CYSTM1",genes)] <- "C5orf32"
genes[grep("MYDGF",genes)] <- "IL27";genes <- genes[-grep("SNHG25",genes)]
genes[grep("ERO1A",genes)] <- "ERO1L";genes[grep("SLIRP",genes)] <- "C14orf156"
genes[grep("PLPP1",genes)] <- "PPAP2A";genes[grep("LAMTOR5",genes)] <- "HBXIP"
genes[grep("SELENOW",genes)] <- "SEPW1"
genes[grep("FAM213A",genes)] <- "C10orf58";genes[grep("SELENOH",genes)] <- "C11orf31"
genes[grep("MSMO1",genes)] <- "SC4MOL";genes[grep("CERS2",genes)] <- "LASS2"
genes[grep("SELENOS",genes)] <- "SELS";genes[grep("AK4",genes)] <- "AK3"
genes[grep("MIEN1",genes)] <- "C17orf37";genes[grep("SELENOT",genes)] <- "SELT"
genes[grep("TMEM258",genes)] <- "C11orf10";genes[grep("MPLKIP",genes)] <- "C7orf11"
genes[grep("MZT2B",genes)] <- "FAM128B";genes[grep("FUOM",genes)] <- "C10orf125"
genes[grep("COA6",genes)] <- "C1orf31";genes[grep("EMC4",genes)] <- "TMEM85"
genes[grep("CCDC167",genes)] <- "C6orf129";genes[grep("CNOT9",genes)] <- "RQCD1"
genes[grep("COX20",genes)] <- "FAM36A";genes[grep("ECI2",genes)] <- "PECI"
genes[grep("MINOS1",genes)] <- "C1orf151";genes[grep("COA4",genes)] <- "CHCHD8"
genes[grep("SELENOF",genes)] <- "SEP15";genes[grep("COPRS",genes)] <- "C17orf79"
genes[grep("LAMTOR2",genes)] <- "ROBLD3";genes <- genes[-grep("MTRNR2L1",genes)]

expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)
pwd = "results/Th1/corr/TCGA/"

for(gene in colnames(expr.norm)[1:8]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

colnames(expr.norm) <- gsub("-",".",colnames(expr.norm))

for(i in 1:8){
  for(j in 9:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x = colnames(expr.norm)[i], y = colnames(expr.norm)[j])) +
      geom_point(size = 1.5, color = '#F9B208',alpha=.7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,colnames(expr.norm)[i],"/",colnames(expr.norm)[j],"_vs_",colnames(expr.norm)[i],".pdf"),p2,width = 6,height = 6)
  }
}
