rm(list = ls());gc()
library(Seurat)
library(ggplot2)
library(dplyr)

scRNA_integrated_without_PVTT_MLN <- readRDS("rds/without_PVTT_MLN/total/1.scRNA_integrated_without_PVTT_MLN.rds")

scRNA_integrated_without_PVTT_MLN$Tissue <- ifelse(
  scRNA_integrated_without_PVTT_MLN$Sample=="HCC03N" | scRNA_integrated_without_PVTT_MLN$Sample=="HCC04N" |
    scRNA_integrated_without_PVTT_MLN$Sample=="HCC05N" | scRNA_integrated_without_PVTT_MLN$Sample=="HCC06N" |
    scRNA_integrated_without_PVTT_MLN$Sample=="HCC07N" | scRNA_integrated_without_PVTT_MLN$Sample=="HCC08N" |
    scRNA_integrated_without_PVTT_MLN$Sample=="HCC09N" | scRNA_integrated_without_PVTT_MLN$Sample=="HCC10N",
  "normal", "tumor"
)

tumor_cells <- subset(scRNA_integrated_without_PVTT_MLN@meta.data, Tissue == "tumor")
tumor <- subset(scRNA_integrated_without_PVTT_MLN, cells = rownames(tumor_cells))

DefaultAssay(tumor) <- "RNA"
# 降维聚类
tumor <- FindVariableFeatures(tumor,selection.method = "vst",nfeatures = 2000)
DefaultAssay(tumor) <- "integrated"
tumor <- ScaleData(tumor,features = rownames(tumor))
tumor <- RunPCA(tumor, features = VariableFeatures(tumor))
p1 <- DimPlot(tumor, reduction = "pca", group.by = "orig.ident",cols=rainbow(length(levels(as.factor(tumor$orig.ident)))))
ggsave("result/without_PVTT_MLN/tumor/1.tumor_beforeQC_pca.pdf",p1,device="pdf",width=8,height=6)
p2 <- ElbowPlot(tumor,ndims=30,reduction="pca")
ggsave("result/without_PVTT_MLN/tumor/2.tumor_beforeQC_elbow.pdf",p2,device="pdf",width=6,height=4)

pc.num = 1:24
# 聚类
tumor <- FindNeighbors(tumor, dims = pc.num) %>% FindClusters(resolution = 0.2)
tumor <- RunUMAP(tumor,reduction="pca",dims=pc.num) %>% RunTSNE(dims=pc.num)

p1 <- DimPlot(tumor,reduction="tsne",cols=rainbow(length(levels(as.factor(tumor$seurat_clusters)))))
p2 <- DimPlot(tumor,reduction="tsne",group.by="orig.ident",cols=rainbow(length(levels(as.factor(tumor$orig.ident))))) + theme(plot.title = element_blank())
pc <- p1 + p2
ggsave("result/without_PVTT_MLN/tumor/3.tumor_beforeQC_tsne.pdf",pc,device = "pdf",width = 14,height = 4.5)

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
ggsave("result/without_PVTT_MLN/tumor/4.tumor_vlnplot_beforeQC.pdf",violin,device = "pdf",width = 8,height = 6)

# QC
minGene=200
maxGene=8000
pctMT=25

tumor <- subset(tumor, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT)

violin <- VlnPlot(tumor, group.by = "project",
                  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                  pt.size = 0.00001, ncol = 3) + NoLegend()

ggsave("result/without_PVTT_MLN/tumor/4.tumor_vlnplot_afterQC.pdf",violin,device = "pdf",width = 8,height = 6)

p1 <- DimPlot(tumor,reduction="tsne",cols=rainbow(length(levels(as.factor(tumor$seurat_clusters)))))
p2 <- DimPlot(tumor,reduction="tsne",group.by="orig.ident",cols=rainbow(length(levels(as.factor(tumor$orig.ident))))) + theme(plot.title=element_blank())
pc <- p1 + p2
ggsave("result/without_PVTT_MLN/tumor/3.tumor_afterQC_tsne.pdf",pc,device = "pdf",width = 14,height = 4.5)
saveRDS(tumor,file = "rds/without_PVTT_MLN/tumor/1.tumor.rds")

Idents(tumor) <- tumor$Sample
DefaultAssay(tumor) <- "RNA"
table(tumor$seurat_clusters,tumor$Type)

# cell type assignment (preset)
current_cluster <- c(0:21)
new_cluster <- c("TAM","T cell","Malignant cell","Malignant cell", "Malignant cell",
                 "Malignant cell","TEC","Malignant cell","CAF","Malignant cell",
                 "Malignant cell","B cell","T cell","Malignant cell","TAM",
                 "T cell","unclassified","B cell","TEC","CAF","unclassified","unclassified")
tumor$CellType <- plyr::mapvalues(tumor$seurat_clusters,from = current_cluster,to = new_cluster)
saveRDS(tumor,file = "rds/without_PVTT_MLN/tumor/2.tumor_celltypeAssigned.rds")

p1 <- DimPlot(tumor,reduction = "tsne",group.by = "seurat_clusters",cols = rainbow(length(levels(as.factor(tumor$seurat_clusters)))))
p2 <- DimPlot(tumor,reduction = "tsne",group.by = "Sample",cols = rainbow(length(levels(as.factor(tumor$Sample)))))
p3 <- DimPlot(tumor,reduction = "tsne",group.by = "CellType",cols = rainbow(length(levels(as.factor(tumor$CellType)))))
pc <- p1 + p2 + p3
ggsave("result/without_PVTT_MLN/tumor/5.tumor_tsne.pdf",pc,device = "pdf",width = 16,height = 4)

# CD3: CD3D, CD3E, CD3G
# CD45: PTPRC
# CD8: CD8A, CD8B
# CD11b: ITGAM
# B220: PTPRC
# CD11C: ITGAX
# NK1.1: KLRB1
# gdT: CD27

markers <- c("CD3D","CD3E","CD3G","PTPRC","CD8A","CD8B","ITGAM","ITGAX","KLRB1","CD27")
Tcell <- DimPlot(tumor,group.by = "CellType",reduction = "tsne",cells.highlight = rownames(subset(tumor@meta.data,CellType=="T cell"))) + ggtitle("T cell") +NoLegend()
TAM <- DimPlot(tumor,group.by = "CellType",reduction = "tsne",cells.highlight = rownames(subset(tumor@meta.data,CellType=="TAM"))) + ggtitle("TAM") +NoLegend()
TEC <- DimPlot(tumor,group.by = "CellType",reduction = "tsne",cells.highlight = rownames(subset(tumor@meta.data,CellType=="TEC"))) + ggtitle("TEC") +NoLegend()
CAF <- DimPlot(tumor,group.by = "CellType",reduction = "tsne",cells.highlight = rownames(subset(tumor@meta.data,CellType=="CAF"))) + ggtitle("CAF") +NoLegend()
Bcell <- DimPlot(tumor,group.by = "CellType",reduction = "tsne",cells.highlight = rownames(subset(tumor@meta.data,CellType=="B cell"))) + ggtitle("B cell") +NoLegend()
MC <- DimPlot(tumor,group.by = "CellType",reduction = "tsne",cells.highlight = rownames(subset(tumor@meta.data,CellType=="Malignant cell"))) + ggtitle("Malignant cell") +NoLegend()
unclassified <- DimPlot(tumor,group.by = "CellType",reduction = "tsne",cells.highlight = rownames(subset(tumor@meta.data,CellType=="unclassified"))) + ggtitle("unclassified") +NoLegend()
ggsave("result/without_PVTT_MLN/tumor/PARTI/preprocessing/Tcell.highlight.pdf",Tcell,width = 6,height = 5.5)
ggsave("result/without_PVTT_MLN/tumor/PARTI/preprocessing/TAM.highlight.pdf",TAM,width = 6,height = 5.5)
ggsave("result/without_PVTT_MLN/tumor/PARTI/preprocessing/TEC.highlight.pdf",TEC,width = 6,height = 5.5)
ggsave("result/without_PVTT_MLN/tumor/PARTI/preprocessing/Bcell.highlight.pdf",Bcell,width = 6,height = 5.5)
ggsave("result/without_PVTT_MLN/tumor/PARTI/preprocessing/CAF.highlight.pdf",CAF,width = 6,height = 5.5)
ggsave("result/without_PVTT_MLN/tumor/PARTI/preprocessing/Malignant cell.highlight.pdf",MC,width = 6,height = 5.5)
ggsave("result/without_PVTT_MLN/tumor/PARTI/preprocessing/unclassified.highlight.pdf",unclassified,width = 6,height = 5.5)


##====================== T subtype ============================================================

Tsub_cells <- subset(tumor@meta.data, CellType == "T cell")
Tsub <- subset(tumor, cells = rownames(Tsub_cells))

Tsub <- FindVariableFeatures(Tsub,selection.method = "vst",nfeatures = 2000)
DefaultAssay(Tsub) <- "integrated"
Tsub <- ScaleData(Tsub,features = rownames(Tsub))
Tsub <- RunPCA(Tsub, features = VariableFeatures(Tsub))
p1 <- DimPlot(Tsub, reduction = "pca", group.by = "orig.ident",cols=rainbow(length(levels(as.factor(Tsub$orig.ident)))))
ggsave("result/without_PVTT_MLN/tumor/sub_T/1.Tsub_pca.pdf",p1,device="pdf",width=8,height=6)
p2 <- ElbowPlot(Tsub,ndims=30,reduction="pca")
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.Tsub_elbow.pdf",p2,device="pdf",width=6,height=4)

pc.num = 1:24
# 聚类
Tsub <- FindNeighbors(Tsub, dims = pc.num) %>% FindClusters(resolution = 0.2)
Tsub <- RunUMAP(Tsub,reduction="pca",dims=pc.num) %>% RunTSNE(dims=pc.num)

pal <- c("#786140","#8C0603","#B2712B","#B48521","#C6EEC4","#DCD3B1","#A7BBA5","#0202BB","#5EC6B9","#4E8A97")

p1 <- DimPlot(Tsub,reduction="tsne",cols=pal)
p2 <- DimPlot(Tsub,reduction="tsne",group.by="orig.ident",cols=rainbow(length(levels(as.factor(Tsub$orig.ident))))) + theme(plot.title = element_blank())
pc <- p1 + p2
ggsave("result/without_PVTT_MLN/tumor/sub_T/3.Tsub_tsne.pdf",pc,device = "pdf",width = 14,height = 4.5)

DefaultAssay(Tsub) <- "RNA"
# percent mt ribo
Tsub[['percent.mt']] <- PercentageFeatureSet(Tsub, pattern="^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(Tsub@assays$RNA)) 
HB.genes <- rownames(Tsub@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)]
Tsub[['percent.hb']] <- PercentageFeatureSet(Tsub, features=HB.genes)
Tsub$project <- "tumor"

violin <- VlnPlot(Tsub, group.by = "project",
                  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                  pt.size = 0.00001, ncol = 3) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/sub_T/4.Tsub_violin.pdf",violin,device = "pdf",width = 8,height = 6)


Idents(Tsub) <- Tsub$Sample
DefaultAssay(Tsub) <- "RNA"
table(Tsub$seurat_clusters,Tsub$Type)

CD3 <- VlnPlot(Tsub, features = "PTPRC",cols = pal,pt.size = 0) +
  coord_flip() + ggtitle("CD3") +
  theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
        panel.border = element_rect(size = 1.5,colour = "black"))
ggsave("result/without_PVTT_MLN/tumor/sub_T/5.Tsub_CD3T.pdf",width = 4,height = 6)

# cell type assignment

current_cluster <- c(0:9)
CD8T <- VlnPlot(Tsub,features = c("CD8A","CD8B"),pt.size = 0,ncol = 1,cols = pal)
ggsave("result/without_PVTT_MLN/tumor/sub_T/6.Tsub_CD8T.pdf",width = 8,height = 4)

# classfied CD3+T cell into CD4+T and CD8+T cells
new_cluster <- c("CD4+T cell","CD4+T cell","CD4+T cell","CD8+T cell","CD4+T cell","CD4+T cell",
                 "CD8+T cell","CD4+T cell","CD4+T cell","CD4+T cell","CD4+T cell","CD4+T cell")
Tsub$subtype <- plyr::mapvalues(Tsub$seurat_clusters,from = current_cluster,to = new_cluster)
saveRDS(Tsub,file = "rds/without_PVTT_MLN/tumor/sub_T/Tsub_celltypeAssigned.rds")

p1 <- DimPlot(Tsub,reduction = "tsne",group.by = "seurat_clusters",cols = pal)
p2 <- DimPlot(Tsub,reduction = "tsne",group.by = "Sample",cols = rainbow(length(levels(as.factor(Tsub$Sample)))))
p3 <- DimPlot(Tsub,reduction = "tsne",group.by = "subtype")
pc <- p1 + p2 + p3
ggsave("result/without_PVTT_MLN/tumor/sub_T/7.Tsub_tsne.pdf",pc,device = "pdf",width = 16,height = 4)

CD8T <- DimPlot(Tsub,reduction = "tsne",group.by = "subtype",cells.highlight = rownames(Tsub@meta.data[Tsub$subtype=="CD8+ T cell",])) + 
  NoLegend() + ggtitle("CD8+T cell") + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/sub_T/8.Tsub_CD8T_tsne.pdf",CD8T,device = "pdf",width = 4,height = 4)

CD4T <- DimPlot(Tsub,reduction = "tsne",group.by = "subtype",cells.highlight = rownames(Tsub@meta.data[Tsub$subtype=="CD4+ T cell",])) + 
  NoLegend() + ggtitle("CD4+T cell") + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/sub_T/9.Tsub_CD4T_tsne.pdf",CD4T,device = "pdf",width = 4,height = 4)

###### Calculating Infiltration Rate
tumor <- readRDS("rds/without_PVTT_MLN/tumor/2.tumor_celltypeAssigned.rds")
Tsub <- readRDS("rds/without_PVTT_MLN/tumor/sub_T/Tsub_celltypeAssigned.rds")

# T cell infiltration rate
total <- cbind(table(tumor$Sample,tumor$CellType),table(Tsub$Sample,Tsub$subtype))
colnames(total)[2] <- "CD3+ T cell"
total <- t(total)
CD3T.IR <- NULL
CD4T.IR <- NULL
CD8T.IR <- NULL
## CD3+T
for(i in seq_len(ncol(total))){CD3T.IR[i] <- total[2,i]/colSums(total[1:7,])[i]}
## CD4+T
for(i in seq_len(ncol(total))){CD4T.IR[i] <- total[8,i]/colSums(total[1:7,])[i]}
## CD8+T
for(i in seq_len(ncol(total))){CD8T.IR[i] <- total[9,i]/colSums(total[1:7,])[i]}
total <- rbind(total,CD3T.IR,CD4T.IR,CD8T.IR)
Sys.setenv("JAVA_HOME" = "C:/Program Files/Java/jre1.8.0_251")
xlsx::write.xlsx(total,file = "result/without_PVTT_MLN/tumor/sub_T/Infiltration_rate.xlsx")


# CD3+: CD3D, CD3E, CD3G
# CD45+: PTPRC
# CD8+: CD8A, CD8B
# CD11b-: ITGAM
# B220: PTPRC
# CD11C-: ITGAX
# NK1.1-: KLRB1
# gdT-: CD27

T.IR <- data.frame(t(total[10:12,]))
T.IR$Patient <- rownames(T.IR)
save(T.IR,file = "rds/without_PVTT_MLN/tumor/sub_T/T.IR.RData")

# add infiltration rate to tumor metadata
tumor$CD3T.IR <- plyr::mapvalues(tumor$Sample, from = T.IR$Patient, to = T.IR$CD3T.IR)
tumor$CD4T.IR <- plyr::mapvalues(tumor$Sample, from = T.IR$Patient, to = T.IR$CD4T.IR)
tumor$CD8T.IR <- plyr::mapvalues(tumor$Sample, from = T.IR$Patient, to = T.IR$CD8T.IR)
# saveRDS(tumor,file = "rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR.rds")


library(ggpubr)
T.IR <- T.IR[order(T.IR[,1],decreasing = T),]
CD3T_IR <- ggbarplot(T.IR,x = "Patient", y = "CD3T.IR",fill = "Patient") + 
  NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 14)) + 
  geom_hline(aes(yintercept=mean(T.IR$CD3T.IR)),color="#BF0032",linetype="dashed") +
  geom_hline(aes(yintercept=median(T.IR$CD3T.IR)),color="#0066A5",linetype="dashed") +
  geom_hline(aes(yintercept=quantile(T.IR$CD3T.IR,0.75)),color="#222222",linetype="dashed") +
  geom_hline(aes(yintercept=quantile(T.IR$CD3T.IR,0.25)),color="#222222",linetype="dashed") +
  annotate("text",x=20.5,y=mean(T.IR$CD3T.IR)+0.02,label="Mean:0.24",color="#BF0032",size=3) +
  annotate("text",x=20.5,y=median(T.IR$CD3T.IR)-0.03,label="Median:0.21",color="#0066A5",size=3) +
  annotate("text",x=20,y=quantile(T.IR$CD3T.IR,0.75)+0.03,label="Upper Quantile:0.30",color="#222222",size=3) +
  annotate("text",x=20,y=quantile(T.IR$CD3T.IR,0.25)+0.03,label="Lower Quantile:0.03",color="#222222",size=3)
ggsave("result/without_PVTT_MLN/tumor/sub_T/CD3T_IR.pdf",CD3T_IR,width = 8,height = 4)

T.IR <- T.IR[order(T.IR[,2],decreasing = T),]
CD4T_IR <- ggbarplot(T.IR,x = "Patient", y = "CD4T.IR",fill = "Patient") + 
  NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 14)) + 
  geom_hline(aes(yintercept=mean(T.IR$CD4T.IR)),color="#BF0032",linetype="dashed") +
  geom_hline(aes(yintercept=median(T.IR$CD4T.IR)),color="#0066A5",linetype="dashed") +
  geom_hline(aes(yintercept=quantile(T.IR$CD4T.IR,0.75)),color="#222222",linetype="dashed") +
  geom_hline(aes(yintercept=quantile(T.IR$CD4T.IR,0.25)),color="#222222",linetype="dashed") +
  annotate("text",x=20.5,y=mean(T.IR$CD4T.IR)+0.02,label="Mean:0.20",color="#BF0032",size=3) +
  annotate("text",x=20.5,y=median(T.IR$CD4T.IR)-0.03,label="Median:0.17",color="#0066A5",size=3) +
  annotate("text",x=20,y=quantile(T.IR$CD4T.IR,0.75)+0.03,label="Upper Quantile:0.24",color="#222222",size=3) +
  annotate("text",x=20,y=quantile(T.IR$CD4T.IR,0.25)+0.05,label="Lower Quantile:0.01",color="#222222",size=3)

ggsave("result/without_PVTT_MLN/tumor/sub_T/CD4T_IR.pdf",CD4T_IR,width = 8,height = 4)

T.IR <- T.IR[order(T.IR[,3],decreasing = T),]

CD8T_IR <- ggbarplot(T.IR,x = "Patient", y = "CD8T.IR",fill = "Patient") + 
  NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 14)) + 
  geom_hline(aes(yintercept=mean(T.IR$CD8T.IR)),color="#BF0032",linetype="dashed") +
  geom_hline(aes(yintercept=median(T.IR$CD8T.IR)),color="#0066A5",linetype="dashed") +
  geom_hline(aes(yintercept=quantile(T.IR$CD8T.IR,0.75)),color="#222222",linetype="dashed") +
  geom_hline(aes(yintercept=quantile(T.IR$CD8T.IR,0.25)),color="#222222",linetype="dashed") +
  annotate("text",x=20.5,y=mean(T.IR$CD8T.IR)-0.005,label="Mean:0.05",color="#BF0032",size=3) +
  annotate("text",x=20.5,y=median(T.IR$CD8T.IR)+0.006,label="Median:0.02",color="#0066A5",size=3) +
  annotate("text",x=20,y=quantile(T.IR$CD8T.IR,0.75)+0.01,label="Upper Quantile:0.06",color="#222222",size=3) +
  annotate("text",x=20,y=quantile(T.IR$CD8T.IR,0.25)+0.005,label="Lower Quantile:0.01",color="#222222",size=3)
ggsave("result/without_PVTT_MLN/tumor/sub_T/CD8T_IR.pdf",CD8T_IR,width = 8,height = 4)


tumor$CellType <- as.character(tumor$CellType)

# 设定阈值，区分High和low组
if(F){
  CD3_IR_threshold = 0.15
  CD4_IR_threshold = 0.10
  CD8_IR_threshold = 0.02
}

T.IR <- as.data.frame(readxl::read_excel("IR.xlsx",col_names = TRUE),stringAsFactors=F)

# add infiltration rate to tumor metadata
tumor$CD3T.IR <- plyr::mapvalues(tumor$Sample, from = T.IR$Patient, to = T.IR$CD3_IR)
tumor$CD4T.IR <- plyr::mapvalues(tumor$Sample, from = T.IR$Patient, to = T.IR$CD4.IR)
tumor$CD8T.IR <- plyr::mapvalues(tumor$Sample, from = T.IR$Patient, to = T.IR$CD8.IR)

CD3_IR_threshold = 0.014
CD4_IR_threshold = 0.012
CD8_IR_threshold = 0.0013

# 高亮阈值
T.IR <- T.IR[order(T.IR[,2],decreasing = T),]
CD3T_IR <- ggbarplot(T.IR, x = "Patient", y = "CD3_IR", fill = "Patient") +
  NoLegend() +
  geom_hline(aes(yintercept=CD3_IR_threshold), color="#BF0032",linetype="dashed") +
  annotate("text", x=20,y=CD3_IR_threshold+0.1,label="threshold: 0.014",color="#BF0032",size=3,fontface="bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 14))
ggsave("../new/result/without_PVTT_MLN/tumor/sub_T/CD3T_IR_threshold.pdf",CD3T_IR,width = 8,height = 2)

T.IR <- T.IR[order(T.IR[,3],decreasing = T),]
CD4T_IR <- ggbarplot(T.IR, x = "Patient", y = "CD4.IR", fill = "Patient") +
  NoLegend() +
  geom_hline(aes(yintercept=CD4_IR_threshold), color="#BF0032",linetype="dashed") +
  annotate("text", x=20,y=CD4_IR_threshold+0.1,label="threshold: 0.012",color="#BF0032",size=3,fontface="bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 14))
ggsave("../new/result/without_PVTT_MLN/tumor/sub_T/CD4T_IR_threshold.pdf",CD4T_IR,width = 8,height = 2)

T.IR <- T.IR[order(T.IR[,4],decreasing = T),]
CD8T_IR <- ggbarplot(T.IR, x = "Patient", y = "CD8.IR", fill = "Patient") +
  NoLegend() +
  geom_hline(aes(yintercept=CD8_IR_threshold), color="#BF0032",linetype="dashed") +
  annotate("text", x=20,y=CD8_IR_threshold+0.02,label="threshold: 0.0013",color="#BF0032",size=3,fontface="bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 14))
ggsave("../new/result/without_PVTT_MLN/tumor/sub_T/CD8T_IR_threshold.pdf",CD8T_IR,width = 8,height = 2)

tumor$Group1 <- ifelse(tumor$CD3T.IR > CD3_IR_threshold, "CD3+T high", "CD3+T low")
tumor$Group2 <- ifelse(tumor$CD4T.IR > CD4_IR_threshold, "CD4+T high", "CD4+T low")
tumor$Group3 <- ifelse(tumor$CD8T.IR > CD8_IR_threshold, "CD8+T high", "CD8+T low")

saveRDS(tumor,file = "rds/tumor_combine_T.IR.rds")

# liver stromal cell marker:VCAM1
p1 <- VlnPlot(tumor,features = "VCAM1",group.by = "CellType",pt.size = 0) + NoLegend()

current_ids <- c("B cell","CAF","Malignant cell","T cell","TAM","TEC","unclassified")
new_ids <- c("B cell","CAF","Malignant cell","T cell","TAM","TEC","Stromal cell")
tumor$CellType <- plyr::mapvalues(tumor$CellType, from = current_ids, to = new_ids)
current_ids <- c("B cell","CAF","CD4+ T cell","CD8+ T cell","Malignant cell","TAM","TEC","unclassified")
new_ids <- c("B cell","CAF","CD4+ T cell","CD8+ T cell","Malignant cell","TAM","TEC","Stromal cell")
tumor$CellType2 <- plyr::mapvalues(tumor$CellType2, from = current_ids, to = new_ids)

saveRDS(tumor,file = "rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR_2.rds")

pal <- c("#0066A5","#46A040","#00AF99","#875691","#98D9E9","#F6313E","#FFA300")

tumor <- readRDS("rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR_2.rds")

p1 <- DimPlot(tumor,reduction = "tsne",group.by = "CellType", cols = pal)

p_tsne_tumor <- DimPlot(tumor,reduction = "tsne",group.by = "CellType",cols = pal,
                      label = T,label.size = 5,label.color = "black") +
  theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.text = element_text(size = 14)) + NoLegend() + ggtitle("CellType")
ggsave("result/without_PVTT_MLN/tumor/CD3T/before_CD45_removal/0.CD3T_before_CD45_removal_celltype.pdf",p_tsne_tumor,width = 6,height = 6)
ggsave("result/without_PVTT_MLN/tumor/CD4T/before_CD45_removal/0.CD4T_before_CD45_removal_celltype.pdf",p_tsne_tumor,width = 6,height = 6)
ggsave("result/without_PVTT_MLN/tumor/CD8T/before_CD45_removal/0.CD8T_before_CD45_removal_celltype.pdf",p_tsne_tumor,width = 6,height = 6)

# CD3
p_CD3_IR <- DimPlot(tumor,reduction = "tsne",group.by = "Group1", cols = c("#BF0032","#0066A5"),pt.size = 0.8) + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        legend.text = element_text(size = 14))
ggsave("result/without_PVTT_MLN/tumor/CD3T/before_CD45_removal/1.CD3.high_vs_low.pdf",p_CD3_IR,width = 6,height = 6)

# CD4 T
p_CD4_IR <- DimPlot(tumor,reduction = "tsne",group.by = "Group2", cols = c("#BF0032","#0066A5"),pt.size = 0.8) + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        legend.text = element_text(size = 14))
ggsave("result/without_PVTT_MLN/tumor/CD4T/before_CD45_removal/1.CD4.high_vs_low.pdf",p_CD4_IR,width = 6,height = 6)

# CD8 T
p_CD8_IR <- DimPlot(tumor,reduction = "tsne",group.by = "Group3", cols = c("#BF0032","#0066A5"),pt.size = 0.8) + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        legend.text = element_text(size = 14))
ggsave("result/without_PVTT_MLN/tumor/CD8T/before_CD45_removal/1.CD8.high_vs_low.pdf",p_CD8_IR,width = 6,height = 6)


CD45 <- VlnPlot(tumor, features = "PTPRC", pt.size = 0, group.by = "seurat_clusters",cols = rainbow(length(levels(as.factor(tumor$seurat_clusters))))) +
  NoLegend() + coord_flip() + ggtitle("CD45")

ggsave("result/without_PVTT_MLN/tumor/CD3T/before_CD45_removal/2.CD45.pdf",CD45,width = 3,height = 8)
ggsave("result/without_PVTT_MLN/tumor/CD4T/before_CD45_removal/2.CD45.pdf",CD45,width = 3,height = 8)
ggsave("result/without_PVTT_MLN/tumor/CD8T/before_CD45_removal/2.CD45.pdf",CD45,width = 3,height = 8)



# 去除CD45高表达的细胞群
tumor_sub_cells <- rownames(subset(tumor@meta.data, seurat_clusters == 2 | seurat_clusters == 3 | seurat_clusters == 4 | seurat_clusters == 5 |
                                     seurat_clusters == 6 | seurat_clusters == 7 | seurat_clusters == 8 | seurat_clusters == 9 |
                                     seurat_clusters == 10 | seurat_clusters == 13 | seurat_clusters == 18))
tumor_sub <- subset(tumor, cells = tumor_sub_cells)

tumor_sub <- FindVariableFeatures(tumor_sub,selection.method = "vst",nfeatures = 2000)
DefaultAssay(tumor_sub) <- "integrated"
tumor_sub <- ScaleData(tumor_sub,features = rownames(tumor_sub))
tumor_sub <- RunPCA(tumor_sub, features = VariableFeatures(tumor_sub))
p1 <- DimPlot(tumor_sub, reduction = "pca", group.by = "orig.ident",cols=rainbow(length(levels(as.factor(tumor_sub$orig.ident)))))
ggsave("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/1.tumor_sub_pca.pdf",p1,device="pdf",width=8,height=6)
ggsave("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/1.tumor_sub_pca.pdf",p1,device="pdf",width=8,height=6)
ggsave("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/1.tumor_sub_pca.pdf",p1,device="pdf",width=8,height=6)
p2 <- ElbowPlot(tumor_sub,ndims=30,reduction="pca")
ggsave("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/2.tumor_sub_elbow.pdf",p2,device="pdf",width=6,height=4)
ggsave("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/2.tumor_sub_elbow.pdf",p2,device="pdf",width=6,height=4)
ggsave("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/2.tumor_sub_elbow.pdf",p2,device="pdf",width=6,height=4)

pc.num = 1:24
# 聚类
tumor_sub <- FindNeighbors(tumor_sub, dims = pc.num) %>% FindClusters(resolution = 0.3)
tumor_sub <- RunUMAP(tumor_sub,reduction="pca",dims=pc.num) %>% RunTSNE(dims=pc.num)

tumor_sub$seurat_clusters <- tumor_sub$integrated_snn_res.0.2
Idents(tumor_sub) <- tumor_sub$seurat_clusters


p_tsne_tumor_sub <- DimPlot(tumor_sub,reduction = "tsne",group.by = "CellType",
                        label = T,label.size = 5,label.color = "black") +
  theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.text = element_text(size = 14)) + NoLegend() + ggtitle("CellType")

ggsave("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/3.CD45_removal.celltype.pdf",p_tsne_tumor_sub, width = 6,height = 6)
ggsave("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/3.CD45_removal.celltype.pdf",p_tsne_tumor_sub, width = 6,height = 6)
ggsave("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/3.CD45_removal.celltype.pdf",p_tsne_tumor_sub, width = 6,height = 6)

p_tsne_CAF <- DimPlot(tumor_sub,reduction = "tsne",group.by = "CellType",
                            label = T,label.size = 5,label.color = "black",
                      cells.highlight = rownames(subset(tumor_sub@meta.data, CellType == "CAF"))) +
  theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.text = element_text(size = 14)) + NoLegend() + ggtitle("CellType")
ggsave("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/3.1.CD45_removal.CAF.pdf",p_tsne_CAF, width = 6,height = 6)
ggsave("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/3.1.CD45_removal.CAF.pdf",p_tsne_CAF, width = 6,height = 6)
ggsave("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/3.1.CD45_removal.CAF.pdf",p_tsne_CAF, width = 6,height = 6)

p_tsne_TEC <- DimPlot(tumor_sub,reduction = "tsne",group.by = "CellType",
                      label = T,label.size = 5,label.color = "black",
                      cells.highlight = rownames(subset(tumor_sub@meta.data, CellType == "TEC"))) +
  theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.text = element_text(size = 14)) + NoLegend() + ggtitle("CellType")
ggsave("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/3.2.CD45_removal.TEC.pdf",p_tsne_TEC, width = 6,height = 6)
ggsave("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/3.2.CD45_removal.TEC.pdf",p_tsne_TEC, width = 6,height = 6)
ggsave("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/3.2.CD45_removal.TEC.pdf",p_tsne_TEC, width = 6,height = 6)

p_tsne_MLG <- DimPlot(tumor_sub,reduction = "tsne",group.by = "CellType",
                      label = T,label.size = 5,label.color = "black",
                      cells.highlight = rownames(subset(tumor_sub@meta.data, CellType == "Malignant cell"))) +
  theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.text = element_text(size = 14)) + NoLegend() + ggtitle("CellType")
ggsave("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/3.3.CD45_removal.Malignant.cell.pdf",p_tsne_MLG, width = 6,height = 6)
ggsave("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/3.3.CD45_removal.Malignant.cell.pdf",p_tsne_MLG, width = 6,height = 6)
ggsave("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/3.3.CD45_removal.Malignant.cell.pdf",p_tsne_MLG, width = 6,height = 6)


# CD3+T
p_CD3 <- DimPlot(tumor_sub,reduction = "tsne",group.by = "Group1", cols = c("#BF0032","#0066A5")) + 
  ggtitle("CD3") + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

p_CD3_mod <- p_CD3 + 
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.text = element_text(size = 14))
ggsave("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/4.CD45_removal.CD3_high_low.pdf",p_CD3_mod,width = 4,height = 4)

# CD3+T high
p_CD3_high <- DimPlot(tumor_sub, reduction = "tsne",cols.highlight = "#BF0032",
                      cells.highlight = rownames(subset(tumor_sub@meta.data,Group1 == "CD3+T high")))+ 
  NoLegend() + ggtitle("CD3+T high") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/4.1.CD45_removal.CD3T_high.pdf",p_CD3_high,width = 4,height = 4)
# CD3+T low
p_CD3_low <- DimPlot(tumor_sub, reduction = "tsne",cols.highlight = "#0066A5",
                      cells.highlight = rownames(subset(tumor_sub@meta.data,Group1 == "CD3+T low")))+ 
  NoLegend() + ggtitle("CD3+T low") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/4.2.CD45_removal.CD3T_low.pdf",p_CD3_low,width = 4,height = 3.5)

# CD4+T
p_CD4 <- DimPlot(tumor_sub,reduction = "tsne",group.by = "Group2", cols = c("#BF0032","#0066A5")) + 
  ggtitle("CD4") + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

p_CD4_mod <- p_CD4 + 
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.text = element_text(size = 14))
ggsave("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/4.CD45_removal.CD4_high_low.pdf",p_CD4_mod,width = 4,height = 4)

# CD4+T high
p_CD4_high <- DimPlot(tumor_sub, reduction = "tsne",cols.highlight = "#BF0032",
                      cells.highlight = rownames(subset(tumor_sub@meta.data,Group2 == "CD4+T high")))+ 
  NoLegend() + ggtitle("CD4+T high") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/4.1.CD45_removal.CD4T_high.pdf",p_CD4_high,width = 4,height = 4)
# CD4+T low
p_CD4_low <- DimPlot(tumor_sub, reduction = "tsne",cols.highlight = "#0066A5",
                     cells.highlight = rownames(subset(tumor_sub@meta.data,Group2 == "CD4+T low")))+ 
  NoLegend() + ggtitle("CD4+T low") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/4.2.CD45_removal.CD4T_low.pdf",p_CD4_low,width = 4,height = 3.5)

# CD8+T
p_CD8 <- DimPlot(tumor_sub,reduction = "tsne",group.by = "Group3", cols = c("#BF0032","#0066A5")) + 
  ggtitle("CD8") + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

p_CD8_mod <- p_CD8 + 
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.text = element_text(size = 14))
ggsave("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/4.CD45_removal.CD8_high_low.pdf",p_CD8_mod,width = 4,height = 4)

# CD8+T high
p_CD8_high <- DimPlot(tumor_sub, reduction = "tsne",cols.highlight = "#BF0032",
                      cells.highlight = rownames(subset(tumor_sub@meta.data,Group3 == "CD8+T high")))+ 
  NoLegend() + ggtitle("CD8+T high") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/4.1.CD45_removal.CD8T_high.pdf",p_CD8_high,width = 4,height = 4)
# CD8+T low
p_CD8_low <- DimPlot(tumor_sub, reduction = "tsne",cols.highlight = "#0066A5",
                     cells.highlight = rownames(subset(tumor_sub@meta.data,Group3 == "CD8+T low")))+ 
  NoLegend() + ggtitle("CD8+T low") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/4.2.CD45_removal.CD8T_low.pdf",p_CD8_low,width = 4,height = 3.5)


saveRDS(tumor_sub,file = "rds/without_PVTT_MLN/tumor/CD45_removal_T.IR.2.rds")


## CD3+ T cell
if(F){
  tumor_CD3_high_cells <- subset(tumor@meta.data, Group1 == "CD3+T high")
  tumor_CD3_high <- subset(tumor, cells = rownames(tumor_CD3_high_cells))
  tumor_CD3_low_cells <- subset(tumor@meta.data, Group1 == "CD3+T low")
  tumor_CD3_low <- subset(tumor, cells = rownames(tumor_CD3_low_cells))
}


##  classification of tumor and non-tumor cells

if(F){
  pal <- c("#00441B","#46A040","#00AF99","#FFC179","#98D9E9","#F6313E","#FFA300")
  
  # CD3+T high
  table(tumor_CD3_high$CellType)
  p_tsne_CD3 <- DimPlot(tumor,reduction = "tsne",group.by = "CellType",cols = pal,
                        label = T,label.size = 5,label.color = "black") +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.text = element_text(size = 14)) + NoLegend() + ggtitle("CellType")
  ggsave("result/without_PVTT_MLN/tumor/CD3T/CD3T_tsne_celltype.pdf",p_tsne_CD3,width = 6,height = 6)
  
  p_tsne_CD3_high <- DimPlot(tumor_CD3_high,reduction = "tsne",group.by = "CellType",cols = pal,
                             label = T,label.size = 5,label.color = "black") +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.text = element_text(size = 14),
          plot.title = element_blank()) + NoLegend()
  ggsave("result/without_PVTT_MLN/tumor/CD3T/1.1.CD3T_high_tsne.pdf",p_tsne_CD3_high,width = 6,height = 6)
  
  p1 <- VlnPlot(tumor_CD3_high, features = "PTPRC",cols = pal,pt.size = 0,group.by = "CellType") +
    coord_flip() + ggtitle("CD45") + NoLegend() +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 14))
  ggsave("result/without_PVTT_MLN/tumor/CD3T/2.CD45_CD3_high.pdf", p1, device = "pdf", width = 4,height = 6)
  
  tumor_CD3_high$Class <- ifelse(
    tumor_CD3_high$CellType == "TEC" | tumor_CD3_high$CellType == "CAF" | tumor_CD3_high$CellType == "Malignant cell",
    "Tumor", "Nontumor"
  )
  
  p_CD3_high <- DimPlot(tumor_CD3_high,reduction = "tsne",group.by = "Class", cols = c("#0066A5","#BF0032")) + 
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_blank(),
          legend.position = "bottom",
          legend.justification = "center",
          legend.text = element_text(size = 14))
  ggsave("result/without_PVTT_MLN/tumor/CD3T/3.CD3_high_tumor_nontumor.pdf",p_CD3_high,width = 4,height = 4)
  
  # keep tumor cells
  tumor_CD3_high_nontumor_removal_cells <- subset(tumor_CD3_high@meta.data, Class == "Tumor")
  tumor_CD3_high_nontumor_removal <- subset(tumor_CD3_high, cells = rownames(tumor_CD3_high_nontumor_removal_cells))
  
  # CD3+T low
  table(tumor_CD3_low$CellType)
  
  p_tsne_CD3_low <- DimPlot(tumor_CD3_low,reduction = "tsne",group.by = "CellType",cols = pal,
                            label = T,label.size = 5,label.color = "black") +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.text = element_text(size = 14),
          plot.title = element_blank()) + NoLegend()
  ggsave("result/without_PVTT_MLN/tumor/CD3T/1.2.CD3T_low_tsne.pdf",p_tsne_CD3_low,width = 6,height = 6)
  
  p2 <- VlnPlot(tumor_CD3_low, features = "PTPRC",cols = pal,pt.size = 0,group.by = "CellType") +
    coord_flip() + ggtitle("CD45") + NoLegend() +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 14))
  ggsave("result/without_PVTT_MLN/tumor/CD3T/2.CD45_CD3_low.pdf",p2,device = "pdf", width = 4,height = 6)
  
  tumor_CD3_low$Class <- ifelse(
    tumor_CD3_low$CellType == "TEC" | tumor_CD3_low$CellType == "CAF" | tumor_CD3_low$CellType == "Malignant cell",
    "Tumor", "Nontumor"
  )
  
  p_CD3_low <- DimPlot(tumor_CD3_low,reduction = "tsne",group.by = "Class", cols = c("#0066A5","#BF0032")) + 
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_blank(),
          legend.position = "bottom",
          legend.justification = "center",
          legend.text = element_text(size = 14))
  ggsave("result/without_PVTT_MLN/tumor/CD3T/3.CD3_low_tumor_nontumor.pdf",p_CD3_low,width = 4,height = 4)
  
  # keep tumor cells
  tumor_CD3_low_nontumor_removal_cells <- subset(tumor_CD3_low@meta.data, Class == "Tumor")
  tumor_CD3_low_nontumor_removal <- subset(tumor_CD3_low, cells = rownames(tumor_CD3_low_nontumor_removal_cells))
  
  ## CD4+ T cell
  tumor_CD4_high_cells <- subset(tumor@meta.data, Group2 == "CD4+T high")
  tumor_CD4_high <- subset(tumor, cells = rownames(tumor_CD4_high_cells))
  tumor_CD4_low_cells <- subset(tumor@meta.data, Group2 == "CD4+T low")
  tumor_CD4_low <- subset(tumor, cells = rownames(tumor_CD4_low_cells))
  
  table(tumor_CD4_high$CellType)
  # CD4+T
  p_tsne_CD4 <- DimPlot(tumor,reduction = "tsne",group.by = "CellType",cols = pal,
                        label = T,label.size = 5,label.color = "black") +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.text = element_text(size = 14)) + NoLegend() + ggtitle("CellType")
  ggsave("result/without_PVTT_MLN/tumor/CD4T/CD4T_tsne_celltype.pdf",p_tsne_CD3,width = 6,height = 6)
  
  # CD4+ T high
  p_tsne_CD4_high <- DimPlot(tumor_CD4_high,reduction = "tsne",group.by = "CellType",cols = pal,
                             label = T,label.size = 5,label.color = "black") +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.text = element_text(size = 14),
          plot.title = element_blank()) + NoLegend()
  ggsave("result/without_PVTT_MLN/tumor/CD4T/1.1.CD4T_high_tsne.pdf",p_tsne_CD4_high,width = 6,height = 6)
  
  p3 <- VlnPlot(tumor_CD4_high, features = "PTPRC",cols = pal,pt.size = 0,group.by = "CellType") +
    coord_flip() + ggtitle("CD45") + NoLegend() +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 14))
  ggsave("result/without_PVTT_MLN/tumor/CD4T/2.CD45_CD4_high.pdf", p3, device = "pdf", width = 4,height = 6)
  
  tumor_CD4_high$Class <- ifelse(
    tumor_CD4_high$CellType == "TEC" | tumor_CD4_high$CellType == "CAF" | tumor_CD4_high$CellType == "Malignant cell",
    "Tumor", "Nontumor"
  )
  
  p_CD4_high <- DimPlot(tumor_CD4_high,reduction = "tsne",group.by = "Class", cols = c("#0066A5","#BF0032")) + 
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_blank(),
          legend.position = "bottom",
          legend.justification = "center",
          legend.text = element_text(size = 14))
  ggsave("result/without_PVTT_MLN/tumor/CD4T/3.CD4_high_tumor_nontumor.pdf",p_CD4_high,width = 4,height = 4)
  
  # keep tumor cells
  tumor_CD4_high_nontumor_removal_cells <- subset(tumor_CD4_high@meta.data, Class == "Tumor")
  tumor_CD4_high_nontumor_removal <- subset(tumor_CD4_high, cells = rownames(tumor_CD4_high_nontumor_removal_cells))
  
  # CD4+ T low
  table(tumor_CD4_low$CellType)
  
  p_tsne_CD4_low <- DimPlot(tumor_CD4_low,reduction = "tsne",group.by = "CellType",cols = pal,
                            label = T,label.size = 5,label.color = "black") +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.text = element_text(size = 14),
          plot.title = element_blank()) + NoLegend()
  ggsave("result/without_PVTT_MLN/tumor/CD4T/1.2.CD4T_low_tsne.pdf",p_tsne_CD4_low,width = 6,height = 6)
  
  p4 <- VlnPlot(tumor_CD4_low, features = "PTPRC",cols = pal,pt.size = 0,group.by = "CellType") +
    coord_flip() + ggtitle("CD45") + NoLegend() +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 14))
  ggsave("result/without_PVTT_MLN/tumor/CD4T/2.CD45_CD4_low.pdf", p4, device = "pdf", width = 4,height = 6)
  
  tumor_CD4_low$Class <- ifelse(
    tumor_CD4_low$CellType == "TEC" | tumor_CD4_low$CellType == "CAF" | tumor_CD4_low$CellType == "Malignant cell",
    "Tumor", "Nontumor"
  )
  
  p_CD4_low <- DimPlot(tumor_CD4_low,reduction = "tsne",group.by = "Class", cols = c("#0066A5","#BF0032")) + 
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_blank(),
          legend.position = "bottom",
          legend.justification = "center",
          legend.text = element_text(size = 14))
  ggsave("result/without_PVTT_MLN/tumor/CD4T/3.CD4_low_tumor_nontumor.pdf",p_CD4_low,width = 4,height = 4)
  
  # keep tumor cells
  tumor_CD4_low_nontumor_removal_cells <- subset(tumor_CD4_low@meta.data, Class == "Tumor")
  tumor_CD4_low_nontumor_removal <- subset(tumor_CD4_low, cells = rownames(tumor_CD4_low_nontumor_removal_cells))
  
  ## CD8+T cell
  tumor_CD8_high_cells <- subset(tumor@meta.data, Group3 == "CD8+T high")
  tumor_CD8_high <- subset(tumor, cells = rownames(tumor_CD8_high_cells))
  tumor_CD8_low_cells <- subset(tumor@meta.data, Group3 == "CD8+T low")
  tumor_CD8_low <- subset(tumor, cells = rownames(tumor_CD8_low_cells))
  
  # classify tumor and nontumor cells
  table(tumor_CD8_high$CellType)
  # CD8+T
  p_tsne_CD8 <- DimPlot(tumor,reduction = "tsne",group.by = "CellType",cols = pal,
                        label = T,label.size = 5,label.color = "black") +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.text = element_text(size = 14)) + NoLegend() + ggtitle("CellType")
  ggsave("result/without_PVTT_MLN/tumor/CD8T/CD8T_tsne_celltype.pdf",p_tsne_CD8,width = 6,height = 6)
  # CD8+T high
  p_tsne_CD8_high <- DimPlot(tumor_CD8_high,reduction = "tsne",group.by = "CellType",cols = pal,
                             label = T,label.size = 5,label.color = "black") +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.text = element_text(size = 14),
          plot.title = element_blank()) + NoLegend()
  ggsave("result/without_PVTT_MLN/tumor/CD8T/1.1.CD8T_high_tsne.pdf",p_tsne_CD8_high,width = 6,height = 6)
  
  p5 <- VlnPlot(tumor_CD8_high, features = "PTPRC",cols = pal,pt.size = 0,group.by = "CellType") +
    coord_flip() + ggtitle("CD45") + NoLegend() +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 14))
  ggsave("result/without_PVTT_MLN/tumor/CD8T/2.CD45_CD8_high.pdf", p5, device = "pdf", width = 4,height = 6)
  
  tumor_CD8_high$Class <- ifelse(
    tumor_CD8_high$CellType == "TEC" | tumor_CD8_high$CellType == "CAF" | tumor_CD8_high$CellType == "Malignant cell",
    "Tumor", "Nontumor"
  )
  
  p_CD8_high <- DimPlot(tumor_CD8_high,reduction = "tsne",group.by = "Class", cols = c("#0066A5","#BF0032")) + 
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_blank(),
          legend.position = "bottom",
          legend.justification = "center",
          legend.text = element_text(size = 14))
  ggsave("result/without_PVTT_MLN/tumor/CD8T/3.CD8_high_tumor_nontumor.pdf",p_CD8_high,width = 4,height = 4)
  
  # keep tumor cells
  tumor_CD8_high_nontumor_removal_cells <- subset(tumor_CD8_high@meta.data, Class == "Tumor")
  tumor_CD8_high_nontumor_removal <- subset(tumor_CD8_high, cells = rownames(tumor_CD8_high_nontumor_removal_cells))
  
  # classify tumor and nontumor cells
  table(tumor_CD8_low$CellType)
  
  p_tsne_CD8_low <- DimPlot(tumor_CD8_low,reduction = "tsne",group.by = "CellType",cols = pal,
                            label = T,label.size = 5,label.color = "black") +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.text = element_text(size = 14),
          plot.title = element_blank()) + NoLegend()
  ggsave("result/without_PVTT_MLN/tumor/CD8T/1.2.CD8T_low_tsne.pdf",p_tsne_CD8_low,width = 6,height = 6)
  
  p6 <- VlnPlot(tumor_CD8_low, features = "PTPRC",cols = pal,pt.size = 0,group.by = "CellType") +
    coord_flip() + ggtitle("CD45") + NoLegend() +
    theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
          panel.border = element_rect(size = 1.5,colour = "black"),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 14))
  ggsave("result/without_PVTT_MLN/tumor/CD8T/2.CD45_CD8_low.pdf", p6, device = "pdf", width = 4,height = 6)
  
  tumor_CD8_low$Class <- ifelse(
    tumor_CD8_low$CellType == "TEC" | tumor_CD8_low$CellType == "CAF" | tumor_CD8_low$CellType == "Malignant cell",
    "Tumor", "Nontumor"
  )
  
  p_CD8_low <- DimPlot(tumor_CD8_low,reduction = "tsne",group.by = "Class", cols = c("#0066A5","#BF0032")) + 
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_blank(),
          legend.position = "bottom",
          legend.justification = "center",
          legend.text = element_text(size = 14))
  ggsave("result/without_PVTT_MLN/tumor/CD8T/3.CD8_low_tumor_nontumor.pdf",p_CD8_low,width = 4,height = 4)
  
  # keep tumor cells
  tumor_CD8_low_nontumor_removal_cells <- subset(tumor_CD8_low@meta.data, Class == "Tumor")
  tumor_CD8_low_nontumor_removal <- subset(tumor_CD8_low, cells = rownames(tumor_CD8_low_nontumor_removal_cells))
  
  save(tumor_CD3_high, file = "rds/without_PVTT_MLN/tumor/CD3T/tumor_CD3_high.RData")
  save(tumor_CD3_low, file = "rds/without_PVTT_MLN/tumor/CD3T/tumor_CD3_low.RData")
  save(tumor_CD4_high, file = "rds/without_PVTT_MLN/tumor/CD4T/tumor_CD4_high.RData")
  save(tumor_CD4_low, file = "rds/without_PVTT_MLN/tumor/CD4T/tumor_CD4_low.RData")
  save(tumor_CD8_high, file = "rds/without_PVTT_MLN/tumor/CD8T/tumor_CD8_high.RData")
  save(tumor_CD8_low, file = "rds/without_PVTT_MLN/tumor/CD8T/tumor_CD8_low.RData")
  
  #save in one RData file
  save(
    tumor_CD3_high,
    tumor_CD3_low,
    tumor_CD4_high,
    tumor_CD4_low,
    tumor_CD8_high,
    tumor_CD8_low, 
    file = "rds/without_PVTT_MLN/tumor/3.tumor_Tcell_tumor_marked.RData"
  )
  
  save(
    tumor_CD3_high_nontumor_removal,
    tumor_CD3_low_nontumor_removal,
    tumor_CD4_high_nontumor_removal,
    tumor_CD4_low_nontumor_removal,
    tumor_CD8_high_nontumor_removal,
    tumor_CD8_low_nontumor_removal,
    file = "rds/without_PVTT_MLN/tumor/4.tumor_Tcell_nontumor_removal.RData"
  )
}


## >>>>>>>>>>>>>>>>>>>>>>>>> Differential Expression Analysis <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

library(Seurat)
rm(list = ls());gc()

tumor_sub <- readRDS("rds/without_PVTT_MLN/tumor/CD45_removal_T.IR.2.rds")
DefaultAssay(tumor_sub) <- "RNA"
Idents(tumor_sub) <- tumor_sub$CellType

##>>>>>> CD3 high组和CD3 low组 <<<<<<<

# CAF细胞亚群差异
CD3T.CAF.DE <- FindMarkers(tumor_sub,ident.1 = "CD3+T high", ident.2 = "CD3+T low", group.by = "Group1", subset.ident = "CAF",
                       only.pos=FALSE,min.diff.pct = 0.1,logfc.threshold = 0)
# TEC细胞亚群差异
CD3T.TEC.DE <- FindMarkers(tumor_sub,ident.1 = "CD3+T high", ident.2 = "CD3+T low", group.by = "Group1", subset.ident = "TEC",
                           only.pos=FALSE,min.diff.pct = 0.1,logfc.threshold = 0)
# Malignant cell细胞亚群差异
CD3T.MLG.DE <- FindMarkers(tumor_sub,ident.1 = "CD3+T high", ident.2 = "CD3+T low", group.by = "Group1", subset.ident = "Malignant cell",
                           only.pos=FALSE,min.diff.pct = 0.1,logfc.threshold = 0.25)

### TO EXCEL

## CAF
CD3T.CAF.DE.UP <- CD3T.CAF.DE[CD3T.CAF.DE$avg_logFC > 0,]
CD3T.CAF.DE.UP <- CD3T.CAF.DE.UP[order(CD3T.CAF.DE.UP[, "avg_logFC"], decreasing = T), ]
xlsx::write.xlsx(CD3T.CAF.DE.UP,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF.CD3T.hi_vs_lo.xlsx", sheetName = "UP")
CD3T.CAF.DE.Down <- CD3T.CAF.DE[CD3T.CAF.DE$avg_logFC < 0,]
CD3T.CAF.DE.Down <- CD3T.CAF.DE.Down[order(CD3T.CAF.DE.Down[, "avg_logFC"], decreasing = F), ]
xlsx::write.xlsx(CD3T.CAF.DE.Down,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF.CD3T.hi_vs_lo.xlsx", sheetName = "Down",append = T)

## TEC
CD3T.TEC.DE.UP <- CD3T.TEC.DE[CD3T.TEC.DE$avg_logFC > 0,]
CD3T.TEC.DE.UP <- CD3T.TEC.DE.UP[order(CD3T.TEC.DE.UP[, "avg_logFC"], decreasing = T), ]
xlsx::write.xlsx(CD3T.TEC.DE.UP,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC.CD3T.hi_vs_lo.xlsx", sheetName = "UP")
CD3T.TEC.DE.Down <- CD3T.TEC.DE[CD3T.TEC.DE$avg_logFC < 0,]
CD3T.TEC.DE.Down <- CD3T.TEC.DE.Down[order(CD3T.TEC.DE.Down[, "avg_logFC"], decreasing = F), ]
xlsx::write.xlsx(CD3T.TEC.DE.Down,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC.CD3T.hi_vs_lo.xlsx", sheetName = "Down",append = T)

## Malignant cell
CD3T.MLG.DE.UP <- CD3T.MLG.DE[CD3T.MLG.DE$avg_logFC > 0,]
CD3T.MLG.DE.UP <- CD3T.MLG.DE.UP[order(CD3T.MLG.DE.UP[, "avg_logFC"], decreasing = T), ]
xlsx::write.xlsx(CD3T.MLG.DE.UP,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Maligant.Cell.CD3T.hi_vs_lo.xlsx", sheetName = "UP")
CD3T.MLG.DE.Down <- CD3T.MLG.DE[CD3T.MLG.DE$avg_logFC < 0,]
CD3T.MLG.DE.Down <- CD3T.MLG.DE.Down[order(CD3T.MLG.DE.Down[, "avg_logFC"], decreasing = F), ]
xlsx::write.xlsx(CD3T.MLG.DE.Down,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Maligant.Cell.CD3T.hi_vs_lo.xlsx", sheetName = "Down",append = T)


##>>>>>> CD4 high组和CD4 low组 <<<<<<<

# CAF细胞亚群差异
CD4T.CAF.DE <- FindMarkers(tumor_sub,ident.1 = "CD4+T high", ident.2 = "CD4+T low", group.by = "Group2", subset.ident = "CAF",
                           only.pos=FALSE,min.diff.pct = 0.1,logfc.threshold = 0)
# TEC细胞亚群差异
CD4T.TEC.DE <- FindMarkers(tumor_sub,ident.1 = "CD4+T high", ident.2 = "CD4+T low", group.by = "Group2", subset.ident = "TEC",
                           only.pos=FALSE,min.diff.pct = 0.1,logfc.threshold = 0)
# Malignant cell细胞亚群差异
CD4T.MLG.DE <- FindMarkers(tumor_sub,ident.1 = "CD4+T high", ident.2 = "CD4+T low", group.by = "Group2", subset.ident = "Malignant cell",
                           only.pos=FALSE,min.diff.pct = 0.1,logfc.threshold = 0.25)

### TO EXCEL

## CAF
CD4T.CAF.DE.UP <- CD4T.CAF.DE[CD4T.CAF.DE$avg_logFC > 0,]
CD4T.CAF.DE.UP <- CD4T.CAF.DE.UP[order(CD4T.CAF.DE.UP[, "avg_logFC"], decreasing = T), ]
xlsx::write.xlsx(CD4T.CAF.DE.UP,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF.CD4T.hi_vs_lo.xlsx", sheetName = "UP")
CD4T.CAF.DE.Down <- CD4T.CAF.DE[CD4T.CAF.DE$avg_logFC < 0,]
CD4T.CAF.DE.Down <- CD4T.CAF.DE.Down[order(CD4T.CAF.DE.Down[, "avg_logFC"], decreasing = F), ]
xlsx::write.xlsx(CD4T.CAF.DE.Down,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF.CD4T.hi_vs_lo.xlsx", sheetName = "Down",append = T)

## TEC
CD4T.TEC.DE.UP <- CD4T.TEC.DE[CD4T.TEC.DE$avg_logFC > 0,]
CD4T.TEC.DE.UP <- CD4T.TEC.DE.UP[order(CD4T.TEC.DE.UP[, "avg_logFC"], decreasing = T), ]
xlsx::write.xlsx(CD4T.TEC.DE.UP,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC.CD4T.hi_vs_lo.xlsx", sheetName = "UP")
CD4T.TEC.DE.Down <- CD4T.TEC.DE[CD4T.TEC.DE$avg_logFC < 0,]
CD4T.TEC.DE.Down <- CD4T.TEC.DE.Down[order(CD4T.TEC.DE.Down[, "avg_logFC"], decreasing = F), ]
xlsx::write.xlsx(CD4T.TEC.DE.Down,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC.CD4T.hi_vs_lo.xlsx", sheetName = "Down",append = T)

## Malignant cell
CD4T.MLG.DE.UP <- CD4T.MLG.DE[CD4T.MLG.DE$avg_logFC > 0,]
CD4T.MLG.DE.UP <- CD4T.MLG.DE.UP[order(CD4T.MLG.DE.UP[, "avg_logFC"], decreasing = T), ]
xlsx::write.xlsx(CD4T.MLG.DE.UP,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Maligant.Cell.CD4T.hi_vs_lo.xlsx", sheetName = "UP")
CD4T.MLG.DE.Down <- CD4T.MLG.DE[CD4T.MLG.DE$avg_logFC < 0,]
CD4T.MLG.DE.Down <- CD4T.MLG.DE.Down[order(CD4T.MLG.DE.Down[, "avg_logFC"], decreasing = F), ]
xlsx::write.xlsx(CD4T.MLG.DE.Down,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Maligant.Cell.CD4T.hi_vs_lo.xlsx", sheetName = "Down",append = T)


##>>>>>> CD8 high组和CD8 low组 <<<<<<<

# CAF细胞亚群差异
CD8T.CAF.DE <- FindMarkers(tumor_sub,ident.1 = "CD8+T high", ident.2 = "CD8+T low", group.by = "Group3", subset.ident = "CAF",
                           only.pos=FALSE,min.diff.pct = 0.1,logfc.threshold = 0)
# TEC细胞亚群差异
CD8T.TEC.DE <- FindMarkers(tumor_sub,ident.1 = "CD8+T high", ident.2 = "CD8+T low", group.by = "Group3", subset.ident = "TEC",
                           only.pos=FALSE,min.diff.pct = 0.1,logfc.threshold = 0.1)
# Malignant cell细胞亚群差异
CD8T.MLG.DE <- FindMarkers(tumor_sub,ident.1 = "CD8+T high", ident.2 = "CD8+T low", group.by = "Group3", subset.ident = "Malignant cell",
                           only.pos=FALSE,min.diff.pct = 0.1,logfc.threshold = 0.25)

### TO EXCEL

## CAF
CD8T.CAF.DE.UP <- CD8T.CAF.DE[CD8T.CAF.DE$avg_logFC > 0,]
CD8T.CAF.DE.UP <- CD8T.CAF.DE.UP[order(CD8T.CAF.DE.UP[, "avg_logFC"], decreasing = T), ]
xlsx::write.xlsx(CD8T.CAF.DE.UP,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF.CD8T.hi_vs_lo.xlsx", sheetName = "UP")
CD8T.CAF.DE.Down <- CD8T.CAF.DE[CD8T.CAF.DE$avg_logFC < 0,]
CD8T.CAF.DE.Down <- CD8T.CAF.DE.Down[order(CD8T.CAF.DE.Down[, "avg_logFC"], decreasing = F), ]
xlsx::write.xlsx(CD8T.CAF.DE.Down,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF.CD8T.hi_vs_lo.xlsx", sheetName = "Down",append = T)

## TEC
CD8T.TEC.DE.UP <- CD8T.TEC.DE[CD8T.TEC.DE$avg_logFC > 0,]
CD8T.TEC.DE.UP <- CD8T.TEC.DE.UP[order(CD8T.TEC.DE.UP[, "avg_logFC"], decreasing = T), ]
xlsx::write.xlsx(CD8T.TEC.DE.UP,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC.CD8T.hi_vs_lo.xlsx", sheetName = "UP")
CD8T.TEC.DE.Down <- CD8T.TEC.DE[CD8T.TEC.DE$avg_logFC < 0,]
CD8T.TEC.DE.Down <- CD8T.TEC.DE.Down[order(CD8T.TEC.DE.Down[, "avg_logFC"], decreasing = F), ]
xlsx::write.xlsx(CD8T.TEC.DE.Down,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC.CD8T.hi_vs_lo.xlsx", sheetName = "Down",append = T)

## Malignant cell
CD8T.MLG.DE.UP <- CD8T.MLG.DE[CD8T.MLG.DE$avg_logFC > 0,]
CD8T.MLG.DE.UP <- CD8T.MLG.DE.UP[order(CD8T.MLG.DE.UP[, "avg_logFC"], decreasing = T), ]
xlsx::write.xlsx(CD8T.MLG.DE.UP,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Maligant.Cell.CD8T.hi_vs_lo.xlsx", sheetName = "UP")
CD8T.MLG.DE.Down <- CD8T.MLG.DE[CD8T.MLG.DE$avg_logFC < 0,]
CD8T.MLG.DE.Down <- CD8T.MLG.DE.Down[order(CD8T.MLG.DE.Down[, "avg_logFC"], decreasing = F), ]
xlsx::write.xlsx(CD8T.MLG.DE.Down,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Maligant.Cell.CD8T.hi_vs_lo.xlsx", sheetName = "Down",append = T)

## Volcano

library(EnhancedVolcano)




# CD3+T
pv <- EnhancedVolcano(CD3T.CAF.DE, lab = rownames(CD3T.CAF.DE), x = 'avg_logFC', y = 'p_val_adj',
                      subtitle = "CAF",title = "CD3+T high vs CD3+T low",labSize = 5) + theme(legend.position = "none") + coord_flip()
ggsave("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF.volcano.pdf",pv,width = 6,height = 9)

pv <- EnhancedVolcano(CD3T.TEC.DE, lab = rownames(CD3T.TEC.DE), x = 'avg_logFC', y = 'p_val_adj',
                      subtitle = "TEC",title = "CD3+T high vs CD3+T low",labSize = 5) + theme(legend.position = "none") + coord_flip()
ggsave("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC.volcano.pdf",pv,width = 6,height = 9)

pv <- EnhancedVolcano(CD3T.MLG.DE, lab = rownames(CD3T.MLG.DE), x = 'avg_logFC', y = 'p_val_adj',
                      subtitle = "Malignant cell",title = "CD3+T high vs CD3+T low",labSize = 5) + 
  theme(legend.position = "none")
ggsave("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant.cell.volcano.pdf",pv,width = 7,height = 9)

# CD4+T
pv <- EnhancedVolcano(CD4T.CAF.DE, lab = rownames(CD4T.CAF.DE), x = 'avg_logFC', y = 'p_val_adj',
                       subtitle = "CAF",title = "CD4+T high vs CD4+T low",labSize = 5) + theme(legend.position = "none") + coord_flip()
ggsave("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF.volcano.pdf",pv,width = 6,height = 9)

pv <- EnhancedVolcano(CD4T.TEC.DE, lab = rownames(CD4T.TEC.DE), x = 'avg_logFC', y = 'p_val_adj',
                       subtitle = "TEC",title = "CD4+T high vs CD4+T low",labSize = 5) + 
  theme(legend.position = "none") + coord_flip()
ggsave("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC.volcano.pdf",pv,width = 6,height = 9)

pv <- EnhancedVolcano(CD4T.MLG.DE, lab = rownames(CD4T.MLG.DE), x = 'avg_logFC', y = 'p_val_adj',
                       subtitle = "Malignant cell",title = "CD4+T high vs CD4+T low",labSize = 5) + 
  theme(legend.position = "none")
ggsave("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant.cell.volcano.pdf",pv,width = 7,height = 9)

# CD8+T
pv <- EnhancedVolcano(CD8T.CAF.DE, lab = rownames(CD8T.CAF.DE), x = 'avg_logFC', y = 'p_val_adj',
                      subtitle = "CAF",title = "CD8+T high vs CD8+T low",labSize = 5) + 
  theme(legend.position = "none") + coord_flip()
ggsave("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF.volcano.pdf",pv,width = 6,height = 9)

pv <- EnhancedVolcano(CD8T.TEC.DE, lab = rownames(CD8T.TEC.DE), x = 'avg_logFC', y = 'p_val_adj',
                      subtitle = "TEC",title = "CD8+T high vs CD8+T low",labSize = 5) + 
  theme(legend.position = "none") + coord_flip()
ggsave("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC.volcano.pdf",pv,width = 6,height = 9)

pv <- EnhancedVolcano(CD8T.MLG.DE, lab = rownames(CD8T.MLG.DE), x = 'avg_logFC', y = 'p_val_adj',
                      subtitle = "Malignant cell",title = "CD8+T high vs CD8+T low",labSize = 5) + 
  theme(legend.position = "none") 
ggsave("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant.cell.volcano.pdf",pv,width = 7,height = 9)



if(F){
  # CD3+T high vs CD3+T low
  CD3T_markers <- FindMarkers(tumor,ident.1="CD3+T high",ident.2="CD3+T low",group.by="Group1",only.pos=FALSE,min.diff.pct=0.25,logfc.threshold=0.25)
  CD3T_markers$Cluster <- ifelse(CD3T_markers$avg_logFC > 0, "CD3+T high", "CD3+T low")
  # CD4+T high vs CD4+T low
  CD4T_markers <- FindMarkers(tumor,ident.1="CD4+T high",ident.2="CD4+T low",group.by="Group2",only.pos=FALSE,min.diff.pct=0.25,logfc.threshold=0.25)
  CD4T_markers$Cluster <- ifelse(CD4T_markers$avg_logFC > 0, "CD4+T high", "CD4+T low")
  # CD8+T high vs CD8+T low
  CD8T_markers <- FindMarkers(tumor,ident.1="CD8+T high",ident.2="CD8+T low",group.by="Group3",only.pos=FALSE,min.diff.pct=0.25,logfc.threshold=0.25)
  CD8T_markers$Cluster <- ifelse(CD8T_markers$avg_logFC > 0, "CD8+T high", "CD8+T low")
  save(CD3T_markers,CD4T_markers,CD8T_markers,file = "rds/without_PVTT_MLN/tumor/sub_T/1.CDT_high_low_DEmarkers.RData")
  
  ### write to EXCEL
  
  # CD3+T cell
  CD3.markers.up <- CD3T_markers[CD3T_markers$avg_logFC > 0,]
  CD3.markers.up <- CD3.markers.up[order(CD3.markers.up[,"avg_logFC"],decreasing = T),]
  CD3.markers.down <- CD3T_markers[CD3T_markers$avg_logFC < 0, ]
  CD3.markers.down <- CD3.markers.down[order(CD3.markers.down[,"avg_logFC"],decreasing = F),]
  xlsx::write.xlsx(CD3.markers.up,file = "result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",sheetName = "UP")
  xlsx::write.xlsx(CD3.markers.down,file = "result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",sheetName = "Down",append = T)
  
  # CD4+T cell
  CD4.markers.up <- CD4T_markers[CD4T_markers$avg_logFC > 0, ]
  CD4.markers.up <- CD4.markers.up[order(CD4.markers.up[,"avg_logFC"],decreasing = T),]
  CD4.markers.down <- CD4T_markers[CD4T_markers$avg_logFC < 0, ]
  CD4.markers.down <- CD4.markers.down[order(CD4.markers.down[,"avg_logFC"],decreasing = F),]
  xlsx::write.xlsx(CD4.markers.up,file = "result/without_PVTT_MLN/tumor/CD4T/CD4.high_vs_low.DEmarkers.xlsx",sheetName = "UP")
  xlsx::write.xlsx(CD4.markers.down,file = "result/without_PVTT_MLN/tumor/CD4T/CD4.high_vs_low.DEmarkers.xlsx",sheetName = "Down",append = T)
  
  # CD8+T cell
  CD8.markers.up <- CD8T_markers[CD8T_markers$avg_logFC > 0, ]
  CD8.markers.up <- CD8.markers.up[order(CD8.markers.up[,"avg_logFC"],decreasing = T),]
  CD8.markers.down <- CD8T_markers[CD8T_markers$avg_logFC < 0, ]
  CD8.markers.down <- CD8.markers.down[order(CD8.markers.down[,"avg_logFC"],decreasing = F),]
  xlsx::write.xlsx(CD8.markers.up,file = "result/without_PVTT_MLN/tumor/CD8T/CD8.high_vs_low.DEmarkers.xlsx",sheetName = "UP")
  xlsx::write.xlsx(CD8.markers.down,file = "result/without_PVTT_MLN/tumor/CD8T/CD8.high_vs_low.DEmarkers.xlsx",sheetName = "Down",append = T)
  
  if(F){
    save(
      CD3.markers.up, CD3.markers.down,
      CD4.markers.up, CD4.markers.down,
      CD8.markers.up, CD8.markers.down,
      file = "rds/without_PVTT_MLN/tumor/sub_T/CDT_high_vs_low_up_down_markers.RData"
    )
  }
}



##>>>>>> CD4high组和CD4low组患者中CD4细胞亚群 <<<<<<<

rm(list = ls());gc()
library(Seurat)

tumor <- readRDS("rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR_2.rds")

## CD4+T cell over CD4+T high vs. CD4+T low
Idents(tumor) <- tumor$CellType2

CD4T.DE <- FindMarkers(tumor,ident.1 = "CD4+T high", ident.2 = "CD4+T low", group.by = "Group2", subset.ident = "CD4+ T cell",
                       only.pos=FALSE,min.diff.pct=0.1,logfc.threshold=0.25)
CD4T.DE.UP <- CD4T.DE[CD4T.DE$avg_logFC > 0,]
CD4T.DE.UP <- CD4T.DE.UP[order(CD4T.DE.UP[,"avg_logFC"], decreasing = T),]
xlsx::write.xlsx(CD4T.DE.UP,file = "result/without_PVTT_MLN/tumor/DE/CD4T/CD4T_DEgenes.CD4T.hi_vs_lo.xlsx", sheetName = "UP")
CD4T.DE.Down <- CD4T.DE[CD4T.DE$avg_logFC < 0,]
CD4T.DE.Down <- CD4T.DE.Down[order(CD4T.DE.Down[,'avg_logFC'], decreasing = F),]
xlsx::write.xlsx(CD4T.DE.Down, file = "result/without_PVTT_MLN/tumor/DE/CD4T/CD4T_DEgenes.CD4T.hi_vs_lo.xlsx", sheetName = "Down",append = T)

### Volcano 
library(EnhancedVolcano)
pv1 <- EnhancedVolcano(CD4T.DE, lab = rownames(CD4T.DE), x = 'avg_logFC', y = 'p_val_adj',
                       subtitle = "",title = "",labSize = 4) + theme(legend.position = "none")
ggsave("result/without_PVTT_MLN/tumor/DE/CD4T/CD4T.volcano.pdf",pv1,width = 10,height = 8)

# CD4+T
p_CD4 <- DimPlot(tumor,reduction = "tsne",group.by = "Group2", cols = c("#BF0032","#0066A5")) + 
  ggtitle("CD4") + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        legend.text = element_text(size = 14))

ggsave("result/without_PVTT_MLN/tumor/DE/CD4T/CD4_high_low.pdf",p_CD4,width = 4,height = 4)

pal <- c("#00441B","#46A040","#00AF99","#FFC179","#98D9E9","#F6313E","#FFA300","#C390D4")

p_CD4T <- DimPlot(tumor, reduction = "tsne",cols.highlight = "#BF0032",pt.size = 0.8,
                     cells.highlight = rownames(subset(tumor@meta.data,CellType2 == "CD4+ T cell")))+ 
  NoLegend() + ggtitle("CD4+T cell") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("result/without_PVTT_MLN/tumor/DE/CD4T/CD4T_highlight.pdf",p_CD4T,width = 4,height = 4)

##>>>>>> CD8high组和CD8low组患者中CD8细胞亚群 <<<<<<<

# CD8+T cell over CD8+T high vs. CD8+T low
CD8T.DE <- FindMarkers(tumor,ident.1 = "CD8+T high", ident.2 = "CD8+T low", group.by = "Group3", subset.ident = "CD8+ T cell",
                       only.pos=FALSE,min.diff.pct=0.1,logfc.threshold=0.25)
CD8T.DE.UP <- CD8T.DE[CD8T.DE$avg_logFC > 0,]
CD8T.DE.UP <- CD8T.DE.UP[order(CD8T.DE.UP[,"avg_logFC"], decreasing = T),]
xlsx::write.xlsx(CD8T.DE.UP,file = "result/without_PVTT_MLN/tumor/DE/CD8T/CD8T_DEgenes.CD8T.hi_vs_lo.xlsx", sheetName = "UP")
CD8T.DE.Down <- CD8T.DE[CD8T.DE$avg_logFC < 0,]
CD8T.DE.Down <- CD8T.DE.Down[order(CD8T.DE.Down[,'avg_logFC'], decreasing = F),]
xlsx::write.xlsx(CD8T.DE.Down, file = "result/without_PVTT_MLN/tumor/DE/CD8T/CD8T_DEgenes.CD8T.hi_vs_lo.xlsx", sheetName = "Down",append = T)

### Volcano 
library(EnhancedVolcano)
pv2 <- EnhancedVolcano(CD8T.DE, lab = rownames(CD8T.DE), x = 'avg_logFC', y = 'p_val_adj',
                       subtitle = "",title = "",labSize = 4) + theme(legend.position = "none")
ggsave("result/without_PVTT_MLN/tumor/DE/CD8T/CD8T.volcano.pdf",pv2,width = 10,height = 8)


# CD8+T
p_CD8 <- DimPlot(tumor,reduction = "tsne",group.by = "Group3", cols = c("#BF0032","#0066A5")) + 
  ggtitle("CD8") + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        legend.text = element_text(size = 14))

ggsave("result/without_PVTT_MLN/tumor/DE/CD8T/CD8_high_low.pdf",p_CD8,width = 4,height = 4)

p_CD8T <- DimPlot(tumor, reduction = "tsne",cols.highlight = "#BF0032",pt.size = 0.8,
                  cells.highlight = rownames(subset(tumor@meta.data,CellType2 == "CD8+ T cell")))+ 
  NoLegend() + ggtitle("CD8+T cell") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave("result/without_PVTT_MLN/tumor/DE/CD8T/CD8T_highlight.pdf",p_CD8T,width = 4,height = 4)


##>>>>>>>> unused <<<<<<<<<<<<<<<<<<<<<<<<
if(FALSE){
  rm(list = ls());gc()
  load("rds/without_PVTT_MLN/tumor/3.tumor_Tcell_tumor_marked.RData")
  
  ### Tumor vs Non tumor
  
  # CD8+T low
  CD8_low.markers <- FindMarkers(tumor_CD8_low,ident.1="Tumor",ident.2="Nontumor",group.by="Class",only.pos=FALSE,min.diff.pct=0.25,logfc.threshold=0.25)
  CD8_low.markers$Cluster <- ifelse(CD8_low.markers$avg_logFC > 0, "Tumor", "Nontumor")
  # CD8+T high
  CD8_high.markers <- FindMarkers(tumor_CD8_high,ident.1="Tumor",ident.2="Nontumor",group.by="Class",only.pos=FALSE,min.diff.pct=0.25,logfc.threshold=0.25)
  CD8_high.markers$Cluster <- ifelse(CD8_high.markers$avg_logFC > 0, "Tumor", "Nontumor")
  # CD3+T low
  CD3_low.markers <- FindMarkers(tumor_CD3_low,ident.1="Tumor",ident.2="Nontumor",group.by="Class",only.pos=FALSE,min.diff.pct=0.25,logfc.threshold=0.25)
  CD3_low.markers$Cluster <- ifelse(CD3_low.markers$avg_logFC > 0, "Tumor", "Nontumor")
  # CD3+T high
  CD3_high.markers <- FindMarkers(tumor_CD3_high,ident.1="Tumor",ident.2="Nontumor",group.by="Class",only.pos=FALSE,min.diff.pct=0.25,logfc.threshold=0.25)
  CD3_high.markers$Cluster <- ifelse(CD3_high.markers$avg_logFC > 0, "Tumor", "Nontumor")
  # CD4+T low
  CD4_low.markers <- FindMarkers(tumor_CD4_low,ident.1="Tumor",ident.2="Nontumor",group.by="Class",only.pos=FALSE,min.diff.pct=0.25,logfc.threshold=0.25)
  CD4_low.markers$Cluster <- ifelse(CD4_low.markers$avg_logFC > 0, "Tumor", "Nontumor")
  # CD4+T high
  CD4_high.markers <- FindMarkers(tumor_CD4_high,ident.1="Tumor",ident.2="Nontumor",group.by="Class",only.pos=FALSE,min.diff.pct=0.25, logfc.threshold=0.25)
  CD4_high.markers$Cluster <- ifelse(CD4_high.markers$avg_logFC > 0, "Tumor", "Nontumor")
  
  save(
    CD8_low.markers = CD8_low.markers,
    CD8_high.markers = CD8_high.markers,
    CD3_low.markers = CD3_low.markers,
    CD3_high.markers = CD3_high.markers,
    CD4_low.markers = CD4_low.markers,
    CD4_high.markers = CD4_high.markers, 
    file = "rds/without_PVTT_MLN/tumor/sub_T/2.tumor_Tumor_vs_nonTumor_DEmarkers.RData"
  )
  
  # save results to EXCEL
  CD3_high.up <- CD3_high.markers[CD3_high.markers$avg_logFC > 0,]
  CD3_high.up <- CD3_high.up[order(CD3_high.up[,"avg_logFC"],decreasing = T),]
  CD3_high.down <- CD3_high.markers[CD3_high.markers$avg_logFC < 0, ]
  CD3_high.down <- CD3_high.down[order(CD3_high.down[,"avg_logFC"],decreasing = F),]
  xlsx::write.xlsx(CD3_high.up,file = "result/without_PVTT_MLN/tumor/sub_T/DEmarkers/CD3.high.tumor_vs_nonTumor.DEmarkers.xlsx",sheetName = "UP")
  xlsx::write.xlsx(CD3_high.down,file = "result/without_PVTT_MLN/tumor/CD3T/CD3.high.tumor_vs_nonTumor.DEmarkers.xlsx",sheetName = "Down",append = T)
  
  CD3_low.up <- CD3_low.markers[CD3_low.markers$avg_logFC > 0, ]
  CD3_low.up <- CD3_low.up[order(CD3_low.up[,"avg_logFC"],decreasing = T),]
  CD3_low.down <- CD3_low.markers[CD3_low.markers$avg_logFC < 0, ]
  CD3_low.down <- CD3_low.down[order(CD3_low.down[,"avg_logFC"],decreasing = F),]
  xlsx::write.xlsx(CD3_low.up, file = "result/without_PVTT_MLN/tumor/CD3T/CD3.low.tumor_vs_nonTumor.DEmarkers.xlsx",sheetName = "UP")
  xlsx::write.xlsx(CD3_low.down, file = "result/without_PVTT_MLN/tumor/CD3T/CD3.low.tumor_vs_nonTumor.DEmarkers.xlsx",sheetName = "Down",append = T)
  
  CD4_high.up <- CD4_high.markers[CD4_high.markers$avg_logFC > 0,]
  CD4_high.up <- CD4_high.up[order(CD4_high.up[,"avg_logFC"],decreasing = T),]
  CD4_high.down <- CD4_high.markers[CD4_high.markers$avg_logFC < 0,]
  CD4_high.down <- CD4_high.down[order(CD4_high.down[,"avg_logFC"],decreasing = F),]
  xlsx::write.xlsx(CD4_high.up,file = "result/without_PVTT_MLN/tumor/CD4T/CD4.high.tumor_vs_nonTumor.DEmarkers.xlsx",sheetName = "UP")
  xlsx::write.xlsx(CD4_high.down,file = "result/without_PVTT_MLN/tumor/CD4T/CD4.high.tumor_vs_nonTumor.DEmarkers.xlsx",sheetName = "Down",append = T)
  
  CD4_low.up <- CD4_low.markers[CD4_low.markers$avg_logFC > 0, ]
  CD4_low.up <- CD4_low.up[order(CD4_low.up[,"avg_logFC"],decreasing = T),]
  CD4_low.down <- CD4_low.markers[CD4_low.markers$avg_logFC < 0, ]
  CD4_low.down <- CD4_low.down[order(CD4_low.down[,"avg_logFC"],decreasing = F),]
  xlsx::write.xlsx(CD4_low.up, file = "result/without_PVTT_MLN/tumor/CD4T/CD4.low.tumor_vs_nonTumor.DEmarkers.xlsx",sheetName = "UP")
  xlsx::write.xlsx(CD4_low.down, file = "result/without_PVTT_MLN/tumor/CD4T/CD4.low.tumor_vs_nonTumor.DEmarkers.xlsx",sheetName = "Down",append = T)
  
  CD8_high.up <- CD8_high.markers[CD8_high.markers$avg_logFC > 0,]
  CD8_high.up <- CD8_high.up[order(CD8_high.up[,"avg_logFC"],decreasing = T),]
  CD8_high.down <- CD8_high.markers[CD8_high.markers$avg_logFC < 0,]
  CD8_high.down <- CD8_high.down[order(CD8_high.down[,"avg_logFC"],decreasing = F),]
  xlsx::write.xlsx(CD8_high.up,file = "result/without_PVTT_MLN/tumor/CD8T/CD8.high.tumor_vs_nonTumor.DEmarkers.xlsx",sheetName = "UP")
  xlsx::write.xlsx(CD8_high.down,file = "result/without_PVTT_MLN/tumor/CD8T/CD8.high.tumor_vs_nonTumor.DEmarkers.xlsx",sheetName = "Down",append = T)
  
  CD8_low.up <- CD8_low.markers[CD8_low.markers$avg_logFC > 0, ]
  CD8_low.up <- CD8_low.up[order(CD8_low.up[,"avg_logFC"],decreasing = T),]
  CD8_low.down <- CD8_low.markers[CD8_low.markers$avg_logFC < 0, ]
  CD8_low.down <- CD8_low.down[order(CD8_low.down[,"avg_logFC"],decreasing = F),]
  xlsx::write.xlsx(CD8_low.up, file = "result/without_PVTT_MLN/tumor/CD8T/CD8.low.tumor_vs_nonTumor.DEmarkers.xlsx",sheetName = "UP")
  xlsx::write.xlsx(CD8_low.down, file = "result/without_PVTT_MLN/tumor/CD8T/CD8.low.tumor_vs_nonTumor.DEmarkers.xlsx",sheetName = "Down",append = T)
  
  # save in one RData file
  save(
    CD3_high.up, CD3_high.down, CD3_low.up, CD3_low.down,
    CD4_high.up, CD4_high.down, CD4_low.up, CD4_low.down,
    CD8_high.up, CD8_high.down, CD8_low.up, CD8_low.down,
    file = "rds/without_PVTT_MLN/tumor/sub_T/2.Tumor_vs_nonTumor_up_down_Markers.RData"
  )
}


if(F){
  load("rds/without_PVTT_MLN/tumor/sub_T/4.tumor_Tcell_nontumor_removal.RData")
  # 1.比较不同肿瘤细胞间差异markers
  Idents(tumor_CD8_low_nontumor_removal) <- tumor_CD8_low_nontumor_removal$CellType
  CD8.low.nonTumor.celltype.markers <- FindAllMarkers(tumor_CD8_low_nontumor_removal,only.pos = FALSE, min.diff.pct = 0.25, logfc.threshold = 0.25)
  Idents(tumor_CD8_high_nontumor_removal) <- tumor_CD8_high_nontumor_removal$CellType
  CD8.high.nonTumor.celltype.markers <- FindAllMarkers(tumor_CD8_high_nontumor_removal,only.pos = FALSE, min.diff.pct = 0.25, logfc.threshold = 0.25)
  Idents(tumor_CD3_low_nontumor_removal) <- tumor_CD3_low_nontumor_removal$CellType
  CD3.low.nonTumor.celltype.markers <- FindAllMarkers(tumor_CD3_low_nontumor_removal,only.pos = FALSE, min.diff.pct = 0.25, logfc.threshold = 0.25)
  Idents(tumor_CD3_high_nontumor_removal) <- tumor_CD3_high_nontumor_removal$CellType
  CD3.high.nonTumor.celltype.markers <- FindAllMarkers(tumor_CD3_high_nontumor_removal,only.pos = FALSE, min.diff.pct = 0.25, logfc.threshold = 0.25)
  Idents(tumor_CD4_low_nontumor_removal) <- tumor_CD4_low_nontumor_removal$CellType
  CD4.low.nonTumor.celltype.markers <- FindAllMarkers(tumor_CD4_low_nontumor_removal,only.pos = FALSE, min.diff.pct = 0.25, logfc.threshold = 0.25)
  Idents(tumor_CD4_high_nontumor_removal) <- tumor_CD4_high_nontumor_removal$CellType
  CD4.high.nonTumor.celltype.markers <- FindAllMarkers(tumor_CD4_high_nontumor_removal,only.pos = FALSE, min.diff.pct = 0.25, logfc.threshold = 0.25)
  save(
    CD8.low.nonTumor.celltype.markers,
    CD8.high.nonTumor.celltype.markers,
    CD3.low.nonTumor.celltype.markers,
    CD3.high.nonTumor.celltype.markers,
    CD4.low.nonTumor.celltype.markers,
    CD4.high.nonTumor.celltype.markers,
    file = "rds/without_PVTT_MLN/tumor/sub_T/tumor_nonTumor_celltype_DEmarkers.RData"
  )
  # 2.不同病人样本间差异Markers
  Idents(tumor_CD8_low_nontumor_removal) <- tumor_CD8_low_nontumor_removal$Sample
  CD8.low.nonTumor.patient.markers <- FindAllMarkers(tumor_CD8_low_nontumor_removal,only.pos = FALSE, min.diff.pct = 0.25, logfc.threshold = 0.25)
  Idents(tumor_CD8_high_nontumor_removal) <- tumor_CD8_high_nontumor_removal$Sample
  CD8.high.nonTumor.patient.markers <- FindAllMarkers(tumor_CD8_high_nontumor_removal,only.pos = FALSE, min.diff.pct = 0.25, logfc.threshold = 0.25)
  Idents(tumor_CD3_low_nontumor_removal) <- tumor_CD3_low_nontumor_removal$Sample
  CD3.low.nonTumor.patient.markers <- FindAllMarkers(tumor_CD3_low_nontumor_removal,only.pos = FALSE, min.diff.pct = 0.25, logfc.threshold = 0.25)
  Idents(tumor_CD3_high_nontumor_removal) <- tumor_CD3_high_nontumor_removal$Sample
  CD3.high.nonTumor.patient.markers <- FindAllMarkers(tumor_CD3_high_nontumor_removal,only.pos = FALSE, min.diff.pct = 0.25, logfc.threshold = 0.25)
  Idents(tumor_CD4_low_nontumor_removal) <- tumor_CD4_low_nontumor_removal$Sample
  CD4.low.nonTumor.patient.markers <- FindAllMarkers(tumor_CD4_low_nontumor_removal,only.pos = FALSE, min.diff.pct = 0.25, logfc.threshold = 0.25)
  Idents(tumor_CD4_high_nontumor_removal) <- tumor_CD4_high_nontumor_removal$Sample
  CD4.high.nonTumor.patient.markers <- FindAllMarkers(tumor_CD4_high_nontumor_removal,only.pos = FALSE, min.diff.pct = 0.25, logfc.threshold = 0.25)
  save(
    CD8.low.nonTumor.patient.markers,
    CD8.high.nonTumor.patient.markers,
    CD3.low.nonTumor.patient.markers,
    CD3.high.nonTumor.patient.markers,
    CD4.low.nonTumor.patient.markers,
    CD4.high.nonTumor.patient.markers,
    file = "rds/without_PVTT_MLN/tumor/sub_T/tumor_nonTume_patient.DEmarkers.RData"
  )
}

##>>>>>>>>>>>>>>>>>>>> DE gene GSEA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if(FALSE){
  rm(list = ls());gc()
  tumor <- readRDS("rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR.rds")
  load("rds/without_PVTT_MLN/tumor/sub_T/1.CDT_high_low_DEmarkers.RData")
  # CD3T high vs CD3T low UP
  expr <- GetAssayData(tumor, slot = "counts")[rownames(CD3.markers.up),]
  expr <- data.frame(NAME = rownames(expr), Description = rep("na",nrow(expr)),expr,stringsAsFactors = F)
  write("#1.2","result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.up.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.up.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.up.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor@meta.data$Group1)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.up.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.up.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.up.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  # CD3T high vs CD3T low Down
  expr <- GetAssayData(tumor, slot = "counts")[rownames(CD3.markers.down),]
  expr <- data.frame(NAME = rownames(expr), Description = rep("na",nrow(expr)),expr,stringsAsFactors = F)
  write("#1.2","result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.down.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.down.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.down.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor@meta.data$Group1)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.down.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.down.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.down.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  # CD4T high vs CD4T low UP
  expr <- GetAssayData(tumor, slot = "counts")[rownames(CD4.markers.up),]
  expr <- data.frame(NAME = rownames(expr), Description = rep("na",nrow(expr)),expr,stringsAsFactors = F)
  write("#1.2","result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.up.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.up.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.up.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor@meta.data$Group2)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.up.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.up.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.up.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  # CD4T high vs CD4T low Down
  expr <- GetAssayData(tumor, slot = "counts")[rownames(CD4.markers.down),]
  expr <- data.frame(NAME = rownames(expr), Description = rep("na",nrow(expr)),expr,stringsAsFactors = F)
  write("#1.2","result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.down.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.down.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.down.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor@meta.data$Group2)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.down.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.down.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.down.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  # CD8T high vs CD8T low UP
  expr <- GetAssayData(tumor, slot = "counts")[rownames(CD8.markers.up),]
  expr <- data.frame(NAME = rownames(expr), Description = rep("na",nrow(expr)),expr,stringsAsFactors = F)
  write("#1.2","result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.up.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.up.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.up.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor@meta.data$Group3)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.up.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.up.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.up.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  # CD8T high vs CD8T low Down
  expr <- GetAssayData(tumor, slot = "counts")[rownames(CD8.markers.down),]
  expr <- data.frame(NAME = rownames(expr), Description = rep("na",nrow(expr)),expr,stringsAsFactors = F)
  write("#1.2","result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.down.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.down.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.down.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor@meta.data$Group3)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.down.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.down.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.down.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  
  ## Tumor vs Nontumor inner each CD high/low
  rm(list = ls());gc()
  load("rds/without_PVTT_MLN/tumor/3.tumor_Tcell_tumor_marked.RData")
  load("rds/without_PVTT_MLN/tumor/sub_T/2.Tumor_vs_nonTumor_up_down_Markers.RData")
  
  ## CD8T_low.up
  expr <- GetAssayData(tumor_CD8_low,slot = "counts")[rownames(CD8_low.up),]
  expr <- data.frame(NAME = rownames(expr),Description = rep('na',nrow(expr)),expr,stringsAsFactors = F)
  write('#1.2',"result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.up.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.up.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.up.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor_CD8_low@meta.data$Class)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.up.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.up.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.up.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  ## CD8T_low.down
  expr <- GetAssayData(tumor_CD8_low,slot = "counts")[rownames(CD8_low.down),]
  expr <- data.frame(NAME = rownames(expr),Description = rep('na',nrow(expr)),expr,stringsAsFactors = F)
  write('#1.2',"result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.down.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.down.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.down.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor_CD8_low@meta.data$Class)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.down.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.down.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.down.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  ## CD8T_high.up
  expr <- GetAssayData(tumor_CD8_high,slot = "counts")[rownames(CD8_high.up),]
  expr <- data.frame(NAME = rownames(expr),Description = rep('na',nrow(expr)),expr,stringsAsFactors = F)
  write('#1.2',"result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.up.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.up.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.up.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor_CD8_high@meta.data$Class)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.up.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.up.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.up.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  ## CD8T_high.down
  expr <- GetAssayData(tumor_CD8_high,slot = "counts")[rownames(CD8_high.down),]
  expr <- data.frame(NAME = rownames(expr),Description = rep('na',nrow(expr)),expr,stringsAsFactors = F)
  write('#1.2',"result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.down.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.down.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.down.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor_CD8_high@meta.data$Class)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.down.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.down.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.down.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  
  ## CD4T_low.up
  expr <- GetAssayData(tumor_CD4_low,slot = "counts")[rownames(CD4_low.up),]
  expr <- data.frame(NAME = rownames(expr),Description = rep('na',nrow(expr)),expr,stringsAsFactors = F)
  write('#1.2',"result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.up.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.up.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.up.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor_CD4_low@meta.data$Class)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.up.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.up.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.up.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  ## CD4T_low.down
  expr <- GetAssayData(tumor_CD4_low,slot = "counts")[rownames(CD4_low.down),]
  expr <- data.frame(NAME = rownames(expr),Description = rep('na',nrow(expr)),expr,stringsAsFactors = F)
  write('#1.2',"result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.down.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.down.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.down.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor_CD4_low@meta.data$Class)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.down.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.down.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.down.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  ## CD4T_high.up
  expr <- GetAssayData(tumor_CD4_high,slot = "counts")[rownames(CD4_high.up),]
  expr <- data.frame(NAME = rownames(expr),Description = rep('na',nrow(expr)),expr,stringsAsFactors = F)
  write('#1.2',"result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.up.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.up.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.up.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor_CD4_high@meta.data$Class)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.up.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.up.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.up.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  ## CD4T_high.down
  expr <- GetAssayData(tumor_CD4_high,slot = "counts")[rownames(CD4_high.down),]
  expr <- data.frame(NAME = rownames(expr),Description = rep('na',nrow(expr)),expr,stringsAsFactors = F)
  write('#1.2',"result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.down.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.down.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.down.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor_CD4_high@meta.data$Class)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.down.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.down.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.down.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  
  ## CD3T_low.up
  expr <- GetAssayData(tumor_CD3_low,slot = "counts")[rownames(CD3_low.up),]
  expr <- data.frame(NAME = rownames(expr),Description = rep('na',nrow(expr)),expr,stringsAsFactors = F)
  write('#1.2',"result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.up.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.up.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.up.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor_CD3_low@meta.data$Class)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.up.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.up.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.up.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  ## CD3T_low.down
  expr <- GetAssayData(tumor_CD3_low,slot = "counts")[rownames(CD3_low.down),]
  expr <- data.frame(NAME = rownames(expr),Description = rep('na',nrow(expr)),expr,stringsAsFactors = F)
  write('#1.2',"result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.down.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.down.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.down.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor_CD3_low@meta.data$Class)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.down.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.down.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.down.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  ## CD3T_high.up
  expr <- GetAssayData(tumor_CD3_high,slot = "counts")[rownames(CD3_high.up),]
  expr <- data.frame(NAME = rownames(expr),Description = rep('na',nrow(expr)),expr,stringsAsFactors = F)
  write('#1.2',"result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.up.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.up.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.up.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor_CD3_high@meta.data$Class)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.up.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.up.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.up.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  ## CD3T_high.down
  expr <- GetAssayData(tumor_CD3_high,slot = "counts")[rownames(CD3_high.down),]
  expr <- data.frame(NAME = rownames(expr),Description = rep('na',nrow(expr)),expr,stringsAsFactors = F)
  write('#1.2',"result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.down.expr.gct",ncolumns = 1)
  write(c(nrow(expr),(ncol(expr)-2)), "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.down.expr.gct",ncolumns = 2,append = T,sep = "\t")
  write.table(expr, "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.down.expr.gct",row.names = F,sep = "\t",append = T,quote = F)
  line.1 <- c((ncol(expr) - 2), 2, 1)
  tmp <- table(tumor_CD3_high@meta.data$Class)
  line.2 <- c("#", names(tmp))
  line.3 <- c(rep(names(tmp)[1],tmp[1]), rep(names(tmp)[2], tmp[2]))
  write(line.1, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.down.group.cls',ncolumns = length(line.1),append = T,sep = "\t")
  write(line.2, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.down.group.cls',ncolumns = length(line.2),append = T,sep = "\t")
  write(line.3, 'result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.down.group.cls',ncolumns = length(line.3),append = T,sep = "\t")
  
  
  library(future)
  plan("multicore",workers = 50)
  # CD8+T.low.up
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.up.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.up.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD8T/output/low/up/")
  # CD8+T.low.down
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.down.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_low.down.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD8T/output/low/down/")
  # CD8+T.high.up
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.up.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.up.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD8T/output/high/up/")
  # CD8+T.high.down
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.down.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8_high.down.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD8T/output/high/down/")
  # CD3+T.low.up
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.up.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.up.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD3T/output/low/up/")
  # CD3+T.low.down
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.down.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_low.down.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD3T/output/low/down/")
  # CD3+T.high.up
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.up.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.up.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD3T/output/high/up/")
  # CD3+T.high.down
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.down.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3_high.down.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD3T/output/high/down/")
  # CD4+T.low.up
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.up.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.up.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD4T/output/low/up/")
  # CD4+T.low.down
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.down.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_low.down.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD4T/output/low/down/")
  # CD4+T.high.up
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.up.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.up.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD4T/output/high/up/")
  # CD4+T.high.down
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.down.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4_high.down.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD4T/output/high/down/")
  
  # CD3+T UP
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.up.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.up.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD3T/output/up/")
  # CD3+T Down
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.down.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD3T/input/CD3.down.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD3T/output/down/")
  
  # CD4+T UP
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.up.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.up.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD4T/output/up/")
  # CD4+T Down
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.down.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD4T/input/CD4.down.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD4T/output/down/")
  
  # CD8+T UP
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.up.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.up.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD8T/output/up/")
  # CD8+T Down
  GSEA::GSEA("result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.down.expr.gct",
             "result/without_PVTT_MLN/tumor/GSEA/CD8T/input/CD8.down.group.cls",
             gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "result/without_PVTT_MLN/tumor/GSEA/CD8T/output/down/")
}


# 批量做GSEA分析

if(FALSE){
  genelist <- CD8_low.up$avg_logFC
  names(genelist) <- toupper(rownames(CD8_low.up))
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GSEABase)
  geneset <- read.gmt("gsea/h.all.v7.2.symbols.gmt")
  egmt <- GSEA(genelist, TERM2GENE = geneset, minGSSize = 1, pvalueCutoff = 0.99, verbose = F)
  {
    genelist <- CD8_low.up$avg_logFC
    names(genelist) <- toupper(rownames(CD8_low.up))
    genelist <- sort(genelist, decreasing = T)
    library(ggplot2)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(GSEABase)
    gmtfile <- "gsea/h.all.v7.2.symbols.gmt"
    geneset <- read.gmt(gmtfile)
    egmt <- GSEA(genelist, TERM2GENE = geneset, minGSSize = 1, pvalueCutoff = 0.99, verbose = F)
    gsea_result_df <- egmt@result
    write.csv(gsea_result_df, file = "result/without_PVTT_MLN/tumor/CD8T/CD8_low.gsea.results.csv")
    library(enrichplot)
    # gseaplot2(egmt, geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',pvalue_table = T)
    # gseaplot2(egmt, geneSetID = 'HALLMARK_MTORC1_SIGNALING',pvalue_table = T)
  }
  
  # 把全部的基因集的GSEA分析结果批量出图
  # down_kegg <- gsea_result_df[gsea_result_df$pvalue<0.05 & gsea_result_df$enrichmentScore < 0,];down_kegg$group <- -1
  # up_kegg <- gsea_result_df[gsea_result_df$pvalue<0.05 & gsea_result_df$enrichmentScore > 0.3,];up_kegg$group <- 1
  
  down_kegg <- gsea_result_df[gsea_result_df$enrichmentScore < 0,];down_kegg$group <- -1
  up_kegg <- gsea_result_df[gsea_result_df$enrichmentScore > 0,];up_kegg$group <- 1
  
  library(enrichplot)
  lapply(1:nrow(down_kegg), function(i){
    gseaplot2(egmt, down_kegg$ID[i],title = down_kegg$Description[i],pvalue_table = T,color = "red")
    ggsave(paste0("down_kegg_",gsub("/","-",down_kegg$Description[i]),".pdf"),width = 10,height = 6,path = "result/without_PVTT_MLN/tumor/CD8T/")
  })
  lapply(1:nrow(up_kegg), function(i){
    gseaplot2(egmt,up_kegg$ID[i],title = up_kegg$Description[i],pvalue_table = T,color = "red")
    ggsave(paste0("up_kegg_",gsub("/","-",up_kegg$Description[i]),".pdf"),width = 10,height = 6,path = "gsea/test")
  })
}


##====================== MDSC ============================================================

# CD45+:   PTPRC
# CD3-:  CD3D, CD3E, CD3G
# B220-:  PTPRC
# NK1.1-: KLRB1
# gdTCR-: 
# HLA-DR-
# CD11b+: ITGAM
# CD14-: CD14
# CD15+: FUT4
# CD33+: CD33
# CD66b+: CEACAM8

rm(list = ls());gc()
tumor <- readRDS("rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR.rds")

# 首先根据CD45+, CD3-, B220-, NK1.1-, gdTCR-, HLA-DR-, CD11b+ CD33+将MDSC细胞区分出来
# 实际用PTPRC,CD3D,CD3E,CD3G,KLRB1,CD33,
markers <- c("PTPRC","CD3D","CD3E","CD3G","KLRB1","CD33")
source("script/StackedVlnPlot.R")

pal <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
                  '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#F3C300','#222222','#828281')
p1 <- DimPlot(tumor,reduction = "tsne",  group.by = "seurat_clusters",cols = pal)
ggsave("result/without_PVTT_MLN/tumor/MDSC/0.1.tumor_tsne.pdf",p1,width = 6,height = 5)

p2 <- StackedVlnPlot(obj = tumor, features = markers, pt.size = 0.0001, group.by = "seurat_clusters")
ggsave("result/without_PVTT_MLN/tumor/MDSC/0.2.CD45_CD3_NK1.1_CD33.pdf",p2,width = 8,height = 10)

tumor$CellType3 <- ifelse(
  tumor$seurat_clusters == "0" | tumor$seurat_clusters == "16" | tumor$seurat_clusters == "20",
  "MDSC","Others"
)
p3 <- DimPlot(tumor,reduction = "tsne",group.by = "CellType3", cols = c("#BF0032","#0066A5")) + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        legend.text = element_text(size = 14))
ggsave("result/without_PVTT_MLN/tumor/MDSC/0.3.MDSC_vs_others_tsne.pdf",p3,width = 5,height = 5)

p3.1 <- DimPlot(tumor,reduction = "tsne", group.by = "CellType3", cells.highlight = rownames(subset(tumor@meta.data,CellType3=="MDSC")),
                cols.highlight = "#BF0032") + ggtitle("MDSC") +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        legend.text = element_text(size = 14))
ggsave("result/without_PVTT_MLN/tumor/MDSC/0.3.MDSC_vs_others_tsne_mod.pdf",p3.1,width = 5,height = 5)

MDSC_cells <- subset(tumor@meta.data, CellType3 == "MDSC")
MDSC <- subset(tumor, cells = rownames(MDSC_cells))

MDSC <- FindVariableFeatures(MDSC, selection.method = "vst",nfeatures = 2000)
DefaultAssay(MDSC) <- "integrated"
MDSC <- ScaleData(MDSC,features = rownames(MDSC))
MDSC <- RunPCA(MDSC, features = VariableFeatures(MDSC))
p1 <- DimPlot(MDSC, reduction = "pca", group.by = "orig.ident",cols=rainbow(length(levels(as.factor(MDSC$orig.ident)))))
ggsave("result/without_PVTT_MLN/tumor/MDSC/1.MDSC_pca.pdf",p1,device="pdf",width=8,height=6)
p2 <- ElbowPlot(MDSC,ndims=30,reduction="pca")
ggsave("result/without_PVTT_MLN/tumor/MDSC/2.MDSC_elbow.pdf",p2,device="pdf",width=6,height=4)

pc.num = 1:24
# 聚类
MDSC <- FindNeighbors(MDSC, dims = pc.num) %>% FindClusters(resolution = 0.2)
MDSC <- RunUMAP(MDSC,reduction="pca",dims=pc.num) %>% RunTSNE(dims=pc.num)

pal <- c("#786140","#8C0603","#B2712B","#B48521","#C6EEC4","#DCD3B1","#A7BBA5","#0202BB","#5EC6B9","#4E8A97","#726DBA")

p1 <- DimPlot(MDSC,reduction="tsne",cols=pal)
p2 <- DimPlot(MDSC,reduction="tsne",group.by="orig.ident",cols=rainbow(length(levels(as.factor(MDSC$orig.ident))))) + theme(plot.title = element_blank())
pc <- p1 + p2
ggsave("result/without_PVTT_MLN/tumor/MDSC/3.MDSC_tsne.pdf",pc,device = "pdf",width = 14,height = 4.5)

DefaultAssay(MDSC) <- "RNA"
# percent mt ribo
MDSC[['percent.mt']] <- PercentageFeatureSet(MDSC, pattern="^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(MDSC@assays$RNA)) 
HB.genes <- rownames(MDSC@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)]
MDSC[['percent.hb']] <- PercentageFeatureSet(MDSC, features=HB.genes)
MDSC$project <- "MDSC"

violin <- VlnPlot(MDSC, group.by = "project",
                  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                  pt.size = 0.00001, ncol = 3) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/MDSC/4.MDSC_violin.pdf",violin,device = "pdf",width = 8,height = 6)


Idents(MDSC) <- MDSC$Sample
DefaultAssay(MDSC) <- "RNA"
table(MDSC$seurat_clusters,MDSC$Type)

# G-MDSC: CD14– CD15+ CD66b+
# M-MSDC: CD14+ CD15– CD66b–

# CD14: CD14, CD15: FUT4, CD66b: CEACAM8
CD14 <- VlnPlot(MDSC, features = "CD14",cols = pal,pt.size = 0,group.by = "seurat_clusters") +
  coord_flip() + ggtitle("CD14") +
  theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
        panel.border = element_rect(size = 1.5,colour = "black"))
ggsave("result/without_PVTT_MLN/tumor/MDSC/5.MDSC_CD14.pdf",width = 4,height = 6)

# cell type assignment

current_cluster <- c(0:10)

# classfied CD3+T cell into CD4+T and CD8+T cells
new_cluster <- c("M-MDSC","M-MDSC","M-MDSC","M-MDSC","M-MDSC","M-MDSC",
                 "M-MDSC","M-MDSC","M-MDSC","G-MDSC","G-MDSC")
MDSC$subtype <- plyr::mapvalues(MDSC$seurat_clusters,from = current_cluster,to = new_cluster)
saveRDS(MDSC,file = "rds/without_PVTT_MLN/tumor/MDSC/MDSC_celltypeAssigned.rds")

p1 <- DimPlot(MDSC,reduction = "tsne",group.by = "seurat_clusters",cols = pal)
p2 <- DimPlot(MDSC,reduction = "tsne",group.by = "Sample",cols = rainbow(length(levels(as.factor(MDSC$Sample)))))
p3 <- DimPlot(MDSC,reduction = "tsne",group.by = "subtype")

p3_mod <- DimPlot(MDSC,reduction = "tsne",group.by = "subtype",cols = c("#BF0032","#0066A5")) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        legend.text = element_text(size = 14))
ggsave("result/without_PVTT_MLN/tumor/MDSC/6.G.MDSC_M.MDSC.tsne.pdf",p3_mod,device = "pdf",width = 4,height = 4)
pc <- p1 + p2 + p3
ggsave("result/without_PVTT_MLN/tumor/MDSC/6.MDSC_tsne.pdf",pc,device = "pdf",width = 16,height = 4)


M_MDSC <- DimPlot(MDSC,reduction = "tsne",group.by = "subtype",cells.highlight = rownames(MDSC@meta.data[MDSC$subtype=="M-MDSC",])) + 
  NoLegend() + ggtitle("M-MDSC") + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/MDSC/7.M_MDSC_tsne.pdf",M_MDSC,device = "pdf",width = 4,height = 4)

G_MDSC <- DimPlot(MDSC,reduction = "tsne",group.by = "subtype",cells.highlight = rownames(MDSC@meta.data[MDSC$subtype=="G-MDSC",])) + 
  NoLegend() + ggtitle("G-MDSC") + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/MDSC/8.G-MDSC_tsne.pdf",G_MDSC,device = "pdf",width = 4,height = 4)


###### Calculating Infiltration Rate
rm(list = ls());gc()
tumor <- readRDS("rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR.rds")
MDSC <- readRDS(MDSC,file = "rds/without_PVTT_MLN/tumor/MDSC/MDSC_celltypeAssigned.rds")

# T cell infiltration rate
MDSC_total <- cbind(table(tumor$Sample,tumor$CellType),table(MDSC$Sample,MDSC$subtype))
MDSC_total <- t(MDSC_total)
M.MDSC.IR <- NULL
G.MDSC.IR <- NULL
## M.MDSC
for(i in seq_len(ncol(MDSC_total))){M.MDSC.IR[i] <- MDSC_total[2,i]/colSums(MDSC_total[1:7,])[i]}
## G.MDSC
for(i in seq_len(ncol(MDSC_total))){G.MDSC.IR[i] <- MDSC_total[8,i]/colSums(MDSC_total[1:7,])[i]}

total <- rbind(MDSC_total,M.MDSC.IR,G.MDSC.IR)
xlsx::write.xlsx(total,file = "result/without_PVTT_MLN/tumor/MDSC/MDSC_Infiltration_rate.xlsx")

MDSC.IR <- data.frame(t(total[10:11,]))
MDSC.IR$Patient <- rownames(MDSC.IR)
save(MDSC.IR,file = "rds/without_PVTT_MLN/tumor/MDSC/MDSC.IR.RData")

# add infiltration rate to tumor metadata
tumor$M.MDSC.IR <- plyr::mapvalues(tumor$Sample, from = MDSC.IR$Patient, to = MDSC.IR$M.MDSC.IR)
tumor$G.MDSC.IR <- plyr::mapvalues(tumor$Sample, from = MDSC.IR$Patient, to = MDSC.IR$G.MDSC.IR)
saveRDS(tumor,file = "rds/without_PVTT_MLN/tumor/4.tumor_combine_T.IR_MDSC.IR.rds")


library(ggpubr)
MDSC.IR <- MDSC.IR[order(MDSC.IR[,1],decreasing = T),]
M.MDSC_IR <- ggbarplot(MDSC.IR,x = "Patient", y = "M.MDSC.IR",fill = "Patient") + 
  NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 14)) + 
  geom_hline(aes(yintercept=mean(MDSC.IR$M.MDSC.IR)),color="#BF0032",linetype="dashed") +
  geom_hline(aes(yintercept=median(MDSC.IR$M.MDSC.IR)),color="#0066A5",linetype="dashed") +
  geom_hline(aes(yintercept=quantile(MDSC.IR$M.MDSC.IR,0.75)),color="#222222",linetype="dashed") +
  geom_hline(aes(yintercept=quantile(MDSC.IR$M.MDSC.IR,0.25)),color="#222222",linetype="dashed") +
  annotate("text",x=20.5,y=mean(MDSC.IR$M.MDSC.IR)+0.01,label="Mean:0.07",color="#BF0032",size=3) +
  annotate("text",x=18,y=median(MDSC.IR$M.MDSC.IR)+0.01,label="Median:0.02",color="#0066A5",size=3) +
  annotate("text",x=20,y=quantile(MDSC.IR$M.MDSC.IR,0.75)-0.01,label="Upper Quantile:0.06",color="#222222",size=3) +
  annotate("text",x=20,y=quantile(MDSC.IR$M.MDSC.IR,0.01)+0.016,label="Lower Quantile:0.01",color="#222222",size=3)
ggsave("result/without_PVTT_MLN/tumor/MDSC/9.M.MDSC.IR.pdf",M.MDSC_IR,width = 8,height = 4)

MDSC.IR <- MDSC.IR[order(MDSC.IR[,2],decreasing = T),]
G.MDSC_IR <- ggbarplot(MDSC.IR,x = "Patient", y = "G.MDSC.IR",fill = "Patient") + 
  NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 14)) + 
  geom_hline(aes(yintercept=mean(MDSC.IR$G.MDSC.IR)),color="#BF0032",linetype="dashed") +
  geom_hline(aes(yintercept=median(MDSC.IR$G.MDSC.IR)),color="#0066A5",linetype="dashed") +
  geom_hline(aes(yintercept=quantile(MDSC.IR$G.MDSC.IR,0.75)),color="#222222",linetype="dashed") +
  geom_hline(aes(yintercept=quantile(MDSC.IR$G.MDSC.IR,0.25)),color="#222222",linetype="dashed") +
  annotate("text",x=20.5,y=mean(MDSC.IR$G.MDSC.IR)+0.02,label="Mean:0.15",color="#BF0032",size=3) +
  annotate("text",x=20.5,y=median(MDSC.IR$G.MDSC.IR)-0.03,label="Median:0.13",color="#0066A5",size=3) +
  annotate("text",x=20,y=quantile(MDSC.IR$G.MDSC.IR,0.75)+0.03,label="Upper Quantile:0.18",color="#222222",size=3) +
  annotate("text",x=20,y=quantile(MDSC.IR$G.MDSC.IR,0.25)+0.02,label="Lower Quantile:0.03",color="#222222",size=3)

ggsave("result/without_PVTT_MLN/tumor/MDSC/10.G.MDSC.IR.pdf",G.MDSC_IR,width = 8,height = 4)


tumor$CellType <- as.character(tumor$CellType)

# 设定阈值，区分High和low组
M.MDSC_IR_threshold = NA
G.MDSC_IR_threshold = NA

tumor$Group4 <- ifelse(tumor$M.MDSC.IR > M.MDSC_IR_threshold, "M-MDSC high", "M-MDSC low")
tumor$Group5 <- ifelse(tumor$G.MDSC.IR > G.MDSC_IR_threshold, "G-MDSC high", "G-MDSC high")

## M-MDSC cell
tumor_M.MDSC_high_cells <- subset(tumor@meta.data, Group4 == "M-MDSC high")
tumor_M.MDSC_high <- subset(tumor, cells = rownames(tumor_M.MDSC_high_cells))
tumor_M.MDSC_low_cells <- subset(tumor@meta.data, Group4 == "M-MDSC low")
tumor_M.MDSC_low <- subset(tumor, cells = rownames(tumor_M.MDSC_low_cells))

## G-MDSC cell
tumor_G.MDSC_high_cells <- subset(tumor@meta.data, Group5 == "G-MDSC high")
tumor_G.MDSC_high <- subset(tumor, cells = rownames(tumor_G.MDSC_high_cells))
tumor_G.MDSC_low_cells <- subset(tumor@meta.data, Group5 == "G-MDSC low")
tumor_G.MDSC_low <- subset(tumor, cells = rownames(tumor_G.MDSC_low_cells))

# classify tumor and nontumor cells
table(tumor_G.MDSC_high$CellType)
p5 <- VlnPlot(tumor_G.MDSC_high, features = "PTPRC",cols = pal,pt.size = 0,group.by = "CellType") +
  coord_flip() + ggtitle("CD45") + NoLegend() +
  theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14))
ggsave("result/without_PVTT_MLN/tumor/MDSC/CD45_tumor_G.MDSC_high.pdf", p5, device = "pdf", width = 4,height = 6)

tumor_G.MDSC_high$Class <- ifelse(
  tumor_G.MDSC_high$CellType == "TEC" | tumor_G.MDSC_high$CellType == "CAF" | tumor_G.MDSC_high$CellType == "Malignant cell",
  "Tumor", "Nontumor"
)


table(tumor_G.MDSC_low$CellType)
p6 <- VlnPlot(tumor_G.MDSC_low, features = "PTPRC",cols = pal,pt.size = 0,group.by = "CellType") +
  coord_flip() + ggtitle("CD45") + NoLegend() +
  theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14))
ggsave("result/without_PVTT_MLN/tumor/MDSC/CD45_tumor_G.MDSC_low.pdf", p6, device = "pdf", width = 4,height = 6)

tumor_G.MDSC_low$Class <- ifelse(
  tumor_G.MDSC_low$CellType == "TEC" | tumor_G.MDSC_low$CellType == "CAF" | tumor_G.MDSC_low$CellType == "Malignant cell",
  "Tumor", "Nontumor"
)


table(tumor_M.MDSC_high$CellType)
p5 <- VlnPlot(tumor_M.MDSC_high, features = "PTPRC",cols = pal,pt.size = 0,group.by = "CellType") +
  coord_flip() + ggtitle("CD45") + NoLegend() +
  theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14))
ggsave("result/without_PVTT_MLN/tumor/MDSC/CD45_tumor_M.MDSC_high.pdf", p5, device = "pdf", width = 4,height = 6)

tumor_M.MDSC_high$Class <- ifelse(
  tumor_M.MDSC_high$CellType == "TEC" | tumor_M.MDSC_high$CellType == "CAF" | tumor_M.MDSC_high$CellType == "Malignant cell",
  "Tumor", "Nontumor"
)

table(tumor_M.MDSC_low$CellType)
p6 <- VlnPlot(tumor_M.MDSC_low, features = "PTPRC",cols = pal,pt.size = 0,group.by = "CellType") +
  coord_flip() + ggtitle("CD45") + NoLegend() +
  theme(axis.text.x = element_text(angle = 0,hjust = 1,vjust = 1),
        panel.border = element_rect(size = 1.5,colour = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 14))
ggsave("result/without_PVTT_MLN/tumor/MDSC/CD45_tumor_M.MDSC_low.pdf", p6, device = "pdf", width = 4,height = 6)

tumor_M.MDSC_low$Class <- ifelse(
  tumor_M.MDSC_low$CellType == "TEC" | tumor_M.MDSC_low$CellType == "CAF" | tumor_M.MDSC_low$CellType == "Malignant cell",
  "Tumor", "Nontumor"
)

save(list(
  tumor_G.MDSC_high = tumor_G.MDSC_high,
  tumor_G.MDSC_low = tumor_G.MDSC_low,
  tumor_M.MDSC_high = tumor_M.MDSC_high,
  tumor_M.MDSC_low = tumor_M.MDSC_low
),file = "rds/without_PVTT_MLN/MDSC/tumor_MDSC_tumor_marked.RData")



## DE genes