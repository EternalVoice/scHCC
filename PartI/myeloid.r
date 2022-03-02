library(Seurat)
library(ggplot2)
library(dplyr)

tumor <- readRDS("rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR_2.rds")

source("script/PartI/StackedVlnPlot.R")
### Myeloid cell

# CD45+, CD11B+, CD3-, CD19-, CD56-, CD16-

markers <- c("PTPRC","ITGAM","CD3D","CD3E","CD3G","CD19","NCAM1","FCGR3A")

p_CD45 <- VlnPlot(tumor,features = "PTPRC",group.by = "seurat_clusters") + NoLegend() + 
  theme(axis.title = element_blank())
p_CD11B <- VlnPlot(tumor,features = "ITGAM",group.by = "seurat_clusters") + NoLegend() +
  theme(axis.title = element_blank())
p_CD3D <- VlnPlot(tumor,features = "CD3D",group.by = "seurat_clusters",pt.size = 0) + NoLegend() +
  theme(axis.title = element_blank())
p_CD3E <- VlnPlot(tumor,features = "CD3E",group.by = "seurat_clusters",pt.size = 0) + NoLegend() +
  theme(axis.title = element_blank())
p_CD3G <- VlnPlot(tumor,features = "CD3G",group.by = "seurat_clusters",pt.size = 0) + NoLegend() +
  theme(axis.title = element_blank())
p_CD19 <- VlnPlot(tumor,features = "CD19",group.by = "seurat_clusters",pt.size = 0) + NoLegend() +
  theme(axis.title = element_blank())
p_CD56 <- VlnPlot(tumor,features = "NCAM1",group.by = "seurat_clusters",pt.size = 0) + NoLegend() +
  theme(axis.title = element_blank())
p_CD16 <- VlnPlot(tumor,features = "FCGR3A",group.by = "seurat_clusters",pt.size = 0) + NoLegend() +
  theme(axis.title = element_blank())

ggsave("result/without_PVTT_MLN/tumor/myeloid/CD45.pdf",p_CD45,width = 8,height = 2.5)
ggsave("result/without_PVTT_MLN/tumor/myeloid/CD11B.pdf",p_CD11B,width = 8,height = 2.5)
ggsave("result/without_PVTT_MLN/tumor/myeloid/CD3D.pdf",p_CD3D,width = 8,height = 2.5)
ggsave("result/without_PVTT_MLN/tumor/myeloid/CD3E.pdf",p_CD3E,width = 8,height = 2.5)
ggsave("result/without_PVTT_MLN/tumor/myeloid/CD3G.pdf",p_CD3G,width = 8,height = 2.5)
ggsave("result/without_PVTT_MLN/tumor/myeloid/CD19.pdf",p_CD19,width = 8,height = 2.5)
ggsave("result/without_PVTT_MLN/tumor/myeloid/CD56.pdf",p_CD56,width = 8,height = 2.5)
ggsave("result/without_PVTT_MLN/tumor/myeloid/CD16.pdf",p_CD16,width = 8,height = 2.5)

# cluster 2-11,13
cells <- NULL
for(i in c(2:11)){
  cell <- rownames(subset(tumor@meta.data,seurat_clusters == i))
  cells <- c(cells,cell)
}
cells <- rownames(subset(tumor@meta.data,seurat_clusters == c(2:11) | seurat_clusters == 13))
cell <- rownames(subset(tumor@meta.data,seurat_clusters == 13))
cells <- c(cells,cell)

p_cluster <- DimPlot(tumor,reduction = "tsne",group.by = "seurat_clusters",label = TRUE,label.size = 5) + NoLegend() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 1.5,colour = "black"))
ggsave("result/without_PVTT_MLN/tumor/myeloid/seurat_clusters.pdf",p_cluster,width = 5.6,height = 5.7)

p_mye <- DimPlot(tumor,cells.highlight = cells,reduction = "tsne",pt.size = 0.7) + 
  NoLegend() + ggtitle("myeloid") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(size = 1.5,colour = "black"))
ggsave("result/without_PVTT_MLN/tumor/myeloid/mye.highlight.pdf",p_mye,width = 5.6,height = 5.7)

mye <- subset(tumor,cells = cells)
mye@meta.data <- mye@meta.data[,-c(6:27)]
saveRDS(mye,file = "rds/without_PVTT_MLN/tumor/mye.rds")

#
rm(list = ls());gc()
mye <- readRDS("rds/without_PVTT_MLN/tumor/mye.rds")

DefaultAssay(mye) <- "RNA"
# 降维聚类
mye <- FindVariableFeatures(mye,selection.method = "vst",nfeatures = 2000)
DefaultAssay(mye) <- "integrated"
mye <- ScaleData(mye,features = rownames(mye))
mye <- RunPCA(mye, features = VariableFeatures(mye))
p2 <- ElbowPlot(mye,ndims=30,reduction="pca")
mye$project <- "tumor"
pc.num = 1:26
# 聚类
mye <- FindNeighbors(mye, dims = pc.num) %>% FindClusters(resolution = 0.2)
mye <- RunUMAP(mye,reduction="pca",dims=pc.num) %>% RunTSNE(dims=pc.num)

p1 <- DimPlot(mye,reduction="tsne",cols=rainbow(length(levels(as.factor(mye$seurat_clusters)))))
p2 <- DimPlot(mye,reduction="tsne",group.by="orig.ident",cols=rainbow(length(levels(as.factor(mye$orig.ident))))) + theme(plot.title = element_blank())
pc <- p1 + p2
ggsave("result/without_PVTT_MLN/tumor/myeloid/clusters.pdf",pc,device = "pdf",width = 14,height = 4.5)

DefaultAssay(mye) <- "RNA"
# percent mt ribo
mye[['percent.mt']] <- PercentageFeatureSet(mye, pattern="^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(mye@assays$RNA)) 
HB.genes <- rownames(mye@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)]
mye[['percent.hb']] <- PercentageFeatureSet(mye, features=HB.genes)

violin <- VlnPlot(mye, group.by = "project",
                  features = c("nFeature_RNA","nCount_RNA","percent.mt"),
                  pt.size = 0.00001, ncol = 3) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/myeloid/myeloid.threeStd.pdf",violin,device = "pdf",width = 8,height = 6)

# Find Markers
markers <- FindAllMarkers(mye,only.pos = TRUE,logfc.threshold = 0.25)

pwd = "result/without_PVTT_MLN/tumor/PARTV/PARTV.2/DEmarkers/"
for(i in c(0:12)){
  if (!file.exists(paste0(pwd,paste0("c",i)))) {
    dir.create(paste0(pwd,paste0("c",i)))
  }
  markers <- FindMarkers(mye,ident.1 = i, only.pos = TRUE)
  for(j in seq_along(rownames(markers))){
    p <- FeaturePlot(mye,features = rownames(markers)[j],pt.size = 0.78,cols = c("lightgrey","#008956"),reduction = "tsne") + NoLegend() +
      theme(panel.border = element_rect(size = 1,colour = "grey"),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank())
    ggsave(paste0(pwd,paste0("c",i),"/",j,"_",rownames(markers)[j],".pdf"),p,width = 4,height = 4)
  }
}

saveRDS(mye,file = "rds/without_PVTT_MLN/tumor/mye.DimRed.rds")

mye <- readRDS("rds/without_PVTT_MLN/tumor/mye.DimRed.rds")

oldCluster <- c(0:12)
novelCluster <- c("C0_PAGE2B_mye","C1_FGF19_mye","C2_CA9_mye","C3_BEX1_mye","C4_ECSCR_mye",
                  "C5_SEPHS2_mye","C6_ALDH3A1_mye","C7_FOXS1_mye","C8_ARHGEF35_mye",
                  "C9_CD79A_mye","C10_SLC22A1_mye","C11_CD2_mye","C12_FGL2_mye")
mye$novelClusters <- plyr::mapvalues(mye$seurat_clusters, from = oldCluster, to = novelCluster)


ClusterVisualization <- function(object,reduction,group) {
  p1 <- DimPlot(object, reduction = reduction, group.by = group)
  p2 <- DimPlot(object, reduction = reduction, group.by = "novelClusters")
  p <- p1 + p2
  return(p)
}

CD3 <- ClusterVisualization(mye,reduction = "tsne",group = "Group1")
ggsave("result/without_PVTT_MLN/tumor/PARTV/PARTV.2/CD3_novelcluster.pdf",CD3,width=14,height=6)
CD4 <- ClusterVisualization(mye,reduction = "tsne",group = "Group2")
ggsave("result/without_PVTT_MLN/tumor/PARTV/PARTV.2/CD4_novelcluster.pdf",CD4,width=14,height=6)
CD8 <- ClusterVisualization(mye,reduction = "tsne",group = "Group3")
ggsave("result/without_PVTT_MLN/tumor/PARTV/PARTV.2/CD8_novelcluster.pdf",CD8,width=14,height=6)

ImmDistinguish <- function(object, group, feature,fileName) {
  p1 <- DimPlot(object, group.by = group, reduction = "tsne")
  p2 <- DimPlot(object, group.by = "novelClusters", reduction = "tsne")
  p3 <- FeaturePlot(object, features = feature, reduction = "tsne")
  p <- p1 + p2 + p3
  ggsave(paste0("result/without_PVTT_MLN/tumor/myeloid/",fileName),p,device = "pdf",width = 20,height = 5.6)
}

##>>>>>>> CD3 <<<<<<<<<<<

# C0_PAGE2B_mye
ImmDistinguish(mye,"Group1","PAGE2B","CD3/CD3_noveltype_PAGE2B.pdf")
# C1_FGF19_mye
ImmDistinguish(mye,"Group1","FGF19","CD3/CD3_novelType_FGF19.pdf")
# C2_CA9_mye
ImmDistinguish(mye,"Group1","CA9","CD3/CD3_novelType_CA9.pdf")
# C3_BEX1_mye
ImmDistinguish(mye,"Group1","BEX1","CD3/CD3_novelType_BEX1.pdf")
# C4_ECSCR_mye
ImmDistinguish(mye,"Group1","ECSCR","CD3/CD3_novelType_ECSCR.pdf")
# C5_SEPHS2_mye
ImmDistinguish(mye,"Group1","SEPHS2","CD3/CD3_novelType_SEPHS2.pdf")
# C6_ALDH3A1_mye
ImmDistinguish(mye,"Group1","ALDH3A1","CD3/CD3_novelType_ALDH3A1.pdf")
# C7_FOXS1_mye
ImmDistinguish(mye,"Group1","FOXS1","CD3/CD3_novelType_FOXS1.pdf")
# C8_ARHGEF35_mye
ImmDistinguish(mye,"Group1","ARHGEF35","CD3/CD3_novelType_ARHGEF35.pdf")
# C9_CD79A_mye
ImmDistinguish(mye,"Group1","CD79A","CD3/CD3_novelType_CD79A.pdf")
# C10_SLC22A1_mye
ImmDistinguish(mye,"Group1","SLC22A1","CD3/CD3_novelType_SLC22A1.pdf")
# C11_CD2_mye
ImmDistinguish(mye,"Group1","CD2","CD3/CD3_novelType_CD2.pdf")
# C12_FGL2_mye
ImmDistinguish(mye,"Group1","FGL2","CD3/CD3_novelType_FGL2.pdf")

# >>>>>>>> CD3+T HIGH Group <<<<<<<<<<<

# C0_PAGE2B_mye,C2_CA9_mye,C3_BEX1_mye,C7_FOXS1_mye,C12_FGL2_mye

# >>>>>>>> CD3+T LOW Group <<<<<<<<<<<

# C1_FGF19_mye,C4_ECSCR_mye,C5_SEPHS2_mye,C6_ALDH3A1_mye,C8_ARHGEF35_mye,C9_CD79A_mye,
# C10_SLC22A1_mye,C11_CD2_mye


## >>>>>>>>> CD4 <<<<<<<<<<< 

# C0_PAGE2B_mye
ImmDistinguish(mye,"Group2","PAGE2B","CD4/CD4_noveltype_PAGE2B.pdf")
# C1_FGF19_mye
ImmDistinguish(mye,"Group2","FGF19","CD4/CD4_novelType_FGF19.pdf")
# C2_CA9_mye
ImmDistinguish(mye,"Group2","CA9","CD4/CD4_novelType_CA9.pdf")
# C3_BEX1_mye
ImmDistinguish(mye,"Group2","BEX1","CD4/CD4_novelType_BEX1.pdf")
# C4_ECSCR_mye
ImmDistinguish(mye,"Group2","ECSCR","CD4/CD4_novelType_ECSCR.pdf")
# C5_SEPHS2_mye
ImmDistinguish(mye,"Group2","SEPHS2","CD4/CD4_novelType_SEPHS2.pdf")
# C6_ALDH3A1_mye
ImmDistinguish(mye,"Group2","ALDH3A1","CD4/CD4_novelType_ALDH3A1.pdf")
# C7_FOXS1_mye
ImmDistinguish(mye,"Group2","FOXS1","CD4/CD4_novelType_FOXS1.pdf")
# C8_ARHGEF35_mye
ImmDistinguish(mye,"Group2","ARHGEF35","CD4/CD4_novelType_ARHGEF35.pdf")
# C9_CD79A_mye
ImmDistinguish(mye,"Group2","CD79A","CD4/CD4_novelType_CD79A.pdf")
# C10_SLC22A1_mye
ImmDistinguish(mye,"Group2","SLC22A1","CD4/CD4_novelType_SLC22A1.pdf")
# C11_CD2_mye
ImmDistinguish(mye,"Group2","CD2","CD4/CD4_novelType_CD2.pdf")
# C12_FGL2_mye
ImmDistinguish(mye,"Group2","FGL2","CD4/CD4_novelType_FGL2.pdf")

# >>>>> CD4+T HIGH Group <<<<<<<<<<<

# C0_PAGE2B_mye,C2_CA9_mye,C3_BEX1_mye,C7_FOXS1_mye,C12_FGL2_mye

# >>>>> CD4+T LOW Group <<<<<<<<<<<

# C1_FGF19_mye,C4_ECSCR_mye,C5_SEPHS2_mye,C6_ALDH3A1_mye,C8_ARHGEF35_mye,C9_CD79A_mye,
# C10_SLC22A1_mye,C11_CD2_mye

## >>>>>>> CD8 <<<<<<<<<<< 

# C0_PAGE2B_mye
ImmDistinguish(mye,"Group3","PAGE2B","CD8/CD8_noveltype_PAGE2B.pdf")
# C1_FGF19_mye
ImmDistinguish(mye,"Group3","FGF19","CD8/CD8_novelType_FGF19.pdf")
# C2_CA9_mye
ImmDistinguish(mye,"Group3","CA9","CD8/CD8_novelType_CA9.pdf")
# C3_BEX1_mye
ImmDistinguish(mye,"Group3","BEX1","CD8/CD8_novelType_BEX1.pdf")
# C4_ECSCR_mye
ImmDistinguish(mye,"Group3","ECSCR","CD8/CD8_novelType_ECSCR.pdf")
# C5_SEPHS2_mye
ImmDistinguish(mye,"Group3","SEPHS2","CD8/CD8_novelType_SEPHS2.pdf")
# C6_ALDH3A1_mye
ImmDistinguish(mye,"Group3","ALDH3A1","CD8/CD8_novelType_ALDH3A1.pdf")
# C7_FOXS1_mye
ImmDistinguish(mye,"Group3","FOXS1","CD8/CD8_novelType_FOXS1.pdf")
# C8_ARHGEF35_mye
ImmDistinguish(mye,"Group3","ARHGEF35","CD8/CD8_novelType_ARHGEF35.pdf")
# C9_CD79A_mye
ImmDistinguish(mye,"Group3","CD79A","CD8/CD8_novelType_CD79A.pdf")
# C10_SLC22A1_mye
ImmDistinguish(mye,"Group3","SLC22A1","CD8/CD8_novelType_SLC22A1.pdf")
# C11_CD2_mye
ImmDistinguish(mye,"Group3","CD2","CD8/CD8_novelType_CD2.pdf")
# C12_FGL2_mye
ImmDistinguish(mye,"Group3","FGL2","CD8/CD8_novelType_FGL2.pdf")


# >>>>> CD8+T HIGH Group <<<<<<<<<<<

# C0_PAGE2B_mye,C2_CA9_mye,C7_FOXS1_mye,C8_ARHGEF35_mye,C9_CD79A_mye,C11_CD2_mye

# >>>>> CD8+T LOW Group <<<<<<<<<<<

# C1_FGF19_mye,C3_BEX1_mye,C4_ECSCR_mye,C5_SEPHS2_mye,C6_ALDH3A1_mye,C10_SLC22A1_mye,C12_FGL2_mye


#>>>>>> Correlation Analysis <<<<<<<<<<

# GMDSC: CD45+, CD3-, B220-, NK1.1-, gdTCR-, CD11b+, CD14-, CD15+, CD33+, CD66b+
# M-MDSC:CD45+,CD3-, B220-, NK1.1-, gdTCR-, CD11b+, CD14+, CD15-, CD33+, CD66b-
MDSCFeatures <- c("PTPRC","CD3G","KLRB1","ITGAM","CD14","FUT4","CD33","CEACAM8")
sigs <- c("PAGE2B","FGF19","CA9","BEX1","ECSCR","SEPHS2","ALDH3A1","FOXS1","ARHGEF35","CD79A","SLC22A1","CD2","FGL2")
genes <- c(MDSCFeatures,sigs)

expr <- GetAssayData(mye,assay = "RNA",slot = "data")[genes,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(mye), to = mye$Sample)
colnames(expr.T) <- c(genes,"Sample")

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
colnames(patient) <- genes
xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/myeloid/corr/GMDSC/avgExpPerPatient.xlsx")

pwd = "result/without_PVTT_MLN/tumor/myeloid/corr/GMDSC/singlecell/"

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
rm(list = ls());gc()
load("refdb/TCGA_LIHC.rnaseq.RData")

genesTBL <- xlsx::read.xlsx("result/without_PVTT_MLN/tumor/myeloid/corr/GMDSC/avgExpPerPatient.xlsx",sheetIndex = 1)
genes <- colnames(genesTBL)[-1]
setdiff(genes, colnames(LIHC.rnaseq))

expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)
pwd = "result/without_PVTT_MLN/tumor/myeloid/corr/GMDSC/TCGA/"

for(gene in colnames(expr.norm)[1:8]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

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


# TAM:CD45+CD3e-CD19-CD56-CD16-CD206+CD163+
TAMFeatures <- c("PTPRC","CD3E","CD19","NCAM1","FCGR3A","MRC1","CD163")
sigs <- c("PAGE2B","FGF19","CA9","BEX1","ECSCR","SEPHS2","ALDH3A1","FOXS1","ARHGEF35","CD79A","SLC22A1","CD2","FGL2")
genes <- c(TAMFeatures,sigs)

expr <- GetAssayData(mye,assay = "RNA",slot = "data")[genes,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(mye), to = mye$Sample)
colnames(expr.T) <- c(genes,"Sample")

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
colnames(patient) <- genes
xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/myeloid/corr/TAM/avgExpPerPatient.xlsx")

pwd = "result/without_PVTT_MLN/tumor/myeloid/corr/TAM/singlecell/"

for(gene in colnames(patient)[1:7]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

for(i in 1:7){
  for(j in 8:ncol(patient)){
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
rm(list = ls());gc()
load("refdb/TCGA_LIHC.rnaseq.RData")

genesTBL <- xlsx::read.xlsx("result/without_PVTT_MLN/tumor/myeloid/corr/TAM/avgExpPerPatient.xlsx",sheetIndex = 1)
genes <- colnames(genesTBL)[-1]
setdiff(genes, colnames(LIHC.rnaseq))
expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)
pwd = "result/without_PVTT_MLN/tumor/myeloid/corr/TAM/TCGA/"

for(gene in colnames(expr.norm)[1:7]){
  if(!file.exists(paste0(pwd,gene))){
    dir.create(paste0(pwd,gene))
  }
}

for(i in 1:7){
  for(j in 8:ncol(expr.norm)){
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
