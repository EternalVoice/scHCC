
library(Seurat)
library(ggplot2)

## >>>>> Requirement I: TAM

# 分析TAM中高表达的基因如LAIR1，并分析它们的表达在我们的数据库以及在TCGA 肝癌数据库中与CD163，CD206, CD8，IFNg的相关性

tumor <- readRDS("rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR_2.rds")

# pal1 <- c("#BF0032","#0066A5","#008956","#E35822","#875691","#892C16","#B3446C")

cellEmbeddings <- function(object, cellCluster){
  cells <- rownames(subset(object@meta.data, novelCluster == cellCluster))
  cells.emb <- object@reductions$umap@cell.embeddings
  cells.emb <- as.data.frame(cells.emb, stringAsFactors = F)
  
  target = NULL
  for(cell in cells){
    pos <- cells.emb[rownames(cells.emb) == cell, ]
    target <- rbind(target, pos)
  }
  return(target)
}

target <- cellEmbeddings(object = tumor, cellCluster = "TAM")

target <- target[target$UMAP_1 > 0 & target$UMAP_2 > -9 & target$UMAP_2 < 3,]

TAM <- DimPlot(tumor,group.by = "CellType",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
               cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())

ggsave("result/without_PVTT_MLN/tumor/TAM/TAM.highlight.pdf",TAM,device = "pdf",width = 6,height = 6)

Idents(tumor) <- tumor$CellType

markers <- FindMarkers(tumor,only.pos = TRUE,ident.1 = "TAM")
genes <- rownames(markers)
for(i in seq_along(genes)){
  p <- FeaturePlot(tumor,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/TAM/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}

## 分析TAM中高表达的基因如LAIR1，并分析它们的表达在数据集以及在TCGA中与
## CD163、CD206、CD8、IFNg的相关性。

## VALIDATION

# CD163 MRC1 CD8A CD8B IFNG

gene <- c("CD163","MRC1","CD8A","CD8B","IFNG",rownames(markers))

# singlecell
expr <- GetAssayData(tumor,assay = "RNA",slot = "data")[gene,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(tumor), to = tumor$Sample)
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
xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/TAM/avgExpPerPatient.xlsx")

colnames(patient) <- gsub("-",".",colnames(patient))
patient <- patient[-c(75,138)]
pwd = "result/without_PVTT_MLN/tumor/TAM/"
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
load("refdb/TCGA_36cancers.rnaseq.RData")
LIHC.rnaseq <- TCGA_36cancers.rnaseq$LIHC.rnaseq
save(LIHC.rnaseq,file = "refdb/TCGA_LIHC.rnaseq.RData")

rm(TCGA_36cancers.rnaseq);gc()

genesTBL <- xlsx::read.xlsx("result/without_PVTT_MLN/tumor/TAM/avgExpPerPatient.xlsx",sheetIndex = 1)

genes <- gsub("\\.","-",colnames(genesTBL)[-1]);genes <- genes[-c(75,138,142)]

setdiff(genes, colnames(LIHC.rnaseq))

# CXCL8 ==> IL8; CTSL ==> CTSL1; TMSB4X ==> NULL; LINC01272 ==> NULL
# THEMIS2 ==> C1orf38; MOB1A ==> MOBKL1B; JAML ==> AMICA1; TMIGD3 ==> NULL
# SCIMP ==> C17orf87; ADGRE2 ==> EMR2; MIR3945HG ==> NULL; NABP1 ==> OBFC2A
# LINC00936 ==> NULL; HACD4 ==> PTPLAD2; YBX3 ==> CSDA; PKM ==> PKM2
# SRGAP2B ==> NULL; PPP1R18 ==> KIAA1949; MTPN ==> NULL; ERO1A ==> ERO1L

genes[grep("CXCL8",genes)] <- "IL8";genes[grep("CTSL",genes)] <- "CTSL1"
genes <- genes[-grep("TMSB4X",genes)];genes <- genes[-grep("LINC01272",genes)]
genes[grep("THEMIS2",genes)] <- "C1orf38";genes[grep("MOB1A",genes)] <- "MOBKL1B"
genes[grep("JAML",genes)] <- "AMICA1";genes <- genes[-grep("TMIGD3",genes)]
genes[grep("SCIMP",genes)] <- "C17orf87";genes[grep("ADGRE2",genes)] <- "EMR2"
genes <- genes[-grep("MIR3945HG",genes)];genes[grep("NABP1",genes)] <- "OBFC2A"
genes <- genes[-grep("LINC00936",genes)];genes[grep("HACD4",genes)] <- "PTPLAD2"
genes[grep("YBX3",genes)] <- "CSDA";genes[grep("PKM",genes)] <- "PKM2"
genes <- genes[-grep("SRGAP2B",genes)];genes[grep("PPP1R18",genes)] <- "KIAA1949"
genes <- genes[-grep("MTPN",genes)];genes[grep("ERO1A",genes)] <- "ERO1L"

expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)
pwd = "result/without_PVTT_MLN/tumor/TAM/corr/TCGA/"

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



## >>>>> Requirement I: MDSC

rm(list = ls());gc()

load("refdb/TCGA_LIHC.rnaseq.RData")
tumor <- readRDS("RDS/without_PVTT_MLN/tumor/3.tumor_combine_T.IR_2.rds")

# CD45+ CD3- B220- NK1.1- gdTCR- CD11b+ CD14- CD15+ CD33+ CD66b+
# PTPRC+, CD3G-, KLRB1-, ITGAM+, CD14-, FUT4+, CD33+, CEACAM8+
markers <- c("PTPRC","CD3G","KLRB1","ITGAM","CD14","FUT4","CD33","CEACAM8")
p1 <- VlnPlot(tumor,features = markers[1],pt.size = 0,group.by = "seurat_clusters") + NoLegend() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0,size = 14),
        axis.title.y = element_blank())
p2 <- VlnPlot(tumor,features = markers[2],pt.size = 0,group.by = "seurat_clusters") + NoLegend() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
p3 <- VlnPlot(tumor,features = markers[3],pt.size = 0,group.by = "seurat_clusters") + NoLegend() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
p4 <- VlnPlot(tumor,features = markers[4],group.by = "seurat_clusters",slot = "data") + NoLegend() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
p5 <- VlnPlot(tumor,features = markers[5],pt.size = 0,group.by = "seurat_clusters") + NoLegend() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
p6 <- VlnPlot(tumor,features = markers[6],group.by = "seurat_clusters") + NoLegend() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
p7 <- VlnPlot(tumor,features = markers[7],group.by = "seurat_clusters") + NoLegend() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
p8 <- VlnPlot(tumor,features = markers[8],group.by = "seurat_clusters") + NoLegend() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
pc <- p8/p7/p6/p5/p4/p3/p2/p1

ggsave("result/without_PVTT_MLN/tumor/MDSC/markers.pdf",pc,width = 8,height = 16)

p1 <- DimPlot(tumor,group.by = "seurat_clusters",label = T,label.size = 5)+NoLegend()+theme(plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/MDSC/seurat_clusters.pdf",p1,width = 6,height = 6)

# seurat_cluster 17: G-MDSC

GMDSC.cells <- rownames(subset(tumor@meta.data,seurat_clusters == "17"))
p2 <- DimPlot(tumor,cells.highlight = GMDSC.cells) + NoLegend() + ggtitle("G-MDSC") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black",size = 1))
ggsave("result/without_PVTT_MLN/tumor/MDSC/G-MDSC_highlight.pdf",p2,width = 6,height = 6)

# seurat_cluster 14,16,20: M-MDSC
MMDSC.cells <- rownames(subset(tumor@meta.data, seurat_clusters=="14" | seurat_clusters=="16" | seurat_clusters=="20"))
p3 <- DimPlot(tumor,cells.highlight = MMDSC.cells) + NoLegend() + ggtitle("M-MDSC") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black",size = 1))
ggsave("result/without_PVTT_MLN/tumor/MDSC/M-MDSC_highlight.pdf",p3,width = 6,height = 6)

## 分析G-MDSC中高表达的基因，并分析它们的表达在数据集以及在TCGA中与CD33,CD66b,CD11B,CD8,IFNg的相关性。
Idents(tumor) <- tumor$seurat_clusters
markers <- FindMarkers(tumor,only.pos = TRUE,ident.1 = "17")
genes <- rownames(markers)

for(i in seq_along(genes)){
  p <- FeaturePlot(tumor,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/MDSC/G-MDSC/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}

# CD163 MRC1 CD8A CD8B IFNG

gene <- c("CD33","CEACAM8","ITGAM","CD8A","CD8B","IFNG",rownames(markers))
expr <- GetAssayData(tumor,assay = "RNA",slot = "data")[gene,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(tumor), to = tumor$Sample)
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
xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/MDSC/G-MDSC.avgExpPerPatient.xlsx")


## validation

# singecell
pwd = "result/without_PVTT_MLN/tumor/MDSC/G-MDSC/"
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
rm(list = ls());gc()
load("refdb/TCGA_LIHC.rnaseq.RData")
genesTBL <- xlsx::read.xlsx("result/without_PVTT_MLN/tumor/MDSC/G-MDSC.avgExpPerPatient.xlsx",sheetIndex = 1)
genes <- gsub("\\.","-",colnames(genesTBL)[-1])
setdiff(genes, colnames(LIHC.rnaseq))

genes <- genes[-grep("RP5-887A10-1",genes)];genes <- genes[-grep("AL928768-3",genes)]
genes <- genes[-grep("IGHD",genes)];genes <- genes[-grep("LINC00926",genes)]
genes <- genes[-grep("RP11-693J15-5",genes)];genes <- genes[-grep("AC079767-4",genes)]
genes <- genes[-grep("RP11-231C14-7",genes)];genes[grep("FCMR",genes)] <- "FAIM3"
genes[grep("KIAA0226L",genes)] <- "C13orf18";genes[grep("SMIM14",genes)] <- "C4orf34"
genes[grep("TMEM243",genes)] <- "C7orf23";genes[grep("ZFAS1",genes)] <- "C20orf199"
genes <- genes[-grep("IGHM",genes)];genes[grep("KMT2E",genes)] <- "MLL5"
genes[grep("SCIMP",genes)] <- "C17orf87";genes[grep("SRSF5",genes)] <- "SFRS5"
genes <- genes[-grep("RP11-347P5-1",genes)];genes <- genes[-grep("AC090498-1",genes)]
genes <- genes[-grep("RP5-1171I10-5",genes)];genes <- genes[-grep("CTD-3252C9-4",genes)]
genes[grep("SRSF7",genes)] <- "SFRS7";genes <- genes[-grep("C1orf132",genes)]
genes[grep("PNISR",genes)] <- "SFRS18";genes[grep("HNRNPDL",genes)] <- "HNRPDL"
genes <- genes[-grep("RP11-138A9-2",genes)];genes[grep("SMDT1",genes)] <- "C22orf32"
genes[grep("RSRP1",genes)] <- "C1orf63";genes <- genes[-grep("AC016831-7",genes)]
genes <- genes[-grep("RP11-138A9-1",genes)];genes[grep("ADGRE5",genes)] <- "CD97"
genes[grep("SARAF",genes)] <- "TMEM66";genes[grep("LDLRAD4",genes)] <- "C18orf1"
genes[grep("THEMIS2",genes)] <- "C1orf38";genes[grep("CNTRL",genes)] <- "CEP110"
genes[grep("MOB3A",genes)] <- "MOBKL2A";genes <- genes[-grep("RP1-313I6-12",genes)]
genes <- genes[-grep("LINC-PINT",genes)];genes[grep("PCNXL4",genes)] <- "C14orf135"
genes[grep("ELMSAN1",genes)] <- "C14orf43";genes[grep("MOB1A",genes)] <- "MOBKL1B"
genes[grep("MIS18BP1",genes)] <- "C14orf106";genes[grep("PRMT9",genes)] <- "PRMT10"
genes[grep("CYB561A3",genes)] <- "CYBASC3";genes[grep("KIAA1551",genes)] <- "C12orf35"
genes[grep("BOD1L1",genes)] <- "BOD1L";genes[grep("METTL21A",genes)] <- "FAM119A"

expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)

pwd = "result/without_PVTT_MLN/tumor/MDSC/G-MDSC/corr/TCGA/"
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


## 分析M-MDSC中高表达的基因，并分析它们的表达在数据集以及在TCGA中与CD33,CD11B,CD8的相关性。

Idents(tumor) <- tumor$seurat_clusters
markers <- FindMarkers(tumor,only.pos = TRUE,ident.1 = c(14,16,20))
genes <- rownames(markers)

for(i in seq_along(genes)){
  p <- FeaturePlot(tumor,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/MDSC/M-MDSC/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}

# CD33 CD11B CD8
gene <- c("CD33","ITGAM","CD8A","CD8B",rownames(markers))
expr <- GetAssayData(tumor,assay = "RNA",slot = "data")[gene,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(tumor), to = tumor$Sample)
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
xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/MDSC/M-MDSC.avgExpPerPatient.xlsx")

## validation

# singecell
pwd = "result/without_PVTT_MLN/tumor/MDSC/M-MDSC/corr/singlecell/"
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
rm(list = ls());gc()
load("refdb/TCGA_LIHC.rnaseq.RData")
genesTBL <- xlsx::read.xlsx("result/without_PVTT_MLN/tumor/MDSC/M-MDSC.avgExpPerPatient.xlsx",sheetIndex = 1)
genes <- gsub("\\.","-",colnames(genesTBL)[-1])
genes <- genes[-304]
setdiff(genes, colnames(LIHC.rnaseq))

genes[grep("CXCL8",genes)] <- "IL8";genes[grep("CTSL",genes)] <- "CTSL1"
genes[grep("THEMIS2",genes)] <- "C1orf38";genes <- genes[-grep("TMSB4X",genes)]
genes <- genes[-grep("TMIGD3",genes)];genes[grep("SCIMP",genes)] <- "C17orf87"
genes[grep("IGFLR1",genes)] <- "U2AF1L4";genes[grep("MOB1A",genes)] <- "MOBKL1B"
genes[grep("DNAAF1",genes)] <- "LRRC50";genes <- genes[-grep("MIR3945HG",genes)]
genes <- genes[-grep("LINC00998",genes)];genes <- genes[-grep("LINC01272",genes)]
genes <- genes[-grep("APOC4-APOC2",genes)];genes[grep("ABRACL",genes)] <- "C6orf115"
genes[grep("VIMP",genes)] <- "SELS";genes[grep("SLC18B1",genes)] <- "C6orf192"

expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)

pwd = "result/without_PVTT_MLN/tumor/MDSC/M-MDSC/corr/TCGA/"
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


## >>>>> Requirement I: Tex CD8+T

# 分析Tex CD8+T中高表达的基因，并分析它们的表达在我们的数据集以及在TCGA数据集中与
# CD8,T-bet,IFNg,GranzymeB,IL-2,PD-1,LAG3,TIM3,TOX的相关性

rm(list = ls());gc()
tumor <- readRDS("rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR_2.rds")

CD8T <- readRDS("rds/without_PVTT_MLN/tumor/sub_T/CD8T.Tsub.assigned.rds")
## map new CD8T subclusters
curr_ids <- c(0:9)
new_ids <- c("Unknown","Unknown","CD8_Tex","CD8_CTLs","Unknown","Unknown","CD8_Tcm","CD8_Teff","CD8_Ttm","Unknown")
CD8T$subcluster <- plyr::mapvalues(CD8T$seurat_clusters,from = curr_ids,to = new_ids)
CD8.meta <- data.frame(row.names = colnames(CD8T), subcluster = as.character(CD8T$subcluster))
CD8.meta$subcluster <- as.character(CD8.meta$subcluster)

CD4T <- readRDS("rds/without_PVTT_MLN/tumor/sub_T/CD4T.Tsub.assigned.rds")
## map new CD4T subclusters
curr_ids <- c(0:14)
new_ids <- c("Unknown","Unknown","Unknown","CD4_Tem","CD4_Tem","CD4_Tem","CD4_Tem","CD4_Tscm","CD4_Treg","CD4_Tcm",
             "CD4_Tex","CD4_Temra","CD4_NKT","CD4_Teff","Unknown")
CD4T$subcluster <- plyr::mapvalues(CD4T$seurat_clusters,from = curr_ids, to = new_ids)
CD4.meta <- data.frame(row.names = colnames(CD4T), subcluster = as.character(CD4T$subcluster))
CD4.meta$subcluster <- as.character(CD4.meta$subcluster)

CDT.meta <- rbind(CD8.meta,CD4.meta)
CDT.meta$barcode <- rownames(CDT.meta)

meta <- data.frame(row.names = colnames(tumor), CellType2 = as.character(tumor$CellType2))
meta$CellType2 <- as.character(meta$CellType2)
meta$barcode <- rownames(meta)
tmp <- meta[-grep("^CD",meta$CellType2),]
colnames(tmp)[1] <- "subcluster"

meta.T <- rbind(CDT.meta,tmp)

tumor$CellType3 <- plyr::mapvalues(colnames(tmp),from = rownames(meta.T), to = meta.T$subcluster)

cellEmbeddings <- function(object, cellCluster){
  cells <- rownames(subset(object@meta.data, CellType3 == cellCluster))
  cells.emb <- object@reductions$umap@cell.embeddings
  cells.emb <- as.data.frame(cells.emb, stringAsFactors = F)
  
  target = NULL
  for(cell in cells){
    pos <- cells.emb[rownames(cells.emb) == cell, ]
    target <- rbind(target, pos)
  }
  return(target)
}

target <- cellEmbeddings(object = tumor, cellCluster = "CD8_Tex")

CD8_Tex <- DimPlot(tumor,group.by = "CellType3",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
               cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())

ggsave("result/without_PVTT_MLN/tumor/CD8_Tex/CD8_Tex.highlight.pdf",CD8_Tex,width = 8,height = 8)

saveRDS(tumor,file = "rds/without_PVTT_MLN/tumor/tumor.celltypes.rds")

Idents(tumor) <- tumor$CellType3
markers <- FindMarkers(tumor, only.pos = TRUE, ident.1 = "CD8_Tex")
genes <- rownames(markers)

for(i in seq_along(genes)){
  p <- FeaturePlot(tumor,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD8_Tex/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}


# CD8,T-bet,IFNg,GranzymeB,IL-2,PD-1,LAG3,TIM3,TOX
gene <- c("CD8A","CD8B","TBX21","IFNG","GZMB","IL2","PDCD1","LAG3","HAVCR2","TOX",rownames(markers))
gene <- gene[-c(13,14,83,27,15,485,29)]
expr <- GetAssayData(tumor,assay = "RNA",slot = "data")[gene,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(tumor), to = tumor$Sample)
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
xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD8_Tex/CD8_Tex.avgExpPerPatient.xlsx")


## validation

# singecell
pwd = "result/without_PVTT_MLN/tumor/CD8_Tex/corr/singlecell/"
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
rm(list = ls());gc()
load("refdb/TCGA_LIHC.rnaseq.RData")
genesTBL <- xlsx::read.xlsx("result/without_PVTT_MLN/tumor/CD8_Tex/CD8_Tex.avgExpPerPatient.xlsx",sheetIndex = 1)
genes <- gsub("\\.","-",colnames(genesTBL)[-1])
setdiff(genes, colnames(LIHC.rnaseq))

# make gene symbols in DE genes compatible with TCGA
genes <- genes[-grep("RP11-291B21-2",genes)];genes <- genes[-grep("TRBC2",genes)]
genes <- genes[-grep("TRAC",genes)];genes <- genes[-grep("TRGC2",genes)]
genes[grep("CNTRL",genes)] <- "CEP110";genes[grep("CENPU",genes)] <- "MLF1IP"
genes <- genes[-grep("AC092580-4",genes)];genes[grep("DDX39A",genes)] <- "DDX39"
genes[grep("ORC6",genes)] <- "ORC6L";genes <- genes[-grep("LINC01480",genes)]
genes[grep("MIS18BP1",genes)] <- "C14orf106";genes[grep("LRR1",genes)] <- "PPIL5"
genes[grep("MZT2B",genes)] <- "FAM128B";genes[grep("MZT2A",genes)] <- "FAM128A"
genes[grep("PKM",genes)] <- "PKM2";genes[grep("KNL1",genes)] <- "CASC5"
genes[grep("PPP1R18",genes)] <- "KIAA1949";genes[grep("SRSF2",genes)] <- "SFRS2"
genes <- genes[-grep("AC133644-2",genes)];genes[grep("ABRACL",genes)] <- "C6orf115"
genes[grep("CCDC167",genes)] <- "C6orf129";genes[grep("SLF1",genes)] <- "ANKRD32"
genes[grep("ANAPC15",genes)] <- "C11orf51";genes[grep("ANAPC15",genes)] <- "C11orf51"
genes[grep("SRSF3",genes)] <- "SFRS3";genes[grep("RABL6",genes)] <- "C9orf86"
genes[grep("SRSF7",genes)] <- "SFRS7";genes[grep("ABHD17A",genes)] <- "FAM108A1"
genes[grep("KIAA1551",genes)] <- "C12orf35";genes[grep("HNRNPDL",genes)] <- "HNRPDL"
genes <- genes[-grep("PSMB8-AS1",genes)];genes[grep("SRSF10",genes)] <- "SFRS13A"
genes <- genes[-grep("TRBC1",genes)];genes[grep("RGCC",genes)] <- "C13orf15"
genes[grep("TUBB4B",genes)] <- "TUBB2C";genes <- genes[-grep("LINC00152",genes)]
genes[grep("CEP128",genes)] <- "C14orf145";genes[grep("PRRC2C",genes)] <- "BAT2L2"
genes <- genes[-grep("TMSB4X",genes)];genes[grep("NCBP3",genes)] <- "C17orf85"
genes[grep("ADGRE5",genes)] <- "CD97";genes[grep("SRSF4",genes)] <- "SFRS4"
genes[grep("NABP2",genes)] <- "OBFC2B";genes[grep("SPDL1",genes)] <- "CCDC99"
genes[grep("HNRNPLL",genes)] <- "HNRPLL";genes[grep("PCLAF",genes)] <- "KIAA0101"
genes <- genes[-grep("X7-Sep",genes)];genes[grep("SRSF9",genes)] <- "SFRS9"
genes[grep("SRSF9",genes)] <- "SFRS9";genes[grep("LSM8",genes)] <- "NAA38"
genes[grep("SRSF11",genes)] <- "SFRS11";genes[grep("CMC2",genes)] <- "C16orf61"
genes[grep("TMA7",genes)] <- "CCDC72";genes[grep("CTDNEP1",genes)] <- "DULLARD"
genes[grep("ANKRD36C",genes)] <- "CCDC72";genes[grep("KMT2A",genes)] <- "MLL"
genes[grep("CEP295",genes)] <- "KIAA1731";genes[grep("ELOB",genes)] <- "TCEB2"
genes[grep("SNU13",genes)] <- "NHP2L1";genes[grep("PNISR",genes)] <- "SFRS18"
genes[grep("COPS9",genes)] <- "MYEOV2";genes[grep("NRDC",genes)] <- "NRD1"
genes[grep("RACK1",genes)] <- "GNB2L1";genes[grep("CHTOP",genes)] <- "C1orf77"
genes[grep("BOD1L1",genes)] <- "BOD1L";genes[grep("EMC9",genes)] <- "FAM158A"
genes <- genes[-grep("RNR1",genes)];genes <- genes[-grep("RNR2",genes)]


expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)

pwd = "result/without_PVTT_MLN/tumor/CD8_Tex/corr/TCGA/"
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


## >>>>> Requirement I: Tex CD4+T

# 分析Tex CD4+T中高表达的基因，并分析它们的表达在我们的数据集以及在TCGA数据集中与
# CD4,T-bet,IFNg,GranzymeB,IL-2,PD-1,LAG3,TIM3,TOX的相关性
rm(list = ls());gc()

tumor <- readRDS("rds/without_PVTT_MLN/tumor/tumor.celltypes.rds")

target <- cellEmbeddings(object = tumor, cellCluster = "CD4_Tex")

CD4_Tex <- DimPlot(tumor,group.by = "CellType3",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
                   cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/CD4_Tex/CD4_Tex.highlight.pdf",CD4_Tex,width = 8,height = 8)

Idents(tumor) <- tumor$CellType3
markers <- FindMarkers(tumor, only.pos = TRUE, ident.1 = "CD4_Tex")
genes <- rownames(markers)

for(i in seq_along(genes)){
  p <- FeaturePlot(tumor,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD4_Tex/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}


# CD4,T-bet,IFNg,GranzymeB,IL-2,PD-1,LAG3,TIM3,TOX
gene <- c("CD4","TBX21","IFNG","GZMB","IL2","PDCD1","LAG3","HAVCR2","TOX",rownames(markers))
gene <- gene[-c(330,97)]
expr <- GetAssayData(tumor,assay = "RNA",slot = "data")[gene,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(tumor), to = tumor$Sample)
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
xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD4_Tex/CD4_Tex.avgExpPerPatient.xlsx")


## validation

# singecell
pwd = "result/without_PVTT_MLN/tumor/CD4_Tex/corr/singlecell/"
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
rm(list = ls());gc()
load("refdb/TCGA_LIHC.rnaseq.RData")
genesTBL <- xlsx::read.xlsx("result/without_PVTT_MLN/tumor/CD4_Tex/CD4_Tex.avgExpPerPatient.xlsx",sheetIndex = 1)
genes <- gsub("\\.","-",colnames(genesTBL)[-1])
setdiff(genes, colnames(LIHC.rnaseq))

genes <- genes[-grep("TRAC",genes)];genes <- genes[-grep("LINC01281",genes)]
genes <- genes[-grep("SNHG25",genes)];genes <- genes[-grep("TRBC2",genes)]
genes <- genes[-grep("TRBC1",genes)];genes <- genes[-grep("AC133644-2",genes)]
genes <- genes[-grep("RP5-1028K7-2",genes)];genes <- genes[-grep("AC104820-2",genes)]
genes <- genes[-grep("AC017002-1",genes)];genes <- genes[-grep("MT-ND1",genes)]
genes[grep("FCMR",genes)] <- "FAIM3";genes <- genes[-grep("TMSB4X",genes)]
genes <- genes[-grep("MT-ND2",genes)];genes <- genes[-grep("AC092580-4",genes)]
genes[grep("TMA7",genes)] <- "CCDC72";genes[grep("PCED1B",genes)] <- "FAM113B"
genes[grep("MZT2A",genes)] <- "FAM128A";genes <- genes[-grep("RP11-796E2-4",genes)]
genes[grep("PKM",genes)] <- "PKM2";genes[grep("PRMT9",genes)] <- "PRMT10"
genes[grep("MZT2B",genes)] <- "FAM128B";genes <- genes[-grep("CTD-3252C9-4",genes)]
genes[grep("PNISR",genes)] <- "SFRS18";genes[grep("PPP1R18",genes)] <- "KIAA1949"
genes[grep("ZFAS1",genes)] <- "C20orf199";genes <- genes[-grep("DLGAP1-AS1",genes)]
genes <- genes[-grep("RP11-51J9-5",genes)];genes <- genes[-grep("LINC01420",genes)]
genes[grep("ABRACL",genes)] <- "C6orf115";genes <- genes[-grep("LINC00152",genes)]
genes[grep("DNPH1",genes)] <- "C6orf108";genes <- genes[-grep("MT-ND6",genes)]
genes[grep("TMEM261",genes)] <- "C9orf123";genes[grep("CCDC167",genes)] <- "C6orf129"
genes <- genes[-grep("GABPB1-AS1",genes)];genes <- genes[-grep("SMIM10L1",genes)]
genes[grep("HNRNPLL",genes)] <- "HNRPLL"

expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)

pwd = "result/without_PVTT_MLN/tumor/CD4_Tex/corr/TCGA/"
for(gene in colnames(expr.norm)[1:10]){
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


## >>>>> CTL

# 分析Treg中高表达的基因,并分析它们的表达在我们的数据库以及在TCGA数据库中与 
#  FoxP3,Helios,CD8,CD11B,CD33基因相关性

target <- cellEmbeddings(object = tumor, cellCluster = "CD4_Treg")

CD4_Treg <- DimPlot(tumor,group.by = "CellType3",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
                   cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/CD4_Treg/CD4_Treg.highlight.pdf",CD4_Treg,width = 8,height = 8)

Idents(tumor) <- tumor$CellType3
markers <- FindMarkers(tumor, only.pos = TRUE, ident.1 = "CD4_Treg")
genes <- rownames(markers)

for(i in seq_along(genes)){
  p <- FeaturePlot(tumor,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD4_Treg/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}


# FoxP3,Helios,CD8,CD11B,CD33
gene <- c("FOXP3","IKZF2","CD8A","CD8B","ITGAM","CD33",rownames(markers))
gene <- gene[-c(22,30)]
expr <- GetAssayData(tumor,assay = "RNA",slot = "data")[gene,]
expr <- as.data.frame(expr,stringsAsFactors = F)
expr.T <- data.frame(t(expr))
expr.T$Patient <- plyr::mapvalues(rownames(expr.T),from = colnames(tumor), to = tumor$Sample)
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
xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD4_Treg/CD4_Treg.avgExpPerPatient.xlsx")


## validation

# singecell
pwd = "result/without_PVTT_MLN/tumor/CD4_Treg/corr/singlecell/"
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
rm(list = ls());gc()
load("refdb/TCGA_LIHC.rnaseq.RData")
genesTBL <- xlsx::read.xlsx("result/without_PVTT_MLN/tumor/CD4_Treg/CD4_Treg.avgExpPerPatient.xlsx",sheetIndex = 1)
genes <- gsub("\\.","-",colnames(genesTBL)[-1])
setdiff(genes, colnames(LIHC.rnaseq))

genes <- genes[-grep("TRAC",genes)];genes <- genes[-grep("AC133644-2",genes)]
genes <- genes[-grep("AC002331-1",genes)];genes <- genes[-grep("TRBC2",genes)]
genes <- genes[-grep("AC017002-1",genes)];genes <- genes[-grep("AC145110-1",genes)]
genes <- genes[-grep("RP11-1399P15-1",genes)];genes <- genes[-grep("TRBC1",genes)]
genes[grep("FCMR",genes)] <- "FAIM3";genes <- genes[-grep("RP11-347P5-1",genes)]
genes[grep("SARAF",genes)] <- "TMEM66";genes <- genes[-grep("RP11-214O1-3",genes)]
genes[grep("PCED1B",genes)] <- "FAM113B";genes[grep("PCED1B",genes)] <- "FAM113B"
genes[grep("HACD1",genes)] <- "PTPLA";genes[grep("RGCC",genes)] <- "C13orf15"
genes <- genes[-grep("TMSB4X",genes)];genes[grep("SRSF7",genes)] <- "SFRS7"
genes <- genes[-grep("LINC01588",genes)];genes <- genes[-grep("MT-ND3",genes)]
genes[grep("SRSF5",genes)] <- "SFRS5";genes[grep("DNPH1",genes)] <- "C6orf108"
genes <- genes[-grep("AC074289-1",genes)];genes <- genes[-grep("RP11-796E2-4",genes)]
genes <- genes[-grep("LINC00152",genes)];genes[grep("RSRP1",genes)] <- "C1orf63"
genes[grep("PKM",genes)] <- "PKM2";genes <- genes[-grep("MIR4435-2HG",genes)]
genes <- genes[-grep("AC016831-7",genes)];genes <- genes[-grep("LINC-PINT",genes)]
genes[grep("PPP1R18",genes)] <- "KIAA1949";genes[grep("KMT2E",genes)] <- "MLL5"
genes[grep("NSMCE3",genes)] <- "NDNL2";genes[grep("NABP1",genes)] <- "OBFC2A"
genes <- genes[-grep("RP11-138A9-2",genes)];genes[grep("CENPC",genes)] <- "CENPC1"
genes[grep("CCDC186",genes)] <- "C10orf118";genes[grep("ABRACL",genes)] <- "C6orf115"
genes[grep("PRRC2C",genes)] <- "BAT2L2";genes[grep("IGFLR1",genes)] <- "U2AF1L4"
genes <- genes[-grep("PET100",genes)];genes[grep("ABHD17A",genes)] <- "FAM108A1"
genes[grep("KIAA1551",genes)] <- "C12orf35";genes[grep("CCSER2",genes)] <- "FAM190B"
genes <- genes[-grep("GABPB1-AS1",genes)];genes[grep("LSM8",genes)] <- "NAA38"
genes[grep("HNRNPLL",genes)] <- "HNRPLL";genes <- genes[-grep("DLGAP1-AS1",genes)]
genes <- genes[-grep("RP11-138A9-1",genes)];genes[grep("BOD1L1",genes)] <- "BOD1L"
genes[grep("PNISR",genes)] <- "SFRS18";genes[grep("SAYSD1",genes)] <- "C6orf64"
genes[grep("CFAP20",genes)] <- "C16orf80";genes[grep("SCAF11",genes)] <- "SFRS2IP"
genes[grep("CREBRF",genes)] <- "C5orf41";genes[grep("SRPRA",genes)] <- "SRPR"
genes[grep("NT5C3A",genes)] <- "NT5C3";genes[grep("SC5D",genes)] <- "SC5DL"

expr.LIHC <- LIHC.rnaseq[,genes]
expr.norm <- log2(expr.LIHC + 1)

pwd = "result/without_PVTT_MLN/tumor/CD4_Treg/corr/TCGA/"
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
