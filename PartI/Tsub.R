
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

if(!file.exists("rds/without_PVTT_MLN/tumor/sub_T/Tcell.rds")){
  ## obtain total T cell population
  tumor <- readRDS("rds/without_PVTT_MLN/tumor/4.tumor_combine_T.IR_MDSC.IR.rds")
  Tcell <- subset(tumor, cells = rownames(subset(tumor@meta.data, CellType == "T cell")))
} else {
  Tcell <- readRDS("rds/without_PVTT_MLN/tumor/sub_T/Tcell.rds")
}

DefaultAssay(Tcell) <- "RNA"
# 降维聚类
Tcell <- FindVariableFeatures(Tcell,selection.method = "vst",nfeatures = 2000)
DefaultAssay(Tcell) <- "integrated"
Tcell <- ScaleData(Tcell,features = rownames(Tcell))
Tcell <- RunPCA(Tcell, features = VariableFeatures(Tcell))
p1 <- DimPlot(Tcell, reduction = "pca", group.by = "orig.ident",cols=rainbow(length(levels(as.factor(Tcell$orig.ident)))))
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/1.Tsub_pca.pdf",p1,device="pdf",width=8,height=6)
p2 <- ElbowPlot(Tcell,ndims=30,reduction="pca")
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/2.Tsub_elbow.pdf",p2,device="pdf",width=6,height=4)

pc.num = 1:24

Tcell <- FindNeighbors(Tcell, dims = pc.num) %>% FindClusters(resolution = 0.5)
Tcell <- RunUMAP(Tcell,reduction="pca",dims=pc.num) %>% RunTSNE(dims=pc.num)

pal <- c("#786140","#8C0603","#B2712B","#B48521","#C6EEC4","#DCD3B1","#A7BBA5","#0202BB",
         "#5EC6B9","#4E8A97","#726DBA","#090795","#AD66B6","#B71CB4","#A01153","#008000")
pal2 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
  '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#F3C300','#222222')
p1 <- DimPlot(Tcell,reduction="tsne",cols=pal,pt.size = 0.9)
p2 <- DimPlot(Tcell,reduction="tsne",group.by="orig.ident",cols=pal2,pt.size = 0.9) + 
  theme(plot.title = element_blank())
p3 <- DimPlot(Tcell, reduction = "tsne", group.by = "CellType2",pt.size = 0.9, cols = c("#BF0032","#0066A5")) + 
  theme(plot.title = element_blank())
pc <- p1 + p2 +p3
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/3.Tcell_DimRed.pdf",pc, width = 18,height = 4.5)

# CD11B- CD11C-
DefaultAssay(Tcell) <- "RNA"
p <- VlnPlot(Tcell,features = c("ITGAM","ITGAX"),pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/4.CD11B_CD11C.pdf",p,width = 6,height = 4)

saveRDS(Tcell,file = "rds/without_PVTT_MLN/tumor/sub_T/Tcell.DimRed.rds")

rm(list = ls());gc()

###>>>>>>>>>>>>>>>>> Cluster-specific Genes   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

table(Tcell$seurat_clusters,Tcell$CellType2)

# cell control

CellFilter <- function(seurat.obj,seurat.cluster,CellType){
  cluster.cell <- colnames(seurat.obj)[seurat.obj$seurat_clusters == seurat.cluster & seurat.obj$CellType2 == CellType]
  clust.cell <- colnames(seurat.obj)
  for(i in cluster.cell) { clust.cell = clust.cell[-which(clust.cell == i)]}
  seurat.obj <- subset(seurat.obj, cells = clust.cell)
  return(seurat.obj)
}

# cluster 0
Tcell_new <- CellFilter(Tcell,0,"CD8+ T cell")
# cluster 1
Tcell_new <- CellFilter(Tcell_new,1,"CD8+ T cell")
# cluster 2
Tcell_new <- CellFilter(Tcell_new,2,"CD8+ T cell")
# cluster 3
Tcell_new <- CellFilter(Tcell_new,3,"CD4+ T cell")
# cluster 5
Tcell_new <- CellFilter(Tcell_new,5,"CD4+ T cell")
# cluster 7
Tcell_new <- CellFilter(Tcell_new,7,"CD8+ T cell")
# cluster 8
Tcell_new <- CellFilter(Tcell_new,8,"CD8+ T cell")
# cluster 9
Tcell_new <- CellFilter(Tcell_new,9,"CD8+ T cell")
# cluster 10
Tcell_new <- CellFilter(Tcell_new,10,"CD8+ T cell")
# cluster 11
Tcell_new <- CellFilter(Tcell_new,11,"CD8+ T cell")
# cluster 12
Tcell_new <- CellFilter(Tcell_new,12,"CD4+ T cell")
# cluster 15
Tcell_new <- CellFilter(Tcell_new,15,"CD8+ T cell")

saveRDS(Tcell_new,file = "rds/without_PVTT_MLN/tumor/sub_T/Tcell.new.DimRed.rds")

rm(list = ls());gc()
Tcell_new <- readRDS("rds/without_PVTT_MLN/tumor/sub_T/Tcell.new.DimRed.rds")

p <- DimPlot(Tcell_new)
ggsave("result/without_PVTT_MLN/tumor/PARTIII/2.sub_T/3.novel_Clusters/clusters.pdf",p,width = 4.5,height = 4)

# cluster 0: CD4_XIST
f0 <- "XIST"
p0 <- VlnPlot(Tcell_new,features = f0,group.by = "seurat_clusters",pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/PARTIII/2.sub_T/3.novel_Clusters/markers/CD4_XIST.pdf",p0,width = 6,height = 2) 
# Cluster 2: CD4_CHGA
f2 <- "CHGA"
p2 <- VlnPlot(Tcell_new,features = f2,group.by = "seurat_clusters",pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/PARTIII/2.sub_T/3.novel_Clusters/markers/CD4_CHGA.pdf",p2,width = 6,height = 2) 
# Cluster 4: CD4_RBP7
f4 <- "RBP7"
p4 <- VlnPlot(Tcell_new,features = f4,group.by = "seurat_clusters",pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/PARTIII/2.sub_T/3.novel_Clusters/markers/CD4_RBP7.pdf",p4,width = 6,height = 2) 
# Cluster 5: CD8_CXCL13
f5 <- "CXCL13"
p5 <- VlnPlot(Tcell_new,features = f5,group.by = "seurat_clusters",pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/PARTIII/2.sub_T/3.novel_Clusters/markers/XIST.pdf",p5,width = 6,height = 2) 
# Cluster 7:CD4_VEGFB
f7 <- "VEGFB"
p7 <- VlnPlot(Tcell_new,features = f7,group.by = "seurat_clusters",pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/PARTIII/2.sub_T/3.novel_Clusters/markers/CD4_VEGFB.pdf",p7,width = 6,height = 2) 
# Cluster 8: CD4_VCX3A
f8 <- "VCX3A"
p8 <- VlnPlot(Tcell_new,features = f8,group.by = "seurat_clusters",pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/PARTIII/2.sub_T/3.novel_Clusters/markers/CD4_VCX3A.pdf",p8,width = 6,height = 2) 
# Cluster 9: CD4_PLCG2
f9 <- "PLCG2"
p9 <- VlnPlot(Tcell_new,features = f9,group.by = "seurat_clusters",pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/PARTIII/2.sub_T/3.novel_Clusters/markers/CD4_PLCG2.pdf",p9,width = 6,height = 2) 
# Cluster 10: CD4_VCX3A
f10 <- "TSPAN13"
p10 <- VlnPlot(Tcell_new,features = f10,group.by = "seurat_clusters",pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/PARTIII/2.sub_T/3.novel_Clusters/markers/CD4_VCX3A.pdf",p10,width = 6,height = 2) 
# Cluster 11: CD4_IL1R2
f11 <- "IL1R2"
p11 <- VlnPlot(Tcell_new,features = f11,group.by = "seurat_clusters",pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/PARTIII/2.sub_T/3.novel_Clusters/markers/CD4_IL1R2.pdf",p11,width = 6,height = 2) 
# Cluster 12: CD8_FCGR3A
f12 <- "FCGR3A"
p12 <- VlnPlot(Tcell_new,features = f12,group.by = "seurat_clusters",pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/PARTIII/2.sub_T/3.novel_Clusters/markers/CD8_FCGR3A.pdf",p12,width = 6,height = 2) 
# Cluster 13:B3GNT7
# Cluster 13:LAT2
# cluster 13:IRF8
# Cluster 13:CD4_FGR
f13 <- "FGR"
p13 <- VlnPlot(Tcell_new,features = f13,group.by = "seurat_clusters",pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/PARTIII/2.sub_T/3.novel_Clusters/markers/CD4_FGR.pdf",p13,width = 6,height = 2) 

# cluster 14: CD4_AKR1C4
# cluster 14: ADH1B
# cluster 14: GPX2
# Cluster 14: KNG1
# Cluster 14: CHI3L1
# Cluster 14: LEAP2
f14 <- "AKR1C4"
p14 <- VlnPlot(Tcell_new,features = f14,group.by = "seurat_clusters",pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/PARTIII/2.sub_T/3.novel_Clusters/markers/CD4_AKR1C4.pdf",p14,width = 6,height = 2) 
# cluster 15: CD4_CRP
f15 <- "CRP"
p15 <- VlnPlot(Tcell_new,features = f15,group.by = "seurat_clusters",pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/PARTIII/2.sub_T/3.novel_Clusters/markers/CD4_CRP.pdf",p15,width = 6,height = 2) 

current_id <- c(0:15)
new_id <- c("CD4_XIST","Unknown","CD4_CHGA","Unknown","CD4_RBP7","CD8_CXCL13","Unknown","CD4_VEGFB","CD4_VCX3A",
            "CD4_PLCG2","CD4_VCX3A","CD4_IL1R2","CD8_FCGR3A","CD4_FGR","CD4_AKR1C4","CD4_CRP")

Tcell_new$novelCluster <- plyr::mapvalues(Tcell_new$seurat_clusters,from = current_id,to = new_id)

pal <- c("#BF0032","#828281","#008956","#892C16","#875691","#2A3C26",
         "#0066A5","#644422","#F3C300","#222222","#E68FAC","#E35822","#C3B381")
p <- DimPlot(Tcell_new,cols = pal,group.by = "novelCluster") + theme(plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/sub_T/Tcell_novel_subset.pdf",p,device = "pdf",width = 6,height = 4)


# markers <- c("XIST","CHGA","RBP7","CXCL13","VEGFB","VCX3A","PLCG2","VCX3A","IL1R2","FCGR3A","FGR","AKR1C4","CRP")
# p2 <- StackedVlnPlot(obj = Tcell, features = markers)
# ggsave("newClusters.pdf",p2,height = 13,width = 7)

saveRDS(Tcell_new, file = "rds/without_PVTT_MLN/tumor/sub_T/Tcell.new.DimRed.rds")


#############################   CD4+T CD8+T ###############################
CD4T <- subset(Tcell, cells = rownames(subset(Tcell@meta.data, CellType2 == "CD4+ T cell")))
CD8T <- subset(Tcell, cells = rownames(subset(Tcell@meta.data, CellType2 == "CD8+ T cell")))
saveRDS(CD4T,file = "rds/without_PVTT_MLN/tumor/sub_T/CD4T.rds")
saveRDS(CD8T,file = "rds/without_PVTT_MLN/tumor/sub_T/CD8T.rds")


CD4T <- readRDS("rds/without_PVTT_MLN/tumor/sub_T/CD4T.rds")

DefaultAssay(CD4T) <- "RNA"
# 降维聚类
CD4T <- FindVariableFeatures(CD4T,selection.method = "vst",nfeatures = 2000)
DefaultAssay(CD4T) <- "integrated"
CD4T <- ScaleData(CD4T,features = rownames(CD4T))
CD4T <- RunPCA(CD4T, features = VariableFeatures(CD4T))
p1 <- DimPlot(CD4T, reduction = "pca", group.by = "orig.ident",cols=rainbow(length(levels(as.factor(CD4T$orig.ident)))))
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/1.CD4T_pca.pdf",p1,device="pdf",width=8,height=6)
p2 <- ElbowPlot(CD4T,ndims=30,reduction="pca")
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/2.CD4T_elbow.pdf",p2,device="pdf",width=6,height=4)

pc.num = 1:27

CD4T <- FindNeighbors(CD4T, dims = pc.num) %>% FindClusters(resolution = 0.5)
CD4T <- RunUMAP(CD4T,reduction="pca",dims=pc.num) %>% RunTSNE(dims=pc.num)

pal <- c("#786140","#8C0603","#B2712B","#B48521","#C6EEC4","#DCD3B1","#A7BBA5","#0202BB",
         "#5EC6B9","#4E8A97","#726DBA","#090795","#AD66B6","#B71CB4","#A01153")
pal2 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#F3C300','#222222')
p1 <- DimPlot(CD4T,reduction="tsne",cols=pal,pt.size = 0.9)
p2 <- DimPlot(CD4T,reduction="tsne",group.by="orig.ident",cols=pal2,pt.size = 0.9) + 
  theme(plot.title = element_blank())
pc <- p1 + p2
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/3.CD4T_DimRed.pdf",pc, width = 12,height = 5.5)

saveRDS(CD4T, file = "rds/without_PVTT_MLN/tumor/sub_T/CD4T.DimRed.rds")


# CD11B- CD11C-
DefaultAssay(CD4T) <- "RNA"
p <- VlnPlot(CD4T,features = c("ITGAM","ITGAX"),pt.size = 0,ncol = 1) + NoLegend()
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/4.CD11B_CD11C.pdf",p,width = 6,height = 4)


CD8T <- readRDS("rds/without_PVTT_MLN/tumor/sub_T/CD8T.rds")

DefaultAssay(CD8T) <- "RNA"
# 降维聚类
CD8T <- FindVariableFeatures(CD8T,selection.method = "vst",nfeatures = 2000)
DefaultAssay(CD8T) <- "integrated"
CD8T <- ScaleData(CD8T,features = rownames(CD8T))
CD8T <- RunPCA(CD8T, features = VariableFeatures(CD8T))
p1 <- DimPlot(CD8T, reduction = "pca", group.by = "orig.ident",cols=rainbow(length(levels(as.factor(CD8T$orig.ident)))))
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/1.CD8T_pca.pdf",p1,device="pdf",width=8,height=6)
p2 <- ElbowPlot(CD8T,ndims=30,reduction="pca")
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/2.CD8T_elbow.pdf",p2,device="pdf",width=6,height=4)

pc.num = 1:27

CD8T <- FindNeighbors(CD8T, dims = pc.num) %>% FindClusters(resolution = 0.5)
CD8T <- RunUMAP(CD8T,reduction="pca",dims=pc.num) %>% RunTSNE(dims=pc.num)

pal <- c("#786140","#8C0603","#B2712B","#B48521","#C6EEC4","#DCD3B1","#A7BBA5","#0202BB","#5EC6B9","#4E8A97")
pal2 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#F3C300','#222222')
p1 <- DimPlot(CD8T,reduction="tsne",cols=pal,pt.size = 1)
p2 <- DimPlot(CD8T,reduction="tsne",group.by="orig.ident",cols=pal2,pt.size = 1) + 
  theme(plot.title = element_blank())
pc <- p1 + p2
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/3.CD8T_DimRed.pdf",pc, width = 12,height = 5.5)

saveRDS(CD8T, file = "rds/without_PVTT_MLN/tumor/sub_T/CD8T.DimRed.rds")

rm(list = ls());gc()

#*************************************************************

###        Stacked Violin Plot

#*************************************************************

library(Seurat)
library(patchwork)
library(ggplot2)

## remove the x-axis text and tick
## plot margin to adjust the white space between each plot
## pass any arguments to VlnPlot in Seurat

modify_vlnplot <- function(obj,feature,pt.size=0,plot.margin=unit(c(-0.75,0,-0.75,0),"cm"),...){
  p <- VlnPlot(obj,features = feature, pt.size = pt.size,...) +
    xlab("") + ylab(feature) + 
    theme(legend.position = "none",
          plot.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0),
          axis.text.y = element_text(size = rel(1)),
          plot.margin = plot.margin)
  return(p)
}

## extract the max value of the y axis
extract_max <- function(p){
  ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

# main function
StackedVlnPlot <- function(obj,features,pt.size=0,plot.margin=unit(c(-0.75,0,-0.75,0),"cm"),...){
  plot_list <- purrr::map(features,function(x) modify_vlnplot(obj = obj,feature = x,...))
  
  # Add back x-axis title to bottom plot.
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(),
          axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value
  ymaxs <- purrr::map_dbl(plot_list, extract_max)
  plot_list <- purrr::map2(plot_list,ymaxs,function(x,y) x+scale_y_continuous(breaks = c(y)) + expand_limits(y=y))
  p <- patchwork::wrap_plots(plot_list,ncol = 1)
  return(p)
}

CD4T <- readRDS("rds/without_PVTT_MLN/tumor/sub_T/CD4T.DimRed.rds")
CD8T <- readRDS("rds/without_PVTT_MLN/tumor/sub_T/CD8T.DimRed.rds")


# CD8 MARKERS
DefaultAssay(CD8T) <- "RNA"

pal <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC")

wkd <- "result/without_PVTT_MLN/tumor/sub_T/2.new/"
#>>> Tn, naive T cell

# CD45+, CD3+, CD8+, CD27+, CD28+, CD45RA+, CD45RO-, CD57-, CD62L+(L-Selectin),
# CD69-, CD95-(FasR), CD127+, CD197+ 

# PTPRC,CD3D/E/G, CD8A/B, CD27,CD28,B3GAT1,SELL,CD69,FAS,IL7R,CCR7

Tn_markers <- c("PTPRC","CD3D","CD3E","CD3G","CD8A","CD8B","CD27","CD28","B3GAT1","SELL","CD69","FAS","IL7R","CCR7")
p1 <- StackedVlnPlot(obj = CD8T, features = Tn_markers,cols = pal)
ggsave(paste(wkd,"CD8_Tn_markers.pdf"),p1,width = 8,height = 14)

#>>> Tscm, Stem Cell Memory T Cell

# CD45+,CD3+,CD8+,CD27+,CD28+,CD45RA+,CD45RO-,CD57-,CD62L+(L-Selectin),
# CD95+(FasR), CD127+, CD197+(CCR7)

# PTPRC,CD3D/E/G,CD8A/B,CD27,CD28,B3GAT1,SELL,FAS,IL7R,CCR7

Tscm_markers <- c("PTPRC","CD3D","CD3E","CD3G","CD8A","CD8B","CD27","CD28","B3GAT1","SELL","FAS","IL7R","CCR7")
p2 <- StackedVlnPlot(obj = CD8T, features = Tscm_markers,cols = pal)
ggsave(paste(wkd,"CD8_Tscm_markers.pdf"),p2,width = 8,height = 14)

#>>> Tcm, central memory T cell

# CD45+,CD3+,CD8+,CD27+,CD28,CD45RA-,CD45R0+,CD57-,CD62L+(L-Selectin),
# CD95+(FasR),CD127+,CD197+(CCR7)

# PTPRC,CD3D/E/G,CD8A/B,CD27,CD28,B3GAT1,SELL,FAS,IL7R,CCR7

Tcm_markers <- c("PTPRC","CD3D","CD3E","CD3G","CD8A","CD8B","CD27","CD28","B3GAT1","SELL","FAS","IL7R","CCR7")
p3 <- StackedVlnPlot(obj = CD8T, features = Tcm_markers,cols=pal)
ggsave(paste(wkd,"CD8_Tcm_markers.pdf"),p3,width = 8,height = 14)

#>>> Teff, effector T cell

# CD45+,CD8+,CD3+,CD25+(IL2RA),CD69+,CD95+(FasR),CD134+(OX40),CD137+(4-1BB),
# CD154+(CD40L),Ki-67+,KLRG1+,CD62L-,CCR7+

# PTPRC,CD8A/B,CD3D/E/G,IL2RA,CD69,FAS,TNFRSF4,TNFRSF9,CD40LG,MKI67,KLRG1,SELL,CCR7

Teff_markers <- c("PTPRC","CD3D","CD3E","CD3G","CD8A","CD8B","IL2RA","CD69","FAS","TNFRSF4","TNFRSF9","CD40LG","MKI67","KLRG1","SELL","CCR7")
p4 <- StackedVlnPlot(obj = CD8T, features = Teff_markers,cols=pal)
ggsave(paste(wkd,"CD8_Teff_markers.pdf"),p4,width = 8,height = 14)


#>>> Ttm

# CD3+,CD8+,CD27+,CD28+,CD45RA-,CD45RO+,CD57+,CD62L-,CD95+,CD127+,CD197-

# CD3D/E/G,CD8A/B,CD27,CD28,B3GAT1,SELL,FAS,IL7R,CCR7

Ttm_markers <- c("CD3D","CD3E","CD3G","CD8A","CD8B","CD27","CD28","B3GAT1","SELL","FAS","IL7R","CCR7")
p5 <- StackedVlnPlot(obj = CD8T, features = Ttm_markers,cols=pal)
ggsave(paste(wkd,"CD8_Ttm_markers.pdf"),p5,width = 8,height = 14)

#>>> Tem, effector memory cell

# CD3+,CD8+,CD27+/-,CD28-,CD45RA-,CD45RO+,CD57+,CD62L-,CD95+,CD127+/-,CD197-

# CD3D/E/G,CD8A/B,CD27,CD28,B3GAT1,SELL,FAS,IL7R,CCR7

Tem_markers <- c("CD3D","CD3E","CD3G","CD8A","CD8B","CD27","CD28","B3GAT1","SELL","FAS","IL7R","CCR7")
p6 <- StackedVlnPlot(obj = CD8T, features = Tem_markers,cols=pal)
ggsave(paste(wkd,"CD8_Tem_markers.pdf"),p6,width = 8,height = 14)

#>>> Trm, Tissue-resident memory cell

# CD45+,CD8+,CD3+,CD69+,CD103+,CTLA4+

# PTPRC,CD8A/B,CD3D/E/G,CD69,ITGAE,CTLA4

Trm_markers <- c("PTPRC","CD3D","CD3E","CD3G","CD8A","CD8B","CD69","ITGAE","CTLA4")
p7 <- StackedVlnPlot(obj = CD8T, features = Trm_markers,cols=pal)
ggsave(paste(wkd,"CD8_Trm_markers.pdf"),p7,width = 8,height = 14)

#>>> Temra

# CD3+,CD8+,CD27-,CD28-,CD45RA+,CD45RO-,CD57+,CD62L-,CD95+,CD127-,CD197-

# CD3D/E/G,CD8A/B,CD27,CD28,B3GAT1,SELL,FAS,IL7R,CCR7

Temra_markers <- c("CD3D","CD3E","CD3G","CD8A","CD8B","CD27","CD28","B3GAT1","SELL","FAS","IL7R","CCR7")
p8 <- StackedVlnPlot(obj = CD8T, features = Temra_markers,cols=pal)
ggsave(paste(wkd,"CD8_Temra_markers.pdf"),p8,width = 8,height = 14)

#>>> Tex, exhausted T cell

# CD45+,CD8+,CD3+,CD96+(TACTILE),CD152+(CTLA4),CD160+(NK1),CD223+(LAG-3),
# CD244+(2B4),CD272+(BTLA),CD279+(PD1),CD366+(TIM-3),EOMES+,ICOS+,TIGIT+
# VISTA+,Tox+

# PTPRC,CD8A/B,CD3D/E/G,CD96,CTLA4,CD160,LAG3,CD244,BTLA,PDCD1,HAVCR2,EOMES,
# ICOS,TIGIT,C10orf54(VISTA/VSIR),TOX

Tex_markers <- c("PTPRC","CD3D","CD3E","CD3G","CD8A","CD8B","CD96","CTLA4","CD160","LAG3","CD244","BTLA","PDCD1",
                 "HAVCR2","EOMES","ICOS","TIGIT","C10orf54","TOX")
p9 <- StackedVlnPlot(obj = CD8T, features = Tex_markers,cols=pal)
ggsave(paste(wkd,"CD8_Tex_markers.pdf"),p9,width = 8,height = 16)

#>>> CTLs

# CD45+,CD3+,CD8+,CD107a+(LAMP-1),CD178+(FasL),GranzymeB+,IFN-γ+,Perforin+

# PTPRC,CD3D/E/G,CD8A/B,LAMP1,FASLG,GZMB,IFNG,PRF1
CTL_markers <- c("PTPRC","CD3D","CD3E","CD3G","CD8A","CD8B","LAMP1","FASLG","GZMB","IFNG","PRF1")
p10 <- StackedVlnPlot(obj = CD8T, features = CTL_markers,cols=pal)
ggsave(paste(wkd,"CD8_CTLs_markers.pdf"),p10,width = 8,height = 14)


## map new T subclusters
curr_ids <- c(0:9)
new_ids <- c("Unknown","Unknown","Tex","CTLs","Unknown","Unknown","Tcm","Teff","Ttm","Unknown")
CD8T$subcluster <- plyr::mapvalues(CD8T$seurat_clusters,from = curr_ids,to = new_ids)

pal <- c("#828281","#BF0032","#0066A5","#008956","#E35822","#604E97")

pCD8 <- DimPlot(CD8T,group.by = "subcluster",cols = pal,pt.size = 1) + ggtitle("CD8T subclusters")
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/CD8T_subclusters.pdf",pCD8,width = 6,height = 4.6)

saveRDS(CD8T,file = "rds/without_PVTT_MLN/tumor/sub_T/CD8T.Tsub.assigned.rds")


# >>>>>>>>>>>>>>>> CD4 MARKERS <<<<<<<<<<<<<<<<<<<<<<<<
DefaultAssay(CD4T) <- "RNA"
wkd <- "result/without_PVTT_MLN/tumor/sub_T/2.new/"
pal <- c("#786140","#8C0603","#B2712B","#B48521","#C6EEC4","#DCD3B1","#A7BBA5","#0202BB","#5EC6B9","#4E8A97",
         "#726DBA","#090795","#AD66B6","#B71CB4","#A01153","#008000")

#>>> Tn, naive T cell

# CD45+,CD3+,CD4+,CD27+,CD28+,CD45RA+,CD45RO-,CD57-,CD62L+(L-Selectin),CD69-,CD95-(FasR),
# CD127+,CD197+(CCR7)

# PTPRC,CD3D/E/G,CD4,CD27,CD28,B3GAT1,SELL,CD69,FAS,IL7R,CCR7

Tn_markers <- c("PTPRC","CD3D","CD3E","CD3G","CD4","CD27","CD28","B3GAT1","SELL","CD69","FAS","IL7R","CCR7")
p1 <- StackedVlnPlot(obj = CD4T, features = Tn_markers,cols=pal)
ggsave(paste(wkd,"CD4_Tn_markers.pdf"),p1,width = 8,height = 14)


#>>> Tscm, stem cell memory T cell

# CD45+,CD3+,CD4+,CD27+,CD28+,CD45RA+,CD45RO-,CD57-,CD62L+(L-Selectin),CD95+(FasR),
# CD127+,CD197+(CCR7)

# PTPRC,CD3D/E/G,CD4,CD27,CD28,B3GAT1,SELL,FAS,IL7R,CCR7

Tscm_markers <- c("PTPRC","CD3D","CD3E","CD3G","CD4","CD27","CD28","B3GAT1","SELL","FAS","IL7R","CCR7")
p2 <- StackedVlnPlot(obj = CD4T, features = Tscm_markers,cols=pal)
ggsave(paste(wkd,"CD4_Tscm_markers.pdf"),p2,width = 8,height = 14)

#>>> Tcm, central memory T cell

# CD45+,CD3+,CD4+,CD27+,CD28+,CD45RA-,CD45RO+,CD57-,CD62+(L-Selectin),CD95+(FasR),
# CD127+,CD197+(CCR7)

# PTPRC,CD3D/E/G,CD4,CD27,CD28,B3GAT1,SELL,FAS,IL7R,CCR7

Tcm_markers <- c("PTPRC","CD3D","CD3E","CD3G","CD4","CD27","CD28","B3GAT1","SELL","FAS","IL7R","CCR7")
p3 <- StackedVlnPlot(obj = CD4T, features = Tcm_markers,cols=pal)
ggsave(paste(wkd,"CD4_Tcm_markers.pdf"),p3,width = 8,height = 14)


#>>> Teff, effector T cell

# CD45+,CD4+,CD3+,CD25+(IL2RA),CD69+,CD95+(FasR),CD134+(OX40),CD137+(4-1BB),
# CD154+(CD40L),Ki-67+,KLRG1+,CD62L-,CCR7+

# PTPRC,CD4,CD3D/E/G,IL2RA,CD69,FAS,TNFRSF4,TNFRSF9,CD40LG,MKI67,KLRG1,SELL,CCR7

Teff_markers <- c("PTPRC","CD4","CD3D","CD3E","CD3G","IL2RA","CD69","FAS","TNFRSF4","TNFRSF9","CD40LG","MKI67","KLRG1","SELL","CCR7")
p4 <- StackedVlnPlot(obj = CD4T, features = Teff_markers,cols=pal)
ggsave(paste(wkd,"CD4_Teff_markers.pdf"),p4,width = 8,height = 14)


#>>> Ttm

# CD45+,CD3+,CD4+,CD27+,CD28+,CD45RA-,CD45RO+,CD57+,CD62L-(L-Selectin),
# CD95+(FasR),CD127+,CD197-(CCR7)

# PTPRC,CD3D/E/G,CD4,CD27,CD28,B3GAT1,SELL,FAS,IL7R,CCR7

Ttm_markers <- c("PTPRC","CD3D","CD3E","CD3G","CD4","CD27","CD28","B3GAT1","SELL","FAS","IL7R","CCR7")
p5 <- StackedVlnPlot(obj = CD4T, features = Ttm_markers,cols=pal)
ggsave(paste(wkd,"CD4_Ttm_markers.pdf"),p5,width = 8,height = 14)


#>>> Tem, effector memory cell

# CD45+,CD3+,CD4+,CD27+/-,CD28-,CD45RA-,CD45RO+,CD57+,CD62L-,CD95+,CD127+/-,CD197-

# PTPRC,CD3D/E/G,CD4,CD27,CD28,B3GAT1,SELL,FAS,IL7R,CCR7

Tem_markers <- c("PTPRC","CD3D","CD3E","CD3G","CD4","CD27","CD28","B3GAT1","SELL","FAS","IL7R","CCR7")
p6 <- StackedVlnPlot(obj = CD4T, features = Tem_markers,cols=pal)
ggsave(paste(wkd,"CD4_Tem_markers.pdf"),p6,width = 8,height = 14)


#>>> Trm, Tissue-resident memory cell

# CD45+,CD4+,CD3+,CD69+,CD103+,CTLA4+

# PTPRC,CD4,CD3D/E/G,CD69,ITGAE,CTLA4

Trm_markers <- c("PTPRC","CD4","CD3D","CD3E","CD3G","CD69","ITGAE","CTLA4")
p7 <- StackedVlnPlot(obj = CD4T, features = Trm_markers,cols=pal)
ggsave(paste(wkd,"CD4_Trm_markers.pdf"),p7,width = 8,height = 14)

#>>> Tex, exhausted T cell

# CD45+,CD4+,CD3+,CD96+,CD152+,CD160+,CD223+,CD244+,CD272+,CD279+,
# CD366+,EOMES+,ICOS+,TIGIT,VISTA+,Tox+

# PTPRC,CD4,CD3D/E/G,CD96,CTLA4,CD160,LAG3,CD244,BTLA,PDCD1,HAVCR2,
# EOMES,ICOS,TIGIT,C10orf54(VISTA/VSIR),TOX

Tex_markers <- c("PTPRC","CD4","CD3D","CD3E","CD3G","CD96","CTLA4","CD160","LAG3","CD244","BTLA","PDCD1","HAVCR2",
                 "EOMES","ICOS","TIGIT","C10orf54","TOX")
p8 <- StackedVlnPlot(obj = CD4T, features = Tex_markers,cols=pal)
ggsave(paste(wkd,"CD4_Tex_markers.pdf"),p8,width = 8,height = 16)


#>>> Temra

# CD45+,CD3+,CD4+,CD27-,CD28-,CD45RA+,CD45RO-,CD57+,CD62L-,CD95+,CD127-,CD197-

# PTPRC,CD3D/E/G,CD4,CD27,CD28,B3GAT1,SELL,FAS,IL7R,CCR7

Temra_markers <- c("PTPRC","CD3D","CD3E","CD3G","CD4","CD27","CD28","B3GAT1","SELL","FAS","IL7R","CCR7")
p9 <- StackedVlnPlot(obj = CD4T, features = Temra_markers,cols=pal)
ggsave(paste(wkd,"CD4_Temra_markers.pdf"),p9,width = 8,height = 14)

#>>> Treg

# CD45+,CD4+,CD3+,Foxp3+,CD25+,Helicos+,CD127-,CTLA4+

# PTPRC,CD4,CD3D/E/G,FOXP3,IL2RA,S100A8,IL7R,CTLA4

Treg_markers <- c("PTPRC","CD4","CD3D","CD3E","CD3G","FOXP3","IL2RA","S100A8","IL7R","CTLA4")
p10 <- StackedVlnPlot(obj = CD4T, features = Treg_markers,cols=pal)
ggsave(paste(wkd,"CD4_Treg_markers.pdf"),p10,width = 8,height = 14)

#>>> NKT

# CD45+,CD3+,CD4-,CD8-,CD16+,CD56+(NCAM),CD57+,CD161+

# PTPRC,CD3D/E/G,CD4,CD8A/B,FCGR3A,NCAM1,B3GAT1,KLRB1

NKT_markers <- c("PTPRC","CD3D","CD3E","CD3G","CD4","CD8A","CD8B","FCGR3A","NCAM1","B3GAT1","KLRB1")
p11 <- StackedVlnPlot(obj = CD4T, features = NKT_markers,cols=pal)
ggsave(paste(wkd,"CD4_NKT_markers.pdf"),p11,width = 8,height = 12)


## map new T subclusters
curr_ids <- c(0:14)
new_ids <- c("Unknown","Unknown","Unknown","Tem","Tem","Tem","Tem","Tscm","Treg","Tcm","Tex","Temra","NKT","Teff","Unknown")
CD4T$subcluster <- plyr::mapvalues(CD4T$seurat_clusters,from = curr_ids, to = new_ids)
pal <- c("#828281","#00441B","#46A040","#00AF99","#FFC179","#98D9E9","#F6313E","#FFA300","#C390D4")
pCD4 <- DimPlot(CD4T,group.by = "subcluster",cols = pal,pt.size = 1,reduction = "tsne") + ggtitle("CD4T subclusters")
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/CD4T_subclusters.pdf",pCD4,width = 6,height = 4.6)

saveRDS(CD4T, file = "rds/without_PVTT_MLN/tumor/sub_T/CD4T.Tsub.assigned.rds")


## >>>>> Requirement II:

## 在CD3+T高浸润和低浸润的患者中分析哪些新亚群富集在CD3+T高浸润肿瘤中，哪些新亚群富集在CD3+T低浸润肿瘤中

Tcell <- readRDS("rds/without_PVTT_MLN/tumor/sub_T/Tcell.new.DimRed.rds")

p1 <- DimPlot(Tcell, group.by = "Group1",cols = c("#BF0032","#0066A5")) + 
  theme(plot.title = element_blank(),
        legend.position = "top")

pal <- c("#BF0032","#828281","#008956","#892C16","#875691","#2A3C26",
         "#0066A5","#644422","#F3C300","#222222","#E68FAC","#E35822","#C3B381")
p2 <- DimPlot(Tcell,cols = pal,group.by = "novelCluster") + theme(plot.title = element_blank())

pc <- p1 | p2
ggsave("result/without_PVTT_MLN/tumor/sub_T/IR_vs_novelClusters.pdf",pc,width = 13,height = 6)

CHGA <- FeaturePlot(Tcell, features = "CHGA")
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/CHGA.pdf",CHGA,width = 5,height = 4.5)

CXCL13 <- FeaturePlot(Tcell, features = "CXCL13")
ggsave("result/without_PVTT_MLN/tumor/sub_T/2.new/CXCL13.pdf",CXCL13,width = 5,height = 4.5)


#1.对于在CD3+T高浸润/低浸润肿瘤中富集的新CD4+T细胞亚型，分析该亚型的特异性基因表达差异(相对其他15种T细胞），并分析该群细胞的
#  特征基因与Th1，Th2，Treg，Th17，Tfh，和CD4+Teff, CD4+Tem，CD4+Tcm，CD4+Tex细胞的特征基因相似性

rm(list = ls());gc()
Tcell <- readRDS("rds/without_PVTT_MLN/tumor/sub_T/Tcell.DimRed.rds")
# CD8T
CD8T <- readRDS("rds/without_PVTT_MLN/tumor/sub_T/CD8T.Tsub.assigned.rds")
## map new T subclusters
curr_ids <- c(0:9)
new_ids <- c("Unknown","Unknown","CD8_Tex","CD8_CTLs","Unknown","Unknown","CD8_Tcm","CD8_Teff","CD8_Ttm","Unknown")
CD8T$subcluster <- plyr::mapvalues(CD8T$seurat_clusters,from = curr_ids,to = new_ids)

# CD4T
CD4T <- readRDS("rds/without_PVTT_MLN/tumor/sub_T/CD4T.Tsub.assigned.rds")
curr_ids <- c(0:14)
new_ids <- c("Unknown","Unknown","Unknown","CD4_Tem","CD4_Tem","CD4_Tem","CD4_Tem","CD4_Tscm","CD4_Treg","CD4_Tcm","CD4_Tex",
             "CD4_Temra","CD4_NKT","CD4_Teff","Unknown")
CD4T$subcluster <- plyr::mapvalues(CD4T$seurat_clusters,from = curr_ids, to = new_ids)

# combine CD4T subtypes and CD8T subtypes into T cell seurat object
CD8.meta <- data.frame(row.names = colnames(CD8T), subtype = as.character(CD8T$subcluster))
CD4.meta <- data.frame(row.names = colnames(CD4T), subtype = as.character(CD4T$subcluster))

CDT.meta <- rbind(CD8.meta, CD4.meta)
CDT.meta$subtype <- as.character(CDT.meta$subtype)

Tcell$subtype <- plyr::mapvalues(rownames(Tcell@meta.data), from = rownames(CDT.meta), to = CDT.meta$subtype)

p <- DimPlot(Tcell, group.by = "subtype") + theme(plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/sub_T/Tcell.subtypes.pdf",p,width = 6,height = 5)

saveRDS(Tcell, file = "rds/without_PVTT_MLN/tumor/sub_T/Tcell.subtypes.rds")

## CD8+T细胞特征基因与CTL，Teff, Tem，Tcm，Tex细胞的特征基因相似性

### CD8+T cell:

# CD8_CTLs: CD45+CD3+ CD8+ CD107a+ (LAMP-1) CD178+ (FasL) Granzyme B+ IFN-γ+ Perforin+
CD8_CTLs <- c("PTPRC","CD3D","CD3E","CD3G","LAMP1","FAS","GZMB","IFNG","PRF1")
# CD8_Teff: CD45+CD8+CD3+CD25+ (IL2RA) CD69+ CD95+ (FasR) CD134+ (OX40) CD137+ (4-1BB) CD154+ (CD40L) Ki-67+ KLRG1+CD62L- CCR7+
CD8_Teff <- c("PTPRC","CD8A","CD8B","IL2RA","CD69","FAS","TNFRSF4","TNFRSF9","CD40LG","MKI67","KLRG1","SELL","CCR7")
# CD8_Tem: CD3+ CD8+ CD27+/– CD28– CD45RA– CD45RO+ CD57+ CD62L– (L-Selectin) CD95+ (FasR) CD127+/– CD197– (CCR7)
CD8_Tem <- c("CD3D","CD3E","CD3G","CD27","CD8A","CD8B","CD28","B3GAT1","SELL","FAS","IL7R","CCR7")
# CD8_Tcm: CD45+CD3+ CD8+ CD27+ CD28+ CD45RA– CD45RO+ CD57– CD62L+ (L-Selectin) CD95+ (FasR) CD127+ CD197+ (CCR7)
CD8_Tcm <- c("PTPRC","CD3D","CD3E","CD3G","CD8A","CD8B","CD27","CD28","B3GAT1","SELL","FAS","IL7R","CCR7")
# CD8_Tex: CD45+CD8+CD3+CD96+ (TACTILE) CD152+ (CTLA-4) CD160+ (NK1) CD223+ (LAG-3) CD244+ (2B4) CD272+ (BTLA) CD279+ (PD1) CD366+ (TIM-3) EOMES+ ICOS+ TIGIT+ VISTA+Tox+
CD8_Tex <- c("PTPRC","CD8A","CD8B","CD3D","CD3E","CD3G","CD96","CTLA4","CD160","LAG3","CD244","BTLA","PDCD1","HAVCR2","EOMES","ICOS","TIGIT","C10orf54","TOX")

refCells_CD8 <- list(CD8_CTLs,CD8_Teff,CD8_Tem,CD8_Tcm,CD8_Tex)
names(refCells_CD8) <- c("CD8_CTLs","CD8_Teff","CD8_Tem","CD8_Tcm","CD8_Tex")

## CD4+T细胞特征基因与Th1，Th2，Treg，Th17，Tfh，和CD4+Teff, CD4+Tem，CD4+Tcm，CD4+Tex细胞的特征基因相似性

### CD4+T cell:

##1. CD4_Th1: CD45+CD3+ CD4+Tbet+ IFN-γ+CXCR3+
# genes: PTPRC, CD3D, CD3E, CD3G, CD4, TBX21, IFNG, CXCR3
CD4_Th1 <- c("PTPRC","CD3D","CD3E","CD3G","CD4","TBX21","IFNG","CXCR3")
##2. CD4_Th2: CRTH2+,CD4+,CRs+,CCR3+,CCR8+
# genes: PTGDR2,CD4,CR1,CCR3,CCR8
CD4_Th2 <- c("GPR44","CD4","CR1","CCR3","CCR8")
CD4_Th2_sc <- c("PTGDR2","CD4","CR1","CCR3","CCR8")
##3. CD4_Treg: CD45+CD4+CD3+Foxp3+CD25+Helios+CD127-CTLA4+
# genes: PTPRC,CD4,CD3D,CD3E,CD3G,FOXP3,IL2RA,IKZF2,IL7R,CTLA4
CD4_Treg <- c("PTPRC","CD4","CD3D","CD3E","CD3G","FOXP3","IL2RA","IKZF2","IL7R","CTLA4")
##4. CD4_Th17: CD45+,CD4+,CCR4,CCR7
# genes: PTPRC, CD4, CCR4, CCR7
CD4_Th17 <- c("PTPRC","CD4","CCR4","CCR7")
##5. CD4_Tfh: CD45+,CXCR5,BCL6+,Tbet-,GATA3-,RORγt-
# genes: PTPRC, CXCR5, BCL6, TBX21, GATA3, RORC
CD4_Tfh <- c("PTPRC","CXCR5","BCL6","TBX21","GATA3","RORC")
##6. CD4_Teff: CD45+CD4+CD3+CD25+ (IL2RA) CD69+ CD95+ (FasR) CD134+ (OX40) CD137+ (4-1BB) CD154+ (CD40L) Ki-67+ KLRG1+CD62L- CCR7+
# genes: PTPRC,CD4,CD3D,CD3E,CD3G,IL2RA,CD69,FAS,TNFRSF4,TNFRSF9,CD154,MKI67,KLRG1,SELL,CCR7
CD4_Teff <- c("PTPRC","CD4","CD3D","CD3E","CD3G","IL2RA","CD69","FAS","TNFRSF4","TNFRSF9","CD40LG","MKI67","KLRG1","SELL","CCR7")
##7. CD4_Tem: CD45+CD3+ CD4+ CD27+/– CD28– CD45RA– CD45RO+ CD57+ CD62L– (L-Selectin) CD95+ (FasR) CD127+/– CD197– (CCR7)
# genes: PTPRC,CD3D,CD3E,CD3G,CD4,CD27,CD28,B3GAT1,SELL,FAS,IL7R,CCR7
CD4_Tem <- c("PTPRC","CD3D","CD3E","CD3G","CD4","CD27","CD28","B3GAT1","SELL","FAS","IL7R","CCR7")
##8. CD4_Tcm: CD45+CD3+ CD4+ CD27+ CD28+ CD45RA– CD45RO+ CD57– CD62L+ (L-Selectin) CD95+ (FasR) CD127+ CD197+ (CCR7)
# genes: PTPRC,CD3D,CD3E,CD3G,CD4,CD27,CD28,B3GAT1,SELL,FAS,IL7R,CCR7
CD4_Tcm <- c("PTPRC","CD3D","CD3E","CD3G","CD4","CD27","CD28","B3GAT1","SELL","FAS","IL7R","CCR7")
##9. CD4_Tex: CD45+CD4+CD3+CD96+ (TACTILE) CD152+ (CTLA-4) CD160+ (NK1) CD223+ (LAG-3) CD244+ (2B4) CD272+ (BTLA) CD279+ (PD1) CD366+ (TIM-3) EOMES+ ICOS+ TIGIT+ VISTA+Tox+
# genes: PTPRC,CD4,CD3D,CD3E,CD3G,CD96,CTLA4,CD160,LAG3,CD244,BTLA,PDCD1,HAVCR2,EOMES,ICOS,TIGIT,VSIR,TOX
CD4_Tex <- c("PTPRC","CD4","CD3D","CD3E","CD3G","CD96","CTLA4","CD160","LAG3","CD244","BTLA","PDCD1","HAVCR2","EOMES","ICOS","TIGIT","C10orf54","TOX")

refCells <- list(CD4_Th1,CD4_Th2,CD4_Treg,CD4_Th17,CD4_Tfh,CD4_Teff,CD4_Tem,CD4_Tcm,CD4_Tex)
refCells_sc <- list(CD4_Th1,CD4_Th2_sc,CD4_Treg,CD4_Th17,CD4_Tfh,CD4_Teff,CD4_Tem,CD4_Tcm,CD4_Tex)
names(refCells) <- c("CD4_Th1","CD4_Th2","CD4_Treg","CD4_Th17","CD4_Tfh","CD4_Teff","CD4_Tem","CD4_Tcm","CD4_Tex")
names(refCells_sc) <- c("CD4_Th1","CD4_Th2","CD4_Treg","CD4_Th17","CD4_Tfh","CD4_Teff","CD4_Tem","CD4_Tcm","CD4_Tex")

library(Seurat)
# Load in dataset
Tcell <- readRDS("rds/without_PVTT_MLN/tumor/sub_T/Tcell.new.DimRed.rds")

##>>>>  CD3 High Group: CD4_XIST,CD4_RBP7,CD8_CXCL13,CD4_VEGFB,CD4_VCX3A,CD4_PLCG2,CD4_IL1R2,CD8_FCGR3A,CD4_FGR,CD4_CRP

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

# >>>> CD4_XIST

target <- cellEmbeddings(object = Tcell, cellCluster = "CD4_XIST")
target <- target[target$UMAP_1 > 5 & target$UMAP_2 < -1,]
# cells highlight
p_CD4_XIST <- DimPlot(Tcell,group.by = "novelCluster",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
               cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())

ggsave("result/without_PVTT_MLN/tumor/CD3_High/CD4_XIST/CD4_XIST.highlight.pdf",p_CD4_XIST,device = "pdf",width = 6,height = 6)

# FIND MARKERS
Idents(Tcell) <- Tcell$novelCluster
markers <- FindMarkers(Tcell,only.pos = TRUE,ident.1 = "CD4_XIST")
genes <- rownames(markers)
for(i in seq_along(genes)){
  p <- FeaturePlot(Tcell,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3_High/CD4_XIST/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}


# signature Correlation Analysis

load("refdb/TCGA_LIHC.rnaseq.RData")

for(i in seq_along(refCells)){
  genes <- c("XIST",refCells_sc[[i]])
  avgExpPerSample <- function(object, gene){
    expr <- GetAssayData(object, assay = "RNA", slot = "data")[gene,]
    expr <- as.data.frame(expr, stringAsFactors = F)
    expr.T <- data.frame(t(expr))
    expr.T$Patient <- plyr::mapvalues(rownames(expr.T), from = colnames(object), to = object$Sample)
    colnames(expr.T) <- c(gene, "Sample")
    
    patient <- NULL
    for(j in 1:(ncol(expr.T) - 1)){
      patient[j] <- data.frame(tapply(expr.T[,j], expr.T$Sample, sum))
    }
    names(patient) <- colnames(expr.T)[-ncol(expr.T)]
    patient <- data.frame(patient)
    rownames(patient) <- names(tapply(expr.T[,1], expr.T$Sample, sum))
    Patient_sample <- as.data.frame(table(expr.T$Sample), stringAsFactors = F)
    patient$Sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- as.data.frame(patient, stringAsFactors = F)
    patient$Sample <- as.numeric(patient$Sample)
    for(k in 1:(ncol(patient) - 1)){
      for(l in 1:nrow(patient)){
        patient[l,k] <- patient[l,k]/patient$Sample[l]
      }
    }
    patient$Sample <- NULL
    colnames(patient) <- gene
    return(patient)
  }
  patient <- avgExpPerSample(Tcell,genes)

  if (i == 1) {
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD4_XIST/avgExpPerPatient.xlsx",sheetName = names(refCells)[1])
  }else{
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD4_XIST/avgExpPerPatient.xlsx",sheetName = names(refCells)[i],append = TRUE)
  }
  
  # singlecell
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD4_XIST/corr/singlecell/"

  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  for(m in 2:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[1], y = colnames(patient)[m])) +
      geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(patient)[m],"_vs_",colnames(patient)[1],".pdf"),p2,width = 6,height = 6)
  }
  
  # TCGA
  colnames(patient) <- c("XIST",refCells[[i]])
  expr.LIHC <- LIHC.rnaseq[,colnames(patient)]
  expr.norm <- log2(expr.LIHC + 1)
  
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD4_XIST/corr/TCGA/"
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  
  for(n in 2:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x=colnames(expr.norm)[1], y=colnames(expr.norm)[n])) +
      geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(expr.norm)[n],"_vs_",colnames(expr.norm)[1],".pdf"),p2,width = 6,height = 6)
  }
}

# >>>> CD4_RBP7

target <- cellEmbeddings(object = Tcell, cellCluster = "CD4_RBP7")
# cells highlight
p_CD4_RBP7 <- DimPlot(Tcell,group.by = "novelCluster",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
                      cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/CD3_High/CD4_RBP7/CD4_RBP7.highlight.pdf",p_CD4_RBP7,device = "pdf",width = 6,height = 6)

Idents(Tcell) <- "novelCluster"
RBP7_markers <- FindMarkers(Tcell, ident.1 = "CD4_RBP7",only.pos = TRUE)

genes <- rownames(markers)
for(i in seq_along(genes)){
  p <- FeaturePlot(Tcell,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3_High/CD4_RBP7/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}


# signature Correlation Analysis

for(i in seq_along(refCells)){
  genes <- c("RBP7",refCells_sc[[i]])
  avgExpPerSample <- function(object, gene){
    expr <- GetAssayData(object, assay = "RNA", slot = "data")[gene,]
    expr <- as.data.frame(expr, stringAsFactors = F)
    expr.T <- data.frame(t(expr))
    expr.T$Patient <- plyr::mapvalues(rownames(expr.T), from = colnames(object), to = object$Sample)
    colnames(expr.T) <- c(gene, "Sample")
    
    patient <- NULL
    for(j in 1:(ncol(expr.T) - 1)){
      patient[j] <- data.frame(tapply(expr.T[,j], expr.T$Sample, sum))
    }
    names(patient) <- colnames(expr.T)[-ncol(expr.T)]
    patient <- data.frame(patient)
    rownames(patient) <- names(tapply(expr.T[,1], expr.T$Sample, sum))
    Patient_sample <- as.data.frame(table(expr.T$Sample), stringAsFactors = F)
    patient$Sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- as.data.frame(patient, stringAsFactors = F)
    patient$Sample <- as.numeric(patient$Sample)
    for(k in 1:(ncol(patient) - 1)){
      for(l in 1:nrow(patient)){
        patient[l,k] <- patient[l,k]/patient$Sample[l]
      }
    }
    patient$Sample <- NULL
    colnames(patient) <- gene
    return(patient)
  }
  patient <- avgExpPerSample(Tcell,genes)
  
  if (i == 1) {
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD4_XIST/avgExpPerPatient.xlsx",sheetName = names(refCells)[1])
  }else{
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD4_RBP7/avgExpPerPatient.xlsx",sheetName = names(refCells)[i],append = TRUE)
  }
  
  # singlecell
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD4_RBP7/corr/singlecell/"
  
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  for(m in 2:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[1], y = colnames(patient)[m])) +
      geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(patient)[m],"_vs_",colnames(patient)[1],".pdf"),p2,width = 6,height = 6)
  }
  
  # TCGA
  colnames(patient) <- c("RBP7",refCells[[i]])
  expr.LIHC <- LIHC.rnaseq[,colnames(patient)]
  expr.norm <- log2(expr.LIHC + 1)
  
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD4_RBP7/corr/TCGA/"
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  
  for(n in 2:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x=colnames(expr.norm)[1], y=colnames(expr.norm)[n])) +
      geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(expr.norm)[n],"_vs_",colnames(expr.norm)[1],".pdf"),p2,width = 6,height = 6)
  }
}


# >>>> CD8_CXCL13

target <- cellEmbeddings(object = Tcell, cellCluster = "CD8_CXCL13")
target <- target[target$UMAP_2 > 3,]
# cells highlight
p_CD8_CXCL13 <- DimPlot(Tcell,group.by = "novelCluster",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
                      cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/CD3_High/CD8_CXCL13/CD8_CXCL13.highlight.pdf",p_CD8_CXCL13,device = "pdf",width = 6,height = 6)

# markers
Idents(Tcell) <- "novelCluster"
markers <- FindMarkers(Tcell, ident.1 = "CD8_CXCL13",only.pos = TRUE)

genes <- rownames(markers)
for(i in seq_along(genes)){
  p <- FeaturePlot(Tcell,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3_High/CD8_CXCL13/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}

# signature Correlation Analysis

for(i in seq_along(refCells_CD8)){
  genes <- c("CXCL13",refCells_CD8[[i]])
  avgExpPerSample <- function(object, gene){
    expr <- GetAssayData(object, assay = "RNA", slot = "data")[gene,]
    expr <- as.data.frame(expr, stringAsFactors = F)
    expr.T <- data.frame(t(expr))
    expr.T$Patient <- plyr::mapvalues(rownames(expr.T), from = colnames(object), to = object$Sample)
    colnames(expr.T) <- c(gene, "Sample")
    
    patient <- NULL
    for(j in 1:(ncol(expr.T) - 1)){
      patient[j] <- data.frame(tapply(expr.T[,j], expr.T$Sample, sum))
    }
    names(patient) <- colnames(expr.T)[-ncol(expr.T)]
    patient <- data.frame(patient)
    rownames(patient) <- names(tapply(expr.T[,1], expr.T$Sample, sum))
    Patient_sample <- as.data.frame(table(expr.T$Sample), stringAsFactors = F)
    patient$Sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- as.data.frame(patient, stringAsFactors = F)
    patient$Sample <- as.numeric(patient$Sample)
    for(k in 1:(ncol(patient) - 1)){
      for(l in 1:nrow(patient)){
        patient[l,k] <- patient[l,k]/patient$Sample[l]
      }
    }
    patient$Sample <- NULL
    colnames(patient) <- gene
    return(patient)
  }
  patient <- avgExpPerSample(Tcell,genes)
  
  if (i == 1) {
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD8_CXCL13/avgExpPerPatient.xlsx",sheetName = names(refCells_CD8)[1])
  }else{
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD8_CXCL13/avgExpPerPatient.xlsx",sheetName = names(refCells_CD8)[i],append = TRUE)
  }
  
  # singlecell
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD8_CXCL13/corr/singlecell/"
  
  if(!file.exists(paste0(pwd,names(refCells_CD8[i])))){
    dir.create(paste0(pwd,names(refCells_CD8)[i]))
  }
  for(m in 2:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[1], y = colnames(patient)[m])) +
      geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells_CD8)[i],"/",colnames(patient)[m],"_vs_",colnames(patient)[1],".pdf"),p2,width = 6,height = 6)
  }
  
  # TCGA
  colnames(patient) <- c("CXCL13",refCells_CD8[[i]])
  expr.LIHC <- LIHC.rnaseq[,colnames(patient)]
  expr.norm <- log2(expr.LIHC + 1)
  
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD8_CXCL13/corr/TCGA/"
  if(!file.exists(paste0(pwd,names(refCells_CD8[i])))){
    dir.create(paste0(pwd,names(refCells_CD8)[i]))
  }
  
  for(n in 2:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x=colnames(expr.norm)[1], y=colnames(expr.norm)[n])) +
      geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells_CD8)[i],"/",colnames(expr.norm)[n],"_vs_",colnames(expr.norm)[1],".pdf"),p2,width = 6,height = 6)
  }
}

#>>>>>> CD4_VEGFB

target <- cellEmbeddings(object = Tcell, cellCluster = "CD4_VEGFB")
target <- target[target$UMAP_1 > 0 & target$UMAP_1 < 5 & target$UMAP_2 < 1.5,]
# cells highlight
p_CD4_VEGFB <- DimPlot(Tcell,group.by = "novelCluster",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
                        cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/CD3_High/CD4_VEGFB/CD4_VEGFB.highlight.pdf",p_CD4_VEGFB,device = "pdf",width = 6,height = 6)

# markers
Idents(Tcell) <- "novelCluster"
markers <- FindMarkers(Tcell, ident.1 = "CD4_VEGFB",only.pos = TRUE)

genes <- rownames(markers)
for(i in seq_along(genes)){
  p <- FeaturePlot(Tcell,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3_High/CD4_VEGFB/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}

# signature Correlation Analysis

for(i in seq_along(refCells)){
  genes <- c("VEGFB",refCells_sc[[i]])
  avgExpPerSample <- function(object, gene){
    expr <- GetAssayData(object, assay = "RNA", slot = "data")[gene,]
    expr <- as.data.frame(expr, stringAsFactors = F)
    expr.T <- data.frame(t(expr))
    expr.T$Patient <- plyr::mapvalues(rownames(expr.T), from = colnames(object), to = object$Sample)
    colnames(expr.T) <- c(gene, "Sample")
    
    patient <- NULL
    for(j in 1:(ncol(expr.T) - 1)){
      patient[j] <- data.frame(tapply(expr.T[,j], expr.T$Sample, sum))
    }
    names(patient) <- colnames(expr.T)[-ncol(expr.T)]
    patient <- data.frame(patient)
    rownames(patient) <- names(tapply(expr.T[,1], expr.T$Sample, sum))
    Patient_sample <- as.data.frame(table(expr.T$Sample), stringAsFactors = F)
    patient$Sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- as.data.frame(patient, stringAsFactors = F)
    patient$Sample <- as.numeric(patient$Sample)
    for(k in 1:(ncol(patient) - 1)){
      for(l in 1:nrow(patient)){
        patient[l,k] <- patient[l,k]/patient$Sample[l]
      }
    }
    patient$Sample <- NULL
    colnames(patient) <- gene
    return(patient)
  }
  patient <- avgExpPerSample(Tcell,genes)
  
  if (i == 1) {
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/PARTV.1/CD3_High/CD4_VEGFB/avgExpPerPatient.xlsx",sheetName = names(refCells)[1])
  }else{
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/PARTV.1/CD3_High/CD4_VEGFB/avgExpPerPatient.xlsx",sheetName = names(refCells)[i],append = TRUE)
  }
  
  # singlecell
  pwd = "result/without_PVTT_MLN/tumor/PARTV.1/CD3_High/CD4_VEGFB/corr/singlecell/"
  
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  for(m in 2:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[1], y = colnames(patient)[m])) +
      geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(patient)[m],"_vs_",colnames(patient)[1],".pdf"),p2,width = 6,height = 6)
  }
  
  # TCGA
  colnames(patient) <- c("VEGFB",refCells[[i]])
  expr.LIHC <- LIHC.rnaseq[,colnames(patient)]
  expr.norm <- log2(expr.LIHC + 1)
  
  pwd = "result/without_PVTT_MLN/tumor/PARTV.1/CD3_High/CD4_VEGFB/corr/TCGA/"
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  
  for(n in 2:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x=colnames(expr.norm)[1], y=colnames(expr.norm)[n])) +
      geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(expr.norm)[n],"_vs_",colnames(expr.norm)[1],".pdf"),p2,width = 6,height = 6)
  }
}


#>>>>>>> CD4_VCX3A

target <- cellEmbeddings(object = Tcell, cellCluster = "CD4_VCX3A")
target <- target[target$UMAP_1 < -6,]
# cells highlight
p_CD4_VCX3A <- DimPlot(Tcell,group.by = "novelCluster",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
                        cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/CD3_High/CD4_VCX3A/CD4_VCX3A.highlight.pdf",p_CD4_VCX3A,device = "pdf",width = 6,height = 6)

# markers
Idents(Tcell) <- "novelCluster"
markers <- FindMarkers(Tcell, ident.1 = "CD4_VCX3A",only.pos = TRUE)

genes <- rownames(markers)
for(i in seq_along(genes)){
  p <- FeaturePlot(Tcell,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3_High/CD4_VCX3A/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}

# signature Correlation Analysis

for(i in seq_along(refCells)){
  genes <- c("VCX3A",refCells_sc[[i]])
  avgExpPerSample <- function(object, gene){
    expr <- GetAssayData(object, assay = "RNA", slot = "data")[gene,]
    expr <- as.data.frame(expr, stringAsFactors = F)
    expr.T <- data.frame(t(expr))
    expr.T$Patient <- plyr::mapvalues(rownames(expr.T), from = colnames(object), to = object$Sample)
    colnames(expr.T) <- c(gene, "Sample")
    
    patient <- NULL
    for(j in 1:(ncol(expr.T) - 1)){
      patient[j] <- data.frame(tapply(expr.T[,j], expr.T$Sample, sum))
    }
    names(patient) <- colnames(expr.T)[-ncol(expr.T)]
    patient <- data.frame(patient)
    rownames(patient) <- names(tapply(expr.T[,1], expr.T$Sample, sum))
    Patient_sample <- as.data.frame(table(expr.T$Sample), stringAsFactors = F)
    patient$Sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- as.data.frame(patient, stringAsFactors = F)
    patient$Sample <- as.numeric(patient$Sample)
    for(k in 1:(ncol(patient) - 1)){
      for(l in 1:nrow(patient)){
        patient[l,k] <- patient[l,k]/patient$Sample[l]
      }
    }
    patient$Sample <- NULL
    colnames(patient) <- gene
    return(patient)
  }
  patient <- avgExpPerSample(Tcell,genes)
  
  if (i == 1) {
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD4_VCX3A/avgExpPerPatient.xlsx",sheetName = names(refCells)[1])
  }else{
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD4_VCX3A/avgExpPerPatient.xlsx",sheetName = names(refCells)[i],append = TRUE)
  }
  
  # singlecell
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD4_VCX3A/corr/singlecell/"
  
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  for(m in 2:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[1], y = colnames(patient)[m])) +
      geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(patient)[m],"_vs_",colnames(patient)[1],".pdf"),p2,width = 6,height = 6)
  }
  
  # TCGA
  colnames(patient) <- c("VCX3A",refCells[[i]])
  expr.LIHC <- LIHC.rnaseq[,colnames(patient)]
  expr.norm <- log2(expr.LIHC + 1)
  
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD4_VCX3A/corr/TCGA/"
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  
  for(n in 2:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x=colnames(expr.norm)[1], y=colnames(expr.norm)[n])) +
      geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(expr.norm)[n],"_vs_",colnames(expr.norm)[1],".pdf"),p2,width = 6,height = 6)
  }
}


#>>>>>>>> CD4_PLCG2

target <- cellEmbeddings(object = Tcell, cellCluster = "CD4_PLCG2")
target <- target[target$UMAP_1 > -0.6 & target$UMAP_1 < 5 & target$UMAP_2 < 0,]
# cells highlight
p_CD4_PLCG2 <- DimPlot(Tcell,group.by = "novelCluster",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
                        cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/CD3_High/CD4_PLCG2/CD4_PLCG2.highlight.pdf",p_CD4_PLCG2,device = "pdf",width = 6,height = 6)

# markers
Idents(Tcell) <- "novelCluster"
markers <- FindMarkers(Tcell, ident.1 = "CD4_PLCG2",only.pos = TRUE)

genes <- rownames(markers)
for(i in seq_along(genes)){
  p <- FeaturePlot(Tcell,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3_High/CD4_PLCG2/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}

# signature Correlation Analysis

for(i in seq_along(refCells)){
  genes <- c("PLCG2",refCells_sc[[i]])
  avgExpPerSample <- function(object, gene){
    expr <- GetAssayData(object, assay = "RNA", slot = "data")[gene,]
    expr <- as.data.frame(expr, stringAsFactors = F)
    expr.T <- data.frame(t(expr))
    expr.T$Patient <- plyr::mapvalues(rownames(expr.T), from = colnames(object), to = object$Sample)
    colnames(expr.T) <- c(gene, "Sample")
    
    patient <- NULL
    for(j in 1:(ncol(expr.T) - 1)){
      patient[j] <- data.frame(tapply(expr.T[,j], expr.T$Sample, sum))
    }
    names(patient) <- colnames(expr.T)[-ncol(expr.T)]
    patient <- data.frame(patient)
    rownames(patient) <- names(tapply(expr.T[,1], expr.T$Sample, sum))
    Patient_sample <- as.data.frame(table(expr.T$Sample), stringAsFactors = F)
    patient$Sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- as.data.frame(patient, stringAsFactors = F)
    patient$Sample <- as.numeric(patient$Sample)
    for(k in 1:(ncol(patient) - 1)){
      for(l in 1:nrow(patient)){
        patient[l,k] <- patient[l,k]/patient$Sample[l]
      }
    }
    patient$Sample <- NULL
    colnames(patient) <- gene
    return(patient)
  }
  patient <- avgExpPerSample(Tcell,genes)
  
  if (i == 1) {
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD4_PLCG2/avgExpPerPatient.xlsx",sheetName = names(refCells)[1])
  }else{
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD4_PLCG2/avgExpPerPatient.xlsx",sheetName = names(refCells)[i],append = TRUE)
  }
  
  # singlecell
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD4_PLCG2/corr/singlecell/"
  
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  for(m in 2:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[1], y = colnames(patient)[m])) +
      geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(patient)[m],"_vs_",colnames(patient)[1],".pdf"),p2,width = 6,height = 6)
  }
  
  # TCGA
  colnames(patient) <- c("PLCG2",refCells[[i]])
  expr.LIHC <- LIHC.rnaseq[,colnames(patient)]
  expr.norm <- log2(expr.LIHC + 1)
  
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD4_PLCG2/corr/TCGA/"
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  
  for(n in 2:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x=colnames(expr.norm)[1], y=colnames(expr.norm)[n])) +
      geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(expr.norm)[n],"_vs_",colnames(expr.norm)[1],".pdf"),p2,width = 6,height = 6)
  }
}


#>>>>> CD4_IL1R2

target <- cellEmbeddings(object = Tcell, cellCluster = "CD4_IL1R2")
target <- target[target$UMAP_2 < 0,]
# cells highlight
p_CD4_PLCG2 <- DimPlot(Tcell,group.by = "novelCluster",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
                       cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/CD3_High/CD4_IL1R2/CD4_IL1R2.highlight.pdf",p_CD4_PLCG2,device = "pdf",width = 6,height = 6)

# markers
Idents(Tcell) <- "novelCluster"
markers <- FindMarkers(Tcell, ident.1 = "CD4_IL1R2",only.pos = TRUE)

genes <- rownames(markers)
for(i in seq_along(genes)){
  p <- FeaturePlot(Tcell,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3_High/CD4_IL1R2/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}

# signature Correlation Analysis

for(i in seq_along(refCells)){
  genes <- c("IL1R2",refCells_sc[[i]])
  avgExpPerSample <- function(object, gene){
    expr <- GetAssayData(object, assay = "RNA", slot = "data")[gene,]
    expr <- as.data.frame(expr, stringAsFactors = F)
    expr.T <- data.frame(t(expr))
    expr.T$Patient <- plyr::mapvalues(rownames(expr.T), from = colnames(object), to = object$Sample)
    colnames(expr.T) <- c(gene, "Sample")
    
    patient <- NULL
    for(j in 1:(ncol(expr.T) - 1)){
      patient[j] <- data.frame(tapply(expr.T[,j], expr.T$Sample, sum))
    }
    names(patient) <- colnames(expr.T)[-ncol(expr.T)]
    patient <- data.frame(patient)
    rownames(patient) <- names(tapply(expr.T[,1], expr.T$Sample, sum))
    Patient_sample <- as.data.frame(table(expr.T$Sample), stringAsFactors = F)
    patient$Sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- as.data.frame(patient, stringAsFactors = F)
    patient$Sample <- as.numeric(patient$Sample)
    for(k in 1:(ncol(patient) - 1)){
      for(l in 1:nrow(patient)){
        patient[l,k] <- patient[l,k]/patient$Sample[l]
      }
    }
    patient$Sample <- NULL
    colnames(patient) <- gene
    return(patient)
  }
  patient <- avgExpPerSample(Tcell,genes)
  
  if (i == 1) {
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD4_IL1R2/avgExpPerPatient.xlsx",sheetName = names(refCells)[1])
  }else{
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD4_IL1R2/avgExpPerPatient.xlsx",sheetName = names(refCells)[i],append = TRUE)
  }
  
  # singlecell
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD4_IL1R2/corr/singlecell/"
  
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  for(m in 2:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[1], y = colnames(patient)[m])) +
      geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(patient)[m],"_vs_",colnames(patient)[1],".pdf"),p2,width = 6,height = 6)
  }
  
  # TCGA
  colnames(patient) <- c("IL1R2",refCells[[i]])
  expr.LIHC <- LIHC.rnaseq[,colnames(patient)]
  expr.norm <- log2(expr.LIHC + 1)
  
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD4_IL1R2/corr/TCGA/"
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  
  for(n in 2:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x=colnames(expr.norm)[1], y=colnames(expr.norm)[n])) +
      geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(expr.norm)[n],"_vs_",colnames(expr.norm)[1],".pdf"),p2,width = 6,height = 6)
  }
}


#>>>>>> CD8_FCGR3A

target <- cellEmbeddings(object = Tcell, cellCluster = "CD8_FCGR3A")
target <- target[target$UMAP_1 > 0  & target$UMAP_2 > 0,]
# cells highlight
p_CD4_PLCG2 <- DimPlot(Tcell,group.by = "novelCluster",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
                       cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/CD3_High/CD8_FCGR3A/CD8_FCGR3A.highlight.pdf",p_CD4_PLCG2,device = "pdf",width = 6,height = 6)

# markers
Idents(Tcell) <- "novelCluster"
markers <- FindMarkers(Tcell, ident.1 = "CD8_FCGR3A",only.pos = TRUE)

genes <- rownames(markers)
for(i in seq_along(genes)){
  p <- FeaturePlot(Tcell,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3_High/CD8_FCGR3A/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}

# signature Correlation Analysis

for(i in seq_along(refCells_CD8)){
  genes <- c("FCGR3A",refCells_CD8[[i]])
  avgExpPerSample <- function(object, gene){
    expr <- GetAssayData(object, assay = "RNA", slot = "data")[gene,]
    expr <- as.data.frame(expr, stringAsFactors = F)
    expr.T <- data.frame(t(expr))
    expr.T$Patient <- plyr::mapvalues(rownames(expr.T), from = colnames(object), to = object$Sample)
    colnames(expr.T) <- c(gene, "Sample")
    
    patient <- NULL
    for(j in 1:(ncol(expr.T) - 1)){
      patient[j] <- data.frame(tapply(expr.T[,j], expr.T$Sample, sum))
    }
    names(patient) <- colnames(expr.T)[-ncol(expr.T)]
    patient <- data.frame(patient)
    rownames(patient) <- names(tapply(expr.T[,1], expr.T$Sample, sum))
    Patient_sample <- as.data.frame(table(expr.T$Sample), stringAsFactors = F)
    patient$Sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- as.data.frame(patient, stringAsFactors = F)
    patient$Sample <- as.numeric(patient$Sample)
    for(k in 1:(ncol(patient) - 1)){
      for(l in 1:nrow(patient)){
        patient[l,k] <- patient[l,k]/patient$Sample[l]
      }
    }
    patient$Sample <- NULL
    colnames(patient) <- gene
    return(patient)
  }
  patient <- avgExpPerSample(Tcell,genes)
  
  if (i == 1) {
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD8_FCGR3A/avgExpPerPatient.xlsx",sheetName = names(refCells_CD8)[1])
  }else{
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD8_FCGR3A/avgExpPerPatient.xlsx",sheetName = names(refCells_CD8)[i],append = TRUE)
  }
  
  # singlecell
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD8_FCGR3A/corr/singlecell/"
  
  if(!file.exists(paste0(pwd,names(refCells_CD8[i])))){
    dir.create(paste0(pwd,names(refCells_CD8)[i]))
  }
  for(m in 2:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[1], y = colnames(patient)[m])) +
      geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells_CD8)[i],"/",colnames(patient)[m],"_vs_",colnames(patient)[1],".pdf"),p2,width = 6,height = 6)
  }
  
  # TCGA
  colnames(patient) <- c("FCGR3A",refCells_CD8[[i]])
  expr.LIHC <- LIHC.rnaseq[,colnames(patient)]
  expr.norm <- log2(expr.LIHC + 1)
  
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD8_FCGR3A/corr/TCGA/"
  if(!file.exists(paste0(pwd,names(refCells_CD8[i])))){
    dir.create(paste0(pwd,names(refCells_CD8)[i]))
  }
  
  for(n in 2:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x=colnames(expr.norm)[1], y=colnames(expr.norm)[n])) +
      geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells_CD8)[i],"/",colnames(expr.norm)[n],"_vs_",colnames(expr.norm)[1],".pdf"),p2,width = 6,height = 6)
  }
}


#>>>>>>> CD4_FGR

target <- cellEmbeddings(object = Tcell, cellCluster = "CD4_FGR")
# cells highlight
p_CD4_FGR <- DimPlot(Tcell,group.by = "novelCluster",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
                       cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/CD3_High/CD4_FGR/CD4_FGR.highlight.pdf",p_CD4_FGR,device = "pdf",width = 6,height = 6)

# markers
Idents(Tcell) <- "novelCluster"
markers <- FindMarkers(Tcell, ident.1 = "CD4_FGR",only.pos = TRUE)

genes <- rownames(markers)
for(i in seq_along(genes)){
  p <- FeaturePlot(Tcell,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3_High/CD4_FGR/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}

# signature Correlation Analysis

for(i in seq_along(refCells)){
  genes <- c("FGR",refCells_sc[[i]])
  avgExpPerSample <- function(object, gene){
    expr <- GetAssayData(object, assay = "RNA", slot = "data")[gene,]
    expr <- as.data.frame(expr, stringAsFactors = F)
    expr.T <- data.frame(t(expr))
    expr.T$Patient <- plyr::mapvalues(rownames(expr.T), from = colnames(object), to = object$Sample)
    colnames(expr.T) <- c(gene, "Sample")
    
    patient <- NULL
    for(j in 1:(ncol(expr.T) - 1)){
      patient[j] <- data.frame(tapply(expr.T[,j], expr.T$Sample, sum))
    }
    names(patient) <- colnames(expr.T)[-ncol(expr.T)]
    patient <- data.frame(patient)
    rownames(patient) <- names(tapply(expr.T[,1], expr.T$Sample, sum))
    Patient_sample <- as.data.frame(table(expr.T$Sample), stringAsFactors = F)
    patient$Sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- as.data.frame(patient, stringAsFactors = F)
    patient$Sample <- as.numeric(patient$Sample)
    for(k in 1:(ncol(patient) - 1)){
      for(l in 1:nrow(patient)){
        patient[l,k] <- patient[l,k]/patient$Sample[l]
      }
    }
    patient$Sample <- NULL
    colnames(patient) <- gene
    return(patient)
  }
  patient <- avgExpPerSample(Tcell,genes)
  
  if (i == 1) {
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD4_FGR/avgExpPerPatient.xlsx",sheetName = names(refCells)[1])
  }else{
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_High/CD4_FGR/avgExpPerPatient.xlsx",sheetName = names(refCells)[i],append = TRUE)
  }
  
  # singlecell
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD4_FGR/corr/singlecell/"
  
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  for(m in 2:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[1], y = colnames(patient)[m])) +
      geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(patient)[m],"_vs_",colnames(patient)[1],".pdf"),p2,width = 6,height = 6)
  }
  
  # TCGA
  colnames(patient) <- c("FGR",refCells[[i]])
  expr.LIHC <- LIHC.rnaseq[,colnames(patient)]
  expr.norm <- log2(expr.LIHC + 1)
  
  pwd = "result/without_PVTT_MLN/tumor/CD3_High/CD4_FGR/corr/TCGA/"
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  
  for(n in 2:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x=colnames(expr.norm)[1], y=colnames(expr.norm)[n])) +
      geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(expr.norm)[n],"_vs_",colnames(expr.norm)[1],".pdf"),p2,width = 6,height = 6)
  }
}

# >>>> CD4_CRP

target <- cellEmbeddings(object = Tcell, cellCluster = "CD4_CRP")
target <- target[target$UMAP_2 > -3.5,]
# cells highlight
p_CD4_CRP <- DimPlot(Tcell,group.by = "novelCluster",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
                      cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())

ggsave("result/without_PVTT_MLN/tumor/PARTV.1/CD3_High/CD4_CRP/CD4_CRP.highlight.pdf",p_CD4_CRP,device = "pdf",width = 6,height = 6)

# FIND MARKERS
Idents(Tcell) <- Tcell$novelCluster
markers <- FindMarkers(Tcell,only.pos = TRUE,ident.1 = "CD4_CRP")
genes <- rownames(markers)
for(i in seq_along(genes)){
  p <- FeaturePlot(Tcell,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/PARTV.1/CD3_High/CD4_CRP/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}


# signature Correlation Analysis

load("refdb/TCGA_LIHC.rnaseq.RData")

for(i in seq_along(refCells)){
  genes <- c("CRP",refCells_sc[[i]])
  avgExpPerSample <- function(object, gene){
    expr <- GetAssayData(object, assay = "RNA", slot = "data")[gene,]
    expr <- as.data.frame(expr, stringAsFactors = F)
    expr.T <- data.frame(t(expr))
    expr.T$Patient <- plyr::mapvalues(rownames(expr.T), from = colnames(object), to = object$Sample)
    colnames(expr.T) <- c(gene, "Sample")
    
    patient <- NULL
    for(j in 1:(ncol(expr.T) - 1)){
      patient[j] <- data.frame(tapply(expr.T[,j], expr.T$Sample, sum))
    }
    names(patient) <- colnames(expr.T)[-ncol(expr.T)]
    patient <- data.frame(patient)
    rownames(patient) <- names(tapply(expr.T[,1], expr.T$Sample, sum))
    Patient_sample <- as.data.frame(table(expr.T$Sample), stringAsFactors = F)
    patient$Sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- as.data.frame(patient, stringAsFactors = F)
    patient$Sample <- as.numeric(patient$Sample)
    for(k in 1:(ncol(patient) - 1)){
      for(l in 1:nrow(patient)){
        patient[l,k] <- patient[l,k]/patient$Sample[l]
      }
    }
    patient$Sample <- NULL
    colnames(patient) <- gene
    return(patient)
  }
  patient <- avgExpPerSample(Tcell,genes)
  
  if (i == 1) {
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/PARTV.1/CD3_High/CD4_CRP/avgExpPerPatient.xlsx",sheetName = names(refCells)[1])
  }else{
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/PARTV.1/CD3_High/CD4_CRP/avgExpPerPatient.xlsx",sheetName = names(refCells)[i],append = TRUE)
  }
  
  # singlecell
  pwd = "result/without_PVTT_MLN/tumor/PARTV.1/CD3_High/CD4_CRP/corr/singlecell/"
  
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  for(m in 2:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[1], y = colnames(patient)[m])) +
      geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(patient)[m],"_vs_",colnames(patient)[1],".pdf"),p2,width = 6,height = 6)
  }
  
  # TCGA
  colnames(patient) <- c("CRP",refCells[[i]])
  expr.LIHC <- LIHC.rnaseq[,colnames(patient)]
  expr.norm <- log2(expr.LIHC + 1)
  
  pwd = "result/without_PVTT_MLN/tumor/PARTV.1/CD3_High/CD4_CRP/corr/TCGA/"
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  
  for(n in 2:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x=colnames(expr.norm)[1], y=colnames(expr.norm)[n])) +
      geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(expr.norm)[n],"_vs_",colnames(expr.norm)[1],".pdf"),p2,width = 6,height = 6)
  }
}




## CD3+ low group: 

# CD4_CHGA, CD4_AKR1C4

#>>>>> CD4_CHGA

target <- cellEmbeddings(object = Tcell, cellCluster = "CD4_CHGA")
target <- target[target$UMAP_1 < 0 & target$UMAP_1 > -8 & target$UMAP_2 < 2 & target$UMAP_2 > -3,]
# cells highlight
p_CD4_CHGA <- DimPlot(Tcell,group.by = "novelCluster",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
                     cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/CD3_Low/CD4_CHGA/CD4_CHGA.highlight.pdf",p_CD4_CHGA,device = "pdf",width = 6,height = 6)

# markers
Idents(Tcell) <- "novelCluster"
markers <- FindMarkers(Tcell, ident.1 = "CD4_CHGA",only.pos = TRUE)

genes <- rownames(markers)
for(i in seq_along(genes)){
  p <- FeaturePlot(Tcell,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3_Low/CD4_CHGA/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}

# signature Correlation Analysis

for(i in seq_along(refCells)){
  genes <- c("CHGA",refCells_sc[[i]])
  avgExpPerSample <- function(object, gene){
    expr <- GetAssayData(object, assay = "RNA", slot = "data")[gene,]
    expr <- as.data.frame(expr, stringAsFactors = F)
    expr.T <- data.frame(t(expr))
    expr.T$Patient <- plyr::mapvalues(rownames(expr.T), from = colnames(object), to = object$Sample)
    colnames(expr.T) <- c(gene, "Sample")
    
    patient <- NULL
    for(j in 1:(ncol(expr.T) - 1)){
      patient[j] <- data.frame(tapply(expr.T[,j], expr.T$Sample, sum))
    }
    names(patient) <- colnames(expr.T)[-ncol(expr.T)]
    patient <- data.frame(patient)
    rownames(patient) <- names(tapply(expr.T[,1], expr.T$Sample, sum))
    Patient_sample <- as.data.frame(table(expr.T$Sample), stringAsFactors = F)
    patient$Sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- as.data.frame(patient, stringAsFactors = F)
    patient$Sample <- as.numeric(patient$Sample)
    for(k in 1:(ncol(patient) - 1)){
      for(l in 1:nrow(patient)){
        patient[l,k] <- patient[l,k]/patient$Sample[l]
      }
    }
    patient$Sample <- NULL
    colnames(patient) <- gene
    return(patient)
  }
  patient <- avgExpPerSample(Tcell,genes)
  
  if (i == 1) {
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_Low/CD4_CHGA/avgExpPerPatient.xlsx",sheetName = names(refCells)[1])
  }else{
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_Low/CD4_CHGA/avgExpPerPatient.xlsx",sheetName = names(refCells)[i],append = TRUE)
  }
  
  # singlecell
  pwd = "result/without_PVTT_MLN/tumor/CD3_Low/CD4_CHGA/corr/singlecell/"
  
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  for(m in 2:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[1], y = colnames(patient)[m])) +
      geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(patient)[m],"_vs_",colnames(patient)[1],".pdf"),p2,width = 6,height = 6)
  }
  
  # TCGA
  colnames(patient) <- c("CHGA",refCells[[i]])
  expr.LIHC <- LIHC.rnaseq[,colnames(patient)]
  expr.norm <- log2(expr.LIHC + 1)
  
  pwd = "result/without_PVTT_MLN/tumor/CD3_Low/CD4_CHGA/corr/TCGA/"
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  
  for(n in 2:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x=colnames(expr.norm)[1], y=colnames(expr.norm)[n])) +
      geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(expr.norm)[n],"_vs_",colnames(expr.norm)[1],".pdf"),p2,width = 6,height = 6)
  }
}


#>>>>> CD4_AKR1C4

target <- cellEmbeddings(object = Tcell, cellCluster = "CD4_AKR1C4")
target <- target[target$UMAP_1 < 0 & target$UMAP_1 > -8 & target$UMAP_2 < 2 & target$UMAP_2 > -3,]
# cells highlight
p_CD4_AKR1C4 <- DimPlot(Tcell,group.by = "novelCluster",pt.size = 0.9,label = TRUE,label.box = TRUE,label.color = "white",label.size = 5,
                      cells.highlight = rownames(target)) + NoLegend() +
  theme(panel.border = element_rect(colour = "black",size = 1),
        plot.title = element_blank())
ggsave("result/without_PVTT_MLN/tumor/CD3_Low/CD4_AKR1C4/CD4_AKR1C4.highlight.pdf",p_CD4_AKR1C4,device = "pdf",width = 6,height = 6)

# markers
Idents(Tcell) <- "novelCluster"
markers <- FindMarkers(Tcell, ident.1 = "CD4_AKR1C4",only.pos = TRUE)

genes <- rownames(markers)
for(i in seq_along(genes)){
  p <- FeaturePlot(Tcell,features = genes[i],pt.size = 0.98,cols = c("lightgrey","#008956")) + NoLegend() +
    theme(panel.border = element_rect(size = 1,colour = "grey"),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3_Low/CD4_AKR1C4/markers/",i,"_",genes[i],".pdf"),p,width = 3,height = 3)
}

# signature Correlation Analysis

for(i in seq_along(refCells)){
  genes <- c("AKR1C4",refCells_sc[[i]])
  avgExpPerSample <- function(object, gene){
    expr <- GetAssayData(object, assay = "RNA", slot = "data")[gene,]
    expr <- as.data.frame(expr, stringAsFactors = F)
    expr.T <- data.frame(t(expr))
    expr.T$Patient <- plyr::mapvalues(rownames(expr.T), from = colnames(object), to = object$Sample)
    colnames(expr.T) <- c(gene, "Sample")
    
    patient <- NULL
    for(j in 1:(ncol(expr.T) - 1)){
      patient[j] <- data.frame(tapply(expr.T[,j], expr.T$Sample, sum))
    }
    names(patient) <- colnames(expr.T)[-ncol(expr.T)]
    patient <- data.frame(patient)
    rownames(patient) <- names(tapply(expr.T[,1], expr.T$Sample, sum))
    Patient_sample <- as.data.frame(table(expr.T$Sample), stringAsFactors = F)
    patient$Sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- as.data.frame(patient, stringAsFactors = F)
    patient$Sample <- as.numeric(patient$Sample)
    for(k in 1:(ncol(patient) - 1)){
      for(l in 1:nrow(patient)){
        patient[l,k] <- patient[l,k]/patient$Sample[l]
      }
    }
    patient$Sample <- NULL
    colnames(patient) <- gene
    return(patient)
  }
  patient <- avgExpPerSample(Tcell,genes)
  
  if (i == 1) {
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_Low/CD4_AKR1C4/avgExpPerPatient.xlsx",sheetName = names(refCells)[1])
  }else{
    xlsx::write.xlsx(patient,file = "result/without_PVTT_MLN/tumor/CD3_Low/CD4_AKR1C4/avgExpPerPatient.xlsx",sheetName = names(refCells)[i],append = TRUE)
  }
  
  # singlecell
  pwd = "result/without_PVTT_MLN/tumor/CD3_Low/CD4_AKR1C4/corr/singlecell/"
  
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  for(m in 2:ncol(patient)){
    p1 <- ggplot(patient, aes_string(x = colnames(patient)[1], y = colnames(patient)[m])) +
      geom_point(size = 1.5, color = "#F9B208", alpha = .7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25, 'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(patient)[m],"_vs_",colnames(patient)[1],".pdf"),p2,width = 6,height = 6)
  }
  
  # TCGA
  colnames(patient) <- c("AKR1C4",refCells[[i]])
  expr.LIHC <- LIHC.rnaseq[,colnames(patient)]
  expr.norm <- log2(expr.LIHC + 1)
  
  pwd = "result/without_PVTT_MLN/tumor/CD3_Low/CD4_AKR1C4/corr/TCGA/"
  if(!file.exists(paste0(pwd,names(refCells[i])))){
    dir.create(paste0(pwd,names(refCells)[i]))
  }
  
  for(n in 2:ncol(expr.norm)){
    p1 <- ggplot(expr.norm, aes_string(x=colnames(expr.norm)[1], y=colnames(expr.norm)[n])) +
      geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
      ) +
      geom_smooth(method = "lm",se = T, color = '#BF0032', size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3, size = 6)
    p2 <- ggExtra::ggMarginal(p1,type = "density",
                              xparams = list(binwidth = 0.1, fill = '#B3E283',size=.7),
                              yparams = list(binwidth = 0.1, fill = '#8AB6D6',size=.7))
    ggsave(paste0(pwd,names(refCells)[i],"/",colnames(expr.norm)[n],"_vs_",colnames(expr.norm)[1],".pdf"),p2,width = 6,height = 6)
  }
}
