
##================== 1.Dealing with Overlapped Genes =====================

## 1.1 Dealing with Overlapped DE genes between CD4+T & CD8+T
###================================================
tumor_sub <- readRDS("rds/without_PVTT_MLN/tumor/CD45_removal_T.IR.2.rds")

## >>>>>>>>>>>>>>>>>>>>>> CAF  <<<<<<<<<<<<<<<<<<<<<<<<<<<

## CAF.UP (CD 4+T high_vs_CD8+T high)
##============================================================
wkd = "result/without_PVTT_MLN/tumor/overlap/CD4T_vs_CD8T/without_CD45/UP/CAF/"

CAF.UP.genes <- read.table(paste0(wkd,"CAF.common.28.txt"),header = F,stringsAsFactors = F)
CAF.UP.genes$V1 <- gsub("\\.","-",CAF.UP.genes$V1)
# obtain expression matrix
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF.UP.genes$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient

if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}


if(F){
  patient <- NULL
  for (i in 1:nrow(CAF.UP.genes)) {
    patient[[i]] <- AverageExpression(tumor_sub,slot = "data",features = CAF.UP.genes$V1[i],assays = "RNA")
  }
}

# final expression matrix
CAF.UP.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF.UP.exp,file = paste0(wkd,"CAF.exp.xlsx"),sheetName = "normalized")

normEXP <- CAF.UP.exp[rowSums(CAF.UP.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 2) scaleEXP[scaleEXP > 2] = 2
if(min(scaleEXP) < -2) scaleEXP[scaleEXP < -2] = -2

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                   filename = paste0(wkd,ncol(scaleEXP),".CAF.genes.heatmap.pdf"))

### barplot
library(ggpubr)
CAF.UP.exp.dt <- data.frame(t(CAF.UP.exp))
CAF.UP.exp.dt$Patient <- rownames(CAF.UP.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF.UP.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF.UP.exp.dt,x = "Patient",y=colnames(CAF.UP.exp.dt)[i],
                        fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF.UP.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"barplot/",colnames(CAF.UP.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###==============================================================


##>>>> CAF.Down (CD4+T low_vs_CD8+T low)
###==============================================================
wkd = "result/without_PVTT_MLN/tumor/overlap/CD4T_vs_CD8T/without_CD45/Down/CAF/"

CAF.Down.genes <- read.table(paste0(wkd,"CAF.common.147.txt"),header = F,stringsAsFactors = F)
# obtain expression matrix
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF.Down.genes$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

if(F){
  # calculate mean expression values of Overlapped DE genes in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}

# final expression matrix
CAF.Down.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF.Down.exp,file = paste0(wkd,"CAF.exp.xlsx"),sheetName = "normalized")

normEXP <- CAF.Down.exp[rowSums(CAF.Down.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                   filename = paste0(wkd,ncol(scaleEXP),".CAF.genes.heatmap.pdf"))

### barplot
library(ggpubr)
CAF.Down.exp.dt <- data.frame(t(CAF.Down.exp))
CAF.Down.exp.dt$Patient <- rownames(CAF.Down.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF.Down.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF.Down.exp.dt,x = "Patient",y=colnames(CAF.Down.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF.Down.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"barplot/",colnames(CAF.Down.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###==============================================================


## >>>>>>>>>>>>>>>>>>>>>> TEC  <<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>> TEC.UP (CD 4+T high_vs_CD8+T high)
###==============================================================
wkd = "result/without_PVTT_MLN/tumor/overlap/CD4T_vs_CD8T/without_CD45/UP/TEC/"
TEC.UP.genes <- read.table(paste0(wkd,"TEC.common.23.txt"),header = F,stringsAsFactors = F)
# obtain expression matrix
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[TEC.UP.genes$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient

if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}


if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}

# final expression matrix
TEC.UP.exp <- data.frame(t(patient))
xlsx::write.xlsx(TEC.UP.exp,file = paste0(wkd,"TEC.exp.xlsx"),sheetName = "normalized")

normEXP <- TEC.UP.exp[rowSums(TEC.UP.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                   filename = paste0(wkd,ncol(scaleEXP),".TEC.genes.heatmap.pdf"))

### barplot
library(ggpubr)
TEC.UP.exp.dt <- data.frame(t(TEC.UP.exp))
TEC.UP.exp.dt$Patient <- rownames(TEC.UP.exp.dt)
plt <- list()
for(i in 1:(ncol(TEC.UP.exp.dt)-1)){
  plt[[i]] <- ggbarplot(TEC.UP.exp.dt,x = "Patient",y=colnames(TEC.UP.exp.dt)[i],
                        fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(TEC.UP.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"barplot/",colnames(TEC.UP.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###==============================================================

##>>>> TEC.Down (CD4+T low_vs_CD8+T low)
###==============================================================
wkd = "result/without_PVTT_MLN/tumor/overlap/CD4T_vs_CD8T/without_CD45/Down/TEC/"
TEC.Down.genes <- read.table(paste0(wkd,"TEC.common.266.txt"),header = F,stringsAsFactors = F)
# obtain expression matrix
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[TEC.Down.genes$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}

# final expression matrix
TEC.Down.exp <- data.frame(t(patient))
xlsx::write.xlsx(TEC.Down.exp,file = paste0(wkd,"TEC.exp.xlsx"),sheetName = "normalized")

normEXP <- TEC.Down.exp[rowSums(TEC.Down.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                   filename = paste0(wkd,ncol(scaleEXP),".TEC.genes.heatmap.pdf"))

### barplot
library(ggpubr)
TEC.Down.exp.dt <- data.frame(t(TEC.Down.exp))
TEC.Down.exp.dt$Patient <- rownames(TEC.Down.exp.dt)
plt <- list()
for(i in 1:(ncol(TEC.Down.exp.dt)-1)){
  plt[[i]] <- ggbarplot(TEC.Down.exp.dt,x = "Patient",y=colnames(TEC.Down.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(TEC.Down.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"barplot/",colnames(TEC.Down.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###==============================================================


## >>>>>>>>>>>>>>>>>>>>>> Malignant Cell  <<<<<<<<<<<<<<<<<<<<<<<<<<<

##>>>> Malignant.Cell.UP (CD 4+T high_vs_CD8+T high)
###==============================================================
wkd = "result/without_PVTT_MLN/tumor/overlap/CD4T_vs_CD8T/without_CD45/UP/Malignant.Cell/"
MLG.UP.genes <- read.table(paste0(wkd,"Malignant.Cell.common.32.txt"),header = F,stringsAsFactors = F)
# obtain expression matrix
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[MLG.UP.genes$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient

if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}


# final expression matrix
MLG.UP.exp <- data.frame(t(patient))
xlsx::write.xlsx(MLG.UP.exp,file = paste0(wkd,"Malignant.Cell.exp.xlsx"),sheetName = "normalized")

normEXP <- MLG.UP.exp[rowSums(MLG.UP.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                   filename = paste0(wkd,ncol(scaleEXP),".Malignant.Cell.genes.heatmap.pdf"))

### barplot
library(ggpubr)
MLG.UP.exp.dt <- data.frame(t(MLG.UP.exp))
MLG.UP.exp.dt$Patient <- rownames(MLG.UP.exp.dt)
plt <- list()
for(i in 1:(ncol(MLG.UP.exp.dt)-1)){
  plt[[i]] <- ggbarplot(MLG.UP.exp.dt,x = "Patient",y=colnames(MLG.UP.exp.dt)[i],
                        fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(MLG.UP.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"barplot/",colnames(MLG.UP.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###==============================================================

##>>>> Malignant Cell.Down (CD4+T low_vs_CD8+T low)
###==============================================================
wkd = "result/without_PVTT_MLN/tumor/overlap/CD4T_vs_CD8T/without_CD45/Down/Malignant.Cell/"
MLG.Down.genes <- read.table(paste0(wkd,"Malignant.cell.common.167.txt"),header = F,stringsAsFactors = F)
MLG.Down.genes$V1 <- sub("\\.", "-", MLG.Down.genes$V1)
# obtain expression matrix
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[MLG.Down.genes$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}

# final expression matrix
MLG.Down.exp <- data.frame(t(patient))
xlsx::write.xlsx(MLG.Down.exp,file = paste0(wkd,"Malignant.Cell.exp.xlsx"),sheetName = "normalized")

normEXP <- MLG.Down.exp[rowSums(MLG.Down.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                   filename = paste0(wkd,ncol(scaleEXP),".Malignant.Cell.genes.heatmap.pdf"))

### barplot
library(ggpubr)
MLG.Down.exp.dt <- data.frame(t(MLG.Down.exp))
MLG.Down.exp.dt$Patient <- rownames(MLG.Down.exp.dt)
plt <- list()
for(i in 1:(ncol(MLG.Down.exp.dt)-1)){
  plt[[i]] <- ggbarplot(MLG.Down.exp.dt,x = "Patient",y=colnames(MLG.Down.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(MLG.Down.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"barplot/",colnames(MLG.Down.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###=====================================================================================================


##>>>>>>>> 1.2 Dealing with Overlapped DE genes in each T cells <<<<<<<<<<
###================================================
tumor_sub <- readRDS("rds/without_PVTT_MLN/tumor/CD45_removal_T.IR.2.rds")

##>>>> CD3+T UP (CAF, TEC, Malignant cell)
###==============================================================
wkd = "result/without_PVTT_MLN/tumor/overlap/CD3T/UP/"

##===== 1.CAF TEC MLG ========
CAF_TEC_MLG <- read.table(paste0(wkd,"CAF_TEC_Malignant.cell/","CAF_TEC_Malignant.cell.common.8.txt"),header = F,stringsAsFactors = F)

expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_TEC_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}

# final expression matrix
CAF_TEC_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_TEC_MLG.exp,file = paste0(wkd,"CAF_TEC_Malignant.cell/","norm.exp.xlsx"))

normEXP <- CAF_TEC_MLG.exp[rowSums(CAF_TEC_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 6,height = 6,
                   filename = paste0(wkd,"CAF_TEC_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_TEC_MLG.exp.dt <- data.frame(t(CAF_TEC_MLG.exp))
CAF_TEC_MLG.exp.dt$Patient <- rownames(CAF_TEC_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_TEC_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_TEC_MLG.exp.dt,x = "Patient",y=colnames(CAF_TEC_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_TEC_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_TEC_Malignant.cell/","barplot/",colnames(CAF_TEC_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 2.CAF TEC ========
CAF_TEC <- read.table(paste0(wkd,"CAF_TEC/","CAF_TEC.common.42.txt"),header = F,stringsAsFactors = F)

expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_TEC$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_TEC.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_TEC.exp,file = paste0(wkd,"CAF_TEC/","norm.exp.xlsx"))

normEXP <- CAF_TEC.exp[rowSums(CAF_TEC.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                   filename = paste0(wkd,"CAF_TEC/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_TEC.exp.dt <- data.frame(t(CAF_TEC.exp))
CAF_TEC.exp.dt$Patient <- rownames(CAF_TEC.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_TEC.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_TEC.exp.dt,x = "Patient",y=colnames(CAF_TEC.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_TEC.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_TEC/","barplot/",colnames(CAF_TEC.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 3.CAF Malignant cell ========
CAF_MLG <- read.table(paste0(wkd,"CAF_Malignant.cell/","CAF_Malignant cell.common.10.txt"),header = F,stringsAsFactors = F)

expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}

# final expression matrix
CAF_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_MLG.exp,file = paste0(wkd,"CAF_Malignant.cell/","norm.exp.xlsx"))

normEXP <- CAF_MLG.exp[rowSums(CAF_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 6,height = 6,
                   filename = paste0(wkd,"CAF_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_MLG.exp.dt <- data.frame(t(CAF_MLG.exp))
CAF_MLG.exp.dt$Patient <- rownames(CAF_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_MLG.exp.dt,x = "Patient",y=colnames(CAF_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_Malignant.cell/","barplot/",colnames(CAF_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 4.TEC Malignant cell ========
TEC_MLG <- read.table(paste0(wkd,"TEC_Malignant.cell/","TEC_Malignant cell.common.11.txt"),header = F,stringsAsFactors = F)

expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[TEC_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient

if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
TEC_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(TEC_MLG.exp,file = paste0(wkd,"TEC_Malignant.cell/","norm.exp.xlsx"))

normEXP <- TEC_MLG.exp[rowSums(TEC_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 6,height = 6,
                   filename = paste0(wkd,"TEC_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
TEC_MLG.exp.dt <- data.frame(t(TEC_MLG.exp))
TEC_MLG.exp.dt$Patient <- rownames(TEC_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(TEC_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(TEC_MLG.exp.dt,x = "Patient",y=colnames(TEC_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(TEC_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"TEC_Malignant.cell/","barplot/",colnames(TEC_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================

##>>>> CD3+T Down (CAF, TEC, Malignant cell)
###==============================================================
wkd = "result/without_PVTT_MLN/tumor/overlap/CD3T/Down/"

##===== 1.CAF TEC MLG ========
CAF_TEC_MLG <- read.table(paste0(wkd,"CAF_TEC_Malignant.cell/","CAF_TEC_Malignant cell.common.59.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_TEC_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_TEC_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_TEC_MLG.exp,file = paste0(wkd,"CAF_TEC_Malignant.cell/","norm.exp.xlsx"))

normEXP <- CAF_TEC_MLG.exp[rowSums(CAF_TEC_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                   filename = paste0(wkd,"CAF_TEC_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_TEC_MLG.exp.dt <- data.frame(t(CAF_TEC_MLG.exp))
CAF_TEC_MLG.exp.dt$Patient <- rownames(CAF_TEC_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_TEC_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_TEC_MLG.exp.dt,x = "Patient",y=colnames(CAF_TEC_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_TEC_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_TEC_Malignant.cell/","barplot/",colnames(CAF_TEC_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 2.CAF TEC ========
CAF_TEC <- read.table(paste0(wkd,"CAF_TEC/","CAF_TEC.common.69.txt"),header = F,stringsAsFactors = F)

expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_TEC$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_TEC.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_TEC.exp,file = paste0(wkd,"CAF_TEC/","norm.exp.xlsx"))

normEXP <- CAF_TEC.exp[rowSums(CAF_TEC.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 22,height = 6,
                   filename = paste0(wkd,"CAF_TEC/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_TEC.exp.dt <- data.frame(t(CAF_TEC.exp))
CAF_TEC.exp.dt$Patient <- rownames(CAF_TEC.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_TEC.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_TEC.exp.dt,x = "Patient",y=colnames(CAF_TEC.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_TEC.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_TEC/","barplot/",colnames(CAF_TEC.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 3.CAF Malignant cell ========
CAF_MLG <- read.table(paste0(wkd,"CAF_Malignant.cell/","CAF_Malignant cell.common.62.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}


if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_MLG.exp,file = paste0(wkd,"CAF_Malignant.cell/","norm.exp.xlsx"))

normEXP <- CAF_MLG.exp[rowSums(CAF_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                   filename = paste0(wkd,"CAF_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_MLG.exp.dt <- data.frame(t(CAF_MLG.exp))
CAF_MLG.exp.dt$Patient <- rownames(CAF_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_MLG.exp.dt,x = "Patient",y=colnames(CAF_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_Malignant.cell/","barplot/",colnames(CAF_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 4.TEC Malignant cell ========
TEC_MLG <- read.table(paste0(wkd,"TEC_Malignant.cell/","TEC_Malignant cell.common.11.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[TEC_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}


if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
TEC_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(TEC_MLG.exp,file = paste0(wkd,"TEC_Malignant.cell/","norm.exp.xlsx"))

normEXP <- TEC_MLG.exp[rowSums(TEC_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 6,height = 6,
                   filename = paste0(wkd,"TEC_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
TEC_MLG.exp.dt <- data.frame(t(TEC_MLG.exp))
TEC_MLG.exp.dt$Patient <- rownames(TEC_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(TEC_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(TEC_MLG.exp.dt,x = "Patient",y=colnames(TEC_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(TEC_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"TEC_Malignant.cell/","barplot/",colnames(TEC_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================

##>>>> CD4+T UP (CAF, TEC, Malignant cell)
###==============================================================
wkd = "result/without_PVTT_MLN/tumor/overlap/CD4T/UP/"

##===== 1.CAF TEC MLG ========
CAF_TEC_MLG <- read.table(paste0(wkd,"CAF_TEC_Malignant.cell/","CAF_TEC_Malignant cell.common.8.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_TEC_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}


if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_TEC_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_TEC_MLG.exp,file = paste0(wkd,"CAF_TEC_Malignant.cell/","norm.exp.xlsx"))

normEXP <- CAF_TEC_MLG.exp[rowSums(CAF_TEC_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 6,height = 6,
                   filename = paste0(wkd,"CAF_TEC_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_TEC_MLG.exp.dt <- data.frame(t(CAF_TEC_MLG.exp))
CAF_TEC_MLG.exp.dt$Patient <- rownames(CAF_TEC_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_TEC_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_TEC_MLG.exp.dt,x = "Patient",y=colnames(CAF_TEC_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_TEC_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_TEC_Malignant.cell/","barplot/",colnames(CAF_TEC_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 2.CAF TEC ========
CAF_TEC <- read.table(paste0(wkd,"CAF_TEC/","CAF_TEC.common.42.txt"),header = F,stringsAsFactors = F)

expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_TEC$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}


if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_TEC.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_TEC.exp,file = paste0(wkd,"CAF_TEC/","norm.exp.xlsx"))

normEXP <- CAF_TEC.exp[rowSums(CAF_TEC.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 18,height = 6,
                   filename = paste0(wkd,"CAF_TEC/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_TEC.exp.dt <- data.frame(t(CAF_TEC.exp))
CAF_TEC.exp.dt$Patient <- rownames(CAF_TEC.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_TEC.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_TEC.exp.dt,x = "Patient",y=colnames(CAF_TEC.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_TEC.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_TEC/","barplot/",colnames(CAF_TEC.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 3.CAF Malignant cell ========
CAF_MLG <- read.table(paste0(wkd,"CAF_Malignant.cell/","CAF_Malignant cell.common.10.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}


if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_MLG.exp,file = paste0(wkd,"CAF_Malignant.cell/","norm.exp.xlsx"))

normEXP <- CAF_MLG.exp[rowSums(CAF_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 6,height = 6,
                   filename = paste0(wkd,"CAF_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_MLG.exp.dt <- data.frame(t(CAF_MLG.exp))
CAF_MLG.exp.dt$Patient <- rownames(CAF_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_MLG.exp.dt,x = "Patient",y=colnames(CAF_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_Malignant.cell/","barplot/",colnames(CAF_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 4.TEC Malignant cell ========
TEC_MLG <- read.table(paste0(wkd,"TEC_Malignant.cell/","TEC_Malignant cell.common.11.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[TEC_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
TEC_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(TEC_MLG.exp,file = paste0(wkd,"TEC_Malignant.cell/","norm.exp.xlsx"))

normEXP <- TEC_MLG.exp[rowSums(TEC_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 6,height = 6,
                   filename = paste0(wkd,"TEC_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
TEC_MLG.exp.dt <- data.frame(t(TEC_MLG.exp))
TEC_MLG.exp.dt$Patient <- rownames(TEC_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(TEC_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(TEC_MLG.exp.dt,x = "Patient",y=colnames(TEC_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(TEC_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"TEC_Malignant.cell/","barplot/",colnames(TEC_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================


##>>>> CD4+T Down (CAF, TEC, Malignant cell)
###==============================================================
wkd = "result/without_PVTT_MLN/tumor/overlap/CD4T/Down/"

##===== 1.CAF TEC MLG ========
CAF_TEC_MLG <- read.table(paste0(wkd,"CAF_TEC_Malignant.cell/","CAF_TEC_Malignant cell.common.59.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_TEC_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_TEC_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_TEC_MLG.exp,file = paste0(wkd,"CAF_TEC_Malignant.cell/","norm.exp.xlsx"))

normEXP <- CAF_TEC_MLG.exp[rowSums(CAF_TEC_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                   filename = paste0(wkd,"CAF_TEC_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_TEC_MLG.exp.dt <- data.frame(t(CAF_TEC_MLG.exp))
CAF_TEC_MLG.exp.dt$Patient <- rownames(CAF_TEC_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_TEC_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_TEC_MLG.exp.dt,x = "Patient",y=colnames(CAF_TEC_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_TEC_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_TEC_Malignant.cell/","barplot/",colnames(CAF_TEC_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 2.CAF TEC ========
CAF_TEC <- read.table(paste0(wkd,"CAF_TEC/","CAF_TEC.common.69.txt"),header = F,stringsAsFactors = F)

expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_TEC$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}


if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_TEC.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_TEC.exp,file = paste0(wkd,"CAF_TEC/","norm.exp.xlsx"))

normEXP <- CAF_TEC.exp[rowSums(CAF_TEC.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 22,height = 6,
                   filename = paste0(wkd,"CAF_TEC/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_TEC.exp.dt <- data.frame(t(CAF_TEC.exp))
CAF_TEC.exp.dt$Patient <- rownames(CAF_TEC.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_TEC.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_TEC.exp.dt,x = "Patient",y=colnames(CAF_TEC.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_TEC.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_TEC/","barplot/",colnames(CAF_TEC.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 3.CAF Malignant cell ========
CAF_MLG <- read.table(paste0(wkd,"CAF_Malignant.cell/","CAF_Malignant cell.common.62.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_MLG.exp,file = paste0(wkd,"CAF_Malignant.cell/","norm.exp.xlsx"))

normEXP <- CAF_MLG.exp[rowSums(CAF_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                   filename = paste0(wkd,"CAF_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_MLG.exp.dt <- data.frame(t(CAF_MLG.exp))
CAF_MLG.exp.dt$Patient <- rownames(CAF_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_MLG.exp.dt,x = "Patient",y=colnames(CAF_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_Malignant.cell/","barplot/",colnames(CAF_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 4.TEC Malignant cell ========
TEC_MLG <- read.table(paste0(wkd,"TEC_Malignant.cell/","TEC_Malignant cell.common.11.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[TEC_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
TEC_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(TEC_MLG.exp,file = paste0(wkd,"TEC_Malignant.cell/","norm.exp.xlsx"))

normEXP <- TEC_MLG.exp[rowSums(TEC_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 6,height = 6,
                   filename = paste0(wkd,"TEC_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
TEC_MLG.exp.dt <- data.frame(t(TEC_MLG.exp))
TEC_MLG.exp.dt$Patient <- rownames(TEC_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(TEC_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(TEC_MLG.exp.dt,x = "Patient",y=colnames(TEC_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(TEC_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"TEC_Malignant.cell/","barplot/",colnames(TEC_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================


##>>>> CD8+T UP (CAF, TEC, Malignant cell)
###==============================================================
wkd = "result/without_PVTT_MLN/tumor/overlap/CD8T/UP/"

##===== 1.CAF TEC MLG ========
CAF_TEC_MLG <- read.table(paste0(wkd,"CAF_TEC_Malignant.cell/","CAF_TEC_Malignant cell.common.3.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_TEC_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample) 
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_TEC_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_TEC_MLG.exp,file = paste0(wkd,"CAF_TEC_Malignant.cell/","norm.exp.xlsx"))

normEXP <- CAF_TEC_MLG.exp[rowSums(CAF_TEC_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 3,height = 6,
                   filename = paste0(wkd,"CAF_TEC_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_TEC_MLG.exp.dt <- data.frame(t(CAF_TEC_MLG.exp))
CAF_TEC_MLG.exp.dt$Patient <- rownames(CAF_TEC_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_TEC_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_TEC_MLG.exp.dt,x = "Patient",y=colnames(CAF_TEC_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_TEC_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_TEC_Malignant.cell/","barplot/",colnames(CAF_TEC_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 2.CAF TEC ========
CAF_TEC <- read.table(paste0(wkd,"CAF_TEC/","CAF_TEC.common.7.txt"),header = F,stringsAsFactors = F)

expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_TEC$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_TEC.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_TEC.exp,file = paste0(wkd,"CAF_TEC/","nomr.exp.xlsx"))

normEXP <- CAF_TEC.exp[rowSums(CAF_TEC.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 6,height = 6,
                   filename = paste0(wkd,"CAF_TEC/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_TEC.exp.dt <- data.frame(t(CAF_TEC.exp))
CAF_TEC.exp.dt$Patient <- rownames(CAF_TEC.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_TEC.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_TEC.exp.dt,x = "Patient",y=colnames(CAF_TEC.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_TEC.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_TEC/","barplot/",colnames(CAF_TEC.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 3.CAF Malignant cell ========
CAF_MLG <- read.table(paste0(wkd,"CAF_Malignant.cell/","CAF_Malignant cell.common.7.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_MLG.exp,file = paste0(wkd,"CAF_Malignant.cell/","norm.exp.xlsx"))

normEXP <- CAF_MLG.exp[rowSums(CAF_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 6,height = 6,
                   filename = paste0(wkd,"CAF_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_MLG.exp.dt <- data.frame(t(CAF_MLG.exp))
CAF_MLG.exp.dt$Patient <- rownames(CAF_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_MLG.exp.dt,x = "Patient",y=colnames(CAF_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_Malignant.cell/","barplot/",colnames(CAF_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 4.TEC Malignant cell ========
TEC_MLG <- read.table(paste0(wkd,"TEC_Malignant.cell/","TEC_Malignant cell.common.3.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[TEC_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
TEC_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(TEC_MLG.exp,file = paste0(wkd,"TEC_Malignant.cell/","norm.exp.xlsx"))

normEXP <- TEC_MLG.exp[rowSums(TEC_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 3,height = 6,
                   filename = paste0(wkd,"TEC_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
TEC_MLG.exp.dt <- data.frame(t(TEC_MLG.exp))
TEC_MLG.exp.dt$Patient <- rownames(TEC_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(TEC_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(TEC_MLG.exp.dt,x = "Patient",y=colnames(TEC_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(TEC_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"TEC_Malignant.cell/","barplot/",colnames(TEC_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================


##>>>> CD8+T Down (CAF, TEC, Malignant cell)
###==============================================================
wkd = "result/without_PVTT_MLN/tumor/overlap/CD8T/Down/"

##===== 1.CAF TEC MLG ========
CAF_TEC_MLG <- read.table(paste0(wkd,"CAF_TEC_Malignant.cell/","CAF_TEC_Malignant cell.common.54.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_TEC_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_TEC_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_TEC_MLG.exp,file = paste0(wkd,"CAF_TEC_Malignant.cell/","norm.exp.xlsx"))

normEXP <- CAF_TEC_MLG.exp[rowSums(CAF_TEC_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                   filename = paste0(wkd,"CAF_TEC_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_TEC_MLG.exp.dt <- data.frame(t(CAF_TEC_MLG.exp))
CAF_TEC_MLG.exp.dt$Patient <- rownames(CAF_TEC_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_TEC_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_TEC_MLG.exp.dt,x = "Patient",y=colnames(CAF_TEC_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_TEC_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_TEC_Malignant.cell/","barplot/",colnames(CAF_TEC_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 2.CAF TEC ========
CAF_TEC <- read.table(paste0(wkd,"CAF_TEC/","CAF_TEC.common.100.txt"),header = F,stringsAsFactors = F)

expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_TEC$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_TEC.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_TEC.exp,file = paste0(wkd,"CAF_TEC/","norm.exp.xlsx"))

normEXP <- CAF_TEC.exp[rowSums(CAF_TEC.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 30,height = 6,
                   filename = paste0(wkd,"CAF_TEC/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_TEC.exp.dt <- data.frame(t(CAF_TEC.exp))
CAF_TEC.exp.dt$Patient <- rownames(CAF_TEC.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_TEC.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_TEC.exp.dt,x = "Patient",y=colnames(CAF_TEC.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_TEC.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_TEC/","barplot/",colnames(CAF_TEC.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 3.CAF Malignant cell ========
CAF_MLG <- read.table(paste0(wkd,"CAF_Malignant.cell/","CAF_Malignant cell.common.27.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[CAF_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
CAF_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(CAF_MLG.exp,file = paste0(wkd,"CAF_Malignant.cell/","norm.exp.xlsx"))

normEXP <- CAF_MLG.exp[rowSums(CAF_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 16,height = 6,
                   filename = paste0(wkd,"CAF_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
CAF_MLG.exp.dt <- data.frame(t(CAF_MLG.exp))
CAF_MLG.exp.dt$Patient <- rownames(CAF_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(CAF_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(CAF_MLG.exp.dt,x = "Patient",y=colnames(CAF_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(CAF_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"CAF_Malignant.cell/","barplot/",colnames(CAF_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================
##===== 4.TEC Malignant cell ========
TEC_MLG <- read.table(paste0(wkd,"TEC_Malignant.cell/","TEC_Malignant cell.common.60.txt"),header = F,stringsAsFactors = F)
expr <- GetAssayData(tumor_sub, slot = "data", assay = "RNA")[TEC_MLG$V1, ]
expr <- data.frame(expr)
# replace "." to "-" to match barcodes in metadata
colnames(expr) <- gsub("\\.", "-", colnames(expr))
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr), stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)

# calculate mean expression values of Overlapped DE genes in each patient
if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr)) - 1 )){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
}

if(T){
  patient <- NULL
  for(i in 1:(length(colnames(expr)) - 1)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
  }
  names(patient) <- colnames(expr)[-ncol(expr)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  Patient_sample <- data.frame(table(expr$Sample))
  patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
  patient <- data.frame(patient)
  patient$sample <- as.numeric(patient$sample)
  for(i in 1:(ncol(patient)-1)){
    for(j in 1:nrow(patient)){
      patient[j,i] <- patient[j,i]/patient$sample[j]
    }
  }
  patient$sample <- NULL
}
# final expression matrix
TEC_MLG.exp <- data.frame(t(patient))
xlsx::write.xlsx(TEC_MLG.exp,file = paste0(wkd,"TEC_Malignant.cell/","norm.exp.xlsx"))

normEXP <- TEC_MLG.exp[rowSums(TEC_MLG.exp) > 0, ]
scaleEXP <- apply(normEXP, 1, scale)
rownames(scaleEXP) <- colnames(normEXP)
scaleEXP <- t(scaleEXP)

if(max(scaleEXP) > 3) scaleEXP[scaleEXP > 3] = 3
if(min(scaleEXP) < -3) scaleEXP[scaleEXP < -3] = -3

scaleEXP <- t(scaleEXP)

pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)

pheatmap::pheatmap(scaleEXP,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                   filename = paste0(wkd,"TEC_Malignant.cell/",ncol(scaleEXP),".GeneExp.heatmap.pdf"))

### barplot
library(ggpubr)
TEC_MLG.exp.dt <- data.frame(t(TEC_MLG.exp))
TEC_MLG.exp.dt$Patient <- rownames(TEC_MLG.exp.dt)
plt <- list()
for(i in 1:(ncol(TEC_MLG.exp.dt)-1)){
  plt[[i]] <- ggbarplot(TEC_MLG.exp.dt,x = "Patient",y=colnames(TEC_MLG.exp.dt)[i],fill = "Patient") + NoLegend() +
    labs(y = "Normalized Expression") +
    ggtitle(colnames(TEC_MLG.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0(wkd,"TEC_Malignant.cell/","barplot/",colnames(TEC_MLG.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}
###===========================

###====== ADD ENSEMBL ID ========

load("refdb/GeneSymbol2ENSEMBL.RData")

### CD3 T
wkd = "result/without_PVTT_MLN/tumor/overlap/CD3T/"
## UP
# CAF_TEC_MLG
CAF_TEC_MLG <- readxl::read_excel(paste0(wkd,"UP/","CAF_TEC_Malignant.cell/","exp.xlsx"))
colnames(CAF_TEC_MLG)[1] <- "SYMBOL"
CAF_TEC_MLG$ENSEMBL <- plyr::mapvalues(CAF_TEC_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_TEC_MLG),file = paste0(wkd,"UP/","CAF_TEC_Malignant.cell/","norm.exp.xlsx"))
# CAF TEC
CAF_TEC <- readxl::read_excel(paste0(wkd,"UP/","CAF_TEC/","exp.xlsx"))
colnames(CAF_TEC)[1] <- "SYMBOL"
CAF_TEC$ENSEMBL <- plyr::mapvalues(CAF_TEC$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_TEC),file = paste0(wkd,"UP/","CAF_TEC/","norm.exp.xlsx"))
# CAF MLG
CAF_MLG <- readxl::read_excel(paste0(wkd,"UP/","CAF_Malignant.cell/","exp.xlsx"))
colnames(CAF_MLG)[1] <- "SYMBOL"
CAF_MLG$ENSEMBL <- plyr::mapvalues(CAF_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_MLG),file = paste0(wkd,"UP/","CAF_Malignant.cell/","norm.exp.xlsx"))
# TEC MLG
TEC_MLG <- readxl::read_excel(paste0(wkd,"UP/","TEC_Malignant.cell/","exp.xlsx"))
colnames(TEC_MLG)[1] <- "SYMBOL"
TEC_MLG$ENSEMBL <- plyr::mapvalues(TEC_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
TEC_MLG$ENSEMBL[3] <- "ENSG00000228253"
TEC_MLG$ENSEMBL[11] <- "ENSG00000232527"
openxlsx::write.xlsx(data.frame(TEC_MLG),file = paste0(wkd,"UP/","TEC_Malignant.cell/","norm.exp.xlsx"))

## Down
# CAF_TEC_MLG
CAF_TEC_MLG <- readxl::read_excel(paste0(wkd,"Down/","CAF_TEC_Malignant.cell/","exp.xlsx"))
colnames(CAF_TEC_MLG)[1] <- "SYMBOL"
CAF_TEC_MLG$ENSEMBL <- plyr::mapvalues(CAF_TEC_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_TEC_MLG),file = paste0(wkd,"Down/","CAF_TEC_Malignant.cell/","norm.exp.xlsx"))
# CAF TEC
CAF_TEC <- readxl::read_excel(paste0(wkd,"Down/","CAF_TEC/","exp.xlsx"))
colnames(CAF_TEC)[1] <- "SYMBOL"
CAF_TEC$ENSEMBL <- plyr::mapvalues(CAF_TEC$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_TEC),file = paste0(wkd,"Down/","CAF_TEC/","norm.exp.xlsx"))
# CAF MLG
CAF_MLG <- readxl::read_excel(paste0(wkd,"Down/","CAF_Malignant.cell/","exp.xlsx"))
colnames(CAF_MLG)[1] <- "SYMBOL"
CAF_MLG$ENSEMBL <- plyr::mapvalues(CAF_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_MLG),file = paste0(wkd,"Down/","CAF_Malignant.cell/","norm.exp.xlsx"))
# TEC MLG
TEC_MLG <- readxl::read_excel(paste0(wkd,"Down/","TEC_Malignant.cell/","exp.xlsx"))
colnames(TEC_MLG)[1] <- "SYMBOL"
TEC_MLG$ENSEMBL <- plyr::mapvalues(TEC_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(TEC_MLG),file = paste0(wkd,"Down/","TEC_Malignant.cell/","norm.exp.xlsx"))


### CD4 T
wkd = "result/without_PVTT_MLN/tumor/overlap/CD4T/"

## UP
# CAF_TEC_MLG
CAF_TEC_MLG <- readxl::read_excel(paste0(wkd,"UP/","CAF_TEC_Malignant.cell/","exp.xlsx"))
colnames(CAF_TEC_MLG)[1] <- "SYMBOL"
CAF_TEC_MLG$ENSEMBL <- plyr::mapvalues(CAF_TEC_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_TEC_MLG),file = paste0(wkd,"UP/","CAF_TEC_Malignant.cell/","norm.exp.xlsx"))
# CAF TEC
CAF_TEC <- readxl::read_excel(paste0(wkd,"UP/","CAF_TEC/","exp.xlsx"))
colnames(CAF_TEC)[1] <- "SYMBOL"
CAF_TEC$ENSEMBL <- plyr::mapvalues(CAF_TEC$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_TEC),file = paste0(wkd,"UP/","CAF_TEC/","norm.exp.xlsx"))
# CAF MLG
CAF_MLG <- readxl::read_excel(paste0(wkd,"UP/","CAF_Malignant.cell/","exp.xlsx"))
colnames(CAF_MLG)[1] <- "SYMBOL"
CAF_MLG$ENSEMBL <- plyr::mapvalues(CAF_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_MLG),file = paste0(wkd,"UP/","CAF_Malignant.cell/","norm.exp.xlsx"))
# TEC MLG
TEC_MLG <- readxl::read_excel(paste0(wkd,"UP/","TEC_Malignant.cell/","exp.xlsx"))
colnames(TEC_MLG)[1] <- "SYMBOL"
TEC_MLG$ENSEMBL <- plyr::mapvalues(TEC_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
TEC_MLG$ENSEMBL[3] <- "ENSG00000228253"
TEC_MLG$ENSEMBL[11] <- "ENSG00000232527"
openxlsx::write.xlsx(data.frame(TEC_MLG),file = paste0(wkd,"UP/","TEC_Malignant.cell/","norm.exp.xlsx"))

## Down
# CAF_TEC_MLG
CAF_TEC_MLG <- readxl::read_excel(paste0(wkd,"Down/","CAF_TEC_Malignant.cell/","exp.xlsx"))
colnames(CAF_TEC_MLG)[1] <- "SYMBOL"
CAF_TEC_MLG$ENSEMBL <- plyr::mapvalues(CAF_TEC_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_TEC_MLG),file = paste0(wkd,"Down/","CAF_TEC_Malignant.cell/","norm.exp.xlsx"))
# CAF TEC
CAF_TEC <- readxl::read_excel(paste0(wkd,"Down/","CAF_TEC/","exp.xlsx"))
colnames(CAF_TEC)[1] <- "SYMBOL"
CAF_TEC$ENSEMBL <- plyr::mapvalues(CAF_TEC$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_TEC),file = paste0(wkd,"Down/","CAF_TEC/","norm.exp.xlsx"))
# CAF MLG
CAF_MLG <- readxl::read_excel(paste0(wkd,"Down/","CAF_Malignant.cell/","exp.xlsx"))
colnames(CAF_MLG)[1] <- "SYMBOL"
CAF_MLG$ENSEMBL <- plyr::mapvalues(CAF_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_MLG),file = paste0(wkd,"Down/","CAF_Malignant.cell/","norm.exp.xlsx"))
# TEC MLG
TEC_MLG <- readxl::read_excel(paste0(wkd,"Down/","TEC_Malignant.cell/","exp.xlsx"))
colnames(TEC_MLG)[1] <- "SYMBOL"
TEC_MLG$ENSEMBL <- plyr::mapvalues(TEC_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(TEC_MLG),file = paste0(wkd,"Down/","TEC_Malignant.cell/","norm.exp.xlsx"))


### CD8 T
wkd = "result/without_PVTT_MLN/tumor/overlap/CD8T/"

## UP
# CAF_TEC_MLG
CAF_TEC_MLG <- readxl::read_excel(paste0(wkd,"UP/","CAF_TEC_Malignant.cell/","exp.xlsx"))
colnames(CAF_TEC_MLG)[1] <- "SYMBOL"
CAF_TEC_MLG$ENSEMBL <- plyr::mapvalues(CAF_TEC_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_TEC_MLG),file = paste0(wkd,"UP/","CAF_TEC_Malignant.cell/","norm.exp.xlsx"))
# CAF TEC
CAF_TEC <- readxl::read_excel(paste0(wkd,"UP/","CAF_TEC/","exp.xlsx"))
colnames(CAF_TEC)[1] <- "SYMBOL"
CAF_TEC$ENSEMBL <- plyr::mapvalues(CAF_TEC$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_TEC),file = paste0(wkd,"UP/","CAF_TEC/","norm.exp.xlsx"))
# CAF MLG
CAF_MLG <- readxl::read_excel(paste0(wkd,"UP/","CAF_Malignant.cell/","exp.xlsx"))
colnames(CAF_MLG)[1] <- "SYMBOL"
CAF_MLG$ENSEMBL <- plyr::mapvalues(CAF_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_MLG),file = paste0(wkd,"UP/","CAF_Malignant.cell/","norm.exp.xlsx"))
# TEC MLG
TEC_MLG <- readxl::read_excel(paste0(wkd,"UP/","TEC_Malignant.cell/","exp.xlsx"))
colnames(TEC_MLG)[1] <- "SYMBOL"
TEC_MLG$ENSEMBL <- plyr::mapvalues(TEC_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(TEC_MLG),file = paste0(wkd,"UP/","TEC_Malignant.cell/","norm.exp.xlsx"))

## Down
# CAF_TEC_MLG
CAF_TEC_MLG <- readxl::read_excel(paste0(wkd,"Down/","CAF_TEC_Malignant.cell/","exp.xlsx"))
colnames(CAF_TEC_MLG)[1] <- "SYMBOL"
CAF_TEC_MLG$ENSEMBL <- plyr::mapvalues(CAF_TEC_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_TEC_MLG),file = paste0(wkd,"Down/","CAF_TEC_Malignant.cell/","norm.exp.xlsx"))
# CAF TEC
CAF_TEC <- readxl::read_excel(paste0(wkd,"Down/","CAF_TEC/","exp.xlsx"))
colnames(CAF_TEC)[1] <- "SYMBOL"
CAF_TEC$ENSEMBL <- plyr::mapvalues(CAF_TEC$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
CAF_TEC$ENSEMBL[10] <- "ENSG00000196126"
CAF_TEC$ENSEMBL[85] <- "ENSG00000172965"
openxlsx::write.xlsx(data.frame(CAF_TEC),file = paste0(wkd,"Down/","CAF_TEC/","norm.exp.xlsx"))
# CAF MLG
CAF_MLG <- readxl::read_excel(paste0(wkd,"Down/","CAF_Malignant.cell/","exp.xlsx"))
colnames(CAF_MLG)[1] <- "SYMBOL"
CAF_MLG$ENSEMBL <- plyr::mapvalues(CAF_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(CAF_MLG),file = paste0(wkd,"Down/","CAF_Malignant.cell/","norm.exp.xlsx"))
# TEC MLG
TEC_MLG <- readxl::read_excel(paste0(wkd,"Down/","TEC_Malignant.cell/","exp.xlsx"))
colnames(TEC_MLG)[1] <- "SYMBOL"
TEC_MLG$ENSEMBL <- plyr::mapvalues(TEC_MLG$SYMBOL, from = mapped$symbol, to = mapped$ENSEMBL)
openxlsx::write.xlsx(data.frame(TEC_MLG),file = paste0(wkd,"Down/","TEC_Malignant.cell/","norm.exp.xlsx"))

###===========================================================


##========================== 2. Obtain expression matrix from UP/Down group ==========================

if(FALSE){
  tumor <- readRDS("rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR.rds")
  
  ## CD4+T low vs CD+8 T high
  high.genes <- read.table("result/without_PVTT_MLN/tumor/overlapped/high/CD4T_CD8T_high.CD4+T high.and.CD8+T high.common.65.txt",header = F,stringsAsFactors = F)
  
  # obtain expression matrix (raw count)
  expr <- GetAssayData(tumor,slot = "counts",assay = "RNA")[high.genes$V1, ]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.integer(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-66]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  ol.high.exp <- data.frame(t(tmp))
  xlsx::write.xlsx(ol.high.exp,file = "result/without_PVTT_MLN/tumor/overlapped/high/CD4T_CD8T_high.geneExp.xlsx",sheetName = "raw")
  
  normDt <- log2(ol.high.exp + 1)
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  xlsx::write.xlsx(scaleDt,file = "result/without_PVTT_MLN/tumor/overlapped/high/CD4T_CD8T_high.geneExp.xlsx",sheetName = "normalized",append = T)
  
  scaleDt[scaleDt > 2] = 2
  scaleDt[scaleDt < -2 ] = -2
  scaleDt <- t(scaleDt)
  pheatmap(scaleDt,border=F,color = pal1,fontsize = 14,angle_col = 45,width = 20,height = 6,filename = "result/without_PVTT_MLN/tumor/overlapped/high/ol.65.heatmap.png")
  
  ## CD4+T low vs CD+8 T low
  low.genes <- read.table("result/without_PVTT_MLN/tumor/overlapped/low/CD4T_CD8T_low.CD4+T low.and.CD8+T low.common.204.txt",header = F,stringsAsFactors = F)
  expr <- GetAssayData(tumor,assay = "RNA",slot = "counts")[low.genes$V1,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.integer(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-205]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  ol.low.exp <- data.frame(t(tmp))
  xlsx::write.xlsx(ol.low.exp,file = "result/without_PVTT_MLN/tumor/overlapped/low/CD4T_CD8T_low.geneExp.xlsx",sheetName = "raw")
  
  normDt <- log2(ol.low.exp + 1)
  scaleDt <- apply(normDt,1,scale)
  scaleDt <- t(scaleDt)
  colnames(scaleDt) <- colnames(normDt)
  xlsx::write.xlsx(scaleDt,file = "result/without_PVTT_MLN/tumor/overlapped/low/CD4T_CD8T_low.geneExp.xlsx",sheetName = "normalized",append = T)
  
  scaleDt[scaleDt > 2] = 2
  scaleDt[scaleDt < -2 ] = -2
  scaleDt <- t(scaleDt)
  pheatmap(scaleDt,border=F,color = pal1,angle_col = 45,width = 30,height = 6,filename = "result/without_PVTT_MLN/tumor/overlapped/low/ol.204.heatmap.png")
  
  ### barplot
  library(ggpubr)
  
  # CD4+T vs CD8+T high
  ol.high.exp.dt <- data.frame(t(ol.high.exp))
  ol.high.exp.dt$Patient <- rownames(ol.high.exp.dt)
  plt <- list()
  for(i in 1:(ncol(ol.high.exp.dt)-1)){
    ol.high.exp.dt <- ol.high.exp.dt[order(ol.high.exp.dt[,i],decreasing = T),]
    plt[[i]] <- ggbarplot(ol.high.exp.dt,x = "Patient",y=colnames(ol.high.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(ol.high.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/overlapped/high/barplot/",
                  colnames(ol.high.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 4)
  }
  
  # CD4+T vs CD8+T low
  ol.low.exp.dt <- data.frame(t(ol.low.exp))
  ol.low.exp.dt$Patient <- rownames(ol.low.exp.dt)
  plt <- list()
  for(i in 1:(ncol(ol.low.exp.dt)-1)){
    ol.low.exp.dt <- ol.low.exp.dt[order(ol.low.exp.dt[,i],decreasing = T),]
    plt[[i]] <- ggbarplot(ol.low.exp.dt, x = "Patient", y = colnames(ol.low.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(ol.low.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/overlapped/low/barplot/",
                  colnames(ol.low.exp.dt)[i],".pdf"),plt[[i]], width = 8, height = 4)
  }
}

## >>>>>>>>>>>>>>>>>>>>>> CD3+T  <<<<<<<<<<<<<<<<<<<<<<<<<<<

tumor_sub <- readRDS("rds/without_PVTT_MLN/tumor/CD45_removal_T.IR.2.rds")

## >>>>>> CAF <<<<<<<

if(FALSE){
  ## UP
  CD3T.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",sheet = 1)
  colnames(CD3T.CAF.up)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD3T.CAF.up$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD3T.CAF.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD3T.CAF.up.exp,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.up.exp.xlsx",sheetName = "raw")
  
  normDt <- log2(CD3T.CAF.up.exp + 1)
  normDt <- normDt[rowSums(normDt) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  xlsx::write.xlsx(scaleDt,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.up.exp.xlsx",sheetName = "normalized",append = T)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/",
                                       ncol(scaleDt),".CAF.up.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD3T.CAF.up.exp.dt <- data.frame(t(CD3T.CAF.up.exp))
  CD3T.CAF.up.exp.dt$Patient <- rownames(CD3T.CAF.up.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD3T.CAF.up.exp.dt)-1)){
    # CD3T.CAF.up.exp.dt <- CD3T.CAF.up.exp.dt[order(CD3T.CAF.up.exp.dt[,i],decreasing = T),]
    plt[[i]] <- ggbarplot(CD3T.CAF.up.exp.dt,x = "Patient",y=colnames(CD3T.CAF.up.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD3T.CAF.up.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/barplot/UP/",
                  colnames(CD3T.CAF.up.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

########################################################################################################

## UP
CD3T.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD3T.CAF.up)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD3T.CAF.up$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group1)
# only calculate expression values in CAF cell
expr <- expr[expr$CellType == "CAF", ]
expr$CellType <- NULL
# calculate mean expression values of DE gens in each patient

## byGroup
# High: Sum(T1+T2+T3+..Tn)/Sum(N1+N2+N3+..+Nn) = T1/Sum(N1+...+Nn) + T2/Sum(N1+..Nn) + ... + Tn/Sum(N1+...Nn)
# Low: Sum(T1+T2+T3+..Tn)/Sum(N1+N2+N3+..+Nn) = T1/Sum(N1+...+Nn) + T2/Sum(N1+..Nn) + ... + Tn/Sum(N1+...Nn)

# CAF Malignant cell   TEC
# CD3+T high  1911           7491  1805
# CD3+T low    598          16686  1474

## bySample
# T1/N1 + T2/N2 + ... Tn/Nn


if(F){
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-2)){
    for(j in seq_along(names(table(expr$Sample)))){
      if(expr$Group1 == "CD3+T high"){
        patient[[j]] <- sum(as.numeric(expr[1:table(expr$Sample)[j],i]))/table(expr$Group1)[1]
      }else if (expr$Group1 == "CD3+T low"){
        patient[[j]] <- sum(as.numeric(expr[1:table(expr$Sample)[j],i]))/table(expr$Group1)[2]
      }
      # patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
}

#### Normalized gene expression & visualization ####
library(ggpubr)

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group1)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD3+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c("H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC05T","HCC9","HCC03T","HCC1")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.norm.exp_byGroup.xlsx"), sheetName = "UP")
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.norm.exp_bySample.xlsx"), sheetName = "UP")
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD3T = rep(c("CD3+T high","CD3+T low"), c(12,9)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.up.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.up.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.up.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.up.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/barplot/UP/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/barplot/UP/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")


if(F){
  tmp <- data.frame(t(patient))
  HCC <- data.frame(HCC2 = rep(0,nrow(tmp)), HCC5 = rep(0,nrow(tmp)), HCC9 = rep(0, nrow(tmp)), row.names = rownames(tmp))
  tmp <- cbind(tmp,HCC)
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(4,7,3,13,6,8,20,18,15,1,16,14,5,19,11,9,2,12,21,10,17)
  
  CD3T.CAF.up.exp <- tmp[,labs]
  CD3T.CAF.up.exp <- CD3T.CAF.up.exp[rowSums(CD3T.CAF.up.exp) > 0, ]
  CD3T.CAF.up.exp <- round(CD3T.CAF.up.exp, digits = 3)
  xlsx::write.xlsx(CD3T.CAF.up.exp,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.norm.exp.xlsx",sheetName = "UP")
  
  scaleDt <- apply(CD3T.CAF.up.exp,1,scale)
  rownames(scaleDt) <- colnames(CD3T.CAF.up.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  anno_col <- data.frame(CD3T = rep(c("CD3+T high","CD3+T low"), c(12,9)),row.names = rownames(scaleDt))
  
  pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                     show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                     annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/",
                                       ncol(scaleDt),".CAF.up.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD3T.CAF.up.exp.dt <- data.frame(t(CD3T.CAF.up.exp))
  CD3T.CAF.up.exp.dt$Patient <- rownames(CD3T.CAF.up.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD3T.CAF.up.exp.dt)-1)){
    # CD3T.CAF.up.exp.dt <- CD3T.CAF.up.exp.dt[order(CD3T.CAF.up.exp.dt[,i],decreasing = T),]
    plt[[i]] <- ggbarplot(CD3T.CAF.up.exp.dt,x = "Patient",y=colnames(CD3T.CAF.up.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD3T.CAF.up.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/barplot/UP/",
                  colnames(CD3T.CAF.up.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

### FeaturePlot
library(ggplot2)
library(Seurat)
library(patchwork)

DefaultAssay(tumor_sub) <- "RNA"

CD3T.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD3T.CAF.up)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
        '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group1",cols = pal3)

for(i in seq_along(CD3T.CAF.up$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD3T.CAF.up$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD3T.CAF.UP]", "Finished FeaturePlot for ------>",i,CD3T.CAF.up$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/FeaturePlot/UP/",CD3T.CAF.up$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}

########################################################################################################

if(FALSE){
  ## DOWN
  CD3T.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",sheet = 2)
  colnames(CD3T.CAF.down)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD3T.CAF.down$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD3T.CAF.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD3T.CAF.down.exp,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.down.exp.xlsx",sheetName = "raw")
  
  normDt <- log2(CD3T.CAF.down.exp + 1)
  normDt <- normDt[rowSums(normDt) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  xlsx::write.xlsx(scaleDt,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.down.exp.xlsx",sheetName = "normalized",append = T)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/",
                                       ncol(scaleDt),".CAF.down.genes.heatmap.pdf"))
}

########################################################################################################

## DOWN
CD3T.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD3T.CAF.down)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD3T.CAF.down$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group1)
expr <- expr[expr$CellType == "CAF",]
expr$CellType = NULL

# calculate mean expression values of DE gens in each patient

if(F){
  if(F){
    patient <- NULL
    Patient <- list()
    for(i in 1:(length(colnames(expr))-2)){
      for(j in 1:(length(names(table(expr$Sample))))){
        if(j == 1){
          if(expr$Group1[table(expr$Sample)[j]] == "CD3+T high"){
            patient[[j]] <- sum(as.numeric(expr[1:table(expr$Sample)[j],i]))/table(expr$Group1)[1]
          }else {
            patient[[j]] <- sum(as.numeric(expr[1:table(expr$Sample)[j],i]))/table(expr$Group1)[2]
          }
        }else {
          if(expr$Group1[table(expr$Sample)[j]] == "CD3+T high"){
            patient[[j]] <- sum(as.numeric(expr[(table(expr$Sample)[j-1]+1):(table(expr$Sample)[j-1]+table(expr$Sample)[j]), i]))/table(expr$Group1)[1]
          }else {
            patient[[j]] <- sum(as.numeric(expr[(table(expr$Sample)[j-1]+1):(table(expr$Sample)[j-1]+table(expr$Sample)[j]), i]))/table(expr$Group1)[2]
          }
        }
        Patient[[i]] <- patient
      }
    }
  }
  
  patient <- NULL
  for(i in 1:(length(colnames(expr))-2)){
    patient[i] <- data.frame(tapply(expr[,i], expr$Sample, sum))
    names(patient[i]) <- colnames(expr)[i]
  }
  
  names(patient) <- colnames(expr)[-c(ncol(expr),ncol(expr)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(expr[,1], expr$Sample, sum))
  
  patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group1)
  
  for(i in 1:(ncol(patient)-1)){
    for (j in 1:nrow(patient)) {
      if(patient$Group[j] == "CD+3 high"){
        patient[j,i] <- patient[j,i]/table(expr$Group)[1]
      }else {
        patient[j,i] <- patient[j,i]/table(expr$Group)[2]
      }
    } 
  }
  patient$Group <- NULL
  
  if(F){
    Patient_sample <- data.frame(table(expr$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    for(i in 1:(ncol(patient)-1)){
      for (j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  
  tmp <- data.frame(t(patient))
  HCC <- data.frame(HCC2 = rep(0,nrow(tmp)), HCC5 = rep(0,nrow(tmp)), HCC9 = rep(0, nrow(tmp)), row.names = rownames(tmp))
  tmp <- cbind(tmp,HCC)
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(4,7,3,13,6,8,20,18,15,1,16,14,5,19,11,9,2,12,21,10,17)
  CD3T.CAF.down.exp <- tmp[,labs]
  CD3T.CAF.down.exp <- CD3T.CAF.down.exp[rowSums(CD3T.CAF.down.exp) > 0,]
  CD3T.CAF.down.exp <- round(CD3T.CAF.down.exp, digits = 3)
  xlsx::write.xlsx(CD3T.CAF.down.exp,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.norm.exp.xlsx",sheetName = "Down",append = T)
  
  scaleDt <- apply(CD3T.CAF.down.exp,1,scale)
  rownames(scaleDt) <- colnames(CD3T.CAF.down.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD3T = rep(c("CD3+T high","CD3+T low"), c(12,9)),row.names = rownames(scaleDt))
  pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                     show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                     annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/",
                                       ncol(scaleDt),".CAF.down.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD3T.CAF.down.exp.dt <- data.frame(t(CD3T.CAF.down.exp))
  CD3T.CAF.down.exp.dt$Patient <- rownames(CD3T.CAF.down.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD3T.CAF.down.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD3T.CAF.down.exp.dt,x = "Patient",y=colnames(CD3T.CAF.down.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD3T.CAF.down.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/barplot/Down/",
                  colnames(CD3T.CAF.down.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group1)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD3+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70",
            "HCC2","HCC04T","HCC02T","H38","HCC05T","HCC9","HCC03T","HCC1")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.norm.exp_byGroup.xlsx"), sheetName = "Down",append = TRUE)
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.norm.exp_bySample.xlsx"), sheetName = "Down",append = TRUE)
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD3T = rep(c("CD3+T high","CD3+T low"), c(12,9)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.down.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.down.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.down.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.down.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/barplot/Down/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/barplot/Down/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD3T.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD3T.CAF.down)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group1",cols = pal3)

for(i in seq_along(CD3T.CAF.down$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD3T.CAF.down$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD3T.CAF.DOWN]", "Finished FeaturePlot for ------>",i,CD3T.CAF.down$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/FeaturePlot/Down/",CD3T.CAF.down$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}


## >>>>>> TEC <<<<<<<

########################################################################################################
if(FALSE){
  ## UP
  CD3T.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",sheet = 1)
  colnames(CD3T.TEC.up)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD3T.TEC.up$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD3T.TEC.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD3T.TEC.up.exp,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.up.exp.xlsx",sheetName = "normalized")
}

## UP
CD3T.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD3T.TEC.up)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD3T.TEC.up$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group1)
expr <- expr[expr$CellType == "TEC",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD3T.TEC.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD3T.TEC.up.exp,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.norm.exp.xlsx",sheetName = "UP")
  
  
  normDt <- CD3T.TEC.up.exp[rowSums(CD3T.TEC.up.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/",
                                       ncol(scaleDt),".TEC.up.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD3T.TEC.up.exp.dt <- data.frame(t(CD3T.TEC.up.exp))
  CD3T.TEC.up.exp.dt$Patient <- rownames(CD3T.TEC.up.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD3T.TEC.up.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD3T.TEC.up.exp.dt,x = "Patient",y=colnames(CD3T.TEC.up.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD3T.TEC.up.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/barplot/UP/",
                  colnames(CD3T.TEC.up.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group1)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD3+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70",
            "HCC2","HCC04T","HCC02T","H38","HCC05T","HCC9","HCC03T","HCC1")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.norm.exp_byGroup.xlsx"), sheetName = "UP")
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.norm.exp_bySample.xlsx"), sheetName = "UP")
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD3T = rep(c("CD3+T high","CD3+T low"), c(12,9)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.up.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.up.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.up.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.up.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/barplot/UP/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/barplot/UP/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD3T.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD3T.TEC.up)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group1",cols = pal3)

for(i in seq_along(CD3T.TEC.up$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD3T.TEC.up$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD3T.TEC.UP]", "Finished FeaturePlot for ------>",i,CD3T.TEC.up$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/FeaturePlot/UP/",CD3T.TEC.up$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}



########################################################################################################
## DOWN
if(FALSE){
  CD3T.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",sheet = 2)
  colnames(CD3T.TEC.down)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD3T.TEC.down$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD3T.TEC.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD3T.TEC.down.exp,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.down.exp.xlsx",sheetName = "raw")
}

CD3T.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD3T.TEC.down)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD3T.TEC.down$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group1)
expr <- expr[expr$CellType == "TEC",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD3T.TEC.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD3T.TEC.down.exp,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.norm.exp.xlsx",sheetName = "Down",append = TRUE)
  
  normDt <- CD3T.TEC.down.exp[rowSums(CD3T.TEC.down.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/",
                                       ncol(scaleDt),".TEC.down.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD3T.TEC.down.exp.dt <- data.frame(t(CD3T.TEC.down.exp))
  CD3T.TEC.down.exp.dt$Patient <- rownames(CD3T.TEC.down.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD3T.TEC.down.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD3T.TEC.down.exp.dt,x = "Patient",y=colnames(CD3T.TEC.down.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD3T.TEC.down.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/barplot/Down/",
                  colnames(CD3T.TEC.down.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group1)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD3+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70",
            "HCC2","HCC04T","HCC02T","H38","HCC05T","HCC9","HCC03T","HCC1")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.norm.exp_byGroup.xlsx"), sheetName = "Down",append = TRUE)
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.norm.exp_bySample.xlsx"), sheetName = "Down",append = TRUE)
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD3T = rep(c("CD3+T high","CD3+T low"), c(12,9)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.down.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.down.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.down.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.down.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/barplot/Down/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/barplot/Down/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD3T.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD3T.TEC.down)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group1",cols = pal3)

for(i in seq_along(CD3T.TEC.down$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD3T.TEC.down$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD3T.TEC.DOWN]", "Finished FeaturePlot for ------>",i,CD3T.TEC.down$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/FeaturePlot/Down/",CD3T.TEC.down$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}


########################################################################################################

## >>>>>> Malignant Cell <<<<<<<

if(FALSE){
  ## UP
  CD3T.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Maligant.Cell.CD3T.hi_vs_lo.xlsx",sheet = 1)
  colnames(CD3T.MLG.up)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD3T.MLG.up$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD3T.MLG.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD3T.MLG.up.exp,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.cell.up.exp.xlsx",sheetName = "raw")
}
## UP
CD3T.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD3T.MLG.up)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD3T.MLG.up$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
# extract malignant cells
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group1)
expr <- expr[expr$CellType == "Malignant cell",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr_mlg))-1)){
    for(j in seq_along(names(table(expr_mlg$Sample)))){
      patient[j] <- mean(as.numeric(expr_mlg[1:table(expr_mlg$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr_mlg)[-ncol(expr_mlg)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- names(table(expr_mlg$Sample))
  # final expression matrix
  tmp <- data.frame(t(tmp))
  H21 <- data.frame(H21 = rep(0,nrow(tmp)),row.names = rownames(tmp))
  tmp <- cbind(H21,tmp)
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(4,7,3,13,6,8,20,18,15,1,16,14,5,19,11,9,2,12,21,10,17)
  CD3T.MLG.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD3T.MLG.up.exp,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.cell.norm.exp.xlsx",sheetName = "UP")
  
  normDt <- CD3T.MLG.up.exp[rowSums(CD3T.MLG.up.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/",
                                       ncol(scaleDt),".Malignant.cell.up.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD3T.MLG.up.exp.dt <- data.frame(t(CD3T.MLG.up.exp))
  CD3T.MLG.up.exp.dt$Patient <- rownames(CD3T.MLG.up.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD3T.MLG.up.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD3T.MLG.up.exp.dt,x = "Patient",y=colnames(CD3T.MLG.up.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD3T.MLG.up.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/barplot/UP/",
                  colnames(CD3T.MLG.up.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group1)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD3+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70",
            "HCC2","HCC04T","HCC02T","H38","HCC05T","HCC9","HCC03T","HCC1")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant_cell.norm.exp_byGroup.xlsx"), sheetName = "UP")
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant_cell.norm.exp_bySample.xlsx"), sheetName = "UP")
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD3T = rep(c("CD3+T high","CD3+T low"), c(12,9)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.up.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.up.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.up.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.up.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/barplot/UP/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/barplot/UP/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD3T.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD3T.MLG.up)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group1",cols = pal3)

for(i in seq_along(CD3T.MLG.up$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD3T.MLG.up$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD3T.MLG.UP]", "Finished FeaturePlot for ------>",i,CD3T.MLG.up$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/FeaturePlot/UP/",CD3T.MLG.up$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}


########################################################################################################

## DOWN
if(FALSE){
  CD3T.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Maligant.Cell.CD3T.hi_vs_lo.xlsx",sheet = 2)
  colnames(CD3T.MLG.down)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD3T.MLG.down$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD3T.MLG.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD3T.MLG.down.exp,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.cell.norm.exp.xlsx",sheetName = "UP")
}

CD3T.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD3T.MLG.down)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD3T.MLG.down$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
# extract malignant cells
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group1)
expr <- expr[expr$CellType == "Malignant cell",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD3T.MLG.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD3T.MLG.down.exp,file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.cell.norm.exp.xlsx",sheetName = "Down",append = TRUE)
  
  
  normDt <- CD3T.MLG.down.exp[rowSums(CD3T.MLG.down.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/",
                                       ncol(scaleDt),".Malignant.cell.down.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD3T.MLG.down.exp.dt <- data.frame(t(CD3T.MLG.down.exp))
  CD3T.MLG.down.exp.dt$Patient <- rownames(CD3T.MLG.down.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD3T.MLG.down.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD3T.MLG.down.exp.dt,x = "Patient",y=colnames(CD3T.MLG.down.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD3T.MLG.down.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/barplot/Down/",
                  colnames(CD3T.MLG.down.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group1)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD3+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70",
            "HCC2","HCC04T","HCC02T","H38","HCC05T","HCC9","HCC03T","HCC1")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant_cell.norm.exp_byGroup.xlsx"), sheetName = "Down",append = TRUE)
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant_cell.norm.exp_bySample.xlsx"), sheetName = "Down",append = TRUE)
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD3T = rep(c("CD3+T high","CD3+T low"), c(12,9)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.down.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.down.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.down.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.down.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/barplot/Down/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/barplot/Down/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD3T.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD3T.MLG.down)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group1",cols = pal3)

for(i in seq_along(CD3T.MLG.down$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD3T.MLG.down$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD3T.MLG.DOWN]", "Finished FeaturePlot for ------>",i,CD3T.MLG.down$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/FeaturePlot/Down/",CD3T.MLG.down$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}


## >>>>>>>>>>>>>>>>>>>>>> CD4+T  <<<<<<<<<<<<<<<<<<<<<<<<<<<

########################################################################################################

## >>>>>> CAF <<<<<<<

if(FALSE){
  ## UP
  CD4T.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 1)
  colnames(CD4T.CAF.up)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD4T.CAF.up$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","HCC01T","HCC08T","HCC10T","H21","H72","HCC09T","HCC07T","HCC5","H70","HCC2","HCC04T","H38","HCC9","HCC1","HCC05T","HCC02T","HCC03T"
  labs <- c(20,19,21,16,7,8,14,1,18,13,15,5,17,4,11,2,6,3,9,10,12)
  CD4T.CAF.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD4T.CAF.up.exp,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.up.exp.xlsx",sheetName = "raw")
}

## UP
CD4T.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD4T.CAF.up)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD4T.CAF.up$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group2)
expr <- expr[expr$CellType == "CAF",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","HCC01T","HCC08T","HCC10T","H21","H72","HCC09T","HCC07T","HCC5","H70","HCC2","HCC04T","H38","HCC9","HCC1","HCC05T","HCC02T","HCC03T"
  labs <- c(20,19,21,16,7,8,14,1,18,13,15,5,17,4,11,2,6,3,9,10,12)
  CD4T.CAF.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD4T.CAF.up.exp,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.norm.exp.xlsx",sheetName = "UP")
  
  normDt <- CD4T.CAF.up.exp[rowSums(CD4T.CAF.up.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/",
                                       ncol(scaleDt),".CAF.up.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD4T.CAF.up.exp.dt <- data.frame(t(CD4T.CAF.up.exp))
  CD4T.CAF.up.exp.dt$Patient <- rownames(CD4T.CAF.up.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD4T.CAF.up.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD4T.CAF.up.exp.dt,x = "Patient",y=colnames(CD4T.CAF.up.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD4T.CAF.up.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/barplot/UP/",
                  colnames(CD4T.CAF.up.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group2)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD4+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("H62","H75","H41","HCC06T","HCC01T","HCC08T","HCC10T","H21","H72","HCC09T","HCC07T","HCC5","H70",
            "HCC2","HCC04T","H38","HCC9","HCC1","HCC05T","HCC02T","HCC03T")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.norm.exp_byGroup.xlsx"), sheetName = "UP")
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.norm.exp_bySample.xlsx"), sheetName = "UP")
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD4T = rep(c("CD4+T high","CD4+T low"), c(12,9)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.up.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.up.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.up.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.up.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/barplot/UP/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/barplot/UP/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD4T.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD4T.CAF.up)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group2",cols = pal3)

for(i in seq_along(CD4T.CAF.up$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD4T.CAF.up$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD4T.CAF.UP]", "Finished FeaturePlot for ------>",i,CD4T.CAF.up$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/FeaturePlot/UP/",CD4T.CAF.up$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}


########################################################################################################

## DOWN

if(FALSE){
  CD4T.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 2)
  colnames(CD4T.CAF.down)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD4T.CAF.down$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","HCC01T","HCC08T","HCC10T","H21","H72","HCC09T","HCC07T","HCC5","H70","HCC2","HCC04T","H38","HCC9","HCC1","HCC05T","HCC02T","HCC03T"
  labs <- c(20,19,21,16,7,8,14,1,18,13,15,5,17,4,11,2,6,3,9,10,12)
  CD4T.CAF.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD4T.CAF.down.exp,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.down.exp.xlsx",sheetName = "raw")
}

CD4T.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD4T.CAF.down)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD4T.CAF.down$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group2)
expr <- expr[expr$CellType == "CAF",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","HCC01T","HCC08T","HCC10T","H21","H72","HCC09T","HCC07T","HCC5","H70","HCC2","HCC04T","H38","HCC9","HCC1","HCC05T","HCC02T","HCC03T"
  labs <- c(20,19,21,16,7,8,14,1,18,13,15,5,17,4,11,2,6,3,9,10,12)
  CD4T.CAF.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD4T.CAF.down.exp,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.norm.exp.xlsx",sheetName = "Down",append = TRUE)
  
  normDt <- CD4T.CAF.down.exp[rowSums(CD4T.CAF.down.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/",
                                       ncol(scaleDt),".CAF.down.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD4T.CAF.down.exp.dt <- data.frame(t(CD4T.CAF.down.exp))
  CD4T.CAF.down.exp.dt$Patient <- rownames(CD4T.CAF.down.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD4T.CAF.down.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD4T.CAF.down.exp.dt,x = "Patient",y=colnames(CD4T.CAF.down.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD4T.CAF.down.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/barplot/Down/",
                  colnames(CD4T.CAF.down.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group2)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD4+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("H62","H75","H41","HCC06T","HCC01T","HCC08T","HCC10T","H21","H72","HCC09T","HCC07T","HCC5","H70",
            "HCC2","HCC04T","H38","HCC9","HCC1","HCC05T","HCC02T","HCC03T")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.norm.exp_byGroup.xlsx"), sheetName = "Down",append = TRUE)
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.norm.exp_bySample.xlsx"), sheetName = "Down",append = TRUE)
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD4T = rep(c("CD4+T high","CD4+T low"), c(12,9)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.down.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.down.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.down.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.down.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/barplot/Down/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/barplot/Down/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD4T.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD4T.CAF.down)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group2",cols = pal3)

for(i in seq_along(CD4T.CAF.down$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD4T.CAF.down$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD4T.CAF.DOWN]","Finished FeaturePlot for ------>",i,CD4T.CAF.down$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/FeaturePlot/Down/",CD4T.CAF.down$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}


## >>>>>> TEC <<<<<<<

########################################################################################################

if(FALSE){
  ## UP
  CD4T.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 1)
  colnames(CD4T.TEC.up)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD4T.TEC.up$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD4T.TEC.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD4T.TEC.up.exp,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.up.exp.xlsx",sheetName = "raw")
}

## UP
CD4T.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD4T.TEC.up)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD4T.TEC.up$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group2)
expr <- expr[expr$CellType == "TEC",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD4T.TEC.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD4T.TEC.up.exp,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.norm.exp.xlsx",sheetName = "UP")
  
  normDt <- CD4T.TEC.up.exp[rowSums(CD4T.TEC.up.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/",
                                       ncol(scaleDt),".TEC.up.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD4T.TEC.up.exp.dt <- data.frame(t(CD4T.TEC.up.exp))
  CD4T.TEC.up.exp.dt$Patient <- rownames(CD4T.TEC.up.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD4T.TEC.up.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD4T.TEC.up.exp.dt,x = "Patient",y=colnames(CD4T.TEC.up.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD4T.TEC.up.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/barplot/UP/",
                  colnames(CD4T.TEC.up.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group2)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD4+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("H62","H75","H41","HCC06T","HCC01T","HCC08T","HCC10T","H21","H72","HCC09T","HCC07T","HCC5","H70",
            "HCC2","HCC04T","H38","HCC9","HCC1","HCC05T","HCC02T","HCC03T")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.norm.exp_byGroup.xlsx"), sheetName = "UP")
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.norm.exp_bySample.xlsx"), sheetName = "UP")
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD4T = rep(c("CD4+T high","CD4+T low"), c(12,9)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.up.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.up.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.up.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.up.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/barplot/UP/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/barplot/UP/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD4T.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD4T.TEC.up)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group2",cols = pal3)

for(i in seq_along(CD4T.TEC.up$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD4T.TEC.up$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD4T.TEC.UP]","Finished FeaturePlot for ------>",i,CD4T.TEC.up$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/FeaturePlot/UP/",CD4T.TEC.up$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}

########################################################################################################
## DOWN

if(FALSE){
  CD4T.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 2)
  colnames(CD4T.TEC.down)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD4T.TEC.down$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD4T.TEC.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD4T.TEC.down.exp,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.down.exp.xlsx",sheetName = "raw")
}

CD4T.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD4T.TEC.down)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD4T.TEC.down$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group2)
expr <- expr[expr$CellType == "TEC",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD4T.TEC.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD4T.TEC.down.exp,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.norm.exp.xlsx",sheetName = "Down",append = TRUE)
  
  normDt <- CD4T.TEC.down.exp[rowSums(CD4T.TEC.down.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/",
                                       ncol(scaleDt),".TEC.down.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD4T.TEC.down.exp.dt <- data.frame(t(CD4T.TEC.down.exp))
  CD4T.TEC.down.exp.dt$Patient <- rownames(CD4T.TEC.down.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD4T.TEC.down.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD4T.TEC.down.exp.dt,x = "Patient",y=colnames(CD4T.TEC.down.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD4T.TEC.down.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/barplot/Down/",
                  colnames(CD4T.TEC.down.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group2)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD4+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("H62","H75","H41","HCC06T","HCC01T","HCC08T","HCC10T","H21","H72","HCC09T","HCC07T","HCC5","H70",
            "HCC2","HCC04T","H38","HCC9","HCC1","HCC05T","HCC02T","HCC03T")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.norm.exp_byGroup.xlsx"), sheetName = "Down",append = TRUE)
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.norm.exp_bySample.xlsx"), sheetName = "Down",append = TRUE)
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD4T = rep(c("CD4+T high","CD4+T low"), c(12,9)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.down.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.down.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.down.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.down.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/barplot/Down/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/barplot/Down/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD4T.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD4T.TEC.down)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group2",cols = pal3)

for(i in seq_along(CD4T.TEC.up$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD4T.TEC.down$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD4T.TEC.DOWN]","Finished FeaturePlot for ------>",i,CD4T.TEC.down$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/FeaturePlot/Down/",CD4T.TEC.down$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}


## >>>>>> Malignant Cell <<<<<<<

########################################################################################################
if(FALSE){
  ## UP
  CD4T.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Maligant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 1)
  colnames(CD4T.MLG.up)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD4T.MLG.up$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD4T.MLG.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD4T.MLG.up.exp,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.cell.up.exp.xlsx",sheetName = "raw")
}

## UP
if(FALSE){
  CD4T.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Maligant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 1)
  colnames(CD4T.MLG.up)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD4T.MLG.up$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD4T.MLG.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD4T.MLG.up.exp,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.cell.up.exp.xlsx",sheetName = "normalized")
}

CD4T.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD4T.MLG.up)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD4T.MLG.up$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group2)
expr <- expr[expr$CellType == "Malignant cell",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD4T.MLG.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD4T.MLG.up.exp,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.cell.norm.exp.xlsx",sheetName = "UP")
  
  normDt <- CD4T.MLG.up.exp[rowSums(CD4T.MLG.up.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/",
                                       ncol(scaleDt),".Malignant.cell.up.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD4T.MLG.up.exp.dt <- data.frame(t(CD4T.MLG.up.exp))
  CD4T.MLG.up.exp.dt$Patient <- rownames(CD4T.MLG.up.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD4T.MLG.up.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD4T.MLG.up.exp.dt,x = "Patient",y=colnames(CD4T.MLG.up.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD4T.MLG.up.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/barplot/UP/",
                  colnames(CD4T.MLG.up.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group2)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD4+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("H62","H75","H41","HCC06T","HCC01T","HCC08T","HCC10T","H21","H72","HCC09T","HCC07T","HCC5","H70",
            "HCC2","HCC04T","H38","HCC9","HCC1","HCC05T","HCC02T","HCC03T")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant_cell.norm.exp_byGroup.xlsx"), sheetName = "UP")
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant_cell.norm.exp_bySample.xlsx"), sheetName = "UP")
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD4T = rep(c("CD4+T high","CD4+T low"), c(12,9)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.up.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.up.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.up.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.up.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/barplot/UP/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/barplot/UP/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD4T.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD4T.MLG.up)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group2",cols = pal3)

for(i in seq_along(CD4T.MLG.up$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD4T.MLG.up$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD4T.MLG.UP]","Finished FeaturePlot for ------>",i,CD4T.MLG.up$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/FeaturePlot/UP/",CD4T.MLG.up$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}

########################################################################################################
## DOWN
if(FALSE){
  CD4T.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Maligant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 2)
  colnames(CD4T.MLG.down)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD4T.MLG.down$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD4T.MLG.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD4T.MLG.down.exp,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.cell.down.exp.xlsx",sheetName = "raw")
}

CD4T.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD4T.MLG.down)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD4T.MLG.down$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group2)
expr <- expr[expr$CellType == "Malignant cell",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
  labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)
  CD4T.MLG.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD4T.MLG.down.exp,file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.cell.norm.exp.xlsx",sheetName = "Down",append = T)
  
  normDt <- CD4T.MLG.down.exp[rowSums(CD4T.MLG.down.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/",
                                       ncol(scaleDt),".Malignant.cell.down.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD4T.MLG.down.exp.dt <- data.frame(t(CD4T.MLG.down.exp))
  CD4T.MLG.down.exp.dt$Patient <- rownames(CD4T.MLG.down.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD4T.MLG.down.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD4T.MLG.down.exp.dt,x = "Patient",y=colnames(CD4T.MLG.down.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD4T.MLG.down.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/barplot/Down/",
                  colnames(CD4T.MLG.down.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group2)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD4+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("H62","H75","H41","HCC06T","HCC01T","HCC08T","HCC10T","H21","H72","HCC09T","HCC07T","HCC5","H70",
            "HCC2","HCC04T","H38","HCC9","HCC1","HCC05T","HCC02T","HCC03T")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant_cell.norm.exp_byGroup.xlsx"), sheetName = "Down",append = TRUE)
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant_cell.norm.exp_bySample.xlsx"), sheetName = "Down",append = TRUE)
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD4T = rep(c("CD4+T high","CD4+T low"), c(12,9)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.down.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.down.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.down.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/",
                                         ncol(scaleDt),".Malignant_cell.down.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/barplot/Down/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/barplot/Down/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD4T.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD4T.MLG.down)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group2",cols = pal3)

for(i in seq_along(CD4T.MLG.down$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD4T.MLG.down$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD4T.MLG.DOWN]","Finished FeaturePlot for ------>",i, CD4T.MLG.down$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/FeaturePlot/Down/",CD4T.MLG.down$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}

## >>>>>>>>>>>>>>>>>>>>>> CD8+T  <<<<<<<<<<<<<<<<<<<<<<<<<<<

## >>>>>> CAF <<<<<<<

########################################################################################################

if(FALSE){
  ## UP
  CD8T.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 1)
  colnames(CD8T.CAF.up)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD8T.CAF.up$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T","HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T"
  labs <- c(5,21,18,19,20,17,10,4,16,14,1,13,8,9,2,15,12,6,7,3,11)
  CD8T.CAF.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD8T.CAF.up.exp,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.up.exp.xlsx",sheetName = "raw")
}

## UP
CD8T.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD8T.CAF.up)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD8T.CAF.up$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group3)
expr <- expr[expr$CellType == "CAF",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T","HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T"
  labs <- c(5,21,18,19,20,17,10,4,16,14,1,13,8,9,2,15,12,6,7,3,11)
  CD8T.CAF.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD8T.CAF.up.exp,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.norm.exp.xlsx",sheetName = "UP")
  
  
  normDt <- CD8T.CAF.up.exp[rowSums(CD8T.CAF.up.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/",
                                       ncol(scaleDt),".CAF.up.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD8T.CAF.up.exp.dt <- data.frame(t(CD8T.CAF.up.exp))
  CD8T.CAF.up.exp.dt$Patient <- rownames(CD8T.CAF.up.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD8T.CAF.up.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD8T.CAF.up.exp.dt,x = "Patient",y=colnames(CD8T.CAF.up.exp.dt)[i],fill = "Patient") + 
      NoLegend() + labs(y = "expression") +
      ggtitle(colnames(CD8T.CAF.up.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/barplot/UP/",
                  colnames(CD8T.CAF.up.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group3)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD8+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T",
            "HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.norm.exp_byGroup.xlsx"), sheetName = "UP")
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.norm.exp_bySample.xlsx"), sheetName = "UP")
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD8T = rep(c("CD8+T high","CD8+T low"), c(11,10)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.up.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.up.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.up.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.up.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/barplot/UP/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/barplot/UP/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD8T.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD8T.CAF.up)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group3",cols = pal3)

for(i in seq_along(CD8T.CAF.up$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD8T.CAF.up$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD8T.CAF.UP]","Finished FeaturePlot for ------>",i,CD8T.CAF.up$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/FeaturePlot/UP/",CD8T.CAF.up$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}

#######################################################################################################

## DOWN

if(FALSE){
  CD8T.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 2)
  colnames(CD8T.CAF.down)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD8T.CAF.down$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T","HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T"
  labs <- c(5,21,18,19,20,17,10,4,16,14,1,13,8,9,2,15,12,6,7,3,11)
  CD8T.CAF.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD8T.CAF.down.exp,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.down.exp.xlsx",sheetName = "raw")
}

CD8T.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD8T.CAF.down)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD8T.CAF.down$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group3)
expr <- expr[expr$CellType == "CAF",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T","HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T"
  labs <- c(5,21,18,19,20,17,10,4,16,14,1,13,8,9,2,15,12,6,7,3,11)
  CD8T.CAF.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD8T.CAF.down.exp,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.norm.exp.xlsx",sheetName = "Down",append = TRUE)
  
  normDt <- CD8T.CAF.down.exp[rowSums(CD8T.CAF.down.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/",
                                       ncol(scaleDt),".CAF.down.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD8T.CAF.down.exp.dt <- data.frame(t(CD8T.CAF.down.exp))
  CD8T.CAF.down.exp.dt$Patient <- rownames(CD8T.CAF.down.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD8T.CAF.down.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD8T.CAF.down.exp.dt,x = "Patient",y=colnames(CD8T.CAF.down.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD8T.CAF.down.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/barplot/Down/",
                  colnames(CD8T.CAF.down.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group3)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD8+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T",
            "HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.norm.exp_byGroup.xlsx"),sheetName="Down",append=TRUE)
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.norm.exp_bySample.xlsx"),sheetName="Down",append=TRUE)
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD8T = rep(c("CD8+T high","CD8+T low"), c(11,10)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.down.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.down.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.down.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/",
                                         ncol(scaleDt),".CAF.down.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/barplot/Down/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/barplot/Down/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD8T.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD8T.CAF.down)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group3",cols = pal3)

for(i in seq_along(CD8T.CAF.down$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD8T.CAF.down$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD8T.CAF.DOWN]","Finished FeaturePlot for ------>",i,CD8T.CAF.down$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/FeaturePlot/Down/",CD8T.CAF.down$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}

## >>>>>> TEC <<<<<<<

########################################################################################################

if(FALSE){
  ## UP
  CD8T.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 1)
  colnames(CD8T.TEC.up)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD8T.TEC.up$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T","HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T"
  labs <- c(5,21,18,19,20,17,10,4,16,14,1,13,8,9,2,15,12,6,7,3,11)
  CD8T.TEC.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD8T.TEC.up.exp,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.up.exp.xlsx",sheetName = "raw")
}

## UP
CD8T.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD8T.TEC.up)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD8T.TEC.up$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group3)
expr <- expr[expr$CellType == "TEC",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T","HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T"
  labs <- c(5,21,18,19,20,17,10,4,16,14,1,13,8,9,2,15,12,6,7,3,11)
  CD8T.TEC.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD8T.TEC.up.exp,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.norm.exp.xlsx",sheetName = "UP")
  
  normDt <- CD8T.TEC.up.exp[rowSums(CD8T.TEC.up.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/",
                                       ncol(scaleDt),".TEC.up.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD8T.TEC.up.exp.dt <- data.frame(t(CD8T.TEC.up.exp))
  CD8T.TEC.up.exp.dt$Patient <- rownames(CD8T.TEC.up.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD8T.TEC.up.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD8T.TEC.up.exp.dt,x = "Patient",y=colnames(CD8T.TEC.up.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD8T.TEC.up.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/barplot/UP/",
                  colnames(CD8T.TEC.up.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group3)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD8+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T",
            "HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.norm.exp_byGroup.xlsx"), sheetName = "UP")
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.norm.exp_bySample.xlsx"), sheetName = "UP")
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD8T = rep(c("CD8+T high","CD8+T low"), c(11,10)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.up.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.up.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.up.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.up.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/barplot/UP/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/barplot/UP/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD8T.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD8T.TEC.up)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group3",cols = pal3)

for(i in seq_along(CD8T.TEC.up$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD8T.TEC.up$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD8T.TEC.UP]","Finished FeaturePlot for ------>",i,CD8T.TEC.up$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/FeaturePlot/UP/",CD8T.TEC.up$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}

########################################################################################################
## DOWN
if(FALSE){
  CD8T.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 2)
  colnames(CD8T.TEC.down)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD8T.TEC.down$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T","HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T"
  labs <- c(5,21,18,19,20,17,10,4,16,14,1,13,8,9,2,15,12,6,7,3,11)
  CD8T.TEC.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD8T.TEC.down.exp,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.down.exp.xlsx",sheetName = "raw")
  
  normDt <- log2(CD8T.TEC.down.exp + 1)
  normDt <- normDt[rowSums(normDt) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  xlsx::write.xlsx(scaleDt,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.down.exp.xlsx",sheetName = "normalized",append = T)
}

CD8T.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD8T.TEC.down)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD8T.TEC.down$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group3)
expr <- expr[expr$CellType == "TEC",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T","HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T"
  labs <- c(5,21,18,19,20,17,10,4,16,14,1,13,8,9,2,15,12,6,7,3,11)
  CD8T.TEC.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD8T.TEC.down.exp,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.norm.exp.xlsx",sheetName = "Down",append = T)
  
  normDt <- CD8T.TEC.down.exp[rowSums(CD8T.TEC.down.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/",
                                       ncol(scaleDt),".TEC.down.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD8T.TEC.down.exp.dt <- data.frame(t(CD8T.TEC.down.exp))
  CD8T.TEC.down.exp.dt$Patient <- rownames(CD8T.TEC.down.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD8T.TEC.down.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD8T.TEC.down.exp.dt,x = "Patient",y=colnames(CD8T.TEC.down.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD8T.TEC.down.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/barplot/Down/",
                  colnames(CD8T.TEC.down.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group3)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD8+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T",
            "HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T")
  a <- data.frame(rep(0,nrow(tmp)),rep(0,nrow(tmp)),rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.norm.exp_byGroup.xlsx"),sheetName="Down",append=TRUE)
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.norm.exp_bySample.xlsx"),sheetName="Down",append=TRUE)
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD8T = rep(c("CD8+T high","CD8+T low"), c(11,10)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.down.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.down.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.down.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/",
                                         ncol(scaleDt),".TEC.down.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/barplot/Down/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/barplot/Down/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD8T.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD8T.TEC.down)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group3",cols = pal3)

for(i in seq_along(CD8T.TEC.down$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD8T.TEC.down$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD8T.TEC.DOWN]","Finished FeaturePlot for ------>",i,CD8T.TEC.down$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/FeaturePlot/Down/",CD8T.TEC.down$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}

########################################################################################################

## >>>>>> Malignant Cell <<<<<<<

if(FALSE){
  ## UP
  CD8T.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Maligant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 1)
  colnames(CD8T.MLG.up)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD8T.MLG.up$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T","HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T"
  labs <- c(5,21,18,19,20,17,10,4,16,14,1,13,8,9,2,15,12,6,7,3,11)
  CD8T.MLG.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD8T.MLG.up.exp,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.cell.up.exp.xlsx",sheetName = "raw")
}

## UP
CD8T.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD8T.MLG.up)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD8T.MLG.up$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group3)
expr <- expr[expr$CellType == "Malignant cell",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T","HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T"
  labs <- c(5,21,18,19,20,17,10,4,16,14,1,13,8,9,2,15,12,6,7,3,11)
  CD8T.MLG.up.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD8T.MLG.up.exp,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.cell.norm.exp.xlsx",sheetName = "UP")
  
  normDt <- CD8T.MLG.up.exp[rowSums(CD8T.MLG.up.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/",
                                       ncol(scaleDt),".Malignant.cell.up.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD8T.MLG.up.exp.dt <- data.frame(t(CD8T.MLG.up.exp))
  CD8T.MLG.up.exp.dt$Patient <- rownames(CD8T.MLG.up.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD8T.MLG.up.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD8T.MLG.up.exp.dt,x = "Patient",y=colnames(CD8T.MLG.up.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD8T.MLG.up.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/barplot/UP/",
                  colnames(CD8T.MLG.up.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group3)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD8+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T",
            "HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T")
  a <- data.frame(rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/Malignant_Cell.norm.exp_byGroup.xlsx"), sheetName = "UP")
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/Malignant_Cell.norm.exp_bySample.xlsx"), sheetName = "UP")
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD8T = rep(c("CD8+T high","CD8+T low"), c(11,10)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/",
                                         ncol(scaleDt),".Malignant_Cell.up.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/",
                                         ncol(scaleDt),".Malignant_Cell.up.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/",
                                         ncol(scaleDt),".Malignant_Cell.up.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/",
                                         ncol(scaleDt),".Malignant_Cell.up.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/barplot/UP/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/barplot/UP/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD8T.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD8T.MLG.up)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group3",cols = pal3)

for(i in seq_along(CD8T.MLG.up$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD8T.MLG.up$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD8T.MLG.UP]","Finished FeaturePlot for ------>",i,CD8T.MLG.up$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/FeaturePlot/UP/",CD8T.MLG.up$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}

########################################################################################################

## DOWN
if(FALSE){
  CD8T.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Maligant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 2)
  colnames(CD8T.MLG.down)[1] <- "Gene"
  expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "counts")[CD8T.MLG.down$Gene,]
  expr <- data.frame(expr)
  # replace . to - to match barcodes in metadata
  colnames(expr) <- gsub("\\.","-",colnames(expr))
  # obtain metadata
  metadata <- tumor_sub@meta.data
  # convert expression matrix to dataFrame format
  expr <- data.frame(t(expr),stringsAsFactors = F)
  # map patient info in metadata to expression matrix
  expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T","HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T"
  labs <- c(5,21,18,19,20,17,10,4,16,14,1,13,8,9,2,15,12,6,7,3,11)
  CD8T.MLG.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD8T.MLG.down.exp,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.cell.down.exp.xlsx",sheetName = "raw")
}

CD8T.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD8T.MLG.down)[1] <- "Gene"
expr <- GetAssayData(tumor_sub,assay = "RNA",slot = "data")[CD8T.MLG.down$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor_sub@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
expr$CellType <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$CellType)
expr$Group <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Group3)
expr <- expr[expr$CellType == "Malignant cell",]
expr$CellType = NULL

if(F){
  # calculate mean expression values of DE gens in each patient
  patient <- NULL
  Patient <- list()
  for(i in 1:(length(colnames(expr))-1)){
    for(j in seq_along(unique(names(table(expr$Sample))))){
      patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
      Patient[[i]] <- patient
    }
  }
  # add gene info
  names(Patient) <- colnames(expr)[-ncol(expr)]
  # add patient info
  tmp <- data.frame(Patient)
  rownames(tmp) <- unique(expr$Sample)
  # final expression matrix
  tmp <- data.frame(t(tmp))
  
  # "HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T","HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T"
  labs <- c(5,21,18,19,20,17,10,4,16,14,1,13,8,9,2,15,12,6,7,3,11)
  CD8T.MLG.down.exp <- tmp[,labs]
  
  xlsx::write.xlsx(CD8T.MLG.down.exp,file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.cell.norm.exp.xlsx",sheetName = "Down",append = T)
  
  normDt <- CD8T.MLG.down.exp[rowSums(CD8T.MLG.down.exp) > 0,]
  scaleDt <- apply(normDt,1,scale)
  rownames(scaleDt) <- colnames(normDt)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 3) scaleDt[scaleDt > 3] = 3
  if(min(scaleDt) < -3) scaleDt[scaleDt < -3] = -3
  
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  
  pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                     filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/",
                                       ncol(scaleDt),".Malignant.cell.down.genes.heatmap.pdf"))
  
  ### barplot
  library(ggpubr)
  CD8T.MLG.down.exp.dt <- data.frame(t(CD8T.MLG.down.exp))
  CD8T.MLG.down.exp.dt$Patient <- rownames(CD8T.MLG.down.exp.dt)
  plt <- list()
  for(i in 1:(ncol(CD8T.MLG.down.exp.dt)-1)){
    plt[[i]] <- ggbarplot(CD8T.MLG.down.exp.dt,x = "Patient",y=colnames(CD8T.MLG.down.exp.dt)[i],
                          fill = "Patient") + NoLegend() +
      labs(y = "expression") +
      ggtitle(colnames(CD8T.MLG.down.exp.dt)[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5))
    ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/barplot/Down/",
                  colnames(CD8T.MLG.down.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
  }
}

DeExpStats <- function(x, option){
  
  patient <- NULL
  for(i in 1:(length(colnames(x)) - 2)){
    patient[i] <- data.frame(tapply(x[,i], x$Sample, sum))
    names(patient[i]) <- colnames(x)[i]
  }
  names(patient) <- colnames(x)[-c(ncol(x),ncol(x)-1)]
  patient <- data.frame(patient)
  rownames(patient) <- names(tapply(x[,1], x$Sample, sum))
  
  if(option == "byGroup"){
    #### Option I: Normalized gene expression based on Groups ####
    patient$Group <- plyr::mapvalues(rownames(patient), from = metadata$Sample, to = metadata$Group3)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        if(patient$Group[j] == "CD8+T high"){
          patient[j,i] <- patient[j,i]/table(x$Group)[1]
        }else {
          patient[j,i] <- patient[j,i]/table(x$Group)[2]
        }
      }
    }
    patient$Group = NULL
  }else if(option == "bySample"){
    ### normalized gene expression based on patient samples #####
    Patient_sample <- data.frame(table(x$Sample))
    patient$sample <- plyr::mapvalues(rownames(patient), from = Patient_sample$Var1, to = Patient_sample$Freq)
    patient <- data.frame(patient)
    patient$sample <- as.numeric(patient$sample)
    for(i in 1:(ncol(patient)-1)){
      for(j in 1:nrow(patient)){
        patient[j,i] <- patient[j,i]/patient$sample[j]
      }
    }
    patient$sample <- NULL
  }
  tmp <- data.frame(t(patient))
  labs <- c("HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T",
            "HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T")
  a <- data.frame(rep(0,nrow(tmp)),row.names = rownames(tmp))
  colnames(a) <- setdiff(labs,names(table(x$Sample)))
  tmp <- cbind(tmp,a)
  
  final.exp <- tmp[,labs]
  final.exp <- final.exp[rowSums(final.exp) > 0, ]
  final.exp <- round(final.exp, digits = 3)
  if(option == "byGroup"){
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/Malignant_Cell.norm.exp_byGroup.xlsx"), sheetName = "Down",append = T)
  }else {
    xlsx::write.xlsx(final.exp, file = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/Malignant_Cell.norm.exp_bySample.xlsx"), sheetName = "Down",append = T)
  }
  
  scaleDt <- apply(final.exp,1,scale)
  rownames(scaleDt) <- colnames(final.exp)
  scaleDt <- t(scaleDt)
  
  if(max(scaleDt) > 2) scaleDt[scaleDt > 2] = 2
  if(min(scaleDt) < -1) scaleDt[scaleDt < -1] = -1
  scaleDt <- t(scaleDt)
  
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(10,'RdBu')))(256)
  anno_col <- data.frame(CD8T = rep(c("CD8+T high","CD8+T low"), c(11,10)),row.names = rownames(scaleDt))
  
  if(option == "byGroup"){
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/",
                                         ncol(scaleDt),".Malignant_Cell.down.genes.heatmap.byGroup_cluster.pdf"))
    
    # unclustered
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/",
                                         ncol(scaleDt),".Malignant_Cell.down.genes.heatmap.byGroup_uncluster.pdf"))
  }else {
    # cluster groups
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/",
                                         ncol(scaleDt),".Malignant_Cell.down.genes.heatmap.bySample_cluster.pdf"))
    
    # uncluster groups (sort by infiltration ratio)
    pheatmap::pheatmap(scaleDt,border_color = F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,
                       show_colnames = F,cluster_rows = F,annotation_row = anno_col,
                       annotation_names_row = F,annotation_colors = ,gaps_row = 12,
                       filename = paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/",
                                         ncol(scaleDt),".Malignant_Cell.down.genes.heatmap.bySample_uncluster.pdf"))
  }
  
  final.exp.dt <- data.frame(t(final.exp))
  final.exp.dt$Patient <- rownames(final.exp.dt)
  plt <- list()
  
  if(option == "byGroup"){
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/barplot/Down/byGroup/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }else {
    for(i in 1:(ncol(final.exp.dt)-1)){
      plt[[i]] <- ggbarplot(final.exp.dt,x = "Patient",y=colnames(final.exp.dt)[i],
                            fill = "Patient") + NoLegend() +
        labs(y = "expression") +
        ggtitle(colnames(final.exp.dt)[i]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              axis.title = element_text(size = 14),
              plot.title = element_text(hjust = 0.5))
      ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/barplot/Down/bySample/",
                    colnames(final.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
    }
  }
}

# byGroup
DeExpStats(x = expr, option = "byGroup")
# bySample
DeExpStats(x = expr, option = "bySample")

# FeaturePlot
CD8T.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD8T.MLG.down)[1] <- "Gene"

pal1 <- c('#2A3C26','#E35822','#644422','#8DB700','#892C16','#DBD200','#B3446C','#F6A600','#604E97','#F99378','#0066A5',
          '#E68FAC','#008956','#F3F3F4','#C3B381','#BF0032','#A0C8EF','#F38500','#875691','#222222','#828281')
pal2 <- c("#BF0032","#008956","#0066A5")
pal3 <- c("#BF0032","#0066A5")
p1 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Sample", cols = pal1, pt.size = 0.8)
p2 <- DimPlot(tumor_sub,reduction = "umap",group.by = "CellType",cols = pal2,pt.size = 0.8)
p3 <- DimPlot(tumor_sub,reduction = "umap",group.by = "Group3",cols = pal3)

for(i in seq_along(CD8T.MLG.down$Gene)){
  p4 <- FeaturePlot(tumor_sub,features = CD8T.MLG.down$Gene[i], reduction = "umap", pt.size = 0.9,cols = c("#F3F3F4","#BF0032"))
  message(paste(Sys.time(),"[CD8T.MLG.DOWN]","Finished FeaturePlot for ------>",i,CD8T.MLG.down$Gene[i]))
  pc <- p4 | p3 | p1 | p2
  ggsave(paste0("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/FeaturePlot/Down/",CD8T.MLG.down$Gene[i],".pdf"),
         pc, width = 18, height = 4)
}


####===================================================================================

tumor <- readRDS("rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR_2.rds")

## >>>>>>>>>>>>>>>>> CD3+T UP <<<<<<<<<<<<<<<<<<<<<

CD3T.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/UP/CD3T/CD3.high_vs_low.UP.xlsx")
colnames(CD3T.up)[1] <- "Gene"
expr <- GetAssayData(tumor,assay = "RNA",slot = "counts")[CD3T.up$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
# calculate mean expression values of DE gens in each patient
patient <- NULL
Patient <- list()
for(i in 1:(length(colnames(expr))-1)){
  for(j in seq_along(unique(names(table(expr$Sample))))){
    patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
    Patient[[i]] <- patient
  }
}
# add gene info
names(Patient) <- colnames(expr)[-534]
# add patient info
tmp <- data.frame(Patient)
rownames(tmp) <- unique(expr$Sample)
# final expression matrix
tmp <- data.frame(t(tmp))

# "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)

CD3.up.exp <- tmp[,labs]

xlsx::write.xlsx(CD3.up.exp,file = "result/without_PVTT_MLN/tumor/UP/CD3T/CD3T.high_vs_CD3T.low.UP.geneExp.xlsx",sheetName = "raw")

normDt <- log2(CD3.up.exp + 1)
normDt <- normDt[rowSums(normDt) > 0,]
scaleDt <- apply(normDt,1,scale)
rownames(scaleDt) <- colnames(normDt)
scaleDt <- t(scaleDt)
xlsx::write.xlsx(scaleDt,file = "result/without_PVTT_MLN/tumor/UP/CD3T/CD3T.high_vs_CD3T.low.UP.geneExp.xlsx",sheetName = "normalized",append = T)

scaleDt[scaleDt > 3] = 3
scaleDt[scaleDt < -3 ] = -3
scaleDt <- t(scaleDt)

pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                   filename = "result/without_PVTT_MLN/tumor/Down/CD3T/679.UP.genes.heatmap.png")

### barplot
library(ggpubr)
# CD3+T high vs CD3+T low DOWN

CD3.up.exp.dt <- data.frame(t(CD3.up.exp))
CD3.up.exp.dt$Patient <- rownames(CD3.up.exp.dt)
plt <- list()
for(i in 1:(ncol(CD3.up.exp.dt)-1)){
  # CD3.up.exp.dt <- CD3.up.exp.dt[order(CD3.up.exp.dt[,i],decreasing = T),]
  plt[[i]] <- ggbarplot(CD3.up.exp.dt,x = "Patient",y=colnames(CD3.up.exp.dt)[i],
                        fill = "Patient") + NoLegend() +
    labs(y = "expression") +
    ggtitle(colnames(CD3.up.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0("result/without_PVTT_MLN/tumor/UP/CD3T/barplot/",
                colnames(CD3.up.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}


## >>>>>>>>>>>>>>>>> CD3+T Down <<<<<<<<<<<<<<<<<<<<<

if(F){
  CD3T.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",sheet = 2)
  xlsx::write.xlsx(CD3T.down,file = "result/without_PVTT_MLN/tumor/Down/CD3T/CD3.high_vs_low.Down.xlsx")
}

CD3T.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/Down/CD3T/CD3.high_vs_low.Down.xlsx")
colnames(CD3T.down)[1] <- "Gene"
expr <- GetAssayData(tumor,assay = "RNA",slot = "counts")[CD3T.down$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
# calculate mean expression values of DE gens in each patient
patient <- NULL
Patient <- list()
for(i in 1:(length(colnames(expr))-1)){
  for(j in seq_along(unique(names(table(expr$Sample))))){
    patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
    Patient[[i]] <- patient
  }
}
# add gene info
names(Patient) <- colnames(expr)[-534]
# add patient info
tmp <- data.frame(Patient)
rownames(tmp) <- unique(expr$Sample)
# final expression matrix
tmp <- data.frame(t(tmp))

# "H62","H75","H41","HCC06T","H72","HCC01T","HCC5","HCC10T","HCC08T","H21","HCC09T","HCC07T","H70","HCC2","HCC04T","HCC02T","H38","HCC005T","HCC9","HCC03T","HCC1"
labs <- c(20,19,21,16,18,7,5,14,8,1,13,15,17,4,11,10,2,9,6,12,3)

CD3.down.exp <- tmp[,labs]

xlsx::write.xlsx(CD3.down.exp,file = "result/without_PVTT_MLN/tumor/Down/CD3T/CD3T.high_vs_CD3T.low.Down.geneExp.xlsx",sheetName = "raw")

normDt <- log2(CD3.down.exp + 1)
normDt <- normDt[rowSums(normDt) > 0,]
scaleDt <- apply(normDt,1,scale)
rownames(scaleDt) <- colnames(normDt)
scaleDt <- t(scaleDt)
xlsx::write.xlsx(scaleDt,file = "result/without_PVTT_MLN/tumor/Down/CD3T/CD3T.high_vs_CD3T.low.Down.geneExp.xlsx",sheetName = "normalized",append = T)

scaleDt[scaleDt > 3] = 3
scaleDt[scaleDt < -3 ] = -3
scaleDt <- t(scaleDt)

pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
         filename = "result/without_PVTT_MLN/tumor/Down/CD3T/533.Down.genes.heatmap.png")

### barplot
library(ggpubr)
# CD3+T high vs CD3+T low DOWN

CD3.down.exp.dt <- data.frame(t(CD3.down.exp))
CD3.down.exp.dt$Patient <- rownames(CD3.down.exp.dt)
plt <- list()
for(i in 1:(ncol(CD3.down.exp.dt)-1)){
  # CD3.down.exp.dt <- CD3.down.exp.dt[order(CD3.down.exp.dt[,i],decreasing = T),]
  plt[[i]] <- ggbarplot(CD3.down.exp.dt,x = "Patient",y=colnames(CD3.down.exp.dt)[i],
                        fill = "Patient") + NoLegend() +
    labs(y = "expression") +
    ggtitle(colnames(CD3.down.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0("result/without_PVTT_MLN/tumor/Down/CD3T/barplot/",
                colnames(CD3.down.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}


## >>>>>>>>>>>>>>>>>>>>>> CD4+T  <<<<<<<<<<<<<<<<<<<<<<<<<<<

## >>>>>>> CD4+T UP <<<<<<<<<<<<<<<<

CD4T.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/Up/CD4T/CD4.high_vs_low.UP.xlsx")
colnames(CD4T.up)[1] <- "Gene"
expr <- GetAssayData(tumor,assay = "RNA",slot = "counts")[CD4T.up$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
# calculate mean expression values of DE gens in each patient
patient <- NULL
Patient <- list()
for(i in 1:(length(colnames(expr))-1)){
  for(j in seq_along(unique(names(table(expr$Sample))))){
    patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
    Patient[[i]] <- patient
  }
}
# add gene info
names(Patient) <- colnames(expr)[-ncol(expr)]
# add patient info
tmp <- data.frame(Patient)
rownames(tmp) <- unique(expr$Sample)
# final expression matrix
tmp <- data.frame(t(tmp))

# "H62","H75","H41","HCC06T","HCC01T","HCC08T","HCC10T","H21","H72","HCC09T","HCC07T","HCC5","H70","HCC2","HCC04T","H38","HCC9","HCC1","HCC05T","HCC02T","HCC03T"
labs <- c(20,19,21,16,7,8,14,1,18,13,15,5,17,4,11,2,6,3,9,10,12)
CD4.up.exp <- tmp[,labs]

xlsx::write.xlsx(CD4.up.exp,file = "result/without_PVTT_MLN/tumor/UP/CD4T/CD4T.high_vs_CD4T.low.UP.geneExp.xlsx",sheetName = "raw")

normDt <- log2(CD4.up.exp + 1)
normDt <- normDt[rowSums(normDt) > 0,]
scaleDt <- apply(normDt,1,scale)
rownames(scaleDt) <- colnames(normDt)
scaleDt <- t(scaleDt)
xlsx::write.xlsx(scaleDt,file = "result/without_PVTT_MLN/tumor/UP/CD4T/CD4T.high_vs_CD4T.low.UP.geneExp.xlsx",sheetName = "normalized",append = T)

scaleDt[scaleDt > 3] = 3
scaleDt[scaleDt < -3 ] = -3
scaleDt <- t(scaleDt)

pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                   filename = "result/without_PVTT_MLN/tumor/UP/CD4T/537.UP.genes.heatmap.png")

### barplot

CD4.up.exp.dt <- data.frame(t(CD4.up.exp))
CD4.up.exp.dt$Patient <- rownames(CD4.up.exp.dt)
plt <- list()
for(i in 1:(ncol(CD4.up.exp.dt)-1)){
  # CD4.up.exp.dt <- CD4.up.exp.dt[order(CD4.up.exp.dt[,i],decreasing = T),]
  plt[[i]] <- ggbarplot(CD4.up.exp.dt,x = "Patient",y=colnames(CD4.up.exp.dt)[i],
                        fill = "Patient") + NoLegend() +
    labs(y = "expression") +
    ggtitle(colnames(CD4.up.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0("result/without_PVTT_MLN/tumor/UP/CD4T/barplot/",
                colnames(CD4.up.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}


## >>>>>>>>>>>> CD4+T Down <<<<<<<<<<<<<<<<<<<<

CD4T.DEgenes <- readxl::read_excel("result/without_PVTT_MLN/tumor/Down/CD4T/CD4.high_vs_low.Down.xlsx")
colnames(CD4T.DEgenes)[1] <- "Gene"
expr <- GetAssayData(tumor,assay = "RNA",slot = "counts")[CD4T.DEgenes$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
# calculate mean expression values of DE gens in each patient
patient <- NULL
Patient <- list()
for(i in 1:(length(colnames(expr))-1)){
  for(j in seq_along(unique(names(table(expr$Sample))))){
    patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
    Patient[[i]] <- patient
  }
}
# add gene info
names(Patient) <- colnames(expr)[-534]
# add patient info
tmp <- data.frame(Patient)
rownames(tmp) <- unique(expr$Sample)
# final expression matrix
tmp <- data.frame(t(tmp))

# "H62","H75","H41","HCC06T","HCC01T","HCC08T","HCC10T","H21","H72","HCC09T","HCC07T","HCC5","H70","HCC2","HCC04T","H38","HCC9","HCC1","HCC05T","HCC02T","HCC03T"
labs <- c(20,19,21,16,7,8,14,1,18,13,15,5,17,4,11,2,6,3,9,10,12)
CD4.down.exp <- tmp[,labs]

xlsx::write.xlsx(CD4.down.exp,file = "result/without_PVTT_MLN/tumor/Down/CD4T/CD4T.high_vs_CD4T.low.Down.geneExp.xlsx",sheetName = "raw")

normDt <- log2(CD4.down.exp + 1)
normDt <- normDt[rowSums(normDt) > 0,]
scaleDt <- apply(normDt,1,scale)
rownames(scaleDt) <- colnames(normDt)
scaleDt <- t(scaleDt)
xlsx::write.xlsx(scaleDt,file = "result/without_PVTT_MLN/tumor/Down/CD4T/CD4T.high_vs_CD4T.low.Down.geneExp.xlsx",sheetName = "normalized",append = T)

scaleDt[scaleDt > 3] = 3
scaleDt[scaleDt < -3 ] = -3
scaleDt <- t(scaleDt)
pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
         filename = "result/without_PVTT_MLN/tumor/Down/CD4T/533.Down.genes.heatmap.png")

### barplot
library(ggpubr)

CD4.down.exp.dt <- data.frame(t(CD4.down.exp))
CD4.down.exp.dt$Patient <- rownames(CD4.down.exp.dt)
plt <- list()
for(i in 1:(ncol(CD4.down.exp.dt)-1)){
  # CD4.down.exp.dt <- CD4.down.exp.dt[order(CD4.down.exp.dt[,i],decreasing = T),]
  plt[[i]] <- ggbarplot(CD4.down.exp.dt, x = "Patient", y = colnames(CD4.down.exp.dt)[i],
                        fill = "Patient") + NoLegend() +
    labs(y = "expression") +
    ggtitle(colnames(CD4.down.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0("result/without_PVTT_MLN/tumor/Down/CD4T/barplot/",
                colnames(CD4.down.exp.dt)[i],".pdf"),plt[[i]], width = 8, height = 3)
}


## >>>>>>>>>>>>>>>>>>>>>> CD8+T  <<<<<<<<<<<<<<<<<<<<<<<<<<<

## >>>>>>>>> CD8+T Up  <<<<<<<<

CD8T.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/Up/CD8T/CD8.high_vs_low.UP.xlsx")
colnames(CD8T.up)[1] <- "Gene"
expr <- GetAssayData(tumor,assay = "RNA",slot = "counts")[CD8T.up$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
# calculate mean expression values of DE gens in each patient
patient <- NULL
Patient <- list()
for(i in 1:(length(colnames(expr))-1)){
  for(j in seq_along(unique(names(table(expr$Sample))))){
    patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
    Patient[[i]] <- patient
  }
}
# add gene info
names(Patient) <- colnames(expr)[-ncol(expr)]
# add patient info
tmp <- data.frame(Patient)
rownames(tmp) <- unique(expr$Sample)
# final expression matrix
tmp <- data.frame(t(tmp))

# "HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T","HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T"
labs <- c(5,21,18,19,20,17,10,4,16,14,1,13,8,9,2,15,12,6,7,3,11)
CD8.up.exp <- tmp[,labs]

xlsx::write.xlsx(CD8.up.exp,file = "result/without_PVTT_MLN/tumor/UP/CD8T/CD8T.high_vs_CD8T.low.UP.geneExp.xlsx",sheetName = "raw")

normDt <- log2(CD8.up.exp + 1)
normDt <- normDt[rowSums(normDt) > 0,]
scaleDt <- apply(normDt,1,scale)
rownames(scaleDt) <- colnames(normDt)
scaleDt <- t(scaleDt)
xlsx::write.xlsx(scaleDt,file = "result/without_PVTT_MLN/tumor/UP/CD8T/CD8T.high_vs_CD8T.low.UP.geneExp.xlsx",sheetName = "normalized",append = T)

scaleDt[scaleDt > 3] = 3
scaleDt[scaleDt < -3 ] = -3
scaleDt <- t(scaleDt)

pheatmap::pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
                   filename = "result/without_PVTT_MLN/tumor/UP/CD8T/176.UP.genes.heatmap.png")

### barplot

CD8.up.exp.dt <- data.frame(t(CD8.up.exp))
CD8.up.exp.dt$Patient <- rownames(CD8.up.exp.dt)
plt <- list()
for(i in 1:(ncol(CD8.up.exp.dt)-1)){
  # CD8.up.exp.dt <- CD8.up.exp.dt[order(CD8.up.exp.dt[,i],decreasing = T),]
  plt[[i]] <- ggbarplot(CD8.up.exp.dt,x = "Patient",y=colnames(CD8.up.exp.dt)[i],
                        fill = "Patient") + NoLegend() +
    labs(y = "expression") +
    ggtitle(colnames(CD8.up.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0("result/without_PVTT_MLN/tumor/UP/CD8T/barplot/",
                colnames(CD8.up.exp.dt)[i],".pdf"),plt[[i]],width = 8,height = 3)
}


## >>>>>>>>>>>>>>> CD8+T Down <<<<<<<<<<<<<<<<<<<<<<<<<

CD8T.DEgenes <- readxl::read_excel("result/without_PVTT_MLN/tumor/Down/CD8T/CD8.high_vs_low.Down.xlsx")
colnames(CD8T.DEgenes)[1] <- "Gene"
expr <- GetAssayData(tumor,assay = "RNA",slot = "counts")[CD8T.DEgenes$Gene,]
expr <- data.frame(expr)
# replace . to - to match barcodes in metadata
colnames(expr) <- gsub("\\.","-",colnames(expr))
# obtain metadata
metadata <- tumor@meta.data
# convert expression matrix to dataFrame format
expr <- data.frame(t(expr),stringsAsFactors = F)
# map patient info in metadata to expression matrix
expr$Sample <- plyr::mapvalues(rownames(expr), from = rownames(metadata), to = metadata$Sample)
# calculate mean expression values of DE gens in each patient
patient <- NULL
Patient <- list()
for(i in 1:(length(colnames(expr))-1)){
  for(j in seq_along(unique(names(table(expr$Sample))))){
    patient[j] <- mean(as.numeric(expr[1:table(expr$Sample)[j],i]))
    Patient[[i]] <- patient
  }
}
# add gene info
names(Patient) <- colnames(expr)[-379]
# add patient info
tmp <- data.frame(Patient)
rownames(tmp) <- unique(expr$Sample)
# final expression matrix
tmp <- data.frame(t(tmp))

# "HCC5","H41","H72","H75","H62","H70","HCC02T","HCC2","HCC06T","HCC10T","H21","HCC09T","HCC08T","HCC05T","H38","HCC07T","HCC03T","HCC9","HCC01T","HCC1","HCC04T"
labs <- c(5,21,18,19,20,17,10,4,16,14,1,13,8,9,2,15,12,6,7,3,11)

CD8.down.exp <- tmp[,labs]

xlsx::write.xlsx(CD8.down.exp,file = "result/without_PVTT_MLN/tumor/Down/CD8T/CD8T.high_vs_CD8T.low.Down.geneExp.xlsx",sheetName = "raw")

normDt <- log2(CD8.down.exp + 1)
normDt <- normDt[rowSums(normDt) > 0,]
scaleDt <- apply(normDt,1,scale)
rownames(scaleDt) <- colnames(normDt)
scaleDt <- t(scaleDt)
xlsx::write.xlsx(scaleDt,file = "result/without_PVTT_MLN/tumor/Down/CD8T/CD8T.high_vs_CD8T.low.Down.geneExp.xlsx",sheetName = "normalized",append = T)

scaleDt[scaleDt > 3] = 4
scaleDt[scaleDt < -3 ] = -3
scaleDt <- t(scaleDt)
pheatmap(scaleDt,border=F,color = pal,fontsize = 14,angle_col = 45,width = 20,height = 6,show_colnames = F,
         filename = "result/without_PVTT_MLN/tumor/Down/CD8T/378.Down.genes.heatmap.png")

### barplot
library(ggpubr)

# CD8+T high vs CD8+T low DOWN
CD8.down.exp.dt <- data.frame(t(CD8.down.exp))
CD8.down.exp.dt$Patient <- rownames(CD8.down.exp.dt)
plt <- list()
for(i in 1:(ncol(CD8.down.exp.dt)-1)){
  # CD8.down.exp.dt <- CD8.down.exp.dt[order(CD8.down.exp.dt[,i],decreasing = T),]
  plt[[i]] <- ggbarplot(CD8.down.exp.dt, x = "Patient", y = colnames(CD8.down.exp.dt)[i],
                        fill = "Patient") + NoLegend() +
    labs(y = "expression") +
    ggtitle(colnames(CD8.down.exp.dt)[i]) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0("result/without_PVTT_MLN/tumor/Down/CD8T/barplot/",
                colnames(CD8.down.exp.dt)[i],".pdf"),plt[[i]], width = 8, height = 3)
}
