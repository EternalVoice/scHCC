###  1.createObj_merge_into_a_list.r

setwd("/home/tongqiang/project/shengzhou/liver_cancer")

# create directories
if(!file.exists("results")) dir.create("results")
if(!file.exists("results/GSE149614")) dir.create("results/GSE149614")
if(!file.exists("results/GSE146115")) dir.create("results/GSE146115")
if(!file.exists("results/GSE151530")) dir.create("results/GSE151530")
if(!file.exists("results/GSE103867")) dir.create("results/GSE103867")
if(!file.exists("results/GSE125449_Set1")) dir.create("results/GSE125449_Set1")
if(!file.exists("results/GSE125449_Set2")) dir.create("results/GSE125449_Set2")

# obtain working directories
wkd <- "/home/tongqiang/project/shengzhou/liver_cancer/results"
samples <- list.files(wkd, "^GSE")
filePath <- sapply(samples,function(x) {paste(wkd,x, sep = "/")})

# load required pkgs
if(!requireNamespace("Seurat")) install.packages("Seurat")
library(Seurat)

# GSE149614
GSE149614_matrix <- as.data.frame(data.table::fread("GSE149614/GSE149614_HCC.scRNAseq.S71915.count.txt.gz",
                                                    stringsAsFactors = F,encoding = "UTF-8"))
GSE149614_matrix <- GSE149614_matrix[!duplicated(GSE149614_matrix[,1]),]
rownames(GSE149614_matrix) <- GSE149614_matrix[,1]
GSE149614_matrix <- GSE149614_matrix[,-1]
rownames(GSE149614_matrix) <- gsub("_", "-", rownames(GSE149614_matrix))
# create seurat object
GSE149614_obj <- CreateSeuratObject(GSE149614_matrix,min.cells = 3)
GSE149614_obj$GEO <- "GSE149614"

# GSE146115
GSE146115_matrix <- as.data.frame(data.table::fread("GSE146115/GSE146115_HCC1-2-5-9_count_with_ERCC.txt.gz",
                                             stringsAsFactors = F,encoding = "UTF-8"))
GSE146115_matrix <- GSE146115_matrix[!duplicated(GSE146115_matrix[,1]),]
rownames(GSE146115_matrix) <- GSE146115_matrix[,1]
GSE146115_matrix <- GSE146115_matrix[,-1]
rownames(GSE146115_matrix) <- gsub("_", "-", rownames(GSE146115_matrix))
# create seurat object
GSE146115_obj <- CreateSeuratObject(GSE146115_matrix,min.cells = 3)
GSE146115_obj$GEO <- "GSE146115"

# GSE151530
GSE151530_matrix <- Read10X("GSE151530/feature_bc_matrix/")
rownames(GSE151530_matrix) <- gsub("_","-",rownames(GSE151530_matrix))
GSE151530_obj <- CreateSeuratObject(GSE151530_matrix, min.cells = 3)
GSE151530_obj$GEO <- "GSE151530"

# GSE103867
GSE103867_matrix <- Read10X("GSE103867/feature_bc_matrix/")
rownames(GSE103867_matrix) <- gsub("_","-",rownames(GSE103867_matrix))
GSE103867_obj <- CreateSeuratObject(GSE103867_matrix, min.cells = 3)
GSE103867_obj$GEO <- "GSE103867"

# GSE125449_Set1
GSE125449_Set1_matrix <- Read10X("GSE125449/Set1/")
rownames(GSE125449_Set1_matrix) <- gsub("_","-",rownames(GSE125449_Set1_matrix))
GSE125449_Set1_obj <- CreateSeuratObject(GSE125449_Set1_matrix, min.cells = 3)
GSE125449_Set1_obj$GEO <- "GSE125449_Set1"

# GSE125449_Set2
GSE125449_Set2_matrix <- Read10X("GSE125449/Set2/")
rownames(GSE125449_Set2_matrix) <- gsub("_","-",rownames(GSE125449_Set2_matrix))
GSE125449_Set2_obj <- CreateSeuratObject(GSE125449_Set2_matrix, min.cells = 3)
GSE125449_Set2_obj$GEO <- "GSE125449_Set2"

# combine all seurat objects into a list

# GSE103867
# GSE125449_Set1
# GSE125449_Set2
# GSE146115
# GSE149614
# GSE151530
scRNAlist <- list(GSE103867 = GSE103867_obj,GSE125449_Set1 = GSE125449_Set1_obj,
                  GSE125449_Set2 = GSE125449_Set2_obj,GSE146115 = GSE146115_obj,
                  GSE149614 = GSE149614_obj,GSE151530 = GSE151530_obj)

### Add metadata to SeuratObject

# GSE151530
GSE151530_info <- read.table("metadata/GSE151530_Info.txt.gz",header = T,sep = "\t")
GSE151530_info$Type <- gsub("cells","cell",GSE151530_info$Type)
GSE151530_info$Type <- gsub("CAFs","CAF",GSE151530_info$Type)
GSE151530_info$Type <- gsub("TAMs","TAM",GSE151530_info$Type)
GSE151530_info$Type <- gsub("TECs","TEC",GSE151530_info$Type)

for(i in 1:nrow(GSE151530_info)){scRNAlist$GSE151530@meta.data[colnames(scRNAlist$GSE151530) == GSE151530_info$Cell[i],"Sample"] <- GSE151530_info$Sample[i]}
for(i in 1:nrow(GSE151530_info)){scRNAlist$GSE151530@meta.data[colnames(scRNAlist$GSE151530) == GSE151530_info$Cell[i],"S_ID"] <- GSE151530_info$S_ID[i]}
for(i in 1:nrow(GSE151530_info)){scRNAlist$GSE151530@meta.data[colnames(scRNAlist$GSE151530) == GSE151530_info$Cell[i],"Type"] <- GSE151530_info$Type[i]}

# GSE125449_Set1
GSE125449_Set1_info <- read.table("metadata/GSE125449_Set1_samples.txt.gz",header = T,sep = "\t")

for(i in 1:nrow(GSE125449_Set1_info)){scRNAlist$GSE125449_Set1@meta.data[colnames(scRNAlist$GSE125449_Set1) == GSE125449_Set1_info$'Cell.Barcode'[i],"Sample"] <- GSE125449_Set1_info$Sample[i]}
for(i in 1:nrow(GSE125449_Set1_info)){scRNAlist$GSE125449_Set1@meta.data[colnames(scRNAlist$GSE125449_Set1) == GSE125449_Set1_info$'Cell.Barcode'[i],"Type"] <- GSE125449_Set1_info$Type[i]}
# GSE125449_Set2
GSE125449_Set2_info <- read.table("metadata/GSE125449_Set2_samples.txt.gz",header = T,sep = "\t")
for(i in 1:nrow(GSE125449_Set2_info)){scRNAlist$GSE125449_Set2@meta.data[colnames(scRNAlist$GSE125449_Set2) == GSE125449_Set2_info$'Cell.Barcode'[i],"Sample"] <- GSE125449_Set2_info$Sample[i]}
for(i in 1:nrow(GSE125449_Set2_info)){scRNAlist$GSE125449_Set2@meta.data[colnames(scRNAlist$GSE125449_Set2) == GSE125449_Set2_info$'Cell.Barcode'[i],"Type"] <- GSE125449_Set2_info$Type[i]}


# Add clinical information to SeuratObject

# GSE125449_Set1
GSE125449_clin <- readxl::read_excel("metadata/GSE125449_clinical.xlsx")

for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set1@meta.data[which(scRNAlist$GSE125449_Set1$Sample == GSE125449_clin$Sample[i]),"sample"] <- GSE125449_clin$ID[i]}
scRNAlist$GSE125449_Set1$Sample <- scRNAlist$GSE125449_Set1$sample
scRNAlist$GSE125449_Set1$sample <- NULL

for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set1@meta.data[which(scRNAlist$GSE125449_Set1$Sample == GSE125449_clin$ID[i]),"Sex"] <- GSE125449_clin$Sex[i]}

# c("Age","Race","Diagnosis","Stage","Etiology","BiopsyTiming","Treatment","Action")
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set1@meta.data[which(scRNAlist$GSE125449_Set1$Sample == GSE125449_clin$ID[i]),'Age'] <- GSE125449_clin$Age[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set1@meta.data[which(scRNAlist$GSE125449_Set1$Sample == GSE125449_clin$ID[i]),'Race'] <- GSE125449_clin$Race[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set1@meta.data[which(scRNAlist$GSE125449_Set1$Sample == GSE125449_clin$ID[i]),'Diagnosis'] <- GSE125449_clin$Diagnosis[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set1@meta.data[which(scRNAlist$GSE125449_Set1$Sample == GSE125449_clin$ID[i]),'Stage'] <- GSE125449_clin$Stage[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set1@meta.data[which(scRNAlist$GSE125449_Set1$Sample == GSE125449_clin$ID[i]),'Etiology'] <- GSE125449_clin$Etiology[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set1@meta.data[which(scRNAlist$GSE125449_Set1$Sample == GSE125449_clin$ID[i]),'BiopsyTiming'] <- GSE125449_clin$BiopsyTiming[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set1@meta.data[which(scRNAlist$GSE125449_Set1$Sample == GSE125449_clin$ID[i]),'Treatment'] <- GSE125449_clin$Treatment[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set1@meta.data[which(scRNAlist$GSE125449_Set1$Sample == GSE125449_clin$ID[i]),'Action'] <- GSE125449_clin$Action[i]}

# GSE125449_Set2
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set2@meta.data[which(scRNAlist$GSE125449_Set2$Sample == GSE125449_clin$Sample[i]),"sample"] <- GSE125449_clin$ID[i]}
scRNAlist$GSE125449_Set2$Sample <- scRNAlist$GSE125449_Set2$sample
scRNAlist$GSE125449_Set2$sample <- NULL

# c("Sex",Age","Race","Diagnosis","Stage","Etiology","BiopsyTiming","Treatment","Action")
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set2@meta.data[which(scRNAlist$GSE125449_Set2$Sample == GSE125449_clin$ID[i]),"Sex"] <- GSE125449_clin$Sex[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set2@meta.data[which(scRNAlist$GSE125449_Set2$Sample == GSE125449_clin$ID[i]),'Age'] <- GSE125449_clin$Age[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set2@meta.data[which(scRNAlist$GSE125449_Set2$Sample == GSE125449_clin$ID[i]),'Race'] <- GSE125449_clin$Race[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set2@meta.data[which(scRNAlist$GSE125449_Set2$Sample == GSE125449_clin$ID[i]),'Diagnosis'] <- GSE125449_clin$Diagnosis[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set2@meta.data[which(scRNAlist$GSE125449_Set2$Sample == GSE125449_clin$ID[i]),'Stage'] <- GSE125449_clin$Stage[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set2@meta.data[which(scRNAlist$GSE125449_Set2$Sample == GSE125449_clin$ID[i]),'Etiology'] <- GSE125449_clin$Etiology[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set2@meta.data[which(scRNAlist$GSE125449_Set2$Sample == GSE125449_clin$ID[i]),'BiopsyTiming'] <- GSE125449_clin$BiopsyTiming[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set2@meta.data[which(scRNAlist$GSE125449_Set2$Sample == GSE125449_clin$ID[i]),'Treatment'] <- GSE125449_clin$Treatment[i]}
for(i in 1:nrow(GSE125449_clin)){scRNAlist$GSE125449_Set2@meta.data[which(scRNAlist$GSE125449_Set2$Sample == GSE125449_clin$ID[i]),'Action'] <- GSE125449_clin$Action[i]}

# rename GSE151530 sample name
levels(scRNAlist$GSE151530$Sample) <- gsub("C26a","C26",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("C26b","C26",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("C46a","C46",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("C46b","C46",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("H34a","H34",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("H34a","H34",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("H34b","H34",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("H34c","H34",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("H49a","H49",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("H49b","H49",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("H58a","H58",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("H58b","H58",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("H58c","H58",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("H68a","H68",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("H68b","H68",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("H73a","H73",levels(scRNAlist$GSE151530$Sample))
levels(scRNAlist$GSE151530$Sample) <- gsub("H73b","H73",levels(scRNAlist$GSE151530$Sample))

GSE151530_clin <- readxl::read_excel("metadata/GSE151530_clinical.xlsx")

# c("Sex","Age","Race","Diagnosis","Stage","Etiology","BiopsyTiming","Treatment","Action")
for(i in 1:nrow(GSE151530_clin)){scRNAlist$GSE151530@meta.data[which(scRNAlist$GSE151530$Sample == GSE151530_clin$ID[i]),"Sex"] <- GSE151530_clin$Sex[i]}
for(i in 1:nrow(GSE151530_clin)){scRNAlist$GSE151530@meta.data[which(scRNAlist$GSE151530$Sample == GSE151530_clin$ID[i]),"Age"] <- GSE151530_clin$Age[i]}
for(i in 1:nrow(GSE151530_clin)){scRNAlist$GSE151530@meta.data[which(scRNAlist$GSE151530$Sample == GSE151530_clin$ID[i]),"Race"] <- GSE151530_clin$Race[i]}
for(i in 1:nrow(GSE151530_clin)){scRNAlist$GSE151530@meta.data[which(scRNAlist$GSE151530$Sample == GSE151530_clin$ID[i]),"Diagnosis"] <- GSE151530_clin$Diagnosis[i]}
for(i in 1:nrow(GSE151530_clin)){scRNAlist$GSE151530@meta.data[which(scRNAlist$GSE151530$Sample == GSE151530_clin$ID[i]),"Stage"] <- GSE151530_clin$Stage[i]}
for(i in 1:nrow(GSE151530_clin)){scRNAlist$GSE151530@meta.data[which(scRNAlist$GSE151530$Sample == GSE151530_clin$ID[i]),"Etiology"] <- GSE151530_clin$Etiology[i]}
for(i in 1:nrow(GSE151530_clin)){scRNAlist$GSE151530@meta.data[which(scRNAlist$GSE151530$Sample == GSE151530_clin$ID[i]),"Treatment"] <- GSE151530_clin$Treatment[i]}
for(i in 1:nrow(GSE151530_clin)){scRNAlist$GSE151530@meta.data[which(scRNAlist$GSE151530$Sample == GSE151530_clin$ID[i]),"Action"] <- GSE151530_clin$Action[i]}

# 
scRNAlist$GSE146115$Action <- "NA"
scRNAlist$GSE146115$Diagnosis <- "HCC"
scRNAlist$GSE149614$Action <- "NA"

GSE149614_clin <- readxl::read_excel("metadata/GSE149614_clinical.xlsx")
scRNAlist$GSE149614$Sample <- scRNAlist$GSE149614$orig.ident
for(i in 1:nrow(GSE149614_clin)){
  scRNAlist$GSE149614@meta.data[which(scRNAlist$GSE149614$Sample == GSE149614_clin$ID[i]),"Diagnosis"] <- GSE149614_clin$Diagnosis[i]
}
scRNAlist$GSE103867 <- NULL
scRNAlist$GSE146115$Sample <- scRNAlist$GSE146115$orig.ident

for(i in seq_along(scRNAlist)){
  scRNAlist[[i]]$orig.ident <- scRNAlist[[i]]$Sample
  Idents(scRNAlist[[i]]) <- scRNAlist[[i]]$Sample
}

scRNAlist$GSE146115$Type <- "unclassified"
scRNAlist$GSE149614$Type <- "unclassified"

saveRDS(scRNAlist,file = paste(wkd, "rds","scRNAlist.raw.rds",sep = "/"))

rm(list = ls());gc()

## Selection criteria: 1.Normal tissue && HCC patients without drug treatment

scRNAlist <- readRDS(paste(wkd, "rds","scRNAlist.raw.rds",sep = "/"))
# Diagnosis = "HCC" & Action = "NA"
GSE125449_1_HCC_cells <- subset(scRNAlist$GSE125449_Set1@meta.data, Diagnosis=="HCC" & Action == "NA")
GSE125449_1_HCC <- subset(scRNAlist$GSE125449_Set1, cells = rownames(GSE125449_1_HCC_cells))

# Diagnosis = "HCC" & Action = "NA"
GSE125449_2_HCC_cells <- subset(scRNAlist$GSE125449_Set2@meta.data, Diagnosis=="HCC" & Action == "NA")

# Action = "NA"
scRNAlist$GSE151530$Sample <- as.character(scRNAlist$GSE151530$Sample)
GSE151530_NA_cells <- subset(scRNAlist$GSE151530@meta.data, Diagnosis=="HCC" & Action == "NA")
GSE151530_NA <- subset(scRNAlist$GSE151530, cells = rownames(GSE151530_NA_cells))

# combine all samples into a list
scRNAlist_sub <- list(
  GSE125449 = GSE125449_1_HCC,
  GSE146115 = scRNAlist$GSE146115,
  GSE149614 = scRNAlist$GSE149614,
  GSE151530 = GSE151530_NA
)

# The following analysis are all based on this data
# with PVTT,MLN
saveRDS(scRNAlist_sub, file = "rds/with_PVTT_MLN/1.scRNAlist_sub_with_PVTT_MLN.rds")
# remove PVTT, MLN
GSE149614_sub_cells <- subset(scRNAlist_sub$GSE149614@meta.data, Diagnosis=="HCC" | Diagnosis=="Normal")
GSE149614_sub <- subset(scRNAlist_sub$GSE149614, cells = rownames(GSE149614_sub_cells))
scRNAlist_sub$GSE149614 <- GSE149614_sub
saveRDS(scRNAlist_sub, file = "rds/without_PVTT_MLN/1.scRNAlist_sub_without_PVTT_MLN.rds")

##########################################################################################################