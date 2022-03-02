

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  TCGA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>> 1. Pan-cancer  <<<<<<<<<<<<<<<<<<<<<<<<<<

library(RTCGA)
library(stringr)

# Prepare TCGA gene expression matrix
if(!file.exists("refdb/TCGA_36cancers.rnaseq.RData")){
  library(RTCGA.rnaseq)
  # cohorts <- infoTCGA() %>% rownames() %>% sub("-counts", ".rnaseq", x=.)
  
  # [1] "ACC.rnaseq"      "BLCA.rnaseq"     "BRCA.rnaseq"     "CESC.rnaseq"     "CHOL.rnaseq"    "COAD.rnaseq"    
  # [7] "COADREAD.rnaseq" "DLBC.rnaseq"     "ESCA.rnaseq"     "GBM.rnaseq"      "GBMLGG.rnaseq"  "HNSC.rnaseq"
  # [13] "KICH.rnaseq"    "KIPAN.rnaseq"    "KIRC.rnaseq"     "KIRP.rnaseq"     "LAML.rnaseq"    "LGG.rnaseq"
  # [19] "LIHC.rnaseq"     "LUAD.rnaseq"    "LUSC.rnaseq"     "OV.rnaseq"       "PAAD.rnaseq"    "PCPG.rnaseq"
  # [25] "PRAD.rnaseq"     "READ.rnaseq"    "SARC.rnaseq"     "SKCM.rnaseq"     "STAD.rnaseq"    "STES.rnaseq"
  # [31] "TGCT.rnaseq"     "THCA.rnaseq"    "THYM.rnaseq"     "UCEC.rnaseq"     "UCS.rnaseq"     "UVM.rnaseq"
  
  TCGA_list <- list(
    ACC.rnaseq = ACC.rnaseq, BLCA.rnaseq = BLCA.rnaseq, BRCA.rnaseq = BRCA.rnaseq,
    CESC.rnaseq = CESC.rnaseq, CHOL.rnaseq = CHOL.rnaseq, COAD.rnaseq = COAD.rnaseq,
    COADREAD.rnaseq = COADREAD.rnaseq, DLBC.rnaseq = DLBC.rnaseq, ESCA.rnaseq = ESCA.rnaseq,
    GBM.rnaseq = GBM.rnaseq, GBMLGG.rnaseq = GBMLGG.rnaseq, HNSC.rnaseq = HNSC.rnaseq,
    KICH.rnaseq = KICH.rnaseq, KIPAN.rnaseq = KIPAN.rnaseq, KIRC.rnaseq = KIRC.rnaseq,
    KIRP.rnaseq = KIRP.rnaseq, LAML.rnaseq = LAML.rnaseq, LGG.rnaseq = LGG.rnaseq,
    LIHC.rnaseq = LIHC.rnaseq, LUAD.rnaseq = LUAD.rnaseq, LUSC.rnaseq = LUSC.rnaseq,
    OV.rnaseq = OV.rnaseq, PAAD.rnaseq = PAAD.rnaseq, PCPG.rnaseq = PCPG.rnaseq, 
    PRAD.rnaseq = PRAD.rnaseq, READ.rnaseq = READ.rnaseq, SARC.rnaseq = SARC.rnaseq,
    SKCM.rnaseq = SKCM.rnaseq, STAD.rnaseq = STAD.rnaseq, STES.rnaseq = STES.rnaseq, 
    TGCT.rnaseq = TGCT.rnaseq, THCA.rnaseq = THCA.rnaseq, THYM.rnaseq = THYM.rnaseq, 
    UCEC.rnaseq = UCEC.rnaseq, UCS.rnaseq = UCS.rnaseq, UVM.rnaseq = UVM.rnaseq
  )
  TCGA_36cancers.rnaseq <- lapply(TCGA_list, function(x){
    for(i in 2:30){colnames(x)[i] <- gsub("\\?\\|","",colnames(x)[i])}
    colnames(x) <- str_split(colnames(x),"\\|",simplify = T)[,1]
    return(x)
  })
  save(TCGA_36cancers.rnaseq, file = "refdb/TCGA_36cancers.rnaseq.RData")
}else {
  load("refdb/TCGA_36cancers.rnaseq.RData")
}

# define a pan-cacancer function
pan_cancer_val <- function(DEgenes, type, wkd){
  if(type == "UP"){
    de.genes <- readxl::read_excel(DEgenes,sheet = 1)
  }else if(type == "DOWN"){
    de.genes <- readxl::read_excel(DEgenes, sheet = 2)
  }else {
    stop("Error: Unrecognized Differential Gene Types !!!")
  }
  
  colnames(de.genes)[1] <- "Genes"
  genes <- de.genes$Genes[which(de.genes$Genes %in% colnames(TCGA_36cancers.rnaseq$ACC.rnaseq))]
  expr <- expressionsTCGA(
    TCGA_36cancers.rnaseq$ACC.rnaseq, TCGA_36cancers.rnaseq$BLCA.rnaseq,
    TCGA_36cancers.rnaseq$BRCA.rnaseq, TCGA_36cancers.rnaseq$CESC.rnaseq,
    TCGA_36cancers.rnaseq$CHOL.rnaseq, TCGA_36cancers.rnaseq$COAD.rnaseq,
    TCGA_36cancers.rnaseq$COADREAD.rnaseq, TCGA_36cancers.rnaseq$DLBC.rnaseq,
    TCGA_36cancers.rnaseq$ESCA.rnaseq, TCGA_36cancers.rnaseq$GBM.rnaseq,
    TCGA_36cancers.rnaseq$GBMLGG.rnaseq, TCGA_36cancers.rnaseq$HNSC.rnaseq,
    TCGA_36cancers.rnaseq$KICH.rnaseq, TCGA_36cancers.rnaseq$KIPAN.rnaseq,
    TCGA_36cancers.rnaseq$KIRC.rnaseq, TCGA_36cancers.rnaseq$KIRP.rnaseq,
    TCGA_36cancers.rnaseq$LAML.rnaseq, TCGA_36cancers.rnaseq$LGG.rnaseq,
    TCGA_36cancers.rnaseq$LIHC.rnaseq, TCGA_36cancers.rnaseq$LUAD.rnaseq,
    TCGA_36cancers.rnaseq$LUSC.rnaseq, TCGA_36cancers.rnaseq$OV.rnaseq,
    TCGA_36cancers.rnaseq$PAAD.rnaseq, TCGA_36cancers.rnaseq$PCPG.rnaseq,
    TCGA_36cancers.rnaseq$PRAD.rnaseq, TCGA_36cancers.rnaseq$READ.rnaseq,
    TCGA_36cancers.rnaseq$SARC.rnaseq, TCGA_36cancers.rnaseq$SKCM.rnaseq,
    TCGA_36cancers.rnaseq$STAD.rnaseq, TCGA_36cancers.rnaseq$STES.rnaseq,
    TCGA_36cancers.rnaseq$TGCT.rnaseq, TCGA_36cancers.rnaseq$THCA.rnaseq,
    TCGA_36cancers.rnaseq$THYM.rnaseq, TCGA_36cancers.rnaseq$UCEC.rnaseq,
    TCGA_36cancers.rnaseq$UCS.rnaseq, TCGA_36cancers.rnaseq$UVM.rnaseq,
    extract.cols = genes
  )
  expr$dataset <- gsub(".rnaseq","",str_split(expr$dataset, "\\$", simplify = T)[,2])
  Type <- ifelse(substr(expr$bcr_patient_barcode,14,15) == "11", "Normal", "Tumor")
  expr <- data.frame(tibble::add_column(expr, Type, .before = 2))
  tmp <- t(log2(t(expr[,4:ncol(expr)])+1))
  norm.exp <- cbind(expr[,1:3],tmp)
  colnames(norm.exp) <- gsub("\\.","-",colnames(norm.exp))
  for(item in genes){
    p <- ggpubr::ggboxplot(norm.exp, x = "dataset", y = item, title = item, ylab = "Expression",
                           color = "dataset", palette = rainbow(36),bxp.errorbar = TRUE) +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                     legend.position = "none", axis.title.x = element_blank())
    ggsave(paste0(wkd,item,".pdf"), p, device = "pdf", width = 12,height = 5)
  }
}

#======= Validation for all Differentially Expressed Genes ====================

if(FALSE){
  # CD3T Up
  pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",type = "UP",
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD3T/up/")
  # CD3T Down
  pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",type = "DOWN",
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD3T/down/")
  #CD4T Up
  pan_cancer_val("result/without_PVTT_MLN/tumor/CD4T/CD4.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "UP",
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD4T/up/")
  # CD4T Down
  pan_cancer_val("result/without_PVTT_MLN/tumor/CD4T/CD4.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "DOWN",
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD4T/down/")
  #CD8T Up
  pan_cancer_val("result/without_PVTT_MLN/tumor/CD8T/CD8.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "UP",
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD8T/up/")
  # CD8T Down
  pan_cancer_val("result/without_PVTT_MLN/tumor/CD8T/CD8.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "DOWN",
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD8T/down/")
}


#======= Validation for all Differentially Expressed Genes (CD45 removed) ====================

## >>>>>>>>>>>> CD3+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD3T/without_CD45/CAF/up/")
### CAF Down
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD3T/without_CD45/CAF/down/")

### TEC UP
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD3T/without_CD45/TEC/up/")

### TEC Down
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD3T/without_CD45/TEC/down/")

### Malignant_cell UP
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Maligant.Cell.CD3T.hi_vs_lo.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD3T/without_CD45/Malignant_cell/up/")

### Malignant_cell Down
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Maligant.Cell.CD3T.hi_vs_lo.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD3T/without_CD45/Malignant_cell/down/")


## >>>>>>>>>>>> CD4+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD4T/without_CD45/CAF/up/")
### CAF Down
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD4T/without_CD45/CAF/down/")

### TEC UP
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD4T/without_CD45/TEC/up/")

### TEC Down
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD4T/without_CD45/TEC/down/")

### Malignant_cell UP
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD4T/without_CD45/Malignant_cell/up/")

### Malignant_cell Down
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD4T/without_CD45/Malignant_cell/down/")


## >>>>>>>>>>>> CD8+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD8T/without_CD45/CAF/up/")
### CAF Down
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD8T/without_CD45/CAF/down/")

### TEC UP
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD8T/without_CD45/TEC/up/")

### TEC Down
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD8T/without_CD45/TEC/down/")

### Malignant_cell UP
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD8T/without_CD45/Malignant_cell/up/")

### Malignant_cell Down
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/CD8T/without_CD45/Malignant_cell/down/")


### >>>>>>>>>>>> OVERLAPS between CD4+T & CD8+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/CAF.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/overlap_CD4T_CD8T/CAF/up/")
### CAF Down
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/CAF.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/overlap_CD4T_CD8T/CAF/down/")

### TEC UP
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/TEC.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/overlap_CD4T_CD8T/TEC/up/")

### TEC Down
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/TEC.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/overlap_CD4T_CD8T/TEC/down/")

### Malignant_cell UP
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/Malignant_cell.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/overlap_CD4T_CD8T/Malignant_cell/up/")

### Malignant_cell Down
pan_cancer_val(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/Malignant_cell.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/pan-cancer/overlap_CD4T_CD8T/Malignant_cell/down/")



#######################################################################################

# >>>>>>>>>>>>>> 2. Expression: Tumor_vs_Normal  <<<<<<<<<<<<<<<<<<<<<<<<<<

## >>> PAN-CANCER <<< ##
DEgeneExp_pancancer <- function(DEgenes, type, wkd){
  if(type == "UP"){
    de.genes <- readxl::read_excel(DEgenes,sheet = 1)
  }else if(type == "DOWN"){
    de.genes <- readxl::read_excel(DEgenes, sheet = 2)
  }else {
    stop("Error: Unrecognized Differential Gene Types !!!")
  }
  
  colnames(de.genes)[1] <- "Genes"
  genes <- de.genes$Genes[which(de.genes$Genes %in% colnames(TCGA_36cancers.rnaseq$ACC.rnaseq))]
  expr <- expressionsTCGA(
    TCGA_36cancers.rnaseq$ACC.rnaseq, TCGA_36cancers.rnaseq$BLCA.rnaseq,
    TCGA_36cancers.rnaseq$BRCA.rnaseq, TCGA_36cancers.rnaseq$CESC.rnaseq,
    TCGA_36cancers.rnaseq$CHOL.rnaseq, TCGA_36cancers.rnaseq$COAD.rnaseq,
    TCGA_36cancers.rnaseq$COADREAD.rnaseq, TCGA_36cancers.rnaseq$DLBC.rnaseq,
    TCGA_36cancers.rnaseq$ESCA.rnaseq, TCGA_36cancers.rnaseq$GBM.rnaseq,
    TCGA_36cancers.rnaseq$GBMLGG.rnaseq, TCGA_36cancers.rnaseq$HNSC.rnaseq,
    TCGA_36cancers.rnaseq$KICH.rnaseq, TCGA_36cancers.rnaseq$KIPAN.rnaseq,
    TCGA_36cancers.rnaseq$KIRC.rnaseq, TCGA_36cancers.rnaseq$KIRP.rnaseq,
    TCGA_36cancers.rnaseq$LAML.rnaseq, TCGA_36cancers.rnaseq$LGG.rnaseq,
    TCGA_36cancers.rnaseq$LIHC.rnaseq, TCGA_36cancers.rnaseq$LUAD.rnaseq,
    TCGA_36cancers.rnaseq$LUSC.rnaseq, TCGA_36cancers.rnaseq$OV.rnaseq,
    TCGA_36cancers.rnaseq$PAAD.rnaseq, TCGA_36cancers.rnaseq$PCPG.rnaseq,
    TCGA_36cancers.rnaseq$PRAD.rnaseq, TCGA_36cancers.rnaseq$READ.rnaseq,
    TCGA_36cancers.rnaseq$SARC.rnaseq, TCGA_36cancers.rnaseq$SKCM.rnaseq,
    TCGA_36cancers.rnaseq$STAD.rnaseq, TCGA_36cancers.rnaseq$STES.rnaseq,
    TCGA_36cancers.rnaseq$TGCT.rnaseq, TCGA_36cancers.rnaseq$THCA.rnaseq,
    TCGA_36cancers.rnaseq$THYM.rnaseq, TCGA_36cancers.rnaseq$UCEC.rnaseq,
    TCGA_36cancers.rnaseq$UCS.rnaseq, TCGA_36cancers.rnaseq$UVM.rnaseq,
    extract.cols = genes
  )
  expr$dataset <- gsub(".rnaseq","",str_split(expr$dataset, "\\$", simplify = T)[,2])
  Type <- ifelse(substr(expr$bcr_patient_barcode,14,15) == "11", "Normal", "Tumor")
  expr <- data.frame(tibble::add_column(expr, Type, .before = 2))
  tmp <- t(log2(t(expr[,4:ncol(expr)])+1))
  norm.exp <- cbind(expr[,1:3],tmp)
  colnames(norm.exp) <- gsub("\\.","-",colnames(norm.exp))
  for(item in genes){
    p <- ggpubr::ggboxplot(norm.exp, x = "Type", y = item, title = item, ylab = "Expression",
                           color = "Type", palette = "jco",bxp.errorbar = TRUE) +
    ggpubr::stat_compare_means(comparisons = list(c("Tumor","Normal"))) +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                     legend.position = "none", axis.title.x = element_blank())
    ggsave(paste0(wkd,item,".pdf"), p, device = "pdf", width = 5,height = 5)
  }
}

#======= Validation for all Differentially Expressed Genes ====================

if(F){
  # CD3T Up
  DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",type = "UP",
                      wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD3T/with_CD45/up/")
  # CD3T Down
  DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",type = "DOWN",
                      wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD3T/with_CD45/down/")
  #CD4T Up
  DEgeneExp_pancancer("result/without_PVTT_MLN/tumor/CD4T/CD4.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "UP",
                      wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD4T/with_CD45/up/")
  # CD4T Down
  DEgeneExp_pancancer("result/without_PVTT_MLN/tumor/CD4T/CD4.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "DOWN",
                      wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD4T/with_CD45/down/")
  #CD8T Up
  DEgeneExp_pancancer("result/without_PVTT_MLN/tumor/CD8T/CD8.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "UP",
                      wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD8T/with_CD45/up/")
  # CD8T Down
  DEgeneExp_pancancer("result/without_PVTT_MLN/tumor/CD8T/CD8.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "DOWN",
                      wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD8T/with_CD45/down/")
}

#======= Validation for all Differentially Expressed Genes (CD45 removed) ====================


## >>>>>>>>>>>> CD3+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",type = "UP",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD3T/without_CD45/CAF/up/")
### CAF Down
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",type = "DOWN",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD3T/without_CD45/CAF/down/")
### TEC UP
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",type = "UP",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD3T/without_CD45/TEC/up/")
### TEC Down
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",type = "DOWN",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD3T/without_CD45/TEC/down/")
### Malignant_Cell UP
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",type = "UP",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD3T/without_CD45/Malignant_cell/up/")
### Malignant_Cell Down
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",type = "DOWN",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD3T/without_CD45/Malignant_cell/down/")

## >>>>>>>>>>>> CD4+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",type = "UP",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD4T/without_CD45/CAF/up/")
### CAF Down
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",type = "DOWN",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD4T/without_CD45/CAF/down/")
### TEC UP
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",type = "UP",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD4T/without_CD45/TEC/up/")
### TEC Down
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",type = "DOWN",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD4T/without_CD45/TEC/down/")
### Malignant_Cell UP
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",type = "UP",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD4T/without_CD45/Malignant_cell/up/")
### Malignant_Cell Down
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",type = "DOWN",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD4T/without_CD45/Malignant_cell/down/")


## >>>>>>>>>>>> CD8+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",type = "UP",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD8T/without_CD45/CAF/up/")
### CAF Down
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",type = "DOWN",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD8T/without_CD45/CAF/down/")
### TEC UP
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",type = "UP",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD8T/without_CD45/TEC/up/")
### TEC Down
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",type = "DOWN",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD8T/without_CD45/TEC/down/")
### Malignant_Cell UP
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",type = "UP",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD8T/without_CD45/Malignant_cell/up/")
### Malignant_Cell Down
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",type = "DOWN",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/CD8T/without_CD45/Malignant_cell/down/")


### >>>>>>>>>>>> OVERLAPS between CD4+T & CD8+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/CAF.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/overlap_CD4T_CD8T/CAF/up/")
### CAF Down
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/CAF.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/overlap_CD4T_CD8T/CAF/down/")

### TEC UP
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/TEC.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/overlap_CD4T_CD8T/TEC/up/")

### TEC Down
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/TEC.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/overlap_CD4T_CD8T/TEC/down/")

### Malignant_cell UP
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/Malignant_cell.xlsx", type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/overlap_CD4T_CD8T/Malignant_cell/up/")

### Malignant_cell Down
DEgeneExp_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/Malignant_cell.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/pan-cancer/overlap_CD4T_CD8T/Malignant_cell/down/")



## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> LIHC <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

DEgeneExp_LIHC <- function(DEgenes, type, wkd){
  if(type == "UP"){
    de.genes <- readxl::read_excel(DEgenes,sheet = 1)
  }else if(type == "DOWN"){
    de.genes <- readxl::read_excel(DEgenes, sheet = 2)
  }else {
    stop("Error: Unrecognized Differential Gene Types !!!")
  }
  
  colnames(de.genes)[1] <- "Genes"
  genes <- de.genes$Genes[which(de.genes$Genes %in% colnames(TCGA_36cancers.rnaseq$ACC.rnaseq))]
  expr <- expressionsTCGA(
    TCGA_36cancers.rnaseq$ACC.rnaseq, TCGA_36cancers.rnaseq$BLCA.rnaseq,
    TCGA_36cancers.rnaseq$BRCA.rnaseq, TCGA_36cancers.rnaseq$CESC.rnaseq,
    TCGA_36cancers.rnaseq$CHOL.rnaseq, TCGA_36cancers.rnaseq$COAD.rnaseq,
    TCGA_36cancers.rnaseq$COADREAD.rnaseq, TCGA_36cancers.rnaseq$DLBC.rnaseq,
    TCGA_36cancers.rnaseq$ESCA.rnaseq, TCGA_36cancers.rnaseq$GBM.rnaseq,
    TCGA_36cancers.rnaseq$GBMLGG.rnaseq, TCGA_36cancers.rnaseq$HNSC.rnaseq,
    TCGA_36cancers.rnaseq$KICH.rnaseq, TCGA_36cancers.rnaseq$KIPAN.rnaseq,
    TCGA_36cancers.rnaseq$KIRC.rnaseq, TCGA_36cancers.rnaseq$KIRP.rnaseq,
    TCGA_36cancers.rnaseq$LAML.rnaseq, TCGA_36cancers.rnaseq$LGG.rnaseq,
    TCGA_36cancers.rnaseq$LIHC.rnaseq, TCGA_36cancers.rnaseq$LUAD.rnaseq,
    TCGA_36cancers.rnaseq$LUSC.rnaseq, TCGA_36cancers.rnaseq$OV.rnaseq,
    TCGA_36cancers.rnaseq$PAAD.rnaseq, TCGA_36cancers.rnaseq$PCPG.rnaseq,
    TCGA_36cancers.rnaseq$PRAD.rnaseq, TCGA_36cancers.rnaseq$READ.rnaseq,
    TCGA_36cancers.rnaseq$SARC.rnaseq, TCGA_36cancers.rnaseq$SKCM.rnaseq,
    TCGA_36cancers.rnaseq$STAD.rnaseq, TCGA_36cancers.rnaseq$STES.rnaseq,
    TCGA_36cancers.rnaseq$TGCT.rnaseq, TCGA_36cancers.rnaseq$THCA.rnaseq,
    TCGA_36cancers.rnaseq$THYM.rnaseq, TCGA_36cancers.rnaseq$UCEC.rnaseq,
    TCGA_36cancers.rnaseq$UCS.rnaseq, TCGA_36cancers.rnaseq$UVM.rnaseq,
    extract.cols = genes
  )
  expr$dataset <- gsub(".rnaseq","",str_split(expr$dataset, "\\$", simplify = T)[,2])
  Type <- ifelse(substr(expr$bcr_patient_barcode,14,15) == "11", "Normal", "Tumor")
  expr <- data.frame(tibble::add_column(expr, Type, .before = 2))
  tmp <- t(log2(t(expr[,4:ncol(expr)])+1))
  norm.exp <- cbind(expr[,1:3],tmp)
  colnames(norm.exp) <- gsub("\\.","-",colnames(norm.exp))
  LIHC <- norm.exp[norm.exp$dataset == "LIHC", ]
  for(item in genes){
    p <- ggpubr::ggboxplot(LIHC, x = "Type", y = item, title = item, ylab = "Expression",
                           color = "Type", palette = "jco",bxp.errorbar = TRUE) +
      ggpubr::stat_compare_means(comparisons = list(c("Tumor","Normal"))) +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                     legend.position = "none", axis.title.x = element_blank())
    ggsave(paste0(wkd,item,".pdf"), p, device = "pdf", width = 5,height = 5)
  }
}

#======= Validation for all Differentially Expressed Genes ====================

if(F){
  # CD3T Up
  DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",type = "UP",
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD3T/up/")
  # CD3T Down
  DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",type = "DOWN",
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD3T/down/")
  #CD4T Up
  DEgeneExp_LIHC("result/without_PVTT_MLN/tumor/CD4T/CD4.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "UP",
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD4T/up/")
  # CD4T Down
  DEgeneExp_LIHC("result/without_PVTT_MLN/tumor/CD4T/CD4.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "DOWN",
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD4T/down/")
  #CD8T Up
  DEgeneExp_LIHC("result/without_PVTT_MLN/tumor/CD8T/CD8.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "UP",
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD8T/up/")
  # CD8T Down
  DEgeneExp_LIHC("result/without_PVTT_MLN/tumor/CD8T/CD8.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "DOWN",
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD8T/down/")
}

#======= Validation for all Differentially Expressed Genes (CD45 removed) ====================

## >>>>>>>>>>>> CD3+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD3T/without_CD45/CAF/up/")
### CAF Down
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD3T/without_CD45/CAF/down/")

### TEC UP
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD3T/without_CD45/TEC/up/")
### TEC Down
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD3T/without_CD45/TEC/down/")

### Malignant_Cell UP
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD3T/without_CD45/Malignant_cell/up/")
### Malignant_Cell Down
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD3T/without_CD45/Malignant_cell/down/")


## >>>>>>>>>>>> CD4+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD4T/without_CD45/CAF/up/")
### CAF Down
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD4T/without_CD45/CAF/down/")

### TEC UP
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD4T/without_CD45/TEC/up/")
### TEC Down
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD4T/without_CD45/TEC/down/")

### Malignant_Cell UP
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD4T/without_CD45/Malignant_cell/up/")
### Malignant_Cell Down
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD4T/without_CD45/Malignant_cell/down/")


## >>>>>>>>>>>> CD8+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD8T/without_CD45/CAF/up/")
### CAF Down
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD8T/without_CD45/CAF/down/")

### TEC UP
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD8T/without_CD45/TEC/up/")
### TEC Down
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD8T/without_CD45/TEC/down/")

### Malignant_Cell UP
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",type = "UP",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD8T/without_CD45/Malignant_cell/up/")
### Malignant_Cell Down
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/CD8T/without_CD45/Malignant_cell/down/")


### >>>>>>>>>>>> OVERLAPS between CD4+T & CD8+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/CAF.xlsx", type = "UP",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/overlap_CD4T_CD8T/CAF/up/")
### CAF Down
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/CAF.xlsx", type = "DOWN",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/overlap_CD4T_CD8T/CAF/down/")

### TEC UP
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/TEC.xlsx", type = "UP",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/overlap_CD4T_CD8T/TEC/up/")

### TEC Down
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/TEC.xlsx", type = "DOWN",
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/overlap_CD4T_CD8T/TEC/down/")

### Malignant_cell UP
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/Malignant_cell.xlsx", type = "UP",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/overlap_CD4T_CD8T/Malignant_cell/up/")

### Malignant_cell Down
DEgeneExp_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/Malignant_cell.xlsx", type = "DOWN",
                    wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/expression/LIHC/overlap_CD4T_CD8T/Malignant_cell/down/")



#######################################################################################

# >>>>>>>>>>>>>> 3. Survival Analysis  <<<<<<<<<<<<<<<<<<<<<<<<<<

library(RTCGA.clinical)

## >>> PAN-CANCER <<< ##

clin <- survivalTCGA(
  ACC.clinical, BLCA.clinical, BRCA.clinical, CESC.clinical, CHOL.clinical, COAD.clinical,
  COADREAD.clinical, DLBC.clinical, ESCA.clinical, GBM.clinical, GBMLGG.clinical, HNSC.clinical,
  KICH.clinical, KIPAN.clinical, KIRC.clinical, KIRP.clinical, LAML.clinical, LGG.clinical,
  LIHC.clinical, LUAD.clinical, LUSC.clinical, OV.clinical, PAAD.clinical, PCPG.clinical,
  PRAD.clinical, READ.clinical, SARC.clinical, SKCM.clinical, STAD.clinical, STES.clinical,
  TGCT.clinical, THCA.clinical, THYM.clinical, UCEC.clinical, UCS.clinical, UVM.clinical
)

exp.TCGA <- expressionsTCGA(
  TCGA_36cancers.rnaseq$ACC.rnaseq, TCGA_36cancers.rnaseq$BLCA.rnaseq,
  TCGA_36cancers.rnaseq$BRCA.rnaseq, TCGA_36cancers.rnaseq$CESC.rnaseq,
  TCGA_36cancers.rnaseq$CHOL.rnaseq, TCGA_36cancers.rnaseq$COAD.rnaseq,
  TCGA_36cancers.rnaseq$COADREAD.rnaseq, TCGA_36cancers.rnaseq$DLBC.rnaseq,
  TCGA_36cancers.rnaseq$ESCA.rnaseq, TCGA_36cancers.rnaseq$GBM.rnaseq,
  TCGA_36cancers.rnaseq$GBMLGG.rnaseq, TCGA_36cancers.rnaseq$HNSC.rnaseq,
  TCGA_36cancers.rnaseq$KICH.rnaseq, TCGA_36cancers.rnaseq$KIPAN.rnaseq,
  TCGA_36cancers.rnaseq$KIRC.rnaseq, TCGA_36cancers.rnaseq$KIRP.rnaseq,
  TCGA_36cancers.rnaseq$LAML.rnaseq, TCGA_36cancers.rnaseq$LGG.rnaseq,
  TCGA_36cancers.rnaseq$LIHC.rnaseq, TCGA_36cancers.rnaseq$LUAD.rnaseq,
  TCGA_36cancers.rnaseq$LUSC.rnaseq, TCGA_36cancers.rnaseq$OV.rnaseq,
  TCGA_36cancers.rnaseq$PAAD.rnaseq, TCGA_36cancers.rnaseq$PCPG.rnaseq,
  TCGA_36cancers.rnaseq$PRAD.rnaseq, TCGA_36cancers.rnaseq$READ.rnaseq,
  TCGA_36cancers.rnaseq$SARC.rnaseq, TCGA_36cancers.rnaseq$SKCM.rnaseq,
  TCGA_36cancers.rnaseq$STAD.rnaseq, TCGA_36cancers.rnaseq$STES.rnaseq,
  TCGA_36cancers.rnaseq$TGCT.rnaseq, TCGA_36cancers.rnaseq$THCA.rnaseq,
  TCGA_36cancers.rnaseq$THYM.rnaseq, TCGA_36cancers.rnaseq$UCEC.rnaseq,
  TCGA_36cancers.rnaseq$UCS.rnaseq, TCGA_36cancers.rnaseq$UVM.rnaseq)

surv_pancancer <- function(DEgenes,type,expr,clindata, wkd){
  if(type == "UP"){
    de.genes <- readxl::read_excel(DEgenes,sheet = 1)
  }else if(type == "DOWN"){
    de.genes <- readxl::read_excel(DEgenes, sheet = 2)
  }else {
    stop("Error: Unrecognized Differential Gene Types !!!")
  }
  colnames(de.genes)[1] <- "Genes"
  genes <- de.genes$Genes[which(de.genes$Genes %in% colnames(TCGA_36cancers.rnaseq$ACC.rnaseq))]
  expr <- expressionsTCGA(expr,extract.cols = genes)
  expr$dataset <- gsub(".rnaseq","",str_split(expr$dataset, "\\$", simplify = T)[,2])
  Type <- ifelse(substr(expr$bcr_patient_barcode,14,15) == "11", "Normal", "Tumor")
  expr <- data.frame(tibble::add_column(expr, Type, .before = 2))
  tmp <- t(log2(t(expr[,4:ncol(expr)])+1))
  norm.exp <- cbind(expr[,1:3],tmp)
  colnames(norm.exp) <- gsub("\\.","-",colnames(norm.exp))
  norm.exp$bcr_patient_barcode <- substr(norm.exp$bcr_patient_barcode, 1, 12)
  survdata <- dplyr::inner_join(clindata,norm.exp,by="bcr_patient_barcode")
  my.surv <- survival::Surv(survdata$times,survdata$patient.vital_status)
  for(item in genes){
    # group
    group <- ifelse(survdata[,item] > median(na.omit(survdata[,item])),'high','low')
    sfit <- survival::survfit(my.surv~group, data = survdata)
    survp <- survminer::ggsurvplot(sfit, conf.int = FALSE, pval = FALSE)
    p1 <- survp$plot + ggplot2::annotate("text",label = paste0("p = ",signif(as.numeric(survminer::surv_pvalue(sfit)[2]),digits = 3)),
                                         x = 2000,y=0.10,size=5)
    ggsave(paste0(wkd,item,".pdf"),p1, width = 5,height = 5)
  }
}

#======= Validation for all Differentially Expressed Genes ====================

if(F){
  # CD3T Up
  surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",type = "UP",
                 expr = exp.TCGA, clindata = clin,
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD3T/up/")
  # CD3T Down
  surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",type = "DOWN",
                 expr = exp.TCGA, clindata = clin,
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD3T/down/")
  #CD4T Up
  surv_pancancer("result/without_PVTT_MLN/tumor/CD4T/CD4.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "UP",
                 expr = exp.TCGA, clindata = clin,
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD4T/up/")
  # CD4T Down
  surv_pancancer("result/without_PVTT_MLN/tumor/CD4T/CD4.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "DOWN",
                 expr = exp.TCGA, clindata = clin,
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD4T/down/")
  #CD8T Up
  surv_pancancer("result/without_PVTT_MLN/tumor/CD8T/CD8.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "UP",
                 expr = exp.TCGA, clindata = clin,
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD8T/up/")
  # CD8T Down
  surv_pancancer("result/without_PVTT_MLN/tumor/CD8T/CD8.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "DOWN",
                 expr = exp.TCGA, clindata = clin,
                 wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD8T/down/")
}

#======= Validation for all Differentially Expressed Genes (CD45 removed) ====================

## >>>>>>>>>>>> CD3+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD3T/without_CD45/CAF/up/")
### CAF Down
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD3T/without_CD45/CAF/down/")
### TEC UP
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD3T/without_CD45/TEC/up/")
### TEC Down
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD3T/without_CD45/TEC/down/")
### Malignant_cell UP
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD3T/without_CD45/Malignant_cell/up/")
### Malignant_cell Down
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD3T/without_CD45/Malignant_cell/down/")


## >>>>>>>>>>>> CD4+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD4T/without_CD45/CAF/up/")
### CAF Down
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD4T/without_CD45/CAF/down/")
### TEC UP
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD4T/without_CD45/TEC/up/")
### TEC Down
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD4T/without_CD45/TEC/down/")
### Malignant_cell UP
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD4T/without_CD45/Malignant_cell/up/")
### Malignant_cell Down
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD4T/without_CD45/Malignant_cell/down/")


## >>>>>>>>>>>> CD8+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD8T/without_CD45/CAF/up/")
### CAF Down
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD8T/without_CD45/CAF/down/")
### TEC UP
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD8T/without_CD45/TEC/up/")
### TEC Down
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD8T/without_CD45/TEC/down/")
### Malignant_cell UP
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD8T/without_CD45/Malignant_cell/up/")
### Malignant_cell Down
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/CD8T/without_CD45/Malignant_cell/down/")


### >>>>>>>>>>>> OVERLAPS between CD4+T & CD8+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/CAF.xlsx", type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/overlap_CD4T_CD8T/CAF/up/")
### CAF Down
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/CAF.xlsx", type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/overlap_CD4T_CD8T/CAF/down/")
### TEC UP
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/TEC.xlsx", type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/overlap_CD4T_CD8T/TEC/up/")
### TEC Down
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/TEC.xlsx", type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/overlap_CD4T_CD8T/TEC/down/")
### Malignant_cell UP
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/Malignant_cell.xlsx", type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/overlap_CD4T_CD8T/Malignant_cell/up/")
### Malignant_cell Down
surv_pancancer(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/Malignant_cell.xlsx", type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/pan-cancer/overlap_CD4T_CD8T/Malignant_cell/down/")



## >>>>>>>> LIHC <<<<<<<<< ##

surv_LIHC <- function(DEgenes,type,expr,clindata, wkd){
  if(type == "UP"){
    de.genes <- readxl::read_excel(DEgenes,sheet = 1)
  }else if(type == "DOWN"){
    de.genes <- readxl::read_excel(DEgenes, sheet = 2)
  }else {
    stop("Error: Unrecognized Differential Gene Types !!!")
  }
  colnames(de.genes)[1] <- "Genes"
  genes <- de.genes$Genes[which(de.genes$Genes %in% colnames(TCGA_36cancers.rnaseq$ACC.rnaseq))]
  expr <- expressionsTCGA(expr,extract.cols = genes)
  expr$dataset <- gsub(".rnaseq","",str_split(expr$dataset, "\\$", simplify = T)[,2])
  Type <- ifelse(substr(expr$bcr_patient_barcode,14,15) == "11", "Normal", "Tumor")
  expr <- data.frame(tibble::add_column(expr, Type, .before = 2))
  tmp <- t(log2(t(expr[,4:ncol(expr)])+1))
  norm.exp <- cbind(expr[,1:3],tmp)
  colnames(norm.exp) <- gsub("\\.","-",colnames(norm.exp))
  norm.exp$bcr_patient_barcode <- substr(norm.exp$bcr_patient_barcode, 1, 12)
  LIHC <- norm.exp[norm.exp$dataset == "LIHC", ]
  survdata <- dplyr::inner_join(clindata,LIHC,by="bcr_patient_barcode")
  my.surv <- survival::Surv(survdata$times,survdata$patient.vital_status)
  for(item in genes){
    # group
    group <- ifelse(survdata[,item] > median(na.omit(survdata[,item])),'high','low')
    sfit <- survival::survfit(my.surv~group, data = survdata)
    survp <- survminer::ggsurvplot(sfit, conf.int = FALSE, pval = FALSE)
    p1 <- survp$plot + ggplot2::annotate("text",label = paste0("p = ",signif(as.numeric(survminer::surv_pvalue(sfit)[2]),digits = 3)),
                                         x = 2000,y=0.10,size=5)
    ggsave(paste0(wkd,item,".pdf"),p1, width = 5,height = 5)
  }
}

#======= Validation for all Differentially Expressed Genes ====================

if(F){
  # CD3T Up
  surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",type = "UP",
            expr = exp.TCGA, clindata = clin,
            wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD3T/up/")
  # CD3T Down
  surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD3.high_vs_low.DEmarkers.xlsx",type = "DOWN",
            expr = exp.TCGA, clindata = clin,
            wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD3T/down/")
  #CD4T Up
  surv_LIHC("result/without_PVTT_MLN/tumor/CD4T/CD4.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "UP",
            expr = exp.TCGA, clindata = clin,
            wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD4T/up/")
  # CD4T Down
  surv_LIHC("result/without_PVTT_MLN/tumor/CD4T/CD4.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "DOWN",
            expr = exp.TCGA, clindata = clin,
            wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD4T/down/")
  #CD8T Up
  surv_LIHC("result/without_PVTT_MLN/tumor/CD8T/CD8.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "UP",
            expr = exp.TCGA, clindata = clin,
            wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD8T/up/")
  # CD8T Down
  surv_LIHC("result/without_PVTT_MLN/tumor/CD8T/CD8.high.tumor_vs_nonTumor.DEmarkers.xlsx",type = "DOWN",
            expr = exp.TCGA, clindata = clin,
            wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD8T/down/")
}

#======= Validation for all Differentially Expressed Genes (CD45 removed) ====================

## >>>>>>>>>>>> CD3+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD3T/without_CD45/CAF/up/")
### CAF Down
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD3T/without_CD45/CAF/down/")
### TEC UP
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD3T/without_CD45/TEC/up/")
### TEC Down
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD3T/without_CD45/TEC/down/")
### Malignant_cell UP
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD3T/without_CD45/Malignant_cell/up/")
### Malignant_cell Down
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD3T/without_CD45/Malignant_cell/down/")


## >>>>>>>>>>>>>>>>>> CD4+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD4T/without_CD45/CAF/up/")
### CAF Down
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD4T/without_CD45/CAF/down/")
### TEC UP
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD4T/without_CD45/TEC/up/")
### TEC Down
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD4T/without_CD45/TEC/down/")
### Malignant_cell UP
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD4T/without_CD45/Malignant_cell/up/")
### Malignant_cell Down
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD4T/without_CD45/Malignant_cell/down/")


## >>>>>>>>>>>> CD8+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD8T/without_CD45/CAF/up/")
### CAF Down
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD8T/without_CD45/CAF/down/")
### TEC UP
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD8T/without_CD45/TEC/up/")
### TEC Down
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD8T/without_CD45/TEC/down/")
### Malignant_cell UP
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",type = "UP",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD8T/without_CD45/Malignant_cell/up/")
### Malignant_cell Down
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",type = "DOWN",
               expr = exp.TCGA, clindata = clin,
               wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/CD8T/without_CD45/Malignant_cell/down/")


### >>>>>>>>>>>> OVERLAPS between CD4+T & CD8+T <<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/CAF.xlsx", type = "UP",
          expr = exp.TCGA, clindata = clin,
          wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/overlap_CD4T_CD8T/CAF/up/")
### CAF Down
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/CAF.xlsx", type = "DOWN",
          expr = exp.TCGA, clindata = clin,
          wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/overlap_CD4T_CD8T/CAF/down/")
### TEC UP
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/TEC.xlsx", type = "UP",
          expr = exp.TCGA, clindata = clin,
          wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/overlap_CD4T_CD8T/TEC/up/")
### TEC Down
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/TEC.xlsx", type = "DOWN",
          expr = exp.TCGA, clindata = clin,
          wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/overlap_CD4T_CD8T/TEC/down/")
### Malignant_cell UP
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/Malignant_cell.xlsx", type = "UP",
          expr = exp.TCGA, clindata = clin,
          wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/overlap_CD4T_CD8T/Malignant_cell/up/")
### Malignant_cell Down
surv_LIHC(DEgenes = "result/without_PVTT_MLN/tumor/overlap.CD4T_vs_CD8T/without_CD45/Malignant_cell.xlsx", type = "DOWN",
          expr = exp.TCGA, clindata = clin,
          wkd = "result/without_PVTT_MLN/tumor/validation/TCGA/survival/LIHC/overlap_CD4T_CD8T/Malignant_cell/down/")



#######################################################################################


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  CancerSEA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>> Prepare CancerGEA reference Database <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#================= Protein Coding Genes ====================================

#>>>> BeforeQC

if(!file.exists("refdb/CancerSEA/beforeQC/CancerSEA.PCG.metadata.rds")){
  system("for i in `seq 1 9`;do axel -n 30 http://biocc.hrbmu.edu.cn/CancerSEA/download/Expression/EXP000${i}_PCG.zip")
  system("for i in `seq 10 72`;do axel -n 30 http://biocc.hrbmu.edu.cn/CancerSEA/download/Expression/EXP00${i}_PCG.zip")
  system("for i in `seq 1 9`;do axel -n 30 http://biocc.hrbmu.edu.cn/CancerSEA/download/Expression/EXP000${i}_lnc.zip")
  system("for i in `seq 10 72`;do axel -n 30 http://biocc.hrbmu.edu.cn/CancerSEA/download/Expression/EXP00${i}_lnc.zip")
  
  Cancers <- rep(
    c("Alveolar rhabdomyosarcoma","Bronchoalveolar carcinoma","Acute lymphoblastic leukemia",
      "Chronic myelogenous leukemia","Colorectal cancer","Breast cancer","Cervix cancer",
      "Glioblastoma","Glioma","Neuroblastoma","Esophageal squamous cell carcinoma",
      "Hepatocellular carcinoma","Lung adenocarcinoma","Non-small cell lung cancer",
      "Pancreatic ductal adenocarcinoma","Prostate cancer","Melanoma","Myxoid liposarcoma",
      "Hepatocellular carcinoma","Acute myeloid leukemia","Chronic myelogenous leukemia",
      "Colorectal cancer","Breast cancer","Astrocytoma","Glioblastoma","Glioma",
      "High-grade glioma","Oligodendroglioma","Head and neck cancer","Renal cell carcinoma",
      "Lung adenocarcinoma","Non-small cell lung cancer","Ovarian carcinoma","Prostate cancer","Melanoma"
    ),
    times = c(1,1,1,7,2,4,2,2,1,1,1,2,7,4,1,2,3,1,1,5,1,1,4,1,2,2,1,1,1,2,2,1,1,1,2)
  )
  
  # PCG beforeQC
  pcg.exp.files <- list.files("refdb/CancerSEA/beforeQC/")
  meta <- list()
  pcg.exp_list <- list()
  for(i in seq_along(pcg.exp.files)){
    pcg.exp_list[[i]] <- read.table(paste0("refdb/CancerSEA/beforeQC/",pcg.exp.files[i]),header = T,sep = "\t",row.names = 1)
    meta[[i]] <- data.frame(
      CellNames = colnames(pcg.exp_list[[i]]),
      Samples = as.character(unlist(pcg.exp_list[[i]][1,],use.names = F)), 
      CanceType = Cancers[i]
    )
  }
  # metadata
  metadata <- Reduce(function(x,y) rbind(x,y), meta, accumulate = F)
  
  # gene expression matrix
  pcg.exp_list <- lapply(pcg.exp_list,function(x){x[-1,]})
  pcg.exp <- lapply(pcg.exp_list, function(x){tibble::add_column(x, ENSEMBLID = rownames(x), .before = 1)})
  library(dplyr)
  pcg.exp <- Reduce(function(x,y) full_join(x,y, by = "ENSEMBLID"),pcg.exp,accumulate = FALSE)
  rownames(pcg.exp) <- pcg.exp$ENSEMBLID
  pcg.exp$ENSEMBLID <- NULL
  pcg.exp.new <- apply(pcg.exp,2,as.numeric)
  rownames(pcg.exp.new) <- rownames(pcg.exp)
  pcg.exp.new[is.na(pcg.exp.new)] <- 0
  
  rownames(metadata) <- colnames(pcg.exp.new)
  metadata$CellNames <- NULL
  saveRDS(metadata, file = "refdb/CancerSEA/beforeQC/CancerSEA.PCG.metadata.rds")
  saveRDS(pcg.exp.new,file = "refdb/CancerSEA/beforeQC/CancerSEA.exp.rds")
} else {
  pcg.metadata <- readRDS("refdb/CancerSEA/beforeQC/CancerSEA.PCG.metadata.rds")
  pcg.exp <- readRDS("refdb/CancerSEA/beforeQC/CancerSEA.exp.rds")
}


#>>>> AfterQC

if(!file.exists("refdb/CancerSEA/log2_afterQC/CancerSEA.PCG.norm.metadata.rds")){
  pcg.norm.exp.files <- list.files("refdb/CancerSEA/log2_afterQC/")
  meta <- list()
  pcg.norm.exp_list <- list()
  for(i in seq_along(pcg.norm.exp.files)){
    pcg.norm.exp_list[[i]] <- read.table(paste0("refdb/CancerSEA/log2_afterQC/",pcg.norm.exp.files[i]),header = T,sep = "\t",row.names = 1)
    meta[[i]] <- data.frame(
      CellNames = colnames(pcg.norm.exp_list[[i]]),
      Samples = as.character(unlist(pcg.norm.exp_list[[i]][1, ], use.names = F)),
      CancerType = Cancers[i]
    )
  }
  metadata <- Reduce(function(x,y) rbind(x,y), meta, accumulate = F)
  
  pcg.norm.exp_list <- lapply(pcg.norm.exp_list, function(x) x[-1, ])
  pcg.norm.exp_list <- lapply(pcg.norm.exp_list, function(x){tibble::add_column(x, ENSEMBLID = rownames(x), .before = 1)})
  pcg.norm.exp <- Reduce(function(x,y) full_join(x,y, by = "ENSEMBLID"),pcg.norm.exp_list,accumulate = FALSE)
  rownames(pcg.norm.exp) <- pcg.norm.exp$ENSEMBLID
  pcg.norm.exp$ENSEMBLID <- NULL
  pcg.norm.exp.new <- apply(pcg.norm.exp,2,as.numeric)
  rownames(pcg.norm.exp.new) <- rownames(pcg.norm.exp)
  pcg.norm.exp.new[is.na(pcg.norm.exp.new)] <- 0
  rownames(metadata) <- colnames(pcg.norm.exp.new)
  metadata$CellNames <- NULL
  saveRDS(metadata, file = "refdb/CancerSEA/log2_afterQC/CancerSEA.PCG.norm.metadata.rds")
  saveRDS(pcg.norm.exp.new,file = "refdb/CancerSEA/log2_afterQC/CancerSEA.PCG.norm.exp.rds")
}else {
  pcg.norm.metadata <- readRDS("refdb/CancerSEA/log2_afterQC/CancerSEA.PCG.norm.metadata.rds")
  pcg.norm.exp <- readRDS("refdb/CancerSEA/log2_afterQC/CancerSEA.PCG.norm.exp.rds")
}


min.cells = 0
min.features = 0

pcg.sce <- CreateSeuratObject(counts = pcg.exp, meta.data = pcg.metadata, min.cells = min.cells, min.features = min.features)

saveRDS(pcg.sce,file = "refdb/CancerSEA/CancerSEA.sce.rds")

#>>>>>> lncRNA
if(F){
  Cancers <- c(
    "Alveolar rhabdomyosarcoma","Bronchoalveolar carcinoma","Acute lymphoblastic leukemia",
    "Chronic myelogenous leukemia","Chronic myelogenous leukemia","Chronic myelogenous leukemia",
    "Chronic myelogenous leukemia","Colorectal cancer","Colorectal cancer","Breast cancer",
    "Breast cancer","Cervix cancer","Cervix cancer","Glioblastoma","Glioblastoma","Glioma",
    "Neuroblastoma","Esophageal squamous cell carcinoma","Hepatocellular carcinoma",
    "Hepatocellular carcinoma","Lung adenocarcinoma","Lung adenocarcinoma","Lung adenocarcinoma",
    "Lung adenocarcinoma","Lung adenocarcinoma","Lung adenocarcinoma","Lung adenocarcinoma",
    "Non-small cell lung cancer","Non-small cell lung cancer","Non-small cell lung cancer",
    "Non-small cell lung cancer","Pancreatic ductal adenocarcinoma","Prostate cancer",
    "Prostate cancer","Melanoma","Melanoma","Melanoma","Myxoid liposarcoma","Hepatocellular carcinoma",
    "Acute myeloid leukemia","Acute lymphoblastic leukemia","Acute lymphoblastic leukemia",
    "Acute lymphoblastic leukemia","Acute lymphoblastic leukemia","Chronic myelogenous leukemia",
    "Colorectal cancer","Breast cancer","Breast cancer","Breast cancer","Breast cancer","Astrocytoma",
    "Glioblastoma","Glioblastoma","High-grade glioma","Renal cell carcinoma","Renal cell carcinoma",
    "Lung adenocarcinoma","Non-small cell lung cancer","Melanoma"
  )
  
  
  lnc.exp.files <- list.files("/home/tongqiang/reference/CancerSEA/lncRNA/beforeQC/")
  
  meta <- list()
  lnc.exp.list <- list()
  for(i in seq_along(lnc.exp.files)){
    lnc.exp.list[[i]] <- read.table(paste0("/home/tongqiang/reference/CancerSEA/lncRNA/beforeQC/",lnc.exp.files[i]),header = T,sep = "\t")
    lnc.exp.list[[i]] <- lnc.exp.list[[i]][!duplicated(lnc.exp.list[[i]][,1]),]
    rownames(lnc.exp.list[[i]]) <- lnc.exp.list[[i]][,1]
    meta[[i]] <- data.frame(
      CellNames = colnames(lnc.exp.list[[i]]),
      Samples = as.character(unlist(lnc.exp.list[[i]][1,],use.names = F)),
      CanceType = Cancers[i]
    )
  }
  # metadata
  metadta <- Reduce(function(x,y) rbind(x,y), meta, accumulate = F)
  saveRDS(metadata, file = "/home/tongqiang/reference/CancerSEA/lncRNA/beforeQC/CancerSEA.PCG.metadata.rds")
  # gene expression matrix
  lnc.exp.list <- lapply(lnc.exp.list,function(x){x[-1,]})
  lnc.exp <- lapply(lnc.exp.list, function(x){tibble::add_column(x, ENSEMBLID = rownames(x), .before = 1)})
  
  library(dplyr)
  lnc.exp <- Reduce(function(x,y) full_join(x,y, by = "ENSEMBLID"),lnc.exp,accumulate = FALSE)
  rownames(lnc.exp) <- lnc.exp$ENSEMBLID
  lnc.exp$ENSEMBLID <- NULL
  lnc.exp.new <- apply(lnc.exp,2,as.numeric)
  rownames(lnc.exp.new) <- rownames(lnc.exp)
  lnc.exp.new[is.na(lnc.exp.new)] <- 0
  
  saveRDS(lnc.exp.new,file = "/home/tongqiang/reference/CancerSEA/lncRNA/beforeQC/CancerSEA.lnc.metadata.rds")
  
  # lncRNA  afterQC
  lnc.norm.exp.files <- list.files("/home/tongqiang/reference/CancerSEA/lncRNA/log2_afterQC/")
  meta <- list()
  lnc.norm.exp_list <- list()
  for(i in seq_along(lnc.norm.exp.files)){
    lnc.norm.exp_list[[i]] <- read.table(paste0("/home/tongqiang/reference/CancerSEA/lncRNA/log2_afterQC/",lnc.norm.exp.files[i]),header = T,sep = "\t")
    lnc.norm.exp_list[[i]] <- lnc.norm.exp_list[[i]][!duplicated(lnc.norm.exp_list[[i]][,1]),]
    rownames(lnc.norm.exp_list[[i]]) <- lnc.norm.exp_list[[i]][,1]
    meta[[i]] <- data.frame(
      CellNames = colnames(lnc.norm.exp_list[[i]]),
      Samples = as.character(unlist(lnc.norm.exp_list[[i]][1, ], use.names = F)),
      CancerType = Cancers[i]
    )
  }
  metadta <- Reduce(function(x,y) rbind(x,y), meta, accumulate = F)
  saveRDS(metadata, file = "/home/tongqiang/reference/CancerSEA/log2_afterQC/CancerSEA.lnc.norm.metadata.rds")
  
  lnc.norm.exp_list <- lapply(lnc.norm.exp_list, function(x) x[-1, ])
  lnc.norm.exp_list <- lapply(lnc.norm.exp_list, function(x){tibble::add_column(x, ENSEMBLID = rownames(x), .before = 1)})
  lnc.norm.exp <- Reduce(function(x,y) full_join(x,y, by = "ENSEMBLID"),lnc.norm.exp,accumulate = FALSE)
  rownames(lnc.norm.exp) <- lnc.norm.exp$ENSEMBLID
  lnc.norm.exp$ENSEMBLID <- NULL
  lnc.norm.exp.new <- apply(lnc.norm.exp,2,as.numeric)
  rownames(lnc.norm.exp.new) <- rownames(lnc.norm.exp)
  lnc.norm.exp.new[is.na(lnc.norm.exp.new)] <- 0
  
  saveRDS(lnc.norm.exp.new,file = "/home/tongqiang/reference/CancerSEA/log2_afterQC/CancerSEA.lnc.norm.exp.rds")
}


if(!file.exists("refdb/GeneSymbol2ENSEMBL.rds")){
  library(org.Hs.eg.db)
  eg2symbol=toTable(org.Hs.egSYMBOL)
  eg2ensembl = toTable(org.Hs.egENSEMBL)
  eg2symbol$ENSEMBL <- plyr::mapvalues(eg2symbol$gene_id, from = eg2ensembl$gene_id, to = eg2ensembl$ensembl_id)
  # mapped genes
  mapped <- eg2symbol[grepl("^ENSG", eg2symbol$ENSEMBL),]
  # unmapped gene id
  ummapped <- eg2symbol[!grepl("^ENSG", eg2symbol$ENSEMBL),]
  write.table(ummapped$ENSEMBL, file = "tmp/uampped_id1.txt",col.names = F,row.names = F,quote = F)
  # round 1
  new_id <- read.table("bioDBnet_db2db_210712015101_1684696328.txt",header = T,sep = "\t",stringsAsFactors = F)
  ummapped$ENSEMBL <- plyr::mapvalues(ummapped$gene_id,from = new_id$Gene.ID, to = new_id$Ensembl.Gene.ID)
  mapped2 <- ummapped[grepl("^ENSG",ummapped$ENSEMBL),]
  ummapped2 <- ummapped[!grepl("^ENSG",ummapped$ENSEMBL),]
  mapped_new <- rbind(mapped,mapped2)
  write.table(ummapped2, file = "tmp/uampped_id2.txt",col.names = F,row.names = F,quote = F,sep = "\t")
  # round 2
  new_id <- read.table("tmp/biodb.txt",header = T,sep = "\t",stringsAsFactors = F)
  ummapped2$ENSEMBL <- plyr::mapvalues(ummapped2$symbol, from = new_id$Gene.ID, to = new_id$Ensembl.Gene.ID)
  mapped3 <- ummapped2[grepl("^ENSG",ummapped2$ENSEMBL),]
  mapped_new <- rbind(mapped_new, mapped3)
  ummapped3 <- ummapped2[!grepl("ENSG",ummapped2$ENSEMBL),]
  ummapped3$ENSEMBL <- NULL
  write.table(ummapped3, file = "tmp/uampped_id3.txt",col.names = F,row.names = F,quote = F,sep = "\t")
  
  saveRDS(mapped_new,file = "refdb/GeneSymbol2ENSEMBL.rds")
  saveRDS(ummapped3, file = "refdb/unTransformedGeneSymbol.rds")
}else {
  mapped <- readRDS("refdb/GeneSymbol2ENSEMBL.rds")
}


#======= Validation for all Differentially Expressed Genes (CD45 removed) ====================

#======= ADD ENSEMBL ID ====================

####>>>>>>>>>>>>>> CD3 T <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
CD3.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",sheet = 1) 
colnames(CD3.CAF.up)[1] <- "Genes"
CD3.CAF.up$ENSEMBL <- plyr::mapvalues(CD3.CAF.up$Genes, from = mapped$symbol, to = mapped$ENSEMBL)
### CAF Down
CD3.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD3.CAF.down)[1] <- "Genes"
CD3.CAF.down$ENSEMBL <- plyr::mapvalues(CD3.CAF.down$Genes, from = mapped$symbol, to = mapped$ENSEMBL)

xlsx::write.xlsx(data.frame(CD3.CAF.up),file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",sheetName = "UP",row.names = FALSE)
xlsx::write.xlsx(data.frame(CD3.CAF.down),file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/CAF/CAF.CD3T.hi_vs_lo.xlsx",sheetName = "Down",append = TRUE,row.names = FALSE)

### TEC UP
CD3.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",sheet = 1) 
colnames(CD3.TEC.up)[1] <- "Genes"
CD3.TEC.up$ENSEMBL <- plyr::mapvalues(CD3.TEC.up$Genes, from = mapped$symbol, to = mapped$ENSEMBL)
### TEC Down
CD3.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",sheet = 2) 
colnames(CD3.TEC.down)[1] <- "Genes"
CD3.TEC.down$ENSEMBL <- plyr::mapvalues(CD3.TEC.down$Genes, from = mapped$symbol, to = mapped$ENSEMBL)

xlsx::write.xlsx(data.frame(CD3.TEC.up),file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",sheetName = "UP",row.names = FALSE)
xlsx::write.xlsx(data.frame(CD3.TEC.down),file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/TEC/TEC.CD3T.hi_vs_lo.xlsx",sheetName = "Down",append = TRUE,row.names = FALSE)

### Malignant.Cell UP
CD3.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD3.MLG.up)[1] <- "Genes"
CD3.MLG.up$ENSEMBL <- plyr::mapvalues(CD3.MLG.up$Genes, from = mapped$symbol, to = mapped$ENSEMBL)
### Malignant.Cell Down
CD3.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD3.MLG.down)[1] <- "Genes"
CD3.MLG.down$ENSEMBL <- plyr::mapvalues(CD3.MLG.down$Genes, from = mapped$symbol, to = mapped$ENSEMBL)

xlsx::write.xlsx(data.frame(CD3.MLG.up),file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",sheetName = "UP",row.names = F)
xlsx::write.xlsx(data.frame(CD3.MLG.down),file = "result/without_PVTT_MLN/tumor/CD3T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD3T.hi_vs_lo.xlsx",sheetName = "Down",row.names = F,append = T)

####>>>>>>>>>>>>>> CD4 T <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
CD4.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 1) 
colnames(CD4.CAF.up)[1] <- "Genes"
CD4.CAF.up$ENSEMBL <- plyr::mapvalues(CD4.CAF.up$Genes, from = mapped$symbol, to = mapped$ENSEMBL)
### CAF Down
CD4.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD4.CAF.down)[1] <- "Genes"
CD4.CAF.down$ENSEMBL <- plyr::mapvalues(CD4.CAF.down$Genes, from = mapped$symbol, to = mapped$ENSEMBL)

xlsx::write.xlsx(data.frame(CD4.CAF.up),file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheetName = "UP",row.names = FALSE)
xlsx::write.xlsx(data.frame(CD4.CAF.down),file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheetName = "Down",append = TRUE,row.names = FALSE)

### TEC UP
CD4.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 1) 
colnames(CD4.TEC.up)[1] <- "Genes"
CD4.TEC.up$ENSEMBL <- plyr::mapvalues(CD4.TEC.up$Genes, from = mapped$symbol, to = mapped$ENSEMBL)
### TEC Down
CD4.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 2) 
colnames(CD4.TEC.down)[1] <- "Genes"
CD4.TEC.down$ENSEMBL <- plyr::mapvalues(CD4.TEC.down$Genes, from = mapped$symbol, to = mapped$ENSEMBL)

xlsx::write.xlsx(data.frame(CD4.TEC.up),file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheetName = "UP",row.names = FALSE)
xlsx::write.xlsx(data.frame(CD4.TEC.down),file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheetName = "Down",append = TRUE,row.names = FALSE)

### Malignant.Cell UP
CD4.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD4.MLG.up)[1] <- "Genes"
CD4.MLG.up$ENSEMBL <- plyr::mapvalues(CD4.MLG.up$Genes, from = mapped$symbol, to = mapped$ENSEMBL)
### Malignant.Cell Down
CD4.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD4.MLG.down)[1] <- "Genes"
CD4.MLG.down$ENSEMBL <- plyr::mapvalues(CD4.MLG.down$Genes, from = mapped$symbol, to = mapped$ENSEMBL)

xlsx::write.xlsx(data.frame(CD4.MLG.up),file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheetName = "UP",row.names = F)
xlsx::write.xlsx(data.frame(CD4.MLG.down),file = "result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheetName = "Down",row.names = F,append = T)


####>>>>>>>>>>>>>> CD8 T <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

### CAF UP
CD8.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 1) 
colnames(CD8.CAF.up)[1] <- "Genes"
CD8.CAF.up$ENSEMBL <- plyr::mapvalues(CD8.CAF.up$Genes, from = mapped$symbol, to = mapped$ENSEMBL)
### CAF Down
CD8.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD8.CAF.down)[1] <- "Genes"
CD8.CAF.down$ENSEMBL <- plyr::mapvalues(CD8.CAF.down$Genes, from = mapped$symbol, to = mapped$ENSEMBL)

xlsx::write.xlsx(data.frame(CD8.CAF.up),file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheetName = "UP",row.names = FALSE)
xlsx::write.xlsx(data.frame(CD8.CAF.down),file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheetName = "Down",append = TRUE,row.names = FALSE)

### TEC UP
CD8.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 1) 
colnames(CD8.TEC.up)[1] <- "Genes"
CD8.TEC.up$ENSEMBL <- plyr::mapvalues(CD8.TEC.up$Genes, from = mapped$symbol, to = mapped$ENSEMBL)
### TEC Down
CD8.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 2) 
colnames(CD8.TEC.down)[1] <- "Genes"
CD8.TEC.down$ENSEMBL <- plyr::mapvalues(CD8.TEC.down$Genes, from = mapped$symbol, to = mapped$ENSEMBL)

xlsx::write.xlsx(data.frame(CD8.TEC.up),file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheetName = "UP",row.names = FALSE)
xlsx::write.xlsx(data.frame(CD8.TEC.down),file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheetName = "Down",append = TRUE,row.names = FALSE)

### Malignant.Cell UP
CD8.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 1)
colnames(CD8.MLG.up)[1] <- "Genes"
CD8.MLG.up$ENSEMBL <- plyr::mapvalues(CD8.MLG.up$Genes, from = mapped$symbol, to = mapped$ENSEMBL)
### Malignant.Cell Down
CD8.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 2)
colnames(CD8.MLG.down)[1] <- "Genes"
CD8.MLG.down$ENSEMBL <- plyr::mapvalues(CD8.MLG.down$Genes, from = mapped$symbol, to = mapped$ENSEMBL)

xlsx::write.xlsx(data.frame(CD8.MLG.up),file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheetName = "UP",row.names = F)
xlsx::write.xlsx(data.frame(CD8.MLG.down),file = "result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheetName = "Down",row.names = F,append = T)



##======= Correlation between signature gens and DE genes =======

library(ggplot2)

##>>>>>>>>>>>>>>>>>>>>>>> 1.CD8+T <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 1.CD3+CD8+;
# 2.CD3+CD8+PD-1+Tim3+Lag3+；
# 3.CD3+CD8+IFNg；
# 4.CD3+CD8+GranzymeB；
# 5.CD3+CD8+TNFa；

##======= 1.CD3+CD8+ =======

geneSigCor <- function(pancancer,genes,sig.genes,workdir){
  if(!file.exists(paste(workdir,genes,sep="/"))){
    dir.create(paste(workdir,genes,sep="/"))
  }
  cancer <- gsub(".rnaseq","",names(pancancer))
  for(i in seq_along(cancer)){
    if(!file.exists(paste(workdir,genes,cancer[i],sep = "/"))){
      dir.create(paste(workdir,genes,cancer[i],sep="/"))
    }
  }
  dt <- lapply(pancancer,function(x){
    log2(x[,c(sig.genes,genes)] + 1)
  })
  r <- lapply(dt,function(x){
    r1 <- cor(x[,length(sig.genes)+1],y = x[,1],method = "pearson")
    r2 <- cor(x[,length(sig.genes)+1],y = x[,2],method = "pearson")
    r3 <- cor(x[,length(sig.genes)+1],y = x[,3],method = "pearson")
    r4 <- cor(x[,length(sig.genes)+1],y = x[,4],method = "pearson")
    r5 <- cor(x[,length(sig.genes)+1],y = x[,5],method = "pearson")
    c(r1,r2,r3,r4,r5)
  })
  cor_result <- t(data.frame(r))
  colnames(cor_result) <- sig.genes
  rownames(cor_result) <- gsub(".rnaseq","",rownames(cor_result))
  n_patient <- lapply(pancancer,function(x) return(nrow(x)))
  for(i in seq_along(n_patient)){
    rownames(cor_result)[i] <- paste0(rownames(cor_result)[i],"(n=",n_patient[i],")")
  }
  pheatmap::pheatmap(cor_result,cluster_rows = F,cluster_cols = F,display_numbers = T,legend = F,
                     fontsize = 10,border_color = NA,gaps_row = c(1:35),gaps_col = c(1:4),
                     angle_col = 45,fontsize_row = 14,fontsize_col = 14,main = paste0(genes),
                     filename = paste(workdir,genes,"cor.pancancer.pdf",sep = "/"),height = 10,width = 5)
  
  for(i in seq_along(dt)){
    for(gene in sig.genes){
      p1 <- ggplot(dt[[i]], aes_string(x = genes, y = gene)) +
        geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
        theme_bw() +
        theme(axis.title = element_text(size = 16),
              axis.text = element_text(size = 14),
              axis.ticks.length = unit(0.25,'cm'),
              axis.ticks = element_line(size = 1),
              panel.border = element_rect(size = 1.5),
              panel.grid = element_blank()) +
        geom_smooth(method = 'lm', se = T, color = "#BF0032", size = 1.5, fill = '#D6D6D6') +
        ggpubr::stat_cor(method = "pearson",digits = 3,size = 6)
      
      p2 <- ggExtra::ggMarginal(p1,type = "density",
                                xparams = list(binwidth = 0.1, fill = '#B3E283',size = .7),
                                yparams = list(binwidth = 0.1, fill = '#8AB6D6', size = .7))
      ggsave(paste(workdir,genes,cancer[i],paste0(cancer[i],".",genes,"_vs_",gene,".pdf"),sep = "/"),
             plot=p2,device='pdf',width=6,height=6)
    }
  }
}

s1 <- c("CD3D","CD3E","CD3G","CD8A","CD8B")

# CAF UP
CD8.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.CAF.up.gene <- CD8.CAF.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.CAF.up.gene)])
in_use <- in_use[-39]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s1,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8/CAF/up")
})

# CAF Down
CD8.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.CAF.down.gene <- CD8.CAF.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.CAF.down.gene)])
in_use <- in_use[-c(194,195)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s1,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8/CAF/down")
})


# TEC UP
CD8.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.TEC.up.gene <- CD8.TEC.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.TEC.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s1,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8/TEC/up")
})

# TEC Down
CD8.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.TEC.down.gene <- CD8.TEC.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.TEC.down.gene)])
in_use <- in_use[-c(269,270,271)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s1,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8/TEC/down")
})


# MLG UP
CD8.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.MLG.up.gene <- CD8.MLG.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.MLG.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s1,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8/Malignant.cell/up")
})

# MLG Down
CD8.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.MLG.down.gene <- CD8.MLG.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.MLG.down.gene)])
in_use <- in_use[-c(171:176)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s1,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8/Malignant.cell/down")
})


##======= 2.CD3+CD8+PD-1+Tim3+Lag3+ =======

geneSigCor <- function(pancancer,genes,sig.genes,workdir){
  if(!file.exists(paste(workdir,genes,sep="/"))){
    dir.create(paste(workdir,genes,sep="/"))
  }
  cancer <- gsub(".rnaseq","",names(pancancer))
  for(i in seq_along(cancer)){
    if(!file.exists(paste(workdir,genes,cancer[i],sep = "/"))){
      dir.create(paste(workdir,genes,cancer[i],sep="/"))
    }
  }
  dt <- lapply(pancancer,function(x){
    log2(x[,c(sig.genes,genes)] + 1)
  })
  r <- lapply(dt,function(x){
    r1 <- cor(x[,length(sig.genes)+1],y = x[,1],method = "pearson")
    r2 <- cor(x[,length(sig.genes)+1],y = x[,2],method = "pearson")
    r3 <- cor(x[,length(sig.genes)+1],y = x[,3],method = "pearson")
    r4 <- cor(x[,length(sig.genes)+1],y = x[,4],method = "pearson")
    r5 <- cor(x[,length(sig.genes)+1],y = x[,5],method = "pearson")
    r6 <- cor(x[,length(sig.genes)+1],y = x[,6],method = "pearson")
    r7 <- cor(x[,length(sig.genes)+1],y = x[,7],method = "pearson")
    r8 <- cor(x[,length(sig.genes)+1],y = x[,8],method = "pearson")
    c(r1,r2,r3,r4,r5,r6,r7,r8)
  })
  cor_result <- t(data.frame(r))
  colnames(cor_result) <- sig.genes
  rownames(cor_result) <- gsub(".rnaseq","",rownames(cor_result))
  n_patient <- lapply(pancancer,function(x) return(nrow(x)))
  for(i in seq_along(n_patient)){
    rownames(cor_result)[i] <- paste0(rownames(cor_result)[i],"(n=",n_patient[i],")")
  }
  pheatmap::pheatmap(cor_result,cluster_rows = F,cluster_cols = F,display_numbers = T,legend = F,
                     fontsize = 10,border_color = NA,gaps_row = c(1:35),gaps_col = c(1:4),
                     angle_col = 45,fontsize_row = 14,fontsize_col = 14,main = paste0(genes),
                     filename = paste(workdir,genes,"cor.pancancer.pdf",sep = "/"),height = 10,width = 5)
  
  for(i in seq_along(dt)){
    for(gene in sig.genes){
      p1 <- ggplot(dt[[i]], aes_string(x = genes, y = gene)) +
        geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
        theme_bw() +
        theme(axis.title = element_text(size = 16),
              axis.text = element_text(size = 14),
              axis.ticks.length = unit(0.25,'cm'),
              axis.ticks = element_line(size = 1),
              panel.border = element_rect(size = 1.5),
              panel.grid = element_blank()) +
        geom_smooth(method = 'lm', se = T, color = "#BF0032", size = 1.5, fill = '#D6D6D6') +
        ggpubr::stat_cor(method = "pearson",digits = 3,size = 6)
      
      p2 <- ggExtra::ggMarginal(p1,type = "density",
                                xparams = list(binwidth = 0.1, fill = '#B3E283',size = .7),
                                yparams = list(binwidth = 0.1, fill = '#8AB6D6', size = .7))
      ggsave(paste(workdir,genes,cancer[i],paste0(cancer[i],".",genes,"_vs_",gene,".pdf"),sep = "/"),
             plot=p2,device='pdf',width=6,height=6)
    }
  }
}

s2 <- c("CD3D","CD3E","CD3G","CD8A","CD8B","PDCD1","HAVCR2","LAG3")

# CAF UP
CD8.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.CAF.up.gene <- CD8.CAF.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.CAF.up.gene)])
in_use <- in_use[-39]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s2,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_PD1_TIM3_LAG3/CAF/up")
})

# CAF Down
CD8.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.CAF.down.gene <- CD8.CAF.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.CAF.down.gene)])
in_use <- in_use[-c(194,195)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s2,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_PD1_TIM3_LAG3/CAF/down")
})


# TEC UP
CD8.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.TEC.up.gene <- CD8.TEC.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.TEC.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s2,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_PD1_TIM3_LAG3/TEC/up")
})

# TEC Down
CD8.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.TEC.down.gene <- CD8.TEC.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.TEC.down.gene)])
in_use <- in_use[-c(269,270,271)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s2,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_PD1_TIM3_LAG3/TEC/down")
})


# MLG UP
CD8.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.MLG.up.gene <- CD8.MLG.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.MLG.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s2,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_PD1_TIM3_LAG3/Malignant.cell/up")
})

# MLG Down
CD8.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.MLG.down.gene <- CD8.MLG.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.MLG.down.gene)])
in_use <- in_use[-c(171:176)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s2,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_PD1_TIM3_LAG3/Malignant.cell/down")
})


##======= 3.CD3+CD8+IFNγ =======

geneSigCor <- function(pancancer,genes,sig.genes,workdir){
  if(!file.exists(paste(workdir,genes,sep="/"))){
    dir.create(paste(workdir,genes,sep="/"))
  }
  cancer <- gsub(".rnaseq","",names(pancancer))
  for(i in seq_along(cancer)){
    if(!file.exists(paste(workdir,genes,cancer[i],sep = "/"))){
      dir.create(paste(workdir,genes,cancer[i],sep="/"))
    }
  }
  dt <- lapply(pancancer,function(x){
    log2(x[,c(sig.genes,genes)] + 1)
  })
  r <- lapply(dt,function(x){
    r1 <- cor(x[,length(sig.genes)+1],y = x[,1],method = "pearson")
    r2 <- cor(x[,length(sig.genes)+1],y = x[,2],method = "pearson")
    r3 <- cor(x[,length(sig.genes)+1],y = x[,3],method = "pearson")
    r4 <- cor(x[,length(sig.genes)+1],y = x[,4],method = "pearson")
    r5 <- cor(x[,length(sig.genes)+1],y = x[,5],method = "pearson")
    r6 <- cor(x[,length(sig.genes)+1],y = x[,6],method = "pearson")
    c(r1,r2,r3,r4,r5,r6)
  })
  cor_result <- t(data.frame(r))
  colnames(cor_result) <- sig.genes
  rownames(cor_result) <- gsub(".rnaseq","",rownames(cor_result))
  n_patient <- lapply(pancancer,function(x) return(nrow(x)))
  for(i in seq_along(n_patient)){
    rownames(cor_result)[i] <- paste0(rownames(cor_result)[i],"(n=",n_patient[i],")")
  }
  pheatmap::pheatmap(cor_result,cluster_rows = F,cluster_cols = F,display_numbers = T,legend = F,
                     fontsize = 10,border_color = NA,gaps_row = c(1:35),gaps_col = c(1:4),
                     angle_col = 45,fontsize_row = 14,fontsize_col = 14,main = paste0(genes),
                     filename = paste(workdir,genes,"cor.pancancer.pdf",sep = "/"),height = 10,width = 5)
  
  for(i in seq_along(dt)){
    for(gene in sig.genes){
      p1 <- ggplot(dt[[i]], aes_string(x = genes, y = gene)) +
        geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
        theme_bw() +
        theme(axis.title = element_text(size = 16),
              axis.text = element_text(size = 14),
              axis.ticks.length = unit(0.25,'cm'),
              axis.ticks = element_line(size = 1),
              panel.border = element_rect(size = 1.5),
              panel.grid = element_blank()) +
        geom_smooth(method = 'lm', se = T, color = "#BF0032", size = 1.5, fill = '#D6D6D6') +
        ggpubr::stat_cor(method = "pearson",digits = 3,size = 6)
      
      p2 <- ggExtra::ggMarginal(p1,type = "density",
                                xparams = list(binwidth = 0.1, fill = '#B3E283',size = .7),
                                yparams = list(binwidth = 0.1, fill = '#8AB6D6', size = .7))
      ggsave(paste(workdir,genes,cancer[i],paste0(cancer[i],".",genes,"_vs_",gene,".pdf"),sep = "/"),
             plot=p2,device='pdf',width=6,height=6)
    }
  }
}

s3 <- c("CD3D","CD3E","CD3G","CD8A","CD8B","INFG")

# CAF UP
CD8.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.CAF.up.gene <- CD8.CAF.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.CAF.up.gene)])
in_use <- in_use[-39]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s3,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_IFNγ/CAF/up")
})

# CAF Down
CD8.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.CAF.down.gene <- CD8.CAF.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.CAF.down.gene)])
in_use <- in_use[-c(194,195)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s3,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_IFNγ/CAF/down")
})


# TEC UP
CD8.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.TEC.up.gene <- CD8.TEC.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.TEC.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s3,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_IFNγ/TEC/up")
})

# TEC Down
CD8.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.TEC.down.gene <- CD8.TEC.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.TEC.down.gene)])
in_use <- in_use[-c(269,270,271)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s3,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_IFNγ/TEC/down")
})


# MLG UP
CD8.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.MLG.up.gene <- CD8.MLG.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.MLG.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s3,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_IFNγ/Malignant.cell/up")
})

# MLG Down
CD8.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.MLG.down.gene <- CD8.MLG.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.MLG.down.gene)])
in_use <- in_use[-c(171:176)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s3,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_IFNγ/Malignant.cell/down")
})



##======= 4.CD3+CD8+GranzymeB =======

geneSigCor <- function(pancancer,genes,sig.genes,workdir){
  if(!file.exists(paste(workdir,genes,sep="/"))){
    dir.create(paste(workdir,genes,sep="/"))
  }
  cancer <- gsub(".rnaseq","",names(pancancer))
  for(i in seq_along(cancer)){
    if(!file.exists(paste(workdir,genes,cancer[i],sep = "/"))){
      dir.create(paste(workdir,genes,cancer[i],sep="/"))
    }
  }
  dt <- lapply(pancancer,function(x){
    log2(x[,c(sig.genes,genes)] + 1)
  })
  r <- lapply(dt,function(x){
    r1 <- cor(x[,length(sig.genes)+1],y = x[,1],method = "pearson")
    r2 <- cor(x[,length(sig.genes)+1],y = x[,2],method = "pearson")
    r3 <- cor(x[,length(sig.genes)+1],y = x[,3],method = "pearson")
    r4 <- cor(x[,length(sig.genes)+1],y = x[,4],method = "pearson")
    r5 <- cor(x[,length(sig.genes)+1],y = x[,5],method = "pearson")
    r6 <- cor(x[,length(sig.genes)+1],y = x[,6],method = "pearson")
    c(r1,r2,r3,r4,r5,r6)
  })
  cor_result <- t(data.frame(r))
  colnames(cor_result) <- sig.genes
  rownames(cor_result) <- gsub(".rnaseq","",rownames(cor_result))
  n_patient <- lapply(pancancer,function(x) return(nrow(x)))
  for(i in seq_along(n_patient)){
    rownames(cor_result)[i] <- paste0(rownames(cor_result)[i],"(n=",n_patient[i],")")
  }
  pheatmap::pheatmap(cor_result,cluster_rows = F,cluster_cols = F,display_numbers = T,legend = F,
                     fontsize = 10,border_color = NA,gaps_row = c(1:35),gaps_col = c(1:4),
                     angle_col = 45,fontsize_row = 14,fontsize_col = 14,main = paste0(genes),
                     filename = paste(workdir,genes,"cor.pancancer.pdf",sep = "/"),height = 10,width = 5)
  
  for(i in seq_along(dt)){
    for(gene in sig.genes){
      p1 <- ggplot(dt[[i]], aes_string(x = genes, y = gene)) +
        geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
        theme_bw() +
        theme(axis.title = element_text(size = 16),
              axis.text = element_text(size = 14),
              axis.ticks.length = unit(0.25,'cm'),
              axis.ticks = element_line(size = 1),
              panel.border = element_rect(size = 1.5),
              panel.grid = element_blank()) +
        geom_smooth(method = 'lm', se = T, color = "#BF0032", size = 1.5, fill = '#D6D6D6') +
        ggpubr::stat_cor(method = "pearson",digits = 3,size = 6)
      
      p2 <- ggExtra::ggMarginal(p1,type = "density",
                                xparams = list(binwidth = 0.1, fill = '#B3E283',size = .7),
                                yparams = list(binwidth = 0.1, fill = '#8AB6D6', size = .7))
      ggsave(paste(workdir,genes,cancer[i],paste0(cancer[i],".",genes,"_vs_",gene,".pdf"),sep = "/"),
             plot=p2,device='pdf',width=6,height=6)
    }
  }
}

s4 <- c("CD3D","CD3E","CD3G","CD8A","CD8B","GZMB")

# CAF UP
CD8.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.CAF.up.gene <- CD8.CAF.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.CAF.up.gene)])
in_use <- in_use[-39]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s4,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_GranzymeB/CAF/up")
})

# CAF Down
CD8.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.CAF.down.gene <- CD8.CAF.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.CAF.down.gene)])
in_use <- in_use[-c(194,195)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s4,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_GranzymeB/CAF/down")
})


# TEC UP
CD8.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.TEC.up.gene <- CD8.TEC.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.TEC.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s4,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_GranzymeB/TEC/up")
})

# TEC Down
CD8.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.TEC.down.gene <- CD8.TEC.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.TEC.down.gene)])
in_use <- in_use[-c(269,270,271)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s4,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_GranzymeB/TEC/down")
})


# MLG UP
CD8.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.MLG.up.gene <- CD8.MLG.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.MLG.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s4,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_GranzymeB/Malignant.cell/up")
})

# MLG Down
CD8.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.MLG.down.gene <- CD8.MLG.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.MLG.down.gene)])
in_use <- in_use[-c(171:176)]

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s4,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_GranzymeB/Malignant.cell/down")
})


##======= 5.CD3+CD8+TNFα =======

geneSigCor <- function(pancancer,genes,sig.genes,workdir){
  if(!file.exists(paste(workdir,genes,sep="/"))){
    dir.create(paste(workdir,genes,sep="/"))
  }
  cancer <- gsub(".rnaseq","",names(pancancer))
  for(i in seq_along(cancer)){
    if(!file.exists(paste(workdir,genes,cancer[i],sep = "/"))){
      dir.create(paste(workdir,genes,cancer[i],sep="/"))
    }
  }
  dt <- lapply(pancancer,function(x){
    log2(x[,c(sig.genes,genes)] + 1)
  })
  r <- lapply(dt,function(x){
    r1 <- cor(x[,length(sig.genes)+1],y = x[,1],method = "pearson")
    r2 <- cor(x[,length(sig.genes)+1],y = x[,2],method = "pearson")
    r3 <- cor(x[,length(sig.genes)+1],y = x[,3],method = "pearson")
    r4 <- cor(x[,length(sig.genes)+1],y = x[,4],method = "pearson")
    r5 <- cor(x[,length(sig.genes)+1],y = x[,5],method = "pearson")
    r6 <- cor(x[,length(sig.genes)+1],y = x[,6],method = "pearson")
    c(r1,r2,r3,r4,r5,r6)
  })
  cor_result <- t(data.frame(r))
  colnames(cor_result) <- sig.genes
  rownames(cor_result) <- gsub(".rnaseq","",rownames(cor_result))
  n_patient <- lapply(pancancer,function(x) return(nrow(x)))
  for(i in seq_along(n_patient)){
    rownames(cor_result)[i] <- paste0(rownames(cor_result)[i],"(n=",n_patient[i],")")
  }
  pheatmap::pheatmap(cor_result,cluster_rows = F,cluster_cols = F,display_numbers = T,legend = F,
                     fontsize = 10,border_color = NA,gaps_row = c(1:35),gaps_col = c(1:4),
                     angle_col = 45,fontsize_row = 14,fontsize_col = 14,main = paste0(genes),
                     filename = paste(workdir,genes,"cor.pancancer.pdf",sep = "/"),height = 10,width = 5)
  
  for(i in seq_along(dt)){
    for(gene in sig.genes){
      p1 <- ggplot(dt[[i]], aes_string(x = genes, y = gene)) +
        geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
        theme_bw() +
        theme(axis.title = element_text(size = 16),
              axis.text = element_text(size = 14),
              axis.ticks.length = unit(0.25,'cm'),
              axis.ticks = element_line(size = 1),
              panel.border = element_rect(size = 1.5),
              panel.grid = element_blank()) +
        geom_smooth(method = 'lm', se = T, color = "#BF0032", size = 1.5, fill = '#D6D6D6') +
        ggpubr::stat_cor(method = "pearson",digits = 3,size = 6)
      
      p2 <- ggExtra::ggMarginal(p1,type = "density",
                                xparams = list(binwidth = 0.1, fill = '#B3E283',size = .7),
                                yparams = list(binwidth = 0.1, fill = '#8AB6D6', size = .7))
      ggsave(paste(workdir,genes,cancer[i],paste0(cancer[i],".",genes,"_vs_",gene,".pdf"),sep = "/"),
             plot=p2,device='pdf',width=6,height=6)
    }
  }
}

s5 <- c("CD3D","CD3E","CD3G","CD8A","CD8B","TNF")

# CAF UP
CD8.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.CAF.up.gene <- CD8.CAF.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.CAF.up.gene)])
in_use <- in_use[-39]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s5,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_TNFα/CAF/up")
})

# CAF Down
CD8.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.CAF.down.gene <- CD8.CAF.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.CAF.down.gene)])
in_use <- in_use[-c(194,195)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s5,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_TNFα/CAF/down")
})


# TEC UP
CD8.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.TEC.up.gene <- CD8.TEC.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.TEC.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s5,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_TNFα/TEC/up")
})

# TEC Down
CD8.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/TEC/TEC.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.TEC.down.gene <- CD8.TEC.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.TEC.down.gene)])
in_use <- in_use[-c(269,270,271)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s5,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_TNFα/TEC/down")
})


# MLG UP
CD8.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 1)
CD8.MLG.up.gene <- CD8.MLG.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.MLG.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s5,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_TNFα/Malignant.cell/up")
})

# MLG Down
CD8.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD8T.hi_vs_lo.xlsx",sheet = 2)
CD8.MLG.down.gene <- CD8.MLG.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD8.MLG.down.gene)])
in_use <- in_use[-c(171:176)]

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s5,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD8/CD3_CD8_TNFα/Malignant.cell/down")
})


##>>>>>>>>>>>>>>>>>>>>>>> 2.CD4+T <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 1.CD3+CD4+;
# 2.CD3+CD4+PD-1+Tim3+Lag3+；
# 3.CD3+CD4+IFNg；
# 4.CD3+CD4+GranzymeB；
# 5.CD3+CD4+TNFa；

##======= 1.CD3+CD4+ ======= 

geneSigCor <- function(pancancer,genes,sig.genes,workdir){
  if(!file.exists(paste(workdir,genes,sep="/"))){
    dir.create(paste(workdir,genes,sep="/"))
  }
  cancer <- gsub(".rnaseq","",names(pancancer))
  for(i in seq_along(cancer)){
    if(!file.exists(paste(workdir,genes,cancer[i],sep = "/"))){
      dir.create(paste(workdir,genes,cancer[i],sep="/"))
    }
  }
  dt <- lapply(pancancer,function(x){
    log2(x[,c(sig.genes,genes)] + 1)
  })
  r <- lapply(dt,function(x){
    r1 <- cor(x[,length(sig.genes)+1],y = x[,1],method = "pearson")
    r2 <- cor(x[,length(sig.genes)+1],y = x[,2],method = "pearson")
    r3 <- cor(x[,length(sig.genes)+1],y = x[,3],method = "pearson")
    r4 <- cor(x[,length(sig.genes)+1],y = x[,4],method = "pearson")
    c(r1,r2,r3,r4)
  })
  cor_result <- t(data.frame(r))
  colnames(cor_result) <- sig.genes
  rownames(cor_result) <- gsub(".rnaseq","",rownames(cor_result))
  n_patient <- lapply(pancancer,function(x) return(nrow(x)))
  for(i in seq_along(n_patient)){
    rownames(cor_result)[i] <- paste0(rownames(cor_result)[i],"(n=",n_patient[i],")")
  }
  pheatmap::pheatmap(cor_result,cluster_rows = F,cluster_cols = F,display_numbers = T,legend = F,
                     fontsize = 10,border_color = NA,gaps_row = c(1:35),gaps_col = c(1:4),
                     angle_col = 45,fontsize_row = 14,fontsize_col = 14,main = paste0(genes),
                     filename = paste(workdir,genes,"cor.pancancer.pdf",sep = "/"),height = 10,width = 5)
  
  for(i in seq_along(dt)){
    for(gene in sig.genes){
      p1 <- ggplot(dt[[i]], aes_string(x = genes, y = gene)) +
        geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
        theme_bw() +
        theme(axis.title = element_text(size = 16),
              axis.text = element_text(size = 14),
              axis.ticks.length = unit(0.25,'cm'),
              axis.ticks = element_line(size = 1),
              panel.border = element_rect(size = 1.5),
              panel.grid = element_blank()) +
        geom_smooth(method = 'lm', se = T, color = "#BF0032", size = 1.5, fill = '#D6D6D6') +
        ggpubr::stat_cor(method = "pearson",digits = 3,size = 6)
      
      p2 <- ggExtra::ggMarginal(p1,type = "density",
                                xparams = list(binwidth = 0.1, fill = '#B3E283',size = .7),
                                yparams = list(binwidth = 0.1, fill = '#8AB6D6', size = .7))
      ggsave(paste(workdir,genes,cancer[i],paste0(cancer[i],".",genes,"_vs_",gene,".pdf"),sep = "/"),
             plot=p2,device='pdf',width=6,height=6)
    }
  }
}

s1 <- c("CD3D","CD3E","CD3G","CD4")

# CAF UP
CD4.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.CAF.up.gene <- CD4.CAF.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.CAF.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s1,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4/CAF/up")
})

# CAF Down
CD4.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.CAF.down.gene <- CD4.CAF.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.CAF.down.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s1,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4/CAF/down")
})


# TEC UP
CD4.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.TEC.up.gene <- CD4.TEC.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.TEC.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s1,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4/TEC/up")
})

# TEC Down
CD4.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.TEC.down.gene <- CD4.TEC.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.TEC.down.gene)])
in_use <- in_use[-c(159,160,161)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s1,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4/TEC/down")
})


# MLG UP
CD4.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.MLG.up.gene <- CD4.MLG.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.MLG.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s1,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4/Malignant.cell/up")
})

# MLG Down
CD4.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.MLG.down.gene <- CD4.MLG.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.MLG.down.gene)])
in_use <- in_use[-c(220:223)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s1,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4/Malignant.cell/down")
})


##======= 2.CD3+CD4+PD-1+Tim3+Lag3+ =======

geneSigCor <- function(pancancer,genes,sig.genes,workdir){
  if(!file.exists(paste(workdir,genes,sep="/"))){
    dir.create(paste(workdir,genes,sep="/"))
  }
  cancer <- gsub(".rnaseq","",names(pancancer))
  for(i in seq_along(cancer)){
    if(!file.exists(paste(workdir,genes,cancer[i],sep = "/"))){
      dir.create(paste(workdir,genes,cancer[i],sep="/"))
    }
  }
  dt <- lapply(pancancer,function(x){
    log2(x[,c(sig.genes,genes)] + 1)
  })
  r <- lapply(dt,function(x){
    r1 <- cor(x[,length(sig.genes)+1],y = x[,1],method = "pearson")
    r2 <- cor(x[,length(sig.genes)+1],y = x[,2],method = "pearson")
    r3 <- cor(x[,length(sig.genes)+1],y = x[,3],method = "pearson")
    r4 <- cor(x[,length(sig.genes)+1],y = x[,4],method = "pearson")
    r5 <- cor(x[,length(sig.genes)+1],y = x[,5],method = "pearson")
    r6 <- cor(x[,length(sig.genes)+1],y = x[,6],method = "pearson")
    r7 <- cor(x[,length(sig.genes)+1],y = x[,7],method = "pearson")
    c(r1,r2,r3,r4,r5,r6,r7)
  })
  cor_result <- t(data.frame(r))
  colnames(cor_result) <- sig.genes
  rownames(cor_result) <- gsub(".rnaseq","",rownames(cor_result))
  n_patient <- lapply(pancancer,function(x) return(nrow(x)))
  for(i in seq_along(n_patient)){
    rownames(cor_result)[i] <- paste0(rownames(cor_result)[i],"(n=",n_patient[i],")")
  }
  pheatmap::pheatmap(cor_result,cluster_rows = F,cluster_cols = F,display_numbers = T,legend = F,
                     fontsize = 10,border_color = NA,gaps_row = c(1:35),gaps_col = c(1:4),
                     angle_col = 45,fontsize_row = 14,fontsize_col = 14,main = paste0(genes),
                     filename = paste(workdir,genes,"cor.pancancer.pdf",sep = "/"),height = 10,width = 5)
  
  for(i in seq_along(dt)){
    for(gene in sig.genes){
      p1 <- ggplot(dt[[i]], aes_string(x = genes, y = gene)) +
        geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
        theme_bw() +
        theme(axis.title = element_text(size = 16),
              axis.text = element_text(size = 14),
              axis.ticks.length = unit(0.25,'cm'),
              axis.ticks = element_line(size = 1),
              panel.border = element_rect(size = 1.5),
              panel.grid = element_blank()) +
        geom_smooth(method = 'lm', se = T, color = "#BF0032", size = 1.5, fill = '#D6D6D6') +
        ggpubr::stat_cor(method = "pearson",digits = 3,size = 6)
      
      p2 <- ggExtra::ggMarginal(p1,type = "density",
                                xparams = list(binwidth = 0.1, fill = '#B3E283',size = .7),
                                yparams = list(binwidth = 0.1, fill = '#8AB6D6', size = .7))
      ggsave(paste(workdir,genes,cancer[i],paste0(cancer[i],".",genes,"_vs_",gene,".pdf"),sep = "/"),
             plot=p2,device='pdf',width=6,height=6)
    }
  }
}

s2 <- c("CD3D","CD3E","CD3G","CD4","PDCD1","HAVCR2","LAG3")

# CAF UP
CD4.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.CAF.up.gene <- CD4.CAF.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.CAF.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s2,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_PD1_TIM3_LAG3/CAF/up/")
})

# CAF Down
CD4.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.CAF.down.gene <- CD4.CAF.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.CAF.down.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s2,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_PD1_TIM3_LAG3/CAF/down")
})


# TEC UP
CD4.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.TEC.up.gene <- CD4.TEC.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.TEC.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s2,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_PD1_TIM3_LAG3/TEC/up")
})

# TEC Down
CD4.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.TEC.down.gene <- CD4.TEC.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.TEC.down.gene)])
in_use <- in_use[-c(159,160,161)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s2,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_PD1_TIM3_LAG3/TEC/down")
})


# MLG UP
CD4.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.MLG.up.gene <- CD4.MLG.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.MLG.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s2,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_PD1_TIM3_LAG3/Malignant.cell/up")
})

# MLG Down
CD4.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.MLG.down.gene <- CD4.MLG.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.MLG.down.gene)])
in_use <- in_use[-c(220:223)]

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s2,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_PD1_TIM3_LAG3/Malignant.cell/down")
})


##======= 3.CD3+CD4+IFNγ =======

geneSigCor <- function(pancancer,genes,sig.genes,workdir){
  if(!file.exists(paste(workdir,genes,sep="/"))){
    dir.create(paste(workdir,genes,sep="/"))
  }
  cancer <- gsub(".rnaseq","",names(pancancer))
  for(i in seq_along(cancer)){
    if(!file.exists(paste(workdir,genes,cancer[i],sep = "/"))){
      dir.create(paste(workdir,genes,cancer[i],sep="/"))
    }
  }
  dt <- lapply(pancancer,function(x){
    log2(x[,c(sig.genes,genes)] + 1)
  })
  r <- lapply(dt,function(x){
    r1 <- cor(x[,length(sig.genes)+1],y = x[,1],method = "pearson")
    r2 <- cor(x[,length(sig.genes)+1],y = x[,2],method = "pearson")
    r3 <- cor(x[,length(sig.genes)+1],y = x[,3],method = "pearson")
    r4 <- cor(x[,length(sig.genes)+1],y = x[,4],method = "pearson")
    r5 <- cor(x[,length(sig.genes)+1],y = x[,5],method = "pearson")
    c(r1,r2,r3,r4,r5)
  })
  cor_result <- t(data.frame(r))
  colnames(cor_result) <- sig.genes
  rownames(cor_result) <- gsub(".rnaseq","",rownames(cor_result))
  n_patient <- lapply(pancancer,function(x) return(nrow(x)))
  for(i in seq_along(n_patient)){
    rownames(cor_result)[i] <- paste0(rownames(cor_result)[i],"(n=",n_patient[i],")")
  }
  pheatmap::pheatmap(cor_result,cluster_rows = F,cluster_cols = F,display_numbers = T,legend = F,
                     fontsize = 10,border_color = NA,gaps_row = c(1:35),gaps_col = c(1:4),
                     angle_col = 45,fontsize_row = 14,fontsize_col = 14,main = paste0(genes),
                     filename = paste(workdir,genes,"cor.pancancer.pdf",sep = "/"),height = 10,width = 5)
  
  for(i in seq_along(dt)){
    for(gene in sig.genes){
      p1 <- ggplot(dt[[i]], aes_string(x = genes, y = gene)) +
        geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
        theme_bw() +
        theme(axis.title = element_text(size = 16),
              axis.text = element_text(size = 14),
              axis.ticks.length = unit(0.25,'cm'),
              axis.ticks = element_line(size = 1),
              panel.border = element_rect(size = 1.5),
              panel.grid = element_blank()) +
        geom_smooth(method = 'lm', se = T, color = "#BF0032", size = 1.5, fill = '#D6D6D6') +
        ggpubr::stat_cor(method = "pearson",digits = 3,size = 6)
      
      p2 <- ggExtra::ggMarginal(p1,type = "density",
                                xparams = list(binwidth = 0.1, fill = '#B3E283',size = .7),
                                yparams = list(binwidth = 0.1, fill = '#8AB6D6', size = .7))
      ggsave(paste(workdir,genes,cancer[i],paste0(cancer[i],".",genes,"_vs_",gene,".pdf"),sep = "/"),
             plot=p2,device='pdf',width=6,height=6)
    }
  }
}

s3 <- c("CD3D","CD3E","CD3G","CD4","INFG")

# CAF UP
CD4.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.CAF.up.gene <- CD4.CAF.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.CAF.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s3,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_IFNγ/CAF/up/")
})

# CAF Down
CD4.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.CAF.down.gene <- CD4.CAF.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.CAF.down.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s3,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_IFNγ/CAF/down")
})


# TEC UP
CD4.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.TEC.up.gene <- CD4.TEC.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.TEC.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s3,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_IFNγ/TEC/up")
})

# TEC Down
CD4.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.TEC.down.gene <- CD4.TEC.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.TEC.down.gene)])
in_use <- in_use[-c(159,160,161)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s3,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_IFNγ/TEC/down")
})


# MLG UP
CD4.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.MLG.up.gene <- CD4.MLG.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.MLG.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s3,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_IFNγ/Malignant.cell/up")
})

# MLG Down
CD4.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.MLG.down.gene <- CD4.MLG.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.MLG.down.gene)])
in_use <- in_use[-c(220:223)]

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s3,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_IFNγ/Malignant.cell/down")
})


##======= 4.CD3+CD4+GranzymeB =======

geneSigCor <- function(pancancer,genes,sig.genes,workdir){
  if(!file.exists(paste(workdir,genes,sep="/"))){
    dir.create(paste(workdir,genes,sep="/"))
  }
  cancer <- gsub(".rnaseq","",names(pancancer))
  for(i in seq_along(cancer)){
    if(!file.exists(paste(workdir,genes,cancer[i],sep = "/"))){
      dir.create(paste(workdir,genes,cancer[i],sep="/"))
    }
  }
  dt <- lapply(pancancer,function(x){
    log2(x[,c(sig.genes,genes)] + 1)
  })
  r <- lapply(dt,function(x){
    r1 <- cor(x[,length(sig.genes)+1],y = x[,1],method = "pearson")
    r2 <- cor(x[,length(sig.genes)+1],y = x[,2],method = "pearson")
    r3 <- cor(x[,length(sig.genes)+1],y = x[,3],method = "pearson")
    r4 <- cor(x[,length(sig.genes)+1],y = x[,4],method = "pearson")
    r5 <- cor(x[,length(sig.genes)+1],y = x[,5],method = "pearson")
    c(r1,r2,r3,r4,r5)
  })
  cor_result <- t(data.frame(r))
  colnames(cor_result) <- sig.genes
  rownames(cor_result) <- gsub(".rnaseq","",rownames(cor_result))
  n_patient <- lapply(pancancer,function(x) return(nrow(x)))
  for(i in seq_along(n_patient)){
    rownames(cor_result)[i] <- paste0(rownames(cor_result)[i],"(n=",n_patient[i],")")
  }
  pheatmap::pheatmap(cor_result,cluster_rows = F,cluster_cols = F,display_numbers = T,legend = F,
                     fontsize = 10,border_color = NA,gaps_row = c(1:35),gaps_col = c(1:4),
                     angle_col = 45,fontsize_row = 14,fontsize_col = 14,main = paste0(genes),
                     filename = paste(workdir,genes,"cor.pancancer.pdf",sep = "/"),height = 10,width = 5)
  
  for(i in seq_along(dt)){
    for(gene in sig.genes){
      p1 <- ggplot(dt[[i]], aes_string(x = genes, y = gene)) +
        geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
        theme_bw() +
        theme(axis.title = element_text(size = 16),
              axis.text = element_text(size = 14),
              axis.ticks.length = unit(0.25,'cm'),
              axis.ticks = element_line(size = 1),
              panel.border = element_rect(size = 1.5),
              panel.grid = element_blank()) +
        geom_smooth(method = 'lm', se = T, color = "#BF0032", size = 1.5, fill = '#D6D6D6') +
        ggpubr::stat_cor(method = "pearson",digits = 3,size = 6)
      
      p2 <- ggExtra::ggMarginal(p1,type = "density",
                                xparams = list(binwidth = 0.1, fill = '#B3E283',size = .7),
                                yparams = list(binwidth = 0.1, fill = '#8AB6D6', size = .7))
      ggsave(paste(workdir,genes,cancer[i],paste0(cancer[i],".",genes,"_vs_",gene,".pdf"),sep = "/"),
             plot=p2,device='pdf',width=6,height=6)
    }
  }
}

s4 <- c("CD3D","CD3E","CD3G","CD4","GZMB")

# CAF UP
CD4.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.CAF.up.gene <- CD4.CAF.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.CAF.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s4,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_GranzymeB/CAF/up/")
})

# CAF Down
CD4.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.CAF.down.gene <- CD4.CAF.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.CAF.down.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s4,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_GranzymeB/CAF/down")
})


# TEC UP
CD4.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.TEC.up.gene <- CD4.TEC.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.TEC.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s4,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_GranzymeB/TEC/up")
})

# TEC Down
CD4.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.TEC.down.gene <- CD4.TEC.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.TEC.down.gene)])
in_use <- in_use[-c(159,160,161)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s4,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_GranzymeB/TEC/down")
})


# MLG UP
CD4.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.MLG.up.gene <- CD4.MLG.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.MLG.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s4,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_GranzymeB/Malignant.cell/up")
})

# MLG Down
CD4.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.MLG.down.gene <- CD4.MLG.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.MLG.down.gene)])
in_use <- in_use[-c(220:223)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s4,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_GranzymeB/Malignant.cell/down")
})


##======= 5.CD3+CD4+TNFα =======

geneSigCor <- function(pancancer,genes,sig.genes,workdir){
  if(!file.exists(paste(workdir,genes,sep="/"))){
    dir.create(paste(workdir,genes,sep="/"))
  }
  cancer <- gsub(".rnaseq","",names(pancancer))
  for(i in seq_along(cancer)){
    if(!file.exists(paste(workdir,genes,cancer[i],sep = "/"))){
      dir.create(paste(workdir,genes,cancer[i],sep="/"))
    }
  }
  dt <- lapply(pancancer,function(x){
    log2(x[,c(sig.genes,genes)] + 1)
  })
  r <- lapply(dt,function(x){
    r1 <- cor(x[,length(sig.genes)+1],y = x[,1],method = "pearson")
    r2 <- cor(x[,length(sig.genes)+1],y = x[,2],method = "pearson")
    r3 <- cor(x[,length(sig.genes)+1],y = x[,3],method = "pearson")
    r4 <- cor(x[,length(sig.genes)+1],y = x[,4],method = "pearson")
    r5 <- cor(x[,length(sig.genes)+1],y = x[,5],method = "pearson")
    c(r1,r2,r3,r4,r5)
  })
  cor_result <- t(data.frame(r))
  colnames(cor_result) <- sig.genes
  rownames(cor_result) <- gsub(".rnaseq","",rownames(cor_result))
  n_patient <- lapply(pancancer,function(x) return(nrow(x)))
  for(i in seq_along(n_patient)){
    rownames(cor_result)[i] <- paste0(rownames(cor_result)[i],"(n=",n_patient[i],")")
  }
  pheatmap::pheatmap(cor_result,cluster_rows = F,cluster_cols = F,display_numbers = T,legend = F,
                     fontsize = 10,border_color = NA,gaps_row = c(1:35),gaps_col = c(1:4),
                     angle_col = 45,fontsize_row = 14,fontsize_col = 14,main = paste0(genes),
                     filename = paste(workdir,genes,"cor.pancancer.pdf",sep = "/"),height = 10,width = 5)
  
  for(i in seq_along(dt)){
    for(gene in sig.genes){
      p1 <- ggplot(dt[[i]], aes_string(x = genes, y = gene)) +
        geom_point(size = 1.5, color = '#F9B208', alpha=.7) +
        theme_bw() +
        theme(axis.title = element_text(size = 16),
              axis.text = element_text(size = 14),
              axis.ticks.length = unit(0.25,'cm'),
              axis.ticks = element_line(size = 1),
              panel.border = element_rect(size = 1.5),
              panel.grid = element_blank()) +
        geom_smooth(method = 'lm', se = T, color = "#BF0032", size = 1.5, fill = '#D6D6D6') +
        ggpubr::stat_cor(method = "pearson",digits = 3,size = 6)
      
      p2 <- ggExtra::ggMarginal(p1,type = "density",
                                xparams = list(binwidth = 0.1, fill = '#B3E283',size = .7),
                                yparams = list(binwidth = 0.1, fill = '#8AB6D6', size = .7))
      ggsave(paste(workdir,genes,cancer[i],paste0(cancer[i],".",genes,"_vs_",gene,".pdf"),sep = "/"),
             plot=p2,device='pdf',width=6,height=6)
    }
  }
}

s5 <- c("CD3D","CD3E","CD3G","CD4","TNF")

# CAF UP
CD4.CAF.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.CAF.up.gene <- CD4.CAF.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.CAF.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s5,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_TNFα/CAF/up/")
})

# CAF Down
CD4.CAF.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/CAF/CAF.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.CAF.down.gene <- CD4.CAF.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.CAF.down.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s5,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_TNFα/CAF/down")
})


# TEC UP
CD4.TEC.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.TEC.up.gene <- CD4.TEC.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.TEC.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s5,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_TNFα/TEC/up")
})

# TEC Down
CD4.TEC.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/TEC/TEC.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.TEC.down.gene <- CD4.TEC.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.TEC.down.gene)])
in_use <- in_use[-c(159,160,161)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s5,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_TNFα/TEC/down")
})


# MLG UP
CD4.MLG.up <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 1)
CD4.MLG.up.gene <- CD4.MLG.up$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.MLG.up.gene)])

sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s5,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_TNFα/Malignant.cell/up")
})

# MLG Down
CD4.MLG.down <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD4T/CD45_removal/DE/Malignant_Cell/Malignant.Cell.CD4T.hi_vs_lo.xlsx",sheet = 2)
CD4.MLG.down.gene <- CD4.MLG.down$Genes
in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% CD4.MLG.down.gene)])
in_use <- in_use[-c(220:223)]
sapply(in_use, function(x){
  geneSigCor(pancancer = TCGA_36cancers.rnaseq,genes = x,sig.genes = s5,
             workdir = "result/without_PVTT_MLN/tumor/validation/TCGA/signature/CD4/CD3_CD4_TNFα/Malignant.cell/down")
})


## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> For singlecell <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#======== Correlation between CD8 T cell & CD8 Down genes

if(F){
  tumor <- readRDS("rds/without_PVTT_MLN/tumor/3.tumor_combine_T.IR_2.rds")
  DefaultAssay(tumor) <- "RNA"
  
  all.exp <- GetAssayData(tumor,assay = "RNA",slot = "data")
  
  all.exp.transfer <- t(as.matrix(all.exp))
  
  ##========= CD8 T 
  
  CAF <- readxl::read_excel("result/without_PVTT_MLN/tumor/CD8T/CD45_removal/DE/CAF/CAF.CD8T.hi_vs_lo.xlsx",sheet = 2)
  # 一般选用normalized之后的数据
  CD8.exp <- GetAssayData(tumor,assay = "RNA",slot = "data")[c("CD8A","CD8B"),]
  CD8.exp <- t(data.frame(CD8.exp))
  CAF.gene.exp <- GetAssayData(tumor,assay = "RNA",slot = "data")[CAF$Genes,]
  CAF.gene.exp <- t(data.frame(CAF.gene.exp))
  
  all.exp <- cbind(CD8.exp,CAF.gene.exp)
  
  all.exp <- data.frame(all.exp)
  
  colnames(all.exp) <- gsub("-","\\.",colnames(all.exp))
  
  target_gene <- "CD8A"
  
  cor_list <- list()
  for(i in colnames(all.exp)){
    # 提取目的基因的表达量
    tar <- all.exp[,target_gene]
    # 相关性检验
    cor_res <- cor(x = tar, y = all.exp[,i],method = "pearson")
    # pvalue
    cor_pval <- cor.test(x = tar, y = all.exp[,i])$p.value
    # 合并结果
    final_res <- data.frame(targ_genename = target_gene,gene_name = i,
                            cor_results = cor_res, cor_pvalue = cor_pval)
    cor_list[[i]] <- final_res
  }
  gene_corres <- do.call('rbind',cor_list)
  
  pos_cor <- gene_corres[gene_corres$cor_results > 0,]
  neg_cor <- gene_corres[gene_corres$cor_results < 0,]
  
  pos_gene <- pos_cor$gene_name
  pos_gene <- gsub("-","\\.",pos_gene)
  
  plotlist <- list()
  
  for(i in pos_gene){
    p <- ggplot(all.exp,aes_string(x=i,y=target_gene)) +
      geom_point(size=2,color='#EC0101',alpha=0.5) + theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = 'lm',se = T,color='#F9B208',size=1.5,fill='#FEA82F') +
      stat_cor(method = "pearson",digits = 3,size=6)
    
    p1 <- ggExtra::ggMarginal(p, type = "densigram",
                              xparams = list(binwidth = 0.1, fill="#B3E283",size=.7),
                              yparams = list(binwidth=0.1,fill='#8AB6D6',size=.7))
    plotlist[[i]] <- p1
  }
  
  for(n in seq(3,84,3)){
    print(n)
    nstart = n - 2
    tmp_list <- plotlist[nstart:n]
    tmp_plot <- cowplot::plot_grid(plotlist = tmp_list,ncol = 3,align = 'hv')
    ggsave(tmp_plot,filename = paste0('result/without_PVTT_MLN/tumor/validation/Self/CD8T/CAF/up/CD8A/',nstart,n,'_corr.pdf'),
           width = 12,height = 4)
  }
  
  p <- ggplot(all.exp,aes_string(x="MIR4435.2HG",y=target_gene)) +
    geom_point(size=2,color='#EC0101',alpha=0.5) + theme_bw() +
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.ticks.length = unit(0.25,'cm'),
          axis.ticks = element_line(size = 1),
          panel.border = element_rect(size = 1.5),
          panel.grid = element_blank()) +
    geom_smooth(method = 'lm',se = T,color='#F9B208',size=1.5,fill='#FEA82F') +
    stat_cor(method = "pearson",digits = 3,size=6)
  
  p1 <- ggExtra::ggMarginal(p, type = "densigram",
                            xparams = list(binwidth = 0.1, fill="#B3E283",size=.7),
                            yparams = list(binwidth=0.1,fill='#8AB6D6',size=.7))
}

library(Seurat)
library(ggplot2)
library(dplyr)

cancerSEA.sce <- readRDS("refdb/CancerSEA/CancerSEA.sce.rds")

cancerSEA.sce <- FindVariableFeatures(cancerSEA.sce,selection.method = "vst",nfeatures = 2000)
cancerSEA.sce <- ScaleData(cancerSEA.sce,features = rownames(cancerSEA.sce))
cancerSEA.sce <- RunPCA(cancerSEA.sce,features = VariableFeatures(cancerSEA.sce))

p1 <- DimPlot(cancerSEA.sce,reduction = "pca",group.by = "CancerType",cols = rainbow(length(levels(as.factor(cancerSEA.sce$CancerType)))))
ggsave("cancerSEA_pca.pdf",p1,device = "pdf",width = 8,height = 6)
p2 <- ElbowPlot(cancerSEA.sce,ndims = 50,reduction = "pca")
ggsave("cancerSEA_elbow.pdf",p2,device = "pdf",width = 6,height = 4)

pc.num = 1:50

cancerSEA.sce <- FindNeighbors(cancerSEA.sce,dims = pc.num) %>%
  FindClusters(resolution = 0.5)

cancerSEA.sce <- RunUMAP(cancerSEA.sce,reduction = "pca",dims = pc.num)

p1 <- DimPlot(cancerSEA.sce,pt.size = 0.9)
p2 <- DimPlot(cancerSEA.sce,group.by="CancerType",pt.size = 0.9) + 
  theme(plot.title = element_blank())
p3 <- DimPlot(cancerSEA.sce, group.by = "Samples",pt.size = 0.9) + 
  theme(plot.title = element_blank())
ggsave("p1.pdf",p1,width = 12,height = 6)
ggsave("p2.pdf",p2,width = 12,height = 6)
ggsave("p3.pdf",p3,width = 12,height = 6)

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EXTENSION: correlation bewteen All human genes and CD8A <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

genes <- readxl::read_excel("result/without_PVTT_MLN/tumor/validation/TCGA/signature/Human-protein-coding gene-19195.xlsx")

genes <- genes$symbol

load("refdb/TCGA_36cancers.rnaseq.RData")

in_use <- colnames(TCGA_36cancers.rnaseq[[1]][which(colnames(TCGA_36cancers.rnaseq[[1]]) %in% genes)])

cancer <- gsub(".rnaseq","",names(TCGA_36cancers.rnaseq))
for(i in seq_along(cancer)){
  if(!file.exists(paste0("AllhumanGenes/",cancer[i]))){
    dir.create(paste0("AllhumanGenes/",cancer[i]))
  }
}


dt <- lapply(TCGA_36cancers.rnaseq, function(x){
  log2(x[,c('CD8A',in_use)] + 1)
})


for(i in seq_along(dt)){
  for(j in 2:16505){
    p1 <- ggplot(dt[[i]],aes_string(x = colnames(dt[[i]])[j], y =  colnames(dt[[i]])[1])) +
      geom_point(size = 1.5, color = '#F9B208', alpha = .7) +
      theme_bw() +
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 14),
            axis.ticks.length = unit(0.25,'cm'),
            axis.ticks = element_line(size = 1),
            panel.border = element_rect(size = 1.5),
            panel.grid = element_blank()) +
      geom_smooth(method = 'lm', se = T, color = "#BF0032", size = 1.5, fill = '#D6D6D6') +
      ggpubr::stat_cor(method = "pearson",digits = 3,size = 6)
    p2 <-  p2 <- ggExtra::ggMarginal(p1,type = "density",
                                     xparams = list(binwidth = 0.1, fill = '#B3E283',size = .7),
                                     yparams = list(binwidth = 0.1, fill = '#8AB6D6', size = .7))
    ggsave(paste("AllhumanGenes",cancer[i],paste0(colnames(dt[[i]])[j],"_vs_CD8A",".pdf"),sep = "/"),
           plot=p2,device='pdf',width=6,height=6)
  }
}