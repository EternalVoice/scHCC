
wkd="/home/tongqiang/project/shengzhou/liver_cancer/results/new"
setwd(wkd)
if(!file.exists("GSEA")) {dir.create("GSEA")}
if(!file.exists("GSEA/CD3T")) {dir.create("GSEA/CD3T")}
if(!file.exists("GSEA/CD4T")) {dir.create("GSEA/CD4T")}
if(!file.exists("GSEA/CD8T")) {dir.create("GSEA/CD8T")}
if(!file.exists("GSEA/CD3T/input")) {dir.create("GSEA/CD3T/input")}
if(!file.exists("GSEA/CD3T/output")) {dir.create("GSEA/CD3T/output")}
if(!file.exists("GSEA/CD3T/output/high")) {dir.create("GSEA/CD3T/output/high")}
if(!file.exists("GSEA/CD3T/output/high/up")) {dir.create("GSEA/CD3T/output/high/up")}
if(!file.exists("GSEA/CD3T/output/high/down")) {dir.create("GSEA/CD3T/output/high/down")}
if(!file.exists("GSEA/CD3T/output/low")) {dir.create("GSEA/CD3T/output/low")}
if(!file.exists("GSEA/CD3T/output/low/up")) {dir.create("GSEA/CD3T/output/low/up")}
if(!file.exists("GSEA/CD3T/output/low/down")) {dir.create("GSEA/CD3T/output/low/down")}
if(!file.exists("GSEA/CD4T/input")) {dir.create("GSEA/CD4T/input")}
if(!file.exists("GSEA/CD4T/output")) {dir.create("GSEA/CD4T/output")}
if(!file.exists("GSEA/CD4T/output/high")) {dir.create("GSEA/CD4T/output/high")}
if(!file.exists("GSEA/CD4T/output/high/up")) {dir.create("GSEA/CD4T/output/high/up")}
if(!file.exists("GSEA/CD4T/output/high/down")) {dir.create("GSEA/CD4T/output/high/down")}
if(!file.exists("GSEA/CD4T/output/low")) {dir.create("GSEA/CD4T/output/low")}
if(!file.exists("GSEA/CD4T/output/low/up")) {dir.create("GSEA/CD4T/output/low/up")}
if(!file.exists("GSEA/CD4T/output/low/down")) {dir.create("GSEA/CD4T/output/low/down")}
if(!file.exists("GSEA/CD8T/input")) {dir.create("GSEA/CD8T/input")}
if(!file.exists("GSEA/CD8T/output")) {dir.create("GSEA/CD8T/output")}
if(!file.exists("GSEA/CD8T/output/high")) {dir.create("GSEA/CD8T/output/high")}
if(!file.exists("GSEA/CD8T/output/high/up")) {dir.create("GSEA/CD8T/output/high/up")}
if(!file.exists("GSEA/CD8T/output/high/down")) {dir.create("GSEA/CD8T/output/high/down")}
if(!file.exists("GSEA/CD8T/output/low")) {dir.create("GSEA/CD8T/output/low")}
if(!file.exists("GSEA/CD8T/output/low/up")) {dir.create("GSEA/CD8T/output/low/up")}
if(!file.exists("GSEA/CD8T/output/low/down")) {dir.create("GSEA/CD8T/output/low/down")}

library(future)
plan("multiprocess", workers = 32)

library(Seurat)

load("RData/3.tumor_Tcell_tumor_marked.RData")
load("RData/2.Tumor_vs_nonTumor_up_down_Markers.RData")

# prepare GSEA inputs
# CD8T_low.up
expr <- GetAssayData(tumor_CD8_low, slot = "counts")[rownames(CD8_low.up),]
expr <- data.frame(NAME = rownames(expr), Description = rep('na',nrow(expr)), expr, stringAsFactors = F)
write('#1.2', 'GSEA/CD8T/input/CD8_low.up.expr.gct', ncolumns =1)
write(c(nrow(expr), (ncol(expr) - 2)), 'GSEA/CD8T/input/CD8_low.up.expr.gct', ncolumns = 2, append = T, sep = "\t")
write.table(expr, 'GSEA/CD8T/input/CD8_low.up.expr.gct', append=T, sep="\t", quote = F, row.names = F)
line.1 <- c((ncol(expr) - 2), 2, 1)
tmp <- table(tumor_CD8_low@meta.data$Class)
line.2 <- c('#', names(tmp))
line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
write(line.1, "GSEA/CD8T/input/CD8_low.up.group.cls", ncolumns = length(line.1), append = T, sep = "\t")
write(line.2, "GSEA/CD8T/input/CD8_low.up.group.cls", ncolumns = length(line.2), append = T, sep = "\t")
write(line.3, "GSEA/CD8T/input/CD8_low.up.group.cls", ncolumns = length(line.3), append = T, sep = "\t")

# CD8T_low down
expr <- GetAssayData(tumor_CD8_low, slot = "counts")[rownames(CD8_low.down),]
expr <- data.frame(NAME = rownames(expr), Description = rep('na',nrow(expr)), expr, stringAsFactors = F)
write('#1.2', 'GSEA/CD8T/input/CD8_low.down.expr.gct', ncolumns =1)
write(c(nrow(expr), (ncol(expr) - 2)), 'GSEA/CD8T/input/CD8_low.down.expr.gct', ncolumns = 2, append = T, sep = "\t")
write.table(expr, 'GSEA/CD8T/input/CD8_low.down.expr.gct', append=T, sep="\t", quote = F, row.names = F)
line.1 <- c((ncol(expr) - 2), 2, 1)
tmp <- table(tumor_CD8_low@meta.data$Class)
line.2 <- c('#', names(tmp))
line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
write(line.1, "GSEA/CD8T/input/CD8_low.down.group.cls", ncolumns = length(line.1), append = T, sep = "\t")
write(line.2, "GSEA/CD8T/input/CD8_low.down.group.cls", ncolumns = length(line.2), append = T, sep = "\t")
write(line.3, "GSEA/CD8T/input/CD8_low.down.group.cls", ncolumns = length(line.3), append = T, sep = "\t")

# CD8T_high.up
expr <- GetAssayData(tumor_CD8_high, slot = "counts")[rownames(CD8_high.up),]
expr <- data.frame(NAME = rownames(expr), Description = rep('na',nrow(expr)), expr, stringAsFactors = F)
write('#1.2', 'GSEA/CD8T/input/CD8_high.up.expr.gct', ncolumns =1)
write(c(nrow(expr), (ncol(expr) - 2)), 'GSEA/CD8T/input/CD8_high.up.expr.gct', ncolumns = 2, append = T, sep = "\t")
write.table(expr, 'GSEA/CD8T/input/CD8_high.up.expr.gct', append=T, sep="\t", quote = F, row.names = F)
line.1 <- c((ncol(expr) - 2), 2, 1)
tmp <- table(tumor_CD8_high@meta.data$Class)
line.2 <- c('#', names(tmp))
line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
write(line.1, "GSEA/CD8T/input/CD8_high.up.group.cls", ncolumns = length(line.1), append = T, sep = "\t")
write(line.2, "GSEA/CD8T/input/CD8_high.up.group.cls", ncolumns = length(line.2), append = T, sep = "\t")
write(line.3, "GSEA/CD8T/input/CD8_high.up.group.cls", ncolumns = length(line.3), append = T, sep = "\t")

# CD8T_high.down
expr <- GetAssayData(tumor_CD8_high, slot = "counts")[rownames(CD8_high.down),]
expr <- data.frame(NAME = rownames(expr), Description = rep('na',nrow(expr)), expr, stringAsFactors = F)
write('#1.2', 'GSEA/CD8T/input/CD8_high.down.expr.gct', ncolumns =1)
write(c(nrow(expr), (ncol(expr) - 2)), 'GSEA/CD8T/input/CD8_high.down.expr.gct', ncolumns = 2, append = T, sep = "\t")
write.table(expr, 'GSEA/CD8T/input/CD8_high.down.expr.gct', append=T, sep="\t", quote = F, row.names = F)
line.1 <- c((ncol(expr) - 2), 2, 1)
tmp <- table(tumor_CD8_high@meta.data$Class)
line.2 <- c('#', names(tmp))
line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
write(line.1, "GSEA/CD8T/input/CD8_high.down.group.cls", ncolumns = length(line.1), append = T, sep = "\t")
write(line.2, "GSEA/CD8T/input/CD8_high.down.group.cls", ncolumns = length(line.2), append = T, sep = "\t")
write(line.3, "GSEA/CD8T/input/CD8_high.down.group.cls", ncolumns = length(line.3), append = T, sep = "\t")

# CD4T_low.up
expr <- GetAssayData(tumor_CD4_low, slot = "counts")[rownames(CD4_low.up),]
expr <- data.frame(NAME = rownames(expr), Description = rep('na',nrow(expr)), expr, stringAsFactors = F)
write('#1.2', 'GSEA/CD4T/input/CD4_low.up.expr.gct', ncolumns =1)
write(c(nrow(expr), (ncol(expr) - 2)), 'GSEA/CD4T/input/CD4_low.up.expr.gct', ncolumns = 2, append = T, sep = "\t")
write.table(expr, 'GSEA/CD4T/input/CD4_low.up.expr.gct', append=T, sep="\t", quote = F, row.names = F)
line.1 <- c((ncol(expr) - 2), 2, 1)
tmp <- table(tumor_CD4_low@meta.data$Class)
line.2 <- c('#', names(tmp))
line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
write(line.1, "GSEA/CD4T/input/CD4_low.up.group.cls", ncolumns = length(line.1), append = T, sep = "\t")
write(line.2, "GSEA/CD4T/input/CD4_low.up.group.cls", ncolumns = length(line.2), append = T, sep = "\t")
write(line.3, "GSEA/CD4T/input/CD4_low.up.group.cls", ncolumns = length(line.3), append = T, sep = "\t")

# CD4T_low.down
expr <- GetAssayData(tumor_CD4_low, slot = "counts")[rownames(CD4_low.down),]
expr <- data.frame(NAME = rownames(expr), Description = rep('na',nrow(expr)), expr, stringAsFactors = F)
write('#1.2', 'GSEA/CD4T/input/CD4_low.down.expr.gct', ncolumns =1)
write(c(nrow(expr), (ncol(expr) - 2)), 'GSEA/CD4T/input/CD4_low.down.expr.gct', ncolumns = 2, append = T, sep = "\t")
write.table(expr, 'GSEA/CD4T/input/CD4_low.down.expr.gct', append=T, sep="\t", quote = F, row.names = F)
line.1 <- c((ncol(expr) - 2), 2, 1)
tmp <- table(tumor_CD4_low@meta.data$Class)
line.2 <- c('#', names(tmp))
line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
write(line.1, "GSEA/CD4T/input/CD4_low.down.group.cls", ncolumns = length(line.1), append = T, sep = "\t")
write(line.2, "GSEA/CD4T/input/CD4_low.down.group.cls", ncolumns = length(line.2), append = T, sep = "\t")
write(line.3, "GSEA/CD4T/input/CD4_low.down.group.cls", ncolumns = length(line.3), append = T, sep = "\t")

# CD4T_high.up
expr <- GetAssayData(tumor_CD4_high, slot = "counts")[rownames(CD4_high.up),]
expr <- data.frame(NAME = rownames(expr), Description = rep('na',nrow(expr)), expr, stringAsFactors = F)
write('#1.2', 'GSEA/CD4T/input/CD4_high.up.expr.gct', ncolumns =1)
write(c(nrow(expr), (ncol(expr) - 2)), 'GSEA/CD4T/input/CD4_high.up.expr.gct', ncolumns = 2, append = T, sep = "\t")
write.table(expr, 'GSEA/CD4T/input/CD4_high.up.expr.gct', append=T, sep="\t", quote = F, row.names = F)
line.1 <- c((ncol(expr) - 2), 2, 1)
tmp <- table(tumor_CD4_high@meta.data$Class)
line.2 <- c('#', names(tmp))
line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
write(line.1, "GSEA/CD4T/input/CD4_high.up.group.cls", ncolumns = length(line.1), append = T, sep = "\t")
write(line.2, "GSEA/CD4T/input/CD4_high.up.group.cls", ncolumns = length(line.2), append = T, sep = "\t")
write(line.3, "GSEA/CD4T/input/CD4_high.up.group.cls", ncolumns = length(line.3), append = T, sep = "\t")

# CD4T_high.down
expr <- GetAssayData(tumor_CD4_high, slot = "counts")[rownames(CD4_high.down),]
expr <- data.frame(NAME = rownames(expr), Description = rep('na',nrow(expr)), expr, stringAsFactors = F)
write('#1.2', 'GSEA/CD4T/input/CD4_high.down.expr.gct', ncolumns =1)
write(c(nrow(expr), (ncol(expr) - 2)), 'GSEA/CD4T/input/CD4_high.down.expr.gct', ncolumns = 2, append = T, sep = "\t")
write.table(expr, 'GSEA/CD4T/input/CD4_high.down.expr.gct', append=T, sep="\t", quote = F, row.names = F)
line.1 <- c((ncol(expr) - 2), 2, 1)
tmp <- table(tumor_CD4_high@meta.data$Class)
line.2 <- c('#', names(tmp))
line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
write(line.1, "GSEA/CD4T/input/CD4_high.down.group.cls", ncolumns = length(line.1), append = T, sep = "\t")
write(line.2, "GSEA/CD4T/input/CD4_high.down.group.cls", ncolumns = length(line.2), append = T, sep = "\t")
write(line.3, "GSEA/CD4T/input/CD4_high.down.group.cls", ncolumns = length(line.3), append = T, sep = "\t")

# CD3T_low.up
expr <- GetAssayData(tumor_CD3_low, slot = "counts")[rownames(CD3_low.up),]
expr <- data.frame(NAME = rownames(expr), Description = rep('na',nrow(expr)), expr, stringAsFactors = F)
write('#1.2', 'GSEA/CD3T/input/CD3_low.up.expr.gct', ncolumns =1)
write(c(nrow(expr), (ncol(expr) - 2)), 'GSEA/CD3T/input/CD3_low.up.expr.gct', ncolumns = 2, append = T, sep = "\t")
write.table(expr, 'GSEA/CD3T/input/CD3_low.up.expr.gct', append=T, sep="\t", quote = F, row.names = F)
line.1 <- c((ncol(expr) - 2), 2, 1)
tmp <- table(tumor_CD3_low@meta.data$Class)
line.2 <- c('#', names(tmp))
line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
write(line.1, "GSEA/CD3T/input/CD3_low.up.group.cls", ncolumns = length(line.1), append = T, sep = "\t")
write(line.2, "GSEA/CD3T/input/CD3_low.up.group.cls", ncolumns = length(line.2), append = T, sep = "\t")
write(line.3, "GSEA/CD3T/input/CD3_low.up.group.cls", ncolumns = length(line.3), append = T, sep = "\t")

# CD3T_low.down
expr <- GetAssayData(tumor_CD3_low, slot = "counts")[rownames(CD3_low.down),]
expr <- data.frame(NAME = rownames(expr), Description = rep('na',nrow(expr)), expr, stringAsFactors = F)
write('#1.2', 'GSEA/CD3T/input/CD3_low.down.expr.gct', ncolumns =1)
write(c(nrow(expr), (ncol(expr) - 2)), 'GSEA/CD3T/input/CD3_low.down.expr.gct', ncolumns = 2, append = T, sep = "\t")
write.table(expr, 'GSEA/CD3T/input/CD3_low.down.expr.gct', append=T, sep="\t", quote = F, row.names = F)
line.1 <- c((ncol(expr) - 2), 2, 1)
tmp <- table(tumor_CD3_low@meta.data$Class)
line.2 <- c('#', names(tmp))
line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
write(line.1, "GSEA/CD3T/input/CD3_low.down.group.cls", ncolumns = length(line.1), append = T, sep = "\t")
write(line.2, "GSEA/CD3T/input/CD3_low.down.group.cls", ncolumns = length(line.2), append = T, sep = "\t")
write(line.3, "GSEA/CD3T/input/CD3_low.down.group.cls", ncolumns = length(line.3), append = T, sep = "\t")

# CD3T_high.up
expr <- GetAssayData(tumor_CD3_high, slot = "counts")[rownames(CD3_high.up),]
expr <- data.frame(NAME = rownames(expr), Description = rep('na',nrow(expr)), expr, stringAsFactors = F)
write('#1.2', 'GSEA/CD3T/input/CD3_high.up.expr.gct', ncolumns =1)
write(c(nrow(expr), (ncol(expr) - 2)), 'GSEA/CD3T/input/CD3_high.up.expr.gct', ncolumns = 2, append = T, sep = "\t")
write.table(expr, 'GSEA/CD3T/input/CD3_high.up.expr.gct', append=T, sep="\t", quote = F, row.names = F)
line.1 <- c((ncol(expr) - 2), 2, 1)
tmp <- table(tumor_CD3_high@meta.data$Class)
line.2 <- c('#', names(tmp))
line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
write(line.1, "GSEA/CD3T/input/CD3_high.up.group.cls", ncolumns = length(line.1), append = T, sep = "\t")
write(line.2, "GSEA/CD3T/input/CD3_high.up.group.cls", ncolumns = length(line.2), append = T, sep = "\t")
write(line.3, "GSEA/CD3T/input/CD3_high.up.group.cls", ncolumns = length(line.3), append = T, sep = "\t")

# CD3T_high.down
expr <- GetAssayData(tumor_CD3_high, slot = "counts")[rownames(CD3_high.down),]
expr <- data.frame(NAME = rownames(expr), Description = rep('na',nrow(expr)), expr, stringAsFactors = F)
write('#1.2', 'GSEA/CD3T/input/CD3_high.down.expr.gct', ncolumns =1)
write(c(nrow(expr), (ncol(expr) - 2)), 'GSEA/CD3T/input/CD3_high.down.expr.gct', ncolumns = 2, append = T, sep = "\t")
write.table(expr, 'GSEA/CD3T/input/CD3_high.down.expr.gct', append=T, sep="\t", quote = F, row.names = F)
line.1 <- c((ncol(expr) - 2), 2, 1)
tmp <- table(tumor_CD3_high@meta.data$Class)
line.2 <- c('#', names(tmp))
line.3 <- c(rep(names(tmp)[1], tmp[1]), rep(names(tmp)[2], tmp[2]))
write(line.1, "GSEA/CD3T/input/CD3_high.down.group.cls", ncolumns = length(line.1), append = T, sep = "\t")
write(line.2, "GSEA/CD3T/input/CD3_high.down.group.cls", ncolumns = length(line.2), append = T, sep = "\t")
write(line.3, "GSEA/CD3T/input/CD3_high.down.group.cls", ncolumns = length(line.3), append = T, sep = "\t")

GSEA::GSEA("GSEA/CD8T/input/CD8_low.up.expr.gct","GSEA/CD8T/input/CD8_low.up.group.cls",gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "GSEA/CD8T/output/low/up")
GSEA::GSEA("GSEA/CD8T/input/CD8_low.down.expr.gct","GSEA/CD8T/input/CD8_low.down.group.cls",gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "GSEA/CD8T/output/low/down")
GSEA::GSEA("GSEA/CD8T/input/CD8_high.up.expr.gct","GSEA/CD8T/input/CD8_high.up.group.cls",gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "GSEA/CD8T/output/high/up")
GSEA::GSEA("GSEA/CD8T/input/CD8_high.down.expr.gct","GSEA/CD8T/input/CD8_high.down.group.cls",gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "GSEA/CD8T/output/high/down")

GSEA::GSEA("GSEA/CD4T/input/CD4_low.up.expr.gct","GSEA/CD4T/input/CD4_low.up.group.cls",gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "GSEA/CD4T/output/low/up")
GSEA::GSEA("GSEA/CD4T/input/CD4_low.down.expr.gct","GSEA/CD4T/input/CD4_low.down.group.cls",gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "GSEA/CD4T/output/low/down")
GSEA::GSEA("GSEA/CD4T/input/CD4_high.up.expr.gct","GSEA/CD4T/input/CD4_high.up.group.cls",gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "GSEA/CD4T/output/high/up")
GSEA::GSEA("GSEA/CD4T/input/CD4_high.down.expr.gct","GSEA/CD4T/input/CD4_high.down.group.cls",gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "GSEA/CD4T/output/high/down")

GSEA::GSEA("GSEA/CD3T/input/CD3_low.up.expr.gct","GSEA/CD3T/input/CD3_low.up.group.cls",gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "GSEA/CD3T/output/low/up")
GSEA::GSEA("GSEA/CD3T/input/CD3_low.down.expr.gct","GSEA/CD3T/input/CD3_low.down.group.cls",gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "GSEA/CD3T/output/low/down")
GSEA::GSEA("GSEA/CD3T/input/CD3_high.up.expr.gct","GSEA/CD3T/input/CD3_high.up.group.cls",gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "GSEA/CD3T/output/high/up")
GSEA::GSEA("GSEA/CD3T/input/CD3_high.down.expr.gct","GSEA/CD3T/input/CD3_high.down.group.cls",gs.db = "refdb/h.all.v7.4.symbols.gmt",output.directory = "GSEA/CD3T/output/high/down")