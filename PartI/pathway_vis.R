
# pathway visualization

library(stringr)
library(ggpubr)

for(i in 1:8){
  pathway <- readxl::read_excel("result/without_PVTT_MLN/tumor/overlap/CD3T/select_pathway.xlsx",skip = 1,sheet = i)
  p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
    theme(axis.text = element_text(size = 14),
          axis.title.y = element_blank(),
          axis.title = element_text(size = 18),
          panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
    scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
  ggsave(paste0(i,".pdf"),p,width = 8,height = 10)
}

##===== CD3T =========

wkd <- "result/without_PVTT_MLN/tumor/overlap/CD3T/"

### CD3+T UP

# CAF_MLG
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 1)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"UP/","CAF_Malignant.cell/","pathway.pdf"),p,width = 5,height = 6)
# CAF_TEC
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 2)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"UP/","CAF_TEC/","pathway.pdf"),p,width = 7,height = 6)

# TEC_MLG
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 4)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"UP/","TEC_Malignant.cell/","pathway.pdf"),p,width = 5,height = 2)

### CD3+T Down

# CAF_MLG
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 5)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","CAF_Malignant.cell/","pathway.pdf"),p,width = 8,height = 12)
# CAF_TEC
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 6)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","CAF_TEC/","pathway.pdf"),p,width = 8,height = 12)

# CAF_TEC_Malignant.cell
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 7)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","CAF_TEC_Malignant.cell/","pathway.pdf"),p,width = 8,height = 8)

# TEC_MLG
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 8)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","TEC_Malignant.cell/","pathway.pdf"),p,width = 6,height = 3)

##===== CD4T =========

wkd <- "result/without_PVTT_MLN/tumor/overlap/CD4T/"

### CD4+T UP

# CAF_MLG
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 1)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"UP/","CAF_Malignant.cell/","pathway.pdf"),p,width = 6,height = 3)
# CAF_TEC
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 2)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"UP/","CAF_TEC/","pathway.pdf"),p,width = 8,height = 8)

# TEC_MLG
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 4)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"UP/","TEC_Malignant.cell/","pathway.pdf"),p,width = 6,height = 3)

### CD4+T Down

# CAF_MLG
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 5)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","CAF_Malignant.cell/","pathway.pdf"),p,width = 9,height = 12)
# CAF_TEC
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 6)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","CAF_TEC/","pathway.pdf"),p,width = 8,height = 16)

# CAF_TEC_Malignant.cell
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 7)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","CAF_TEC_Malignant.cell/","pathway.pdf"),p,width = 8,height = 12)

# TEC_MLG
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 8)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","TEC_Malignant.cell/","pathway.pdf"),p,width = 8,height = 6)


##===== CD8T =========

wkd <- "result/without_PVTT_MLN/tumor/overlap/CD8T/"

### CD8+T UP

# CAF_TEC
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 2)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"UP/","CAF_TEC/","pathway.pdf"),p,width = 7,height = 4)


### CD4+T Down

# CAF_MLG
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 5)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","CAF_Malignant.cell/","pathway.pdf"),p,width = 8,height = 10)
# CAF_TEC
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 6)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","CAF_TEC/","pathway.pdf"),p,width = 8,height = 16)

# CAF_TEC_Malignant.cell
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 7)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","CAF_TEC_Malignant.cell/","pathway.pdf"),p,width = 8,height = 10)

# TEC_MLG
pathway <- readxl::read_excel(paste0(wkd,"select_pathway.xlsx"),skip = 1,sheet = 8)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","TEC_Malignant.cell/","pathway.pdf"),p,width = 7,height = 8)
##======================================================================================

##================= overlap between CD4T & CD8T ==============================

wkd = "result/without_PVTT_MLN/tumor/overlap/CD4T_vs_CD8T/without_CD45/"

## UP

# CAF
pathway <- readxl::read_excel(paste0(wkd,"UP/","select_pathway.xlsx"),skip = 1,sheet = 1)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"UP/","CAF/","pathway.pdf"),p,width = 9,height = 3)

# TEC
pathway <- readxl::read_excel(paste0(wkd,"UP/","select_pathway.xlsx"),skip = 1,sheet = 2)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"UP/","TEC/","pathway.pdf"),p,width = 9,height = 3)

# MLG
pathway <- readxl::read_excel(paste0(wkd,"UP/","select_pathway.xlsx"),skip = 1,sheet = 3)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"UP/","Malignant.Cell/","pathway.pdf"),p,width = 9,height = 10)

## Down

# CAF
pathway <- readxl::read_excel(paste0(wkd,"Down/","select_pathway.xlsx"),skip = 1,sheet = 1)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","CAF/","pathway.pdf"),p,width = 6,height = 9)

# TEC
pathway <- readxl::read_excel(paste0(wkd,"Down/","select_pathway.xlsx"),skip = 1,sheet = 2)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","TEC/","pathway.pdf"),p,width = 6,height = 9)

# MLG
pathway <- readxl::read_excel(paste0(wkd,"Down/","select_pathway.xlsx"),skip = 1,sheet = 3)
p <- ggbarplot(pathway,x="Description",y="-LogP",fill = rainbow(length(pathway$Description))) + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 18),
        panel.border = element_rect(size = 1.5,colour = "black",fill = NA)) +
  scale_x_discrete(labels = function(x) str_wrap(x,width = 60))
ggsave(paste0(wkd,"Down/","Malignant.Cell/","pathway.pdf"),p,width = 7,height = 12)
##======================================================================================