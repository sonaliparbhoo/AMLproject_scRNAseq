#This script does the following:
#1) Pick the resolution for the clustering
#2) Using scGate (https://github.com/carmonalab/scGate) creates a gating model to subset the different cell subsets
# This can be done also in Seurat but I prefer scGate because you can also exclude negative margers to make the model 
#more precise and subset the object on the output 
#3) According to the expression and the DGE (FindAllMarkers) labels the clusters and then creates a seurat object for 
#each subset

rm(list = ls())

#Load packages
library(Seurat)
library(scGate)
library(ggplot2)
library(magrittr)

#library(ComplexHeatmap)
#library(cowplot)
#library(dplyr)
#library(muscat)
#library(purrr)
#library(RColorBrewer)
#library(viridis)
#library(scran)
#library(SingleCellExperiment)
#library(scater)
#library(tibble)
#library(patchwork)
#library(scGate)

setwd("~/Documents/AML_project/scRNA_AMLproj")
load("scripts/scRNAseq_step2.rds")

# set cluster IDs to resolution 1 clustering
sobj <- SetIdent(sobj, value = "integrated_snn_res.0.8")

DefaultAssay(sobj) <- "RNA" #Put "RNA" as default assay

#UMAP (Visualize clusters)
DimPlot(sobj, reduction = "umap", label = TRUE)

#Plot different phenotypes
#CD3
CD3 <- paste0("CD3D$|CD3G$")
CD3 <- rownames(sobj)[grep(CD3, rownames(sobj))]
CD3 <- gating_model(name = "CD3", signature = CD3)

CD3 <- scGate(sobj, model = CD3)

tiff("./plots/CD3.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
DimPlot(CD3, group.by = "is.pure", cols = c(list(Impure = "blue", Pure = "red"))) + ggtitle("CD3+ T cells")
dev.off()

#CD8
CD8 <- paste0("CD8A$|CD8B$")
CD8 <- rownames(sobj)[grep(CD8, rownames(sobj))]
CD8 <- gating_model(name = "CD8", signature = CD8)
CD8 <- scGate(sobj, model = CD8)

tiff("./plots/CD8.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
DimPlot(CD8, group.by = "is.pure", cols = c(list(Impure = "blue", Pure = "red"))) + ggtitle("CD8+ T cells")
dev.off()

#CD4
CD4 <- paste0("[.]CD4$")
CD4 <- rownames(sobj)[grep(CD4, rownames(sobj))]
CD4 <- gating_model(name = "CD4", signature = CD4)
CD4 <- scGate(sobj, model = CD4)

tiff("./plots/CD4.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
DimPlot(CD4, group.by = "is.pure", cols = c(list(Impure = "blue", Pure = "red"))) + ggtitle("CD4+ T cells")
dev.off()

#NK
NK <- paste0("NCAM1$|KLRD1$|KLRG1$")
NK <- rownames(sobj)[grep(NK, rownames(sobj))]
CD3 <- paste0("CD3D$|CD3G$")
CD3 <- rownames(sobj)[grep(CD3, rownames(sobj))]
NK <- gating_model(level=1, name = "NK", signature = NK)
NK <- gating_model(model=NK, level=1, name = "CD3", signature = CD3, negative=TRUE)
NK <- scGate(sobj, model = NK)

tiff("./plots/NK.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
DimPlot(NK, group.by = "is.pure", cols = c(list(Impure = "blue", Pure = "red"))) + ggtitle("NK cells")
dev.off()

#NKT
NKT <- paste0("CD3D$|CD3G$|NCAM1$")
NKT <- rownames(sobj)[grep(NKT, rownames(sobj))]
CD4 <- paste0("[.]CD4$")
CD4 <- rownames(sobj)[grep(CD4, rownames(sobj))]
CD8 <- paste0("CD8A$|CD8B$")
CD8 <- rownames(sobj)[grep(CD8, rownames(sobj))]
NKT <- gating_model(level=1, name = "NKT", signature = NKT)
NKT <- gating_model(model=NKT, level=1, name = "CD8", signature = CD8, negative=TRUE)
NKT <- gating_model(model=NKT, level=1, name = "CD4", signature = CD4, negative=TRUE)
NKT <- scGate(sobj, model = NKT)

tiff("./plots/GD.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
DimPlot(NKT, group.by = "is.pure", cols = c(list(Impure = "blue", Pure = "red"))) + ggtitle("NKT cells")
dev.off()

#GD
GD <- paste0("CD3D$|TRGC1$|TRDV1$|TRDV2$|TRDV3")
GD <- rownames(sobj)[grep(GD, rownames(sobj))]
CD4 <- paste0("[.]CD4$")
CD4 <- rownames(sobj)[grep(CD4, rownames(sobj))]
CD8 <- paste0("CD8A$|CD8B$")
CD8 <- rownames(sobj)[grep(CD8, rownames(sobj))]
GD <- gating_model(level=1, name="GD", signature=GD)
GD <- gating_model(model=GD, level=1, name = "CD8", signature = CD8, negative=TRUE)
GD <- gating_model(model=GD, level=1, name = "CD4", signature = CD4, negative=TRUE)
GD <- scGate(sobj, model = GD)

tiff("./plots/GD.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
DimPlot(GD, group.by = "is.pure", cols = c(list(Impure = "blue", Pure = "red"))) + ggtitle("GD cells")
dev.off()

#MAIT
MAIT <- paste0("TRAV1-2$|SLC4A10$")
MAIT <- rownames(sobj)[grep(MAIT, rownames(sobj))]
MAIT <- gating_model(name = "MAIT", signature = MAIT)

MAIT <- scGate(sobj, model = MAIT)

tiff("./plots/MAIT.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
DimPlot(MAIT, group.by = "is.pure", cols = c(list(Impure = "blue", Pure = "red"))) + ggtitle("MAIT cells")
dev.off()

tiff("./plots/sobj.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
DimPlot(sobj, label = TRUE)
dev.off()

#Remove non-T, non-NK cells (cluster 13, 15, 16, 17)
plot <- DimPlot(sobj, reduction = "umap")
T_NK <- CellSelector(plot=plot)
sobj <- subset(sobj, cells = T_NK)

#Look at the top 50 markers per cluster
DefaultAssay(sobj) <- "integrated"
mark <- FindAllMarkers(sobj)
top100 <- mark %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
df <- data.frame(top100$cluster, top100$gene)
df <- reshape(transform(df, indx = ave(as.character(top100.cluster), top100.cluster, FUN = seq)), 
              idvar = "indx", timevar = "top100.cluster", direction = "wide")
View(df)
writexl::write_xlsx(df, "df_all.xlsx")

#Label subsets
sobj <- RenameIdents(object = sobj,
                     "0" = "CD4+ T cells",
                     "1" = "CD4+ T cells", 
                     "4" = "CD4+ T cells",
                     "9" = "CD4+ T cells",
                     "10" = "CD4+ T cells",
                     "9" = "CD4+ T cells",
                     "12" = "MAIT", 
                     "11" = "GD T cells",
                     "5" = "NK cells",
                     "6" = "NKT cells",
                     "17" = "NK cells",
                     "14" = "CD8+ T cells",
                     "3" = "CD8+ T cells",
                     "7" = "CD8+ T cells",
                     "8" = "CD8+ T cells",
                     "2" = "CD8+ T cells")


tiff("./plots/sobj_ann.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
DimPlot(sobj, label = TRUE)
dev.off()

#Subsets
#CD8 (take all CD8+ excluding CD4, MAIT)
CD8 <- paste0("CD8A$|CD8B$")
CD8 <- rownames(sobj)[grep(CD8, rownames(sobj))]
CD4 <- paste0("[.]CD4$")
CD4 <- rownames(sobj)[grep(CD4, rownames(sobj))]
mmCD8 <- scGate::gating_model(level=1, name="CD8T", signature = CD8)
mmCD8 <- scGate::gating_model(model=mmCD8, level=1, name="MAIT", signature = c("ENSG00000144290.SLC4A10", "ENSG00000256553.TRAV1-2"), negative=TRUE)
mmCD8 <- scGate::gating_model(model=mmCD8, level=1, name="CD4T", signature = CD4, negative=TRUE)
CD8 <- scGate(sobj, model = mmCD8)
CD8sub <- CD8
CD8sub <- subset(CD8, subset = `is.pure.level1` == "Pure")

dim(CD8sub)

tiff("./plots/CD8sub.tiff", width = 5*300, height = 5*300, res = 300, pointsize = 5)     
DimPlot(CD8sub)
dev.off()

saveRDS(CD8sub, file = "scripts/CD8sub.rds")

#CD4
CD4 <- paste0("[.]CD4$")
CD4 <- rownames(sobj)[grep(CD4, rownames(sobj))]
mmCD4 <- gating_model(name = "CD4", signature = CD4)
CD4 <- scGate(sobj, model = mmCD4)
CD4sub <- CD4
CD4sub <- subset(CD4, subset = `is.pure.level1` == "Pure")
tiff("./plots/CD4sub.tiff", width = 5*300, height = 5*300, res = 300, pointsize = 5)     
DimPlot(CD4sub)
dev.off()

dim(CD4sub)

saveRDS(CD4sub, file = "scripts/CD4sub.rds")

#NK
NK <- paste0("NCAM1$|KLRD1$|KLRG1$")
NK <- rownames(sobj)[grep(NK, rownames(sobj))]
CD3 <- paste0("CD3D$|CD3G$")
CD3 <- rownames(sobj)[grep(CD3, rownames(sobj))]
mmNK <- gating_model(level=1, name = "NK", signature = NK)
mmNK <- gating_model(model=mmNK, level=1, name = "CD3", signature = CD3, negative=TRUE)
NK <- scGate(sobj, model=mmNK)
NKsub <- NK
NKsub <- subset(NK, subset = `is.pure.level1` == "Pure")

tiff("./plots/NKsub.tiff", width = 5*300, height = 5*300, res = 300, pointsize = 5)     
DimPlot(NKsub)
dev.off()

dim(NKsub)

saveRDS(NKsub, file = "scripts/NKsub.rds")

#MAIT
MAIT <- paste0("TRAV1-2$|SLC4A10$")
MAIT <- rownames(sobj)[grep(MAIT, rownames(sobj))]
mmMAIT <- gating_model(name = "MAIT", signature = MAIT)
MAIT <- scGate(sobj, model = mmMAIT)
MAITsub <- MAIT
MAITsub <- subset(MAITsub, subset = `is.pure.level1` == "Pure")
tiff("./plots/NKsub.tiff", width = 5*300, height = 5*300, res = 300, pointsize = 5)     
DimPlot(MAITsub)
dim(MAITsub)

saveRDS(MAITsub, file = "scripts/MAITsub.rds")

#GD
GD <- paste0("CD3D$|TRGC1$|TRDV1$|TRDV2$|TRDV3")
GD <- rownames(sobj)[grep(GD, rownames(sobj))]
CD4 <- paste0("[.]CD4$")
CD4 <- rownames(sobj)[grep(CD4, rownames(sobj))]
CD8 <- paste0("CD8A$|CD8B$")
CD8 <- rownames(sobj)[grep(CD8, rownames(sobj))]
mmGD <- gating_model(level=1, name="GD", signature=GD)
mmGD <- gating_model(model=mmGD, level=1, name = "CD8", signature = CD8, negative=TRUE)
mmGD <- gating_model(model=mmGD, level=1, name = "CD4", signature = CD4, negative=TRUE)
GD <- scGate(sobj, model = mmGD)
GDsub <- GD
GDsub <- subset(GDsub, subset = `is.pure.level1` == "Pure")

tiff("./plots/GDsub.tiff", width = 5*300, height = 5*300, res = 300, pointsize = 5)     
DimPlot(GDsub)
dim(GDsub)

saveRDS(GDsub, file = "scripts/GDsub.rds")

#Pathway analysis (Need to be reviewed!!!)
library(stringr)
library(SCPA)

sobj$cluster_id <- sobj@active.ident

CD8_NR <- seurat_extract(sobj,
                         meta1 = "cluster_id", value_meta1 = "CD8+ T cells",
                         meta2 = "group_id", value_meta2 = "NonRes")

CD8_CR <- seurat_extract(sobj,
                         meta1 = "cluster_id", value_meta1 = "CD8+ T cells",
                         meta2 = "group_id", value_meta2 = "Res")

ss1 <- strsplit(rownames(CD8_NR), ".", fixed=TRUE)
rownames(CD8_NR) <- sapply(ss1, .subset, 2)

ss2 <- strsplit(rownames(CD8_CR), ".", fixed = TRUE)
rownames(CD8_CR) <- sapply(ss2, .subset, 2)

pathways <- msigdbr::msigdbr("Homo sapiens", "H") %>%
  format_pathways()

CD8_CR_NR <- compare_pathways(samples = list(CD8_CR, CD8_NR), 
                              pathways = pathways)

plot_rank(scpa_out = CD8_CR_NR, 
          pathway = "HALLMARK_INFLAMMATORY_RESPONSE", 
          base_point_size = 2, 
          highlight_point_size = 3)

pathways <- "combined_metabolic_pathways.csv"

CD8_CR_NR <- CD8_CR_NR %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'seagreen2',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

aa <- CD8_CR_NR %>% 
  filter(grepl(pattern = "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES", ignore.case = T, x = Pathway))
View(CD8_CR_NR)
ggplot(CD8_CR_NR, aes(FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = CD8_CR_NR$color, stroke = 0.3) +
  geom_point(data = aa, shape = 21, cex = 2.8, fill = "orangered2", color = "black", stroke = 0.3) +
  xlim(-20, 80) +
  ylim(0, 11) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)
