##STEP4: CD8
#This script does the followings:
#Loads in the CD8 subset, find the most variable features, scale the data, run dimensionality reduction, run clustering
#Plots single genes and scores to identify different subsets
#Labels clusters and look at clusters abundancies

###In progress
#Trajectory inference
#single cell pathway analysis
#DGE along and across trajectories (TradeSeq)
#DGE between conditions (condiments)
#velocity inference
#clonotype analysis

rm(list = ls())

library(Seurat)
library(scGate)
library(slingshot)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(magrittr)
library(SCPA)
library(dyno)
library(scRepertoire)
library(SeuratWrappers)
library(monocle3)
library(hdf5r)
library(tibble)
library(msigdbr)
library(miloR)
library(statmod)
library(scater)
library(SCPA)
library(stringr)
library(scuttle)
library(destiny)

setwd("~/Documents/AML_project/scRNA_AMLproj")

CD8 <-readRDS("scripts/CD8sub.rds")

#Subclustering the integrated assay as suggested here https://github.com/satijalab/seurat/issues/2087 by timoast
# Identify the 2000 most variable genes in the RNA assay 
DefaultAssay(CD8) <- "RNA"
CD8 <- FindVariableFeatures(CD8, selection.method = "vst", 
                            nfeatures = 2000, 
                            verbose = FALSE)

# Scale the counts in the integrated assay
DefaultAssay(CD8) <- "integrated"
CD8 <- ScaleData(CD8)

# Run PCA and UMAP
CD8 <- RunPCA(CD8)
CD8 <- RunUMAP(CD8, dims = 1:30, verbose = FALSE)
#CD8 <- RunTSNE(CD8, dims = 1:30, verbose = FALSE)

# Look at the batches
DimPlot(CD8, reduction = "pca", group.by = "batch")
DimPlot(CD8, reduction = "umap", group.by = "batch")

# Plot the elbow plot
ElbowPlot(object = CD8, ndims = 30)

# Determine the K-nearest neighbor graph
CD8 <- FindNeighbors(object = CD8, dims = 1:20)

# Determine the clusters for various resolutions                                
CD8 <- FindClusters(object = CD8, resolution = 1.2)

# set cluster IDs to resolution 1 clustering
CD8 <- SetIdent(CD8, value = "integrated_snn_res.1.2")

DefaultAssay(CD8) <- "RNA" #Put "RNA" as default assay

#Visualize the new subcluster
p <- DimPlot(object = CD8, 
             reduction = "umap", 
             label = TRUE,
             label.size = 3,
             repel = TRUE)
p

###Visualize single markers and scores distributions
#Check they are CD8
T_CD8 <- paste0("CD8A$|CD8B$")
T_CD8 <- rownames(CD8)[grep(T_CD8, rownames(CD8))]
FeaturePlot(CD8, T_CD8) 

#Check there are no CD4
T_CD4 <- paste0("[.]CD4$")
T_CD4 <- rownames(CD8)[grep(T_CD4, rownames(CD8))]
FeaturePlot(CD8, T_CD4) 

#Single Markers
TCF7 <- "TCF7$"
TCF7 <- rownames(CD8)[grep(TCF7, rownames(CD8))]
FeaturePlot(CD8, TCF7,min.cutoff = "q10", max.cutoff = "q90") + ggtitle("TCF7")

CX3CR1 <- "CX3CR1$"
CX3CR1 <- rownames(CD8)[grep(CX3CR1, rownames(CD8))]
FeaturePlot(CD8, CX3CR1, min.cutoff = "q10", max.cutoff = "q90") + ggtitle("CX3CR1")

tiff("./plots/CX3CR1.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, CX3CR1, pt.size = 0.1, order = T, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.8, 1.9), label = c("0", "Maximum")) +
  ggtitle("Effectors CX3CR1+")
dev.off()

LEF1 <- "LEF1$"
LEF1 <- rownames(CD8)[grep(LEF1, rownames(CD8))]
FeaturePlot(CD8, LEF1, min.cutoff = "q10", max.cutoff = "q90") + ggtitle("LEF1")

CD69 <- "CD69$"
CD69 <- rownames(CD8)[grep(CD69, rownames(CD8))]

tiff("./plots/CD69.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, CD69, min.cutoff = "q10", max.cutoff = "q90") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.05, 0.68), label = c("0", "Maximum")) +
  ggtitle("CD69")
dev.off()

PD1 <- "PDCD1$"
PD1 <- rownames(CD8)[grep(PD1, rownames(CD8))]
FeaturePlot(CD8, PD1, min.cutoff = "q10", max.cutoff = "q90") + ggtitle("PD1")

IL7R <- "IL7R$"
IL7R <- rownames(CD8)[grep(IL7R, rownames(CD8))]
tiff("./plots/IL7R.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, IL7R, min.cutoff = "q10", max.cutoff = "q90") + ggtitle("IL7R") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.05, 0.68), label = c("0", "Maximum")) +
  ggtitle("IL7R")
dev.off()

GZMK <- "GZMK$"
GZMK <- rownames(CD8)[grep(GZMK, rownames(CD8))]
tiff("./plots/GZMK.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, GZMK, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.05, 0.68), label = c("0", "Maximum")) +
  ggtitle("GZMK")
dev.off()

EOMES <- "EOMES$"
EOMES <- rownames(CD8)[grep(EOMES, rownames(CD8))]
FeaturePlot(CD8, EOMES, min.cutoff = "q10", max.cutoff = "q90") + ggtitle("EOMES")

TBX21 <- "TBX21$"
TBX21 <- rownames(CD8)[grep(TBX21, rownames(CD8))]
FeaturePlot(CD8, TBX21, min.cutoff = "q10", max.cutoff = "q90") + ggtitle("TBX21")

GZMB <- "GZMB$"
GZMB <- rownames(CD8)[grep(GZMB, rownames(CD8))]
tiff("./plots/GZMB.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, GZMB, min.cutoff = "q10", max.cutoff = "q90") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.05, 0.68), label = c("0", "Maximum")) +
  ggtitle("GZMB")
dev.off()

GZMA <- "GZMA$"
GZMA <- rownames(CD8)[grep(GZMA, rownames(CD8))]
FeaturePlot(CD8, GZMA, min.cutoff = "q10", max.cutoff = "q90") + ggtitle("GZMA")

GNLY <- "GNLY$"
GNLY <- rownames(CD8)[grep(GNLY, rownames(CD8))]
tiff("./plots/GNLY.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, GNLY, min.cutoff = "q10", max.cutoff = "q90") +   
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")), breaks=c(0.05, 0.68), label = c("0", "Maximum")) +
  ggtitle("GNLY")
dev.off()

FCGR3A <- "FCGR3A$"
FCGR3A <- rownames(CD8)[grep(FCGR3A, rownames(CD8))]
FeaturePlot(CD8, FCGR3A, min.cutoff = "q10", max.cutoff = "q90") 

CRTAM <- "CRTAM$"
CRTAM <- rownames(CD8)[grep(CRTAM, rownames(CD8))]
FeaturePlot(CD8, CRTAM, min.cutoff = "q10", max.cutoff = "q90") 

KLF2 <- "KLF2$"
KLF2 <- rownames(CD8)[grep(KLF2, rownames(CD8))]
FeaturePlot(CD8, KLF2, min.cutoff = "q10", max.cutoff = "q90") 

CXCR6 <- "CXCR6$"
CXCR6 <- rownames(CD8)[grep(CXCR6, rownames(CD8))]
FeaturePlot(CD8, CXCR6, min.cutoff = "q10", max.cutoff = "q90") 

#Scores
naive_genes <- "TCF7$|CCR7|SELL|LEF1|GZMK"
naive <- list(rownames(CD8)[grep(naive_genes, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = naive, name = "naive", random.seed = 1)
FeaturePlot(CD8, "naive1", min.cutoff = "q10", max.cutoff = "q90") + ggtitle("Tnaive")

Tcm_genes <- "CCR7|SELL|IL7R"
Tcm <- list(rownames(CD8)[grep(Tcm_genes, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = Tcm, name = "Tcm", random.seed = 1)
FeaturePlot(CD8, "Tcm1", min.cutoff = "q10", max.cutoff = "q90") + ggtitle("Tcm")

Teff_genes <- "CX3CR1|FGFBP2|FCGR3A"
Teff <- list(rownames(CD8)[grep(Teff_genes, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = Teff, name = "Teff", random.seed = 1)
FeaturePlot(CD8, "Teff1") + ggtitle("Teff")

Tpex <- "GZMK|IL7R"
Tpex <- list(rownames(CD8)[grep(Tpex, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = Tpex, name = "Tpex", random.seed = 1)
FeaturePlot(CD8, "Tpex1") + ggtitle("Tpex")

tiff("./plots/stem-like.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, "Tpex1", pt.size = 0.1, order = T, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu"))) + NoLegend() +
  ggtitle("stem-like IL7R+GZMK+")
dev.off()

# signatures
#Scores
#####Change --> all AddModuleScore or all AddModuleScore_UCell
#signatures
sig <- readxl::read_xlsx("signatures/sig.xlsx", sheet = 1)
sig_naive <- as.list(sig$`Naive (from Szabo et al. 2019)`[!is.na(sig$`Naive (from Szabo et al. 2019)`)])
sig_naive <- lapply(sig_naive, function(x){paste0("[.]", x, "$")})
sig_naive <- paste0(sig_naive, collapse = "|")
sig_naive <- list(rownames(CD8)[grep(sig_naive, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = sig_naive, name = "sig_naive", random.seed = 1)

tiff("./plots/naive.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, "sig_naive1", pt.size = 0.1, order = T, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.02, 0.16), label = c("0", "Maximum")) +
  ggtitle("Naive (from Szabo et al. 2019)")
dev.off()

sig_stem <- as.list(sig$`Stemness (from Pace et al. 2018)`[!is.na(sig$`Stemness (from Pace et al. 2018)`)])
sig_stem <- lapply(sig_stem, function(x){paste0("[.]", x, "$")})
sig_stem <- paste0(sig_stem, collapse = "|")
sig_stem <- list(rownames(CD8)[grep(sig_stem, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = sig_stem, name = "sig_stem", random.seed = 1)

tiff("./plots/stem.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, "sig_stem1", pt.size = 0.1, order = T, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.015, 0.125), label = c("0", "Maximum")) +
  ggtitle("Stemness (from Pace et al. 2018)")
dev.off()

Tex <- as.list(sig$`Terminal exhausted (from Miller et al. 2019)`[!is.na(sig$`Terminal exhausted (from Miller et al. 2019)`)])
Tex <- lapply(Tex, function(x){paste0("[.]",x, "$")})
Tex <- paste0(Tex, collapse = "|")
Tex <- list(rownames(CD8)[grep(Tex, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = Tex, name = "Tex", random.seed = 1)


tiff("./plots/Tex.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, "Tex1", pt.size = 0.1, order = T, min.cutoff = "q10", max.cutoff = "q90") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),  
                         breaks=c(0.06, 0.32), label = c("0", "Maximum")) +
  ggtitle("Terminally exhausted (from Miller et al. 2019)")
dev.off()

#Bulk senescence signature
genes <- readxl::read_xlsx("signatures/Model_5_results.xlsx", sheet = 2)
sen <- as.list(genes$gene_name[1:100])
sen <- lapply(sen, function(x){paste0("[.]",x, "$")})
sen <- paste0(sen, collapse = "|")
sen <- list(rownames(CD8)[grep(sen, rownames(CD8))])
CD8 <- AddModuleScore(CD8, features = sen, name = "sen", random.seed = 1)
tiff("./plots/sen_group.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
FeaturePlot(CD8, features = "sen1", pt.size = 0.1, order = T,  min.cutoff = "q10", max.cutoff = "q90") + labs(col = "EXPRESSION") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdBu")),
                         breaks=c(0.09, 0.42), label = c("0", "Maximum")) + 
  ggtitle("")
dev.off()

tiff("./plots/UMAP.tiff", width = 8*250, height = 5*300, res = 300, pointsize = 5)     
DimPlot(CD8, label = TRUE)
dev.off()

#Remove patient dependent subset
CD8 <- subset(CD8, idents = '10', invert = TRUE)

#Clean workspace
rm(list=setdiff(ls(), "CD8"))

CD8 <- RenameIdents(CD8,
                    "0" = "Naive",
                    "1" = "Naive",
                    "2" = "Int",
                    "3" = "Stem-like",
                    "4" = "Senescent-like",
                    "5" = "Tex", 
                    "6" = "Stem-like", 
                    "7" = "Naive", 
                    "8" = "Int", 
                    "9" = "Int-GZMK+",
                    "11" = "Senescent-like", 
                    "12" = "Senescent-like",
                    "13" = "Naive",
                    "14" = "Senescent-like",
                    "15" = "Stem-like",
                    "16" = "Stem-like",
                    "17" = "Stem-like",
                    "18" = "Stem-like",
                    "19" = "Tex",
                    "20" = "Int",
                    "21" = "Int")

tiff("./plots/UMAP_ann.tiff", width = 5*300, height = 5*200, res = 300, pointsize = 5)     
p <- DimPlot(CD8, label = TRUE)
p
dev.off()

tiff("./plots/UMAP_groupId.tiff", width = 5*600, height = 5*250, res = 300, pointsize = 5)     
p <- DimPlot(CD8, label = TRUE, split.by = "group_id")
p
dev.off()

tiff("./plots/UMAP_density.tiff", width = 5*600, height = 5*400, res = 300, pointsize = 5)     
a <- DimPlot(CD8, reduction = 'umap', split.by = "RespTmp")  + NoLegend() + NoAxes() 
a +  geom_density_2d_filled(a$data, mapping = aes(x = a$data[,"UMAP_1"], y = a$data[,"UMAP_2"]), contour_var = "ndensity") + 
  facet_wrap(vars(RespTmp))
dev.off()

mark <- FindAllMarkers(CD8)
mark %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
df <- data.frame(top20$cluster, top20$gene)
df <- reshape(transform(df, indx = ave(as.character(top20.cluster), top20.cluster, FUN = seq)), 
              idvar = "indx", timevar = "top20.cluster", direction = "wide")
View(top20)
writexl::write_xlsx(df, "./genes.xlsx")

tiff("./plots/DoHeat.tiff", width = 5*400, height = 5*500, res = 150, pointsize = 5)     
DoHeatmap(object = CD8, features = top20$gene, size = 3, label = FALSE) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 6)) + ggtitle("Top 20")
dev.off()

#Relative abundance of the each clusters per condition (How to check for "significance"? What test should we apply?)
freq_table <- table(CD8@active.ident, CD8$group_id)
fqs <- prop.table(table(CD8@active.ident, CD8$sample_id), 2)
df <- set_colnames(reshape2::melt(fqs), c("cluster_id", "sample_id", "frequency"))
df$group_id <- CD8$group_id[match(df$sample_id, CD8$sample_id)]

tiff("./plots/box.tiff", width = 5*400, height = 5*60, res = 150, pointsize = 5)     
ggplot(df, aes(x = group_id, y = frequency, color = group_id)) +
  labs(x = NULL, y = "Proportion [%]") +
  theme_bw() + theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.key.height  =  unit(0.8, "lines")) +
  geom_boxplot(aes_string(color = "group_id", fill = "group_id"), position = position_dodge(), alpha = 0.2, 
               outlier.color = NA, show.legend = FALSE) + 
  geom_point(aes_string(x = "group_id", col = "group_id"), position = position_jitter(width = 0.2)) +
  facet_wrap(~ cluster_id, scales = "free_y", ncol = 6) +
  theme_classic()
dev.off()

#DE between int and int-GZMK using DESeq2
df1 <- FindMarkers(CD8, ident.1 = "Int", ident.2 = "Int-GZMK+", test.use = "DESeq2")
df1 %>%
  top_n(n = 100, wt = avg_log2FC) -> top100

DefaultAssay(CD8) <- "integrated"

#Plot heatmap with top100 genes
DoHeatmap(CD8, features = top100$gene, label = TRUE)

#DE between Tex and senescent-like using DESeq2
df2 <- FindMarkers(CD8, ident.1 = "Senescent-like", ident.2 = "Tex", test.use = "DESeq2")
df2 %>%
  top_n(n = 100, wt = avg_log2FC) -> top100

#Plot heatmap with top100 genes
DoHeatmap(CD8, features = top100$gene, label = TRUE)

#Single cell pathway analysis
#See if there are metabolic changes at the trajectory split
#(https://jackbibby1.github.io/SCPA/articles/quick_start.html)
CD8$cluster_id <- CD8@active.ident
sen <- seurat_extract(CD8, meta1 = "cluster_id", value_meta1 = "Senescent-like",)
tex <- seurat_extract(CD8, meta1 = "cluster_id", value_meta1 = "Tex")

ss1 <- strsplit(rownames(sen), ".", fixed=TRUE)
rownames(sen) <- sapply(ss1, .subset, 2)

ss2 <- strsplit(rownames(tex), ".", fixed = TRUE)
rownames(tex) <- sapply(ss2, .subset, 2)

pathways <- msigdbr::msigdbr("Homo sapiens", "H") %>%
  format_pathways()

pathways <- "/home/francesco/Documents/AML_project/scRNA_AMLproj/combined_metabolic_pathways.csv"

sen_tex <- compare_pathways(samples = list(sen, tex), 
                            pathways = pathways)

plot_rank(scpa_out = int_intGZMK, 
          pathway = "HALLMARK_INFLAMMATORY_RESPONSE", 
          base_point_size = 2, 
          highlight_point_size = 3)

sen_tex <- sen_tex %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

aa_path <- sen_tex %>% 
  filter(grepl(pattern = "reactome_respiratory", ignore.case = T, x = Pathway))

ggplot(aa_path, aes(-FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = sen_tex$color, stroke = 0.3) +
  geom_point(data = aa_path, shape = 21, cex = 2.8, fill = "orangered2", color = "black", stroke = 0.3) +
  xlim(-20, 80) +
  ylim(0, 11) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)

CD8$cluster_id <- CD8@active.ident
cell_types <- unique(CD8$cluster_id)
CD8_cond <- SplitObject(CD8, split.by = "group_id")
pathways <- "h_k_r_go_pid_reg_wik.csv"

scpa_out <- list()
for (i in cell_types) {
  HD <- seurat_extract(CD8_cond$HD, 
                       meta1 = "cluster_id", value_meta1 = i)
  NonRes <- seurat_extract(CD8_cond$NonRes, 
                           meta1 = "cluster_id", value_meta1 = i)
  Res <-  seurat_extract(CD8_cond$Res, 
                         meta1 = "cluster_id", value_meta1 = i)
}

ls.cond <- list(HD, NonRes, Res)
ls.cond <- lapply(ls.cond, function(x) {
  ss <- strsplit(rownames(x), ".", fixed = TRUE)
  rownames(x) <- sapply(ss, .subset, 2)
  return(x)})

for(i in cell_types) {
  scpa_out[[i]] <- compare_pathways(ls.cond, pathways) %>%
    select(Pathway, qval) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_")))
}

# DiffusionMap
logcounts <- GetAssayData(CD8, "data")

# transpose matrix (genes as columns, cells as rows)
input_matrix <- t(logcounts[VariableFeatures(CD8), ])
rm(list=setdiff(ls(), "CD8")
   # generate a diffusion map
   dm <- DiffusionMap(as.matrix(input_matrix))
   
   # store the diffusion map as a custom dimensional reduction in the Seurat object
   CD8[["DM"]] <- CreateDimReducObject(embeddings = dm@eigenvectors, key = "DM_", assay = DefaultAssay(CD8))
   
   # plot the diffusion map
   DimPlot(CD8, reduction = "DM") + ggtitle("Diffusion Map")
   DimPlot(CD8, reduction = "DM", split.by = "group_id") + ggtitle("Diffusion Map")
   
   # calculate the diffusion pseudotime (DPT)
   dpt <- DPT(dm)
   
   # color single cells on diffusion map plot according to DPT
   p1 <- plot(dpt, 1:2) + ggtitle("Diffusion Pseudotime (DPT)")
   p1
   
   # create data.frame for easy plotting of DPT
   #tmp <- data.frame(DC1 = dm$DC1,
   #DC2 = dm$DC2,
   #timepoint = cell_type,
   #dpt = dpt$DPT1)
   
   # color single cells on diffusion map plot according to timepoint
   #p2 <- ggplot(tmp, aes(x = DC1, y = DC2, colour = timepoint)) +
   #geom_point() +
   #theme_classic() +
   #ggtitle("Timepoint")
   
   #plot_grid(p1, p2, align = "hv", ncol = 2, rel_widths = c(1, 1))
   
   # plot timepoint vs. DPT
   #ggplot(tmp, aes(timepoint, dpt, colour = timepoint)) +
   #geom_point() + 
   #geom_jitter() +
   #ggtitle ("Timepoint vs. DPT")
#TRAJECTORY INFERENCE using Slingshot
sds <- slingshot(Embeddings(CD8, "umap"), clusterLabels = CD8@active.ident, 
                allow.breaks = TRUE, stretch = 2, reducedDim = "umap", start.clus = "Naive") #Calcualting the trajectory
sds <- SlingshotDataSet(sds)
   
#Change palette
cell_pal <- function(cell_vars, pal_fun,...) {
 if (is.numeric(cell_vars)) {
   pal <- pal_fun(100, ...)
   return(pal[cut(cell_vars, breaks = 100)])
 } else {
   categories <- sort(unique(cell_vars))
   pal <- setNames(pal_fun(length(categories), ...), categories)
   return(pal[cell_vars])
 }
}
   
#Create dataframe pseudotimes
pt <- as.data.frame(slingPseudotime(sds))
colnames(pt) <- c("pseudoT1", "pseudoT2")
names <- rownames(pt)

#Add pseudotime values to Seurat 
CD8$pt1 <- pt$pseudoT1
CD8$pt2 <- pt$pseudoT2
 
#Extract the curves coordinates
curve1 <- slingCurves(sds)[[1]]
curve2 <- slingCurves(sds)[[2]]
 
dev.off()
pdf("./DataAnalysis/CellType/CD8.Tcells/slingshot/Trajectory2.pdf", height=4, width=4)
plot(reducedDim(sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(sds, lwd = 2, col = 'black')
dev.off()
   
plot(reducedDim(sds), col = cell_colors_clust, pch = 16, cex = 0.25)
lines(sds, lwd = 2, type = 'lineages', col = 'black')
#UMAP coordinates for the 1st lineage
p1 <- FeaturePlot(CD8, c("pt1")) + ggtitle("Lineage 1")
UMAP_1 <- p1$data$UMAP_1
UMAP_2 <- p1$data$UMAP_2
   
#Plot 1st lineage
tiff("./plots/traj1.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
p1 <- p1 + geom_path(aes(x = UMAP_1, y = UMAP_2), data = curve1$s[curve1$ord, ] %>% as.data.frame(),
                    col = "black", size = 1, arrow = arrow(), lineend = "round") + scale_color_viridis_c() +  labs(color = "Pseudotime")
p1
dev.off()
   
#UMAP coordinates for the second lineage
p2 <- FeaturePlot(CD8, c("pt2")) + ggtitle("Lineage 2") 
UMAP_1 <- p2$data$UMAP_1
UMAP_2 <- p2$data$UMAP_2
   
tiff("./plots/traj2.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
p2 <- p2 +   geom_path(aes(x = UMAP_1, y = UMAP_2), data = curve2$s[curve2$ord, ] %>% as.data.frame(),
                      col = "black", size = 1, arrow = arrow(), lineend = "round") + scale_color_viridis_c() +  labs(color = "Pseudotime")
p
dev.off()
   
tiff("./plots/traj.tiff", width = 8*130, height = 5*300, res = 150, pointsize = 5)  
plot_grid(p1, p2, ncol = 1, align = "h")
dev.off()
   
   
## Identifying differentially expressed genes along a trajectory
# select the ptime values 
sceCD8$pt1 <- pt$pseudoT1
sceCD8$pt2 <- pt$pseudoT2
 
# get cells in that lineage
lineage_cells <- colnames(sceCD8)[!is.na(sceCD8$pt1)]
 
# remove values for cells not in the lineage
pt1 <- sceCD8$pt1[!is.na(sceCD8$pt1)]
   
# just test variable genes to save some time
genes_to_test <- VariableFeatures(CD8)[1:1000]
   
# get log normalized data to test
cnts <- logcounts(sceCD8)[genes_to_test, lineage_cells]
   
# fit a GAM with a loess term for pseudotime
gam.pval <- apply(cnts, 1, function(z){
 d <- data.frame(z = z, ptime = pt1)
 tmp <- suppressWarnings(gam(z ~ lo(pt1), data=d))
 p <- summary(tmp)[4][[1]][1, 5]
 p
})
   
# adjust pvalues 
res <- tibble(
 id = names(gam.pval),
 pvals = gam.pval,
 qval = p.adjust(gam.pval, method = "fdr")) %>% 
 arrange(qval)
   
head(res)
   
# get log normalized counts 
to_plot <- as.matrix(logcounts(sceCD8)[res$id[1:100], lineage_cells])
   
# arrange cells by pseudotime
ptime_order <- colnames(to_plot)[order(pt1)]
   
# add useful annotations
annotations <- colData(sceCD8)[lineage_cells, 
                              c("pt1", 
                                "ident", 
                                "group_id")] %>% 
 as.data.frame()
ha <- HeatmapAnnotation(df = annotations)
   
Heatmap(to_plot,
       column_order = ptime_order,
       show_column_names = FALSE,
       show_row_names = FALSE,
       top_annotation = ha)

#Fit negative binomial model
set.seed(5)
icMat <- evaluateK(counts = counts(sceCD8), sds = sds, k = 3:10, 
                  nGenes = 200, verbose = T)
 
set.seed(7)
counts <- CD8[["RNA"]]@counts
###Try to run on all the genes on the server 
sce <- fitGAM(counts = counts, sds = sds, nknots = 6, genes = 1:2000, verbose = TRUE, sce=TRUE)
   
# plot our Slingshot lineage trajectories, this time illustrating the new tradeSeq knots
tiff("./plots/traj.tiff", width = 5*500, height = 5*300, res = 300, pointsize = 5)     
plotGeneCount(curve = sds, counts = counts,
             clusters = CD8@active.ident,
             models = sce)
dev.off()
   
#Association test
assoRes <- associationTest(sce)
head(assoRes)
   
#Progenitor markers
startRes <- startVsEndTest(sce)
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[10]]
plotSmoothers(sce, counts, gene = sigGeneStart)
   
### Discovering differentiated cell type markers
# discover marker genes for the differentiated cell types
endRes <- diffEndTest(sce, pairwise = TRUE)
head(endRes)
   
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[1]]
plotSmoothers(sceCD8, counts, sigGene)
   
plotGeneCount(sds, counts, gene = sigGene)
   
earlyDERes <- earlyDETest(sceCD8, knots = c(4, 6))
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
head(rownames(earlyDERes)[oEarly])
plotSmoothers(sce, counts, gene = rownames(earlyDERes)[oEarly][1])
plotGeneCount(sds, counts, gene = rownames(earlyDERes)[oEarly][1])
   
#Finish analysis whole dataset
   
#Trajectory by condition
   
   
library(devtools)
install_github("epurdom/clusterExperiment")
   
nPointsClus <- 20
clusPat <- clusterExpressionPatterns(sce, nPoints = nPointsClus,
                                    genes = rownames(counts)[1:100])
   
df <- bind_cols(
 as.data.frame(reducedDims(sce)$UMAP),
 as.data.frame(colData(sce))
) %>%
 sample_frac(1)
   
   
scores <- condiments::imbalance_score(
 Object = df %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
 conditions = df$group_id,
 k = 20, smooth = 40)
   
df$scores <- scores$scaled_scores
p3 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scores)) +
 geom_point(size = .7) +
 scale_color_viridis_c(option = "C") +
 labs(col = "Scores")
p3
   
#Analyze two condition only
sceCD8 <- as.SingleCellExperiment(CD8, assay = "RNA")
sceCD8 <- slingshot(sceCD8, reducedDim = 'UMAP', clusterLabels = colData(sceCD8)$ident,
                   start.clus = 'Tnaive', approx_points = 150)
df <- bind_cols(
 as.data.frame(reducedDims(sce)$UMAP),
 as.data.frame(colData(sce))
) %>%
 sample_frac(1)
   
#topologyTest(SlingshotDataSet(sceCD8), condition = sceCD8$group_id, rep = 100, methods = "Classifier", threshs = .01)
slingPseudotime_1 <- sce$slingPseudotime_1
slingPseudotime_2 <- sce$slingPseudotime_2
df <- cbind(df, slingPseudotime_1, slingPseudotime_2)
 
curve <- slingCurves(sceCD8)[[1]]
   
p4 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = slingPseudotime_1)) +
 geom_point(size = .7) +
 scale_color_viridis_c() +
 labs(col = "Pseudotime") +
 geom_path(data = curve$s[curve$ord, ] %>% as.data.frame(),
           col = "black", size = 1.5)
p4
   
#Try Clustering using RSEC, clusterExperiment from TradeSeq workflow
#Condiments
df <- bind_cols(
 as.data.frame(reducedDims(sceCD8)$UMAP),
 as.data.frame(colData(sceCD8))
) %>%
 sample_frac(1)
   
scores <- condiments::imbalance_score(
 Object = df %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
 conditions = df$group_id,
 k = 20, smooth = 40)
df$scores <- scores$scaled_scores
p3 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scores)) +
 geom_point(size = .7) +
 scale_color_viridis_c(option = "C") +
 labs(col = "Scores")
p3
   
sceCD8 <- slingshot(sceCD8, reducedDim = 'UMAP',
                   clusterLabels = sceCD8$cluster_id,
                   start.clus = 'Naive', approx_points = 100)
 
set.seed(821)
#topologyTest(SlingshotDataSet(sceCD8), sceCD8$group_id, rep = 100,
#methods = "KS_mean", threshs = .01)
   
curve <- slingCurves(sceCD8)[[1]]
p4 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = pt1)) +
 geom_point(size = .7) +
 scale_color_viridis_c() +
 labs(col = "Pseudotime") +
 geom_path(data = curve$s[curve$ord, ] %>% as.data.frame(),
           col = "black", size = 1.5)
p4
   
p5 <- ggplot(df, aes(x = pt1)) +
 geom_density(alpha = .8, aes(fill = group_id), col = "transparent") +
 geom_density(aes(col = group_id), fill = "transparent",
              guide = FALSE, size = 1.5) +
 labs(x = "Pseudotime", fill = "group_id") +
 guides(col = FALSE, fill = guide_legend(
   override.aes = list(size = 1.5, col = c("#7FC97F", "#BEAED4", "#AA4371"))
 )) +
 scale_fill_brewer(palette = "Accent") +
 scale_color_brewer(palette = "Accent")
   
p5
   
p6 <- ggplot(df, aes(x = pt2)) +
 geom_density(alpha = .8, aes(fill = group_id), col = "transparent") +
 geom_density(aes(col = group_id), fill = "transparent",
              guide = FALSE, size = 1.5) +
 labs(x = "Pseudotime", fill = "group_id") +
 guides(col = FALSE, fill = guide_legend(
   override.aes = list(size = 1.5, col = c("#7FC97F", "#BEAED4", "#AA4371"))
 )) +
 scale_fill_brewer(palette = "Accent") +
 scale_color_brewer(palette = "Accent")
   
p6

progressionTest(SlingshotDataSet(sceCD8), conditions = sceCD8$group_id)
sceCD8 <- fitGAM(counts = as.matrix(assays(sceCD8)$counts),
                sds = sds, genes = 1:1000,
                conditions = factor(colData(sceCD8)$group_id),
                nknots = 6, verbose = TRUE)
   
sceCD8 <- readRDS("scripts/sceCD8_AM.rds")
mean(rowData(sceCD8)$tradeSeq$converged)
   
rowData(sceCD8)$assocRes <- associationTest(sceCD8, lineages = TRUE, l2fc = log2(2))
assocRes <- rowData(sceCD8)$assocRes

NRGenes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionNonRes, "fdr") <= 0.05)
]
RGenes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionRes, "fdr") <= 0.05)
]

HDGenes <-  rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_lineage1_conditionHD, "fdr") <= 0.05)
]

length(RGenes)
length(NRGenes)
length(HDGenes)
UpSetR::upset(fromList(list(Res = RGenes, NonRes = NRGenes, HD = HDGenes)))

### based on mean smoother
yhatSmooth <- predictSmooth(sceCD8, gene = RGenes, nPoints = 50, tidy = FALSE)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                       cluster_cols = FALSE,
                       show_rownames = FALSE, 
                       show_colnames = FALSE)

## the hierarchical trees constructed here, can also be used for 
## clustering of the genes according to their average expression pattern.
cl <- sort(cutree(heatSmooth$tree_row, k = 6))
table(cl)

######Retrieve sceCD8 with UMAP dimred

####To be improved!!


## C7 category is according to gene ontology grouping: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707969/pdf/nihms-743907.pdf
geneSets <- msigdbr(species = "Homo sapiens", category = "C7")
### filter background to only include genes that we assessed.
geneSets$gene_symbol <- toupper(geneSets$gene_symbol)

genenames <- gsub(".*[.]", "", names(sceCD8))
geneSets <- geneSets[geneSets$gene_symbol %in% genenames,]
m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
stats <- assocRes$waldStat_lineage1_conditionRes
names(stats) <- gsub(".*[.]", "", rownames(assocRes))
eaRes <- fgsea(pathways = m_list, stats = stats, nperm = 5e4, minSize = 10)

ooEA <- order(eaRes$pval, decreasing = FALSE)
kable(head(eaRes[ooEA, 1:3], n = 20))

#Differential expression between conditions
condRes <- conditionTest(sceCD8, l2fc = log2(2))
condRes$padj <- p.adjust(condRes$pvalue, "fdr")
mean(condRes$padj <= 0.05, na.rm = TRUE)
sum(condRes$padj <= 0.05, na.rm = TRUE)

conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
conditionGenes <- conditionGenes[!is.na(conditionGenes)]

scales <- brewer.pal(3, "Accent")[1:6]

# plot genes
oo <- order(condRes$waldStat, decreasing = TRUE)

# most significant gene
p6 <- plotSmoothers(sceCD8, assays(sceCD8)$counts,
                    gene = rownames(assays(sceCD8)$counts)[oo[1]],
                    alpha = 1, border = TRUE, curvesCols = scales) +
  scale_color_manual(values = scales) +
  ggtitle(rownames(assays(sceCD8)$counts)[oo[3]])
p6
### based on mean smoother
yhatSmooth <- predictSmooth(sceCD8, gene = conditionGenes, nPoints = 50, tidy = FALSE) %>%
  log1p()

yhatSmoothScaled <- t(apply(yhatSmooth,1, scales::rescale))
heatSmooth <- pheatmap(yhatSmoothScaled[, 51:100],
                       cluster_cols = FALSE,
                       show_rownames = FALSE, show_colnames = FALSE, legend = FALSE,
                       silent = TRUE
)
matchingHeatmap <- pheatmap(yhatSmoothScaled[heatSmooth$tree_row$order, 1:50],
                            cluster_cols = FALSE, cluster_rows = FALSE,
                            show_rownames = FALSE, show_colnames = FALSE, main = "",
                            legend = FALSE, silent = TRUE 
)

p9 <- plot_grid(heatSmooth_TGF[[4]], matchingHeatmap_mock[[4]], ncol = 2)
p9

#Velocity using velocyto.R
#Load in the files
ld.list <- list.files("velocity", pattern = "loom", all.files = TRUE)
ld <- lapply(ld.list, function(x){
  ReadVelocity(paste0("velocity/", x))
})

ld.list <- gsub(".loom", "_", ld.list)

#Rename colnames for each object according to the colnames of the integrated Seurat object
so <- lapply(ld, as.Seurat)
so <- lapply(so, function(x) RenameCells(x, new.names = gsub("x", "-1", colnames(x)), for.merge = T))
so <- lapply(so, function(x) RenameCells(x, new.names = gsub(".*:", "", colnames(x)), for.merge = T))
so <-  lapply(seq_along(ld.list), function(y) RenameCells(so[[y]], new.names = paste0(ld.list[[y]], colnames(so[[y]]))))

#Merge the objects
som <- merge(so[[1]], so[-(1)], merge.data = TRUE)

# Extract only the intersection
som <- som[, intersect(colnames(CD8), colnames(som[, colnames(CD8)]))] # This works

#Create the specific assays
spliced <- CreateAssayObject(GetAssayData(som, assay = "spliced"))
unspliced <- CreateAssayObject(GetAssayData(som, assay = "unspliced"))
ambiguous <- CreateAssayObject(GetAssayData(som, assay = "ambiguous"))

#Import the assays in the original Seurat object
CD8[["spliced"]] <- spliced
CD8[["unspliced"]] <- unspliced
CD8[["ambiguous"]] <- ambiguous

CD8 <- RunVelocity(CD8, deltaT = 1, kCells = 25, fit.quantile = 0.02)
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = CD8)))
names(x = ident.colors) <- levels(x = CD8)
cell.colors <- ident.colors[Idents(object = CD8)]
names(x = cell.colors) <- colnames(x = CD8)
par(mar=c(1,1,1,1))
show.velocity.on.embedding.cor(emb = Embeddings(object = CD8, reduction = "umap"), vel = Tool(object = CD8, 
                                                                                              slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)

#Running velocity using scvelo
#bash
#pip install -U scvelo

# bash
#conda activate r-velocity
#python
#>>> import scvelo as scv

use_condaenv("r-velocity", required = TRUE)
scv <- import("scvelo")

#Assign the spliced assay to the RNA assay and put it as default
CD8[['RNA']] <- CD8[["spliced"]]
DefaultAssay(CD8) <- "RNA"

#save seurat and convert into h5ad
SaveH5Seurat(CD8, filename = "CD8scv.h5Seurat")
Convert("CD8scv.h5Seurat", dest = "h5ad")

#read the h5ad file
adata <- scv$read("CD8scv.h5ad")
adata #check the adata file

## get embedding
emb <- data.frame(adata$obsm['X_umap'])
clusters <- CD8@active.ident
rownames(emb) <- names(clusters)

## get jackbibby1.github.io/SCPA/clusters, convert to colors
col <- rainbow(length(levels(clusters)), s=0.8, v=0.8)
cell.cols <- col[clusters]
names(cell.cols) <- names(clusters)

## simple plot
par(mar=c(1,1,1,1))
plot(emb, col=cell.cols, pch=16,
     xlab='UMAP X', ylab='UMAP Y')
legend(x=4, y=-3.4, 
       legend=levels(clusters),
       col=col, 
       pch=16)

## run scvelo dynamic model
scv$pp$filter_and_normalize(adata, min_shared_counts=20)
scv$pp$moments(adata, method='hnsw') ## normalize and compute moments
scv$tl$recover_dynamics(adata) ## model
scv$tl$velocity(adata, mode = 'dynamical')
scv$tl$velocity_graph(adata)

#adata$write('data/scvelo.h5ad', compression='gzip')
#adata = scv$read('data/scvelo.h5ad')

scv$pl$velocity_embedding_stream(adata, basis="umap", color = "seurat_clusters")

#calculate latent time
scv$tl$latent_time(adata)
scv$pl$scatter(adata, color='latent_time', color_map='gnuplot', size=80)


#Dyno trajectories metabolism
df <- as.matrix(CD8[["RNA"]]@data)
var_genes <- names(sort(apply(df, 1, var), decreasing = TRUE))[1:1000]

counts <- Matrix::t(as(as.matrix(CD8@assays$RNA@counts[var_genes,]), 'sparseMatrix'))
expression <- Matrix::t(as(as.matrix(CD8@assays$RNA@data[var_genes,]), 'sparseMatrix'))

ds <- wrap_expression(expression = expression,
                      counts = counts)
model <- infer_trajectory(ds, dimred = "umap", method = dynmethods::ti_slingshot(),  verbose = T)

model <- model %>% 
  add_dimred(dimred = as.matrix(CD8@reductions$umap@cell.embeddings),
             expression_source = ds$expression)

plot_dimred(model, expression_source = ds$expression, )


plot_dimred(model,
            "pseudotime", 
            expression_source = ds$expression,
            pseudotime = calculate_pseudotime(model), 
            hex_cells = F,
            plot_trajectory = T, 
            size_cells = 1, alpha_cells = 0.8) + 
  theme(aspect.ratio = 1)

plot_dimred(model, 
            expression_source = ds$expression,
            grouping = group_onto_nearest_milestones(model), 
            hex_cells = F,
            plot_trajectory = T, 
            size_cells = 1, alpha_cells = 0.8) + 
  theme(aspect.ratio = 1)

mile_group <- data.frame(group_onto_nearest_milestones(model)) %>%
  set_colnames("milestone") %>%
  rownames_to_column("cell")
CD8$milestone <- mile_group$milestone
CD8_pseudo <- list()
for (i in 1:max(mile_group$milestone)) {
  CD8_pseudo[[i]] <- seurat_extract(CD8, meta1 = "milestone", value_meta1 = i)
}

pathways <- "combined_metabolic_pathways.csv"

CD8_meta <- compare_pathways(samples = CD8_pseudo, 
                             pathways = pathways)

#Clonotype analysis
#Delete HD (we have VDJ only for responders and nonresponders)
CD8 <- subset(CD8, group_id == c("Res", "NonRes"))

CD8@meta.data$cloneType <- factor(CD8@meta.data$cloneType, levels = c("Hyperexpanded (100 < X <= 500)", "Large (20 < X <= 100)", "Medium (5 < X <= 20)", "Small (1 < X <= 5)", "Single (0 < X <= 1)", NA))
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
DimPlot(CD8, group.by = "cloneType", split.by = "group_id") + scale_color_manual(values = c(colorblind_vector(5)), na.value="grey")
TCD8 <- expression2List(CD8, group = "cluster")
Samples <- expression2List(CD8, group = "RespTmp")

p229 <- subset(CD8, orig.ident == "p229")
p229 <- expression2List(p229, group = "cluster")

p219 <- subset(CD8, orig.ident == "p219")
p219 <- expression2List(p219, group = "cluster")

p264 <- subset(CD8, orig.ident == "p264")
p264 <- expression2List(p264, group = "cluster")

pV03 <- subset(CD8, orig.ident == "pV03")
pV03 <- expression2List(pV03, group = "cluster")

pV09 <- subset(CD8, orig.ident == "pV09")
pV09 <- expression2List(pV09, group = "cluster")

clonalProportion(p229, cloneCall = "gene+nt")
clonalProportion(p219, cloneCall = "gene+nt")
clonalProportion(p264, cloneCall = "gene+nt")
clonalProportion(pV03, cloneCall = "gene+nt")
clonalProportion(pV09, cloneCall = "gene+nt")

tiff("./pipelines/plots/clono.tiff", width = 5*10, height = 5*10, res = 300, pointsize = 5)     
compareClonotypes(p229, numbers = 15, 
                  cloneCall="gene+nt", graph = "alluvial")
compareClonotypes(p219, numbers = 15, 
                  cloneCall="gene+nt", graph = "alluvial")
compareClonotypes(p264, numbers = 15, 
                  cloneCall="gene+nt", graph = "alluvial")
compareClonotypes(pV03, numbers = 15, 
                  cloneCall="gene+nt", graph = "alluvial")
compareClonotypes(pV09, numbers = 15, 
                  cloneCall="gene+nt", graph = "alluvial")

quantContig(TCD8, cloneCall="gene+nt", scale = T) + guides(fill = F)
quantContig(Samples, cloneCall="gene+nt", scale = TRUE) + guides(fill = F)

clonalHomeostasis(TCD8, cloneCall = "gene+nt")
clonalProportion(TCD8, cloneCall = "gene+nt")

clonalOverlap(TCD8, cloneCall = "gene+nt", method = "morisita") + 
  scale_fill_gradientn(colors = rev(colorblind_vector(11)), na.value = "white")

clonalOverlap(TCD8, cloneCall = "gene+nt", method = "overlap") + 
  scale_fill_gradientn(colors = rev(colorblind_vector(11)), na.value = "white")

clonalDiversity(CD8, cloneCall = "gene+nt", group = "cluster")

clonalDiversity(CD8, cloneCall = "gene+nt", group = "RespTmp")

tiff("./plots/clonoDivClus.tiff", width = 5*1500, height = 5*800, res = 300, pointsize = 5)     
clonalDiversity(CD8, cloneCall = "gene+nt", group.by = "RespTmp",
                x.axis = "cluster", n.boots = 100)
dev.off()

alluvialClonotypes(CD8, cloneCall = "gene", 
                   y.axes = c("cluster", "group_id"), 
                   color = "cluster") 

CD8@meta.data$cluster_id <- CD8@active.ident



tiff("./plots/clonoUMAP.tiff", width = 5*500, height = 5*300, res = 300, pointsize = 5)     
p <- clonalOverlay(CD8, reduction = "umap", 
                   freq.cutpoint = 30, bins = 10, facet = "group_id") + 
  guides(color = FALSE) 
p
dev.off()


library(RColorBrewer)
cols<- brewer.pal(5, "Reds")
pal <- colorRampPalette(cols)

slot(CD8, "meta.data")$cloneType <- factor(slot(CD8, "meta.data")$cloneType, 
                                           levels = c(NA, "Single (0 < X <= 1)", 
                                                      "Small (1 < X <= 5)",
                                                      "Medium (5 < X <= 20)",
                                                      "Large (20 < X <= 100)",
                                                      "Hyperexpanded (100 < X <= 500)"))

table(CD8@meta.data$cloneType, CD8@active.ident)
meta <- melt(table(CD8@meta.data[!is.na(CD8@meta.data$Frequency),
                                 c("functional.cluster", "cloneType")]), varnames = c("functional.cluster", "cloneType"))
ggplot(data = meta, aes(x = CD8@active.ident, y = value, fill = cloneType)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_fill_manual(values = c(palette(5)),
                                                                                            na.value = "grey")


CD8_split <- SplitObject(CD8, split.by = "group_id")
CD8_R <- CD8_split$Res
CD8_NR <- CD8_split$NonRes



mark_res <- FindAllMarkers(CD8_R)
mark_nr <- FindAllMarkers(CD8_NR)
mark_res %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

mark %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

tiff("./plots/DoHeat1.tiff", width = 5*400, height = 5*500, res = 150, pointsize = 5)     
DoHeatmap(object = CD8_R, features = top20$gene, size = 8, label = FALSE) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 14)) + ggtitle("Top 20")
dev.off()

tiff("./plots/DoHeatNR.tiff", width = 5*400, height = 5*500, res = 150, pointsize = 5)     
DoHeatmap(object = CD8_NR,  group.by = "cloneType", features = top20$gene, size = 8, label = FALSE) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) + 
  theme(text = element_text(size = 14)) + ggtitle("Top 20")
dev.off()

# get log normalized values for top 50 variable genes
to_plot <- GetAssayData(CD8_R, slot = "data")[top20$gene, ]

cluster_anno<- CD8@meta.data$cloneType


Heatmap(as.matrix(to_plot),
        show_column_names = TRUE)


comb2 <- expression2List(CD8, split.by = "cluster")

length(comb2) 

tiff("./plots/diversity.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
clonalDiversity(comb2, cloneCall = "nt")
dev.off()
library(scRepertoire)

CD8_list <- SplitObject(CD8, split.by = "group_id")

tiff("./plots/occupiedResp.tiff", width = 5*450, height = 5*150, res = 300, pointsize = 5)     
occupiedscRepertoire(CD8_list$Res, x.axis = "cluster") 
dev.off()

tiff("./plots/occupiedNonResp.tiff", width = 5*450, height = 5*150, res = 300, pointsize = 5)     
occupiedscRepertoire(CD8_list$NonRes, x.axis = "cluster") 
dev.off()

#ggraph needs to be loaded due to issues with ggplot
library(ggraph)

clonalHomeostasis(comb2, cloneCall = "nt")
clonalProportion(comb2, cloneCall = "nt")

tiff("./plots/overlap.tiff", width = 5*350, height = 5*300, res = 300, pointsize = 5)     
clonalOverlap(comb2, cloneCall="aa", method="overlap")
dev.off()

save(list = ls(), file ="scRNAseq_CD8_step4.rds")
load("scRNAseq_CD8_step4.rds")


save(CD8, model, file = "CD8.rds")
load("CD8.rds")




