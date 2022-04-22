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
library(tidyr)
library(mclust)
library(SCPA)
library(dyno)
library(scRepertoire)
library(SeuratWrappers)
library(monocle3)
library(hdf5r)
library(circlize)
library(ComplexHeatmap)
library(tibble)
library(msigdbr)
library(statmod)
library(scater)
library(SCPA)
library(mgcv)
library(gam)
library(tradeSeq)
library(condiments)
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

pathways <- "combined_metabolic_pathways.csv"

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

scpa_out <- scpa_out %>% 
   reduce(full_join, by = "Pathway") %>% 
   set_colnames(gsub(colnames(.), pattern = " ", replacement = "_")) %>%
   select(c("Pathway", grep("_qval", colnames(.)))) %>%
   filter_all(any_vars(. > 2)) %>%
   column_to_rownames("Pathway")

blood_paths <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
                 "HALLMARK_COMPLEMENT", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                 "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                 "HALLMARK_MTORC1_SIGNALING", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                 "HALLMARK_MYC_TARGETS_V1", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")

position <- which(rownames(scpa_out) %in% blood_paths)
row_an <- rowAnnotation(Genes = anno_mark(at = which(rownames(scpa_out) %in% blood_paths),
                                          labels = rownames(scpa_out)[position],
                                          labels_gp = gpar(fontsize = 7),
                                          link_width = unit(2.5, "mm"),
                                          padding = unit(1, "mm"),
                                          link_gp = gpar(lwd = 0.5)))
col_hm <- colorRamp2(colors = c("blue", "white", "red"), breaks = c(0, 3, 6))


Heatmap(scpa_out,
        col = col_hm,
        name = "Qval",
        show_row_names = F,
        right_annotation = row_an,
        column_names_gp = gpar(fontsize = 8),
        border = T,
        column_km = 3,
        row_km = 3, 
        column_labels = c("HD", "NonRes", "Res"))

#TRAJECTORY INFERENCE using Slingshot
clusterLabels <- CD8@active.ident
sceCD8 <- as.SingleCellExperiment(CD8, assay = "RNA")
sds <- slingshot(sceCD8, clusterLabels = clusterLabels, 
                 allow.breaks = TRUE, stretch = 2, reducedDim = "UMAP", start.clus = "Naive") #Calcualting the trajectory
sds <- SlingshotDataSet(sds)
saveRDS(sds, "sds")
#Create dataframe pseudotimes
pt <- as.data.frame(slingPseudotime(sds))
colnames(pt) <- c("pseudoT1", "pseudoT2")
names <- rownames(pt)
   
#Extract the curves coordinates
curve1 <- slingCurves(sds)[[1]]
curve2 <- slingCurves(sds)[[2]]
   
#Create dataframe with UMAP coordinates
umap_df <- reducedDim(sceCD8, "UMAP")
   
#Create dataframe with UMAP coordinates and pseudotimes
df <- cbind(umap_df, pt)
   
#Plot according to pseudotime values  (pseudoT1)
p1 <- ggplot(df, aes(UMAP_1, UMAP_2)) +
   geom_point(aes_string(color = df$pseudoT1),
              alpha = 0.5) +
   scale_colour_viridis_c() +
   theme_minimal() + labs(colour = "Pseudotime") 
   
   
#Plot 1st lineage
tiff("./plots/traj1.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
p1 <- p1 + geom_path(aes(x = UMAP_1, y = UMAP_2), data = curve1$s[curve1$ord, ] %>% as.data.frame(),
                     col = "black", size = 1, arrow = arrow(), lineend = "round") + scale_color_viridis_c() +  labs(color = "Pseudotime") + ggtitle("Senescence")
dev.off()
   
#Plot according to pseudotime (pseudoT2)
p2 <- ggplot(df, aes(UMAP_1, UMAP_2)) +
   geom_point(aes_string(color = df$pseudoT2),
              alpha = 0.5) +
   scale_colour_viridis_c() +
   theme_minimal() + labs(colour = "Pseudotime")
   
tiff("./plots/traj2.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
p2 <- p2 +   geom_path(aes(x = UMAP_1, y = UMAP_2), data = curve2$s[curve2$ord, ] %>% as.data.frame(),
                       col = "black", size = 1, arrow = arrow(), lineend = "round") + scale_color_viridis_c() +  labs(color = "Pseudotime") + ggtitle("Exhaustion")
dev.off()
   
tiff("./plots/traj.tiff", width = 8*130, height = 5*300, res = 150, pointsize = 5)  
plot_grid(p1, p2, ncol = 1, align = "h")
dev.off()

#One dimensional clustering with Mclust
pt1 <- pt %>% select(pseudoT1) %>% drop_na()
pt1.mc <- Mclust(pt1, G=1:4)$classification %>% 
   as.data.frame() %>% 
   set_colnames("milestone") %>% rownames_to_column("cell")

df_mile <- cbind(pt, pt1.mc[, "milestone"][match(rownames(pt), pt1.mc$cell)]) %>% set_colnames(c("pt1", "pt2", "milestone1"))
CD8$milestone1 <- df_mile$milestone1

CD8_pseudo1 <- list()
for (i in 1:max(df_mile$milestone1, na.rm = TRUE)) {
   CD8_pseudo1[[i]] <- seurat_extract(CD8, meta1 = "milestone1", value_meta1 = i)
}

pathways <- "combined_metabolic_pathways.csv"

CD8_pseudo1 <- lapply(CD8_pseudo1, function(x) {
   ss <- strsplit(rownames(x), ".", fixed = TRUE)
   rownames(x) <- sapply(ss, .subset, 2)
   return(x)})

CD8.meta.pt1 <- compare_pathways(samples = CD8_pseudo1, 
                                   pathways = pathways)

CD8.meta.pt1 <- CD8.meta.pt1 %>%
   data.frame() %>%
   select(Pathway, qval) %>% remove_rownames() %>% 
   column_to_rownames("Pathway")
View(CD8.meta.pt1)
col_hm <- colorRamp2(colors = c("white", "red"), breaks = c(0, max(CD8.meta.pt1$qval)))

paths <- c("REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
           "REACTOME_METABOLISM_OF_CARBOHYDRATES", "REACTOME_FATTY_ACID_METABOLISM", "HALLMARK_FATTY_ACID_METABOLISM",
           "REACTOME_SYNTHESIS_OF_VERY_LONG_CHAIN_FATTY_ACYL_COAS", "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE")
position <- which(rownames(CD8.meta.pt1) %in% paths)

row_an <- HeatmapAnnotation(Genes = anno_mark(at = which(rownames(CD8.meta.pt1) %in% paths),
                                          labels = rownames(CD8.meta.pt1)[position],
                                          labels_gp = gpar(fontsize = 5),
                                          link_width = unit(0, "mm")))
                                          #link_gp = gpar(lwd = 0.5)))

#Annotation to be added(!!!)
Heatmap(t(CD8.meta.pt1),
        name = "Qvalue",
        col = col_hm,
        bottom_annotation = row_an,
        border = T,
        rect_gp = gpar(col = "white", lwd = 0.1),
        heatmap_height = unit(10, "cm"),
        show_column_dend = F,
        show_row_names = F,
        show_column_names = F,
        column_title = "Metabolic pathway (senescence trajectory") 

pt2 <- pt %>% select(pseudoT2) %>% drop_na()
pt2.mc <- Mclust(pt2, G=1:4)$classification %>% 
   as.data.frame() %>% 
   set_colnames("milestone") %>% rownames_to_column("cell")

df_mile2 <- cbind(pt, pt2.mc[, "milestone"][match(rownames(pt), pt2.mc$cell)]) %>% set_colnames(c("pt1", "pt2", "milestone2"))
CD8$milestone2 <- df_mile2$milestone2

CD8_pseudo2 <- list()
for (i in 1:max(df_mile2$milestone2, na.rm = TRUE)) {
   CD8_pseudo2[[i]] <- seurat_extract(CD8, meta1 = "milestone2", value_meta1 = i)
}

pathways <- "combined_metabolic_pathways.csv"

CD8_pseudo2 <- lapply(CD8_pseudo2, function(x) {
   ss <- strsplit(rownames(x), ".", fixed = TRUE)
   rownames(x) <- sapply(ss, .subset, 2)
   return(x)})

CD8.meta.pt2 <- compare_pathways(samples = CD8_pseudo2, 
                                 pathways = pathways)

CD8.meta.pt2 <- CD8.meta.pt2 %>%
   data.frame() %>%
   select(Pathway, qval) %>% remove_rownames() %>% 
   column_to_rownames("Pathway")

col_hm <- colorRamp2(colors = c("white", "red"), breaks = c(0, max(CD8.meta.pt2$qval)))

paths <- c("REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
           "REACTOME_METABOLISM_OF_CARBOHYDRATES", "REACTOME_FATTY_ACID_METABOLISM", "HALLMARK_FATTY_ACID_METABOLISM",
           "REACTOME_SYNTHESIS_OF_VERY_LONG_CHAIN_FATTY_ACYL_COAS", "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE")
position <- which(rownames(CD8.meta.pt2) %in% paths)

row_an <- HeatmapAnnotation(Genes = anno_mark(at = which(rownames(CD8.meta.pt2) %in% paths),
                                              labels = rownames(CD8.meta.pt2)[position],
                                              labels_gp = gpar(fontsize = 5),
                                              link_width = unit(0, "mm")))


Heatmap(t(CD8.meta.pt2),
        name = "Qvalue",
        col = col_hm,
        bottom_annotation = row_an,
        border = T,
        rect_gp = gpar(col = "white", lwd = 0.1),
        heatmap_height = unit(10, "cm"),
        show_column_dend = F,
        show_row_names = F,
        show_column_names = F,
        column_title = "Metabolic pathway (exhaustion trajectory)")

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
   
# fit a GAM with a loess term for pseudotime pt1
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
par(mar=c(1,1,1,1))
icMat <- evaluateK(counts = counts(sceCD8), sds = sds, k = 3:7, 
                  nGenes = 100, verbose = T, plot = TRUE)
print(icMat[1:2,])
set.seed(7)

###Parallel computinh
BPPARAM <- BiocParallel::bpparam()
BPPARAM # lists current options
BPPARAM$workers <- 2

#fitGAM
sce <- fitGAM(counts = counts(sceCD8), sds = sds, nknots = 5, verbose = TRUE, BPPARAM = BPPARAM)
   
# plot our Slingshot lineage trajectories, this time illustrating the new tradeSeq knots
tiff("./plots/traj.tiff", width = 5*500, height = 5*300, res = 300, pointsize = 5)     
plotGeneCount(curve = sds, counts = counts,
             clusters = CD8@active.ident,
             models = sce)
dev.off()
   
#Association test
assoRes <- associationTest(sce)
head(assoRes)
   
### Discovering differentiated cell type markers
# discover marker genes for the differentiated cell types
endRes <- diffEndTest(sce) #Nothing interesting with this analysis
head(endRes)
   
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[2]]
plotSmoothers(sce, counts(sce), sigGene)
   
plotGeneCount(sds, counts, gene = sigGene)
   
earlyDERes <- earlyDETest(sce, knots = c(3, 4))
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
head(rownames(earlyDERes)[oEarly])
plotSmoothers(sce, counts(sce), gene = rownames(earlyDERes)[oEarly][3])
plotGeneCount(sds, counts(sce), gene = rownames(earlyDERes)[oEarly][3])

#Genes with different expression patterns
patternRes <- patternTest(sce)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
rownames(patternRes)[oPat][1]

tiff("./plots/DEtrajGZMK.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
plotSmoothers(sce, counts(sce), gene = rownames(patternRes)[oPat][233]) + ggtitle ("GZMK")
dev.off()

tiff("./plots/DEtrajGZMKumap.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
plotGeneCount(sds, counts(sce), gene = rownames(patternRes)[oPat][1]) + ggtitle ("GZMK")
dev.off()

tiff("./plots/DEtrajGNLY.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
plotSmoothers(sce, counts(sce), gene = rownames(patternRes)[oPat][7]) + ggtitle ("GNLY")
dev.off()

tiff("./plots/DEtrajGNLYumap.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
plotGeneCount(sds, counts(sce), gene = rownames(patternRes)[oPat][7]) + ggtitle ("GNLY")
dev.off()

tiff("./plots/DEtrajPRF1.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
plotSmoothers(sce, counts(sce), gene = rownames(patternRes)[oPat][8]) + ggtitle ("PRF1")
dev.off()

tiff("./plots/DEtrajPRF1umap.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
plotGeneCount(sds, counts(sce), gene = rownames(patternRes)[oPat][8]) + ggtitle ("PRF1")
dev.off()

tiff("./plots/DEtrajGZMB.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
plotSmoothers(sce, counts(sce), gene = rownames(patternRes)[oPat][9]) + ggtitle ("GZMB")
dev.off()

tiff("./plots/DEtrajGZMBumap.tiff", width = 8*100, height = 5*100, res = 150, pointsize = 5)  
plotGeneCount(sds, counts(sce), gene = rownames(patternRes)[oPat][9]) + ggtitle ("GZMB")
dev.off()

#Finish analysis whole dataset
sce.cond <- subset(sce,,group_id %in% c("Res", "NonRes"))
sds <- slingshot(sce.cond, reducedDim = 'UMAP', clusterLabels = colData(sce.cond)$ident,
                    start.clus = 'Tnaive', approx_points = 150)
sds <- SlingshotDataSet(sds)

#Trajectory by condition (compare Res and NonRes)
df <- bind_cols(
 as.data.frame(reducedDim(sce.cond, "UMAP")),
 as.data.frame(colData(sce.cond))
) %>%
 sample_frac(1)

scores <- condiments::imbalance_score(
 Object = df %>% select(UMAP_1, UMAP_2) %>% as.matrix(), 
 conditions = df$group_id)
   
df$scores <- scores$scores
p3 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scores)) +
 geom_point(size = .7) +
 scale_color_viridis_c(option = "C") +
 labs(col = "Scores")
p3

df$scaled_scores <- scores$scaled_scores
p4 <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = scaled_scores)) +
  geom_point(size = .7) +
  scale_color_viridis_c(option = "C") +
  labs(col = "Scores")
p4

topologyTest(SlingshotDataSet(sds), sce.cond$group_id, rep = 100,
             methods = "KS_mean", threshs = .01)
knitr::kable(top_res)

psts <- slingPseudotime(sds) %>%
  as.data.frame() %>%
  mutate(cells = rownames(.),
         conditions = df$group_id) %>%
  pivot_longer(starts_with("Lineage"), values_to = "pseudotime", names_to = "lineages")

ggplot(psts, aes(x = pseudotime, fill = conditions)) +
  geom_density(alpha = .5) +
  scale_fill_brewer(type = "qual") +
  facet_wrap(~lineages) +
  theme(legend.position = "bottom")

prog_res <- progressionTest(sds, conditions = df$group_id, global = TRUE, lineages = TRUE)
knitr::kable(prog_res)

df$weight_1 <- slingCurveWeights(sds, as.probs = TRUE)[, 1]
ggplot(df, aes(x = weight_1, fill = group_id)) +
  geom_density(alpha = .5) +
  scale_fill_brewer(type = "qual") +
  labs(x = "Curve weight for the first lineage")

dif_res <- differentiationTest(sds, conditions = df$group_id, global = FALSE, pairwise = TRUE)
knitr::kable(dif_res)

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

View(CD8_pseudo)
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




