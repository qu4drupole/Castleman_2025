library(AUCell)
library(dplyr)
library(reshape2)
library(Seurat)
library(org.Hs.eg.db)
library(cowplot)
library(reactome.db)
library(here)
library(ggplot2)
library(here)

setwd("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration")

seu.stromal.sub <- readRDS("processed/02-annotation/stromal_subtypes.RDS")

i_am("code/02-annotation_DGS/Stromal_subtyping.R")
# Heatmap showing gene expression for Seurat-defined clusters of known markers for
# 
# Fibroblastic stromal cells (PDGFRA, PDGFRB, CXCL13, APOE, CCL21, CCL19, and PDPN), 
# blood endothelial cells (CDH5, ENG, CD34, PECAM1), and 
# lymphatic endothelial cells markers (PROX1, PECAM1, PDPN). 
# 
# Heatmap showing gene expression of known markers for 
# fibroblastic stromal cell subsets, including 
# ACTA2+ perivascular reticular cells (ACTA2, TAGLN, TPM2, PDGFRB), 
# vascular smooth muscle cells (ACTA2, MYH11, MCAM), 
# CCL19hi T-zone fibroblastic reticular cells (CCL19, CCL21, CXCL12, CXCL9), 
# CCL19lo T-zone fibroblastic reticular cells (LUM, DCN, PDPN, PDGFRA), 
# follicular dendritic cells (CXCL13, CLU, FDCSP, DES). 

seu <- readRDS("processed/02-annotation/symphony_transfer_integrated_DGS_edited.RDS")
Idents(seu) <- "pred_anno"
unique(seu$pred_anno)
seu.stromal <- subset(seu, idents=c("Follicular dendritic cells","Fibroblasts","Stromal cells",
                                    "Endothelial cells","Lymphatics"))

DefaultAssay(seu.stromal) <- "RNA"
seu.stromal <- FindVariableFeatures(seu.stromal)
seu.stromal <- ScaleData(seu.stromal)
seu.stromal <- RunPCA(seu.stromal, verbose = FALSE)
seu.stromal <- FindNeighbors(seu.stromal, dims = 1:10)
seu.stromal <- FindClusters(seu.stromal, resolution = 0.5)
seu.stromal <- RunUMAP(seu.stromal, reduction = "pca", dims = 1:30, verbose = FALSE)
DimPlot(seu.stromal, label = T, label.box = T, repel = T)
DimPlot(seu.stromal, group.by = "disease")
DimPlot(seu.stromal, split.by = "disease", ncol = 2)

markers <- unique(c('PDGFRA', 'PDGFRB', 'CXCL13', 'APOE', 'CCL21', 'CCL19', 'PDPN',
             'CDH5', 'ENG', 'CD34', 'PECAM1',
             'PROX1', 'PECAM1', 'PDPN',
             'ACTA2', 'TAGLN', 'TPM2', 'PDGFRB',
             'ACTA2', 'MYH11', 'MCAM',
             'CCL19', 'CCL21', 'CXCL12', 'CXCL9',
             'LUM', 'DCN', 'PDPN', 'PDGFRA',
             'CXCL13', 'CLU', 'FDCSP', 'DES'))

marker_hit <- sapply(markers, function(m) grep(m, rownames(seu.stromal), ignore.case = T, value = T))
# they all have exact hits

DotPlot(seu.stromal, features = markers, assay = "SCT") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

DimPlot(seu.stromal, group.by = "pred_anno")
FeaturePlot(seu.stromal, "nFeature_SCT")

# 0. Fibroblastic stromal cell, type 1 (HVCD, R)
# 1. Fibroblastic stromal cell, type 1
# 2. ACTA2+ perivascular reticular cells
# 3. lymphatic endothelial cell
# 4. Fibroblastic stromal cell, type 2 (MCD)
# 5. Fibroblastic stromal cell, type 3 (MCD, PDGFRAhi, CXCL12hi)
# 6. blood endothelial cell
# 7. Fibroblastic stromal cell, type 2 (strange spread...)
# 8. FDC
# 9. low quality?
# 10. FDC (HVCD)
# 11. low quality
# 12. Fibroblastic stromal cell, type 1 (KFD)

new_cell_type <- c("Fibroblastic stromal cell, type 1", "Fibroblastic stromal cell, type 1", "ACTA2+ perivascular reticular cells",
                   "lymphatic endothelial cell", "Fibroblastic stromal cell, type 2", "Fibroblastic stromal cell, type 2",
                   "blood endothelial cell", "Fibroblastic stromal cell, type 2", "FDC", "drop", "FDC", "drop",
                   "Fibroblastic stromal cell, type 1")
names(new_cell_type) <- seq(0,12)
seu.stromal$pred_anno_2 <- new_cell_type[seu.stromal$seurat_clusters]
seu.stromal.sub <- subset(seu.stromal, subset = pred_anno_2 != "drop")

marker_list <- c("IL6", "IL6R", "IL6ST", "VEGFA","VEGFB","VEGFC","FLT1","IL1B","OSMR","LIFR")

DimPlot(seu.stromal.sub, group.by = "pred_anno_2")

DefaultAssay(seu.stromal.sub) <- "SCT"

pdf(
  here("code/02-annotation_DGS/outs/stromal_subtype/UMAP_features_of_interest_celltype_split.pdf"),
  height = 8,
  width = 11
)

for(m in marker_list){
  p <- FeaturePlot(seu.stromal.sub, m, split.by = "pred_anno_2", cols = c("lightgrey", "red"),
                   pt.size = 0.5, combine = F)
  print(plot_grid(plotlist = p, ncol = 3))
}

dev.off()

pdf(
  here("code/02-annotation_DGS/outs/stromal_subtype/UMAP_features_of_interest_disease_split.pdf"),
  height = 8,
  width = 11
)

for(m in marker_list){
  p <- FeaturePlot(seu.stromal.sub, m, split.by = "disease", cols = c("lightgrey", "red"),
                   pt.size = 0.5, combine = F)
  print(plot_grid(plotlist = p, ncol = 2))
}

dev.off()

saveRDS(seu.stromal.sub, "processed/02-annotation/stromal_subtypes.RDS")

### more plots

pdf(
  here("code/02-annotation_DGS/outs/stromal_subtype/violin_markerlist.pdf"),
  height = 8,
  width = 11
)

for(m in marker_list){
  print(VlnPlot(seu.stromal.sub, features = m, group.by = "pred_anno_2", assay = "RNA"))
}

dev.off()

### drop KFD
Idents(seu.stromal.sub)<-"pred_anno_2"
seu.stromal.cd <- subset(seu.stromal.sub, subset=disease=="KFD",invert=T)
DefaultAssay(seu.stromal.cd)<-"RNA"
seu.stromal.cd <- FindVariableFeatures(seu.stromal.cd)
# seu.stromal.cd <- ScaleData(seu.stromal.cd)
seu.stromal.cd <- RunPCA(seu.stromal.cd, verbose = FALSE)
seu.stromal.cd <- FindNeighbors(seu.stromal.cd, dims = 1:10)
seu.stromal.cd <- FindClusters(seu.stromal.cd, resolution = 0.5)
seu.stromal.cd <- RunUMAP(seu.stromal.cd, reduction = "pca", dims = 1:15, verbose = FALSE)
DimPlot(seu.stromal.cd, label = T, label.box = T, repel = T)
DimPlot(seu.stromal.cd, group.by = "pred_anno_2")
DimPlot(seu.stromal.cd, group.by = "disease")


###
# new dotplot with abbreviations
unique(seu.stromal.sub$pred_anno_2)
seu.stromal.sub$pred_anno_3 <- as.character(seu.stromal.sub$pred_anno_2)

new_names <- c("FSC, type 1","ACTA2+ PRC","LEC",
               "FDC","BEC","FSC, type 2")
names(new_names) <- unique(unique(seu.stromal.sub$pred_anno_2))
seu.stromal.sub$pred_anno_3 <- new_names[seu.stromal.sub$pred_anno_3]
Idents(seu.stromal.sub) <- "pred_anno_3"
DotPlot(seu.stromal.sub, features = markers, assay = "SCT", scale = F) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
