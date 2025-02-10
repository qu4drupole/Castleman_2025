# adding some additional editing and cleaning of the data
setwd("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration")
library(Seurat)
library(tidyverse)
library(cowplot)
Sys.setenv(RETICULATE_PYTHON = "/mnt/isilon/cscb/software/anaconda3/envs/scanpy-basic/bin/python3")
library(reticulate)
library(future)
library(here)
plan(multisession, workers = 4)
options(future.globals.maxSize= 50000*1024^2)

here::i_am("code/02-annotation_DGS/Additional_editing.R")

# if coming from symphony...
seu <- gex.integrated

# else...
seu <- readRDS("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration/processed/02-annotation/symphony_transfer_integrated_DGS.RDS")

#####
# Removing the mysterious "B cell subset"

DimPlot(seu, group.by = "pred_anno")
seu <- subset(seu, subset = pred_anno != "B cell subset")

#####
# Renaming some stuff
seu$pred_anno <- as.character(seu$pred_anno)
seu$pred_anno[seu$pred_anno=="Classic Dendritic cells 1"] <- "Classic Dendritic cells"
seu$pred_anno[seu$pred_anno=="T follicular helper cells"] <- "T Follicular Helper cells"
seu$pred_anno[seu$pred_anno=="Classic Dendritic cells 2 and 3"] <- "Monocytes"

#####
# Should I remove high ribo cells?
# enforce naming within clusters?
FeaturePlot(seu, features = c("percent.rps", "nFeature_SCT"))
quantile(seu$percent.rps, c(0.9,0.95,0.99,1))
quantile(seu$nFeature_SCT, c(0,0.01,0.05,0.1))
ncol(seu)
seu <- subset(seu, subset = percent.rps < 13 & nFeature_SCT > 550)


DimPlot(seu, group.by = "pred_anno", label = T, label.box = T, repel = T) + theme(legend.position = "none")
DimPlot(seu, group.by = "seurat_clusters", label = T, label.box = T, repel = T) + theme(legend.position = "none")
table(seu$seurat_clusters[seu$pred_anno == "Stromal cells"])
table(seu$pred_anno[seu$seurat_clusters==4])

#####
# Re-clustering
ElbowPlot(seu)

seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, algorithm=4, method = "igraph", resolution=1)
# DimPlot(seu, group.by = "seurat_clusters")


#####
# More metadata
seu$edit.ident <- seu$orig.ident
seu$edit.ident[seu$edit.ident == "KFD1_3seq"] <- "KFD2_3seq"

cell_names <- seu$orig.ident
chem <- sapply(cell_names, function(x) strsplit(x, "_")[[1]][2])
seu$seq_chem <- chem

disease_rep <- sapply(cell_names, function(x) strsplit(x, "_")[[1]][1])
disease <- sapply(disease_rep, function(s) sub("\\d", "", s))
seu$sample_name <- disease_rep
seu$disease <- disease

saveRDS(seu, "processed/02-annotation/symphony_transfer_integrated_DGS_edited.RDS")

#####
# Additional editing, 5/30/23
#   -These changes may be best left in. Instead restrict DE to gene sets excluding mito, ribo, etc. genes.

seu <- readRDS("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration/processed/02-annotation/symphony_transfer_integrated_DGS_edited.RDS")

DefaultAssay(seu) <- "SCT"

# Dropping ribo and mito genes
drop_genes <- grep("^MT-", rownames(seu), value = T)
drop_genes <- c(drop_genes, grep("^RP[LS]|^MRPL", rownames(seu), value = T))
keep <- rownames(seu)[!(rownames(seu) %in% drop_genes)]

seu <- subset(seu, features = keep)
saveRDS(seu, "processed/02-annotation/symphony_transfer_integrated_DGS_edited.RDS")



#####
# Additional editing, 7/28/2023
#   -new UMAP with just KFD and R
seu <- readRDS("processed/02-annotation/symphony_transfer_integrated_DGS_edited.RDS")
Idents(seu) <- 'disease'
seu.sub <- subset(seu, idents = c("R", "KFD"))
seu.sub$grade <- "R"
seu.sub$grade[seu.sub$sample_name=="KFD2"] <- "KFD_peak"
seu.sub$grade[seu.sub$sample_name=="KFD3"] <- "KFD_early"

dat_list <- SplitObject(seu.sub, split.by = "edit.ident")
gex_list <- list()
for(s in seq_along(dat_list)){
  DefaultAssay(dat_list[[s]]) <- "RNA"
  gex_list[[s]] <- DietSeurat(dat_list[[s]], assays="RNA")
  gex_list[[s]] <- SCTransform(gex_list[[s]], vst.flavor = "v2") %>% RunPCA(npcs = 30, verbose = FALSE)
}

features <- SelectIntegrationFeatures(object.list = gex_list, nfeatures = 2000)
gex_list <- PrepSCTIntegration(object.list = gex_list, anchor.features = features)
gex.anchors <- FindIntegrationAnchors(object.list = gex_list, normalization.method = "SCT",
                                         anchor.features = features, reduction = 'rpca')
gex.combined.sct <- IntegrateData(anchorset = gex.anchors, normalization.method = "SCT")

gex.combined.sct <- RunPCA(gex.combined.sct, verbose = FALSE)
gex.combined.sct <- RunUMAP(gex.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
gex.combined.sct <- RunTSNE(gex.combined.sct, dims = 1:30)


###############
# UMAP

pdf(
  here("outs/KFD/Dimension_reductions/UMAPs.pdf"),
  height = 8,
  width = 11
)

p1 <- DimPlot(gex.combined.sct, reduction = "umap", group.by = "edit.ident")
p2 <- DimPlot(gex.combined.sct, reduction = "umap", group.by = "pred_anno", 
        label = T, label.box = T, label.size = 2, repel = T) + theme(legend.position = "none")
p3 <- FeaturePlot(gex.combined.sct, features = 'nCount_SCT', reduction = "umap")
p4 <- FeaturePlot(gex.combined.sct, features = 'percent.mt', reduction = "umap")
p5 <- FeaturePlot(gex.combined.sct, features = 'percent.rps', reduction = "umap")

main_umaps <- plot_grid(p1,p2, ncol = 2)
bottom_umaps <- plot_grid(p3,p4,p5, nrow = 1)
umaps <- plot_grid(main_umaps, bottom_umaps, nrow = 2, rel_heights = c(1,0.5))
ptitle <- ggdraw() + 
  draw_label("KFD, R: UMAPs",
             fontface = 'bold',x = 0,hjust = 0
  ) + theme(plot.margin = margin(0, 0, 0, 7))
plot_grid(ptitle, umaps, ncol = 1, rel_heights = c(0.1, 1))

dev.off()

###############
# t-SNE

pdf(
  here("outs/KFD/Dimension_reductions/TSNEs.pdf"),
  height = 8,
  width = 11
)

p1 <- DimPlot(gex.combined.sct, reduction = "tsne", group.by = "edit.ident")
p2 <- DimPlot(gex.combined.sct, reduction = "tsne", group.by = "pred_anno", 
              label = T, label.box = T, label.size = 2, repel = T) + theme(legend.position = "none")
p3 <- FeaturePlot(gex.combined.sct, features = 'nCount_SCT', reduction = "tsne")
p4 <- FeaturePlot(gex.combined.sct, features = 'percent.mt', reduction = "tsne")
p5 <- FeaturePlot(gex.combined.sct, features = 'percent.rps', reduction = "tsne")

main_umaps <- plot_grid(p1,p2, ncol = 2)
bottom_umaps <- plot_grid(p3,p4,p5, nrow = 1)
umaps <- plot_grid(main_umaps, bottom_umaps, nrow = 2, rel_heights = c(1,0.5))
ptitle <- ggdraw() + 
  draw_label("KFD, R: UMAPs",
             fontface = 'bold',x = 0,hjust = 0
  ) + theme(plot.margin = margin(0, 0, 0, 7))
plot_grid(ptitle, umaps, ncol = 1, rel_heights = c(0.1, 1))

dev.off()

saveRDS(gex.combined.sct, "/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration/processed/02-annotation/symphony_transfer_integrated_DGS_edited_onlyKFD-R.RDS")


################
# Making some barcodes to match scVelo
seu <- readRDS("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration/processed/02-annotation/symphony_transfer_integrated_DGS_edited.RDS")
###
# apparently MCD2_5seq is mislabelled
seu$edit.ident[seu$edit.ident=="MCD2_5seq"] <- "MCD3_5seq"
saveRDS(seu, "/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration/processed/02-annotation/symphony_transfer_integrated_DGS_edited.RDS")

meta <- seu@meta.data
raw_bc <- sapply(rownames(meta), function(s) strsplit(s,"-")[[1]][1])
meta$velo_bc <- paste(meta$edit.ident, raw_bc, sep = "_")
rownames(meta) <- meta$velo_bc
write.csv(meta, "code/03-trajectory/seu_meta.csv")

umap_coord <- seu@reductions$umap@cell.embeddings
rownames(umap_coord) <- meta$velo_bc
write.csv(umap_coord, "code/03-trajectory/seu_integrated_umap.csv")

###
# investigating some count discrepancies between velocyto and CellRanger, the 5' samples seem very poor
seu_list <- SplitObject(seu, split.by = "edit.ident")

VlnPlot(seu_list$MCD4_5seq, features = "nCount_RNA", group.by = "edit.ident")

bc <- sapply(meta$velo_bc, function(s) strsplit(s,"_")[[1]][3])
meta$raw_bc <- bc

table(meta$edit.ident)
sum(meta$raw_bc[meta$edit.ident=="KFD2_5seq"] %in% meta$raw_bc[meta$edit.ident!="KFD2_5seq"])
tail(meta$raw_bc[meta$edit.ident=="KFD2_5seq"])


