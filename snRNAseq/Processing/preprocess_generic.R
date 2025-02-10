# /mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration

library(reshape2)
library(stringr)
library(SoupX)
library(DoubletFinder)
library(Seurat)
library(SeuratDisk)
Sys.setenv(RETICULATE_PYTHON = "/mnt/isilon/cscb/software/anaconda3/envs/scanpy-basic/bin/python3")
library(reticulate)

library(here)
library(ggplot2)
library(cowplot)
# library(cscb.tools)
# too lazy to reinstall my edits from github. source cscb.tools from my local library (i.e. copy/paste)
library(dplyr)
library(parallel)

# # turning off future for integration
library(future)
plan("multisession", workers = 4)
options(future.globals.maxSize= 500000*1024^2)
#########################################
# Restarting the entire processing and integration for Vinodh's data
# Things were too messy and confusing
#
# NOTE: parallel and future do not work well together. Had to adjust 'future' plan throughout script 
#########################################

setwd("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration")

here::i_am("code/01-process_DGS/preprocess_generic.R")
source("../../SCTC-VP-22/vinodh-5-3-integration/code/01-5_3-comp/SoupX_load10x.R")

# need to track down the CR outputs from the various sample preparations
# SCTC-VP-3
# "/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-3/cscb-lymph-vinodh/code/02-counts/"
VP_3_pathnames <- c("KD1_wIntrons", "MCD1-redo_wIntrons-force", "MCD2_wIntrons-force", "MCD3_wIntrons", "NVCD1_wIntrons-force", "R1_wIntrons-force")
VP_3_snames <- c("KFD1_3seq", "MCD1_3seq", "MCD2_3seq", "MCD3_3seq", "HVCD1_3seq", "R1_3seq")

# SCTC-VP-16
# There appear to be two batches? 
# /mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/code/02-multi/runs
#    -KFD1 is not listed in the sample name manifest...oh, skipping bc it's very low quality
#    -R2 looks pretty bad, but will keep anyway
VP_16_1_pathnames <- c("HVCD2_VDJ_GEX","HVCD3_VDJ_GEX","KFD3_VDJ_GEX","MCD4_VDJ_GEX","R2_VDJ_GEX")
VP_16_1_snames <- c("HVCD2_5seq", "HVCD3_5seq", "KFD3_5seq", "MCD4_5seq", "R2_5seq")

# /mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/code/05-multi-pt2/runs
VP_16_2_pathnames <- c("KFD2_VDJ_GEX","MCD2_VDJ_GEX")
VP_16_2_snames <- c("KFD2_5seq", "MCD2_5seq")

# cr_paths <- paste0("data/",snames,"/filtered_feature_bc_matrix")
VP_3_paths <- paste0("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-3/cscb-lymph-vinodh/code/02-counts/",
                     VP_3_pathnames,
                     "/outs")
VP_16_1_paths <- paste0("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/code/02-multi/runs/",
                        VP_16_1_pathnames,
                        "/outs")
VP_16_2_paths <- paste0("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/code/05-multi-pt2/runs/",
                        VP_16_2_pathnames,
                        "/outs")

all_paths <- c(VP_3_paths, VP_16_1_paths, VP_16_2_paths)
all_snames <- c(VP_3_snames, VP_16_1_snames, VP_16_2_snames)
# VP_3_cluster_paths <- paste0(VP_3_paths,"/analysis/clustering/graphclust/clusters.csv")

soupx_obj <- mclapply(all_paths, load10X_dev, mc.cores = 6)

adjust_wrapper <- function(soup){
  soup <- autoEstCont(soup, soupQuantile = 0.8, tfidfMin = 0.8)
  soup <-  adjustCounts(soup)
  return(soup)
}
soupx_obj_adj <- mclapply(soupx_obj, adjust_wrapper, mc.cores = 12)

gex_obj <- list()
for (i in seq_along(all_snames)){
  # raw_obj[[i]] <- Read10X(mtx_folder[i])
  gex_obj[[i]] <- CreateSeuratObject(soupx_obj_adj[[i]], 
                                     project = all_snames[i],
                                     min.cells = 10,
                                     min.features = 200)
}

names(gex_obj) <- all_snames

# adding raw UMI counts
# for (i in 1:length(raw_obj)) {
#   raw_assay <- CreateAssayObject(counts = raw_obj[[i]])
#   
#   # add this assay to the previously created Seurat object
#   gex_obj[[i]][["RAW"]] <- raw_assay
#   
# }

########################
## GEX preprocessing####
########################

for (i in seq_along(all_snames)) {
  gex_obj[[i]] <- get_percent(gex_obj[[i]])
}

pdf(
  here("code/01-process_DGS/outs/sample_summary_stats.pdf"),
  height = 8,
  width = 11
)

for (i in seq_along(all_snames)){
  p1 <- VlnPlot(gex_obj[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rps"), ncol = 4)
  p2 <- FeatureScatter(gex_obj[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p3 <- FeatureScatter(gex_obj[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  p4 <- FeatureScatter(gex_obj[[i]], feature1 = "nCount_RNA", feature2 = "percent.rps")
  bottomrow <- plot_grid(p2,p3,p4, ncol = 3)
  print(plot_grid(p1, bottomrow, ncol = 1))
}
dev.off()

count_cutoffs <- c(40000, 20000, 20000, 40000, 50000, 50000, 15000, 40000, 40000, 25000, 20000, 25000, 20000)
mt_cutoffs <- c(2, 2, 2, 5, 2.5, 5, 15, 15, 15, 15, 15, 5, 20)
rps_cutoffs <- c(15, 4, 3, 5, 20, 20, 9, 20, 15, 12, 10, 12, 10)

# https://github.com/satijalab/seurat/issues/1227
# gex_obj_proc <- mclapply(1:length(gex_obj), function(i) seurat_process(gex_obj[[i]],
#                                                                        counts_name = snames[i],
#                                                                        count.cutoff = count_cutoffs[i],
#                                                                        mito.cutoff = mt_cutoffs[i],
#                                                                        rps.cutoff = rps_cutoffs[i],
#                                                                        cc_adjust = TRUE),
#                           mc.cores = 2)
gex_obj_proc <- lapply(1:length(gex_obj), function(i) seurat_process(gex_obj[[i]],
                                                                       counts_name = snames[i],
                                                                       count.cutoff = count_cutoffs[i],
                                                                       mito.cutoff = mt_cutoffs[i],
                                                                       rps.cutoff = rps_cutoffs[i],
                                                                       cc_adjust = TRUE))

names(gex_obj_proc) <- all_snames

pdf(
  here("code/01-process_DGS/outs/sample_UMAPs.pdf"),
  height = 8,
  width = 11
)

for (seu in gex_obj_proc){
  p1 <- DimPlot(seu)
  p2 <- FeaturePlot(seu, features = "nCount_SCT")
  p3 <- FeaturePlot(seu, features = "percent.mt")
  p4 <- FeaturePlot(seu, features = "percent.rps")
  umaps <- plot_grid(p1,p2,p3,p4, ncol = 2)
  ptitle <- ggdraw() + 
    draw_label(
      seu@project.name,fontface = 'bold',x = 0,hjust = 0
    ) + theme(plot.margin = margin(0, 0, 0, 7))
  print(plot_grid(ptitle, umaps, ncol = 1, rel_heights = c(0.1, 1)))
}
dev.off()

gex_pk <- mclapply(gex_obj_proc, sum_sweep, mc.cores = 12)
# gex_pk <- lapply(gex_obj_proc, sum_sweep)

gex_obj_proc_db <- mclapply(1:length(gex_obj_proc), 
                            function(i) doubFinder(gex_obj_proc[[i]], sweep.stats.sample = gex_pk[[i]]),
                            mc.cores = 12)
names(gex_obj_proc_db) <- all_snames

# for the most part, pretty clean cells
lapply(gex_obj_proc_db, function(x) sum(x$doublet=="Singlet")/ncol(x))

gex_obj_proc_db <- lapply(gex_obj_proc_db, subset, subset = doublet == "Singlet")

saveRDS(gex_obj_proc_db, "processed/01-5_3-comp/processed_samples_list.RDS")

########################
#### INTEGRATION    ####
########################

count_means <- mclapply(gex_obj_proc_db, rowMeans, mc.cores = 6)
lapply(count_means, function(x) quantile(x, seq(0,1,0.2)))
# I guess the GEM well counts are fairly different...

###
# We have 3 "batches" of GEM well preparation. 2 batches of libraries. How does data look, pre-integration?
#    - going to treat prep chemistry as batch
#    - Check if I need to redo sample normalization...probably not?
#!!      -- I should rerun SCTransform on each GEM using scale_factor = min of GEM median counts, then merge
# merging 3'

vp3_dat <- gex_obj_proc_db[grepl("3seq$", names(gex_obj_proc_db))]
vp3_names <- all_snames[grepl("3seq$", all_snames)]

# vp3.merged <- merge(vp3_dat[[1]], vp3_dat[2:length(vp3_dat)], add.cell.ids = vp3_names, merge.data=T)
# vp3.merged <- RunPCA(vp3.merged, features = features)
# ElbowPlot(vp3.merged)
# vp3.merged <- FindNeighbors(vp3.merged, dims = 1:15)
# vp3.merged <- RunUMAP(vp3.merged, dims = 1:15)
# DimPlot(vp3.merged, group.by = "orig.ident")
# 
# pdf(
#   here("code/01-process_DGS/outs/VP3_pre-integration_UMAPs.pdf"),
#   height = 8,
#   width = 11
# )
# 
# p1 <- DimPlot(vp3.merged, group.by = "orig.ident")
# p2 <- FeaturePlot(vp3.merged, features = "nCount_SCT")
# p4 <- FeaturePlot(vp3.merged, features = "nFeature_RNA")
# p5 <- FeaturePlot(vp3.merged, features = "percent.rps")
# p6 <- FeaturePlot(vp3.merged, features = "percent.mt")
# main_umaps <- plot_grid(p1,p2, ncol = 2)
# bottom_umaps <- plot_grid(p4,p5,p6, nrow = 1)
# umaps <- plot_grid(main_umaps, bottom_umaps, nrow = 2, rel_heights = c(1,0.5))
# ptitle <- ggdraw() + 
#   draw_label("Pre-integration UMAPs",
#     fontface = 'bold',x = 0,hjust = 0
#   ) + theme(plot.margin = margin(0, 0, 0, 7))
# plot_grid(ptitle, umaps, ncol = 1, rel_heights = c(0.1, 1))
# 
# dev.off()
# 
# # do I drop the low expressing cells??
# quantile(vp3.merged$nCount_SCT, seq(0,1,0.1))
# vp3.merged <- subset(vp3.merged, subset = nCount_SCT > 1000)
# vp3.merged@meta.data <- vp3.merged@meta.data[,!grepl("^[pANN|DF]", colnames(vp3.merged@meta.data))]

###
# Should we adjust for sequencing depths?? Yes
for(i in seq_along(vp3_dat)){
  DefaultAssay(vp3_dat[[i]]) <- "RNA"
}

# FOR SCT NORMS
min_median <- min(sapply(vp3_dat, function(x) median(colSums(x))))
vp3_dat <- lapply(vp3_dat, SCTransform, scale_factor = min_median, method = "glmGamPoi")
vp3.merged <- merge(vp3_dat[[1]], vp3_dat[2:length(vp3_dat)], add.cell.ids = vp3_names, merge.data=T)

# ...if you want to review merge
features <- SelectIntegrationFeatures(vp3_dat)
vp3.merged <- RunPCA(vp3.merged, features = features)
ElbowPlot(vp3.merged)
vp3.merged <- FindNeighbors(vp3.merged, dims = 1:15)
vp3.merged <- RunUMAP(vp3.merged, dims = 1:15)

pdf(
  here("code/01-process_DGS/outs/VP3_pre-integration_UMAPs_standardNorm.pdf"),
  height = 8,
  width = 11
)

p1 <- DimPlot(vp3.merged, group.by = "orig.ident")
p2 <- FeaturePlot(vp3.merged, features = "nCount_SCT")
p4 <- FeaturePlot(vp3.merged, features = "nFeature_RNA")
p5 <- FeaturePlot(vp3.merged, features = "percent.rps")
p6 <- FeaturePlot(vp3.merged, features = "percent.mt")
main_umaps <- plot_grid(p1,p2, ncol = 2)
bottom_umaps <- plot_grid(p4,p5,p6, nrow = 1)
umaps <- plot_grid(main_umaps, bottom_umaps, nrow = 2, rel_heights = c(1,0.5))
ptitle <- ggdraw() +
  draw_label("Pre-integration UMAPs",
             fontface = 'bold',x = 0,hjust = 0
  ) + theme(plot.margin = margin(0, 0, 0, 7))
plot_grid(ptitle, umaps, ncol = 1, rel_heights = c(0.1, 1))

dev.off()

quantile(vp3.merged$nCount_SCT, seq(0,1,0.1))
vp3.merged <- subset(vp3.merged, subset = nCount_SCT > 500)
vp3.merged@meta.data <- vp3.merged@meta.data[,!grepl("^[pANN|DF]", colnames(vp3.merged@meta.data))]


# FOR STANDARD NORMALIZATION
for(i in seq_along(vp3_dat)){
  DefaultAssay(seu) <- "RNA"
}
vp3.merged <- merge(vp3_dat[[1]], vp3_dat[2:length(vp3_dat)], add.cell.ids = vp3_names)

#############
# merging 5'
vp16_dat <- gex_obj_proc_db[grepl("5seq$", names(gex_obj_proc_db))]
vp16_names <- all_snames[grepl("5seq$", all_snames)]
# 
# vp16.merged <- merge(vp16_dat[[1]], vp16_dat[2:length(vp16_dat)], add.cell.ids = vp16_names, merge.data=T)
# vp16.merged <- RunPCA(vp16.merged, features = features)
# ElbowPlot(vp16.merged)
# vp16.merged <- FindNeighbors(vp16.merged, dims = 1:15)
# vp16.merged <- RunUMAP(vp16.merged, dims = 1:15)
# 
# pdf(
#   here("code/01-process_DGS/outs/VP16_pre-integration_UMAPs.pdf"),
#   height = 8,
#   width = 11
# )
# 
# p1 <- DimPlot(vp16.merged, group.by = "orig.ident")
# p2 <- FeaturePlot(vp16.merged, features = "nCount_SCT")
# p4 <- FeaturePlot(vp16.merged, features = "nFeature_RNA")
# p5 <- FeaturePlot(vp16.merged, features = "percent.rps")
# p6 <- FeaturePlot(vp16.merged, features = "percent.mt")
# main_umaps <- plot_grid(p1,p2, ncol = 2)
# bottom_umaps <- plot_grid(p4,p5,p6, nrow = 1)
# umaps <- plot_grid(main_umaps, bottom_umaps, nrow = 2, rel_heights = c(1,0.5))
# ptitle <- ggdraw() + 
#   draw_label("Pre-integration UMAPs",
#              fontface = 'bold',x = 0,hjust = 0
#   ) + theme(plot.margin = margin(0, 0, 0, 7))
# plot_grid(ptitle, umaps, ncol = 1, rel_heights = c(0.1, 1))
# 
# dev.off()
# 
# quantile(vp16.merged$nCount_SCT, seq(0,1,0.1))
# vp16.merged <- subset(vp16.merged, subset = nCount_SCT > 1000)
# vp16.merged@meta.data <- vp16.merged@meta.data[,!grepl("^[pANN|DF]", colnames(vp16.merged@meta.data))]

# rescaling for integration
for(i in seq_along(vp16_dat)){
  DefaultAssay(vp16_dat[[i]]) <- "RNA"
}
min_median <- min(sapply(vp16_dat, function(x) median(colSums(x))))
vp16_dat <- lapply(vp16_dat, SCTransform, scale_factor = min_median, method = "glmGamPoi")
vp16.merged <- merge(vp16_dat[[1]], vp16_dat[2:length(vp16_dat)], add.cell.ids = vp16_names, merge.data=T)

#...if you want to review merge
features <- SelectIntegrationFeatures(vp16_dat)
vp16.merged <- RunPCA(vp16.merged, features = features)
ElbowPlot(vp16.merged)
vp16.merged <- FindNeighbors(vp16.merged, dims = 1:15)
vp16.merged <- RunUMAP(vp16.merged, dims = 1:15)

pdf(
  here("code/01-process_DGS/outs/VP16_pre-integration_UMAPs_reSCT.pdf"),
  height = 8,
  width = 11
)

p1 <- DimPlot(vp16.merged, group.by = "orig.ident")
p2 <- FeaturePlot(vp16.merged, features = "nCount_SCT")
p4 <- FeaturePlot(vp16.merged, features = "nFeature_RNA")
p5 <- FeaturePlot(vp16.merged, features = "percent.rps")
p6 <- FeaturePlot(vp16.merged, features = "percent.mt")
main_umaps <- plot_grid(p1,p2, ncol = 2)
bottom_umaps <- plot_grid(p4,p5,p6, nrow = 1)
umaps <- plot_grid(main_umaps, bottom_umaps, nrow = 2, rel_heights = c(1,0.5))
ptitle <- ggdraw() +
  draw_label("Pre-integration UMAPs",
             fontface = 'bold',x = 0,hjust = 0
  ) + theme(plot.margin = margin(0, 0, 0, 7))
plot_grid(ptitle, umaps, ncol = 1, rel_heights = c(0.1, 1))

dev.off()

quantile(vp16.merged$nCount_SCT, seq(0,1,0.1))
vp16.merged <- subset(vp16.merged, subset = nCount_SCT > 600)
vp16.merged@meta.data <- vp16.merged@meta.data[,!grepl("^[pANN|DF]", colnames(vp16.merged@meta.data))]

# FOR STANDARD NORMALIZATION
for(i in seq_along(vp16_dat)){
  DefaultAssay(vp16_dat[[i]]) <- "RNA"
}
vp16.merged <- merge(vp16_dat[[1]], vp16_dat[2:length(vp16_dat)], add.cell.ids = vp16_names)


####################################
# Running RPCA integration    ######
####################################
batch_list <- list(vp3.merged, vp16.merged)

####
# THE SCT METHOD, BATCH SEPARATED
features <- SelectIntegrationFeatures(gex_obj_proc_db, nfeatures = 2000)
# doesn't seem to like "AP003327.2"
features <- features[features != "AP003327.2"]
batch_list <- PrepSCTIntegration(batch_list, anchor.features = features)
batch_list <- lapply(batch_list, RunPCA, features = features)
sample.anchors <- FindIntegrationAnchors(batch_list, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)

# saveRDS(sample.anchors, "processed/01-5_3-comp/integration_checkpoint.RDS")
gex.integrated <- IntegrateData(sample.anchors, normalization.method = "SCT", dims = 1:30)
saveRDS(gex.integrated, "processed/01-5_3-comp/samples_integrated_SCT_seqMerge.RDS")

# Run the standard workflow for visualization and clustering
gex.integrated <- RunPCA(gex.integrated, verbose = FALSE)
gex.integrated <- RunUMAP(gex.integrated, reduction = "pca", dims = 1:30)

####
# THE STANDARD METHOD, BATCH SEPARATED
batch_list <- lapply(X = batch_list, FUN = function(x) {
  DefaultAssay(x) <- 'RNA'
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = batch_list)
batch_list <- lapply(X = batch_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

sample.anchors <- FindIntegrationAnchors(object.list = batch_list, anchor.features = features, reduction = "rpca")
samples.combined <- IntegrateData(anchorset = sample.anchors)
samples.combined <- ScaleData(samples.combined, verbose = FALSE)
samples.combined <- RunPCA(samples.combined, npcs = 30, verbose = FALSE)
samples.combined <- RunUMAP(samples.combined, reduction = "pca", dims = 1:30)
saveRDS(samples.combined, "processed/01-5_3-comp/samples_integrated_RNA.RDS")

DimPlot(samples.combined, group.by = 'orig.ident')

####
# THE SCT METHOD, ALL GEMS
features <- SelectIntegrationFeatures(gex_obj_proc_db, nfeatures = 2000)
# doesn't seem to like "AP003327.2"
features <- features[features != "AP003327.2"]
gex_obj_proc_db <- PrepSCTIntegration(gex_obj_proc_db, anchor.features = features)
# gex_obj_proc_db <- lapply(gex_obj_proc_db, RunPCA, features = features)
sample.anchors <- FindIntegrationAnchors(gex_obj_proc_db, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)

# saveRDS(sample.anchors, "processed/01-5_3-comp/integration_checkpoint.RDS")
gex.integrated <- IntegrateData(sample.anchors, normalization.method = "SCT", dims = 1:15)
gex.integrated <- RunPCA(gex.integrated)
gex.integrated <- RunUMAP(gex.integrated, reduction = "pca", dims = 1:15)

saveRDS(gex.integrated, "processed/01-5_3-comp/samples_integrated_SCT_allGEM.RDS")


quantile(gex.integrated$nCount_SCT, seq(0,1,0.1))
gex.integrated <- subset(gex.integrated, subset = nCount_SCT > 1000)
gex.integrated@meta.data <- gex.integrated@meta.data[,!grepl("^[pANN|DF]", colnames(gex.integrated@meta.data))]
gex.integrated <- get_percent(gex.integrated)

gex.integrated <- RunPCA(gex.integrated, verbose = FALSE)
gex.integrated <- RunUMAP(gex.integrated, reduction = "pca", dims = 1:30, 
                          spread = 1, min.dist = 0.2, n.neighbors = 20)

####
# Some plotting

pdf(
  here("code/01-process_DGS/outs/post-integration_UMAPs_refined.pdf"),
  height = 8,
  width = 11
)

# DimPlot(gex.integrated, group.by = "orig.ident")
# dev.off()

p1 <- DimPlot(gex.integrated, group.by = "orig.ident")
p2 <- FeaturePlot(gex.integrated, features = "nCount_SCT")
p4 <- FeaturePlot(gex.integrated, features = "nFeature_RNA")
p5 <- FeaturePlot(gex.integrated, features = "percent.rps")
p6 <- FeaturePlot(gex.integrated, features = "percent.mt")
main_umaps <- plot_grid(p1,p2, ncol = 2)
bottom_umaps <- plot_grid(p4,p5,p6, nrow = 1)
umaps <- plot_grid(main_umaps, bottom_umaps, nrow = 2, rel_heights = c(1,0.5))
ptitle <- ggdraw() + 
  draw_label("Pre-integration UMAPs",
             fontface = 'bold',x = 0,hjust = 0
  ) + theme(plot.margin = margin(0, 0, 0, 7))
plot_grid(ptitle, umaps, ncol = 1, rel_heights = c(0.1, 1))

dev.off()

saveRDS(gex.integrated, "processed/01-5_3-comp/samples_integrated_SCT_allGEM_refined.RDS")


## UMAP ####
# 
# # looking at cluster composition
# pdata <- lapply(levels(dat.integrated$seurat_clusters), function(x) table(dat.integrated$orig.ident[dat.integrated$seurat_clusters==x]))
# pdata <- Reduce(cbind, pdata)
# colnames(pdata) <- levels(dat.integrated$seurat_clusters)
# pdata_sums <- apply(pdata, 2, sum)
# pdata_rel <- sweep(pdata, 2, pdata_sums, FUN = "/")
# pdata <- melt(pdata) 
# ggplot(pdata, aes(x=Var2, y=value, fill=Var1))+
#   geom_bar(position="stack", stat="identity")
# 
# 
# pdata_rel <- melt(pdata_rel)
# ggplot(pdata_rel, aes(x=Var2, y=value, fill=Var1))+
#   geom_bar(position="stack", stat="identity")
