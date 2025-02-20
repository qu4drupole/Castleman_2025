# investigating cluster composition and differences between experiment groups

library(Seurat)
library(ggplot2)
library(dplyr)
library(parallel)
library(here)
library(RColorBrewer)
library(tibble)
library(data.table)
library(future)
plan("multisession", workers = 8)
options(future.globals.maxSize= 20000*1024^2)
library(SingleCellExperiment)
library("miloR")

setwd(
  "/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration/"
)

seu <- readRDS("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration/processed/02-annotation/symphony_transfer_integrated_DGS_edited_onlyKFD-R.RDS")

###
# Add cell cycle score


###
# MiloR

# from checkpoint
dat_milo <- readRDS("code/02-clustering/milo_checkpoint.RDS")

# sce <- as.SingleCellExperiment(seu)
sce <- readRDS("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration/KFD_cellFilter_RNA_checkpoint.RDS") #
pca_dat <- seu@reductions$pca@cell.embeddings
reducedDim(sce, "pca") <- pca_dat[colnames(sce),]
umap_dat <- seu@reductions$umap@cell.embeddings
reducedDim(sce, "UMAP") <- umap_dat[colnames(sce),]
dat_milo <- Milo(sce)
# dat_milo@graph$graph <- buildFromAdjacency(seu.merged@graphs$SCT_nn, is.binary = T)
dat_milo <- buildGraph(dat_milo, k = 30, d = 30, reduced.dim = "pca")
dat_milo <- makeNhoods(dat_milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "pca")
# plotNhoodSizeHist(dat_milo)
dat_milo <- countCells(dat_milo, meta.data = as.data.frame(colData(dat_milo)), sample="edit.ident")

exp_design <- data.frame(colData(dat_milo))[,c("edit.ident", "disease", "seq_chem")]
exp_design$disease <- as.factor(exp_design$disease)
exp_design$seq_chem <- as.factor(exp_design$seq_chem)
exp_design <- distinct(exp_design)
rownames(exp_design) <- exp_design$edit.ident

dat_milo <- calcNhoodDistance(dat_milo, d=30, reduced.dim = "pca")
dat_milo <- buildNhoodGraph(dat_milo)
# contrast.1 <- c("groupDENV - groupDHF")
# da_results <- testNhoods(dat_milo, design = ~ 0 + disease + seq_chem, design.df = exp_design, model.contrasts = contrast.1,
#                          fdr.weighting="graph-overlap", norm.method="TMM")
da_results <- testNhoods(dat_milo, design = ~ seq_chem + disease, design.df = exp_design,
                         fdr.weighting="graph-overlap", norm.method="TMM", reduced.dim = "pca")
# da_results <- testNhoods(dat_milo, design = ~ disease, design.df = exp_design,
#                          fdr.weighting="k-distance", norm.method="RLE", reduced.dim = "pca")

nh_graph_pl <- plotNhoodGraphDA(dat_milo, da_results, layout="UMAP",alpha=0.3) 
nh_graph_pl

da_results <- annotateNhoods(dat_milo, da_results, coldata_col = "pred_anno")
ggplot(da_results, aes(pred_anno_fraction)) + geom_histogram(bins=50)
da_results$celltype <- ifelse(da_results$pred_anno_fraction < 0.6, "Mixed", da_results$pred_anno)
plotDAbeeswarm(da_results, group.by = "celltype", alpha = 0.5)

###
# saving
saveRDS(dat_milo, "code/02-clustering/milo_checkpoint.RDS")
