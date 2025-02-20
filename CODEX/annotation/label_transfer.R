# Using Symphony for automated label transfer
# reference: https://github.com/immunogenomics/symphony/blob/main/vignettes/pbmcs_tutorial.ipynb/
# Sys.setenv(RETICULATE_PYTHON = "/home/smithd30/software/anaconda3/envs/scanpy-basic/bin/python")
Sys.setenv(RETICULATE_PYTHON = "/mnt/isilon/cscb/software/anaconda3/envs/scanpy-basic/bin/python")
library(reticulate)
# Sys.setenv(BASILISK_EXTERNAL_CONDA = '/mnt/isilon/cscb/software/anaconda3')
# 
Sys.setenv(BASILISK_EXTERNAL_DIR  = '/mnt/isilon/cscb/software/R-basilisk')
library(basilisk)
library(zellkonverter)
# library(SeuratDisk)
library(Seurat)
library(harmony)
library(symphony)
library(irlba)
library(ggplot2)
library(dplyr)
library(parallel)
library(pbmcapply)
# sc <- import("scanpy")
setwd('/mnt/isilon/cscb/codex/pillaiv/SCTC-VP-15/')
# r1 <- readH5AD('outs/R1_region1_CellSeg_cluster-checkpoint-mod1.h5ad', reader='python', raw=T)
r1 <- readH5AD('outs/cleaned_ref_data/R1_v2.h5ad', reader='python', raw=T)
r1 <- as.Seurat(r1, counts = 'X', data = 'X', project = "r1")
# k2 <- readH5AD('outs/K2_region1_CellSeg_cluster-checkpoint.h5ad', reader='python', raw=T)
k2 <- readH5AD('outs/cleaned_ref_data/K2_v2.h5ad', reader='python', raw=T)
k2 <- as.Seurat(k2, counts = 'X', data = 'X', project = 'k2')
# mcd1 <- readH5AD('outs/...', reader='python', raw=T)
mcd1 <- readH5AD('outs/cleaned_ref_data/MCD1_v2.h5ad', reader='python', raw=T)
mcd1 <- as.Seurat(mcd1, counts = 'X', data = 'X', project = 'mcd1')
# hvcd1 <- readH5AD('outs/HVCD1_region2_CellSeg_cluster-checkpoint.h5ad', reader='python', raw=T)
hvcd1 <- readH5AD('outs/cleaned_ref_data/HVCD1_v2.h5ad', reader='python', raw=T)
hvcd1 <- as.Seurat(hvcd1, counts = 'X', data = 'X', project = 'hvcd1')

# Trim to common markers
markers <- Reduce(intersect, list(rownames(hvcd1), rownames(r1), rownames(mcd1), rownames(k2)))
r1 <- r1[markers,]
k2 <- k2[markers,]
mcd1 <- mcd1[markers,]
hvcd1 <- hvcd1[markers,]

# Check cell type name consistency
Reduce(union, list(unique(hvcd1$new_cell_type), 
                   unique(r1$new_cell_type), 
                   unique(mcd1$new_cell_type), 
                   unique(k2$new_cell_type)))

# Need to change Macrophage_MPO_pos, which should only be K2?
# k2$cell_type <- recode(k2$cell_type, 'Macrophage_MPO_pos' = 'Macrophage_MPO+')

# raw.ref <- merge(r1, k2, add.cell.ids = c("r1","k2"))
raw.ref <- merge(r1, k2)
raw.ref <- merge(raw.ref, mcd1)
raw.ref <- merge(raw.ref, hvcd1)
raw.ref$ID <- Idents(raw.ref)

# # scale.factor = 1 bc it should already be normalized (but not log-transformed?)
# VlnPlot(k2, features = c("nFeature_originalexp", "nCount_originalexp"), ncol = 2)
# VlnPlot(raw.ref, features = c("nFeature_originalexp", "nCount_originalexp"), ncol = 2)
# #  looks like it's the scale data...
# plot(density(k2@assays$originalexp@data[2,]))
# min(k2@assays$originalexp@data[5,])
# min(k2@assays$originalexp@counts[5,])

# raw.ref@assays$originalexp@scale.data <- raw.ref@assays$originalexp@data
# raw.ref <- RunPCA(raw.ref, features = rownames(raw.ref), assay = "originalexp")


vargenes_means_sds = data.frame(symbol = rownames(raw.ref), mean = Matrix::rowMeans(raw.ref@assays$raw@data))
vargenes_means_sds$stddev = singlecellmethods::rowSDs(raw.ref@assays$raw@data, vargenes_means_sds$mean)


set.seed(59)
s = irlba(raw.ref@assays$originalexp@data, nv = 15)
Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
loadings = s$u

# V <- cell_lines$scaled_pcs
# meta_data <- cell_lines$meta_data

ref_harmObj = HarmonyMatrix(
  data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
  meta_data = raw.ref@meta.data, ## dataframe with cell labels
  theta = c(2),             ## cluster diversity enforcement
  vars_use = c('ID'),    ## variable to integrate out
  nclust = 100,             ## number of clusters in Harmony model
  max.iter.harmony = 20,
  return_object = TRUE,     ## return the full Harmony model object
  do_pca = FALSE            ## don't recompute PCs
)


# Harmony part...
reference = buildReferenceFromHarmonyObj(
  ref_harmObj,            # output object from HarmonyMatrix()
  raw.ref@meta.data,           # reference cell metadata
  vargenes_means_sds,     # gene names, means, and std devs for scaling
  loadings,               # genes x PCs matrix
  verbose = TRUE,         # verbose output
  do_umap = TRUE,
  save_uwot_path = './model_generic_cell_types')

saveRDS(reference, 'outs/R1_K2_MCD1_HVCD1_reference_generic_cell_type.rds')

umap_labels = cbind(raw.ref@meta.data, reference$umap$embedding)
pdata <- umap_labels[sample(1:nrow(umap_labels), 50000), ]

ggplot(pdata, aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(color = ID), size=0.3)

ggplot(pdata, aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(color = cell_type))

# TRYING TO RUN A QUERY

#   if not already loaded...
reference <- readRDS('outs/R1_K2_MCD1_HVCD1_reference.rds')

query <- readH5AD('outs/simple_h5ad/R1_reg2_statistics_growth5_comp.h5ad', reader='python', raw=T)
query <- as.Seurat(query, counts = 'X', data = 'X', project = "r1_s2")

query = mapQuery(query@assays$raw@data,             # query gene expression (genes x cells)
                 query@meta.data,        # query metadata (cells x attributes)
                 reference,             # Symphony reference object
                 do_normalize = FALSE,  # perform log(CP10k+1) normalization on query
                 do_umap = TRUE)  

query = knnPredict(query, reference, reference$meta_data$cell_type, k = 5)

# Visualize!
# Sync the column names for both data frames
reference$meta_data$ref_query = 'reference'
query$meta_data$ref_query = 'query'
query$meta_data$ID = 'query'

# Add the UMAP coordinates to the metadata
keep_cols <- intersect(colnames(query$meta_data), colnames(reference$meta_data))
meta_data_combined = rbind(query$meta_data[, keep_cols], reference$meta_data[,keep_cols])
umap_combined = rbind(query$umap, reference$umap$embedding)
umap_combined_labels = cbind(meta_data_combined, umap_combined)

pdata <- umap_combined_labels[sample(1:nrow(umap_combined_labels), 50000), ]
ggplot(pdata, aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(color = ref_query), size=0.3)
ggplot(pdata, aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(color = ID), size=0.3)

# EXPORTING
write.csv(query$meta_data, "outs/simple_h5ad/R1_reg2_pred.csv")


#~~~~~~~~~~~~~~~~~~~~~
# Iterating predictions over all files
reference <- readRDS('outs/R1_K2_MCD1_HVCD1_reference.rds')
in_path <- 'outs/simple_h5ad'
s_files <- list.files(in_path, pattern = '*.h5ad')
s_files
# for(f in s_files[26:29]){ #job timed out...
#   # browser()
#   print(f)
#   in_file = file.path(in_path, f)
#   pname = paste(unlist(strsplit(f, '_'))[1:2], collapse = "")
#   query <- readH5AD(in_file, reader='python', raw=T)
#   query <- as.Seurat(query, counts = 'X', data = 'X', project = pname)
#   
#   query = mapQuery(query@assays$raw@data,             # query gene expression (genes x cells)
#                    query@meta.data,        # query metadata (cells x attributes)
#                    reference,             # Symphony reference object
#                    do_normalize = FALSE,  # perform log(CP10k+1) normalization on query
#                    do_umap = TRUE)  
#   
#   query = knnPredict(query, reference, reference$meta_data$cell_type, k = 5)
#   
#   out_name <- paste(unlist(strsplit(f, '\\.'))[1], 'pred.csv', sep = '_')
#   out_path <- file.path(in_path, out_name)
#   write.csv(query$meta_data, out_path)
# }

# parallelizing it
predict_cell <- function(f, ref = reference){
  in_file = file.path(in_path, f)
  pname = paste(unlist(strsplit(f, '_'))[1:2], collapse = "")
  query <- readH5AD(in_file, reader='python', raw=T)
  query <- as.Seurat(query, counts = 'X', data = 'X', project = pname)
  
  query = mapQuery(query@assays$raw@data,             # query gene expression (genes x cells)
                   query@meta.data,        # query metadata (cells x attributes)
                   ref,             # Symphony reference object
                   do_normalize = FALSE,  # perform log(CP10k+1) normalization on query
                   do_umap = TRUE)  
  
  query = knnPredict(query, ref, ref$meta_data$new_cell_type, k = 5)
  
  out_name <- paste(unlist(strsplit(f, '\\.'))[1], 'pred.csv', sep = '_')
  out_path <- file.path(in_path, out_name)
  write.csv(query$meta_data, out_path)
}

pbmclapply(s_files, predict_cell, mc.cores = 8)
