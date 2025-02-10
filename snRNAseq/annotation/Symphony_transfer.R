##########################################
# Running Harmony/Symphony label transfer
# 
# There are a few options:
#    1. Use the original labels from VP3, transfer to integrated data
#    2. Match the previous integrated labels to new integration, predict the remaining cells
#    3. Match the previous integrated labels to new integration, manually set the rest
##########################################

library(Seurat)
library(SeuratDisk)
library(symphony)
library(ggplot2)
library(tidyverse)
library(future)
plan("multisession", workers = 4)
options(future.globals.maxSize= 500000*1024^2)
library(irlba)
library(here)

####
# Are the high ribo cells garbage? 
#   (coming from gex.integrated in preprocess_generic.R)
FeatureScatter(gex.integrated, group.by = "orig.ident", feature1 = "nFeature_RNA", feature2 = "percent.rps")
quantile(gex.integrated$percent.rps, seq(0,1,0.1))


####
# STARTING THE DATA PREP
setwd("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration")
gex.integrated <- readRDS("processed/01-5_3-comp/samples_integrated_SCT_allGEM_refined.RDS")

# previous 3', 5' integration labels
# SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/test-cell-labels-4XGSVF2B.csv
anno.integrated <- read.csv("../../SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/test-cell-labels-4XGSVF2B.csv",
                            skip=2,
                            header = T)
# ...Then there's these annotations...
#    "processed/02-annotation/raw_seurat/cellxgene_metadata.csv"
old.meta <- read.csv("processed/02-annotation/raw_seurat/cellxgene_metadata.csv")

# VP3 annotations
anno.vp3 <- read.csv("../../SCTC-VP-3/cscb-lymph-vinodh/processed-data/05-annotation/vinodh_labels.csv",
                     header=T)

#############################
# MATCHING BARCODES         #
#############################
# since there are fewer cells in the annotations, coercing those barcode names to match integrated data

bc_match <- function(bc, ref){
  res <- rep(0,13)
  # bc_trim <- strsplit(bc,"_")[[1]][1]
  hit <- grep(bc, ref, value = T)
  # increment BC hit
  for(s in hit){
    res_i <- as.integer(strsplit(s,"_")[[1]][2])
    res[res_i] <- res[res_i] + 1
  }
  return(res)
}


# get some test barcodes from annotations
split.bc <- Reduce(rbind, strsplit(anno.integrated$index, "_"))
split.bc <- data.frame(split.bc)
colnames(split.bc) <- c("base_bc", "sample")
anno.integrated <- cbind(anno.integrated, split.bc)

test_bcs <- lapply(unique(anno.integrated$sample), function(x) sample(anno.integrated$base_bc[anno.integrated$sample==x], 100))

# match the sample IDs
bc_survey <- lapply(test_bcs, function(l) which.max(rowSums(sapply(l, bc_match, ref=Cells(gex.integrated)))))
anno.integrated$corrected_bc <- unlist(bc_survey)[as.integer(anno.integrated$sample)]
anno.integrated$new_bc <- paste(anno.integrated$base_bc, anno.integrated$corrected_bc, sep = "_")

# ...how many intersect? 59169, 85%
length(intersect(Cells(gex.integrated), anno.integrated$new_bc))
bc_intersect <- intersect(Cells(gex.integrated), anno.integrated$new_bc)

# there appear to be a few duplicated BCs...but, they're not in the intersection
dup_bc <- anno.integrated$new_bc[duplicated(anno.integrated$new_bc)]
View(anno.integrated[anno.integrated$new_bc %in% dup_bc,])
grep("AAGGAGCTCTGCGGCA", Cells(gex.integrated), value = T)
# it's a mix up between 9 and 10
table(anno.integrated$corrected_bc[anno.integrated$sample=="9"])
table(anno.integrated$corrected_bc[anno.integrated$sample=="10"])
table(anno.integrated$corrected_bc)
table(anno.integrated$sample)
# summary
# 1 -> 5
# 10 - no direct match
# 11 -> 10
# 12 -> 11
# 2 -> 6
# 3 -> 2
# 4 -> 3
# 5 -> 4
# 6 -> 1
# 7 -> 7
# 8 -> 8
# 9 - no direct match
#
# from new integration: 12 and 13 are entirely missing...but I expected one to be missing
# old 9 and 10 all go to new 9 now...
# ...is old 9 (1011 BCs) all duplicates of 10 (5861 BCs)?
sum(anno.integrated$base_bc[anno.integrated$sample=="9"] %in% anno.integrated$base_bc[anno.integrated$sample=="10"])
# no, only 8 BCs overlap, so it would seem one of the new sample labels is composed of 2 different "samples"????

new_sample_ids <- sapply(Cells(gex.integrated), function(x) strsplit(x,"_")[[1]][2])
table(new_sample_ids)

# so, what are these samples...
# old 9 = "KFD1"
# old 10 = "KFD3"
# new 9 = "KFD3"
# new 12 = "KFD2" (this should match old "KFD1_5seq"...)
#     but, old "KFD1_5seq" is old sample ID 9, which matched new "KFD3_5seq"
#     and, new "KFD2_5seq" does not match anything from old data
# new 13 = "MCD2" (not included in original data)
sapply(head(anno.integrated$base_bc[anno.integrated$sample=="10"]), function(s) grep(s, old.meta$orig.bc, value=T))

head(gex.integrated$orig.ident[new_sample_ids=="12"])

unique(old.meta$comp.ident)
unique(gex.integrated$orig.ident)
unique(gex.integrated$orig.ident)[!(unique(gex.integrated$orig.ident) %in% unique(old.meta$comp.ident))]
grep("KFD", unique(gex.integrated$orig.ident), value=T)
grep("KFD", unique(old.meta$comp.ident), value=T)

# !!!!!!!!!!!!!!!!!!!!!!!!! 
# Looking back through the notes, "KFD1_3seq" is actually "KFD2", which pairs with "KFD2_5seq"
# !!!!!!!!!!!!!!!!!!!!!!!!! 

bc_test <- sample(old.meta$orig.bc[old.meta$comp.ident=="KFD1_5seq"], 6)
bc_test <- sapply(bc_test, function(s) strsplit(s,"_")[[1]][2], simplify = T)
sapply(bc_test, function(s) grep(s, anno.integrated$index, value=T))

bc_test <- sample(Cells(gex.integrated)[gex.integrated$orig.ident=="KFD2_5seq"], 6)
bc_test <- sapply(bc_test, function(s) strsplit(s,"_")[[1]][1], simplify = T)
sapply(bc_test, function(s) grep(s, old.meta$orig.bc, value=T))

# ....I just don't know...I'll stick with the new names
anno.integrated <- anno.integrated[!duplicated(anno.integrated$new_bc),]
rownames(anno.integrated) <- anno.integrated$new_bc
# add missing cells to annotation
new.anno <- anno.integrated[intersect(anno.integrated$new_bc, Cells(gex.integrated)),]
missing_bc <- Cells(gex.integrated)[!(Cells(gex.integrated) %in% rownames(anno.integrated))]
filler <- data.frame(matrix(NA, nrow = length(missing_bc), ncol = ncol(anno.integrated)),
                     row.names = missing_bc)
names(filler) <- names(anno.integrated)
new.anno <- rbind(new.anno,filler)
new.anno <- new.anno[Cells(gex.integrated),]
# ...finally, how do they look...
gex.integrated$old_anno <- new.anno$Pillai.manual.annotation
DimPlot(gex.integrated, group.by = "old_anno", 
        label = T, label.size = 2, label.box = T, repel = T) + 
  theme(legend.position = "none")



#############################
# RUNNING HARMONY           #
#############################
gex_obj_proc_db <- readRDS("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration/processed/01-5_3-comp/processed_samples_list.RDS")
# Will use the SCT values
#     - should I standardize the scale factor?
#       -- I'll build a harmony reference using the batch separated SCT

sample.chem <- sapply(gex.integrated$orig.ident, function(s) strsplit(s,"_")[[1]][2])
gex.integrated$chem <- sample.chem
batch.list <- SplitObject(gex.integrated, split.by = "chem")

for(i in seq_along(batch.list)){
  DefaultAssay(batch.list[[i]]) <- "RNA"
}

# FOR SCT NORMS
# VP3...
vp3_dat <- SplitObject(batch.list$`3seq`, split.by = "orig.ident")
min_median <- min(sapply(vp3_dat, function(x) median(colSums(x))))
vp3_dat <- lapply(vp3_dat, SCTransform, scale_factor = min_median, method = "glmGamPoi")
vp3.merged <- merge(vp3_dat[[1]], vp3_dat[2:length(vp3_dat)], add.cell.ids = names(vp3_dat), merge.data=T)

# VP16...
vp16_dat <- SplitObject(batch.list$`5seq`, split.by = "orig.ident")
min_median <- min(sapply(vp16_dat, function(x) median(colSums(x))))
vp16_dat <- lapply(vp16_dat, SCTransform, scale_factor = min_median, method = "glmGamPoi")
vp16.merged <- merge(vp16_dat[[1]], vp16_dat[2:length(vp16_dat)], add.cell.ids = names(vp16_dat), merge.data=T)

# Now merge these two together...
gex.merged <- merge(vp3.merged, vp16.merged, merge.data = T)
saveRDS(gex.merged, "processed/01-5_3-comp/all_samples_seqScaledSCT_merged_DGS.RDS")
# and then picking up from, http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/docs/SeuratV3.html
#   or, this: https://github.com/immunogenomics/symphony/blob/main/vignettes/pbmcs_tutorial.ipynb/
#     but I think the Harmony instructions with RunHarmony is just a wrapper for HarmonyMatrix called in the Symphony vignette

features <- SelectIntegrationFeatures(gex_obj_proc_db, nfeatures = 3000)
# ...subset merged data to only include cells with known labels
keep_bc <- Cells(gex.integrated)[!is.na(gex.integrated$old_anno)]

base_bc <- sapply(Cells(gex.merged), function(s) strsplit(s,"_")[[1]][3])
base_id <- sapply(Cells(gex.merged), function(s) strsplit(s,"_")[[1]][4])
gex.merged$orig_bc <- paste(base_bc, base_id, sep="_")
gex.merged.ref <- subset(gex.merged, cells = Cells(gex.merged)[gex.merged$orig_bc %in% keep_bc])
gex.merged.ref <- RunPCA(gex.merged.ref, features = features, npcs = 20, verbose = FALSE)
# gex.merged$chem <- factor(gex.merged$chem)
# gex.merged.ref <- harmony::RunHarmony(gex.merged.ref, "chem", plot_convergence=T)

harmony_meta <- gex.merged.ref@meta.data
# harmony_ref <- gex.merged.ref@reductions$harmony
ex <- GetAssayData(gex.merged.ref)
ex <- ex[features,]
vargenes_means_sds = tibble(symbol = features, mean = Matrix::rowMeans(ex))
vargenes_means_sds$stddev <- rowSDs(ex, row_means = vargenes_means_sds$mean)
ref_exp_scaled = scaleDataWithStats(ex, vargenes_means_sds$mean, vargenes_means_sds$stddev, 1)
set.seed(0)
s = irlba(ref_exp_scaled, nv = 20)
Z_pca_ref = diag(s$d) %*% t(s$v)
loadings = s$u
harmony_ref <- harmony::HarmonyMatrix(t(Z_pca_ref), 
                                      gex.merged.ref$chem, 
                                      do_pca = F, 
                                      return_object = T)


reference = buildReferenceFromHarmonyObj(
  harmony_ref,            # output object from HarmonyMatrix()
  harmony_meta,           # reference cell metadata
  vargenes_means_sds,     # gene names, means, and std devs for scaling
  loadings,               # genes x PCs matrix
  verbose = TRUE,         # verbose output
  do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
  save_uwot_path = './reference_prev_vp3_vp16_integration')

saveRDS(reference, 'processed/02-annotation/symphony_vp3_vp16.rds')
reference <- readRDS('processed/02-annotation/symphony_vp3_vp16.rds')
# Now to map query, which will be all the data
query_exp <- GetAssayData(gex.merged, slot = 'data')
# query_exp <- query_exp[features,]
query_metadata <- gex.merged@meta.data
query = mapQuery(query_exp,             # query gene expression (genes x cells)
                 query_metadata,        # query metadata (cells x attributes)
                 reference,             # Symphony reference object
                 do_normalize = FALSE,  # perform log(CP10k+1) normalization on query
                 do_umap = T)        # project query cells into reference UMAP

query = knnPredict(query, reference, reference$meta_data$old_anno, k = 10)



pdf(
  here("code/02-annotation_DGS/outs/Symphony_integrations_seqSCTscaled.pdf"),
  height = 8,
  width = 11
)

pdata <- data.frame(reference$umap)
ggplot(pdata, aes(x=embedding.UMAP1, y=embedding.UMAP2))+
  geom_point(aes(colour=reference$meta_data$old_anno), size=0.1)+
  ggtitle("reference, old labels")

ggplot(pdata, aes(x=embedding.UMAP1, y=embedding.UMAP2))+
  geom_point(aes(colour=reference$meta_data$orig.ident), size=0.1)+
  ggtitle("reference, sample ID")

ggplot(pdata, aes(x=embedding.UMAP1, y=embedding.UMAP2))+
  geom_point(aes(colour=reference$meta_data$chem), size=0.1)+
  ggtitle("reference, seq. chemistry")

pdata <- data.frame(query$umap)
ggplot(pdata, aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(colour=query$meta_data$cell_type_pred_knn), size=0.1)+
  ggtitle("query, predicted labels")

ggplot(pdata, aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(colour=query$meta_data$orig.ident), size=0.1)+
  ggtitle("query, sample ID")

ggplot(pdata, aes(x=UMAP1, y=UMAP2))+
  geom_point(aes(colour=query$meta_data$chem), size=0.1)+
  ggtitle("query, seq. chemistry")

dev.off()


query_meta <- query$meta_data
rownames(query_meta) <- query_meta$orig_bc
query_meta <- query_meta[Cells(gex.integrated),]

gex.integrated$pred_anno <- query_meta$cell_type_pred_knn
gex.integrated$pred_anno_prob <- query_meta$cell_type_pred_knn_prob
# FeaturePlot(gex.integrated, features = "nCount_SCT")

pdf(
  here("code/02-annotation_DGS/outs/cell_type_predictions.pdf"),
  height = 8,
  width = 11
)

DimPlot(gex.integrated, group.by = "old_anno", 
        label = T, label.size = 2, label.box = T, repel = T) + 
  theme(legend.position = "none") +
  ggtitle("Original matched labels")

DimPlot(gex.integrated, group.by = "pred_anno", 
        label = T, label.size = 2, label.box = T, repel = T) + 
  theme(legend.position = "none") +
  ggtitle("Predicted, all cells")

dev.off()

pdf(
  here("code/02-annotation_DGS/outs/random_UMAPs.pdf"),
  height = 8,
  width = 11
)

DimPlot(seu, group.by = "pred_anno", 
        label = T, label.size = 2, label.box = T, repel = T) + 
  theme(legend.position = "none") +
  ggtitle("Predicited cell type")

DimPlot(seu, group.by = "edit.ident")

DimPlot(seu, group.by = "chem")

dev.off()

saveRDS(gex.integrated, "processed/02-annotation/symphony_transfer_integrated_DGS.RDS")


##################
# EXPORTING      #
##################
DefaultAssay(gex.integrated) <- "SCT"

str(gex.integrated@meta.data)
gex.integrated$pred_anno <- as.character(gex.integrated$pred_anno)

SaveH5Seurat(gex.integrated, filename = "processed/02-annotation/symphony_transfer_integrated_DGS.H5Seurat", overwrite = T)
Convert("processed/02-annotation/symphony_transfer_integrated_DGS.H5Seurat", dest = "h5ad", overwrite = T)

