library(Seurat)
library(tidyverse)
library(cowplot)
library(dplyr)
library(ggplot2)
# library(Platypus) # not available on CRAN for R 4.4
library(RColorBrewer)
library(here)
library(ggsci)
library(viridis)
library(future)
plan(multisession, workers = 4)
options(future.globals.maxSize= 20000*1024^2)
setwd("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration")
here::i_am("code/making_figures.R")
seu <- readRDS("processed/02-annotation/symphony_transfer_integrated_DGS_edited_CD-R.RDS")

sort(unique(seu$pred_anno))
ct_colors <- c('#33cc33', '#FFFF66', '#CCCC66', '#CC3366', '#333399', '#663399', '#9999cc', '#006633', '#993366',
               '#9966CC', '#FFCC33', '#ff6633', '#996633', '#009933', '#FF6666', '#ff6699', '#993399', '#FF9933',
               '#cc6699', '#339999', '#CC99CC', '#FF9966')

sort(unique(seu$sample_name))
# donor_color <- c('#9999cc', '#666699', '#663399', '#ffcc66', '#ff9933', '#ff6633', '#cc3333', '#66cc66', '#339933')
donor_color <- c('#C197C6', '#854C74', '#963C6D', '#FBB14B', '#F58535', '#F15C34', '#D42027', '#60A9C3', '#364B8B')
names(donor_color) <- sort(unique(seu$sample_name))

sort(unique(seu$seq_chem))
seq_color <- c('#0099ff', '#ff0099')

##########################################
# Figure 1/2 -- interfollicular distance
##########################################
setwd("/mnt/isilon/cscb/Projects/codex/pillaiv/SCTC-VP-15/")
cc_dist <- read.csv("code/CD_analysis/figures/Figure2/interfollicular/distance_distributions/all_interfoll/dist_200/stats.csv")
source_c <- t(sapply(cc_dist$X, function(s) strsplit(s,"--")[[1]])) %>% data.frame
source_c <- apply(source_c,2, trimws) %>% data.frame
cc_dist$source <- source_c$X1
cc_dist$neighbor <- source_c$X2

hdata <- cc_dist %>% 
  mutate(diff_mcd = avg_MCD_distance - avg_R_distance, diff_ucd = avg_HVCD_distance - avg_R_distance) %>%
  dplyr::select(source, neighbor, diff_mcd, diff_ucd) %>%
  melt(id.vars=c("source","neighbor")) %>%
  mutate(scaled_dist = scale(value)[,1])

ggplot(hdata, aes(x=source, y=neighbor, fill=scaled_dist)) +
  geom_tile()+
  # scale_fill_viridis(option="magma", )+
  scale_fill_gradient2(low = "purple", mid="white", high = "orange")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~variable, ncol = 2)

ggplot(hdata, aes(x=source, y=neighbor, fill=value)) +
  geom_tile()+
  # scale_fill_viridis(option="magma", )+
  scale_fill_gradient2(low = "purple", mid="white", high = "orange")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~variable, ncol = 2)

pdata <- hdata[hdata$variable == "diff_mcd",] %>% 
  pivot_wider(id_cols = source, names_from = neighbor, values_from = value, values_fill = NA) %>%
  remove_rownames %>% column_to_rownames(var="source") 

pdf(
  here("figures/Figure2/heatmap_interfollicular_distance_targetScaled.pdf"),
  height = 4,
  width = 6
)

pdata <- hdata[hdata$variable == "diff_mcd",] %>% 
  pivot_wider(id_cols = source, names_from = neighbor, values_from = value, values_fill = NA) %>%
  remove_rownames %>% column_to_rownames(var="source") 

print(pheatmap(t(as.matrix(pdata)), scale = "row", cluster_cols = F, 
               main = "MCD"))

pdata <- hdata[hdata$variable == "diff_ucd",] %>% 
  pivot_wider(id_cols = source, names_from = neighbor, values_from = value, values_fill = NA) %>%
  remove_rownames %>% column_to_rownames(var="source") 

print(pheatmap(t(as.matrix(pdata)), scale = "row", cluster_cols = F, 
               main = "UCD"))

dev.off()

##########################################
# Figure 3 -- snRNA-seq
##########################################
#####
# UMAPs
seu$disease[seu$disease=="HVCD"] <- "UCD"

DimPlot(seu, group.by="pred_anno", raster = T, cols = ct_colors, pt.size = 2)

DimPlot(seu, group.by="pred_anno_2", raster = T, cols = ct_colors)

DimPlot(seu, group.by="pred_anno", split.by="sample_name", 
        raster = T, cols = ct_colors, ncol = 4, pt.size = 3) + theme(legend.position = "none")

DimPlot(seu, group.by="pred_anno", split.by="disease", 
        raster = T, cols = ct_colors) + theme(legend.position = "none")

DimPlot(seu, group.by="sample_name", raster = T, cols = donor_color)

DimPlot(seu, group.by="seq_chem", raster = T, cols = seq_color)


DimPlot(seu, group.by="pred_anno", split.by = "disease", 
        cols = ct_colors, pt.size = 0.2, ncol = 1) + theme(legend.position = "none")

DimPlot(seu, group.by="pred_anno", 
        cols = ct_colors, pt.size = 0.2,
        label = T, label.size = 5, label.box = T, repel = T) + theme(legend.position = "none")


#####
# cell type bar charts

cell_comp <- seu@meta.data %>% dplyr::group_by(sample_name, pred_anno) %>%
  dplyr::summarise("cell_count"=length(pred_anno))

cell_comp_sums <-table(seu$pred_anno)
cell_comp$rel_count <- apply(cell_comp, 1, function(x) as.numeric(x[3]) / cell_comp_sums[x[2]])
sample_comp_sums <-table(seu$sample_name)
cell_comp$rel_count_sample <- apply(cell_comp, 1, function(x) as.numeric(x[3]) / sample_comp_sums[x[1]])

cell_comp$pred_anno <- factor(cell_comp$pred_anno, 
                              levels = c("Naive B cells", "Activated and memory B cells", "Germinal center B cells", "Plasma cells",
                                         "Proliferating plasma cells", "Naive CD4 and Naive CD8 T cells", "Cytotoxic CD8 T cells", "Memory or effector T cells",
                                         "Regulatory T cells", "T Follicular Helper cells", "Monocytes", "Classic Dendritic cells", "Activated and migratory cDC",
                                         "Macrophages", "Plasmacytoid dendritic cells", "Granulocytes", "NK cells", "Stromal cells", "Follicular dendritic cells",
                                         "Lymphatics", "Fibroblasts", "Endothelial cells"),
                              ordered = T)

ggplot(cell_comp, aes(fill=sample_name, y=rel_count, x=pred_anno)) + 
  geom_bar(position="stack", stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5,
                                   hjust=1,
  ))

# 
cell_comp <- seu@meta.data %>% dplyr::group_by(pred_anno) %>%
  dplyr::summarise("cell_count"=length(pred_anno))

cell_comp$pred_anno <- factor(cell_comp$pred_anno, 
                              levels = c("Naive B cells", "Activated and memory B cells", "Germinal center B cells", "Plasma cells",
                                         "Proliferating plasma cells", "Naive CD4 and Naive CD8 T cells", "Cytotoxic CD8 T cells", "Memory or effector T cells",
                                         "Regulatory T cells", "T Follicular Helper cells", "Monocytes", "Classic Dendritic cells", "Activated and migratory cDC",
                                         "Macrophages", "Plasmacytoid dendritic cells", "Granulocytes", "NK cells", "Stromal cells", "Follicular dendritic cells",
                                         "Lymphatics", "Fibroblasts", "Endothelial cells"),
                              ordered = T)

cell_comp$pct.comp <- cell_comp$cell_count/sum(cell_comp$cell_count)
ggplot(cell_comp, aes(y=pct.comp, x=pred_anno)) + 
  geom_bar(position="stack", stat="identity") +
  theme_minimal()


###
# Immune profiling stuff
all_vgm <- readRDS("code/08-Immune-profiling/vdj_dat_checkpoint.RDS")
b_cell_meta <- all_vgm$b
seu$raw_bc <- sapply(Cells(seu), function(s) strsplit(s,"-")[[1]][1])
seu$bc_convert <- paste(seu$raw_bc, seu$edit.ident, sep = "_")
b_cell_meta$raw_bc <- sapply(b_cell_meta$barcode, function(s) strsplit(s,"_")[[1]][2])
b_cell_meta$bc_convert <- paste(b_cell_meta$raw_bc, b_cell_meta$group_id, "5seq", sep = "_")
b_cell_meta <- b_cell_meta[b_cell_meta$bc_convert %in% seu$bc_convert, ]
rownames(b_cell_meta) <- b_cell_meta$bc_convert

seu$BCRlabel <- "unknown"
seu$BCRlabel_simple <- seu$BCRlabel
b_cell_meta$VDJ_cgene[b_cell_meta$VDJ_cgene==""] <- "unknown"
seu_bc_inorder <- seu$bc_convert[seu$bc_convert %in% b_cell_meta$bc_convert]
b_cell_meta <- b_cell_meta[seu_bc_inorder, ]
b_cell_meta$VDJ_cgene_simple <- sapply(b_cell_meta$VDJ_cgene, function(s) gsub("\\d$","",s))
seu$BCRlabel[seu$bc_convert %in% b_cell_meta$bc_convert] <- b_cell_meta$VDJ_cgene
seu$BCRlabel_simple[seu$bc_convert %in% b_cell_meta$bc_convert] <- b_cell_meta$VDJ_cgene_simple

# filter out non-b cell bcr labels
seu$BCRlabel[!(seu$pred_anno_2 %in% c("Plasma cells","Naive B cells","Activated and memory B cells",
                                      "Germinal center B cells"))] <- "unknown"
b_cell_highlight <- unique(seu$BCRlabel)
b_cell_highlight <- b_cell_highlight[-1] # drop "unknown"
b_cell_highlight_ <- sapply(b_cell_highlight, function(x) Cells(seu)[seu$BCRlabel==x])
DimPlot(seu, group.by = "BCRlabel",
        cells.highlight = b_cell_highlight_, sizes.highlight = 0.4, 
        cols.highlight = brewer.pal(length(b_cell_highlight_), "Set1"))

###
# The SHM plot
VDJ_mixcr_out_all <- readRDS("code/08-Immune-profiling/allSamples_SHM_checkpoint.RDS")
sample_names_b_cd <- c("HVCD2", "HVCD3",  "MCD4", "MCD3", "R2")
VDJ_mixcr_out_all <- VDJ_mixcr_out_all[VDJ_mixcr_out_all$sample_id %in% sample_names_b_cd,]
plist <- VDJ_plot_SHM(VDJ_mixcr_out_all)
new_donor_color <- donor_color[sort(unique(VDJ_mixcr_out_all$sample_id))]

plist[[1]] +
  scale_color_manual(values=unname(new_donor_color)) +
  scale_fill_manual("blue") +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5,
                                   hjust=1,
  ))

##########################################
# Figure 4 -- Cell type DEG, enrichments
##########################################

interesting_cells <- c("Stromal cells","Lymphatics","Endothelial cells","Fibroblasts",
                       "Follicular dendritic cells","Activated and migratory cDC","Plasmacytoid dendritic cells",
                       "Activated and memory B cells")

###
# VEGF violin plots
Idents(seu) <- "pred_anno_2"

pdf(
  here("figures/Figure4/Violin_VEGFA_byDisease_RNA.pdf"),
  height = 12,
  width = 12
)

plist <- list()
for(ct in unique(Idents(seu))){
  plist[[ct]] <- VlnPlot(seu, features = "VEGFA", idents = ct,
                assay = "RNA", layer="data", 
                group.by = "disease", pt.size = 0) +
          ggtitle(ct)+
    ylim(0,4)
}

plot_grid(plotlist = plist, ncol = 4)

dev.off()

###
# IL-6 module violins
DefaultAssay(seu) <- "RNA"
Idents(seu) <- "pred_anno_2"
il6_mod <- list(c("IL6","LIFR","OSMR","IL6ST"),
                c("IL6","LIFR","OSMR","IL6R","IL6ST"))
# il6_mod[[1]] %in% rownames(seu)

seu <- AddModuleScore(seu, features = il6_mod, assay = "SCT", name = "IL6_mod_sct")
ct<-"Follicular dendritic cells"
VlnPlot(seu, features = "IL6_mod_sct1", idents = ct, group.by = "edit.ident")

pdf(
  here("figures/Figure4/Violin_IL6module_bySample_RNA.pdf"),
  height = 24,
  width = 12
)

plist <- list()
for(ct in unique(Idents(seu))){
  plist[[ct]] <- VlnPlot(seu, features = "IL6_mod_rna2", idents = ct,
                         group.by = "edit.ident", pt.size = 0) +
    ggtitle(ct)+
    ylim(0,2)
}

plot_grid(plotlist = plist, ncol = 4)

dev.off()

###
# VEGF AUC
#   see "AUCell_analysis_CD.R"
pdata <- t(cells_AUC@assays@data$AUC)
pdata <- as.data.frame(pdata)
pdata$pred_anno_2 <- seu_cd$pred_anno_2
pdata$disease <- seu_cd$disease
pdata <- reshape2::melt(pdata, id.vars = c("pred_anno_2","disease"))
pdata_path <- pdata[pdata$variable=="VEGF ligand-receptor interactions",]
pdata_path <- pdata_path[pdata_path$pred_anno_2 %in% interesting_cells,]

ggplot(pdata_path, aes(x=disease, y=value, fill=disease))+
          # geom_violin(linewidth=0.25, position=position_dodge(width = 0.75),
          #             width=4, alpha=0.8, trim=T, draw_quantiles = c(0.5),
          #             outlier.colour="transparent", scale="count",
          #             adjust=200, bw=0.00005)+
          # geom_violin(linewidth=0.25, position=position_dodge(width = 0.75), trim=T, draw_quantiles = c(0.5),
          #             outlier.colour="transparent",
          #             adjust=200, bw=0.00005)+
          geom_violin(scale="width", bw=0.03, position=position_dodge(width = 1)) +
          # geom_boxplot()+
          ggtitle(ct)+
          theme_minimal()+
  facet_wrap(~pred_anno_2, ncol = 3)


###
# plasma cell abundance vs. module scores
plasma_comp <- cell_comp[cell_comp$pred_anno == "Plasma cells",]
plasma_comp$disease <- sapply(plasma_comp$sample_name, function(s) str_match(s,"\\D*")[1])
plasma_comp$il6.fdc <- sapply(plasma_comp$sample_name, 
                              function(x) sum(seu$IL6_mod_sct2[(seu$pred_anno_2 == "Follicular dendritic cells") & (seu$sample_name==x)]))
plasma_comp$il6.fibro <- sapply(plasma_comp$sample_name, 
                              function(x) sum(seu$IL6_mod_sct2[(seu$pred_anno_2 == "Fibroblasts") & (seu$sample_name==x)]))
plasma_comp$il6.combined <- plasma_comp$il6.fdc + plasma_comp$il6.fibro
ggplot(plasma_comp, aes(x=il6.combined, y=log(cell_count), color=disease))+
  geom_point()+
  theme_bw()

pdata <- melt(plasma_comp, id.vars = c("sample_name", "cell_count", "rel_count_sample"), measure.vars = c("il6.fdc","il6.fibro"))

ggplot(pdata, aes(x=value, y=log(cell_count), color=disease))+
  geom_point()+
  theme_bw()+
  facet_wrap(~variable, scales = "free_x")

###
# t test

x1 <- seu$IL6_mod_rna2[(seu$pred_anno_2 == "Follicular dendritic cells") & (seu$disease == "R")]
x2 <- seu$IL6_mod_rna2[(seu$pred_anno_2 == "Follicular dendritic cells") & (seu$disease == "MCD")]

t.test(x1, x2)

###
# enrichment bar plots
#   ...taken from middle of code chunk...

plts <- list()
for(ct in ctoi){
  hdata_ <- read.csv(paste0("figures/rnaseq_enrich/ct_enrichments/",ct,"_aggRes.csv"), row.names = 1)
  hdata_ <- hdata_[hdata_$select.term == "y",]
  pdata <- data.frame(p_value = c(hdata_$p_value.UCD, hdata_$p_value.MCD), 
                      f_score = c(hdata_$f_score.UCD, hdata_$f_score.MCD),
                      term = rep(rownames(hdata_),2),
                      disease = c(rep("hvcd", nrow(hdata_)), rep("mcd", nrow(hdata_))))
  pdata$term <- factor(pdata$term, levels = rownames(hdata_)[order(hdata_$diff)])
  pdata$logp <- -1*log(pdata$p_value)
  pdata$logp[pdata$logp>16] <- 16
  pdata$logp[is.na(pdata$logp)] <- 0
  plts[[ct]] <- ggplot(pdata, aes(x=f_score, y=term, fill=logp))+
    geom_bar(position="dodge2", stat="identity")+
    scale_fill_viridis_c(option='D', direction = 1, end = max(pdata$logp)/16) +
    theme_minimal()+
    ggtitle(ct)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size=5))
}

plot_grid(plotlist = plts, col=3)

##########################################
# Figure 5 -- Stromal subtype
##########################################

###
# Dot plot

seu.stromal.sub <- readRDS("processed/02-annotation/stromal_subtypes.RDS")
Idents(seu.stromal.sub) <- "disease"
seu.stromal.sub <- subset(seu.stromal.sub, idents = "KFD", invert=T)

Idents(seu.stromal.sub) <- "pred_anno_2"

VlnPlot(seu.stromal.sub, features = "DLL4", 
        idents = "Fibroblastic stromal cell, type 2", group.by = "disease",
        assay = "RNA")
coi <- Cells(seu.stromal.sub)[seu.stromal.sub$pred_anno_2=="Fibroblastic stromal cell, type 2"]
small_seu <- subset(seu, cells = coi)
VlnPlot(small_seu, features = "DLL4", group.by = "disease")

sort(unique(seu.stromal.sub$pred_anno_2))
stromal_col_pal <- pal_nejm()
stromal_colors <- stromal_col_pal(length(unique(seu.stromal.sub$pred_anno_2)))
seu.stromal.sub$pred_anno_2 <- factor(seu.stromal.sub$pred_anno_2, levels = c("Fibroblastic stromal cell, type 1",
                                                                  "ACTA2+ perivascular reticular cells",
                                                                  "lymphatic endothelial cell",
                                                                  "FDC",
                                                                  "blood endothelial cell",
                                                                  "Fibroblastic stromal cell, type 2"))

markers <- unique(c('PDGFRA', 'PDGFRB', 'CXCL13', 'APOE', 'CCL21', 'CCL19', 'PDPN',
                    'CDH5', 'ENG', 'CD34', 'PECAM1',
                    'PROX1', 'PECAM1', 'PDPN',
                    'ACTA2', 'TAGLN', 'TPM2', 'PDGFRB',
                    'ACTA2', 'MYH11', 'MCAM',
                    'CCL19', 'CCL21', 'CXCL12', 'CXCL9',
                    'LUM', 'DCN', 'PDPN', 'PDGFRA',
                    'CXCL13', 'CLU', 'FDCSP', 'DES'))

DotPlot(seu.stromal.sub, features = markers, assay = "SCT", cols = c("gray","red")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

###
# UMAP
DefaultAssay(seu.stromal.sub) <- "SCT"
VariableFeatures(seu.stromal.sub) <- rownames(seu.stromal.sub@assays$SCT$scale.data)
seu.stromal.sub <- RunPCA(seu.stromal.sub)
ElbowPlot(seu.stromal.sub)
seu.stromal.sub <- FindNeighbors(seu.stromal.sub, dims = 1:10)
seu.stromal.sub <- FindClusters(seu.stromal.sub, resolution = 0.5)
seu.stromal.sub <- RunUMAP(seu.stromal.sub, dims = 1:6, min.dist = 0.3, spread = 0.5, 
                           repulsion.strength = 2, negative.sample.rate = 15, 
                           local.connectivity = 10, n.neighbors = 20)

DimPlot(seu.stromal.sub)
DimPlot(seu.stromal.sub , group.by = "pred_anno_2", cols=stromal_colors,
        pt.size = 4, raster = T)
DimPlot(seu.stromal.sub, group.by = "pred_anno_2", cols=stromal_colors, split.by = "disease", ncol = 3, 
        raster = T, pt.size = 4)+
  theme(legend.position = "none")

###
# compositions
cell_comp <- seu.stromal.sub@meta.data %>% dplyr::group_by(sample_name, pred_anno_2) %>%
  dplyr::summarise("cell_count"=length(pred_anno_2))

cell_comp_sums <-table(seu.stromal.sub$pred_anno_2)
cell_comp$rel_count <- apply(cell_comp, 1, function(x) as.numeric(x[3]) / cell_comp_sums[x[2]])

cell_comp$pred_anno_2 <- factor(cell_comp$pred_anno_2, 
                              levels = c("lymphatic endothelial cell", "blood endothelial cell", "Fibroblastic stromal cell, type 2", "Fibroblastic stromal cell, type 1",
                                         "ACTA2+ perivascular reticular cells", "FDC"),
                              ordered = T)

ggplot(cell_comp, aes(fill=sample_name, y=rel_count, x=pred_anno_2)) + 
  geom_bar(position="stack", stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5,
                                   hjust=1,
  ))

cell_comp <- seu.stromal.sub@meta.data %>% dplyr::group_by(pred_anno_2) %>%
  dplyr::summarise("cell_count"=length(pred_anno_2))
cell_comp$pred_anno_2 <- factor(cell_comp$pred_anno_2, 
                                levels = c("lymphatic endothelial cell", "blood endothelial cell", "Fibroblastic stromal cell, type 2", "Fibroblastic stromal cell, type 1",
                                           "ACTA2+ perivascular reticular cells", "FDC"),
                                ordered = T)
cell_comp$pct.comp <- cell_comp$cell_count/sum(cell_comp$cell_count)
ggplot(cell_comp, aes(y=pct.comp, x=pred_anno_2)) + 
  geom_bar(position="stack", stat="identity") +
  theme_minimal()

###
# gene expression/violin plots
DefaultAssay(seu.stromal.sub) <- "RNA"
VlnPlot(seu.stromal.sub, features = "VEGFA", 
        assay = "SCT", layer="data", 
        group.by = "pred_anno_2", pt.size = 0)

FeaturePlot(seu.stromal.sub, features = "VEGFA", split.by = "disease",
            cols=c("gray","red"), order = T)


# IL6 module
seu.stromal.sub <- AddModuleScore(seu.stromal.sub, features = il6_mod, 
                                  assay = "SCT", name = "IL6_mod_sct",
                                  nbin=10)

ct<-"Fibroblastic stromal cell, type 2"
VlnPlot(seu.stromal.sub, features = "IL6_mod_sct1", 
        group.by = "pred_anno_2",
        pt.size = 0)


VlnPlot(seu.stromal.sub, features = c("IL6_mod_sct1"),
        group.by = "sample_name", idents = ct,
        y.max = 2, same.y.lims = T,
        assay = "SCT", raster = F, pt.size = 0, alpha = 0.5,
        cols = donor_color)

FeaturePlot(seu.stromal.sub, features = "IL6_mod_sct1", split.by = "disease",
            cols=c("gray","red"), order = T)

###
# adding veg-f and IL-6 to dotplot?

DotPlot(seu.stromal.sub, features = c(markers,"VEGFA","IL6_mod_sct1"), assay = "SCT", cols = c("gray","red")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


Idents(seu.stromal.sub) <- "pred_anno_2"

pdf(
  here("figures/Figure5/Violin_VEGF_bySample_SCT_raster.pdf"),
  height = 10,
  width = 10
)

plist <- list()
for(ct in unique(seu.stromal.sub$pred_anno_2)){
  p <- VlnPlot(seu.stromal.sub, features = c("VEGFA"),
               group.by = "sample_name", idents = ct,
               y.max = 2, same.y.lims = T,
               assay = "SCT", raster = T, pt.size = 0.001, alpha = 0.5,
               cols = donor_color) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(ct)
  plist[[ct]] <- p
}
plot_grid(plotlist = plist, ncol=2)

dev.off()

##########################################
# Figure 6 -- Ligand-receptor interactions
##########################################
library(liana)
library(magrittr)
library(circlize)
library(pheatmap)
library(ComplexHeatmap)
library(viridis)
library(org.Hs.eg.db)

loading_pal <- colorRamp2(breaks = c(0,1), colors = c("yellow","blue"))
# loading_pal <- colorRamp2(breaks = c(0,1), colors = c("white", "purple"))

sce_c2c <- readRDS("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration/code/07-LIANA/CD_checkpoint.RDS")

sce_liana <- readRDS("code/07-LIANA/CD_SCT_stromasubset_checkpoint.RDS")

lr_agg <- lapply(sce_liana@metadata$liana_res, liana_aggregate)

for(i in 1:length(lr_agg)){
  lr_agg[[i]]$full_lr_pair <- paste(lr_agg[[i]]$source, 
                                    lr_agg[[i]]$target, 
                                    lr_agg[[i]]$ligand.complex, 
                                    lr_agg[[i]]$receptor.complex, sep="_")
}

lr_sig_r <- lr_agg$R1_3seq$full_lr_pair[lr_agg$R1_3seq$aggregate_rank < 0.1]

###
# Re-ordering LIANA dot plots
#   We plot HVCD1, MCD2, MCD4
#   categories:
#     cytokine/inflamm -- IL, etc.
#     stromal -- collagenase, integrins
#     plasma -- SDC1, CD138

pdata <- lr_agg$MCD4_5seq %>% filter(aggregate_rank < 0.05, !(full_lr_pair %in% lr_sig_r))
source_groups = c("FDC", "FSC, type 1", "FSC, type 2", "PRC")
target_groups = c("Macrophages", "Naive B cells", "Plasma cells", "BEC", "LEC")
pdata_filter <- pdata[(pdata$source %in% source_groups) & (pdata$target %in% target_groups),]

common_sources <- c()
for(x in c("MCD4_5seq", "MCD2_3seq", "HVCD1_3seq")){
  pdata <- lr_agg[[x]] %>% filter(aggregate_rank < 0.05, !(full_lr_pair %in% lr_sig_r))
  source_groups = c("FDC", "FSC, type 1", "FSC, type 2", "PRC")
  target_groups = c("Macrophages", "Naive B cells", "Plasma cells", "BEC", "LEC")
  pdata_filter <- pdata[(pdata$source %in% source_groups) & (pdata$target %in% target_groups),]
  common_sources <- c(common_sources, pdata_filter$ligand.complex[1:30])
}
common_sources <- unique(common_sources)
ligand_go <- mapIds(org.Hs.eg.db, common_sources, column = "GO", keytype = "SYMBOL")

# bring the annotations back in
GO_category <- read.csv("code/07-LIANA/GO_clustering_simple.csv", row.names = 1)

pdf(
  here("figures/Figure6/dotplot_LIANA_GOsorted.pdf"),
  height = 10,
  width = 12
)

plist <- list()
for(x in c("MCD4_5seq", "MCD2_3seq", "HVCD1_3seq")){
  pdata <- lr_agg[[x]] %>% filter(aggregate_rank < 0.05, !(full_lr_pair %in% lr_sig_r))
  source_groups = c("FDC", "FSC, type 1", "FSC, type 2", "PRC")
  target_groups = c("Macrophages", "Naive B cells", "Plasma cells", "BEC", "LEC")
  pdata_filter <- pdata[(pdata$source %in% source_groups) & (pdata$target %in% target_groups),]
  pdata_sort <- pdata_filter[1:30,]
  pdata_sort$cat <- GO_category[pdata_sort$ligand.complex,"final"]
  pdata_sort <- pdata_sort[order(pdata_sort$cat),]
  plist[[x]] <- liana_dotplot(pdata_sort,
                source_groups = c("FDC", "FSC, type 1", "FSC, type 2", "PRC"),
                target_groups = c("Macrophages", "Naive B cells", "Plasma cells", "BEC", "LEC"),
                ntop = 30)+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size=6),
          strip.text = element_text(angle = 45, vjust = 0.5, hjust=0.5, size=8),
          axis.text.y = element_text(size = 8))+
    ggtitle(x)
}

plot_grid(plotlist = plist, ncol=2)

dev.off()

###
# Cell2Cell stuff
factors <- get_c2c_factors(sce_c2c, group_col = "disease", sample_col = "edit.ident")

plot_context_heat(sce_c2c, group_col = "disease", sample_col="edit.ident", 
                  col = inferno(6))

progeny <- decoupleR::get_progeny(organism = 'human', top=5000) %>%
  select(-p_value)

# convert to LR sets
progeny_lr <- generate_lr_geneset(sce_c2c,
                                  resource = progeny)

mat <- factors$interactions %>%
  column_to_rownames("lr") %>%
  as.matrix()

# run enrichment analysis with decoupler
# (we fit a univariate linear model for each gene set)
# We don't consider genesets with minsize < 10
res <- decoupleR::run_ulm(mat = mat,
                          network = progeny_lr,
                          .source = "set",
                          .target = "lr",
                          minsize=10) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

res$condition <- factor(res$condition, levels = rev(c("Factor.1","Factor.9","Factor.8","Factor.12","Factor.4",
                                                  "Factor.2","Factor.6","Factor.11","Factor.10","Factor.3",
                                                  "Factor.7","Factor.5")),
                        ordered = T)


res %>% # sig/isnig flag
  # mutate(significant = if_else(p_adj <= 0.05, "signif.", "not")) %>%
  filter(p_adj <= 0.05) %>%
  ggplot(aes(x=source, y=condition, 
             colour=score, size=-log10(p_value+1e-36))) +
  geom_point() +
  # scale_colour_gradient2(high = "red", low="blue") +
  scale_colour_viridis_c()+
  scale_size_continuous(range = c(3, 12)) +
  # scale_shape_manual(values=c(21, 16)) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x="Pathway",
       y="Factor",
       colour="Activity"
  )



###
# scatter plots of enrichment weights

lrs <-  factors$interactions %>%
  left_join(progeny_lr, by="lr") %>%
  filter(set=="JAK-STAT") %>%
  select(lr, set, mor, loading = Factor.9) %>%
  mutate(lr = gsub(as.character(str_glue("\\^")), " -> ", lr)) %>%
  mutate(weight = if_else(mor >= 0, "positive", "negative"))
pf9 <-lrs %>%
  # only label those that are > x
  mutate(lr = if_else(loading>=0.001 & abs(mor) > 2, lr, "")) %>%
  ggplot(aes(x=mor, y=loading, colour=weight)) +
  # label only top 20
  stat_smooth(method = "lm", col = "red") +
  geom_point(alpha = 0.5) +
  ggrepel::geom_label_repel(aes(label = lr)) +
  theme_bw(base_size = 15) +
  scale_colour_manual(values = c("royalblue3", "red")) +
  labs(x="Pathway Weight", y="LR Loading")

###
# making sender-receiver loadings dotplot
factor_order <- c("Factor.1","Factor.9","Factor.8","Factor.12","Factor.4","Factor.2",
                  "Factor.6","Factor.11","Factor.10","Factor.3","Factor.7","Factor.5")

intersting_cells <- c("Lymphatics", "Macrophages", "Endothelial cells", "Monocytes", "Stromal cells", "Plasmacytoid dendritic cells", "Cytotoxic CD8 T cells", "Fibroblasts", "Activated and migratory cDC", "T Follicular Helper cells", "Follicular dendritic cells")
intersting_cells2 <- c("Lymphatics", "Macrophages", "Endothelial cells", "Monocytes", "Stromal cells", "Plasmacytoid dendritic cells", "Cytotoxic CD8 T cells", "Fibroblasts", "Activated and migratory cDC", "T Follicular Helper cells", "Follicular dendritic cells",
                       "Activated and memory B cells", "Naive B cells", "Germinal center B cells", "Plasma cells")

hdata <- data.frame(sce_c2c@metadata$tensor_res$senders) %>% column_to_rownames("celltype")
hdata <- hdata[intersting_cells,]
# pheatmap(t(hdata), scale = "none", cluster_rows = F)
# Heatmap(t(hdata), col=brewer.pal(9, "YlGnBu"), row_order = factor_order)
Heatmap(t(hdata), col=inferno(8), row_order = factor_order)

hdata <- data.frame(sce_c2c@metadata$tensor_res$receivers) %>% column_to_rownames("celltype")
hdata <- hdata[intersting_cells2,]
# pheatmap(t(hdata), scale = "none", cluster_rows = F)
Heatmap(t(hdata), col=inferno(6), row_order = factor_order)

##########################################
# Figure N -- Xenium integration
#   focus on factors 2 and 4
##########################################

xen_genes <- read.csv("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration/code/11-Xenium_integration/xen_seq_gene_intersection.csv",
                      row.names = 1)

c2c_intxn <- sce_c2c@metadata$tensor_res$interactions
intxn_genes <- sapply(as.character(c2c_intxn$lr), function(s) strsplit(s, "\\^"))
# keep <- unlist(lapply(intxn_genes, function(x) any(x %in% xen_genes$X0)))
keep <- unlist(lapply(intxn_genes, function(x) all(x %in% xen_genes$X0)))
c2c_intxn <- c2c_intxn[keep,]

hdata <- c2c_intxn %>% slice_max(order_by = Factor.2, n = 20) %>% tibble::column_to_rownames("lr")
pheatmap(hdata, cluster_cols = F, scale = "row")

Idents(seu) <- "pred_anno_2"
unique(seu$pred_anno_2)
VlnPlot(seu, features = "DLL4", idents = "Fibroblasts", group.by = "disease")
VlnPlot(seu, features = "NOTCH1", idents = "Memory or effector T cells", group.by = "disease")
