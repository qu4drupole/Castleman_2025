library(ggplot2)
library(reshape2)
library(here)
library(circlize)
setwd("/mnt/isilon/cscb/Projects/10x/pillaiv/SCTC-VP-22/vinodh-5-3-integration")

here::i_am("code/07-LIANA/CD_LIANA_plot_assist.R")

DimPlot(integration.combined, group.by = "Man_anno", label = T) + theme(legend.position = "none")

table(integration.combined$Man_anno[integration.combined$orig.ident=="HVCD1"])

source_cells <- c("FDC", "FSC, type 1", "FSC, type 2", "PRC")

pdata <- lr_agg$HVCD3_5seq %>% filter(aggregate_rank < 0.1, !(full_lr_pair %in% lr_sig_r))

pdf(
  here("code/07-LIANA/outs/CD/LIANA/LIANA_MCD1_stromaSubtype.pdf"),
  height = 11,
  width = 11
)

heat_freq(pdata)

chord_freq(pdata, 
           source_groups = source_cells,
           target_groups = c("Macrophages", "Naive B cells", "Plasma cells", "BEC", "LEC"),
           scale=F,
           cex=0.5,
           adj = c(0,0.5))

chord_freq(pdata, 
           source_groups = source_cells,
           target_groups = unique(sce$pred_anno_2)[!unique(sce$pred_anno_2) %in% source_cells],
           scale=T,
           cex=0.5,
           adj = c(0,0.5))


liana_dotplot(pdata,
              source_groups = c("FDC", "FSC, type 1", "FSC, type 2", "PRC"),
              target_groups = c("Macrophages", "Naive B cells", "Plasma cells", "BEC", "LEC"),
              ntop = 30)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size=6),
        strip.text = element_text(angle = 45, vjust = 0.5, hjust=0.5, size=8),
        axis.text.y = element_text(size = 8))

for(ct in source_cells){
  p <- liana_dotplot(pdata,
                source_groups = ct,
                target_groups = c("Macrophages", "Naive B cells", "Plasma cells", "BEC", "LEC"),
                ntop = 30)+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size=6),
          strip.text = element_text(angle = 45, vjust = 0.5, hjust=0.5, size=8),
          axis.text.y = element_text(size = 8))
  if(nrow(p$data)>0){
    print(p)
  }
}

dev.off()

# [1] "Lymphatics"                      "Follicular dendritic cells"      "Fibroblasts"                     "Classic Dendritic cells 2 and 3"
# [5] "NK cells"                        "Classic Dendritic cells 1"       "Activated and memory B cells"    "Granulocytes"                   
# [9] "Cytotoxic CD8 T cells"           "Plasmacytoid dendritic cells"    "Macrophages"                     "B cell subset"                  
# [13] "Monocytes"                       "Regulatory T cells"              "Endothelial cells"               "Naive CD4 and Naive CD8 T cells"
# [17] "Naive B cells"                   "Memory or effector T cells"      "Activated and migratory cDC"     "Plasma cells"                   
# [21] "Stromal cells"                   "Proliferating plasma cells"      "Germinal center B cells"         "T follicular helper cells" 

ct <- "Plasmacytoid dendritic cells"
pdata <- all_res[all_res$sourceKFD == ct,]
pdata <- pdata[order(pdata$d_rank, decreasing = T),]

# how many of common interactions are "upregulated"/"downregulated"
hdata_R <- all_res %>% filter(d_rank < 0, abs(d_rank) > 1) %>%
  pivot_wider(id_cols = sourceKFD, names_from = targetKFD, values_from = targetKFD, values_fn = length) %>%
  replace(is.na(.), 0) %>%
  melt
hdata_R$sourceKFD <- factor(hdata_R$sourceKFD, levels = unique(all_res$sourceKFD), ordered = T)
hdata_R$variable <- factor(hdata_R$variable, levels = unique(all_res$targetKFD), ordered = T)

ggplot(hdata_R, aes(x=variable, y=sourceKFD, fill=value))+
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8))

hdata_KFD <- all_res %>% filter(d_rank > 0, abs(d_rank) > 1) %>%
  pivot_wider(id_cols = sourceKFD, names_from = targetKFD, values_from = targetKFD, values_fn = length) %>%
  replace(is.na(.), 0) %>%
  melt
hdata_KFD$sourceKFD <- factor(hdata_KFD$sourceKFD, levels = unique(all_res$sourceKFD), ordered = T)
hdata_KFD$variable <- factor(hdata_KFD$variable, levels = unique(all_res$targetKFD), ordered = T)

hdata_KFD <- hdata_KFD[order(hdata_KFD$sourceKFD, hdata_KFD$variable),]
hdata_R <- hdata_R[order(hdata_R$sourceKFD, hdata_R$variable),]

hdata_KFD$diff_diff <- hdata_KFD$value - hdata_R$value
ggplot(hdata_KFD, aes(x=variable, y=sourceKFD, fill=diff_diff))+
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=8)) +
  scale_fill_gradient2()
