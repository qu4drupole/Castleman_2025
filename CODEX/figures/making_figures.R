library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)
library(reshape2)
# Follicle/CD21 EDA
# setwd("Z:/Projects/codex/pillaiv/SCTC-VP-15/code/")
setwd("/mnt/isilon/cscb/Projects/codex/pillaiv/SCTC-VP-15/code/")

# dat <- read_xlsx("CD_analysis/aSMA_CD21_analysis/aggregated_region_stats_formatted_method2.xlsx", sheet = "for_R") %>% data.frame()
# dat <- read.csv("CD_analysis/aSMA_CD21_analysis/aggregated_region_stats_method2.csv")
dat <- read.csv("CD_analysis/aSMA_CD21_analysis/aggregated_region_stats_method3.csv")
dat$donor <- sapply(dat$X, function(s) strsplit(s,"_")[[1]][1]) # no formatting
# dat$donor <- sapply(dat$Sample, function(s) strsplit(s,"_")[[1]][1])
dat$disease <- sapply(dat$donor, function(s) gsub("\\d","",s))
dat$donor <- sapply(dat$donor, function(s) gsub("R","RLN",s))
dat <- dat[dat$avg..CD21.brightness..cluster.!=0,]
dat$sq_area <- sqrt(dat$avg..cluster.area*(0.377**2))


donor_colors = c("#66cc66", "#339933", "#9999cc", "#666699", "#ff9933", "#ff6633", "#cc3333")
names(donor_colors) <- c("RLN1","RLN2","HVCD1","HVCD2","MCD1","MCD3","MCD4")
dat$color <- donor_colors[dat$donor]

ggplot(dat, aes(dat$avg..CD21.cell.signal..cluster., dat$avg..entropy.threshold..cluster.))+
  geom_point(aes(color=disease))

dat$entropy_ratio <- dat$avg..entropy.threshold..cluster./dat$avg..CD21.cell.signal..cluster.

pdata <- melt(dat, id.vars = c("X","donor","disease","color"))

ggplot(pdata, aes(x=disease, y=value, group=disease))+
  geom_boxplot(fill="white", position="dodge", outlier.shape = NA)+
  geom_jitter(aes(color=donor), width = 0.1)+
  scale_color_manual(values = donor_colors[sort(names(donor_colors))])+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~variable, nrow = 3, scales = 'free')


id_df <- dat %>% select(is.character)
num_df <- dat %>% select(is.numeric)
num_df <- scale(num_df)

pdata <- cbind(id_df, num_df)
pdata <- melt(pdata, id.vars = c("Sample","donor","disease"))

ggplot(pdata, aes(x=variable, y=value))+
  geom_boxplot(aes(fill=disease), position="dodge", outlier.shape = NA)+
  geom_jitter(aes(color=donor), width = 0.1)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

###
# t.tests

a <- dat$avg..entropy.threshold..cluster.[dat$disease=="MCD"]
b <- dat$avg..entropy.threshold..cluster.[dat$disease=="R"]
t.test(a,b)




# Random

FeaturePlot(seu, features = c("FTH1", "FTL"))
Idents(seu) <- "pred_anno_2"
VlnPlot(seu, features = c("FTH1", "FTL"), 
        group.by = "disease", idents = "Activated and memory B cells",
        assay = "RNA")

plot(dat_test@assays$combat@data["FTL",], dat_test@assays$RNA@data["FTL",])

tgenes <- rownames(dat_test@assays$combat)
a <- rowSums(dat_test@assays$RNA@data[tgenes,])
b <- rowSums(dat_test@assays$combat@data[tgenes,])
assay_diff <- abs(a-b); assay_diff <- assay_diff[order(assay_diff, decreasing = T)]

plot(dat_test@assays$SCT@counts["IL6ST",], dat_test@assays$RNA@data["IL6ST",])

