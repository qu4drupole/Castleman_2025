require(Seurat)
require(ggplot2)
require(DoubletFinder)
require(SoupX)
require(sva)
require(Matrix)

#' Load Cell Ranger Multi Results into SoupX
#'
#' Takes in 10x Cell Ranger-generated `multi` results and produces a modified table of counts, with background contamination removed.
#' Learn more here \url{https://cran.r-project.org/web/packages/SoupX/SoupX.pdf} or here \url{https://github.com/constantAmateur/SoupX}
#'
#' @param raw file path of raw multi results; found in outs/multi/count/raw_feature_bc_matrix/
#' @param filtered file path of filtered multi results; found in count/sample_filtered_feature_bc_matrix/
#' @param clusters file path of clusters.csv; found in count/analysis/clustering/gene_expression_graphclust/clusters.csv
#'
#' @return SoupX object adjusted for ambient RNA contamination
#' @export
#'
#' @examples
#' gex_multi_soup_raw <-
#'    R.utils::getAbsolutePath(Sys.glob(
#'   file.path(
#'     "../../SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/code/02-multi/runs/*_VDJ_GEX/outs/multi/count/raw_feature_bc_matrix/"
#'   )
#' ))
#' gex_multi_soup_filtered <-
#'   R.utils::getAbsolutePath(Sys.glob(
#'     file.path(
#'       "../../SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/code/02-multi/runs/*_VDJ_GEX/outs/per_sample_outs/*_VDJ_GEX/count/sample_filtered_feature_bc_matrix/"
#'     )
#'   ))
#' gex_multi_soup_clusters <-
#'   R.utils::getAbsolutePath(Sys.glob(
#'     file.path(
#'       "../../SCTC-VP-16/Data_Analysis/vinodh-vdj-gex/code/02-multi/runs/*_VDJ_GEX/outs/per_sample_outs/*_VDJ_GEX/count/analysis/clustering/gene_expression_graphclust/clusters.csv"
#'     )
#'   ))
#'
#' names(gex_multi_soup_raw) <- NULL
#' names(gex_multi_soup_filtered) <- NULL
#' names(gex_multi_soup_clusters) <- NULL
#'
#' gex_obj <- list()
#'
#' for (i in 1:length(gex_multi_soup_raw)) {
#'   gex_obj[[i]] <-
#'     soupx_load_multi(gex_multi_soup_raw[i],
#'                      gex_multi_soup_filtered[i],
#'                      gex_multi_soup_clusters[i])
#' }
soupx_load_multi <- function(raw, filtered, clusters) {
  drops <- SoupX::Read10X(raw)
  dat.filtered <-
    SoupX::Read10X(filtered)
  
  sc = SoupX::SoupChannel(drops, dat.filtered)
  
  cr.clusters <-
    read.csv(clusters)
  clust <- cr.clusters$Cluster
  names(clust) <- cr.clusters$Barcode
  sc = SoupX::setClusters(sc, clust)
  
  # depending on the read depth, these last two steps may need to be manually adjusted
  #sc = autoEstCont(sc, soupQuantile = 0.5)
  #out = adjustCounts(sc)
  return(sc)
}


#' Load Cell Ranger count data
#'
#' @param dir root of CellRanger count output
#' @param s list of sample names/directories output from CellRanger count
#' @param extra file path modification if multiple CellRanger count instances exist
#'
#' @return Seurat object
#' @export
#'
#' @examples
#' snames <- c('2hr_1', '2hr_2', '2hr_3', '2hr_ctrl', '2hr_untag', '12hr_1', '12hr_2', '12hr_3')
#' seq2.dat <- mclapply(snames, function(x) load.seur.data(x, "rerun"), mc.cores = 4)
#' snames <- paste(snames, "seq2", sep="_")
#' names(seq2.dat) <- snames
load_seur_data <- function(dir, s, extra = "") {
  dpath <- paste(dir,
                 extra,
                 s,
                 'outs/filtered_feature_bc_matrix/',
                 sep = '/')
  dat.raw <- Seurat::Read10X(dpath)
  dat <-
    Seurat::CreateSeuratObject(
      counts = dat.raw$`Gene Expression`,
      # I guess this assumes multi...
      project = s,
      min.cells = 200,
      min.features = 200
    )
  dat[["percent.mt"]] <- Seurat::PercentageFeatureSet(dat, pattern = "^MT-")
  dat[["percent.ribo"]] <-
    Seurat::PercentageFeatureSet(dat, pattern = "^RP[LS]|^MRPL")
  return(dat)
}

#' Adds layers to Seurat object, percent mitochondrial and ribosomal genes
#'
#' Matches string prefixes to gene symbols ("MT" and "RPS") to find percent mitochondrial and ribosomal genes per cell
#'
#' @param seuratObj A Seurat object with RNA assay
#'
#' @return A Seurat object with two new layers
#' @export
#'
#'
#' @examples
#' for (i in 1:length(gex_obj)) {
#'  gex_obj[[i]] <- get.percent(gex_obj[[i]])
#' }
get_percent <- function(seuratObj) {
  seuratObj[["percent.mt"]] <-
    Seurat::PercentageFeatureSet(object = seuratObj,
                                 pattern = "^MT-",
                                 assay = "RNA")
  
  seuratObj[["percent.rps"]] <-
    Seurat::PercentageFeatureSet(object = seuratObj,
                                 pattern = "^RP[LS]|^MRPL",
                                 assay = "RNA")
  
  return(seuratObj)
}


#' Plots mitochondrial and ribosomal genes with ggplot2
#'
#' @param seuratObj A Seurat object with percent mitochondrial and ribosomal genes calculated with `get.percent()`
#' @param mito.cutoff Cutoff of percent mitochondrial genes, default 10
#' @param rps.cutoff Cutoff of percent ribosomal genes, default 10
#' @param title Title of ggplot, default NULL
#'
#' @return A dual-faceted ggplot printed to active device
#' @export
#'
#' @examples
#'
#' plot_mito(gex_obj[[1]], title = gex_names[[1]])
#' plot_mito(gex_obj[[2]], title = gex_names[[2]])
#' plot_mito(gex_obj[[3]], title = gex_names[[3]])
#' plot_mito(gex_obj[[4]], title = gex_names[[4]])
#' plot_mito(gex_obj[[6]], title = gex_names[[6]])
#' plot_mito(gex_obj[[5]], title = gex_names[[5]])
plot_mito <-
  function(seuratObj,
           mito.cutoff = 10,
           rps.cutoff = 10,
           title = NULL) {
    count_mt <-
      Seurat::FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
      ggplot2::geom_hline(yintercept = 0.05,
                          linetype = "dashed",
                          color = "black") +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ggtitle(title) + ggplot2::geom_hline(yintercept = mito.cutoff) +
      ggplot2::geom_point(aes(colour = cut(percent.mt, c(0, mito.cutoff))))
    
    count_rps <-
      Seurat::FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.rps") +
      ggplot2::geom_hline(yintercept = 0.05,
                          linetype = "dashed",
                          color = "black") +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ggtitle(NULL) + ggplot2::geom_hline(yintercept = rps.cutoff) +
      ggplot2::geom_point(aes(colour = cut(percent.rps, c(0, rps.cutoff))))
    
    
    return((count_mt + count_rps))
  }

#' Generate feature plots of transcript counts vs gene counts
#'
#' @param seuratObj A Seurat object
#' @param count.cutoff A cutoff indicating max counts per cell; default 10k
#' @param downsample A number representing number of cells to subset
#'
#' @return A FeatureScatter printed to device
#' @export
#'
#' @examples
#'
#' seurat.counts(gex_obj[[1]], count.cutoff = 22000) + ggtitle(gex_names[[1]])
#' seurat.counts(gex_obj[[2]], count.cutoff = 52000) + ggtitle(gex_names[[2]])
#' seurat.counts(gex_obj[[3]], count.cutoff = 10000) + ggtitle(gex_names[[3]])
#' seurat.counts(gex_obj[[4]], count.cutoff = 70000) + ggtitle(gex_names[[4]])
#' seurat.counts(gex_obj[[5]], count.cutoff = 65000) + ggtitle(gex_names[[5]])
#' seurat.counts(gex_obj[[6]], count.cutoff = 60000) + ggtitle(gex_names[[6]])
#'
seurat_counts <-
  function(seuratObj,
           count.cutoff = 10000,
           downsample = 0) {
    if (downsample > 0) {
      seuratCounts <- FeatureScatter(
        subset(seuratObj, downsample = downsample),
        feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA"
      ) +
        theme(legend.position = "none") +
        ggtitle(NULL) +
        geom_vline(xintercept = count.cutoff) +
        geom_point(aes(colour = cut(nCount_RNA, c(
          -Inf, count.cutoff
        ))))
      return(seuratCounts)
    }
    else {
      seuratCounts <-
        FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
        theme(legend.position = "none") +
        ggtitle(NULL) +
        geom_vline(xintercept = count.cutoff) +
        geom_point(aes(colour = cut(nCount_RNA, c(
          -Inf, count.cutoff
        ))))
      return(seuratCounts)
    }
  }

#' Preprocess Seurat object
#'
#' Prunes cells from Seurat object that do not meet cutoff criteria. Stores parameter values in `seuratObj@assays$RNA@misc` as a data frame.
#'
#' @param counts Seurat object
#' @param counts_name Name of sample corresponding to Seurat object
#' @param count.cutoff Cutoff of counts used in seurat_counts, default 10k
#' @param mito.cutoff Cutoff of percent mitochondrial genes per cell, default 10
#' @param rps.cutoff Cutoff of percent ribosomal genes per cell, default 10
#' @param cc_adjust Boolean, whether to regress out cell cycle scores
#'
#' @return Seurat object that has been QC'd, SCTransformed, PCA'd, Clustered and UMAP.
#' @export
#'
#' @examples
#' for (i in 1:length(gex_obj)) {
#' gex_obj[[i]] <-
#'   seurat.process(gex_obj[[i]], counts_name = gex_names[[i]])
#' }
#' # with future
#' combined.dat <- future.apply::future_lapply(combined.dat, seurat.process, future.seed = TRUE)
seurat_process_v2 <-
  function(counts,
           counts_name,
           feature.cutoff = c(0,50000),
           count.cutoff,
           # can we store these options in the seurat object?
           mito.cutoff = 10,
           rps.cutoff = 10,
           cc_adjust = FALSE) {
    if (class(counts)[1] != "Seurat") {
      seuratObj <- CreateSeuratObject(
        counts = counts,
        project = counts_name,
        min.cells = 10,
        min.features = 200
      )
    } else {
      # if `counts` is already a Seurat object:
      seuratObj <- counts
      DefaultAssay(seuratObj) <- "RNA"
    }
    
    # Absolutely imperative that these thresholds are set PER SAMPLE
    seuratObj <-
      subset(
        seuratObj,
        subset = nFeature_RNA > 200 &
          nCount_RNA < count.cutoff &
          nFeature_RNA > feature.cutoff[1] & nFeature_RNA < feature.cutoff[2] &
          percent.mt < mito.cutoff &
          percent.rps < rps.cutoff
      )
    
    seuratObj@assays$RNA@misc <-
      data.frame (
        parameter  = c("Count Cutoff", "Mito Cutoff", "Ribo Cutoff", "CC Adjust?"),
        value = c(count.cutoff, mito.cutoff, rps.cutoff, cc_adjust)
      )
    
    # normalize
    
    # TODO expanding limit to prevent some arbitrary error msg
    # options(future.globals.maxSize = 20971520000)
    if (cc_adjust) {
      seuratObj <-
        CellCycleScoring(
          seuratObj,
          s.features = cc.genes$s.genes,
          g2m.features = cc.genes$g2m.genes,
          set.ident = TRUE
        )
      seuratObj <-
        FindVariableFeatures(seuratObj,
                             selection.method = "vst",
                             nfeatures = 2000)
      if (dim(seuratObj)[2] < 10000) {
        seuratObj <-
          SCTransform(
            seuratObj,
            ncells = 2000,
            vars.to.regress = c("S.Score", "G2M.Score")
          )
      }
      else{
        seuratObj <-
          SCTransform(
            seuratObj,
            ncells = 0.25 * ncol(seuratObj),
            vars.to.regress = c("S.Score", "G2M.Score")
          )
      }
    }
    else {
      seuratObj <- SCTransform(seuratObj, ncells = 2000)
    }
    seuratObj <- RunPCA(seuratObj)
    
    eb <- ElbowPlot(seuratObj, reduction = "pca")
    
    seuratObj <-
      FindNeighbors(seuratObj, dims = 1:10) # we pick this based on elbow plot
    seuratObj <- FindClusters(seuratObj, resolution = 0.5)
    seuratObj <- RunUMAP(seuratObj, dims = 1:10)
    return(seuratObj)
  }



#' Running paramSweep and summarizeSweep from DoubletFinder
#'
#' @param sample Seurat object
#' @param cores Number of cores, default 1
#'
#' @return Summary states from summarizeSweep
#' @export
#'
#' @examples
#' for (i in 1:length(gex_obj)) {
#'   message(crayon::red$underline$bold(
#'   paste0("Finding Parameter Sweep with DoubletFinder for ", gex_names[[i]])
#' ))
#' gex_pk[[i]] <-
#'   sum.sweep(gex_obj[[i]])
#'   message(crayon::red$underline$bold(paste0("Saving ", gex_names[[i]])))
#'   save(gex_pk, file = here("gex_obj_postDoublet.RData"))
#' }
sum_sweep <- function(sample, cores = 1) {
  message("Parameter sweep...")
  
  system.time({
    sweep.res.sample <-
      paramSweep_v3(sample,
                    PCs = 1:10,
                    sct = TRUE,
                    num.cores = cores)
  })
  
  sweep.stats.sample <-
    summarizeSweep(sweep.res.sample, GT = FALSE)
  
  return(sweep.stats.sample)
}


#' Plots result from find.pk
#'
#' @param sweep.stats.sample Summary stats generated by sum.sweep
#' @param title Plot title, default NA
#'
#' @return R plot printed to active device
#' @export
#'
#' @examples
#'
#' pdf(here("processed/01-5_3-comp/ParamSweepPeaks.pdf"))
#'
#' plot.pk(gex_pk[[1]], title = gex_names[[1]])
#' plot.pk(gex_pk[[2]], title = gex_names[[2]])
#' plot.pk(gex_pk[[3]], title = gex_names[[3]])
#' plot.pk(gex_pk[[4]], title = gex_names[[4]])
#' plot.pk(gex_pk[[5]], title = gex_names[[5]])
#' plot.pk(gex_pk[[6]], title = gex_names[[6]])
#'
#' dev.off()
plot_pk <- function(sweep.stats.sample,  title = NA) {
  bcmvn_sample <- find.pK(sweep.stats.sample)
  x <-
    plot(
      x = bcmvn_sample$ParamID,
      y = bcmvn_sample$MeanAUC,
      pch = 18,
      col = "black",
      cex = 0.75,
      xlab = NA,
      ylab = NA,
      main = title
    )
  x <-
    lines(
      x = bcmvn_sample$ParamID,
      y = bcmvn_sample$MeanAUC,
      col = "black",
      lty = 2
    )
  par(new = TRUE)
  x <-
    plot(
      x = bcmvn_sample$ParamID,
      y = bcmvn_sample$BCmetric,
      pch = 16,
      col = "#41b6c4",
      cex = 0.75
    )
  axis(side = 4)
  x <- lines(x = bcmvn_sample$ParamID,
             y = bcmvn_sample$BCmetric,
             col = "#41b6c4")
  
  return(x)
  
}

#' Predict and remove doublets from Seurat object
#'
#' @param sample Seurat object
#' @param sweep.stats.sample Summary stats from sum.sweep()
#' @param rm.doubs logical; default FALSE; removes doublets from `sample` if TRUE
#'
#' @return A Seurat object with identified doublets
#' @export
#'
#' @examples
#'
#' for (i in 1:length(gex_obj)) {
#'   message(crayon::red$underline$bold(paste0(
#'     "Finding Doublets with DoubletFinder for ", gex_names[[i]]
#'   )))
#'   gex_doubs[[i]] <-
#'     doubFinder(gex_obj[[i]], sweep.stats.sample = gex_pk[[i]])
#'   message(crayon::red$underline$bold(paste0("Saving ", gex_names[[i]])))
#'   save(gex_doubs, file = here("gex_obj_postDoublet_actual.RData"))
#' }
#' # with future
#' combined.dat <- future.apply::future_lapply(combined.dat, doubFinder, future.seed = TRUE)
doubFinder <-
  function(sample, sweep.stats.sample, rm.doubs = FALSE) {
    bcmvn_sample <- find.pK(sweep.stats.sample) # save this plot
    
    message("Selecting pK Value...")
    # alternative method
    maxima <- which(diff(sign(diff(
      bcmvn_sample$BCmetric
    ))) == 2)
    if (length(maxima) == 0) {
      maxima <- which.max(bcmvn_sample$BCmetric)
    } else {
      maxima <-
        maxima[which(bcmvn_sample$BCmetric[maxima] > min(bcmvn_sample$BCmetric))]
    }
    pK_val <- as.numeric(levels(bcmvn_sample$pK)[maxima[1]])
    
    
    # peaks <-
    #   bcmvn_sample[bcmvn_sample$BCmetric > (max(bcmvn_sample$BCmetric) / 2),]
    
    # # consider adding trycatch
    # # https://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r
    
    # # pay special attention to the different types of boolean operators in R
    # # https://medium.com/biosyntax/single-or-double-and-operator-and-or-operator-in-r-442f00332d5b
    # for (v in 1:nrow(peaks)) {
    #   if (nrow(peaks) == 1) {
    #     pK_val <- as.numeric(as.vector(peaks[v, ]$pK))
    #     break
    #   }
    #   else if (peaks[v, ]$BCmetric > peaks[v + 1, ]$BCmetric |
    #            is.na(peaks[v + 1, ]$BCmetric)) {
    #     pK_val <- as.numeric(as.vector(peaks[v, ]$pK))
    #     break
    #   }
    #   else {
    #     next
    #   }
    # }
    
    message("Homotypic Doublet Proportion Estimation")
    homotypic.prop <-
      modelHomotypic(sample@meta.data$SCT_snn_res.0.5) ## ex: annotations <- sample@meta.data$ClusteringResults
    assumed_db <- (0.0008*ncol(sample) - 0.067)/100  
    nExp_poi <-
      round(assumed_db * nrow(sample@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    
    ## Run DoubletFinder with varying classification stringencies
    ## Remember to update with appropriate params
    # options(future.globals.maxSize = 2097152000)
    
    message("Running doubletFinder...")
    sample <-
      doubletFinder_v3(
        sample,
        PCs = 1:10,
        pN = 0.25,
        pK = pK_val,
        nExp = nExp_poi,
        reuse.pANN = FALSE,
        sct = TRUE
      )
    
    # grep for the column name that has pANN at the start
    pANN_col <-
      colnames(sample@meta.data)[grep("pANN*", colnames(sample@meta.data))]
    
    sample <-
      doubletFinder_v3(
        sample,
        PCs = 1:10,
        pN = 0.25,
        pK = pK_val,
        nExp = nExp_poi.adj,
        reuse.pANN = pANN_col,
        sct = TRUE
      )
    
    class_col <-
      colnames(sample@meta.data)[grep("*classification*", colnames(sample@meta.data))]
    
    class_col <- class_col[length(class_col)]
    
    colnames(sample@meta.data)[ncol(sample@meta.data)] <- "doublet"
    
    if (rm.doubs) {
      sample <- subset(sample, subset = doublet == "Singlet")
    }
    
    return(sample)
    
  }

#' Visualize doublets on UMAP embedding w/ RNA counts for reference
#'
#' @param seuratObj Seurat object with predicted doublets
#' @param title Sample name
#'
#' @return Printed ggplot2 to active device
#' @export
#'
#' @examples
#' pdf(
#'   here("processed/01-5_3-comp/GEX_doubletClusters.pdf"),
#'     height = 8,
#'     width = 14
#'   )
#'
#' viz.doubs(gex_doubs[[1]], title = gex_names[[1]])
#' viz.doubs(gex_doubs[[2]], title = gex_names[[2]])
#' viz.doubs(gex_doubs[[3]], title = gex_names[[3]])
#' viz.doubs(gex_doubs[[4]], title = gex_names[[4]])
#' viz.doubs(gex_doubs[[5]], title = gex_names[[5]])
#' viz.doubs(gex_doubs[[6]], title = gex_names[[6]])
#'
#' dev.off()
viz_doubs <- function(seuratObj, title = "Sample") {
  print(
    DimPlot(seuratObj, reduction = 'umap', group.by = 'doublet') +
      ggtitle(paste0("Doublet Clusters for ", title)) + FeaturePlot(seuratObj,
                                                                    reduction = "umap",
                                                                    features = "nCount_RNA")
  )
}

#' ComBat Batch-Correct SCT Counts in Seurat Object
#'
#' @param seuratObj A Seurat object with an SCT assay
#' @param seuratBatch A Vector of batch labels for each
#' @param merge Logical indicating whether to merge; if true, supply vector of Seurat objects as seuratObj
#'
#' @return Seurat object with combatBatch Assay
#' @export
#'
#' @examples
#' seuratObj <-
#'   merge(
#'     rna_seq[[1]],
#'     y = c(
#'       rna_seq[[2]],
#'       rna_seq[[3]],
#'       rna_seq[[4]],
#'       rna_seq[[5]],
#'       rna_seq[[6]],
#'       gex_doubs[[1]],
#'       gex_doubs[[2]],
#'       gex_doubs[[3]],
#'       gex_doubs[[4]],
#'       gex_doubs[[5]],
#'       gex_doubs[[6]]
#'     )
#'   )
#'
#' combat_proccess <-
#'   process_combat(seuratObj, seuratObj$comp.ident)
process_combat <- function(seuratObj, seuratBatch, merge = FALSE) {
  if (merge == TRUE) {
    seuratObj <- merge(seuratObj)
  }
  
  
  pheno <- seuratObj@meta.data
  
  batch <- seuratBatch
  
  edata <-
    as.matrix(GetAssayData(seuratObj, assay = "SCT", slot = "data"))
  
  mod0 <- model.matrix( ~ 1, data = pheno)
  
  combat_edata = ComBat(
    dat = edata,
    batch = batch,
    mod = mod0,
    par.prior = TRUE,
    prior.plots = FALSE
  )
  
  seuratObj[["combatBatch"]] <-
    CreateAssayObject(data = combat_edata)
  
  DefaultAssay(seuratObj) <- "combatBatch"
  
  seuratObj <-
    ScaleData(seuratObj, verbose = TRUE, assay = "combatBatch")
  
  return(seuratObj)
}


#' Delete genes from Seurat object with prefix
#'
#' Based on [this GitHub issue.](https://github.com/satijalab/seurat/issues/2610#issuecomment-585935810)
#'
#' @param seuratObj Input Seurat object
#' @param gene_name Prefix of genes to filter out (delete); also takes vector of strings
#' @param assay Assay, defaults to RNA
#'
#' @return Seurat object with deleted genes
#' @export
#'
#' @examples
#' seuratObj <- del_genes(seuratObj, "HLA")
#' # multiple genes
#' gene_name_prefixes <- c("HLA", "RPS", "RPL", "IG")
#' seuratObj <- del_genes(seuratObj, gene_name_prefixes, assay = "SCT")

# There's an easy way to do this. for example:
# DefaultAssay(integration.combined) <- "RNA"
# integration.combined <- subset(integration.combined, 
#                                  features = rownames(integration.combined)[!grepl("^(RPS|IG|HLA|MT-|RPL)", rownames(integration.combined))])
# subsetting will apply to all assays, so if you run the subset on the largest expression matrix it will take care of the smaller assays like integrated too

del_genes <- function(seuratObj, gene_name_prefix, assay = "RNA") {
  counts <- GetAssayData(seuratObj, assay = seuratObj@active.assay)
  counts <-
    counts[-(which(rownames(counts) %in% gene_name_prefix)),]
  seuratObj <- subset(seuratObj, features = rownames(counts))
  return(seuratObj)
}


#' Fix h5ad to work with cellxgene
#'
#' For some reason, Seurat objects which have been converted to h5ad using SeuratDisk have some sort of indexing issue
#' which prevent them from being used with cellxgene or being rewritten after importing into Python. This function
#' uses Reticulate to call a few lines in Python to fix this issue. Learn more [here](https://github.com/theislab/scvelo/issues/255#issuecomment-739995301).
#'
#' @param file Path to h5ad converted from h5Seurat by SeuratDisk
#' @param file_out Output path, overwrites if blank
#'
#' @return In-place function, no return
#' @export
#'
#' @examples
#' SaveH5Seurat(
#'   samples.combined,
#'   filename = "processed/02-annotation/cellxgene_out_UNFILTERED_no4k_doublets.h5Seurat",
#'   verbose = TRUE,
#'   overwrite = TRUE
#' )
#'
#' Convert(
#'   source = "processed/02-annotation/cellxgene_out_UNFILTERED_no4k_doublets.h5Seurat",
#'   dest = "h5ad",
#'   assay = "combatBatch",
#'   verbose = TRUE,
#'   overwrite = TRUE
#' )
#'
#' fix_cellxgene(file = "processed/02-annotation/cellxgene_out_UNFILTERED_no4k_doublets.h5ad")
fix_cellxgene_py <- function (file, file_out = NULL)
{
  reticulate::source_python(system.file("python", "fix_cellxgene.py", package = "cscb.tools"))
  if (is.null(file_out))
  {
    file_out = file
  }
  fix_cellxgene(file = file, file_out = file_out)
  message(paste0(
    file,
    " has been re-written to ",
    file_out,
    " with indexing issues fixed."
  ))
}


#' Remove Bottom or Top n Cells by nFeature_RNA or nCount_RNA
#'
#' @param seuratObj A Seurat object
#' @param n Number of desired cells to remove
#' @param nFeature Whether to sort by nFeature or nCount
#' @param top Top n cells by feature or bottom
#'
#' @return A Seurat object sorted by feature and with n cells removed
#' @export
#'
#' @examples
#' # to remove the bottom 4k cells by nFeature_RNA:
#' seuratObj <- remove_n_cells(seuratObj, 4000, top = FALSE)
remove_n_cells <-
  function(seuratObj,
           n,
           nFeature = TRUE,
           top = TRUE) {
    if (nFeature & top) {
      seuratObj@meta.data <-
        seuratObj@meta.data[order(seuratObj$nFeature_RNA, decreasing = TRUE) ,]
      
      seuratObj <-
        seuratObj[, order(seuratObj$nFeature_RNA)]
      
      cells_to_remove <- colnames(seuratObj)[c(1:n)]
      
      seuratObj <-
        seuratObj[,!colnames(seuratObj) %in% cells_to_remove]
    } else if (nFeature & top == FALSE) {
      seuratObj@meta.data <-
        seuratObj@meta.data[order(seuratObj$nFeature_RNA, decreasing = FALSE) ,]
      
      seuratObj <-
        seuratObj[, order(seuratObj$nFeature_RNA)]
      
      cells_to_remove <- colnames(seuratObj)[c(1:n)]
      
      seuratObj <-
        seuratObj[,!colnames(seuratObj) %in% cells_to_remove]
    } else if (nFeature == FALSE & top) {
      seuratObj@meta.data <-
        seuratObj@meta.data[order(seuratObj$nCount_RNA, decreasing = TRUE) ,]
      
      seuratObj <-
        seuratObj[, order(seuratObj$nCount_RNA)]
      
      cells_to_remove <- colnames(seuratObj)[c(1:n)]
      
      seuratObj <-
        seuratObj[,!colnames(seuratObj) %in% cells_to_remove]
    } else if (nFeature == FALSE & top == FALSE) {
      seuratObj@meta.data <-
        seuratObj@meta.data[order(seuratObj$nCount_RNA, decreasing = FALSE) ,]
      
      seuratObj <-
        seuratObj[, order(seuratObj$nCount_RNA)]
      
      cells_to_remove <- colnames(seuratObj)[c(1:n)]
      
      seuratObj <-
        seuratObj[,!colnames(seuratObj) %in% cells_to_remove]
    }
    return(seuratObj)
  }


#' Prepare ranked gene list for GSEA
#'
#' @param DESeq_res DESeq2 Results
#' @param gene_counts Raw counts file
#' @param rank Rank by either "log2FC" or "padj"
#'
#' @return A named vector with gene names and either log2FC or padj in descending order
#' @export
#'
#' @examples
#' geneList_res_WTvKO_cond1_input <- prepGeneList(res_WTvKO_cond1, gene_counts)
#'
#' gsea_gene_WTvKO_cond1_sig_genes_tx <- gseGO(
#'    geneList     = geneList_res_WTvKO_cond1_input,
#'    OrgDb        = org.Mm.eg.db,
#'    keyType = "ENSEMBL",
#'    pAdjustMethod = "BH",
#'    ont = "BP",
#'    minGSSize    = 100,
#'    maxGSSize    = 500,
#'    pvalueCutoff = 0.10,
#'    verbose      = TRUE
#' )
prepGeneList <- function(DESeq_res, gene_counts, rank = "log2") {
  geneList <-
    DESeq_res %>% as.data.frame() %>% rownames_to_column("gene_id")
  
  # geneList_res_WTvKO_cond2 <- geneList_res_WTvKO_cond2[geneList_res_WTvKO_cond2$gene_id %in% gene_counts$gene_name,]
  
  gene_count_sub <-
    gene_counts[gene_counts$gene_name %in% geneList$gene_id, ][, c("gene_name", "gene_id")] %>% as.data.frame()
  
  geneList <-
    merge(gene_count_sub,
          geneList,
          by.x = "gene_name",
          by.y = "gene_id")
  
  if (rank == "padj") {
    geneList <-
      geneList[, c("gene_id", "padj")]
    geneList_input <-
      geneList$padj
  } else if (rank == "log2") {
    geneList <-
      geneList[, c("gene_id", "log2FoldChange")]
    geneList_input <-
      geneList$log2FoldChange
  }
  
  names(geneList_input) <-
    geneList$gene_id
  
  geneList_input <-
    geneList_input[order(geneList_input, decreasing = TRUE)]
  
  return(geneList_input)
}

#' Writes Seurat object to file
#'
#' A "manual" solution to writing a Seurat object to file
#'
#' @param seuratObj a Seurat object with UMAP coordinates
#' @param assay_write Assay to write to file, default RNA
#' @param dir_out output directory, default cwd
#'
#' @return Nothing, writes to file
#' @export
#'
#' @examples
#' writeSeurat(gex_obj[[i]], assay = "RNA")
writeSeurat <- function(seuratObj,
                        assay = "RNA",
                        dir_out = "./") {
  seuratObj$barcode <- colnames(seuratObj)
  seuratObj$UMAP_1 <- seuratObj@reductions$umap@cell.embeddings[, 1]
  seuratObj$UMAP_2 <- seuratObj@reductions$umap@cell.embeddings[, 2]
  
  write.csv(
    seuratObj@meta.data,
    file = paste0(
      dir_out,
      "/cellxgene_metadata.csv",
      quote = F,
      row.names = F
    )
  )
  
  cmatrix <-
    Seurat::GetAssayData(seuratObj, assay = assay, slot = "counts")
  
  writeMM(cmatrix,
          file = paste0(dir_out, "/counts.mtx"))
  
  # write dimesnionality reduction matrix, in this example case pca matrix
  write.csv(
    seuratObj@reductions$pca@cell.embeddings,
    file = paste0(dir_out, "/pca.csv",
                  quote =
                    F,
                  row.names = F)
  )
  
  # write gene names
  write.table(
    data.frame('gene' = rownames(cmatrix)),
    file = paste0(dir_out, "/gene_names.csv"),
    quote = F,
    row.names = F,
    col.names = F
  )
}
