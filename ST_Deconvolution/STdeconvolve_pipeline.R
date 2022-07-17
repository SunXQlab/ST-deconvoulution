STdeconvolve_pipeline <- function(sc_data, st_data, num_ct){
  
  library(STdeconvolve)
  library(pryr)
  library(Matrix)
  library(slam)
  
  #load ST data
  ptm <-proc.time()
  
  counts_st <- st_data@assays$Spatial@counts
  spatial_loc <- st_data@meta.data 
  
  #data pre-processing
  genes_0_st <- which(! rowSums(as.matrix(counts_st) == 0) == ncol(counts_st))
  counts_st <- counts_st[genes_0_st, ] 
  keep_st <- which(rowSums(counts_st > 0) >= round(0.05 * ncol(counts_st)))
  counts_st <- counts_st[keep_st,] 
  counts_st <- as.matrix(counts_st)
  
  #remove pixels with too few genes
  counts_st <- cleanCounts(counts_st, min.lib.size = 100)
  counts_st <- round(counts_st)
  
  #feature select for genes，因为merfish dataset的gene太少了，所以不进行gene select
  #corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05) #232 genes
  
  #choose optimal number of cell-types
  ldas <- fitLDA(t(as.matrix(counts_st)), 
                 Ks = seq(2, 20, by = 1),
                 perc.rare.thresh = 0.05,
                 plot=TRUE,
                 verbose=TRUE)
  
  #get best model results
  optLDA <- optimalModel(models = ldas, opt = num_ct) # or extract the model for any K that was used
  
  # extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
  st_results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
  deconProp <- st_results$theta #proportion matrix: 3071 pixels * 9 celltype
  deconGexp <- st_results$beta #gene expression matrix: 9 celltypes * 135genes
  
  
  sc_meta <- sc_seurat@meta.data
  counts_sc <- t(sc_seurat@assays$RNA@counts)
  celltype <- unique(sc_meta$Cell_class)
  
  #compute ground truth gene expression from single cell count matrix
  sc_gexp <- NULL
  for (i in celltype){
    pos <- which(sc_meta$Cell_class==i)
    counts_ct <- apply(counts_sc[pos,],2,mean)
    sc_gexp <- cbind(sc_gexp,counts_ct)
  }
  colnames(sc_gexp) <- celltype #135genes * 9celltypes
  
  #strategy2: GSEA 
  Markers <- list()
  
  ## make the tissue layers the rows and genes the columns
  gexp <- t(as.matrix(sc_gexp)) #9 celltypes * 135 genes
  
  for (i in seq(length(rownames(sc_gexp)))){
    celltype <- i
    ## log2FC relative to other cell-types
    ## highly expressed in cell-type of interest
    highgexp <- names(which(sc_gexp[celltype,] > 4 ))
    ## high log2(fold-change) compared to other deconvolved cell-types and limit to top 200
    log2fc <- sort(log2(gexp[celltype,highgexp]/colMeans(gexp[-celltype,highgexp])), decreasing=TRUE)
    
    ## for gene set of the ground truth cell-type, get the genes
    ## with log2FC > 1 (so FC > 2 over the mean exp of the other cell-types)
    markers <- names(log2fc[log2fc > 1])
    Markers[[ rownames(gexp)[celltype] ]] <- markers
  }
  
  celltype_annotations <- annotateCellTypesGSEA(beta = st_results$beta, gset = Markers, qval = 0.05)
  celltype_annotations$results$`2`
  pred_ct <- celltype_annotations$predictions
  
  colnames(deconProp) <- pred_ct
  runtime <- proc.time() - ptm
  memory <- mem_used()
  
  cat(paste0("STdeconvolve's running time is:", round(runtime[3],3), "min"),
      paste0("STdeconvolve's memory usage is:", round(memory/1024/1024/1024,3), "GB"),
      sep = "\n")
  
  return(list(decon_mtrx = deconProp, pred_ct = pred_st))
}



