Seurat_pipeline <- function(sc_data, st_data){
  #Seurat deconvolution
  
  library(Seurat) 
  library(dplyr)
  library(slam)
  library(Matrix)
  library(ggsci)
  library(ggpubr)
  library(ggplot2)
  
  #output_path <- getwd()
  #load scRNA-seq data and ST data
  ptm <-proc.time()
  
  #scRNA-seq data pre-processing
  sc_count <- sc_data@assays$RNA@counts
  
  genes_0_sc <- which(! rowSums(as.matrix(sc_count) == 0) == ncol(sc_count))
  counts_sc <- sc_count[genes_0_sc, ] 
  keep_sc <- which(Matrix::rowSums(counts_sc > 0) >= round(0.05 * ncol(counts_sc)))
  counts_sc <- counts_sc[keep_sc,] 
  
  se_sc <- CreateSeuratObject(counts = counts_sc,
                              meta.data = sc_data@meta.data,
                              assay = "RNA")
  se_sc <- SCTransform(se_sc, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>% 
    RunUMAP(dims = 1:30)
  
  #ST data pre-processing
  counts_st <- st_data@assays$Spatial@counts
  spatial_loc <- st_data@meta.data
  
  #data pre-processing
  genes_0_st <- which(! rowSums(as.matrix(counts_st) == 0) == ncol(counts_st))
  counts_st <- counts_st[genes_0_st, ] 
  keep_st <- which(Matrix::rowSums(counts_st > 0) >= round(0.05 * ncol(counts_st)))
  counts_st <- counts_st[keep_st,] 
  counts_st <- as.matrix(counts_st) 
  counts_st <- as.data.frame(counts_st)
  
  se_st <- CreateSeuratObject(counts = counts_st,
                              meta.data = spatial_loc,
                              assay = "Spatial")
  se_st <- SCTransform(se_st, assay = "Spatial",  verbose = FALSE)
  se_st <- FindVariableFeatures(se_st, verbose = FALSE)
  se_st <- RunPCA(se_st, assay = "SCT", verbose = FALSE)
  se_st <- FindNeighbors(se_st, reduction = "pca", dims = 1:30)
  se_st <- FindClusters(se_st, verbose = FALSE)
  se_st <- RunUMAP(se_st, reduction = "pca", dims = 1:30)
  
  anchors <- FindTransferAnchors(reference = se_sc, query = se_st, 
                                 normalization.method = "SCT")
  
  predictions.assay <- TransferData(anchorset = anchors, 
                                    refdata = se_st$Cell_class, 
                                    prediction.assay = TRUE, 
                                    weight.reduction = merfish_st[["pca"]],
                                    dims = 1:30)
  se_st[["predictions"]] <- predictions.assay
  
  decon_mtrx <- se_st@assays[["predictions"]]@data
  decon_mtrx <- t(decon_mtrx[-10, ])
  
  # write.table(decon_mtrx, 
  #           paste0(output_path, '/Seurat_result.csv'),
  #           row.names = FALSE, col.names = TRUE, sep=",")
  
  runtime <- proc.time() - ptm
  memory <- mem_used()
  
  cat(paste0("Running time is:", round(runtime[3],3), "min"),
      paste0("Memory usage is:", round(memory/1024/1024/1024,3), "GB"),
      sep = "\n")
  return(decon_mtrx)
}



