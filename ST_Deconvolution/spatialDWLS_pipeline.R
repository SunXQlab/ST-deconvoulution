spatialDWLS_pipeline <- function(sc_obj, st_obj){
  
  #spatialDWLS deconvolution
  
  library(Giotto)
  library(ggplot2)
  library(scatterpie)
  library(data.table)
  library(slam)
  library(pryr)
  library(dplyr)
  library(tibble)
  library(philentropy)
  
  ptm <- proc.time()
  
  #data pre-processing
  sc_meta <- sc_obj@meta.data
  sc_count  <- sc_obj@assays$RNA@counts
  
  genes_0_sc <- which(! rowSums(as.matrix(sc_count ) == 0) == ncol(sc_count ))
  count_sc <- sc_count[genes_0_sc, ] 
  keep_sc <- which(Matrix::rowSums(count_sc > 0) >= round(0.05 * ncol(count_sc)))
  count_sc <- count_sc[keep_sc,] 
  count_sc <- as.data.frame(count_sc)
  
  #Identification of marker gene expression in single cell
  my_python_path = "/usr/bin/python3"
  instrs = createGiottoInstructions(python_path = my_python_path)
  sc_data <- createGiottoObject(raw_exprs = count_sc,instructions = instrs)
  sc_data <- normalizeGiotto(gobject = sc_data, scalefactor = 6000, verbose = T)
  sc_data <- calculateHVG(gobject = sc_data)
  gene_metadata = fDataDT(sc_data) 
  sc_featgenes = gene_metadata[hvg == 'yes']$gene_ID 
  
  sc_data <- runPCA(gobject = sc_data, genes_to_use = sc_featgenes, scale_unit = F)
  
  #calculate Sig for deconvolution, This step use DEG function implemented in Giotto
  sc_data@cell_metadata$leiden_clus <- sc_meta$Cell_class 
  scran_markers_subclusters = findMarkers_one_vs_all(gobject = sc_data,
                                                     method = 'scran',
                                                     expression_values = 'normalized',
                                                     cluster_column = 'leiden_clus')
  Sig_scran <- unique(scran_markers_subclusters$genes[which(scran_markers_subclusters$ranking <= 100)])#1155 signature gene
  
  #Calculate median expression value of signature genes in each cell type
  norm_exp <- 2^(sc_data@norm_expr)-1
  id <- sc_data@cell_metadata$leiden_clus 
  ExprSubset <- norm_exp[Sig_scran,] 
  Sig_exp <- NULL 
  for (i in unique(id)){
    Sig_exp <- cbind(Sig_exp,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(Sig_exp) <- unique(id) 
  
  #Spatial transcriptomic data analysis
  count_st <- st_obj@assays$Spatial@counts
  spatial_loc <- st_obj@meta.data 
  
  #data pre-processing
  genes_0_st <- which(! rowSums(as.matrix(count_st) == 0) == ncol(count_st))
  count_st <- count_st[genes_0_st, ] 
  keep_st <- which(rowSums(count_st > 0) >= round(0.05 * ncol(count_st)))
  count_st <- count_st[keep_st,] 
  count_st <- as.data.frame(count_st)
  
  ##Generate Giotto objects and cluster spots
  st_data <- createGiottoObject(raw_exprs = count_st,
                                spatial_locs = spatial_loc,
                                instructions = instrs)
  
  #st_data <- filterGiotto(gobject = st_data,
  #                        expression_threshold = 1,
  #                        gene_det_in_min_cells = 10,
  #                        min_det_genes_per_cell = 10,
  #                        expression_values = c('raw'),
  #                        verbose = T) 
  
  st_data <- normalizeGiotto(gobject = st_data)
  st_data <- calculateHVG(gobject = st_data)
  gene_metadata = fDataDT(st_data) 
  st_featgenes = gene_metadata[hvg == 'yes']$gene_ID 
  
  st_data <- runPCA(gobject = st_data, genes_to_use = st_featgenes, scale_unit = F)
  st_data <- createNearestNetwork(gobject = st_data, dimensions_to_use = 1:10, k = 10)
  st_data <- doLeidenCluster(gobject = st_data, resolution = 0.4, n_iterations = 1000)
  st_data <- runDWLSDeconv(gobject = st_data, sign_matrix = Sig_exp)
  
  decon_mtrx <- st_data@spatial_enrichment$DWLS
  decon_mtrx <- column_to_rownames(decon_mtrx,var="cell_ID")
  decon_mtrx <- as.matrix(decon_mtrx)
  
  # write.table(st_data@spatial_enrichment$DWLS, 
  #           paste0(output_path, '/SpatialDWLS_result.csv'),
  #           row.names = FALSE, col.names = TRUE, sep=",")
  
  runtime <- (proc.time() - ptm)/60
  memory <- mem_used()
  
  cat(paste0("Running time is:", round(runtime[3],3), "min"),
      paste0("Memory usage is:", round(memory/1024/1024/1024,3), "GB"),
      sep = "\n")
  return(decon_mtrx)
}





