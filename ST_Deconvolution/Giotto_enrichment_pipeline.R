Giotto_enrichment_pipeline <- function(sc_obj, st_obj, enrich_method) {
  #Giotto/Hypergeometric/RANK/PAGE enrichment
  #Note that the deconvoluted result of enrich_method="hypergeometric" on the MPOA (MERFISH) dataset
  #contains 'inf'.
  library(Giotto)
  library(ggplot2)
  library(scatterpie)
  library(data.table)
  library(slam)
  library(tibble)
  library(pryr)
  
  #output_path <- getwd()
  #load scRNA-seq data and ST data
  ptm <-proc.time()
  
  #scRNA-seq data pre-processing
  #sc_meta <- annot.table_new
  #filter_count_matrix <- t(counts_new[rownames(sc_meta), ])
  sc_count <- sc_obj@assays$RNA@counts
  
  genes_0_sc <- which(! rowSums(as.matrix(sc_count) == 0) == ncol(sc_count))
  counts_sc <- sc_count[genes_0_sc, ] 
  keep_sc <- which(Matrix::rowSums(counts_sc > 0) >= round(0.05 * ncol(counts_sc)))
  counts_sc <- counts_sc[keep_sc,]
  
  my_python_path = "/usr/bin/python3"
  instrs = createGiottoInstructions(python_path = my_python_path)
  sc_data <- createGiottoObject(raw_exprs = counts_sc,instructions = instrs)
  sc_data <- normalizeGiotto(gobject = sc_data, scalefactor = 6000, verbose = T)
  sc_data <- calculateHVG(gobject = sc_data)
  gene_metadata <- fDataDT(sc_data) 
  sc_featgenes <- gene_metadata[hvg == 'yes']$gene_ID 
  sc_data <- runPCA(gobject = sc_data, genes_to_use = sc_featgenes, scale_unit = F)
  
  #calculate Sig for deconvolution, This step use DEG function implemented in Giotto
  sc_data@cell_metadata$leiden_clus <- sc_obj@meta.data$Cell_class 
  scran_markers_subclusters = findMarkers_one_vs_all(gobject = sc_data,
                                                     method = 'scran',
                                                     expression_values = 'normalized',
                                                     cluster_column = 'leiden_clus')
  Sig_scran <- unique(scran_markers_subclusters$genes[which(scran_markers_subclusters$ranking <= 100)])
  
  #Calculate median expression value of signature genes in each cell type
  norm_exp <- 2^(sc_data@norm_expr)-1
  id <- sc_data@cell_metadata$leiden_clus 
  ExprSubset <- norm_exp[Sig_scran,]
  Sig_exp <- NULL 
  for (i in unique(id)){
    Sig_exp <- cbind(Sig_exp,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(Sig_exp) <- unique(id) 
  
  #ST data pre-processing
  counts_st <- st_obj@assays$Spatial@counts
  spatial_loc <- st_obj@meta.data 
  
  genes_0_st <- which(! rowSums(as.matrix(counts_st) == 0) == ncol(counts_st))
  counts_st <- counts_st[genes_0_st, ] 
  keep_st <- which(Matrix::rowSums(counts_st > 0) >= round(0.05 * ncol(counts_st)))
  counts_st <- counts_st[keep_st,] 
  counts_st <- as.matrix(counts_st) 
  counts_st <- as.data.frame(counts_st)
  
  #Generate Giotto objects and cluster spots
  st_data <- createGiottoObject(raw_exprs = counts_st,spatial_locs = spatial_loc,
                                instructions = instrs)
  
  st_data <- normalizeGiotto(gobject = st_data)
  st_data <- calculateHVG(gobject = st_data)
  gene_metadata = fDataDT(st_data) 
  st_featgenes = gene_metadata[hvg == 'yes']$gene_ID 
  st_data <- runPCA(gobject = st_data, genes_to_use = st_featgenes, scale_unit = F)
  st_data <- createNearestNetwork(gobject = st_data, dimensions_to_use = 1:10, k = 10)
  st_data <- doLeidenCluster(gobject = st_data, resolution = 0.4, n_iterations = 1000)
  
  #signature matrix for PAGE and Hypergeomtric
  logbase=2
  cluster_column="leiden_clus"
  expr_values = st_data@norm_expr
  nolog_expr = logbase^(expr_values) - 1
  cell_metadata = pDataDT(st_data)
  if (!cluster_column %in% colnames(cell_metadata)) {
    stop("\n cluster column not found \n")
  }
  cluster = cell_metadata[[cluster_column]]
  Sig_exp = as.matrix(Sig_exp)
  intersect_gene = intersect(rownames(Sig_exp), rownames(nolog_expr))
  filter_Sig = Sig_exp[intersect_gene, ]
  filter_Sig = filter_Sig[rowSums(filter_Sig) > 0, ]
  enrich_matrix = matrix(0, nrow = dim(filter_Sig)[1], ncol = dim(filter_Sig)[2])
  rowmax_col = Rfast::rowMaxs(filter_Sig)
  for (i in 1:length(rowmax_col)) {
    enrich_matrix[i, rowmax_col[i]] = 1
  }
  colsum_ct_binary <- colSums(enrich_matrix)
  for (i in 1:length(colsum_ct_binary)) {
    if (colsum_ct_binary[i] <= 2) {
      rank = rank(-ct_exp[, i])
      enrich_matrix[rank <= 2, i] = 1
    }
  }
  rownames(enrich_matrix) = rownames(filter_Sig)
  colnames(enrich_matrix) = colnames(filter_Sig)
  
  #rank matrix for rank 
  single_cell_matrix <- sc_data@norm_expr
  cell_metadata_sc <- pDataDT(sc_data)
  cell_annotations <- cell_metadata_sc$leiden_clus
  rank_matrix <- makeSignMatrixRank(sc_matrix = single_cell_matrix, 
                                    sc_cluster_ids = cell_annotations)
  
  #enrichment analysis
  if (!enrich_method == "rank"){
    st_data = runSpatialEnrich(gobject = st_data,
                               sign_matrix = enrich_matrix,
                               enrich_method = enrich_method,#"rank','hypergeometric'
                               top_percentage = 5,
                               min_overlap_genes = 2) #default = 'PAGE'
  }else if(enrich_method == "rank"){
    st_data = runSpatialEnrich(gobject = st_data,
                               sign_matrix = rank_matrix,
                               enrich_method = enrich_method,#"rank','hypergeometric'
                               top_percentage = 5,
                               min_overlap_genes = 2)
  }

  
  if (enrich_method =='hypergeometric'){
    enrich_result <- st_data@spatial_enrichment$hypergeometric
    enrich_result <- column_to_rownames(enrich_result,var="cell_ID")
    enrich_result <- as.matrix(enrich_result) + 1e-12
    norm_result <- t(apply(enrich_result,1,function(x) x/sum(x)))
    # write.table(norm_result,paste0(output_path, '/Giotto_hypergeometric_result.csv'),
    #           row.names = FALSE, col.names = TRUE, sep=",")
  } else if (enrich_method =='PAGE'){
    enrich_result <- st_data@spatial_enrichment$PAGE
    enrich_result <- column_to_rownames(enrich_result,var="cell_ID")
    enrich_result <- as.matrix(enrich_result)
    enrich_result <- apply(enrich_result,2,function(x) ifelse(x < 0, 0, x)) +1e-12
    norm_result <- t(apply(enrich_result,1,function(x) x/sum(x)))
    # write.table(norm_result,paste0(output_path, '/Giotto_PAGE_result.csv'),
    #           row.names = FALSE, col.names = TRUE, sep=",")
  } else if (enrich_method =='rank'){
    enrich_result <- st_data@spatial_enrichment$rank
    enrich_result <- column_to_rownames(enrich_result,var="cell_ID")
    enrich_result <- as.matrix(enrich_result)
    norm_result <- t(apply(enrich_result,1,function(x) x/sum(x)))
    # write.table(norm_result,paste0(output_path, '/Giotto_RANK_result.csv'),
    #           row.names = FALSE, col.names = TRUE, sep=",")
  }
  
  runtime <- (proc.time() - ptm)/60
  memory <- mem_used()
  
  cat(paste0("Running time is:", round(runtime[3],3), "min"),
      paste0("Memory usage is:", round(memory/1024/1024/1024,3), "GB"),
      sep = "\n")
  
  return(norm_result)
}
