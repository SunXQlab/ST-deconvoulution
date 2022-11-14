cell2loc_process <- function(sc_data,st_data){
  
  output_path <- getwd()
  
  # single cell annoation information
  
  Cell_class <- sc_data@meta.data$Cell_class 
  names(Cell_class) <- rownames(sc_data@meta.data) 
  sc_meta <- as.data.frame(Cell_class)
  
  write.table(sc_meta, file=paste0(output_path, '/single_cell_metadata.csv'),
              sep = ",",row.names = TRUE,col.names = TRUE)
  
  # single cell count data
  sc_count_raw <- sc_data@assays$RNA@counts %>% as.matrix(.)
  genes_0_sc <- which(! rowSums(as.matrix(sc_count_raw) == 0) == ncol(sc_count_raw))
  sc_count <- sc_count_raw[genes_0_sc, ] 
  keep_sc <- which(Matrix::rowSums(sc_count > 0) >= round(0.05 * ncol(sc_count)))
  sc_count <- sc_count[keep_sc,] 
  sc_count <- t(as.matrix(sc_count)) %>% as.data.frame(.)
  sc_count <- round(sc_count)
  dim(sc_count)
  
  write.table(sc_count,file=paste0(output_path, '/single_cell_countdata.csv'),
              sep = ",",row.names = TRUE,col.names = TRUE)
  
  # spatial count data
  st_count_raw <- st_data@assays$Spatial@counts %>% as.matrix()
  colnames(st_count_raw ) <- gsub(",","_",colnames(st_count_raw ))
  genes_0_st <- which(! rowSums(as.matrix(st_count_raw) == 0) == ncol(st_count_raw))
  st_count <- st_count_raw[genes_0_st, ] 
  keep_st <- which(Matrix::rowSums(st_count > 0) >= round(0.05 * ncol(st_count)))
  st_count <- st_count[keep_st,] 
  st_count <- t(as.matrix(st_count)) %>% as.data.frame(.)
  st_count <- round(st_count)
  dim(st_count)
  
  write.table(st_count,file=paste0(output_path, '/spatial_countdata.csv'),
              sep = ",",row.names = TRUE,col.names = TRUE)
  
  sc_meta_path = paste0(output_path, '/single_cell_metadata.csv')
  sc_count_path = paste0(output_path, '/single_cell_countdata.csv')
  st_count_path = paste0(output_path, '/spatial_countdata.csv')
  
  path_ls = list(sc_meta_path,sc_count_path,st_count_path)
  return(path_ls)
}
