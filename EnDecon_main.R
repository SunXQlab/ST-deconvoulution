EnDecon_main <- function(sc_data, st_data, python_path){
  
  output_path <- getwd()
  ## EnDecon three methods
  
  DeconResults <- list()
  
  RCTD_result <- RCTD_pipeline(sc_data = sc_data,
                               st_data = st_data)
  
  
  spatialDWLS_result <- spatialDWLS_pipeline(sc_obj = sc_data,
                                             st_obj = st_data,
                                             python_path = python_path)
  
  # run cell2location
  #data process
  
  path_ls <- cell2loc_process(sc_data,st_data)
  
  st_cnt_path = path_ls[[3]]
  sc_cnt_path = path_ls[[2]]
  sc_mta_path = path_ls[[1]]
 
  use_python(python_path)
  
  WorkDir <- paste0(output_path,"/deconvolution_results/decon_", "cell2location")
  dir.create(WorkDir, recursive = TRUE, showWarnings = F)
  cat(paste0("WorkDir: ", WorkDir, "\n"))
  means_cell_abun <- cell2location_main(sc_cnt_dir = sc_cnt_path,
                                        sc_mta_dir = sc_mta_path,
                                        st_cnt_dir = st_cnt_path)
  
  means_cell_abun <- as.matrix(means_cell_abun)
  colnames(means_cell_abun) <- gsub("meanscell_abundance_w_sf_","",colnames(means_cell_abun))
  cell2loc_result <- t(apply(means_cell_abun,1,as.numeric))
  colnames(cell2loc_result) <- colnames(means_cell_abun)
  
  write.table(cell2loc_result,
              paste0(WorkDir, '/decon_result.csv'),
              row.names = TRUE, col.names = TRUE, sep=",")
  
  cell2loc_result <- cell2loc_result[rownames(RCTD_result),]
  spatialDWLS_result <- spatialDWLS_result[rownames(RCTD_result),]
  
  DeconResults[[1]] <- cell2loc_result
  DeconResults[[2]] <- RCTD_result
  DeconResults[[3]] <- spatialDWLS_result
  
  # EnDecon model
  
  w <- data.frame(cell2location = 0.651577272, RCTD = 0.111148121, spatialDWLS = 0.237274608)
  EnDecon_result <- w[1,1]*DeconResults[[1]] + w[1,2]*DeconResults[[2]] + w[1,3]*DeconResults[[3]]
  
  WorkDir <- paste0(output_path,"/deconvolution_results/decon_", "EnDecon")
  dir.create(WorkDir, recursive = TRUE, showWarnings = F)
  cat(paste0("WorkDir: ", WorkDir, "\n"))
  write.table(EnDecon_result,
              paste0(WorkDir, '/decon_result.csv'),
              row.names = TRUE, col.names = TRUE, sep=",")
  
  DeconResults[[4]] <- EnDecon_result
  
  return(DeconResults)
}
