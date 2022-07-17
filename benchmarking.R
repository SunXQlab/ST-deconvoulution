benchmark <- function (sc_data, st_data, number_ct, methods, ground_truth){
  
  
  # Check variables
  if (!is.character(methods)) stop("ERROR: method must be a character!")
  
  dir=paste0(getwd(),'/Deconvolution/')
  list.files(dir,pattern = ".R$")
  for (i in list.files(path=dir, pattern = ".R$")){
    source(paste0(getwd(),'/Deconvolution/R_deconvolution/',i))
  }
  
  source(paste0(getwd(),'/Evaluation/benchmark_performance_fun.R'))
  
  WorkDir <- paste0("./deconvolution_results/decon_", methods)
  dir.create(WorkDir, recursive = TRUE, showWarnings = F)
  cat(paste0("WorkDir: ", WorkDir, "\n"))
  
  if (methods == "RCTD"){
    decon_result <- RCTD_pipeline(sc_data = sc_data,
                                  st_data = st_data)
    write.table(decon_result, 
                paste0(WorkDir, '/decon_result.csv'),
                row.names = TRUE, col.names = TRUE, sep=",")
    
    result <- benchmark_performance(ground_truth_mtrx = ground_truth,
                                    deconvoluted_mtrx = decon_result)
    
  }else if (methods == "spatialDWLS"){
    decon_result <- spatialDWLS_pipeline(sc_obj = sc_data,
                                         st_obj = st_data)
    
    write.table(decon_result, 
                paste0(WorkDir, '/decon_result.csv'),
                row.names = TRUE, col.names = TRUE, sep=",")
    
    result <- benchmark_performance(ground_truth_mtrx = ground_truth,
                                    deconvoluted_mtrx = decon_result)
    
  }else if (methods == "Seurat"){
    decon_result <- Seurat_pipeline(sc_data = sc_data,
                                    st_data = st_data)
    
    write.table(decon_result, 
                paste0(WorkDir, '/decon_result.csv'),
                row.names = TRUE, col.names = TRUE, sep=",")
    
    result <- benchmark_performance(ground_truth_mtrx = ground_truth,
                                    deconvoluted_mtrx = decon_result)
    
  }else if (methods == "PAGE"){
    decon_result <- Giotto_enrichment_pipeline(sc_obj = sc_data,
                                    st_obj = st_data,
                                    enrich_method = methods)
    write.table(decon_result, 
                paste0(WorkDir, '/decon_result.csv'),
                row.names = FALSE, col.names = TRUE, sep=",")
    
    result <- benchmark_performance(ground_truth_mtrx = ground_truth,
                                    deconvoluted_mtrx = decon_result)
    
  }else if (methods == "STdeconvolve"){
    decon_result <- STdeconvolve_pipeline(sc_data = sc_data,
                                          st_data = st_data,
                                          num_ct = number_ct)
    
    write.table(decon_result, 
                paste0(WorkDir, '/decon_result.csv'),
                row.names = TRUE, col.names = TRUE, sep=",")
    # result <- benchmark_performance(ground_truth_mtrx = ground_truth,
    #                                 deconvoluted_mtrx = decon_result)
    
  }else if (methods == "SPOTlight"){
    decon_result <- SPOTlight_pipeline(sc_obj = sc_data,
                                       st_obj = st_data)
    
    write.table(decon_result, 
                paste0(WorkDir, '/decon_result.csv'),
                row.names = TRUE, col.names = TRUE, sep=",")
    
    result <- benchmark_performance(ground_truth_mtrx = ground_truth,
                                    deconvoluted_mtrx = decon_result)
    
  }else if (methods == "rank"){
    decon_result <- Giotto_enrichment_pipeline(sc_obj = sc_data,
                                               st_obj = st_data,
                                               enrich_method = methods)
    write.table(decon_result, 
                paste0(WorkDir, '/decon_result.csv'),
                row.names = TRUE, col.names = TRUE, sep=",")
    
    result <- benchmark_performance(ground_truth_mtrx = ground_truth,
                                    deconvoluted_mtrx = decon_result)
    
  }else if (methods == "hypergeometric"){
    decon_result <- Giotto_enrichment_pipeline(sc_obj = sc_data,
                                               st_obj = st_data,
                                               enrich_method = methods)
    write.table(decon_result, 
                paste0(WorkDir, '/decon_result.csv'),
                row.names = TRUE, col.names = TRUE, sep=",")
    
    result <- benchmark_performance(ground_truth_mtrx = ground_truth,
                                    deconvoluted_mtrx = decon_result)
  }
  
  return(list(decon_result = decon_result, bench_result = result))
}



