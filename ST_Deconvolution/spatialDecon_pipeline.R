#spatialDecon deconvolution

library(SpatialDecon)
library(Seurat)
library(dplyr)
library(pryr)

source(file="D:/AA-luluyan-phd/code/ST deconvoluting method/Deconvolution methods/Regression-based/SpatialDecon-master/R/create_profile_matrix.R")
source(file="D:/AA-luluyan-phd/code/ST deconvoluting method/Deconvolution methods/Regression-based/SpatialDecon-master/R/runspatialdecon.R")
source(file="D:/AA-luluyan-phd/code/ST deconvoluting method/Deconvolution methods/Regression-based/SpatialDecon-master/R/runCollapseCellTypes.R")

output_path <- getwd()
#load data
ptm <-proc.time()
load("/home/yll/spatial_decon_methods/STdeconvolve_ReferenceFree/MERFISH_singlecell_dataset_new.RData")
load("~/spatial_decon_methods/spatialDWLS_dataset-main/codes_yll/MERFISH_spatialspot_dataset.RData")

#scRNA-seq data pre-processing
sc_meta <- annot.table_new
filter_count_matrix <- round(t(as.matrix(counts_new)))

genes_0_sc <- which(! rowSums(as.matrix(filter_count_matrix) == 0) == ncol(filter_count_matrix))
counts_sc <- filter_count_matrix[genes_0_sc, ] 
keep_sc <- which(Matrix::rowSums(counts_sc > 0) >= round(0.05 * ncol(counts_sc)))
counts_sc <- counts_sc[keep_sc,] 

ref_seur_obj <- CreateSeuratObject(counts = counts_sc,
                                   meta.data = sc_meta)

#normalizetion
ref_seur_obj <- Seurat::SCTransform(ref_seur_obj, verbose=FALSE) 

ref_cell_type <- data.frame(Cell_class = sc_meta$Cell_class,Neuron_cluster_ID = rownames(counts_sc))
#names(ref_cell_type) <- colnames(ref_seur_obj@assays$SCT@counts)
temp_custom_ref <- create_profile_matrix(mtx = ref_seur_obj@assays$SCT@counts,
                                         cellAnnots = ref_cell_type,
                                         cellTypeCol = "Cell_class",
                                         cellNameCol = "Neuron_cluster_ID",
                                         minGenes = 0)

#ST data pre-processing
counts_st <- t(as.matrix(simFN7$sim))
counts_st <- round(counts_st) 
spatial_loc <- simFN7$cellCounts 

genes_0_st <- which(! rowSums(as.matrix(counts_st) == 0) == ncol(counts_st))
counts_st <- counts_st[genes_0_st, ] 
keep_st <- which(rowSums(counts_st > 0) >= round(0.05 * ncol(counts_st)))
counts_st <- counts_st[keep_st,] 
counts_st <- as.matrix(counts_st) 

spatial_metadata <- mutate(spatial_loc,spot_ID = colnames(counts_st))
rownames(spatial_metadata) <- colnames(counts_st)
st_seur_obj <- CreateSeuratObject(counts = counts_st, 
                                  meta.data = spatial_metadata,
                                  assay = "Spatial")

temp_spatialdecon_result <- runspatialdecon(object = st_seur_obj, # Seurat object
                                            X = temp_custom_ref, # safeTME matrix, used by default
                                            bg = 0.01)

temp_beta <- t(temp_spatialdecon_result$beta)
observed_beta <- t(apply(temp_beta,1,function(x) x/sum(x))) #proportion matrix

write.csv(observed_beta, 
          paste0(output_path, '/spatialDecon_result.csv'),
          row.names = FALSE, col.names = TRUE, sep=",")

runtime <- proc.time() - ptm
memory <- mem_used()

cat(paste0("Running time is:", round(runtime[3],3), "min"),
    paste0("Memory usage is:", round(memory/1024/1024/1024,3), "GB"),
    sep = "\n")



