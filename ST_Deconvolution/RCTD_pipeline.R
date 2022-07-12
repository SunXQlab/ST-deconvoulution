#RCTD deconvolution

library(spacexr)
library(pryr)
library(slam)
library(Matrix)
library(Seurat)

output_path <- getwd()
#load scRNA-seq data and ST data
ptm <-proc.time()
load("/home/yll/spatial_decon_methods/STdeconvolve_ReferenceFree/MERFISH_singlecell_dataset_new.RData")
load("~/spatial_decon_methods/spatialDWLS_dataset-main/codes_yll/MERFISH_spatialspot_dataset.RData")

#scRNA-seq data pre-processing
cell_metadata <- annot.table_new
filter_count_matrix <- round(t(as.matrix(counts_new)))

genes_0_sc <- which(! rowSums(as.matrix(filter_count_matrix) == 0) == ncol(filter_count_matrix))
counts_sc <- filter_count_matrix[genes_0_sc, ] 
keep_sc <- which(Matrix::rowSums(counts_sc > 0) >= round(0.05 * ncol(counts_sc)))
counts_sc <- counts_sc[keep_sc,] 

cell_type <- factor(cell_metadata$Cell_class) 
names(cell_type) <- colnames(counts_sc)
reference <- Reference(counts=counts_sc, cell_types=cell_type, nUMI=NULL) 

#ST data pre-processing
counts_st <- t(as.matrix(simFN7$sim))
counts_st <- round(counts_st) 
spatial_loc <- simFN7$cellCounts 

genes_0_st <- which(! rowSums(as.matrix(counts_st) == 0) == ncol(counts_st))
counts_st <- counts_st[genes_0_st, ] 
keep_st <- which(rowSums(counts_st > 0) >= round(0.05 * ncol(counts_st)))
counts_st <- counts_st[keep_st,] 
counts_st <- as.matrix(counts_st) 

#creat spatialRNA object
spot_coords <- data.frame(x=spatial_loc$x,
                          y=spatial_loc$y,
                          row.names = colnames(counts_st))
spaceRNA <- SpatialRNA(coords=spot_coords, counts=counts_st, nUMI=NULL)

#creat RCTD object
myRCTD <- create.RCTD(spaceRNA, reference)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

decon_mtrx <- as.matrix(myRCTD@results$weights) #3071 spots
norm_mtrx <- normalize_weights(decon_mtrx)

write.csv(norm_mtrx, 
          paste0(output_path, '/RCTD_result.csv'),
          row.names = FALSE, col.names = TRUE, sep=",")

runtime <- (proc.time() - ptm)/60
memory <- mem_used()

cat(paste0("Running time is:", round(runtime[3],3), "min"),
    paste0("Memory usage is:", round(memory/1024/1024/1024,3), "GB"),
    sep = "\n")
