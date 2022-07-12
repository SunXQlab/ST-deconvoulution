#cell2location data process

library(readr)

output_path <- getwd()

## MPOA (MERFISH) data process for cell2location
load("/home/yll/spatial_decon_methods/STdeconvolve_ReferenceFree/MERFISH_singlecell_dataset_new.RData")
load("~/spatial_decon_methods/spatialDWLS_dataset-main/codes_yll/MERFISH_spatialspot_dataset.RData")

#single cell annoation information
cell_metadata <- simFN7$annotDf
cell_metadata <- cell_metadata[,-6]

write.table(cell_metadata, file=paste0(output_path, '/merfish_cell_metadata.csv'),
            sep = ",",row.names = TRUE,col.names = TRUE)

#single cell count data
filter_count_matrix <- counts_new[rownames(cell_metadata), ]
filter_count_matrix <- as.data.frame(filter_count_matrix)
filter_count_matrix <- round(filter_count_matrix)
dim(filter_count_matrix)
# genes_0_sc <- which(! rowSums(as.matrix(filter_count_matrix) == 0) == ncol(filter_count_matrix))
# counts_sc <- filter_count_matrix[genes_0_sc, ] 
# keep_sc <- which(Matrix::rowSums(counts_sc > 0) >= round(0.05 * ncol(counts_sc)))
# counts_sc <- counts_sc[keep_sc,]

write.table(filter_count_matrix,file=paste0(output_path, '/merfish_cell_countdata.csv'),
            sep = ",",row.names = TRUE,col.names = TRUE)

#spatial count data
counts_st <- as.matrix(simFN7$sim)
counts_st <- as.data.frame(counts_st) #3072*135
counts_st <- round(counts_st)
dim(counts_st)

write.table(counts_st,file=paste0(output_path, '/merfish_spatial_countdata.csv'),
            sep = ",",row.names = TRUE,col.names = TRUE)

## embryo (sci-Space) data process for cell2location
load("~/spatial_decon_methods/spatialDWLS_dataset-main/codes_yll/embryo_singlecell_remove_rare_ct.RData")
load("~/spatial_decon_methods/spatialDWLS_dataset-main/codes_yll/embryo_spatial_remove_rare_ct.RData")
#single cell annoation information
sc_meta <- data.frame(cluster=filter_cell_metadata$cluster,
                    celltype=filter_cell_metadata$final_cluster_label)
write.table(sc_meta, file=paste0(output_path, '/embryo_cell_metadata.csv'),
            sep = ",",row.names = TRUE,col.names = TRUE)

#data pre-processing
genes_0_sc <- which(! rowSums(as.matrix(filter_count_matrix) == 0) == ncol(filter_count_matrix))
counts_sc <- filter_count_matrix[genes_0_sc, ] 
keep_sc <- which(Matrix::rowSums(counts_sc > 0) >= round(0.05 * ncol(counts_sc)))
counts_sc <- counts_sc[keep_sc,] 
counts_sc <- t(as.matrix(counts_sc))
counts_sc <- as.data.frame(counts_sc)
dim(counts_sc)

write.table(counts_st,file=paste0(output_path, '/embryo_cell_countdata.csv'),
            sep = ",",row.names = TRUE,col.names = TRUE)

#spatial count data
colnames(slide14_spatial_count) <- gsub(",","_",colnames(slide14_spatial_count))
genes_0_st <- which(! rowSums(as.matrix(slide14_spatial_count) == 0) == ncol(slide14_spatial_count))
counts_st <- slide14_spatial_count[genes_0_st, ] 
keep_st <- which(Matrix::rowSums(counts_st > 0) >= round(0.05 * ncol(counts_st)))
counts_st <- counts_st[keep_st,] 
counts_st <- t(as.matrix(counts_st))
counts_st <- as.data.frame(counts_st)
dim(counts_st)

write.table(counts_st,file=paste0(output_path,'/embryo_spatial_countdata.csv'),
            sep = ",",row.names = TRUE,col.names = TRUE)
            
## brain (cellTrek) data process for cell2location
load("~/scRNA-seq_Mapping_to_STdata/brain_celltrek_usedcells_data.RData")
load("~/SequencingDepth/brain_spatial_SeqDepth0.01_data.RData")

#single cell annoation information
brain_celltrek <- sim_brain_celltrek
cell_metadata <- brain_celltrek@meta.data
write.table(cell_metadata,file=paste0(output_path,'/brain_cell_metadata.csv'),
            sep = ",",row.names = TRUE,col.names = TRUE)

#data pre-processing
sc_count <- as.matrix(brain_celltrek@assays$RNA@counts)
genes_0_sc <- which(! rowSums(sc_count == 0) == ncol(sc_count))
sc_count <- sc_count[genes_0_sc, ] 
keep_sc <- which(Matrix::rowSums(sc_count > 0) >= round(0.05 * ncol(sc_count)))
sc_count <- sc_count[keep_sc,] 
sc_count <- as.data.frame(t(sc_count))
dim(sc_count)

write.table(sc_count ,file=paste0(output_path,'/brain_cell_countdata.csv'),
            sep = ",",row.names = TRUE,col.names = TRUE)

#spatial data pre-processing
st_count <- as.matrix(sim_brain_st@assays$spatial@counts)
genes_0_st <- which(! rowSums(st_count == 0) == ncol(st_count))
st_count <- st_count[genes_0_st, ] 
keep_st <- which(Matrix::rowSums(st_count > 0) >= round(0.05 * ncol(st_count)))
st_count <- st_count[keep_st,] 
st_count <- as.data.frame(t(st_count))
dim(st_count)

write.table(spatial_loc,file=paste0(output_path,'/brain_spatial_seq0.01_countdata.csv'),
            sep = ",",row.names = TRUE,col.names = TRUE)
