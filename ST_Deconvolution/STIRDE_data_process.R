#STRIDE data procss

library(slam)
library(Matrix)
library(tibble)

output_path <- getwd()

## MPOA (MERFISH) data process
load("/home/yll/spatial_decon_methods/STdeconvolve_ReferenceFree/MERFISH_singlecell_dataset_new.RData")
load("~/spatial_decon_methods/spatialDWLS_dataset-main/codes_yll/MERFISH_spatialspot_dataset.RData")

#save single cell annoation infor as .txt format
cell_metadata <- annot.table_new
cell_metadata$Cell_class <- gsub(" ", "\\.",cell_metadata$Cell_class)
cell_lineage <- data.frame(v1 = rownames(cell_metadata),v2 = cell_metadata$final_cluster_label)
cell_lineage <- data.frame(v1 = colnames(filter_count_matrix),v2 = cell_metadata$final_cluster_label)

write.table(mer_lineage,file=paste0(output_path, '/merfish_sc_lineage.txt'),
            sep = "\t",quote = FALSE, eol = "\n",row.names = FALSE,col.names = FALSE)

#save single cell count as .txt format
filter_count_matrix <- t(counts_new[rownames(cell_metadata), ])
genes_0_sc <- which(! rowSums(as.matrix(filter_count_matrix) == 0) == ncol(filter_count_matrix))
counts_sc <- filter_count_matrix[genes_0_sc, ] 
keep_sc <- which(Matrix::rowSums(counts_sc > 0) >= round(0.05 * ncol(counts_sc)))
counts_sc <- counts_sc[keep_sc,]

write.table(counts_sc,file=paste0(output_path, '/merfish_sc_gene_count.txt'),
            sep = "\t",quote = FALSE, eol = "\n", row.names = TRUE,col.names = TRUE)

#save spatial count as .txt format
counts_st <- t(as.matrix(simFN7$sim))
counts_st <- round(counts_st)
genes_0_st <- which(! rowSums(as.matrix(counts_st) == 0) == ncol(counts_st))
counts_st <- counts_st[genes_0_st, ] 
keep_st <- which(Matrix::rowSums(counts_st > 0) >= round(0.05 * ncol(counts_st)))
counts_st <- counts_st[keep_st,] 
counts_st <- as.matrix(counts_st) 
counts_st <- as.data.frame(counts_st)

write.table(st_cnt,file= paste0(output_path, '/merfish_st_gene_count.txt'),
            sep = "\t",quote = FALSE, eol = "\n",row.names = TRUE,col.names = TRUE)


## embryo (sci-Space) data process
load(file = "/home/yll/spatial_decon_methods/RCTD/embryo_singlecell_remove_rare_ct.RData")
load(file="/home/yll/spatial_decon_methods/RCTD/embryo_spatial_remove_rare_ct.RData")

cell_metadata <- filter_cell_metadata
cell_metadata$final_cluster_label <- gsub(" ", "\\.",cell_metadata$final_cluster_label)
sc_lineage <- data.frame(v1 = rownames(cell_metadata),v2 = cell_metadata$final_cluster_label)
sc_lineage <- data.frame(v1 = colnames(filter_count_matrix),v2 = cell_metadata$final_cluster_label)

write.table(sc_lineage,file = paste0(output_path, '/embryo_sc_lineage.txt'),
            sep = "\t",quote = FALSE, eol = "\n",row.names = FALSE,col.names = FALSE)

#embryo single cell count data
filter_count_matrix <- as.matrix(filter_count_matrix)
filter_count_matrix <- as.data.frame(filter_count_matrix)
genes_0_sc <- which(! rowSums(as.matrix(filter_count_matrix) == 0) == ncol(filter_count_matrix))
counts_sc <- filter_count_matrix[genes_0_sc, ] #35148*17309
keep_sc <- which(Matrix::rowSums(counts_sc > 0) >= round(0.05 * ncol(counts_sc)))
counts_sc <- counts_sc[keep_sc,] #9081 genes*17309 cells

write.table(counts_sc,file=paste0(output_path, '/embryo_sc_gene_count.txt'),
            sep = "\t",quote = FALSE, eol = "\n", row.names = TRUE,col.names = TRUE)

#spatial data process
genes_0_st <- which(! rowSums(as.matrix(slide14_spatial_count) == 0) == ncol(slide14_spatial_count))
counts_st <- slide14_spatial_count[genes_0_st, ] #35147*1393
keep_st <- which(Matrix::rowSums(counts_st > 0) >= round(0.05 * ncol(counts_st)))
counts_st <- counts_st[keep_st,] #16774 genes*1393 spots
counts_st <- as.matrix(counts_st)#16774*1393 spots
counts_st <- as.data.frame(counts_st)

write.table(counts_st,file=paste0(output_path, '/embryo_st_gene_count.txt'),
            sep = "\t",quote = FALSE, eol = "\n",row.names = TRUE,col.names = TRUE)

## mouse brain (cellTrek) data process
load("~/scRNA-seq_Mapping_to_STdata/brain_celltrek_data.RData")
load("~/scRNA-seq_Mapping_to_STdata/simulated_brain_spatialdata.RData")

sc_meta_data <- brain_celltrek@meta.data # 34617 genes * 5944 cells

sc_meta_data$subclass <- gsub(" ", "\\.",sc_meta_data$subclass)
sc_meta_data$subclass <- gsub("/", "\\.",sc_meta_data$subclass)
rownames(sc_meta_data) <- gsub("\\.","_", rownames(sc_meta_data))

bra_lineage <- data.frame(v1 = rownames(sc_meta_data),v2 = sc_meta_data$subclass)

write.table(bra_lineage,file=paste0(output_path,'/brain_sc_lineage.txt'),
            sep = "\t",quote = FALSE, eol = "\n",row.names = FALSE,col.names = FALSE)

#brain single cell count data
sc_count <- as.matrix(brain_celltrek@assays$RNA@counts)
filter_count_matrix1 <- as.matrix(sc_count)
filter_count_matrix1 <- as.data.frame(filter_count_matrix1)
colnames(filter_count_matrix1) <- gsub("\\.","_", colnames(filter_count_matrix1))

genes_0_sc <- which(! rowSums(as.matrix(filter_count_matrix1) == 0) == ncol(filter_count_matrix1))
bra_counts_sc <- filter_count_matrix1[genes_0_sc, ] 
keep_sc <- which(Matrix::rowSums(bra_counts_sc > 0) >= round(0.05 * ncol(bra_counts_sc)))
bra_counts_sc <- bra_counts_sc[keep_sc,] 

write.table(counts_sc,file=paste0(output_path,'/brain_sc_gene_count.txt'),
            sep = "\t",quote = FALSE, eol = "\n", row.names = TRUE,col.names = TRUE)

#spatial data process
st_count <- as.matrix(sim_brain_st@assays$spatial@counts)

#data pre-processing
genes_0_st <- which(! rowSums(st_count == 0) == ncol(st_count))
st_count <- st_count[genes_0_st, ] #34617 genes *739 spots
keep_st <- which(Matrix::rowSums(st_count > 0) >= round(0.05 * ncol(st_count)))
st_count <- st_count[keep_st,] #25246 genes*739 spots
st_count <- as.data.frame(st_count)
colnames(st_count) <- paste0("spot_",seq(dim(st_count)[2]))

write.table(st_count,file=paste0(output_path,'/brain_st_gene_count.txt'),
            sep = "\t",quote = FALSE, eol = "\n",row.names = TRUE,col.names = TRUE)
