#stereoscope data process
library(tibble)
library(slam)

output_path = getwd()
## MPOA (MERFISH) data process
load("/home/yll/spatial_decon_methods/STdeconvolve_ReferenceFree/MERFISH_singlecell_dataset_new.RData")
load("~/spatial_decon_methods/spatialDWLS_dataset-main/codes_yll/MERFISH_spatialspot_dataset.RData")

#spatial spot data
st_count <- as.matrix(simFN7$sim)
st_count <- as.data.frame(st_count)
rownames(st_count) <- gsub("_","x",rownames(st_count))
write.table(st_count,file=paste0(output_path,'/merfish_spatial_cnt.tsv'),
                                 row.names = TRUE,col.names = TRUE,sep="\t")

#select single cell data, the cell number of each cell type =250
annotDf_N7 <-simFN7$annotDf
uni_ct <- unique(annotDf_N7$Cell_class)
new_cnt <- NULL
new_mta <-NULL
for (i in uni_ct){
  mta_ct <- annotDf_N7[annotDf_N7$Cell_class==i,]#choose "celltype i" meta data
  cnt_ct <- counts_[rownames(mta_ct), ]#according to rownames(mta_ct), choose "celltype i" count data
  set.seed(1234)
  zidx <- sample(x = rownames(mta_ct),size=250)#random samlping 250 sample
  pos <- which(rownames(cnt_ct) %in%  zidx)
  new_ct <- cnt_ct[pos,]
  new_meta <- mta_ct[pos, ]
  new_cnt <- rbind(new_cnt,new_ct)
  new_mta <- rbind(new_mta,new_meta)
}

#save selected merfish single cell count data
new_sc_count <- as.data.frame(new_cnt)
new_sc_count <- round(new_sc_count)
new_sc_count <- rownames_to_column(new_sc_count,var = "cell")
write.table(new_sc_count,file=paste0(output_path, '/merfish_selected_singlecell_cnt.tsv'),
            row.names = FALSE,col.names = TRUE,sep="\t")

new_mer_meta <- data.frame(cell = new_sc_count$cell, bio_celltype = new_mta$Cell_class)
new_mer_meta$bio_celltype <- gsub(" ","\\.",new_mer_meta$bio_celltype)
write.table(new_mer_meta,file=paste0(output_path,'/merfish_selected_singlecell_mta.tsv'),
            row.names = FALSE,col.names = TRUE,sep="\t")

#single cell stats data
new_stats <- table(new_mer_meta$bio_celltype)
new_stats <- as.data.frame(new_stats)
new_stats$Var1 <- as.character(new_stats$Var1)
colnames(new_stats) <- c("cell","members")
new_stats$cell <- gsub(" ","\\.",new_stats$cell)
write.table(new_stats,file=paste0(output_path,'/merfish_selected_singlecell_stats.tsv'),
            row.names = FALSE,col.names = TRUE,sep="\t")


## embryo (sci-Space) data process
load("~/spatial_decon_methods/RCTD/embryo_spatial_remove_rare_ct.RData")
load("~/spatial_decon_methods/RCTD/embryo_singlecell_remove_rare_ct.RData")

# data pre-processing
genes_0_st <- which(! rowSums(as.matrix(slide14_spatial_count) == 0) == ncol(slide14_spatial_count))
counts_st <- slide14_spatial_count[genes_0_st, ] 
keep_st <- which(Matrix::rowSums(counts_st > 0) >= round(0.05 * ncol(counts_st)))
counts_st <- counts_st[keep_st,] 
counts_st <- t(counts_st)
dim(counts_st)
counts_st <- as.data.frame(counts_st) 
rownames(counts_st) <- gsub(",","x",rownames(counts_st))
rownames(counts_st) <- gsub("spot_","",rownames(counts_st))

write.table(emb_st_count,file=pastes0(output_path,'/embryo_spatial_cnt.tsv'),
            row.names = TRUE,col.names = TRUE,sep="\t")

#single cell data
filter_count_matrix<-as.matrix(filter_count_matrix)
genes_0_sc <- which(! rowSums(as.matrix(filter_count_matrix) == 0) == ncol(filter_count_matrix))
counts_sc <- filter_count_matrix[genes_0_sc, ] #35147*17301
keep_sc <- which(Matrix::rowSums(counts_sc > 0) >= round(0.05 * ncol(counts_sc)))
counts_sc <- counts_sc[keep_sc,] #9081 genes*17301 cells
counts_sc <- t(counts_sc)
counts_sc <- as.data.frame(counts_sc)

#select single cell data, the cell number of each cell type =250
uni_ct <- unique(filter_cell_metadata$final_cluster_label)
new_cnt <- NULL
new_mta <-NULL
new_sc_data <- lapply(uni_ct,function(n){
  ct_num <- length(which(filter_cell_metadata$final_cluster_label==n))
  if (ct_num >= 250){
    mta_ct <- filter_cell_metadata[filter_cell_metadata$final_cluster_label==n,]#choose "celltype i" meta data
    cnt_ct <- counts_sc[rownames(mta_ct), ]#according to rownames(mta_ct), choose "celltype i" count data
    set.seed(1234)
    zidx <- sample(x = rownames(mta_ct),size=250)#random samlping 250 sample
    pos <- which(rownames(cnt_ct) %in%  zidx)
    new_ct <- cnt_ct[pos,]
    new_meta <- mta_ct[pos, ]
  }else{
    new_meta <- filter_cell_metadata[filter_cell_metadata$final_cluster_label==n,]#choose "celltype i" meta data
    new_ct <- counts_sc[rownames(new_meta), ]#according to rownames(mta_ct), choose "celltype i" count data
  }
  new_cnt <- rbind(new_cnt,new_ct)
  new_mta <- rbind(new_mta,new_meta)
  return(list(new_cnt,new_mta))
})

new_cnt <- NULL
new_mta <- NULL
for (i in seq_len(18)){
  new_cnt <- dplyr::bind_rows(new_cnt,new_sc_data[[i]][[1]])
  new_mta <- dplyr::bind_rows(new_mta,new_sc_data[[i]][[2]])
}

#save selected single cell count data
new_sc_count <- as.data.frame(new_cnt)
new_sc_count <- rownames_to_column(new_sc_count,var = "cell")
write.table(new_sc_count,file=paste0(output_path, '/embryo_selected_singlecell_cnt.tsv'),
            row.names = FALSE,col.names = TRUE,sep="\t")

new_emb_meta <- data.frame(cell = new_sc_count$cell, bio_celltype = new_mta$final_cluster_label)
new_emb_meta$bio_celltype <- gsub(" ","\\.",new_emb_meta$bio_celltype)
write.table(new_emb_meta,file=paste0(output_path, '/embryo_selected_singlecell_mta.tsv'),
            row.names = FALSE,col.names = TRUE,sep="\t")

#single cell stats data
new_stats <- table(new_emb_meta$bio_celltype)
new_stats <- as.data.frame(new_stats)
new_stats$Var1 <- as.character(new_stats$Var1)
colnames(new_stats) <- c("cell","members")
new_stats$cell <- gsub(" ","\\.",new_stats$cell)
write.table(new_stats,file=paste0(output_path, '/embryo_selected_singlecell_stats.tsv'),
            row.names = FALSE,col.names = TRUE,sep="\t")

## brain (cellTrek) data process
load("~/scRNA-seq_Mapping_to_STdata/simulated_brain_spatialdata.RData")#load ST data
load("~/scRNA-seq_Mapping_to_STdata/brain_celltrek_data.RData")#load singcell data
#spatial spot data
st_count <- as.matrix(sim_brain_st@assays$spatial@counts)
genes_0_st <- which(! rowSums(as.matrix(st_count) == 0) == ncol(st_count))
counts_st <- st_count[genes_0_st, ] #35147*1393
keep_st <- which(Matrix::rowSums(counts_st > 0) >= round(0.05 * ncol(counts_st)))
counts_st <- counts_st[keep_st,] #[1] 25267   739
counts_st <- as.data.frame(counts_st)
counts_st <- t(counts_st) #[1]   739 25267,matrix
bra_st_count <- as.data.frame(counts_st) # 739spots * 25267genes, data.frame
rownames(bra_st_count) <- gsub("_","x",rownames(bra_st_count))
write.table(bra_st_count,file=paste0(output_path,'/brain_spatial_cnt.tsv'),
            row.names = TRUE,col.names = TRUE,sep="\t")

#single cell data
brain_celltrek <- sim_brain_celltrek
counts_sc <- as.matrix(brain_celltrek@assays$RNA@counts)
counts_sc <- as.data.frame(counts_sc)
colnames(counts_sc) <- gsub("\\.","_", colnames(counts_sc))
genes_0_sc <- which(! rowSums(as.matrix(counts_sc) == 0) == ncol(counts_sc))
counts_sc <- counts_sc[genes_0_sc, ] #35148*17309
keep_sc <- which(Matrix::rowSums(counts_sc > 0) >= round(0.05 * ncol(counts_sc)))
counts_sc <- counts_sc[keep_sc,] #9081 genes*17309 cells
counts_sc <- t(counts_sc)
counts_sc <- as.data.frame(counts_sc)

#select single cell data, the cell number of each cell type =250
bra_meta <- brain_celltrek@meta.data
rownames(bra_meta) <- rownames(counts_sc)
bra_meta$subclass <- gsub(" ","\\.",bra_meta$subclass)
bra_meta$subclass <- gsub("/","\\.",bra_meta$subclass)
uni_ct <- unique(bra_meta$subclass)
new_cnt <- NULL
new_mta <- NULL
new_sc_data <- lapply(uni_ct,function(n){
  ct_num <- length(which(bra_meta$subclass==n))
  if (ct_num >= 250){
    mta_ct <- bra_meta[bra_meta$subclass==n,]#choose "celltype i" meta data
    cnt_ct <- counts_sc[rownames(mta_ct), ]#according to rownames(mta_ct), choose "celltype i" count data
    set.seed(1234)
    zidx <- sample(x = rownames(mta_ct),size=250)#random samlping 250 sample
    #pos <- which( rownames(cnt_ct)  %in% zidx)
    new_ct <- cnt_ct[zidx,]
    new_meta <- mta_ct[zidx, ]
  }else{
    new_meta <- bra_meta[bra_meta$subclass==n,]#choose "celltype i" meta data
    new_ct <- counts_sc[rownames(new_meta), ]#according to rownames(mta_ct), choose "celltype i" count data
  }
  new_cnt <- rbind(new_cnt,new_ct)
  new_mta <- rbind(new_mta,new_meta)
  return(list(new_cnt,new_mta))
})

new_cnt <- NULL
new_mta <- NULL
for (i in seq_len(15)){
  new_cnt <- rbind(new_cnt,new_sc_data[[i]][[1]])
  new_mta <- rbind(new_mta,new_sc_data[[i]][[2]])
}

new_sc_count <- as.data.frame(new_cnt)
new_sc_count <- rownames_to_column(new_sc_count,var = "cell")
write.table(new_sc_count,file=paste0(output_path,'/brain_selected_singlecell_cnt.tsv'),
            row.names = FALSE,col.names = TRUE,sep="\t")

new_bra_meta <- data.frame(cell = new_sc_count$cell, bio_celltype = new_mta$subclass)
new_bra_meta$bio_celltype <- gsub(" ","\\.",new_bra_meta$bio_celltype)
new_bra_meta$bio_celltype <- gsub("/","\\.",new_bra_meta$bio_celltype)
write.table(new_bra_meta,file=paste0(output_path,'/brain_selected_singlecell_mta.tsv'),
            row.names = FALSE,col.names = TRUE,sep="\t")

#single cell stats data
new_stats <- table(new_bra_meta$bio_celltype)
new_stats <- as.data.frame(new_stats)
new_stats$Var1 <- as.character(new_stats$Var1)
colnames(new_stats) <- c("cell","members")
new_stats$cell <- gsub(" ","\\.",new_stats$cell)
write.table(new_stats,file=paste0(output_path,'/brain_selected_singlecell_stats.tsv'),
            row.names = FALSE,col.names = TRUE,sep="\t")





                       

            
            

