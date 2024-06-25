setwd("/home/yll/benchmark_github/Dataset/embryo_sci-Space_raw")
library(readr)
library(Matrix)
library(dplyr)
library(vec2dtransf)
library(tidyr)
library(tibble)

#load raw data
cell_meta <- read_tsv('GSE166692_sciSpace_cell_metadata_new.tsv')
gene_meta <- read_tsv('GSE166692_sciSpace_gene_metadata.tsv')
count_mat <- Matrix::readMM('GSE166692_sciSpace_count_matrix.mtx')
table(cell_meta$slide_id)

#remove cell type with cell number < 25 
pos_OPCs <- which(cell_meta$final_cluster_label == "OPCs")
pos_Testis <- which(cell_meta$final_cluster_label == "Testis Cells")
pos <- c(pos_OPCs,pos_Testis)
cell_meta <- cell_meta[-pos,]
count_mat <- count_mat[,-pos]
colnames(count_mat) <- cell_meta$Cell

#selected silde14 data
slide14_cell_meta <- subset.data.frame(cell_meta,
                                       cell_meta$slide_id=="Slide 14")
rownames(slide14_cell_meta) <- NULL
pos_slide14 <- which(cell_meta$slide_id=="Slide 14")
slide14_count_mat <- count_mat[,pos_slide14]

slide14_cell_meta <- subset.data.frame(cell_meta,cell_meta$slide_id=="Slide 14")
rownames(slide14_cell_meta) <- NULL
pos_slide14 <- which(cell_meta$slide_id=="Slide 14")
slide14_count_mat <- count_mat[,pos_slide14]

#obtain gene name
n <- length(gene_meta$gene_short_name)
gene <- lapply(1:n,function(x){
  strsplit(gene_meta$gene_short_name[x],"\\t")[[1]][2]
})
gene <- unlist(gene) 
rownames(slide14_count_mat) <- gene
colnames(slide14_count_mat) <- slide14_cell_meta$Cell

#obtain spatial coordinate
x_slide14 <- slide14_cell_meta$Row 
y_slide14 <- slide14_cell_meta$Col 
x_slide14 <- as.data.frame(x_slide14)
y_slide14 <- as.data.frame(y_slide14)
cell_coords <- mutate(x_slide14,y_slide14) 
colnames(cell_coords) <- c("Row","Col")

unique_coords <- unique(paste(cell_coords$Row,cell_coords$Col,sep=","))
celltype_list <- list()  
local_info_list <- list() 
for (k in unique_coords){
  
  pos_cell = which(cell_coords$Row==strsplit(k,split=",")[[1]][1] & cell_coords$Col==strsplit(k,split=",")[[1]][2])
  local_celltype = cell_meta$final_cluster_label[pos_cell]
  local_info = as.data.frame(cbind(pos_cell,local_celltype))
  local_info_list[[k]] = local_info  
  
  celltype_number = as.data.frame(table(cell_meta$final_cluster_label[pos_cell]))
  cell_prop = celltype_number$Freq/sum(celltype_number$Freq)
  celltype_info = as.data.frame(cbind(celltype_number,cell_prop))
  celltype_list[[k]] = celltype_info
}

#generate simulated spot by synthesising cells with the same coordinates
ds_spot<-lapply(unique_coords,function(n){
  
  
  if (dim(celltype_list[[n]])[1]==1){
    
    pos_spot<-as.numeric(unlist(local_info_list[[n]][1]))
    prop<-celltype_list[[n]][3]
    syn_spot_count<-slide14_count_mat[,pos_spot]
    name_simp <- paste("spot_", n, sep = "")
    
    tmp_ds<-celltype_list[[n]] 
    colnames(tmp_ds)<-c("cluster","Freq","cell_prop")
    spot_ds <- tmp_ds %>%
      dplyr::select(cluster,cell_prop) %>%
      tidyr::pivot_wider(names_from = cluster,
                         values_from = cell_prop) %>%
      dplyr::mutate(name = name_simp)
    
    syn_spot <- rowSums(as.matrix(syn_spot_count)); sum(syn_spot)
    
  }else{
    
    #celltype_select<-celltype_list[[n]][which(celltype_list[[n]]$Freq>=2),]
    #celltype_name_spot<-as.character(celltype_select$Var1) 
    #pos_select<-local_info_list[[n]][which(local_info_list[[n]]$local_celltype %in% celltype_name_spot),]
    
    pos_spot<-as.numeric(unlist(local_info_list[[n]][1]))
    prop<-celltype_list[[n]][3]
    #tmp_ds<-cell_metadata_new[pos_spot,] %>% mutate(prop)
    syn_spot_count<-slide14_count_mat[,pos_spot]
    name_simp <- paste("spot_", n, sep = "")
    
    tmp_ds<-celltype_list[[n]] 
    colnames(tmp_ds)<-c("cluster","Freq","cell_prop")
    spot_ds <- tmp_ds %>%
      dplyr::select(cluster,cell_prop) %>%
      tidyr::pivot_wider(names_from = cluster,
                         values_from = cell_prop) %>%
      dplyr::mutate(name = name_simp)
    
    syn_spot <- rowSums(as.matrix(syn_spot_count)); sum(syn_spot)
    
  }
  
  return(list(syn_spot, spot_ds))
})

#generate simulated ST metadata
ds_spot_metadata <- purrr::map(ds_spot, 2) %>% 
  dplyr::bind_rows() %>%
  data.frame()
ds_spot_metadata[is.na(ds_spot_metadata)] <- 0

ds_spot_metadata <- select(ds_spot_metadata,name, Connective.Tissue.Progenitors, everything())

ds_spot_count <- purrr::map(ds_spot,1) %>%
  dplyr::bind_cols() %>%
  data.frame()

#generate spot name
ds_spot_name<-list()
for (k in unique_coords){
  name = paste("spot_", k, sep = "")
  ds_spot_name[[k]] = name
}
ds_spot_name <- unlist(ds_spot_name)
colnames(ds_spot_count) <- ds_spot_name

#remove repeated gene
ds_count <- mutate(ds_spot_count,GENE=gene)
ds_count_unique <- aggregate(ds_count,by=list(ds_count$GENE),FUN=mean,na.rm=TRUE)
ds_spot_unique <- ds_count_unique[,-1395]
ds_spot_unique <- column_to_rownames(ds_count_unique,var="Group.1")
ds_spot_unique <- ds_spot_unique[,-1394]

#spatial count data and meta data
slide14_spatial_count <- ds_spot_unique

#Genearte scRNA-seq Seurat object
sc_data <- CreateSeuratObject(counts = slide14_count_mat,
                              meta.data = slide14_cell_meta,
                              assay = "RNA")

#Generate ST Seurat object
rownames(ds_spot_metadata) <- colnames(ds_spot_unique)
st_data <- CreateSeuratObject(counts = ds_spot_unique,
                              meta.data = ds_spot_metadata,
                              assay = "Spatial")

spatial_coords <- read_csv(file=paste0(getwd(),'/Dataset/embryo_sci-Space_raw/spatial_coordinate.csv'))
st_data@images$coordinate <- spatial_coords

saveRDS(sc_data,file="/home/yll/benchmark_github/synthetic_st_dataset/embryo_singlecell_dataset.rds")
saveRDS(st_data,file="/home/yll/benchmark_github/synthetic_st_dataset/embryo_spatialspot_dataset.rds")








