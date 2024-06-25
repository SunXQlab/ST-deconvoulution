setwd("/home/yll/scRNA-seq_Mapping_to_STdata")
# library(devtools)
# install_github("navinlabcode/CellTrek")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ConsensusClusterPlus")

library(CellTrek)
library(dplyr)
library(Seurat)
library(viridis)
library(ConsensusClusterPlus)
library(tidyr)

output_path = getwd()
## using cellTek mapping single cells to spatial location
#load raw data
brain_st_cortex <- readRDS("brain_st_cortex.rds")#31053gene, 1075 spots
brain_sc <- readRDS("brain_sc.rds")#34617genes, 4785 cells

## Rename the cells/spots with syntactically valid names
brain_st_cortex <- RenameCells(brain_st_cortex, new.names=make.names(Cells(brain_st_cortex)))
brain_sc <- RenameCells(brain_sc, new.names=make.names(Cells(brain_sc)))

## co-embed ST and scRNA-seq datasets using traint
brain_traint <- CellTrek::traint(st_data=brain_st_cortex, 
                                 sc_data=brain_sc, 
                                 sc_assay='RNA', 
                                 cell_names='cell_type')

brain_celltrek <- CellTrek::celltrek(st_sc_int=brain_traint, int_assay='traint', 
                                     sc_data=brain_sc, sc_assay = 'RNA', 
                                     reduction='pca', intp=T, intp_pnt=5000, 
                                     intp_lin=F, nPCs=30, ntree=1000, 
                                     dist_thresh=0.55, top_spot=5, 
                                     spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$celltrek

brain_celltrek$cell_type <- factor(brain_celltrek$cell_type, levels=sort(unique(brain_celltrek$cell_type)))
id_raw <- brain_celltrek_test@meta.data$id_raw
counts <- brain_sc@assays$RNA@counts[,id_raw]
brain_celltrek@assays$RNA@counts <- counts

## generate synthetic ST data 
#set spot size
patch_size <- 105
sc_coord <- brain_celltrek@images$anterior1@coordinates

x_edges <- seq(min(sc_coord$imagerow), max(sc_coord$imagerow), patch_size)
inner_x_edges <- x_edges[2:(length(x_edges) - 1)]
y_edges <- seq(min(sc_coord$imagecol), max(sc_coord$imagecol), patch_size)
inner_y_edges <- y_edges[2:(length(y_edges) - 1)]

sc_coord$patch_id <- character(length(rownames(sc_coord)))
for (x in inner_x_edges) {
  for (y in inner_y_edges) {
    patch_id <- paste0(as.character(x), "_", as.character(y))
    patch <- sc_coord[which((sc_coord$imagerow > x) & (sc_coord$imagerow < x + patch_size) & 
                              (sc_coord$imagecol > y) & (sc_coord$imagecol < y + patch_size)), ]
    if (length(rownames(patch)) > 0) {
      sc_coord[rownames(patch), ]$patch_id <- patch_id
    }
  }
}
sc_coord$cell_class <- as.vector(brain_celltrek$cell_type)
patch_cell_types <- sc_coord[which(sc_coord$patch_id != ""), c("patch_id", "cell_class")]
cellTypeTable <- table(patch_cell_types[])#patch_id*cell_class,每个patch中对应的cell_class数目
cellType_stat <- as.data.frame(cellTypeTable)
cellType_stat$cell_class <- factor(cellType_stat$cell_class)
cellType_stat <- spread(data = cellType_stat,
                        key = cell_class,
                        value = Freq)
A <- column_to_rownames(cellType_stat,var="patch_id")
cellType_stat$cellnumber <- apply(A,1,sum)

#write.table(cellTypeTable,file="/home/yll/scRNA-seq_Mapping_to_STdata/patch_size150_celltype.txt")
counts <- t(as.matrix(brain_celltrek@assays$RNA@counts))
patchTotalCells <- rowSums(cellTypeTable)#每个patch种总共的cell number
cellTypeCount <- c()
for (i in seq_len(length(rownames(cellTypeTable)))) {
  patch_num_cell_types <- length(which(cellTypeTable[i, ] != 0))#每个patch中celltype的number
  cellTypeCount <- append(cellTypeCount, patch_num_cell_types)
}
patches <- unique(sc_coord$patch_id[which(!sc_coord$patch_id == "")])
patchGexp <- do.call(rbind, lapply(patches, function(patch) {
  cells <- rownames(sc_coord[which(sc_coord$patch_id == patch), ])
  mat <- as.matrix(counts[cells, ])
  if (length(cells) == 1) {
    patch_counts <- as.vector(mat)
  }
  else if (length(cells) > 1) {
    patch_counts <- colSums(mat)
  }
  else if (length(cells) == 0) {
    cat("WARNING:", bregma, "patch", patch, "had no cells in `counts` to use for simulated gene expression", 
        "\n")
  }
  patch_counts
}))
rownames(patchGexp) <- patches
#cell type proportion
cellType_stat <- column_to_rownames(cellType_stat,var="patch_id")
cellTypeProp <- as.matrix(cellType_stat[,-16])
cellTypeProp <- t(apply(cellTypeProp,1,function(x){x/sum(x)}))
#each spot coordinate
patches <- unique(sc_coord$patch_id[which(!sc_coord$patch_id == "")])
spot_coord.x <- NULL  
spot_coord.y <- NULL 

for (k in patches){
  
  x <- strsplit(k,split="_")[[1]][1] %>% as.numeric() %>% round(3)
  y <- strsplit(k,split="_")[[1]][2] %>% as.numeric() %>% round(3)
  spot_coord.x <- rbind(spot_coord.x,x)
  spot_coord.y <- rbind(spot_coord.y,y)
  
}

new_coord <- data.frame(imagerow=spot_coord.x,imagecol=spot_coord.y)
rownames(new_coord) <- patches
st_counts <- patchGexp

# #downsamping SequencingDepth 
# set.seed(100)
# #st_counts <- sim_brain_st@assays$spatial@counts
# new.st_cnt <- downsampleMatrix(t(st_counts), prop=0.01)
# nUMI_stnew <- colSums(new.st_cnt)
# max(nUMI_stnew)
# st_count_new <- as.matrix(new.st_cnt)

#creat spatial seurat object
cellTypeProp <- as.data.frame(cellTypeProp)
sim_brain_st <- CreateSeuratObject(counts = st_count_new,
                                   meta.data = cellTypeProp,
                                   assay="spatial")
sim_brain_st@images$coordinate <- new_coord

#creat single cell seurat object
sc_count_use <- counts[rownames(patch_cell_types), ]
sc_meta_data_use <- brain_celltrek@meta.data[rownames(patch_cell_types), ]
sc_coord_use <- sc_coord[rownames(patch_cell_types), ]
sim_brain_celltrek <- CreateSeuratObject(counts = t(sc_count_use),
                                         meta.data = sc_meta_data_use,
                                         assay="RNA")
sim_brain_celltrek@images$coordinate <- sc_coord_use

saveRDS(brain_celltrek,file=paste0(output_path,'/brain_celltrek_data.Rds'))
save(sim_brain_st,file=paste0(output_path,'/brain_spatialspot_data.Rds'))


