#SPOTlight deconvolution

library(Matrix)
library(data.table)
library(Seurat)
library("dplyr",warn.conflicts = FALSE)
library(SPOTlight)
library("igraph",warn.conflicts = FALSE)
library(RColorBrewer)
library(pryr)
library(ggplot2)
library(data.table)
library(tibble)

output_path <- getwd()
#load data
ptm <- proc.time()
load("/home/yll/spatial_decon_methods/STdeconvolve_ReferenceFree/MERFISH_singlecell_dataset_new.RData")
load("~/spatial_decon_methods/spatialDWLS_dataset-main/codes_yll/MERFISH_spatialspot_dataset.RData")

#reference data pre-processing
sc_meta <- annot.table_new
filter_count_matrix <- t(as.matrix(counts_new))

genes_0_sc <- which(! rowSums(as.matrix(filter_count_matrix) == 0) == ncol(filter_count_matrix))
count_sc <- filter_count_matrix[genes_0_sc, ] 
keep_sc <- which(Matrix::rowSums(count_sc > 0) >= round(0.05 * ncol(count_sc)))
count_sc <- count_sc[keep_sc,] 
count_sc <- as.data.frame(count_sc)

se_sc <- CreateSeuratObject(counts = count_sc,
                               meta.data = sc_meta,
                               assay = "RNA")
#spatial data pre-processing
count_st <- t(as.matrix(simFN7$sim))
spatial_loc <- simFN7$cellCounts 
rownames(spatial_loc) <- colnames(count_st)

genes_0_st <- which(! rowSums(as.matrix(count_st) == 0) == ncol(count_st))
count_st <- count_st[genes_0_st, ] 
keep_st <- which(rowSums(count_st > 0) >= round(0.05 * ncol(count_st)))
count_st <- count_st[keep_st,] 
count_st <- as.data.frame(count_st)

st_meta <- mutate(spatial_loc,spot_ID = colnames(count_st))
rownames(st_meta) <- colnames(count_st)
se_st <- CreateSeuratObject(counts = count_st,
                               meta.data = st_meta,
                               assay = "Spatial")

#pre-processing
set.seed(123)
se_sc <- Seurat::SCTransform(se_sc, verbose=FALSE) %>%
  Seurat::RunPCA(., verbose=FALSE)%>%
  Seurat::RunUMAP(., dims=1:30, verbose=FALSE)

#Compute marker gene
Seurat::Idents(object=se_sc) <- se_sc@meta.data$celltype
cluster_markers_all <- Seurat::FindAllMarkers(object = se_sc,
                                            assay = "SCT",
                                            slot = "data",
                                            verbose=TRUE,
                                            only.pos=TRUE)

saveRDS(object = cluster_markers_all,
        file = here::here("inst/markers_embryo_sc.RDS"))

#Spotlight decomposition
set.seed(123)
spotlight_ls <- spotlight_deconvolution(
  se_sc = se_sc,
  counts_spatial = se_st@assays$Spatial@counts,
  clust_vr = "celltype", #variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all,#Dataframe with the marker genes
  cl_n = 100,#number of cells per cell type to use
  hvg = 3000,#number of HVG to use
  ntop = NULL,#how many of the marker genes to use (by default all)
  transf = "uv",#perform unit-variance scaling per cell and spot prior to factorization and NLS
  method = "nsNMF",#factorization method
  min_cont = 0.01 #remove those cells contributing to a spot below a certain threshold
)

decon_mtrx <- spotlight_ls[[2]]
decon_mtrx <- decon_mtrx[,-16]

write.csv(decon_mtrx, 
          paste0(output_path, '/SPOTlight_result.csv'),
          row.names = FALSE, col.names = TRUE, sep=",")

runtime_spotlight <- (proc.time() - ptm)/60
memory_spotlight <- mem_used()

cat(paste0("Running time is:", round(runtime[3],3), "min"),
    paste0("Memory usage is:", round(memory/1024/1024/1024,3), "GB"),
    sep = "\n")
