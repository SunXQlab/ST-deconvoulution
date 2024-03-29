setwd("/home/yll/benchmark_github/Dataset/MPOA_MERFISH_raw")
#synthetic MPOA (MERFISH) ST data
library(readr)
library(Matrix)
library(dplyr)
library(vec2dtransf)
library(tidyr)
library(tibble)

## Moffit et al. 2018 raw data downloaded from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248/
moffit <- read.csv2(file = "Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv",
                    sep = ",")

## split into metadata
annot.table_ <- moffit[,c(1:9)]
## and gene counts
counts_ <- moffit[,-c(1:9)]

## convert gene counts to numerics
counts_ <- as.matrix(data.frame(apply(counts_, 2, function(x) as.numeric(as.character(x)))))

rownames(annot.table_) <- annot.table_[,"Cell_ID"]
rownames(counts_) <- annot.table_[,"Cell_ID"]

## Moffit et al. 2018 Supplemental Table S6 contains the genes profiled
merfish_genes <- readxl::read_excel(path = "aau5324_moffitt_table-s6.xlsx",
                                    range = "A2:D163")
## find the smFISH genes
smFISH_genes <- dplyr::filter(merfish_genes, Barcode == "Sequential stain")

## filter for MERFISH profiled genes
counts_ <- counts_[, which(!colnames(counts_) %in% smFISH_genes$`Gene name`)]
## remove Blank measurements
counts_ <- counts_[, !grepl("Blank", colnames(counts_))] #135genes

## collect the pertinent metadata for the cells in the tissue sections
annot.table_ <- annot.table_[annot.table_$Animal_ID == 2,
                             c('Centroid_X', 'Centroid_Y', 'Bregma', "Cell_class", "Neuron_cluster_ID")]

## collapse the Oligodendrocytes and Endothelial cell-types
annot.table_[grep(pattern = "OD Mature",
                  x = annot.table_$Cell_class),]$Cell_class <- "OD Mature"
annot.table_[grep(pattern = "OD Immature",
                  x = annot.table_$Cell_class),]$Cell_class <- "OD Immature"
annot.table_[grep(pattern = "Endothelial",
                  x = annot.table_$Cell_class),]$Cell_class <- "Endothelial"

## remove "Ambiguous" cell-types
annot.table_ <- annot.table_[annot.table_$Cell_class != "Ambiguous",]
dim(annot.table_)
## [1] 59651     5

annot.table_ <- na.omit(annot.table_) 
## change the type of some of the metadata columns to numerics
annot.table_$Centroid_X <- as.numeric(as.character(annot.table_$Centroid_X))
annot.table_$Centroid_Y <- as.numeric(as.character(annot.table_$Centroid_Y))

## filter the gene counts matrix to only include the filtered cells
counts_ <- counts_[rownames(annot.table_),]
dim(counts_)
## [1] 59651   135

#####################
# generate row data #
#####################
library(STdeconvolve)

## get list of major cell-types of each cell
majorCellTypes <- annot.table_$Cell_class
length(majorCellTypes) #59651 singlecell

## get a similar list but with neuronal cell-types expanded to subtypes
neuronalCellsubtypes <- unlist(lapply(rownames(annot.table_), function(cell){
  class <- annot.table_[cell,]$Cell_class
  neuron <- annot.table_[cell,]$Neuron_cluster_ID
  if (neuron != ""){
    i <- neuron
  } else {
    i <- class
  }
  i
}))
length(neuronalCellsubtypes)

## Using the gene counts of the individual cells, generate a ground truth gene expression profile
## for each of the cell-types

## major cell-types:
cellTypes <- majorCellTypes
cells <- rownames(annot.table_)
mat <- counts_[cells,] #single cell count data
mm <- model.matrix(~ 0 + factor(cellTypes)) #mm=zeros[row=59651,col=9],if celltype(cell1)=A,mm[1,1]=1;else mm[1,1]=0
colnames(mm) <- levels(factor(cellTypes))
gtGexpCellTypes <- t(t(as.matrix(mat)) %*% mm)
gtGexpCellTypes <- gtGexpCellTypes/rowSums(gtGexpCellTypes)
dim(gtGexpCellTypes)

## expanded neuronal cell-types:
cellTypes <- neuronalCellsubtypes
cells <- rownames(annot.table_)
mat <- counts_[cells,]
mm <- model.matrix(~ 0 + factor(cellTypes))
colnames(mm) <- levels(factor(cellTypes))
gtGexpNeuronalCellTypes <- t(t(as.matrix(mat)) %*% mm)
gtGexpNeuronalCellTypes <- gtGexpNeuronalCellTypes/rowSums(gtGexpNeuronalCellTypes)
dim(gtGexpNeuronalCellTypes)

#----------100um2 simulated pixels for 9 major cell-types--------------
## The function takes into account the cell-types in the "Cell_class" column of the metadata
## so make sure this is set to the cell-types of interest
annot.table_$Cell_class <- majorCellTypes

FN7_hash <- simulateBregmaSpots(annot.table_,
                                counts = counts_, # gene counts matrix
                                patch_size = 100) # size of the simulated pixels in um2 (units of the Centroids)

# table of number of cell types in each spot for entire FN7
# spot IDs for all spots in FN7
FN7_spotIDs <- unlist(lapply(hash::keys(FN7_hash), function(ix){
  # table to df
  rownames(FN7_hash[[ix]]$cellTypeTable)
}))
# combine all the cellTypeTables of each bregma in FN7 (counts of each cell type in each spot)
FN7_cellTypeTable <- lapply(hash::keys(FN7_hash), function(ix){
  # table to df
  as.data.frame.matrix(FN7_hash[[ix]]$cellTypeTable)
})
# combine into single df, and because some bregmas may be missing cell types,
# use rbindlist to keep all columns and add NAs to spots for cell types
# they are missing
FN7_cellTypeTable <- data.table::rbindlist(FN7_cellTypeTable, fill = TRUE)
# replace NAs with 0s
FN7_cellTypeTable[is.na(FN7_cellTypeTable)] <- 0
# spot IDs as row names
FN7_cellTypeTable <- as.matrix(FN7_cellTypeTable)

rownames(FN7_cellTypeTable) <- FN7_spotIDs
simBregmasFN7 <- lapply(hash::keys(FN7_hash), function(ix){
  bregma <- buildBregmaCorpus(hashTable = FN7_hash, 
                              bregmaID = ix)
  print(bregma$sim)
  print(dim(bregma$gtSpotTopics))
  print(dim(bregma$gtCtGenes))
  bregma
})

names(simBregmasFN7) <- hash::keys(FN7_hash)

# 1. sim
# combine the sim slam matrices to make the corpus for all spots across all bregmas
sim_N7 <- slam::as.simple_triplet_matrix(do.call(rbind, lapply(names(simBregmasFN7), function(ix){
  m <- as.matrix(simBregmasFN7[[ix]]$sim)
  m
})))
sim_N7 #3072 spots *135 genes

# 2. gtSpotTopics
# each gtSpotTopics ref can have different numbers of cell types that are present in each bregma,
# at least for the neuro. So need a way to combine and have all 75 columns,
# and set 0 for patches that do not have any of one ct
gtSpotTopics_N7 <- lapply(names(simBregmasFN7), function(ix){
  simBregmasFN7[[ix]]$gtSpotTopics
})
names(gtSpotTopics_N7) <- names(simBregmasFN7)
gtSpotTopics_N7 <- data.table::rbindlist(gtSpotTopics_N7, fill = TRUE)
dim(gtSpotTopics_N7)  #ground truth

gtSpotTopics_N7 <- as.matrix(gtSpotTopics_N7)
rownames(gtSpotTopics_N7) <- rownames(sim_N7)
# Cts not present in a bregma but columns added here have NA values
# replace NAs with 0's
gtSpotTopics_N7[is.na(gtSpotTopics_N7)] <- 0
gtSpotTopics_N7[1:3,]

# 3. gtCtGenes
# the `gtCtGenes` is the beta of the average gene expression for each cell cluster,
# in this case using all of cells across the bregma in the animal
dim(gtGexpCellTypes)

# 4. cellCounts
# counts of cells in each simulated spot, but also has spot coordinates for easy plotting
cellCounts_N7 <- do.call(rbind, lapply(names(simBregmasFN7), function(ix){
  m <- simBregmasFN7[[ix]]$cellCounts
  m
}))
dim(cellCounts_N7) #spot coordinates

cellCounts_N7[1:10,] 

# 5. annotDf
# recall this is meta data data frame includes information for only the cells that were kept in simulated spots
annotDf_N7 <- do.call(rbind, lapply(names(simBregmasFN7), function(ix){
  m <- simBregmasFN7[[ix]]$annotDf
  m
}))
dim(annotDf_N7)

annotDf_N7[1:10,]

# construct the list similar to the output of `buildBregmaCorpus()`
simFN7 <- list(sim = sim_N7,
               gtSpotTopics = gtSpotTopics_N7,
               gtCtGenes = gtGexpCellTypes,
               cellCounts = cellCounts_N7,
               # classColors = classColors,
               annotDf = annotDf_N7)

#Genearte scRNA-seq Seurat object
annotDf_N7 <- simFN7$annotDf
sc_data <- CreateSeuratObject(counts = t(counts_[rownames(annotDf_N7), ]),
                              meta.data = annotDf_N7,
                              assay = "RNA")

#Generate ST Seurat object
counts_st <- t(as.matrix(simFN7$sim))
spatial_loc <- simFN7$cellCounts #spatial locations
groundtruth <- as.data.frame(simFN7$gtSpotTopics)
sp_meta <- cbind(spatial_loc, groundtruth)
st_data <- CreateSeuratObject(counts = t(as.matrix(simFN7$sim)),
                                 meta.data = sp_meta,
                                 assay = "Spatial")


saveRDS(sc_data,file="/home/yll/benchmark_github/synthetic_st_dataset/MERFISH_singlecell_dataset.rds")
saveRDS(st_data,file="/home/yll/benchmark_github/synthetic_st_dataset/MERFISH_spatialspot_dataset.rds")












