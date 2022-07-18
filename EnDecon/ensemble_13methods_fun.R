ls()
rm(list = ls())
library(philentropy)
library(readr)

results <- read_csv("~/spatial_decon_methods/ensemble_methods_results/ensemble_results_13methods.csv")
emb_results <- results[which(results$Dataset=="embryo(sci-Space)"), ]
mer_results <- results[which(results$Dataset=="MPOA(MERFISH)"), ]
bra_results <- results[which(results$Dataset=="mouse brain(Mapping)"), ]

#reorder RMSE,Pearson,quant_jsd_50%, the rank of lowest RMSE and lowest JSD is top,
#the rank of highest Pearson is top.
## The rank of different methods in embryo (sci-Space) dataset
emb_results <- emb_results[order(emb_results$RMSE), ]
rank_RMSE <- 14 - rank(emb_results$RMSE)
rank_Pearson <- rank(emb_results$Pearson)
rank_JSD <- 14 - rank(emb_results$`quants_jsd_50%`)

emb_rank <- data.frame(Methods = emb_results$Method,Dataset = emb_results$Dataset,
                       rank_RMSE = rank_RMSE, rank_Pearson = rank_Pearson, rank_JSD = rank_JSD)
emb_results$AS <- rowMeans(emb_rank[3:5])
emb_rank$rank_AS <- 14 - rank(emb_results$AS)#(Ave_score)

## The rank of different methods in MPOA (MERFISH) dataset
mer_results <- mer_results[order(mer_results$RMSE), ]
rank_RMSE <- 14 - rank(mer_results$RMSE)
rank_Pearson <- rank(mer_results$Pearson)
rank_JSD <- 14 - rank(mer_results$`quants_jsd_50%`)

mer_rank <- data.frame(Methods = mer_results$Method,Dataset = mer_results$Dataset,
                       rank_RMSE = rank_RMSE, rank_Pearson = rank_Pearson, rank_JSD = rank_JSD)
mer_results$AS <- rowMeans(mer_rank[3:5])
mer_rank$rank_AS <- 14 - rank(mer_results$AS)

## The rank of different methods in mouse brain (Mapping) dataset
bra_results <- bra_results[order(bra_results$RMSE), ]
rank_RMSE <- 14 - rank(bra_results$RMSE)
rank_Pearson <- rank(bra_results$Pearson)
rank_JSD <- 14 - rank(bra_results$`quants_jsd_50%`)

bra_rank <- data.frame(Methods = bra_results$Method,Dataset = bra_results$Dataset,
                       rank_RMSE = rank_RMSE, rank_Pearson = rank_Pearson, rank_JSD = rank_JSD)
bra_results$AS <- rowMeans(bra_rank[3:5])
bra_rank$rank_AS <- 14 - rank(bra_results$AS)

loc1 <- match(emb_rank$Methods, mer_rank$Methods)
mer_rank <- mer_rank[loc1,]
loc2 <- match(emb_rank$Methods, bra_rank$Methods)
bra_rank <- bra_rank[loc2,]
Ave_rank <- data.frame(Methods = mer_rank$Methods, 
                       bra_rank = bra_rank$rank_AS,
                       mer_rank = mer_rank$rank_AS, 
                      emb_rank = emb_rank$rank_AS)
Ave_rank$Ave_score <- apply(Ave_rank[2:4], 1, mean)

# loc3 <- match(Ave_rank$Methods, rank_results$Methods)
# rank_results <- rank_results[loc3, ]

Ave_rank <- Ave_rank[order(Ave_rank$Ave_score), ]

#methods2: weighted = a/sum(a), a = 1/(1+exp(rank)) 
a <- 1/(1+exp(Ave_rank$Ave_score)) 
Ave_rank$weight <- a/sum(a)

#load 13 different deconvlution results
#read deconvolution results using R: RCTD, spatialDWLS, SPOTlight, Giotto/RAGE/rank/hypergeometric,
#Seurat, spatialDecon, STdeconvolve
decon_result <- read_csv("deconvolution_results/decon_RCTD/decon_result.csv")
inputDir <- paste0(getwd(),'/deconvolution_results/')
files <- list.files(inputDir)
diff_deconlist <- list()
for (f in files) {
  
  decon_mtrx = read_csv(paste0(getwd(),'/deconvolution_results/',f,"/decon_result.csv"))
  diff_deconlist[[f]] = decon_mtrx
  
}

#read deconvolution results using python: cell2location, DestVI, Tangram, STRIDE
#cell2location
decon_cell2location <- read_csv("~/cell2loction/cell2location_means_result.csv")
decon_cell2location <- column_to_rownames(decon_cell2location,var = "...1")
colnames(decon_cell2location) <- gsub("meanscell_abundance_w_sf_","",colnames(decon_cell2location))
Colname <- colnames(decon_cell2location)
decon_cell2location <- t(apply(decon_cell2location,1,as.numeric))
colnames(decon_cell2location) <- Colname 

#Tangram
decon_Tangram <- read.csv2(file = "embryo_Tangram_proportion.csv",sep = ",")
decon_Tangram <- column_to_rownames(decon_Tangram,var="X")

#DestVI
decon_DestVI <- read.csv2(file = "merfish_UV_proportion.csv",sep = ",")
decon_DestVI <- column_to_rownames(decon_DestVI,var="X")
decon_DestVI <- decon_DestVI[,-10]

#STRIDE
decon_STRIDE <- read.table(file="/home/yll/STRIDE/Results/MERFISH_spot_celltype_frac.txt",
                           header = TRUE)

all(rownames(decon_cell2location) == rownames(decon_Hyper))
all(rownames(decon_Hyper) == rownames(decon_PAGE))
all(colnames(decon_PAGE) == colnames(decon_rank))
all(colnames(decon_rank) == colnames(decon_spatialDWLS))
all(colnames(decon_spatialDWLS) == colnames(decon_Seurat))
all(colnames(decon_Seurat) == colnames(decon_spatialDecon))
all(colnames(decon_spatialDecon) == colnames(decon_STRIDE))
all(colnames(decon_STRIDE) == colnames(decon_SPOTlight))
all(colnames(decon_SPOTlight) == colnames(decon_STdeconvolve))

pos <- which(!rownames(decon_cell2location) %in% rownames(decon_RCTD))

decon_cell2location <- decon_cell2location[-pos, ]
decon_cell2location <- decon_cell2location[rownames(expected_result),]
colnames(decon_cell2location) <- gsub("\\.", " ",colnames(decon_cell2location))
decon_cell2location <- decon_cell2location[,colnames(expected_result)]

decon_PAGE <- decon_PAGE[rownames(decon_Hyper), ]
decon_PAGE <- decon_PAGE[-pos, ]
decon_PAGE <- decon_PAGE[rownames(expected_result),]
colnames(decon_PAGE) <- gsub("\\.", " ",colnames(decon_PAGE))
decon_PAGE <- decon_PAGE[,colnames(expected_result)]

decon_Hyper <- decon_Hyper[-pos, ]
decon_Hyper <- decon_Hyper[rownames(expected_result),]
colnames(decon_Hyper) <- gsub("\\.", " ",colnames(decon_Hyper))
decon_Hyper <- decon_Hyper[,colnames(expected_result)]

decon_rank <- decon_rank[-pos, ]
decon_rank <- decon_rank[rownames(expected_result),]
colnames(decon_rank) <- gsub("\\.", " ",colnames(decon_rank))
decon_rank <- decon_rank[,colnames(expected_result)]

decon_Seurat <- decon_Seurat[-pos, ]
decon_Seurat <- decon_Seurat[rownames(expected_result),]
colnames(decon_Seurat) <- gsub("\\.", " ",colnames(decon_Seurat))
decon_Seurat <- decon_Seurat[,colnames(expected_result)]

decon_spatialDecon <- decon_spatialDecon[-pos, ]
decon_spatialDecon <- decon_spatialDecon[rownames(expected_result),]
colnames(decon_spatialDecon) <- gsub("\\.", " ",colnames(decon_spatialDecon))
decon_spatialDecon <- decon_spatialDecon[,colnames(expected_result)]

decon_spatialDWLS <- decon_spatialDWLS[-pos, ]
decon_spatialDWLS <- decon_spatialDWLS[rownames(expected_result),]
colnames(decon_spatialDWLS) <- gsub("\\.", " ",colnames(decon_spatialDWLS))
decon_spatialDWLS <- decon_spatialDWLS[,colnames(expected_result)]

decon_SPOTlight <- decon_SPOTlight[-pos, ]
decon_SPOTlight <- decon_SPOTlight[rownames(expected_result),]
colnames(decon_SPOTlight) <- gsub("\\.", " ",colnames(decon_SPOTlight))
decon_SPOTlight <- decon_SPOTlight[,colnames(expected_result)]

decon_STdeconvolve <- decon_STdeconvolve[-pos, ]
decon_STdeconvolve <- decon_STdeconvolve[rownames(expected_result),]
colnames(decon_STdeconvolve) <- gsub("\\.", " ",colnames(decon_STdeconvolve))
decon_STdeconvolve <- decon_STdeconvolve[,colnames(expected_result)]
decon_STdeconvolve <- decon_STdeconvolve - (1e-12)

decon_STRIDE <- decon_STRIDE[-pos, ]
decon_STRIDE <- decon_STRIDE[rownames(expected_result),]
colnames(decon_STRIDE) <- gsub("\\.", " ",colnames(decon_STRIDE))
decon_STRIDE <- decon_STRIDE[,colnames(expected_result)]

decon_Tangram <- decon_Tangram[-pos, ]
decon_Tangram <- decon_Tangram[rownames(expected_result),]
colnames(decon_Tangram) <- gsub("\\.", " ",colnames(decon_Tangram))
decon_Tangram <- decon_Tangram[,colnames(expected_result)]

decon_DestVI <- decon_DestVI[-pos, ]
decon_DestVI <- decon_DestVI[rownames(expected_result),]
colnames(decon_DestVI) <- gsub("\\.", " ",colnames(decon_DestVI))
decon_DestVI <- decon_DestVI[,colnames(expected_result)]

#Intergation 13 methods by experional weight
w <- Ave_rank$weight
decon_ensemble <- w[1]*decon_cell2location+w[2]*decon_RCTD+w[3]*decon_spatialDWLS+
  w[4]*decon_spatialDecon + w[5]*decon_Seurat +  w[6]*decon_Tangram+
  w[7]*decon_STRIDE + w[8]*decon_Hyper + w[9]*decon_PAGE+
  w[10]*decon_SPOTlight+w[11]*decon_rank+w[12]*decon_STdeconvolve+w[13]*decon_DestVI

result <- benchmark_performance(ground_truth_mtrx = as.matrix(st_seurat@meta.data[7:15]),
                                deconvoluted_mtrx = decon_ensemble)







