setwd("/home/yll/spatial_decon_methods/ensemble_methods_results")
ls()
rm(list = ls())

library(readr)
results <- read_csv("~/spatial_decon_methods/ensemble_methods_results/ensemble_results_13methods.csv")
emb_results <- results[which(results$Dataset=="embryo(sci-Space)"), ]
mer_results <- results[which(results$Dataset=="MPOA(MERFISH)"), ]
bra_results <- results[which(results$Dataset=="mouse brain(Mapping)"), ]

#reorder RMSE,Pearson,quant_jsd_50%, the rank of lowest RMSE and lowest JSD is top,
#the rank of highest Pearson is top.
## The rank of different methods in embryo (sci-Space) dataset
emb_results <- emb_results[order(emb_results$RMSE), ]
rank_RMSE <- 14 - rank(emb_results$RMSE)#rank_RMSE值越大表达排名越靠前，性能越好
rank_Pearson <- rank(emb_results$Pearson)#rank_RMSE值越大表达排名越靠前，性能越好
rank_JSD <- 14 - rank(emb_results$`quants_jsd_50%`)#rank_RMSE值越大表达排名越靠前，性能越好

emb_rank <- data.frame(Methods = emb_results$Method,Dataset = emb_results$Dataset,
                       rank_RMSE = rank_RMSE, rank_Pearson = rank_Pearson, rank_JSD = rank_JSD)
emb_results$AS <- rowMeans(emb_rank[3:5])#值越大表明排名越靠前，性能越好
emb_rank$rank_AS <- 14 - rank(emb_results$AS)#值越小表明排名越靠前，性能越好

## The rank of different methods in MPOA (MERFISH) dataset
mer_results <- mer_results[order(mer_results$RMSE), ]
rank_RMSE <- 14 - rank(mer_results$RMSE)#rank_RMSE值越大表达排名越靠前，性能越好
rank_Pearson <- rank(mer_results$Pearson)#rank_RMSE值越大表达排名越靠前，性能越好
rank_JSD <- 14 - rank(mer_results$`quants_jsd_50%`)#rank_RMSE值越大表达排名越靠前，性能越好

mer_rank <- data.frame(Methods = mer_results$Method,Dataset = mer_results$Dataset,
                       rank_RMSE = rank_RMSE, rank_Pearson = rank_Pearson, rank_JSD = rank_JSD)
mer_results$AS <- rowMeans(mer_rank[3:5])#值越大表明排名越靠前，性能越好
mer_rank$rank_AS <- 14 - rank(mer_results$AS)#值越小表明排名越靠前，性能越好

## The rank of different methods in mouse brain (Mapping) dataset
bra_results <- bra_results[order(bra_results$RMSE), ]
rank_RMSE <- 14 - rank(bra_results$RMSE)#rank_RMSE值越大表达排名越靠前，性能越好
rank_Pearson <- rank(bra_results$Pearson)#rank_RMSE值越大表达排名越靠前，性能越好
rank_JSD <- 14 - rank(bra_results$`quants_jsd_50%`)#rank_RMSE值越大表达排名越靠前，性能越好

bra_rank <- data.frame(Methods = bra_results$Method,Dataset = bra_results$Dataset,
                       rank_RMSE = rank_RMSE, rank_Pearson = rank_Pearson, rank_JSD = rank_JSD)
bra_results$AS <- rowMeans(bra_rank[3:5])#值越大表明排名越靠前，性能越好
bra_rank$rank_AS <- 14 - rank(bra_results$AS)#值越小表明排名越靠前，性能越好

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

load("/home/yll/spatial_decon_methods/ensemble_methods_results/brain_11methods_results/brain_11methods_decon_mtrix.RData")
load("~/spatial_decon_methods/ensemble_methods_results/brain_11methods_results/brain_groundtruth.RData")

expected_result <- as.matrix(slide14_ground_truth)
colnames(expected_result) <- gsub("\\.", " ",colnames(expected_result))
pos <- which(!rownames(expected_result) %in% rownames(decon_STdeconvolve))
expected_result <- expected_result[-pos, ]

w <- Ave_rank$weight
decon_ensemble <- w[1]*decon_cell2location+w[2]*decon_RCTD+w[3]*decon_spatialDWLS+
  w[4]*decon_spatialDecon + w[5]*decon_Seurat +  w[6]*decon_Tangram+
  w[7]*decon_STRIDE + w[8]*decon_Hyper + w[9]*decon_PAGE+
  w[10]*decon_SPOTlight+w[11]*decon_rank+w[12]*decon_STdeconvolve+w[13]*decon_DestVI

observed_values <- unlist(lapply(seq_len(ncol(decon_ensemble)), function(i) decon_ensemble[,i]))
expected_values <- unlist(lapply(seq_len(ncol(expected_result)), function(i) expected_result[,i]))

#total RMSE and Pearson
RMSE <- sqrt(mean((observed_values - expected_values)^2)) %>% round(.,4)
Pearson <- cor(observed_values,expected_values) %>% round(.,4)

x <- rbind(observed_values,expected_values)
JSD <- JSD(x = x, unit = "log2",est.prob = "empirical") %>% round(.,4)

#RMSE each type
integrated_RMSE_ct <- NULL
for (i in seq_len(dim(expected_result)[2])){
  rmse<- sqrt(mean((decon_ensemble[,i] - expected_result[,i])^2))
  integrated_RMSE_ct <- c(integrated_RMSE_ct,rmse)
}

##### Get TRUE JSD between real-predicted proportions #####
true_jsd_mtrx <- matrix(nrow = nrow(decon_ensemble), ncol = 1)
for (i in seq_len(nrow(decon_ensemble))) {
  
  # Create matrix to feed to JSD
  x <- rbind(decon_ensemble[i, ],
             expected_result[i, ])
  
  # Calculate JSD and save it in true_JSD_mtrx
  if(sum(expected_result[i, ]) > 0) {
    true_jsd_mtrx[i, 1] <- suppressMessages(JSD(x = x, unit = "log2",
                                                est.prob = "empirical"))
  } else {
    true_jsd_mtrx[i, 1] <- 1
  }
}; rm(i)

quants_jsd <- round(quantile(matrixStats::rowMins(true_jsd_mtrx,
                                                  na.rm = TRUE),
                             c(0.25, 0.5, 0.75)), 5)

cat(sprintf("The following summary statistics are obtained:
              RMSE: %s,
              Pearson: %s,
              JSD: %s,
              JSD quantiles: %s[%s-%s]",
            RMSE, Pearson, JSD,quants_jsd[[2]], quants_jsd[[1]], quants_jsd[[3]]), sep = "\n")


Ave_rank <- Ave_rank[order(Ave_rank$Ave_score), ]
Ave_rank$Methods <- factor(Ave_rank$Methods,
                               levels = c("cell2location","RCTD","spatialDWLS","spatialDecon","Seurat 3.0","Tangram",
                                          "STRIDE","Giotto/Hypergeometric","Giotto/PAGE","SPOTlight","Giotto/rank",
                                          "STdeconvolve","DestVI"))
Ave_rank$Methods <- factor(Ave_rank$Methods,
                          levels = c("DestVI","STdeconvolve","Giotto/rank","SPOTlight","Giotto/PAGE","Giotto/Hypergeometric",
                                     "STRIDE","Tangram","Seurat 3.0","spatialDecon","spatialDWLS","RCTD",
                                      "cell2location"))
mycolor <- c("DestVI" = "#E6550DFF",
             "Giotto/rank" = "#636363FF" ,
             "STdeconvolve"="#9ECAE1FF",
             "SPOTlight"="#969696FF",
             "Giotto/PAGE" = "#756BB1FF",
             "spatialDecon"="#74C476FF",
             "Giotto/Hypergeometric"="#31A354FF",
             "Seurat 3.0"="#FD8D3CFF",
             "STRIDE" = "#A1D99BFF",
             "Tangram" = "#BCBDDCFF",
             "RCTD" = "#6BAED6FF",
             "spatialDWLS"="#9E9AC8FF",
             "cell2location"="#3182BDFF")

ggplot(data=Ave_rank, aes(x=Methods,y=Ave_score,fill=Methods)) +
  geom_bar(stat = "identity",position=position_dodge(0.75),width = 0.6) +
  scale_fill_manual(values =  mycolor) + 
  #facet_wrap(~Dataset) + 
  coord_flip() + 
  #scale_y_reverse() +
  theme_bw() + #去除背景色
  theme(panel.grid = element_blank()) +
  theme(axis.text.x=element_text(size=8,face="bold",vjust=1, hjust = 1,angle=0)) +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size=8)) +
  labs(x="Methods", y="Rank",fill="Methods") +
  guides(fill="none")

ggsave(file="/home/yll/spatial_decon_methods/evaluation_metric/figure/ensembel_13methods_rank.pdf",
       width=8, height=6)



