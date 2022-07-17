#' If you wanna test the performance of your model on synthetic generated test spots you can use this function to benchmark and get a sense of the model's performance.
#'
#' @param test_spots_metadata Object of class matrix containing the ground truth composition of each spot as obtained from the function syn_spot_comb_topic_fun.R, 2nd element.
#' @param spot_composition_mtrx Object of class matrix with the predicted topic probability distributions for each spot.
#' @return This function returns a list with TP, TN, FP, FN and the Jensen-Shannon Divergence index.
#' If you wanna test the performance of your model on synthetic generated test spots you can use this function to benchmark and get a sense of the model's performance.
#'
#' @param ground_truth_mtrx Object of class matrix containing the ground truth composition of each spot as obtained from the function syn_spot_comb_topic_fun.R, 2nd element.
#' @param deconvoluted_mtrx Object of class matrix with the predicted topic probability distributions for each spot.
#' @return This function returns a list with TP, TN, FP, FN and the Jensen-Shannon Divergence index.
#' @export
#' @examples
#'

benchmark_performance <- function(ground_truth_mtrx,
                                  deconvoluted_mtrx) {
  # Check variables
  if (!is.matrix(ground_truth_mtrx)) stop("ERROR: ground_truth_mtrx must be a matrix object!")
  if (!is.matrix(deconvoluted_mtrx)) stop("ERROR: deconvoluted_mtrx must be a matrix object!")
  
  if (!dim(ground_truth_mtrx)[1] == dim(deconvoluted_mtrx)[1]){
    pos <- which(!rownames(ground_truth_mtrx) %in% rownames(deconvoluted_mtrx))
    ground_truth_mtrx <- ground_truth_mtrx[-pos, ]
  }
  
  colnames(deconvoluted_mtrx) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                x = colnames(deconvoluted_mtrx),
                                perl = TRUE)
  colnames(ground_truth_mtrx) <- gsub(pattern = "[[:punct:]]|[[:blank:]]", ".",
                                   x = colnames(ground_truth_mtrx),
                                   perl = TRUE)
  #load required packages
  suppressMessages(require(philentropy))
  suppressMessages(require(tibble))
  suppressMessages(require(dplyr))
  
  deconvoluted_mtrx <- deconvoluted_mtrx[ ,colnames(ground_truth_mtrx)]
  deconvoluted_mtrx <- deconvoluted_mtrx[rownames(ground_truth_mtrx), ]
  
  observed_values <- unlist(lapply(seq_len(ncol(deconvoluted_mtrx)), function(i) deconvoluted_mtrx[,i]))
  expected_values <- unlist(lapply(seq_len(ncol(ground_truth_mtrx)), function(i) ground_truth_mtrx[,i]))
  
  #total RMSE and Pearson
  RMSE <- sqrt(mean((observed_values - expected_values)^2)) %>% round(.,4)
  Pearson <- cor(observed_values,expected_values) %>% round(.,4)
  
  #RMSE for each type
  RMSE_ct <- NULL
  for (i in seq_len(dim(ground_truth_mtrx)[2])){
    rmse<- sqrt(mean((deconvoluted_mtrx[,i] - ground_truth_mtrx[,i])^2))
    RMSE_ct <- c(RMSE_ct,rmse)
  }
  
  #RMSE for each type
   PCC_ct <- NULL
   for (i in seq_len(dim(ground_truth_mtrx)[2])){
     a = deconvoluted_mtrx[,i]
     b = ground_truth_mtrx[,i]
     pcc <- cor(a,b)%>% round(.,4)
     PCC_ct <- c(PCC_ct,pcc)
   }
  ##### JSD between real-predicted proportions #####
  true_jsd_mtrx <- matrix(nrow = nrow(ground_truth_mtrx), ncol = 1)
  for (i in seq_len(nrow(ground_truth_mtrx))) {

    # Create matrix to feed to JSD
    x <- rbind(ground_truth_mtrx[i, ],
               deconvoluted_mtrx[i, ])

    # Calculate JSD and save it in true_JSD_mtrx
    if(sum(deconvoluted_mtrx[i, ]) > 0) {
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
              JSD quantiles: %s[%s-%s]",
              RMSE, Pearson,quants_jsd[[2]], quants_jsd[[1]], quants_jsd[[3]]), sep = "\n")

  return(list(RMSE = RMSE, RMSE_perct = RMSE_ct, 
              PCC = Pearson, PCC_perct = PCC_ct, JSD = quants_jsd))
}
