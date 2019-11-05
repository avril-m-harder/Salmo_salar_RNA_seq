## Using output from permutations.R, calculate means and 83% or 95% quantiles for
## slope and intercept values

setwd("/scratch/snyder/h/harder/ASperm_permutations/log_haz_values_incl_y_val_cals/")
data <- read.csv("permutation_output_coeffs_and_pred_ys.csv", sep=",")

colnames(data)[1] <- "gene"

sink("perm_coeffs_and_pred_ys_quantiles_95.csv", type="output")

for(jean in (unique(data$gene))) {
  temp <- data[which(data$gene==jean),]
  cat(paste0(jean,","))
  # cat(apply(temp[2:ncol(temp)], 2, quantile, probs=c(0.085, 0.915)), sep=",") ## for 83% quantiles
  cat(apply(temp[2:ncol(temp)], 2, quantile, probs=c(0.025, 0.975)), sep=",") ## for 95% quantiles
  cat("\n")
}

sink()