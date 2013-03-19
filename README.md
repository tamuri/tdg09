



Loading the output in R:

library(yaml)
tdg_output <- yaml.load_file('tdg_output.txt')

summary(tdg_output)

lrt_results_df <- as.data.frame(matrix(unlist(tdg_output$LrtResults), ncol=5, byrow=T))
names(lrt_results_df) <- c("site", "deltaLnL", "dof", "lrt", "fdr")


full_results_df <- as.data.frame(matrix(unlist(x$FullResults), ncol=9, byrow=T))
names(full_results_df) <- c("site", "ssfParams", "ssfLnL", "lssfParams", "lssfLnL", "deltaLnL", "dof", "lrt", "fdr")



