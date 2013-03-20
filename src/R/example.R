# An example of how to load output from tdg09 program in R and inspect the results
# Author: Asif Tamuri
# Date: 2013-03-20

# Install the R package 'yaml' to parse the tdg09 output file
install.packages("yaml")
library(yaml)

# Load the tdg09 output file
out <- yaml.load_file(input='/Users/Tester/Documents/tdg09/H1_out.txt')

# Here's what we have
summary(out)

# Prepare the LRT results table
lrt_results <- as.data.frame(matrix(unlist(out$LrtResults), ncol=5, byrow=T))
names(lrt_results) <- c("site", "deltaLnL", "dof", "lrt", "fdr")

# Here's what we have
head(lrt_results)

# How many sites did we find with false discovery rate < 0.05?
sum(lrt_results$fdr <= 0.05)

# Prepare the full results table
full_results <- as.data.frame(matrix(unlist(out$FullResults), ncol=9, byrow=T))
names(full_results) <- c("site", "ssfParams", "ssfLnL", "lssfParams", "lssfLnL", "deltaLnL", "dof", "lrt", "fdr")

# Here's what we have
head(full_results)

# Make a plot of FDR across the sequence alignment
fdr <- as.numeric(levels(full_results$fdr)[full_results$fdr]) # FDR column should be numeric
fdr[is.na(fdr)] <- 1.0 # conserved locations implicitly have no evidence of non-homogeneity

# Split the sequence into 5 chunks
sites <- out$Alignment$SiteCount
plot_ranges <- split(seq(1, sites), cut(seq(1, sites), 5))

# Draw the plot
par(mfrow=c(5,1), mar=c(2.0,0.5,0.5,0.5))
for (p in 1:5) {
    plot(1 - fdr,  xlim=c(plot_ranges[[p]][1], tail(plot_ranges[[p]], n=1)), ty='h', lwd=1, main='', xlab='', ylab='', yaxt='n', col="#1B9E77")
    lines(which(fdr <= 0.20), 1 - fdr[fdr <= 0.20], col="#D95F02", ty='h')
    abline(h=0.95, lty='dashed')
    points(which(fdr <= 0.05), 1 - fdr[fdr <= 0.05], pch=20, col="#DE2D26")
}


