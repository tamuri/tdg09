# An example of how to load output from tdg09 program in R and inspect the results
# @author Asif Tamuri
# @version 1.1

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

# We can use the phangorn and ape packages to generate simulated sequences
library(phangorn)
library(ape)

simtree <- read.tree(text = out$LabelledTree)    # Load the tree

# Get the site amino acid equilibrium frequencies for the homogeneous model
homog_f <- as.data.frame(matrix(unlist(out$HomogeneousFrequencies), ncol=21, byrow=T))
names(homog_f) <- c("Site", strsplit("ARNDCQEGHILKMFPSTWYV", split="")[[1]])

# Get the amino acid frequencies for a particular site (say, 204)
site_pos <- which(homog_f$Site == 204)
freqs <- as.numeric(homog_f[site_pos, 2:21])     # Site 204 is row 98
freqs[freqs == 0] <- 0.000000000001       # Add pseudocounts for unobserved amino acids
freqs <- freqs / sum(freqs)               # Normalise

# Synthesize 1000 sites using the WAG+ssF (site-specific frequencies/homogeneous) model
sim_seqs <- simSeq(tree, l=1000, bf=freqs, type="AA", model="WAG", rate=out$HomogeneousRates[site_pos])

# Save the sequence (there is probably a much easier way than this)
write.dna(file='/Users/Tester/Documents/tdg09/sim_site204.txt', x=toupper(as.character(sim_seqs)), format='sequential', nbcol= -1, colsep='')

# We can analyse the simulated data (1000 sites simulated under site 204's homogeneous model)
# and use the Cox test to calculate statistical significance of the non-homogeneous model for the real data.
# 
# For details, see Goldman, N. (1993). Statistical tests of models of DNA substitution. Journal of Molecular Evolution, 36(2), 182–198.
#
# For example, given the original tree in this example (H1.tree) and the simulated data for site 204 (sim_site204.txt),
# we can analyse the simulated data set using tdg09:
#
# java -cp dist/tdg09.jar tdg09.Analyse -tree etc/H1.tree -alignment sim_site204.txt -groups Av Hu -threads 2 > sim204_out.txt
# 
# Once complete, read the results of this analysis into R and calculate the Monte Carlo P-value
 
sim_out <- yaml.load_file(input = '/Users/Tester/Documents/tdg09/sim204_out.txt')         # Load results of simulation
sim_lrt_results <- as.data.frame(matrix(unlist(sim_out$LrtResults), ncol=5, byrow=T))	  # Prepare the likelihood ratio test table (contains delta lnL)
names(sim_lrt_results) <- c("site", "deltaLnL", "dof", "lrt", "fdr")
summary(sim_lrt_results$deltaLnL)                                                         # Have a look at the distribution of delta lnL

# What is the Monte Carlo p-value? = (number of simulated delta lnL < delta lnL estimated from real data) / (1 + number of replicates)
real_data_delta = lrt_results[lrt_results$site == '204', "deltaLnL"]
pvalue <- 1 - sum(sim_lrt_results$deltaLnL < real_data_delta ) / (1 + length(sim_lrt_results$deltaLnL))

# Plot the Monte Carlo distribution of Delta for the Cox test applied to site 204 of the flu virus H1 protein
plot(density(sim_lrt_results$deltaLnL), xlim=c(0,25), yaxs='i', ylim=c(0,0.7), xlab=expression(Delta), ylab='', yaxt='n', main='')
abline(v = real_data_delta, col='red')
text(x = real_data_delta + 1.0, y = 0.5, pos = 3, labels=bquote(paste(delta, " = ", .(real_data_delta))))



