


These files can be used to recreate some of the results reported in

Tamuri AU, dos Reis M, Hay AJ, Goldstein RA (2009) Identifying Changes in
Selective Constraints: Host Shifts in Influenza. PLoS Comput Biol 5(11): e1000564
doi:10.1371/journal.pcbi.1000564

Please cite this work if you use any of this material in your own publications.

You'll need a recent version of the JRE to run the program. The JDK and the
'ant' build tool if you want to use the build.xml file to build files from
source. However, a binary is included.

Description of directories:
build - Used by ant for compiled source files
dist - The compiled tdg09.jar library
etc - Protein sequence alignments and tree files for influenza genes
lib - Libraries used by tdg09
src - Java source files for tdg09

To run the program, select the appropriate (matching) alignment and tree files
and pass these as command-line arguments to the run.sh program. For example,

./run.sh etc/H1.faa etc/H1.tree

for the H1 gene. Or

./run.sh etc/MP.faa etc/MP.tree

for the MP gene.





Loading the output in R:

library(yaml)
tdg_output <- yaml.load_file('tdg_output.txt')

summary(tdg_output)

lrt_results_df <- as.data.frame(matrix(unlist(tdg_output$LrtResults), ncol=5, byrow=T))
names(lrt_results_df) <- c("site", "deltaLnL", "dof", "lrt", "fdr")


full_results_df <- as.data.frame(matrix(unlist(tdg_output$FullResults), ncol=9, byrow=T))
names(full_results_df) <- c("site", "ssfParams", "ssfLnL", "lssfParams", "lssfLnL", "deltaLnL", "dof", "lrt", "fdr")


# plotting full results
x <- as.numeric(levels(full_results_df$fdr)[full_results_df$fdr])
x[which(is.na(x))] <- 1
plot(1 - x, ty='h', lwd=2)
abline(h = 1 - 0.05, lty='dashed')
points(which(x < 0.05), 1 - x[x < 0.05], pch=20)

1 115 229 343 457
par(mfrow=c(5,1), mar=c(0,0,0,0))
plot(1 - x, ty='h', lwd=2, xlim=c(1, 114), main='', xlab='')
plot(1 - x, ty='h', lwd=2, xlim=c(115, 228), main='', xlab='')
plot(1 - x, ty='h', lwd=2, xlim=c(229, 342), main='', xlab='')
plot(1 - x, ty='h', lwd=2, xlim=c(343, 456), main='', xlab='')
plot(1 - x, ty='h', lwd=2, xlim=c(457, 566), main='', xlab='')

