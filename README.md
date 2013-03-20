# TdG09 program: identifying changes in selective constraints

This program is an implementation of the model described in:

Tamuri AU, dos Reis M, Hay AJ, Goldstein RA (2009) Identifying Changes in
Selective Constraints: Host Shifts in Influenza. _PLoS Comput Biol_ 5(11): e1000564. [doi:10.1371/journal.pcbi.1000564](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000564)

This phylogenetic model uses site-specific amino acid frequencies to distinguish patterns of substitution between two (or more) lineages or groups of  taxa. We first estimate the site-specific amino acid frequencies at a given location in a protein alignment assuming that the pattern of substitution is homogeneous across all branches. We then estimate multiple sets of site-specific frequencies, allowing them to differ among different branches, producing a non-homogeneous model of evolutionary change. Using statistical tests, we then see whether the non-homogeneous model provides a significantly better fit to the data than the homogenous model.

## Tutorial

1. **Install Java:** The program requires a recent version of the Java Runtime (JRE), which, if not already installed, can be downloaded from [Oracle](http://www.oracle.com/technetwork/java/javase/downloads/index.html). Linux packages are usually available in the distribution's repository (e.g. `sudo apt-get install openjdk-7-jre` for Debian, Ubuntu etc. distributions).

2. **Download the program:** The latest version of the program is available for download from the [repository](https://github.com/tamuri/tdg09/tags). The download includes a compiled binary as well as a `build.xml` to compile from sources using [ant](http://ant.apache.org/). Unzip the download and check that the program works:
 
    ```
    $ ls -F
    README.md  build.xml  dist/  etc/  lib/  src/
    
    $ java -cp dist/tdg09.jar tdg09.Analyse
    Usage: java -cp tdg09.jar tdg09.Analyse [options]
     Options:
     * -alignment
          Alignment in PHYLIP format
     * -groups
          Group labels to partition tree e.g. Av Hu
       -threads
          Number of threads to use
          Default: 1
     * -tree
          Tree in NEWICK format
    ```
    
3. **Preparing your data:** The program requires a protein sequence alignment in [PHYLIP format](http://www.bioperl.org/wiki/PHYLIP_multiple_alignment_format) and the corresponding tree in [Newick format](http://en.wikipedia.org/wiki/Newick_format). The tree branch lengths should be optimised for an amino acid model, such as WAG, using a program such as [RAxML](http://sco.h-its.org/exelixis/software.html) or [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html). Most importantly, sequence names must be prefixed by a two-letter identifier that is used to indicate it's lineage/grouping. For example, to identify changes in selective constraints between avian and human flu viral proteins, we can prefix every sequence name with 'Av' or 'Hu' to indicate its lineage:

    ```
    $ cat etc/H1.faa
    434 566
    Hu_HA_AAX56530_H1N1      MKAKLLVLLCAFTATYADTI...
    Hu_HA_AAY78939_H1N2      MKVKLLILLCTFTATYADTI...
    Av_HA_ABB19507_H1N6      MEAKLFVLFCTFTVLKADTI...
    Av_HA_ABB19518_H1N1      MEAKLFVLFCTFTALKADTI...
    ...
    ```
Both the sequence alignment and the tree must follow this convention. 

    ```
    $ cat etc/H1.tree
    ((Av_HA_ABB19607_H1N1:0.0852880,Av_HA_ABG88212_H1N1:0.1036700):
    0.0248320,((Av_HA_ABB19618_H1N1:...
    ```
This is the only way by which the tdg09 program determines which of the different non-homogeneous models a particular tree branch should use. Example data sets of flu viral proteins (used in Tamuri *et al.* 2009) are included in the `etc/` directory.

4. **Running the program:** The command-line options for the program are:
    + -alignment : the sequence alignment file in PHYLIP format e.g. `etc/H1.faa`
    + -tree : the tree file in Newisk format e.g. `etc/H1.tree`
    + -groups : the two-letter identifiers used to partition the sequences e.g. `Av Hu`
    + -threads : specify the numbers of CPU cores/threads to utilise e.g. `2`
    
    The program prints messages to standard out, so this should be captured using `>` or `tee`. We are now ready to run the program:
    
    ```
    $ java -cp dist/tdg09.jar tdg09.Analyse -alignment etc/H1.faa \
    -tree etc/H1.tree -groups Av Hu -threads 2 > H1_out.txt
    ```
    or if you have 'tee' installed:
    
    ```
    $ java -cp dist/tdg09.jar tdg09.Analyse -alignment etc/H1.faa \
    -tree etc/H1.tree -groups Av Hu -threads 2 | tee H1_out.txt
    ```

5. **Inspect the results:** In this example, program output is captured in `H1_out.txt`:

    ```
    $ cat H1_out.txt
    StartTime: 2013-03-19 13:44:10.129
    WorkingDirectory: /Users/Tester/Documents/tdg09
    Options: -alignment etc/H1.faa -tree etc/H1.tree -groups Av Hu -threads 2 
    
    TreeFile: /Users/Tester/Documents/tdg09/etc/H1.tree
    AlignmentFile: /Users/Tester/Documents/tdg09/etc/H1.faa

    Alignment:
      SequenceCount: 434
      SiteCount: 566

    Groups: [Av, Hu]
    
    # The internal nodes of the tree are not labelled. Labelling...
    # Node 432 from [Av, Hu] resolved
    # Node 432 from [Av, Hu] resolved
    # Assuming that root of tree is in group [Av]
    # Switching from group [Av] to [Hu] at branch 432..431
    
    LabelledTree: >
        (((((((((((Av_HA_ABB19607_H1N1:0.0852880,Av_HA_ABG88212_H1N1:0.1036700)
        Av:0.0248320,((Av_HA_ABB19618_H1N1:0.0910450,Av_HA_ABG88201_H1N1:0.0810480)
        ...     
		Hu:0.0618410)Hu:0.0978830)Hu:0.0436020)Hu:0.0482480)Hu:0.0255390)Hu:0.0362960)
		Hu:0.0696670)Hu:0.0620790)Hu:0.0268645);   
    ```   
The output contains a "LabelledTree" that shows the inferred lineage for each ancestral node. This should be checked in a tree viewing program (such as [Dendroscope](http://www-ab.informatik.uni-tuebingen.de/software/dendroscope)) to make sure that the lineages are correct. If not, they can be modified in and the analysis can be rerun with the new, custom-labelled, tree. The output continues:

    ```    
    # 2013-03-19 22:44:10.33 - site 1 complete.
    # 2013-03-19 22:44:11.366 - site 3 complete.
    # 2013-03-19 22:44:11.366 - site 4 complete.    
    
    ...
    
    LrtResults:
    #   Site,  delta lnL,  dof, LRT,       FDR
    - [  204,  21.850690,  3,   0.0000000, 0.0000003 ]
    - [  169,  11.528011,  1,   0.0000016, 0.0001542 ]
	- [  289,  10.733096,  2,   0.0000218, 0.0008550 ]
	- [  252,  10.927903,  2,   0.0000180, 0.0008796 ]
	- [    9,   8.225516,  1,   0.0000499, 0.0008895 ]
	- [  300,  15.027774,  5,   0.0000144, 0.0009396 ]
	- [   62,   8.261257,  1,   0.0000481, 0.0009423 ]
	- [  303,   8.261836,  1,   0.0000480, 0.0010463 ]
	- [  239,   9.647869,  2,   0.0000646, 0.0010545 ]
	- [  315,   9.483944,  2,   0.0000761, 0.0011468 ]
	- [  253,   8.262439,  1,   0.0000480, 0.0011764 ]

    ...
    
    FullResults:
    # Site, WAG+ssF params, WAG+ssF lnL, WAG+lssF params, WAG+lssF params, delta lnL, dof, LRT, FDR
	- [    1,  NA,         NA, NA,          NA,         NA, NA,        NA,        NA ]
	- [    2,  3,  -21.488149,  5,  -13.378828,   8.109321,  2, 0.0003007, 0.0025627 ]
	- [    3,  2,  -48.561153,  3,  -47.763333,   0.797820,  1, 0.2065223, 0.3489515 ]
	- [    4,  2,  -24.137478,  3,  -23.407924,   0.729554,  1, 0.2270721, 0.3708845 ]
	- [    5,  2,  -12.611075,  3,  -12.456865,   0.154209,  1, 0.5786522, 0.6593943 ]
	- [    6,  3,  -23.856049,  5,  -22.761440,   1.094609,  2, 0.3346706, 0.5084918 ]
	- [    7,  2,  -44.839152,  3,  -43.761797,   1.077355,  1, 0.1421333, 0.2509740 ]
	- [    8,  NA,         NA, NA,          NA,         NA, NA,        NA,        NA ]   

    ...
  
 	```
 	Of interest are the "LrtResults" and "FullResults" tables. 
 	
 	The LrtResults lists polymorphic sites (on which the non-homogeneous model was estimated) and orders them by the false discovery rate (a correction on the likelihood ratio test P-value required by multiple hypothesis testing). These are sites for which the non-homogeneous model provides a statistically significant improvement over the homogenous model, indicating that the patterns of substitution are different between the two groups/lineages. 
 	
 	The FullResults table lists further results from all sites. This includes the log-likelihood for the WAG+ssF (site-specific frequencies or *homogeneous model*) and WAG+lssF (lineage and site-specific frequencies or *non-homogeneous model*). Conserved locations are not analysed, so their entries are 'NA'.
 	
6. **Analysing the results**: The output file is in [YAML](http://www.yaml.org/) format, which means it can be read by any other programming language that has a YAML parsing library. Here we show an example of analysing the results using the programming language [R](http://www.r-project.org/). The code is available in the `etc/` directory of the download. Start the R console, install and load the *yaml* library, then load the output file using the `yaml.load_file` function:

    ```
    $ R
	R version 2.15.3 (2013-03-01) -- "Security Blanket"
    ...
    > install.packages("yaml")
    Installing package(s) into ‘/Users/Tester/Library/R/2.15/library’
    (as ‘lib’ is unspecified)
    ...
    > library(yaml)
    > out <- yaml.load_file(input='/Users/Tester/Documents/tdg09/H1_out.txt')
    > summary(out)
	                   Length Class  Mode
	StartTime            1    -none- character
	WorkingDirectory     1    -none- character
	Options              1    -none- character
	TreeFile             1    -none- character
	AlignmentFile        1    -none- character
	Alignment            2    -none- list
	Groups               2    -none- character
	LabelledTree         1    -none- character
	LrtResults         196    -none- list
	FullResults        566    -none- list
	ConservedPositions   2    -none- list
	SiteResults        566    -none- list
	EndTime              1    -none- character
    ```
    
    To convert the R list objects into a type easier to work with, we convert the LrtResults and FullResults objects into data.frames:
    
    ```
    > lrt_results <- as.data.frame(matrix(unlist(out$LrtResults), ncol=5, byrow=T))
	> names(lrt_results) <- c("site", "deltaLnL", "dof", "lrt", "fdr")
	> head(lrt_results)
	  site  deltaLnL dof      lrt       fdr
	1  204 21.850690   3 0.00e+00 0.0000003
	2  169 11.528011   1 1.60e-06 0.0001542
	3  289 10.733096   2 2.18e-05 0.0008550
	4  252 10.927903   2 1.80e-05 0.0008796
	5    9  8.225516   1 4.99e-05 0.0008895
	6  300 15.027774   5 1.44e-05 0.0009396
	> sum(lrt_results$fdr <= 0.05) # how many sites identified with FDR < 0.05?
	[1] 55
	> full_results <- as.data.frame(matrix(unlist(out$FullResults), ncol=9, byrow=T))
	> names(full_results) <- c("site", "ssfParams", "ssfLnL", "lssfParams", "lssfLnL", "deltaLnL", "dof", "lrt", "fdr")
	> head(full_results)
	  site ssfParams     ssfLnL lssfParams    lssfLnL deltaLnL dof       lrt       fdr
	1    1        NA         NA         NA         NA       NA  NA        NA        NA
	2    2         3 -21.488149          5 -13.378828 8.109321   2 0.0003007 0.0025627
	3    3         2 -48.561153          3 -47.763333  0.79782   1 0.2065223 0.3489515
	4    4         2 -24.137478          3 -23.407924 0.729554   1 0.2270721 0.3708845
	5    5         2 -12.611075          3 -12.456865 0.154209   1 0.5786522 0.6593943
	6    6         3 -23.856049          5  -22.76144 1.094609   2 0.3346706 0.5084918
	```
	To produce a plot of FDR values by site:
	
	```
	> fdr <- as.numeric(levels(full_results$fdr)[full_results$fdr]) # FDR column should be numeric
	Warning message:
	NAs introduced by coercion
	> fdr[is.na(fdr)] <- 1.0 # conserved locations implicitly have no evidence of non-homogeneity
	> sites <- out$Alignment$SiteCount
	> plot_ranges <- split(seq(1, sites), cut(seq(1, sites), 5)) # split sites in plot into 5 rows
	> par(mfrow=c(5,1), mar=c(2.0,0.5,0.5,0.5))
	for (p in 1:5) {
		plot(1 - fdr,  
		xlim=c(plot_ranges[[p]][1], tail(plot_ranges[[p]], n=1)), 
		ty='h', lwd=1, main='', xlab='', ylab='', yaxt='n', col="#1B9E77")
		lines(which(fdr <= 0.20), 1 - fdr[fdr <= 0.20], col="#D95F02", ty='h')
		abline(h=0.95, lty='dashed')
		points(which(fdr <= 0.05), 1 - fdr[fdr <= 0.05], pch=20, col="#DE2D26")
	}
	```
	This produces the plot show below. Bars drawn in orange indicate locations with FDR < 0.20, and those locations with FDR < 0.05 have a red dot at their value (we drew the plot with "1 - fdr" values so that smaller FDRs are taller). We can see from this plot a cluster of identified sites at locations 200-210.
	
	![image](./example_plot.png)

