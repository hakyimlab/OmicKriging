OmicKriging
===========


This package provides functions to generate a correlation matrix from a genetic dataset and to use this matrix to predict the phenotype of an individual by using the phenotypes of the remaining individuals through kriging. Kriging is a geostatistical method for optimal prediction or best unbiased linear prediction. It consists of predicting the value of a variable at an unobserved location as a weighted sum of the variable at observed locations. Intuitively, it works as a reverse linear regression: instead of computing correlation (univariate regression coefficients are simply scaled correlation) between a dependent variable Y and independent variables X, it uses known correlation between X and Y to predict Y.






## Performance Improvements ##

### Performance of GDS-based input ###

Using WTCC T1D data set consisting of 4900 individuals and 388548 SNPs.

Basic code:

    f <- function() {
      library(gdsfmt)
      library(SNPRelate)
      ## convert plink file to GDS
      gdsFile <- "wtcct1d.gds"
      snpgdsBED2GDS("data/T1DCC.bed", "data/T1DCC.fam", "data/T1DCC.bim", gdsFile)
      ## open GDS file
      genofile <- openfn.gds(gdsFile)
      g <- read.gdsn(index.gdsn(genofile, "genotype"))
    }

_Results are extremely favorable:_

    > system.time(f())
    SNPRelate: 0.9.16
    Supported by Streaming SIMD Extensions 2 (SSE2).
    Start snpgdsBED2GDS ...
    	open /home/vasya/coxlab_projects/turbo_krigr/data/T1DCC.bed in the SNP-major mode
    	open /home/vasya/coxlab_projects/turbo_krigr/data/T1DCC.fam DONE.
    	open /home/vasya/coxlab_projects/turbo_krigr/data/T1DCC.bim DONE.
    Wed Sep 25 15:43:32 2013 	store sample id, snp id, position, and chromosome.
    	start writing: 4900 samples, 388548 SNPs ...
     	Wed Sep 25 15:43:32 2013	0%
     	Wed Sep 25 15:43:34 2013	100%
    Wed Sep 25 15:43:39 2013 	Done.
       user  system elapsed 
      9.396   2.092  16.432 
      
     
     
### Performance of GDS-based GRM calculation ###
     
Smaller dataset for testing: 4900 samples, 80340 SNPs

    > grm <- make_grm_gds(gdsFile = gdsFile, n.core = ncore)
    > system.time(grm <- make_grm_gds(gdsFile = gdsFile, n.core = ncore))
        user  system elapsed
        952.075   0.509  98.331 



### Performance of irlba-based SVD calculation ###

This particular algorithm computes partial SVDs efficiently. The time depends on the number of top eigenvectors you request, but seems to scale quite favorably.

    > dim(grm)
    [1] 4900 4900
    > system.time(pca <- make_PCs_irlba(grm, n.top = 2))
        user  system elapsed 
        5.013   0.000   5.027
        
    > system.time(pca <- make_PCs_irlba(grm, n.top = 10))
        user  system ela psed 
        7.002   0.014   7.036 

