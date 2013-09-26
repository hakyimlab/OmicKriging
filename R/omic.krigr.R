### Refactor the omic kriging using the GDS format

## dependencies
library(gdsfmt)
library(SNPRelate)
library(Rcpp)
library(RcppEigen)
library(inline)
library(doMC)

## globals to be factored out into high level config
gdsFile <- "test.gds"
bedFile <- "data/T1DCC.subset.bed"
bimFile <- "data/T1DCC.subset.bim"
famFile <-"data/T1DCC.subset.fam"
ncore <- 12

source('src/computeGeneRelMat.R')
grm <- make_grm(bedfile = bedFile, bimFile = bimFile, famFile = famFile, gdsFile = gdsFile)

## n-fold parallel cross validation
source('src/krigrCrossValidation.R')
result <- krigr_cross_validation(n.cores = ncore,
            corlist = list(grm),
            pheno.df = pheno.df,
            pheno.name = "test")

## summarize results: prediction and performance functions



