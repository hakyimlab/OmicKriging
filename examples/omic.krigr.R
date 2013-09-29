 ## Refactor the omic kriging using the GDS format

## dependencies
library(gdsfmt)
library(SNPRelate)
library(Rcpp)
library(RcppEigen)
library(inline)
library(doMC)

## globals (can be factored out into high level config)
gdsFile <- "examples/test.gds"
bedFile <- "data/T1DCC.subset.bed"
bimFile <- "data/T1DCC.subset.bim"
famFile <-"data/T1DCC.subset.fam"
grmDataFile <- "data/T1DCC.subset.GRM.Rdata"
phenoFile <- "data/T1DCC.pheno"
pheno.name <- "PHENO"
ncore <- 4

## load some functions for testing -- ordinarily you will simple load the package
source('R/computeGeneRelMat.R')
source('R/dataInput.R')
source('R/computePCA.R')
source('R/krigrCrossValidation.R')

## load genetic data
load_gene_data(bedFile, bimFile, famFile, gdsFile)

## load phenotype data
pheno <- load_sample_data(phenoFile, main.pheno = pheno.name)

## calculate the genetic relatedness matrix
grm <- make_grm(gdsFile = gdsFile)

## calculate principal components
pca <- make_PCs(gdsFile, n.core = ncore, n.top = 2)


## n-fold parallel cross validation
result <- krigr_cross_validation(n.cores = ncore,
            corlist = list(grm),
            pheno.df = pheno,
            pheno.name = pheno.name)


## summarize results: prediction and performance functions from ROCR


