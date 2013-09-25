### Refactor the omic kriging using the GDS format

## all dependencies
library(gdsfmt)
library(SNPRelate)
library(Rcpp)
library(RcppEigen)
library(inline)

## functions


## globals to be factored out into high level config
gdsFile <- "wtcct1d.gds"
bedFile <- "data/T1DCC.bed"
bimFile <- "data/T1DCC.bim"
famFile <-"data/T1DCC.fam"

## INPUTS 
## convert plink file to GDS
snpgdsBED2GDS(bedFile, famFile, bimFile, gdsFile)
## open GDS file and pull genotype matrix
genofile <- openfn.gds(gdsFile)
X <- read.gdsn(index.gdsn(genofile, "genotype"))
Y <- read.gdsn(index.gdsn(genofile, "sample.annot"))$phenotype

## compute GRM matrix (or read in?)
source('src/rcppcormat.r')
grm <- rcppcormat(t(X))

## cross validation


