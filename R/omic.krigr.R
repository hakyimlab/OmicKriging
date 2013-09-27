### Refactor the omic kriging using the GDS format

## dependencies
library(gdsfmt)
library(SNPRelate)
library(Rcpp)
library(RcppEigen)
library(inline)
library(doMC)

## globals (can be factored out into high level config)
gdsFile <- "test.gds"
bedFile <- "data/T1DCC.subset.bed"
bimFile <- "data/T1DCC.subset.bim"
famFile <-"data/T1DCC.subset.fam"
ncore <- 4
pheno.name <- "phenotype"


## load genetic data
source('R/computeGeneRelMat.R')
grm <- make_grm(bedFile = bedFile,
        bimFile = bimFile,
        famFile = famFile,
        gdsFile = gdsFile)

## load phenotype data
genofile <- openfn.gds(gdsFile)
sample.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
pheno.val <- read.gdsn(index.gdsn(genofile, "sample.annot"))$phenotype
pheno.df <- data.frame(IID=sample.ids)
pheno.df[,pheno.name] <- pheno.val
rownames(pheno.df) <- sample.ids

## n-fold parallel cross validation
source('R/krigrCrossValidation.R')
result <- krigr_cross_validation(n.cores = ncore,
            corlist = list(grm),
            pheno.df = pheno.df,
            pheno.name = pheno.name)


## summarize results: prediction and performance functions from ROCR


