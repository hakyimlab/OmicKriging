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
phenoFile <- "data/T1DCC.pheno"
pheno.name <- "PHENO"
ncore <- 4


## load genetic data
source('R/computeGeneRelMat.R')
grm <- make_grm(bedFile = bedFile,
        bimFile = bimFile,
        famFile = famFile,
        gdsFile = gdsFile)

## load phenotype data
pheno.df <- read.table(phenoFile, header=T)
## subset to individuals in GRM and match the sample order
sub.pheno.df <- subset( pheno.df, pheno.df$IID %in% colnames(grm) )
pheno.df <- sub.pheno.df
pheno.df[pheno.name] <- as.numeric(unlist(pheno.df[pheno.name]))
rownames(pheno.df) <- pheno.df$IID

## n-fold parallel cross validation
source('R/krigrCrossValidation.R')
result <- krigr_cross_validation(n.cores = ncore,
            corlist = list(grm),
            pheno.df = pheno.df,
            pheno.name = pheno.name)


## summarize results: prediction and performance functions from ROCR


