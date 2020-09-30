Update
==================

Code and documentation for a updated release to the Omic Kriging R package. Focused primarily on improving efficiency, improving ease of use, and reducing dependencies.

[Link to latest CRAN release](https://cran.r-project.org/package=OmicKriging)

## Description

This package provides functions to generate a correlation matrix from a genetic dataset and to use this matrix to predict the phenotype of an individual by using the phenotypes of the remaining individuals through kriging. Kriging is a geostatistical method for optimal prediction or best unbiased linear prediction. It consists of predicting the value of a variable at an unobserved location as a weighted sum of the variable at observed locations. Intuitively, it works as a reverse linear regression: instead of computing correlation (univariate regression coefficients are simply scaled correlation) between a dependent variable Y and independent variables X, it uses known correlation between X and Y to predict Y.
More updated versions can be found [here](https://github.com/hakyimlab/OmicKriging)

## Install from Github

```
install.packages("devtools")

library(devtools)

install_github("hakyimlab/omickriging")
```

## Tutorial

Find tutorial [here](docs/Tutorial-OmicKriging.pdf)


### Authors and Contributors: ###
  * Hae Kyung Im
  * Heather E. Wheeler
  * Keston Aquino Michaels
  * Vassily Trubetskoy

Please cite:
H. E. Wheeler, K. Aquino-Michaels, E. R. Gamazon, V. V. Trubetskoy, M. E. Dolan, R. S. Huang, N. J. Cox, and H. K. Im, “Poly-Omic Prediction of Complex Traits: OmicKriging.,” Genetic epidemiology, vol. 38, no. 5, pp. 402–415, May 2014.
http://onlinelibrary.wiley.com/doi/10.1002/gepi.21808/abstract
