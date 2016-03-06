Update
==================

Code and documentation for a updated release to the Omic Kriging R package. Focused primarily on improving efficiency, improving ease of use, and reducing dependencies.

[Link to latest CRAN release](http://cran.r-project.org/web/packages/OmicKriging/index.html)
[Link to github version](https://github.com/hakyimlab/OmicKriging)

## Description

This package provides functions to generate a correlation matrix from a genetic dataset and to use this matrix to predict the phenotype of an individual by using the phenotypes of the remaining individuals through kriging. Kriging is a geostatistical method for optimal prediction or best unbiased linear prediction. It consists of predicting the value of a variable at an unobserved location as a weighted sum of the variable at observed locations. Intuitively, it works as a reverse linear regression: instead of computing correlation (univariate regression coefficients are simply scaled correlation) between a dependent variable Y and independent variables X, it uses known correlation between X and Y to predict Y.

### Authors and Contributors: ###
  * Hae Kyung Im
  * Heather E. Wheeler
  * Keston Aquino Michaels
  * Vassily Trubetskoy
