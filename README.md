
why should they use how should they use how do they get it

<!-- README.md is generated from README.Rmd. Please edit that file -->

# dynUGENE

<!-- badges: start -->

<!-- badges: end -->

dynUGENE is an R package for the inference, simulation, and
visualization of gene regulatory network dynamics from time-series
expression data.

## Description

Purpose, how it add/improves current work, uniqueness

## Installation

<!-- You can install the released version of dynUGENE from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("dynUGENE") -->

<!-- ``` -->

Install the development version with:

``` r
install.packages("devtools")
devtools::install_github("tianyu-lu/dynUGENE", build_vignettes = TRUE)
library("dynUGENE")
```

## Overview

``` r
#ls(package:dynUGENE)
#data(package="dynUGENE")
```

For detailed tutorials, see the vignette here:

``` r
browseVignettes("dynUGENE")
#> No vignettes found by browseVignettes("dynUGENE")
```

## Contributions

describe them, how it builds off of existing functions, which parts are
yours

## References

## Acknowledgements

This package was developed as part of an assessment for 2020 BCB410H:
Applied Bioinformatics, University of Toronto, Toronto, CANADA.

## Example

Load and do inference with repressilator data:

``` r
library(dynUGENE)
## load repressilator dataset

## infer network

## visualize inferred network

## simulate the inferred network dynamics
```
