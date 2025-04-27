medsenssurv
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# medsenssurv

<!-- badges: start -->
<!-- badges: end -->

The goal of medsenssurv is to perform sensitivity analysis for
unmeasured confounder in survival mediation analysis

## Installation

You can install the development version of medsenssurv from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Guo-yi-y/medsenssurv")
```

## Example

This is a basic example which shows you how to perform this sensitivity
analysis on a simulated dataset. This will generate a list containing
data frame of sensitivity analysis and 3D scatter plot

``` r
library(medsenssurv)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(coxed)
```

    ## Warning: package 'coxed' was built under R version 4.4.2

    ## Loading required package: rms

    ## Warning: package 'rms' was built under R version 4.4.2

    ## Loading required package: Hmisc

    ## Warning: package 'Hmisc' was built under R version 4.4.2

    ## 
    ## Attaching package: 'Hmisc'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     src, summarize

    ## The following objects are masked from 'package:base':
    ## 
    ##     format.pval, units

    ## Loading required package: survival

    ## Loading required package: mgcv

    ## Loading required package: nlme

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

    ## This is mgcv 1.9-1. For overview type 'help("mgcv-package")'.

``` r
library(htmlwidgets)
```

    ## Warning: package 'htmlwidgets' was built under R version 4.4.1

``` r
set.seed(123)
data_pas = data.frame(covx1 = rnorm(100), covx2 = rnorm(100), U = rbinom(100, 1, 0.5)) %>%
  mutate(
    px = 1/(1+exp(-(0.1*covx1+(-0.1*covx2)+ 0.5*U+0.05))),
    x = unlist(sapply( px , function(p) rbinom(1, 1, p)  )),
    pm = 1/(1+exp(-(0.05*covx1+0.3*covx2+ 1.2*x+ 0.2*U+0.1))),
    m = unlist(sapply( pm , function(p) rbinom(1, 1, p)  ))
  )

simdata <- sim.survdata(N=100, X = data_pas[, c("U","x","m", "covx1","covx2")], T =365*3,  
                        beta=c(0.3, 0.2, 0.5, 0.1, 0.2), fixed.hazard = T)


surv_data = simdata[[1]]
surv_data$failed = as.numeric(surv_data$failed)
survsens_res <- medsurvsens(surv_data, outcome = "failed", time = "y", m = "m", t = "x", covariates = c("covx1", "covx2"),  scale.m = "binary", scale.t = "binary",
                       range.b.y = c(-1, 1), range.b.m = c(-1, 1), range.b.t = c(-1, 1), parallel = T, B=100)
```

Generated 3D scatter plot can’t be viewed in R. You can save it as a
html object and open it with browser on your computer

``` r
saveWidget(survsens_res$p_nie, file = "nie.html", selfcontained = TRUE)
saveWidget(survsens_res$p_nde, file = "nde.html", selfcontained = TRUE)
```

You can also save the data frame into your computer and upload it to
“<http://127.0.0.1:5000/test>”

Reference
