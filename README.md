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
library(coxed)
library(htmlwidgets)
library(plotly)
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

In the 3D plot, other covariates can be used as benchmarks to help
assess the strength of U. Here, we provide an example including covx1,
covx2, and twice their corresponding values. In practical applications,
you can add benchmark points according to your own needs.

``` r
covar_df = data.frame(covar = c("covx1", "covx2"), x_coef = c(0.1, 0.1), m_coef = c(0.05, 0.3), y_coef = c(0.1, 0.2))

covar_df2 = data.frame(covar = c( "2*covx1", "2*covx2"), x_coef = c(0.2, 0.2), m_coef = c( 0.1, 0.6), y_coef = c(0.2, 0.4))

p_nie_with_covar = survsens_res$p_nie %>%
  add_trace(
    data = covar_df,
    x = ~m_coef, y = ~y_coef, z = ~x_coef,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3, color = "black"),
    showlegend = FALSE
  ) %>%
  add_trace(
    data = covar_df,
    x = ~m_coef, y = ~y_coef, z = ~x_coef,
    type = "scatter3d",
    mode = "text",
    text = ~covar,
    textfont = list(color = "black", size = 13),
    showlegend = FALSE,
    inherit = FALSE
  )%>%
  add_trace(
    data = covar_df2,
    x = ~m_coef, y = ~y_coef, z = ~x_coef,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3, color = "blue"),
    showlegend = FALSE
  ) %>%
  add_trace(
    data = covar_df2,
    x = ~m_coef, y = ~y_coef, z = ~x_coef,
    type = "scatter3d",
    mode = "text",
    text = ~covar,
    textfont = list(color = "black", size = 13),
    showlegend = FALSE,
    inherit = FALSE
  )


p_nde_with_covar = survsens_res$p_nde %>%
  add_trace(
    data = covar_df,
    x = ~m_coef, y = ~y_coef, z = ~x_coef,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3, color = "blue"),
    showlegend = FALSE
  ) %>%
  add_trace(
    data = covar_df,
    x = ~m_coef, y = ~y_coef, z = ~x_coef,
    type = "scatter3d",
    mode = "text",
    text = ~covar,
    textfont = list(color = "black", size = 13),
    showlegend = FALSE,
    inherit = FALSE
  )%>%
  add_trace(
    data = covar_df,
    x = ~m_coef, y = ~y_coef, z = ~x_coef,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3, color = "black"),
    showlegend = FALSE
  ) %>%
  add_trace(
    data = covar_df,
    x = ~m_coef, y = ~y_coef, z = ~x_coef,
    type = "scatter3d",
    mode = "text",
    text = ~covar,
    textfont = list(color = "black", size = 13),
    showlegend = FALSE,
    inherit = FALSE
  )%>%
  add_trace(
    data = covar_df2,
    x = ~m_coef, y = ~y_coef, z = ~x_coef,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 3, color = "blue"),
    showlegend = FALSE
  ) %>%
  add_trace(
    data = covar_df2,
    x = ~m_coef, y = ~y_coef, z = ~x_coef,
    type = "scatter3d",
    mode = "text",
    text = ~covar,
    textfont = list(color = "black", size = 13),
    showlegend = FALSE,
    inherit = FALSE
  )


saveWidget(p_nie_with_covar, file = "nie_with_covar.html", selfcontained = TRUE)
saveWidget(p_nde_with_covar, file = "nde_with_covar.html", selfcontained = TRUE)
```

## Reference

Guo Y, Chen D, Xu X, Zhang Z, Wen Y, Zheng X, Wu Z, Qin X. Sensitivity
Analysis for Unmeasured Confounding in Causal Mediation Analysis With
Survival Outcome. Stat Med. 2026 Apr;45(8-9):e70548.
