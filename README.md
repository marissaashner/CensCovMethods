
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CensCovMethods: Consistent parameter estimation with censored covariates in R

Marissa C. Ashner and Tanya P. Garcia

## Installation

You can install the development version of CensCovMethods from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
#devtools::install_github("marissaashner/CensCovMethods")
```

## Data Generation

To illustrate the functions in CensCovMethods, let’s simulate a dataset
for a sample of *n* = 500 observations under ( ∼ 40%) random right
censoring. Specifically,

-   $X $
-   $C $
-   $Z $
-   *ϵ* ∼ *N*(0, 1)
-   *Y* = *X* + *Z* + *ϵ*
-   *W* = min (*X*, *C*)
-   *Δ* = *I*(*X* ≤ *C*).

``` r
library(CensCovMethods)
#> Warning: replacing previous import 'numDeriv::hessian' by 'rootSolve::hessian'
#> when loading 'CensCovMethods'

### Generate Data 
n = 500
r = 0.4
X = rnorm(n, 0, 1)
C = rnorm(n, 0, 1)
Z = rnorm(n, 0, 1)
e = rnorm(n, 0, 1)
Y = X + Z + e
W = pmin(X, C)
D = ifelse(X <= C, 1, 0)
data_sample = data.frame(Y, W, D, Z)
```

## Complete Case Analysis

``` r
# cc_estimate = cc_censored(formula = Y ~ beta[1]*W + beta[2] + beta[3]*Z,
#                           cens_ind = "D", 
#                           data = data_sample, 
#                           starting_vals = c(0,0,0), 
#                           par_vec = "beta",
#                           sandwich_se = TRUE)
```

## Inverse Probability Weighting

``` r
# ipw_estimate = ipw_censored(formula = Y ~ beta[1]*W + beta[2] + beta[3]*Z,
#                           cens_ind = "D", 
#                           cens_name = "W",
#                           data = data_sample, 
#                           starting_vals = c(0,0,0), 
#                           par_vec = "beta",
#                           sandwich_se = TRUE,
#                           weight_opt = "MVN", 
#                           weights_cov = "Z", 
#                           weights_threshold = 25,
#                           weight_stabilize = "KM")
```

## Augmented Inverse Probability Weighting

``` r
# aipw_estimate = aipw_censored()
```

## Maximum Likelihood

``` r
# mle_estimate = mle_censored()
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
