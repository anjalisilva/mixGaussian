
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mixGaussian

Mixtures of Multivariate Gaussian Model for Clustering Data

<!-- badges: start -->
<!-- https://www.codefactor.io/repository/github/anjalisilva/testingpackage/issues -->

[![CodeFactor](https://www.codefactor.io/repository/github/anjalisilva/mixgaussian/badge)](https://www.codefactor.io/repository/github/anjalisilva/testingpackage)
[![GitHub
issues](https://img.shields.io/github/issues/anjalisilva/mixGaussian)](https://github.com/anjalisilva/mixGaussian/issues)
[![License](https://img.shields.io/badge/license-MIT-green)](./LICENSE)
![GitHub language
count](https://img.shields.io/github/languages/count/anjalisilva/mixGaussian)
![GitHub commit activity
(branch)](https://img.shields.io/github/commit-activity/y/anjalisilva/mixGaussian/master)
<!-- https://shields.io/category/license --> <!-- badges: end -->

## Description

`mixGaussian` is a simple R package for performing clustering using
mixtures of multivariate Gaussian distributions. Main function
***mixGaussianEM*** permit to carry out model-based clustering.
Information criteria (AIC, BIC, AIC3 and ICL) are offered for model
selection. The shiny implementation of *mixGaussian* is available as
***runMixGaussian***. For more information, see details below.

## Installation

To install the latest version of the package:

``` r
require("devtools")
install_github("anjalisilva/mixGaussian", build_vignettes = TRUE)
library("mixGaussian")
```

To run the Shiny app:

``` r
mixGaussian::runMixGaussian()
```

## Overview

To list all functions available in the package:

``` r
ls("package:mixGaussian")
```

`mixGaussian` contains 6 functions.

1.  ***mixGaussianEM*** for carrying out clustering of data using
    mixtures of multivariate Gaussian distributions via
    expectation-maximization (EM).

2.  ***AICFunction*** for model selection.

3.  ***BICFunction*** for model selection.

4.  ***AIC3Function*** for model selection.

5.  ***ICLFunction*** for model selection.

6.  ***runMixGaussian*** for the shiny implementation of *mixGaussian*.

<div style="text-align:center">

<img src="inst/extdata/ShinyLinePlot.png" alt="ShinyLinePlot" width="750" height="550"/>

Figure: Shiny app for mixGaussian package showing cluster results.

<div style="text-align:left">
<div style="text-align:left">
&#10;

## Details

<div style="text-align:left">

<img src="inst/extdata/MixtureGaussian.png" alt="MixtureGaussian" width="750" height="550"/>

<div style="text-align:left">
<div style="text-align:left">

For more details, see vignette.

## `mixGaussian` Specifics

In `mixGaussian` the parameter and group membership estimation is
carried out using the EM algorithm, because the complete-data consists
of the unobserved group membership labels. To check the convergence of
EM algorithm, a modified version of Aitken’s acceleration criterion as
outlined by B¨ohning et al., 1994 is used. Model selection is performed
using AIC, BIC, AIC3 and ICL. For more details, see vignette.

## Tutorials

For tutorials and plot interpretation, refer to the vignette:

``` r
browseVignettes("mixGaussian")
```

## Citation for Package

``` r
citation("mixGaussian")
```

## References for Package

- [Aitken, A. C. (1926). A series formula for the roots of algebraic and
  transcendental equations. *Proceedings of the Royal Society of
  Edinburgh*](https://www.cambridge.org/core/journals/proceedings-of-the-royal-society-of-edinburgh/article/iiia-series-formula-for-the-roots-of-algebraic-and-transcendental-equations/0CC96A97C8B634E2730F5208E506E6A9)

- [B¨ohning, D., E. Dietz, R. Schaub, P. Schlattmann, and B. Lindsay
  (1994). The distribution of the likelihood ratio for mixtures of
  densities from the one-parameter exponential family. *Annals of the
  Institute of Statistical
  Mathematics*](https://link.springer.com/article/10.1007/BF01720593)

- [Dempster, A. P., N. M. Laird, and D. B. Rubin (1977). Maximum
  likelihood from incomplete data via the EM algorithm. *Journal of the
  Royal Statistical Society: Series
  B*](https://www.ece.iastate.edu/~namrata/EE527_Spring08/Dempster77.pdf)

- For others, refer to help page of inidividual functions via
  `?function` or `help(function)`.

## Maintainer

- Anjali Silva (<anjali@alumni.uoguelph.ca>).

## Contributions

`mixGaussian` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/anjalisilva/mixGaussian/issues).
