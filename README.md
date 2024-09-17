
<!-- README.md is generated from README.Rmd. Please edit that file -->

# statgenGxE

[![](https://www.r-pkg.org/badges/version/statgenGxE)](https://www.r-pkg.org/pkg/statgenGxE)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/statgenGxE)](https://www.r-pkg.org/pkg/statgenGxE)
[![R-CMD-check](https://github.com/Biometris/statgenGxE/actions/workflows/check.yaml/badge.svg)](https://github.com/Biometris/statgenGxE/actions/workflows/check.yaml)
[![codecov](https://codecov.io/gh/Biometris/statgenGxE/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Biometris/statgenGxE)

**statgenGxE** is an R package providing functions for Genotype by
Environment (GxE) analysis for data of plant breeding experiments.

The following types of analysis can be done using statgenGxE:

- Mixed model analysis of GxE table of means
- Finlay-Wilkinson Analysis
- AMMI Analysis
- GGE Analysis
- Identifying mega environments
- Stability measures
- Modeling of heterogeneity of genetic variances and correlations

statgenGxE has extensive options for summarizing and visualizing the
results of the analyses. For a full overview of all options it is best
to read the
[**vignette**](https://biometris.github.io/statgenGxE/articles/statgenGxE.html)

## Installation

- Install from CRAN:

``` r
install.packages("statgenGxE")
```

- Install latest development version from GitHub (requires
  [remotes](https://github.com/r-lib/remotes) package):

``` r
remotes::install_github("Biometris/statgenGxE", ref = "develop", dependencies = TRUE)
```
