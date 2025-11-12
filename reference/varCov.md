# S3 class varCov

Function for creating objects of S3 class varCov.  
[`print`](https://rdrr.io/r/base/print.html),
[`summary`](https://rdrr.io/r/base/summary.html),
[`plot`](https://rdrr.io/r/graphics/plot.default.html) and
[`report`](https://rdrr.io/pkg/statgenSTA/man/report.html) methods are
available.

## Usage

``` r
createVarCov(STA, choice, summary, vcov, criterion, engine, dat, trait)
```

## Arguments

- STA:

  An object of class STA, the best fitted model.

- choice:

  A character string indicating the best fitted model.

- summary:

  A data.frame with a summary of the fitted models.

- vcov:

  The covariance matrix of the best fitted model.

- criterion:

  A character string indicating the goodness-of-fit criterion used for
  determinening the best model.

- engine:

  A character string containing the engine used for the analysis.
