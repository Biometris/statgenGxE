# S3 class varComp

Function for creating objects of S3 class varComp.  
[`print`](https://rdrr.io/r/base/print.html),
[`summary`](https://rdrr.io/r/base/summary.html), and
[`plot`](https://rdrr.io/r/graphics/plot.default.html) methods are
available.

## Usage

``` r
createVarComp(
  fitMod,
  modDat,
  trait,
  nestingFactor,
  useLocYear,
  useRegionLocYear,
  fullRandVC,
  aovFullFixedMod,
  engine,
  confoundVars,
  diagTabs
)
```

## Arguments

- fitMod:

  A fitted variance components model.

- modDat:

  A data.frame containing the data used in fitting the model.
