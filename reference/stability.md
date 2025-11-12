# S3 class stability

Function for creating objects of S3 class stability.  
[`print`](https://rdrr.io/r/base/print.html),
[`summary`](https://rdrr.io/r/base/summary.html),
[`plot`](https://rdrr.io/r/graphics/plot.default.html) and
[`report`](https://rdrr.io/pkg/statgenSTA/man/report.html) methods are
available.

## Usage

``` r
createStability(superiority = NULL, static = NULL, wricke = NULL, TD, trait)
```

## Arguments

- superiority:

  A data.frame containing values for the cultivar-superiority measure of
  Lin and Binns.

- static:

  A data.frame containing values for Shukla's stability variance.

- wricke:

  A data.frame containing values for Wricke's ecovalence.

- trait:

  A character string indicating the trait that has been analyzed.
