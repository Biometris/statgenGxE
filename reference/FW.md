# S3 class FW

Function for creating objects of S3 class FW (Finlay-Wilkinson).  
[`print`](https://rdrr.io/r/base/print.html),
[`summary`](https://rdrr.io/r/base/summary.html),
[`plot`](https://rdrr.io/r/graphics/plot.default.html) and
[`report`](https://rdrr.io/pkg/statgenSTA/man/report.html) methods are
available.

## Usage

``` r
createFW(
  estimates,
  anova,
  envEffs,
  trait,
  nGeno,
  nEnv,
  TD,
  fittedGeno,
  tol,
  iter
)
```

## Arguments

- estimates:

  A data.frame containing the estimated values.

- anova:

  A data.frame containing anova scores of the FW analysis.

- envEffs:

  A data.frame containing the environmental effects.

- trait:

  A character value indicating the analysed trait.

- nGeno:

  A numerical value containing the number of genotypes in the analysis.

- nEnv:

  A numerical value containing the number of environments in the
  analysis.

- TD:

  The object of class [`TD`](https://rdrr.io/pkg/statgenSTA/man/TD.html)
  on which the analysis was performed.

- fittedGeno:

  The fitted values for the genotypes.

- tol:

  A numerical value containing the tolerance used during the analysis.

- iter:

  A numerical value containing the number of iterations for the analysis
  to converge.

## See also

[`plot.FW`](plot.FW.md), [`report.FW`](report.FW.md)
