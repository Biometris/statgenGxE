# Form mega environments based on fitted values from an AMMI model

This function fits an AMMI model and then using the fitted values
produces a new factor clustering the trials. This factor is added as a
column megaEnv to the input data. If a column megaEnv already exists
this column is overwritten with a warning.  
  
Mega environments are created by grouping environments based on their
best performing genotype; i.e. environments that share the same best
genotype belong to the same mega environment.

## Usage

``` r
gxeMegaEnv(
  TD,
  trials = names(TD),
  trait,
  method = c("max", "min"),
  byYear = FALSE
)
```

## Arguments

- TD:

  An object of class [`TD`](https://rdrr.io/pkg/statgenSTA/man/TD.html).

- trials:

  A character string specifying the trials to be analyzed. If not
  supplied, all trials are used in the analysis.

- trait:

  A character string specifying the trait to be analyzed.

- method:

  A character string indicating the criterion to determine the best
  genotype per environment, either `"max"` or `"min"`.

- byYear:

  Should the analysis be done by year? If `TRUE` the data is split by
  the variable year, analysis is performed and the results are merged
  together and returned.

## Value

An object of class megaEnv, a list consisting of

- TD:

  An object of class TD, the TD object used as input to the function
  with an extra column megaEnv.

- summTab:

  A data.frame, a summary table containing information on the trials in
  each mega environment.

- trait:

  The trait used for calculating the mega environments.

## References

Atlin, G. N., R. J. Baker, K. B. McRae, and X. Lu. 2000. Selection
Response in Subdivided Target Regions. Crop Sci. 40:7-13.
[doi:10.2135/cropsci2000.4017](https://doi.org/10.2135/cropsci2000.4017)

## See also

Other mega environments: [`plot.megaEnv()`](plot.megaEnv.md),
[`predict.megaEnv()`](predict.megaEnv.md)

## Examples

``` r
## Calculate mega environments for TDMaize.
gemegaEnv <- gxeMegaEnv(TD = TDMaize, trait = "yld")

## Calculate new mega environments based on the genotypes with the lowest
## value per environment.
gemegaEnv2 <- gxeMegaEnv(TD = TDMaize, trait = "yld", method = "min")
```
