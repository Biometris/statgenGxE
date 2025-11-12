# Extract fitted values.

Extract the fitted values for an object of class varCov.

## Usage

``` r
# S3 method for class 'varCov'
fitted(object, ...)
```

## Arguments

- object:

  An object of class varCov

- ...:

  Not used.

## Value

A data.frame with fitted values.

## See also

Other varCov: [`gxeVarCov()`](gxeVarCov.md),
[`plot.varCov()`](plot.varCov.md),
[`report.varCov()`](report.varCov.md),
[`residuals.varCov()`](residuals.varCov.md)

## Examples

``` r
# \donttest{
if (requireNamespace("asreml", quietly = TRUE)) {
  ## Select the best variance-covariance model using asreml for modeling.
  geVarCov <- gxeVarCov(TD = TDMaize, trait = "yld", engine = "asreml")
  ## Extract fitted values from the model.

  fitVarCov <- fitted(geVarCov)
  head(fitVarCov)
  }
# }
```
