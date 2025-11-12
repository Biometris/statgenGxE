# Extract residuals.

Extract the residuals for the best model.

## Usage

``` r
# S3 method for class 'varCov'
residuals(object, ...)
```

## Arguments

- object:

  An object of class varCov

- ...:

  Not used.

## Value

A data.frame with residuals.

## See also

Other varCov: [`fitted.varCov()`](fitted.varCov.md),
[`gxeVarCov()`](gxeVarCov.md), [`plot.varCov()`](plot.varCov.md),
[`report.varCov()`](report.varCov.md)

## Examples

``` r
# \donttest{
## Select the best variance-covariance model using asreml for modeling.
if (requireNamespace("asreml", quietly = TRUE)) {
  geVarCov <- gxeVarCov(TD = TDMaize, trait = "yld", engine = "asreml")

  ## Extract residuals from the model.
  residVarCov <- residuals(geVarCov)
  head(residVarCov)
  }
# }
```
