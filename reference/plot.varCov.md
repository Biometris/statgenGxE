# Plot function for class varCov

Function for plotting a heatmap of the correlation matrix for objects of
class varCov.

## Usage

``` r
# S3 method for class 'varCov'
plot(x, title = paste("Heatmap for model:", x$choice), ..., output = TRUE)
```

## Arguments

- x:

  An object of class varCov

- title:

  A character string used a title for the plot.

- ...:

  Not used.

- output:

  Should the plot be output to the current device? If `FALSE` only a
  ggplot object is invisibly returned.

## See also

Other varCov: [`fitted.varCov()`](fitted.varCov.md),
[`gxeVarCov()`](gxeVarCov.md), [`report.varCov()`](report.varCov.md),
[`residuals.varCov()`](residuals.varCov.md)

## Examples

``` r
# \donttest{
if (requireNamespace("asreml", quietly = TRUE)) {
  ## Select the best variance-covariance model using asreml for modeling.
  geVarCov <- gxeVarCov(TD = TDMaize, trait = "yld", engine = "asreml")

  ## Create a heatmap of the correlation matrix for the best model.
  plot(geVarCov)
  }
# }
```
