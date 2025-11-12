# Plot function for class varComp

A plot is created of either the standard deviations of each of the terms
in the fitted model or the percentage of variance explained by each of
the terms in the fitted model. Also the degrees of freedom for each of
the terms is shown in the plot.

## Usage

``` r
# S3 method for class 'varComp'
plot(x, ..., plotType = c("sd", "percVar"), title = NULL, output = TRUE)
```

## Arguments

- x:

  An object of class varComp

- ...:

  Not used.

- plotType:

  A character string. Either "sd" to plot the standard deviation of the
  variance components, or "percVar" to plot the percentage of variance
  explained by each variance component.

- title:

  A character string used a title for the plot.

- output:

  Should the plot be output to the current device? If `FALSE` only a
  ggplot object is invisibly returned.

## Value

A ggplot object is invisibly returned.

## See also

Other Mixed model analysis: [`CRDR()`](CRDR.md),
[`correlations()`](correlations.md), [`diagnostics()`](diagnostics.md),
[`gxeVarComp()`](gxeVarComp.md), [`herit()`](herit.md),
[`predict.varComp()`](predict.varComp.md), [`vc()`](vc.md)

## Examples

``` r
## Fit a mixed model.
geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")

## Plot the standard deviations.
plot(geVarComp)

## Plot the percentage of variance explained.
plot(geVarComp, plotType = "percVar")

```
