# Extract variance components

Extract variance components from an object of class varComp.

## Usage

``` r
vc(varComp)
```

## Arguments

- varComp:

  An object of class varComp.

## Value

A data.frame with variance components and standard errors for the random
components in the fitted model.

## See also

Other Mixed model analysis: [`CRDR()`](CRDR.md),
[`correlations()`](correlations.md), [`diagnostics()`](diagnostics.md),
[`gxeVarComp()`](gxeVarComp.md), [`herit()`](herit.md),
[`plot.varComp()`](plot.varComp.md),
[`predict.varComp()`](predict.varComp.md)

## Examples

``` r
## Fit a mixed model.
geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")

## Extract variance components.
vc(geVarComp)
#>           Component
#> genotype    6670.93
#> residuals  12446.94
```
