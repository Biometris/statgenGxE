# Get diagnostics for an object of class varComp

Get the diagnostics for the model fitted. This will print a table of
combinations missing in the data. For each random factor in the model a
table is printed.

## Usage

``` r
diagnostics(varComp)
```

## Arguments

- varComp:

  An object of class varComp.

## Value

A list of tables is invisibly returned.

## See also

Other Mixed model analysis: [`CRDR()`](CRDR.md),
[`correlations()`](correlations.md), [`gxeVarComp()`](gxeVarComp.md),
[`herit()`](herit.md), [`plot.varComp()`](plot.varComp.md),
[`predict.varComp()`](predict.varComp.md), [`vc()`](vc.md)

## Examples

``` r
## Fit a mixed model.
geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")

## Display diagnostics.
diagnostics(geVarComp)
#> No missing combinations for genotype x trial.
#> 
```
