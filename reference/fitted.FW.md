# Extract fitted values.

Extract the fitted values for a fitted Finlay-Wilkinson model.

## Usage

``` r
# S3 method for class 'FW'
fitted(object, ...)
```

## Arguments

- object:

  An object of class FW

- ...:

  Not used.

## Value

A data.frame with fitted values.

## See also

Other Finlay-Wilkinson: [`gxeFw()`](gxeFw.md),
[`plot.FW()`](plot.FW.md), [`report.FW()`](report.FW.md),
[`residuals.FW()`](residuals.FW.md)

## Examples

``` r
## Run Finlay-Wilkinson analysis.
geFW <- gxeFw(TD = TDMaize, trait = "yld")
#> Warning: ANOVA F-tests on an essentially perfect fit are unreliable

## Extract fitted values.
fitFW <- fitted(geFW)
head(fitFW)
#>   trial genotype fittedValue seFittedValue
#> 1 HN96b     G001    542.0873      36.14588
#> 2 IS92a     G001    735.8210      43.17104
#> 3 IS94a     G001    467.0806      36.28094
#> 4 LN96a     G001    175.2519      50.51782
#> 5 LN96b     G001     59.5040      59.74929
#> 6 NS92a     G001   1246.5479      85.77817
```
