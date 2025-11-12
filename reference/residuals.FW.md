# Extract residuals.

Extract the residuals for a fitted Finlay-Wilkinson model.

## Usage

``` r
# S3 method for class 'FW'
residuals(object, ...)
```

## Arguments

- object:

  An object of class FW

- ...:

  Not used.

## Value

A data.frame with residuals.

## See also

Other Finlay-Wilkinson: [`fitted.FW()`](fitted.FW.md),
[`gxeFw()`](gxeFw.md), [`plot.FW()`](plot.FW.md),
[`report.FW()`](report.FW.md)

## Examples

``` r
## Run Finlay-Wilkinson analysis.
geFW <- gxeFw(TD = TDMaize, trait = "yld")
#> Warning: ANOVA F-tests on an essentially perfect fit are unreliable

## Extract residuals.
residFW <- residuals(geFW)
head(residFW)
#>   trial genotype     residual
#> 1 HN96b     G001 -114.9126828
#> 2 HN96b     G002  106.8592121
#> 3 HN96b     G003  -70.5527316
#> 4 HN96b     G004   -0.9363679
#> 5 HN96b     G005   -7.1086436
#> 6 HN96b     G006  110.9634426
```
