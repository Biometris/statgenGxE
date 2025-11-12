# Extract fitted values.

Extract the fitted values for an object of class AMMI.

## Usage

``` r
# S3 method for class 'AMMI'
fitted(object, ...)
```

## Arguments

- object:

  An object of class AMMI

- ...:

  Not used.

## Value

A data.frame with fitted values.

## See also

Other AMMI: [`gxeAmmi()`](gxeAmmi.md), [`plot.AMMI()`](plot.AMMI.md),
[`report.AMMI()`](report.AMMI.md),
[`residuals.AMMI()`](residuals.AMMI.md)

## Examples

``` r
## Run AMMI analysis on TDMaize.
geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")

## Extract fitted values.
fitAmmi <- fitted(geAmmi)
head(fitAmmi)
#>   genotype trial fittedValue
#> 1     G001 HN96b    602.9904
#> 2     G002 HN96b    386.4668
#> 3     G003 HN96b    557.0535
#> 4     G004 HN96b    335.6527
#> 5     G005 HN96b    493.7241
#> 6     G006 HN96b    433.7480
```
