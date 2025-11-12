# Extract residuals.

Extract the residuals for the fitted AMMI model.

## Usage

``` r
# S3 method for class 'AMMI'
residuals(object, ...)
```

## Arguments

- object:

  An object of class AMMI

- ...:

  Not used.

## Value

A data.frame with residuals.

## See also

Other AMMI: [`fitted.AMMI()`](fitted.AMMI.md),
[`gxeAmmi()`](gxeAmmi.md), [`plot.AMMI()`](plot.AMMI.md),
[`report.AMMI()`](report.AMMI.md)

## Examples

``` r
## Run AMMI analysis on TDMaize.
geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")

## Extract residuals.
residAmmi <- residuals(geAmmi)
head(residAmmi)
#>   trial genotype   residual
#> 1 HN96b     G001 -54.009646
#> 2 HN96b     G002 -20.533206
#> 3 HN96b     G003 -16.946488
#> 4 HN96b     G004  -7.347288
#> 5 HN96b     G005  -2.275921
#> 6 HN96b     G006  97.747967
```
