# Predictions based on a fitted varComp model.

Predictions are made based on the fitted model in the varComp object.
These predictions can be at genotype level, at genotype x trial level or
at the level of genotype x nestingFactor. If the model was fitted with
trial as year x location then genotype x trial level becomes genotype x
year x location.

## Usage

``` r
# S3 method for class 'varComp'
predict(object, ..., predictLevel = "genotype")
```

## Arguments

- object:

  An object of class varComp.

- ...:

  Not used.

- predictLevel:

  A character string, the level at which prediction should be made.
  Either "genotype" for prediction at genotype level, "trial" for
  predictions at genotype x trial level, the variable used as nesting
  factor for predictions at the level of genotype x nestingFactor level,
  or one or more of the extra terms used in the model. E.g. c("region",
  "year") for a model fitted with `regionLocationYear = TRUE`.

## Value

A data.frame with predictions.

## See also

Other Mixed model analysis: [`CRDR()`](CRDR.md),
[`correlations()`](correlations.md), [`diagnostics()`](diagnostics.md),
[`gxeVarComp()`](gxeVarComp.md), [`herit()`](herit.md),
[`plot.varComp()`](plot.varComp.md), [`vc()`](vc.md)

## Examples

``` r
## Fit a mixed model.
geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")

## Predictions at genotype level.
predGeno <- predict(geVarComp)
head(predGeno)
#>   genotype predictedValue
#> 1     G001       500.1561
#> 2     G002       479.8132
#> 3     G003       471.3800
#> 4     G004       346.9103
#> 5     G005       466.8594
#> 6     G006       425.9303

## Predictions at genotype x trial level.
predGenoTrial <- predict(geVarComp, predictLevel = "trial")
head(predGenoTrial)
#>   genotype trial predictedValue
#> 1     G001 HN96b       529.1555
#> 2     G001 LN96a       227.7905
#> 3     G001 IS92a       684.1152
#> 4     G001 LN96b       133.7147
#> 5     G001 SS94a       457.4688
#> 6     G001 NS92a      1093.0337
```
