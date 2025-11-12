# Mixed model analysis of GxE table of means

This function fits a mixed model best fitting to the data in a TD
object. The exact model fitted is determined by both the structure of
the genotype by environment table of observations and the chosen
parameters.  
  
Six different types of models can be fitted depending on the structure
of the environments in the data. These models are described in the table
below, together with the function parameters used in `gxeVarComp` to fit
the model.

|                                                            |                                                                                                                                                                                                                     |                              |
|------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------|
| Structure of environments                                  | Model                                                                                                                                                                                                               | Function parameters          |
| Environments correspond to trials                          | **trait** = trial + **genotype + genotype:trial**                                                                                                                                                                   |                              |
| Trials form a factorial structure of locations x years     | **trait** = year + location + year:location + **genotype + genotype:year + genotype:location + genotype:year:location**                                                                                             | `locationYear = TRUE`        |
| Trials are nested within year                              | **trait** = year + year:trial + **genotype + genotype:year + genotype:year:trial**                                                                                                                                  | `nestingFactor = "year"`     |
| Trials are nested within locations                         | **trait** = location + location:trial + **genotype + genotype:location + genotype:location:trial**                                                                                                                  | `nestingFactor = "loc"`      |
| Trials correspond to locations within regions across years | **trait** = region + region:location + year + region:year + region:location:year + **genotype + genotype:region + genotype:region:location + genotype:year + genotype:region:year + genotype:region:location:year** | `regionLocationYear = TRUE`  |
| Trials are nested within scenarios                         | **trait** = scenario + scenario:trial + **genotype + genotype:scenario + genotype:scenario:trial**                                                                                                                  | `nestingFactor = "scenario"` |

In the models above the random part of the model is printed bold.  
For data in the form of GxE means, the last random term in all models
above will become a residual term. If the GxE means are provided
together with weights, then a residual term will be added to the models
above.  
  
The function first fits a model where all model terms are included as
fixed terms. Based on the ANOVA table of this model, terms in the fixed
part of the model that are likely to give a problem when fitting the
mixed model are removed because of the reduced connectivity and number
of available observations to estimate that model term. Also a warning is
printed if the mean sum of squares for a model term points to a possible
zero variance component in the mixed model.  
  
Then a model is fitted where all model terms are included as random
terms. Based on the variance components in this model the percentage of
variance explained by each of the model components is determined. The
percentages of variance are printed in the model summary, together with
the variance components. The latter are presented on a standard
deviation scale.  
  
Finally a mixed model is fitted as specified in the overview above.
Based on this mixed model variance components can be computed using
[`vc`](vc.md), heritabilies can be computed using [`herit`](herit.md)
and predictions can be made using
[`predict.varComp`](predict.varComp.md). Predictions of genotypic
performance can be made at the level of individual trials, or for groups
of trials by using `predictLevel`.

## Usage

``` r
gxeVarComp(
  TD,
  trials = names(TD),
  trait,
  engine = c("lme4", "asreml"),
  locationYear = FALSE,
  nestingFactor = NULL,
  regionLocationYear = FALSE,
  useWt = FALSE,
  diagnostics = FALSE
)
```

## Arguments

- TD:

  An object of class [`TD`](https://rdrr.io/pkg/statgenSTA/man/TD.html).

- trials:

  A character string specifying the trials to be analyzed. If not
  supplied, all trials are used in the analysis.

- trait:

  A character string specifying the trait to be analyzed.

- engine:

  A character string specifying the engine used for modeling. Either
  "lme4" or "asreml".

- locationYear:

  Should a model be fitted assuming a factorial structure of locations x
  years?

- nestingFactor:

  A character string specifying a column in TD specifying the nesting
  structure of the trials.

- regionLocationYear:

  Should a model be fitted assuming locations within regions across
  years?

- useWt:

  Should the model be fitted using weights? Doing so requires a column
  wt in the data. If `useWt = FALSE`, the default, and the data contains
  no replicates, the last model term will be dropped and used as
  homogeneous residual.

- diagnostics:

  Should diagnostics on missing combinations of model variables be
  printed?

## Value

An object of class `varComp`, a list containing:

- fitMod:

  The fitted model.

- modDat:

  A data.frame containing the data used when fitting the model.

- nestingFactor:

  A name of the variable used as nesting variable in the model.

- useLocYear:

  A boolean specifying if a model containing location x year interaction
  was fitted.

- fullRandVC:

  A data.frame containing the variance components for the fully random
  model.

- aovFullMixedMod:

  A data.frame containing the ANOVA table for the fully fixed model.

- engine:

  The engine used for fitting the model.

- diagTabs:

  A list of data.frame, one for each random model term, containing the
  missing combinations in the data for that term.

## See also

Other Mixed model analysis: [`CRDR()`](CRDR.md),
[`correlations()`](correlations.md), [`diagnostics()`](diagnostics.md),
[`herit()`](herit.md), [`plot.varComp()`](plot.varComp.md),
[`predict.varComp()`](predict.varComp.md), [`vc()`](vc.md)

## Examples

``` r
## Fit a mixed model.
geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")

## Summarize results.
summary(geVarComp)
#> Fitted model formula final mixed model
#> 
#>  yld ~ trial + (1 | genotype) 
#> 
#> Sources of variation for fully random model:
#>  yld ~ (1 | trial) + (1 | genotype) 
#> 
#>           Component % Variance expl.
#> trial      86448.61          81.89 %
#> genotype    6670.93           6.32 %
#> residuals  12446.94          11.79 %
#> 
#> Analysis of Variance Table for fully fixed model:
#>  yld ~ trial + genotype 
#> 
#>             Df    Sum Sq  Mean Sq   F value    Pr(>F)    
#> trial        7 127771687 18253098 1466.4731 < 2.2e-16 ***
#> genotype   210  13821018    65814    5.2876 < 2.2e-16 ***
#> residuals 1470  18296997    12447                        
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Plot the standard deviations.
plot(geVarComp)


## Generate predictions
pred <- predict(geVarComp, predictLevel = "trial")
head(pred)
#>   genotype trial predictedValue
#> 1     G001 HN96b       529.1555
#> 2     G001 LN96a       227.7905
#> 3     G001 IS92a       684.1152
#> 4     G001 LN96b       133.7147
#> 5     G001 SS94a       457.4688
#> 6     G001 NS92a      1093.0337
```
