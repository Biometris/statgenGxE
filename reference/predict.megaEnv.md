# Compute BLUPS based on a set of mega environments

This function calculates Best Linear Unbiased Predictors (BLUPS) and
associated standard errors based on a set of mega environments.

## Usage

``` r
# S3 method for class 'megaEnv'
predict(
  object,
  ...,
  trials = names(object$TD),
  useYear = FALSE,
  engine = c("lme4", "asreml")
)
```

## Arguments

- object:

  An object of class megaEnv.

- ...:

  Further parameters passed to either `asreml` or `lmer`.

- trials:

  A character string specifying the trials to be analyzed. If not
  supplied, all trials are used in the analysis.

- useYear:

  Should year be used for modeling (as years within trials). If `TRUE`,
  `TD` should contain a column "year".

- engine:

  A character string specifying the engine used for modeling.

## Value

A list consisting of two data.frames, `predictedValue` containing BLUPs
per genotype per mega environment and `standardError` containing
standard errors for those BLUPs.

## See also

Other mega environments: [`gxeMegaEnv()`](gxeMegaEnv.md),
[`plot.megaEnv()`](plot.megaEnv.md)

## Examples

``` r
## Compute mega environments for TDMaize.
geMegaEnv <- gxeMegaEnv(TD = TDMaize, trait = "yld")

## Compute BLUPS and standard errors for those mega environments.
megaEnvPred <- predict(geMegaEnv)
#> Warning: One should be cautious with the interpretation of predictions for mega environments that are based on less than 10 trials.
#> boundary (singular) fit: see help('isSingular')
head(megaEnvPred$predictedValue)
#>      megaEnv_1 megaEnv_2 megaEnv_3
#> G001  489.0515  483.5317  497.8785
#> G002  472.3096  491.5258  498.1123
#> G003  468.2102  457.0879  462.5265
#> G004  363.6653  361.1352  321.1292
#> G005  462.8374  465.0800  467.9159
#> G006  435.4026  419.0985  410.4387
head(megaEnvPred$standardError)
#>      megaEnv_1 megaEnv_2 megaEnv_3
#> G001  32.09261  38.36256  45.40008
#> G002  32.09261  38.36256  45.40008
#> G003  32.09261  38.36256  45.40008
#> G004  32.09261  38.36256  45.40008
#> G005  32.09261  38.36256  45.40008
#> G006  32.09261  38.36256  45.40008
```
