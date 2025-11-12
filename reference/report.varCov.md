# Report method for class varCov

A pdf report will be created containing a summary of an object of class
varCov. Simultaneously the same report will be created as a tex file.

## Usage

``` r
# S3 method for class 'varCov'
report(x, ..., outfile = NULL)
```

## Arguments

- x:

  An object of class varCov.

- ...:

  Not used.

- outfile:

  A character string, the name and location of the output .pdf and .tex
  file for the report. If `NULL`, a report with a default name will be
  created in the current working directory.

## Value

A pdf and tex report.

## See also

Other varCov: [`fitted.varCov()`](fitted.varCov.md),
[`gxeVarCov()`](gxeVarCov.md), [`plot.varCov()`](plot.varCov.md),
[`residuals.varCov()`](residuals.varCov.md)

## Examples

``` r
# \donttest{
## Select the best variance-covariance model using asreml for modeling.
if (requireNamespace("asreml", quietly = TRUE)) {
  geVarCov <- gxeVarCov(TD = TDMaize, trait = "yld", engine = "asreml")

  ## Create a pdf report summarizing the results.
  report(geVarCov, outfile = tempfile(fileext = ".pdf"))
  }
# }
```
