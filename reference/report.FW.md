# Report method for class FW

A pdf report will be created containing a summary of an FW object.
Simultaneously the same report will be created as a tex file.

## Usage

``` r
# S3 method for class 'FW'
report(x, sortBy = c("sens", "genMean", "mse"), ..., outfile = NULL)
```

## Arguments

- x:

  An object of class FW.

- sortBy:

  A character string indicating by which variable the estimates should
  be sorted. Either `sens`(itivity), `genMean` (genotypic Mean) or `mse`
  (mean squared error).

- ...:

  Not used.

- outfile:

  A character string, the name and location of the output .pdf and .tex
  file for the report. If `NULL`, a report with a default name will be
  created in the current working directory.

## Value

A pdf and tex report.

## See also

Other Finlay-Wilkinson: [`fitted.FW()`](fitted.FW.md),
[`gxeFw()`](gxeFw.md), [`plot.FW()`](plot.FW.md),
[`residuals.FW()`](residuals.FW.md)

## Examples

``` r
## Run Finlay-Wilkinson analysis on TDMaize.
geFW <- gxeFw(TDMaize, trait = "yld")
#> Warning: ANOVA F-tests on an essentially perfect fit are unreliable

# \donttest{
## Create a report summarizing the results.
report(geFW, outfile = tempfile(fileext = ".pdf"))
#> Error in report.FW(geFW, outfile = tempfile(fileext = ".pdf")): An installation of LaTeX is required to create a pdf report.
# }
```
