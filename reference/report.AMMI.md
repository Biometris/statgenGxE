# Report method for class AMMI

A pdf report will be created containing a summary of an AMMI object.
Simultaneously the same report will be created as a tex file.

## Usage

``` r
# S3 method for class 'AMMI'
report(x, ..., outfile = NULL)
```

## Arguments

- x:

  An object of class AMMI.

- ...:

  Not used.

- outfile:

  A character string, the name and location of the output .pdf and .tex
  file for the report. If `NULL`, a report with a default name will be
  created in the current working directory.

## Value

A pdf and tex report.

## See also

Other AMMI: [`fitted.AMMI()`](fitted.AMMI.md),
[`gxeAmmi()`](gxeAmmi.md), [`plot.AMMI()`](plot.AMMI.md),
[`residuals.AMMI()`](residuals.AMMI.md)

## Examples

``` r
## Run AMMI analysis on TDMaize.
geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
# \donttest{
## Create a pdf report summarizing the results.
report(geAmmi, outfile = tempfile(fileext = ".pdf"))
#> Error in report.AMMI(geAmmi, outfile = tempfile(fileext = ".pdf")): An installation of LaTeX is required to create a pdf report.
# }
```
