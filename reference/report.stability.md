# Report method for class stability

A pdf report will be created containing a summary of an object of class
stability. Simultaneously the same report will be created as a tex file.

## Usage

``` r
# S3 method for class 'stability'
report(x, ..., outfile = NULL)
```

## Arguments

- x:

  An object of class stability.

- ...:

  Not used.

- outfile:

  A character string, the name and location of the output .pdf and .tex
  file for the report. If `NULL`, a report with a default name will be
  created in the current working directory.

## Value

A pdf and tex report.

## See also

Other stability: [`gxeStability()`](gxeStability.md),
[`plot.stability()`](plot.stability.md)

## Examples

``` r
## Compute three stability measures for TDMaize.
geStab <- gxeStability(TD = TDMaize, trait = "yld")
# \donttest{
## Create a .pdf report summarizing the stability measures.
report(geStab, outfile = tempfile(fileext = ".pdf"))
#> Error in report.stability(geStab, outfile = tempfile(fileext = ".pdf")): An installation of LaTeX is required to create a pdf report.
# }
```
