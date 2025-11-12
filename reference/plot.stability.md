# Plot function for class stability

Function for creating scatter plots of the square roots of the computed
stability measures against the means.

## Usage

``` r
# S3 method for class 'stability'
plot(
  x,
  ...,
  colorGenoBy = NULL,
  colGeno = NULL,
  title = paste("Stability coefficients for", x$trait),
  output = TRUE
)
```

## Arguments

- x:

  An object of class stability.

- ...:

  Not used.

- colorGenoBy:

  A character string indicating a column in the `TD` used as input for
  the stability analysis by which the genotypes should be colored. If
  `NULL` all genotypes will be colored black.

- colGeno:

  A character vector with plot colors for the genotypes. A single color
  when `colorGenoBy = NULL`, a vector of colors otherwise.

- title:

  A character string used a title for the plot.

- output:

  Should the plot be output to the current device? If `FALSE` only a
  list of ggplot objects is invisibly returned.

## Value

A list of ggplot object is invisibly returned.

## See also

Other stability: [`gxeStability()`](gxeStability.md),
[`report.stability()`](report.stability.md)

## Examples

``` r
## Compute three stability measures for TDMaize.
geStab <- gxeStability(TD = TDMaize, trait = "yld")

## Create scatter plots of the computed stability measures against the means.
plot(geStab)

```
