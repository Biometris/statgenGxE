# Estimate missing values in multivariate data

This function estimates missing values for units in a multivariate
dataset, using an iterative regression technique.

## Usage

``` r
multMissing(Y, maxIter = 10, naStrings = NULL)
```

## Arguments

- Y:

  A matrix, data.frame or vector of multivariate data.

- maxIter:

  An integer specifying the maximum number of iterations.

- naStrings:

  A character vector of strings which are to be interpreted as `NA`
  values.

## Value

An object of the same class as the input `Y` with the missing values
replaced by their estimates.

## Details

Initial estimates of the missing values in each variate are formed from
the variate means using the values for units that have no missing values
for any variate. Estimates of the missing values for each variate are
then recalculated as the fitted values from the multiple regression of
that variate on all the other variates. When all the missing values have
been estimated the variate means are recalculated. If any of the means
differs from the previous mean by more than a tolerance (the initial
standard error divided by 1000) the process is repeated, subject to a
maximum number of repetitions defined by `maxIter` option. The default
maximum number of iterations (10) is usually sufficient when there are
few missing values, say two or three. If there are many more, 20 or so,
it may be necessary to increase the maximum number of iterations to
around 30. The method is similar to that of Orchard & Woodbury (1972),
but does not adjust for bias in the variance-covariance matrix as
suggested by Beale & Little (1975).

## References

Beale, E.M.L. & Little, R.J.A. (1975). Missing values in multivariate
analysis. Journal of the Royal Statistical Society, Series B, 37,
129-145.

Orchard, T. & Woodbury, M.A. (1972). A missing information principle:
theory and applications. In: Proceedings of the 6th Berkeley Symposium
in Mathematical Statistics and Probability, Vol I, 697-715.

## Examples

``` r
M <- matrix(c("1", "2", "3", NA, "b", "5", "6",
              "6", "5", "b", NA, "3", "2", "1"), nrow = 7, ncol = 2)

## Estimate missing values treating "b" as NA.
multMissing(M, naStrings = "b")
#>          [,1]     [,2]
#> [1,] 1.000000 6.000000
#> [2,] 2.000000 5.000000
#> [3,] 3.000000 3.999972
#> [4,] 3.499975 3.499995
#> [5,] 3.999851 3.000000
#> [6,] 5.000000 2.000000
#> [7,] 6.000000 1.000000
```
