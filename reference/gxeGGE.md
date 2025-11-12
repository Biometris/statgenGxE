# GGE analysis

The Genotype plus Genotype by Environment interaction (GGE) model fits a
model with trial as main fixed effect. Then a principal component
analysis is done on the residuals. This results in an interaction
characterized by Interaction Principal Components (IPCA) enabling
simultaneous plotting of genotypes and trials.  
  
The parameter `nPC` is used to indicate the number of principal
components that is used in the principal component analysis (PCA). By
setting this parameter to `NULL` the algorithm determines the best
number of principal components (see Details).  
  
By specifying the parameter `byYear = TRUE`, a separate analysis will be
done for every year in the data. Combining the option with `nPC = NULL`
may result in different numbers of principal components per year. The
GGE estimates will still be returned as a single data.frame, but the
other results will be either lists or arrays.

## Usage

``` r
gxeGGE(
  TD,
  trials = names(TD),
  trait,
  nPC = 2,
  byYear = FALSE,
  center = TRUE,
  excludeGeno = NULL,
  useWt = FALSE
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

- nPC:

  An integer specifying the number of principal components used as
  multiplicative term of genotype-by-trial interaction. If `NULL`, the
  number of principal components is determined by the algorithm using
  forward selection. See details.

- byYear:

  Should the analysis be done by year? If `TRUE` the data is split by
  the variable year, analysis is performed and the results are merged
  together and returned.

- center:

  Should the variables be shifted to be zero centered?

- excludeGeno:

  An optional character vector with names of genotypes to be excluded
  from the analysis. If `NULL`, all genotypes are used.

- useWt:

  Should weighting be used when modeling? Requires a column `wt` in
  `TD`.

## Details

First a linear model \\trait = trial + \epsilon\\ is fitted with trial a
fixed component in the model.  
The residuals from the fitted model are then used in a PCA. If `nPC` is
not `NULL` a single PCA is done using
[`prcomp`](https://rdrr.io/r/stats/prcomp.html) with maximum rank
`nPC`.  
In case `nPC = NULL`, the PCA is first done with one PC. Then using
forward selection one by one the number of PCs is increased as long as
the added component is significant in the analysis.  
GGE estimates are then computed using the results of the PCA.  
