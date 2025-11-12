# Selects the best variance-covariance model for a set of trials

This function selects the best covariance structure for genetic
correlations between trials. It fits a range of variance-covariance
models (identity, compound symmetry (cs), diagonal, simple correlation
with heterogeneous variance (outside), heterogeneous compound symmetry
(hcs), first order factor analytic (fa), second order factor analytic
(fa2) and unstructured), and selects the best one using a
goodness-of-fit criterion. See details for the exact models fitted.

## Usage

``` r
gxeVarCov(
  TD,
  trials = names(TD),
  trait,
  models = c("identity", "cs", "diagonal", "hcs", "outside", "fa", "fa2", "unstructured"),
  engine = c("lme4", "asreml"),
  criterion = c("BIC", "AIC"),
  ...
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

- models:

  A character vector specifying the models to be fitted.

- engine:

  A character string specifying the engine used for modeling. Either
  "lme4" or "asreml".

- criterion:

  A string specifying a goodness-of-fit criterion. Either "AIC" or
  "BIC".

- ...:

  Further arguments to be passed to the modeling engine.

## Value

An object of class [`varCov`](varCov.md), a list object containing:

- STA:

  An object of class STA containing the best fitted model.

- choice:

  A character string indicating the best fitted model.

- summary:

  A data.frame with a summary of the fitted models.

- vcov:

  The covariance matrix of the best fitted model.

- criterion:

  A character string indicating the goodness-of-fit criterion used for
  determining the best model, either "AIC" or "BIC".

- engine:

  A character string containing the engine used for the analysis.

- trait:

  A character string containing the trait analyzed.

- dat:

  A data.frame with the full data set used for the analysis.

## Details

The models fitted are of the form \\y\_{ij} = \mu_j + \epsilon\_{ij}\\,
where \\y\_{ij}\\ is the phenotypic value of genotype \\i\\ in
environment \\j\\, \\\mu_j\\ is the environmental mean, and
\\\epsilon\_{ij}\\ represents mainly genetic variation, although some
non-genetic variation may be included as well. The random term
\\\epsilon\_{ij}\\ is modeled in eight ways as described in the table
below.

|              |                                   |                                                      |                                                           |                      |
|--------------|-----------------------------------|------------------------------------------------------|-----------------------------------------------------------|----------------------|
| Model        | Description                       | var(\\g\_{ij}\\)                                     | cov(\\g\_{ij}\\;\\g\_{ik}\\)                              | Number of parameters |
| identity     | identity                          | \\\sigma_G^2\\                                       | 0                                                         | 1                    |
| cs           | compound symmetry                 | \\\sigma_G^2+\sigma\_{GE}^2\\                        | \\\sigma\_{GE}^2\\                                        | 2                    |
| diagonal     | diagonal matrix (heteroscedastic) | \\\sigma\_{GE_j}^2\\                                 | 0                                                         | \\J\\                |
| hcs          | heterogeneous compound symmetry   | \\\sigma_G^2+\sigma\_{GE_j}^2\\                      | \\\sigma_G^2\\                                            | \\J+1\\              |
| outside      | heterogeneity outside             | \\\sigma\_{G_j}^2\\                                  | \\\theta\\                                                | \\J+1\\              |
| fa           | first order factor analytic       | \\\lambda\_{1j}^2+\sigma\_{GE_j}^2\\                 | \\\lambda\_{1j}\lambda\_{1k}\\                            | \\2J\\               |
| fa2          | second order factor analytic      | \\\lambda\_{1j}^2+\lambda\_{2j}^2+\sigma\_{GE_j}^2\\ | \\\lambda\_{1j}\lambda\_{1k}+\lambda\_{2j}\lambda\_{2k}\\ | \\3J-1\\             |
| unstructured | unstructured                      | \\\sigma\_{G_j}^2\\                                  | \\\sigma\_{G\_{j,k}}^2\\                                  | \\J(J+1)/2\\         |

In this table \\J\\ is the number of environments, \\\sigma_G^2\\ the
variance component for the genotype main effects, \\\sigma\_{GE}^2\\ the
variance component for GxE interactions. \\\sigma\_{G_j}^2\\ and
\\\sigma\_{GE_j}^2\\ are the environment specific variance components
for the genotype main effects and GxE interaction in environment \\j\\.
\\\sigma\_{G\_{j,k}}^2\\ is the genetic covariance between environments
\\j\\ and \\k\\. \\\theta\\ is the common correlation between
environments and \\\lambda\_{1j}\\ and \\\lambda\_{2j}\\ are environment
specific multiplicative parameters.

## Note

If `engine = "lme4"`, only the compound symmetry model can be fitted.

## See also

Other varCov: [`fitted.varCov()`](fitted.varCov.md),
[`plot.varCov()`](plot.varCov.md),
[`report.varCov()`](report.varCov.md),
[`residuals.varCov()`](residuals.varCov.md)

## Examples

``` r
## Select the best variance-covariance model using lme4 for modeling.
geVarCov <- gxeVarCov(TD = TDMaize, trait = "yld")

## Summarize results.
summary(geVarCov)
#> Best model: cs, based on BIC.
#> 
#>         AIC      BIC Deviance NParameters
#> cs 21005.28 21016.14 21001.28           2

# \donttest{
## Create a pdf report summarizing the results.
report(geVarCov, outfile = tempfile(fileext = ".pdf"))
#> Error in report.varCov(geVarCov, outfile = tempfile(fileext = ".pdf")): An installation of LaTeX is required to create a pdf report.
# }

# \donttest{
if (requireNamespace("asreml", quietly = TRUE)) {
  ## Select the best variance-covariance model using asreml for modeling.
  ## Use BIC as a goodness-of-fit criterion.
  geVarCov2 <- gxeVarCov(TD = TDMaize, trait = "yld", engine = "asreml",
                         criterion = "BIC")

  summary(geVarCov2)

  ## Plot a heatmap of the correlation matrix for the best model.
  plot(geVarCov2)
  }
# }
```
