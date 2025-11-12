# Compute different types of correlations.

Compute three types of correlations for models fitted with a nesting
factor.

- correlation between scenarios or environment types: \$\$\sigma_G^2 /
  (\sigma_G^2 + \sigma\_{GS}^2)\$\$

- correlation between trials within scenarios or environment types:
  \$\$(\sigma_G^2 + \sigma\_{GS}^2) / (\sigma_G^2 + \sigma\_{GS}^2 +
  \sigma_E^2)\$\$

- correlation between trials that belong to different
  scenarios/environment types: \$\$\sigma_G^2 / (\sigma_G^2 +
  \sigma\_{GS}^2 + \sigma_E^2)\$\$

In these formulas the \\\sigma\\ terms stand for the standard deviations
of the respective model terms. So \\\sigma_S\\ is the standard deviation
for the scenario term in the model, \\\sigma\_{GS}\\ for the standard
deviation of the genotype by scenario term and \\\sigma_E\\ corresponds
to the residual standard deviation.

## Usage

``` r
correlations(varComp)
```

## Arguments

- varComp:

  An object of class varComp.

## Value

A list with three correlations.

## See also

Other Mixed model analysis: [`CRDR()`](CRDR.md),
[`diagnostics()`](diagnostics.md), [`gxeVarComp()`](gxeVarComp.md),
[`herit()`](herit.md), [`plot.varComp()`](plot.varComp.md),
[`predict.varComp()`](predict.varComp.md), [`vc()`](vc.md)
