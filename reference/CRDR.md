# Calculate the correlated response to selection

Calculate the correlated response to selection (CRDR) based on the
fitted model. The CRDR is calculated as described by Atlin et al. E.g.
for a model with trials nested within scenarios, which has a random part
that looks like this: genotype + genotype:scenario +
genotype:scenario:trial the CRDR is calculated as:  
  
\$\$H1 = \sigma_G^2 / (\sigma_G^2 + \sigma_S^2 / s + \sigma\_{ST}^2 /
st + \sigma_E^2 / str)\$\$ \$\$H2 = (\sigma_G^2 + \sigma_S^2) /
(\sigma_G^2 + \sigma_S^2 + \sigma\_{ST}^2 / st + \sigma_E^2 / str)\$\$
\$\$CRDR = (\sigma_G^2 / (\sigma_G^2 + \sigma_S^2)) \* sqrt(H1 / H2)\$\$
In these formulas the \\\sigma\\ terms stand for the standard deviations
of the respective model terms, and the lower case letters for the number
of levels for the respective model terms. So \\\sigma_S\\ is the
standard deviation for the scenario term in the model and \\s\\ is the
number of scenarios. \\\sigma_E\\ corresponds to the residual standard
deviation and \\r\\ to the number of replicates.

## Usage

``` r
CRDR(varComp)
```

## Arguments

- varComp:

  An object of class varComp.

## References

Atlin, G. N., Baker, R. J., McRae, K. B., & Lu, X. (2000). Selection
response in subdivided target regions. Crop Science, 40(1), 7–13.
[doi:10.2135/cropsci2000.4017](https://doi.org/10.2135/cropsci2000.4017)

## See also

Other Mixed model analysis: [`correlations()`](correlations.md),
[`diagnostics()`](diagnostics.md), [`gxeVarComp()`](gxeVarComp.md),
[`herit()`](herit.md), [`plot.varComp()`](plot.varComp.md),
[`predict.varComp()`](predict.varComp.md), [`vc()`](vc.md)
