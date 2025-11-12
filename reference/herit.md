# Calculate heritability

Calculate the heritability based on the fitted model. For balanced data,
the heritability is calculated as described by Atlin et al. E.g. for a
model with trials nested within locations, which has a random part that
looks like this: genotype + genotype:location + genotype:location:trial
the heritability is computed as  
  
\$\$\sigma_G^2 / (\sigma_G^2 + \sigma_L^2 / l + \sigma\_{LT}^2 / lt +
\sigma_E^2 / ltr)\$\$ In this formula the \\\sigma\\ terms stand for the
standard deviations of the respective model terms, and the lower case
letters for the number of levels for the respective model terms. So
\\\sigma_L\\ is the standard deviation for the location term in the
model and \\l\\ is the number of locations. \\\sigma_E\\ corresponds to
the residual standard deviation and \\r\\ to the number of replicates.  
  
When the data is unbalanced a more general form of this formula is used
as described in Holland et al. Here the numerator \\l\\ is replaced by
the harmonic means of the number of locations across genotypes. The
other numerators are replaced correspondingly. For balanced data this
more general form gives identical results as the form described by Atlin
et al.

## Usage

``` r
herit(varComp)
```

## Arguments

- varComp:

  An object of class varComp.

## References

Atlin, G. N., Baker, R. J., McRae, K. B., & Lu, X. (2000). Selection
response in subdivided target regions. Crop Science, 40(1), 7–13.
[doi:10.2135/cropsci2000.4017](https://doi.org/10.2135/cropsci2000.4017)

Holland, J.B., W.E. Nyquist, and C.T. Cervantes-Martínez. (2003).
Estimating and interpreting heritability for plant breeding: An update.
Plant Breed. Rev. 2003:9–112.
[doi:10.1002/9780470650202.ch2](https://doi.org/10.1002/9780470650202.ch2)

## See also

Other Mixed model analysis: [`CRDR()`](CRDR.md),
[`correlations()`](correlations.md), [`diagnostics()`](diagnostics.md),
[`gxeVarComp()`](gxeVarComp.md), [`plot.varComp()`](plot.varComp.md),
[`predict.varComp()`](predict.varComp.md), [`vc()`](vc.md)

## Examples

``` r
## Fit a mixed model.
geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")

## Compute heritability.
herit(geVarComp)
#> [1] 0.8108781
```
