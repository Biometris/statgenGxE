# Changelog

## statgenGxE 1.0.11

- The method for making predictions from the output of `gxeVarComp` has
  been improved. When `lme4` is used for fitting variance component
  models the predictions for unbalanced data are more accurate now.

## statgenGxE 1.0.10

CRAN release: 2025-06-24

- A bug in `herit` has been fixed. For models were the last model term
  was removed due to confounding with the residual the heritabilitiy
  wasn’t calculated correctly.
- Tests are modified for compatibility with the upcoming version of
  ggplot2.

## statgenGxE 1.0.9

CRAN release: 2024-09-18

- A bug in `gxeFW` that in some cases caused negative sums of squares in
  the resulting anova table has been fixed.
- A bug is fixed that caused the means in the stability plots to be
  plotted as a straight line.

## statgenGxE 1.0.8

CRAN release: 2024-05-06

- Patch release for R 4.4.0. No user visual changes.

## statgenGxE 1.0.7

CRAN release: 2024-03-19

- No user visual changes

## statgenGxE 1.0.6

CRAN release: 2023-12-12

- Functions no longer rely on soft-deprecated ggplot2 functions.
- A small bug is fixed that made plotting of `gxeVarComp` output
  impossible when using asreml for fitting the models.

## statgenGxE 1.0.5

CRAN release: 2022-08-11

- The predict function for `gxeVarComp` output is extended so all
  variables in the fitted model can now be used for making predictions.
- The plot functions for AMMI and GGE analysis now have an argument
  `rotatePC` allowing the specification of a trial that is aligned with
  the positive x-axis in the plot.
- The `gxeVarCov` function now has an argument `models` allowing a
  subset of the available models to be fitted.
- A small bug in `gxeMegaEnv` that sometimes caused NA for predicted
  values is fixed.
- Some minor changes in order and capitalization of outputs.

## statgenGxE 1.0.4

CRAN release: 2021-01-07

- When in a Finlay-Wilkinson analysis one or more genotypes are observed
  in only one trial a user friendly warning message is now shown.
- A wrong classification of one of the environments in the data used in
  the vignette is corrected.
- The degrees of freedom in the anova table for the GGE analysis are
  corrected. The value for the PC terms is lowered by one.

## statgenGxE 1.0.3

CRAN release: 2020-11-09

- No user visible changes

## statgenGxE 1.0.2

CRAN release: 2020-10-01

- Initial CRAN version
