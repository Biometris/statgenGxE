# statgenGxE 1.0.4.1

* The predict function for gxeVarComp output is extended so all variables in the fitted model can now be used for making predictions.
* The plot functions for AMMI and GGE analysis now have an argument `rotatePC` allowing the specification of a trial that is aligned with the positive x-axis in the plot.
* The gxeVarComp function now has an argument `models` allowing a subset of the available models to be fitted.
* A small bug in gxeMegaEnv that sometimes caused NA for predicted values is fixed.
* Some minor changes in order and capitalization of outputs.

# statgenGxE 1.0.4

* When in a Finlay-Wilkinson analysis one or more genotypes are observed in only one trial a user friendly warning message is now shown.
* A wrong classification of one of the environments in the data used in the vignette is corrected.
* The degrees of freedom in the anova table for the GGE analysis are corrected. The value for the PC terms is lowered by one.

# statgenGxE 1.0.3

* No user visible changes

# statgenGxE 1.0.2

* Initial CRAN version
