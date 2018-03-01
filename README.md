# RAP

R Analytical Pipeline for 
* Single trial analysis
* GxE analysis
* QTL mapping

# Installation

For direct installation from gitlab use the following code:

``` r
## Replace the location for public and private key with your own.
creds <- git2r::cred_ssh_key(publickey = "M:\\.ssh\\id_rsa.pub",
                             privatekey = "M:\\.ssh\\id_rsa")
devtools::install_git(url = "git@git.wur.nl:rossu027/RAP.git",
                      credentials = creds)

```

# Implemented functionality

So far the following functionality has been implemented:
* Single trial analysis
 * Single trial analysis using SpATS, asreml or lme4
 * Computation of statistics e.g. BLUEs, BLUPs, SE of BLUEs, SE of BLUPs 
and heritability based on the fitted models.
* GxE analysis
 * Selection of the best variance-covariance model
 * AMMI analysis
 * Finlay-Wilkinson anlysis
 * Computation of mega-environments based on AMMI results
 * Stability measures
* QTL Mapping
 * Quality Control and cleaning
 * QTL Detection and selection of peaks
 * Multi QTL Mapping 
* Plotting and reporting (.pdf and .tex) functions for most analyses

