[![pipeline status](https://git.wur.nl/statistical-genetic-pipeline/statgenGxE/badges/master/pipeline.svg)](https://git.wur.nl/statistical-genetic-pipeline/statgenGxE/commits/master)
[![coverage report](https://git.wur.nl/statistical-genetic-pipeline/statgenGxE/badges/master/coverage.svg)](https://git.wur.nl/statistical-genetic-pipeline/statgenGxE/commits/master)

# statgenGxE

R package for Genotype by Environment (GxE) analysis

# Installation

For direct installation from gitlab use the following code:

``` r
## Replace the location for public and private key with your own.
creds <- git2r::cred_ssh_key(publickey = "M:\\.ssh\\id_rsa.pub",
                             privatekey = "M:\\.ssh\\id_rsa")
devtools::install_git(url = "git@git.wur.nl:statistical-genetic-pipeline/statgenGxE.git",
                      credentials = creds)

```

# Implemented functionality

The following functionality has been implemented:

* GxE analysis
    * Selection of the best variance-covariance model
    * AMMI analysis
    * Finlay-Wilkinson anlysis
    * Computation of mega environments based on AMMI results
    * Stability measures

* Plotting and reporting (.pdf and .tex) functions for most analyses.
