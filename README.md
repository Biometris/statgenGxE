# RAP

R Analytical Pipeline for 
* Single trial analysis
* GxE analysis
* QTL mapping

# Installation

For direct installation from gitlab use the following code:

``` r
## Replace the location for public and privatekey with your own.
creds <- git2r::cred_ssh_key(publickey = "M:\\.ssh\\id_rsa.pub",
                             privatekey = "M:\\.ssh\\id_rsa")
devtools::install_git(url = "git@git.wur.nl:rossu027/RAP.git",
                      credentials = creds)

```


