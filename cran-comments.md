## Initial release.

----

## Test environments

* local Windows 10 install, R 4.0.2
* Ubuntu (on github actions, devel and release)
* macOS (on github actions, release)
* R-hub (devel and release)

----

## R CMD check results

There were no ERRORs or WARNINGs.

There were 2 NOTES:

  * checking CRAN incoming feasibility ... NOTE    
  
    New submission

    Suggests or Enhances not in mainstream repositories: asreml

    - asreml is a commercial R package that is used as one of two alternatives in some functions.

  * checking package dependencies ... NOTE  
    Package suggested but not available for checking: 'asreml'
    
    - asreml is a commercial R package that is used as one of two alternatives in some functions.

## Notes from submitting version 1.0.0

* checking CRAN incoming feasibility ... NOTE  

  Possibly mis-spelled words in DESCRIPTION:
    Biometris (53:53)
    GxE (3:33, 50:39)
    
  - These are spelled correctly.
  
  Found the following (possibly) invalid URLs:
    URL: https://cordis.europa.eu/project/id/244374
      From: inst/doc/statgenGxE.html
      Status: Error
      Message: libcurl error code 35:
      	 schannel: next InitializeSecurityContext failed: SEC_E_ILLEGAL_MESSAGE    (0x80090326) - This error usually occurs when a fatal SSL/TLS alert is received (e.g. handshake failed).
      	 
  - The link works fine from a browser      	 

** running examples for arch 'x64' ... [49s] NOTE
Examples with CPU (user + system) or elapsed time > 10s
         user system elapsed
plot.FW 10.57   0.08   10.69

  - Suppressed testing of some of the examples

