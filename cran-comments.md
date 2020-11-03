## Release to correct test error for noLD by decreasing tolerance for one test.

----

## Test environments

* local Windows 10 install, R 4.0.3
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



  

