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

---

## Comments from submission v.1.0.1

* Please do not start the description with "This package", package name, title or similar.

  - Fixed.

* If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form authors (year) <doi:...> authors (year) <arXiv:...> authors (year, ISBN:...) or if those are not available: <https:...> with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking.
(If you want to add a title as well please put it in quotes: "Title")

  - Reference to article added.

* Please always make sure to reset to user's options(), working directory or par() after you changed it in examples and vignettes and demos.
e.g.:

old <- options()
options(digits = 3)
...
options(old)
e.g.: Vignette

  - Options are now reset at end of vignette.

* Please do not modify the global environment (e.g. by using <<-) in your functions. This is not allowed by the CRAN policies.

  - I did use <<-, but as far as I could see I was not modifying the global environment.

  

