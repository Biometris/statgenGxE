#' Choose either \code{asreml} or \code{lme4} for mixed modelling
#'
#' \code{asreml} is used as the default package. If the package is not installed, use
#' \code{lme4} package instead.
#'
#' @return A character string: either \code{"lme4"} or \code{"asreml"} is returned
#' as \code{engine}.
#'
#' @export
asremlORlme4 <- function () {
  ## choose engine (default asreml if engine = 0)
  engine <- 0
  packages <- as.character(as.data.frame(installed.packages())$Package)
  if (!("asreml" %in% packages)) {
    warning("asreml is not found in R library.\n
         lme4 will be used as an engine for mixed modelling.\n")
    engine <- 1
  } else {
    if (requireNamespace("asreml", quietly = TRUE)) {
      ## check if the license is valid
      ## through a simulated example
      dat <- data.frame(y = rnorm(3), x = seq(1,3))
      license <- suppressPackageStartupMessages(capture.output(
        asreml::asreml(y ~ x, data = dat)$license))
      license <- license[length(license)]
      licenseDate <- unlist(regmatches(x = license,
                                       m = gregexpr(pattern =
                                                      "[[:digit:]]+-+[[:alpha:]]+-+[[:digit:]]{4}",
                                                    text = license)))
      if (length(licenseDate) == 0) {
        plic <- unlist(regmatches(license, gregexpr(pattern = "permanent",
                                                    text = license)))
        if(length(plic) == 0){
          ## use lme4 as the engine
          warning("The asreml license date is not found. \n
                    lme4 will be used as an engine for mixed modelling.\n")
          engine <- 1
        }
      } else {
        ## check if asreml license is valid
        licenseDate <- as.Date(licenseDate, "%d-%b-%Y")
        if (licenseDate < Sys.Date()) {
          warning("The asreml license had expired.\n
                    lme4 will be used as an engine for mixed modelling.\n")
          engine <- 1
        }
      }
    } else {
      engine <- 1
    }
  }
  ## check if the following pakcages is in the library:
  if (engine == 1) {
    ## if not in the library, install them from CRAN
    ## internet connection needed
    if (!"lme4" %in% packages) {
      install.packages("lme4", repos = "http://cran.r-project.org")
    }
    if (requireNamespace("lme4", quietly = TRUE)) {
      stop("No mixed modelling engine is loaded.\n")
    }
  }
  return(ifelse(engine == 0, "asreml", "lme4"))
}
