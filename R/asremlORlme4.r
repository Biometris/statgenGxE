#' Choose either \code{asreml} or \code{lme4} for mixed modelling
#'
#' \code{asreml} is used as the default package. If the package is not installed, use \code{lme4} package instead.
#' @return A character string: either \code{"lme4"} or \code{"asreml"} is returned as \code{engine}.
#'
#' @export
asremlORlme4 = function (){

    # choose engine (default asreml if engine = 0)
    engine <- 0
    PACKages <- as.character(as.data.frame(installed.packages())$Package)
    test0 <-("asreml" %in% PACKages)
    if (!test0){
        warning("asreml is not found in R library.\n
         lme4 will be used as an engine for mixed modelling.\n")
        engine <- 1
    }else{
        flag=try(require(asreml,quietly=T),silent=T)
        if(!inherits(flag, "try-error")){
            flag1 = require(asreml,quietly=T)
            # check if the license is valid
            # through a simulated example
            dat <- data.frame(y=rnorm(20),x=seq(1,20))
            license <- capture.output(asreml::asreml(y ~ x, data=dat)$license)
            license.date <- unlist(regmatches(license,
            gregexpr("[[:digit:]]+-+[[:alpha:]]+-+[[:digit:]]{4}",license)))
            if (length(license.date)==0){
                plic <- unlist(regmatches(license, gregexpr("permanent",license)))
                if(length(plic)==0){
                    # use lme4 as the engine
                    warning("The asreml license date is not found. \n
                    lme4 will be used as an engine for mixed modelling.\n")
                    engine <- 1
                }
            }else{
                # check if asreml license is valid
                license.date <- as.Date(license.date,"%d-%b-%Y")
                if (as.numeric(license.date-Sys.Date())<0){
                    warning("The asreml license expires. \n
                    lme4 will be used as an engine for mixed modelling.\n")
                    engine <- 1
                }
            }
        }else{
            engine <- 1
        }
    }

    # check if the following pakcages is in the library:
    if (engine ==1){
        packs.req <- "lme4"
        test <- (packs.req %in% PACKages)
        # if not in the library, install them from CRAN
        # internet connection needed
        if (!test)
            install.packages(packs.req, repos = "http://cran.r-project.org")
        flag1=suppressPackageStartupMessages(require(lme4))
    }
    if (flag1){
        if (engine ==1) engine="lme4"
        if (engine ==0) engine="asreml"
    }else{
        stop("No mixed modelling engine is loaded.\n")
    }
    engine
}
