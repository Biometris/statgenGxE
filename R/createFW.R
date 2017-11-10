#' S3 class FW
#'
#' Function for creating objects of S3 class FW (Finlay Wilkinson).
#'
#' \code{\link{print}} and \code{\link{plot}} methods are available.
#'
#' @param estimates a data.frame containing the estimated values
#' @param anova a data.frame containing anova scores of the FW analysis
#' @param envEffs a data.frame containing the environmental effects
#' @param data the data.frame on which the analysis was performed
#' @param fittedGeno the fitted values for the genotypes
#' @param trait a character value indicating the analysed trait
#' @param nGeno a numerical value containing the number of genotypes in the analysis
#' @param nEnv a numerical value containing the number of environments in the analysis
#' @param tol a numerical value containing the tolerance used during the analysis
#' @param iter a numberical value containing the number of iterations for the
#' analysis to converge
#' @param x an \code{R} object
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.FW}}
#'
#' @name FW
NULL

#' @rdname FW
#' @export
createFW <- function(estimates,
                     anova,
                     envEffs,
                     trait,
                     nGeno,
                     nEnv,
                     data,
                     fittedGeno,
                     tol,
                     iter) {
  FW <- structure(list(estimates = estimates,
                       anova = anova,
                       envEffs = envEffs,
                       data = data,
                       fittedGeno = fittedGeno,
                       trait = trait,
                       nGeno = nGeno,
                       nEnv = nEnv,
                       tol = tol,
                       iter = iter),
                  class = "FW")
  return(FW)
}

#' @rdname FW
#' @export
is.FW <- function(x) {
  inherits(x, "FW")
}

#' @export
print.FW <- function(x, ...) {
  cat("Environmental effects",
      "\n===================\n")
  print(x$envEffs)
  cat("\nAnova",
      "\n=====\n")
  printCoefmat(x$anova)
  cat("\nEstimates",
      "\n=========\n")
  print(x$estimates, ..., row.names = FALSE)
}


