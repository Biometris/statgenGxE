#' S3 class varComp
#'
#' Function for creating objects of S3 class varComp.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param SSA The best fitted model.
#' @param choice A character string indicating the best fitted model.
#' @param summary A data.frame with a summary of the fitted models.
#' @param vcov The covariance matrix of the best fitted model.
#' @param criterion A character string indicating the criterion used for
#' determinening the best model, either "AIC" or "BIC".
#' @param engine A character string containing the engine used for
#' the analysis.
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.varComp}}, \code{\link{report.varComp}}
#'
#' @name varComp
NULL

#' @rdname varComp
#' @export
createVarComp <- function(SSA,
                          choice,
                          summary,
                          vcov,
                          criterion,
                          engine) {
  varComp <- structure(list(SSA = SSA,
                            choice = choice,
                            summary = summary,
                            vcov = vcov,
                            criterion = criterion,
                            engine = engine),
                       class = "varComp")
  attr(varComp, which = "timestamp") <- Sys.time()
  return(varComp)
}

#' @export
print.varComp <- function(x, ...) {
  x$summary
}

#' @export
summary.varComp <- function(object, ...) {
  print(object, ...)
}

#' Plot Function for Class varComp
#'
#' Function for plotting a heatmap of the correlation matrix for objects of
#' class varComp.
#'
#' @param x An object of class varComp
#' @param ... not unused
#'
#' @examples
#' \dontrun{
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld", engine = "asreml")
#' plot(geVarComp)
#' }
#'
#' @import stats
#' @export
plot.varComp <- function(x, ...) {
  dotArgs <- list(...)
  ## Set arguments for plot
  plotArgs <- list(corMat = cov2cor(x$vcov),
                   main = paste("Heatmap for correlations for model:", x$choice))
  ## Add and overwrite args with custom args from ...
  fixedArgs <- c("corMat")
  plotArgs <- modifyList(plotArgs, dotArgs[!names(dotArgs) %in% fixedArgs])
  do.call(plotCorMat, args = plotArgs)
}

#' Report method for class varComp
#'
#' A pdf report will be created containing a summary of selection of the
#' best variance model.
#' Simultaneously the same report will be created as a tex file.
#'
#' @param x An object of class varComp.
#' @param ... Further arguments passed on from other functions - not used yet.
#' @param outfile A character string, the name and location of the output .pdf and .tex
#' file for the report. If \code{NULL} a report will be created in the current working
#' directory.
#'
#' @export
report.varComp <- function(x,
                           ...,
                           outfile = NULL) {
  createReport(x = x, reportName = "varCompReport.Rnw",
               outfile = outfile, ...)
}



