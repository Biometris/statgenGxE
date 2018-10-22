#' S3 class varComp
#'
#' Function for creating objects of S3 class varComp.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param SSA An object of class SSA, the best fitted model.
#' @param choice A character string indicating the best fitted model.
#' @param summary A data.frame with a summary of the fitted models.
#' @param vcov The covariance matrix of the best fitted model.
#' @param criterion A character string indicating the goodness-of-fit criterion
#' used for determinening the best model.
#' @param engine A character string containing the engine used for the analysis.
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.varComp}}, \code{\link{report.varComp}}
#'
#' @name varComp
NULL

#' @rdname varComp
#' @keywords internal
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
  cat(paste0("Best model: ", x$choice, ", based on ", x$criterion, ".\n\n"))

  x$summary
}

#' @export
summary.varComp <- function(object, ...) {
  print(object, ...)
}

#' Plot function for class varComp
#'
#' Function for plotting a heatmap of the correlation matrix for objects of
#' class varComp.
#'
#' @param x An object of class varComp
#' @param ... Not used
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @examples
#' \dontrun{
#' ## Select the best variance-covariance model using asreml for modeling.
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld", engine = "asreml")
#' ## Create a heatmap of the correlation matrix for the best model.
#' plot(geVarComp)
#' }
#'
#' @import stats
#' @importFrom utils modifyList
#' @export
plot.varComp <- function(x,
                         ...,
                         output = TRUE) {
  dotArgs <- list(...)
  ## Set arguments for plot
  plotArgs <- list(varMat = x$vcov,
                   main = paste("Heatmap for model:", x$choice))
  ## Add and overwrite args with custom args from ...
  fixedArgs <- c("varMat")
  plotArgs <- modifyList(plotArgs, dotArgs[!names(dotArgs) %in% fixedArgs])
  p <- do.call(plotCorMat, args = plotArgs)
  if (output) {
    plot(p)
  }
  invisible(p)
}

#' Report method for class varComp
#'
#' A pdf report will be created containing a summary of an object of class
#' varComp. Simultaneously the same report will be created as a tex
#' file.
#'
#' @inheritParams report.AMMI
#'
#' @param x An object of class varComp.
#'
#' @return A pdf and tex report.
#'
#' @examples
#' \dontrun{
#' ## Select the best variance-covariance model using asreml for modeling.
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld", engine = "asreml")
#' ## Create a pdf report summarizing the results.
#' report(geVarComp, outfile = "./testReports/reportVarComp.pdf")
#' }
#' @export
report.varComp <- function(x,
                           ...,
                           outfile = NULL) {
  createReport(x = x, reportName = "varCompReport.Rnw", outfile = outfile, ...)
}
