#' S3 class varComp
#'
#' Function for creating objects of S3 class varComp.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param model The best fitted model.
#' @param choice A character string indicating the best fitted model.
#' @param summary A data.frame with a summary of the fitted models.
#' @param vcov The covariance matrix of the best fitted model.
#' @param x an \code{R} object
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.VarComp}}, \code{\link{report.VarComp}}
#'
#' @name varComp
NULL

#' @rdname varComp
#' @export
createVarComp <- function(model,
                          choice,
                          summary,
                          vcov,
                          criterion) {
  varComp <- structure(list(model = model,
                            choice = choice,
                            summary = summary,
                            vcov = vcov,
                            criterion = criterion),
                       class = "varComp")
  return(varComp)
}

#' @rdname varComp
#' @export
is.varComp <- function(x) {
  inherits(x, "varComp")
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
#' @param x an object of class varComp
#' @param ... not unused
#'
#' @import stats
#' @export
plot.varComp <- function(x, ...) {
  corMat <- cov2cor(x$vcov)
  heatmap(corMat, Rowv = NA, symm = TRUE)
}

#' Report method for class varComp
#'
#' A pdf report will be created containing a summary of selection of the
#' best variance model.
#' Simultaneously the same report will be created as a tex file.
#'
#' @param x an object of class varComp.
#' @param ... further arguments passed on from other functions - not used yet.
#' @param outfile a character string, the name and location of the output .pdf and .tex
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


