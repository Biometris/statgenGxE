#' S3 class QTLDet
#'
#' Function for creating objects of S3 class QTLDet.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param scores A data.frame containing the lod scores.
#' @param type A character string indicating the type of QTLDetection performed.
#' @param cross An object of class cross in the \code{qtl} package.
#' @param x an \code{R} object
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.QTLDet}}, \code{\link{report.QTLDet}}
#'
#' @name QTLDet
NULL

#' @rdname QTLDet
#' @export
createQTLDet <- function(scores,
                         type,
                         cross) {
  QTLDet <- structure(list(scores = scores,
                           type = type,
                           cross = cross),
                      class = "QTLDet")
  attr(QTLDet, which = "timestamp") <- Sys.time()
  return(QTLDet)
}

#' @rdname QTLDet
#' @export
is.QTLDet <- function(x) {
  inherits(x, "QTLDet")
}

#' @export
print.QTLDet <- function(x, ...) {

}

#' @export
summary.QTLDet <- function(object, ...) {
  print(object, ...)
}

#' Plot Function for Class QTLDet
#'
#' Function for creating scatter, line and trellis plots for objects of class QTLDet.
#'
#' @param x an object of class QTLDet
#' @param ... not unused
#' @param plotType a character vector indicating which plot(s) will be drawn. Possible values
#' "scatter", "line"  and "trellis" for creating a scatter plot of sensitivities, a plot of
#' fitted lines for each genotype and a trellis plot of the individual genotype slopes
#' respectively.
#' @param sortBySens A character string specifying whether the results are to be sorted
#' in an increasing (or decreasing) order of sensitivities.
#' By default, \code{sortBySens = "ascending"}. Other options are "descending" and NA.

#' @return Plots as described in \code{plotType}
#'
#' @import graphics grDevices
#' @export
plot.QTLDet <- function(x,
                        ...,
                        plotType = c("scatter", "line", "trellis"),
                        sortBySens = "ascending") {

}


#' Report method for class QTLDet
#'
#' A pdf report will be created containing a summary of QTLDet analysis.
#' Simultaneously the same report will be created as a tex file.
#'
#' @param x an object of class QTLDet.
#' @param ... further arguments passed on from other functions - not used yet.
#' @param outfile a character string, the name and location of the output .pdf and .tex
#' file for the report. If \code{NULL} a report will be created in the current working
#' directory.
#'
#' @export
report.QTLDet <- function(x,
                          ...,
                          outfile = NULL) {
  createReport(x = x, reportName = "QTLDetReport.Rnw",
               outfile = outfile, ...)
}


