#' S3 class multiQTL
#'
#' Function for creating objects of S3 class multiQTL.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param qtl A fitted multi QTL model.
#' @param x an \code{R} object
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.multiQTL}}, \code{\link{report.multiQTL}}
#'
#' @name multiQTL
NULL

#' @rdname multiQTL
#' @export
createMultiQTL <- function(qtl) {
  multiQTL <- structure(list(qtl = qtl),
                        class = "multiQTL")
  attr(multiQTL, which = "timestamp") <- Sys.time()
  return(multiQTL)
}

#' @rdname multiQTL
#' @export
is.multiQTL <- function(x) {
  inherits(x, "multiQTL")
}

#' @export
print.multiQTL <- function(x, ...) {

}

#' @export
summary.multiQTL <- function(object, ...) {
  print(object, ...)
}

#' Plot Function for Class multiQTL
#'
#' Function for creating scatter, line and trellis plots for objects of class multiQTL.
#'
#' @param x an object of class multiQTL
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
plot.multiQTL <- function(x,
                          ...,
                          plotType = c("scatter", "line", "trellis"),
                          sortBySens = "ascending") {

}


#' Report method for class multiQTL
#'
#' A pdf report will be created containing a summary of multiQTL analysis.
#' Simultaneously the same report will be created as a tex file.
#'
#' @param x an object of class multiQTL.
#' @param ... further arguments passed on from other functions - not used yet.
#' @param outfile a character string, the name and location of the output .pdf and .tex
#' file for the report. If \code{NULL} a report will be created in the current working
#' directory.
#'
#' @export
report.multiQTL <- function(x,
                            ...,
                            outfile = NULL) {
  createReport(x = x, reportName = "multiQTLReport.Rnw",
               outfile = outfile, ...)
}


