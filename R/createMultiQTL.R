#' S3 class multiQTL
#'
#' Function for creating objects of S3 class multiQTL.\cr
#' \code{\link{print}}, \code{\link{summary}} and \code{\link{report}}
#' methods are available.
#'
#' @param qtl A fitted multi QTL model.
#' @param QTLDet The object of class \code{\link{QTLDet}} used as base for fitting
#' the QTL model.
#' @param selection A character string indictating the type of selection used for
#' selecting the markers in the final model.
#' @param thr A numerical value indicating the threshold for dropping terms in
#' the backwards elemination process.
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{report.multiQTL}}
#'
#' @name multiQTL
NULL

#' @rdname multiQTL
#' @export
createMultiQTL <- function(qtl,
                           QTLDet,
                           selection,
                           thr) {
  multiQTL <- structure(list(qtl = qtl,
                             QTLDet = QTLDet,
                             selection = selection,
                             thr = thr),
                        class = "multiQTL")
  attr(multiQTL, which = "timestamp") <- Sys.time()
  return(multiQTL)
}

#' @export
print.multiQTL <- function(x, ...) {
  summary(x, ...)
}

#' @export
summary.multiQTL <- function(object, ...) {
  summary(object$qtl, ...)
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


