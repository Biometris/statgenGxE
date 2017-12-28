#' S3 class QTLDet
#'
#' Function for creating objects of S3 class QTLDet.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param scores A data.frame containing the lod scores.
#' @param peaks A data.frame containing the peaks found.
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
                         peaks,
                         type,
                         cross) {
  QTLDet <- structure(list(scores = scores,
                           peaks = peaks,
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
#' Function for creating a manhattan plot for objects of class QTLDet.
#'
#' @param x an object of class QTLDet
#' @param ... not unused
#'
#' @import graphics grDevices
#' @export
plot.QTLDet <- function(x,
                        ...) {
  plot(x$scores, ylab = "LOD")
  if (x$type == "CIM") {
    qtl::add.cim.covar(x$scores)
  }
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


