#' S3 class QTLDet
#'
#' Function for creating objects of S3 class QTLDet.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param scores A data.frame containing the lod scores.
#' @param peaks A data.frame containing the peaks found.
#' @param type A character string indicating the type of QTL detection performed.
#' @param cross An object of class cross in the \code{qtl} package.
#' @param trait A character string indicating the trait for which the analysis
#' is done.
#' @param info A list containing information on the settings used for
#' QTL detection, i.e. step, threshold and window.
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
                         cross,
                         trait,
                         info) {
  QTLDet <- structure(list(scores = scores,
                           peaks = peaks,
                           type = type,
                           cross = cross,
                           trait = trait,
                           info = info),
                      class = "QTLDet")
  attr(QTLDet, which = "timestamp") <- Sys.time()
  return(QTLDet)
}

#' @export
print.QTLDet <- function(x, ...) {
  cat("Peaks\n")
  cat("=====\n")
  print(x$peaks)
}

#' @export
summary.QTLDet <- function(object, ...) {
  print(object, ...)
}

#' Plot function for class QTLDet
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
  plot(x$scores, ylab = "LOD", ...)
  if (x$type == "CIM") {
    qtl::add.cim.covar(x$scores)
  }
}

#' Report method for class QTLDet
#'
#' A pdf report will be created containing a summary of a QTLDet analysis.
#' Simultaneously the same report will be created as a tex file.
#'
#' @inheritParams report.AMMI
#'
#' @param x An object of class QTLDet.
#'
#' @return A pdf and tex report.
#'
#' @examples
#' F2 <- qtl::read.cross(format="csv",
#'                       file = system.file("extdata",
#'                                         "F2_maize_practical3_ex2.csv",
#'                                         package = "RAP"),
#'                       genotypes = c("AA", "AB", "BB"),
#'                       alleles = c("A", "B"), estimate.map = FALSE)
#' ## Perform a composite interval mapping for detecting QTLs.
#' QTLDet <- QTLDetect(cross = F2, trait = "trait", type = "CIM")
#' \dontrun{
#' ## Create a pdf report summarizing the results.
#' report(QTLDet, outfile = "./testReports/reportQTLDectection.pdf")
#' }
#'
#' @export
report.QTLDet <- function(x,
                          ...,
                          outfile = NULL) {
  createReport(x = x, reportName = "QTLDetReport.Rnw",
               outfile = outfile, ...)
}


