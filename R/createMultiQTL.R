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
#' A pdf report will be created containing a summary of a multiQTL analysis.
#' Simultaneously the same report will be created as a tex file.
#'
#' @inheritParams report.AMMI
#'
#' @param x An object of class multiQTL.
#'
#' @return A pdf and tex report.
#'
#' @examples
#' ## Read the data
#' F2 <- qtl::read.cross(format="csv",
#'                       file = system.file("extdata",
#'                                         "F2_maize_practical3_ex2.csv",
#'                                         package = "RAP"),
#'                       genotypes = c("AA", "AB", "BB"),
#'                       alleles = c("A", "B"), estimate.map = FALSE)
#' ## Perform QTL detection using simple interval mapping.
#' QTLDet <- QTLDetect(cross = F2, trait = "trait", type = "SIM")
#' ## Fit a multi QTL model.
#' multiFit <- multiQTLFit(QTLDet)
#' \dontrun{
#' ## Create a pdf report summarizing results.
#' report(multiFit, outfile = "./testReports/reportMultiQTLFit.pdf")
#' }
#'
#' @export
report.multiQTL <- function(x,
                            ...,
                            outfile = NULL) {
  createReport(x = x, reportName = "multiQTLReport.Rnw",
               outfile = outfile, ...)
}


