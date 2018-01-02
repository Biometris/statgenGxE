#' Function for doing quality control
#'
#' Function for doing quality control of input before doing a QTL analysis.
#'
#' @param cross An object of class cross created by the qtl package
#'
#' @export
QTLMapQC <- function(cross) {



  ## Re-estimate map based on the observed markers
  new_map <- qtl::est.map(cross, error.prob = 0.001)
  plot(new_map)
  qtl::plot.map(cross, new_map)





}

#' Report method for class cross
#'
#' A pdf report will be created containing a summary of cross analysis.
#' Simultaneously the same report will be created as a tex file.
#'
#' @param x an object of class cross.
#' @param ... further arguments passed on from other functions - not used yet.
#' @param outfile a character string, the name and location of the output .pdf and .tex
#' file for the report. If \code{NULL} a report will be created in the current working
#' directory.
#'
#' @examples
#' ## Read the data
#' F2 <- qtl::read.cross(format="csv",
#'                       file = system.file("extdata", "F2_maize_practical3_ex2.csv", package = "RAP"),
#'                       genotypes = c("AA", "AB", "BB"),
#'                       alleles = c("A", "B"), estimate.map = FALSE)
#' report(F2, outfile = "./testReports/reportCross.pdf")
#'
#' @export
report.cross <- function(x,
                          ...,
                          outfile = NULL) {
  createReport(x = x, reportName = "crossReport.Rnw",
               outfile = outfile, ...)
}





