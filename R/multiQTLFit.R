#' Fit an additive multi QTL model
#'
#' @inheritParams QTLDetect
#'
#' @param qtlDet An object of class \code{\link{QTLDet}}
#'
#' @examples
#' ## Read the data
#' F2 <- qtl::read.cross(format="csv",
#'                       file = system.file("extdata", "F2_maize_practical3_ex2.csv",
#'                       package = "RAP"),
#'                       genotypes = c("AA", "AB", "BB"),
#'                       alleles = c("A", "B"), estimate.map = FALSE)
#' QTLDet <- QTLDetect(F2, "SIM")
#' multiFit <- multiQTLFit(QTLDet)
#' report(multiFit, outfile = "./testReports/reportMultiQTLFit.pdf")
#'
#' @export
multiQTLFit <- function(qtlDet,
                        ...) {
  ## Create an object of class qtl for modelling.
  qtl <- qtl::makeqtl(qtlDet$cross, chr = qtlDet$peaks$chr,
                      pos = qtlDet$peaks$pos, what = "prob")
  ## Construct model formula.
  qtlForm <- paste("y ~ ", paste(qtl$altname, collapse = "+"))
  qtlFinal <- qtl::fitqtl(qtlDet$cross, qtl = qtl, formula = qtlForm,
                          method = "hk", get.ests = TRUE, dropone = TRUE,
                          ...)
  multiQtl <- createMultiQTL(qtl = qtlFinal)
  return(multiQtl)
}
