#' Fit an additive multi QTL model
#'
#' @inheritParams QTLDetect
#'
#' @param qtlDet An object of class \code{\link{QTLDet}}
#'
#' @export
multiQTLFit <- function(qtlDet,
                        thr = 3,
                        window = 15,
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
