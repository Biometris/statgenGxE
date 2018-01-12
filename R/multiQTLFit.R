#' Fit an additive multi QTL model
#'
#' @inheritParams QTLDetect
#'
#' @param qtlDet An object of class \code{\link{QTLDet}}
#' @param selection An integer string indicating whether backward selection should
#' be applied or no selection at all.
#'
#' @return An object of class \code{\link{multiQTL}}
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
                        selection = c("backward", "none"),
                        ...) {
  ## Checks
  if (!is.QTLDet(qtlDet)) {
    stop("qtlDet should be an object of class qtlDET.\n")
  }
  selection <- match.arg(selection)
  ## Create an object of class qtl for modelling.
  qtl <- qtl::makeqtl(qtlDet$cross, chr = qtlDet$peaks$chr,
                      pos = qtlDet$peaks$pos, what = "prob")
  ## Construct model formula.
  qtlForm <- paste("y ~ ", paste(qtl$altname, collapse = "+"))
  ## Fit full model with all markers.
  qtlFit <- qtl::fitqtl(qtlDet$cross, qtl = qtl, formula = qtlForm,
                        method = "hk", get.ests = TRUE, dropone = TRUE,
                        ...)
  if (selection == "backward") {
    ## While there are markers with Pvalue remove the one with the
    ## and refit the model without this marker highest value.
    while (any(qtlFit$result.drop[, "Pvalue(F)"] > 0.05)) {
      ## Drop the marker with the highest P-value.
      qtl <- qtl::dropfromqtl(qtl, qtl.name =
                                names(which.max(qtlFit$result.drop[, "Pvalue(F)"])))
      ## Rebuild the fitting formula for the remaining markers.
      qtlForm <- paste("y ~ ", paste(qtl$altname, collapse = "+"))
      ## Refit the model.
      qtlFit <- qtl::fitqtl(qtlDet$cross, qtl = qtl, formula = qtlForm,
                            method = "hk", get.ests = TRUE,
                            dropone = TRUE, ...)
    }
  }
  multiQtl <- createMultiQTL(qtl = qtlFit)
  return(multiQtl)
}


