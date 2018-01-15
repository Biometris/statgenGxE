#' Fit an additive multi QTL model
#'
#' An additive multi QTL model is fitted based on the peaks in the QTLDet object.
#' Fitting is done using the \code{\link[qtl]{fitqtl}} function in the qtl package.
#' After fitting the model backward elemination is done until all markers in
#' the model have a significant P-value.
#'
#' @inheritParams QTLDetect
#'
#' @param QTLDet An object of class \code{\link{QTLDet}}
#' @param selection An integer string indicating whether backward selection should
#' be applied or no selection at all.
#'
#' @return An object of class \code{\link{multiQTL}}
#'
#' @seealso \code{\link[qtl]{fitqtl}}
#'
#' @references Broman et al. (2003) R/qtl: QTL mapping in experimental crosses.
#' Bioinformatics 19:889-890
#'
#' @examples
#' ## Read the data
#' F2 <- qtl::read.cross(format="csv",
#'                       file = system.file("extdata", "F2_maize_practical3_ex2.csv",
#'                       package = "RAP"),
#'                       genotypes = c("AA", "AB", "BB"),
#'                       alleles = c("A", "B"), estimate.map = FALSE)
#' ## Perform QTL detection using simple interval mapping.
#' QTLDet <- QTLDetect(cross = F2, trait = "trait", type = "SIM")
#' ## Fit the multi QTL model.
#' multiFit <- multiQTLFit(QTLDet)
#' ## Create a report.
#' report(multiFit, outfile = "./testReports/reportMultiQTLFit.pdf")
#'
#' @export
multiQTLFit <- function(QTLDet,
                        selection = c("backward", "none"),
                        ...) {
  ## Checks
  if (!is.QTLDet(QTLDet)) {
    stop("QTLDet should be an object of class QTLDet.\n")
  }
  selection <- match.arg(selection)
  ## Create an object of class qtl for modelling.
  qtl <- qtl::makeqtl(QTLDet$cross, chr = QTLDet$peaks$chr,
                      pos = QTLDet$peaks$pos, what = "prob")
  ## Construct model formula.
  qtlForm <- paste("y ~ ", paste(qtl$altname, collapse = "+"))
  ## Fit full model with all markers.
  qtlFit <- qtl::fitqtl(QTLDet$cross, qtl = qtl, formula = qtlForm,
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
      qtlFit <- qtl::fitqtl(QTLDet$cross, qtl = qtl, formula = qtlForm,
                            method = "hk", get.ests = TRUE,
                            dropone = TRUE, ...)
    }
  }
  ## Rename rows to contain markerNames instead of chr@pos.
  rownames(qtlFit$result.drop) <- qtlPosToName(rownames(qtlFit$result.drop),
                                               cross = QTLDet$cross)$chrNames
  estNames <- qtlPosToName(names(qtlFit$ests$ests)[-1], cross = QTLDet$cross)
  names(qtlFit$ests$ests)[-1] <- paste0(estNames$chrNames, estNames$ext)
  ## Create multiQTL object.
  multiQtl <- createMultiQTL(qtl = qtlFit, QTLDet = QTLDet, selection = selection)
  return(multiQtl)
}


