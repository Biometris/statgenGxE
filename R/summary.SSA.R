#' Summarizing (\code{SSA})Model Fits
#'
#' \code{summary} method for class \code{SSA}.
#'
#' @param object An object of class \code{SSA}.
#' @param digits Integer. The number of significant digits to use when printing.
#' @param nBest Integer. The number of the best genotypes (sorted by either BLUEs or BLUPs)
#' is to print. If \code{NA}, print all of them.
#' @param sortBy A string specifying by which the genotypes will be sorted. The options
#' are \code{"BLUEs"}, \code{"BLUPs"} and \code{NA} (i.e. no sort).
#' @param naLast For controlling the treatment of NAs. If TRUE, missing values in the data
#' are put last; if FALSE, they are put first; if NA, they are removed.
#' @param decreasing Logical. Should the sort order be increasing or decreasing?
#' @param ... Further arguments passed to \code{\link[stats]{printCoefmat}}.
#'
#' @examples
#' mydat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames=c("Env","Genotype","Rep","Subblock","Row","Column"),
#'                      traitNames="yield", env ="Env", rowSelect="HEAT05",
#'                      colSelect=c("Env","Genotype","Rep","Row","Column","yield"))
#' mymodel <- ST.run.model(mydat, design="res.rowcol", trait="yield",
#'                         genotype="Genotype", rep="Rep", row="Row",
#'                         col="Column", tryspatial=NA)
#' summary(mymodel)
#'
#' @export
summary.SSA <- function(object,
                        digits = max(getOption("digits") - 2, 3),
                        nBest = 20,
                        sortBy = "BLUEs",
                        naLast = TRUE,
                        decreasing = TRUE,
                        ...) {
  # get summary stats for raw data
  data <- object$data
  trait <- object$trait
  stats <- summary.TD(object = data, trait = trait)
  stats <- na.omit(stats)
  attr(stats, "na.action") <- NULL
  # get predicted means (BLUEs & BLUPs)
  extr <- ST.extract(object)
  meanTab <- extr$stats[-1]
  if (!is.na(sortBy)) {
    if (sortBy == "BLUEs") {
      oList <- order(meanTab[["predicted (BLUEs)"]], na.last = naLast, decreasing = decreasing)
      meanTab <- meanTab[oList, ]
    } else {
      if (sortBy == "BLUPs") {
        oList <- order(meanTab[["predicted (BLUPs)"]], na.last = naLast, decreasing = decreasing)
        meanTab <- meanTab[oList, ]
      }
    }
  }
  if (!is.na(nBest)) {
    meanTab <- meanTab[1:nBest, ]
  }
  cat("Summary statistics:\n", "===================\n", sep = "")
  printCoefmat(stats, digits = digits, ...)
  cat("\nEstimated heritability\n", "======================\n", sep = "")
  cat("\nHeritability:", extr$heritability, "\n")
  cat("\nPredicted means (BLUEs & BLUPs)\n", "===============================\n", sep = "")
  if (!is.na(nBest)) {
    cat("Best", nBest,"genotypes\n")
  } else {
    cat("\n")
  }
  printCoefmat(meanTab, digits = digits, ...)
  if(object$engine == "asreml" && !is.null(extr$predictionsSed) &&
     !is.null(extr$predictionsLsd)) {
    cat("\nStandard Error of Difference (genotypes modelled as fixed effect)\n",
        "===================================================================\n", sep = "")
    sed <- as.data.frame(extr$predictionsSed)
    names(sed) <- "s.e.d."
    printCoefmat(sed, digits = digits, ...)
    cat("\nLeast Significant Difference (genotypes modelled as fixed effect)\n",
        "===================================================================\n", sep = "")
    lsd  <- as.data.frame(extr$predictionsLsd)
    names(lsd) <- "l.s.d."
    printCoefmat(lsd, digits = digits, ...)
  }
}
