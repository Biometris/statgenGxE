#' Summarizing (\code{SSA})Model Fits
#'
#' \code{summary} method for class \code{SSA}.
#'
#' @param object An object of class \code{SSA}.
#' @param digits Integer. The number of significant digits to use when printing.
#' @param n.best Integer. The number of the best genotypes (sorted by either BLUEs or BLUPs) is to print. If \code{NA}, print all of them.
#' @param sortBy A string specifying by which the genotypes will be sorted. The options are \code{"BLUEs"}, \code{"BLUPs"} and \code{NA} (i.e. no sort).
#' @param na.last For controlling the treatment of NAs. If TRUE, missing values in the data are put last; if FALSE, they are put first; if NA, they are removed.
#' @param decreasing Logical. Should the sort order be increasing or decreasing?
#' @param ... Further arguments passed to \code{\link[stats]{printCoefmat}}.
#' @details summarize the results obtained from a single site analysis.
#' @examples
#' mydat <- ST.read.csv(file.path(path.package("RAP"),"SB_yield.csv"),
#'                      factor.names=c("Env","Genotype","Rep","Subblock","Row","Column"),
#'                      trait.names="yield", env ="Env", rowSelect="HEAT05",
#'                      colSelect=c("Env","Genotype","Rep","Row","Column","yield"))
#' mymodel <- ST.run.model(mydat, design="res.rowcol", trait="yield",
#'                         genotype="Genotype", rep="Rep", row="Row",
#'                         col="Column", tryspatial=NA)
#' ## This crashes
#' ## summary(mymodel)
#'
#' @method summary SSA
#' @export
summary.SSA <- function(object, digits = max(getOption("digits") - 2, 3), n.best = 20,
    sortBy = "BLUEs", na.last=TRUE, decreasing = TRUE, ...)
{
  # get summary stats for raw data
  data <- object$Data
  trait <- attr(object, "Trait")
  stats <- ST.summary.trait(data, trait, printTable = F)
  stats <- na.omit(stats)
  attr(stats,"na.action") <- NULL

  # get predicted means (BLUEs & BLUPs)
  extr <- ST.extract(object)
  meanTab <- extr$Stats[-1]
  if (!is.na(sortBy)){
    if (sortBy == "BLUEs"){
      olist <- order(meanTab[["predicted (BLUEs)"]], na.last=na.last, decreasing=decreasing)
      meanTab <- meanTab[olist,]
    }else{
      if (sortBy == "BLUPs"){
        olist <- order(meanTab[["predicted (BLUPs)"]], na.last=na.last, decreasing=decreasing)
        meanTab <- meanTab[olist,]
      }
    }
  }
  if (!is.na(n.best))
    meanTab <- meanTab[1:n.best,]

  cat("Summary statistics:\n", "===================\n", sep = "")
  printCoefmat(stats, digits = digits,...)
  cat("\nEstimated heritability\n", "======================\n", sep = "")
  cat("\nHeritability:", extr$heritability, "\n")

  cat("\nPredicted means (BLUEs & BLUPs)\n", "===============================\n", sep = "")
  if (!is.na(n.best))
    cat("Best", n.best,"genotypes\n")
  else
    cat("\n")
  printCoefmat(meanTab, digits = digits,...)

  if(attr(object, "Engine") == "asreml" && !is.null(extr$predictions.sed) && !is.null(extr$predictions.lsd)){
    cat("\nStandard Error of Difference (genotypes modelled as fixed effect)\n", "===================================================================\n", sep = "")
    sed  <- as.data.frame(extr$predictions.sed)
    names(sed) <- "s.e.d."
    printCoefmat(sed, digits = digits,...)

    cat("\nLeast Significant Difference (genotypes modelled as fixed effect)\n", "===================================================================\n", sep = "")
    lsd  <- as.data.frame(extr$predictions.lsd)
    names(lsd) <- "l.s.d."
    printCoefmat(lsd, digits = digits,...)
  }
}
