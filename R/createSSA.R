#' S3 class SSA
#'
#' Function for creating objects of S3 class Single Site Analysis (SSA).
#'
#' @param mMix a mixed model created using either asreml or lme4
#' @param mFix a fixed model created using either asreml or lme4
#' @param data an object of class TD containing the data on which mMix and mFix are based.
#' @param trait a character sting indicating the trait for which the analysis is done.
#' @param genotype a character sting indicating genotype column in the data.
#' @param repId a character sting indicating the replicates column in the data.
#' @param design a character string containing the design of the trial.
#' @param engine a character string containing the engine used to do the analysis.
#' @param x \code{R} object
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{summary.SSA}}
#'
#' @name SSA
NULL

#' @rdname SSA
#' @export
createSSA <- function(mMix,
                      mFix,
                      data,
                      trait = NULL,
                      genotype = NULL,
                      repId = NULL,
                      design = NULL,
                      engine = NULL) {
  SSA <- structure(list(mMix = mMix,
                        mFix = mFix,
                        data = data,
                        trait = trait,
                        genotype = genotype,
                        repId = repId,
                        design = design,
                        engine = engine),
                   class = "SSA")
  return(SSA)
}

#' @rdname SSA
#' @export
is.SSA <- function(x) {
  inherits(x, "SSA")
}

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
#' myDat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames = c("Env", "Genotype", "Rep", "Subblock", "Row", "Column"),
#'                      traitNames = "yield", env = "Env", rowSelect = "HEAT05",
#'                      colSelect = c("Env", "Genotype", "Rep", "Row", "Column", "yield"))
#' myTD <- createTD(data = myDat, genotype = "Genotype", env = "Env")
#' myModel <- ST.run.model(TD = myTD, design = "res.rowcol", trait = "yield",
#'                         repId = "Rep", rowId = "Row", colId = "Column")
#' summary(myModel)
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
  if (object$engine == "asreml" && !is.null(extr$predictionsSed) &&
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

#' Diagnostic Plots of Models
#'
#' This function draws four plots, a histogram of residuals, a normal Q-Q plot, a residuals
#' vs fitted values plot and an absolute residuals vs fitted values plot.
#'
#' @param x an object of class SSA.
#' @param ... Other graphical parameters (see \code{\link[lattice]{xyplot}} for details).
#' @param plotType character indiciting whether the fixed (\code{plotType = "fix"}) or mixed
#' (\code{plotType = "fix"}) should be plotted.
#'
#' @seealso \code{\link{createSSA}}
#'
#' @examples
#' myDat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames = c("Env","Genotype","Rep","Row","Column"),
#'                      traitNames = "yield", env = "Env", rowSelect = "HEAT05",
#'                      colSelect = c("Env","Genotype", "Rep", "Row", "Column", "yield"))
#' myTD <- createTD(data = myDat, genotype = "Genotype", env = "Env")
#' myModel <- ST.run.model(TD = myTD, design = "res.rowcol", trait = "yield",
#'                         repId = "Rep", rowId = "Row", colId = "Column",
#'                         tryspatial = "always")
#' plot(myModel, plotType = "fix")
#'
#' @export

plot.SSA <- function(x,
                     ...,
                     plotType = "fix") {
  if (plotType == "fix") {
    model <- x$mFix
  } else if (plotType == "mix") {
    model <- x$mMix
  } else {
    stop("plotType should either be fix or mix.")
  }
  # Diagnostic plots
  if (class(model) == "asreml") {
    resid <- model$residuals
    fitted <- model$fitted.values
  } else {
    resid <- residuals(model)
    fitted <- fitted(model)
  }
  trellisObj <- vector(mode = "list", length = 4)
  names(trellisObj) <- c("histogram", "qq", "residFitted", "absResidFitted")
  # Histogram of residuals
  trellisObj[["histogram"]] <- lattice::histogram(x = ~resid, xlab = "Residuals", ...)
  # Q-Q plot of residuals
  trellisObj[["qq"]] <- lattice::qqmath(~resid, xlab = "Normal quantiles",
                                        ylab = "Residuals", ...)
  # Residuals vs fitted values
  trellisObj[["residFitted"]] <- lattice::xyplot(resid ~ fitted,
                                                 panel = function(x, y, ...) {
                                                   lattice::panel.xyplot(x, y, ...,
                                                                         type = c("p", "g"))
                                                   lattice::panel.abline(h = 0)
                                                   lattice::panel.loess(x, y,
                                                                        col = "red", ...)
                                                 }, ylab = "Residuals",
                                                 xlab = "fitted values", ...)
  # Residuals vs fitted values
  trellisObj[["absResidFitted"]] <- lattice::xyplot(abs(resid) ~ fitted,
                                                    panel = function(x, y, ...) {
                                                      lattice::panel.xyplot(x, y, ...,
                                                                            type = c("p", "g"))
                                                      lattice::panel.loess(x, y,
                                                                           col = "red", ...)
                                                    }, ylab = "|Residuals|",
                                                    xlab = "fitted values", ...)
  adt <- lattice::trellis.par.get("add.text")
  xlb <- lattice::trellis.par.get("par.xlab.text")
  ylb <- lattice::trellis.par.get("par.ylab.text")
  zlb <- lattice::trellis.par.get("par.zlab.text")
  axt <- lattice::trellis.par.get("axis.text")
  syx <- lattice::trellis.par.get("plot.symbol")
  lattice::trellis.par.set("add.text", list(cex = 0.75))
  lattice::trellis.par.set("par.xlab.text", list(cex = 0.75))
  lattice::trellis.par.set("par.ylab.text", list(cex = 0.75))
  lattice::trellis.par.set("par.zlab.text", list(cex = 0.75))
  lattice::trellis.par.set("axis.text", list(cex = 0.75))
  lattice::trellis.par.set("plot.symbol", list(cex = 0.6))
  print(trellisObj[["histogram"]], position = c(0, 0.5, 0.5, 1), more = TRUE)
  print(trellisObj[["qq"]], position = c(0.5, 0.5, 1, 1), more = TRUE)
  suppressWarnings(print(trellisObj[["residFitted"]], position = c(0, 0, 0.5, 0.5),
                         more = TRUE))
  suppressWarnings(print(trellisObj[["absResidFitted"]], position = c(0.5, 0, 1, 0.5)))
  lattice::trellis.par.set("add.text", adt)
  lattice::trellis.par.set("par.xlab.text", xlb)
  lattice::trellis.par.set("par.ylab.text", ylb)
  lattice::trellis.par.set("par.zlab.text", zlb)
  lattice::trellis.par.set("axis.text", axt)
  lattice::trellis.par.set("plot.symbol", syx)
  invisible(trellisObj)
}


