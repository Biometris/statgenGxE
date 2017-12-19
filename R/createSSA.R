#' S3 class SSA
#'
#' Function for creating objects of S3 class Single Site Analysis (SSA).
#'
#' @param mRand a mixed model created using either asreml or lme4
#' @param mFix a fixed model created using either asreml or lme4
#' @param data an object of class TD containing the data on which mRand and mFix are based.
#' @param traits a character vector indicating the traits for which the analysis is done.
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
createSSA <- function(mRand,
                      mFix,
                      data,
                      traits = NULL,
                      design = NULL,
                      spatial = NULL,
                      engine = NULL) {
  SSA <- structure(list(mRand = mRand,
                        mFix = mFix,
                        data = data,
                        traits = traits,
                        design = design,
                        spatial = spatial,
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
#' @param object an object of class \code{SSA}.
#' @param trait a string indicating the trait to summarize over. If
#' \code{trait = NULL} and only one trait is modelled this trait is summarized.
#' @param digits an integer. The number of significant digits to use when
#' printing.
#' @param nBest an integer. The number of the best genotypes (sorted by either
#' BLUEs or BLUPs) is to print. If \code{NA}, print all of them.
#' @param sortBy a string specifying by which the genotypes will be sorted.
#' The options are \code{"BLUEs"}, \code{"BLUPs"} and \code{NA} (i.e. no sort).
#' @param naLast for controlling the treatment of NAs. If TRUE, missing values
#' in the data are put last; if FALSE, they are put first; if NA, they are removed.
#' @param decreasing should the sort order be decreasing?
#' @param ... further arguments passed to \code{\link[stats]{printCoefmat}}.
#'
#' @examples
#' data(TDHeat05)
#' myModel <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield")
#' summary(myModel)
#'
#' @export
summary.SSA <- function(object,
                        trait = NULL,
                        digits = max(getOption("digits") - 2, 3),
                        nBest = 20,
                        sortBy = if (!is.null(object$mFix)) "BLUEs" else "BLUPs",
                        naLast = TRUE,
                        decreasing = TRUE,
                        ...) {
  ## Checks.
  if (is.null(trait) && length(object$traits) > 1) {
    stop("No trait provided but multiple traits found in SSA object.\n")
  }
  if (!is.null(trait) && (!is.character(trait) || length(trait) > 1 ||
                          !trait %in% colnames(object$data))) {
    stop("Trait has to be a single character string defining a column in data.\n")
  }
  ## get summary stats for raw data
  TD <- object$data
  if (is.null(trait)) {
    trait <- object$traits
  }
  stats <- summary.TD(object = TD, traits = trait)
  stats <- na.omit(stats)
  attr(stats, "na.action") <- NULL
  ## get predicted means (BLUEs & BLUPs).
  extr <- STExtract(object)
  joinList <- Filter(f = Negate(f = is.null),
                     x = list(extr$BLUEs, extr$seBLUEs,
                              extr$BLUPs, extr$seBLUPs))
  meanTab <- Reduce(f = function(x, y) {
    dplyr::full_join(x, y, all = TRUE, by = "genotype")
  }, x = joinList)
  colnames(meanTab) <- c("genotype",
                         if (!is.null(extr$BLUEs)) c("predictedBLUEs", "se predictedBLUEs"),
                         if (!is.null(extr$BLUPs)) c("predictedBLUPs", "se predictedBLUPs"))
  rownames(meanTab) <- meanTab$genotype
  if (!is.na(sortBy)) {
    if (sortBy == "BLUEs") {
      oList <- order(meanTab$predictedBLUEs, na.last = naLast, decreasing = decreasing)
      meanTab <- meanTab[oList, ]
    } else {
      if (sortBy == "BLUPs") {
        oList <- order(meanTab$predictedBLUPs, na.last = naLast, decreasing = decreasing)
        meanTab <- meanTab[oList, ]
      }
    }
  }
  if (!is.na(nBest)) {
    meanTab <- meanTab[1:nBest, ]
  }
  cat("Summary statistics:\n", "===================\n", sep = "")
  printCoefmat(stats, digits = digits, ...)
  if (!is.null(object$mRand)) {
    cat("\nEstimated heritability\n", "======================\n", sep = "")
    cat("\nHeritability:", extr$heritability, "\n")
  }
  cat("\nPredicted means (BLUEs & BLUPs)\n", "===============================\n", sep = "")
  if (!is.na(nBest)) {
    cat("Best", nBest,"genotypes\n")
  } else {
    cat("\n")
  }
  printCoefmat(meanTab[, -1], digits = digits, ...)
  if (object$engine == "asreml" && !is.null(extr$sed) &&
      !is.null(extr$lsd)) {
    cat("\nStandard Error of Difference (genotypes modelled as fixed effect)\n",
        "===================================================================\n", sep = "")
    sed <- as.data.frame(extr$lsd)
    names(sed) <- "s.e.d."
    printCoefmat(sed, digits = digits, ...)
    cat("\nLeast Significant Difference (genotypes modelled as fixed effect)\n",
        "===================================================================\n", sep = "")
    lsd  <- as.data.frame(extr$lsd)
    names(lsd) <- "l.s.d."
    printCoefmat(lsd, digits = digits, ...)
  }
  invisible(meanTab)
}

#' Diagnostic Plots of Models
#'
#' This function draws four plots, a histogram of residuals, a normal Q-Q plot, a residuals
#' vs fitted values plot and an absolute residuals vs fitted values plot.
#'
#' @inheritParams summary.SSA
#'
#' @param x an object of class SSA.
#' @param ... Other graphical parameters (see \code{\link[lattice]{xyplot}} for details).
#' @param what character indiciting whether the fixed (\code{what = "fix"}) or mixed
#' (\code{what = "fix"}) model should be plotted.
#' @param plotType not used for now.
#'
#' @seealso \code{\link{createSSA}}
#'
#' @examples
#' data(TDHeat05)
#' myModel <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield",
#'                       trySpatial = "always")
#' plot(myModel, plotType = "fix")
#'
#' @export

plot.SSA <- function(x,
                     ...,
                     trait = NULL,
                     what = "fixed",
                     plotType = "base") {
  ## Checks.
  if (is.null(trait) && length(x$traits) > 1) {
    stop("No trait provided but multiple traits found in SSA x\n")
  }
  if (!is.null(trait) && (!is.character(trait) || length(trait) > 1 ||
                          !trait %in% colnames(x$data))) {
    stop("Trait has to be a single character string defining a column in data.\n")
  }
  if (!is.character(what) || length(what) > 1 ||
      !what %in% c("fixed", "random")) {
    stop("what should be fixed or random.\n")
  }
  if (is.null(trait)) {
    trait <- x$traits
  }
  if (what == "fixed") {
    model <- x$mFix[[trait]]
  } else if (what == "random") {
    model <- x$mRand[[trait]]
  }
  engine <- x$engine
  ## Diagnostic plots.
  resid <- residuals(model)
  fitted <- fitted(model)
  trellisObj <- setNames(vector(mode = "list", length = 4),
                         c("histogram", "qq", "residFitted", "absResidFitted"))
  # Histogram of residuals
  trellisObj[["histogram"]] <- lattice::histogram(x = ~resid, xlab = "Residuals", ...)
  # Q-Q plot of residuals
  trellisObj[["qq"]] <- lattice::qqmath(~resid, xlab = "Normal quantiles",
                                        ylab = "Residuals", ...)
  # Residuals vs fitted values
  trellisObj[["residFitted"]] <-
    lattice::xyplot(resid ~ fitted,
                    panel = function(x, y, ...) {
                      lattice::panel.xyplot(x, y, ...,
                                            type = c("p", "g"))
                      lattice::panel.abline(h = 0)
                      lattice::panel.loess(x, y,
                                           col = "red", ...)
                    }, ylab = "Residuals",
                    xlab = "Fitted values", ...)
  # Residuals vs fitted values
  trellisObj[["absResidFitted"]] <-
    lattice::xyplot(abs(resid) ~ fitted,
                    panel = function(x, y, ...) {
                      lattice::panel.xyplot(x, y, ...,
                                            type = c("p", "g"))
                      lattice::panel.loess(x, y,
                                           col = "red", ...)
                    }, ylab = "|Residuals|",
                    xlab = "Fitted values", ...)
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

#' Create a report with basic results.
#'
#' A pdf report will be created containing a summary of the results of the model.
#' Simultaneously the same report will be created as a tex file.
#'
#' @param x an object of class SSA.
#' @param outfile a character string, the name and location of the output .pdf and .tex
#' file for the report. If \code{NULL} a report will be created in the current working
#' directory.
#'
#' @export
report.SSA <- function(x, ..., outfile = NULL) {
  if (!is.null(outfile)) {
    if (!is.character(outfile) || length(outfile) > 1 ||
        !dir.exists(dirname(outfile)) || file_ext(outfile) != "pdf") {
      stop("invalid output filename provided.\n")
    }
  } else {
    timeStamp <- format(Sys.time(), "%Y%m%d%H%M%S")
    outfile <- paste0(getwd(), "/modelReport_", timeStamp, ".pdf")
  }
  outBase <- substring(basename(outfile), first = 1,
                       last = nchar(basename(outfile)) - 3)
  outTex <- paste0(system.file("latex", package = "RAP"), "/", outBase, "tex")
  reportFile <- system.file("latex", "modelReport.Rnw", package = "RAP")
  knitr::knit(input = reportFile, output = outTex, quiet = FALSE)
  system2(command = Sys.which("pdflatex"),
          args = c(paste0(' -output-directory="', dirname(outfile), '"'),
                   "-interaction=nonstopmode",
                   paste0(' "',  outTex, '"')))
  system2(command = Sys.which("pdflatex"),
          args = c(paste0(' -output-directory="', dirname(outfile), '"'),
                   "-interaction=nonstopmode",
                   paste0(' "',  outTex, '"')))
  for (extension in c("aux", "log", "out", "toc", "xwm")) {
    unlink(paste0(dirname(outfile), "/", outBase, extension))
  }
  invisible(file.rename(from = outTex,
                        to = paste0(dirname(outfile), "/", basename(outTex))))
}

