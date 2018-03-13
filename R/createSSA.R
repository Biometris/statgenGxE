#' S3 class SSA
#'
#' Function for creating objects of S3 class Single Site Analysis (SSA).\cr
#' \code{\link{summary}}, \code{\link{plot}} and \code{\link{report}}
#' methods are available.
#'
#' @param mRand A list of models with fitted with genotype as random effect.
#' @param mFix A list of models fitted with genotype as fixed effect.
#' @param TD An object of class \code{\link{TD}} containing the data on which
#' \code{mRand} and \code{mFix} are based.
#' @param traits A character vector indicating the traits for which the analysis
#' is done.
#' @param design A character string containing the design of the trial.
#' (see \code{\link{STRunModel}} for the possible designs).
#' @param spatial A character string indicating the spatial part of the model.
#' \code{FALSE} if no spatial design has been used.
#' @param engine A character string containing the engine used for the analysis.
#' @param predicted A character string indicating the variable that has been
#' predicted.
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{summary.SSA}}, \code{\link{plot.SSA}}, \code{\link{report.SSA}}
#'
#' @name SSA
NULL

#' @rdname SSA
#' @export
createSSA <- function(mRand,
                      mFix,
                      TD,
                      traits = NULL,
                      design = NULL,
                      spatial = NULL,
                      engine = NULL,
                      predicted = "genotype") {
  SSA <- structure(list(mRand = mRand,
                        mFix = mFix,
                        TD = TD,
                        traits = traits,
                        design = design,
                        spatial = spatial,
                        engine = engine),
                   class = "SSA")
  attr(SSA, which = "timestamp") <- Sys.time()
  attr(SSA, which = "predicted") <- predicted
  return(SSA)
}

#' Summarizing objects of class \code{SSA}
#'
#' \code{summary} method for class \code{SSA}.
#'
#' @param object An object of class \code{SSA}.
#' @param trait A character string indicating the trait to summarize. If
#' \code{trait = NULL} and only one trait is modelled this trait is summarized.
#' @param digits An integer indicating the number of significant digits for
#' printing.
#' @param nBest An integer indicating the number of the best genotypes (sorted
#' by either BLUEs or BLUPs) to print. If \code{NA} all genotypes will be
#' printed.
#' @param sortBy A character string specifying how the genotypes will be sorted.
#' Either \code{"BLUEs"}, \code{"BLUPs"} or \code{NA} (i.e. no sorting).
#' @param naLast Should missing values in the data be put last when sorting?
#' @param decreasing Should the sort order be decreasing?
#' @param ... Further arguments passed to \code{\link[stats]{printCoefmat}}.
#'
#' @examples
#' ## Run a single trait analysis using SpATS.
#' myModel <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield")
#' ## Print a summary of the fitted model.
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
                          !trait %in% colnames(object$TD))) {
    stop("Trait has to be a single character string defining a column in TD.\n")
  }
  ## get summary stats for raw data
  TD <- object$TD
  if (is.null(trait)) {
    trait <- object$traits
  }
  stats <- summary.TD(object = TD, traits = trait)
  ## get predicted means (BLUEs + BLUPs).
  extr <- STExtract(object)
  ## Merge BLUEs, BLUPs + SE.
  joinList <- Filter(f = Negate(f = is.null),
                     x = list(extr$BLUEs, extr$seBLUEs,
                              extr$BLUPs, extr$seBLUPs))
  meanTab <- Reduce(f = function(x, y) {
    dplyr::full_join(x, y, all = TRUE, by = "genotype")
  }, x = joinList) %>%
    ## Move genotype to rowname for proper printing with printCoefMat
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "genotype") %>%
    ## Set colnames. Because of duplicate colname SE no selection on columns can
    ## be done anymore after this.
    setNames(c(if (!is.null(extr$BLUEs)) c("BLUEs", "SE"),
               if (!is.null(extr$BLUPs)) c("BLUPs", "SE")))
  if (!is.na(sortBy)) {
    ## Sort by sortBy with options from input params.
    oList <- order(meanTab[[sortBy]], na.last = naLast, decreasing = decreasing)
    meanTab <- meanTab[oList, ]
  }
  if (!is.na(nBest)) {
    ## Extract the n best genotypes.
    meanTab <- meanTab[1:nBest, ]
  }
  cat("Summary statistics:", "\n===================\n")
  ## Print stats using printCoefMat for a nicer layout.
  printCoefmat(stats, digits = digits, ...)
  if (!is.null(object$mRand)) {
    cat("\nEstimated heritability", "\n======================\n")
    cat("\nHeritability:", extr$heritability, "\n")
  }
  meansTxt <- paste(c(if (!is.null(extr$BLUEs)) "BLUEs",
                      if (!is.null(extr$BLUPs)) "BLUPs"), collapse = " & ")
  cat(paste0("\nPredicted means (", meansTxt, ")"), "\n===============================\n")
  if (!is.na(nBest)) {
    cat("Best", nBest,"genotypes\n")
  } else {
    cat("\n")
  }
  ## Print meanTab using printCoefMat for a nicer layout.
  printCoefmat(meanTab, digits = digits, ...)
  if (object$engine == "asreml" && !is.null(extr$sed) &&
      !is.null(extr$lsd)) {
    cat("\nStandard Error of Difference (genotype modelled as fixed effect)",
        "\n================================================================\n")
    sed <- data.frame("s.e.d" = extr$sed)
    printCoefmat(sed, digits = digits, ...)
    cat("\nLeast Significant Difference (genotype modelled as fixed effect)",
        "\n================================================================\n")
    lsd  <- data.frame("l.s.d." = extr$lsd)
    printCoefmat(lsd, digits = digits, ...)
  }
  invisible(meanTab)
}

#' Plot function for class SSA
#'
#' This function draws either four base plots:
#' \itemize{
#' \item{A histogram of the residuals}
#' \item{A normal Q-Q plot}
#' \item{A residuals vs fitted values plot}
#' \item{An absolute residuals vs fitted values plot}
#' }
#' or five or (in case SpATS is used for modelling) six spatial plots:
#' \itemize{
#' \item{A spatial plot of the raw data}
#' \item{A spatial plot of the fitted data}
#' \item{A spatial plot of the residuals}
#' \item{A spatial plot of the estimated spatial trend (SpATS only)}
#' \item{A spatial plot of the BLUEs or BLUPs}
#' \item{A histogram of the BLUEs or BLUPs}
#' }
#' Spatial plots can only be made if the data contains both row and column
#' information.
#'
#' @param x An object of class SSA.
#' @param ... Further graphical parameters (see \code{\link[lattice]{xyplot}}
#' for details).
#' @param trait a character string indicating the trait to plot. If
#' \code{trait = NULL} and only one trait is modelled this trait is plotted.
#' @param what A character string indicating whether the fitted model with
#' genotype as fixed or genotype as random factor should be plotted.
#' If \code{x} contains only one model this model is chosen automatically.
#' @param plotType character string indicating whether \code{base} plots or
#' \code{spatial} plots should be made.
#'
#' @seealso \code{\link{SSA}}
#'
#' @examples
#' ## Run a single trait analysis using SpATS.
#' myModel <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield")
#' ## Create base plots.
#' plot(myModel, what = "fixed", plotType = "base")
#' ## Create spatial plots.
#' plot(myModel, what = "fixed", plotType = "spatial")
#'
#' @export
plot.SSA <- function(x,
                     ...,
                     trait = NULL,
                     what = ifelse(is.null(x$mFix), "random", "fixed"),
                     plotType = c("base", "spatial")) {
  ## Checks.
  if (is.null(trait) && length(x$traits) > 1) {
    stop("No trait provided but multiple traits found in SSA x\n")
  }
  if (!is.null(trait) && (!is.character(trait) || length(trait) > 1 ||
                          !trait %in% colnames(x$TD))) {
    stop("Trait has to be a single character string defining a column in TD.\n")
  }
  what <- match.arg(what, choices = c("fixed", "random"))
  plotType <- match.arg(arg = plotType)
  ## If no trait is given as input extract it from the SSA object.
  if (is.null(trait)) {
    trait <- x$traits
  }
  ## Extract the model to plot from the SSA object.
  if (what == "fixed") {
    model <- x$mFix[[trait]]
  } else if (what == "random") {
    model <- x$mRand[[trait]]
  }
  ## Extract fitted and predicted values from model.
  fitted <- STExtract(x, what = ifelse(what == "fixed", "fitted", "rMeans"))[[trait]]
  pred <- STExtract(x, what = ifelse(what == "fixed", "BLUEs", "BLUPs"))[[trait]]
  ## Extract raw data and compute residuals.
  response <- x$TD[, trait]
  residuals <- response - fitted
  if (plotType == "base") {
    ## Setup frame for plots.
    trellisObj <- setNames(vector(mode = "list", length = 4),
                           c("histogram", "qq", "residFitted", "absResidFitted"))
    ## Plot histogram of residuals.
    trellisObj[["histogram"]] <- lattice::histogram(x = ~residuals,
                                                    xlab = "Residuals", ...)
    ## Plot Q-Q plot of residuals.
    trellisObj[["qq"]] <- lattice::qqmath(~residuals, xlab = "Normal quantiles",
                                          ylab = "Residuals", ...)
    ## Plot residuals vs fitted values
    trellisObj[["residFitted"]] <-
      lattice::xyplot(residuals ~ fitted,
                      panel = function(x, y, ...) {
                        lattice::panel.xyplot(x, y, ...,
                                              type = c("p", "g"))
                        lattice::panel.abline(h = 0)
                        lattice::panel.loess(x, y,
                                             col = "red", ...)
                      }, ylab = "Residuals",
                      xlab = "Fitted values", ...)
    ## Plot absolute residuals vs fitted values
    trellisObj[["absResidFitted"]] <-
      lattice::xyplot(abs(residuals) ~ fitted,
                      panel = function(x, y, ...) {
                        lattice::panel.xyplot(x, y, ...,
                                              type = c("p", "g"))
                        lattice::panel.loess(x, y,
                                             col = "red", ...)
                      }, ylab = "|Residuals|",
                      xlab = "Fitted values", ...)
    ## Save trellis options to reset when exiting function
    trellPar <- lattice::trellis.par.get()
    on.exit(lattice::trellis.par.set(trellPar))
    ## Change options for current plot.
    changeOpt <- list(add.text = list(cex = 0.75),
                      par.xlab.text = list(cex = 0.75),
                      par.ylab.text = list(cex = 0.75),
                      par.zlab.text = list(cex = 0.75),
                      axis.text = list(cex = 0.75),
                      plot.symbol = list(cex = 0.6))
    lattice::trellis.par.set(changeOpt)
    ## Fill frame with plots.
    print(trellisObj[["histogram"]], position = c(0, 0.5, 0.5, 1), more = TRUE)
    print(trellisObj[["qq"]], position = c(0.5, 0.5, 1, 1), more = TRUE)
    suppressWarnings(print(trellisObj[["residFitted"]], position = c(0, 0, 0.5, 0.5),
                           more = TRUE))
    suppressWarnings(print(trellisObj[["absResidFitted"]], position = c(0.5, 0, 1, 0.5)))
    invisible(trellisObj)
  } else if (plotType == "spatial") {
    if (x$engine == "SpATS") {
      plot(model, main = "")
    } else {
      ## Check whether data contains row/col information
      if (!all(c("rowCoordinates", "colCoordinates") %in% colnames(x$TD))) {
        stop(paste("Data in", substitute(x), "contains no spatial information.\n"))
      }
      ## Code taken from plot.SpATS and simplified.
      ## Set colors and legends.
      colors = topo.colors(100)
      mainLegends <- c("Raw data", "Fitted data", "Residuals",
                       ifelse(what == "fixed", "Genotypic BLUEs", "Genotypic BLUPs"),
                       "Histogram")
      ## Extract spatial coordinates from data.
      colCoord <- x$TD[, "colCoordinates"]
      rowCoord <- x$TD[, "rowCoordinates"]
      ## Order plotcols and rows and fill gaps if needed.
      plotCols <- seq(from = min(colCoord), to = max(colCoord),
                      by = min(diff(sort(unique(colCoord)))))
      plotRows <- seq(from = min(rowCoord), to = max(rowCoord),
                      by = min(diff(sort(unique(rowCoord)))))
      ## Compute range of values in response + fitted data so same scale
      ## can be used over plots.
      range <- range(c(response, fitted), na.rm = TRUE)
      ## Save plot options to reset when exiting function.
      op <- par(mfrow = c(2, 3), oma = c(2, 1, 3, 2),
                mar = c(2.5, 4, 2.5, 2.5), mgp = c(1.7, 0.5, 0))
      on.exit(par(op))
      ## Spatial plot of raw data.
      fields::image.plot(plotCols, plotRows, t(matrix(response, ncol = length(plotCols),
                                                      nrow = length(plotRows))),
                         main = mainLegends[1], col = colors,
                         xlab = "colCoordinates", ylab = "rowCoordinates",
                         zlim = range, graphics.reset = TRUE,
                         ...)
      ## Spatial plot of fitted data.
      fields::image.plot(plotCols, plotRows, t(matrix(fitted, ncol = length(plotCols),
                                                      nrow = length(plotRows))),
                         main = mainLegends[2], col = colors,
                         xlab = "colCoordinates", ylab = "rowCoordinates",
                         zlim = range, graphics.reset = TRUE,
                         ...)
      ## Spatial plot of residuals.
      fields::image.plot(plotCols, plotRows, t(matrix(residuals,
                                                      ncol = length(plotCols),
                                                      nrow = length(plotRows))),
                         main = mainLegends[3], col = colors,
                         xlab = "colCoordinates", ylab = "rowCoordinates",
                         graphics.reset = TRUE, ...)
      ## Spatial plot of BLUEs or BLUPs.
      fields::image.plot(plotCols, plotRows, t(matrix(pred, ncol = length(plotCols),
                                                      nrow = length(plotRows))),
                         main = mainLegends[4], col = colors,
                         xlab = "colCoordinates", ylab = "rowCoordinates",
                         graphics.reset = TRUE, ...)
      ## Histogram of BLUEs or BLUPs.
      suppressWarnings(hist(pred, main = mainLegends[5],
                            xlab = mainLegends[5], ...))
    }
  }
}

#' Report method for class SSA
#'
#' A pdf report will be created containing a summary of the results of the
#' fitted model. Simultaneously the same report will be created as a tex file.
#'
#' @inheritParams report.AMMI
#'
#' @param x An object of class SSA.
#' @param descending Should the trait be ordered in descending order? Set to
#' \code{FALSE} if low values of the trait indicate better performance.
#' @param what A character string indicating whether the model with genotype
#' fixed or genotype random should be reported. Can be omitted if only one
#' model has been fitted.
#'
#' @return A pdf and tex report.
#'
#' @examples
#' ## Fit model using lme4.
#' myModel1 <- STRunModel(TD = TDHeat05, design = "ibd", traits = "yield")
#' \dontrun{
#' ## Create a pdf report summarizing the results for the model with genotype
#' ## as fixed factor.
#' report(myModel1, outfile = "./testReports/reportModelLme4.pdf",
#'       what = "fixed")
#' ## Create a pdf report summarizing the results for the model with genotype
#' ## as random factor. Order the results in ascending order.
#' report(myModel1, outfile = "./testReports/reportModelLme4.pdf",
#'       what = "random", descending = FALSE)
#' }
#'
#' @export
report.SSA <- function(x,
                       ...,
                       descending = TRUE,
                       outfile = NULL,
                       what = if (is.null(x$mFix)) "random" else "fixed") {
  if (length(x$traits) > 1) {
    stop("Model contains models for multiple traits. Reporting can only be done
         for a single trait.\n")
  }
  what <- match.arg(what, choices = c("fixed", "random"))
  if (!is.null(x$mFix) && !is.null(x$mRand)) {
    warning("Model contains both a fitted model with fixed genotype and random
            genotype. Reporting can be done for only one. By default the model with
            genotype fixed is reported. Use option what for changing this.\n",
            call. = FALSE)
  }
  if (what == "fixed") {
    x$mRand <- NULL
  } else {
    x$mFix <- NULL
  }
  createReport(x = x, reportName = "modelReport.Rnw",
               outfile = outfile, ..., descending = descending)
}

#' Convert SSA to Cross
#'
#' Convert an SSA object to a cross object from class qtl. Genotypic information
#' should be available in a .csv file.\cr
#' The only way to create an object of class cross is by importing both the
#' phenotypic and the genotypic data from external files. Therefore the
#' phenotypic data, either the BLUEs or the BLUPs from the fitted model are
#' first written to a temporary file. The genotypic data has to be available in
#' a .csv file in the correct format as well, see \code{genoFile} for a
#' description of this format. These phenotypic and genotypic files are then
#' imported into a cross object using the read.cross function in the qtl
#' package.
#'
#' @param SSA An object of class \code{\link{SSA}}.
#' @param traits A character string containing the traits to be exported.
#' @param what A character string containing the statistics to be exported as
#' phenotype in the cross object. This can be either \code{BLUEs} or
#' \code{BLUPs}.
#' @param genoFile A character string indicating a filename containing
#' phenotypic data. The data should be in the format required by the
#' qtl package. The first column  should contain the individuals, starting
#' from row 4. The following columns contain markers with in the second and
#' third row the chromosome and position on the chromosome and in the
#' following rows the genotypes.
#' @param genotypes A character vector specifying the genotype codes
#' corresponding to AA, AB, BB, not BB and not AA.
#' @param ... Further arguments to be passed to the read.cross function.
#' See \code{\link[qtl]{read.cross}}.
#'
#' @seealso \code{\link[qtl]{read.cross}}
#'
#' @examples
#' ## Run model using SpATS.
#' myModel <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield",
#'                      what = "fixed")
#' ## Create cross object with BLUEs from myModel using genotypic information
#' ## from markers.csv in the package.
#' cross <- SSAtoCross(myModel, genoFile = system.file("extdata", "markers.csv",
#'                                                    package = "RAP"))
#'
#' @export
SSAtoCross <- function(SSA,
                       traits = SSA$traits,
                       what = c("BLUEs", "BLUPs"),
                       genoFile,
                       genotypes = c("A", "H", "B", "D", "C"),
                       ...) {
  ## Checks
  if (!inherits(SSA, "SSA")) {
    stop("SSA is not a valid object of class SSA.\n")
  }
  what <- match.arg(what)
  if (!is.character(genoFile) || length(genoFile) > 1 || !file.exists(genoFile)) {
    stop("genoFile is not a valid filename.\n")
  }
  ## Extract predictions from the model.
  pred <- STExtract(SSA, traits = traits, what = what)
  ## Rename first column to match first column in genoFile.
  colnames(pred)[1] <- colnames(utils::read.csv(genoFile, nrow = 1))[1]
  ## Write predictions to temporary file.
  tmp <- tempfile()
  utils::write.csv(pred, file = tmp, row.names = FALSE)
  ## Read cross from temporary file and supplied genoFile.
  cross <- qtl::read.cross(format = "csvs",
                           phefile = tmp,
                           genfile = genoFile,
                           genotypes = genotypes,
                           ...)
  unlink(tmp)
  return(cross)
}



