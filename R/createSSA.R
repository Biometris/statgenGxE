#' S3 class SSA
#'
#' Function for creating objects of S3 class Single Site Analysis (SSA).\cr
#' \code{\link{summary}}, \code{\link{plot}} and \code{\link{report}}
#' methods are available.
#'
#' @param models A list of trials with for each trial the following elements
#' \itemize{
#' \item{mRand}{A list of models with fitted with genotype as random effect.}
#' \item{mFix}{A list of models fitted with genotype as fixed effect.}
#' \item{TD}{An object of class \code{\link{TD}} containing the data on which
#' \code{mRand} and \code{mFix} are based.}
#' \item{traits}{A character vector indicating the traits for which the analysis
#' is done.}
#' \item{design}{A character string containing the design of the trial.
#' (see \code{\link{STRunModel}} for the possible designs).}
#' \item{spatial}{A character string indicating the spatial part of the model.
#' \code{FALSE} if no spatial design has been used.}
#' \item{engine}{A character string containing the engine used for the
#' analysis.}
#' \item{predicted}{A character string indicating the variable that has been
#' predicted.}
#' }
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{summary.SSA}}, \code{\link{plot.SSA}},
#' \code{\link{report.SSA}}
#'
#' @name SSA
NULL

#' @rdname SSA
#' @keywords internal
createSSA <- function(models) {
  SSA <- structure(models,
                   class = c("SSA", "list"),
                   timestamp = Sys.time())
  return(SSA)
}

#' Summarizing objects of class \code{SSA}
#'
#' \code{summary} method for class \code{SSA}.
#'
#' @param object An object of class \code{SSA}.
#' @param trial A character string indicating the trial to summarize. If
#' \code{trial = NULL} and only one trial is modelled this trial is summarized.
#' @param trait A character string indicating the trait to summarize. If
#' \code{trait = NULL} and only one trait is modelled this trait is summarized.
#' @param nBest An integer indicating the number of the best genotypes (sorted
#' by either BLUEs or BLUPs) to print. If \code{NA} all genotypes will be
#' printed.
#' @param sortBy A character string specifying how the genotypes will be sorted.
#' Either \code{"BLUEs"}, \code{"BLUPs"} or \code{NA} (i.e. no sorting).
#' @param naLast Should missing values in the data be put last when sorting?
#' @param decreasing Should the sort order be decreasing?
#' @param ... Further arguments - not used.
#'
#' @examples
#' ## Run a single trait analysis using SpATS.
#' myModel <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield")
#' ## Print a summary of the fitted model.
#' summary(myModel)
#'
#' @export
summary.SSA <- function(object,
                        trial = NULL,
                        trait = NULL,
                        nBest = 20,
                        sortBy = NULL,
                        naLast = TRUE,
                        decreasing = TRUE,
                        ...) {
  ## Checks.
  if (is.null(trial) && length(object) > 1) {
    stop("No trial provided but multiple trials found in SSA object.\n")
  }
  if (!is.null(trial) && (!is.character(trial) || length(trial) > 1 ||
                          !trial %in% names(object))) {
    stop("Trial has to be a single character string defining a trial in SSA.\n")
  }
  if (is.null(trial)) {
    trial <- names(object)
  }
  if (is.null(trait) && length(object[[trial]]$traits) > 1) {
    stop("No trait provided but multiple traits found.\n")
  }
  if (!is.null(trait) && (!is.character(trait) || length(trait) > 1 ||
                          !trait %in% colnames(object[[trial]]$TD))) {
    stop("Trait has to be a single character string defining a column in TD.\n")
  }
  if (is.null(sortBy)) {
    sortBy <- ifelse(!is.null(object[[trial]]$mFix), "BLUEs", "BLUPs")
  } else {
    sortBy <- match.arg(sortBy)
  }
  ## get summary stats for raw data
  TD <- object[[trial]]$TD
  if (is.null(trait)) {
    trait <- object[[trial]]$traits
  }
  stats <- summary.TD(object = TD, traits = trait)
  ## get predicted means (BLUEs + BLUPs).
  extr <- STExtract(object, trials = trial)[[trial]]
  ## Merge results using a loop to avoid warnings over suffixes caused by
  ## merge when using using Reduce.
  joinList <- Filter(f = Negate(f = is.null),
                     x = extr[c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs")])
  meanTab <- joinList[[1]]
  for (i in 2:length(joinList)) {
    meanTab <- merge(meanTab, joinList[[i]], all = TRUE, by = "genotype",
                     suffixes = c(i, i + 1))
  }
  ## Move genotype to rowname for proper printing with printCoefMat.
  rownames(meanTab) <- meanTab$genotype
  meanTab <- meanTab[colnames(meanTab) != "genotype"]
  ## Set colnames. Because of duplicate colname SE no selection on columns can
  ## be done anymore after this.
  colnames(meanTab) <- c(if (!is.null(extr$BLUEs)) c("BLUEs", "SE"),
                         if (!is.null(extr$BLUPs)) c("BLUPs", "SE"))
  meansTxt <- paste(c(if (!is.null(extr$BLUEs)) "BLUEs",
                      if (!is.null(extr$BLUPs)) "BLUPs"), collapse = " & ")
  attr(x = meanTab, which = "title") <- meansTxt
  if (!is.na(sortBy)) {
    ## Sort by sortBy with options from input params.
    oList <- order(meanTab[[sortBy]], na.last = naLast, decreasing = decreasing)
    meanTab <- meanTab[oList, ]
  }
  if (!is.na(nBest)) {
    ## Set nBest to number of rows in meanTab to prevent printing of NA rows.
    nBest <- min(nrow(meanTab), nBest)
    ## Extract the n best genotypes.
    meanTab <- meanTab[1:nBest, ]
    attr(x = meanTab, which = "nBest") <- nBest
  }
  ## Extract selected spatial model when applicable.
  if (object[[trial]]$engine == "asreml" &&
      is.character(object[[trial]]$spatial[[trait]])) {
    selSpatMod <- object[[trial]]$spatial[[trait]]
  } else {
    selSpatMod <- NULL
  }
  return(structure(list(selSpatMod = selSpatMod, stats = stats,
                        meanTab = meanTab, heritability = extr$heritability,
                        sed = data.frame("s.e.d" = extr$sed),
                        lsd = data.frame("l.s.d." = extr$lsd)),
                   class = c("summary.SSA")))
}

#' Printing summazed objects of class SSA
#'
#' \code{print} method for object of class summary.SSA created by summarizing
#' objects of class SSA.
#'
#' @param x An object of class \code{summary.SSA}
#' @param digits An integer indicating the number of significant digits for
#' printing.
#' @param ... Further arguments passed to \code{\link[stats]{printCoefmat}}.
#'
#' @export
print.summary.SSA <- function(x,
                              digits = max(getOption("digits") - 2, 3),
                              ...) {
  if (!is.null(x$selSpatMod)) {
    cat("Selected spatial model: ", x$selSpatMod, "\n\n")
  }
  cat("Summary statistics",
      "\n==================\n")
  ## Print stats using printCoefMat for a nicer layout.
  printCoefmat(x$stats, digits = digits, ...)
  if (!is.null(x$heritability)) {
    cat("\nEstimated heritability",
        "\n======================\n")
    cat("\nHeritability:", x$heritability, "\n")
  }
  cat(paste0("\nPredicted means (", attr(x = x$meanTab, which = "title"), ")"),
      "\n===============================\n")
  if (!is.null(attr(x = x$meanTab, which = "nBest"))) {
    cat("Best", attr(x = x$meanTab, which = "nBest"), "genotypes\n")
  } else {
    cat("\n")
  }
  printCoefmat(x$meanTab, digits = digits, ...)
  if (nrow(x$lsd) > 0) {
    cat("\nStandard Error of Difference (genotype modelled as fixed effect)",
        "\n================================================================\n")
    printCoefmat(x$sed, digits = digits, ...)
  }
  if (nrow(x$lsd) > 0) {
    cat("\nLeast Significant Difference (genotype modelled as fixed effect)",
        "\n================================================================\n")
    printCoefmat(x$lsd, digits = digits, ...)
  }
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
#' @param ... Further graphical parameters.
#' @param trial A character string indicating the trial to plot. If
#' \code{trial = NULL} and only one trial is modelled this trial is plotted.
#' @param trait A character string indicating the trait to plot. If
#' \code{trait = NULL} and only one trait is modelled this trait is plotted.
#' @param what A character string indicating whether the fitted model with
#' genotype as fixed or genotype as random factor should be plotted.
#' If \code{x} contains only one model this model is chosen automatically.
#' @param plotType A Character string indicating whether \code{base} plots or
#' \code{spatial} plots should be made.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @return A list containing ggplot object for the selected plots.
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
                     trial = NULL,
                     trait = NULL,
                     what = NULL,
                     plotType = c("base", "spatial"),
                     output = TRUE) {
  ## Checks.
  if (is.null(trial) && length(x) > 1) {
    stop("No trial provided but multiple trials found in SSA object.\n")
  }
  if (!is.null(trial) && (!is.character(trial) || length(trial) > 1 ||
                          !trial %in% names(x))) {
    stop("Trial has to be a single character string defining a trial in SSA.\n")
  }
  if (is.null(trial)) {
    trial <- names(x)
  }
  if (is.null(trait) && length(x[[trial]]$traits) > 1) {
    stop("No trait provided but multiple traits found.\n")
  }
  if (!is.null(trait) && (!is.character(trait) || length(trait) > 1 ||
                          !trait %in% colnames(x[[trial]]$TD[[trial]]))) {
    stop("Trait has to be a single character string defining a column in TD.\n")
  }
  if (is.null(what)) {
    what <- ifelse(is.null(x[[trial]]$mFix), "random", "fixed")
  } else {
    what <- match.arg(arg = what, choices = c("fixed", "random"))
  }
  plotType <- match.arg(arg = plotType)
  dotArgs <- list(...)
  ## Check whether data contains row/col information.
  spatCols <- c("colCoord", "rowCoord")
  if (plotType == "spatial" && !all(spatCols %in%
                                    colnames(x[[trial]]$TD[[trial]]))) {
    stop(paste("Data in", substitute(x), "contains no spatial information.\n"))
  }
  ## If no trait is given as input extract it from the SSA object.
  if (is.null(trait)) {
    trait <- x[[trial]]$traits
  }
  ## Extract the model to plot from the SSA object.
  if (what == "fixed") {
    model <- x[[trial]]$mFix[[trait]]
  } else if (what == "random") {
    model <- x[[trial]]$mRand[[trait]]
  }
  if (is.null(model)) {
    stop(paste("No model with genotype", what, "in SSA object.\n"))
  }
  predicted <- x[[trial]]$predicted
  ## Extract fitted and predicted values from model.
  fitted <- STExtract(x, trials = trial, traits = trait,
                      what = ifelse(what == "fixed", "fitted", "rMeans"),
                      keep = if (plotType == "spatial") spatCols else
                        NULL)[[trial]][[ifelse(what == "fixed",
                                               "fitted", "rMeans")]]
  pred <- STExtract(x, trials = trial,
                    what = ifelse(what == "fixed", "BLUEs", "BLUPs"))[[trial]][[ifelse(what == "fixed",
                                                                                       "BLUEs", "BLUPs")]][c(predicted,
                                                                                                             trait)]
  ## Extract raw data and compute residuals.
  response <- x[[trial]]$TD[[trial]][, c(predicted, trait,
                                         if (plotType == "spatial") spatCols)]
  ## Create plot data by merging extracted data together and renaming some
  ## columns.
  plotDat <- merge(response,
                   fitted, by = c(predicted,
                                  if (plotType == "spatial") spatCols))
  plotDat <- merge(plotDat, pred, by = predicted)
  plotDat$response <- plotDat[[paste0(trait, ".x")]]
  plotDat$fitted <- plotDat[[paste0(trait, ".y")]]
  plotDat$pred <- plotDat[[trait]]
  plotDat$pred[is.na(plotDat$fitted)] <- NA
  plotDat$residuals <- plotDat$response - plotDat$fitted
  ## Create empty list for storing plots
  plots <- vector(mode = "list")
  ## Create main plot title.
  plotTitle <- ifelse(!is.null(dotArgs$title), dotArgs$title,
                      paste("Trial:", trial, "Trait:", trait))
  if (plotType == "base") {
    plotDat <- ggplot2::remove_missing(plotDat, na.rm = TRUE)
    ## Plot histogram of residuals.
    plots$p1 <- ggplot2::ggplot(data = plotDat) +
      ggplot2::geom_histogram(ggplot2::aes(x = residuals,
                                           y = (..count..)/sum(..count..)),
                              fill = "cyan", col = "black", bins = 10,
                              boundary = 0) +
      ggplot2::scale_y_continuous(labels = function(x) {paste0(100 * x, "%")}) +
      ggplot2::labs(y = "Percent of Total", x = "Residuals")
    ## Plot Q-Q plot of residuals.
    plots$p2 <- ggplot2::ggplot(data = plotDat,
                                ggplot2::aes_string(sample = "residuals")) +
      ggplot2::stat_qq(col = "blue") +
      ggplot2::labs(y = "Residuals", x = "Normal quantiles")
    ## Plot residuals vs fitted values.
    plots$p3 <- ggplot2::ggplot(data = plotDat,
                                ggplot2::aes_string(x = "fitted",
                                                    y = "residuals")) +
      ggplot2::geom_point(col = "blue", shape = 1) +
      ggplot2::geom_smooth(method = "loess", col = "red") +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::labs(y = "Residuals", x = "Fitted values")
    ## Plot absolute value of residuals vs fitted values.
    plots$p4 <- ggplot2::ggplot(data = plotDat,
                                ggplot2::aes_string(x = "fitted",
                                                    y = "abs(residuals)")) +
      ggplot2::geom_point(col = "blue", shape = 1) +
      ggplot2::geom_smooth(method = "loess", col = "red") +
      ggplot2::labs(y = "|Residuals|", x = "Fitted values")
    if (output) {
      ## do.call is needed since grid.arrange doesn't accept lists as input.
      do.call(gridExtra::grid.arrange,
              args = c(plots, list(ncol = 2, top = plotTitle)))
    }
  } else if (plotType == "spatial") {
    if (x[[trial]]$engine == "SpATS") {
      ## Execute this part first since it needs plotData without missings
      ## removed.
      ## Code mimickes code from SpATS packages but is adapted to create a
      ## data.frame useable by ggplot.
      plotDat <- plotDat[order(plotDat$colCoord, plotDat$rowCoord), ]
      nCol <- length(unique(plotDat$colCoord))
      nRow <- length(unique(plotDat$rowCoord))
      p1 <- 100 %/% nCol + 1
      p2 <- 100 %/% nRow + 1
      ## Get spatial trend from SpATS object.
      spatTr <- SpATS::obtain.spatialtrend(model,
                                           grid = c(nCol * p1, nRow * p2))
      ## spatial trend contains values for all data points, so NA in original
      ## data need to be removed. The kronecker multiplication is needed to
      ## convert the normal row col pattern to the smaller grid extending the
      ## missing values.
      spatTrDat <- kronecker(matrix(data = ifelse(is.na(plotDat$response),
                                                  NA, 1),
                                    ncol = nCol, nrow = nRow),
                             matrix(data = 1, ncol = p1, nrow = p2)) *
        spatTr$fit
      ## Melt to get the data in ggplot shape. Rows and columns in the
      ## spatial trend coming from SpATS are swapped so therefore use t()
      plotDatSpat <- reshape2::melt(t(spatTrDat),
                                    varnames = c("colCoord", "rowCoord"))
      ## Add true values for columns and rows for plotting.
      plotDatSpat$colCoord <- spatTr$col.p
      plotDatSpat$rowCoord <- rep(x = spatTr$row.p, each = p1 * nCol)
      ## Remove missings from data.
      plotDatSpat <- ggplot2::remove_missing(plotDatSpat, na.rm = TRUE)
    }
    plotDat <- ggplot2::remove_missing(plotDat, na.rm = TRUE)
    ## Code taken from plot.SpATS and simplified.
    ## Set colors and legends.
    colors = topo.colors(100)
    legends <- c("Raw data", "Fitted data", "Residuals",
                 "Fitted Spatial Trend",
                 ifelse(what == "fixed", "Genotypic BLUEs",
                        "Genotypic BLUPs"), "Histogram")
    ## Compute range of values in response + fitted data so same scale
    ## can be used over plots.
    zlim <- range(c(plotDat$response, plotDat$fitted), na.rm = TRUE)
    plots$p1 <- fieldPlot(plotDat = plotDat, fillVar = "response",
                          title = legends[1], colors = colors, zlim = zlim)
    plots$p2 <- fieldPlot(plotDat = plotDat, fillVar = "fitted",
                          title = legends[2], colors = colors, zlim = zlim)
    plots$p3 <- fieldPlot(plotDat = plotDat, fillVar = "residuals",
                          title = legends[3], colors = colors)
    if (x[[trial]]$engine == "SpATS") {
      ## Get tickmarks from first plot to be used as ticks.
      ## Spatial plot tends to use different tickmarks by default.
      xTicks <-
        ggplot2::ggplot_build(plots[[1]])$layout$panel_params[[1]]$x.major_source
      plots$p4 <- fieldPlot(plotDat = plotDatSpat, fillVar = "value",
                            title = legends[4], colors = colors,
                            xTicks = xTicks)
    }
    plots$p5 <- fieldPlot(plotDat = plotDat, fillVar = "pred",
                          title = legends[5], colors = colors)
    plots$p6 <- ggplot2::ggplot(data = plotDat) +
      ggplot2::geom_histogram(ggplot2::aes(x = residuals),
                              fill = "white", col = "black", bins = 10,
                              boundary = 0) +
      ## Remove empty space between ticks and actual plot.
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ## No background. Center and resize title. Resize axis labels.
      ggplot2::theme(panel.background = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0.5,
                                                        size = 10),
                     axis.title = ggplot2::element_text(size = 9)) +
      ggplot2::labs(y = "Frequency", x = legends[5]) +
      ggplot2::ggtitle(legends[6])
    if (output) {
      ## do.call is needed since grid.arrange doesn't accept lists as input.
      do.call(gridExtra::grid.arrange,
              args = c(Filter(f = Negate(f = is.null), x = plots),
                       list(ncol = 3, top = plotTitle)))
    }
  }
  invisible(plots)
}

## Helper function for creating field plots.
fieldPlot <- function(plotDat,
                      fillVar,
                      title,
                      colors,
                      zlim = range(plotDat[fillVar]),
                      xTicks = ggplot2::waiver(),
                      ...) {
  p <- ggplot2::ggplot(data = plotDat,
                       ggplot2::aes_string(x = "colCoord", y = "rowCoord",
                                           fill = fillVar)) +
    ggplot2::geom_raster() +
    ## Remove empty space between ticks and actual plot.
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = xTicks) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ## Adjust plot colors.
    ggplot2::scale_fill_gradientn(limits = zlim, colors = colors) +
    ## No background. Center and resize title. Resize axis labels.
    ## Remove legend title and resize legend entries.
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
                   axis.title = ggplot2::element_text(size = 9),
                   legend.title = ggplot2::element_blank(),
                   legend.text =
                     ggplot2::element_text(size = 8,
                                           margin = ggplot2::margin(l = 5))) +
    ggplot2::ggtitle(title)
  return(p)
}

#' Report method for class SSA
#'
#' A pdf report will be created containing a summary of the results of the
#' fitted model. Simultaneously the same report will be created as a tex file.
#'
#' @inheritParams report.AMMI
#'
#' @param x An object of class SSA.
#' @param trial A character string indicating the trial to be reported. If
#' \code{NULL} and \code{SSA} contains only one trial that trial is reported.
#' @param trait A character string indicating the trait to be reported. If
#' \code{NULL} and \code{SSA} contains only one trait that trait is reported.
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
#'        what = "fixed")
#' ## Create a pdf report summarizing the results for the model with genotype
#' ## as random factor. Order the results in ascending order.
#' report(myModel1, outfile = "./testReports/reportModelLme4.pdf",
#'        what = "random", descending = FALSE)
#' }
#'
#' @export
report.SSA <- function(x,
                       ...,
                       trial = NULL,
                       trait = NULL,
                       descending = TRUE,
                       outfile = NULL,
                       what = c("fixed", "random")) {
  if (is.null(trial) && length(x) > 1) {
    stop("No trial provided but multiple trials found in SSA object.\n")
  }
  if (!is.null(trial) && (!is.character(trial) || length(trial) > 1 ||
                          !trial %in% names(x))) {
    stop("Trial has to be a single character string defining a trial in SSA.\n")
  }
  if (is.null(trial)) {
    trial <- names(x)
  }
  if (is.null(trait) && length(x[[trial]]$traits) > 1) {
    stop("No trait provided but multiple traits found.\n")
  }
  if (!is.null(trait) && (!is.character(trait) || length(trait) > 1 ||
                          !trait %in% colnames(x[[trial]]$TD[[trial]]))) {
    stop("Trait has to be a single character string defining a column in TD.\n")
  }
  ## If no trait is given as input extract it from the SSA object.
  if (is.null(trait)) {
    trait <- x[[trial]]$traits
  }
  what <- match.arg(what)
  if (is.null(x[[trial]]$mFix)) {
    what <- "random"
  }
  if (is.null(x[[trial]]$mRand)) {
    what <- "fixed"
  }
  if (is.null(what) && !is.null(x[[trial]]$mFix) &&
      !is.null(x[[trial]]$mRand)) {
    warning("Model contains both a fitted model with fixed genotype and random
            genotype. Reporting can be done for only one. By default the model with
            genotype fixed is reported. Use option 'what' for changing this.\n",
            call. = FALSE)
  }
  if (what == "fixed") {
    x[[trial]]$mRand <- NULL
  } else {
    x[[trial]]$mFix <- NULL
  }
  createReport(x = x, reportName = "modelReport.Rnw",
               outfile = outfile, ..., trial = trial, trait = trait,
               descending = descending)
}

#' Convert SSA to Cross
#'
#' Convert an SSA object to a cross object from package qtl. Genotypic
#' information should be available in a .csv file.\cr
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
#' @param trial A character string indicating the trial to be exported. If
#' \code{NULL} and \code{SSA} contains only one trial that trial is exported.
#' @param traits A character string containing the traits to be exported. If
#' \code{NULL} all traits for the selected trial are exported.
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
#'                       what = "fixed")
#' ## Create cross object with BLUEs from myModel using genotypic information
#' ## from markers.csv in the package.
#' cross <- SSAtoCross(myModel, genoFile = system.file("extdata", "markers.csv",
#'                                                     package = "RAP"))
#'
#' @export
SSAtoCross <- function(SSA,
                       trial = NULL,
                       traits = NULL,
                       what = c("BLUEs", "BLUPs"),
                       genoFile,
                       genotypes = c("A", "H", "B", "D", "C"),
                       ...) {
  ## Checks
  if (!inherits(SSA, "SSA")) {
    stop("SSA is not a valid object of class SSA.\n")
  }
  if (is.null(trial) && length(SSA) > 1) {
    stop("No trial provided but multiple trials found in SSA object.\n")
  }
  if (!is.null(trial) && (!is.character(trial) || length(trial) > 1 ||
                          !trial %in% names(SSA))) {
    stop("Trial has to be a single character string defining a trial in SSA.\n")
  }
  if (is.null(trial)) {
    trial <- names(SSA)
  }
  if (!is.null(traits) && (!is.character(traits) ||
                           !all(traits %in% colnames(SSA[[trial]]$TD)))) {
    stop("Trait has to be a character vector defining columns in TD.\n")
  }
  if (is.null(traits)) {
    traits <- SSA[[trial]]$traits
  }
  what <- match.arg(what)
  if (!is.character(genoFile) || length(genoFile) > 1 || !file.exists(genoFile)) {
    stop("genoFile is not a valid filename.\n")
  }
  ## Extract predictions from the model.
  pred <- STExtract(SSA, traits = traits, what = what)[[trial]][[what]]
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

#' Convert SSA to TD
#'
#' Convert an SSA object to a TD object.\cr
#' To be able to use the output of a single site analysis in Genotype by
#' Environment (GxE) analysis the output first needs to be converted bakc to
#' an TD object. This function does exactly that. It extracts BLUEs, BLUPs and
#' their standard errors from the SSA object and creates a new TD object using
#' these. Also a column wt may also be added. Weights are then calculated as
#' 1/SE BLUEs.
#'
#' Trial information for the trials in the SSA object will be copied from the
#' original TD object on which the modeling was done.
#'
#' @param SSA An object of class \code{\link{SSA}}.
#' @param what A character string containing the statistics to be included as
#' traits in the TD object. Multiple statistics can be included in which case
#' they will appear as \code{statistic_trait} in the output
#' @param traits A character string containing the traits to be included in the
#' TD object. If \code{NULL} all traits are exported.
#' @param keep Columns from the TD object used as input for the SSA model to
#' be copied to the output. see \code{\link{STExtract}} for possible columns to
#' copy. If if it is available in TD the column \code{trial} will always be
#' copied.
#' @param addWt Should a column wt be added to the output? If \code{TRUE}
#' weight is calculated as 1/SE BLUEs. If multiple traits are included in the
#' output multiple weight columns will be added, 1 for each trait. These will
#' be named \code{wt_trait}.
#'
#' @examples
#' ## Run model using SpATS.
#' myModel <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield",
#'                       what = "fixed")
#' ## Create TD object from the fitted model.
#' myTD <- SSAtoTD(myModel)
#'
#' @export
SSAtoTD <- function(SSA,
                    what = c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs"),
                    traits = NULL,
                    keep = NULL,
                    addWt = FALSE) {
  ## Checks
  if (!inherits(SSA, "SSA")) {
    stop("SSA is not a valid object of class SSA.\n")
  }
  if (!is.null(traits) && (!is.character(traits) ||
                           !all(traits %in% colnames(SSA[[1]]$TD)))) {
    stop("Trait has to be a character vector defining columns in TD.\n")
  }
  what <- match.arg(what, several.ok = TRUE)
  if (any(c("BLUEs", "seBLUEs") %in% what) && is.null(SSA[[1]]$mFix)) {
    warning(paste("BLUEs and seBLUEs can only extracted if a model with",
                  "genotype fixed is fitted\nRemoving them from what"),
            call. = FALSE)
    what <- what[!what %in% c("BLUEs", "seBLUEs")]
  }
  if (any(c("BLUPs", "seBLUPs") %in% what) && is.null(SSA[[1]]$mRand)) {
    warning(paste("BLUPs and seBLUPs can only extracted if a model with",
                  "genotype random is fitted\nRemoving them from what"),
            call. = FALSE)
    what <- what[!what %in% c("BLUPs", "seBLUPs")]
  }
  if (is.null(what)) {
    stop("No statistics left to extract.")
  }
  if (addWt && is.null(SSA[[1]]$mFix)) {
    warning(paste("Weights can only be added if a model with genotype fixed is",
                  "fitted.\naddWt set to FALSE"),
            call. = FALSE)
    addWt <- FALSE
  }
  if (addWt && !"seBLUEs" %in% what) {
    warning(paste("Weights can only be added together with seBLUEs.\n",
                  "seBLUEs added to what"), call. = FALSE)
    what <- c(what, "seBLUEs")
  }
  if (is.null(traits)) {
    traits <- SSA[[1]]$traits
  }
  if (!"trial" %in% keep && hasName(x = SSA[[1]]$TD[[1]], name = "trial")) {
    keep <- c(keep, "trial")
  }
  ## Extract predictions from the model.
  pred <- STExtract(SSA, traits = traits,
                    what = what,
                    keep = keep)
  ## Create a list of dataframes with all statistics per trial
  predTrTot <- lapply(X = pred, FUN = function(trial) {
    if (length(what) + addWt > 1) {
      ## Rename columns if more than one column per trait will appear in the
      ## output. Add the name of the statistic as prefix to the traits.
      for (ext in names(trial)) {
        colNames <- colnames(trial[[ext]])
        colnames(trial[[ext]])[colNames %in% traits] <-
          paste0(ext, "_", colNames[colNames %in% traits])
      }
    }
    ## Merge all statistics togethter. Because of the renaming above the is
    ## never a problem with duplicate column and merging is done on all other
    ## columns than the traits.
    predTr <- Reduce(f = merge, x = trial)
    if (addWt && "seBLUEs" %in% what) {
      ## Add a wt column.
      nTr <- length(traits)
      for (trait in traits) {
        wtName <- ifelse(nTr == 1, "wt", paste0("wt_", trait))
        predTr[[wtName]] <- 1 / predTr[[paste0("seBLUEs_", trait)]]
      }
    }
    return(predTr)
  })
  ## Rbind all data together and create a new TD data set.
  predTot <- Reduce(f = rbind, x = predTrTot)
  predTD <- createTD(data = predTot)
  ## Copy meta data from the original TD to the new TD.
  predTD <- setMeta(TD = predTD, meta = getMeta(SSA[[1]]$TD))
  return(predTD)
}






