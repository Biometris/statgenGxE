#' S3 class FW
#'
#' Function for creating objects of S3 class FW (Finlay Wilkinson).\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param estimates a data.frame containing the estimated values
#' @param anova a data.frame containing anova scores of the FW analysis
#' @param envEffs a data.frame containing the environmental effects
#' @param data the data.frame on which the analysis was performed
#' @param fittedGeno the fitted values for the genotypes
#' @param trait a character value indicating the analysed trait
#' @param nGeno a numerical value containing the number of genotypes in the analysis
#' @param nEnv a numerical value containing the number of environments in the analysis
#' @param tol a numerical value containing the tolerance used during the analysis
#' @param iter a numberical value containing the number of iterations for the
#' analysis to converge
#' @param x an \code{R} object
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.FW}}, \code{\link{report.FW}}
#'
#' @name FW
NULL

#' @rdname FW
#' @export
createFW <- function(estimates,
                     anova,
                     envEffs,
                     trait,
                     nGeno,
                     nEnv,
                     data,
                     fittedGeno,
                     tol,
                     iter) {
  FW <- structure(list(estimates = estimates,
                       anova = anova,
                       envEffs = envEffs,
                       data = data,
                       fittedGeno = fittedGeno,
                       trait = trait,
                       nGeno = nGeno,
                       nEnv = nEnv,
                       tol = tol,
                       iter = iter),
                  class = "FW")
  attr(FW, which = "timestamp") <- Sys.time()
  return(FW)
}

#' @rdname FW
#' @export
is.FW <- function(x) {
  inherits(x, "FW")
}

#' @export
print.FW <- function(x, ...) {
  cat("Environmental effects",
      "\n===================\n")
  print(x$envEffs)
  cat("\nAnova",
      "\n=====\n")
  printCoefmat(x$anova)
  cat("\nEstimates",
      "\n=========\n")
  print(x$estimates, ..., row.names = FALSE)
}

#' @export
summary.FW <- function(object, ...) {
  print(object, ...)
}

#' Plot Function for Class FW
#'
#' Three types of plot can be made. A scatter plot for genotypic mean,
#' mse and sensitivity, a line plot with fitted lines for each genotype and
#' a trellis plot with individual slopes per genotype (for max 64 genotypes).
#'
#' @param x An object of class FW
#' @param ... Other graphical parameters passed on to actual plot function.
#' @param plotType A character string indicating which plot should be made.
#' Possible values are "scatter", "line"  and "trellis" for creating a scatter
#' plot of sensitivities, a plot of fitted lines for each genotype and a trellis
#' plot of the individual genotype slopes respectively.
#' @param sorted A character string specifying whether the results are to be sorted
#' in an increasing (or decreasing) order of sensitivities.

#' @return A plot depending on \code{plotType}
#'
#' @examples
#' # Run Finlay-Wilkinson analysis.
#' geFW <- gxeFw(TD = TDMaize, trait = "yld")
#' # Create scatter plot.
#' plot(geFW, plotType = "scatter")
#'
#' @import graphics grDevices
#' @export
plot.FW <- function(x,
                    ...,
                    plotType = c("scatter", "line", "trellis"),
                    sorted = c("ascending", "descending", "none")) {
  plotType <- match.arg(plotType, several.ok = TRUE)
  sorted <- match.arg(sorted)
  dotArgs <- list(...)
  envEffs <- x$envEffs$Effect
  if ("scatter" %in% plotType) {
    selCols = c(1, if (!all(is.na(x$estimates$mse))) 2, 3)
    scatterData <- setNames(x$estimates[, c("genMean", "mse", "sens")[selCols]],
                            c("Mean", "m.s.deviation", "Sensitivity")[selCols])
    ## Set arguments for plot.
    plotArgs <- list(x = scatterData, upper.panel = NULL,
                     main = paste0("Finlay & Wilkinson analysis for ", x$trait))
    ## Add and overwrite args with custom args from ...
    fixedArgs <- c("x")
    plotArgs <- modifyList(plotArgs, dotArgs[-which(names(dotArgs) %in% fixedArgs)])
    do.call(ifelse(!all(is.na(x$estimates$mse)), pairs, plot), args = plotArgs)
  } else if ("line" %in% plotType) {
    fVal <- tapply(X = x$fittedGeno, INDEX = x$data[, c("env", "genotype")],
                   FUN = mean, na.rm = TRUE)
    if (sorted == "none") {
      orderEnv <- 1:length(envEffs)
    } else {
      orderEnv <- order(envEffs, decreasing = (sorted == "descending"))
    }
    ## Set arguments for plot.
    plotArgs <- list(x = envEffs[orderEnv],
                     y = fVal[orderEnv, ], type = "l",
                     main = paste0("Finlay & Wilkinson analysis for ", x$trait),
                     ylab = x$trait, xlab = "Environment",
                     xlim = range(envEffs, na.rm = TRUE),
                     ylim = range(x$fittedGeno, na.rm = TRUE),
                     xaxt = "n", col = 1:ncol(fVal))
    ## Add and overwrite args with custom args from ...
    fixedArgs <- c("x", "y", "xlim", "ylim", "xaxt")
    plotArgs <- modifyList(plotArgs, dotArgs[-which(names(dotArgs) %in% fixedArgs)])
    do.call(matplot, args = plotArgs)
    ## Add environments as ticks on axis.
    axis(side = 1, at = envEffs, labels = levels(x$envEffs$Environment),
         las = 2, cex.axis = .75)
  } else if ("trellis" %in% plotType) {
    trellisData <- data.frame(genotype = x$data$genotype,
                              trait = x$data[[x$trait]],
                              fitted = x$fittedGen,
                              xEff = rep(envEffs, x$nGeno))
    if (x$nGeno > 64) {
      ## Select first 64 genotypes for plotting.
      first64 <- x$data$genotype %in% levels(x$estimates$genotype)[1:64]
      trellisData <- droplevels(trellisData[first64, ])
    }
    ## Define panelfunction for xy plot.
    panelFunc <- function(x, y, subscripts) {
      lattice::panel.xyplot(x, y)
      lattice::panel.lines(trellisData$xEff[subscripts],
                           trellisData$fitted[subscripts])
    }
    ## Set arguments for plot.
    plotArgs <- list(x = trait + fitted ~ xEff | genotype, data = trellisData,
      panel = panelFunc, as.table = TRUE, subscripts = TRUE,
      xlab = "Environment", ylab = x$trait,
      main = paste0("Finlay & Wilkinson analysis for ", x$trait))
    fixedArgs <- c("x", "data", "panel")
    plotArgs <- modifyList(plotArgs, dotArgs[-which(names(dotArgs) %in% fixedArgs)])
    ## Add and overwrite args with custom args from ...
    do.call(lattice::xyplot, args = plotArgs)
  }
}

#' Report method for class FW
#'
#' A pdf report will be created containing a summary of FW analysis.
#' Simultaneously the same report will be created as a tex file.
#'
#' @param x an object of class FW.
#' @param sortBy character string indicating by which variable the estimates
#' should be sorted. Either \code{sens}(itivity), \code{genMean} (genotypic Mean) or
#' \code{mse} (mean squared error).
#' @param ... further arguments passed on from other functions - not used yet.
#' @param outfile a character string, the name and location of the output .pdf and .tex
#' file for the report. If \code{NULL} a report will be created in the current working
#' directory.
#'
#' @export
report.FW <- function(x,
                      sortBy = "sens",
                      ...,
                      outfile = NULL) {
  sortBy <- match.arg(arg = sortBy, choices = c("sens", "genMean", "mse"))
  createReport(x = x,
               reportName = "FWReport.Rnw",
               outfile = outfile,
               ...,
               sortBy = sortBy)
}


