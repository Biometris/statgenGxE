#' S3 class FW
#'
#' Function for creating objects of S3 class FW (Finlay-Wilkinson).\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param estimates A data.frame containing the estimated values.
#' @param anova A data.frame containing anova scores of the FW analysis.
#' @param envEffs A data.frame containing the environmental effects.
#' @param TD The object of class \code{\link{TD}} on which the analysis was
#' performed.
#' @param fittedGeno The fitted values for the genotypes.
#' @param trait A character value indicating the analysed trait.
#' @param nGeno A numerical value containing the number of genotypes in the
#' analysis.
#' @param nEnv A numerical value containing the number of environments in the
#' analysis.
#' @param tol A numerical value containing the tolerance used during the
#' analysis.
#' @param iter A numerical value containing the number of iterations for the
#' analysis to converge.
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
                     TD,
                     fittedGeno,
                     tol,
                     iter) {
  FW <- structure(list(estimates = estimates,
                       anova = anova,
                       envEffs = envEffs,
                       TD = TD,
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

#' @export
print.FW <- function(x, ...) {
  cat("Environmental effects",
      "\n===================\n")
  print(x$envEffs)
  cat("\nAnova",
      "\n=====\n")
  printCoefmat(x$anova, na.print = "")
  cat("\nEstimates",
      "\n=========\n")
  print(x$estimates, ..., row.names = FALSE)
}

#' @export
summary.FW <- function(object, ...) {
  print(object, ...)
}

#' Plot function for class FW
#'
#' Three types of plot can be made. A scatter plot for genotypic mean,
#' mse and sensitivity, a line plot with fitted lines for each genotype and
#' a trellis plot with individual slopes per genotype (for max 64 genotypes).
#' If there are more than 64 genotypes only the first 64 are plotted in the
#' trellis plot.
#'
#' @param x An object of class FW.
#' @param ... Further graphical parameters passed on to actual plot function.
#' @param plotType A character string indicating which plot should be made.
#' Either "scatter", "line" or "trellis" for creating a scatter
#' plot of genotypic means, mse and sensitivities, a plot of fitted lines for
#' each genotype or a trellis plot of the individual genotype slopes
#' respectively.
#' @param sorted A character string specifying whether the results should be
#' sorted in an increasing (or decreasing) order of sensitivities.

#' @return A plot depending on \code{plotType}.
#'
#' @examples
#' ## Run Finlay-Wilkinson analysis.
#' geFW <- gxeFw(TD = TDMaize, trait = "yld")
#' ## Create a scatter plot.
#' plot(geFW)
#' ## Create a line plot.
#' plot(geFW, plotType = "line")
#' ## Create a line plot.
#' plot(geFW, plotType = "trellis")
#'
#' @import graphics grDevices
#' @importFrom utils modifyList
#' @export
plot.FW <- function(x,
                    ...,
                    plotType = c("scatter", "line", "trellis"),
                    sorted = c("ascending", "descending", "none")) {
  plotType <- match.arg(plotType, several.ok = TRUE)
  sorted <- match.arg(sorted)
  dotArgs <- list(...)
  envEffs <- x$envEffs$effect
  TDTot <- Reduce(f = rbind, x = x$TD)
  if ("scatter" %in% plotType) {
    selCols = c(1, if (!all(is.na(x$estimates$mse))) 2, 3)
    scatterData <- setNames(x$estimates[, c("genMean", "mse", "sens")[selCols]],
                            c("Mean", "m.s.deviation", "Sensitivity")[selCols])
    ## Set arguments for plot.
    plotArgs <- list(x = scatterData, upper.panel = NULL,
                     main = paste0("Finlay & Wilkinson analysis for ", x$trait))
    ## Add and overwrite args with custom args from ...
    fixedArgs <- c("x", "title")
    plotArgs <- modifyList(plotArgs, dotArgs[!names(dotArgs) %in% fixedArgs])
    do.call(ifelse(!all(is.na(x$estimates$mse)), pairs, plot), args = plotArgs)
  } else if ("line" %in% plotType) {
    fVal <- tapply(X = x$fittedGeno, INDEX = TDTot[, c("trial", "genotype")],
                   FUN = mean, na.rm = TRUE)
    if (sorted == "none") {
      orderEnv <- 1:length(envEffs)
    } else {
      orderEnv <- order(envEffs, decreasing = (sorted == "descending"))
    }
    lineDat <- reshape2::melt(fVal)
    ## Set arguments for plot aesthetics.
    aesArgs <- list(x = rep(envEffs, nlevels(lineDat$genotype)),
                    y = "value", color = "genotype")
    fixedArgs <- c("x", "y", "color", "title")
    ## Add and overwrite args with custom args from ...
    aesArgs <- modifyList(aesArgs, dotArgs[!names(dotArgs) %in% fixedArgs])
    ## Create plot.
    ggplot2::ggplot(data = lineDat,
                    do.call(ggplot2::aes_string, args = aesArgs)) +
      ggplot2::geom_point() + ggplot2::geom_line(size = 0.5, alpha = 0.7) +
      ggplot2::scale_x_continuous(breaks = envEffs,
                                  labels = levels(lineDat$trial)) +
      ggplot2::theme(legend.position = "none",
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::ggtitle(ifelse(!is.null(dotArgs$title), dotArgs$title,
                              paste0("Finlay & Wilkinson analysis for ",
                                     x$trait))) +
      ggplot2::labs(x = "Environment", y = x$trait)
  } else if ("trellis" %in% plotType) {
    trellisData <- data.frame(genotype = TDTot$genotype,
                              trait = TDTot[[x$trait]],
                              fitted = x$fittedGen,
                              xEff = rep(envEffs, x$nGeno))
    if (x$nGeno > 64) {
      ## Select first 64 genotypes for plotting.
      first64 <- TDTot$genotype %in% levels(x$estimates$genotype)[1:64]
      trellisData <- trellisData[first64, ]
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
    ## Add and overwrite args with custom args from ...
    fixedArgs <- c("x", "data", "panel")
    plotArgs <- modifyList(plotArgs, dotArgs[!names(dotArgs) %in% fixedArgs])
    do.call(lattice::xyplot, args = plotArgs)
  }
}

#' Report method for class FW
#'
#' A pdf report will be created containing a summary of a Finlay-Wilkinson
#' analysis. Simultaneously the same report will be created as a tex file.
#'
#' @inheritParams report.AMMI
#'
#' @param x An object of class FW.
#' @param sortBy A character string indicating by which variable the estimates
#' should be sorted. Either \code{sens}(itivity), \code{genMean} (genotypic
#' Mean) or \code{mse} (mean squared error).
#'
#' @return A pdf and tex report.
#'
#' @examples
#' ## Run Finlay-Wilkinson analysis on TDMaize.
#' geFW <- gxeFw(TDMaize, trait = "yld")
#' \dontrun{
#' ## Create a report summarizing the results.
#' report(geFW, outfile = "./testReports/reportFW.pdf")
#' }
#'
#' @export
report.FW <- function(x,
                      sortBy = c("sens", "genMean", "mse"),
                      ...,
                      outfile = NULL) {
  sortBy <- match.arg(arg = sortBy)
  createReport(x = x, reportName = "FWReport.Rnw", outfile = outfile, ...,
               sortBy = sortBy)
}


