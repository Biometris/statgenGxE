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
#' @seealso \code{\link{plot.FW}}, \code{\link{report.FW}}
#'
#' @name FW
NULL

#' @rdname FW
#' @keywords internal
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

#' @importFrom utils head
#' @export
print.FW <- function(x, ...) {
  cat("Environmental effects",
      "\n=====================\n")
  print(x$envEffs)
  cat("\nAnova",
      "\n=====\n")
  print(x$anova)
  if (all(x$estimates[["rank"]] == 1:nrow(x$estimates))) {
    cat("\nMost sensitive genotypes")
    cat("\n========================\n")
  } else if (all(x$estimates[["rank"]] == nrow(x$estimates):1)) {
    cat("\nLeast sensitive genotypes")
    cat("\n=========================\n")
  } else {
    cat("\nFirst five genotypes")
    cat("\n====================\n")
  }
  print(head(x$estimates, 5),  ..., row.names = FALSE)
}

#' @export
summary.FW <- function(object, ...) {
  print(object, ...)
}

#' Plot function for class FW
#'
#' Four types of plot can be made. A scatter plot for genotypic mean,
#' mean squared deviation and sensitivity, a line plot with fitted lines for
#' each genotype, a trellis plot with individual slopes per genotype and a
#' scatter plot of fitted values in the worst and best trial.\cr
#' It is possible to select genotypes for the trellis plot using the
#' \code{genotypes} parameter. If there are more than 64 genotypes, only the
#' first 64 are plotted in the trellis plot.
#'
#' @param x An object of class FW.
#' @param ... Further graphical parameters passed on to actual plot function.
#' @param plotType A character string indicating which plot should be made.
#' Either "scatter", "line" or "trellis" for creating a scatter
#' plot of genotypic means, mse and sensitivities, a plot of fitted lines for
#' each genotype or a trellis plot of the individual genotype slopes
#' respectively.
#' @param order A character string specifying whether the results in the line
#' plot should be ordered in an increasing (or decreasing) order of
#' sensitivities. Ignored if \code{plotType} is not "line".
#' @param lineTrait A character string specifying whether in the line plot the
#' "fitted" or the "raw" data should be plotted. Ignored if \code{plotType} is
#' not "line".
#' @param genotypes An optional character string containing the genotypes to
#' be plotted in the trellis plot. If \code{NULL} all genotypes are plotted.
#' If more than 64 genotypes are selected, only the first 64 are plotted.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE}, only a list of ggplot objects is invisibly returned.
#'
#' @return A plot depending on \code{plotType}.
#'
#' @examples
#' ## Run Finlay-Wilkinson analysis.
#' geFW <- gxeFw(TD = TDMaize, trait = "yld")
#'
#' ## Create a scatter plot.
#' plot(geFW)
#'
#' ## Create a line plot.
#' plot(geFW, plotType = "line")
#'
#' ## Create a line plot showing raw data value for genotypes and fitted lines.
#' ## Display trials in descending order.
#' plot(geFW, plotType = "line", order = "descending", lineTrait = "raw")
#'
#' ## Create a trellis plot.
#' plot(geFW, plotType = "trellis")
#'
#' ## Create a scatter plot of fitted values for the worst and best trials.
#' plot(geFW, plotType = "scatterFit")
#'
#' @export
plot.FW <- function(x,
                    ...,
                    plotType = c("scatter", "line", "trellis", "scatterFit"),
                    order = c("ascending", "descending"),
                    lineTrait = c("fitted", "raw"),
                    genotypes = NULL,
                    output = TRUE) {
  plotType <- match.arg(plotType)
  trait <- x$trait
  dotArgs <- list(...)
  envEffs <- x$envEffs[c("trial", "envEff")]
  TDTot <- Reduce(f = rbind, x = x$TD)
  TDTot[["genoMean"]] <- ave(x = TDTot[[trait]], TDTot[["trial"]],
                             TDTot[["genotype"]], FUN = mean)
  genoMeans <- unique(TDTot[c("trial", "genotype", "genoMean")])
  plotTitle <- ifelse(!is.null(dotArgs$title), dotArgs$title,
                      paste0("Finlay & Wilkinson analysis for ", trait))
  if (plotType == "scatter") {
    selCols = c(1:2, if (!all(is.na(x$estimates$MSdeviation))) 3, 4)
    scatterDat <- setNames(x$estimates[, c("genotype", "genMean",
                                           "MSdeviation", "sens")[selCols]],
                           c("genotype", "Mean", "MSDeviation",
                             "Sensitivity")[selCols])
    scatterDat <- remove_missing(scatterDat, na.rm = TRUE)
    ## Create plot of mean x mse. No x axis because of position in grid.
    p1 <- ggplot(data = scatterDat,
                 aes_string(x = "Mean", y = "MSDeviation")) +
      geom_point() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(y = "Mean Squared Deviation")
    ## Create plot of mean x sensitivity.
    p2 <- ggplot(data = scatterDat, aes_string(x = "Mean", y = "Sensitivity")) +
      geom_point()
    ## Create plot of mse x sensitivity. No y axis because of position in grid.
    p3 <- ggplot(data = scatterDat, aes_string(x = "MSDeviation",
                                               y = "Sensitivity")) +
      geom_point() +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      labs(x = "Mean Squared Deviation")
    ## Create empty plot for top right grid position.
    pEmpty <- ggplot() + theme(panel.background = element_blank())
    # Convert to Grobs to make alignment of axis possible.
    p1Gr <- ggplotGrob(p1)
    p2Gr <- ggplotGrob(p2)
    p3Gr <- ggplotGrob(p3)
    pEmpty <- ggplotGrob(pEmpty)
    # Create grid by first binding rows to assure axis alignment and then
    # by columns.
    c1 <- gridExtra::gtable_rbind(p1Gr, p2Gr)
    c2 <- gridExtra::gtable_rbind(pEmpty, p3Gr)
    tot <- gridExtra::gtable_cbind(c1, c2)
    if (output) {
      # grid.arrange automatically plots the results.
      # Assign to variable to avoid double output plot.
      p <- gridExtra::grid.arrange(tot, top = plotTitle)
    }
    invisible(list(p1 = p1, p2 = p2, p3 = p3))
    ## Set arguments for plot.
  } else if (plotType == "line") {
    order <- match.arg(order)
    lineTrait <- match.arg(lineTrait)
    lineDat <- data.frame(genotype = genoMeans[["genotype"]],
                          trait = genoMeans[["genoMean"]],
                          trial = rep(x = envEffs[["trial"]], each = x$nGeno),
                          fitted = x$fittedGeno,
                          envEff = rep(x = envEffs[["envEff"]], each = x$nGeno))
    lineDat <- remove_missing(lineDat, na.rm = TRUE)
    ## Set arguments for plot aesthetics.
    yVar <- ifelse(lineTrait == "raw", "trait", "fitted")
    aesArgs <- list(x = "envEff", y = yVar, color = "genotype")
    fixedArgs <- c("x", "y", "color", "title")
    ## Add and overwrite args with custom args from ...
    aesArgs <- utils::modifyList(aesArgs,
                                 dotArgs[!names(dotArgs) %in% fixedArgs])
    ## Order descending can be achieved by reversing the x-axis.
    if (order == "descending") {
      xTrans <- "reverse"
    } else {
      xTrans <- "identity"
    }
    ## Create plot.
    p <- ggplot(data = lineDat, do.call(aes_string, args = aesArgs)) +
      geom_point() +
      geom_line(aes_string(y = "fitted"), size = 0.5, alpha = 0.7) +
      scale_x_continuous(breaks = envEffs[["envEff"]], minor_breaks = NULL,
                         labels = envEffs[["trial"]], trans = xTrans) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = plotTitle, x = "Environment", y = trait)
    if (output) {
      plot(p)
    }
    invisible(p)
  } else if (plotType == "trellis") {
    if (!is.null(genotypes) && !all(genotypes %in% TDTot[["genotype"]])) {
      stop("All genotypes should be in TD.\n")
    }
    trellisDat <- data.frame(genotype = genoMeans[["genotype"]],
                             trait = genoMeans[["genoMean"]],
                             fitted = x$fittedGeno,
                             envEff = rep(x = envEffs$envEff, each = x$nGeno))
    if (!is.null(genotypes)) {
      trellisDat <- trellisDat[trellisDat[["genotype"]] %in% genotypes, ]
      trellisDat <- droplevels(trellisDat)
    }
    if (nlevels(trellisDat[["genotype"]]) > 64) {
      ## Select first 64 genotypes for plotting.
      trellisDat <- trellisDat[trellisDat[["genotype"]] %in%
                                 levels(trellisDat[["genotype"]])[1:64], ]
    }
    trellisDat <- remove_missing(trellisDat, na.rm = TRUE)
    ## The data needs to be ordered for the lines to be drawn properly.
    trellisDat <- trellisDat[order(trellisDat[["genotype"]],
                                   trellisDat[["envEff"]]), ]
    p <- ggplot(data = trellisDat,
                aes_string(x = "envEff", y = "trait")) +
      geom_point() +
      geom_line(data = trellisDat, aes_string(x = "envEff", y = "fitted")) +
      facet_wrap(facets = "genotype") +
      labs(x = "Environment", y = trait) +
      ggtitle(plotTitle) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            panel.spacing = unit(.2, "cm"),
            axis.text = element_text(size = 6))
    if (output) {
      plot(p)
    }
    invisible(p)
  } else if (plotType == "scatterFit") {
    ## Get worst and best trials.
    trialMin <- as.character(envEffs[which.min(envEffs[["envEff"]]), "trial"])
    trialMax <- as.character(envEffs[which.max(envEffs[["envEff"]]), "trial"])
    ## Construct plot data, fitted values for worst and best trials.
    plotDat <- data.frame(trial = TDTot[["trial"]], fitted = x$fittedGeno)
    plotDat <- data.frame(genotype = TDTot[["genotype"]],
                          trMin = plotDat[plotDat[["trial"]] == trialMin, "fitted"],
                          trMax = plotDat[plotDat[["trial"]] == trialMax, "fitted"])
    ## Create scatter plot of fitted values.
    p <- ggplot(data = plotDat,
                aes_string(x = "trMin", y = "trMax")) +
      geom_point() +
      labs(x = paste("Fitted values for worst trial:", trialMin),
           y = paste("Fitted values for best trial:", trialMax)) +
      ggtitle(plotTitle) +
      theme(plot.title = element_text(hjust = 0.5))
    if (output) {
      plot(p)
    }
    invisible(p)
  }
}

#' Report method for class FW
#'
#' A pdf report will be created containing a summary of an FW object.
#' Simultaneously the same report will be created as a tex file.
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
  ## Checks.
  if (nchar(Sys.which("pdflatex")) == 0) {
    stop("An installation of LaTeX is required to create a pdf report.\n")
  }
  sortBy <- match.arg(arg = sortBy)
  createReport(x = x, reportName = "FWReport.Rnw", reportPackage = "statgenGxE",
               outfile = outfile, ..., sortBy = sortBy)
}
