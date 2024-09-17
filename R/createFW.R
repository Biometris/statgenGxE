#' S3 class FW
#'
#' Function for creating objects of S3 class FW (Finlay-Wilkinson).\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param estimates A data.frame containing the estimated values.
#' @param anova A data.frame containing anova scores of the FW analysis.
#' @param envEffs A data.frame containing the environmental effects.
#' @param TD The object of class \code{\link[statgenSTA]{TD}} on which the
#' analysis was performed.
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
  if (all(x$estimates[["Rank"]] == 1:nrow(x$estimates))) {
    cat("\nMost sensitive genotypes")
    cat("\n========================\n")
  } else if (all(x$estimates[["Rank"]] == nrow(x$estimates):1)) {
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
#' square root of mean squared deviation and sensitivity, a line plot with
#' fitted lines for each genotype, a trellis plot with individual slopes per
#' genotype and a scatter plot of fitted values in the worst and best trial.\cr
#' It is possible to select genotypes for the trellis plot using the
#' \code{genotypes} parameter. If there are more than 64 genotypes, only the
#' first 64 are plotted in the trellis plot.
#'
#' @param x An object of class FW.
#' @param ... Not used.
#' @param plotType A character string indicating which plot should be made.
#' Either "scatter", "line" or "trellis" for creating a scatter
#' plot of genotypic means, mse and sensitivities, a plot of fitted lines for
#' each genotype or a trellis plot of the individual genotype slopes
#' respectively.
#' @param order A character string specifying whether the results in the line
#' plot should be ordered in an increasing (or decreasing) order of
#' sensitivities. Ignored if \code{plotType} is not "line".
#' @param response A character string specifying whether in the line plot the
#' "predicted" or the "observed" data should be plotted. Ignored if
#' \code{plotType} is not "line".
#' @param colorGenoBy A character string indicating a column in the \code{TD}
#' used as input for the Finlay Wilkinson analysis by which the genotypes
#' should be colored. If \code{NULL} all genotypes will be colored differently.
#' @param colGeno A character vector with plot colors for the genotypes. A
#' single color when \code{colorGenoBy = NULL}, a vector of colors otherwise.
#' @param genotypes An optional character string containing the genotypes to
#' be plotted in the trellis plot. If \code{NULL} all genotypes are plotted.
#' If more than 64 genotypes are selected, only the first 64 are plotted.
#' @param title A character string used a title for the plot.
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
#' ## Create a line plot showing observed data value for genotypes and fitted lines.
#' ## Display trials in descending order.
#' plot(geFW, plotType = "line", order = "descending", response = "observed")
#'
#' \donttest{
#' ## Create a trellis plot.
#' plot(geFW, plotType = "trellis")
#'
#' ## Create a scatter plot of fitted values for the worst and best trials.
#' plot(geFW, plotType = "scatterFit")
#' }
#'
#' @family Finlay-Wilkinson
#'
#' @importFrom grDevices topo.colors
#' @export
plot.FW <- function(x,
                    ...,
                    plotType = c("scatter", "line", "trellis", "scatterFit"),
                    order = c("ascending", "descending"),
                    response = c("predicted", "observed"),
                    colorGenoBy = NULL,
                    colGeno = NULL,
                    genotypes = NULL,
                    title = paste("Finlay & Wilkinson analysis for", x$trait),
                    output = TRUE) {
  plotType <- match.arg(plotType)
  chkChar(title, len = 1)
  trait <- x$trait
  envEffs <- x$envEffs[c("Trial", "EnvMean")]
  TDTot <- do.call(rbind, x$TD)
  if (!is.null(colorGenoBy)) {
    chkCol(colorGenoBy, TDTot)
  }
  chkChar(colGeno)
  fittedGeno <- x$fittedGeno
  colnames(fittedGeno)[colnames(fittedGeno) == "fittedValue"] <- "fitted"
  TDTot <- merge(fittedGeno, TDTot[c("trial", "genotype", colorGenoBy, trait)],
                 all.x = TRUE)
  TDTot[["genoMean"]] <- ave(x = TDTot[[trait]], TDTot[["trial"]],
                             TDTot[["genotype"]], FUN = mean)
  genoDat <- unique(TDTot[c("trial", "genotype", "genoMean", "fitted",
                            colorGenoBy)])
  if (!is.null(colorGenoBy)) {
    if (!is.factor(genoDat[[colorGenoBy]])) {
      genoDat[[colorGenoBy]] <- as.factor(genoDat[[colorGenoBy]])
    }
    genoDat <- genoDat[order(genoDat[[colorGenoBy]]), ]
  } else {
    genoDat[[".colorGenoBy"]] <- factor(1)
    colorGenoBy <- ".colorGenoBy"
  }
  ## Get number of colors.
  nColGeno <- nlevels(genoDat[[colorGenoBy]])
  if (length(colGeno) == 0) {
    ## Defaults to black for one color for genotypes.
    ## For more than one colors from statgen.genoColors are used.
    ## Fall back to topo.colors if number of colors in option is too small.
    if (nColGeno == 1) {
      colGeno <- "black"
    } else if (length(getOption("statgen.genoColors")) >= nColGeno) {
      colGeno <- getOption("statgen.genoColors")[1:nColGeno]
    } else {
      colGeno <- topo.colors(nColGeno)
    }
  } else {
    nColGenoArg <- length(colGeno)
    if (nColGenoArg != nColGeno) {
      stop("Number of colors provided doesn't match number of genotype groups:",
           "\n", nColGenoArg, " colors provided, ", nColGeno,
           " groups in data.\n")
    }
  }
  if (plotType == "scatter") {
    selCols = c(1:2, if (!all(is.na(x$estimates$MSdeviation))) 3, 4)
    scatterDat <- x$estimates[, c("Genotype", "GenMean", "MSdeviation",
                                  "Sens")[selCols]]
    scatterDat <- merge(scatterDat, genoDat[, c("genotype", colorGenoBy)],
                        by.x = "Genotype", by.y = "genotype")
    scatterDat <- ggplot2::remove_missing(scatterDat, na.rm = TRUE)
    ## Create plot of mean x mse. No x axis because of position in grid.
    p1 <- ggplot2::ggplot(data = scatterDat,
                          ggplot2::aes(x = .data[["GenMean"]],
                                       y = sqrt(.data[["MSdeviation"]]),
                                       color = .data[[colorGenoBy]])) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = colGeno) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank()) +
      ggplot2::labs(x = "Mean", y = "Square root of\n Mean Squared Deviation")
    if (colorGenoBy != ".colorGenoBy") {
      ## Build plot to extract legend.
      p1Gtable <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(p1))
      legendPos <- sapply(X = p1Gtable$grobs, FUN = `[[`, "name") == "guide-box"
      legend <- p1Gtable$grobs[[which(legendPos)]]
    } else {
      legend <- NULL
    }
    ## Now remove the legend from the plot.
    p1 <- p1 + ggplot2::theme(legend.position = "none")
    ## Create plot of mean x sensitivity.
    p2 <- ggplot2::ggplot(data = scatterDat,
                          ggplot2::aes(x = .data[["GenMean"]],
                                       y = .data[["Sens"]],
                                       color = .data[[colorGenoBy]])) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = colGeno) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(x = "Mean", y = "Sensitivity")
    ## Create plot of mse x sensitivity. No y axis because of position in grid.
    p3 <- ggplot2::ggplot(data = scatterDat,
                          ggplot2::aes(x = sqrt(.data[["MSdeviation"]]),
                                       y = .data[["Sens"]],
                                       color = .data[[colorGenoBy]])) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = colGeno) +
      ggplot2::theme(legend.position = "none",
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank()) +
      ggplot2::labs(x = "Square root of Mean Squared Deviation",
                    y = "Sensitivity")
    ## Create empty plot for top right grid position.
    pEmpty <- ggplot2::ggplot() +
      ggplot2::theme(panel.background = ggplot2::element_blank())
    # Convert to Grobs to make alignment of axis possible.
    p1Gr <- ggplot2::ggplotGrob(p1)
    p2Gr <- ggplot2::ggplotGrob(p2)
    p3Gr <- ggplot2::ggplotGrob(p3)
    pEmpty <- ggplot2::ggplotGrob(pEmpty)
    # Create grid by first binding rows to assure axis alignment and then
    # by columns.
    c1 <- gridExtra::gtable_rbind(p1Gr, p2Gr)
    c2 <- gridExtra::gtable_rbind(pEmpty, p3Gr)
    tot <- gridExtra::gtable_cbind(c1, c2)
    if (output) {
      # grid.arrange automatically plots the results.
      # Assign to variable to avoid double output plot.
      p <- gridExtra::grid.arrange(tot, right = legend, top = title)
    }
    invisible(list(p1 = p1, p2 = p2, p3 = p3))
    ## Set arguments for plot.
  } else if (plotType == "line") {
    order <- match.arg(order)
    response <- match.arg(response)
    lineDat <- merge(genoDat, envEffs, by.x = "trial", by.y = "Trial")
    lineDat <- ggplot2::remove_missing(lineDat, na.rm = TRUE)
    ## Set arguments for plot aesthetics.
    yVar <- ifelse(response == "observed", "genoMean", "fitted")
    ## Order descending can be achieved by reversing the x-axis.
    if (order == "descending") {
      xTrans <- "reverse"
    } else {
      xTrans <- "identity"
    }
    plotLims <- range(c(lineDat[["EnvMean"]], lineDat[[yVar]]))
    colorVar <- if (colorGenoBy == ".colorGenoBy") "genotype" else
      enquote(colorGenoBy)
    ## Create plot.
    p <- ggplot2::ggplot(
      data = lineDat,
      ggplot2::aes(x = .data[["EnvMean"]], y = .data[[yVar]],
                   group = .data[["genotype"]],
                   color = .data[[colorVar]])) +
      ggplot2::geom_point() +
      ggplot2::geom_line(ggplot2::aes(y = .data[["fitted"]]),
                         linewidth = 0.5, alpha = 0.7) +
      ggplot2::scale_x_continuous(trans = xTrans,
                                  sec.axis = ggplot2::dup_axis(name = "Environment",
                                                               breaks = envEffs[["EnvMean"]],
                                                               labels = envEffs[["Trial"]])) +
      ggplot2::geom_vline(xintercept = mean(TDTot[[trait]], na.rm = TRUE),
                          color = "red", linetype = "dashed") +
      ggplot2::coord_equal(xlim = plotLims, ylim = plotLims) +
      ggplot2::theme(legend.position = if (colorGenoBy == ".colorGenoBy") "none" else "right",
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     axis.text.x.top = ggplot2::element_text(angle = 90,
                                                             hjust = 1)) +
      ggplot2::labs(title = title, x = NULL, y = trait)
    if (colorGenoBy != ".colorGenoBy") {
      p <- p + ggplot2::scale_color_manual(values = colGeno)
    }
    if (output) {
      plot(p)
    }
    invisible(p)
  } else if (plotType == "trellis") {
    if (!is.null(genotypes) && !all(genotypes %in% TDTot[["genotype"]])) {
      stop("All genotypes should be in TD.\n")
    }
    trellisDat <- merge(genoDat, envEffs, by.x = "trial", by.y = "Trial")
    if (!is.null(genotypes)) {
      trellisDat <- trellisDat[trellisDat[["genotype"]] %in% genotypes, ]
      trellisDat <- droplevels(trellisDat)
    }
    if (nlevels(trellisDat[["genotype"]]) > 64) {
      ## Select first 64 genotypes for plotting.
      trellisDat <- trellisDat[trellisDat[["genotype"]] %in%
                                 levels(trellisDat[["genotype"]])[1:64], ]
    }
    trellisDat <- ggplot2::remove_missing(trellisDat, na.rm = TRUE)
    ## The data needs to be ordered for the lines to be drawn properly.
    trellisDat <- trellisDat[order(trellisDat[["genotype"]],
                                   trellisDat[["EnvMean"]]), ]
    p <- ggplot2::ggplot(data = trellisDat,
                         ggplot2::aes(x = .data[["EnvMean"]],
                                      y = .data[["genoMean"]])) +
      ggplot2::geom_point() +
      ggplot2::geom_line(data = trellisDat,
                         ggplot2::aes(x = .data[["EnvMean"]],
                                      y = .data[["fitted"]])) +
      ggplot2::facet_wrap(facets = "genotype") +
      ggplot2::labs(x = "Environment", y = trait) +
      ggplot2::ggtitle(title) +
      ggplot2::theme(legend.position = "none",
                     plot.title = ggplot2::element_text(hjust = 0.5),
                     panel.spacing = ggplot2::unit(.2, "cm"),
                     axis.text = ggplot2::element_text(size = 6))
    if (output) {
      plot(p)
    }
    invisible(p)
  } else if (plotType == "scatterFit") {
    ## Get worst and best trials.
    trialMin <- as.character(envEffs[which.min(envEffs[["EnvMean"]]), "Trial"])
    trialMax <- as.character(envEffs[which.max(envEffs[["EnvMean"]]), "Trial"])
    ## Construct plot data, fitted values for worst and best trials.
    plotDat <- TDTot[c("trial", "fitted")]
    plotDat <- data.frame(genotype = levels(TDTot[["genotype"]]),
                          trMin = plotDat[plotDat[["trial"]] == trialMin, "fitted"],
                          trMax = plotDat[plotDat[["trial"]] == trialMax, "fitted"])
    plotDat <- merge(plotDat, genoDat[, c("genotype", colorGenoBy)])
    ## Create scatter plot of fitted values.
    p <- ggplot2::ggplot(data = plotDat,
                         ggplot2::aes(x = .data[["trMin"]],
                                      y = .data[["trMax"]],
                                      color = .data[[colorGenoBy]])) +
      ggplot2::geom_point(na.rm = TRUE,
                          show.legend = colorGenoBy != ".colorGenoBy") +
      ggplot2::scale_color_manual(values = colGeno) +
      ggplot2::labs(x = paste("Fitted values for worst trial:", trialMin),
                    y = paste("Fitted values for best trial:", trialMax)) +
      ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    if (output) {
      plot(p)
    }
    invisible(p)
  }
}

#' Extract fitted values.
#'
#' Extract the fitted values for a fitted Finlay-Wilkinson model.
#'
#' @param object An object of class FW
#' @param ... Not used.
#'
#' @return A data.frame with fitted values.
#'
#' @examples
#' ## Run Finlay-Wilkinson analysis.
#' geFW <- gxeFw(TD = TDMaize, trait = "yld")
#'
#' ## Extract fitted values.
#' fitFW <- fitted(geFW)
#' head(fitFW)
#'
#' @family Finlay-Wilkinson
#'
#' @export
fitted.FW <- function(object,
                      ...) {
  return(object$fittedGeno)
}

#' Extract residuals.
#'
#' Extract the residuals for a fitted Finlay-Wilkinson model.
#'
#' @param object An object of class FW
#' @param ... Not used.
#'
#' @return A data.frame with residuals.
#'
#' @examples
#' ## Run Finlay-Wilkinson analysis.
#' geFW <- gxeFw(TD = TDMaize, trait = "yld")
#'
#' ## Extract residuals.
#' residFW <- residuals(geFW)
#' head(residFW)
#'
#' @family Finlay-Wilkinson
#'
#' @export
residuals.FW <- function(object,
                         ...) {
  trait <- object$trait
  fittedGeno <- object$fittedGeno
  TDTot <- do.call(rbind, object$TD)
  residGeno <- merge(fittedGeno, TDTot, by = c("trial", "genotype"))
  residGeno[["residual"]] <- residGeno[["fittedValue"]] - residGeno[[trait]]
  residGeno <- residGeno[c("trial", "genotype", "residual")]
  return(residGeno)
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
#'
#' \donttest{
#' ## Create a report summarizing the results.
#' report(geFW, outfile = tempfile(fileext = ".pdf"))
#' }
#'
#' @family Finlay-Wilkinson
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
