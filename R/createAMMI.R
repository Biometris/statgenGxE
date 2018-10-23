#' S3 class AMMI
#'
#' Function for creating objects of S3 class AMMI.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param envScores A matrix containing environmental scores.
#' @param genoScores A matrix containing genotypic scores.
#' @param importance A data.frame containing the importance of the principal
#' components.
#' @param anova A data.frame containing anova scores of the AMMI analysis.
#' @param fitted A matrix containing fitted values from the AMMI model.
#' @param trait A character string indicating the analyzed trait.
#' @param envMean A numerical vector containing the environmental means.
#' @param genoMean A numerical vector containing the genotypic means.
#' @param overallMean A numerical value containing the overall mean.
#' @param GGE Has a GGE analysis been performed?
#' @param byYear Has the analysis been performed by year?
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.AMMI}}, \code{\link{report.AMMI}}
#'
#' @name AMMI
NULL

#' @rdname AMMI
#' @keywords internal
createAMMI <- function(envScores,
                       genoScores,
                       importance,
                       anova,
                       fitted,
                       trait,
                       envMean,
                       genoMean,
                       overallMean,
                       dat,
                       GGE,
                       byYear) {
  AMMI <- structure(list(envScores = envScores,
                         genoScores = genoScores,
                         importance = importance,
                         anova = anova,
                         fitted = fitted,
                         trait = trait,
                         envMean = envMean,
                         genoMean = genoMean,
                         overallMean = overallMean,
                         dat = dat,
                         GGE = GGE,
                         byYear = byYear),
                    class = "AMMI",
                    timestamp = Sys.time())
  return(AMMI)
}

#' @export
print.AMMI <- function(x, ...) {
  cat("Principal components",
      "\n====================\n")
  if (x$byYear) {
    years <- names(x$importance)
    for (year in years) {
      cat(paste(year, "\n"))
      print(x$importance[[year]][, 1:ncol(x$envScores[[year]]), drop = FALSE])
      cat("\n")
    }
  } else {
    print(x$importance[, 1:ncol(x$envScores), drop = FALSE])
  }
  cat("\nAnova",
      "\n=====\n")
  if (x$byYear) {
    for (year in years) {
      cat(paste(year, "\n"))
      printCoefmat(x$anova[[year]], na.print = "")
      cat("\n")
    }
  } else {
    printCoefmat(x$anova, na.print = "")
  }
  cat("\nEnvironment scores",
      "\n==================\n")
  if (x$byYear) {
    for (year in years) {
      cat(paste(year, "\n"))
      print(x$envScores[[year]], ...)
      cat("\n")
    }
  } else {
    print(x$envScores, ...)
  }
  cat("\nGenotypic scores",
      "\n================\n")
  if (x$byYear) {
    for (year in years) {
      cat(paste(year, "\n"))
      print(x$genoScores[[year]], ...)
      cat("\n")
    }
  } else {
    print(x$genoScores, ...)
  }
}

#' @export
summary.AMMI <- function(object, ...) {
  print(object, ...)
}

#' Plot function for class AMMI
#'
#' Two types of biplot can be made. A biplot of genotype and environment
#' means vs PC1 (AMMI1) or a biplot of genotypes and environment interaction
#' with PC1 and PC2 (AMMI2).\cr\cr
#' If the AMMI analysis was done by year a separate plot will be made for
#' every year in the data. If for some of the years it is not possible to
#' make the plot because the number of principal components for those year is
#' lower than the number specified as secondary axis those years will be
#' skipped when plotting. If this would be the case for all years the function
#' returns an error.
#'
#' @param x An object of class AMMI
#' @param ... Further graphical parameters passed on to actual plot function.
#' @param plotType A character string indicating which plot should be made.
#' Either "AMMI1" or "AMMI2" for an AMMI1 biplot (genotype and
#' environment means vs PC1) or an AMMI2 biplot (genotypes and environment
#' interaction with PC1 and PC2) respectively.
#' @param scale A numerical value. The variables are scaled by
#' \code{lambda ^ scale} and the observations by \code{lambda ^ (1 - scale)}
#' where \code{lambda} are the singular values computed by
#' \code{\link[stats]{princomp}} in \code{\link{gxeAmmi}}. Normally
#' \code{0 <= scale <= 1}, and a warning will be issued if the specified
#' scale is outside this range.
#' @param col A vector with plot colors for genotype and environment. This can
#' either be named colors or color numbers.
#' @param primAxis A character string indicating the principal component to be
#' plotted on the primary axis of the AMMI2 plot. Has to be given as
#' \code{"PCn"} where n is the number of the principal component.
#' @param secAxis A character string indicating the principal component to be
#' plotted on the secondary axis of the AMMI2 plot. Has to be given as
#' \code{"PCn"} where n is the number of the principal component. n Has to
#' differ from \code{primAxis}.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @return A biplot depending on \code{plotType}.
#'
#' @examples
#' ## Run AMMI analysis.
#' geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
#' ## Create a biplot of genotype and environment means vs PC1.
#' plot(geAmmi)
#' ## Create a biplot of genotypes and environment interaction with PC1 and PC2.
#' plot(geAmmi, plotType = "AMMI2")
#'
#' @import graphics grDevices
#' @importFrom utils modifyList
#' @export
plot.AMMI <- function(x,
                      ...,
                      plotType = c("AMMI1", "AMMI2"),
                      scale = 1,
                      colorBy = NULL,
                      col = c("black", "red"),
                      primAxis = "PC1",
                      secAxis = "PC2",
                      output = TRUE) {
  ## Checks.
  if (!is.numeric(scale) || length(scale) > 1) {
    stop("scale should be a single numerical value.\n")
  }
  if (scale < 0 || scale > 1) {
    warning("Scale is outside [0, 1].\n", call. = FALSE)
  }
  if (length(col) != 2) {
    stop("col should contain exactly two colors.\n")
  }
  plotType <- match.arg(plotType)
  if (!is.null(colorBy) && (!is.character(colorBy) || length(colorBy) > 1)) {
    stop("colorBy should be a single character string.\n")
  }
  if (!is.null(colorBy)) {
    if (x$byYear) {
      for (dat in x$dat) {
        if (!hasName(x = dat, name = colorBy)) {
          stop("colorBy should be a column in data.\n")
        }
        colTab <- unique(dat[c("genotype", colorBy)])
        if (nrow(colTab) != nlevels(dat$genotype)) {
          stop("colorBy should have exactly one value per genotype")
        }
      }
    } else {
      if (!hasName(x = x$dat, name = colorBy)) {
        stop("colorBy should be a column in data.\n")
      }
      colTab <- unique(x$dat[c("genotype", colorBy)])
      if (nrow(colTab) != nlevels(x$dat$genotype)) {
        stop("colorBy should have exactly one value per genotype")
      }
    }
  }
  dotArgs <- list(...)
  if (plotType == "AMMI1") {
    if (x$byYear) {
      ## Create a list of AMMI1 plots.
      p <- sapply(X = names(x$envScores), FUN = function(year) {
        plotAMMI1(loadings = x$envScores[[year]], scores = x$genoScores[[year]],
                  importance = x$importance[[year]],
                  overallMean = x$overallMean[[year]],
                  genoMean = x$genoMean[[year]], envMean = x$envMean[[year]],
                  trait = x$trait, dat = x$dat[[year]], GGE = x$GGE, year = year,
                  scale = scale, col = col, colorBy = colorBy)
      }, simplify = FALSE)
    } else {
      ## Create a single AMMI1 plot.
      p <- plotAMMI1(loadings = x$envScores, scores = x$genoScores,
                     importance = x$importance, overallMean = x$overallMean,
                     genoMean = x$genoMean, envMean = x$envMean,
                     trait = x$trait, dat = x$dat, GGE = x$GGE, scale = scale,
                     col = col, colorBy = colorBy)
    }
  } else if (plotType == "AMMI2") {
    if (!is.character(primAxis) || length(primAxis) > 1 ||
        substring(text = primAxis, first = 1, last = 2) != "PC") {
      stop("primAxis should be a single character string starting with PC.\n")
    }
    nPC1 <- suppressWarnings(as.numeric(substring(text = primAxis, first = 3)))
    if (is.na(nPC1)) {
      stop(paste("Invalid value provided for primAxis Make sure the value is",
                 "of the form primAxis = 'PCn' where n is the principal",
                 "component to plot on the secondary axis.\n"))
    }
    if (!is.character(secAxis) || length(secAxis) > 1 ||
        substring(text = secAxis, first = 1, last = 2) != "PC") {
      stop("secAxis should be a single character string starting with PC.\n")
    }
    nPC2 <- suppressWarnings(as.numeric(substring(text = secAxis, first = 3)))
    if (is.na(nPC2)) {
      stop(paste("Invalid value provided for secAxis. Make sure the value is",
                 "of the form secAxis = 'PCn' where n is the principal",
                 "component to plot on the secondary axis.\n"))
    }
    if (nPC1 == nPC2) {
      stop("primAxis should differ from secAxis.\n")
    }
    if (x$byYear) {
      nPCs <- sapply(X = x$envScores, FUN = ncol)
      maxPC <- max(nPCs)
      if (nPC1 > maxPC) {
        stop(paste0("Highest number of principal components is ", maxPC,
                    ". Plotting of PC", nPC1, " is not possible.\n"))
      }
      if (nPC2 > maxPC) {
        stop(paste0("Highest number of principal components is ", maxPC,
                    ". Plotting of PC", nPC2, " is not possible.\n"))
      }
      p <- sapply(X = names(x$envScores)[nPCs >= max(nPC1, nPC2)],
                  FUN = function(year) {
                    plotAMMI2(loadings = x$envScores[[year]],
                              scores = x$genoScores[[year]],
                              importance = x$importance[[year]],
                              trait = x$trait, dat = x$dat[[year]], GGE = x$GGE,
                              year = year, primAxis = primAxis,
                              secAxis = secAxis, scale = scale, col = col,
                              colorBy = colorBy)
                  }, simplify = FALSE)


    } else {
      if (nPC1 > ncol(x$envScores)) {
        stop(paste0("AMMI was run with ", ncol(x$envScores), " principal ",
                    "components. Plotting of PC", nPC1, " is not possible.\n"))
      }
      if (nPC2 > ncol(x$envScores)) {
        stop(paste0("AMMI was run with ", ncol(x$envScores), " principal ",
                    "components. Plotting of PC", nPC2, " is not possible.\n"))
      }
      ## Create a single AMMI2 plot.
      p <- plotAMMI2(loadings = x$envScores, scores = x$genoScores,
                     importance = x$importance, trait = x$trait, dat = x$dat,
                     GGE = x$GGE, primAxis = primAxis, secAxis = secAxis,
                     scale = scale, col = col, colorBy = colorBy)
    }
  }
  if (output) {
    if (x$byYear) {
      lapply(X = p, FUN = plot)
    } else {
      plot(p)
    }
  }
  invisible(p)
}

#' Helper function for creating AMMI1 plot
#' @keywords internal
plotAMMI1 <- function(loadings,
                      scores,
                      importance,
                      overallMean,
                      genoMean,
                      envMean,
                      trait,
                      dat,
                      GGE,
                      year = "",
                      scale,
                      col,
                      colorBy) {
  percPC1 <- round(importance[2, 1] * 100, 1)
  ## Calculate lambda scale
  lam <- importance[1, 1]
  lam <- lam ^ scale
  ## Put overallMean in variable since using x$ in plot produces an error.
  ovMean <- overallMean
  ## Create dataframes for genotypes and environments.
  genoDat <- data.frame(x = genoMean, y = scores[, 1] / lam)
  if (!is.null(colorBy)) {
    genoDat <- merge(genoDat, dat[c("genotype", colorBy)],
                     by.x = "row.names", by.y = "genotype")
  }
  envDat <- data.frame(x = envMean, y = loadings[, 1] * lam)
  plotRatio <- (max(c(genoMean, envMean)) - min(c(genoMean, envMean))) /
    (max(c(scores[, 1] / lam, loadings[, 1] * lam)) -
       min(c(scores[, 1] / lam, loadings[, 1] * lam)))
  p <- ggplot2::ggplot(genoDat, ggplot2::aes_string(x = "x", y = "y")) +
    ## Plot genotypes as points.
    ggplot2::geom_point(ggplot2::aes_string(color = colorBy)) +
    ## Needed for a square plot output.
    ggplot2::coord_fixed(ratio = plotRatio, clip = "off") +
    ## Plot environments as texts.
    ggplot2::geom_text(data = envDat,
                       ggplot2::aes_string(x = "x", y = "y",
                                           label = "rownames(envDat)"),
                       size = 3, vjust = 1, color = col[2]) +
    ## Add reference axes.
    ggplot2::geom_vline(ggplot2::aes(xintercept = ovMean),
                        linetype = "dashed", show.legend = FALSE) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = "dashed",
                        show.legend = FALSE) +
    ## Add labeling.
    ggplot2::labs(x = "Main Effects", y = paste0("PC1 (", percPC1, "%)")) +
    ggplot2::ggtitle(paste0(ifelse(GGE, "GGE", "AMMI1"), " biplot for ",
                            trait, " ", year)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

#' Helper function for creating AMMI2 plot
#' @keywords internal
plotAMMI2 <- function(loadings,
                      scores,
                      importance,
                      trait,
                      dat,
                      GGE,
                      year = "",
                      primAxis = "PC1",
                      secAxis = "PC2",
                      scale,
                      col,
                      colorBy) {
  percPC1 <- round(importance[2, primAxis] * 100, 1)
  percPC2 <- round(importance[2, secAxis] * 100, 1)
  if (scale == 1) {
    info <- "environment scaling"
  } else if (scale == 0) {
    info <- "genotype scaling"
  } else if (scale == 0.5) {
    info <- "symmetric scaling"
  } else {
    info <- paste0(round(importance[3, secAxis] * 100, 1), "%")
  }
  ## Calculate lambda scale.
  lam <- as.numeric(importance[1, c(primAxis, secAxis)])
  lam <- lam * sqrt(nrow(scores))
  lam <- lam ^ scale
  ## Create dataframes for genotypes and environments.
  genoDat <- as.data.frame(t(t(scores[, c(primAxis, secAxis)]) / lam))
  if (!is.null(colorBy)) {
    genoDat <- merge(genoDat, dat[c("genotype", colorBy)],
                     by.x = "row.names", by.y = "genotype")
  }
  envDat <- as.data.frame(t(t(loadings[, c(primAxis, secAxis)]) * lam))
  ## Compute multiplication factor for rescaling environmental data.
  mult <- min(
    (max(genoDat[[primAxis]]) - min(genoDat[[primAxis]])) /
      (max(envDat[[primAxis]]) - min(envDat[[primAxis]])),
    (max(genoDat[[secAxis]]) - min(genoDat[[secAxis]])) /
      (max(envDat[[secAxis]]) - min(envDat[[secAxis]]))
  )
  ## Rescale data. 0.6 is more or less random but seems to work well in
  ## practice.
  envDat <- envDat * mult * 0.6
  plotRatio <- (max(c(envDat[[primAxis]], genoDat[[primAxis]])) -
                  min(c(envDat[[primAxis]], genoDat[[primAxis]]))) /
    (max(c(envDat[[secAxis]], genoDat[[secAxis]])) -
       min(c(envDat[[secAxis]], genoDat[[secAxis]])))
  p <- ggplot2::ggplot(genoDat,
                       ggplot2::aes_string(x = primAxis, y = secAxis)) +
    ## Plot genotypes as points.
    ggplot2::geom_point(ggplot2::aes_string(color = colorBy)) +
    ## Needed for a square plot output.
    ggplot2::coord_fixed(clip = "off", ratio = plotRatio) +
    ## Plot environments as texts.
    ggplot2::geom_text(data = envDat,
                       ggplot2::aes_string(x = primAxis, y = secAxis,
                                           label = "rownames(envDat)"),
                       size = 3, vjust = "outward", hjust = "outward",
                       color = col[2]) +
    ## Add arrows from origin to environments.
    ## Adding alpha = for transparency causes the arrows not being plotted
    ## after turning off clipping which is needed since labels may fall off
    ## the plot otherwise.
    ggplot2::geom_segment(data = envDat,
                          ggplot2::aes_string(x = 0, y = 0, xend = primAxis,
                                              yend = secAxis),
                          arrow = ggplot2::arrow(length =
                                                   ggplot2::unit(0.2, "cm")),
                          color = col[2]) +
    ## Add labeling.
    ggplot2::labs(x = paste0(primAxis, " (", percPC1, "%)"),
                  y = paste0(secAxis, " (", percPC2, "%)")) +
    ggplot2::ggtitle(paste0(ifelse(GGE, "GGE", "AMMI2"), " biplot for ", trait,
                            " (", info, ") ", year)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

#' Report method for class AMMI
#'
#' A pdf report will be created containing a summary of an AMMI model.
#' Simultaneously the same report will be created as a tex file.
#'
#' @param x An object of class AMMI.
#' @param ... Further arguments passed on from other functions - not used yet.
#' @param outfile A character string, the name and location of the output .pdf
#' and .tex file for the report. If \code{NULL} a report with a default name
#' will be created in the current working directory.
#'
#' @return A pdf and tex report.
#'
#' @examples
#' ## Run AMMI analysis on TDMaize.
#' geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
#' \dontrun{
#' ## Create a pdf report summarizing the results.
#' report(geAmmi, outfile = "./testReports/reportAmmi.pdf")
#' }
#'
#' @export
report.AMMI <- function(x,
                        ...,
                        outfile = NULL) {
  createReport(x = x, reportName = "ammiReport.Rnw", outfile = outfile, ...)
}



