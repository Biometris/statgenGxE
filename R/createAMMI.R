#' S3 class AMMI
#'
#' Function for creating objects of S3 class AMMI.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param envScores a matrix containing environment scores
#' @param genoScores a matrix containing genotypic scores
#' @param importance a data.frame containing the importance of the principal components
#' @param anova a data.frame containing anova scores of the AMMI analysis
#' @param fitted a matrix containing fitted values from the AMMI model
#' @param trait a character value indicating the analysed trait
#' @param envMean a numerical vector containing the means per environment
#' @param genoMean a numerical vector containing the means per genotype
#' @param overallMean a numerical value containing the overall mean
#' @param x an \code{R} object
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.AMMI}}, \code{\link{report.AMMI}}
#'
#' @name AMMI
NULL

#' @rdname AMMI
#' @export
createAMMI <- function(envScores,
                       genoScores,
                       importance,
                       anova,
                       fitted,
                       trait,
                       envMean,
                       genoMean,
                       overallMean) {
  AMMI <- structure(list(envScores = envScores,
                         genoScores = genoScores,
                         importance = importance,
                         anova = anova,
                         fitted = fitted,
                         trait = trait,
                         envMean = envMean,
                         genoMean = genoMean,
                         overallMean = overallMean),
                    class = "AMMI")
  attr(AMMI, which = "timestamp") <- Sys.time()
  return(AMMI)
}

#' @rdname AMMI
#' @export
is.AMMI <- function(x) {
  inherits(x, "AMMI")
}

#' @export
print.AMMI <- function(x, ...) {
  cat("Principal components",
      "\n====================\n")
  print(x$importance[, 1:ncol(x$envScores)])
  cat("\nAnova",
      "\n=====\n")
  printCoefmat(x$anova)
  cat("\nEnvironment scores",
      "\n==================\n")
  print(x$envScores, ...)
  cat("\nGenotypic scores",
      "\n================\n")
  print(x$genoScores, ..., max.print = 50)
}

#' @export
summary.AMMI <- function(object, ...) {
  print(object, ...)
}

#' Plot Function for Class AMMI
#'
#' Two types of biplot can be made. A biplot of genotype and environment
#' means vs PC1 (AMMI1) or a biplot of genotypes and environment interaction
#' with PC1 and PC2 (AMMI2).
#'
#' @param x An object of class AMMI
#' @param ... Other graphical parameters passed on to actual plot function.
#' @param plotType A character string indicating which plot should be made.
#' Either "AMMI1" or "AMMI2" for an AMMI1 biplot (genotype and
#' environment means vs PC1) or an AMMI2 biplot respectively.
#' @param scale A numeric value. The variables are scaled by
#' \code{lambda ^ scale} and the observations are by \code{lambda ^ (1 - scale)}
#' where \code{lambda} are the singular values computed by
#' \code{\link[stats]{princomp}} in \code{\link{gxeAmmi}}. Normally
#' \code{0 <= scale <= 1}, and a warning will be issued if the specified
#' scale is outside this range.
#' @param col A character vector with plot colors for genotype and environment.
#'
#' @return A biplot depending on \code{plotType}
#'
#' @examples
#' # Run AMMI.
#' geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
#' # Create AMMI2 biplot.
#' plot(geAmmi, plotType = "AMMI2")
#'
#' @import graphics grDevices
#' @export
plot.AMMI <- function(x,
                      ...,
                      plotType = c("AMMI1", "AMMI2"),
                      scale = 1,
                      col = c("orange3", "navyblue")) {
  ## Checks.
  if (length(col) != 2) {
    stop("col should contain exactly two colors.\n")
  }
  if (!is.numeric(scale) || length(scale) > 1) {
    stop("scale should be a single numeric value.\n")
  }
  if (scale < 0 || scale > 1) {
    warning("Scale is outside [0, 1].\n", call. = FALSE)
  }
  plotType <- match.arg(plotType)
  dotArgs <- list(...)
  loadings <- x$envScores
  scores <- x$genoScores
  percPC1 <- round(x$importance[2, 1] * 100, 1)
  percPC2 <- round(x$importance[2, 2] * 100, 1)
  if (plotType == "AMMI1") {
    ## Calculate lambda scale
    lam <- x$importance[1, 1]
    lam <- lam ^ scale
    ## Set arguments for biplot1
    bp1Args <- list(x = 1, type = 'n', main = paste0("AMMI1 biplot for ", x$trait),
                    xlab = "Main Effects",
                    ylab = paste0("PC1 (", percPC1, "%)"),
                    xlim = range(c(x$envMean, x$genoMean)),
                    ylim = range(c(loadings[, 1] * lam, scores[, 1] / lam)))
    ## Add and overwrite args with custom args from ...
    fixedArgs <- c("x", "type", "xlim", "ylim")
    bp1Args <- modifyList(bp1Args, dotArgs[!names(dotArgs) %in% fixedArgs])
    ## Setup plot frame.
    do.call(plot, args = bp1Args)
    ## Add genotypes to empty plot.
    text(x = x$genoMean, y = scores[, 1] / lam, labels = row.names(x$genoMean),
         adj = c(0.5, 0.5), col = col[1])
    ## Add environments to empty plot
    text(x = x$envMean, y = loadings[, 1] * lam, labels = row.names(x$envMean),
         adj = c(0.5, 0.5), col = col[2])
    abline(h = 0, v = x$overallMean, lty = 5)
  } else if (plotType == "AMMI2") {
    if (scale == 1) {
      info <- "environment scaling"
    } else if (scale == 0) {
      info <- "genotype scaling"
    } else if (scale == 0.5) {
      info <- "symmetric scaling"
    } else {
      info <- paste0(round(x$importance[3, 2] * 100, 1), "%")
    }
    ## Calculate lambda scale.
    lam <- as.numeric(x$importance[1, 1:2])
    lam <- lam * sqrt(nrow(scores))
    lam <- lam ^ scale
    ## Set arguments for biplot2.
    bp2Args <- list(x = t(t(scores[, 1:2]) / lam),
                    y = t(t(loadings[, 1:2]) * lam), col = col,
                    cex = c(par("cex") / 2, par("cex")),
                    xlabs = rep("o", nrow(scores)),
                    main = paste0("AMMI2 biplot for ", x$trait, " (", info, ")"),
                    xlab = paste0("PC1 (", percPC1, "%)"),
                    ylab = paste0("PC2 (", percPC2, "%)"),
                    xpd = TRUE)
    ## Add and overwrite args with custom args from ...
    fixedArgs <- c("x", "type", "xlim", "ylim")
    bp2Args <- modifyList(bp2Args, dotArgs[!names(dotArgs) %in% fixedArgs])
    do.call(biplot, args = bp2Args)
  }
}

#' Report method for class AMMI
#'
#' A pdf report will be created containing a summary of AMMI model.
#' Simultaneously the same report will be created as a tex file.
#'
#' @param x An object of class AMMI.
#' @param ... Further arguments passed on from other functions - not used yet.
#' @param outfile A character string, the name and location of the output .pdf and .tex
#' file for the report. If \code{NULL} a report will be created in the current working
#' directory.
#'
#' @export
report.AMMI <- function(x,
                        ...,
                        outfile = NULL) {
  createReport(x = x, reportName = "ammiReport.Rnw",
               outfile = outfile, ...)
}



