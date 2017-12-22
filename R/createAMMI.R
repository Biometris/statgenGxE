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
  print(x$importance)
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
#' Function for creating AMMI1 and AMMI2 biplots for objects of class AMMI.
#'
#' @param x an object of class AMMI
#' @param ... not unused
#' @param plotType a character vector indicating which plot(s) will be drawn. Possible values
#' "AMMI1" and "AMMI2" for creating an AMMI1 biplot (genotypes and
#' environments means vs PC1) or an AMMI2 biplot respectively.
#' @param scaleAMMI1 For AMMI1 biplot the variables are scaled by \code{lambda ^ scale} and
#' the observations are by \code{lambda ^ (1 - scale)} where \code{lambda} are the singular values
#' computed by \code{\link[stats]{princomp}} in \code{\link{gxeAmmi}}. Normally
#' \code{0 <= scale <= 1}, and a warning will be issued if the specified scale is
#' outside this range.
#' @param scaleAMMI2 Similar to \code{scaleAMMI1}, for AMMI2 biplot.
#'
#' @return One or two biplots as described in \code{plotType}
#'
#' @import graphics grDevices
#' @export

plot.AMMI <- function(x,
                      ...,
                      plotType = c("AMMI1", "AMMI2"),
                      scaleAMMI1 = 1,
                      scaleAMMI2 = 1) {
  oldPar <- par(xpd = NA)
  loadings <- x$envScores
  scores <- x$genoScores
  propPC1 <- x$importance[2, 1]
  propPC2 <- x$importance[2, 2]
  if ("AMMI1" %in% plotType) {
    if (scaleAMMI1 < 0 || scaleAMMI1 > 1) {
      warning("'scale' is outside [0, 1]")
    }
    # calculate lambda scale
    lam <- x$importance[1, 1]
    if (scaleAMMI1 != 0) {
      lam <- lam ^ scaleAMMI1
    } else {
      lam <- 1
    }
    plot(x = 1, type = 'n', xlim = range(c(x$envMean, x$genoMean)),
         ylim = range(c(loadings[, 1] * lam, scores[, 1] / lam)),
         xlab = "Main Effects", ylab = paste0("PC1 (", round(propPC1 * 100, 1), "%)"),
         main = paste0("AMMI1 biplot for ", x$trait))
    points(x = x$envMean, y = loadings[, 1] * lam, type = "n", col = "navyblue", lwd = 5)
    text(x = x$envMean, y = loadings[, 1] * lam, labels = row.names(x$envMean),
         adj = c(0.5, 0.5), col = "navyblue")
    points(x = x$genoMean, y = scores[, 1] / lam, type = "n", col = "orange3", lwd = 5)
    text(x = x$genoMean, y = scores[, 1] / lam, labels = row.names(x$genoMean),
         adj = c(0.5, 0.5), col = "orange3")
    abline(h = 0, v = x$overallMean, lty = 5)
  }
  if ("AMMI2" %in% plotType) {
    if (scaleAMMI2 < 0 || scaleAMMI2 > 1) {
      warning("'scale' is outside [0, 1]")
    }
    if (scaleAMMI2 == 1) {
      info <- "environment scaling"
    } else if (scaleAMMI2 == 0) {
      info <- "genotype scaling"
    } else if (scaleAMMI2 == 0.5) {
      info <- "symmetric scaling"
    } else {
      info <- paste0(round(x$importance[3, 2] * 100, 1), "%")
    }
    # calculate lambda scale
    lam <- as.numeric(x$importance[1, 1:2])
    lam <- lam * sqrt(nrow(scores))
    if (scaleAMMI2 != 0) {
      lam <- lam ^ scaleAMMI2
    } else {
      lam <- 1
    }
    biplot(x = t(t(scores[, 1:2]) / lam),
           y = t(t(loadings[, 1:2]) * lam),
           col = c("orange3", "navyblue"),
           cex = c(par("cex") / 2, par("cex")),
           xlabs = rep("o", nrow(scores)),
           main = paste0("AMMI2 biplot for ", x$trait, " (", info, ")"),
           xlab = paste0("PC1 (", round(propPC1 * 100, 1), "%)"),
           ylab = paste0("PC2 (", round(propPC2 * 100, 1), "%)"))
  }
  par <- oldPar
}

#' Report method for class AMMi
#'
#' A pdf report will be created containing a summary of AMMI model.
#' Simultaneously the same report will be created as a tex file.
#'
#' @param x an object of class AMMI.
#' @param ... further arguments passed on from other functions - not used yet.
#' @param outfile a character string, the name and location of the output .pdf and .tex
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



