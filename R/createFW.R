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
#' Function for creating scatter, line and trellis plots for objects of class FW.
#'
#' @param x an object of class FW
#' @param ... not unused
#' @param plotType a character vector indicating which plot(s) will be drawn. Possible values
#' "scatter", "line"  and "trellis" for creating a scatter plot of sensitivities, a plot of
#' fitted lines for each genotype and a trellis plot of the individual genotype slopes
#' respectively.
#' @param sortBySens A character string specifying whether the results are to be sorted
#' in an increasing (or decreasing) order of sensitivities.
#' By default, \code{sortBySens = "ascending"}. Other options are "descending" and NA.

#' @return Plots as described in \code{plotType}
#'
#' @import graphics grDevices
#' @export
plot.FW <- function(x,
                    ...,
                    plotType = c("scatter", "line", "trellis"),
                    sortBySens = "ascending") {
  mse <- x$estimates$mse
  genMean <- x$estimates$genMean
  sens <- x$estimates$sens
  envEffs <- x$envEffs$Effect
  fVal <- tapply(X = x$fittedGeno, INDEX = x$data[, c("genotype", "env")], FUN = function(x) {
    mean(x, na.rm = TRUE)
  })
  if ("scatter" %in% plotType) {
    if (!all(is.na(mse))) {
      scatterData <- cbind(genMean, mse, sens)
      colnames(scatterData) <- c("Mean", "m.s.deviation", "Sensitivity")
      pairs(x = scatterData, upper.panel = NULL,
            main = "Finlay & Wilkinson analysis")
    } else {
      plot(x = genMean, y = sens, xlab = "Mean", ylab = "Sensitivity",
           main = "Finlay & Wilkinson analysis")
    }
  }
  if ("line" %in% plotType) {
    minFVal <- min(x$fittedGeno, na.rm = TRUE)
    maxFVal <- max(x$fittedGeno, na.rm = TRUE)
    minXEff <- min(envEffs, na.rm = TRUE)
    maxXEff <- max(envEffs, na.rm = TRUE)
    plot(x = NA, xlim = c(minXEff, maxXEff), ylim = c(minFVal, maxFVal),
         ylab = x$trait, xlab = "Environment", xaxt = "n")
    axis(side = 1, envEffs, levels(x$envEffs$Environment), las = 2, cex.axis = .75)
    color <- 1
    for (i in 1:x$nGeno) {
      if (!is.na(sortBySens)) {
        xfVal <- fVal[names(sens[i]), ]
      } else {
        xfVal <- fVal[i, ]
      }
      lines(envEffs[order(envEffs)], xfVal[order(envEffs)], col = color)
      color <- color + 1
    }
  }
  if ("trellis" %in% plotType) {
    trellisdata <- data.frame(genotype = x$data[["genotype"]], trait = x$data[[x$trait]],
                              fitted = x$fittedGen, xEff = rep(envEffs, x$nGeno))
    if (x$nGeno > 64) {
      first64 <- levels(x$estimates$genotype)[1:64]
      first64 <- x$data[["genotype"]] %in% first64
      trellisdata <- droplevels(trellisdata[first64, ])
    }
    print(lattice::xyplot(trait + fitted ~ xEff | genotype, data = trellisdata,
                          panel = function(x, y, subscripts) {
                            lattice::panel.xyplot(x, y)
                            lattice::panel.lines(trellisdata$xEff[subscripts],
                                                 trellisdata$fitted[subscripts])
                          }, as.table = TRUE, subscripts = TRUE,
                          xlab = "Environment", ylab = x$trait,
                          main = paste0("Finlay & Wilkinson analysis for ", x$trait)))
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


