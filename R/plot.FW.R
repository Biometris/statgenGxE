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
    #fittedGen, , Y, env, , fVal, orderSens) {
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
      first64 <- levels(x$estimates$G)[1:64]
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

