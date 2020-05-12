#' S3 class varComp
#'
#' Function for creating objects of S3 class varComp.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param fitMod A fitted variance components model.
#' @param modDat A data.frame containing the data used in fitting the model.
#'
#' @seealso \code{\link{plot.varComp}}, \code{\link{report.varComp}}
#'
#' @keywords internal
createVarComp <- function(fitMod,
                          modDat,
                          nesting,
                          useLocYear,
                          fullRandVC,
                          aovFullFixedMod,
                          engine,
                          diagTabs) {
  varComp <- structure(list(fitMod = fitMod,
                            modDat = modDat,
                            nesting = nesting,
                            useLocYear = useLocYear,
                            fullRandVC = fullRandVC,
                            aovFullFixedMod = aovFullFixedMod,
                            engine = engine,
                            diagTabs = diagTabs),
                       class = "varComp")
  attr(varComp, which = "timestamp") <- Sys.time()
  return(varComp)
}

#' @export
print.varComp <- function(x,
                          ...) {
  summary(x$fitmod)
}

#' @export
summary.varComp <- function(object,
                            ...) {
  if (object$engine == "lme4") {
    fitModCall <- deparse(formula(getCall(object$fitMod)), width.cutoff = 500)

  } else if (object$engine == "asreml") {
    fitModCallFixed <- deparse(object$fitMod$call$fixed, width.cutoff = 500)
    fitModCallRandTerms <- attr(x = terms(object$fitMod$call$random),
                                which = "term.labels")
    fitModCallRand <- paste0("(1 | ", fitModCallRandTerms, ")",
                             collapse = " + ")
    fitModCall <- paste(fitModCallFixed, "+", fitModCallRand)
  }
  ## Print source of variation as percentage.
  fullRandVC <- object$fullRandVC
  fullRandVC[["vcov"]] <- sprintf("%1.2f", fullRandVC[["vcov"]])
  fullRandVC[["vcovPerc"]] <- sprintf("%1.2f %%", 100 * fullRandVC[["vcovPerc"]])
  colnames(fullRandVC) <- c("component", "% variance expl.")
  ## Print output
  cat("Fitted model formula\n")
  cat(fitModCall, "\n\n")
  cat("Sources of variation\n")
  print(fullRandVC)
}

#' Plot function for class varComp
#'
#' Plot function for class varComp.
#'
#' @param x An object of class varComp
#' @param plotType A character string. Either "sd" to plot the standard
#' deviation of the variance components, or percVar to plot the percentage of
#' variance explained by each variance component.
#' @param ... Not used
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @export
plot.varComp <- function(x,
                         plotType = c("sd", "percVar"),
                         ...,
                         output = TRUE) {
  plotType <- match.arg(plotType)
  plotVar <- if (plotType == "sd") "sd" else "vcovPerc"
  fullRandVC <- x$fullRandVC
  aovFullFixedMod <- x$aovFullFixedMod
  fullRandVC[["vars"]] <- sapply(X = strsplit(x = rownames(fullRandVC),
                                              split = ":"),
                                 FUN = function(var) {
                                   paste0(sort(var), collapse = "_")
                                 })
  aovFullFixedMod[["vars"]] <- sapply(X = strsplit(x = rownames(aovFullFixedMod),
                                                   split = ":"),
                                      FUN = function(var) {
                                        paste0(sort(var), collapse = "_")
                                      })
  fullRandVC[["Df"]] <- aovFullFixedMod[["Df"]][match(aovFullFixedMod[["vars"]],
                                                      fullRandVC[["vars"]])]
  fullRandVC[["term"]] <- paste(rownames(fullRandVC), "\t\t", fullRandVC[["Df"]])
  fullRandVC[["term"]] <- factor(fullRandVC[["term"]],
                                 levels = rev(fullRandVC[["term"]]))
  fullRandVC[["sd"]] <- sqrt(fullRandVC[["vcov"]])
  annoPosX <- -max(fullRandVC[[plotVar]]) / 5e5
  p <- ggplot(fullRandVC, aes_string(x = plotVar, y = "term")) +
    geom_point(na.rm = TRUE, size = 2) +
    geom_segment(aes_string(xend = plotVar, yend = "term"), x = 0) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(xlim = c(0, NA), clip = "off") +
    theme(panel.background = element_blank(),
          panel.grid.major.x = element_line(color = "grey50"),
          axis.line = element_line(),
          axis.title.y = element_blank(),
          axis.ticks.length.y = grid::unit(0, "mm"),
          axis.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5)) +
    annotation_custom(grid::textGrob("source \t\t df", just = "right",
                                     gp = grid::gpar(size = 14)),
                      xmin = annoPosX, xmax = annoPosX,
                      ymin = Inf, ymax = Inf)
  if (plotType == "sd") {
    p <- p + labs(title = "Standard deviations",
                  x = "Standard deviation estimate")
  } else if (plotType == "percVar") {
    p <- p + labs(title = "Percentage of variance explained",
                  x = "Percentage of variance explained")
  }
  if (output) {
    plot(p)
  }
  invisible(p)
}

#' Report method for class varComp
#'
#' A pdf report will be created containing a summary of an object of class
#' varComp. Simultaneously the same report will be created as a tex
#' file.
#'
#' @inheritParams report.AMMI
#'
#' @param x An object of class varComp.
#'
#' @return A pdf and tex report.
#'
#' @export
report.varComp <- function(x,
                           ...,
                           outfile = NULL) {
  ## Checks.
  if (nchar(Sys.which("pdflatex")) == 0) {
    stop("An installation of LaTeX is required to create a pdf report.\n")
  }
  createReport(x = x, reportName = "varCompReport.Rnw", outfile = outfile,
               reportPackage = "statgenGxE", ...)
}

## @importFrom stats predict

#' @export
predict.varComp <- function(object,
                            ...,
                            predictLevel = c("genotype", "trial", "nesting")) {
  fitMod <- object$fitMod
  modDat <- object$modDat
  predictLevel <- match.arg(predictLevel)
  if (object$useLocYear) {
    envVars <- c("loc", "year")
  } else {
    envVars <- "trial"
  }
  predLevels <- "genotype"
  if (predictLevel == "trial") {
    predLevels <- c(predLevels, envVars)
  } else if (predictLevel == "nesting") {
    predLevels <- c("genotype", object$nesting)
  }
  if (object$engine == "lme4") {
    modDat[["preds"]] <- predict(fitMod)
    preds <- aggregate(x = modDat[["preds"]], by = modDat[predLevels],
                       FUN = mean, na.rm = TRUE)
    colnames(preds)[ncol(preds)] <- "predictedValue"
  } else if (object$engine == "asreml") {
    classForm <- paste0(predLevels, collapse = ":")
    presVars <- union(rownames(attr(terms(update(fitMod$call$fixed, "NULL ~ .")),
                                          "factors")),
                      rownames(attr(terms(fitMod$call$random), "factors")))
    preds <- predictAsreml(model = fitMod, classify = classForm,
                           TD = modDat, present = presVars,
                           vcov = FALSE)$pvals
    preds <- preds[preds[["status"]] == "Estimable", ]
    preds <- as.data.frame(preds[, 1:(ncol(preds) - 2)])
    colnames(preds)[ncol(preds)] <- "predictedValue"
  }
  return(preds)
}

#' Extract variance components
#'
#' Extract variance components from an object of class varComp.
#'
#' @param varComp An object of class varComp.
#'
#' @return A data.frame with variance components and standard errors for
#' the random components in the fitted model.#'
#'
#' @export
vc <- function(varComp) {
  if (varComp$engine == "lme4") {
    varcomps <- as.data.frame(lme4::VarCorr(varComp$fitMod))
    rownames(varcomps) <- varcomps[["grp"]]
    varcomps <- varcomps[c((nrow(varcomps)-1):1, nrow(varcomps)),
                         "vcov", drop = FALSE]
    colnames(varcomps) <- "component"
  } else if (varComp$engine == "asreml") {
    varcomps <- summary(varComp$fitMod)$varcomp
    rownames(varcomps)[nrow(varcomps)] <- "Residual"
    varcomps <- varcomps[, "component", drop = FALSE]
  }
  return(varcomps)
}

#' Compute heritability
#'
#' Compute the heritability based on the fitted model.
#'
#' @param varComp An object of class varComp.
#'
#' @export
herit <- function(varComp) {
  fitMod <- varComp$fitMod
  modDat <- varComp$modDat
  varcomps <- vc(varComp)
  sigmaG <- varcomps["genotype", "component"]
  sigmaRes <- varcomps["Residual", "component"]
  numerator <- sigmaG
  modTerms <- rownames(varcomps)
  if (varComp$engine == "lme4") {
    modVars <- rownames(attr(x = terms(fitMod, random.only = TRUE),
                             which = "factors"))[-c(1, 2)]

  } else if (varComp$engine == "asreml") {
    modVars <- rownames(attr(x = terms(fitMod$call$random), which = "factors"))[-1]
  }
  for (term in modTerms[-c(1, length(modTerms))]) {
    sigmaTerm <- varcomps[term, "component"]
    termVars <- unlist(strsplit(x = term, split = ":"))[-1]
    numerator <- numerator + sigmaTerm /
      prod(sapply(X = termVars, FUN = function(termVar) {
        nlevels(modDat[[termVar]])}))
  }
  if (length(modVars) > 0) {
    numerator <- numerator + sigmaRes /
      prod(sapply(X = modVars, FUN = function(modVar) {
        nlevels(modDat[[modVar]])}))
  } else {
    numerator <- numerator + sigmaRes
  }
  return(sigmaG / numerator)
}

#' Get diagnostics for an object of class varComp
#'
#' Get the diagnostics for the model fitted. This will print a table of
#' combinations missing in the data. For each random factor in the model a
#' table is printed.
#'
#' @param varComp An object of class varComp.
#'
#' @return A list of tables is invisibly returned
#'
#' @export
diagnostics <- function(varComp) {
  if (!inherits(varComp, "varComp")) {
    stop(varComp, "should be an object of class varComp.\n")
  }
  diagTabs <- varComp$diagTabs
  for (diagTab in diagTabs) {
    if (nrow(diagTab) > 0) {
      cat(nrow(diagTab), " missing combinations for ",
          paste(colnames(diagTab), collapse = " x "), ".\n", sep = "")
      print(diagTab, row.names = FALSE)
      cat("\n\n")
    } else {
      cat("No missing combinations for ",
          paste(colnames(diagTab), collapse = " x "), ".\n\n", sep = "")
    }
  }
  invisible(diagTabs)
}
