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
                          trialGroup,
                          useLocYear,
                          fullRandVC,
                          engine,
                          diagTabs) {
  varComp <- structure(list(fitMod = fitMod,
                            modDat = modDat,
                            trialGroup = trialGroup,
                            useLocYear = useLocYear,
                            fullRandVC = fullRandVC,
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
  object$fullRandVC[["vcovPerc"]] <-
    sprintf("%1.2f %%", 100 * object$fullRandVC[["vcovPerc"]])
  ## Print output
  cat("Fitted model formula\n")
  cat(fitModCall, "\n\n")
  cat("Sources of variation\n")
  print(setNames(object$fullRandVC[, "vcovPerc", drop = FALSE], NULL))
}

#' Plot function for class varComp
#'
#' Plot function for class varComp.
#'
#' @param x An object of class varComp
#' @param ... Not used
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @export
plot.varComp <- function(x,
                         ...,
                         output = TRUE) {
  fullRandVC <- x$fullRandVC
  fullRandVC[["term"]] <- factor(rownames(fullRandVC),
                                 levels = rev(rownames(fullRandVC)))
  p <- ggplot(fullRandVC, aes_string(x = "vcov", y = "term")) +
    geom_point()
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
                            predictLevel = c("genotype", "trial", "trialGroup")) {
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
  } else if (predictLevel == "trialGroup") {
    predLevels <- c("genotype", object$trialGroup)
  }
  if (object$engine == "lme4") {
    gridLevels <- sapply(X = predLevels, FUN = function(predLevel) {
      levels(modDat[[predLevel]])
    }, simplify = FALSE)
    newDat <- do.call(expand.grid, gridLevels)
    predicted.value <- predict(fitMod, newdata = newDat)
    preds <- cbind(newDat, predicted.value)
  } else if (object$engine == "asreml") {
    classForm <- paste0(predLevels, collapse = ":")
    preds <- predictAsreml(model = fitMod, classify = classForm,
                           TD = object$modDat,
                           aliased = TRUE, vcov = FALSE)$pvals
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
