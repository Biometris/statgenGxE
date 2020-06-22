#' S3 class varComp
#'
#' Function for creating objects of S3 class varComp.\cr
#' \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} methods
#' are available.
#'
#' @param fitMod A fitted variance components model.
#' @param modDat A data.frame containing the data used in fitting the model.
#'
#' @seealso \code{\link{plot.varComp}}
#'
#' @keywords internal
createVarComp <- function(fitMod,
                          modDat,
                          trait,
                          nestingFactor,
                          useLocYear,
                          useRegionLocYear,
                          fullRandVC,
                          aovFullFixedMod,
                          engine,
                          diagTabs) {
  varComp <- structure(list(fitMod = fitMod,
                            modDat = modDat,
                            trait = trait,
                            nestingFactor = nestingFactor,
                            useLocYear = useLocYear,
                            useRegionLocYear = useRegionLocYear,
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
  summary(x)
}

#' @export
summary.varComp <- function(object,
                            ...) {
  if (object$engine == "lme4") {
    ## Display model formula in text form.
    ## This might cut off the formula if it is longer than 500 character.
    ## This is however highly unlikely given the fixed structure of the models.
    fitModCall <- deparse(formula(getCall(object$fitMod)), width.cutoff = 500)

  } else if (object$engine == "asreml") {
    ## For asreml, rewrite the formula to match the display of formulas
    ## from lme4,
    ## No changes needed in the fixed part.
    fitModCallFixed <- deparse(object$fitMod$call$fixed, width.cutoff = 500)
    ## The terms in the random part need to be surrounded by (1 | ...).
    fitModCallRandTerms <- attr(x = terms(object$fitMod$call$random),
                                which = "term.labels")
    fitModCallRand <- paste0("(1 | ", fitModCallRandTerms, ")",
                             collapse = " + ")
    ## Concatenate fixed and random part to form the full formula.
    fitModCall <- paste(fitModCallFixed, "+", fitModCallRand)
  }
  ## Print source of variation as percentage.
  fullRandVC <- object$fullRandVC
  ## Prevent scientific notation in variance component.
  fullRandVC[["vcov"]] <- sprintf("%1.2f", fullRandVC[["vcov"]])
  fullRandVC[["vcovPerc"]] <- sprintf("%1.2f %%", 100 * fullRandVC[["vcovPerc"]])
  colnames(fullRandVC) <- c("component", "% variance expl.")
  ## Print ANOVA for fully fixed model with alternative header.
  aovFullFixedMod <- object$aovFullFixedMod
  attr(x = aovFullFixedMod, which = "heading") <-
    "Analysis of Variance Table for fully fixed model"
  ## Print output.
  cat("Fitted model formula\n")
  cat(fitModCall, "\n\n")
  cat("Sources of variation\n")
  print(fullRandVC)
  cat("\n")
  print(aovFullFixedMod)
}

#' Plot function for class varComp
#'
#' A plot is created of either the standard deviations of each of the terms in
#' the fitted model or the percentage of variance explained by each of the
#' terms in the fitted model. Also the degrees of freedom for each of the
#' terms is shown in the plot.
#'
#' @param x An object of class varComp
#' @param plotType A character string. Either "sd" to plot the standard
#' deviation of the variance components, or "percVar" to plot the percentage of
#' variance explained by each variance component.
#' @param ... Not used
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a ggplot object is invisibly returned.
#'
#' @return A ggplot object is invisibly returned.
#'
#' @examples
#' ## Fit a mixed model.
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")
#'
#' ## Plot the standard deviations.
#' plot(geVarComp)
#' ## Plot the percentage of variance explained.
#' plot(geVarComp, plotType = "percVar")
#'
#' @export
plot.varComp <- function(x,
                         plotType = c("sd", "percVar"),
                         ...,
                         output = TRUE) {
  plotType <- match.arg(plotType)
  ## Extract mu from the fitted model.
  mu <- round(mean(fitted(x$fitMod)))
  ## The actual variable to plot depends on the plotType.
  plotVar <- if (plotType == "sd") "sd" else "vcovPerc"
  ## Extract var comps for random model and anova for fixed model.
  fullRandVC <- x$fullRandVC
  aovFullFixedMod <- x$aovFullFixedMod
  ## Degrees of freedom in the plot come from the anova for the fixed model.
  ## These need to be merged to the random model.
  ## lm (used for the anova) orders interaction terms alphabetically.
  ## To merge these interactions need to be split into their terms and
  ## matching can be done on the terms.
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
  ## term will be used as y-axis label. It consists of term + df
  spaces <- sapply(X = 6 - nchar(fullRandVC[["Df"]]),
                   FUN = function(i) {
                     paste0(rep(" ", times = i), collapse = "")
                   })
  fullRandVC[["term"]] <- paste(rownames(fullRandVC), spaces,
                                fullRandVC[["Df"]])
  ## Revert levels term to get a nice ordering on the y-axis.
  fullRandVC[["term"]] <- factor(fullRandVC[["term"]],
                                 levels = rev(fullRandVC[["term"]]))
  ## Compute standard deviations.
  fullRandVC[["sd"]] <- sqrt(fullRandVC[["vcov"]])
  ## Compute position for annotation. 5e5 is found empirically.
  ## The idea is to annotate the y-axis just left of the x = 0.
  annoPosX <- -max(fullRandVC[[plotVar]]) / 5e5
  p <- ggplot(fullRandVC, aes_string(x = plotVar, y = "term")) +
    geom_point(na.rm = TRUE, size = 2) +
    ## Add line from y-axis to points.
    geom_segment(aes_string(xend = plotVar, yend = "term"), x = 0) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    ## Set lower xlim to 0. This assures 0 is always displayed on the x-axis
    ## even if the lowest variance component is e.g. 1e-8.
    coord_cartesian(xlim = c(0, NA), clip = "off") +
    theme(panel.background = element_blank(),
          panel.grid.major.x = element_line(color = "grey50"),
          axis.line = element_line(),
          axis.title.y = element_blank(),
          axis.ticks.length.y = grid::unit(0, "mm"),
          axis.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5)) +
    annotation_custom(grid::textGrob("Source      df ", just = "right",
                                     gp = grid::gpar(size = 14)),
                      xmin = annoPosX, xmax = annoPosX,
                      ymin = Inf, ymax = Inf)
  if (plotType == "sd") {
    p <- p + labs(title = paste0("Standard deviations (general mean = ",
                                 mu, ")"),
                  x = "Square root of variance estimate")
  } else if (plotType == "percVar") {
    p <- p + labs(title = paste0("Percentage of variance explained ",
                                 "(general mean = ", mu, ")"),
                  x = "Percentage of variance explained")
  }
  if (output) {
    plot(p)
  }
  invisible(p)
}

#' Predictions based on a fitted varComp model.
#'
#' Predictions are made based on the fitted model in the varComp object.
#' These predictions can be at genotype level, at genotype x trial level or at
#' the level of genotype x nestingFactor. If the model was fitted with trial as
#' year x location then genotype x trial level becomes genotype x year x location.
#'
#' @param object An object of class varComp.
#' @param ... Not used.
#' @param predictLevel A character string, the level at which prediction should
#' be made. Either "genotype" for prediction at genotype level, "trial" for
#' predictions at genotype x trial level or the variable used as nesting factor
#' for predictions at the level of genotype x nestingFactor level.
#'
#' @return A data.frame with predictions.
#'
#' @examples
#' ## Fit a mixed model.
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")
#'
#' ## Predictions at genotype level.
#' predGeno <- predict(geVarComp)
#' ## Predictions at genotype x trial level.
#' predGenoTrial <- predict(geVarComp, predictLevel = "trial")
#'
#' @importFrom stats predict
#'
#' @export
predict.varComp <- function(object,
                            ...,
                            predictLevel = c("genotype", "trial",
                                             object$nestingFactor)) {
  predictLevel <- match.arg(predictLevel)
  ## Extract fitted model and model data from object.
  fitMod <- object$fitMod
  modDat <- object$modDat
  ## Variables for environment depend on the fitted model.
  ## Either trial or location x year.
  if (object$useLocYear) {
    envVars <- c("loc", "year")
  } else if (object$useRegionLocYear) {
    envVars <- c("region", "loc", "year")
  } else {
    envVars <- "trial"
  }
  ## Construct vector of levels on which predictions should be made.
  ## Always include genotype.
  predLevels <- "genotype"
  if (predictLevel == "trial") {
    ## For predictLevel trial predict genotype x envVars.
    predLevels <- c(predLevels, envVars)
  } else if (!is.null(object$nestingFactor) &&
             predictLevel == object$nestingFactor) {
    ## For predictLevel nestingFactor predict genotype x nestingFactor variable.
    predLevels <- c(predLevels, object$nestingFactor)
  }
  if (object$engine == "lme4") {
    ## Make predictions for all observations in the data.
    modDat[!is.na(modDat[[object$trait]]), "preds"] <- predict(fitMod)
    ## Compute means per predict level.
    preds <- aggregate(x = modDat[["preds"]], by = modDat[predLevels],
                       FUN = mean, na.rm = TRUE)
    ## Rename column to match asreml output.
    colnames(preds)[ncol(preds)] <- "predictedValue"
  } else if (object$engine == "asreml") {
    ## Construct formula for classify used in predict.
    classForm <- paste0(predLevels, collapse = ":")
    ## Only use observations that where present in the input data for making
    ## predictions. All variables used in the model need to be included here.
    presVars <- union(rownames(attr(terms(update(fitMod$call$fixed, "NULL ~ .")),
                                    "factors")),
                      rownames(attr(terms(fitMod$call$random), "factors")))
    preds <- predictAsreml(model = fitMod, classify = classForm,
                           TD = modDat, present = presVars,
                           vcov = FALSE)$pvals
    ## In some cases NA predictions are returned for Estimable combinations.
    ## Remove those.
    preds <- preds[preds[["status"]] == "Estimable", ]
    ## Convert to data.frame to get rid of attributes.
    ## Remove status column.
    preds <- as.data.frame(preds[, 1:(ncol(preds) - 1)])
    ## Rename column to match lme4 output.
    colnames(preds)[colnames(preds) %in% c("predicted.value", "std.error")] <-
      c("predictedValue", "stdError")
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
#' the random components in the fitted model.
#'
#' @export
vc <- function(varComp) {
  if (!inherits(varComp, "varComp")) {
    stop(varComp, " should be an object of class varComp.\n")
  }
  ## Extract variance component and rename so rows/columns to assure
  ## matching outputs for lme4/asreml.
  if (varComp$engine == "lme4") {
    varcomps <- as.data.frame(lme4::VarCorr(varComp$fitMod))
    rownames(varcomps) <- varcomps[["grp"]]
    rownames(varcomps)[nrow(varcomps)] <- "residual"
    varcomps <- varcomps[c((nrow(varcomps)-1):1, nrow(varcomps)),
                         "vcov", drop = FALSE]
    colnames(varcomps) <- "component"
  } else if (varComp$engine == "asreml") {
    varcomps <- summary(varComp$fitMod)$varcomp
    rownames(varcomps)[nrow(varcomps)] <- "residual"
    varcomps <- varcomps[, "component", drop = FALSE]
  }
  return(varcomps)
}

#' Calculate heritability
#'
#' Calculate the heritability based on the fitted model. The heritability is
#' calculated as described by Atlin et al. E.g. for a model with trials nested
#' within locations, which has a random part that looks like this: genotype +
#' genotype:location + genotype:location:trial the heritability is computed
#' as\cr\cr
#' \deqn{\sigma_G^2 / (\sigma_G^2 + \sigma_L^2 / l + \sigma_{LT}^2 / lt +
#' \sigma_E^2 / ltr)}
#' In this formula the \eqn{\sigma} terms stand for the standard deviations of the
#' respective model terms, and the lower case letters for the number of levels
#' for the respective model terms. So \eqn{\sigma_L} is the standard deviation for
#' the location term in the model and \eqn{l} is the number of locations.
#' \eqn{\sigma_E} corresponds to the residual standard deviation and \eqn{r} to the
#' number of replicates.
#'
#' @param varComp An object of class varComp.
#'
#' @references Atlin, G. N., Baker, R. J., McRae, K. B., & Lu, X. (2000).
#' Selection response in subdivided target regions. Crop Science, 40(1), 7–13.
#' \url{https://doi.org/10.2135/cropsci2000.4017}
#'
#' @export
herit <- function(varComp) {
  if (!inherits(varComp, "varComp")) {
    stop(varComp, " should be an object of class varComp.\n")
  }
  ## Extract fitted model and model data.
  fitMod <- varComp$fitMod
  modDat <- varComp$modDat
  ## Compute variance components.
  varcomps <- vc(varComp)
  ## Extract variance components for genotype and residual.
  sigmaG <- varcomps["genotype", "component"]
  sigmaRes <- varcomps["residual", "component"]
  ## Numerator is constructed by looping over all random model terms and
  ## Adding their share. It always includes sigmaG.
  numerator <- sigmaG
  ## Get the terms used in the random part of the model.
  modTerms <- rownames(varcomps)
  ## Extract all variables used in the random part of the model.
  ## The are needed for computing the contribution of the residual variance.
  if (varComp$engine == "lme4") {
    modVars <- rownames(attr(x = terms(fitMod, random.only = TRUE),
                             which = "factors"))[-c(1, 2)]

  } else if (varComp$engine == "asreml") {
    modVars <- rownames(attr(x = terms(fitMod$call$random),
                             which = "factors"))[-1]
  }
  for (term in modTerms[-c(1, length(modTerms))]) {
    ## Get variance for current term.
    sigmaTerm <- varcomps[term, "component"]
    ## Get variables in current term, exclude gentype (always the first var).
    termVars <- unlist(strsplit(x = term, split = ":"))[-1]
    ## Divide variance by product of #levels for all variables in current term.
    ## Add that to numberator.
    numerator <- numerator + sigmaTerm /
      prod(sapply(X = termVars, FUN = function(termVar) {
        nlevels(modDat[[termVar]])}))
  }
  if (length(modVars) > 0) {
    ## Contribution for residual variance is computed by dividing sigmaRes by
    ## product of #levels of all variables in random part of model.
    numerator <- numerator + sigmaRes /
      prod(sapply(X = modVars, FUN = function(modVar) {
        nlevels(modDat[[modVar]])}))
  } else {
    ## No other variables in random part. Just add sigmaRes to numerator.
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
#' @return A list of tables is invisibly returned.
#'
#' @export
diagnostics <- function(varComp) {
  if (!inherits(varComp, "varComp")) {
    stop(varComp, " should be an object of class varComp.\n")
  }
  ## Get diagTabs from varComp.
  diagTabs <- varComp$diagTabs
  ## For each diagTab print either its content or display a message that
  ## it has no content.
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
