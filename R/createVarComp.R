#' S3 class varComp
#'
#' Function for creating objects of S3 class varComp.\cr
#' \code{\link{print}}, \code{\link{summary}}, and \code{\link{plot}} methods
#' are available.
#'
#' @param fitMod A fitted variance components model.
#' @param modDat A data.frame containing the data used in fitting the model.
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
  engine <- object$engine
  trait <- object$trait
  if (engine == "lme4") {
    ## Display model formula in text form.
    ## This might cut off the formula if it is longer than 500 character.
    ## This is however highly unlikely given the fixed structure of the models.
    fitModCall <- deparse(formula(getCall(object$fitMod)), width.cutoff = 500)
  } else if (engine == "asreml") {
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
  if (engine == "asreml") {
    fullRandVC[["stdError"]] <- sprintf("%1.3f", fullRandVC[["stdError"]])
  }
  fullRandVC[["vcovPerc"]] <- sprintf("%1.2f %%", 100 * fullRandVC[["vcovPerc"]])
  colnames(fullRandVC) <- c("Component", if (engine == "asreml") "SE",
                            "% Variance expl.")
  aovFullFixedMod <- object$aovFullFixedMod
  ## Construct fully fixed and fully random model formula from ANOVA.
  modTerms <- rownames(aovFullFixedMod)[-nrow(aovFullFixedMod)]
  fixModCall <- paste(trait, "~", paste0(modTerms, collapse = " + "))
  randModCall <- paste(trait, "~",
                       paste0("(1 | ", modTerms, ")", collapse = " + "))
  ## Print ANOVA for fully fixed model with alternative header.
  attr(x = aovFullFixedMod, which = "heading") <-
    paste("Analysis of Variance Table for fully fixed model:\n",
          fixModCall, "\n")
  ## Print output.
  cat("Fitted model formula final mixed model\n\n")
  cat("", fitModCall, "\n\n")
  cat("Sources of variation for fully random model:\n")
  cat("", randModCall, "\n\n")
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
#' @param ... Not used.
#' @param plotType A character string. Either "sd" to plot the standard
#' deviation of the variance components, or "percVar" to plot the percentage of
#' variance explained by each variance component.
#' @param title A character string used a title for the plot.
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
#' @family Mixed model analysis
#'
#' @export
plot.varComp <- function(x,
                         ...,
                         plotType = c("sd", "percVar"),
                         title = NULL,
                         output = TRUE) {
  plotType <- match.arg(plotType)
  chkChar(title, len = 1)
  ## Extract mu from the fitted model.
  mu <- round(mean(fitted(x$fitMod)))
  if (is.null(title)) {
    title <- paste0(ifelse(plotType == "sd", "Standard deviations",
                           "Percentage of variance explained"),
                    " (general mean = ", mu, ")")
  }
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
  p <- ggplot2::ggplot(fullRandVC,
                       ggplot2::aes_string(x = plotVar, y = "term")) +
    ggplot2::geom_point(na.rm = TRUE, size = 2) +
    ## Add line from y-axis to points.
    ggplot2::geom_segment(ggplot2::aes_string(xend = plotVar, yend = "term"),
                          x = 0) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ## Set lower xlim to 0. This assures 0 is always displayed on the x-axis
    ## even if the lowest variance component is e.g. 1e-8.
    ggplot2::coord_cartesian(xlim = c(0, NA), clip = "off") +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(color = "grey50"),
                   axis.line = ggplot2::element_line(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.ticks.length.y = grid::unit(0, "mm"),
                   axis.text = ggplot2::element_text(size = 12),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::annotation_custom(grid::textGrob("Source      df ", just = "right",
                                              gp = grid::gpar(size = 14)),
                               xmin = annoPosX, xmax = annoPosX,
                               ymin = Inf, ymax = Inf) +
    ggplot2::ggtitle(title)
  if (plotType == "sd") {
    p <- p + ggplot2::labs(x = "Square root of variance estimate")
  } else if (plotType == "percVar") {
    p <- p + ggplot2::labs(x = "Percentage of variance explained")
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
#' predictions at genotype x trial level, the variable used as nesting factor
#' for predictions at the level of genotype x nestingFactor level, or one or
#' more of the extra terms used in the model. E.g. c("region", "year") for a
#' model fitted with \code{regionLocationYear = TRUE}.
#'
#' @return A data.frame with predictions.
#'
#' @examples
#' ## Fit a mixed model.
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")
#'
#' ## Predictions at genotype level.
#' predGeno <- predict(geVarComp)
#' head(predGeno)
#'
#' ## Predictions at genotype x trial level.
#' predGenoTrial <- predict(geVarComp, predictLevel = "trial")
#' head(predGenoTrial)
#'
#' @importFrom stats predict
#'
#' @family Mixed model analysis
#'
#' @export
predict.varComp <- function(object,
                            ...,
                            predictLevel = "genotype") {
  ## Extract fitted model and model data from object.
  fitMod <- object$fitMod
  modDat <- object$modDat
  ## Variables for environment depend on the fitted model.
  ## Either trial or location x year.
  if (object$useLocYear) {
    predVars <- c("genotype", "trial", "loc", "year")
  } else if (object$useRegionLocYear) {
    predVars <- c("genotype", "trial", "region", "loc", "year")
  } else {
    predVars <- c("genotype", "trial", object$nestingFactor)
  }
  predLevels <- match.arg(predictLevel, choices = predVars, several.ok = TRUE)
  ## Always include genotype.
  predLevels <- unique(c("genotype", predLevels))
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
#' @examples
#' ## Fit a mixed model.
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")
#'
#' ## Extract variance components.
#' vc(geVarComp)
#'
#' @family Mixed model analysis
#'
#' @export
vc <- function(varComp) {
  if (!inherits(varComp, "varComp")) {
    stop(varComp, " should be an object of class varComp.\n")
  }
  fitMod <- varComp$fitMod
  ## Extract variance component and rename so rows/columns to assure
  ## matching outputs for lme4/asreml.
  if (varComp$engine == "lme4") {
    modTerms <- colnames(attr(x = terms(getCall(fitMod)$formula,
                                        keep.order = TRUE), which = "factors"))
    modTerms <- gsub(pattern = "1 | ", replacement = "", x = modTerms,
                     fixed = TRUE)
    varcomps <- as.data.frame(lme4::VarCorr(fitMod))
    rownames(varcomps) <- varcomps[["grp"]]
    rownames(varcomps)[nrow(varcomps)] <- "residuals"
    modTermsRand <- modTerms[modTerms %in% rownames(varcomps)]
    varcomps <- varcomps[c(modTermsRand, "residuals"), "vcov", drop = FALSE]
    colnames(varcomps) <- "Component"
  } else if (varComp$engine == "asreml") {
    modTerms <- colnames(attr(x = terms(fitMod$call$random, keep.order = TRUE),
                              which = "factors"))
    varcomps <- summary(fitMod)$varcomp
    rownames(varcomps)[nrow(varcomps)] <- "residuals"
    varcomps <- varcomps[c(modTerms, "residuals"), c("component", "std.error")]
    colnames(varcomps) <- c("Component", "SE")
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
#' In this formula the \eqn{\sigma} terms stand for the standard deviations of
#' the respective model terms, and the lower case letters for the number of
#' levels for the respective model terms. So \eqn{\sigma_L} is the standard
#' deviation for the location term in the model and \eqn{l} is the number of
#' locations. \eqn{\sigma_E} corresponds to the residual standard deviation and
#' \eqn{r} to the number of replicates.
#'
#' @param varComp An object of class varComp.
#'
#' @examples
#' ## Fit a mixed model.
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")
#'
#' ## Compute heritability.
#' herit(geVarComp)
#'
#' @family Mixed model analysis
#'
#' @references Atlin, G. N., Baker, R. J., McRae, K. B., & Lu, X. (2000).
#' Selection response in subdivided target regions. Crop Science, 40(1), 7–13.
#' \doi{10.2135/cropsci2000.4017}
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
  sigmaG <- varcomps["genotype", "Component"]
  sigmaRes <- varcomps["residual", "Component"]
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
  ## Get median number for times genotypes are tested within modVars.
  nLevModVars <- sapply(X = modVars, FUN = function(modVar) {
    median(rowSums(table(modDat[["genotype"]], modDat[[modVar]]) > 0))
  })
  for (term in modTerms[-c(1, length(modTerms))]) {
    ## Get variance for current term.
    sigmaTerm <- varcomps[term, "Component"]
    ## Get variables in current term, exclude genotype (always the first var).
    termVars <- unlist(strsplit(x = term, split = ":"))[-1]
    ## Divide variance by product of #levels for all variables in current term.
    ## Add that to numerator.
    numerator <- numerator + sigmaTerm / prod(nLevModVars[termVars])
  }
  nReps <- median(table(modDat[["genotype"]], modDat[["trial"]]))
  if (length(modVars) > 0) {
    ## Contribution for residual variance is computed by dividing sigmaRes by
    ## product of #levels of all variables in random part of model and
    ## #replicates.
    numerator <- numerator + sigmaRes / prod(nLevModVars, nReps)
  } else {
    ## No other variables in random part.
    ## Just divide sigmaRes by #replicates.
    numerator <- numerator + sigmaRes / nReps
  }
  return(sigmaG / numerator)
}

#' Calculate the correlated response to selection
#'
#' Calculate the correlated response to selection (CRDR) based on the fitted
#' model. The CRDR is calculated as described by Atlin et al. E.g. for a model
#' with trials nested within scenarios, which has a random part that looks like
#' this: genotype + genotype:scenario + genotype:scenario:trial the CRDR is
#' calculated as:\cr\cr
#' \deqn{H1 = \sigma_G^2 / (\sigma_G^2 + \sigma_S^2 / s + \sigma_{ST}^2 / st +
#' \sigma_E^2 / str)}
#' \deqn{H2 = (\sigma_G^2 + \sigma_S^2) / (\sigma_G^2 + \sigma_S^2 +
#' \sigma_{ST}^2 / st + \sigma_E^2 / str)}
#' \deqn{CRDR = (\sigma_G^2 / (\sigma_G^2 + \sigma_S^2)) * sqrt(H1 / H2)}
#' In these formulas the \eqn{\sigma} terms stand for the standard deviations of
#' the respective model terms, and the lower case letters for the number of
#' levels for the respective model terms. So \eqn{\sigma_S} is the standard
#' deviation for the scenario term in the model and \eqn{s} is the number of
#' scenarios. \eqn{\sigma_E} corresponds to the residual standard deviation and
#' \eqn{r} to the number of replicates.
#'
#' @inheritParams herit
#'
#' @references Atlin, G. N., Baker, R. J., McRae, K. B., & Lu, X. (2000).
#' Selection response in subdivided target regions. Crop Science, 40(1), 7–13.
#' \doi{10.2135/cropsci2000.4017}
#'
#' @family Mixed model analysis
#'
#' @export
CRDR <- function(varComp) {
  if (!inherits(varComp, "varComp")) {
    stop(varComp, " should be an object of class varComp.\n")
  }
  if (is.null(varComp$nestingFactor) && isFALSE(varComp$useRegionLocYear)) {
    stop("CRDR can only be computed when a model is fitted with a nesting ",
         "structure or when regions are included in the model.\n")
  }
  H1 <- herit(varComp)
  ## Extract fitted model and model data.
  fitMod <- varComp$fitMod
  modDat <- varComp$modDat
  ## Get factor for computing H2
  if (!is.null(varComp$nestingFactor)) {
    H2factor <- paste0("genotype:", varComp$nestingFactor)
  } else if (varComp$useRegionLocYear) {
    H2factor <- "genotype:region"
  }
  ## Compute variance components.
  varcomps <- vc(varComp)
  ## Extract variance components for genotype, H2factor and residual.
  sigmaG <- varcomps["genotype", "component"]
  sigmaH2 <- varcomps[H2factor, "component"]
  sigmaRes <- varcomps["residual", "component"]
  ## Numerator is constructed by looping over all random model terms and
  ## Adding their share. It always includes sigmaG and sigmaH2.
  numerator <- sigmaG + sigmaH2
  ## Get the terms used in the random part of the model.
  modTerms <- rownames(varcomps)
  ## Extract all variables used in the random part of the model.
  ## They are needed for computing the contribution of the residual variance.
  if (varComp$engine == "lme4") {
    modVars <- rownames(attr(x = terms(fitMod, random.only = TRUE),
                             which = "factors"))[-c(1, 2)]

  } else if (varComp$engine == "asreml") {
    modVars <- rownames(attr(x = terms(fitMod$call$random),
                             which = "factors"))[-1]
  }
  ## Get median number for times genotypes are tested within modVars.
  nLevModVars <- sapply(X = modVars, FUN = function(modVar) {
    median(rowSums(table(modDat[["genotype"]], modDat[[modVar]]) > 0))
  })
  H2factorPos <- which(modTerms == H2factor)
  for (term in modTerms[-c(1, H2factorPos, length(modTerms))]) {
    ## Get variance for current term.
    sigmaTerm <- varcomps[term, "component"]
    ## Get variables in current term, exclude genotype (always the first var).
    termVars <- unlist(strsplit(x = term, split = ":"))[-1]
    ## Divide variance by product of #levels for all variables in current term.
    ## Add that to numerator.
    numerator <- numerator + sigmaTerm / prod(nLevModVars[termVars])
  }
  nReps <- median(table(modDat[["genotype"]], modDat[["trial"]]))
  if (length(modVars) > 0) {
    ## Contribution for residual variance is computed by dividing sigmaRes by
    ## product of #levels of all variables in random part of model and
    ## #replicates.
    numerator <- numerator + sigmaRes / prod(nLevModVars, nReps)
  } else {
    ## No other variables in random part.
    ## Just divide sigmaRes by #replicates.
    numerator <- numerator + sigmaRes / nReps
  }
  H2 <- (sigmaG + sigmaH2) / numerator
  r <- sigmaG / (sigmaG + sigmaH2)
  return(r * sqrt(H1 / H2))
}

#' Compute different types of correlations.
#'
#' Compute three types of correlations for models fitted with a nesting factor.
#' \itemize{
#' \item{correlation between scenarios or environment types:
#' \deqn{\sigma_G^2 / (\sigma_G^2 + \sigma_{GS}^2)}
#' }
#' \item{correlation between trials within scenarios or environment types:
#' \deqn{(\sigma_G^2 + \sigma_{GS}^2) / (\sigma_G^2 + \sigma_{GS}^2 +
#' \sigma_E^2)}
#' }
#' \item{correlation between trials that belong to different
#' scenarios/environment types:
#' \deqn{\sigma_G^2 / (\sigma_G^2 + \sigma_{GS}^2 + \sigma_E^2)}
#' }
#' }
#' In these formulas the \eqn{\sigma} terms stand for the standard deviations of
#' the respective model terms. So \eqn{\sigma_S} is the standard deviation for
#' the scenario term in the model, \eqn{\sigma_{GS}} for the standard deviation
#' of the genotype by scenario term and \eqn{\sigma_E} corresponds to the
#' residual standard deviation.
#'
#' @inheritParams herit
#'
#' @return A list with three correlations.
#'
#' @family Mixed model analysis
#'
#' @export
correlations <- function(varComp) {
  if (!inherits(varComp, "varComp")) {
    stop(varComp, " should be an object of class varComp.\n")
  }
  ## Get factor for computing correlations.
  if (!is.null(varComp$nestingFactor)) {
    corFactor <- paste0("genotype:", varComp$nestingFactor)
  } else {
    stop("correlations can only be computed when a model is fitted with a ",
         "nesting structure.\n")
  }
  ## Compute variance components.
  varComps <- vc(varComp)
  varGeno <- varComps["genotype", "component"]
  varCorFactor <- varComps[corFactor, "component"]
  varRes <- varComps["residuals", "component"]
  ## Compute correlation between scenarios.
  rScen <- varGeno / (varGeno + varCorFactor)
  ## Compute correlation between trials within scenarios.
  rTrScen <- (varGeno + varCorFactor) / (varGeno + varCorFactor + varRes)
  ## Compute correlattion between trials that belong to different scenarios.
  rTrDiffScen <- varGeno / (varGeno + varCorFactor + varRes)
  return(list(rScen = rScen, rTrScen = rTrScen, rTrDiffScen = rTrDiffScen))
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
#' @examples
#' ## Fit a mixed model.
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")
#'
#' ## Display diagnostics.
#' diagnostics(geVarComp)
#'
#' @family Mixed model analysis
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
