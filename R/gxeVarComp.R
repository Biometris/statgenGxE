#' Selects the best variance-covariance model for a set of trials
#'
#' This function selects the best covariance structure for genetic correlations
#' between trials. It fits a range of variance-covariance models (identity,
#' compound symmetry (cs), diagonal, simple correlation with heterogeneous
#' variance (outside), heterogeneous compound symmetry (hcs),
#' first order factor analytic (fa), second order factor analytic (fa2) and
#' unstructured), and selects the best one using a goodness-of-fit criterion.
#'
#' @inheritParams gxeAmmi
#'
#' @param engine A character string specifying the engine used for modeling.
#' Either "lme4" or "asreml".
#' @param trialGroup A character string specifying a column in TD.......
#' @param ... Further arguments to be passed to \code{asreml}.
#'
#' @note If \code{engine = "lme4"}, only the compound symmetry model can be
#' fitted.
#'
#' @return An object of class \code{\link{varComp}}, a list object containing:
#' \item{STA}{An object of class STA containing the best fitted model.}
#' \item{choice}{A character string indicating the best fitted model.}
#' \item{summary}{A data.frame with a summary of the fitted models.}
#' \item{engine}{A character string containing the engine used for
#' the analysis.}
#'
#' @examples
#' ## Select the best variance-covariance model using lme4 for modeling.
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")
#'
#' ## Summarize results.
#' summary(geVarComp)
#'
#' \dontrun{
#' ## Create a pdf report summarizing the results.
#' report(geVarComp, outfile = "./testReports/reportVarComp.pdf")
#' }
#'
#' \dontrun{
#' ## Select the best variance-covariance model using asreml for modeling.
#' ## Use BIC as a goodness-of-fit criterion.
#' geVarComp2 <- gxeVarComp(TD = TDMaize, trait = "yld", engine = "asreml",
#'                         criterion = "BIC")
#'
#' summary(geVarComp2)
#'
#' ## Plot a heatmap of the correlation matrix for the best model.
#' plot(geVarComp2)
#' }
#' @export
gxeVarComp <- function(TD,
                       trials = names(TD),
                       trait,
                       engine = c("lme4", "asreml"),
                       #trialGroup = NULL,
                       locationYear = FALSE,
                       nesting = NULL,
                       regionLocationYear = FALSE,
                       useWt = FALSE,
                       ...) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  trials <- chkTrials(trials, TD)
  TDTot <- Reduce(f = rbind, x = TD[trials])
  chkCol(trait, TDTot)
  chkCol("trial", TDTot)
  chkCol("genotype", TDTot)
  if (locationYear) {
    chkCol("loc", TDTot)
    chkCol("year", TDTot)
    ## Check whether year x locations are crossed or nested.
    locYearTab <- table(TDTot[["year"]], TDTot[["loc"]])
    locYearCrossed <- all(locYearTab > 0)
  }
  if (!is.null(nesting)) {
    chkCol(nesting, TDTot)
  }
  if (regionLocationYear) {
    chkCol("loc", TDTot)
    chkCol("year", TDTot)
    chkCol("region", TDTot)
  }
  if (useWt) {
    chkCol("wt", TDTot)
  } else {
    TDTot[["wt"]] <- 1
  }
  engine <- match.arg(engine)
  TDTot <- droplevels(TDTot)
  #useLocYear <- hasName(TDTot, "year") & hasName(TDTot, "loc")
  ## Increase maximum number of iterations for asreml.
  maxIter <- 200
  ## Add combinations of trial and genotype currently not in TD to TD.
  ## No missing combinations are allowed when fitting asreml models.
  ## Asreml can't handle missing weights, so set them to 0 when missing.
  TDTot[is.na(TDTot[["wt"]]), "wt"] <- 0
  ## Check if the trial is nested within the trialGroup.
  #hasGroup <- !is.null(trialGroup)
  ## Check if the data contains replicates.
  repTab <- table(TDTot[["trial"]], TDTot[["genotype"]])
  hasReps <- any(repTab > 1)
  ## Construct formula for fixed part - first as text.
  ## Trying to fit this in something 'smart' actually makes it unreadable.
  ## First create a vector with the separate terms.
  ## This avoids difficult constructions to get the +-es correct.
  # fixedTerms <- c(if (!useLocYear && !hasGroup) "trial",
  #                 if (hasGroup) trialGroup,
  #                 if (!useLocYear && hasGroup) paste0(trialGroup, ":trial"),
  #                 if (useLocYear) "year",
  #                 if (useLocYear && !hasGroup) "loc",
  #                 if (useLocYear && !hasGroup) "loc:year",
  #                 if (useLocYear && hasGroup) paste0(trialGroup, ":year"),
  #                 if (useLocYear && hasGroup) paste0(trialGroup, ":loc:year"))

  fixedTerms <- c(if (!locationYear && is.null(nesting) &&
                      !regionLocationYear) "trial",
                  if (locationYear) c("year", if (locYearCrossed) "loc",
                                      "year:loc"),
                  if (!is.null(nesting)) c(nesting, paste0(nesting, ":trial")),
                  if (regionLocationYear) c("region", "region:loc", "year",
                                            "region:year", "region:loc",
                                            "region:loc:year"))

  ## Construct formula for random part in a similar way.
  # randTerms <- c("genotype",
  #                if (!useLocYear && !hasGroup && (hasReps || useWt)) "genotype:trial",
  #                if (hasGroup) paste0("genotype:", trialGroup),
  #                if (!useLocYear && hasGroup && (hasReps || useWt))
  #                  paste0("genotype:", trialGroup, ":trial"),
  #                if (useLocYear) "genotype:year",
  #                if (useLocYear && !hasGroup) "genotype:loc",
  #                if (useLocYear && !hasGroup && (hasReps || useWt)) "genotype:loc:year",
  #                if (useLocYear && hasGroup) paste0("genotype:", trialGroup, ":year"),
  #                if (useLocYear && hasGroup  && (hasReps || useWt))
  #                  paste0("genotype:", trialGroup, ":loc:year"))

  if (hasReps || useWt) {
    randTermIncl <- fixedTerms
  } else {
    randTermIncl <- fixedTerms[-length(fixedTerms)]
  }
  randTerms <- c("genotype",
                 if (length(randTermIncl) > 0) paste0("genotype:", randTermIncl))

  ## First fit a model with all terms fixed to determine:
  ## - should all terms in the fixed part really be present.
  ## - Predict which terms in the random part of the model will probably
  ##   have a zero variance component.
  fullFixedTxt <- paste0("`", trait, "`~",
                         paste(c(fixedTerms, randTerms), collapse = "+"))

  ## Construct input for full random model.
  if (engine == "lme4") {
    fullRandTxt <- paste0("`", trait, "`~",
                          paste(paste0("(1|", c(fixedTerms, randTerms), ")"),
                                collapse = "+"))
  } else if (engine == "asreml") {
    fullRandTxt <- paste("~", paste(c(fixedTerms, randTerms), collapse = "+"))
  }


  fullFixedMod <- lm(formula(fullFixedTxt), data = TDTot)
  aovFullFixedMod <- anova(fullFixedMod)
  ## Get all model terms as used by lm (might involve reordered terms).
  fullFixedLabs <- attr(x = terms(fullFixedMod), which = "term.labels")
  if (!all(fullFixedLabs %in% rownames(aovFullFixedMod))) {
    ## At least one terms missing from anova.
    ## If this is a fixed term remove it from fixed.
    missTerms <- fullFixedLabs[!fullFixedLabs %in% rownames(aovFullFixedMod)]
    for (missTerm in missTerms) {
      fixedTermSets <- strsplit(x = fixedTerms, split = ":")
      missTermSet <- unlist(strsplit(x = missTerm, split = ":"))
      remPos <- !sapply(X = fixedTermSets, FUN = setequal, missTermSet)
      fixedTerms <- fixedTerms[remPos]
    }
  }
  ## Get rand terms that indicate zero variance components.
  ## lm reorders variables in model terms.
  ## use sets of variables in terms to compare them.
  aovTermSets <- strsplit(x = rownames(aovFullFixedMod), split = ":")
  for (randTerm in randTerms) {
    ## Convert term to set of variables in term.
    randTermSet <- unlist(strsplit(x = randTerm, split = ":"))
    ## Get position of term in anova table by comparing sets.
    randTermPos <- sapply(X = aovTermSets, FUN = setequal, randTermSet)
    ## Get MSS for current term.
    MSSRandTerm <- aovFullFixedMod[randTermPos, "Mean Sq"]
    ## For all other terms in the anova table that have the current term
    ## as a subset the MSS cannot be higher.
    ## If it is the corresponding variance component is possibly zero.
    for (i in seq_along(aovTermSets)) {
      if ((all(randTermSet %in% aovTermSets[i]) ||
           ## Always include the residual term for comparison.
           i == nrow(aovFullFixedMod)) &&
          aovFullFixedMod[i, "Mean Sq"] > MSSRandTerm) {
        warning("Mean Sum of Squares for ", randTerm, " smaller than Mean ",
                "Sum of Squares for ", rownames(aovFullFixedMod)[i], ".\n",
                "Possible zero variance components.\n", call. = FALSE)
      }
    }
  }

  if (engine == "lme4") {
    fullRandMod <- lme4::lmer(formula(fullRandTxt), data = TDTot)
  } else if (engine == "asreml") {
    fullRandMod <- tryCatchExt(asreml::asreml(fixed = formula(paste0("`", trait, "`~ 1")),
                                  random = formula(fullRandTxt), data = TDTot,
                                  trace = TRUE))
    if (!is.null(fullRandMod$warning)) {
      ## Check if param 1% increase is significant. Remove warning if not.
      fullRandMod <- chkLastIter(fullRandMod)
    }
    if (length(fullRandMod$warning) != 0) {
      warning(fullRandMod$warning, "\n", call. = FALSE)
    }
    if (!is.null(fullRandMod$error)) {
      stop(fullRandMod$error)
    } else {
      fullRandMod <- fullRandMod$value
    }
  }


  ## Create the full fixed part of the model as a character.
  ## This is identical for lme4 and asreml so only needs to be done once.
  fixedTxt <- paste0("`", trait, "`~", paste(fixedTerms, collapse = "+"))
  if (engine == "lme4") {
    randTxt <- paste(paste0("(1|", randTerms, ")"), collapse = "+")
    formTxt <- paste(fixedTxt, "+", randTxt)
    ## Fit the actual model.
    mr <- lme4::lmer(formula(formTxt), data = TDTot, weights = TDTot[["wt"]],
                     ...)
    ## Construct STA object.
  } else if (engine == "asreml") {
    if (requireNamespace("asreml", quietly = TRUE)) {
      randTxt <- paste("~ ", paste(randTerms, collapse = "+"))
      ## Put arguments for models in a list to make it easier to switch
      ## between asreml3 and asreml4. Usually only one or two arguments differ.
      ## Also some arguments are identical for all models
      modArgs0 <- list(fixed = formula(fixedTxt), random = formula(randTxt),
                       data = TDTot, weights = "wt", maxiter = maxIter,
                       trace = TRUE)
      modArgs <- modArgs0
      ## Fit the actual model.
      mr <- tryCatchExt(do.call(asreml::asreml, modArgs))
      if (!is.null(mr$warning)) {
        ## Check if param 1% increase is significant. Remove warning if not.
        mr <- chkLastIter(mr)
      }
      if (length(mr$warning) != 0) {
        warning("Asreml gave the following warning:\n", mr$warning, "\n",
                call. = FALSE)
      }
      if (!is.null(mr$error)) {
        warning("Asreml gave the following error:\n", mr$error, call. = FALSE)
        mr <- list(loglik = -Inf)
      } else {
        mr <- mr$value
        mr$call$fixed <- eval(mr$call$fixed)
        mr$call$random <- eval(mr$call$random)
        mr$call$rcov <- eval(mr$call$rcov)
        mr$call$G.param <- eval(mr$call$G.param)
        mr$call$R.param <- eval(mr$call$R.param)
        if (!mr$converge) {
          warning("No convergence.\n", call. = FALSE)
          mr$loglik <- -Inf
        }
      }
    } else {
      stop("Failed to load 'asreml'.\n")
    }
  }
  ## Create output.
  res <- createVarComp(fitMod = mr, modDat = TDTot, trialGroup = nesting,
                       useLocYear = locationYear, engine = engine)
  return(res)
}


