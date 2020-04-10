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
                       trialGroup = NULL,
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
  if (!is.null(trialGroup)) {
    chkCol(trialGroup, TDTot)
  }
  if (useWt) {
    chkCol("wt", TDTot)
  } else {
    TDTot[["wt"]] <- 1
  }
  engine <- match.arg(engine)
  TDTot <- droplevels(TDTot)
  hasYear <- hasName(TDTot, "year")
  useLocYear <- hasYear & hasName(TDTot, "loc")
  envVar <- ifelse(useLocYear, "loc", "trial")


  ## Set to FALSE for developing purposes.
  hasYear <- FALSE


  ## Increase maximum number of iterations for asreml. Needed for more complex
  ## designs to converge.
  maxIter <- 200
  ## Add combinations of trial and genotype currently not in TD to TD.
  ## No missing combinations are allowed when fitting asreml models.
  # fullModLevs <- list(genotype = levels(TDTot[["genotype"]]))
  # fullModLevs[[envVar]] <- levels(TDTot[[envVar]])
  # if (useLocYear) fullModLevs[["year"]] <- levels(TDTot[["year"]])
  # TD0 <- expand.grid(fullModLevs)
  # TDTot <- merge(TD0, TDTot, all.x = TRUE)
  ## Asreml can't handle missing weights, so set them to 0 for added combinations.
  TDTot[is.na(TDTot[["wt"]]), "wt"] <- 0
  ## Check if the trial is nested within the trialGroup.
  hasGroup <- !is.null(trialGroup)
  if (hasGroup) {
    groupTab <- table(TDTot[[envVar]], TDTot[[trialGroup]])
    isNestedTrialGroup <- sum(groupTab > 0) == nlevels(TDTot[[envVar]])
  } else {
    isNestedTrialGroup <- FALSE
  }
  isNestedTrialGroup <- FALSE
  ## Check if the data contains replicates.
  repTab <- table(TDTot[["trial"]], TDTot[["genotype"]])
  hasReps <- any(repTab > 1)
  ## Main procedure to fit mixed models.
  modCols <- c(trialGroup, if (hasYear) "year", envVar)

  ## Construct formula for fixed part - first as text.
  ## Trying to fit this in something 'smart' actually makes it unreadable.
  ## First create a vector with the separate terms.
  ## This avoids difficult constructions to get the +-es correct.
  # fixedTerms <- c(if (!useLocYear && (!hasGroup || (hasReps && useWt))) "trial",
  #                 if (useLocYear) "loc + year + year:loc",
  #                 if (hasGroup && !useLocYear) paste0(trialGroup, "+", envVar, ":", trialGroup),
  #                 if (hasGroup && useLocYear) (trialGroup))

  fixedTerms <- c(if (!useLocYear && !hasGroup) "trial",
                  if (hasGroup) trialGroup,
                  if (!useLocYear && hasGroup) paste0(trialGroup, ":trial"),
                  if (useLocYear && !hasGroup) "loc",
                  if (useLocYear) "year",
                  if (useLocYear && !hasGroup) "loc:year",
                  if (useLocYear && hasGroup) paste0(trialGroup, ":year"),
                  if (useLocYear && hasGroup) paste0(trialGroup, ":loc:year"))



  fixedTxt <- paste0("`", trait, "`~",
                     paste(fixedTerms, collapse = "+"))
  ## Construct formula for random part in a similar way.
  # randTerms <- c("genotype",
  #                if (hasGroup) paste0("genotype:", trialGroup),
  #                if (useLocYear || hasReps) paste0("genotype:", envVar),
  #                if (useLocYear) "genotype:year",
  #                # if (hasGroup && !isNestedTrialGroup && (hasReps || useWt))
  #                  # paste0("genotype:", trialGroup, ":", envVar),
  #                if (useLocYear && (hasReps || useWt)) "genotype:year:loc")

  randTerms <- c("genotype",
                 if (!useLocYear && !hasGroup) "genotype:trial",
                 if (hasGroup) paste0("genotype:", trialGroup),
                 if (!useLocYear && hasGroup && (hasReps || useWt))
                   paste0("genotype:", trialGroup, ":trial"),
                 if (useLocYear && !hasGroup) "genotype:loc",
                 if (useLocYear) "genotype:year",
                 if (useLocYear && !hasGroup) "genotype:loc:year",
                 if (useLocYear && hasGroup) paste0("genotype:", trialGroup, ":year"),
                 if (useLocYear && hasGroup  && (hasReps || useWt))
                   paste0("genotype:", trialGroup, ":loc:year"))


  if (engine == "lme4") {
    randTxt <- paste(paste0("(1|", randTerms, ")"), collapse = "+")
    formTxt <- paste(fixedTxt, "+", randTxt)
    ## Fit the actual model.
    mr <- lme4::lmer(formula(formTxt), data = TDTot,  weights = TDTot[["wt"]],
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
  res <- createVarComp(fitMod = mr, modDat = TDTot, trialGroup = trialGroup,
                       useLocYear = useLocYear, engine = engine)
  return(res)
}


