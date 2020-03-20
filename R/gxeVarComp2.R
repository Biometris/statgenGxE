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
#' @param group A character string specifying a column in TD.......
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
gxeVarComp2 <- function(TD,
                        trials = names(TD),
                        trait,
                        engine = c("lme4", "asreml"),
                        group = NULL,
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
  if (!is.null(group)) {
    chkCol(group, TDTot)
  }
  if (useWt) {
    chkCol("wt", TDTot)
  } else {
    TDTot[["wt"]] <- 1
  }
  engine <- match.arg(engine)
  TDTot <- droplevels(TDTot)
  hasYear <- hasName(TDTot, "year")


  ## Set to FALSE for developing purposes.
  hasYear <- FALSE


  ## Increase maximum number of iterations for asreml. Needed for more complex
  ## designs to converge.
  maxIter <- 200
  ## Add combinations of trial and genotype currently not in TD to TD.
  ## No missing combinations are allowed when fitting asreml models.
  TD0 <- expand.grid(genotype = levels(TDTot[["genotype"]]),
                     trial = levels(TDTot[["trial"]]))
  TDTot <- merge(TD0, TDTot, all.x = TRUE)
  ## Asreml can't handle missing weights, so set them to 0 for added combinations.
  TDTot[is.na(TDTot[["wt"]]), "wt"] <- 0
  ## Check if the trial is nested within the group.
  hasGroup <- !is.null(group)
  if (hasGroup) {
    groupTab <- table(TDTot[["trial"]], TDTot[[group]])
    isNestedTrialGroup <- sum(groupTab > 0) == nlevels(TDTot[["trial"]])
  } else {
    isNestedTrialGroup <- FALSE
  }
  ## Check if the data contains replicates.
  repTab <- table(TDTot[["trial"]], TDTot[["genotype"]])
  hasReps <- any(repTab > 1)
  ## Main procedure to fit mixed models.
  modCols <- c(group, if (hasYear) "year", "trial")
  if (engine == "lme4") {
    ## Just the basic model.
    mr <- lme4::lmer(formula(paste0("`", trait, "`~ (",
                                    paste(modCols, collapse = "+"), ")^2",
                                    "+ (1|genotype)",
                                    if (!is.null(group)) {
                                      paste(" + (", group, "|genotype)")
                                    })),
                    data = TDTot, weights = TDTot[["wt"]], ...)
    ## Construct STA object.
  } else if (engine == "asreml") {
    if (requireNamespace("asreml", quietly = TRUE)) {
      ## Construct formula for fixed past - first as text.
      ## Trying to fit this in something 'smart' actually makes it unreadable.
      ## First create a vector with the separate terms.
      ## This avoids difficult constructions to get the +-es correct.
      fixedTerms <- c(if (!isNestedTrialGroup) "trial",
                      if (hasYear) "year",
                      if (hasGroup) paste0(group, "+ trial:", group))
      fixedTxt <- paste0("`", trait, "`~",
                         paste(fixedTerms, collapse = "+"))
      ## Construct formula for random part in a similar way.
      randTerms <- c("genotype",
                     if (hasGroup) paste0("genotype:", group),
                     if (hasReps) "genotype:trial")
      randTxt <- paste("~ ", paste(randTerms, collapse = "+"))
      ## Put arguments for models in a list to make it easier to switch
      ## between asreml3 and asreml4. Usually only one or two arguments differ.
      ## Also some arguments are identical for all models
      modArgs0 <- list(fixed = formula(fixedTxt), random = formula(randTxt),
                       data = TDTot, weights = "wt", maxiter = maxIter,
                       trace = TRUE)
      modArgs <- modArgs0
      # modArgs[[ifelse(asreml4(), "residual", "rcov")]] <-
      #   formula("~genotype:trial")
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
  # model <- setNames(list(list(mRand = NULL,
  #                             mFix = setNames(list(mr), trait),
  #                             TD = createTD(TDTot), traits = trait,
  #                             engine = engine, predicted = "trial")),
  #                   rownames(bestTab)[1])
  # STA <- createSTA(models = model)
  # res <- createVarComp(STA = STA, choice = rownames(bestTab)[1],
  #                      summary = bestTab, vcov = vcovBest,
  #                      criterion = criterion, engine = engine)
  return(mr)
}



