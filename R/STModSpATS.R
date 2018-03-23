#' Fit Single Trial Model using SpATS
#'
#' Fit Single Trial Model using SpATS
#'
#' @inheritParams STRunModel
#'
#' @seealso \code{\link{STRunModel}}
#'
#' @keywords internal
STModSpATS <- function(TD,
                       trials = names(TD),
                       traits,
                       what = c("fixed", "random"),
                       covariates = NULL,
                       useCheckId = FALSE,
                       trySpatial = FALSE,
                       design = "rowcol",
                       control = NULL,
                       checks = TRUE,
                       ...) {
  if (checks) {
    ## Checks.
    checkOut <- modelChecks(TD = TD, trials = trials, design = design,
                            traits = traits, what = what,
                            covariates = covariates, trySpatial = trySpatial,
                            engine = "SpATS", useCheckId = useCheckId,
                            control = control)
    ## Convert output to variables.
    list2env(x = checkOut, envir = environment())
  }
  ## Set default value for nestDiv
  nestDiv <- 2
  ## If valid values for nSeg are provided in control use these instead.
  if ("nestDiv" %in% names(control)) {
    nestDivCt <- control$nestDiv
    if (length(nestDivCt) == 1) {
      nestDivCt <- rep(x = nestDivCt, times = 2)
    }
    if (is.numeric(nestDivCt) && length(nestDivCt) <= 2 && all(nestDivCt >= 1)) {
      nestDiv <- nestDivCt
    } else {
      warning(paste("Invalid value for control parameter nestDiv.",
                    "Using default values instead.\n"))
    }
  }
  models <- sapply(X = trials, FUN = function(trial) {
    TDTr <- TD[[trial]]
    ## Should repId be used as fixed effect in the model.
    useRepIdFix <- design %in% c("res.ibd", "res.rowcol", "rcbd")
    ## Indicate extra random effects.
    if (design %in% c("ibd", "res.ibd")) {
      randEff <- "subBlock"
    } else if (design %in% c("rowcol", "res.rowcol")) {
      randEff <- c("rowId", "colId")
    } else if (design == "rcbd") {
      randEff <- character()
    }
    ## Compute number of segments.
    nSeg <- c(ceiling(nlevels(TDTr$colId) / 2), ceiling(nlevels(TDTr$rowId) / 2))
    ## If valid values for nSeg are provided in control use these instead.
    if ("nSeg" %in% names(control)) {
      nSegCt <- control$nSeg
      if (length(nSegCt) == 1) {
        nSegCt <- rep(x = nSegCt, times = 2)
      }
      if (is.numeric(nSegCt) && length(nSegCt) <= 2 && all(nSegCt >= 1) &&
          all(nSegCt <= c(nlevels(TDTr$colId), nlevels(TDTr$rowId)))) {
        nSeg <- nSegCt
      } else {
        warning(paste("Invalid value for control parameter nSeg.",
                      "Using default values instead.\n"))
      }
    }
    ## Construct formula for fixed part.
    fixedForm <- as.formula(paste("~",
                                  if (useRepIdFix) "repId" else "1",
                                  if (useCheckId) "+ checkId",
                                  if (!is.null(covariates)) paste(c("", covariates),
                                                                  collapse = "+")))
    ## Construct formula for random part. Include repId depending on design.
    if (length(randEff) != 0) {
      randomForm <- as.formula(paste0("~", if (useRepIdFix) "repId:",
                                      "(", paste(randEff, collapse = "+"), ")"))
    } else {
      randomForm <- NULL
    }
    if ("random" %in% what) {
      mr <- sapply(X = traits, FUN = function(trait) {
        ## Fit model with genotype random.
        SpATS::SpATS(response = trait, genotype = "genotype",
                     genotype.as.random = TRUE,
                     spatial = ~ SpATS::PSANOVA(colCoordinates, rowCoordinates,
                                                nseg = nSeg, nest.div = nestDiv),
                     fixed = fixedForm,
                     random = randomForm,
                     data = TDTr, control = list(monitoring = 0), ...)
      }, simplify = FALSE)
    } else {
      mr <- NULL
    }
    if ("fixed" %in% what) {
      mf <- sapply(X = traits, FUN = function(trait) {
        ## Fit model with genotype fixed.
        SpATS::SpATS(response = trait, genotype = "genotype",
                     genotype.as.random = FALSE,
                     spatial = ~ SpATS::PSANOVA(colCoordinates, rowCoordinates,
                                                nseg = nSeg, nest.div = nestDiv),
                     fixed = fixedForm,
                     random = randomForm,
                     data = TDTr, control = list(monitoring = 0), ...)
      }, simplify = FALSE)
    } else {
      mf <- NULL
    }

    return(list(mRand = mr, mFix = mf, TD = TD[trial], traits = traits,
                  design = design, spatial = "2 dimensional P-splines",
                  engine = "SpATS", predicted = "genotype"))
  }, simplify = FALSE)
  ## Construct and return SSA object.
  return(createSSA(models = models))
}
