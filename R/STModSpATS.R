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
                       trial = NULL,
                       traits,
                       what = c("fixed", "random"),
                       covariates = NULL,
                       useCheckId = FALSE,
                       trySpatial = FALSE,
                       design = "rowcol",
                       control = NULL,
                       checks = TRUE,
                       ...) {
  ## Base check.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (checks) {
    ## Checks.
    checkOut <- modelChecks(TD = TD, trial = trial, design = design,
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
  TDTr <- droplevels(TD[[trial]])
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
      modTrR <- tryCatchExt(
        SpATS::SpATS(response = trait, genotype = "genotype",
                     genotype.as.random = TRUE,
                     spatial = ~ SpATS::PSANOVA(colCoord, rowCoord,
                                                nseg = nSeg,
                                                nest.div = nestDiv),
                     fixed = fixedForm, random = randomForm, data = TDTr,
                     control = list(monitoring = 0), ...)
      )
      if (length(modTrR$warning) != 0) {
        modTrR <- wrnToErr(modTrR)
      }
      if (length(modTrR$warning) != 0) {
        warning(paste0("Warning in SpATS for genotype random, trait ", trait,
                       " in trial ", trial, ":\n", modTrR$warning, "\n"),
                call. = FALSE)
      }
      if (is.null(modTrR$error)) {
        return(modTrR$value)
      } else {
        warning(paste0("Error in SpATS for genotype random, trait ", trait,
                       " in trial ", trial, ":\n", modTrR$error, "\n"),
                call. = FALSE)
        return(NULL)
      }
    }, simplify = FALSE)
  } else {
    mr <- NULL
  }
  if ("fixed" %in% what) {
    mf <- sapply(X = traits, FUN = function(trait) {
      ## Fit model with genotype fixed.
      modTrF <- tryCatchExt(
        SpATS::SpATS(response = trait, genotype = "genotype",
                     genotype.as.random = FALSE,
                     spatial = ~ SpATS::PSANOVA(colCoord, rowCoord,
                                                nseg = nSeg,
                                                nest.div = nestDiv),
                     fixed = fixedForm, random = randomForm, data = TDTr,
                     control = list(monitoring = 0), ...)
      )
      if (length(modTrF$warning) != 0) {
        modTrF <- wrnToErr(modTrF)
      }
      if (length(modTrF$warning) != 0) {
        warning(paste0("Warning in SpATS for genotype fixed, trait ", trait,
                       " in trial ", trial, ":\n", modTrF$warning, "\n"),
                call. = FALSE)
      }
      if (is.null(modTrF$error)) {
        return(modTrF$value)
      } else {
        warning(paste0("Error in SpATS for genotype fixed, trait ", trait,
                       " in trial ", trial, ":\n", modTrF$error, "\n"),
                call. = FALSE)
        return(NULL)
      }
    }, simplify = FALSE)
  } else {
    mf <- NULL
  }
  spatial <- setNames(rep("2 dimensional P-splines", times = length(traits)),
                      traits)
  sumTab <- setNames(vector(mode = "list", length = length(traits)), traits)
  return(list(mRand = mr, mFix = mf, TD = TD[trial], traits = traits,
              design = design, spatial = spatial, engine = "SpATS",
              predicted = "genotype", sumTab = sumTab))
}
