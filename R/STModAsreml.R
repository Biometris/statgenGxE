#' Fit Single Trial Model using asreml
#'
#' Fit Single Trial Model using asreml
#'
#' @inheritParams STRunModel
#'
#' @seealso \code{\link{STRunModel}}
#'
#' @keywords internal
STModAsreml <- function(TD,
                        trials = names(TD),
                        traits,
                        what = c("fixed", "random"),
                        covariates = NULL,
                        useCheckId = FALSE,
                        design = "rowcol",
                        trySpatial = FALSE,
                        control = NULL,
                        checks = TRUE,
                        ...) {
  if (checks) {
    ## Checks.
    checkOut <- modelChecks(TD = TD, trials = trials, design = design,
                            traits = traits, what = what,
                            covariates = covariates, trySpatial = trySpatial,
                            engine = "asreml", useCheckId = useCheckId,
                            control = control)
    ## Convert output to variables.
    list2env(x = checkOut, envir = environment())
  }
  models <- sapply(X = trials, FUN = function(trial) {
    TDTr <- TD[[trial]]
    ## Should repId be used as fixed effect in the model.
    useRepIdFix <- design %in% c("res.ibd", "res.rowcol", "rcbd")
    # Check if spatial models can be fitted.
    if (trySpatial) {
      ## Set default value for criterion
      criterion = "AIC"
      if ("criterion" %in% names(control)) {
        critCt <- control$criterion
        if (critCt %in% c("AIC", "BIC")) {
          criterion <- critCt
        } else {
          warning(paste("Invalid value for control parameter criterion.",
                        "Using default value instead.\n"))
        }
      }
      if (useRepIdFix) {
        repTab <- table(TDTr$repId, TDTr$rowId, TDTr$colId)
      } else {
        repTab <- table(TDTr$rowId, TDTr$colId)
      }
      if (min(repTab) > 1) {
        warning(paste("There should only be one plot at each combination of",
                      if (useRepIdFix) "replicate", "row and column.\n",
                      "Spatial models will not be tried"))
        trySpatial <- FALSE
      }
    }
    if (!trySpatial) {
      ## Indicate extra random effects.
      if (design %in% c("ibd", "res.ibd")) {
        randEff <- "subBlock"
      } else if (design %in% c("rowcol", "res.rowcol")) {
        randEff <- c("rowId", "colId")
      } else if (design == "rcbd") {
        randEff <- character()
      }
      ## Create tempfile to suppress asreml output messages.
      tmp <- tempfile()
      ## Construct formula for fixed part.
      fixedForm <- paste("~",
                         if (useRepIdFix) "repId" else "1",
                         if (useCheckId) "+ checkId",
                         if (!is.null(covariates)) paste(c("", covariates),
                                                         collapse = "+"))
      ## Construct formula for random part. Include repId depending on design.
      if (length(randEff) != 0) {
        randomForm <- paste0(if (useRepIdFix) "repId:",
                             paste(randEff,
                                   collapse = paste("+",
                                                    if (useRepIdFix) "repId:")))
      } else {
        randomForm <- character()
      }
      ## Create empty base lists.
      mr <- mf <- setNames(vector(mode = "list", length = length(traits)),
                           traits)
      for (trait in traits) {
        if ("random" %in% what) {
          ## Fit model with genotype random.
          sink(file = tmp)
          mrTrait <- asreml::asreml(fixed = as.formula(paste(trait, fixedForm)),
                                    random = as.formula(paste("~", randomForm,
                                                              if (length(randomForm) != 0) "+",
                                                              "genotype")),
                                    rcov = ~ units, aom = TRUE, data = TDTr, ...)
          sink()
          if ("fixed" %in% what) {
            ## Constrain variance of the variance components to be fixed as the values in mr.
            GParamTmp <- mrTrait$G.param
            for (randEf in randEff) {
              ## When there are no replicates the structure is [[randEf]][[randEf]]
              ## otherwise it is [[repId:randEf]][[repId]]
              GParamTmp[[paste0(ifelse(useRepIdFix, "repId:", ""),
                                randEf)]][[ifelse(useRepIdFix, "repId", randEf)]]$con <- "F"
            }
          }
          ## evaluate call terms in mr and mfTrait so predict can be run.
          mrTrait$call$fixed <- eval(mrTrait$call$fixed)
          mrTrait$call$random <- eval(mrTrait$call$random)
          mrTrait$call$rcov <- eval(mrTrait$call$rcov)
          # Run predict.
          mrTrait <- predictAsreml(mrTrait, TD = TDTr)
          mrTrait$call$data <- substitute(TDTr)
          mr[[trait]] <- mrTrait
        }
        if ("fixed" %in% what) {
          ## Fit model with genotype fixed.
          if (!"random" %in% what) {
            GParamTmp <- NULL
          }
          sink(file = tmp)
          if (length(randomForm) != 0) {
            mfTrait <- asreml::asreml(fixed = as.formula(paste(trait, fixedForm,
                                                               "+ genotype")),
                                      random = as.formula(paste("~",
                                                                randomForm)),
                                      rcov = ~ units, G.param = GParamTmp,
                                      aom = TRUE, data = TDTr, ...)
          } else {
            mfTrait <- asreml::asreml(fixed = as.formula(paste(trait, fixedForm,
                                                               "+ genotype")),
                                      rcov = ~ units, G.param = GParamTmp,
                                      aom = TRUE, data = TDTr, ...)
          }
          sink()
          mfTrait$call$fixed <- eval(mfTrait$call$fixed)
          mfTrait$call$random <- eval(mfTrait$call$random)
          mfTrait$call$rcov <- eval(mfTrait$call$rcov)
          ## Construct assocForm for use in associate in predict.
          if (useCheckId) {
            assocForm <- as.formula("~ checkId:genotype")
          } else {
            assocForm <- as.formula("~ NULL")
          }
          ## Run predict.
          mfTrait <- predictAsreml(mfTrait, TD = TDTr, associate = assocForm)
          mfTrait$call$data <- substitute(TDTr)
          mf[[trait]] <- mfTrait
        }
      }
      unlink(tmp)
      ## Construct SSA object.
      return(list(mRand = if ("random" %in% what) mr else NULL,
                  mFix = if ("fixed" %in% what) mf else NULL, TD = TDTr,
                  traits = traits, design = design, engine = "asreml",
                  predicted = "genotype"))
    } else {# trySpatial
      regular <- min(repTab) == 1 && max(repTab) == 1
      return(bestSpatMod(TD = TDTr, traits = traits, what = what,
                         regular = regular, criterion = criterion,
                         useCheckId = useCheckId, design = design,
                         covariates = covariates, ...))
    }
  }, simplify = FALSE)
  ## Construct and return SSA object.
  return(createSSA(models = models))
}

#' Helper function for calculating best spatial model using asreml.
#' @keywords internal
bestSpatMod <- function(TD,
                        traits,
                        what = c("fixed", "random"),
                        regular = TRUE,
                        criterion = "AIC",
                        useCheckId = FALSE,
                        design = "rowcol",
                        covariates = NULL,
                        ...) {
  ## Create tempfile to suppress asreml output messages.
  tmp <- tempfile()
  ## TD needs to be sorted by row and column to prevent asreml from crashing.
  TD <- TD[order(TD$rowId, TD$colId), ]
  useRepIdFix <- design == "res.rowcol"
  ## Define random terms of models to try.
  randomTerm <- c(rep(x = "NULL", times = 3),
                  "repId:rowId", "repId:colId", "repId:rowId + repId:colId")
  if (!useRepIdFix) {
    ## If no repId remove this from randomTerm
    randomTerm <- gsub(pattern = "repId:", replacement = "", x = randomTerm)
  }
  if (regular) {
    ## Define spatial terms of models to try.
    spatialChoice <- rep(x = c("exp(x)id", "id(x)exp",
                               "isotropic exponential"), times = 2)
    spatialTerm <- rep(x = c("exp(rowCoordinates):colCoordinates",
                             "rowCoordinates:exp(colCoordinates)",
                             "iexp(rowCoordinates,colCoordinates)"),
                       times = 2)
  } else {
    spatialChoice <- rep(x = c("AR1(x)id", "id(x)AR1", "AR1(x)AR1"), times = 2)
    spatialTerm <- rep(x = c("ar1(rowId):colId",
                             "rowId:ar1(colId)",
                             "ar1(rowId):ar1(colId)"),
                       times = 2)
  }
  ## Create empty base lists.
  mr <- mf <- spatial <- setNames(vector(mode = "list", length = length(traits)),
                                  traits)
  for (trait in traits) {
    ## Create formula for the fixed part.
    fixedFormR <- as.formula(paste(trait, "~",
                                   if (useRepIdFix) "repId" else "1",
                                   if (useCheckId) "+ checkId",
                                   if (!is.null(covariates)) paste(c("", covariates),
                                                                   collapse = "+")))
    ## Fit model with genotype random for all different random/spatial terms.
    for (i in 1:length(randomTerm)) {
      sink(file = tmp)
      mrTrait <- asreml::asreml(fixed = fixedFormR,
                                random = as.formula(paste("~ genotype +",
                                                          randomTerm[i])),
                                rcov = as.formula(paste("~", spatialTerm[i])),
                                aom = TRUE, data = TD, ...)
      sink()
      ## If current model is better than best so far based on chosen criterion
      ## define best model as current model.
      if (i == 1) {
        bestModelTrait <- mrTrait
        bestLoc <- 1
      } else {
        if (criterion == "AIC") {
          criterionCur  <- -2 * mrTrait$loglik + 2 * length(mrTrait$gammas)
          criterionPrev <- -2 * bestModelTrait$loglik +
            2 * length(bestModelTrait$gammas)
        } else {
          criterionCur  <- -2 * mrTrait$loglik +
            log(length(mrTrait$fitted.values)) * length(mrTrait$gammas)
          criterionPrev <- -2 * bestModelTrait$loglik +
            log(length(bestModelTrait$fitted.values)) * length(bestModelTrait$gammas)
        }
        if (criterionCur < criterionPrev) {
          bestModelTrait <- mrTrait
          bestLoc <- i
        }
      }
    }
    fixedFormfTrait <- as.formula(paste(deparse(fixedFormR), "+ genotype"))
    ## Constrain variance of the variance components to be fixed as the values in the best model.
    GParamTmp <- bestModelTrait$G.param
    for (randEf in c("rowId", "colId")) {
      ## When there are no replicates the structure is [[randEf]][[randEf]]
      ## otherwise it is [[repId:randEf]][[repId]]
      GParamTmp[[paste0(ifelse(useRepIdFix, "repId:", ""),
                        randEf)]][[ifelse(useRepIdFix, "repId", randEf)]]$con <- "F"
    }
    sink(file = tmp)
    ## Fit the model with genotype fixed only for the best model.
    mfTrait <- asreml::asreml(fixed = fixedFormfTrait,
                              random = as.formula(paste("~", randomTerm[bestLoc])),
                              rcov = as.formula(paste("~", spatialTerm[bestLoc])),
                              G.param = GParamTmp, aom = TRUE, data = TD, ...)
    sink()
    ## evaluate call terms in bestNidek and mfTrait so predict can be run.
    bestModelTrait$call$fixed <- eval(bestModelTrait$call$fixed)
    bestModelTrait$call$random <- eval(bestModelTrait$call$random)
    bestModelTrait$call$rcov <- eval(bestModelTrait$call$rcov)
    # Run predict.
    bestModelTrait <- predictAsreml(bestModelTrait, TD = TD)
    mfTrait$call$fixed <- eval(mfTrait$call$fixed)
    mfTrait$call$random <- eval(mfTrait$call$random)
    mfTrait$call$rcov <- eval(mfTrait$call$rcov)
    ## Construct assocForm for use in associate in predict.
    if (useCheckId) {
      assocForm <- as.formula("~ checkId:genotype")
    } else {
      assocForm <- as.formula("~ NULL")
    }
    ## Run predict.
    mfTrait <- predictAsreml(mfTrait, TD = TD, associate = assocForm)
    mr[[trait]] <- bestModelTrait
    mf[[trait]] <- mfTrait
    spatial[[trait]] <- spatialChoice[bestLoc]
  }
  unlink(tmp)
  return(list(mRand = if ("random" %in% what) mr else NULL,
              mFix = if ("fixed" %in% what) mf else NULL, TD = TD,
              traits = trait, design = design, spatial = spatial,
              engine = "asreml", predicted = "genotype"))
}




