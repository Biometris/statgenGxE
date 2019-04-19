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
                        trial = NULL,
                        traits,
                        what = c("fixed", "random"),
                        covariates = NULL,
                        useCheckId = FALSE,
                        design = "rowcol",
                        trySpatial = FALSE,
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
                            engine = "asreml", useCheckId = useCheckId,
                            control = control)
    ## Convert output to variables.
    list2env(x = checkOut, envir = environment())
  }
  TDTr <- droplevels(TD[[trial]])
  ## Should repId be used as fixed effect in the model.
  repIdFix <- design %in% c("res.ibd", "res.rowcol", "rcbd")
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
    if (repIdFix) {
      repTab <- table(TDTr$repId, TDTr$rowId, TDTr$colId)
    } else {
      repTab <- table(TDTr$rowId, TDTr$colId)
    }
    if (min(repTab) > 1) {
      warning(paste("There should only be one plot at each combination of",
                    if (repIdFix) "replicate", "row and column.\n",
                    "Spatial models will not be tried"))
      trySpatial <- FALSE
    }
  }
  ## Indicate extra random effects.
  if (design %in% c("ibd", "res.ibd")) {
    randEff <- "subBlock"
  } else if (design %in% c("rowcol", "res.rowcol")) {
    randEff <- c("rowId", "colId")
  } else if (design == "rcbd") {
    randEff <- character()
  }
  ## Construct formula for fixed part.
  fixedForm <- paste("~",
                     if (repIdFix) "repId" else "1",
                     if (useCheckId) "+ checkId",
                     if (!is.null(covariates)) paste(c("", covariates),
                                                     collapse = "+"))
  ## Construct formula for random part. Include repId depending on design.
  if (length(randEff) != 0) {
    randomForm <- paste0(if (repIdFix) "repId:",
                         paste(randEff,
                               collapse = paste("+",
                                                if (repIdFix) "repId:")))
  } else {
    randomForm <- character()
  }
  if (!trySpatial) {
    ## Create tempfile to suppress asreml output messages.
    tmp <- tempfile()
    ## Increase max number of iterations for asreml.
    maxIter <- 200
    ## Create empty base lists.
    mr <- mf <- spatial <- sumTab <- setNames(vector(mode = "list",
                                                     length = length(traits)),
                                              traits)
    for (trait in traits) {
      if ("random" %in% what) {
        ## Fit model with genotype random.
        sink(file = tmp)
        mrTrait <- tryCatchExt(
          asreml::asreml(fixed = formula(paste(trait, fixedForm)),
                         random = formula(paste("~", randomForm,
                                                if (length(randomForm) != 0) "+",
                                                "genotype")),
                         aom = TRUE, data = TDTr, maxiter = maxIter, ...))
        sink()
        if (!is.null(mrTrait$warning)) {
          mrTrait <- chkLastIter(mrTrait)
          mrTrait <- wrnToErr(mrTrait)
        }
        if (length(mrTrait$warning) != 0) {
          warning(paste0("Warning in asreml for genotype random, trait ", trait,
                         " in trial ", trial, ":\n", mrTrait$warning, "\n"),
                  call. = FALSE)
        }
        if (is.null(mrTrait$error)) {
          mrTrait <- mrTrait$value
        } else {
          warning(paste0("Error in asreml for genotype random, trait ", trait,
                         " in trial ", trial, ":\n", mrTrait$error, "\n"),
                  call. = FALSE)
          mrTrait <- NULL
        }
        if (!is.null(mrTrait)) {
          if ("fixed" %in% what) {
            ## Constrain variance of the variance components to be fixed as
            ## the values in mr.
            GParamTmp <- mrTrait$G.param
            for (randEf in randEff) {
              ## When there are no replicates the structure is
              ## [[randEf]][[randEf]] otherwise it is [[repId:randEf]][[repId]]
              GParamTmp[[paste0(ifelse(repIdFix, "repId:", ""),
                                randEf)]][[ifelse(repIdFix,
                                                  "repId", randEf)]]$con <- "F"
            }
          }
          ## evaluate call terms in mr and mfTrait so predict can be run.
          mrTrait$call$fixed <- eval(mrTrait$call$fixed)
          mrTrait$call$random <- eval(mrTrait$call$random)
          mrTrait$call$rcov <- eval(mrTrait$call$rcov)
          mrTrait$call$data <- substitute(TDTr)
        }
        mr[[trait]] <- mrTrait
      } # End random.
      if ("fixed" %in% what) {
        ## Fit model with genotype fixed.
        if (!"random" %in% what) {
          GParamTmp <- NULL
        }
        sink(file = tmp)
        ## There is no way to specify random as NULL or NA so therefore split
        ## cases on whether there is or there is not a random term.
        if (length(randomForm) != 0) {
          mfTrait <- tryCatchExt(
            asreml::asreml(fixed = formula(paste(trait, fixedForm,
                                                 "+ genotype")),
                           random = formula(paste("~", randomForm)),
                           G.param = GParamTmp, aom = TRUE, data = TDTr,
                           maxiter = maxIter, ...))
        } else {
          mfTrait <- tryCatchExt(
            asreml::asreml(fixed = formula(paste(trait, fixedForm,
                                                 "+ genotype")),
                           G.param = GParamTmp, aom = TRUE, data = TDTr,
                           maxiter = maxIter, ...))
        }
        sink()
        if (!is.null(mfTrait$warning)) {
          mfTrait <- chkLastIter(mfTrait)
          mfTrait <- wrnToErr(mfTrait)
        }
        if (length(mfTrait$warning) != 0) {
          warning(paste0("Warning in asreml for genotype fixed, trait ", trait,
                         " in trial ", trial, ":\n", mfTrait$warning, "\n"),
                  call. = FALSE)
        }
        if (is.null(mfTrait$error)) {
          mfTrait <- mfTrait$value
        } else {
          warning(paste0("Error in asreml for genotype fixed, trait ", trait,
                         " in trial ", trial, ":\n", mfTrait$error, "\n"),
                  call. = FALSE)
          mfTrait <- NULL
        }
        if (!is.null(mfTrait)) {
          mfTrait$call$fixed <- eval(mfTrait$call$fixed)
          mfTrait$call$random <- eval(mfTrait$call$random)
          mfTrait$call$rcov <- eval(mfTrait$call$rcov)
          mfTrait$call$data <- substitute(TDTr)
          ## Construct assocForm for use in associate in predict.
          if (useCheckId) {
            assocForm <- formula("~ checkId:genotype")
          } else {
            assocForm <- formula("~ NULL")
          }
        }
        mf[[trait]] <- mfTrait
      } # End fixed.
      spatial[trait] <- FALSE
    } # End for traits.
    unlink(tmp)
    ## Construct SSA object.
    return(list(mRand = if ("random" %in% what) mr else NULL,
                mFix = if ("fixed" %in% what) mf else NULL, TD = TD[trial],
                traits = traits, design = design, spatial = spatial,
                engine = "asreml", predicted = "genotype", sumTab = sumTab))
  } else {# trySpatial
    regular <- min(repTab) == 1 && max(repTab) == 1
    return(bestSpatMod(TD = TD[trial], traits = traits, what = what,
                       regular = TRUE, criterion = criterion,
                       #regular = regular, criterion = criterion,
                       useCheckId = useCheckId, design = design,
                       covariates = covariates, fixedForm = fixedForm,
                       randomForm = randomForm, ...))
  }
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
                        fixedForm,
                        randomForm,
                        ...) {
  dotArgs <- list(...)
  ## Create tempfile to suppress asreml output messages.
  tmp <- tempfile()
  ## Increase max number of iterations for asreml.
  maxIter <- 200
  TDTr <- droplevels(TD[[1]])
  ## Add empty observations.
  TDTab <- as.data.frame(table(TDTr$colId, TDTr$rowId))
  TDTab <- TDTab[TDTab$Freq == 0, , drop = FALSE]
  if (nrow(TDTab) > 0) {
    extObs <- setNames(as.data.frame(matrix(nrow = nrow(TDTab),
                                            ncol = ncol(TDTr))),
                       colnames(TDTr))
    extObs$trial <- TDTr$trial[1]
    extObs[, c("colId", "rowId")] <- TDTab[, c("Var1", "Var2")]
    extObs[, c("colCoord", "rowCoord")] <-
      c(as.numeric(levels(TDTab[, "Var1"]))[TDTab[, "Var1"]],
        as.numeric(levels(TDTab[, "Var2"]))[TDTab[, "Var2"]])
    TDTr <- rbind(TDTr, extObs)
  }
  ## TD needs to be sorted by row and column to prevent asreml from crashing.
  TDTr <- TDTr[order(TDTr$rowId, TDTr$colId), ]
  repIdFix <- design == "res.rowcol"
  ## Define random terms of models to try.
  randTerm <- c("NULL", rep(x = c("NULL", "units"), each = 3))
  if (regular) {
    ## Define spatial terms of models to try.
    spatCh <- c("none", rep(x = c("exp(x)id", "id(x)exp",
                                  "isotropic exponential"), times = 2))
    spatTerm <- c(NA, paste("~", rep(x = c("exp(rowCoord):colCoord",
                                           "rowCoord:exp(colCoord)",
                                           "iexp(rowCoord,colCoord)"),
                                     times = 2)))
  } else {
    spatCh <- c("none", rep(x = c("AR1(x)id", "id(x)AR1", "AR1(x)AR1"),
                            times = 2))
    spatTerm <- c(NA, paste("~",
                            rep(x = c("ar1(rowId):colId", "rowId:ar1(colId)",
                                      "ar1(rowId):ar1(colId)"), times = 2)))
  }
  ## Create empty base lists.
  mr <- mf <- spatial <- sumTab <- setNames(vector(mode = "list",
                                                   length = length(traits)),
                                            traits)
  btCols <- c("spatial", "random", "AIC", "BIC", "row", "col", "error",
              "correlated error", "converge")
  for (trait in traits) {
    ## Reset criterion to Inf.
    criterionBest <- Inf
    ## Create data.frame for storing summary for current trait.
    modSum <- as.data.frame(matrix(nrow = length(spatCh), ncol = length(btCols),
                                   dimnames = list(NULL, btCols)))
    ## Create formula for the fixed part.
    fixedFormR <- formula(paste(trait, fixedForm))
    ## Fit model with genotype random for all different random/spatial terms.
    for (i in seq_along(randTerm)) {
      if (length(randomForm) > 0) {
        randFormR <- formula(paste("~ genotype +", randomForm, "+", randTerm[i]))
      } else {
        randFormR <- formula(paste("~ genotype +", randTerm[i]))
      }
      asrArgsR <- c(list(fixed = fixedFormR, random = randFormR, aom = TRUE,
                         data = TDTr, maxiter = maxIter), dotArgs)
      ## In asreml3 na.method.X and na.method.Y are used.
      ## In asreml4 this is replaced by na.action.
      asrArgsR <- c(asrArgsR, if (asreml4()) {
        list(na.action = asreml::na.method(x = "include"))
      } else {
        list(na.method.X = "include")
      })
      if (!is.na(spatTerm[i])) {
        ## In asreml4 rcov is replaced by residual.
        asrArgsR[[ifelse(asreml4(), "residual", "rcov")]] <- formula(spatTerm[i])
      }
      sink(file = tmp)
      mrTrait <- tryCatchExt(do.call(what = asreml::asreml, args = asrArgsR))
      sink()
      if (!is.null(mrTrait$warning)) {
        mrTrait <- chkLastIter(mrTrait)
        mrTrait <- wrnToErr(mrTrait)
      }
      if (length(mrTrait$warning) != 0) {
        warning(paste0("Warning in asreml for model ", spatCh[i],
                       " genotype random, trait ", trait, " in trial ",
                       TDTr$trial[1], ":\n", mrTrait$warning,
                       "\n"), call. = FALSE)
      }
      if (is.null(mrTrait$error)) {
        mrTrait <- mrTrait$value
      } else {
        warning(paste0("Error in asreml for model ", spatCh[i],
                       " genotype random, trait ", trait, " in trial ",
                       TDTr$trial[1], ":\n", mrTrait$error,
                       "\n"), call. = FALSE)
        mrTrait <- NULL
      }
      ## Fill model summary table.
      modSum[i, "spatial"] <- spatCh[i]
      modSum[i, "random"] <- randTerm[i]
      modSum[i, "converge"] <- isTRUE(!is.null(mrTrait) & mrTrait$converge)
      if (!is.null(mrTrait)) {
        summ <- summary(mrTrait)$varcomp["component"]
        modSum[i, "AIC"] <- -2 * mrTrait$loglik + 2 * nrow(summ)
        modSum[i, "BIC"] <- -2 * mrTrait$loglik +
          log(length(fitted(mrTrait))) * nrow(summ)
        ## Row and column output differs for regular/non-regular.
        ## Always max. one of the possibilities is in summary so rowVal and
        ## colVal are always a single value.
        rowVal <- summ[rownames(summ) %in%
                         c("R!rowId.cor", "R!rowCoord.pow", "R!pow",
                           ## New naming for asreml4.
                           "rowId:colId!rowId!cor",
                           "rowCoord:colCoord!rowCoord!pow",
                           "iexp(rowCoord,colCoord)!pow"), ]
        modSum[i, "row"] <- ifelse(length(rowVal) == 0, NA, rowVal)
        colVal <- summ[rownames(summ) %in%
                         c("R!colId.cor", "R!colCoord.pow", "R!pow",
                           ## New naming for asreml4
                           "rowId:colId!colId!cor",
                           "rowCoord:colCoord!colCoord!pow",
                           "iexp(rowCoord,colCoord)!pow"), ]
        modSum[i, "col"] <- ifelse(length(colVal) == 0, NA, colVal)
        modSum[i, "error"] <- ifelse(randTerm[i] == "units",
                                     summ[rownames(summ) %in%
                                            c("units!units.var",
                                              ## New naming for asreml4.
                                              "units"), ],
                                     summ[rownames(summ) %in%
                                            c("R!variance",
                                              ## New naming for asreml4.
                                              "units!R", "rowId:colId!R",
                                              "rowCoord:colCoord!R",
                                              "iexp(rowCoord,colCoord)!R"), ])
        if (randTerm[i] == "units") {
          modSum[i, "correlated error"] <- summ[rownames(summ) %in%
                                                  c("R!variance",
                                                    ## New naming for asreml4.
                                                    "rowId:colId!R",
                                                    "rowCoord:colCoord!R",
                                                    "iexp(rowCoord,colCoord)!R"), ]
        }
        ## If current model is better than best so far based on chosen criterion
        ## define best model as current model.
        if (mrTrait$converge) {
          if (criterion == "AIC") {
            criterionCur <- modSum[i, "AIC"]
          } else {
            criterionCur <- modSum[i, "BIC"]
          }
        }
        if (criterionCur < criterionBest) {
          bestModTr <- mrTrait
          ## Evaluate call terms in bestModTr and mfTrait so predict can be run.
          ## Needs to be called in every iteration to prevent final result
          ## from always having the values of the last iteration.
          bestModTr$call$fixed <- eval(bestModTr$call$fixed)
          bestModTr$call$random <- eval(bestModTr$call$random)
          bestModTr$call$rcov <- eval(bestModTr$call$rcov)
          criterionBest <- criterionCur
          bestMod <- i
        }
      }
    }
    fixedFormfTrait <- update(fixedFormR, ~ . + genotype)
    ## Constrain variance of the variance components to be fixed as the values
    ## in the best model.
    GParamTmp <- bestModTr$G.param
    for (randEf in c("rowId", "colId")) {
      ## When there are no replicates the structure is [[randEf]][[randEf]]
      ## otherwise it is [[repId:randEf]][[repId]]
      GParamTmp[[paste0(ifelse(repIdFix, "repId:", ""),
                        randEf)]][[ifelse(repIdFix, "repId", randEf)]]$con <- "F"
    }
    if (length(randomForm) > 0) {
      randFormF <- formula(paste("~", randomForm, "+", randTerm[bestMod]))
    } else {
      randFormF <- formula(paste("~", randTerm[bestMod]))
    }
    asrArgsF <- c(list(fixed = fixedFormfTrait, random = randFormF,
                       G.param = GParamTmp, aom = TRUE, data = TDTr,
                       maxiter = maxIter), dotArgs)
    ## In asreml3 na.method.X and na.method.Y are used.
    ## In asreml4 this is replaced by na.action.
    asrArgsR <- c(asrArgsR, if (asreml4()) {
      list(na.action = asreml::na.method(x = "include"))
    } else {
      list(na.method.X = "include")
    })
    if (!is.na(spatTerm[bestMod])) {
      ## In asreml4 parameter rcov is replaced by residual.
      asrArgsF[[ifelse(asreml4(), "residual", "rcov")]] <- formula(spatTerm[bestMod])
    }
    sink(file = tmp)
    ## Fit the model with genotype fixed only for the best model.
    mfTrait <- tryCatchExt(do.call(what = asreml::asreml, args = asrArgsF))
    sink()
    if (!is.null(mfTrait$warning)) {
      mfTrait <- chkLastIter(mfTrait)
      mfTrait <- wrnToErr(mfTrait)
    }
    if (length(mfTrait$warning) != 0) {
      warning(paste0("Warning in asreml for model ", spatCh[i],
                     " genotype fixed, trait ", trait, " in trial ",
                     TDTr$trial[1], ":\n", mrTrait$warning,
                     "\n"), call. = FALSE)
    }
    if (is.null(mfTrait$error)) {
      mfTrait <- mfTrait$value
    } else {
      warning(paste0("Error in asreml for model ", spatCh[i],
                     " genotype fixed, trait ", trait, " in trial ",
                     TDTr$trial[1], ":\n", mfTrait$error,
                     "\n"), call. = FALSE)
      mfTrait <- NULL
    }
    if (!is.null(mfTrait)) {
      mfTrait$call$fixed <- eval(mfTrait$call$fixed)
      mfTrait$call$random <- eval(mfTrait$call$random)
      mfTrait$call$rcov <- eval(mfTrait$call$rcov)
    }
    mr[[trait]] <- bestModTr
    mf[[trait]] <- mfTrait
    spatial[[trait]] <- spatCh[bestMod]
    attr(x = modSum, which = "chosen") <- bestMod
    sumTab[[trait]] <- modSum
  } # End for traits.
  unlink(tmp)
  TD[[1]] <- TDTr
  return(list(mRand = if ("random" %in% what) mr else NULL,
              mFix = if ("fixed" %in% what) mf else NULL, TD = TD,
              traits = traits, design = design, spatial = spatial,
              engine = "asreml", predicted = "genotype", sumTab = sumTab))
}
