#' Extract statistics from Fitted Models
#'
#' This function extracts and calculates various results for fitted models such
#' as BLUEs, BLUPs, unit errors and heritabilities. Note that most results can
#' only be calculated if a model is fitted with genotype as fixed or random.
#' This is indicated in the list below with "F" and "R"
#'
#' Possible options for \code{what} are:
#' \describe{
#' \item{F - BLUEs}{Best Lineair Unbiased Estimators.}
#' \item{F - seBLUES}{Standard errors of the BLUEs.}
#' \item{R - BLUPs}{Best Lineair Unbiased Predictors.}
#' \item{R - seBLUPs}{Standard errors of the BLUPs.}
#' \item{F - ue}{Unit errors - only for \code{lme4} and \code{asreml}.}
#' \item{R - heritability}{Heritability.}
#' \item{F - varCompF}{Variance components for model with genotype as fixed
#' component.}
#' \item{R - varCompR}{Variance components for model with genotype as random
#' component.}
#' \item{R - varGen}{Genetic variance component(s).}
#' \item{R - varErr}{Residual variance component(s) - only for \code{lme4}
#' and \code{asreml}.}
#' \item{R - varSpat}{Spatial variance components - only for \code{SpATS}.}
#' \item{F - fitted}{Fitted values for the model with genotype as fixed
#' component.}
#' \item{F - resid}{Residuals for the model with genotype as fixed component.}
#' \item{F - stdRes}{Standardized residuals for the model with genotype as fixed
#' component}
#' \item{R - rMeans}{Fitted values for the model with genotype as random
#' component.}
#' \item{R - ranEf}{Random genetic effects.}
#' \item{F - residR}{Residuals for the model with genotype as random component.}
#' \item{F - stdResR}{Standardized residuals for the model with genotype as
#' random component}
#' \item{F - wald}{Results of the wald test - only for \code{lme4} and
#' \code{asreml}.}
#' \item{F - CV}{Coefficient of variation - only for \code{lme4} and
#' \code{asreml}.}
#' \item{F - rDf}{Residual degrees of freedom for the model with genotype as
#' fixed component.}
#' \item{R - rDfR}{Residual degrees of freedom for the model with genotype as
#' random component.}
#' \item{R - effDim}{Effective dimensions - only for \code{SpATS}.}
#' \item{R - ratEffDim}{Ratio's of the effective dimensions -
#' only for \code{SpATS}.}
#' \item{F - sed}{Standard error of difference - only for \code{asreml}.}
#' \item{F - lsd}{Least significant difference - only for \code{asreml}.}
#' \item{all}{All available statistics.}
#' }
#'
#' @param SSA An object of class SSA.
#' @param trials A character vector of trials for which the statistics should be
#' computed. If not supplied, statistics are computed for all trials that have
#' been modeled.
#' @param traits A character vector of traits for which the statistics should be
#' computed. If not supplied, statistics are computed for all traits that have
#' been modeled.
#' @param what A character vector indicating which statistics should be
#' computed. Most statistics are available for all models, some only for models
#' fitted using a certain engine. If this is the case, this is indicated in the
#' list with options in details.\cr
#' If \code{what = "all"}, all available statistics are computed.
#' @param keep A character vector of column(s) in the object of class
#' \code{\link{TD}} used for modeling. These columns will be kept as output when
#' computing fitted values, residuals, standardized residuals and rMeans.
#' Columns can also be kept when computing (se)BLUEs and (se)BLUPs but only if
#' the column to keep contains unique values for the modeled variables, i.e. a
#' column repId with several different values per genotype cannot be kept.
#' @param restoreColNames Should the original column names be restored in the
#' output of the extracted data?
#'
#' @return A list with, per trial for which statistics have been extracted, a
#' list of those statistics.
#'
#' @seealso \code{\link{fitTD}}
#'
#' @examples
#' ## Fit model using SpATS.
#' myModel <- fitTD(TD = TDHeat05, design = "res.rowcol", traits = "yield")
#' ## Extract all available statistics from the fitted model.
#' extr <- extract(myModel)
#' ## Extract only the BLUEs from the fitted model.
#' BLUEs <- extract(myModel, what = "BLUEs")
#' ## Extract only the BLUEs from the fitted model and keep trial as variable in
#' ## the output.
#' BLUEs2 <- extract(myModel, what = "BLUEs", keep = "trial")
#'
#' @importFrom utils capture.output
#' @export
extract <- function(SSA,
                    trials = names(SSA),
                    traits = NULL,
                    what = "all",
                    keep = NULL,
                    restoreColNames = FALSE) {
  ## Checks.
  if (!inherits(SSA, "SSA")) {
    stop("SSA has to be an object of class SSA.\n")
  }
  if (is.null(trials) || !is.character(trials) ||
      !all(trials %in% names(SSA))) {
    stop("All trials should be in SSA.")
  }
  if (!is.null(traits) && !is.character(traits)) {
    stop("traits should be NULL or a character vector.")
  }
  if (!is.null(keep) && !is.character(keep)) {
    stop("keep should be NULL or a character vector.")
  }
  sapply(X = trials, FUN = function(trial) {
    SSATr <- SSA[[trial]]
    if (is.null(traits)) {
      traitsTr <- SSATr$traits
    } else {
      traitsTr <- traits
    }
    ## Trial specific checks.
    if (!all(traits %in% colnames(SSATr$TD[[trial]]))) {
      stop(paste0("All traits should be columns in ", trial, ".\n"))
    }
    if (!all(keep %in% colnames(SSATr$TD[[trial]]))) {
      stop(paste0("All keep should be columns in ", trial, ".\n"))
    }
    traitsTr <- traitsTr[sapply(X = traitsTr, FUN = function(trait) {
      !is.null(SSATr$mRand[[trait]]) || !is.null(SSATr$mFix[[trait]])})]
    if (length(traitsTr) == 0) {
      return(NULL)
    }
    engine <- SSATr$engine
    ## Set useRepId to TRUE when it is used as fixed effect in the model.
    useRepId <- SSATr$design %in% c("res.ibd", "res.rowcol", "rcbd")
    ## Extract statistics from fitted model.
    result <- do.call(what = paste0("extract", tools::toTitleCase(engine)),
                      args = list(SSA = SSATr, traits = traitsTr, what = what,
                                  useRepId = useRepId, keep = keep,
                                  restore = restoreColNames))
    attr(x = result, which = "traits") <- traitsTr
    attr(x = result, which = "design") <- SSATr$design
    attr(x = result, which = "engine") <- engine
    return(result)
  }, simplify = FALSE)
}

#' Extract statistics from model fitted using SpATS
#'
#' @importFrom SpATS predict.SpATS
#' @keywords internal
extractSpATS <- function(SSA,
                         traits = SSA$traits,
                         what = "all",
                         keep = NULL,
                         useRepId,
                         restore = FALSE) {
  mf <- SSA$mFix[names(SSA$mFix) %in% traits]
  mr <- SSA$mRand[names(SSA$mRand) %in% traits]
  TD <- SSA$TD[[1]]
  renCols <- attr(TD, "renamedCols")
  predicted <- SSA$predicted
  useCheckId <- length(grep(pattern = "checkId",
                            x = deparse(mr[[1]]$model$fixed))) > 0
  whatTot <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "heritability",
               "varCompF", "varCompR", "varGen", "varSpat", "fitted", "resid",
               "stdRes", "rMeans", "ranEf", "residR", "stdResR", "rDf", "rDfR",
               "effDim", "ratEffDim")
  whatMod <- c("F", "F", "R", "R", "R", "F", "R", "R", "R", "F", "F", "F", "R",
               "R", "R", "R", "F", "R", "R", "R")
  whatSSA <- c(if (!is.null(mf)) "F", if (!is.null(mr)) "R")
  whatPred <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "ranEf")
  opts <- Reduce("|", sapply(whatSSA, FUN = grepl, x = whatMod,
                             simplify = FALSE))
  if (what[[1]] == "all") {
    what <- whatTot[opts]
  } else {
    what <- match.arg(arg = what, choices = whatTot[opts], several.ok = TRUE)
  }
  ## Fitted values and residuals are not set to NA by SpATS in case the original
  ## data is NA. Get missing values per trait to set them to NA later.
  naTr <- sapply(X = traits, FUN = function(trait) {
    which(is.na(TD[[trait]]))
  }, simplify = FALSE)
  ## Create baseData and baseDataPred to which further results will be merged.
  base <- createBaseData(TD, predicted, keep, useRepId,
                         bdPred = any(what %in% whatPred))
  baseData <- base$baseData
  baseDataPred <- base$baseDataPred
  ## Create empty result list.
  result <- setNames(vector(mode = "list", length = length(what)),
                     what)
  ## Compute BLUEs and se of BLUEs from fixed model.
  if ("BLUEs" %in% what) {
    predVals <- lapply(X = traits, FUN = function(trait) {
      predVal <- predict(mf[[trait]],
                         which = predicted)[c(predicted, "predicted.values")]
      colnames(predVal) <- c(predicted, trait)
      return(predVal)
    })
    BLUEs <- Reduce(f = function(x, y) merge(x, y, all = TRUE),
                    x = predVals, init = baseDataPred)
    result[["BLUEs"]] <- restoreColNames(renDat = BLUEs, renamedCols = renCols,
                                         restore = restore)
  }
  if ("seBLUEs" %in% what) {
    predErrs <- lapply(X = traits, FUN = function(trait) {
      predErr <- predict(mf[[trait]],
                         which = predicted)[c(predicted, "standard.errors")]
      colnames(predErr) <- c(predicted, trait)
      return(predErr)
    })
    seBLUEs <- Reduce(f = function(x, y) merge(x, y, all = TRUE),
                      x = predErrs, init = baseDataPred)
    result[["seBLUEs"]] <- restoreColNames(renDat = seBLUEs,
                                           renamedCols = renCols,
                                           restore = restore)
  }
  ## Compute BLUPs and se of BLUPs from mixed model.
  if ("BLUPs" %in% what) {
    whichPred <- c(predicted, if (useCheckId) "checkId")
    predVals <- lapply(X = traits, FUN = function(trait) {
      predVal <- predict(mr[[trait]],
                         which = whichPred)[c(whichPred, "predicted.values")]
      colnames(predVal) <- c(whichPred, trait)
      return(predVal)
    })
    BLUPs <- Reduce(f = function(x, y) merge(x, y, all = TRUE),
                    x = predVals, init = baseDataPred)
    result[["BLUPs"]] <- restoreColNames(renDat = BLUPs, renamedCols = renCols,
                                         restore = restore)
  }
  if ("seBLUPs" %in% what) {
    whichPred <- c(predicted, if (useCheckId) "checkId")
    predErrs <- lapply(X = traits, FUN = function(trait) {
      predErr <- predict(mr[[trait]],
                         which = whichPred)[c(whichPred, "standard.errors")]
      colnames(predErr) <- c(whichPred, trait)
      return(predErr)
    })
    seBLUPs <- Reduce(f = function(x, y) merge(x, y, all = TRUE),
                      x = predErrs, init = baseDataPred)
    result[["seBLUPs"]] <- restoreColNames(renDat = seBLUPs,
                                           renamedCols = renCols,
                                           restore = restore)
  }
  ## Compute generalized heritability.
  if ("heritability" %in% what) {
    result[["heritability"]] <- sapply(X = mr, FUN = function(mr0) {
      unname(SpATS::getHeritability(mr0))
    })
  }
  ## Extract variance components for genotype fixed.
  if ("varCompF" %in% what) {
    result[["varCompF"]] <- lapply(X = mf, FUN = extractVarComp,
                                   engine = "SpATS")
  }
  ## Extract variance components for genotype random.
  if ("varCompR" %in% what) {
    result[["varCompR"]] <- lapply(X = mr, FUN = extractVarComp,
                                   engine = "SpATS")
  }
  ## Extract genetic variance.
  if ("varGen" %in% what) {
    result[["varGen"]] <- sapply(X = mr, FUN = function(mr0) {
      unname(mr0$var.comp[predicted])
    })
  }
  ## Extract spatial variance.
  if ("varSpat" %in% what) {
    result[["varSpat"]] <- sapply(X = mr, FUN = function(mr0) {
      mr0$var.comp[grep(pattern = "Coord", x = names(mr0$var.comp))]
    })
  }
  ## Extract fitted values.
  if ("fitted" %in% what) {
    fitVal <- cbind(baseData, sapply(X = traits, FUN = function(trait) {
      fitVals <- fitted(mf[[trait]])
      fitVals[naTr[[trait]]] <- NA
      fitVals
    }))
    result[["fitted"]] <- restoreColNames(renDat = fitVal, renamedCols = renCols,
                                          restore = restore)
  }
  ## Extract residuals.
  if ("resid" %in% what) {
    resVal <- cbind(baseData, sapply(X = traits, FUN = function(trait) {
      resVals <- residuals(mf[[trait]])
      resVals[naTr[[trait]]] <- NA
      resVals
    }))
    result[["resid"]] <- restoreColNames(renDat = resVal, renamedCols = renCols,
                                         restore = restore)
  }
  ## Extract standardized residuals.
  if ("stdRes" %in% what) {
    stdRes <- cbind(baseData, sapply(X = mf, FUN = function(mf0) {
      residuals(mf0) / sd(residuals(mf0), na.rm = TRUE)
    }))
    result[["stdRes"]] <- restoreColNames(renDat = stdRes, renamedCols = renCols,
                                          restore = restore)
  }
  ## Extract rMeans.
  if ("rMeans" %in% what) {
    rMeans <- cbind(baseData, sapply(X = traits, FUN = function(trait) {
      fitVals <- fitted(mr[[trait]])
      fitVals[naTr[[trait]]] <- NA
      fitVals
    }))
    result[["rMeans"]] <- restoreColNames(renDat = rMeans, renamedCols = renCols,
                                          restore = restore)
  }
  ## Extract random effects.
  if ("ranEf" %in% what) {
    ranEffs <- lapply(X = traits, FUN = function(trait) {
      effs <- mr[[trait]]$coeff[names(mr[[trait]]$coeff) %in%
                                  baseDataPred[[predicted]]]
      ranEff <- data.frame(names(effs), effs)
      colnames(ranEff) <- c(predicted, trait)
      return(ranEff)
    })
    ranEf <- Reduce(f = function(x, y) merge(x, y, all = TRUE),
                    x = ranEffs, init = baseDataPred)
    result[["ranEf"]] <- restoreColNames(renDat = ranEf, renamedCols = renCols,
                                         restore = restore)
  }
  ## Extract residuals for genotype random.
  if ("residR" %in% what) {
    resVal <- cbind(baseData, sapply(X = traits, FUN = function(trait) {
      resVals <- residuals(mr[[trait]])
      resVals[naTr[[trait]]] <- NA
      resVals
    }))
    result[["residR"]] <- restoreColNames(renDat = resVal,
                                          renamedCols = renCols,
                                          restore = restore)
  }
  ## Extract standardized residuals for genotype random.
  if ("stdResR" %in% what) {
    stdRes <- cbind(baseData, sapply(X = mr, FUN = function(mr0) {
      residuals(mr0) / sd(residuals(mr0), na.rm = TRUE)
    }))
    result[["stdResR"]] <- restoreColNames(renDat = stdRes,
                                           renamedCols = renCols,
                                           restore = restore)
  }
  ## Extract residual degrees of freedom.
  if ("rDf" %in% what) {
    result[["rDf"]] <- sapply(X = mf, FUN = function(mf0) {
      round(mf0[["nobs"]] - sum(mf0[["eff.dim"]]))
    })
  }
  ## Extract residual degrees of freedom.
  if ("rDfR" %in% what) {
    result[["rDfR"]] <- sapply(X = mr, FUN = function(mr0) {
      round(mr0[["nobs"]] - sum(mr0[["eff.dim"]]))
    })
  }
  ## Extract effective dimensions.
  if ("effDim" %in% what) {
    result[["effDim"]] <- sapply(X = mr, FUN = `[[`, "eff.dim")
  }
  ## Extract ratio's of effective dimensions.
  if ("ratEffDim" %in% what) {
    result[["ratEffDim"]] <- sapply(X = mr, FUN = function(mr0) {
      capture.output(ratTot <- summary(mr0)$p.table.dim[, "Ratio"])
      setNames(as.numeric(ratTot[1:(length(ratTot) - 4)]),
               names(ratTot[1:(length(ratTot) - 4)]))
    })
  }
  return(result)
}

#' Extract statistics from model fitted using lme4
#'
#' @keywords internal
extractLme4 <- function(SSA,
                        traits = SSA$traits,
                        what = "all",
                        keep = NULL,
                        useRepId,
                        restore = FALSE) {
  mf <- SSA$mFix[names(SSA$mFix) %in% traits]
  mr <- SSA$mRand[names(SSA$mRand) %in% traits]
  TD <- SSA$TD[[1]]
  renCols <- attr(TD, "renamedCols")
  predicted = SSA$predicted
  whatTot <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "ue", "heritability",
               "varCompF", "varCompR", "varGen", "varErr", "fitted", "resid",
               "stdRes", "rMeans", "ranEf", "residR", "stdResR", "wald", "CV",
               "rDf", "rDfR")
  whatMod <- c("F", "F", "R", "R", "F", "R", "F", "R", "R", "R", "F", "F", "F",
               "R", "R", "R", "R", "F", "F", "F", "R")
  whatSSA <- c(if (!is.null(mf)) "F", if (!is.null(mr)) "R")
  whatPred <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "ue", "ranEf")
  opts <- Reduce("|", sapply(whatSSA, FUN = grepl, x = whatMod,
                             simplify = FALSE))
  if (what[[1]] == "all") {
    what <- whatTot[opts]
  } else {
    what <- match.arg(arg = what, choices = whatTot[opts], several.ok = TRUE)
  }
  ## Create baseData and baseDataPred to which further results will be merged.
  base <- createBaseData(TD, predicted, keep, useRepId,
                         bdPred = any(what %in% whatPred))
  baseData <- base$baseData
  baseDataPred <- base$baseDataPred
  ## Create empty result list.
  result <- setNames(vector(mode = "list", length = length(what)),
                     what)
  ## Use emmeans to convert mf to emmGrid per trait.
  ## Option lmer.df = "asymptotic" for speed. Default options gives
  ## exact the same results as asreml but is very slow.
  em <- sapply(X = mf, FUN = function(mf0) {
    emmeans::emmeans(mf0, specs = predicted, lmer.df = "asymptotic")
  }, simplify = FALSE)
  ## Summarize emGrid to calculate BLUEs en se per trait.
  emStats <- sapply(X = em, FUN = summary, simplify = FALSE)
  ## Extract BLUEs and se of BLUEs from emStats.
  if ("BLUEs" %in% what) {
    predVals <- lapply(X = traits, FUN = function(trait) {
      setNames(emStats[[trait]][c(predicted, "emmean")],
               c(predicted, trait))
    })
    BLUEs <- Reduce(f = function(x, y) merge(x, y, all = TRUE),
                    x = predVals, init = baseDataPred)
    result[["BLUEs"]] <- restoreColNames(renDat = BLUEs, renamedCols = renCols,
                                         restore = restore)
  }
  if ("seBLUEs" %in% what) {
    predErrs <- lapply(X = traits, FUN = function(trait) {
      setNames(emStats[[trait]][c(predicted, "SE")],
               c(predicted, trait))
    })
    seBLUEs <- Reduce(f = function(x, y) merge(x, y, all = TRUE),
                      x = predErrs, init = baseDataPred)
    result[["seBLUEs"]] <- restoreColNames(renDat = seBLUEs,
                                           renamedCols = renCols,
                                           restore = restore)
  }
  ## Compute BLUPs and se of BLUPs from mixed model.
  if ("BLUPs" %in% what) {
    predVals <- lapply(X = traits, FUN = function(trait) {
      ## Extract fixed effects.
      fixEfMr <- lme4::fixef(mr[[trait]])
      ## Get positions of 'repId' within fixed effects.
      repIdPos <- grep(pattern = "repId",
                       x = names(fixEfMr))
      repIdMix <- fixEfMr[repIdPos]
      ## Compute Blups.
      coefs <- coef(mr[[trait]])[[predicted]]
      predVal <- data.frame(rownames(coefs),
                            coefs[, "(Intercept)"] + mean(c(repIdMix, 0)))
      colnames(predVal) <- c(predicted, trait)
      return(predVal)
    })
    BLUPs <- Reduce(f = function(x, y) merge(x, y, all = TRUE),
                    x = predVals, init = baseDataPred)
    result[["BLUPs"]] <- restoreColNames(renDat = BLUPs, renamedCols = renCols,
                                         restore = restore)
  }
  if ("seBLUPs" %in% what) {
    predErrs <- lapply(X = traits, FUN = function(trait) {
      ranEffs <- lme4::ranef(mr[[trait]], condVar = TRUE)[[predicted]]
      predErr <- data.frame(rownames(ranEffs),
                            as.vector(sqrt(attr(ranEffs, "postVar"))))
      colnames(predErr) <- c(predicted, trait)
      return(predErr)
    })
    seBLUPs <- Reduce(f = function(x, y) merge(x, y, all = TRUE),
                      x = predErrs, init = baseDataPred)
    result[["seBLUPs"]] <- restoreColNames(renDat = seBLUPs,
                                           renamedCols = renCols,
                                           restore = restore)
  }
  ## Compute unit errors.
  if ("ue" %in% what) {
    unitErrs <- lapply(X = traits, FUN = function(trait) {
      ## Extract and invert variance covariance matrix.
      V <- vcov(em[[trait]])
      Vinv <- try(chol2inv(chol(V)), silent = TRUE)
      ## Compute unit errors.
      if (!inherits(Vinv, "try-error")) {
        ue <- 1 / diag(Vinv)
      } else {
        ue <- 1 / diag(solve(V))
      }
      unitErr <- data.frame(levels(em[[trait]]), ue)
      colnames(unitErr) <- c(predicted, trait)
      return(unitErr)
    })
    ue <- Reduce(f = function(x, y) merge(x, y, all = TRUE),
                 x = unitErrs, init = baseDataPred)
    result[["ue"]] <- restoreColNames(renDat = ue, renamedCols = renCols,
                                      restore = restore)
  }
  ## Extract variance components for genotype fixed.
  if ("varCompF" %in% what) {
    result[["varCompF"]] <- lapply(X = mf, FUN = extractVarComp,
                                   engine = "lme4")
  }
  ## Extract variance components for genotype random.
  if ("varCompR" %in% what) {
    result[["varCompR"]] <- lapply(X = mr, FUN = extractVarComp,
                                   engine = "lme4")
  }
  ## Extract variances.
  if (any(c("varGen", "varErr", "heritability") %in% what)) {
    varCor <- lapply(X = mr, FUN = lme4::VarCorr)
    varGen <- sapply(X = varCor, FUN = function(vc) {
      vc[[predicted]][1, 1]
    })
    varErr <- sapply(X = varCor, FUN = "attr", "sc") ^ 2
    if ("varGen" %in% what) {
      result[["varGen"]] <- varGen
    }
    if ("varErr" %in% what) {
      result[["varErr"]] <- varErr
    }
    if ("heritability" %in% what) {
      ## Estimate heritability on a line mean basis.
      if (useRepId) {
        result[["heritability"]] <- varGen /
          (varGen + (varErr / length(unique(TD$repId))))
      } else {
        result[["heritability"]] <- varGen / (varGen + varErr)
      }
    }
  }
  ## Extract fitted values.
  if ("fitted" %in% what) {
    fitVal <- cbind(baseData, sapply(X = mf, FUN = fitted))
    result[["fitted"]] <- restoreColNames(renDat = fitVal, renamedCols = renCols,
                                          restore = restore)
  }
  ## Extract residuals.
  if ("resid" %in% what) {
    resVal <- cbind(baseData, sapply(X = mf, FUN = residuals))
    result[["resid"]] <- restoreColNames(renDat = resVal, renamedCols = renCols,
                                         restore = restore)
  }
  ## Extract standardized residuals.
  if ("stdRes" %in% what) {
    stdRes <- cbind(baseData,
                    sapply(X = mf, FUN = function(mf0) {
                      if (inherits(mf0, "lm")) {
                        stdRes <- rstandard(mf0)
                      } else if (inherits(mf0, "lmerMod")) {
                        stdRes <- residuals(mf0, scaled = TRUE)
                      }
                    }))
    result[["stdRes"]] <- restoreColNames(renDat = stdRes, renamedCols = renCols,
                                          restore = restore)
  }
  ## Compute rMeans.
  ## Use napredict to fill in NAs in data with NAs.
  if ("rMeans" %in% what) {
    rMeans <- cbind(baseData,
                    sapply(X = mr, FUN = function(mr0) {
                      napredict(attr(model.frame(mr0), "na.action"),
                                x = lme4::getME(mr0, "mu"))
                    }))
    result[["rMeans"]] <- restoreColNames(renDat = rMeans, renamedCols = renCols,
                                          restore = restore)
  }
  ## Extract random effects.
  if ("ranEf" %in% what) {
    ranEffs <- lapply(X = traits, FUN = function(trait) {
      effs <- lme4::ranef(mr[[trait]], drop = TRUE)[[predicted]]
      ranEff <- data.frame(names(effs), effs)
      colnames(ranEff) <- c(predicted, trait)
      return(ranEff)
    })
    ranEf <- Reduce(f = function(x, y) merge(x, y, all = TRUE),
                    x = ranEffs, init = baseDataPred)
    result[["ranEf"]] <- restoreColNames(renDat = ranEf, renamedCols = renCols,
                                         restore = restore)
  }
  ## Extract residuals for genotype random.
  if ("residR" %in% what) {
    resVal <- cbind(baseData, sapply(X = mr, FUN = residuals))
    result[["residR"]] <- restoreColNames(renDat = resVal,
                                          renamedCols = renCols,
                                          restore = restore)
  }
  ## Extract standardized residuals for genotype random.
  if ("stdResR" %in% what) {
    stdRes <- cbind(baseData,
                    sapply(X = mr, FUN = function(mr0) {
                      if (inherits(mr0, "lm")) {
                        stdRes <- rstandard(mr0)
                      } else if (inherits(mr0, "lmerMod")) {
                        stdRes <- residuals(mr0, scaled = TRUE)
                      }
                    }))
    result[["stdResR"]] <- restoreColNames(renDat = stdRes,
                                           renamedCols = renCols,
                                           restore = restore)
  }
  ## Compute wald test.
  if ("wald" %in% what) {
    result[["wald"]] <- lapply(X = em, FUN = emmeans::test, joint = TRUE)
  }
  if ("CV" %in% what) {
    ## Compute Coefficient of Variation.
    result[["CV"]] <- sapply(X = mf, FUN = function(mf0) {
      100 * sigma(mf0) / mean(fitted(mf0), na.rm = TRUE)
    })
  }
  if ("rDf" %in% what) {
    result[["rDf"]] <- sapply(X = mf, FUN = df.residual)
  }
  if ("rDfR" %in% what) {
    result[["rDfR"]] <- sapply(X = mr, FUN = df.residual)
  }
  return(result)
}

#' Extract statistics from model fitted using asreml
#'
#' @keywords internal
extractAsreml <- function(SSA,
                          traits = SSA$traits,
                          what = "all",
                          keep = NULL,
                          useRepId,
                          restore = FALSE) {
  if (!requireNamespace("asreml", quietly = TRUE)) {
    stop("asreml cannot be successfully loaded.\n")
  }
  mf <- SSA$mFix[names(SSA$mFix) %in% traits]
  mr <- SSA$mRand[names(SSA$mRand) %in% traits]
  TD <- SSA$TD[[1]]
  renCols <- attr(TD, "renamedCols")
  predicted <- SSA$predicted
  whatTot <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "ue", "heritability",
               "varCompF", "varCompR", "varGen", "varErr", "fitted", "resid",
               "stdRes", "rMeans", "ranEf", "residR", "stdResR", "wald", "CV",
               "rDf", "rDfR", "sed", "lsd")
  whatMod <- c("F", "F", "R", "R", "F", "R", "F", "R", "R", "R", "F", "F", "F",
               "R", "R", "R", "R", "F", "F", "F", "R", "F", "F")
  whatSSA <- c(if (!is.null(mf)) "F", if (!is.null(mr)) "R")
  whatPred <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "ue", "ranEf")
  opts <- Reduce("|", sapply(whatSSA, FUN = grepl, x = whatMod,
                             simplify = FALSE))
  if (what[[1]] == "all") {
    what <- whatTot[opts]
  } else {
    what <- match.arg(arg = what, choices = whatTot[opts], several.ok = TRUE)
  }
  ## Create baseData and baseDataPred to which further results will be merged.
  base <- createBaseData(TD, predicted, keep, useRepId,
                         bdPred = any(what %in% whatPred))
  baseData <- base$baseData
  baseDataPred <- base$baseDataPred
  ## Create empty result list.
  result <- setNames(vector(mode = "list", length = length(what)), what)
  ## Construct assocForm for use in associate in predict.
  if (length(grep(pattern = "+ checkId +", x = getCall(mf[[1]]))) > 0) {
    assocForm <- formula("~ checkId:genotype")
  } else {
    assocForm <- formula("~ NULL")
  }
  ## Extract BLUEs and se of BLUEs from fixed model.
  if ("BLUEs" %in% what) {
    ## asreml3 saves the predictions inside the asreml object.
    ## asreml4 creates a list containing nothing but the predictions.
    predVals <- lapply(X = traits, FUN = function(trait) {
      mfPred <- predictAsreml(mf[[trait]], TD = TD, associate = assocForm)
      setNames(if (asreml4()) {
        mfPred$pvals[c(predicted, "predicted.value")]
      } else {
        mfPred$predictions$pvals[c(predicted, "predicted.value")]
      }, c(predicted, trait))
    })
    BLUEs <- Reduce(f = merge, x = predVals, init = baseDataPred)
    result[["BLUEs"]] <- restoreColNames(renDat = BLUEs, renamedCols = renCols,
                                         restore = restore)
  }
  if ("seBLUEs" %in% what) {
    ## asreml3 saves the predictions inside the asreml object.
    ## asreml4 creates a list containing nothing but the predictions.
    predErrs <- lapply(X = traits, FUN = function(trait) {
      mfPred <- predictAsreml(mf[[trait]], TD = TD, associate = assocForm)
      setNames(if (asreml4()) {
        mfPred$pvals[c(predicted, "std.error")]
      } else {
        mfPred$predictions$pvals[c(predicted, "standard.error")]
      }, c(predicted, trait))
    })
    seBLUEs <- Reduce(f = merge, x = predErrs, init = baseDataPred)
    result[["seBLUEs"]] <- restoreColNames(renDat = seBLUEs,
                                           renamedCols = renCols,
                                           restore = restore)
  }
  ## Extract BLUPs and se of BLUPs from fixed model.
  if ("BLUPs" %in% what) {
    ## asreml3 saves the predictions inside the asreml object.
    ## asreml4 creates a list containing nothing but the predictions.
    predVals <- lapply(X = traits, FUN = function(trait) {
      mrPred <- predictAsreml(mr[[trait]], TD = TD)
      setNames(if (asreml4()) {
        mrPred$pvals[c(predicted, "predicted.value")]
      } else {
        mrPred$predictions$pvals[c(predicted, "predicted.value")]
      }, c(predicted, trait))
    })
    BLUPs <- Reduce(f = merge, x = predVals, init = baseDataPred)
    result[["BLUPs"]] <- restoreColNames(renDat = BLUPs, renamedCols = renCols,
                                         restore = restore)
  }
  if ("seBLUPs" %in% what) {
    ## asreml3 saves the predictions inside the asreml object.
    ## asreml4 creates a list containing nothing but the predictions.
    predErrs <- lapply(X = traits, FUN = function(trait) {
      mrPred <- predictAsreml(mr[[trait]], TD = TD)
      setNames(if (asreml4()) {
        mrPred$pvals[c(predicted, "std.error")]
      } else {
        mrPred$predictions$pvals[c(predicted, "standard.error")]
      }, c(predicted, trait))
    })
    seBLUPs <- Reduce(f = merge, x = predErrs, init = baseDataPred)
    result[["seBLUPs"]] <- restoreColNames(renDat = seBLUPs, renamedCols = renCols,
                                           restore = restore)
  }
  ## Compute unit errors.
  if ("ue" %in% what) {
    ue <- cbind(baseDataPred, sapply(X = mf, FUN = function(mf0) {
      ## asreml3 saves the predictions inside the asreml object.
      ## asreml4 creates a list containing nothing but the predictions.
      mfPred <- predictAsreml(mf0, TD = TD)
      ## Extract V from mf.
      V <- if (asreml4()) {
        mfPred$vcov
      } else {
        mfPred$predictions$vcov
      }
      ## Remove columns and rows containing NA.
      VMiss <- apply(X = V, MARGIN = 2, FUN = anyNA)
      V <- V[!VMiss, !VMiss]
      Vinv <- try(chol2inv(chol(V)), silent = TRUE)
      ## Compute unit errors.
      ue <- rep(x = NA, times = nrow(baseDataPred))
      if (!inherits(Vinv, "try-error")) {
        ue[!VMiss] <- 1 / diag(Vinv)
      } else {
        ue[!VMiss] <- 1 / diag(solve(V))
      }
      return(ue)
    }))
    result[["ue"]] <- restoreColNames(renDat = ue, renamedCols = renCols,
                                      restore = restore)
  }
  ## Extract variance components for genotype fixed.
  if ("varCompF" %in% what) {
    result[["varCompF"]] <- lapply(X = mf, FUN = extractVarComp,
                                   engine = "asreml")
  }
  ## Extract variance components for genotype random.
  if ("varCompR" %in% what) {
    result[["varCompR"]] <- lapply(X = mr, FUN = extractVarComp,
                                   engine = "asreml")
  }
  ## Extract variances.
  ## In asreml3 variances are stored in gammas, in asreml4 in vparameters.
  ## Also the naming is different.
  varGen <- sapply(X = mr, FUN = function(mr0) {
    if (asreml4()) {
      unname(mr0$vparameters[predicted] * mr0$sigma2)
    } else {
      unname(mr0$gammas[pattern = paste0(predicted, "!", predicted, ".var")] *
               mr0$sigma2)
    }
  })
  if ("varGen" %in% what) {
    result[["varGen"]] <- varGen
  }
  if ("varErr" %in% what) {
    result[["varErr"]] <- sapply(X = mr, FUN = function(mr0) {
      if (asreml4()) {
        unname(mr0$vparameters["units!R"] * mr0$sigma2)
      } else {
        unname(mr0$gammas["R!variance"] * mr0$sigma2)
      }
    })
  }
  ## Estimate heritability on a line mean basis.
  ## asreml3 saves the predictions inside the asreml object.
  ## asreml4 creates a list containing nothing but the predictions.
  if ("heritability" %in% what) {
    result[["heritability"]] <- sapply(X = traits, FUN = function(trait) {
      mrPred <- predictAsreml(model = mr[[trait]], classify = predicted,
                              vcov = FALSE, TD = TD, only = predicted,
                              sed = TRUE)
      sedSq <- if (asreml4()) {
        mrPred$sed ^ 2
      } else {
        mrPred$predictions$sed ^ 2
      }
      unname(1 - mean(sedSq[lower.tri(sedSq)]) / (2 * varGen[[trait]]))
    })
  }
  ## Extract fitted values.
  if ("fitted" %in% what) {
    fitVal <- cbind(baseData, sapply(X = mf, FUN = fitted))
    result[["fitted"]] <- restoreColNames(renDat = fitVal, renamedCols = renCols,
                                          restore = restore)
  }
  ## Extract residuals.
  if ("resid" %in% what) {
    resVal <- cbind(baseData, sapply(X = mf, FUN = residuals,
                                     type = "response"))
    result[["resid"]] <- restoreColNames(renDat = resVal, renamedCols = renCols,
                                         restore = restore)
  }
  ## Extract standardized residuals.
  if ("stdRes" %in% what) {
    stdRes <- cbind(baseData, sapply(X = mf, FUN = residuals, type = "stdCond"))
    result[["stdRes"]] <- restoreColNames(renDat = stdRes, renamedCols = renCols,
                                          restore = restore)
  }
  ## Extract rMeans.
  if ("rMeans" %in% what) {
    rMeans <- cbind(baseData, sapply(X = mr, FUN = fitted))
    result[["rMeans"]] <- restoreColNames(renDat = rMeans, renamedCols = renCols,
                                          restore = restore)
  }
  ## Extract random effects.
  if ("ranEf" %in% what) {
    ranEffs <- lapply(X = traits, FUN = function(trait) {
      coefs <- mr[[trait]]$coe$random
      ## In asreml3 coefficients are stored as a vector,
      ## in asreml4 as a data.frame.
      ranEff <- if (asreml4()) {
        data.frame(gsub(pattern = paste0(predicted, "_"), replacement = "",
                        x = rownames(coefs)[grep(pattern = predicted,
                                                 x = rownames(coefs))]),
                   coefs[grep(pattern = predicted, x = rownames(coefs))])
      } else {
        data.frame(gsub(pattern = paste0(predicted, "_"), replacement = "",
                        x = names(coefs)[grep(pattern = predicted,
                                              x = names(coefs))]),
                   coefs[grep(pattern = predicted, x = names(coefs))])
      }
      colnames(ranEff) <- c(predicted, trait)
      return(ranEff)
    })
    ranEf <- Reduce(f = merge, x = ranEffs, init = baseDataPred)
    result[["ranEf"]] <- restoreColNames(renDat = ranEf, renamedCols = renCols,
                                         restore = restore)
  }
  ## Extract residuals for genotype random.
  if ("residR" %in% what) {
    resVal <- cbind(baseData, sapply(X = mr, FUN = residuals,
                                     type = "response"))
    result[["residR"]] <- restoreColNames(renDat = resVal,
                                          renamedCols = renCols,
                                          restore = restore)
  }
  ## Extract standardized residuals.
  if ("stdResR" %in% what) {
    stdRes <- cbind(baseData, sapply(X = mr, FUN = residuals, type = "stdCond"))
    result[["stdResR"]] <- restoreColNames(renDat = stdRes,
                                           renamedCols = renCols,
                                           restore = restore)
  }
  ## Compute wald test.
  if ("wald" %in% what) {
    tmpfile <- tempfile()
    sink(file = tmpfile)
    result[["wald"]] <- lapply(X = mf, function(mf0) {
      wtt <- asreml::wald.asreml(mf0, ssType = "conditional", denDF = "numeric",
                                 maxiter = 200, data = TD)
      pos <- grep(pattern = predicted, x = row.names(wtt$Wald))
      chi2 <- wtt$Wald$F.con[pos] * wtt$Wald$Df[pos]
      prob <- 1 - pchisq(q = chi2, df = wtt$Wald$Df[pos])
      list(chi2 = c(chi2 = chi2, df = wtt$Wald$Df[pos], P = prob),
           Ftest = c(Fstat = wtt$Wald$F.con[pos],
                     df1 = wtt$Wald$Df[pos],
                     df2 = wtt$Wald$denDF[pos],
                     P = wtt$Wald$Pr[pos]))
    })
    sink()
    unlink(tmpfile)
  }
  ## Compute Coefficient of Variation.
  if ("CV" %in% what) {
    result[["CV"]] <- sapply(X = mf, function(mf0) {
      100 * summary(mf0)$sigma / mean(fitted(mf0), na.rm = TRUE)
    })
  }
  ## Extract residual degrees of freedom.
  if ("rDf" %in% what) {
    result[["rDf"]] <- sapply(X = mf, FUN = "[[", "nedf")
  }
  ## Extract residual degrees of freedom for genotype random.
  if ("rDfR" %in% what) {
    result[["rDfR"]] <- sapply(X = mr, FUN = "[[", "nedf")
  }
  ## Extract standard error of difference.
  ## asreml3 saves the predictions inside the asreml object.
  ## asreml4 creates a list containing nothing but the predictions.
  if ("sed" %in% what) {
    result[["sed"]] <- lapply(X = mf, FUN = function(mf0) {
      mfPred <- predictAsreml(mf0, TD = TD)
      if (asreml4()) {
        mfPred$avsed
      } else {
        mfPred$predictions$avsed
      }
    })
  }
  ## Compute lsd with signifcance level 5%.
  ## asreml3 saves the predictions inside the asreml object.
  ## asreml4 creates a list containing nothing but the predictions.
  if ("lsd" %in% what) {
    result[["lsd"]] <- lapply(X = mf, FUN = function(mf0) {
      mfPred <- predictAsreml(mf0, TD = TD)
      if (asreml4()) {
        qt(p = .975, df = mf0$nedf) * mfPred$avsed
      } else {
        qt(p = .975, df = mf0$nedf) * mfPred$predictions$avsed
      }
    })
  }
  return(result)
}

#' Helper function for creating baseData
#'
#' @keywords internal
createBaseData <- function(TD,
                           predicted,
                           keep = NULL,
                           useRepId = FALSE,
                           bdPred = FALSE) {
  ## Create baseData consisting of predicted variable, possibly repId and
  ## selected keep columns.
  baseData <- TD[, colnames(TD) %in% c(predicted, ifelse(useRepId, "repId", ""),
                                       keep), drop = FALSE]
  ## Create baseData for predictions with predicted variable(s).
  if (bdPred) {
    baseDataPred <- unique(TD[predicted])
    ## Add columns in keep one-by-one. Only data that is constant within
    ## predicted is actually kept. Other columns are dropped with a warning.
    for (col in keep) {
      TDKeep <- unique(TD[, c(predicted, col)])
      if (!anyDuplicated(TDKeep[, predicted])) {
        baseDataPred <- merge(baseDataPred, TDKeep, by = predicted)
      } else {
        warning(paste0("Duplicate values for ", deparse(substitute(col)), ". ",
                       "Column dropped.\n"), call. = FALSE)
      }
    }
  } else {
    baseDataPred <- NULL
  }
  return(list(baseData = baseData, baseDataPred = baseDataPred))
}

#' Helper function for adding back original colnames
#'
#' @keywords internal
restoreColNames <- function(renDat,
                            renamedCols,
                            restore = FALSE) {
  if (restore && !is.null(renamedCols)) {
    renCols <- colnames(renDat)
    ## Get original columnnames from renamedCols data.frame.
    origCols <- sapply(X = renCols, FUN = function(renCol) {
      if (renCol %in% renamedCols$new) {
        renamedCols$orig[renCol == renamedCols$new]
      } else {
        ## If no renaming took place keep current name.
        renCol
      }
    })
    colnames(renDat) <- origCols
    ## Columns might be duplicated now because one column was renamed to multiple
    ## new columns. Remove those duplicated columns.
    ## Set column to NULL to prevent attributes from being dropped.
    renDat[duplicated(colnames(renDat))] <- NULL
  }
  return(renDat)
}
