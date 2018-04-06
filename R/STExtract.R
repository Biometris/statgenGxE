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
#' \item{R - varGen}{Genetic variance component.}
#' \item{R - varErr}{Residual variance component - only for \code{lme4}
#' and \code{asreml}.}
#' \item{R - varSpat}{Spatial variance components - only for \code{SpATS}.}
#' \item{F - fitted}{Fitted values for the model with genotype as fixed
#' component.}
#' \item{F - resid}{Residuals for the model with genotype as fixed component.}
#' \item{F - stdRes}{Standardized residuals for the model with genotype as fixed
#' component - only for \code{lme4} and \code{asreml}.}
#' \item{R - rMeans}{Fitted values for the model with genotype as random
#' component.}
#' \item{R - ranEf}{Random genetic effects.}
#' \item{F - wald}{Results of the wald test - only for \code{lme4} and
#' \code{asreml}.}
#' \item{F - CV}{Coefficient of variation - only for \code{lme4} and
#' \code{asreml}.}
#' \item{F - rDf}{Residual degrees of freedom.}
#' \item{R - effDim}{Effective dimensions - only for \code{SpATS}.}
#' \item{F - sed}{Standard error of difference - only for \code{asreml}.}
#' \item{F - lsd}{Least significant difference - only for \code{asreml}.}
#' \item{all}{All available statistics.}
#' }
#'
#' @param SSA An object of class SSA.
#' @param trials A character vector of trials for which the statistics should be
#' computed. If not supplied statistics are computed for all trials that have
#' been modelled.
#' @param traits A character vector of traits for which the statistics should be
#' computed. If not supplied statistics are computed for all traits that have
#' been modelled.
#' @param what A character vector indicating which statistics should be
#' computed. Most statistics are available for all models, some only for models
#' fitted using a certain engine. If this is the case this is indicated in the
#' list with options in details.\cr
#' If \code{what = "all"} all available statistics are computed.
#' @param keep A character vector of column(s) in the object of class
#' \code{\link{TD}} used for modeling. These columns will be kept as output when
#' computing fitted values, residuals, standardized residuals and rMeans.
#' Columns can also be kept when computing (se)BLUEs and (se)BLUPs but only if
#' the column to keep contains unique values for the modeled variables, i.e. a
#' column repId with several different values per genotype cannot be kept.
#'
#' @return A list with per trial for which statistics have been extracted either
#' a list of those statistics or if only one statistic is extracted
#' a single object containing this statistic.
#'
#' @seealso
#' \code{\link{STRunModel}}, \code{\link{STModSpATS}}, \code{\link{STModLme4}}
#' and \code{\link{STModAsreml}}
#'
#' @examples
#' ## Fit model using SpATS.
#' myModel <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield")
#' ## Extract all available statistics from the fitted model.
#' extr <- STExtract(myModel)
#' ## Extract only the BLUEs from the fitted model.
#' BLUEs <- STExtract(myModel, what = "BLUEs")
#' ## Extract only the BLUEs from the fitted model and keep trial as variable in
#' ## the output.
#' BLUEs2 <- STExtract(myModel, what = "BLUEs", keep = "trial")
#'
#' @export
STExtract <- function(SSA,
                      trials = names(SSA),
                      traits = NULL,
                      what = "all",
                      keep = NULL) {
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
      traits <- SSATr$traits
    }
    if (!all(traits %in% colnames(SSATr$TD[[trial]]))) {
      stop(paste0("All traits should be columns in ", trial, ".\n"))
    }
    if (!all(keep %in% colnames(SSATr$TD[[trial]]))) {
      stop(paste0("All keep should be columns in ", trial, ".\n"))
    }
    engine <- SSATr$engine
    ## Set useRepId to TRUE when it is used as fixed effect in the model.
    useRepId <- (SSATr$design %in% c("res.ibd", "res.rowcol", "rcbd"))
    ## Extract statistics from fitted model.
    result <- do.call(what = paste0("extract", tools::toTitleCase(engine)),
                      args = list(SSA = SSATr, traits = traits, what = what,
                                  useRepId = useRepId, keep = keep))
    attr(x = result, which = "traits") <- traits
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
                         useRepId) {
  mf <- SSA$mFix
  mr <- SSA$mRand
  TD <- SSA$TD[[1]]
  predicted <- SSA$predicted
  useCheckId <- length(grep(pattern = "checkId",
                            x = deparse(mr[[1]]$model$fixed))) > 0
  whatTot <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "heritability", "varGen",
               "varSpat", "fitted", "resid", "rMeans", "ranEf", "rDf", "effDim")
  whatMod <- c("F", "F", "R", "R", "R", "R", "R", "F", "F", "R", "R", "F", "F")
  whatSSA <- c(if (!is.null(mf)) "F", if (!is.null(mr)) "R")
  whatPred <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "ranEf")
  if (what[[1]] == "all") {
    what <- whatTot[whatMod %in% whatSSA]
  } else {
    what <- match.arg(arg = what, choices = whatTot[whatMod %in% whatSSA],
                      several.ok = TRUE)
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
    result[["BLUEs"]] <- BLUEs
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
    result[["seBLUEs"]] <- seBLUEs
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
    result[["BLUPs"]] <- BLUPs
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
    result[["seBLUPs"]] <- seBLUPs
  }
  ## Compute generalized heritability.
  if ("heritability" %in% what) {
    result[["heritability"]] <- sapply(X = mr, FUN = function(mr0) {
      unname(SpATS::getHeritability(mr0))
    })
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
      mr0$var.comp[grep(pattern = "Coordinates", x = names(mr0$var.comp))]
    })
  }
  ## Extract fitted values.
  if ("fitted" %in% what) {
    fitVal <- cbind(baseData, sapply(X = traits, FUN = function(trait) {
      fitVals <- fitted(mf[[trait]])
      fitVals[naTr[[trait]]] <- NA
      fitVals
    }))
    result[["fitted"]] <- fitVal
  }
  ## Extract residuals.
  if ("resid" %in% what) {
    resVal <- cbind(baseData, sapply(X = traits, FUN = function(trait) {
      resVals <- residuals(mf[[trait]])
      resVals[naTr[[trait]]] <- NA
      resVals
    }))
    result[["resid"]] <- resVal
  }
  ## Extract rMeans.
  if ("rMeans" %in% what) {
    rMeans <- cbind(baseData, sapply(X = traits, FUN = function(trait) {
      fitVals <- fitted(mr[[trait]])
      fitVals[naTr[[trait]]] <- NA
      fitVals
    }))
    result[["rMeans"]] <- rMeans
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
    result[["ranEf"]] <- ranEf
  }
  ## Extract residual degrees of freedom.
  if ("rDf" %in% what) {
    result[["rDf"]] <- sapply(X = mf, FUN = function(mf0) {
      unname(mf0$dim[predicted])
    })
  }
  ## Extract effective dimensions.
  if ("effDim" %in% what) {
    result[["effDim"]] <- sapply(X = mr, FUN = function(mr0) {
      mr0$eff.dim
    })
  }
  if (length(result) == 1) {
    result <- result[[1]]
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
                        useRepId) {
  mf <- SSA$mFix
  mr <- SSA$mRand
  TD <- SSA$TD[[1]]
  predicted = SSA$predicted
  whatTot <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "ue", "heritability",
               "varGen", "varErr", "fitted", "resid", "stdRes", "rMeans", "ranEf",
               "wald", "CV", "rDf")
  whatMod <- c("F", "F", "R", "R", "F", "R", "R", "R", "F", "F", "F", "R", "R",
               "F", "F", "F")
  whatSSA <- c(if (!is.null(mf)) "F", if (!is.null(mr)) "R")
  whatPred <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "ue", "ranEf")
  if (what[[1]] == "all") {
    what <- whatTot[whatMod %in% whatSSA]
  } else {
    what <- match.arg(arg = what, choices = whatTot[whatMod %in% whatSSA],
                      several.ok = TRUE)
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
    result[["BLUEs"]] <- BLUEs
  }
  if ("seBLUEs" %in% what) {
    predErrs <- lapply(X = traits, FUN = function(trait) {
      setNames(emStats[[trait]][c(predicted, "SE")],
               c(predicted, trait))
    })
    seBLUEs <- Reduce(f = function(x, y) merge(x, y, all = TRUE),
                      x = predErrs, init = baseDataPred)
    result[["seBLUEs"]] <- seBLUEs
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
    result[["BLUPs"]] <- BLUPs
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
    result[["seBLUPs"]] <- seBLUPs
  }
  ## Compute unit errors.
  if ("ue" %in% what) {
    ue <- cbind(baseDataPred,
                lapply(X = em, FUN = function(em0) {
                  ## Extract and invert variance covariance matrix.
                  V <- vcov(em0)
                  Vinv <- try(chol2inv(chol(V)), silent = TRUE)
                  ## Compute unit errors.
                  if (!inherits(Vinv, "try-error")) {
                    ue <- 1 / diag(Vinv)
                  } else {
                    ue <- 1 / diag(solve(V))
                  }
                }))
    result[["ue"]] <- ue
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
      ## Estimatie heritability on a line mean basis.
      if (useRepId) {
        result[["heritability"]] <- varGen /
          (varGen + (varErr / dplyr::n_distinct(TD$repId)))
      } else {
        result[["heritability"]] <- varGen / (varGen + varErr)
      }
    }
  }
  ## Extract fitted values.
  if ("fitted" %in% what) {
    fitVal <- cbind(baseData, sapply(X = mf, FUN = fitted))
    result[["fitted"]] <- fitVal
  }
  ## Extract residuals.
  if ("resid" %in% what) {
    resVal <- cbind(baseData, sapply(X = mf, FUN = residuals))
    result[["resid"]] <- resVal
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
    result[["stdRes"]] <- stdRes
  }
  ## Compute rMeans.
  ## Use napredict to fill in NAs in data with NAs.
  if ("rMeans" %in% what) {
    rMeans <- cbind(baseData,
                    sapply(X = mr, FUN = function(mr0) {
                      napredict(attr(model.frame(mr0), "na.action"),
                                x = lme4::getME(mr0, "mu"))
                    }))
    result[["rMeans"]] <- rMeans
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
    result[["ranEf"]] <- ranEf
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
  if (length(result) == 1) {
    result <- result[[1]]
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
                          useRepId) {
  if (!requireNamespace("asreml", quietly = TRUE)) {
    stop("asreml cannot be successfully loaded.\n")
  }
  mf <- SSA$mFix
  mr <- SSA$mRand
  TD <- SSA$TD[[1]]
  predicted <- SSA$predicted
  whatTot <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "ue", "heritability", "varGen",
               "varErr", "fitted", "resid", "stdRes", "rMeans", "ranEf",
               "wald", "CV", "rDf", "sed", "lsd")
  whatMod <- c("F", "F", "R", "R", "F", "R", "R", "R", "F", "F", "F", "R", "R",
               "F", "F", "F", "F", "F")
  whatSSA <- c(if (!is.null(mf)) "F", if (!is.null(mr)) "R")
  whatPred <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "ue", "ranEf")
  if (what[[1]] == "all") {
    what <- whatTot[whatMod %in% whatSSA]
  } else {
    what <- match.arg(arg = what, choices = whatTot[whatMod %in% whatSSA],
                      several.ok = TRUE)
  }
  ## Create baseData and baseDataPred to which further results will be merged.
  base <- createBaseData(TD, predicted, keep, useRepId,
                         bdPred = any(what %in% whatPred))
  baseData <- base$baseData
  baseDataPred <- base$baseDataPred
  ## Create empty result list.
  result <- setNames(vector(mode = "list", length = length(what)),
                     what)
  ## Extract BLUEs and se of BLUEs from fixed model.
  if ("BLUEs" %in% what) {
    predVals <- lapply(X = traits, FUN = function(trait) {
      setNames(mf[[trait]]$predictions$pvals[c(predicted, "predicted.value")],
               c(predicted, trait))
    })
    BLUEs <- Reduce(f = merge, x = predVals, init = baseDataPred)
    result[["BLUEs"]] <- BLUEs
  }
  if ("seBLUEs" %in% what) {
    predErrs <- lapply(X = traits, FUN = function(trait) {
      setNames(mf[[trait]]$predictions$pvals[c(predicted, "standard.error")],
               c(predicted, trait))
    })
    seBLUEs <- Reduce(f = merge, x = predErrs, init = baseDataPred)
    result[["seBLUEs"]] <- seBLUEs
  }
  ## Extract BLUPs and se of BLUPs from fixed model.
  if ("BLUPs" %in% what) {
    predVals <- lapply(X = traits, FUN = function(trait) {
      setNames(mr[[trait]]$predictions$pvals[c(predicted, "predicted.value")],
               c(predicted, trait))
    })
    BLUPs <- Reduce(f = merge, x = predVals, init = baseDataPred)
    result[["BLUPs"]] <- BLUPs
  }
  if ("seBLUPs" %in% what) {
    predErrs <- lapply(X = traits, FUN = function(trait) {
      setNames(mr[[trait]]$predictions$pvals[c(predicted, "standard.error")],
               c(predicted, trait))
    })
    seBLUPs <- Reduce(f = merge, x = predErrs, init = baseDataPred)
    result[["seBLUPs"]] <- seBLUPs
  }
  ## Compute unit errors.
  if ("ue" %in% what) {
    ue <- cbind(baseDataPred, sapply(X = mf, FUN = function(mf0) {
      ## Extract V from mf.
      V <- mf0$predictions$vcov
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
    result[["ue"]] <- ue
  }
  ## Extract variances
  varGen <- sapply(X = mr, FUN = function(mr0) {
    unname(mr0$gammas[grep(pattern = paste0(predicted, "!", predicted, ".var"),
                           x = names(mr0$gammas))] * mr0$sigma2)
  })
  if ("varGen" %in% what) {
    result[["varGen"]] <- varGen
  }
  if ("varErr" %in% what) {
    result[["varErr"]] <- sapply(X = mr, FUN = function(mr0) {
      unname(mr0$gammas[grep(pattern = "R!variance",
                             x = names(mr0$gammas))] * mr0$sigma2)
    })
  }
  ## Estimate heritability on a line mean basis.
  if ("heritability" %in% what) {
    result[["heritability"]] <- sapply(X = traits, FUN = function(trait) {
      sedSq <- predictAsreml(model = mr[[trait]], classify = predicted,
                             vcov = FALSE, TD = TD, only = predicted,
                             sed = TRUE)$predictions$sed ^ 2
      unname(1 - mean(sedSq[lower.tri(sedSq)]) / (2 * varGen[[trait]]))
    })
  }
  ## Extract fitted values.
  if ("fitted" %in% what) {
    fitVal <- cbind(baseData, sapply(X = mf, FUN = fitted))
    result[["fitted"]] <- fitVal
  }
  ## Extract residuals.
  if ("resid" %in% what) {
    resVal <- cbind(baseData, sapply(X = mf, FUN = residuals,
                                     type = "response"))
    result[["resid"]] <- resVal
  }
  ## Extract standardized residuals.
  if ("stdRes" %in% what) {
    stdRes <- cbind(baseData, sapply(X = mf, FUN = residuals, type = "stdCond"))
    result[["stdRes"]] <- stdRes
  }
  ## Extract rMeans.
  if ("rMeans" %in% what) {
    rMeans <- cbind(baseData, sapply(X = mr, FUN = fitted))
    result[["rMeans"]] <- rMeans
  }
  ## Extract random effects.
  if ("ranEf" %in% what) {
    ranEffs <- lapply(X = traits, FUN = function(trait) {
      coefs <- mr[[trait]]$coe$random
      ranEff <- data.frame(gsub(pattern = paste0(predicted, "_"),
                                replacement = "",
                                x = names(coefs)[grep(pattern = predicted,
                                                      x = names(coefs))]),
                           coefs[grep(pattern = predicted, x = names(coefs))])
      colnames(ranEff) <- c(predicted, trait)
      return(ranEff)
    })
    ranEf <- Reduce(f = merge, x = ranEffs, init = baseDataPred)
    result[["ranEf"]] <- ranEf
  }
  ## Compute wald test.
  if ("wald" %in% what) {
    tmpfile <- tempfile()
    sink(file = tmpfile)
    result[["wald"]] <- lapply(X = mf, function(mf0) {
      wtt <- asreml::wald.asreml(mf0, ssType = "conditional", denDF = "numeric")
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
  ## Extract standard error of difference.
  if ("sed" %in% what) {
    result[["sed"]] <- lapply(X = mf, FUN = function(mf0) {
      mf0$predictions$avsed
    })
  }
  ## Compute lsd with signifcance level 5%.
  if ("lsd" %in% what) {
    result[["lsd"]] <- lapply(X = mf, FUN = function(mf0) {
      qt(p = .975, df = mf0$nedf) * mf0$predictions$avsed
    })
  }
  if (length(result) == 1) {
    result <- result[[1]]
  }
  return(result)
}

#' Helper function for creating baseData
#'
#' @keywords internal
createBaseData <- function(TD, predicted, keep, useRepId, bdPred) {
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




