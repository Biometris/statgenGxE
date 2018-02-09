#' Extracts statistics from Fitted Models
#'
#' This function extracts and calculates various results for fitted models such as
#' BLUEs, BLUPs, unit errors, heritabilities. Note that most results can only
#' calculated if a model is fitted with genotype is fixed or random. This is
#' indicated in the list below with "F" and "R"
#'
#' Possible options for \code{what} are:
#' \describe{
#' \item{F - BLUEs}{Best Lineair Unbiased Estimatiors}
#' \item{F - seBLUES}{Standard errors of the BLUEs}
#' \item{R - BLUPs}{Best Lineair Unbiased Predictors}
#' \item{R - seBLUPs}{Standard errors of the BLUPs}
#' \item{F - ue}{unit errors - only for \code{lme4} and \code{asreml}}
#' \item{R - heritability}{heritability}
#' \item{R - varGen}{genetic variance component}
#' \item{R - varErr}{residual variance component - only for \code{lme4}
#' and \code{asreml}}
#' \item{R - varSpat}{spatial variance components - only for \code{SpATS}}
#' \item{F - fitted}{fitted values for the model with genotype as fixed component}
#' \item{F - resid}{residuals for the model with genotype as fixed component}
#' \item{F - stdRes}{standardized residuals for the model with genotype as fixed
#' component - only for \code{lme4} and \code{asreml}}
#' \item{R - rMeans}{fitted values for the model with genotype as random component}
#' \item{R - ranEf}{random genetic effects}
#' \item{F - wald}{results of the wald test}
#' \item{F - CV}{coefficient of variation - only for \code{lme4} and \code{asreml}}
#' \item{F - rDf}{residual degrees of freedom}
#' \item{R - effDim}{effective dimensions - only for \code{SpATS}}
#' \item{F - sed}{standard error of difference - only for \code{asreml}}
#' \item{F - lsd}{least significant difference - only for \code{asreml}}
#' }
#'
#' @param SSA An object of class SSA.
#' @param traits A character vector of traits for which the statistics should be
#' computed If not supplied statistics are computed for all traits that have
#' been modelled.
#' @param what A character vector indicating which statistics should be computed.
#' Most statistics are available for all models, some only for models fitted using
#' a certain engine. If this is the case this is indicated in the list with options
#' in details.\cr
#' If \code{what = "all"} all available statistics are computed.\cr
#' @param keep A character vector of column(s) in the object of class
#' \code{\link{TD}} used for modeling. These columns will be kept as output when
#' computing fitted values, residuals, standardized residuals and rMeans. Columns
#' can also be kept when computing (se)BLUEs and (se)BLUPs but only if the column
#' to keep contains unique values for the modeled variables, i.e. a column repId
#' with several values per genotype cannot be kept.
#'
#' @return A list of extracted statistics.
#'
#' @seealso
#' \code{\link{STRunModel}}, \code{\link{STModSpATS}}, \code{\link{STModLme4}} and
#' \code{\link{STModAsreml}}
#'
#' @examples
#' ## Fit model using lme4.
#' myModel <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield")
#'
#' ## Extract statistics from fitted model.
#' extr <- STExtract(myModel)
#'
#' @export
STExtract <- function(SSA,
                      traits = SSA$traits,
                      what = "all",
                      keep = NULL) {
  ## Checks.
  if (!is.SSA(SSA)) {
    stop("SSA has to be an object of class SSA.\n")
  }
  if (is.null(traits) || !is.character(traits) || !all(traits %in% colnames(SSA$TD))) {
    stop("All traits have to be columns in TD.\n")
  }
  if (!is.null(keep) && (!is.character(keep) || !all(keep %in% colnames(SSA$TD)))) {
    stop("All items in keep have to be columns in TD.\n")
  }
  engine <- SSA$engine
  ## Set useRepId to TRUE when it is used as fixed effect in the model.
  useRepId <- (SSA$design %in% c("res.ibd", "res.rowcol", "rcbd"))
  ## Extract statistics from fitted model.
  result <- do.call(what = paste0("extract", tools::toTitleCase(engine)),
                    args = list(SSA = SSA, traits = traits, what = what,
                                useRepId = useRepId, keep = keep))
  attr(x = result, which = "traits") <- traits
  attr(x = result, which = "design") <- SSA$design
  attr(x = result, which = "engine") <- engine
  return(result)
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
  TD <- SSA$TD
  predicted <- attr(SSA, "predicted")
  useCheckId <- length(grep(pattern = "checkId",
                            x = deparse(mr[[1]]$model$fixed))) > 0
  whatTot <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "heritability", "varGen",
               "varSpat", "fitted", "resid", "rMeans", "ranEf", "rDf", "effDim")
  whatMod <- c("F", "F", "R", "R", "R", "R", "R", "F", "F", "R", "R", "F", "F")
  whatSSA <- c(if (!is.null(mf)) "F", if (!is.null(mr)) "R")
  if (what[[1]] == "all") {
    what <- whatTot[whatMod %in% whatSSA]
  } else {
    what <- match.arg(arg = what, choices = whatTot[whatMod %in% whatSSA],
                      several.ok = TRUE)
  }
  ## Create baseData and baseDataPred to which further results will be merged.
  base <- createBaseData(TD, predicted, keep)
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
    result[["BLUEs"]] <- createTD(data = Reduce(f = merge,
                                                x = c(list(baseDataPred),
                                                      predVals)))
  }
  if ("seBLUEs" %in% what) {
    predErrs <- lapply(X = traits, FUN = function(trait) {
      predErr <- predict(mf[[trait]],
                         which = predicted)[c(predicted, "standard.errors")]
      colnames(predErr) <- c(predicted, trait)
      return(predErr)
    })
    result[["seBLUEs"]] <- createTD(data = Reduce(f = merge,
                                                  x = c(list(baseDataPred),
                                                        predErrs)))
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
    result[["BLUPs"]] <- createTD(data = Reduce(f = merge,
                                                x = c(list(baseDataPred),
                                                      predVals)))
  }
  if ("seBLUPs" %in% what) {
    whichPred <- c(predicted, if (useCheckId) "checkId")
    predErrs <- lapply(X = traits, FUN = function(trait) {
      predErr <- predict(mr[[trait]],
                         which = whichPred)[c(whichPred, "standard.errors")]
      colnames(predErr) <- c(whichPred, trait)
      return(predErr)
    })
    result[["seBLUPs"]] <- createTD(data = Reduce(f = merge,
                                                  x = c(list(baseDataPred),
                                                        predErrs)))
  }
  ## Compute generalized heritability.
  if ("heritability" %in% what) {
    result[["heritability"]] <- sapply(X = mr, FUN = function(mr0) {
      unname(mr0$eff.dim[predicted] / (mr0$dim[predicted] - 1))
    })
  }
  ## Extract genetic variance.
  if ("varGen" %in% what) {
    result[["varGen"]] <- sapply(X = mr, FUN = function(mr0) {
      unname(mr0$var.comp[predicted])
    })
  }
  ## Extract spartial variance.
  if ("varSpat" %in% what) {
    result[["varSpat"]] <- sapply(X = mr, FUN = function(mr0) {
      mr0$var.comp[grep(pattern = "Coordinates", x = names(mr0$var.comp))]
    })
  }
  ## Extract fitted values.
  if ("fitted" %in% what) {
    result[["fitted"]] <- createTD(data = cbind(baseData,
                                                sapply(X = mf, FUN = fitted)))
  }
  ## Extract residuals.
  if ("resid" %in% what) {
    result[["resid"]] <- createTD(data = cbind(baseData,
                                               sapply(X = mf, FUN = residuals)))
  }
  ## Extract rMeans.
  if ("rMeans" %in% what) {
    result[["rMeans"]] <- createTD(data = cbind(baseData,
                                                sapply(X = mr, FUN = fitted)))
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
    result[["ranEf"]] <- createTD(data = Reduce(f = merge,
                                                x = c(list(baseDataPred),
                                                      ranEffs)))
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
  TD <- SSA$TD
  predicted = attr(SSA, "predicted")
  whatTot <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "ue", "heritability",
               "varGen", "varErr", "fitted", "resid", "stdRes", "rMeans", "ranEf",
               "wald", "CV", "rDf")
  whatMod <- c("F", "F", "R", "R", "F", "R", "R", "R", "F", "F", "F", "R", "R",
               "F", "F", "F")
  whatSSA <- c(if (!is.null(mf)) "F", if (!is.null(mr)) "R")
  if (what[[1]] == "all") {
    what <- whatTot[whatMod %in% whatSSA]
  } else {
    what <- match.arg(arg = what, choices = whatTot[whatMod %in% whatSSA],
                      several.ok = TRUE)
  }
  ## Create baseData and baseDataPred to which further results will be merged.
  base <- createBaseData(TD, predicted, keep)
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
    result[["BLUEs"]] <- createTD(data = Reduce(f = merge,
                                                x = c(list(baseDataPred),
                                                      predVals)))
  }
  if ("seBLUEs" %in% what) {
    predErrs <- lapply(X = traits, FUN = function(trait) {
      setNames(emStats[[trait]][c(predicted, "SE")],
               c(predicted, trait))
    })
    result[["seBLUEs"]] <- createTD(data = Reduce(f = merge,
                                                  x = c(list(baseDataPred),
                                                        predErrs)))
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
    result[["BLUPs"]] <- createTD(data = Reduce(f = merge,
                                                x = c(list(baseDataPred),
                                                      predVals)))
  }
  if ("seBLUPs" %in% what) {
    predErrs <- lapply(X = traits, FUN = function(trait) {
      ranEffs <- lme4::ranef(mr[[trait]], condVar = TRUE)[[predicted]]
      predErr <- data.frame(rownames(ranEffs),
                            as.vector(sqrt(attr(ranEffs, "postVar"))))
      colnames(predErr) <- c(predicted, trait)
      return(predErr)
    })
    result[["seBLUPs"]] <- createTD(data = Reduce(f = merge,
                                                  x = c(list(baseDataPred),
                                                        predErrs)))
  }
  ## Compute unit errors.
  if ("ue" %in% what) {
    result[["ue"]] <- cbind(baseDataPred,
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
    result[["fitted"]] <- createTD(data = cbind(baseData,
                                                sapply(X = mf, FUN = fitted)))
  }
  ## Extract residuals.
  if ("resid" %in% what) {
    result[["resid"]] <- createTD(data = cbind(baseData,
                                               sapply(X = mf, FUN = residuals)))
  }
  ## Extract standardized residuals.
  if ("stdRes" %in% what) {
    result[["stdRes"]] <- createTD(data = cbind(baseData,
                                                sapply(X = mf, FUN = function(mf0) {
                                                  if (inherits(mf0, "lm")) {
                                                    stdRes <- rstandard(mf0)
                                                  } else if (inherits(mf0, "lmerMod")) {
                                                    stdRes <- residuals(mf0, scaled = TRUE)
                                                  }
                                                })))
  }
  ## Compute rMeans.
  ## Use napredict to fill in NAs in data with NAs.
  if ("rMeans" %in% what) {
    result[["rMeans"]] <- createTD(data = cbind(baseData,
                                                sapply(X = mr, FUN = function(mr0) {
                                                  napredict(attr(model.frame(mr0),
                                                                 "na.action"),
                                                            x = lme4::getME(mr0, "mu"))
                                                })))
  }
  ## Extract random effects.
  if ("ranEf" %in% what) {
    ranEffs <- lapply(X = traits, FUN = function(trait) {
      effs <- lme4::ranef(mr[[trait]], drop = TRUE)[[predicted]]
      ranEff <- data.frame(names(effs), effs)
      colnames(ranEff) <- c(predicted, trait)
      return(ranEff)
    })
    result[["ranEf"]] <- createTD(data = Reduce(f = merge,
                                                x = c(list(baseDataPred),
                                                      ranEffs)))
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
  TD <- SSA$TD
  predicted <- attr(SSA, "predicted")
  whatTot <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "ue", "heritability", "varGen",
               "varErr", "fitted", "resid", "stdRes", "rMeans", "ranEf",
               "wald", "CV", "rDf", "sed", "lsd")
  whatMod <- c("F", "F", "R", "R", "F", "R", "R", "R", "F", "F", "F", "R", "R",
               "F", "F", "F", "F", "F")
  whatSSA <- c(if (!is.null(mf)) "F", if (!is.null(mr)) "R")
  if (what[[1]] == "all") {
    what <- whatTot[whatMod %in% whatSSA]
  } else {
    what <- match.arg(arg = what, choices = whatTot[whatMod %in% whatSSA],
                      several.ok = TRUE)
  }
  ## Create baseData and baseDataPred to which further results will be merged.
  base <- createBaseData(TD, predicted, keep)
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
    result[["BLUEs"]] <- createTD(data = Reduce(f = merge,
                                                x = c(list(baseDataPred),
                                                      predVals)))
  }
  if ("seBLUEs" %in% what) {
    predErrs <- lapply(X = traits, FUN = function(trait) {
      setNames(mf[[trait]]$predictions$pvals[c(predicted, "standard.error")],
               c(predicted, trait))
    })
    result[["seBLUEs"]] <- createTD(data = Reduce(f = merge,
                                                  x = c(list(baseDataPred),
                                                        predErrs)))
  }
  ## Extract BLUPs and se of BLUPs from fixed model.
  if ("BLUPs" %in% what) {
    predVals <- lapply(X = traits, FUN = function(trait) {
      setNames(mr[[trait]]$predictions$pvals[c(predicted, "predicted.value")],
               c(predicted, trait))
    })
    result[["BLUPs"]] <- createTD(data = Reduce(f = merge,
                                                x = c(list(baseDataPred),
                                                      predVals)))
  }
  if ("seBLUPs" %in% what) {
    predErrs <- lapply(X = traits, FUN = function(trait) {
      setNames(mr[[trait]]$predictions$pvals[c(predicted, "standard.error")],
               c(predicted, trait))
    })
    result[["seBLUPs"]] <- createTD(data = Reduce(f = merge,
                                                  x = c(list(baseDataPred),
                                                        predErrs)))
  }
  ## Compute unit errors.
  if ("ue" %in% what) {
    result[["ue"]] <- createTD(data = cbind(baseDataPred,
                                            sapply(X = mf, FUN = function(mf0) {
                                              ## Extract V from mf.
                                              V <- mf0$predictions$vcov
                                              ## Remove columns and rows containing NA.
                                              VMiss <- apply(X = V, MARGIN = 2,
                                                             FUN = anyNA)
                                              V <- V[!VMiss, !VMiss]
                                              Vinv <- try(chol2inv(chol(V)),
                                                          silent = TRUE)
                                              ## Compute unit errors.
                                              ue <- rep(x = NA,
                                                        times = nrow(baseDataPred))
                                              if (!inherits(Vinv, "try-error")) {
                                                ue[!VMiss] <- 1 / diag(Vinv)
                                              } else {
                                                ue[!VMiss] <- 1 / diag(solve(V))
                                              }
                                              return(ue)
                                            })))
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
    result[["fitted"]] <- createTD(data = cbind(baseData,
                                                sapply(X = mf, FUN = fitted)))
  }
  ## Extract residuals.
  if ("resid" %in% what) {
    result[["resid"]] <- createTD(data = cbind(baseData,
                                               sapply(X = mf, FUN = residuals,
                                                      type = "response")))
  }
  ## Extract standardized residuals.
  if ("stdRes" %in% what) {
    result[["stdRes"]] <- createTD(data = cbind(baseData,
                                                sapply(X = mf, FUN = residuals,
                                                       type = "stdCond")))
  }
  ## Extract rMeans.
  if ("rMeans" %in% what) {
    result[["rMeans"]] <- createTD(data = cbind(baseData,
                                                sapply(X = mr, FUN = fitted)))
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
      ranEff
    })
    result[["ranEf"]] <- createTD(data = Reduce(f = merge,
                                                x = c(list(baseDataPred),
                                                      ranEffs)))
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

#' Extract statistics from model fitted using sommer
#'
#' @keywords internal
extractSommer <- function(SSA,
                          traits = SSA$traits,
                          what = "all",
                          keep = NULL,
                          useRepId) {
  mf <- SSA$mFix
  mr <- SSA$mRand
  TD <- SSA$TD
  predicted = attr(SSA, "predicted")
  whatTot <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "ue", "heritability",
               "varGen", "varErr", "fitted", "resid", "stdRes", "rMeans",
               "ranEf", "rDf", "CV")
  whatMod <- c("F", "F", "R", "R", "F", "R", "R", "R", "F", "F", "F", "R", "R",
               "F", "F")
  whatSSA <- c(if (!is.null(mf)) "F", if (!is.null(mr)) "R")
  if (what[[1]] == "all") {
    what <- whatTot[whatMod %in% whatSSA]
  } else {
    what <- match.arg(arg = what, choices = whatTot[whatMod %in% whatSSA],
                      several.ok = TRUE)
  }
  ## Create baseData and baseDataPred to which further results will be merged.
  base <- createBaseData(TD, predicted, keep)
  baseData <- base$baseData
  baseDataPred <- base$baseDataPred
  ## Create empty result list.
  result <- setNames(vector(mode = "list", length = length(what)), what)
  # ## Compute BLUEs.
  if ("BLUEs" %in% what) {
    result[["BLUEs"]] <- cbind(baseDataPred,
                               sapply(X = mf, FUN = function(mf0) {
                                 betaTot <- coef(mf0)
                                 beta0 <- betaTot[-grep(pattern = predicted,
                                                        x = rownames(betaTot)), ]
                                 ic <- beta0[1] + (sum(beta0[2:length(beta0)]) /
                                                     length(beta0))
                                 c(ic, betaTot[grep(pattern = predicted,
                                                    x = rownames(betaTot)), ] +
                                     ic)
                               }))
  }
  # if ("seBLUEs" %in% what) {
  #   result[["seBLUEs"]] <- cbind(baseDataPred,
  #                                sapply(X = emStats, FUN = "[[", "SE"))
  # }
  ## Compute BLUPs and se of BLUPs from mixed model.
  if ("BLUPs" %in% what) {
    predVals <- lapply(X = traits, FUN = function(trait) {
      beta <- coef(mr[[trait]])
      ic <- beta[1, ] + sum(beta[2:nrow(beta), ]) / nrow(beta)
      vals <- mr[[trait]]$u.hat[[predicted]]
      predVal <- data.frame(gsub(pattern = predicted, replacement = "",
                                 x = rownames(vals)),
                            vals + ic)
      colnames(predVal) <- c(predicted, trait)
      return(predVal)
    })
    result[["BLUPs"]] <- createTD(data = Reduce(f = merge,
                                                x = c(list(baseDataPred),
                                                      predVals)))
  }
  if ("seBLUPs" %in% what) {
    result[["seBLUPs"]] <- cbind(baseDataPred,
                                 sapply(X = mr, FUN = function(mr0) {
                                   ## Compute se of Blups.
                                   Matrix::diag(mr0$PEV.u.hat[[predicted]][[traits]]) /
                                     sqrt(nrow(baseDataPred))
                                 }))
  }
  ## Compute unit errors.
  if ("ue" %in% what) {
    result[["ue"]] <- cbind(baseDataPred,
                            lapply(X = traits, FUN = function(trait) {
                              ## Extract inverted variance covariance matrix.
                              ues <- data.frame(1 / Matrix::diag(mf[[trait]]$V.inv) *
                                                  (nrow(baseDataPred) /
                                                     nrow(baseData)))
                              colnames(ues) <- trait
                              return(ues)
                            }))
  }
  ## Extract variances.
  varCor <- lapply(X = mr, FUN = "[[", "var.comp")
  if ("varGen" %in% what) {
    result[["varGen"]] <- sapply(X = varCor, FUN = "[[", predicted)
  }
  if ("varErr" %in% what) {
    result[["varErr"]] <- sapply(X = varCor, FUN = "[[", "units")
  }
  if ("heritability" %in% what) {
    ## Estimate heritability on a line mean basis.
    result[["heritability"]] <- lapply(X = mr, FUN = function(mr0) {
      sommer::h2.fun(mr0, TD, gTerm = predicted)[["H2.stdrd"]]
    })
  }
  ## Extract fitted values.
  if ("fitted" %in% what) {
    result[["fitted"]] <- cbind(baseData,
                                sapply(X = mf, FUN = fitted))
  }
  ## Extract residuals.
  if ("resid" %in% what) {
    result[["resid"]] <- cbind(baseData,
                               sapply(X = mf, FUN = residuals))
  }
  ## Extract standardized residuals.
  if ("stdRes" %in% what) {
    result[["stdRes"]] <- cbind(baseData,
                                sapply(X = mf, FUN = function(mf0) {
                                  scale(residuals(mf0))}))
  }
  ## Compute rMeans.
  if ("rMeans" %in% what) {
    result[["rMeans"]] <- cbind(baseData,
                                sapply(X = mr, FUN = fitted))
  }
  ## Extract random effects.
  if ("ranEf" %in% what) {
    ## Use sink to suppress randef default text message.
    tmpfile <- tempfile()
    sink(file = tempfile())
    ranEffs <- lapply(X = traits, FUN = function(trait) {
      effs <- sommer::randef(mr[[trait]])[[predicted]]
      ranEff <- data.frame(gsub(pattern = predicted, replacement = "",
                                 x = rownames(effs)),
                           effs)
      colnames(ranEff) <- c(predicted, trait)
      return(ranEff)
    })
    result[["ranEf"]] <- createTD(data = Reduce(f = merge,
                                                x = c(list(baseDataPred),
                                                      ranEffs)))
    sink()
    unlink(tmpfile)
  }
  if ("CV" %in% what) {
    ## Compute Coefficient of Variation.
    result[["CV"]] <- sapply(X = mf, FUN = function(mf0) {
      unname(100 * sqrt(mf0$sigma["units"]) /
               mean(fitted(mf0), na.rm = TRUE))
    })
  }
  ## Compute residual degrees of freedom.
  if ("rDf" %in% what) {
    result[["rDf"]] <- sapply(X = mf, FUN = function(mf0) {
      length(mf0$Y) - Matrix::rankMatrix(mf0$X)
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
createBaseData <- function(TD, predicted, keep) {
  ## Create baseData consisting of predicted variable, possibly repId and
  ## selected keep columns.
  baseData <- TD[, colnames(TD) %in% c(predicted, "repId", keep), drop = FALSE]
  ## Create baseData for predictions with predicted variable(s).
  baseDataPred <- unique(TD[predicted])
  ## Add columns in keep one-by-one. Only data that is constant within
  ## predicted is actually kept. Other columns are dropped with a warning.
  for (col in keep) {
    TDKeep <- unique(TD[, c(predicted, col)])
    if (!anyDuplicated(TDKeep[, predicted])) {
      baseDataPred <- merge(baseDataPred, TDKeep, by = predicted)
    } else {
      warning(paste0("Duplicate values for ", deparse(substitute(col)),
                     "Column dropped.\n"), call. = FALSE)
    }
  }
  return(list(baseData = baseData, baseDataPred = baseDataPred))
}




