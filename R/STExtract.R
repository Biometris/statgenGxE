#' Extracts statistics from Fitted Models
#'
#' This function extracts and calculates various results for fitted models such as
#' BLUEs, BLUPs, unit errors, heritabilities.
#'
#' @param SSA An object of class SSA.
#'
#' @return A list of extracted statistics.
#'
#' @seealso
#' \code{\link{STRunModel}}, \code{\link{STModSpATS}}, \code{\link{STModLme4}} and
#' \code{\link{STModAsreml}}
#'
#' @examples
#' ## Load data.
#' data(TDHeat05)
#'
#' ## Fit model using lme4.
#' myModel <- STRunModel(TD = TDHeat05, design = "res.rowcol", trait = "yield")
#'
#' ## Extract statistics from fitted model.
#' extr <- STExtract(myModel)
#'
#' @export

STExtract <- function(SSA,
                      traits = SSA$traits,
                      what = "all") {
  engine <- SSA$engine
  ## Set useRepId to TRUE when it is used as fixed effect in the model.
  useRepId <- (SSA$design %in% c("res.ibd", "res.rowcol", "rcbd"))
  ## Extract statistics from fitted model.
  result <- do.call(what = paste0("extract", tools::toTitleCase(engine)),
                    args = list(SSA = SSA, traits = traits, useRepId = useRepId))
  if (useRepId) {
    TD <- SSA$data
    repVec <- tapply(X = TD$repId, INDEX = TD$genotype, FUN = nlevels)
    result$minReps <- min(repVec, na.rm = TRUE)
    result$meanReps <- mean(repVec, na.rm = TRUE)
    result$maxReps <- max(repVec, na.rm = TRUE)
  } else {
    result$minReps <- 1
    result$meanReps <- 1
    result$maxReps <- 1
  }
  result$traits <- SSA$traits
  result$model <- SSA$design
  result$engine <- engine
  return(result)
}


#' Extract statistics from model fitted using SpATS
#'
#' @keywords internal
extractSpATS <- function(SSA,
                         traits = SSA$traits,
                         what = "all",
                         useRepId) {
  if (what == "all") {
    what <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "heritability", "varGen",
              "varSpat", "fitted", "resid", "rMeans", "ranEf", "rDf", "effDim")
  }
  mf <- SSA$mFix
  mr <- SSA$mMix
  TD <- SSA$data
  ## Extract names of genotypes.
  genoNames <- mf[[1]]$terms$geno$geno_names
  ## Create baseData consisting of genotype and possibly repId
  baseData <- TD[, colnames(TD) %in% c("genotype", "repId")]
  ## Create empty result list.
  result <- setNames(vector(mode = "list", length = length(what)),
                     what)
  ## Compute BLUEs and se of BLUEs from fixed model.
  if ("BLUEs" %in% what) {
    result[["BLUEs"]] <- cbind2(data.frame(genotype = genoNames),
                                sapply(X = mf, FUN = function(mf0) {
                                  SpATS::predict.SpATS(mf0,
                                                       which = "genotype")$predicted.values
                                }))
  }
  if ("seBLUEs" %in% what) {
    result[["seBLUEs"]] <- cbind2(data.frame(genotype = genoNames),
                                  sapply(X = mf, FUN = function(mf0) {
                                    SpATS::predict.SpATS(mf0,
                                                         which = "genotype")$standard.errors
                                  }))
  }
  ## Compute BLUPs and se of BLUPs from mixed model.
  if ("BLUPs" %in% what) {
    result[["BLUPs"]] <- cbind2(data.frame(genotype = genoNames),
                                sapply(X = mr, FUN = function(mr0) {
                                  SpATS::predict.SpATS(mr0,
                                                       which = "genotype")$predicted.values
                                }))
  }
  if ("seBLUPs" %in% what) {
    result[["seBLUPs"]] <- cbind2(data.frame(genotype = genoNames),
                                  sapply(X = mr, FUN = function(mr0) {
                                    SpATS::predict.SpATS(mr0,
                                                         which = "genotype")$standard.errors
                                  }))
  }
  ## Compute generalized heritability.
  if ("heritability" %in% what) {
    result[["heritability"]] <- sapply(X = mr, FUN = function(mr0) {
      unname(mr0$eff.dim["genotype"] / (mr0$dim["genotype"]))
    })
  }
  ## Extract genetic variance.
  if ("varGen" %in% what) {
    result[["varGen"]] <- sapply(X = mr, FUN = function(mr0) {
      unname(mr0$var.comp["genotype"])
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
    result[["fitted"]] <- cbind2(baseData,
                                 sapply(X = mf, FUN = fitted))
  }
  ## Extract residuals.
  if ("resid" %in% what) {
    result[["resid"]] <- cbind2(baseData,
                                sapply(X = mf, FUN = residuals))
  }
  ## Extract rMeans.
  if ("rMeans" %in% what) {
    result[["rMeans"]] <- cbind2(baseData,
                                 sapply(X = mr, FUN = fitted))
  }
  ## Extract random effects.
  if ("ranEf" %in% what) {
    result[["ranEf"]] <- cbind2(data.frame(genotype = genoNames),
                                sapply(X = mr, FUN = function(mr0) {
                                  mr0$coeff[names(mr0$coeff) %in% genoNames]
                                }))
  }
  ## Extract residual degrees of freedom.
  if ("rDf" %in% what) {
    result[["rDf"]] <- sapply(X = mf, FUN = function(mf0) {
      unname(mf0$dim["genotype"])
    })
  }
  ## Extract effective dimensions.
  if ("effDim" %in% what) {
    result[["effDim"]] <- sapply(X = mr, FUN = function(mr0) {
      mr0$eff.dim
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
                        useRepId) {
  if (what == "all") {
    what <- c("BLUEs", "seBLUEs", "BLUPs", "seBLUPs", "ue", "heritability", "varGen",
              "varErr", "fitted", "resid", "stdRes", "rMeans", "ranEf",
              "wald", "CV", "rDf")
  }
  mf <- SSA$mFix
  mr <- SSA$mMix
  TD <- SSA$data
  ## Create baseData consisting of genotype and possibly repId
  baseData <- TD[, colnames(TD) %in% c("genotype", "repId")]
  ## Create empty result list.
  result <- setNames(vector(mode = "list", length = length(what)),
                     what)
  ## Use emmeans to convert mf to emmGrid per trait.
  ## Option lmer.df = "asymptotic" for speed. Default options gives
  ## exact the same results as asreml but is very slow.
  em <- sapply(X = mf, FUN = function(mf0) {
    emmeans::emmeans(mf0, specs = "genotype", lmer.df = "asymptotic")
  }, simplify = FALSE)
  ## Summarize emGrid to calculate BLUEs en se per trait.
  emStats <- sapply(X = em, FUN = summary, simplify = FALSE)
  ## Extract genonames from emSTats.
  genoNames <- sapply(X = emStats, FUN = "[[", "genotype", simplify = FALSE)
  ## Extract BLUEs and se of BLUEs from emStats.
  if ("BLUEs" %in% what) {
    result[["BLUEs"]] <- cbind2(data.frame(genotype = genoNames[[1]]),
                                sapply(X = emStats, FUN = "[[", "emmean"))
  }
  if ("seBLUEs" %in% what) {
    result[["seBLUEs"]] <- cbind2(data.frame(genotype = genoNames[[1]]),
                                sapply(X = emStats, FUN = "[[", "SE"))
  }
  ## Compute BLUPs and se of BLUPs from mixed model.
  if ("BLUPs" %in% what) {
    result[["BLUPs"]] <- cbind2(data.frame(genotype = genoNames[[1]]),
                                sapply(X = mr, FUN = function(mr0) {
                                  ## Extract fixed effects.
                                  fixEfMr <- lme4::fixef(mr0)
                                  ## Get positions of 'repId' within fixed effects.
                                  repIdPos <- grep(pattern = "repId",
                                                   x = names(fixEfMr))
                                  repIdMix <- fixEfMr[repIdPos]
                                  ## Compute Blups.
                                  coef(mr0)$genotype[, "(Intercept)"] +
                                    mean(c(repIdMix, 0))
                                }))
  }
  if ("seBLUPs" %in% what) {
    result[["seBLUPs"]] <- cbind2(data.frame(genotype = genoNames[[1]]),
                                  sapply(X = mr, FUN = function(mr0) {
                                    ## Compute se of Blups.
                                    as.vector(sqrt(attr(lme4::ranef(mr0,
                                                                    condVar = TRUE)[["genotype"]],
                                                        "postVar")))
                                  }))
  }
  ## Compute unit errors.
  if ("ue" %in% what) {
    result[["ue"]] <- cbind2(data.frame(genotype = genoNames[[1]]),
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
  varCor <- lapply(X = mr, FUN = lme4::VarCorr)
  varGen <- sapply(X = varCor, FUN = function(vc) {
    vc[["genotype"]][1, 1]
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
        (varGen + (varErr / length(unique(TD$repId))))
    } else {
      result[["heritability"]] <- varGen / (varGen + varErr)
    }
  }
  ## Extract fitted values.
  if ("fitted" %in% what) {
    result[["fitted"]] <- cbind2(baseData,
                                 sapply(X = mf, FUN = fitted))
  }
  ## Extract residuals.
  if ("resid" %in% what) {
    result[["resid"]] <- cbind2(baseData,
                                sapply(X = mf, FUN = residuals))
  }
  ## Extract standardized residuals.
  if ("stdRes" %in% what) {
    result[["stdRes"]] <- cbind2(baseData,
                                 sapply(X = mf, FUN = function(mf0) {
                                   if (class(mf0) == "lm") {
                                     stdRes <- rstandard(mf0)
                                   } else if (class(mf0) == "lmerMod") {
                                     stdRes <- residuals(mf0, scaled = TRUE)
                                   }
                                 }))
  }
  ## Compute rMeans.
  ## Use napredict to fill in NAs in data with NAs.
  if ("rMeans" %in% what) {
    result[["rMeans"]] <- cbind2(baseData,
                                 sapply(X = mr, FUN = function(mr0) {
                                   napredict(attr(model.frame(mr0),
                                                  "na.action"),
                                             x = lme4::getME(mr0, "mu"))
                                 }))
  }
  ## Extract random effects.
  if ("ranEf" %in% what) {
    result[["ranEf"]] <- cbind2(data.frame(genotype = genoNames),
                                sapply(X = mr, FUN = function(mr0) {
                                  lme4::ranef(mr0, drop = TRUE)[["genotype"]]
                                }))
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
  return(result)
}

#' Extract statistics from model fitted using asreml
#'
#' @keywords internal
extractAsreml <- function(SSA,
                          traits = SSA$traits,
                          what = c("BLUEs"),
                          useRepId) {
  mf <- SSA$mFix[[traits[1]]]
  mr <- SSA$mMix[[traits[1]]]
  TD <- SSA$data
  genoNames <- mf$predictions$pvals$genotype
  if (!requireNamespace("asreml", quietly = TRUE)) {
    stop("asreml cannot be successfully loaded.\n")
  }
  ## Extract BLUEs and SE of BLUEs from fixed model.
  predBlues <- setNames(mf$predictions$pvals$predicted.value,
                        genoNames)
  seBlues <- mf$predictions$pvals$standard.error
  ## Extract and invert variance covariance matrix.
  genoMiss <- is.na(predBlues)
  V <- mf$predictions$vcov[!genoMiss, !genoMiss]
  Vinv <- try(chol2inv(chol(V)), silent = TRUE)
  ## Compute unit errors.
  ue <- rep(x = NA, times = length(genoNames))
  if (!inherits(Vinv, "try-error")) {
    ue[!genoMiss] <- 1 / diag(Vinv)
  } else {
    ue[!genoMiss] <- 1 / diag(solve(V))
  }
  ## Extract Blups and se of Blups.
  predBlups <- mr$predictions$pvals$predicted.value
  seBlups <- mr$predictions$pvals$standard.error
  ## Collect predictions in data.frame.
  stats <- data.frame(genotype = names(predBlues),
                      "predictedBLUEs" = predBlues,
                      "seBLUEs" = seBlues,
                      "predictedBLUPs" = predBlups,
                      "seBLUPs" = seBlups,
                      "ue" = ue)
  # Extract variances
  genopos <- grep(pattern = "genotype!genotype.var", x = names(mr$gammas))
  varGen <- unname(mr$gammas[genopos] * mr$sigma2)
  resipos <- grep(pattern = "R!variance", x = names(mr$gammas))
  varErr <- unname(mr$gammas[resipos] * mr$sigma2)
  ## Compute Coefficient of Variation.
  CV <- 100 * summary(mf)$sigma / mean(fitted(mf))
  ## Create tempfile for asreml output.
  tmpfile <- tempfile()
  sink(file = tmpfile)
  ## Calculate wald test for genotype coefficients
  wtt <- asreml::wald.asreml(mf, ssType = "conditional", denDF = "numeric")
  pos <- grep(pattern = "genotype", x = row.names(wtt$Wald))
  chi2 <- wtt$Wald$F.con[pos] * wtt$Wald$Df[pos]
  prob <- 1 - pchisq(q = chi2, df = wtt$Wald$Df[pos])
  resWald <- list(chi2 = c(chi2 = chi2, df = wtt$Wald$Df[pos], P = prob),
                  Ftest = c(Fstat = wtt$Wald$F.con[pos],
                            df1 = wtt$Wald$Df[pos],
                            df2 = wtt$Wald$denDF[pos],
                            P = wtt$Wald$Pr[pos]))
  waldTestGeno <- list(result = resWald)
  ## Estimatie heritability on a line mean basis.
  newSed <- predict(mr, classify = "genotype", only = "genotype",
                    sed = TRUE, data = TD)$predictions$sed
  sink()
  unlink(tmpfile)
  sedSq <- newSed ^ 2
  heritability <- unname(1 - mean(sedSq[lower.tri(sedSq)]) / (varGen * 2))
  ## Calculate LSD; significance level (5%).
  lsd <- qt(p = .975, df = mf$nedf) * mf$predictions$avsed
  ## Combine results.
  genoLong <- TD$genotype
  result <- list(stats = stats, varGen = varGen, varErr = varErr,
                 heritability = heritability,
                 fitted = setNames(fitted(mf), genoLong),
                 resid = setNames(residuals(mf, type = "response"), genoLong),
                 stdRes = setNames(residuals(mf, type = "stdCond"), genoLong),
                 rMeans = setNames(fitted(mr), genoLong),
                 ranEf = setNames(mr$coe$random[grep(pattern = "genotype",
                                                     x = names(mr$coe$random))],
                                  genoNames),
                 waldTestGeno = waldTestGeno, CV = CV,
                 rDf = mf$nedf, predictionsSed = mf$predictions$avsed,
                 predictionsLsd = lsd)
  return(result)
}



