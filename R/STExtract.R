#' Extracts statistics from the fitted model results
#'
#' This function is to extract and calculate various results such as
#' heritabilities, genotypic means, unit errors etc.
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

STExtract <- function(SSA) {
  engine <- SSA$engine
  ## Set useRepId to TRUE when it is used as fixed effect in the model.
  useRepId <- (SSA$design %in% c("res.ibd", "res.rowcol", "rcbd"))
  ## Extract statistics from fitted model.
  result <- do.call(what = paste0("extract", tools::toTitleCase(engine)),
                    args = list(SSA = SSA, useRepId = useRepId))
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
  result$trait <- SSA$trait
  result$model <- SSA$design
  result$engine <- engine
  return(result)
}


#' Extract statistics from model fitted using SpATS
#'
#' @keywords internal
extractSpATS <- function(SSA, useRepId) {
  mf <- SSA$mFix
  mr <- SSA$mMix
  TD <- SSA$data
  genoNames <- mf$terms$geno$geno_names
  ## Compute BLUEs and se of BLUEs from fixed model.
  predBlues <- setNames(SpATS::predict.SpATS(mf, which = "genotype")$predicted.values,
                        genoNames)
  seBlues <- SpATS::predict.SpATS(mf, which = "genotype")$standard.errors
  ## Compute BLUPs and se of BLUPs from fixed model.
  predBlups <- SpATS::predict.SpATS(mr, which = "genotype")$predicted.values
  seBlups <- SpATS::predict.SpATS(mr, which = "genotype")$standard.errors
  ## Collect predictions in data.frame.
  stats <- data.frame(genotype = names(predBlues),
                      "predictedBLUEs" = predBlues,
                      "seBLUEs" = seBlues,
                      "predictedBLUPs" = predBlups,
                      "seBLUPs" = seBlups)
  ## Compute generalized heritability.
  effDimGeno <- mr$eff.dim["genotype"]
  nGeno <- mr$dim["genotype"]
  heritability <- unname(effDimGeno / (nGeno - 1))
  ## Combine results.
  result <- list(stats = stats, varGen = unname(mr$var.comp["genotype"]),
                 varSpat = mr$var.comp[grep(pattern = "Coordinates",
                                            x = names(mr$var.comp))],
                 heritability = heritability,
                 fitted = setNames(fitted(mf), genoNames),
                 resid = setNames(residuals(mf), genoNames),
                 stdRes = NULL,
                 rMeans = setNames(fitted(mr), genoNames),
                 ranEf = mr$coeff[names(mr$coeff) %in% names(predBlues)],
                 rDf = unname(mf$dim["genotype"]),
                 effDims = mr$eff.dim)
  return(result)
}

#' Extract statistics from model fitted using lme4
#'
#' @keywords internal
extractLme4 <- function(SSA, useRepId) {
  mf <- SSA$mFix
  mr <- SSA$mMix
  TD <- SSA$data
  ## Extract BLUEs and se of BLUEs from fixed model using emmeans.
  em <- emmeans::emmeans(mf, specs = "genotype")
  emStats <- summary(em)
  genoNames <- emStats$genotype
  predBlues <- setNames(emStats$emmean, genoNames)
  seBlues <- emStats$SE
  ## Extract and invert variance covariance matrix.
  V <- vcov(em)
  Vinv <- try(chol2inv(chol(V)), silent = TRUE)
  ## Compute unit errors.
  if (!inherits(Vinv, "try-error")) {
    ue <- 1 / diag(Vinv)
  } else {
    ue <- 1 / diag(solve(V))
  }
  ## Extract fixed effects.
  fixEfMr <- lme4::fixef(mr)
  ## Get positions of 'repId', 'genotypes' and NA within fixed effects.
  repIdPos <- grep(pattern = "repId", x = names(fixEfMr))
  repIdMix <- fixEfMr[repIdPos]
  ## Compute Blups.
  blo <- mean(c(repIdMix, 0))
  predBlups <- coef(mr)$genotype[, "(Intercept)"] + blo
  ## Compute se of Blups.
  seBlups <- as.vector(sqrt(attr(lme4::ranef(mr, condVar = TRUE)[["genotype"]],
                                 "postVar")))
  ## Collect predictions in data.frame.
  stats <- data.frame(genotype = names(predBlues),
                      "predictedBLUEs" = predBlues,
                      "seBLUEs" = seBlues,
                      "predictedBLUPs" = predBlups,
                      "seBLUPs" = seBlups,
                      "ue" = ue)
  ## Calculate wald test for genotype coeffients
  waldTestGeno <- emmeans::test(emmeans::emmeans(mf, specs = "genotype"),
                                joint = TRUE)
  ## Compute Coefficient of Variation.
  CV <- 100 * summary(mf)$sigma / mean(fitted(mf))
  ## Compute standardized residuals.
  if (class(mf) == "lm") {
    stdRes <- setNames(rstandard(mf), genoNames)
  } else if (class(mf) == "lmerMod") {
    stdRes <- setNames(resid(mf, scaled = TRUE), genoNames)
  }
  ## Extract variances.
  varCor <- lme4::VarCorr(mr)
  varGen <- varCor[["genotype"]][1, 1]
  varErr <- attr(varCor, "sc") ^ 2
  ## Estimatie heritability on a line mean basis.
  if (useRepId) {
    heritability <- varGen / (varGen + (varErr / length(unique(TD$repId))))
  } else {
    heritability <- varGen / (varGen + varErr)
  }
  ## Combine results.
  result <- list(stats = stats, varGen = varGen, varErr = varErr,
                 heritability = heritability,
                 fitted = setNames(fitted(mf), genoNames),
                 resid = setNames(resid(mf), genoNames),
                 stdRes = stdRes,
                 rMeans = setNames(lme4::getME(mr, "mu"), genoNames),
                 ranEf = lme4::ranef(mr, drop = TRUE)[["genotype"]],
                 waldTestGeno = waldTestGeno, CV = CV,
                 rDf = df.residual(mf))
  return(result)
}

#' Extract statistics from model fitted using asreml
#'
#' @keywords internal
extractAsreml <- function(SSA, useRepId) {
  mf <- SSA$mFix
  mr <- SSA$mMix
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
  V <- mf$predictions$vcov
  Vinv <- try(chol2inv(chol(V)), silent = TRUE)
  ## Compute unit errors.
  if (!inherits(Vinv, "try-error")) {
    ue <- 1 / diag(Vinv)
  } else {
    ue <- 1 / diag(solve(V))
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
  result <- list(stats = stats, varGen = varGen, varErr = varErr,
                 heritability = heritability,
                 fitted = setNames(fitted(mf), genoNames),
                 resid = setNames(residuals(mf, type = "response"), genoNames),
                 stdRes = setNames(residuals(mf, type = "stdCond"), genoNames),
                 rMeans = setNames(fitted(mr), genoNames),
                 ranEf = setNames(mr$coe$random[grep(pattern = "genotype",
                                                     x = names(mr$coe$random))],
                                  genoNames),
                 waldTestGeno = waldTestGeno, CV = CV,
                 rDf = mf$nedf, predictionsSed = mf$predictions$avsed,
                 predictionsLsd = lsd)
  return(result)
}



