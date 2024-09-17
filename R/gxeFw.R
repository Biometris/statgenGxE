#' Finlay-Wilkinson analysis
#'
#' This function performs a Finlay-Wilkinson analysis of data classified by two
#' factors.
#'
#' @inheritParams gxeAmmi
#'
#' @param maxIter An integer specifying the maximum number of iterations in
#' the algorithm.
#' @param tol A positive numerical value specifying convergence tolerance of the
#' algorithm.
#' @param genotypes An optional character string containing the genotypes to
#' which the analysis should be restricted. If \code{NULL}, all genotypes are
#' used.
#' @param sorted A character string specifying the sorting order of the
#' estimated values in the output.
#'
#' @return An object of class \code{\link{FW}}, a list containing:
#' \item{estimates}{A data.frame containing the estimated values, with the
#' following columns:
#' \itemize{
#' \item genotype: The name of the genotype.
#' \item sens: The estimate of the sensitivity.
#' \item se_sens: The standard error of the estimate of the sensitivity.
#' \item genMean: The estimate of the genotypic mean.
#' \item se_genMean: The standard error of the estimate of the genotypic
#' mean.
#' \item MSdeviation: The mean square deviation about the line fitted to
#' each genotype
#' \item rank: The rank of the genotype based on its sensitivity.
#' }
#' }
#' \item{anova}{A data.frame containing anova scores of the FW analysis.}
#' \item{envEffs}{A data.frame containing the environmental effects, with the
#' following columns:
#' \itemize{
#' \item trial: The name of the trial.
#' \item envEff: The estimate of the environment effect.
#' \item se_envEff: The standard error of the estimate of the environment
#' effect.
#' \item envMean: The estimate of the environment mean.
#' \item rank: The rank of the trial based on its mean.
#' }
#' }
#' \item{TD}{The object of class TD on which the analysis was performed.}
#' \item{fittedGeno}{A numerical vector containing the fitted values for the
#' genotypes.}
#' \item{trait}{A character string containing the analyzed trait.}
#' \item{nGeno}{A numerical value containing the number of genotypes in the
#' analysis.}
#' \item{nEnv}{A numerical value containing the number of environments in the
#' analysis.}
#' \item{tol}{A numerical value containing the tolerance used during the
#' analysis.}
#' \item{iter}{A numerical value containing the number of iterations for the
#' analysis to converge.}
#'
#' @references Finlay, K.W. & Wilkinson, G.N. (1963). The analysis of adaptation
#' in a plant-breeding programme. Australian Journal of Agricultural
#' Research, 14, 742-754.
#'
#' @examples
#' ## Run Finlay-Wilkinson analysis on TDMaize.
#' geFW <- gxeFw(TDMaize, trait = "yld")
#'
#' ## Summarize results.
#' summary(geFW)
#'
#' ## Create a scatterplot of the results.
#' plot(geFW, plotType = "scatter")
#'
#' \donttest{
#' ## Create a report summarizing the results.
#' report(geFW, outfile = tempfile(fileext = ".pdf"))
#' }
#'
#' @family Finlay-Wilkinson
#'
#' @export
gxeFw <- function(TD,
                  trials = names(TD),
                  trait,
                  maxIter = 15,
                  tol = 0.001,
                  sorted = c("descending", "ascending", "none"),
                  genotypes = NULL,
                  useWt = FALSE) {
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  trials <- chkTrials(trials, TD)
  TDTot <- do.call(rbind, args = TD[trials])
  chkCol(trait, TDTot)
  chkCol("trial", TDTot)
  chkCol("genotype", TDTot)
  chkNum(maxIter, min = 1, null = FALSE, incl = TRUE)
  chkNum(tol, min = 0, null = FALSE)
  sorted <- match.arg(sorted)
  if (useWt) {
    chkCol("wt", TDTot)
  }
  chkChar(genotypes)
  if (!is.null(genotypes) && !all(genotypes %in% TDTot[["genotype"]])) {
    stop("All genotypes to include should be in TD.\n")
  }
  if (!is.null(genotypes)) {
    TDTot <- TDTot[TDTot[["genotype"]] %in% genotypes, ]
  }
  ## Remove genotypes that contain only NAs.
  allNA <- by(TDTot, TDTot[["genotype"]], FUN = function(x) {
    all(is.na(x[trait]))
  })
  TDTot <- TDTot[!TDTot[["genotype"]] %in% names(allNA[allNA]), ]
  ## Drop levels.
  TDTot <- droplevels(TDTot)
  ## Genotypes that are observed in only one trial will cause warnings when
  ## making predictions based on the fitted model.
  ## Check if the data contains such genotypes and give a user friendly
  ## warning here.
  genoTab <- table(unique(TDTot[c("genotype", "trial")]))
  genoOneObs <- rownames(genoTab)[rowSums(genoTab) == 1]
  if (length(genoOneObs) > 0) {
    warning("The following genotypes are present in only one trial. ",
            "This might give misleading predictions. Please check your data.\n",
            paste(genoOneObs, collapse = ", "), call. = FALSE)
  }
  ## Set wt to 1 if no weighting is used.
  if (!useWt) {
    TDTot[["wt"]] <- 1
  }
  nGeno <- nlevels(TDTot[["genotype"]])
  nEnv <- nlevels(TDTot[["trial"]])
  ## Setup empty vectors for storing rDev and rDF.
  rDev <- rDf <- rep(NA, 5)
  ## Compute total deviance.
  modelA <- lm(as.formula(paste0("`", trait, "`~genotype + trial + trial:genotype")),
               data = TDTot, weights = TDTot[["wt"]], na.action = na.exclude)
  aovA <- anova(modelA)
  rDev[5] <- sum(aovA[["Sum Sq"]])
  rDf[5] <- sum(aovA[["Df"]])
  ## Fit varieties only for first entry in aov.
  modelB <- lm(as.formula(paste0("`", trait, "`~-1 + genotype")), data = TDTot,
               weights = TDTot[["wt"]], na.action = na.exclude)
  aovB <- anova(modelB)
  rDev[1] <- aovB["Residuals", "Sum Sq"]
  rDf[1] <- aovB["Residuals", "Df"]
  ## Set a relative difference to be large.
  maxDiff <- Inf
  ## Set iteration to 0.
  iter <- 0
  TDTot$beta <- 1
  TDTot$envEffs <- NA
  ## Iterate to fit 'y(i,j) = genMean(i)+beta(i)*envEffs(j)'
  while (maxDiff > tol && iter <= maxIter) {
    beta0 <- TDTot[["beta"]]

    ## Fit model with current trial means relevant to each unit.
    model2 <- lm(as.formula(paste0("`", trait, "`~-1 + genotype + trial:beta")),
                 data = TDTot, weights = TDTot[["wt"]], na.action = na.exclude)
    coeffsModel2 <- coefficients(model2)
    ## Store residual ss and degrees of freedom after first iteration - beta = 1.
    if (iter == 0) {
      ## Environments.
      aov1 <- anova(model2)
      rDev[2] <- aov1["Residuals","Sum Sq"]
      rDf[2] <- aov1["Residuals", "Df"]
    }
    ## Update envEffs.
    envEffs <- coeffsModel2[match(paste0("trial", levels(TDTot[["trial"]]), ":beta"),
                                  names(coeffsModel2))]
    naPos <- is.na(envEffs)
    envEffs[naPos] <- 0
    envEffs <- envEffs - mean(envEffs)
    matchPos <- match(x = paste0("trial", TDTot[["trial"]], ":beta"),
                      table = names(envEffs))
    TDTot[["envEffs"]] <- envEffs[matchPos]
    ## Fit model with current genotype sensitivity relevant to each unit.
    model1 <- lm(as.formula(paste0("`", trait,
                                   "`~-1 + genotype + genotype:envEffs")),
                 data = TDTot, weights = TDTot[["wt"]], na.action = na.exclude)
    coeffsModel1 <- coefficients(model1)
    betas <- coeffsModel1[(nGeno + 1):(2 * nGeno)]
    ## Update beta.
    TDTot[["beta"]] <- betas[match(paste0("genotype", TDTot[["genotype"]],
                                          ":envEffs"),
                                   names(betas))]
    ## Compute max difference of sensitivities between the successive iterations.
    maxDiff <- max(abs(TDTot[["beta"]] - beta0), na.rm = TRUE)
    if (iter == maxIter && maxDiff > tol) {
      warning("Convergence not achieved in ", iter, " iterations. Tolerance ",
              tol, ", criterion at last iteration ", signif(maxDiff, 4), ".\n")
    }
    iter <- iter + 1
  }
  ## Store residual ss and degrees of freedom after last iteration - final model.
  aov2 <- anova(model1)
  rDev[4] <- aov2["Residuals","Sum Sq"]
  rDf[4] <- aov2["Residuals", "Df"]
  ## Calculate sums of squares and d.f. - calculations follow those in GenStat.
  rDev[c(3, 2, 1)] <- rDev[c(2, 1, 5)] - rDev[c(4, 2, 1)]
  rDf[c(2, 1)] <- rDf[c(1, 5)] - rDf[c(2, 1)]
  rDf[3] <- rDf[1]
  rDf[4] <- rDf[5] - rDf[1] - rDf[2] - rDf[3]
  ## Calculate mean deviances and F statistics.
  mDev <- rDev / rDf
  devr <- mDev / mDev[4]
  devr[c(4, 5)] <- NA
  devr[rDf == 0] <- NA
  devr[!is.na(devr) & devr < 0] <- NA
  fProb <- pf(q = devr, df1 = rDf, df2 = rDf[4], lower.tail = FALSE)
  ## Construct anova table.
  aovTable <- data.frame("Df" = rDf, "Sum Sq" = rDev, "Mean Sq" = mDev,
                         "F value" = devr, "Pr(>F)" = fProb,
                         row.names = c("Genotype", "Trial", "Sensitivities",
                                       "Residual", "Total"),
                         check.names = FALSE)
  aovTable <- aovTable[c("Trial", "Genotype", "Sensitivities", "Residual",
                         "Total"), ]
  class(aovTable) <- c("anova", "data.frame")
  ## Extract sensitivity beta.
  sens <- as.vector(tapply(X = TDTot[["beta"]], INDEX = TDTot[["genotype"]],
                           FUN = mean, na.rm = TRUE))
  ## Compute the standard errors.
  varE <- sqrt(diag(vcov(model1))[match(paste0("genotype", TDTot[["genotype"]],
                                               ":envEffs"),
                                        names(coeffsModel1))])
  sigmaE <- as.vector(tapply(X = varE, INDEX = TDTot[["genotype"]],
                             FUN = mean, na.rm = TRUE))
  ## Extract the mean for each genotype.
  genMean <- coeffsModel1[1:nGeno]
  ## Genotypic standard error.
  genSigma <- sqrt(diag(vcov(model1))[1:nGeno])
  ## Compute mean squared error (MSE) of the trait means for each genotype.
  mse <- as.vector(tapply(X = residuals(model1), INDEX = TDTot[["genotype"]],
                          FUN = function(x) {
                            checkG <- length(x)
                            if (checkG > 2) {
                              sum(x ^ 2, na.rm = TRUE) / (checkG - 2)
                            } else {
                              NA
                            }
                          }))
  ## Compute sorting order for estimates.
  if (sorted == "none") {
    orderSens <- 1:nGeno
  } else {
    orderSens <- order(sens, decreasing = (sorted == "descending"))
  }
  ## Construct estimate data.frame.
  estimates <- data.frame(Genotype = factor(levels(TDTot[["genotype"]]),
                                            labels = levels(TDTot[["genotype"]])),
                          GenMean = genMean,
                          SE_GenMean = genSigma,
                          Rank = rank(-sens),
                          Sens = sens,
                          SE_Sens = sigmaE,
                          MSdeviation = mse,
                          row.names = 1:length(sens))[orderSens, ]
  ## Construct data.frame with trial effects.
  vcovMod2 <- vcov(model2)
  vcovMod2[is.na(vcovMod2)] <- 0
  vc <- sqrt(diag(vcovMod2))
  vc <- vc[startsWith(x = names(vc), prefix = "trial")]
  vcEnv <- vcovMod2[startsWith(x = colnames(vcovMod2), prefix = "trial"),
                    startsWith(x = colnames(vcovMod2), prefix = "trial")]
  vm <- sum(vcEnv) / (nEnv ^ 2)
  ve <- vc ^ 2 + vm - 2 * rowSums(vcEnv) / nEnv
  seEnvEffs <- sqrt(ve)
  matchPos2 <- match(paste0("trial", levels(TDTot[["trial"]]), ":beta"),
                     names(envEffs))
  ## Create a full grid for making predictions.
  fullDat <- expand.grid(trial = levels(TDTot[["trial"]]),
                         genotype = levels(TDTot[["genotype"]]))
  fullDat[["envEffs"]] <- envEffs
  ## Predict will give a warning if there are genotypes that are present
  ## in only one trial. A more user friendly warning is shown earlier.
  ## Therefore suppress this specific warning if this is the case.
  if (length(genoOneObs) > 0) {
    supprWarn(predGeno <- predict(model1, se.fit = TRUE, newdata = fullDat),
              "has doubtful cases")
  } else {
    predGeno <- predict(model1, se.fit = TRUE, newdata = fullDat)
  }
  fittedGeno <- cbind(fullDat[c("trial", "genotype")],
                      fittedValue = predGeno$fit,
                      seFittedValue = predGeno$se.fit)
  meansFitted <- tapply(X = fittedGeno[["fittedValue"]],
                        INDEX = fittedGeno[["trial"]], FUN = mean, na.rm = TRUE)
  seMeansFitted <- tapply(X = fittedGeno[["seFittedValue"]],
                          INDEX = fittedGeno[["trial"]], FUN = mean, na.rm = TRUE)
  meansFitted <- meansFitted[matchPos2]
  seMeansFitted <- seMeansFitted[matchPos2]
  envEffsSummary <- data.frame(Trial = names(meansFitted),
                               EnvEff = envEffs,
                               SE_EnvEff = seEnvEffs,
                               EnvMean = as.vector(meansFitted),
                               SE_EnvMean = as.vector(seMeansFitted),
                               Rank = rank(-meansFitted), row.names = NULL)
  return(createFW(estimates = estimates, anova = aovTable,
                  envEffs = envEffsSummary, TD = createTD(TDTot),
                  fittedGeno = fittedGeno, trait = trait, nGeno = nGeno,
                  nEnv = nlevels(TDTot[["trial"]]), tol = tol, iter = iter - 1))
}
