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
#' \item{genotype} {The name of the genotype.}
#' \item{sens} {The estimate of the sensitivity.}
#' \item{se_sens} {The standard error of the estimate of the sensitivity.}
#' \item{genMean} {The estimate of the genotypic mean.}
#' \item{se_genMean} {The standard error of the estimate of the genotypic
#' mean.}
#' \item{MSdeviation} {The mean square deviation about the line fitted to
#' each genotype}
#' \item{rank} {The rank of the genotype based on its sensitivity.}
#' }
#' }
#' \item{anova}{A data.frame containing anova scores of the FW analysis.}
#' \item{envEffs}{A data.frame containing the environmental effects, with the
#' following columns:
#' \itemize{
#' \item{trial} {The name of the trial.}
#' \item{envEff} {The estimate of the environment effect.}
#' \item{se_envEff} {The standard error of the estimate of the environment
#' effect.}
#' \item{envMean} {The estimate of the environment mean.}
#' \item{rank} {The rank of the trial based on its mean.}
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
#' @seealso \code{\link{FW}}, \code{\link{plot.FW}}, \code{\link{report.FW}}
#'
#' @examples
#' ## Run Finlay-Wilkinson analysis on TDMaize.
#' geFW <- gxeFw(TDMaize, trait = "yld")
#' ## Summarize results.
#' summary(geFW)
#' ## Create a scatterplot of the results.
#' plot(geFW, plotType = "scatter")
#' \dontrun{
#' ## Create a report summarizing the results.
#' report(geFW, outfile = "./testReports/reportFW.pdf")
#' }
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
  if (!is.character(trials) || !all(trials %in% names(TD))) {
    stop("All trials should be in TD.")
  }
  TDTot <- Reduce(f = rbind, x = TD[trials])
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !hasName(x = TDTot, name = trait)) {
    stop("trait has to be a column in TD.\n")
  }
  if (!"trial" %in% colnames(TDTot)) {
    stop("TD should contain a column trial to be able to run a Finlay
         Wilkinson analysis.\n")
  }
  if (is.null(maxIter) || !is.numeric(maxIter) || length(maxIter) > 1 ||
      maxIter != round(maxIter) || maxIter < 1) {
    stop("maxIter should be a positive integer.\n")
  }
  if (is.null(tol) || !is.numeric(tol) || length(tol) > 1 || tol < 0) {
    stop("tol should be a numerical value > 10^-6.\n")
  }
  if (!is.null(genotypes) && !all(genotypes %in% TDTot[["genotype"]])) {
    stop("All genotypes to include should be in TD.\n")
  }
  if (useWt && !hasName(x = TDTot, name = "wt")) {
    stop("wt has to be a column in TD when using weighting.")
  }
  if (!is.null(genotypes)) {
    TDTot <- TDTot[TDTot[["genotype"]] %in% genotypes, ]
  }
  ## Remove genotypes that contain only NAs
  allNA <- by(TDTot, TDTot$genotype, FUN = function(x) {
    all(is.na(x[trait]))
  })
  TDTot <- TDTot[!TDTot$genotype %in% names(allNA[allNA]), ]
  ## Drop levels to make sure prcomp doesn't crash.
  TDTot$genotype <- droplevels(TDTot$genotype)
  TDTot$trial <- droplevels(TDTot$trial)
  ## Set wt to 1 if no weighting is used.
  if (!useWt) {
    TDTot$wt <- 1
  }
  sorted <- match.arg(sorted)
  nGeno <- nlevels(TDTot$genotype)
  ## Setup empty vectors for storing rDev and rDF
  rDev <- rDf <- rep(NA, 5)
  ## Estimate trial effects with the sensitivity beta = 1.
  model0 <- lm(as.formula(paste0("`", trait, "`~-1 + trial + genotype")),
               data = TDTot, weights = TDTot$wt, na.action = na.exclude)
  aov0 <- anova(model0)
  rDev[2] <- aov0["Residuals", "Sum Sq"]
  rDf[2] <- aov0["Residuals", "Df"]
  coeffsModel0 <- coefficients(model0)
  ## Select coefficients for trials
  coeffsTr <- coeffsModel0[grep(pattern = "trial", x = names(coeffsModel0))]
  ## Center trial effects.
  envEffs0 <- scale(coeffsTr, scale = FALSE)
  ## Remove 'trial' from rownames and add column name.
  rownames(envEffs0) <- substring(rownames(envEffs0), first = 6)
  colnames(envEffs0) <- "envEffs"
  TDTot <- merge(x = TDTot, y = envEffs0, by.x = "trial", by.y = "row.names")
  ## Set initial values for sensitivity beta.
  TDTot$beta <- 1
  ## Set a relative difference to be large.
  maxDiff <- Inf
  ## Set iteration to 1.
  iter <- 1
  ## Iterate to fit 'y(i,j) = genMean(i)+beta(i)*envEffs(j)'
  while (maxDiff > tol && iter <= maxIter) {
    beta0 <- TDTot$beta
    ## Fit model with current genotype sensitivity relevant to each unit.
    model1 <- lm(as.formula(paste0("`", trait, "`~-1 + genotype + genotype:envEffs")),
                 data = TDTot, weights = TDTot$wt, na.action = na.exclude)
    coeffsModel1 <- coefficients(model1)
    ## Update beta.
    TDTot$beta <- coeffsModel1[match(paste0("genotype", TDTot$genotype,
                                            ":envEffs"), names(coeffsModel1))]
    TDTot$beta <- TDTot$beta / mean(TDTot$beta, na.rm = TRUE)
    ## Fit model with current trial means relevant to each unit.
    model2 <- lm(as.formula(paste0("`", trait, "`~-1 + trial:beta")),
                 data = TDTot, weights = TDTot$wt, na.action = na.exclude)
    coeffsModel2 <- coefficients(model2)
    ## Update envEffs.
    TDTot$envEffs <- coeffsModel2[match(paste0("trial", TDTot$trial, ":beta"),
                                     names(coeffsModel2))]
    TDTot[is.na(TDTot$envEffs), "envEffs"] <- 0
    TDTot$envEffs <- TDTot$envEffs - mean(TDTot$envEffs)
    ## Compute max difference of sensitivities between the succesive iterations.
    maxDiff <- max(abs(TDTot$beta - beta0), na.rm = TRUE)
    if (iter == maxIter && maxDiff > tol) {
      warning(paste0("Convergence not achieved in ", iter,
                     " iterations. Tolerance ", tol,
                     ", criterion at last iteration ", signif(maxDiff, 4),
                     ".\n"))
    }
    iter <- iter + 1
  }
  ## Environments.
  aov1 <- anova(model1)
  rDev[4] <- aov1["Residuals","Sum Sq"]
  rDf[4] <- aov1["Residuals", "Df"]
  ## Extract total deviance.
  modelA <- lm(as.formula(paste0("`", trait, "`~ genotype")), data = TDTot,
               weights = TDTot$wt, na.action = na.exclude)
  aovA <- anova(modelA)
  rDev[5] <- sum(aovA[["Sum Sq"]])
  rDf[5] <- sum(aovA[["Df"]])
  ## Fit varieties only for first entry in aov.
  modelB <- lm(as.formula(paste0("`", trait, "`~-1 + genotype")), data = TDTot,
               weights = TDTot$wt, na.action = na.exclude)
  aovB <- anova(modelB)
  rDev[1] <- aovB["Residuals", "Sum Sq"]
  rDf[1] <- aovB["Residuals", "Df"]
  ## Calculate deviances and d.f.
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
  aovTable <- data.frame("Df" = rDf, "Sum Sq" = rDev, "Mean Sq" = mDev,
                         "F value" = devr, "Pr(>F)" = fProb,
                         row.names = c("genotype", "trial", "Sensitivities",
                                       "Residual", "Total"),
                         check.names = FALSE)
  ## Extract sensitivity beta.
  sens <- as.vector(tapply(X = TDTot$beta, INDEX = TDTot$genotype,
                           FUN = mean, na.rm = TRUE))
  ## Compute the standard errors.
  varE <- sqrt(diag(vcov(model1))[match(paste0("genotype", TDTot$genotype,
                                               ":envEffs"),
                                        names(coeffsModel1))])
  sigmaE <- as.vector(tapply(X = varE, INDEX = TDTot$genotype,
                             FUN = mean, na.rm = TRUE))
  ## Extract the mean for each genotype.
  coeffsGen <- coeffsModel1[match(paste0("genotype", TDTot$genotype),
                                  names(coeffsModel1))]
  genMean <- as.vector(tapply(X = coeffsGen, INDEX = TDTot$genotype,
                    FUN = mean, na.rm = TRUE))
  ## Residual standard error.
  varG <- sqrt(diag(vcov(model1))[match(paste0("genotype", TDTot$genotype),
                                         names(coeffsModel1))])
  sigma <- as.vector(tapply(X = varG, INDEX = TDTot$genotype,
                            FUN = mean, na.rm = TRUE))
  ## Compute mean squared error (MSE) of the trait means for each genotype.
  mse <- as.vector(tapply(X = residuals(model1), INDEX = TDTot$genotype,
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
  estimates <- data.frame(genotype = levels(TDTot$genotype), sens,
                          se_sens = sigmaE, genMean, se_genMean = sigma,
                          MSdeviation = mse, rank = rank(-sens),
                          row.names = 1:length(sens))[orderSens, ]
  ## Construct data.frame with trial effects.
  matchPos <- match(paste0("trial", levels(TDTot$trial), ":beta"),
                    names(coeffsModel2))
  envEffs <- coeffsModel2[matchPos]
  naPos <- is.na(envEffs)
  envEffs[naPos] <- 0
  envEffs <- envEffs - mean(envEffs)
  seEnvEffs <- sqrt(diag(vcov(model2)[matchPos[!naPos], matchPos[!naPos]]))
  matchPos2 <- match(paste0("trial", levels(TDTot$trial), ":beta"),
                     names(envEffs))
  if (!is.null(model1$na.action)) {
    meansFitted <- tapply(X = model1$fitted,
                          INDEX = TDTot$trial[-model1$na.action], FUN = mean)
  } else {
    meansFitted <- tapply(X = model1$fitted, INDEX = TDTot$trial, FUN = mean)
  }
  meansFitted <- meansFitted[matchPos2]
  envEffsSummary <- data.frame(trial = names(meansFitted), envEff = envEffs,
                               se_envEff = seEnvEffs,
                               envMean = as.vector(meansFitted),
                               rank = rank(-meansFitted), row.names = NULL)
  return(createFW(estimates = estimates, anova = aovTable,
                  envEffs = envEffsSummary, TD = createTD(TDTot),
                  fittedGeno = unname(fitted(model1)), trait = trait,
                  nGeno = nGeno, nEnv = nlevels(TDTot$trial), tol = tol,
                  iter = iter - 1))
}
