#' AMMI analysis
#'
#' This function fits a model which involves the Additive Main effects (i.e.
#' genotype and trial) along with the Multiplicative Interaction effects
#' of principal component analysis (PCA).
#'
#' @param TD An object of class \code{\link{TD}}.
#' @param trials A character string specifying the trials to be analyzed. If
#' not supplied all trials are used in the analysis.
#' @param trait A character string specifying the trait to be analyzed.
#' @param nPC An integer specifying the number of principal components used
#' as multiplicative term of genotype-by-trial interaction.
#' @param center Should the variables be shifted to be zero centered?
#' @param scale Should the variables be scaled to have unit variance?
#'
#' @return An object of class \code{\link{AMMI}}, a list containing:
#' \item{envScores}{A matrix with environmental scores.}
#' \item{genoScores}{A matrix with genotypic scores.}
#' \item{importance}{A data.frame containing the importance of the principal
#' components.}
#' \item{anova}{A data.frame containing anova scores of the AMMI analysis.}
#' \item{fitted}{A matrix containing fitted values from the AMMI model.}
#' \item{trait}{A character string containing the analyzed trait.}
#' \item{envMean}{A numerical vector containing the environmental means.}
#' \item{genoMean}{A numerical vector containing the genotypic means.}
#' \item{overallMean}{A numerical value containing the overall mean.}
#'
#' @seealso \code{\link{AMMI}}, \code{\link{plot.AMMI}},
#' \code{\link{report.AMMI}}
#'
#' @references Gauch H.G. (1992) Statistical Analysis of Regional Yield Trials:
#' AMMI Analysis of Factorial Designs. Elsevier, Amsterdam.
#'
#' @examples
#' ## Run AMMI analysis on TDMaize.
#' geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
#' ## Summarize results.
#' summary(geAmmi)
#' ## Create a biplot of genotypes and environment interaction with PC1 and PC2.
#' plot(geAmmi, plotType = "AMMI2")
#' ## Create a pdf report summarizing the results.
#' \dontrun{
#' report(geAmmi, outfile = "./testReports/reportAmmi.pdf")
#' }
#'
#' @import stats
#' @export
gxeAmmi <- function(TD,
                    trials = names(TD),
                    trait,
                    nPC = 2,
                    center = TRUE,
                    scale = FALSE) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (!is.character(trials) || !all(trials %in% names(TD))) {
    stop("All trials should be in TD.")
  }
  TDTot <- Reduce(f = rbind, x = TD[trials])
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TDTot)) {
    stop("trait has to be a column in TD.\n")
  }
  if (!"trial" %in% colnames(TDTot)) {
    stop("TD should contain a column trial to be able to run an AMMI analysis.\n")
  }
  ## Dropping levels to make sure prcomp doesn't crash.
  TDTot$genotype <- droplevels(TDTot$genotype)
  TDTot$trial <- droplevels(TDTot$trial)
  ## Count number of genotypes, environments and traits.
  nGeno <- nlevels(TDTot$genotype)
  nEnv <- nlevels(TDTot$trial)
  nTrait <- nrow(TDTot)
  ## At least 3 trials needed.
  if (nEnv < 3) {
    stop("TD should contain at least 3 trials to run the AMMI model.\n")
  }
  ## check if the supplied data contains the genotype by environment means.
  if (nTrait != nGeno * nEnv) {
    stop("TD should contain 1 value per trial per genotype.\n")
  }
  if (!is.numeric(nPC) || length(nPC) > 1 || round(nPC) != nPC || nPC < 0 ||
      nPC > min(nEnv, nGeno)) {
    stop("nPC should be an integer smaller than the number of trials.\n")
  }
  ## Impute missing values
  if (any(is.na(TDTot[[trait]]))) {
    ## Transform data to genotype x trial matrix.
    y0 <- tapply(X = TDTot[[trait]], INDEX = TDTot[, c("genotype", "trial")],
                 FUN = identity)
    yIndex <- tapply(X = 1:nTrait, INDEX = TDTot[, c("genotype", "trial")],
                     FUN = identity)
    ## Actual imputation.
    y1 <- multMissing(y0, maxIter = 10)
    ## Insert imputed values back into original data.
    TDTot[yIndex[is.na(y0)], trait] <- y1[is.na(y0)]
  }
  ## Fit linear model.
  model <- lm(as.formula(paste(trait, "~ genotype + trial")), data = TDTot)
  ## Calculate residuals & fitted values of the linear model.
  resids <- tapply(X = resid(model), INDEX = TDTot[, c("genotype", "trial")],
                   FUN = identity)
  fittedVals <- tapply(X = fitted(model), INDEX = TDTot[, c("genotype", "trial")],
                       FUN = identity)
  # Compute principal components.
  pca <- prcomp(x = resids, retx = TRUE, center = center, scale. = scale,
                rank. = nPC)
  loadings <- pca$rotation
  scores <- pca$x
  ## Compute AMMI-estimates per genotype per trial.
  mTerms <- matrix(data = 0, nrow = nGeno, ncol = nEnv)
  for (i in 1:nPC) {
    mTerms <- mTerms + outer(scores[, i], loadings[, i])
  }
  fitted <- fittedVals + mTerms
  ## Extract ANOVA table for linear model.
  anv <- anova(model)
  rownames(anv)[rownames(anv) == "Residuals"] <- "Interactions"
  ## Create empty base table for extending anova table.
  addTbl <- matrix(data = NA, nrow = nPC + 1, ncol = 5,
                   dimnames = list(c(paste0("PC", 1:nPC), "Residuals"),
                                   colnames(anv)))
  ## Compute degrees of freedom and add to table.
  dfPC <- nGeno + nEnv - 3 - (2 * (1:nPC - 1))
  dfResid <- anv["Interactions", "Df"] - sum(dfPC)
  addTbl[, "Df"] <- c(dfPC, dfResid)
  ## Compute sum of squares for PC and residuals and add to table.
  PCAVar <- pca$sdev ^ 2
  propVar <- PCAVar / sum(PCAVar)
  ssPC <- anv["Interactions", "Sum Sq"] * propVar[1:nPC]
  ssResid <- anv["Interactions", "Sum Sq"] - sum(ssPC)
  addTbl[, "Sum Sq"] <- c(ssPC, ssResid)
  ## Compute mean squares for PC scores and residuals and add to table.
  addTbl[, "Mean Sq"] <- addTbl[, "Sum Sq"] / addTbl[, "Df"]
  ## Convert infinite values to NA.
  addTbl[, "Mean Sq"][is.infinite(addTbl[, "Mean Sq"])] <- NA
  ## Add the F-values for PC scores.
  addTbl[1:nPC, "F value"] <- addTbl[1:nPC, "Mean Sq"] /
    addTbl["Residuals", "Mean Sq"]
  ## Add the p-values for PC scores.
  addTbl[1:nPC, "Pr(>F)"] <- 1 - pf(q = addTbl[1:nPC, "F value"],
                                    df1 = dfPC, df2 = dfResid)
  ## Create complete ANOVA table.
  anv <- rbind(anv, addTbl)
  ## Extract importance from pca object.
  importance <- as.data.frame(summary(pca)$importance)
  colnames(importance) <- paste0("PC", 1:ncol(importance))
  ## Compute means.
  envMean <- tapply(X = TDTot[[trait]], INDEX = TDTot$trial, FUN = mean)
  envMean <- setNames(as.numeric(envMean), names(envMean))
  genoMean <- tapply(X = TDTot[[trait]], INDEX = TDTot$genotype, FUN = mean)
  genoMean <- setNames(as.numeric(genoMean), names(genoMean))
  overallMean <- mean(TDTot[[trait]])
  return(createAMMI(envScores = loadings, genoScores = scores,
                    importance = importance, anova = anv, fitted = fitted,
                    trait = trait, envMean = envMean, genoMean = genoMean,
                    overallMean = overallMean))
}
