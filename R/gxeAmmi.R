#' AMMI analysis
#'
#' This function fits a model which involves the Additive Main effects (i.e.
#' genotype and environment) along with the Multiplicative Interaction effects
#' of principal component analysis (PCA).
#'
#' @param TD an object of class \code{\link{TD}}.
#' @param trait A character string specifying the trait to be analyzed.
#' @param nPC An integer specifying the number of principal components used
#' as multiplicative term of genotype-by-environment interaction.
#' @param center Should the variables be shifted to be zero centered?
#' @param scale Should the variables be scaled to have unit variance before
#' the analysis takes place?
#'
#' @return an object of class \code{\link{AMMI}}, a list containing
#' \item{envScores}{a matrix with environmental scores.}
#' \item{genoScores}{a matrix with genotypic scores.}
#' \item{importance}{a data.frame containing the importance of the principal
#' components.}
#' \item{anova}{a data.frame containing anova scores of the AMMI analysis.}
#' \item{fitted}{a matrix containing fitted values from the AMMI model.}
#' \item{trait}{a character vector containing the analyzed trait.}
#' \item{envMean}{a numerical vector containing the means per environment.}
#' \item{genoMean}{a numerical vector containing the means per genotype.}
#' \item{overallMean}{a numerical value containing the overall mean.}
#'
#' @examples
#' ## Run AMMI
#' geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
#' ## Create report
#' report(geAmmi, outfile = "./testReports/reportAmmi.pdf")
#'
#' @import stats
#' @export

gxeAmmi <- function(TD,
                    trait,
                    nPC = 2,
                    center = TRUE,
                    scale = FALSE) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TD)) {
    stop("trait has to be a column in TD.\n")
  }
  if (!"env" %in% colnames(TD)) {
    stop("TD should contain a column env to be able to run an AMMI analysis.\n")
  }
  ## Count number of genotypes, environments and traits.
  nGeno <- nlevels(droplevels(TD$genotype))
  nEnv <- nlevels(droplevels(TD$env))
  nTrait <- nrow(TD)
  ## At least 3 environments needed.
  if (nEnv < 3) {
    stop("TD should contain at least 3 environments to run the AMMI model.\n")
  }
  ## check if the supplied data contains the genotype by environment means.
  if (nTrait != nGeno * nEnv) {
    stop("TD should contain only 1 value per environment per genotype.\n")
  }
  ## Impute missing values
  if (any(is.na(TD[[trait]]))) {
    ## Transform data to genotype x environment matrix.
    y0 <- tapply(X = TD[[trait]], INDEX = TD[, c("genotype", "env")],
                 FUN = identity)
    yIndex <- tapply(X = 1:nTrait, INDEX = TD[, c("genotype", "env")],
                     FUN = identity)
    ## Actual imputation.
    y1 <- multMissing(y0, maxcycle = 10, na.strings = NA)
    ## Insert imputed values back into original data.
    TD[yIndex[is.na(y0)], trait] <- y1[is.na(y0)]
  }
  ## Fit linear model.
  model <- lm(as.formula(paste(trait, "~ genotype + env")), data = TD)
  ## Calculate residuals & fitted values of the linear model.
  resids <- tapply(X = resid(model), INDEX = TD[, c("genotype", "env")],
                   FUN = identity)
  fittedVals <- tapply(X = fitted(model), INDEX = TD[, c("genotype", "env")],
                       FUN = identity)
  # Compute principal components.
  pca <- prcomp(x = resids, retx = TRUE, center = center, scale. = scale,
                rank. = nPC)
  loadings <- pca$rotation
  scores <- pca$x
  ## Compute AMMI-estimates per genotype per environment.
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
  envMean <- tapply(X = TD[[trait]], INDEX = TD$env, FUN = mean)
  genoMean <- tapply(X = TD[[trait]], INDEX = TD$genotype, FUN = mean)
  overallMean <- mean(TD[[trait]])
  return(createAMMI(envScores = pca$rotation, genoScores = pca$x,
                    importance = importance, anova = anv, fitted = fitted,
                    trait = trait, envMean = envMean, genoMean = genoMean,
                    overallMean = overallMean))
}
