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
#' @param byYear Should the analysis be done by year? If \code{TRUE} the data
#' is split by the variable year, analysis is performed and the results are
#' merged together and returned.
#' @param center Should the variables be shifted to be zero centered?
#' @param scale Should the variables be scaled to have unit variance?
#' @param GGE Should a GGE analysis be performed instead of a regular AMMI
#' analysis. When doing so genotype will be excluded from the model.
#' @param useWt Should weighting be used when modelling? Requires a column
#' \code{wt} in \code{TD}.
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
                    byYear = FALSE,
                    center = TRUE,
                    scale = FALSE,
                    GGE = FALSE,
                    useWt = FALSE) {
  ## Checks.
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
    stop("TD should contain a column trial to be able to run an AMMI analysis.\n")
  }
  if (useWt && !hasName(x = TDTot, name = "wt")) {
    stop("wt has to be a column in TD when using weighting.")
  }
  TDTot$year. <- if (byYear) {
    as.character(TDTot$year)
  } else {
    "0"
  }
  years <- unique(TDTot$year.)
  fitTot <- data.frame(genotype = unique(TDTot$genotype))
  loadTot <- scoreTot <- impTot <- aovTot <- envMeanTot <- genoMeanTot <-
    ovMeanTot <- setNames(vector(mode = "list", length = length(years)), years)
  for (year in years) {
    TDYear <- TDTot[TDTot$year. == year, ]
    ## Remove genotypes that contain only NAs
    allNA <- by(TDYear, TDYear$genotype, FUN = function(x) {
      all(is.na(x[trait]))
    })
    TDYear <- TDYear[!TDYear$genotype %in% names(allNA[allNA]), ]
    ## Drop levels to make sure prcomp doesn't crash.
    TDYear$genotype <- droplevels(TDYear$genotype)
    TDYear$trial <- droplevels(TDYear$trial)
    ## Count number of genotypes, environments and traits.
    nGeno <- nlevels(TDYear$genotype)
    nEnv <- nlevels(TDYear$trial)
    nTrait <- nrow(TDYear)
    ## At least 3 trials needed.
    if (nEnv < 3) {
      stop("TD should contain at least 3 trials to run the AMMI model.\n")
    }
    ## check if the supplied data contains the genotype by environment means.
    if (nTrait > nGeno * nEnv) {
      stop("TD should contain 1 value per trial per genotype.\n")
    }
    # if (!is.numeric(nPC) || length(nPC) > 1 || round(nPC) != nPC || nPC < 0 ||
    #     nPC > min(nEnv, nGeno)) {
    #   stop("nPC should be an integer smaller than the number of trials.\n")
    # }
    ## Add combinations of trial and genotype currently not in TD to TD.
    TDYear <- reshape2::melt(data = reshape2::dcast(data = TDYear,
                                                    formula = trial ~ genotype,
                                                    value.var = trait),
                             id.vars = "trial", variable.name = "genotype",
                             value.name = trait)
    ## Impute missing values
    if (any(is.na(TDYear[[trait]]))) {
      ## Transform data to genotype x trial matrix.
      y0 <- tapply(X = TDYear[[trait]],
                   INDEX = TDYear[, c("genotype", "trial")], FUN = identity)
      ## Actual imputation.
      y1 <- multMissing(y0, maxIter = 50)
      ## Insert imputed values back into original data.
      TDYear[is.na(TDYear[[trait]]), trait] <- y1[is.na(y0)]
    }
    ## Set wt to 1 if no weighting is used.
    if (!useWt) {
      TDYear$wt <- 1
    }
    ## Fit linear model.
    modForm <- formula(paste(trait, "~", if (!GGE) "genotype +", "trial"))
    model <- lm(modForm, data = TDYear, weights = TDYear$wt)
    ## Calculate residuals & fitted values of the linear model.
    resids <- tapply(X = resid(model), INDEX = TDYear[, c("genotype", "trial")],
                     FUN = I)
    fittedVals <- tapply(X = fitted(model),
                         INDEX = TDYear[, c("genotype", "trial")], FUN = I)
    ## Extract ANOVA table for linear model.
    aov <- anova(model)
    rownames(aov)[rownames(aov) == "Residuals"] <- "Interactions"
    ## Compute principal components.
    if (!is.na(nPC)) {
      pca <- prcomp(x = na.omit(resids), retx = TRUE, center = center,
                    scale. = scale, rank. = nPC)
      nPCYear <- nPC
    } else {
      pca <- prcomp(x = na.omit(resids), retx = TRUE, center = center,
                    scale. = scale, rank. = 1)
      for (i in 2:(nEnv - 2)) {
        pcaOrig <- pca
        pca <- prcomp(x = na.omit(resids), retx = TRUE, center = center,
                      scale. = scale, rank. = i)
        pcaAov <- pcaToAov(pca = pca, aov = aov)
        if (pcaAov[i, "Pr(>F)"] > 0.1) {
          pca <- pcaOrig
          break
        }
      }
      nPCYear <- ncol(pca$rotation)
    }
    ## Compute anova part for pca components.
    pcaAov <- pcaToAov(pca = pca, aov = aov)
    ## Create complete ANOVA table.
    aov <- rbind(aov, pcaAov)
    ## Extract loadings and scoress from pca.
    loadings <- pca$rotation
    scores <- pca$x
    ## Compute AMMI-estimates per genotype per trial.
    mTerms <- matrix(data = 0, nrow = nGeno, ncol = nEnv)
    for (i in 1:nPCYear) {
      mTerms <- mTerms + outer(scores[, i], loadings[, i])
    }
    fitted <- fittedVals + mTerms
    ## Extract importance from pca.
    importance <- as.data.frame(summary(pca)$importance)
    colnames(importance) <- paste0("PC", 1:ncol(importance))
    ## Compute means.
    envMean <- tapply(X = TDYear[[trait]], INDEX = TDYear$trial, FUN = mean)
    envMean <- setNames(as.numeric(envMean), names(envMean))
    genoMean <- tapply(X = TDYear[[trait]], INDEX = TDYear$genotype, FUN = mean)
    genoMean <- setNames(as.numeric(genoMean), names(genoMean))
    overallMean <- mean(TDYear[[trait]])
    fitTot <- merge(fitTot, fitted, by.x = "genotype", by.y = "row.names",
                    all.x = TRUE)
    loadTot[[year]] <- loadings
    scoreTot[[year]] <- scores
    impTot[[year]] <- importance
    aovTot[[year]] <- aov
    envMeanTot[[year]] <- envMean
    genoMeanTot[[year]] <- genoMean
    ovMeanTot[[year]] <- overallMean
  }
  rownames(fitTot) <- fitTot$genotype
  fitTot <- as.matrix(fitTot[-1])
  if (!byYear) {
    loadTot <- loadTot[[1]]
    scoreTot <- scoreTot[[1]]
    impTot <- impTot[[1]]
    aovTot <- aovTot[[1]]
    envMeanTot <- envMeanTot[[1]]
    genoMeanTot <- genoMeanTot[[1]]
    ovMeanTot <- ovMeanTot[[1]]
  }
  return(createAMMI(envScores = loadTot, genoScores = scoreTot,
                    importance = impTot, anova = aovTot, fitted = fitTot,
                    trait = trait, envMean = envMeanTot, genoMean = genoMeanTot,
                    overallMean = ovMeanTot))
}

#' @keywords internal
pcaToAov <- function(pca,
                     aov) {
  nPC <- ncol(pca$rotation)
  nEnv <- nrow(pca$rotation)
  nGeno <- nrow(pca$x)
  pcaAov <- matrix(data = NA, nrow = nPC + 1, ncol = 5,
                   dimnames = list(c(paste0("PC", 1:nPC), "Residuals"),
                                   colnames(aov)))
  ## Compute degrees of freedom and add to table.
  dfPC <- nGeno + nEnv - 3 - (2 * (1:nPC - 1))
  dfResid <- aov["Interactions", "Df"] - sum(dfPC)
  pcaAov[, "Df"] <- c(dfPC, dfResid)
  ## Compute sum of squares for PC and residuals and add to table.
  PCAVar <- pca$sdev ^ 2
  propVar <- PCAVar / sum(PCAVar)
  ssPC <- aov["Interactions", "Sum Sq"] * propVar[1:nPC]
  ssResid <- aov["Interactions", "Sum Sq"] - sum(ssPC)
  pcaAov[, "Sum Sq"] <- c(ssPC, ssResid)
  ## Compute mean squares for PC scores and residuals and add to table.
  pcaAov[, "Mean Sq"] <- pcaAov[, "Sum Sq"] / pcaAov[, "Df"]
  ## Convert infinite values to NA.
  pcaAov[, "Mean Sq"][is.infinite(pcaAov[, "Mean Sq"])] <- NA
  ## Add the F-values for PC scores.
  pcaAov[1:nPC, "F value"] <- pcaAov[1:nPC, "Mean Sq"] /
    pcaAov["Residuals", "Mean Sq"]
  ## Add the p-values for PC scores.
  pcaAov[1:nPC, "Pr(>F)"] <- 1 - pf(q = pcaAov[1:nPC, "F value"],
                                    df1 = dfPC, df2 = dfResid)
  return(pcaAov)
}


