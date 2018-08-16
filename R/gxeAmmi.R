#' AMMI analysis
#'
#' The Additive Main Effects and Multiplicative Interaction (AMMI) model fits
#' a model which involves the Additive Main effects (i.e. genotype and trial)
#' along with the Multiplicative Interaction effects. Then a principal component
#' analysis is done on the residuals (multiplicative interaction). This results
#' in an interaction characterized by Interaction Principal Components (IPCA)
#' enabling simultaneous plotting of genotypes and trials.\cr\cr
#' The parameter \code{nPC} is used to indicate the number of principal
#' components that is used in the principal component analysis (PCA). By setting
#' this parameter to \code{NA} the algorithm determines the best number of
#' principal components itself (see Details).\cr\cr
#' When setting the parameter \code{GGE} to \code{TRUE} instead of a regular
#' AMMI analysis a GGE analysis is performed. In this case instead of fitting
#' a model with genotype and trial as main effects only trial is used a a main
#' effect.\cr\cr
#' By specifying the parameter \code{asYear = true} a separate analysis will be
#' done for every year in the data. Combining the option with \code{nPC = NA}
#' may result in different numbers of principal components per year. The AMMI
#' estimates will still be returned as a single data.frame, but the other
#' results will be either lists or arrays.
#'
#' First a linear model \eqn{trait ~ genotype + trial} is fitted with both
#' genotype and trial fixed components in the model.\cr
#' The residuals from the fitted model are then used in a PCA. If \code{nPC} is
#' not \code{NA} a single PCA is done using \code{\link[stats]{prcomp}} with
#' maximum rank \code{nPC}.\cr
#' In case \code{nPC = NA} the PCA is first done with maximum rank 1. Then using
#' forward selection one by one the maximum rank is increased as long as the
#' added component is significant in the analysis.\cr
#' AMMI estimates are then computed using the results of the PCA.\cr
#'
#' @param TD An object of class \code{\link{TD}}.
#' @param trials A character string specifying the trials to be analyzed. If
#' not supplied all trials are used in the analysis.
#' @param trait A character string specifying the trait to be analyzed.
#' @param nPC An integer specifying the number of principal components used
#' as multiplicative term of genotype-by-trial interaction. If \code{NA} the
#' number of principal components is determined by the algorithm using
#' forward selection. See details.
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
#' If \code{byYear} = \code{TRUE} all returned items in the AMMI object except
#' \code{fitted} will consist of a list of results by year.
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
  ## Extract years and define empty objects for output.
  years <- sort(unique(TDTot$year.))
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
    if (!is.null(nPC) && (!is.numeric(nPC) || length(nPC) > 1 ||
        round(nPC) != nPC || nPC < 0 || nPC > min(nEnv, nGeno))) {
     stop("nPC should be an integer smaller than the number of trials.\n")
    }
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
      if (sum(is.na(y0)) / length(y0) > 0.3) {
        stop(ifelse(byYear, paste0("More than 30% missing values for ", year,
                                  ".\n"),
                    "More than 30% missing values.\n"))
      }
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
    rownames(aov)[rownames(aov) == "genotype"] <- "Genotype"
    rownames(aov)[rownames(aov) == "trial"] <- "Environment"
    rownames(aov)[rownames(aov) == "Residuals"] <- "Interactions"
    ## Compute principal components.
    if (!is.null(nPC)) {
      ## nPC is given. Use this in principal components analysis.
      pca <- prcomp(x = na.omit(resids), retx = TRUE, center = center,
                    scale. = scale, rank. = nPC)
      nPCYear <- nPC
    } else {
      ## nPC is not supplied. Do principal component analyses as long as
      ## when adding an extra component this new component is signifacant.
      pca <- prcomp(x = na.omit(resids), retx = TRUE, center = center,
                    scale. = scale, rank. = 1)
      for (i in 2:(nEnv - 2)) {
        pcaOrig <- pca
        pca <- prcomp(x = na.omit(resids), retx = TRUE, center = center,
                      scale. = scale, rank. = i)
        pcaAov <- pcaToAov(pca = pca, aov = aov)
        ## When there are no degrees of freedom left for the residual variance
        ## Pr(>F) will be nan. In this case revert to the previous number of
        ## components as well.
        if (is.nan(pcaAov[i, "Pr(>F)"]) || pcaAov[i, "Pr(>F)"] > 0.001) {
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
  } # End loop over years.
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
                    overallMean = ovMeanTot, byYear = byYear))
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


