#' AMMI analysis
#'
#' The Additive Main Effects and Multiplicative Interaction (AMMI) model fits
#' a model which involves the Additive Main effects (i.e. genotype and trial)
#' along with Multiplicative Interaction effects. The additive effects are the
#' classical ANOVA main effects for genotype and environment, the
#' multiplicative effects follow from a principal component analysis on the
#' interaction residuals (= genotype by environment means after adjustment for
#' additive genotype and environment effects). This results in an interaction
#' characterized by Interaction Principal Components (IPCA) enabling
#' simultaneous plotting of genotypes and trials.\cr\cr
#' The parameter \code{nPC} is used to indicate the number of principal
#' components that is used in the principal component analysis (PCA). By setting
#' this parameter to \code{NULL} the algorithm determines the best number of
#' principal components (see Details).\cr\cr
#' By specifying the parameter \code{byYear = TRUE}, a separate analysis will be
#' done for every year in the data. Combining the option with \code{nPC = NULL}
#' may result in different numbers of principal components per year. The AMMI
#' estimates will still be returned as a single data.frame, but the other
#' results will be either lists or arrays.
#'
#' First a linear model \eqn{trait = genotype + trial + \epsilon} is fitted with
#' both genotype and trial fixed components in the model.\cr
#' The residuals from the fitted model are then used in a PCA. If \code{nPC} is
#' not \code{NULL} a single PCA is done using \code{\link[stats]{prcomp}} with
#' maximum rank \code{nPC}.\cr
#' In case \code{nPC = NULL}, the PCA is first done with one PC. Then using
#' forward selection one by one the number of PCs is increased as long as the
#' added component is significant in the analysis.\cr
#' AMMI estimates are then computed using the results of the PCA.\cr
#'
#' @param TD An object of class \code{\link{TD}}.
#' @param trials A character string specifying the trials to be analyzed. If
#' not supplied, all trials are used in the analysis.
#' @param trait A character string specifying the trait to be analyzed.
#' @param nPC An integer specifying the number of principal components used
#' as multiplicative term of genotype-by-trial interaction. If \code{NULL}, the
#' number of principal components is determined by the algorithm using
#' forward selection. See details.
#' @param byYear Should the analysis be done by year? If \code{TRUE} the data
#' is split by the variable year, analysis is performed and the results are
#' merged together and returned.
#' @param center Should the variables be shifted to be zero centered?
#' @param excludeGeno An optional character vector with names of genotypes to
#' be excluded from the analysis. If \code{NULL}, all genotypes are used.
#' @param useWt Should weighting be used when modeling? Requires a column
#' \code{wt} in \code{TD}.
#'
#' @return An object of class \code{\link{AMMI}}, a list containing:
#' \item{envScores}{A matrix with environmental scores.}
#' \item{genoScores}{A matrix with genotypic scores.}
#' \item{importance}{A data.frame containing the importance of the principal
#' components.}
#' \item{anova}{A data.frame containing anova scores of the AMMI analysis.}
#' \item{fitted}{A data.frame containing fitted values from the AMMI model.}
#' \item{trait}{A character string containing the analyzed trait.}
#' \item{envMean}{A numerical vector containing the environmental means.}
#' \item{genoMean}{A numerical vector containing the genotypic means.}
#' \item{overallMean}{A numerical value containing the overall mean.}
#' If \code{byYear} = \code{TRUE}, all returned items in the AMMI object except
#' \code{fitted} will consist of a list of results by year.
#'
#' @references Gauch H.G. (1992) Statistical Analysis of Regional Yield Trials:
#' AMMI Analysis of Factorial Designs. Elsevier, Amsterdam.
#' @references Yan, W., Kang, M. (2002). GGE Biplot Analysis. Boca Raton: CRC
#' Press.
#'
#' @examples
#' ## Run AMMI analysis on TDMaize.
#' geAmmi <- gxeAmmi(TD = TDMaize, trait = "yld")
#'
#' ## Summarize results.
#' summary(geAmmi)
#'
#' ## Create a biplot of genotypes and environment interaction with PC1 and PC2.
#' plot(geAmmi, plotType = "AMMI2")
#'
#' ## Create a pdf report summarizing the results.
#' \donttest{
#' report(geAmmi, outfile = tempfile(fileext = ".pdf"))
#' }
#'
#' @family AMMI
#'
#' @export
gxeAmmi <- function(TD,
                    trials = names(TD),
                    trait,
                    nPC = 2,
                    byYear = FALSE,
                    center = TRUE,
                    excludeGeno = NULL,
                    useWt = FALSE) {
  return(gxeAmmiHelp(TD = TD, trials = trials, trait = trait, nPC = nPC,
                     byYear = byYear, center = center,
                     excludeGeno = excludeGeno, useWt = useWt, GGE = FALSE))
}

#' GGE analysis
#'
#' The Genotype plus Genotype by Environment interaction (GGE) model fits
#' a model with trial as main fixed effect. Then a principal component
#' analysis is done on the residuals. This results in an interaction
#' characterized by Interaction Principal Components (IPCA)
#' enabling simultaneous plotting of genotypes and trials.\cr\cr
#' The parameter \code{nPC} is used to indicate the number of principal
#' components that is used in the principal component analysis (PCA). By setting
#' this parameter to \code{NULL} the algorithm determines the best number of
#' principal components (see Details).\cr\cr
#' By specifying the parameter \code{byYear = TRUE}, a separate analysis will be
#' done for every year in the data. Combining the option with \code{nPC = NULL}
#' may result in different numbers of principal components per year. The GGE
#' estimates will still be returned as a single data.frame, but the other
#' results will be either lists or arrays.
#'
#' First a linear model \eqn{trait = trial + \epsilon} is fitted with trial a
#' fixed component in the model.\cr
#' The residuals from the fitted model are then used in a PCA. If \code{nPC} is
#' not \code{NULL} a single PCA is done using \code{\link[stats]{prcomp}} with
#' maximum rank \code{nPC}.\cr
#' In case \code{nPC = NULL}, the PCA is first done with one PC. Then using
#' forward selection one by one the number of PCs is increased as long as the
#' added component is significant in the analysis.\cr
#' GGE estimates are then computed using the results of the PCA.\cr
#'
#' @inheritParams gxeAmmi
#'
#' @export
gxeGGE <- function(TD,
                   trials = names(TD),
                   trait,
                   nPC = 2,
                   byYear = FALSE,
                   center = TRUE,
                   excludeGeno = NULL,
                   useWt = FALSE) {
  return(gxeAmmiHelp(TD = TD, trials = trials, trait = trait, nPC = nPC,
                     byYear = byYear, center = center,
                     excludeGeno = excludeGeno, useWt = useWt, GGE = TRUE))
}

gxeAmmiHelp <- function(TD,
                        trials = names(TD),
                        trait,
                        nPC = 2,
                        byYear = FALSE,
                        center = TRUE,
                        GGE = FALSE,
                        excludeGeno = NULL,
                        useWt = FALSE) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  trials <- chkTrials(trials, TD)
  TDTot <- do.call(rbind, args = TD[trials])
  chkCol(trait, TDTot)
  chkCol("genotype", TDTot)
  chkCol("trial", TDTot)
  if (useWt) {
    chkCol("wt", TDTot)
  }
  if (byYear) {
    chkCol("year", TDTot)
  }
  chkNum(nPC, min = 1, incl = TRUE)
  chkChar(excludeGeno)
  if (!is.null(excludeGeno) && !all(excludeGeno %in% TDTot[["genotype"]])) {
    stop("All genotypes to exclude should be in TD.\n")
  }
  TDTot[["year."]] <- if (byYear) {
    as.character(TDTot[["year"]])
  } else {
    "0"
  }
  if (!is.null(excludeGeno)) {
    TDTot <- TDTot[!TDTot[["genotype"]] %in% excludeGeno, ]
    remGeno <- length(unique(TDTot[["genotype"]]))
    if (remGeno < 10) {
      warning("Less than 10 genotypes present in data.\n")
    }
  }
  ## Extract years and define empty objects for output.
  years <- sort(unique(TDTot[["year."]]))
  fitTot <- data.frame(genotype = unique(TDTot[["genotype"]]))
  loadTot <- scoreTot <- impTot <- aovTot <- envMeanTot <- genoMeanTot <-
    ovMeanTot <- datTot <- setNames(vector(mode = "list",
                                           length = length(years)), years)
  for (year in years) {
    TDYear <- TDTot[TDTot[["year."]] == year, ]
    ## Remove genotypes that contain only NAs
    allNA <- by(TDYear, TDYear[["genotype"]], FUN = function(x) {
      all(is.na(x[trait]))
    })
    TDYear <- TDYear[!TDYear[["genotype"]] %in% names(allNA[allNA]), ]
    ## Drop levels to make sure prcomp doesn't crash.
    TDYear <- droplevels(TDYear)
    ## Count number of genotypes, environments and traits.
    nGeno <- nlevels(TDYear[["genotype"]])
    nEnv <- nlevels(TDYear[["trial"]])
    nTrait <- nrow(TDYear)
    ## At least 3 trials needed.
    if (nEnv < 3) {
      if (byYear) {
        warning("There are less than 3 trials for ", year, ".\n",
                "Year ", year, " skipped.\n", call. =  FALSE)
        next
      } else {
        stop("TD should contain at least 3 trials to run the AMMI model.\n")
      }
    }
    ## check if the supplied data contains the genotype by environment means.
    maxTrGeno <- max(table(TDYear[["trial"]], TDYear[["genotype"]]))
    if (maxTrGeno > 1 || nTrait > nGeno * nEnv) {
      if (byYear) {
        warning("More than 1 value per trial per genotype for ", year, ".\n",
                "Year ", year, " skipped.\n", call. =  FALSE)
        next
      } else {
        stop("TD should contain at most 1 value per trial per genotype.\n")
      }
    }
    if (!is.null(nPC) && nPC >= min(nEnv, nGeno)) {
      if (byYear) {
        warning("nPC is larger than the number of trials for ", year, ".\n",
                "Year ", year, " skipped.\n", call. =  FALSE)
        next
      } else {
        stop("nPC should be smaller than the number of trials.\n")
      }
    }
    datTot[[year]] <- TDYear
    ## Add combinations of trial and genotype currently not in TD to TD.
    ## Reshape adds missing combinations to its output.
    ## Reshaping to wide and then back to long retains those missings.
    TDYear <- reshape(reshape(TDYear[c("trial", "genotype", trait)],
                              direction = "wide", timevar = "genotype",
                              idvar = "trial", v.names = trait))
    TDYear <- TDYear[order(TDYear[["genotype"]], TDYear[["trial"]]), ]
    ## Impute missing values.
    if (any(is.na(TDYear[[trait]]))) {
      ## Transform data to genotype x trial matrix.
      y0 <- tapply(X = TDYear[[trait]],
                   INDEX = TDYear[, c("genotype", "trial")], FUN = I)
      if (sum(is.na(y0)) / length(y0) > 0.3) {
        if (byYear) {
          warning("More than 30% missing values for ", year, ".\n",
                  "Year ", year, "skipped.\n", call. =  FALSE)
          next
        } else {
          stop("More than 30% missing values.\n")
        }
      }
      ## Actual imputation.
      y1 <- multMissing(y0, maxIter = 50)
      ## Insert imputed values back into original data.
      TDYear[is.na(TDYear[[trait]]), trait] <- t(y1)[is.na(t(y0))]
    }
    ## Set wt to 1 if no weighting is used.
    if (!useWt) {
      TDYear[["wt"]] <- 1
    }
    ## Fit linear model.
    modForm <- formula(paste0("`", trait, "`~",
                              if (!GGE) "genotype +", "trial"))
    model <- lm(modForm, data = TDYear, weights = TDYear[["wt"]])
    ## Calculate residuals and fitted values of the linear model.
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
                    scale. = FALSE, rank. = nPC)
      nPCYear <- nPC
    } else {
      ## nPC is not supplied. Do principal component analyses as long as
      ## when adding an extra component this new component is significant.
      pca <- prcomp(x = na.omit(resids), retx = TRUE, center = center,
                    scale. = FALSE, rank. = 2)
      if (nEnv > 4) {
        for (i in 3:max(3, (nEnv - 2))) {
          pcaOrig <- pca
          pca <- prcomp(x = na.omit(resids), retx = TRUE, center = center,
                        scale. = FALSE, rank. = i)
          pcaAov <- pcaToAov(pca = pca, aov = aov)
          ## When there are no degrees of freedom left for the residual variance
          ## Pr(>F) will be nan. In this case revert to the previous number of
          ## components as well.
          if (is.nan(pcaAov[i, "Pr(>F)"]) || pcaAov[i, "Pr(>F)"] > 0.001) {
            pca <- pcaOrig
            break
          }
        }
      }
      nPCYear <- ncol(pca$rotation)
    }
    ## Compute anova part for pca components.
    pcaAov <- pcaToAov(pca = pca, aov = aov)
    ## Create complete ANOVA table.
    aov <- rbind(aov, pcaAov)
    ## Extract loadings and scores from pca.
    loadings <- pca$rotation
    scores <- pca$x
    ## For GGE assure loadings for PC1 are positive to assure proper plotting.
    if (GGE && sum(loadings[, "PC1"] < 0) / nrow(loadings) > 0.6) {
      loadings[, "PC1"] <- -loadings[, "PC1"]
      scores[, "PC1"] <- -scores[, "PC1"]
    }
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
    envMean <- tapply(X = TDYear[[trait]], INDEX = TDYear[["trial"]],
                      FUN = mean)
    envMean <- setNames(as.numeric(envMean), names(envMean))
    genoMean <- tapply(X = TDYear[[trait]], INDEX = TDYear[["genotype"]],
                       FUN = mean)
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
  if (ncol(fitTot) > 1) {
    fitTot <- reshape(fitTot, direction = "long",
                      varying = list(genotype = colnames(fitTot)[-1]),
                      ids = fitTot[["genotype"]], idvar = "genotype",
                      times = colnames(fitTot)[-1], timevar = "trial",
                      v.names = "fittedValue")
    rownames(fitTot) <- NULL
  }
  if (!byYear) {
    loadTot <- loadTot[[1]]
    scoreTot <- scoreTot[[1]]
    impTot <- impTot[[1]]
    aovTot <- aovTot[[1]]
    envMeanTot <- envMeanTot[[1]]
    genoMeanTot <- genoMeanTot[[1]]
    ovMeanTot <- ovMeanTot[[1]]
    datTot <- datTot[[1]]
  } else {
    loadTot <- Filter(f = Negate(f = is.null), x = loadTot)
    scoreTot <- Filter(f = Negate(f = is.null), x = scoreTot)
    impTot <- Filter(f = Negate(f = is.null), x = impTot)
    aovTot <- Filter(f = Negate(f = is.null), x = aovTot)
    envMeanTot <- Filter(f = Negate(f = is.null), x = envMeanTot)
    genoMeanTot <- Filter(f = Negate(f = is.null), x = genoMeanTot)
    ovMeanTot <- Filter(f = Negate(f = is.null), x = ovMeanTot)
    datTot <- Filter(f = Negate(f = is.null), x = datTot)
    if (length(loadTot) == 0) {
      stop("All years were skipped.\n")
    }
  }
  return(createAMMI(envScores = loadTot, genoScores = scoreTot,
                    importance = impTot, anova = aovTot, fitted = fitTot,
                    trait = trait, envMean = envMeanTot, genoMean = genoMeanTot,
                    overallMean = ovMeanTot, dat = datTot, GGE = GGE,
                    byYear = byYear))
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

