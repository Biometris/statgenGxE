#' AMMI analysis
#'
#' The Additive Main Effects and Multiplicative Interaction (AMMI) model fits
#' a model which combines the Additive Main effects (i.e. genotype and trial)
#' along with the Multiplicative Interaction effects. Then a principal component
#' analysis is done on the residuals (multiplicative interaction). This results
#' in an interaction characterized by Interaction Principal Components (IPCA)
#' enabling simultaneous plotting of genotypes and trials.\cr\cr
#' The parameter \code{nPC} is used to indicate the number of principal
#' components that is used in the principal component analysis (PCA). By setting
#' this parameter to \code{NA} the algorithm determines the best number of
#' principal components (see Details).\cr\cr
#' By specifying the parameter \code{asYear = true}, a separate analysis will be
#' done for every year in the data. Combining the option with \code{nPC = NA}
#' may result in different numbers of principal components per year. The AMMI
#' estimates will still be returned as a single data.frame, but the other
#' results will be either lists or arrays.
#'
#' First a linear model \eqn{trait = genotype + trial + \epsilon} is fitted with
#' both genotype and trial fixed components in the model.\cr
#' The residuals from the fitted model are then used in a PCA. If \code{nPC} is
#' not \code{NA} a single PCA is done using \code{\link[stats]{prcomp}} with
#' maximum rank \code{nPC}.\cr
#' In case \code{nPC = NA}, the PCA is first done with one PC. Then using
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
#' @param scale Should the variables be scaled to have unit variance?
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
#' \item{fitted}{A matrix containing fitted values from the AMMI model.}
#' \item{trait}{A character string containing the analyzed trait.}
#' \item{envMean}{A numerical vector containing the environmental means.}
#' \item{genoMean}{A numerical vector containing the genotypic means.}
#' \item{overallMean}{A numerical value containing the overall mean.}
#' If \code{byYear} = \code{TRUE}, all returned items in the AMMI object except
#' \code{fitted} will consist of a list of results by year.
#'
#' @seealso \code{\link{AMMI}}, \code{\link{plot.AMMI}},
#' \code{\link{report.AMMI}}
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
#' \dontrun{
#' report(geAmmi, outfile = "./testReports/reportAmmi.pdf")
#' }
#'
#' @export
gxeAmmi <- function(TD,
                    trials = names(TD),
                    trait,
                    nPC = 2,
                    byYear = FALSE,
                    center = TRUE,
                    scale = FALSE,
                    excludeGeno = NULL,
                    useWt = FALSE,
                    maxIter = 10000,
                    tolerance = 1e-9) {
  return(gxeAmmiHelp(TD = TD, trials = trials, trait = trait, nPC = nPC,
                     byYear = byYear, center = center, scale = scale,
                     excludeGeno = excludeGeno, useWt = useWt,
                     maxIter = maxIter, tolerance = tolerance, GGE = FALSE))
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
#' this parameter to \code{NA} the algorithm determines the best number of
#' principal components (see Details).\cr\cr
#' By specifying the parameter \code{asYear = true}, a separate analysis will be
#' done for every year in the data. Combining the option with \code{nPC = NA}
#' may result in different numbers of principal components per year. The GGE
#' estimates will still be returned as a single data.frame, but the other
#' results will be either lists or arrays.
#'
#' First a linear model \eqn{trait = trial + \epsilon} is fitted with trial a
#' fixed component in the model.\cr
#' The residuals from the fitted model are then used in a PCA. If \code{nPC} is
#' not \code{NA} a single PCA is done using \code{\link[stats]{prcomp}} with
#' maximum rank \code{nPC}.\cr
#' In case \code{nPC = NA}, the PCA is first done with one PC. Then using
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
                   scale = FALSE,
                   excludeGeno = NULL,
                   useWt = FALSE,
                   maxIter = 10000,
                   tolerance = 1e-9) {
  return(gxeAmmiHelp(TD = TD, trials = trials, trait = trait, nPC = nPC,
                     byYear = byYear, center = center, scale = scale,
                     excludeGeno = excludeGeno, useWt = useWt,
                     maxIter = maxIter, tolerance = tolerance, GGE = TRUE))
}

gxeAmmiHelp <- function(TD,
                        trials = names(TD),
                        trait,
                        nPC = 2,
                        byYear = FALSE,
                        center = TRUE,
                        scale = FALSE,
                        GGE = FALSE,
                        excludeGeno = NULL,
                        useWt = FALSE,
                        maxIter = 10000,
                        tolerance = 1e-9) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  trials <- chkTrials(trials, TD)
  TDTot <- Reduce(f = rbind, x = TD[trials])
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
    ## Drop levels.
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
    TDYear <- reshape(reshape(TDYear[c("trial", "genotype", trait,
                                       if (useWt) "wt")],
                              direction = "wide", timevar = "genotype",
                              idvar = "trial",
                              v.names = c(trait, if (useWt) "wt")))
    if (useWt) {
      envMeans <- tapply(X = TDYear[[trait]] * TDYear[["wt"]],
                         INDEX = TDYear[["trial"]], FUN = mean, na.rm = TRUE)
      TDYear[is.na(TDYear[["wt"]]), "wt"] <- 0
    } else {
      envMeans <- tapply(X = TDYear[[trait]], INDEX = TDYear[["trial"]],
                         FUN = mean, na.rm = TRUE)
    }
    TDYear <- TDYear[order(TDYear[["trial"]], TDYear[["genotype"]]), ]
    ## Create matrix of observed phenotype.
    thetaHat <- matrix(TDYear[[trait]], nrow = nrow(TDYear))
    ## Construct weight matrix M.
    M <- diag(as.numeric(!is.na(thetaHat)))
    if (useWt) {
      M <- diag(x = TDYear[["wt"]])
    }
    ## Assure all values of M are between 0 and 1.
    M <- M / max(M)
    ## Replace missing values in thetaHat by (weighted) environment mean.
    thetaHat[is.na(thetaHat)] <- envMeans[TDYear[is.na(thetaHat), "trial"]]
    ## Construct design matrices.
    Xm <- matrix(data = 1, nrow = nEnv * nGeno)
    Xg <- matrix(data = 1, nrow = nEnv) %x% diag(x = nGeno)
    Xe <- diag(x = nEnv) %x% matrix(data = 1, nrow = nGeno)
    ## Initialize g, e and vecW.
    g <- matrix(data = 0, nrow = nGeno)
    e <- matrix(data = 0, nrow = nEnv)
    vecW <- matrix(data = 0, nrow = nEnv * nGeno)
    ## Create matrices containing coefficients for linear constraints.
    Pe <- matrix(data = 1, nrow = nEnv)
    Pg <- matrix(data = 1, nrow = nGeno)
    ## Pre compute inverses to speed up computations later on.
    tXgMXgInv <- solve(t(Xg) %*% M %*% Xg)
    tXeMXeInv <- solve(t(Xe) %*% M %*% Xe)
    ## Compute 'base' parts for mu, g and e that don't change across iterations.
    muBase <- diag(M) / sum(diag(M))
    gBase <- (diag(nGeno) - tXgMXgInv %*% Pg %*%
                solve(t(Pg) %*% tXgMXgInv %*% Pg) %*% t(Pg)) %*%
      tXgMXgInv %*% t(Xg) %*% M
    eBase <- (diag(nEnv) - tXeMXeInv %*% Pe %*%
                solve(t(Pe) %*% tXeMXeInv %*% Pe) %*% t(Pe)) %*%
      tXeMXeInv %*% t(Xe) %*% M
    ## Construct square matrices used later in the algorithm.
    fullMat <- diag(x = nEnv * nGeno) - M
    envMat <- matrix(data = 1 / nEnv, nrow = nEnv, ncol = nEnv)
    genoMat <- matrix(data = 1 / nGeno, nrow = nGeno, ncol = nGeno)
    if (!GGE) {
      mu0 <- muBase %*% thetaHat
      g0 <- gBase %*% (thetaHat - Xm %*% mu0)
      e0 <- eBase %*% (thetaHat - Xm %*% mu0)
      theta0 <- Xm %*% mu0 + Xe %*% e0 + Xg %*% g0
    } else {
      mu0 <- muBase %*% thetaHat
      e0 <- eBase %*% (thetaHat - Xm %*% mu0)
      theta0 <- Xm %*% mu0 + Xe %*% e0
    }
    ## Initialize loop parameters.
    i <- 1
    itDiff <- Inf
    thetaIter <- matrix(data = 0, nrow = nEnv * nGeno)
    while (i < maxIter && itDiff > tolerance) {
      if (!GGE) {
        ## Update mu, g and e.
        mu <- muBase %*% (thetaHat - Xe %*% e - Xg %*% g - vecW)
        g <- gBase %*% (thetaHat - Xm %*% mu - Xe %*% e - vecW)
        e <- eBase %*% (thetaHat - Xm %*% mu - Xg %*% g - vecW)
        ## Compute weighted residuals.
        wRes <- matrix(data = M %*% (thetaHat - Xm %*% mu - Xe %*% e - Xg %*% g) +
                        fullMat %*% vecW, nrow = nGeno, ncol = nEnv)
        ## Compute svd of weighted residuals.
        wSvd <- svd(wRes)
        ## Extract and reparameterize U and V from computed svd.
        U <- wSvd$u[, 1:nPC]
        U <- U - genoMat %*% U
        V <- wSvd$v[, 1:nPC]
        V <- V - envMat %*% V
        ## Compute vectorized version of W.
        vecW <- as.vector(U %*% diag(x = wSvd$d[1:nPC]) %*% t(V))
        ## Update theta.
        theta <- Xm %*% mu + Xe %*% e + Xg %*% g + vecW
      } else {
        ## Update mu and e.
        mu <- muBase %*% (thetaHat - Xe %*% e - vecW)
        e <- eBase %*% (thetaHat - Xm %*% mu - vecW)
        ## Compute weighted residuals.
        wRes <- matrix(M %*% (thetaHat - Xm %*% mu - Xe %*% e) +
                         fullMat %*% vecW, nrow = nGeno, ncol = nEnv)
        ## Extract and reparameterize U and V from computed svd.
        wSvd <- svd(wRes)
        U <- wSvd$u[, 1:nPC]
        U <- U - genoMat %*% U
        V <- wSvd$v[, 1:nPC]
        ## Compute vectorized version of W.
        vecW <- as.vector(U %*% diag(x = wSvd$d[1:nPC]) %*% t(V))
        ## Update theta.
        theta <- Xm %*% mu + Xe %*% e + vecW
      }
      itDiff <- sum((thetaIter - theta) ^ 2) / abs(mean(theta))
      cat("iteration", i, "ct =", itDiff , "\n")
      i <- i + 1
      thetaIter <- theta
    }
    ## Convert final vecW to a nGeno x nEnv matrix for a final svd computation.
    W <- matrix(data = vecW, nrow = nGeno, ncol = nEnv)
    ## Compute svd and extract U, V and D.
    wSvd <- svd(W)
    U <- wSvd$u[, 1:nPC]
    D <- wSvd$d[1:nPC]
    V <- wSvd$v[, 1:nPC]
    ## Construct ammi table.
    pcNames <- paste0("PC", 1:nPC)
    geNames <- c("Genotype", "Environment")
    aovAmmi <- data.frame(Df = c(nGeno - 1, nEnv - 1, (nGeno - 1) * (nEnv - 1),
                                 nGeno + nEnv - 1 - 2 * 1:nPC),
                          "Sum Sq" = c(nEnv * sum(g0 ^ 2), nGeno * sum(e0 ^ 2),
                                       sum((thetaHat - theta0) ^ 2), D ^ 2),
                          row.names = c(geNames, "Interactions", pcNames),
                          check.names = FALSE)
    aovAmmi["Residuals", c("Df", "Sum Sq")] <-
      aovAmmi["Interactions", c("Df", "Sum Sq")] -
      colSums(aovAmmi[pcNames, c("Df", "Sum Sq")])
    aovAmmi[["Mean Sq"]] <- aovAmmi[["Sum Sq"]] / aovAmmi[["Df"]]
    ## Convert infinite values to NA.
    aovAmmi[, "Mean Sq"][is.infinite(aovAmmi[, "Mean Sq"])] <- NA
    aovAmmi[geNames, "F value"] <- aovAmmi[geNames, "Mean Sq"] /
      aovAmmi["Interactions", "Mean Sq"]
    aovAmmi[pcNames, "F value"] <- aovAmmi[pcNames, "Mean Sq"] /
      aovAmmi["Residuals", "Mean Sq"]
    aovAmmi[geNames, "Pr(>F)"] <- 1 - pf(q = aovAmmi[geNames, "F value"],
                                         df1 = aovAmmi[geNames, "Df"],
                                         df2 = aovAmmi["Interactions", "Df"])
    aovAmmi[pcNames, "Pr(>F)"] <- 1 - pf(q = aovAmmi[pcNames, "F value"],
                                         df1 = aovAmmi[pcNames, "Df"],
                                         df2 = aovAmmi["Residuals", "Df"])
    class(aovAmmi) <- c("anova", "data.frame")
    ## Construct importance.
    UD <- U %*% diag(x = D)
    importance <- rbind(sqrt(apply(X = UD, MARGIN = 2, FUN = var)),
                        aovAmmi[pcNames, "Sum Sq"] /
                          aovAmmi["Interactions", "Sum Sq"])
    importance <- rbind(importance, cumsum(importance[2, ]))
    colnames(importance) <- pcNames
    rownames(importance) <- c("Standard deviation", "Proportion of Variance",
                              "Cumulative Proportion")
    ## Assign/Compute other return variables.
    envScores <- V
    genoScores <- U
    colnames(envScores) <- colnames(genoScores) <- pcNames
    rownames(envScores) <- levels(TDYear[["trial"]])
    rownames(genoScores) <- levels(TDYear[["genotype"]])
    fitted <- data.frame(trial = TDYear[["trial"]],
                         genotype = TDYear[["genotype"]],
                         fittedValue = theta)





    #     ## Compute anova part for pca components.
    #     pcaAov <- pcaToAov(pca = pca, aov = aov)
    #     ## Create complete ANOVA table.
    #     aov <- rbind(aov, pcaAov)
    #     ## Extract loadings and scores from pca.
    #     loadings <- pca$rotation
    #     scores <- pca$x
    #     ## For GGE assure loadings for PC1 are positive to assure proper plotting.
    #     if (GGE && all(loadings[, "PC1"] < 0)) {
    #       loadings[, "PC1"] <- -loadings[, "PC1"]
    #       scores[, "PC1"] <- -scores[, "PC1"]
    #     }
    #     ## Compute AMMI-estimates per genotype per trial.
    #     mTerms <- matrix(data = 0, nrow = nGeno, ncol = nEnv)
    #     for (i in 1:nPCYear) {
    #       mTerms <- mTerms + outer(scores[, i], loadings[, i])
    #     }
    #     fitted <- fittedVals + mTerms
    #     ## Extract importance from pca.
    #     importance <- as.data.frame(summary(pca)$importance)
    #     colnames(importance) <- paste0("PC", 1:ncol(importance))
    #     ## Compute means.
    #     envMean <- tapply(X = TDYear[[trait]], INDEX = TDYear[["trial"]],
    #                       FUN = mean)
    #     envMean <- setNames(as.numeric(envMean), names(envMean))
    #     genoMean <- tapply(X = TDYear[[trait]], INDEX = TDYear[["genotype"]],
    #                        FUN = mean)
    #     genoMean <- setNames(as.numeric(genoMean), names(genoMean))
    #     overallMean <- mean(TDYear[[trait]])
    #     fitTot <- merge(fitTot, fitted, by.x = "genotype", by.y = "row.names",
    #                     all.x = TRUE)
    #     loadTot[[year]] <- loadings
    #     scoreTot[[year]] <- scores
    #     impTot[[year]] <- importance
    #     aovTot[[year]] <- aov
    #     envMeanTot[[year]] <- envMean
    #     genoMeanTot[[year]] <- genoMean
    #     ovMeanTot[[year]] <- overallMean
    #   } # End loop over years.
    #   if (ncol(fitTot) > 1) {
    #     fitTot <- reshape(fitTot, direction = "long",
    #                       varying = list(genotype = colnames(fitTot)[-1]),
    #                       ids = fitTot[["genotype"]], idvar = "genotype",
    #                       times = colnames(fitTot)[-1], timevar = "trial",
    #                       v.names = "fittedValue")
    #     rownames(fitTot) <- NULL
  }
  #   if (!byYear) {
  #     loadTot <- loadTot[[1]]
  #     scoreTot <- scoreTot[[1]]
  #     impTot <- impTot[[1]]
  #     aovTot <- aovTot[[1]]
  #     envMeanTot <- envMeanTot[[1]]
  #     genoMeanTot <- genoMeanTot[[1]]
  #     ovMeanTot <- ovMeanTot[[1]]
  #     datTot <- datTot[[1]]
  #   } else {
  #     loadTot <- Filter(f = Negate(f = is.null), x = loadTot)
  #     scoreTot <- Filter(f = Negate(f = is.null), x = scoreTot)
  #     impTot <- Filter(f = Negate(f = is.null), x = impTot)
  #     aovTot <- Filter(f = Negate(f = is.null), x = aovTot)
  #     envMeanTot <- Filter(f = Negate(f = is.null), x = envMeanTot)
  #     genoMeanTot <- Filter(f = Negate(f = is.null), x = genoMeanTot)
  #     ovMeanTot <- Filter(f = Negate(f = is.null), x = ovMeanTot)
  #     datTot <- Filter(f = Negate(f = is.null), x = datTot)
  #     if (length(loadTot) == 0) {
  #       stop("All years were skipped.\n")
  #     }
  #   }


  return(createAMMI(envScores = envScores, genoScores = genoScores,
                    importance = importance, anova = aovAmmi, fitted = fitted,
                    trait = trait, envMean = as.numeric(e) + as.numeric(mu),
                    genoMean = as.numeric(g) + as.numeric(mu),
                    overallMean = mu, dat = datTot, GGE = GGE,
                    byYear = byYear))
}



