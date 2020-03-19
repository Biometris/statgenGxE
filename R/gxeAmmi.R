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
#' @param maxIter An integer specifying the maximum number of iterations in the
#' algorithm.
#' @param tolerance A numeric value indicating the tolerance for convergence in
#' the algorithm.
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
    TD0 <- expand.grid(trial = levels(TDYear[["trial"]]),
                       genotype = levels(TDYear[["genotype"]]))
    TDYear <- merge(TD0, TDYear, all.x = TRUE)
    if (useWt) {
      TDYear[is.na(TDYear[[trait]]), "wt"] <- 0
      ## Divide by max value to get all values in 0 to 1 range.
      TDYear[["wt"]] <- TDYear[["wt"]] / max(TDYear[["wt"]])
    } else {
      ## 1 for non-missing trait values, 0 for missing trait values.
      TDYear[["wt"]] <- as.numeric(!is.na(TDYear[[trait]]))
    }

    TDYear$wt <- TDYear$wt

    ## Fit a base additive model and get predictions.
    baseFit <- lm(formula(paste(trait, "~ genotype + trial")), data = TDYear,
                  weights = TDYear[["wt"]])
    ## Create matrix of observed phenotype.
    thetaHat <- matrix(TDYear[[trait]], nrow = nrow(TDYear))
    ## Replace missing values in thetaHat by 0.
    ## The actual replacement value doesn't matter since these observations
    ## get weight 0 and thus won't be used in the algoritm.
    thetaHat[is.na(thetaHat)] <- 0
    ## Construct weight matrix M from wt column in data.
    M <- diag(x = TDYear[["wt"]])
    ## Construct design matrices.
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
      mu0 <- as.numeric(muBase %*% thetaHat)
      g0 <- gBase %*% (thetaHat - mu0)
      e0 <- eBase %*% (thetaHat - mu0)
      theta0 <- mu0 + Xe %*% e0 + Xg %*% g0
    } else {
      mu0 <- as.numeric(muBase %*% thetaHat)
      e0 <- eBase %*% (thetaHat - mu0)
      theta0 <- mu0 + Xe %*% e0
      g0 <- 0
    }
    ## Initialize loop parameters.
    i <- 1
    itDiff <- scoreDiffE <- scoreDiffG <- Inf
    thetaIter <- matrix(data = 0, nrow = nEnv * nGeno)
    eIter <- e
    gIter <- g
    theta <- thetaHat
    nPC0 <- min(nEnv, nGeno)
    decrease <- FALSE
    while ((i < maxIter && itDiff > tolerance &&
            (scoreDiffE > tolerance || scoreDiffG > tolerance)) ||
           is.null(nPC)) {
      if (is.null(nPC) && (nPC0 == 5 ||
                           (nPC0 < 5 && nPC0 == min(nEnv - 1, nGeno - 1)))) {
        nPC0 <- nPC <- testPPB(theta = matrix(theta, nrow = nGeno, ncol = nEnv),
                               GGE = GGE)$K
      }
      if (decrease && (is.null(nPC) || nPC0 > nPC)) {
        nPC0 <- nPC0 - 1
      }
      if (!GGE) {
        ## Update mu, g and e.
        mu <- as.numeric(muBase %*% (thetaHat - Xe %*% e - Xg %*% g - vecW))
        e <- eBase %*% (thetaHat - mu - Xg %*% g - vecW)
        g <- gBase %*% (thetaHat - mu - Xe %*% e - vecW)
        ## Compute weighted residuals.
        wRes <- matrix(data = M %*% (thetaHat - mu - Xe %*% e - Xg %*% g) +
                         fullMat %*% vecW, nrow = nGeno, ncol = nEnv)
        ## Compute svd of weighted residuals.
        wSvd <- La.svd(wRes)
        ## Extract and reparameterize U and V from computed svd.
        U <- wSvd$u[, 1:nPC0]
        U <- U - genoMat %*% U
        Vt <- wSvd$vt[1:nPC0, ]
        Vt <- Vt - Vt %*% envMat
        ## Compute vectorized version of W.
        vecW <- as.vector(U %*% diag(x = wSvd$d[1:nPC0],
                                     nrow = nPC0, ncol = nPC0) %*% Vt)
        ## Update theta.
        theta <- mu + Xe %*% e + Xg %*% g + vecW
      } else {
        ## Update mu and e.
        mu <- as.numeric(muBase %*% (thetaHat - Xe %*% e - vecW))
        e <- eBase %*% (thetaHat - mu - vecW)
        ## Compute weighted residuals.
        wRes <- matrix(M %*% (thetaHat - mu - Xe %*% e) +
                         fullMat %*% vecW, nrow = nGeno, ncol = nEnv)
        ## Extract and reparameterize U and V from computed svd.
        wSvd <- La.svd(wRes)
        U <- wSvd$u[, 1:nPC0]
        U <- U - genoMat %*% U
        Vt <- wSvd$vt[1:nPC0, ]
        ## Compute vectorized version of W.
        vecW <- as.vector(U %*% diag(x = wSvd$d[1:nPC0],
                                     nrow = nPC0, ncol = nPC0) %*% Vt)
        ## Update theta.
        theta <- mu + Xe %*% e + vecW
      }
      itDiff <- sum((thetaIter - theta) ^ 2) / abs(mean(theta))
      scoreDiffE <- sum((eIter - e) ^ 2) / mean(abs(e))
      scoreDiffG <- sum((gIter - g) ^ 2) / mean(abs(g))
      cat("iteration", i, "ct =", itDiff , "scoreDiffE =", scoreDiffE,
          "scoreDiffG =", scoreDiffG, "\n")
      i <- i + 1
      thetaIter <- theta
      eIter <- e
      gIter <- g
      if (!decrease && itDiff < 1e-2) {
        decrease <- TRUE
      }
    }
    ## Convert final vecW to a nGeno x nEnv matrix for a final svd computation.
    W <- matrix(data = vecW, nrow = nGeno, ncol = nEnv)
    ## Compute svd and extract U, V and D.
    wSvd <- svd(W)
    U <- wSvd$u[, 1:nPC, drop = FALSE]
    D <- wSvd$d[1:nPC]
    V <- wSvd$v[, 1:nPC, drop = FALSE]
    ## Construct ammi table.
    pcNames <- paste0("PC", 1:nPC)
    geNames <- c("Genotype", "Environment")
    # aovAmmi <- data.frame(Df = c(nGeno - 1, nEnv - 1, (nGeno - 1) * (nEnv - 1),
    #                              nGeno + nEnv - 1 - 2 * 1:nPC),
    #                       "Sum Sq" = c(nEnv * sum(g ^ 2), nGeno * sum(e ^ 2),
    #                                    sum(residuals(baseFit) ^ 2), D ^ 2),
    #                       row.names = c(geNames, "Interactions", pcNames),
    #                       check.names = FALSE)

    aovAmmi <- anova(baseFit)
    aovAmmi$`Sum Sq` <- aovAmmi$`Sum Sq` * length(diag(M)) / sum(diag(M))
    rownames(aovAmmi) <- c("Genotype", "Environment", "Interactions")
    aovAmmi[pcNames, "Df"] <- nGeno + nEnv - 1 - 2 * 1:nPC
    aovAmmi[pcNames, "Sum Sq"] <- D ^ 2
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
    UD <- U %*% diag(x = D, nrow = nPC, ncol = nPC)
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
    ## For GGE assure loadings for PC1 are positive to assure proper plotting.
    if (GGE && all(envScores[, "PC1"] < 0)) {
      envScores[, "PC1"] <- -envScores[, "PC1"]
      genoScores[, "PC1"] <- -genoScores[, "PC1"]
    }
    # fitTot <- merge(fitTot, fitted, by.x = "genotype", by.y = "row.names",
    #                 all.x = TRUE)
    fitTot <- fitted
    loadTot[[year]] <- envScores
    scoreTot[[year]] <- genoScores
    impTot[[year]] <- importance
    aovTot[[year]] <- aovAmmi
    envMeanTot[[year]] <- as.numeric(e) + mu
    genoMeanTot[[year]] <- as.numeric(g) + mu
    ovMeanTot[[year]] <- mu
  } # End loop over years.
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
                    trait = trait, envMean = envMeanTot,
                    genoMean = genoMeanTot, overallMean = ovMeanTot,
                    dat = datTot, GGE = GGE, byYear = byYear))
}

#' @keywords internal
testPPB <- function(theta,
                    nBoot = 1000,
                    GGE = FALSE)  {
  nGeno <- nrow(theta)
  nEnv <- ncol(theta)
  KMax <- if (!GGE) min(nGeno - 1, nEnv - 1) else min(nGeno, nEnv - 1)
  E <- sweep(x = theta, MARGIN = 1, STATS = rowMeans(theta))
  E <- sweep(x = E, MARGIN = 2, STATS = colMeans(theta))
  E <- E + mean(theta)
  ESvd <- svd(E)
  U <- ESvd$u[, 1:KMax, drop = FALSE]
  V <- ESvd$v[, 1:KMax, drop = FALSE]
  lam <- ESvd$d[1:KMax]
  k <- 0
  pValue <- 0
  while (k < KMax && pValue < 0.05) {
    TObs <- lam[k + 1] ^ 2 / sum(lam[(k + 1):KMax] ^ 2)
    if (k > 0) {
      thetaK <- U[, 1:k, drop = FALSE] %*%
        diag(x = lam[1:k], nrow = k, ncol = k) %*% t(V[, 1:k, drop = FALSE])
    } else {
      thetaK <- matrix(data = 0, nrow = nGeno, ncol = nEnv)
    }
    RB <- U[, (k + 1):KMax, drop = FALSE] %*%
      diag(x = lam[(k+1):KMax], nrow = KMax - k, ncol = KMax - k) %*%
      t(V[, (k + 1):KMax, drop = FALSE])
    TBoot <- sapply(1:nBoot, FUN = function(i) {
      Rb <- sample(RB, nGeno * nEnv, replace = FALSE)
      Eb <- thetaK + Rb
      Ebb <- sweep(x = Eb, MARGIN = 1, STATS = rowMeans(Eb))
      Ebb <- sweep(x = Ebb, MARGIN = 2, STATS = colMeans(Eb))
      Ebb <- Ebb + mean(Eb)
      lamb <- svd(Ebb)$d
      lamb[k + 1] ^ 2 / sum(lamb[(k + 1):KMax] ^ 2)
    })
    pValue <- mean(TBoot > TObs)

    cat(k, pValue, "\n")
    k <- k + 1
  }
  return(list(K = if (pValue < 0.05) k else k - 1, pValue = pValue))
}


