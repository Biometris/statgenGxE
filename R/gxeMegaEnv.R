#' Form mega-environments based on winning genotypes from an AMMI model
#'
#' This function fits an AMMI model and then using the fitted values produces
#' a new factor clustering trials . This factor is added as a column megaEnv to
#' the input data. If a column megaEnv already exists this column is
#' overwritten with a warning.\cr\cr
#' The fitted values from the AMMI model can be used in two ways to determine
#' the mega environments. If \code{useWinGeno = TRUE} then for every trial the
#' genotype with the highest fitted value is extracted and trials with the same
#' genotype are put together in the same mega environment.\cr
#' If \code{useWinGeno = FALSE} instead of only using the best genotype for
#' determining the mega environments a proportion of the best genotypes is used
#' (indicated by \code{cutOff}). Trials are then clustered and the best number
#' of clusters is determined by minimizing the ratio of repeatabilities of the
#' line means of the trials and the mega environments. The number of clusters
#' minimizing this ratio determines the minal number of mega environments.#'
#'
#' @inheritParams gxeAmmi
#'
#' @param useWinGeno Should only the best genotype per trail be used for
#' determining mega environments? If \code{FALSE} clustering of trials is used.
#' @param method A character string indicating the criterion to determine
#' the best genotype per environment, either \code{"max"} or \code{"min"}.
#' @param cutOff A numerical value indicating the proportion of best genotypes
#' per locationh used in the calculation. I.e. a value of 0.8 indicates that
#' the best 80\% genotypes will be used.
#' @param sumTab Should a summary table be printed?
#'
#' @return The input object of class \code{\link{TD}} with an added extra
#' column megaEnv.
#'
#' @examples
#' ## Calculate mega-environments for TDMaize and print a summary of the results.
#' TDmegaEnv <- gxeMegaEnv(TD = TDMaize, trait = "yld")
#' ## Calculate new mega-environments based on the genotypes with the lowest
#' ## value per environment.
#' TDmegaEnv2 <- gxeMegaEnv(TD = TDmegaEnv, trait = "yld", method = "min")
#'
#' @references Atlin, G. N., R. J. Baker, K. B. McRae, and X. Lu. 2000.
#' Selection Response in Subdivided Target Regions. Crop Sci. 40:7-13.
#' https://doi:10.2135/cropsci2000.4017
#' @references Charter, R.A. & Alexander, R.A. Bull. Psychon. Soc. (1993)
#' 31: 123. https://doi.org/10.3758/BF03334158
#'
#' @importFrom methods getFunction
#' @importFrom utils combn
#' @export
gxeMegaEnv <- function(TD,
                       trials = names(TD),
                       trait,
                       useWinGeno = TRUE,
                       method = c("max", "min"),
                       cutOff = 0.8,
                       sumTab = TRUE) {
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
    stop(paste("TD should contain a column trial to be able to run an AMMI",
               "analysis.\n"))
  }
  if ("megaEnv" %in% colnames(TDTot)) {
    warning(paste("TD already contains a column megaEnv. This column will",
                  "be overwritten.\n"), call. = FALSE)
  }
  method <- match.arg(method)
  ## Remove genotypes that contain only NAs
  allNA <- by(TDTot, TDTot$genotype, FUN = function(x) {
    all(is.na(x[trait]))
  })
  TDTot <- TDTot[!TDTot$genotype %in% names(allNA[allNA]), ]
  rmYear <- FALSE
  if (!hasName(x = TDTot, name = "year")) {
    TDTot$year <- 0
    rmYear <- TRUE
  }
  ## Save and then drop factor levels.
  envLevels <- levels(TDTot$trial)
  TDTot$trial <- droplevels(TDTot$trial)
  if (useWinGeno) {
    ## Perform AMMI analysis.
    AMMI <- gxeAmmi(TD = createTD(TDTot), trait = trait, nPC = 2, byYear = TRUE)
    fitted <- AMMI$fitted
    ## Extract position of best genotype per trial.
    winPos <- apply(X = fitted, MARGIN = 2,
                    FUN = getFunction(paste0("which.", method)))
    ## Extract best genotype per trial.
    winGeno <- rownames(fitted)[winPos]
    ## Create factor based on best genotypes.
    megaFactor <- factor(winGeno, labels = "")
    ## Merge factor levels to original data.
    TDTot$megaEnv <- TDTot$trial
    levels(TDTot$megaEnv) <- as.character(megaFactor)
    ## Reapply saved levels to ensure input and output TDTot are identical.
    levels(TDTot$trial) <- envLevels
    ## If year was added, remove if before creating output.
    if (isTRUE(rmYear)) {
      TDTot <- TDTot[-which(colnames(TDTot) == "year")]
    }
    TDOut <- createTD(TDTot)
    ## Create summary table.
    summTab <- data.frame("Mega factor" = megaFactor, Trial = colnames(fitted),
                          "Winning genotype" = winGeno,
                          "AMMI estimates" =
                            fitted[matrix(c(winPos, 1:ncol(fitted)), ncol = 2)],
                          check.names = FALSE)
    summTab <- summTab[order(megaFactor), ]
    attr(TDOut, "sumTab") <- summTab
    if (sumTab) {
      printCoefmat(summTab)
    }
  } else {
    if (!hasName(x = TDTot, name = "loc")) {
      TDTot$loc <- TDTot$trial
    }
    ## Remove locations that appear only in one year.
    locYear <- table(TDTot$loc, TDTot$year)
    rmLocs <- rownames(locYear)[rowSums(locYear > 0) == 1]
    TDTot <- droplevels(TDTot[!TDTot$loc %in% rmLocs, ])
    ## One by one remove locations that don't appear within a year with at
    ## least one other location.
    ## This needs to be done in steps since otherwise in many case all
    ## locations will be removed.
    continue <- TRUE
    while (continue) {
      locs <- as.character(unique(TDTot$loc))
      locYear <- table(TDTot$loc, TDTot$year)
      ## Create a data.frame with the number of observations per combination
      ## of locations.
      locCom <- cbind(t(combn(x = rownames(locYear), m = 2)),
                      combn(x = rownames(locYear), m = 2, FUN = function(x) {
                        sum(pmin(locYear[x[[1]], ], locYear[x[[2]], ]))
                      }))
      ## While there is any 0 in this table continue the removal of locations.
      if (any(locCom[, 3] == 0)) {
        ## Create a table with number of zeros per locations.
        locZeros <- sapply(locs, function(loc) {
          length(locCom[(locCom[, 1] == loc | locCom[, 2] == loc) &
                          locCom[, 3] == 0])
        })
        ## Remove the location with most zeros from the data.
        TDTot <- droplevels(TDTot[TDTot$loc != names(which.max(locZeros)), ])
      } else {
        continue <- FALSE
      }
    }
    ## Perform AMMI analysis.
    AMMI <- gxeAmmi(TD = createTD(TDTot), trait = trait, nPC = NULL,
                    byYear = TRUE)
    ammiRaw <- reshape2::melt(AMMI$fitted, varnames = c("genotype", "trial"),
                              value.name = "ammiPred")
    ammiRaw <- merge(ammiRaw, TDTot[c("genotype", "trial", "loc", "year")])
    ## Compute quantile for determining best genotypes per year per location.
    quant <- tapply(X = ammiRaw$ammiPred, INDEX = list(ammiRaw$year, ammiRaw$loc),
                    FUN = quantile,
                    probs = ifelse(method == "max", cutOff, 1 - cutOff),
                    na.rm = TRUE, type = 1)
    ## Merge quantile to data.
    ammiQnt <- merge(ammiRaw, reshape2::melt(quant), by.x = c("year", "loc"),
                     by.y = c("Var1", "Var2"))
    ## Replace values by 0 when lower than quantile value and 1 otherwise.
    ## Reverse the equation if low values are better for current trait.
    ammiQnt$ammiPred <- if (method == "max") {
      as.numeric(ammiQnt$ammiPred > ammiQnt$value)
    } else {
      as.numeric(ammiQnt$ammiPred < ammiQnt$value)
    }
    ## Reshape data to get locations as header.
    ammiLoc <- reshape2::dcast(ammiQnt, genotype + year ~ loc,
                               value.var = "ammiPred", drop = FALSE)
    ## Compute means and standard deviations.
    ## For sd division by n should be used so this cannot be done by sd().
    Xi = tapply(ammiQnt$ammiPred, list(ammiQnt$year, ammiQnt$loc), FUN = mean)
    SXi = tapply(ammiQnt$ammiPred, list(ammiQnt$year, ammiQnt$loc),
                 FUN = function(x) {
                   sqrt(sum((x - mean(x)) ^ 2) / length(x))
                 })
    ## Compute correlations between locations per year.
    r0 <- by(data = ammiLoc[, c(3:ncol(ammiLoc))], INDICES = ammiLoc$year,
             FUN = cor, use = "pairwise.complete.obs")
    ## Form combinations of locations.
    combs <- combn(locs, m = 2)
    ## Compute correlations across years per genotype.
    ## If numbers of observations are similar per location x year this
    ## has very similar results to just applying cor with pairwise.complete.obs
    ## on the complete data set. However when the numbers start differing so
    ## do the reults.
    combs <- rbind(combs, mapply(FUN = combLocs, combs[1, ], combs[2, ],
                                 MoreArgs = list(ammi = ammiLoc, r0 = r0,
                                                 Xi = Xi, SXi = SXi)))
    ## Put computed correlations in lower half of correlation matrix.
    corMat <- matrix(nrow = length(locs), ncol = length(locs),
                     dimnames = list(locs, locs))
    corMat[lower.tri(corMat)] <- as.numeric(combs[3, ])
    ## Compute distances.
    distMat <- as.dist(sqrt(1 - corMat ^ 2))
    ## Cluster locations.
    tree <- hclust(distMat, method = "ward.D")
    clustRes <- clustGrRes <- data.frame()
    CRDRMin <- Inf
    ## Create tempfile for diverting asreml output
    tmp <- tempfile()
    ## Loop over number of clusters to determine best number of clusters
    ## by minimizing ratio CR/DR.
    for (k in 2:ceiling(length(locs) / 2)) {
      ## Extract cluster groups for k clusters.
      clustGr <- data.frame(megaEnv = cutree(tree, k = k))
      ## Merge cluster groups to data.
      modDat <- merge(TDTot, clustGr, by.x = "loc", by.y = "row.names")
      ## Model with regions, one varcomp.
      sink(file = tmp)
      ## Fit model with regions.
      modReg <- asreml::asreml(
        fixed = formula(paste(trait, "~ 1 + year * loc")),
        random = ~ genotype + genotype:megaEnv + genotype:loc + genotype:year +
          genotype:megaEnv:year, data = modDat, maxiter = 50)
      sink()
      ## Extract variance components.
      vcReg <- summary(modReg)$varcomp$component
      ## Compute number of locations.
      nL <- length(unique(modDat$loc))
      nY <- length(unique(modDat$year))
      ## Compute H2 over locations.
      H2Loc <- vcReg[1] / (vcReg[1] + vcReg[2] / k + vcReg[3] / nY +
                             vcReg[4] / (k * nY) + vcReg[5] / (k * nL))
      ## Compute H2 over regions.
      H2Reg <- (vcReg[1] + vcReg[2]) / (vcReg[1] + vcReg[2] + vcReg[3] / nY +
                                          vcReg[4] / nY + vcReg[5] / nL)
      ## Compute genetic correlation.
      rho <- vcReg[1] / sqrt(vcReg[1] * (vcReg[1] + vcReg[2]))
      ## Compute ratio CD/CR.
      CRDR <- rho * sqrt(H2Loc / H2Reg)
      if (CRDR < CRDRMin) {
        ## If CRDR is smaller than the previous minimum readjust values.
        clustGrRes <- clustGr
        clustRes <- modDat
        CRDRMin <- CRDR
      }
    }
    ## Throw away tempfile.
    unlink(tmp)
    ## If year was added, remove if before creating output.
    if (isTRUE(rmYear)) {
      clustRes <- clustRes[-which(colnames(clustRes) == "year")]
    }
    ## Create TD Output.
    TDOut <- createTD(clustRes)
    ## Attach cluster groups as attribute.
    clustGrRes <- clustGrRes[order(clustGrRes$megaEnv), , drop = FALSE]
    attr(TDOut, "sumTab") <- clustGrRes
    attr(TDOut, "CRDR") <- CRDRMin
    if (sumTab) {
      printCoefmat(clustGrRes)
    }
  }
  return(TDOut)
}

#' Helper function for computing combined correlations for a pair of variables.
#' Following the notation of Charter (1993).
#'
#' @keywords internal
combCor <- function(Xi,
                    Yi,
                    SXi,
                    SYi,
                    ni,
                    ri) {
  N <- sum(ni)
  X <- sum(ni * Xi)
  Y <- sum(ni * Yi)
  X2 <- sum(ni * (Xi ^ 2 + SXi ^ 2))
  Y2 <- sum(ni * (Yi ^ 2 + SYi ^ 2))
  XY <- sum(ni * (ri * SXi * SYi + Xi * Yi))
  r <- (N * XY - X * Y) / sqrt((N * X2 - X ^ 2) * (N * Y2 - Y ^ 2))
  return(r)
}

#' Helper function for computing the combined correlation for 2 locations.
#'
#' @keywords internal
combLocs <- function(l1,
                     l2,
                     ammi,
                     r0,
                     Xi,
                     SXi) {
  l1 <- as.character(l1)
  l2 <- as.character(l2)
  ## Compute number of observations for combination of l1 and l2 per year.
  nl1l2 <- table(ammi[!is.na(ammi[, l1]) & !is.na(ammi[, l2]), ]$year)
  ## Extract correlations for l1 and l2 per year for list of correlations.
  rl1l2 <- sapply(X = r0, FUN = "[", l1, l2)
  ## Determine which observations should be included in calculation of r.
  inclObs <- !is.na(rl1l2)
  ## Compute combined correlation.
  combCor(Xi = Xi[, l1][inclObs], Yi = Xi[, l2][inclObs],
          SXi = SXi[, l1][inclObs], SYi = SXi[, l2][inclObs],
          ni = nl1l2[inclObs], ri = rl1l2[inclObs])
}

