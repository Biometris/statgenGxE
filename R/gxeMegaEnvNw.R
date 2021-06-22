#' Form mega environments based on fitted values from an AMMI model
#'
#' This function fits an AMMI model and then using the fitted values produces
#' a new factor clustering the trials. This factor is added as a column megaEnv
#' to the input data. If a column megaEnv already exists this column is
#' overwritten with a warning.\cr\cr
#' Mega environments can be created by two methods. The first method
#' (\code{useWinGeno = TRUE}) groups
#' environments based on their best performing genotype; i.e. environments that
#' share the same best genotype belong to the same mega environment, regardless
#' whether environments correspond to years or locations.\cr\cr
#' In the second method (\code{useWinGeno = FALSE}),
#' genotypes that are above a certain quantile are used to classify locations
#' into mega environments that are consistent across years. In this method,
#' genotypes are scored according to a \code{cutOff} threshold for the genotypic
#' ranking within each location (one if a genotype is above the cutOff and zero
#' if below). This genotype by location matrix with ones and zeros is used to
#' calculate the correlation between locations. Then, correlations across years
#' are combined using the method by Charter and Alexander (1993). The combined
#' correlations are used to calculate Euclidean distances for hierarchical
#' clustering. The number of mega environments obtained with the hierarchical
#' clustering procedure is chosen to maximize the correlated response to
#' selection within mega environments, as proposed in Atlin et al (2000).
#'
#' @inheritParams gxeAmmi
#'
#' @param useWinGeno Should only the best genotype per trial be used for
#' determining mega environments? If \code{FALSE} hierarchical clustering is
#' applied to classify locations into mega environments that are consistent
#' across years.
#' @param method A character string indicating the criterion to determine
#' the best genotype per environment, either \code{"max"} or \code{"min"}.
#' @param cutOff A numerical value indicating the proportion of best genotypes
#' per location used in the calculation, i.e. a value of 0.8 indicates that
#' the best 80% genotypes will be used.
#' @param sumTab Should a summary table be printed?
#'
#' @return The input object of class \code{\link{TD}} with an added extra
#' column megaEnv.
#'
#' @examples
#' ## Calculate mega environments for TDMaize.
#' gemegaEnv <- gxeMegaEnv(TD = TDMaize, trait = "yld")
#'
#' ## Calculate new mega environments based on the genotypes with the lowest
#' ## value per environment.
#' gemegaEnv2 <- gxeMegaEnv(TD = TDMaize, trait = "yld", method = "min")
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
gxeMegaEnvNw <- function(TD,
                         trials = names(TD),
                         trait,
                         useWinGeno = TRUE,
                         method = c("max", "min"),
                         cutOff = 0.8,
                         byYear = FALSE,
                         sumTab = TRUE) {
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  trials <- chkTrials(trials, TD)
  TDTot <- Reduce(f = rbind, x = TD[trials])
  chkCol(trait, TDTot)
  chkCol("trial", TDTot)
  chkCol("genotype", TDTot)
  if (byYear) {
    chkCol("year", TDTot)
  }
  method <- match.arg(method)
  if (hasName(x = TDTot, name = "megaEnv")) {
    warning("TD already contains a column megaEnv. This column will",
            "be overwritten.\n", call. = FALSE)
  }
  ## Remove genotypes that contain only NAs
  allNA <- by(TDTot, TDTot[["genotype"]], FUN = function(x) {
    all(is.na(x[trait]))
  })
  TDTot <- TDTot[!TDTot[["genotype"]] %in% names(allNA[allNA]), ]
  rmYear <- FALSE
  if (!byYear) {
    TDTot[["year"]] <- 0
    rmYear <- TRUE
  }
  ## Save and then drop factor levels.
  envLevels <- levels(TDTot[["trial"]])[levels(TDTot[["trial"]]) %in% trials]
  TDTot[["trial"]] <- droplevels(TDTot[["trial"]])
  if (useWinGeno) {
    ## Perform AMMI analysis.
    AMMI <- gxeAmmi(TD = createTD(TDTot), trait = trait, nPC = 2, byYear = byYear)
    ## Extract winning genotype per trial.
    winGeno <- by(data = AMMI$fitted, INDICES = AMMI$fitted[["trial"]],
                  FUN = function(trial) {
                    selGeno <- do.call(paste0("which.", method),
                                       args = list(trial[["fittedValue"]]))
                    as.character(trial[["genotype"]])[selGeno]
                  })
    ## Extract values for winning genotype per trial.
    winGenoVal <- by(data = AMMI$fitted, INDICES = AMMI$fitted[["trial"]],
                     FUN = function(trial) {
                       do.call(method, args = list(trial[["fittedValue"]]))
                     })
    ## Create factor based on best genotypes.
    megaFactor <- factor(winGeno,
                         labels = paste0("megaEnv_", seq_along(unique(winGeno))))
    ## Merge factor levels to original data.
    TDTot[["megaEnv"]] <- TDTot[["trial"]]
    levels(TDTot[["megaEnv"]]) <- as.character(megaFactor)
    ## Reapply saved levels to ensure input and output TDTot are identical.
    levels(TDTot[["trial"]]) <- envLevels
    ## If year was added, remove if before creating output.
    if (isTRUE(rmYear)) {
      TDTot <- TDTot[-which(colnames(TDTot) == "year")]
    }
    ## Relevel megaEnv so it is in increasing order.
    TDTot[["megaEnv"]] <- factor(as.character(TDTot[["megaEnv"]]))
    TDOut <- createTD(TDTot)
    ## Create summary table.
    summTab <- data.frame("Mega_factor" = megaFactor,
                          Trial = names(winGeno),
                          "Winning_genotype" = as.character(winGeno),
                          "AMMI_estimates" = as.numeric(winGenoVal))
    summTab <- summTab[order(megaFactor), ]
  } else {
    if (!hasName(x = TDTot, name = "loc")) {
      TDTot[["loc"]] <- TDTot[["trial"]]
    }
    TDTot <- droplevels(TDTot)
    ## Remove locations that appear only in one year.
    locYear <- table(TDTot[["loc"]], TDTot[["year"]])
    # rmLocs <- rownames(locYear)[rowSums(locYear > 0) == 1]
    # TDTot <- droplevels(TDTot[!TDTot[["loc"]] %in% rmLocs, ])
    ## One by one remove locations that don't appear within a year with at
    ## least one other location.
    ## This needs to be done in steps since otherwise in many cases all
    ## locations will be removed.
    continue <- TRUE
    while (continue) {
      locs <- as.character(unique(TDTot[["loc"]]))
      locYear <- table(TDTot[["loc"]], TDTot[["year"]])
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
        TDTot <- droplevels(TDTot[TDTot[["loc"]] != names(which.max(locZeros)), ])
      } else {
        continue <- FALSE
      }
    }
    ## Perform AMMI analysis.
    AMMI <- gxeAmmi(TD = createTD(TDTot), trait = trait, nPC = NULL,
                    byYear = TRUE)
    ammiRaw <- merge(fitted(AMMI), TDTot[c("genotype", "trial", "loc", "year")])
    ## Remove years skipped in AMMI.
    ammiRaw <- droplevels(ammiRaw)
    ## Compute quantile for determining best genotypes per year per location.
    quant <- tapply(X = ammiRaw[["fittedValue"]],
                    INDEX = list(ammiRaw[["year"]], ammiRaw[["loc"]]),
                    FUN = quantile,
                    probs = ifelse(method == "max", cutOff, 1 - cutOff),
                    na.rm = TRUE, type = 1)
    ## Convert quant to data.frame to prevent crash when reshaping.
    quant <- as.data.frame(quant)
    ## Convert quant matrix to long format for merging.
    meltedQuant <- reshape(quant, direction = "long",
                           varying = list(loc = colnames(quant)),
                           ids = rownames(quant), idvar = "year",
                           times = colnames(quant), timevar = "loc",
                           v.names = "qntValue")
    ## Merge quantile to data.
    ammiQnt <- merge(ammiRaw, meltedQuant)
    ## Replace values by 0 when lower than quantile value and 1 otherwise.
    ## Reverse the equation if low values are better for current trait.
    ammiQnt[["fittedValue"]] <- if (method == "max") {
      as.numeric(ammiQnt[["fittedValue"]] > ammiQnt[["qntValue"]])
    } else {
      as.numeric(ammiQnt[["fittedValue"]] < ammiQnt[["qntValue"]])
    }
    ## Reshape data to get locations as header.
    ammiLocTab <- tapply(ammiQnt$fittedValue,
                         list(ammiQnt$genotype, ammiQnt$year, ammiQnt$loc),
                         FUN = I)
    ammiLocMat <- matrix(ammiLocTab, ncol = dim(ammiLocTab)[3],
                         dimnames = list(NULL, dimnames(ammiLocTab)[[3]]))
    nonMissRows <- apply(X = ammiLocMat, MARGIN = 1, FUN = function(x) {
      any(!is.na(x))
    })
    ammiLoc <- cbind(expand.grid(dimnames(ammiLocTab)[1:2]), ammiLocMat)
    ammiLoc <- ammiLoc[nonMissRows, ]
    colnames(ammiLoc)[1:2] <- c("genotype", "year")
    ## Compute means and standard deviations.
    ## For sd division by n should be used so this cannot be done by sd().
    Xi = tapply(X = ammiQnt[["fittedValue"]],
                INDEX = list(ammiQnt[["year"]], ammiQnt[["loc"]]), FUN = mean)
    SXi = tapply(X = ammiQnt[["fittedValue"]],
                 INDEX = list(ammiQnt[["year"]], ammiQnt[["loc"]]),
                 FUN = function(x) {
                   sqrt(sum((x - mean(x)) ^ 2) / length(x))
                 })
    ## Compute correlations between locations per year.
    r0 <- by(data = ammiLoc[, c(3:ncol(ammiLoc))], INDICES = ammiLoc[["year"]],
             FUN = cor, use = "pairwise.complete.obs")
    ## Form combinations of locations.
    locs <- as.character(unique(TDTot[["loc"]]))
    combs <- combn(locs, m = 2)
    ## Compute correlations across years per genotype.
    ## If numbers of observations are similar per location x year this
    ## has very similar results to just applying cor with pairwise.complete.obs
    ## on the complete dataset. However when the numbers start differing so
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
    ## Throw away tempfile.
    on.exit(unlink(tmp), add = TRUE)
    ## Loop over number of clusters to determine best number of clusters
    ## by minimizing ratio CR/DR.
    for (k in 2:ceiling(length(locs) / 2)) {
      ## Extract cluster groups for k clusters.
      clustGr <- data.frame(megaEnv = factor(cutree(tree, k = k)))
      ## Merge cluster groups to data.
      modDat <- merge(TDTot, clustGr, by.x = "loc", by.y = "row.names")
      ## Model with regions, one varcomp.
      sink(file = tmp)
      ## Fit model with regions.
      modReg <- asreml::asreml(
        fixed = formula(paste(trait, "~ 1 + year * loc")),
        random = ~ genotype + genotype:megaEnv + genotype:loc + genotype:year +
          genotype:megaEnv:year,
        data = modDat, maxiter = 50, trace = FALSE)
      sink()
      ## Extract variance components.
      vcReg <- summary(modReg)$varcomp$component
      ## Compute mdian number of locations, years and megaEnv per genotype.
      nR <- median(rowSums(table(modDat[["genotype"]], modDat[["megaEnv"]]) > 0))
      nL <- median(rowSums(table(modDat[["genotype"]], modDat[["loc"]]) > 0))
      nY <- median(rowSums(table(modDat[["genotype"]], modDat[["year"]]) > 0))
      ## Compute H2 over locations.
      H2Loc <- vcReg[1] / (vcReg[1] + (vcReg[2] / nR) + (vcReg[3] / nY) +
                             (vcReg[4] / (nR * nY)) + (vcReg[5] / (nR * nL)) +
                             (vcReg[6] / (nR * nL * nY)))
      ## Compute H2 over regions.
      H2Reg <- (vcReg[1] + vcReg[2]) / (vcReg[1] + vcReg[2] + (vcReg[3] / nY) +
                                          (vcReg[4] / nY) + (vcReg[5] / nL) +
                                          (vcReg[6] / (nL * nY)))
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
    ## If year was added, remove if before creating output.
    if (isTRUE(rmYear)) {
      clustRes <- clustRes[-which(colnames(clustRes) == "year")]
    }
    ## Create TD Output.
    TDOut <- createTD(clustRes)
    ## Attach cluster groups as attribute.
    clustGrRes <- clustGrRes[order(clustGrRes[["megaEnv"]]), , drop = FALSE]
    clustGrRes <- data.frame("Mega factor" = clustGrRes[["megaEnv"]],
                             Trial = rownames(clustGrRes))
    # attr(TDOut, "sumTab") <- clustGrRes
    # attr(TDOut, "CRDR") <- CRDRMin

    summTab <- clustGrRes
    attr(summTab, "CRDR") <- CRDRMin
  }
  return(createMegaEnv(TD = TDOut, summTab = summTab, trait = trait))
}

#' Helper function for computing combined correlations for a pair of variables.
#' Following the notation of Charter (1993).
#'
#' @noRd
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
#' @noRd
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
  nl1l2 <- table(ammi[!is.na(ammi[, l1]) & !is.na(ammi[, l2]), ][["year"]])
  ## Extract correlations for l1 and l2 per year for list of correlations.
  rl1l2 <- sapply(X = r0, FUN = "[", l1, l2)
  ## Determine which observations should be included in calculation of r.
  inclObs <- !is.null(rl1l2) & !is.na(rl1l2)
  ## Compute combined correlation.
  combCor(Xi = Xi[, l1][inclObs], Yi = Xi[, l2][inclObs],
          SXi = SXi[, l1][inclObs], SYi = SXi[, l2][inclObs],
          ni = nl1l2[inclObs], ri = rl1l2[inclObs])
}










