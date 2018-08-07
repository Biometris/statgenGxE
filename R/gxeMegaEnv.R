#' Form mega-environments based on winning genotypes from an AMMI model
#'
#' This function fits an AMMI model and then using the fitted values produces
#' a new factor based on the winning genotype in each environment. This factor
#' is added as a column megaEnv to the input data. If a column megaEnv already
#' exists the existing data is overwritten with a warning.
#'
#' @inheritParams gxeAmmi
#'
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
#' @importFrom methods getFunction
#' @export
gxeMegaEnv <- function(TD,
                       trials = names(TD),
                       trait,
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
  ## Save and then drop factor levels.
  envLevels <- levels(TDTot$trial)
  TDTot$trial <- droplevels(TDTot$trial)
  rmYear <- FALSE
  if (!hasName(x = TDTot, name = "year")) {
    TDTot$year <- 0
    rmYear <- TRUE
  }
  if (!hasName(x = TDTot, name = "loc")) {
    TDTot$loc <- TDTot$trial
  }
  ## Perform AMMI analysis.
  AMMI <- gxeAmmi(TD = createTD(TDTot), trait = trait, nPC = 2, byYear = TRUE)
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
  locs <- unique(ammiQnt$loc)
  combs <- combn(locs, m = 2)
  ## Compute correlations across years per genotype.
  combs <- rbind(combs, mapply(FUN = combLocs, combs[1, ], combs[2, ],
                        MoreArgs = list(ammi = ammiLoc, r0 = r0, Xi = Xi,
                                        SXi = SXi)))
  ## Put computed correlations in lower half of correlation matrix.
  corMat <- matrix(nrow = length(locs), ncol = length(locs),
                   dimnames = list(locs, locs))
  corMat[lower.tri(corMat)] <- combs[3, ]
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
        genotype:megaEnv:year, data = modDat, maxiter = 200)
    sink()
    ## Extract variance components.
    vcReg <- summary(modReg)$varcomp$component
    ## Compute number of locations.
    nL <- length(unique(modDat$loc))
    ## Compute H2 over locations.
    H2Tr <- vcReg[1] / (vcReg[1] + vcReg[2] / k + vcReg[3] / (k * nL))
    ## Compute H2 over regions.
    H2Reg <- (vcReg[1] + vcReg[2]) / (vcReg[1] + vcReg[2] + vcReg[3] / (k * nL))
    ## Compute genetic correlation.
    rho <- vcReg[1] / sqrt(vcReg[1] * (vcReg[1] + vcReg[2]))
    ## Compute ratio CD/CR.
    CRDR <- rho * sqrt(H2Tr / H2Reg)
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
  if (sumTab) {
    printCoefmat(clustGrRes)
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
#' @keywords internal
#'
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

