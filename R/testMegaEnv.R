rawDat <- read.csv('D:/R other projects/testScripts/WheatYieldUK.csv',
                   stringsAsFactors = FALSE)
ammiRaw <- read.csv('D:/R other projects/testScripts/AMMIprediction_UK.csv',
                    quote = "\'", stringsAsFactors = FALSE)

## Helper function for computing combined correlations for a pair of variables.
## Following the notation of Charter (1993).
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

## Helper function for computing the combined correlation for 2 locations.
combLocs <- function(l1,
                     l2,
                     ammi,
                     Xi,
                     SXi) {
  ## Compute number of observations for combination of l1 and l2 per year.
  nl1l2 <- table(ammi[!is.na(ammi[, l1]) & !is.na(ammi[, l2]), ]$year)
  ## Extract correlations for l1 and l2 per year for list of correlations.
  rl1l2 <- sapply(X = r0, FUN = "[", l1, l2)
  ## Determine which observations should be included in calculation of r.
  inclObs <- !is.na(rl1l2)
  ## Compute combined correlation.
  combCor(Xi = Xi[, l1][inclObs], Yi = Xi[, l2][inclObs],
          SXi = SXi[, l1][inclObs], SYi = SXi[, l2][inclObs], ni = nl1l2,
          ri = rl1l2[inclObs])
}

## Compute quantile for determining best genotypes per year per location.
quant <- tapply(X = ammiRaw$ammi_pred, INDEX = list(ammiRaw$year, ammiRaw$loc),
                FUN = quantile, probs = 0.8, na.rm = TRUE, type = 1)
## Merge quantile to data.
ammiQnt <- merge(ammiRaw, reshape2::melt(quant), by.x = c("year", "loc"),
                 by.y = c("Var1", "Var2"))
## Replace values by 0 when lower than quantile value and 1 otherwise.
ammiQnt$ammi_pred <- as.numeric(ammiQnt$ammi_pred > ammiQnt$value)
## Reshape data to get locations as header.
ammiLoc <- reshape2::dcast(ammiQnt, geno + year ~ loc, value.var = "ammi_pred",
                           drop = FALSE)
## Compute means and standard deviations.
## For sd division by n should be used so this cannot be done by sd().
Xi = tapply(ammiQnt$ammi_pred, list(ammiQnt$year, ammiQnt$loc), FUN = mean)
SXi = tapply(ammiQnt$ammi_pred, list(ammiQnt$year, ammiQnt$loc),
             FUN = function(x) {
               sqrt(sum((x - mean(x)) ^ 2) / length(x))
             })
## Compute correlations between environments per year.
r0 <- by(data = ammiLoc[, c(3:ncol(ammiLoc))], INDICES = ammiLoc$year, FUN = cor,
         use = "pairwise.complete.obs")
## Form combinations of locations
locs <- unique(ammiQnt$loc)
combs <- combn(locs, m = 2)
## Compute correlations accross years per genotype.
combs <- rbind(combs, mapply(FUN = combLocs, combs[1, ], combs[2, ],
                             MoreArgs = list(ammi = ammiLoc, Xi = Xi,
                                             SXi = SXi)))
## Put computed correlations in lower half of correlation matrix.
corMat <- matrix(nrow = length(locs), ncol = length(locs),
                 dimnames = list(locs, locs))
corMat[lower.tri(corMat)] <- as.numeric(combs[3, ])
## Compute distances.
distMat <- as.dist(sqrt(1 - corMat ^ 2))
## Cluster locations.
tree <- hclust(distMat, method = "ward.D")
clustRes <- data.frame()
CDCRMin <- Inf
## Create tempfile for diverting asreml output
tmp <- tempfile()
## Loop over number of clusters to determine best number of clusters
## by minimizing ratio CD/CR.
for (k in 2:ceiling(length(locs) / 2)) {
  ## Extract cluster groups for k clusters.
  clustGr <- data.frame(megaEnv = cutree(tree, k = k))
  ## Merge cluster groups and raw data to ammi data.
  ammiTot <- merge(ammiQnt, clustGr, by.x = "loc", by.y = "row.names")
  modDat <- merge(ammiTot, rawDat)

  ## Model without regions, one varcomp.
  # modLoc <- asreml::asreml(fixed = yield ~ 1 + year * loc,
  #                          random = ~ geno + geno:loc + geno:year,
  #                          data = modDat, maxiter = 50)
  # vcLoc <- summary(modLoc)$varcomp$component

  ## Model with regions, one varcomp.
  sink(file = tmp)
  ## Fit model with regions.
  modReg <- asreml::asreml(fixed = yield ~ 1 + year * loc,
                           random = ~ geno + geno:megaEnv + geno:loc + geno:year,
                           data = modDat, maxiter = 200)
  sink()
  ## Extract variance components.
  vcReg <- summary(modReg)$varcomp$component
  ## Compute number of locations.
  nl <- length(unique(modDat$loc))
  ## Compute H2 over locations.
  H2Loc <- vcReg[1] / (vcReg[1] + vcReg[2] / k + vcReg[3] / (k * nl))
  ## Compute H2 over regions.
  H2Reg <- (vcReg[1] + vcReg[2]) / (vcReg[1] + vcReg[2] + vcReg[3] / (k * nl))
  ## Compute genetic correlation.
  rho <- vcReg[1] / sqrt(vcReg[1] * (vcReg[1] + vcReg[2]))
  ## Compute ratio CD/CR.
  CDCR <- rho * sqrt(H2Loc / H2Reg)
  if (CDCR < CDCRMin) {
    ## If CDCR is smaller than the previous minimum readjust values.
    clustRes <- ammiTot
    CDCRMin <- CDCR
  }
}
## Throw away tempfile.
unlink(tmp)
