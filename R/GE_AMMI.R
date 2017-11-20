#' AMMI analysis
#'
#' This function fits a model which involves the Additive Main effects (i.e. genotype and
#' environment) along with the Multiplicative Interaction effects of principal
#' component analysis (PCA).
#'
#' @param TD an object of class TD.
#' @param trait A character string specifying a trait column of the data.
#' @param nPC An integer specifying the number of principal components used.
#' as the multiplicative term of genotype-by-environment. \code{nPC=2} by default.
#' @param center  A logical value indicating whether the variables
#'  should be shifted to be zero centered.
#' @param scale A logical value indicating whether the variables should
#'  be scaled to have unit variance before the analysis takes
#'  place. The default is \code{FALSE}.
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
#' myDat <- GE.read.csv(system.file("extdata", "F2maize_pheno.csv", package = "RAP"),
#'                      env = "env!", genotype = "genotype!", trait = "yld")
#' myTD <- createTD(data = myDat, genotype = "genotype!", env = "env!")
#' GE.AMMI(TD = myTD, trait = "yld", nPC = 2, center = TRUE, scale = FALSE)
#'
#' @import stats
#' @export

GE.AMMI <- function(TD,
                    trait,
                    nPC = 2,
                    center = TRUE,
                    scale = FALSE) {
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (!"env" %in% colnames(TD)) {
   stop("TD should contain a column env to be able to run an AMMI analysis.\n")
  }
  #drop factor levels
  TD$genotype <- droplevels(TD$genotype)
  TD$env <- droplevels(TD$env)
  nGeno <- nlevels(TD$genotype)
  nEnv  <- nlevels(TD$env)
  nTrait <- nrow(TD)
  # requre number of environments >=3
  if (nEnv < 3) {
    stop("Requires number of environments greater and equal than 3 for running the AMMI model.\n")
  }
  #check if the supplied data contains the genotype by environment means
  if (nTrait != nGeno * nEnv) {
    stop("Only allows the genotype by environment means, \ni.e., one trait value per
         genotype per enviroment.\n")
  }
  if (any(is.na(TD[[trait]]))) {
    y0 <- tapply(X = TD[[trait]], INDEX = TD[, c("genotype", "env")], FUN = identity)
    yIndex <- tapply(X = 1:nTrait, INDEX = TD[, c("genotype", "env")], FUN = identity)
    na_yes_no <- is.na(y0)
    # imputation
    y1 <- RAP.multmissing(y0, maxcycle = 10, na.strings = NA)
    replaceVal <- y1[na_yes_no]
    TD[yIndex[na_yes_no], trait] <- replaceVal
  }
  # Descriptive statistics
  envMean <- tapply(X = TD[[trait]], INDEX = TD$env, FUN = mean)
  genoMean <- tapply(X = TD[[trait]], INDEX = TD$genotype, FUN = mean)
  overallMean <- mean(TD[[trait]])
  # Fit the linear model
  model <- lm(as.formula(paste(trait, "~ genotype + env")), data = TD)
  # calculate residuals & fitted values of the linear model
  X <- tapply(X = resid(model), INDEX = TD[, c("genotype", "env")], FUN = identity)
  fittedVals <- tapply(X = fitted(model), INDEX = TD[, c("genotype", "env")],
                       FUN = identity)
  X <- as.matrix(X)
  gNames <- rownames(X)
  if (is.null(gNames)) {
    rownames(X) <- rownames(fittedVals) <- gNames
  }
  # Use R in-built prcomp
  pca <- prcomp(x = X, retx = TRUE, center = center, scale. = scale)
  loadings <- pca$rotation
  scores <- pca$x
  # calculating the AMMI-estimates per genotype per environment
  mTerms <- matrix(data = 0, nrow = nGeno, ncol = nEnv)
  for (ii in 1:nPC) {
    mTerms <- mTerms + outer(scores[, ii], loadings[, ii])
  }
  fitted <- fittedVals + mTerms
  # ANOVA table for linear model
  a1 <- anova(model)
  tNames <- rownames(a1)
  rownames(a1)[which(tNames == "Residuals")] <- "Interactions"
  # Extend the existing ANOVA table
  addTbl <- matrix(NA, nPC + 1, 5)
  colnames(addTbl) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
  tNames <- paste0("PC", 1:nPC)
  rownames(addTbl) <- c(tNames, "Residuals")
  # Add the df for PC scores and residuals
  nGeno <- nlevels(TD$genotype)
  nEnv  <- nlevels(TD$env)
  dfPC <- nGeno + nEnv - 3 - (2 * (1:nPC - 1))
  dfResid <- a1["Interactions", "Df"] - sum(dfPC)
  addTbl[,"Df"] <- c(dfPC, dfResid)
  # Add the sum of squares for PC scores and residuals
  PCAVar <- pca$sdev ^ 2
  totalVar <- sum(PCAVar)
  propVar <- PCAVar / totalVar
  ssPC <- a1["Interactions", "Sum Sq"] * propVar[1:nPC]
  ssresid <- a1["Interactions", "Sum Sq"] - sum(ssPC)
  addTbl[,"Sum Sq"] <- c(ssPC, ssresid)
  # Add the mean squares for PC scores and residuals
  addTbl[,"Mean Sq"] <- addTbl[, "Sum Sq"] / addTbl[, "Df"]
  whichinf <- is.infinite(addTbl[, "Mean Sq"])
  addTbl[, "Mean Sq"][whichinf] <- NA
  # Add the F-values for PC scores
  addTbl[-(nPC + 1), "F value"] <- addTbl[-(nPC + 1), "Mean Sq"] / addTbl["Residuals", "Mean Sq"]
  # Add the p-value for PC scores
  addTbl[-(nPC + 1), "Pr(>F)"] <- 1 - sapply(X = 1:nPC, FUN = function(i) {
    pf(q = addTbl[i, "F value"], df1 = addTbl[i, "Df"], df2 = addTbl["Residuals", "Df"])
  })
  # ANOVA table for AMMI model
  a0 <- rbind(a1, as.data.frame(addTbl))
  return(createAMMI(envScores = pca$rotation, genoScores = pca$x,
                    importance = as.data.frame(summary(pca)$importance),
                    anova = a0, fitted = fitted,
                    trait = trait, envMean = envMean, genoMean = genoMean,
                    overallMean = overallMean))
}

## explicitly use singular value decomposition to compute pca
#  s <- svd(X)
#  U <-s$u
#  L <- s$d
#  V <- s$v
#  LL <- sqrt(diag(L))
#  loadings <- V %*% LL
#  scores <- U %*% LL
#  res <- new.env()
#  res$loadings <- loadings
#  res$scores <- scores
#  as.list(res)
