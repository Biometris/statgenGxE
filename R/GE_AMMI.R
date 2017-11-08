#' AMMI analysis
#'
#' This function fits a model which involves the Additive Main effects (i.e. genotype and
#' environment) along with the Multiplicative Interaction effects of principal
#' component analysis (PCA).
#'
#' @param Y A data frame object.
#' @param trait A character string specifying a trait column of the data.
#' @param genotype A character string specifying a genotype column of the data.
#' @param env A character string specifying an environment column of the data.
#' @param nPC An integer specifying the number of principal components used.
#' as the multiplicative term of genotype-by-environment. \code{nPC=2} by default.
#' @param center  A logical value indicating whether the variables
#'  should be shifted to be zero centered.
#' @param scale A logical value indicating whether the variables should
#'  be scaled to have unit variance before the analysis takes
#'  place. The default is \code{FALSE}.
#' @param AMMI2plot A logical value specifying whether an AMMI2 biplot is drawn.
#' @param scaleAMMI2 The variables are scaled by \code{lambda ^ scale} and the observations
#' are scaled by \code{lambda ^ (1-scale)}
#' where \code{lambda} are the singular values as computed by \code{\link[stats]{princomp}}.
#' Normally \code{0 <= scale <= 1}, and a warning will be issued if the specified scale is
#' outside this range.
#' @param AMMI1plot A logical value determining whether the AMMI1 biplot (genotypes and
#' environments means vs PC1) is drawn.
#' @param scaleAMMI1 as same as described in \code{scaleAMMI2}.
#' @return A list of three objects, a data frame object of environment scores, a data frame
#' object of genotype scores and an object of class \code{\link[stats]{anova}}, and a
#' matrix of the fitted values from the AMMI model.
#'
#' @examples
#' mydat <- GE.read.csv(system.file("extdata", "F2maize_pheno.csv", package = "RAP"),
#'                      env ="env!", genotype ="genotype!", trait = "yld")
#' names(mydat) <- c("env", "genotype","yld")
#' GE.AMMI(Y = mydat, trait = "yld", genotype = "genotype", env = "env", nPC = 2,
#'         center = TRUE, scale = FALSE, AMMI2plot = TRUE, scaleAMMI2 = 1)
#'
#' @import stats graphics grDevices
#' @export

GE.AMMI <- function(Y,
                    trait,
                    genotype,
                    env,
                    nPC = 2,
                    center = TRUE,
                    scale = FALSE,
                    AMMI2plot = TRUE,
                    scaleAMMI2 = 1,
                    AMMI1plot = FALSE,
                    scaleAMMI1 = 1) {
  #drop factor levels
  Y[[genotype]] <- droplevels(Y[[genotype]])
  Y[[env]] <- droplevels(Y[[env]])
  nGeno <- nlevels(Y[[genotype]])
  nEnv  <- nlevels(Y[[env]])
  nTrait <- nrow(Y)
  # requre number of environments >=3
  if (nEnv < 3) {
    stop("Requires number of environments greater and equal than 3 for running the AMMI model.\n")
  }
  #check if the supplied data contains the genotype by environment means
  if (nTrait != nGeno * nEnv) {
    stop("Only allows the genotype by environment means, \ni.e., one trait value per
         genotype per enviroment.\n")
  }
  if (any(is.na(Y[[trait]]))) {
    y0 <- tapply(X = Y[[trait]], INDEX = Y[, c(genotype, env)], FUN = identity)
    yIndex <- tapply(X = 1:nTrait, INDEX = Y[, c(genotype, env)], FUN = identity)
    na_yes_no <- is.na(y0)
    # imputaion
    y1 <- RAP.multmissing(y0, maxcycle = 10, na.strings = NA)
    replaceVal <- y1[na_yes_no]
    Y[yIndex[na_yes_no], trait] <- replaceVal
  }
  # Descriptive statistics
  envMean <- tapply(X = Y[[trait]], INDEX = Y[[env]], FUN = mean)
  genoMean <- tapply(X =Y[[trait]], INDEX = Y[[genotype]], FUN = mean)
  overallMean <- mean(Y[[trait]])
  # Fit the linear model
  model <- lm(as.formula(paste(trait, "~", genotype, "+", env)), data = Y)
  # calculate residuals & fitted values of the linear model
  X <- tapply(X = resid(model), INDEX = Y[, c(genotype, env)], FUN = identity)
  fittedVals <- tapply(X = fitted(model), INDEX = Y[, c(genotype, env)], FUN = identity)
  X <- as.matrix(X)
  gNames <- rownames(X)
  if (is.null(gNames)) {
    rownames(X) <- rownames(fittedVals) <- gNames
  }
  # Use R in-built prcomp
  pca <- prcomp(x = X, retx = TRUE, center = center, scale. = scale)
  cump2 <- summary(pca)$importance[3, 2]
  propPc1 <- summary(pca)$importance[2, 1]
  propPc2 <- summary(pca)$importance[2, 2]
  loadings <- pca$rotation
  scores <- pca$x
  if (AMMI2plot){
    if (scaleAMMI2 == 1) {
      info <- "environment scaling"
    } else if (scaleAMMI2 == 0) {
      info <- "genotype scaling"
    } else if (scaleAMMI2 == 0.5) {
      info <- "symmetric scaling"
    } else {
      info <- paste0(round(cump2 * 100, 1), "%")
    }
    oldPar <- par(xpd = NA)
    biplot(pca, scale = scaleAMMI2, col = c("orange3", "navyblue"),
           main = paste0("AMMI2 biplot for ", trait, " (", info, ")"),
           xlab = paste0("PC1 (", round(propPc1 * 100, 1), "%)"),
           ylab = paste0("PC2 (", round(propPc2 * 100, 1), "%)"))
    par(oldPar)
  }
  if (AMMI1plot){
    # calculate lambda scale
    lam <- pca$sdev[1]
    if (is.null(n <- pca$n.obs)) {
      n <- 1
    }
    lam <- lam * sqrt(n)
    if (scaleAMMI1 < 0 || scaleAMMI1 > 1) {
      warning("'scale' is outside [0, 1]")
    }
    if (scaleAMMI1 != 0) {
      lam <- lam ^ scaleAMMI1
    } else {
      lam <- 1
    }
    dev.new()
    oldPar <- par(xpd = NA)
    plot(x = 1, type = 'n', xlim = range(c(envMean, genoMean)),
         ylim = range(c(loadings[, 1] * lam, scores[, 1] / lam)),
         xlab = "Main Effects", ylab = paste0("PC1 (", round(propPc1 * 100, 1), "%)"),
         main = paste0("AMMI1 biplot for ", trait))
    points(x = envMean, y = loadings[, 1] * lam, type = "n", col = "navyblue", lwd = 5)
    text(x = envMean, y = loadings[, 1] * lam, labels = row.names(envMean),
         adj = c(0.5, 0.5), col = "navyblue")
    points(x = genoMean, y = scores[, 1] / lam, type = "n", col = "orange3", lwd = 5)
    text(x = genoMean, y = scores[, 1] / lam, labels = row.names(genoMean),
         adj = c(0.5, 0.5), col = "orange3")
    abline(h = 0, v = overallMean, lty = 5)
    par(oldPar)
  }
  # calculating the AMMI-estimates per genotype per environment
  mterms <- matrix(data = 0, nrow = nGeno, ncol = nEnv)
  for (ii in 1:nPC) {
    mterms <- mterms + outer(scores[, ii], loadings[, ii])
  }
  fitted <- fittedVals + mterms
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
  nGeno <- nlevels(Y[, genotype])
  nEnv  <- nlevels(Y[, env])
  dfPC <- nGeno + nEnv - 3 - (2 * (1:nPC - 1))
  dfresid <- a1["Interactions", "Df"] - sum(dfPC)
  addTbl[,"Df"] <- c(dfPC, dfresid)
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
  return(list(environmentScores = pca$rotation, genotypeScores = pca$x,
              ANOVA = a0, fitted = fitted))
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
