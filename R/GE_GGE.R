#' GGE biplot
#'
#' This function draws GGE biplots which are useful for assessing the performance of genotypes
#' in different environments.
#'
#' @param Y A data frame object.
#' @param trait A character string specifying a trait column of the data.
#' @param genotype A character string specifying a genotype column of the data.
#' @param env A character string specifying an environment column of the data.
#' @param nPC An integer specifying the number of principal components used
#' as the multiplicative term of genotype by environment. \code{nPC=2} by default.
#' @param center  A logical value indicating whether the variables
#'  should be shifted to be zero centered.
#' @param scale A logical value indicating whether the variables should
#'  be scaled to have unit variance before the analysis takes
#'  place. The default is \code{FALSE}.
#' @param GGE2plot A logical value determining whether the GGE biplot is drawn.
#' @param scaleBiplot The variables are scaled by \code{lambda ^ scale} and the
#' observations are scaled by \code{lambda ^ (1-scale)}
#' where \code{lambda} are the singular values as computed by \code{\link[stats]{princomp}}.
#' Normally \code{0 <= scale <= 1},
#' and a warning will be issued if the specified scale is outside this range.
#' @return A list of two objects, a list with class \code{\link[stats]{prcomp}} and an object of
#' class \code{\link[stats]{anova}}.
#'
#' @examples
#' mydat <- GE.read.csv(system.file("extdata", "F2maize_pheno.csv", package = "RAP"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' names(mydat)=c("env", "genotype","yld")
#' GE.GGE(Y=mydat, trait="yld", genotype="genotype", env="env", scaleBiplot=1)
#'
#' @export
GE.GGE <- function(Y,
                   trait,
                   genotype,
                   env,
                   nPC = 2,
                   center = TRUE,
                   scale = FALSE,
                   GGE2plot = TRUE,
                   scaleBiplot = 0) {
  #drop factor levels
  Y[[genotype]] <- droplevels(Y[[genotype]])
  Y[[env]] <- droplevels(Y[[env]])
  nGeno <- nlevels(Y[[genotype]])
  nEnv  <- nlevels(Y[[env]])
  nTrait <- nrow(Y)
  #check if the supplied data contains the genotype by environment means
  if (nTrait != nGeno * nEnv) {
    stop("Only allows the genotype by environment means, \ni.e., one trait value per
         genotype per enviroment")
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
  # Fit the linear model
  model <- lm(as.formula(paste(trait,"~", env)), data = Y)
  # Use R in-built prcomp
  X <- tapply(X = resid(model), INDEX = Y[, c(genotype, env)], FUN = identity)
  X <- as.matrix(X)
  pca <- prcomp(x = X, center = center, scale. = scale)
  if (GGE2plot) {
    cump2 <- summary(pca)$importance[3, 2]
    propPC1 <- summary(pca)$importance[2, 1]
    propPC2 <- summary(pca)$importance[2, 2]
    if (scaleBiplot == 1) {
      info <- "environment scaling"
    } else if (scaleBiplot == 0) {
      info <- "genotype scaling"
    } else if (scaleBiplot==0.5) {
      info <- "symmetric scaling"
    } else {
      info <- paste0(round(cump2 * 100, 1), "%")
    }
    oldPar <- par(xpd = NA)
    biplot(pca, scale = scaleBiplot, col = c("orange3", "navyblue"),
      main = paste0("GGE2 biplot for ", trait, " (", info, ")"),
      xlab = paste0("PC1 (", round(propPC1 * 100, 1), "%)"),
      ylab = paste0("PC2 (", round(propPC2 * 100, 1), "%)"))
    par(oldPar)
  }
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
  nEnv <- nlevels(Y[, env])
  dfPC <- nGeno + nEnv - 3 - (2 * (1:nPC - 1))
  dfResid <- a1["Interactions", "Df"] - sum(dfPC)
  addTbl[,"Df"] <- c(dfPC, dfResid)
  # Add the sum of squaures for PC scores and residuals
  PCAVar <- pca$sdev ^ 2
  totalVar <- sum(PCAVar)
  propVar <- PCAVar / totalVar
  ssPC <- a1["Interactions", "Sum Sq"] * propVar[1:nPC]
  ssResid <- a1["Interactions", "Sum Sq"] - sum(ssPC)
  addTbl[, "Sum Sq"] <- c(ssPC, ssResid)
  # Add the mean squaures for PC scores and residuals
  addTbl[, "Mean Sq"] <- addTbl[, "Sum Sq"] / addTbl[, "Df"]
  whichinf <- is.infinite(addTbl[, "Mean Sq"])
  addTbl[, "Mean Sq"][whichinf] <- NA
  # Add the F-values for PC scores
  addTbl[-(nPC + 1), "F value"] <- addTbl[-(nPC + 1), "Mean Sq"] /
    addTbl["Residuals", "Mean Sq"]
  # Add the p-value for PC scores
  addTbl[-(nPC + 1),"Pr(>F)"] <- 1 - sapply(X = 1:nPC, FUN = function(i) {
    pf(q = addTbl[i, "F value"], df1 = addTbl[i, "Df"],
       df2 = addTbl["Residuals", "Df"])
    })
  # ANOVA table for AMMI model
  a0 <- rbind(a1, as.data.frame(addTbl))
  return(list(PCA = pca, ANOVA = a0))
}
