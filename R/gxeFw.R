#' Modified joint regression analysis
#'
#' This function performs a modified joint analysis of data classified by two factors.
#'
#' @inheritParams gxeAmmi
#'
#' @param maxCycle An integer specifying the maximum number of iterations to be achieved.
#' By default, \code{maxCycle = 15}.
#' @param tol A small positive numerical value specifying convergence tolerance.
#' By default, \code{tol = 0.001}.
#' @param sorted A character string specifying the sorting order of the results.
#'
#' @return An object of class \code{\link{FW}}, a list containing
#' \item{estimates}{a data.frame containing the estimated values.}
#' \item{anova}{a data.frame containing anova scores of the FW analysis.}
#' \item{envEffs}{a data.frame containing the environmental effects.}
#' \item{data}{the data.frame on which the analysis was performed.}
#' \item{fittedGeno}{a numerical vector containing the fitted values for the genotypes.}
#' \item{trait}{a character vector containing the analyzed trait.}
#' \item{nGeno}{a numerical value containing the number of genotypes in the analysis.}
#' \item{nEnv}{a numerical value containing the number of environments in the analysis.}
#' \item{tol}{a numerical value containing the tolerance used during the analysis.}
#' \item{iter}{a numberical value containing the number of iterations for the
#' analysis to converge.}
#'
#' @references Finlay, K.W. & Wilkinson, G.N. (1963). The analysis of adaptation
#' in a plant-breeding programme. Australian Journal of Agricultural
#' Research, 14, 742-754.
#'
#' @examples
#' ## Run Finlay-Wilkinson analysis
#' geFW <- gxeFw(TDMaize, trait = "yld")
#' report(geFW, outfile = "./testReports/reportFW.pdf")
#'
#' @export
gxeFw <- function(TD,
                  trait,
                  maxCycle = 15,
                  tol = 0.001,
                  sorted = c("ascending", "descending", "none")) {
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (!"env" %in% colnames(TD)) {
    stop("TD should contain a column env to be able to run an Finlay
         Wilkinson analysis.\n")
  }
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TD)) {
    stop("trait has to be a column in TD.\n")
  }
  if (is.null(maxCycle) || !is.numeric(maxCycle) || length(maxCycle) > 1 ||
      maxCycle != round(maxCycle) || maxCycle < 1) {
    stop("maxCycle should be a positive integer.\n")
  }
  if (is.null(tol) || !is.numeric(tol) || length(tol) > 1 || tol < 0) {
    stop("tol should be a numerical value > 10^-6.\n")
  }
  sorted <- match.arg(sorted)
  nGeno <- nlevels(TD$genotype)
  nEnv <- nlevels(TD$env)
  ## Save Sum Sq & Df.
  rDev <- rDf <- rep(NA, 5)
  ## Estimate env effects with the sensitivity beta = 1.
  model0 <- lm(as.formula(paste(trait, "~-1 + env + genotype")),
               data = TD, na.action = na.exclude)
  aov0 <- anova(model0)
  tPos <- rownames(aov0) == "Residuals"
  rDev[2] <- aov0[["Sum Sq"]][tPos]
  rDf[2] <- aov0[["Df"]][tPos]
  coeffsModel0 <- coefficients(model0)
  envEffs0 <- coeffsModel0[grep(pattern = "env", x = names(coeffsModel0))]
  ## Adjust env effects to have mean zero.
  envEffs0 <- envEffs0 - mean(envEffs0, na.rm = TRUE)
  envEffs <- rep(NA, nrow(TD))
  for (ii in names(envEffs0)) {
    tPosition <- sapply(X = TD$env, FUN = function(x) {
      grepl(pattern = x, x = ii)
    })
    envEffs[tPosition] <- envEffs0[ii]
  }
  TD <- cbind(TD, envEffs)
  ## Set initial values for sensitivity beta.
  TD$beta <- 1
  ## Set a relative difference to be large.
  maxDiff <- Inf
  ## Set iteration to be 1.
  iter <- 1
  ## Iterate to fit 'y(i,j) = genMean(i)+beta(i)*envEffs(j)'
  while (maxDiff > tol && iter <= maxCycle) {
    beta0 <- TD$beta
    ## Fit model with current genotype sensitivity relevant to each unit.
    model1 <- lm(as.formula(paste(trait, "~-1 + genotype + genotype:envEffs")),
                 data = TD, na.action = na.exclude)
    coeffsModel1 <- coefficients(model1)
    ## Update beta.
    TD$beta <- coeffsModel1[match(paste0("genotype", TD$genotype, ":envEffs"),
                                  names(coeffsModel1))]
    TD$beta <- TD$beta / mean(TD$beta, na.rm = TRUE)
    ## Fit model with current env means relevant to each unit.
    model2 <- lm(as.formula(paste(trait, "~-1 + genotype + env:beta")),
                 data = TD, na.action = na.exclude)
    coeffsModel2 <- coefficients(model2)
    ## Update envEffs.
    TD$envEffs <- coeffsModel2[match(paste0("env", TD$env, ":beta"),
                                     names(coeffsModel2))]
    TD[is.na(TD$envEffs), "envEffs"] <- 0
    TD$envEffs <- TD$envEffs - mean(TD$envEffs)
    ## Compute maximum difference of sensitivities between the succesive iterations.
    maxDiff <- max(abs(TD$beta - beta0), na.rm = TRUE)
    if (iter == maxCycle && maxDiff > tol) {
      warning(paste0("Convergence not achieved in ", iter," iterations. Tolerance ",
                     tol, ", criterion at last iteration ", signif(maxDiff, 4), ".\n"))
    }
    iter <- iter + 1
  }
  ## Environments.
  aov1 <- anova(model1)
  tPos <- rownames(aov1) == "Residuals"
  rDev[4] <- aov1[["Sum Sq"]][tPos]
  rDf[4] <- aov1[["Df"]][tPos]
  ## Extract total deviance.
  modelA <- lm(as.formula(paste(trait, "~ genotype")), data = TD,
               na.action = na.exclude)
  aovA <- anova(modelA)
  rDev[5] <- sum(aovA[["Sum Sq"]])
  rDf[5] <- sum(aovA[["Df"]])
  ## Fit varieties only for first entry in aov.
  modelB <- lm(as.formula(paste(trait, "~-1 + genotype")), data = TD,
               na.action = na.exclude)
  aovB <- anova(modelB)
  tPos <- rownames(aovB) == "Residuals"
  rDev[1] <- aovB[["Sum Sq"]][tPos]
  rDf[1] <- aovB[["Df"]][tPos]
  ## Calculate deviances and d.f.
  rDev[c(3, 2, 1)] <- rDev[c(2, 1, 5)] - rDev[c(4, 2, 1)]
  rDf[c(2, 1)] <- rDf[c(1, 5)] - rDf[c(2, 1)]
  rDf[3] <- rDf[1]
  rDf[4] <- rDf[5] - rDf[1] - rDf[2] - rDf[3]
  ## Calculate mean deviances and F statistics
  mDev <- rDev / rDf
  rmDev <- mDev[4]
  devr <- mDev / rmDev
  devr[c(4, 5)] <- NA
  devr[rDf == 0] <- NA
  devr[!is.na(devr) & devr < 0] <- NA
  fProb <- pf(devr, rDf, rDf[4], lower.tail = FALSE)
  aovTable <- data.frame("Df" = rDf, "Sum Sq" = rDev, "Mean Sq" = mDev,
                         "F value" = devr, "Pr(>F)" = fProb,
                         row.names = c("genotype", "env", "Sensitivities",
                                       "Residual", "Total"),
                         check.names = FALSE)
  ## Sensitivity beta(i).
  sens <- TD$beta
  ## standard error for env.
  sigmaE <- sqrt(diag(vcov(model1))[match(paste0("genotype", TD$genotype,
                                                 ":envEffs"),
                                          names(coeffsModel1))])
  ## Mean for each genotype genMean(i).
  genMean <- coeffsModel1[match(paste0("genotype", TD$genotype),
                                names(coeffsModel1))]
  ## Residual standard error.
  sigma <- sqrt(diag(vcov(model1))[match(paste0("genotype", TD$genotype),
                                         names(coeffsModel1))])
  fittedGen <- fitted(model1)
  resiGen <- residuals(model1)
  ## mean squared error (MSE) of the trait means for each genotype.
  mse <- tapply(X = resiGen, INDEX = TD$genotype, FUN = function(x) {
    checkG <- length(x)
    if (checkG > 2) {
      sum(x ^ 2, na.rm = TRUE) / (checkG - 2)
    } else {
      rep(NA, checkG)
    }
  })
  ## mean estimates for each genotype.
  G <- levels(TD$genotype)
  sens <- tapply(X = sens, INDEX = TD$genotype, FUN = function(x) {
    mean(x, na.rm = TRUE)
  })
  sigmaE <- tapply(X = sigmaE, INDEX = TD$genotype, FUN = function(x) {
    mean(x, na.rm = TRUE)
  })
  genMean <- tapply(X = genMean, INDEX = TD$genotype, FUN =  function(x) {
    mean(x, na.rm = TRUE)
  })
  sigma <- tapply(X = sigma, INDEX = TD$genotype, FUN = function(x) {
    mean(x, na.rm = TRUE)
  })
  if (sorted == "none") {
    orderSens <- 1:nGeno
  } else {
    orderSens <- order(sens, decreasing = (sorted == "descending"))
  }
  res <- data.frame(genotype = G, sens, sigmaE, genMean, sigma, mse,
                    row.names = 1:length(sens))[orderSens, ]
  ## ANOVA table.
  ## Environment effects.
  matchPos <- match(paste0("env", levels(TD$env), ":beta"), names(coeffsModel2))
  envEffs <- coeffsModel2[matchPos]
  naPos <- is.na(envEffs)
  envEffs[naPos] <- 0
  envEffs <- envEffs - mean(envEffs)
  varEnvEffs <- diag(vcov(model2)[matchPos[!naPos], matchPos[!naPos]])
  seEnvEffs <- rep(NA, length(envEffs))
  names(seEnvEffs) <- names(envEffs)
  seEnvEffs[names(varEnvEffs)] <- sqrt(varEnvEffs)
  matchPos2 <- match(paste0("env", levels(TD$env), ":beta"), names(envEffs))
  if (!is.null(model1$na.action)) {
    meansFitted <- tapply(X = model1$fitted,
                          INDEX = TD$env[-model1$na.action], FUN = mean)
  } else {
    meansFitted <- tapply(X = model1$fitted, INDEX = TD$env, FUN = mean)
  }
  meansFitted <- meansFitted[matchPos2]
  meansRank <- rank(-meansFitted)
  meansNames <- names(meansFitted)
  envEffsSummary <- data.frame(Environment = meansNames, Effect = envEffs,
                               s.e. = seEnvEffs, Mean = meansFitted,
                               Rank = meansRank, row.names = NULL)
  return(createFW(estimates = res, anova = aovTable, envEffs = envEffsSummary,
                  data = TD, fittedGeno = fittedGen,
                  trait = trait, nGeno = nGeno, nEnv = nEnv, tol = tol,
                  iter = iter - 1))
}
