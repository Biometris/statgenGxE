#' Modified joint regression analysis
#'
#' This function performs a modified joint analysis of data classified by two factors.
#'
#' @inheritParams GE.AMMI
#'
#' @param maxcycle A integer specifying the maximum number of iterations to be achieved.
#' By default, \code{maxcycle = 15}.
#' @param tol A small positive numerical value specifying convergence tolerance.
#' By default, \code{tol = 0.001}.
#' @param sortBySens A character string specifying whether the results are to be sorted
#' in an increasing (or decreasing) order of sensitivities.
#' By default, \code{sortBySens = "ascending"}. Other options are "descending" and NA.

#' @examples
#' myDat <- GE.read.csv(system.file("extdata", "F2maize_pheno.csv", package = "RAP"),
#'                      env = "env!", genotype = "genotype!", trait = "yld")
#' myTD <- createTD(data = myDat, genotype = "genotype!", env = "env!")
#' fw.analysis <- GE.FW(myTD, trait = "yld", maxcycle = 15, tol = 0.001, sortBySens = "ascending")
#' fw.analysis
#'
#' @export

GE.FW <- function(TD,
                  trait,
                  maxcycle = 15,
                  tol = 0.001,
                  sortBySens = c("ascending", "descending", NA)) {
  ## ZZ replicates RJOINT procedure in GenStat
  ## Handling missing values ?na.action = na.exclude?
  if (missing(sortBySens)) {
    sortBySens <- "ascending"
  }
  nGeno <- nlevels(TD$genotype)
  nEnv <- nlevels(TD$env)
  ## Save Sum Sq & Df
  rDev <- rDf <- rep(NA, 5)
  ## Initializing...
  ## Estimating env effects with the sensitivity beta=1
  model0 <- lm(as.formula(paste(trait, "~-1 + env + genotype")),
               data = TD, na.action = na.exclude)
  aov0 <- anova(model0)
  tPos <- rownames(aov0) == "Residuals"
  rDev[2] <- aov0[["Sum Sq"]][tPos]
  rDf[2] <- aov0[["Df"]][tPos]
  coeffsModel0 <- coefficients(model0)
  envEffs0 <- coeffsModel0[grep(pattern = "env", x = names(coeffsModel0))]
  # Adjust env effects to have mean zero
  envEffs0 <- envEffs0 - mean(envEffs0, na.rm = TRUE)
  envEffs <- rep(NA, nrow(TD))
  for (ii in names(envEffs0)) {
    tPosition <- sapply(X = TD$env, FUN = function(x) {
      grepl(pattern = x, x = ii)
    })
    envEffs[tPosition] <- envEffs0[ii]
  }
  TD <- cbind(TD, envEffs)
  # Initial values for sensitivity beta
  TD$beta <- 1
  # Set a relative difference to be large
  maxdiff <- 10000
  # Set iteration to be 1
  iter <- 1
  # Iterate to fit 'y(i,j) = genMean(i)+beta(i)*envEffs(j)'
  while (maxdiff > tol && iter <= maxcycle) {
    beta0 <- TD$beta
    # Form variate with current genotype sensitivity relevant to each unit
    model1 <- lm(as.formula(paste(trait, "~-1 + genotype + genotype:envEffs")),
                 data = TD, na.action = na.exclude)
    coeffsModel1 <- coefficients(model1)
    # Update beta
    TD$beta <- coeffsModel1[match(paste0("genotype", TD$genotype, ":envEffs"), names(coeffsModel1))]
    TD$beta <- TD$beta / mean(TD$beta, na.rm = TRUE)
    print(logLik(model1))
    # Form variate with current env means relevant to each unit
    model2 <- lm(as.formula(paste(trait, "~-1 + genotype + env:beta")),
                 data = TD, na.action = na.exclude)
    coeffsModel2 <- coefficients(model2)
    # Update envEffs
    TD$envEffs <- coeffsModel2[match(paste0("env", TD$env, ":beta"), names(coeffsModel2))]
    TD[is.na(TD$envEffs), "envEffs"] <- 0
    TD$envEffs <- TD$envEffs - mean(TD$envEffs)
    # Maximum difference of sensitivities between the succesive iterations
    maxdiff <- max(abs(TD$beta - beta0), na.rm = TRUE)
    if (iter == maxcycle && maxdiff > tol)
      warning(paste0("Convergence not achieved in ", iter," iterations. Tolerance ",
                     tol, ", criterion at last iteration ", signif(maxdiff, 4), ".\n"))
    iter <- iter + 1
  }
  # Environments
  aov1 <- anova(model1)
  tPos <- rownames(aov1) == "Residuals"
  rDev[4] <- aov1[["Sum Sq"]][tPos]
  rDf[4] <- aov1[["Df"]][tPos]
  # Extract total deviance
  modelA <- lm(as.formula(paste(trait, "~ genotype")), data = TD, na.action = na.exclude)
  aovA <- anova(modelA)
  rDev[5] <- sum(aovA[["Sum Sq"]])
  rDf[5] <- sum(aovA[["Df"]])
  # Fit varieties only for first entry in aov
  modelB <- lm(as.formula(paste(trait, "~-1 + genotype")), data = TD, na.action = na.exclude)
  aovB <- anova(modelB)
  tPos <- rownames(aovB) == "Residuals"
  rDev[1] <- aovB[["Sum Sq"]][tPos]
  rDf[1] <- aovB[["Df"]][tPos]
  # Calculate deviances and d.f.
  rDev[c(3, 2, 1)] <- rDev[c(2, 1, 5)] - rDev[c(4, 2, 1)]
  rDf[c(2, 1)] <- rDf[c(1, 5)] - rDf[c(2, 1)]
  rDf[3] <- rDf[1]
  rDf[4] <- rDf[5] - rDf[1] - rDf[2] - rDf[3]
  # Calculate mean deviances and F statistics
  mDev <- rDev / rDf
  rmDev <- mDev[4]
  devr <- mDev / rmDev
  devr[c(4, 5)] <- NA
  devr[rDf == 0] <- NA
  devr[!is.na(devr) & devr < 0] <- NA
  fProb <- pf(devr, rDf, rDf[4], lower.tail = FALSE)
  aovTable <- data.frame("Df" = rDf, "Sum Sq" = rDev, "Mean Sq" = mDev,
                         "F value" = devr, "Pr(>F)" = fProb,
                         row.names = c("genotype", "env", "Sensitivities", "Residual", "Total"),
                         check.names = FALSE)
  # Sensitivity beta(i)
  sens <- TD$beta
  # standard error for env
  sigmaE <- sqrt(diag(vcov(model1))[match(paste0("genotype", TD$genotype, ":envEffs"),
                                          names(coeffsModel1))])
  # Mean for each genotype genMean(i)
  genMean <- coeffsModel1[match(paste0("genotype", TD$genotype), names(coeffsModel1))]
  # residual standard error
  sigma <- sqrt(diag(vcov(model1))[match(paste0("genotype", TD$genotype),
                                         names(coeffsModel1))])
  fittedGen <- fitted(model1)
  resiGen <- residuals(model1)
  # mean squared error (MSE) of the trait means for each genotype
  mse <- tapply(X = resiGen, INDEX = TD$genotype, FUN = function(x) {
    checkG <- length(x)
    if (checkG > 2) {
      sum(x ^ 2, na.rm = TRUE) / (checkG - 2)
    } else {
      rep(NA, checkG)
    }
  })
  # mean estimates for each genotype
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
  if (sortBySens == "ascending") {
    orderSens <- order(sens)
    res <- data.frame(genotype = G, sens, sigmaE, genMean, sigma, mse,
                      row.names = 1:length(sens))[orderSens, ]
  } else if (sortBySens == "descending") {
    orderSens <- order(sens, decreasing = TRUE)
    res <- data.frame(genotype = G, sens, sigmaE, genMean, sigma, mse,
                      row.names = 1:length(sens))[orderSens, ]
  } else {
    res <- data.frame(genotype = G, sens, sigmaE, genMean, sigma, mse,
                      row.names = 1:length(sens))
  }
  # ANOVA table
  # Environment effects
  matchPositions <- match(paste0("env", levels(TD$env), ":beta"), names(coeffsModel2))
  envEffs <- coeffsModel2[matchPositions]
  naPosition <- is.na(envEffs)
  envEffs[naPosition] <- 0
  envEffs <- envEffs - mean(envEffs)
  varEnvEffs <- diag(vcov(model2)[matchPositions[!naPosition], matchPositions[!naPosition]])
  seEnvEffs <- rep(NA, length(envEffs))
  names(seEnvEffs) <- names(envEffs)
  seEnvEffs[names(varEnvEffs)] <- sqrt(varEnvEffs)
  matchPositions2 <- match(paste0("env", levels(TD$env), ":beta"), names(envEffs))
  if (!is.null(model1$na.action)) {
    meansFitted <- tapply(X = model1$fitted, INDEX = TD$env[-model1$na.action], FUN = mean)
  } else {
    meansFitted <- tapply(X = model1$fitted, INDEX = TD$env, FUN = mean)
  }
  meansFitted <- meansFitted[matchPositions2]
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
