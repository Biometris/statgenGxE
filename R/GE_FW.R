#' Modified joint regression analysis
#'
#' This function performs a modified joint analysis of data classified by two factors.
#'
#' @param Y A data frame object.
#' @param trait A character string specifying a trait column of the data.
#' @param genotype A character string specifying a genotype column of the data.
#' @param env A character string specifying the envirnoment column of the data.
#' @param maxcycle A integer specifying the maximum number of iterations to be achieved.
#' By default, \code{maxcycle = 15}.
#' @param tol A small positive numerical value specifying convergence tolerance.
#' By default, \code{tol = 0.001}.
#' @param sortBYsens A character string specifying whether the results are to be sorted
#' in an increasing (or decreasing) order of sensitivities.
#' By default, \code{sortBYsens = "ascending"}. Other options are "descending" and NA.
#' @param scatterplot A logical value specifying if a scatterplot of sensitivities is produced.
#' @param lineplot A logical value specifying if a fitted line for each genotype is produced.
#' @param trellisplot A logical value specifying if the trellis plot of the individual genotype slopes
#' is produced.
#'
#' @examples
#' mydat <- GE.read.csv(file.path(path.package("RAP"),"F2maize_pheno.csv"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' names(mydat)=c("env", "genotype","yld")
#' fw.anlysis <- GE.FW(mydat, trait="yld", genotype="genotype", env="env", maxcycle = 15, tol = 0.001,
#'                sortBYsens = "ascending")
#' fw.anlysis
#'
#' @import graphics grDevices
#' @export

GE.FW <- function(Y,
                  trait,
                  genotype,
                  env,
                  maxcycle = 15,
                  tol = 0.001,
                  sortBYsens = c("ascending", "descending", NA),
                  scatterplot = TRUE,
                  lineplot = FALSE,
                  trellisplot = FALSE) {
  ## ZZ replicates RJOINT procedure in GenStat
  ## Handling missing values ?na.action = na.exclude?
  if (missing(sortBYsens)) {
    sortBYsens <- "ascending"
  }
  nLab <- nlevels(Y[, genotype])
  nEnvs <- nlevels(Y[, env])
  ## Save Sum Sq & Df
  rDev <- rDf <- rep(NA, 5)
  ## Initializing...
  ## Estimating env effects with the sensitivity beta=1
  model0 <- lm(as.formula(paste(trait, "~-1 +", env, "+", genotype)),
               data = Y, na.action = na.exclude)
  aov0 <- anova(model0)
  tPos <- rownames(aov0) == "Residuals"
  rDev[2] <- aov0[["Sum Sq"]][tPos]
  rDf[2] <- aov0[["Df"]][tPos]
  coeffsModel0 <- coefficients(model0)
  envEffs0 <- coeffsModel0[grep(pattern = env, x = names(coeffsModel0))]
  # Adjust env effects to have mean zero
  envEffs0 <- envEffs0 - mean(envEffs0, na.rm = TRUE)
  envEffs <- rep(NA, nrow(Y))
  for (ii in names(envEffs0)) {
    tPosition <- sapply(Y[[env]], function(x) {
      grepl(pattern = x, x = ii)
    })
    envEffs[tPosition] <- envEffs0[ii]
  }
  Y <- cbind(Y, envEffs)
  # Initial values for sensitivity beta
  Y$beta <- rep(1, nLab)
  # Set a relative difference to be large
  maxdiff <- 10000
  # Set iteration to be 1
  iter <- 1
  # Iterate to fit 'y(i,j) = genMean(i)+beta(i)*envEffs(j)'
  while (maxdiff > tol && iter <= maxcycle) {
    beta0 <- Y$beta
    # Form variate with current genotype sensitivity relevant to each unit
    model1 <- lm(as.formula(paste(trait, "~-1 +", genotype, "+", genotype, ":envEffs")),
                 data = Y, na.action = na.exclude)
    coeffsModel1 <- coefficients(model1)
    # Update beta
    Y$beta <- coeffsModel1[match(paste0(genotype, Y[[genotype]], ":envEffs"), names(coeffsModel1))]
    Y$beta <- Y$beta / mean(Y$beta, na.rm = TRUE)
    print(logLik(model1))
    # Form variate with current env means relevant to each unit
    model2 <- lm(as.formula(paste(trait, "~-1 +", genotype, "+", env, ":beta")),
                 data = Y, na.action = na.exclude)
    coeffsModel2 <- coefficients(model2)
    # Update envEffs
    Y$envEffs <- coeffsModel2[match(paste0(env, Y[[env]],":beta"), names(coeffsModel2))]
    Y[is.na(Y$envEffs), "envEffs"] <- 0
    Y$envEffs <- Y$envEffs - mean(Y$envEffs)
    # Maximum difference of sensitivities between the succesive iterations
    maxdiff <- max(abs(Y$beta - beta0), na.rm = TRUE)
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
  modelA <- lm(as.formula(paste(trait, "~", genotype)), data = Y, na.action = na.exclude)
  aovA <- anova(modelA)
  rDev[5] <- sum(aovA[["Sum Sq"]])
  rDf[5] <- sum(aovA[["Df"]])
  # Fit varieties only for first entry in aov
  modelB <- lm(as.formula(paste(trait, "~-1 +", genotype)), data = Y, na.action = na.exclude)
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
  aovtable <- data.frame("Df" = rDf, "Sum Sq" = rDev, "Mean Sq" = mDev,
                         "F value" = devr, "Pr(>F)" = fProb,
                         row.names = c(genotype, env, "Sensitivities", "Residual", "Total"),
                         check.names = FALSE)
  # Sensitivity beta(i)
  sens <- Y$beta
  # standard error for env
  sigmaE <- sqrt(diag(vcov(model1))[match(paste0(genotype, Y[[genotype]], ":envEffs"),
                                          names(coeffsModel1))])
  # Mean for each genotype genMean(i)
  genMean <- coeffsModel1[match(paste0(genotype, Y[[genotype]]), names(coeffsModel1))]
  # residual standard error
  sigma <- sqrt(diag(vcov(model1))[match(paste0(genotype, Y[[genotype]]),
                                         names(coeffsModel1))])
  fittedGen <- fitted(model1)
  resiGen <- residuals(model1)
  fVal <- tapply(X = fittedGen, INDEX = Y[, c(genotype, env)], FUN = function(x) {
    mean(x, na.rm = TRUE)
  })
  # mean squared error (MSE) of the trait means for each genotype
  mse <- tapply(X = resiGen, INDEX = Y[, genotype], FUN = function(x) {
    checkG <- length(x)
    if (checkG > 2) {
      sum(x ^ 2, na.rm = TRUE) / (checkG - 2)
    } else {
      rep(NA, checkG)
    }
  })
  # mean estimates for each genotype
  G <- levels(Y[, genotype])
  sens <- tapply(X = sens, INDEX = Y[, genotype], FUN = function(x) {
    mean(x, na.rm = TRUE)
  })
  sigmaE <- tapply(X = sigmaE, INDEX = Y[, genotype], FUN = function(x) {
    mean(x, na.rm = TRUE)
  })
  genMean <- tapply(X = genMean, INDEX = Y[, genotype],FUN =  function(x) {
    mean(x, na.rm = TRUE)
  })
  sigma <- tapply(X = sigma, INDEX = Y[, genotype], FUN = function(x) {
    mean(x, na.rm = TRUE)
  })
  if (sortBYsens == "ascending") {
    orderSens <- order(sens)
    res <- data.frame(G, sens, sigmaE, genMean, sigma, mse,
                      row.names = 1:length(sens))[orderSens, ]
  } else if (sortBYsens == "descending") {
    orderSens <- order(sens, decreasing = TRUE)
    res <- data.frame(G, sens, sigmaE, genMean, sigma, mse,
                      row.names=1:length(sens))[orderSens, ]
  } else {
    res <- data.frame(G, sens, sigmaE, genMean, sigma, mse, row.names=1:length(sens))
  }
  # ANOVA table
  # Environment effects
  matchPositions <- match(paste0(env, levels(Y[[env]]), ":beta"), names(coeffsModel2))
  envEffs <- coeffsModel2[matchPositions]
  naPosition <- is.na(envEffs)
  envEffs[naPosition] <- 0
  envEffs <- envEffs - mean(envEffs)
  varEnvEffs <- diag(vcov(model2)[matchPositions[!naPosition], matchPositions[!naPosition]])
  seEnvEffs <- rep(NA, length(envEffs))
  names(seEnvEffs) <- names(envEffs)
  seEnvEffs[names(varEnvEffs)] <- sqrt(varEnvEffs)
  matchPositions2 <- match(paste0(env, levels(Y[[env]]), ":beta"), names(envEffs))
  if (!is.null(model1$na.action)) {
    meansFitted <- tapply(X = model1$fitted, INDEX = Y[[env]][-model1$na.action], FUN = mean)
  } else {
    meansFitted <- tapply(X = model1$fitted, INDEX = Y[[env]], FUN = mean)
  }
  meansFitted <- meansFitted[matchPositions2]
  meansRank <- rank(-meansFitted)
  meansNames <- names(meansFitted)
  envEffsSummary <- data.frame(Environment = meansNames, Effect = envEffs,
                               s.e. = seEnvEffs, Mean = meansFitted,
                               Rank = meansRank, row.names = NULL)
  # Draw various plots
  if (scatterplot) {
    if (!all(is.na(mse))) {
      xx <- cbind(genMean, mse, sens)
      colnames(xx) <- c("Mean", "m.s.deviation", "Sensitivity")
      pairs(xx, upper.panel = NULL, main = "Finlay & Wilkinson analysis")
    } else {
      plot(x = genMean, y = sens, xlab = "Mean", ylab = "Sensitivity",
           main = "Finlay & Wilkinson analysis")
    }
  }
  if (lineplot){
    minFVal <- min(fittedGen, na.rm = TRUE)
    maxFVal <- max(fittedGen, na.rm = TRUE)
    xEff <- envEffs
    minXEff <- min(xEff, na.rm = TRUE)
    maxXEff <- max(xEff, na.rm = TRUE)
    dev.new()
    plot(x = NA, xlim = c(minXEff, maxXEff), ylim = c(minFVal, maxFVal),
         ylab = trait, xlab = "Environment", xaxt = "n")
    axis(side = 1, xEff, levels(Y[, env]), las = 2, cex.axis = .75)
    colour <- 1
    for (i in 1:nLab) {
      if (!is.na(sortBYsens)) {
        xfVal <- fVal[orderSens[i], ]
      } else {
        xfVal <- fVal[i, ]
      }
      lines(xEff[order(xEff)], xfVal[order(xEff)], col = colour)
      colour <- colour + 1
    }
  }
  if (trellisplot){
    trellisdata <- data.frame(genotype = Y[[genotype]], trait = Y[[trait]],
                              fitted = fittedGen, xEff = rep(envEffs, nLab))
    if (nLab > 64) {
      first64 <- levels(Y[[genotype]])[1:64]
      first64 <- Y[[genotype]] %in% first64
      trellisdata <- droplevels(trellisdata[first64, ])
    }
    dev.new()
    print(lattice::xyplot(trait + fitted ~ xEff | genotype, data = trellisdata,
                          panel = function(x, y, subscripts) {
                            lattice::panel.xyplot(x,y)
                            lattice::panel.lines(trellisdata$xEff[subscripts],
                                                 trellisdata$fitted[subscripts])
                          }, as.table = TRUE, subscripts = TRUE,
                          xlab = "Environment", ylab = trait,
                          main = paste0("Finlay & Wilkinson analysis for ", trait)))
  }
  info <- data.frame("Response variate" = trait,
                     "Number of genotypes" = nLab,
                     "Number of environments" = nEnvs,
                     "Convergence criterion" = tol,
                     "Number of iterations" = iter - 1,
                     row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
  result <- list(estimates = res, aovTable = aovtable, envEffs = envEffsSummary)
  attr(result, "info") <- info
  return(result)
}
