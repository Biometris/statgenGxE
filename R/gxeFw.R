#' Finlay-Wilkinson analysis
#'
#' This function performs a Finlay-Wilkinson analysis of data classified by two
#' factors.
#'
#' @inheritParams gxeAmmi
#'
#' @param maxIter An integer specifying the maximum number of iterations in
#' the algorithm.
#' @param tol A positive numerical value specifying convergence tolerance of the
#' algorithm.
#' @param sorted A character string specifying the sorting order of the
#' estimated values in the output.
#'
#' @return An object of class \code{\link{FW}}, a list containing:
#' \item{estimates}{A data.frame containing the estimated values.}
#' \item{anova}{A data.frame containing anova scores of the FW analysis.}
#' \item{envEffs}{A data.frame containing the environmental effects.}
#' \item{TD}{The object of class TD on which the analysis was performed.}
#' \item{fittedGeno}{A numerical vector containing the fitted values for the
#' genotypes.}
#' \item{trait}{A character string containing the analyzed trait.}
#' \item{nGeno}{A numerical value containing the number of genotypes in the
#' analysis.}
#' \item{nEnv}{A numerical value containing the number of environments in the
#' analysis.}
#' \item{tol}{A numerical value containing the tolerance used during the
#' analysis.}
#' \item{iter}{A numerical value containing the number of iterations for the
#' analysis to converge.}
#'
#' @references Finlay, K.W. & Wilkinson, G.N. (1963). The analysis of adaptation
#' in a plant-breeding programme. Australian Journal of Agricultural
#' Research, 14, 742-754.
#'
#' @seealso \code{\link{FW}}, \code{\link{plot.FW}}, \code{\link{report.FW}}
#'
#' @examples
#' ## Run Finlay-Wilkinson analysis on TDMaize.
#' geFW <- gxeFw(TDMaize, trait = "yld")
#' ## Summarize results.
#' summary(geFW)
#' ## Create a scatterplot of the results.
#' plot(geFW, plotType = "scatter")
#' \dontrun{
#' ## Create report.
#' report(geFW, outfile = "./testReports/reportFW.pdf")
#' }
#'
#' @export
gxeFw <- function(TD,
                  trait,
                  maxIter = 15,
                  tol = 0.001,
                  sorted = c("ascending", "descending", "none")) {
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (!"env" %in% colnames(TD)) {
    stop("TD should contain a column env to be able to run a Finlay
         Wilkinson analysis.\n")
  }
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TD)) {
    stop("trait has to be a column in TD.\n")
  }
  if (is.null(maxIter) || !is.numeric(maxIter) || length(maxIter) > 1 ||
      maxIter != round(maxIter) || maxIter < 1) {
    stop("maxIter should be a positive integer.\n")
  }
  if (is.null(tol) || !is.numeric(tol) || length(tol) > 1 || tol < 0) {
    stop("tol should be a numerical value > 10^-6.\n")
  }
  sorted <- match.arg(sorted)
  nGeno <- nlevels(TD$genotype)
  ## Setup empty vectors for storing rDev and rDF
  rDev <- rDf <- rep(NA, 5)
  ## Estimate env effects with the sensitivity beta = 1.
  model0 <- lm(as.formula(paste(trait, "~-1 + env + genotype")),
               data = TD, na.action = na.exclude)
  aov0 <- anova(model0)
  rDev[2] <- aov0["Residuals", "Sum Sq"]
  rDf[2] <- aov0["Residuals", "Df"]
  coeffsModel0 <- coefficients(model0)
  envEffs0 <- coeffsModel0[grep(pattern = "env", x = names(coeffsModel0))] %>%
    scale(scale = FALSE) %>% `colnames<-`("envEffs") %>%
    `rownames<-`(substring(rownames(.), first = 4))
  TD <- merge(x = TD, y = envEffs0, by.x = "env", by.y = "row.names")
  ## Set initial values for sensitivity beta.
  TD$beta <- 1
  ## Set a relative difference to be large.
  maxDiff <- Inf
  ## Set iteration to 1.
  iter <- 1
  ## Iterate to fit 'y(i,j) = genMean(i)+beta(i)*envEffs(j)'
  while (maxDiff > tol && iter <= maxIter) {
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
    model2 <- lm(as.formula(paste(trait, "~ -genotype + env:beta")),
                 data = TD, na.action = na.exclude)
    coeffsModel2 <- coefficients(model2)
    ## Update envEffs.
    TD$envEffs <- coeffsModel2[match(paste0("env", TD$env, ":beta"),
                                     names(coeffsModel2))]
    TD[is.na(TD$envEffs), "envEffs"] <- 0
    TD$envEffs <- TD$envEffs - mean(TD$envEffs)
    ## Compute max difference of sensitivities between the succesive iterations.
    maxDiff <- max(abs(TD$beta - beta0), na.rm = TRUE)
    if (iter == maxIter && maxDiff > tol) {
      warning(paste0("Convergence not achieved in ", iter,
                     " iterations. Tolerance ", tol,
                     ", criterion at last iteration ", signif(maxDiff, 4),
                     ".\n"))
    }
    iter <- iter + 1
  }
  ## Environments.
  aov1 <- anova(model1)
  rDev[4] <- aov1["Residuals","Sum Sq"]
  rDf[4] <- aov1["Residuals", "Df"]
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
  rDev[1] <- aovB["Residuals", "Sum Sq"]
  rDf[1] <- aovB["Residuals", "Df"]
  ## Calculate deviances and d.f.
  rDev[c(3, 2, 1)] <- rDev[c(2, 1, 5)] - rDev[c(4, 2, 1)]
  rDf[c(2, 1)] <- rDf[c(1, 5)] - rDf[c(2, 1)]
  rDf[3] <- rDf[1]
  rDf[4] <- rDf[5] - rDf[1] - rDf[2] - rDf[3]
  ## Calculate mean deviances and F statistics.
  mDev <- rDev / rDf
  devr <- mDev / mDev[4]
  devr[c(4, 5)] <- NA
  devr[rDf == 0] <- NA
  devr[!is.na(devr) & devr < 0] <- NA
  fProb <- pf(q = devr, df1 = rDf, df2 = rDf[4], lower.tail = FALSE)
  aovTable <- data.frame("Df" = rDf, "Sum Sq" = rDev, "Mean Sq" = mDev,
                         "F value" = devr, "Pr(>F)" = fProb,
                         row.names = c("genotype", "env", "Sensitivities",
                                       "Residual", "Total"),
                         check.names = FALSE)
  ## Extract sensitivity beta.
  sens <- as.vector(tapply(X = TD$beta, INDEX = TD$genotype,
                           FUN = mean, na.rm = TRUE))
  ## Compute the standard errors.
  sigmaE <- sqrt(diag(vcov(model1))[match(paste0("genotype", TD$genotype,
                                                 ":envEffs"),
                                          names(coeffsModel1))]) %>%
    tapply(INDEX = TD$genotype, FUN = mean, na.rm = TRUE) %>%
    as.vector
  ## Extract the mean for each genotype.
  genMean <- coeffsModel1[match(paste0("genotype", TD$genotype),
                                names(coeffsModel1))] %>%
    tapply(INDEX = TD$genotype, FUN = mean, na.rm = TRUE) %>%
    as.vector
  ## Residual standard error.
  sigma <- sqrt(diag(vcov(model1))[match(paste0("genotype", TD$genotype),
                                         names(coeffsModel1))]) %>%
    tapply(INDEX = TD$genotype, FUN = mean, na.rm = TRUE) %>%
    as.vector
  ## Compute mean squared error (MSE) of the trait means for each genotype.
  mse <- residuals(model1) %>%
    tapply(INDEX = TD$genotype, FUN = function(x) {
      checkG <- length(x)
      if (checkG > 2) {
        sum(x ^ 2, na.rm = TRUE) / (checkG - 2)
      } else {
        rep(NA, checkG)
      }
    }) %>% as.vector
  ## Compute sortting order for estimates.
  if (sorted == "none") {
    orderSens <- 1:nGeno
  } else {
    orderSens <- order(sens, decreasing = (sorted == "descending"))
  }
  ## Construct estimate data.frame.
  estimates <- data.frame(genotype = levels(TD$genotype), sens, sigmaE,
                          genMean, sigma, mse,
                          row.names = 1:length(sens))[orderSens, ]
  ## Construct data.frame with environment effects.
  matchPos <- match(paste0("env", levels(TD$env), ":beta"), names(coeffsModel2))
  envEffs <- coeffsModel2[matchPos]
  naPos <- is.na(envEffs)
  envEffs[naPos] <- 0
  envEffs <- envEffs - mean(envEffs)
  seEnvEffs <- sqrt(diag(vcov(model2)[matchPos[!naPos], matchPos[!naPos]]))
  matchPos2 <- match(paste0("env", levels(TD$env), ":beta"), names(envEffs))
  if (!is.null(model1$na.action)) {
    meansFitted <- tapply(X = model1$fitted,
                          INDEX = TD$env[-model1$na.action], FUN = mean)
  } else {
    meansFitted <- tapply(X = model1$fitted, INDEX = TD$env, FUN = mean)
  }
  meansFitted <- meansFitted[matchPos2]
  envEffsSummary <- data.frame(env = names(meansFitted), effect = envEffs,
                               SE = seEnvEffs, mean = as.vector(meansFitted),
                               rank = rank(-meansFitted), row.names = NULL)
  return(createFW(estimates = estimates, anova = aovTable,
                  envEffs = envEffsSummary, TD = createTD(TD),
                  fittedGeno = unname(fitted(model1)), trait = trait,
                  nGeno = nGeno, nEnv = nlevels(TD$env), tol = tol,
                  iter = iter - 1))
}
