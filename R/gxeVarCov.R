#' Selects the best variance-covariance model for a set of trials
#'
#' This function selects the best covariance structure for genetic correlations
#' between trials. It fits a range of variance-covariance models (identity,
#' compound symmetry (cs), diagonal, simple correlation with heterogeneous
#' variance (outside), heterogeneous compound symmetry (hcs),
#' first order factor analytic (fa), second order factor analytic (fa2) and
#' unstructured), and selects the best one using a goodness-of-fit criterion.
#' See details for the exact models fitted.
#'
#' The models fitted are of the form \eqn{y_{ij} = \mu_j + \epsilon_{ij}}, where
#' \eqn{y_{ij}} is the phenotypic value of genotype \eqn{i} in environment
#' \eqn{j}, \eqn{\mu_j}
#' is the environmental mean, and \eqn{\epsilon_{ij}} represents mainly genetic
#' variation, although some non-genetic variation may be included as well.
#' The random term \eqn{\epsilon_{ij}} is modeled in eight ways as described in
#' the table below.
#'
#' | Model | Description | var(\eqn{g_{ij}}) | cov(\eqn{g_{ij}};\eqn{g_{ik}}) | Number of parameters|
#' |:-----------|:-------------------|:------------|:------------|:-----------------|
#' | identity | identity | \eqn{\sigma_G^2} | 0 | 1 |
#' | cs | compound symmetry | \eqn{\sigma_G^2+\sigma_{GE}^2} | \eqn{\sigma_{GE}^2} | 2 |
#' | diagonal | diagonal matrix (heteroscedastic) | \eqn{\sigma_{GE_j}^2} | 0 | \eqn{J} |
#' | hcs | heterogeneous compound symmetry | \eqn{\sigma_G^2+\sigma_{GE_j}^2} | \eqn{\sigma_G^2} | \eqn{J+1} |
#' | outside | heterogeneity outside | \eqn{\sigma_{G_j}^2} | \eqn{\theta} | \eqn{J+1} |
#' | fa | first order factor analytic | \eqn{\lambda_{1j}^2+\sigma_{GE_j}^2} | \eqn{\lambda_{1j}\lambda_{1k}} | \eqn{2J} |
#' | fa2 | second order factor analytic | \eqn{\lambda_{1j}^2+\lambda_{2j}^2+\sigma_{GE_j}^2} | \eqn{\lambda_{1j}\lambda_{1k}+\lambda_{2j}\lambda_{2k}} | \eqn{3J-1} |
#' | unstructured | unstructured | \eqn{\sigma_{G_j}^2} | \eqn{\sigma_{G_{j,k}}^2} | \eqn{J(J+1)/2} |
#'
#' In this table \eqn{J} is the number of environments, \eqn{\sigma_G^2} the
#' variance component for the genotype main effects, \eqn{\sigma_{GE}^2} the
#' variance component for GxE interactions. \eqn{\sigma_{G_j}^2} and
#' \eqn{\sigma_{GE_j}^2} are the environment specific variance components for
#' the genotype main effects and GxE interaction in environment \eqn{j}.
#' \eqn{\sigma_{G_{j,k}}^2} is the genetic covariance between environments
#' \eqn{j} and \eqn{k}. \eqn{\theta} is the common correlation between
#' environments and \eqn{\lambda_{1j}} and \eqn{\lambda_{2j}} are
#' environment specific multiplicative parameters.
#'
#' @inheritParams gxeAmmi
#'
#' @param engine A character string specifying the engine used for modeling.
#' Either "lme4" or "asreml".
#' @param criterion A string specifying a goodness-of-fit criterion. Either
#' "AIC" or "BIC".
#' @param models A character vector specifying the models to be fitted.
#' @param ... Further arguments to be passed to the modeling engine.
#'
#' @note If \code{engine = "lme4"}, only the compound symmetry model can be
#' fitted.
#'
#' @return An object of class \code{\link{varCov}}, a list object containing:
#' \item{STA}{An object of class STA containing the best fitted model.}
#' \item{choice}{A character string indicating the best fitted model.}
#' \item{summary}{A data.frame with a summary of the fitted models.}
#' \item{vcov}{The covariance matrix of the best fitted model.}
#' \item{criterion}{A character string indicating the goodness-of-fit criterion
#' used for determining the best model, either "AIC" or "BIC".}
#' \item{engine}{A character string containing the engine used for
#' the analysis.}
#' \item{trait}{A character string containing the trait analyzed.}
#' \item{dat}{A data.frame with the full data set used for the analysis.}
#'
#' @examples
#' ## Select the best variance-covariance model using lme4 for modeling.
#' geVarCov <- gxeVarCov(TD = TDMaize, trait = "yld")
#'
#' ## Summarize results.
#' summary(geVarCov)
#'
#' \donttest{
#' ## Create a pdf report summarizing the results.
#' report(geVarCov, outfile = tempfile(fileext = ".pdf"))
#' }
#'
#' \donttest{
#' if (requireNamespace("asreml", quietly = TRUE)) {
#'   ## Select the best variance-covariance model using asreml for modeling.
#'   ## Use BIC as a goodness-of-fit criterion.
#'   geVarCov2 <- gxeVarCov(TD = TDMaize, trait = "yld", engine = "asreml",
#'                          criterion = "BIC")
#'
#'   summary(geVarCov2)
#'
#'   ## Plot a heatmap of the correlation matrix for the best model.
#'   plot(geVarCov2)
#'   }
#' }
#'
#' @family varCov
#'
#' @export
gxeVarCov <- function(TD,
                      trials = names(TD),
                      trait,
                      models = c("identity", "cs", "diagonal", "hcs", "outside",
                                 "fa", "fa2", "unstructured"),
                      engine = c("lme4", "asreml"),
                      criterion = c("BIC", "AIC"),
                      ...) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  trials <- chkTrials(trials, TD)
  TDTot <- Reduce(f = rbind, x = TD[trials])
  chkCol(trait, TDTot)
  chkCol("trial", TDTot)
  chkCol("genotype", TDTot)
  engine <- match.arg(engine)
  criterion <- match.arg(criterion)
  models <- match.arg(models, several.ok = TRUE)

  if (engine == "asreml") {
    choices <- models
  } else {
    if (!"cs" %in% models) {
      stop("With lme4 only compound symmetry (cs) model can be fitted.\n")
    }
    choices <- "cs"
  }
  TDTot <- droplevels(TDTot)
  ## Increase maximum number of iterations for asreml. Needed for more complex
  ## designs to converge.
  maxIter <- 200
  ## Add combinations of trial and genotype currently not in TD to TD.
  ## No missing combinations are allowed when fitting asreml models.
  TD0 <- expand.grid(genotype = levels(TDTot[["genotype"]]),
                     trial = levels(TDTot[["trial"]]))
  TDTot <- merge(TD0, TDTot, all.x = TRUE)
  ## Create matrix for storing results.
  bestTab <- matrix(nrow = length(choices), ncol = 4,
                    dimnames = list(choices,
                                    c("AIC", "BIC", "Deviance", "NParameters")))
  models <- setNames(vector(mode = "list", length = length(choices)), choices)
  ## Main procedure to fit mixed models.
  if (engine == "lme4") {
    ## Compound symmetry ("cs") only.
    mr = lme4::lmer(formula(paste0("`", trait, "`~ trial + (1 | genotype)")),
                    data = TDTot, ...)
    nPar <- 2
    ## Construct STA object.
    models[["cs"]] <- mr
    bestTab["cs", ] <- c(-2 * as.numeric(logLik(mr)) + 2 * nPar,
                         -2 * as.numeric(logLik(mr)) +
                           log(length(fitted(mr))) * nPar,
                         -2 * logLik(mr), nPar)
    bestModel <- models[[rownames(bestTab)[1]]]
    vcovBest <- vcov(emmeans::emmeans(bestModel, specs = "trial",
                                      lmer.df = "asymptotic"))
    colnames(vcovBest) <- rownames(vcovBest) <- levels(TDTot[["trial"]])
  } else if (engine == "asreml") {
    if (requireNamespace("asreml", quietly = TRUE)) {
      nTr <- nlevels(TDTot[["trial"]])
      fixedForm <- formula(paste0("`", trait, "`~ trial"))
      ## Put arguments for models in a list to make it easier to switch
      ## between asreml3 and asreml4. Usually only one or two arguments differ.
      ## Also some arguments are identical for all models
      modArgs0 <- list(fixed = fixedForm, data = TDTot, maxiter = maxIter,
                       trace = FALSE)
      for (choice in choices) {
        if (choice == "identity") {
          modArgs <- modArgs0
          modArgs[[ifelse(asreml4(), "residual", "rcov")]] <-
            formula("~genotype:trial")
          mr <- tryCatchExt(do.call(asreml::asreml, modArgs))
          nPar <- 1
        } else if (choice == "cs") {
          modArgs <- c(modArgs0, list(random = formula("~genotype")))
          modArgs[[ifelse(asreml4(), "residual", "rcov")]] <-
            formula("~ genotype:trial")
          mr <- tryCatchExt(do.call(asreml::asreml, modArgs))
          nPar <- 2
        } else if (choice == "diagonal") {
          modArgs <- modArgs0
          modArgs[[ifelse(asreml4(), "residual", "rcov")]] <-
            formula("~genotype:diag(trial)")
          mr <- tryCatchExt(do.call(asreml::asreml, modArgs))
          nPar <- nTr
        } else if (choice == "hcs") {
          modArgs <- c(modArgs0, list(random = formula("~genotype")))
          modArgs[[ifelse(asreml4(), "residual", "rcov")]] <-
            formula("~ genotype:diag(trial)")
          mr <- tryCatchExt(do.call(asreml::asreml, modArgs))
          nPar <- nTr + 1
        } else if (choice == "outside") {
          modArgs <- c(modArgs0,
                       list(random = formula("~genotype:corh(trial)"),
                            start.values = TRUE))
          startVals <- do.call(asreml::asreml, modArgs)
          tmpValues <- initVals(TD = TDTot, trait = trait, vcmodel = "outside",
                                fixed = fixedForm, ...)
          tmpTable <- startVals[[ifelse(asreml4(), "vparameters.table",
                                        "gammas.table")]]
          tmpTable[, "Value"] <- c(tmpValues$vg,
                                   scale(tmpValues$diag, center = FALSE), 1)
          modArgs <- c(modArgs, list(G.param = tmpTable, R.param = tmpTable,
                                     workspace = 160e6))
          ## Putting start.values to FALSE instead of removing it from the list
          ## crashes asreml.
          modArgs[["start.values"]] <- NULL
          mr <- tryCatchExt(do.call(asreml::asreml, modArgs))
          nPar <- nTr + 1
        } else if (choice == "fa" && nTr > 4) {
          modArgs <- c(modArgs0,
                       list(random = formula("~genotype:fa(trial, 1)"),
                            start.values = TRUE))
          startVals <- do.call(asreml::asreml, modArgs)
          tmpTable <- startVals[[ifelse(asreml4(), "vparameters.table",
                                        "gammas.table")]]
          tmpValues <- initVals(TD = TDTot, trait = trait, vcmodel = "fa",
                                fixed = fixedForm, ...)
          if (!is.null(tmpValues)) {
            tmpTable[, "Value"] <- c(scale(tmpValues$psi, center = FALSE),
                                     tmpValues$gamma, 1)
          }
          modArgs <- c(modArgs, list(G.param = tmpTable, R.param = tmpTable,
                                     workspace = 160e6))
          ## Putting start.values to FALSE instead of removing it from the list
          ## crashes asreml.
          modArgs[["start.values"]] <- NULL
          mr <- tryCatchExt(do.call(asreml::asreml, modArgs))
          nPar <- nTr * 2
        } else if (choice == "fa2" && nTr > 4) {
          modArgs <- c(modArgs0,
                       list(random = formula("~genotype:fa(trial, 2)"),
                            start.values = TRUE))
          startVals <- do.call(asreml::asreml, modArgs)
          tmpTable <- startVals[[ifelse(asreml4(), "vparameters.table",
                                        "gammas.table")]]
          tmpValues <- initVals(TD = TDTot, trait = trait, vcmodel = "fa2",
                                fixed = fixedForm, ...)
          if (!is.null(tmpValues)) {
            ## Keep loadings of factor 2 away from 0.
            tmpValues$gamma[2, tmpValues$gamma[2, ] < 1e-3] <- 1e-3
            ## Make sure that first entry is 0.
            tmpValues$gamma[2, 1] <- 0
            tmpTable[, "Value"] <- c(scale(tmpValues$psi, center = FALSE),
                                     tmpValues$gamma[1, ],
                                     tmpValues$gamma[2, ], 1)
          }
          modArgs <- c(modArgs, list(G.param = tmpTable, R.param = tmpTable,
                                     workspace = 160e6))
          ## Putting start.values to FALSE instead of removing it from the list
          ## crashes asreml.
          modArgs[["start.values"]] <- NULL
          mr <- tryCatchExt(do.call(asreml::asreml, modArgs))
          nPar <- nTr * 3 - 1
        } else if (choice == "unstructured") {
          modArgs <- c(modArgs0,
                       list(random = formula("~genotype:us(trial)"),
                            start.values = TRUE))
          startVals <- do.call(asreml::asreml, modArgs)
          tmpValues <- initVals(TD = TDTot, trait = trait,
                                vcmodel = "unstructured", fixed = fixedForm,
                                ...)
          tmpTable <- startVals[[ifelse(asreml4(), "vparameters.table",
                                        "gammas.table")]]
          tmpValues <- tmpValues$evCov[upper.tri(tmpValues$evCov, diag = TRUE)]
          tmpTable[, "Value"] <- c(scale(tmpValues, center = FALSE),
                                   10 ^ min(floor(log10(min(abs(tmpValues)))), 0))
          ## All off diagonal elements are unconstrained, diagonal elements
          ## should be positive
          tmpTable[-c((1:nTr) * ((1:nTr) + 1) / 2), "Constraint"] <- "U"
          ## Fix residual variance.
          ## Column and component names are different for asreml3/asreml4.
          if (asreml4()) {
            tmpTable[tmpTable[["Component"]] == "units!R", "Constraint"] <- "F"
          } else {
            tmpTable[tmpTable[["Gamma"]] == "R!variance", "Constraint"] <- "F"
          }
          modArgs <- c(modArgs, list(G.param = tmpTable, R.param = tmpTable,
                                     workspace = 160e6))
          ## Putting start.values to FALSE instead of removing it from the list
          ## crashes asreml.
          modArgs[["start.values"]] <- NULL
          mr <- tryCatchExt(do.call(asreml::asreml, modArgs))
          nPar <- nTr * (nTr + 1) / 2
        }
        if (!is.null(mr$warning)) {
          ## Check if param 1% increase is significant. Remove warning if not.
          mr <- chkLastIter(mr)
        }
        if (length(mr$warning) != 0) {
          warning("Asreml gave the following warning for ", choice, ":\n",
                  mr$warning, "\n", call. = FALSE)
        }
        if (!is.null(mr$error)) {
          warning("Asreml gave the following error for ", choice, ":\n",
                  mr$error, call. = FALSE)
          mr <- list(loglik = -Inf)
        } else if (!(nTr <= 4 && choice %in% c("fa", "fa2"))) {
          mr <- mr$value
          mr$call$fixed <- eval(mr$call$fixed)
          mr$call$random <- eval(mr$call$random)
          mr$call$rcov <- eval(mr$call$rcov)
          mr$call$G.param <- eval(mr$call$G.param)
          mr$call$R.param <- eval(mr$call$R.param)
          if (!mr$converge) {
            warning("No convergence for ", choice, ".\n", call. = FALSE)
            mr$loglik <- -Inf
          }
          models[[choice]] <- mr
          bestTab[choice, "AIC"] <- -2 * mr$loglik + 2 * nPar
          bestTab[choice, "BIC"] <- -2 * mr$loglik +
            (log(length(fitted(mr)) - nlevels(TDTot[["genotype"]])) * nPar)
          bestTab[choice, "Deviance"] <- -2 * mr$loglik
          bestTab[choice, "NParameters"] <- nPar
        }
      } # End loop over choices
      bestTab <- bestTab[order(bestTab[, criterion]), , drop = FALSE]
      bestModel <- models[[rownames(bestTab)[1]]]
      bestPred <- predictAsreml(model = bestModel, classify = "trial",
                                TD = TDTot, maxiter = maxIter, ...)
      vcovBest <- if (asreml4()) {
        as.matrix(bestPred$vcov)
      } else {
        bestPred$predictions$vcov
      }
      colnames(vcovBest) <- rownames(vcovBest) <- levels(TDTot[["trial"]])
    } else {
      stop("Failed to load 'asreml'.\n")
    }
  }
  modDat <- TDTot
  modDat[["trial"]] <- rownames(bestTab)[1]
  ## Create output.
  model <- setNames(list(list(mRand = NULL,
                              mFix = setNames(list(bestModel), trait),
                              TD = createTD(modDat), traits = trait,
                              engine = engine, predicted = "trial")),
                    rownames(bestTab)[1])
  STA <- createSTA(models = model)
  res <- createVarCov(STA = STA, choice = rownames(bestTab)[1],
                      summary = bestTab, vcov = vcovBest,
                      criterion = criterion, engine = engine, trait = trait,
                      dat = TDTot)
  return(res)
}

#' Helper function for computing initial values.
#'
#' Replicates '_qvInitial' procedure (S. J. Welham 15/05/09) in GenStat
#'
#' @noRd
#' @keywords internal
initVals <- function(TD,
                     trait,
                     useWt = FALSE,
                     vcmodel = c("identity", "cs", "diagonal", "hcs",
                                 "outside", "fa", "fa2", "unstructured"),
                     fixed = NULL,
                     ...) {
  ## First, form estimate of unstructured matrix
  ## Remove the rows with NA.
  X <- na.omit(TD[, colnames(TD) %in% c(trait, "genotype", "trial")])
  X <- droplevels(X)
  nEnv <- nlevels(X[["trial"]])
  nGeno <- nlevels(X[["genotype"]])
  ## Get fixed df by stealth - in absence of other info, share among trials.
  if (!is.null(fixed)) {
    mr <- asreml::asreml(fixed = fixed, data = X, trace = FALSE, ...)
    P <- length(fitted(mr)) - mr$nedf
    fixedForm <- fixed
  } else {
    P <- 1
    fixedForm <- formula(paste0("`", trait, "`~ trial"))
  }
  ## Get number of effects contributing to each sum of squares.
  nobsEnv <- rowSums(table(X[, c("trial", "genotype")]))
  Rnobs <- matrix(data = nobsEnv, nrow = nEnv, ncol = nEnv, byrow = TRUE)
  Nobs <- ifelse(Rnobs < t(Rnobs), Rnobs, t(Rnobs))
  if (useWt) {
    ## This case is trickier, because of partitioning between two random terms,
    ## but use of the diag structure is better than nothing!
    initValues <- asreml::asreml(fixed = fixedForm,
                                 random = ~ genotype:idh(trial),
                                 weights = "wt", start.values = TRUE,
                                 data = X, trace = FALSE, ...)
    tmpTable <- initValues[[ifelse(asreml4(), "vparameters.table",
                                   "gammas.table")]]
    tmpTable[, "Constraint"] <- as.character(tmpTable[, "Constraint"])
    ## Fix residual variance.
    ## Column and component names are different for asreml3/asreml4.
    if (asreml4()) {
      tmpTable[tmpTable[["Component"]] == "units!R", "Constraint"] <- "F"
    } else {
      tmpTable[tmpTable[["Gamma"]] == "R!variance", "Constraint"] <- "F"
    }
    tmpTable[, "Constraint"] <- as.factor(tmpTable[, "Constraint"])
    mr <- asreml::asreml(fixed = fixedForm, random = ~ genotype:idh(trial),
                         weights = "wt", R.param = tmpTable, data = X,
                         trace = FALSE, ...)
    RMat <- matrix(data = coefficients(mr)$random, nrow = nGeno,
                   ncol = nEnv, byrow = TRUE)
  } else {
    ## This gives correct answers for complete balanced data.
    mr <- asreml::asreml(fixed = fixedForm, data = X, trace = FALSE, ...)
    mr$call$fixed <- fixedForm
    mr$call$data <- X
    res <- residuals(mr, type = "response")
    RMat <- tapply(X = res, INDEX = list(X[["genotype"]], X[["trial"]]),
                   FUN = mean)
    RMat[is.na(RMat)] <- 0
  }
  evCov <- crossprod(RMat) / (Nobs - (P / nEnv))
  ## Get off-diagonal elements of evCov only.
  offDiag <- evCov
  diag(offDiag) <- NA
  ## Get off-diagonal elements of cor(evCov) only.
  offCorr <- cor(evCov)
  diag(offCorr) <- NA
  diagEv <- diag(evCov)
  ## Get initial values for each model.
  if (vcmodel == "identity") {
    vcInitial <- list(vge = mean(diagEv))
  } else if (vcmodel == "cs") {
    vg <- mean(offDiag, na.rm = TRUE)
    vge <- mean(diagEv)
    vge <- ifelse(vge > vg , vge - vg, 0.1 * vge)
    vcInitial <- list(vg = vg, vge = vge)
  } else if (vcmodel == "diagonal") {
    vcInitial <- list(diag = diagEv)
  } else if (vcmodel == "hcs") {
    vg <- mean(offDiag, na.rm = TRUE) / 2
    diag <- ifelse(diagEv > vg, diagEv - vg, 0.1 * diagEv)
    vcInitial <- list(vg = vg, diag = diag)
  } else if (vcmodel == "outside") {
    vg <- mean(offCorr, na.rm = TRUE)
    vcInitial <- list(vg = vg, diag = diagEv)
  } else if (vcmodel == "fa") {
    factorAnalysis <- factanal(factors = 1, covmat = evCov)
    psi <- diagEv * factorAnalysis$uniquenesses
    SMat <- tcrossprod(factorAnalysis$loadings)
    chol <- suppressWarnings(chol(SMat, pivot = TRUE))
    gamma <- chol[1, ]
    vcInitial <- list(gamma = gamma, psi = psi)
  } else if (vcmodel == "fa2") {
    factorAnalysis <- factanal(factors = 2, covmat = evCov)
    psi <- diagEv * factorAnalysis$uniquenesses
    SMat <- tcrossprod(factorAnalysis$loadings[, 1])
    chol <- suppressWarnings(chol(SMat, pivot = TRUE))
    gamma <- chol[1:2, ]
    vcInitial <- list(gamma = gamma, psi = psi)
  } else if (vcmodel == "unstructured") {
    vcInitial <- list(evCov = evCov)
  }
  return(vcInitial)
}

