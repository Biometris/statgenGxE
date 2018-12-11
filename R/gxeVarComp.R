#' Selects the best variance-covariance model for a set of trials
#'
#' This function selects the best covariance structure for genetic correlations
#' between trials. It fits a range of variance-covariance models to
#' compare (identity, compound symmetry, diagonal, heterogeneous compound
#' symmetry, first order factor analysis, second order factor analysis,
#' unstructured), and selects the best one using a goodness-of-fit criterion.
#'
#' @inheritParams gxeAmmi
#'
#' @param engine A character string specifying the engine used for modeling.
#' Either "lme4" or "asreml".
#' @param criterion A string specifying a goodness-of-fit criterion. Either
#' "AIC" or "BIC".
#' @param ... Further arguments to be passed to \code{asreml}.
#'
#' @note If \code{engine = "lme4"}, only the compound symmetry model can be
#' fitted.
#'
#' @return An object of class \code{\link{varComp}}, a list object containing:
#' \item{SSA}{An object of class SSA containing the best fitted model.}
#' \item{choice}{A character string indicating the best fitted model.}
#' \item{summary}{A data.frame with a summary of the fitted models.}
#' \item{vcov}{The covariance matrix of the best fitted model.}
#' \item{criterion}{A character string indicating the goodness-of-fit criterion
#' used for determinening the best model, either "AIC" or "BIC".}
#' \item{engine}{A character string containing the engine used for
#' the analysis.}
#'
#' @examples
#' ## Select the best variance-covariance model using lme4 for modeling.
#' geVarComp <- gxeVarComp(TD = TDMaize, trait = "yld")
#' ## Summarize results.
#' summary(geVarComp)
#' \dontrun{
#' ## Create a pdf report summarizing the results.
#' report(geVarComp, outfile = "./testReports/reportVarComp.pdf")
#' }
#'
#' \dontrun{
#' ## Select the best variance-covariance model using asreml for modeling.
#' ## Use BIC as a goodness-of-fit criterion.
#' geVarComp2 <- gxeVarComp(TD = TDMaize, trait = "yld", engine = "asreml",
#'                         criterion = "BIC")
#' summary(geVarComp2)
#' ## Create a heatmap of the correlation matrix for the best model.
#' plot(geVarComp2)
#' }
#' @export
gxeVarComp <- function(TD,
                       trials = names(TD),
                       trait,
                       engine = c("lme4", "asreml"),
                       criterion = c("BIC", "AIC"),
                       useWt = FALSE,
                       ...) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (!is.character(trials) || !all(trials %in% names(TD))) {
    stop("All trials should be in TD.")
  }
  TDTot <- Reduce(f = rbind, x = TD[trials])
  TDTot$trial <- droplevels(TDTot$trial)
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !hasName(x = TDTot, name = trait)) {
    stop("trait has to be a column in TD.\n")
  }
  if (useWt && !hasName(x = TDTot, name = "wt")) {
    stop("wt has to be a column in TD when using weighting.")
  }
  engine <- match.arg(engine)
  criterion <- match.arg(criterion)
  ## Increase maximum number of iterations for asreml. Needed for more complex
  ## designs to converge.
  maxIter <- 200
  ## Add combinations of trial and genotype currently not in TD to TD.
  ## No missing combinations are allowed when fitting asreml models.
  TD0 <- expand.grid(genotype = levels(TDTot$genotype),
                     trial = levels(TDTot$trial))
  TDTot <- merge(TD0, TDTot, all.x = TRUE)
  if (engine == "asreml") {
    choices <- c("identity", "cs", "diagonal", "hcs", "outside",
                 "fa", "fa2", "unstructured")
  } else {
    choices <- "cs"
  }
  ## Create matrix for storing results.
  bestTab <- matrix(nrow = length(choices), ncol = 4,
                    dimnames = list(choices,
                                    c("AIC", "BIC", "Deviance", "NParameters")))
  models <- setNames(vector(mode = "list", length = length(choices)), choices)
  ## Main procedure to fit mixed models.
  if (engine == "lme4") {
    ## Compound symmetry ("cs") only.
    mr = lme4::lmer(formula(paste(trait, "~ trial + (1 | genotype)")),
                    data = TDTot, ...)
    nPar <- 2
    ## Construct SSA object.
    models[["cs"]] <- mr
    bestTab["cs", ] <- c(-2 * as.numeric(logLik(mr)) + 2 * nPar,
                         -2 * as.numeric(logLik(mr)) +
                           log(length(fitted(mr))) * nPar,
                         -2 * logLik(mr), nPar)
    bestModel <- models[[rownames(bestTab)[1]]]
    vcovBest <- vcov(emmeans::emmeans(bestModel, specs = "trial",
                                      lmer.df = "asymptotic"))
    colnames(vcovBest) <- rownames(vcovBest) <- levels(TDTot$trial)
  } else if (engine == "asreml") {
    if (requireNamespace("asreml", quietly = TRUE)) {
      nTr <- nlevels(TDTot$trial)
      fixedForm <- formula(paste(trait, "~ trial"))
      tmp <- tempfile()
      for (choice in choices) {
        if (choice == "identity") {
          sink(file = tmp)
          mr <- tryCatchExt(asreml::asreml(fixed = fixedForm,
                                           rcov = ~ genotype:trial,
                                           data = TDTot, maxiter = maxIter,
                                           ...))
          sink()
          nPar <- 1
        } else if (choice == "cs") {
          sink(file = tmp)
          mr <- tryCatchExt(asreml::asreml(fixed = fixedForm,
                                           random = ~ genotype,
                                           rcov = ~ genotype:trial,
                                           data = TDTot, maxiter = maxIter,
                                           ...))
          sink()
          nPar <- 2
        } else if (choice == "diagonal") {
          sink(file = tmp)
          mr <- tryCatchExt(asreml::asreml(fixed = fixedForm,
                                           rcov = ~ genotype:diag(trial),
                                           data = TDTot, maxiter = maxIter,
                                           ...))
          sink()
          nPar <- nTr
        } else if (choice == "hcs") {
          sink(file = tmp)
          mr <- tryCatchExt(asreml::asreml(fixed = fixedForm,
                                           random = ~ genotype,
                                           rcov = ~ genotype:diag(trial),
                                           data = TDTot, maxiter = maxIter,
                                           ...))
          sink()
          nPar <- nTr + 1
        } else if (choice == "outside") {
          sink(file = tmp)
          initVals <- asreml::asreml(fixed = fixedForm,
                                     random = ~ genotype:corh(trial),
                                     start.values = TRUE, data = TDTot, ...)
          sink()
          tmpValues <- initVals(TD = TDTot, trait = trait, useWt = useWt,
                                vcmodel = "outside", fixed = fixedForm, ...)
          tmpTable <- initVals$gammas.table
          tmpTable[, "Value"] <- c(tmpValues$vg,
                                   scale(tmpValues$diag, center = FALSE), 1)
          sink(file = tmp)
          if (useWt) {
            mr <- tryCatchExt(
              asreml::asreml(fixed = fixedForm,
                             random = ~ genotype:corh(trial),
                             G.param = tmpTable,
                             R.param = tmpTable, data = TDTot,
                             maxiter = maxIter, workspace = 160e6,
                             weigths = "wt",
                             family = asreml::asreml.gaussian(dispersion = 1),
                             ...))
          } else {
            mr <- tryCatchExt(
              asreml::asreml(fixed = fixedForm,
                             random = ~ genotype:corh(trial),
                             G.param = tmpTable,
                             R.param = tmpTable, data = TDTot,
                             maxiter = maxIter, workspace = 160e6,
                             ...))
          }
          sink()
          nPar <- nTr + 1
        } else if (choice == "fa" && nTr > 4) {
          sink(file = tmp)
          initVals <- asreml::asreml(fixed = fixedForm,
                                     random = ~ genotype:fa(trial, 1),
                                     start.values = TRUE, data = TDTot, ...)
          sink()
          tmpTable <- initVals$gammas.table
          tmpValues <- initVals(TD = TDTot, trait = trait, useWt = useWt,
                                vcmodel = "fa", fixed = fixedForm, ...)
          if (!is.null(tmpValues)) {
            tmpTable[, "Value"] <- c(scale(tmpValues$psi, center = FALSE),
                                     tmpValues$gamma, 1)
          }
          sink(file = tmp)
          if (useWt) {
            mr <- tryCatchExt(
              asreml::asreml(fixed = fixedForm,
                             random = ~ genotype:fa(trial, 1),
                             R.param = tmpTable,
                             G.param = tmpTable, data = TDTot,
                             maxiter = maxIter,
                             workspace = 160e6,
                             weigths = "wt",
                             family = asreml::asreml.gaussian(dispersion = 1),
                             ...))

          } else {
            mr <- tryCatchExt(
              asreml::asreml(fixed = fixedForm,
                             random = ~ genotype:fa(trial, 1),
                             R.param = tmpTable,
                             G.param = tmpTable, data = TDTot,
                             maxiter = maxIter,
                             workspace = 160e6, ...))
          }
          sink()
          nPar <- nTr * 2
        } else if (choice == "fa2" && nTr > 4) {
          sink(file = tmp)
          initVals <- asreml::asreml(fixed = fixedForm,
                                     random = ~ genotype:fa(trial, 2),
                                     start.values = TRUE, data = TDTot, ...)
          sink()
          tmpTable <- initVals$gammas.table
          tmpValues <- initVals(TD = TDTot, trait = trait, useWt = useWt,
                                vcmodel = "fa2", fixed = fixedForm, ...)
          if (!is.null(tmpValues)) {
            ## Keep loadings of factor 2 away from 0.
            tmpValues$gamma[2, tmpValues$gamma[2, ] < 1e-3] <- 1e-3
            ## Make sure that first entry is 0.
            tmpValues$gamma[2, 1] <- 0
            tmpTable[, "Value"] <- c(scale(tmpValues$psi, center = FALSE),
                                     tmpValues$gamma[1, ],
                                     tmpValues$gamma[2, ], 1)
          }
          sink(file = tmp)
          if (useWt) {
            mr <- tryCatchExt(
              asreml::asreml(fixed = fixedForm,
                             random = ~ genotype:fa(trial, 2),
                             R.param = tmpTable,
                             G.param = tmpTable, data = TDTot,
                             maxiter = maxIter,
                             workspace = 160e6,
                             weigths = "wt",
                             family = asreml::asreml.gaussian(dispersion = 1),
                             ...))

          } else {
            mr <- tryCatchExt(
              asreml::asreml(fixed = fixedForm,
                             random = ~ genotype:fa(trial, 2),
                             R.param = tmpTable,
                             G.param = tmpTable, data = TDTot,
                             maxiter = maxIter,
                             workspace = 160e6, ...))
          }
          sink()
          nPar <- nTr * 3 - 1
        } else if (choice == "unstructured") {
          ## Check model.
          sink(file = tmp)
          initVals <- asreml::asreml(fixed = fixedForm,
                                     random = ~ genotype:us(trial),
                                     start.values = TRUE, data = TDTot, ...)
          sink()
          tmpValues <- initVals(TD = TDTot, trait = trait, useWt = useWt,
                                vcmodel = "unstructured", fixed = fixedForm,
                                ...)
          tmpValues <- tmpValues$evCov[upper.tri(tmpValues$evCov, diag = TRUE)]
          tmpTable <- initVals$gammas.table
          tmpTable[, "Value"] <- c(scale(tmpValues, center = FALSE),
                                   10 ^ min(floor(log10(min(abs(tmpValues)))), 0))
          ## All off diagonal elements are unconstrained, diagonal elements
          ## should be positive
          tmpTable[-c((1:nTr) * ((1:nTr) + 1) / 2), "Constraint"] <- "U"
          ## Fix residual variance.
          tmpTable[tmpTable$Gamma == "R!variance", "Constraint"] <- "F"
          sink(file = tmp)
          if (useWt) {
            mr <- tryCatchExt(
              asreml::asreml(fixed = fixedForm,
                             random = ~ genotype:us(trial),
                             G.param = tmpTable,
                             R.param = tmpTable, data = TDTot,
                             maxiter = maxIter, workspace = 160e6,
                             weigths = "wt",
                             family = asreml::asreml.gaussian(dispersion = 1),
                             ...))
          } else {
            mr <- tryCatchExt(
              asreml::asreml(fixed = fixedForm,
                             random = ~ genotype:us(trial),
                             G.param = tmpTable,
                             R.param = tmpTable, data = TDTot,
                             maxiter = maxIter, workspace = 160e6,
                             ...))
          }
          sink()
          nPar <- nTr * (nTr + 1) / 2
        }
        if (!is.null(mr$warning)) {
          ## Check if param 1% increase is significant. Remove warning if not.
          mr <- chkLastIter(mr)
        }
        if (length(mr$warning) != 0) {
          warning(paste0("Asreml gave the following warning for ", choice,
                         ":\n", mr$warning, "\n"), call. = FALSE)
        }
        if (!is.null(mr$error)) {
          warning(paste0("Asreml gave the following error for ", choice,
                         ":\n", mr$error), call. = FALSE)
          mr <- list(loglik = -Inf)
        } else if (!(nTr <= 4 && choice %in% c("fa", "fa2"))) {
          mr <- mr$value
          mr$call$fixed <- eval(mr$call$fixed)
          mr$call$random <- eval(mr$call$random)
          mr$call$rcov <- eval(mr$call$rcov)
          mr$call$G.param <- eval(mr$call$G.param)
          mr$call$R.param <- eval(mr$call$R.param)
          if (!mr$converge) {
            warning(paste0("No convergence for ", choice, ".\n"), call. = FALSE)
            mr$loglik <- -Inf
          }
          models[[choice]] <- mr
          bestTab[choice, "AIC"] <- -2 * mr$loglik + 2 * nPar
          bestTab[choice, "BIC"] <- -2 * mr$loglik +
            (log(length(mr$fitted.values) - nlevels(TDTot$genotype)) * nPar)
          bestTab[choice, "Deviance"] <- -2 * mr$loglik
          bestTab[choice, "NParameters"] <- nPar
        }
      } # End loop over choices
      bestTab <- bestTab[order(bestTab[, criterion]), ]
      bestModel <- models[[rownames(bestTab)[1]]]
      bestModel <- predictAsreml(model = bestModel, classify = "trial",
                                 TD = TDTot, maxiter = maxIter, ...)
      vcovBest <- bestModel$predictions$vcov
      colnames(vcovBest) <- rownames(vcovBest) <- levels(TDTot$trial)
      unlink(tmp)
    } else {
      stop("Failed to load 'asreml'.\n")
    }
  }
  TDTot$trial <- rownames(bestTab)[1]
  ## Create output.
  model <- setNames(list(list(mRand = NULL,
                              mFix = setNames(list(bestModel), trait),
                              TD = createTD(TDTot), traits = trait,
                              engine = engine, predicted = "trial")),
                    rownames(bestTab)[1])
  SSA <- createSSA(models = model)
  res <- createVarComp(SSA = SSA, choice = rownames(bestTab)[1],
                       summary = bestTab, vcov = vcovBest,
                       criterion = criterion, engine = engine)
  return(res)
}

#' Helper function for computing initial values.
#'
#' Replicates '_qvInitial' procedure (S. J. Welham 15/05/09) in GenStat
#'
#' @keywords internal
initVals <- function(TD,
                     trait,
                     useWt = FALSE,
                     vcmodel = c("identity", "cs", "diagonal", "hcs",
                                 "outside", "fa", "fa2", "unstructured"),
                     fixed = NULL,
                     ...) {
  ## First, form estimate of unstructured matrix
  ## Create tempfile for asreml output.
  tmp <- tempfile()
  ## Remove the rows with NA.
  X <- na.omit(TD[, colnames(TD) %in% c(trait, "genotype", "trial", "wt")])
  X <- droplevels(X)
  nEnv <- nlevels(X$trial)
  nGeno <- nlevels(X$genotype)
  ## Get fixed df by stealth - in absence of other info, share among trials.
  if (!is.null(fixed)) {
    sink(file = tmp)
    mr <- asreml::asreml(fixed = fixed, rcov = ~id(units), data = X, ...)
    sink()
    P <- length(mr$fitted.values) - (1 + mr$nedf)
    fixedForm <- fixed
  } else {
    P <- 1
    fixedForm <- formula(paste(trait, "~ trial"))
  }
  ## Get number of effects contributing to each sum of squares.
  nobsEnv <- rowSums(table(X[, c("trial", "genotype")]))
  Rnobs <- matrix(data = nobsEnv, nrow = nEnv, ncol = nEnv, byrow = TRUE)
  Nobs <- ifelse(Rnobs < t(Rnobs), Rnobs, t(Rnobs))
  if (useWt) {
    ## This case is trickier, because of partitioning between two random terms,
    ## but use of the diag structure is better than nothing!
    sink(file = tmp)
    initValues <- asreml::asreml(fixed = fixedForm,
                                 random = ~ genotype:idh(trial),
                                 weights = "wt", start.values = TRUE,
                                 data = X, ...)
    sink()
    tmpTable <- initValues$gammas.table
    tmpTable[, "Constraint"] <- as.character(tmpTable[,"Constraint"])
    tmpTable[tmpTable$Gamma == "R!variance", "Constraint"] <- "F"
    tmpTable[, "Constraint"] <- as.factor(tmpTable[, "Constraint"])
    sink(file = tmp)
    mr <- asreml::asreml(fixed = fixedForm, random = ~ genotype:idh(trial),
                         weights = "wt", R.param = tmpTable, data = X, ...)
    sink()
    RMat <- matrix(data = coefficients(mr)$random, nrow = nGeno,
                   ncol = nEnv, byrow = TRUE)
  } else {
    ## This gives correct answers for complete balanced data.
    sink(file = tmp)
    mr <- asreml::asreml(fixed = fixedForm, data = X, ...)
    sink()
    mr$call$fixed <- fixedForm
    mr$call$data <- X
    res <- residuals(mr, type = "response")
    RMat <- tapply(X = res, INDEX = list(X$genotype, X$trial), FUN = mean)
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
  unlink(tmp)
  return(vcInitial)
}

