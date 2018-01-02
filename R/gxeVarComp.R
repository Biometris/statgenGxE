#' Selects the best variance-covariance model for a set of environments
#'
#' This function selects the best covariance structure for genetic correlations
#' between environments. It fits a range of variance-covariance models to
#' compare (identity, compound symmetry, diagonal, heterogeneous compound
#' symmetry, first order factor analysis, second order factor analysis,
#' unstructured), and selects the best one using a goodness-of-fit criterion.
#'
#' @inheritParams gxeAmmi
#'
#' @param engine A string specifying the name of a mixed modelling engine to be used.
#' @param criterion A string specifying a goodness-of-fit criterion, i.e., "AIC" or "BIC".
#' @param ... Further arguments to be passed to \code{asreml}.
#'
#' @note If \code{engine = "lme4"}, only the compound symmetry model can be fitted.
#'
#' @return A list object consisting of the fitted model objects, a string specifying
#' the best model and its related goodness-of-fit criterion.
#'
#' @examples
#' ## Fit models.
#' geVarComp1 <- gxeVarComp(TD = TDMaize, trait = "yld", engine = "lme4",
#'                         criterion = "BIC")
#' ## Display results.
#' summary(geVarComp1)
#'
#' \dontrun{
#' geVarComp2 <- gxeVarComp(TD = TDMaize, trait = "yld", engine = "asreml",
#'                         criterion = "BIC")
#' report(geVarComp2, outfile = "./testReports/reportVarComp.pdf")
#' }
#'
#' @export

gxeVarComp <- function(TD,
                       trait,
                       engine = "asreml",
                       criterion = "BIC",
                       ...) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TD)) {
    stop("trait has to be a column in TD.\n")
  }
  if (!is.null(engine) && (!is.character(engine) || length(engine) > 1 ||
                           !engine %in% c("asreml", "lme4"))) {
    stop("engine should be asreml or lme4.\n")
  }
  if (!is.null(criterion) && (!is.character(criterion) || length(criterion) > 1 ||
                              !criterion %in% c("AIC", "BIC"))) {
    stop("criterion should be AIC or BIC.\n")
  }
  ## Add combinations of env and genotype currently not in TD to TD.
  TD <- reshape2::melt(data = reshape2::dcast(data = TD,
                                              formula = env ~ genotype,
                                              value.var = trait),
                       id.vars = "env", variable.name = "genotype",
                       value.name = trait)
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
    mr = lme4::lmer(as.formula(paste(trait, "~ env + (1 | genotype)")),
                    data = TD, ...)
    nPar <- 2
    ## Construct SSA object.
    models[["cs"]] <- mr
    bestTab["cs", ] <- c(-2 * as.numeric(logLik(mr)) + 2 * nPar,
                         -2 * as.numeric(logLik(mr)) +
                           log(length(fitted(mr))) * nPar,
                         -2 * logLik(mr), nPar)
    bestModel <- models[[rownames(bestTab)[1]]]
    vcovBest <- vcov(emmeans::emmeans(bestModel, specs = "env",
                                      lmer.df = "asymptotic"))
    colnames(vcovBest) <- rownames(vcovBest) <- levels(TD$env)
  } else if (engine == "asreml") {
    if (requireNamespace("asreml", quietly = TRUE)) {
      tmp <- tempfile()
      for (choice in choices) {
        if (choice == "identity") {
          sink(file = tmp)
          mr <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                               rcov = ~ genotype:env,
                               data = TD, ...)
          sink()
          nPar <- 1
        } else if (choice == "cs") {
          sink(file = tmp)
          mr <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                               random = ~ genotype,
                               rcov = ~ genotype:env,
                               data = TD, ...)
          sink()
          nPar <- 2
        } else if (choice == "diagonal") {
          sink(file = tmp)
          mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                   rcov = as.formula("~ genotype:diag(env)"),
                                   data = TD, ...), silent = TRUE)
          sink()
          nPar <- nlevels(TD$env)
        } else if (choice == "hcs") {
          sink(file = tmp)
          mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                   random = as.formula("~ genotype"),
                                   rcov = as.formula("~ genotype:diag(env)"),
                                   data = TD, ...), silent = TRUE)
          sink()
          nPar <- nlevels(TD$env) + 1
        } else if (choice == "outside") {
          sink(file = tmp)
          initVals <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                     random = as.formula("~ genotype:corh(env)"),
                                     start.values = TRUE, data = TD, ...)
          sink()
          unlink(tmp)
          tmpValues <- qvInitial(TD = TD, trait = trait, unitError = NA,
                                 vcmodel = "outside",
                                 fixed = as.formula(paste(trait, "~ env")),
                                 unitFactor = NA, ...)
          tmpTable <- initVals$gammas.table
          tmpTable[, "Value"] <- c(tmpValues$vg, tmpValues$diag, 1)
          tmpTable[, "Constraint"] <- as.character(tmpTable[, "Constraint"] )
          tmpTable[which(tmpTable[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
          tmpTable[, "Constraint"] <- as.factor(tmpTable[, "Constraint"])
          sink(file = tmp)
          mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                   random = as.formula("~ genotype:corh(env)"),
                                   G.param = tmpTable, R.param = tmpTable, data = TD, ...),
                    silent = TRUE)
          sink()
          nPar <- nlevels(TD$env) + 1
        } else if (choice == "fa" && nlevels(TD$env) > 4) {
          sink(file = tmp)
          initVals <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                     random = as.formula("~ genotype:fa(env, 1)"),
                                     start.values = TRUE, data = TD, ...)
          sink()
          tmpTable <- initVals$gammas.table
          tmpValues <- qvInitial(TD = TD, trait = trait, unitError = NA,
                                 vcmodel = "fa", fixed = as.formula(paste(trait, "~ env")),
                                 unitFactor = NA, ...)
          if (is.null(tmpValues)) {
            sink(file = tmp)
            mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                     random = as.formula("~ genotype:fa(env, 1)"),
                                     data = TD, ...), silent = TRUE)
            sink()
          } else {
            tmpTable[, "Value"] <- c(tmpValues$psi, tmpValues$gamma, 1)
            tmpTable[, "Constraint"] <- as.character(tmpTable[, "Constraint"] )
            tmpTable[which(tmpTable[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
            tmpTable[, "Constraint"] <- as.factor(tmpTable[, "Constraint"])
            sink(file = tmp)
            mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                     random = as.formula("~ genotype:fa(env, 1)"),
                                     R.param = tmpTable, G.param = tmpTable,
                                     data = TD, ...),
                      silent = TRUE)
            sink()
          }
          nPar <- nlevels(TD$env) * 2
        } else if (choice == "fa2" && nlevels(TD$env) > 4) {
          sink(file = tmp)
          initVals <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                     random = as.formula("~ genotype:fa(env, 2)"),
                                     start.values = TRUE, data = TD, ...)
          sink()
          tmpTable <- initVals$gammas.table
          tmpValues <- qvInitial(TD = TD, trait = trait, unitError = NA,
                                 vcmodel = "fa2",
                                 fixed = as.formula(paste(trait, "~ env")),
                                 unitFactor = NA, ...)
          if (is.null(tmpValues)) {
            sink(file = tmp)
            mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                     random = as.formula("~ genotype:fa(env, 2)"),
                                     data = TD, ...), silent = TRUE)
            sink()
          } else {
            ## Keep loadings of factor 2 away from 0.
            tmpValues$gamma[2, which(tmpValues$gamma[2, ] < 1e-3)] <- 1e-3
            ## Make sure that first entry is 0.
            tmpValues$gamma[2, 1] <- 0
            tmpTable[, "Value"] <- c(tmpValues$psi, tmpValues$gamma[1, ],
                                     tmpValues$gamma[2, ], 1)
            tmpTable[, "Constraint"] <- as.character(tmpTable[, "Constraint"])
            tmpTable[which(tmpTable[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
            tmpTable[, "Constraint"] <- as.factor(tmpTable[, "Constraint"])
            sink(file = tmp)
            mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                     random = as.formula("~ genotype:fa(env, 2)"),
                                     R.param = tmpTable, G.param = tmpTable,
                                     data = TD, ...),
                      silent = TRUE)
            sink()
          }
          nPar <- nlevels(TD$env) * 3 - 1
        } else if (choice == "unstructured") {
          ## Check model.
          sink(file = tmp)
          initVals <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                     random = as.formula("~ genotype:us(env)"),
                                     start.values = TRUE, data = TD, ...)
          sink()
          tmpValues <- qvInitial(TD = TD, trait = trait, vcmodel = "unstructured",
                                 fixed = as.formula(paste(trait, "~ env")),
                                 unitFactor = NA, ...)
          tmpValues <- tmpValues$evCov[upper.tri(tmpValues$evCov, diag = TRUE)]
          tmpTable <- initVals$gammas.table
          tmpTable[, "Value"] <- c(tmpValues, 1)
          tmpTable[, "Constraint"] <- as.character(tmpTable[, "Constraint"] )
          tmpTable[which(tmpTable[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
          tmpTable[, "Constraint"] <- as.factor(tmpTable[, "Constraint"])
          sink(file = tmp)
          mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                   random = as.formula("~ genotype:us(env)"),
                                   G.param = tmpTable, R.param = tmpTable,
                                   data = TD, ...),
                    silent = TRUE)
          sink()
          nPar <- nlevels(TD$env) * (nlevels(TD$env) - 1) / 2 +
            nlevels(TD$env)
        }
        if (inherits(mr, "try-error")) {
          mr <- list(loglik = -Inf)
        } else {
          mr$call$fixed <- eval(mr$call$fixed)
          mr$call$random <- eval(mr$call$random)
          mr$call$rcov <- eval(mr$call$rcov)
          mr$call$G.param <- eval(mr$call$G.param)
          mr$call$R.param <- eval(mr$call$R.param)
          models[[choice]] <- mr
          if (!mr$converge) {
            mr$loglik <- -Inf
          }
        }
        if (!(nlevels(TD$env) <= 4 && choice %in% c("fa", "fa2"))) {
          bestTab[choice, "AIC"] <- -2 * mr$loglik + 2 * nPar
          bestTab[choice, "BIC"] <- -2 * mr$loglik +
            log(length(mr$fitted.values)) * nPar
          bestTab[choice, "Deviance"] <- -2 * mr$loglik
          bestTab[choice, "NParameters"] <- nPar
        }
      }
      bestModel <- models[[rownames(bestTab)[1]]]
      bestModel <- predictAsreml(model = bestModel,
                                 classify = "env",
                                 TD = TD)
      bestTab <- bestTab[order(bestTab[, criterion]), ]
      vcovBest <- bestModel$predictions$vcov
      colnames(vcovBest) <- rownames(vcovBest) <- levels(TD$env)
      unlink(tmp)
    } else {
      stop("Failed to load 'asreml'.\n")
    }
  }
  ## Create output.
  model <- createSSA(mRand = NULL, mFix = setNames(list(bestModel), trait),
                     data = TD, traits = trait,
                     engine = engine, predicted = "env")
  res <- createVarComp(model = model, choice = rownames(bestTab)[1],
                       summary = bestTab, vcov = vcovBest,
                       criterion = criterion)
  return(res)
}

#' Helper function
#'
#' Replicates '_qvInitial' procedure (S. J. Welham 15/05/09) in GenStat
#'
#' @keywords internal
qvInitial <- function(TD,
                      trait,
                      unitError = NA,
                      vcmodel = c("identity", "cs", "diagonal", "hcs",
                                  "outside", "fa", "fa2", "unstructured"),
                      fixed = NULL,
                      unitFactor = NA, ...) {

  ## First, form estimate of unstructured matrix
  ## Create tempfile for asreml output.
  tmp <- tempfile()
  ## Remove the rows with NA.
  X <- na.omit(TD[, c(trait, "genotype", "env")])
  nEnv <- nlevels(X$env)
  nGeno <- nlevels(X$genotype)
  ## Get fixed df by stealth - in absence of other info, share among environments
  if (!is.null(fixed)) {
    sink(file = tmp)
    mr <- asreml::asreml(fixed = fixed, rcov = ~id(units), data = X, ...)
    sink()
    P <- length(mr$fitted.values) - (1 + mr$nedf)
    fixedForm <- fixed
  } else {
    P <- 1
    fixedForm <- as.formula(paste(trait, "~ env"))
  }
  ## Get number of effects contributing to each sum of squares.
  nobsEnv <- rowSums(table(X[, c("env", "genotype")]))
  Rnobs <- matrix(data = nobsEnv, nrow = nEnv, ncol = nEnv, byrow = TRUE)
  Cnobs <- t(Rnobs)
  Nobs <- Rnobs * (Rnobs - Cnobs < 0) + Cnobs * (Cnobs - Rnobs <= 0)
  if (!is.na(unitError)) {
    ## This case is trickier, becuase of partitioning between two random terms,
    ## but use of the diag structure is better than nothing!
    weights <- 1 / unitError
    X["weights"] <- weights
    sink(file = tmp)
    initValues <- asreml::asreml(fixed = fixedForm,
                                 random = as.formula("~ genotype:idh(env)"),
                                 weights = weights, start.values = TRUE,
                                 data = X, ...)
    sink()
    tmpTable <- initValues$gammas.table
    tmpTable[, "Constraint"] <- as.character(tmpTable[,"Constraint"])
    tmpTable[which(tmpTable[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
    tmpTable[, "Constraint"] <- as.factor(tmpTable[, "Constraint"])
    sink(file = tmpTable)
    mr <- asreml::asreml(fixed = fixedForm,
                         random = as.formula("~ genotype:idh(env)"),
                         weights = weights, R.param = tmpTable, data = X, ...)
    sink()
    RMat <- matrix(data = coefficients(mr)$random, nrow = nGeno,
                   ncol = nEnv, byrow = TRUE)
  } else {
    ## This gives correct answers for complete balanced data.
    sink(file = tmp)
    mr <- asreml::asreml(fixed = fixedForm, data = X, ...)
    sink()
    mr$call$fixed <- eval(mr$call$fixed)
    mr$call$random <- eval(mr$call$random)
    mr$call$rcov <- eval(mr$call$rcov)
    mr$call$data <- substitute(TD)
    res <- residuals(mr, type = "response")
    RMat <- tapply(X = res, INDEX = list(X$genotype, X$env), FUN = mean)
    RMat[which(is.na(RMat))] <- 0
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
    vge <- (vge > vg) * (vge - vg) + (vge <= vg) * 0.1 * vge
    vcInitial <- list(vg = vg, vge = vge)
  } else if (vcmodel == "diagonal") {
    vcInitial <- list(diag = diagEv)
  } else if (vcmodel == "hcs") {
    vg <- mean(offDiag, na.rm = TRUE) / 2
    diag <- (diagEv > vg) * (diagEv - vg) + (diagEv <= vg) * 0.1 * diagEv
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

