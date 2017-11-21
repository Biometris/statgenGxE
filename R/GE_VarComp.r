#' Selects the best variance-covariance model for a set of environments
#'
#' This function selects the best covariance structure for genetic correlations between
#' environments. It fits a range of variance-covariance models to
#' compare (e.g., identity, compound symmetry, diagonal, heterogeneous compound symmetry,
#' first order factor analysis, second order factor analysis, unstructured),
#' and selects the best one using a goodness-of-fit criterion.
#'
#' @inheritParams GE.AMMI
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
#' myDat <- GE.read.csv(system.file("extdata", "F2maize_pheno.csv", package = "RAP"),
#'                      env = "env!", genotype = "genotype!", trait = "yld")
#' myTD <- createTD(data = myDat, genotype = "genotype!", env = "env!")
#' model1 <- GE.VarComp(TD = myTD, trait = "yld", engine = "lme4", #engine = "asreml",
#'                      criterion = "BIC")
#' model1$BIC
#' model1$choice
#' summary(model1$model[[model1$choice]])
#'
#' @export
GE.VarComp <- function(TD,
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
  ## keep NA values.
  TD <- droplevels(TD[, c("env", "genotype", trait)])
  naEntries <- which(is.na(tapply(X = rep(1, nrow(TD)),
                                  INDEX = TD[, c("genotype", "env")],
                                  FUN = identity)),
                     arr.ind = TRUE)
  if (length(naEntries)) {
    TDExt <- data.frame(TD[naEntries$env, "env"],
                        TD[naEntries$genotype, "genotype"], NA)
    names(TDExt) <- names(TD)
    TD <- rbind(TD, TDExt)
  }
  ## sort data in order of genotype, env.
  TD <- TD[order(TD$genotype, TD$env), ]
  qvInitial <- function(TD,
                        trait,
                        unitError = NA,
                        vcmodel = c("identity", "cs", "diagonal", "hcs",
                                    "outside", "fa", "fa2", "unstructured"),
                        fixed = NULL,
                        unitFactor = NA, ...) {
    # Replicates '_qvInitial' procedure (S. J. Welham 15/05/09) in GenStat
    # TODO: factanal() ? fa()
    # First, form estimate of unstructured matrix
    # exclude the rows including NA
    X <- na.omit(TD[, c(trait, "genotype", "env")])
    nEnv <- nlevels(X$env)
    nGeno <- nlevels(X$genotype)
    ## Get fixed df by stealth - in absence of other info, share among environments
    if (!is.null(fixed)) {
      mr <- asreml::asreml(fixed = fixed, rcov = ~id(units), data = X, ...)
      P <- length(mr$fitted.values) - (1 + mr$nedf)
    } else {
      P <- 1
    }
    ## Get number of effects contributing to each sum of squares.
    tmpTable <- table(X[, c("env", "genotype")])
    nobsEnv <- rowSums(tmpTable)
    Rnobs <- matrix(data = nobsEnv, nrow = nEnv, ncol = nEnv, byrow = TRUE)
    Cnobs <- t(Rnobs)
    Nobs <- Rnobs * (Rnobs - Cnobs < 0) + Cnobs * (Cnobs - Rnobs <= 0)
    if (!is.na(unitError)) {
      ## This case is trickier, becuase of partitioning between two random terms,
      ## but use of the diag structure is better than nothing!
      weights <- 1 / unitError
      X["weights"] <- weights
      if (!is.null(fixed)) {
        initValues <- asreml::asreml(fixed = fixed,
                                     random = as.formula("~ genotype:idh(env)"),
                                     weights = weights, start.values = TRUE,
                                     data = X, ...)
        tmp <- initValues$gammas.table
        tmp[, "Constraint"] <- as.character(tmp[,"Constraint"])
        tmp[which(tmp[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
        tmp[, "Constraint"] <- as.factor(tmp[, "Constraint"])
        mr <- asreml::asreml(fixed = fixed,
                             random = as.formula("~ genotype:idh(env)"),
                             weights = weights, R.param = tmp, data = X, ...)
        Ores <- matrix(data = mr$coeff$random, nrow = nGeno, ncol = nEnv, byrow = TRUE)
      } else {
        initValues <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                     random = as.formula("~ genotype:idh(env)"),
                                     weights = weights, start.values = TRUE,
                                     data = X, ...)
        tmp <- initValues$gammas.table
        tmp[, "Constraint"] <- as.character(tmp[, "Constraint"] )
        tmp[which(tmp[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
        tmp[, "Constraint"] <- as.factor(tmp[, "Constraint"])
        mr <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                             random = as.formula("~ genotype:idh(env)"),
                             weights = weights, R.param = tmp, data = X, ...)
        Ores <- matrix(data = mr$coeff$random, nrow = nGeno, ncol = nEnv, byrow = TRUE)
      }
    } else {
      ## This gives correct answers for complete balanced data.
      if (!is.null(fixed)) {
        mr <- asreml::asreml(fixed = fixed, data = X, ...)
        residuals <- mr$residuals
      } else {
        mr <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")), data = X, ...)
        residuals <- mr$residuals
      }
      Ores <- tapply(X = residuals, INDEX = list(X$genotype, X$env), FUN = mean)
      if (sum(is.na(Ores)) != 0) {
        Ores[which(is.na(Ores))] <- 0
      }
    }
    Rmat <- Ores
    evCov <- crossprod(Rmat) / (Nobs - (P / nEnv))
    ## Get off-diagonal elements of evCov only.
    offdiag <- evCov
    diag(offdiag) <- NA
    ## Get off-diagonal elements of cor(evCov) only.
    offCorr <- cor(evCov)
    diag(offCorr) <- NA
    ## Get initial values for each model.
    if (vcmodel == "identity") {
      vcInitial <- list(vge = mean(diag(evCov)))
    } else if (vcmodel == "cs") {
      vg <- mean(offdiag, na.rm = TRUE)
      vge <- mean(diag(evCov))
      vge <- (vge > vg) * (vge - vg) + (vge <= vg) * 0.1 * vge
      vcInitial <- list(vg = vg, vge = vge)
    } else if (vcmodel == "diagonal") {
      vcInitial <- list(diag = diag(evCov))
    } else if (vcmodel == "hcs") {
      vg <- mean(offdiag, na.rm = TRUE) / 2
      diag <- diag(evCov)
      diag <- (diag > vg) * (diag - vg) + (diag <= vg) * 0.1 * diag
      vcInitial <- list(vg = vg, diag = diag)
    } else if (vcmodel == "outside") {
      vg <- mean(offCorr, na.rm = TRUE)
      diag <- diag(evCov)
      vcInitial <- list(vg = vg, diag = diag)
    } else if (vcmodel == "fa") {
      if (requireNamespace("psych", quietly = TRUE) &&
          requireNamespace("GPArotation", quietly = TRUE)) {
        factorAnalysis <- try(psych::fa(r = evCov, nfactors = 1, fm = "mle"), silent = TRUE)
        if (inherits(factorAnalysis, "try-error")) {
          factorAnalysis <- psych::fa(r = evCov, nfactors = 1, fm = "minres")
        }
        loading <- as.vector(factorAnalysis$loadings)
        comm <- factorAnalysis$communalities
        var <- diag(evCov)
        psi <- var * (1 - comm)
        Smat <- tcrossprod(loading)
        chol <- suppressWarnings(chol(Smat, pivot = TRUE))
        gamma <- chol[1, ]
        vcInitial <- list(gamma = gamma, psi = psi)
      } else {
        vcInitial <- NULL
        warning("psych package is required but failed to load.\n")
      }
    } else if (vcmodel == "fa2") {
      if (requireNamespace("psych", quietly = TRUE) &&
          requireNamespace("GPArotation", quietly = TRUE)) {
        factorAnalysis <- try(psych::fa(r = evCov, nfactors = 2, fm = "mle"),
                              silent = TRUE)
        if (inherits(factorAnalysis, "try-error")) {
          factorAnalysis <- psych::fa(r = evCov, nfactors = 2, fm = "minres")
        }
        loading <- as.vector(factorAnalysis$loadings[, 1])
        comm <- factorAnalysis$communalities
        var <- diag(evCov)
        psi <- var * (1 - comm)
        Smat <- tcrossprod(loading)
        chol <- suppressWarnings(chol(Smat, pivot = TRUE))
        gamma <- chol[1:2, ]
        vcInitial <- list(gamma = gamma, psi = psi)
      } else {
        vcInitial <- NULL
        warning("psych and GPArotation packages are required but failed to load.\n")
      }
    } else if (vcmodel == "unstructured") {
      vcInitial <- list(evCov = evCov)
    }
    return(vcInitial)
  }
  res <- vector(mode = "list")
  ## Main procedure to fit mixed models.
  if (engine == "lme4") {
    ## Compound symmetry ("cs") only.
    if (requireNamespace("lme4", quietly = TRUE)) {
      mr = lme4::lmer(as.formula(paste(trait, "~ env + (1 | genotype)")),
                      data = TD, ...)
      nPar <- 2
      ## Outputs.
      res$model$cs <- mr
      res$choice <- "cs"
      if (criterion == "AIC") {
        res$AICBest <- -2 * as.numeric(logLik(mr)) + 2 * nPar
      } else {
        res$BICBest <- -2 * as.numeric(logLik(mr)) + log(length(fitted(mr))) * nPar
      }
    } else {
      stop("Failed to load 'lme4'.\n")
    }
  } else if (engine == "asreml") {
    if (requireNamespace("asreml", quietly = TRUE)) {
      choices <- c("identity", "cs", "diagonal", "hcs", "outside",
                   "fa", "fa2", "unstructured")
      bestTab <- matrix(nrow = 8, ncol = 4)
      colnames(bestTab) <- c("AIC", "BIC", "Deviance", "NParameters")
      rownames(bestTab) <- choices
      for (i in 1:length(choices)) {
        if (choices[i] == "identity") {
          mr <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                               rcov = as.formula("~ genotype:env"), data = TD, ...)
          mr$call$fixed <- eval(mr$call$fixed)
          mr$call$rcov <- eval(mr$call$rcov)
          res$model[["identity"]] <- mr
          if (!mr$converge) {
            mr$loglik <- -Inf
          }
          nPar <- 1
        } else if (choices[i] == "cs") {
          mr <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                               random = as.formula("~ genotype"),
                               rcov = as.formula("~ genotype:env"), data = TD, ...)
          mr$call$fixed <- eval(mr$call$fixed)
          mr$call$random <- eval(mr$call$random)
          mr$call$rcov <- eval(mr$call$rcov)
          res$model[["cs"]] <- mr
          if (!mr$converge) {
            mr$loglik <- -Inf
          }
          nPar <- 2
        } else if (choices[i] == "diagonal") {
          mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                   rcov = as.formula("~ genotype:diag(env)"),
                                   data = TD, ...), silent = TRUE)
          if (inherits(mr, "try-error")) {
            mr <- list(loglik = -Inf)
          } else {
            mr$call$fixed <- eval(mr$call$fixed)
            mr$call$rcov <- eval(mr$call$rcov)
            res$model[["diagonal"]] <- mr
            if (!mr$converge) {
              mr$loglik <- -Inf
            }
          }
          nPar <- nlevels(TD$env)
        } else if (choices[i] == "hcs") {
          mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                   random = as.formula("~ genotype"),
                                   rcov = as.formula("~ genotype:diag(env)"),
                                   data = TD, ...), silent = TRUE)
          if (inherits(mr, "try-error")) {
            mr <- list(loglik = -Inf)
          } else {
            mr$call$fixed <- eval(mr$call$fixed)
            mr$call$random <- eval(mr$call$random)
            mr$call$rcov <- eval(mr$call$rcov)
            res$model[["hcs"]] <- mr
            if (!mr$converge) {
              mr$loglik <- -Inf
            }
          }
          nPar <- nlevels(TD$env) + 1
        } else if (choices[i] == "outside") {
          initVals <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                     random = as.formula("~ genotype:corh(env)"),
                                     start.values = TRUE, data = TD, ...)
          tmpTable <- initVals$gammas.table
          tmp <- tempfile()
          sink(file = tmp)
          tmpValues <- qvInitial(TD = TD, trait = trait, unitError = NA, vcmodel = "outside",
                                 fixed = as.formula(paste(trait, "~ env")), unitFactor = NA, ...)
          sink()
          unlink(tmp)
          tmpTable[, "Value"] <- c(c(tmpValues$vg, tmpValues$diag), 1)
          tmpTable[, "Constraint"] <- as.character(tmpTable[, "Constraint"] )
          tmpTable[which(tmpTable[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
          tmpTable[, "Constraint"] <- as.factor(tmpTable[, "Constraint"])
          mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                   random = as.formula("~ genotype:corh(env)"),
                                   G.param = tmpTable, R.param = tmpTable, data = TD, ...),
                    silent = TRUE)
          if (inherits(mr, "try-error")) {
            mr <- list(loglik = -Inf)
          } else {
            mr$call$fixed <- eval(mr$call$fixed)
            mr$call$random <- eval(mr$call$random)
            mr$call$rcov <- eval(mr$call$rcov)
            mr$call$G.param <- eval(mr$call$G.param)
            mr$call$R.param <- eval(mr$call$R.param)
            res$model[["outside"]] <- mr
            if (!mr$converge) {
              mr$loglik <- -Inf
            }
          }
          nPar <- nlevels(TD$env) + 1
        }
        if (nlevels(TD$env) > 4) {
          if (choices[i] == "fa") {
            initVals <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                       random = as.formula("~ genotype:fa(env, 1)"),
                                       start.values = TRUE, data = TD, ...)
            tmpTable <- initVals$gammas.table
            tmp <- tempfile()
            sink(file = tmp)
            tmpValues <- qvInitial(TD = TD, trait = trait, unitError = NA,
                                   vcmodel = "fa", fixed = as.formula(paste(trait, "~ env")),
                                   unitFactor = NA, ...)
            sink()
            unlink(tmp)
            if (is.null(tmpValues)) {
              mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                       random = as.formula("~ genotype:fa(env, 1)"),
                                       data = TD, ...), silent = TRUE)
            } else {
              tmpTable[, "Value"] <- c(tmpValues$psi, tmpValues$gamma, 1)
              tmpTable[, "Constraint"] <- as.character(tmpTable[, "Constraint"] )
              tmpTable[which(tmpTable[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
              tmpTable[, "Constraint"] <- as.factor(tmpTable[, "Constraint"])
              mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                       random = as.formula("~ genotype:fa(env, 1)"),
                                       R.param = tmpTable, G.param = tmpTable,
                                       data = TD, ...),
                        silent = TRUE)
            }
            if (inherits(mr, "try-error")) {
              mr <- list(loglik = -Inf)
            } else {
              mr$call$fixed <- eval(mr$call$fixed)
              mr$call$random <- eval(mr$call$random)
              mr$call$rcov <- eval(mr$call$rcov)
              mr$call$G.param <- eval(mr$call$G.param)
              mr$call$R.param <- eval(mr$call$R.param)
              res$model[["fa"]] <- mr
              if (!mr$converge) {
                mr$loglik <- -Inf
              }
            }
            nPar <- nlevels(TD$env) * 2
          }
          if (choices[i] == "fa2") {
            initVals <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                       random = as.formula("~ genotype:fa(env, 2)"),
                                       start.values = TRUE, data = TD, ...)
            tmpTable <- initVals$gammas.table
            tmp <- tempfile()
            sink(file = tmp)
            tmpValues <- qvInitial(TD = TD, trait = trait, unitError = NA,
                                   vcmodel = "fa2",
                                   fixed = as.formula(paste(trait, "~ env")),
                                   unitFactor = NA, ...)
            sink()
            unlink(tmp)
            if (is.null(tmpValues)) {
              mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                       random = as.formula("~ genotype:fa(env, 2)"),
                                       data = TD, ...), silent = TRUE)
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
              mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                       random = as.formula("~ genotype:fa(env, 2)"),
                                       R.param = tmpTable, G.param = tmpTable,
                                       data = TD, ...),
                        silent = TRUE)
            }
            if (inherits(mr, "try-error")) {
              mr <- list(loglik = -Inf)
            } else{
              mr$call$fixed <- eval(mr$call$fixed)
              mr$call$random <- eval(mr$call$random)
              mr$call$rcov <- eval(mr$call$rcov)
              mr$call$G.param <- eval(mr$call$G.param)
              mr$call$R.param <- eval(mr$call$R.param)
              res$model[["fa2"]] <- mr
              if (!mr$converge) {
                mr$loglik <- -Inf
              }
            }
            nPar <- nlevels(TD$env) * 3 - 1
          }
        }
        ## Check model.
        if (choices[i] == "unstructured") {
          initVals <- asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                     random = as.formula("~ genotype:us(env)"),
                                     start.values = TRUE, data = TD, ...)
          tmpTable <- initVals$gammas.table
          tmp <- tempfile()
          sink(file = tmp)
          tmpValues <- qvInitial(TD = TD, trait = trait, vcmodel = "unstructured",
                                 fixed = as.formula(paste(trait, "~ env")),
                                 unitFactor = NA, ...)
          sink()
          unlink(tmp)
          tmpValues <- tmpValues$evCov[upper.tri(tmpValues$evCov, diag = TRUE)]
          tmpTable[, "Value"] <- c(tmpValues, 1)
          tmpTable[, "Constraint"] <- as.character(tmpTable[, "Constraint"] )
          tmpTable[which(tmpTable[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
          tmpTable[, "Constraint"] <- as.factor(tmpTable[, "Constraint"])
          mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                   random = as.formula("~ genotype:us(env)"),
                                   G.param = tmpTable, R.param = tmpTable,
                                   data = TD, ...),
                    silent = TRUE)
          if (inherits(mr, "try-error")) {
            mr <- list(loglik = -Inf)
          } else {
            mr$call$fixed <- eval(mr$call$fixed)
            mr$call$random <- eval(mr$call$random)
            mr$call$rcov <- eval(mr$call$rcov)
            mr$call$G.param <- eval(mr$call$G.param)
            mr$call$R.param <- eval(mr$call$R.param)
            res$model[["unstructured"]] <- mr
            if (!mr$converge) {
              mr$loglik <- -Inf
            }
          }
          nPar <- nlevels(TD$env) * (nlevels(TD$env) - 1) / 2 +
            nlevels(TD$env)
        }
        if (!(nlevels(TD$env) <= 4 && choices[i] %in% c("fa", "fa2"))) {
          bestTab[choices[i], "AIC"] <- -2 * mr$loglik + 2 * nPar
          bestTab[choices[i], "BIC"] <- -2 * mr$loglik +
            log(length(mr$fitted.values)) * nPar
          bestTab[choices[i], "Deviance"] <- -2 * mr$loglik
          bestTab[choices[i], "NParameters"] <- nPar
          if (i == 1) {
            bestModel <- mr
            bestChoice <- choices[i]
            nBest <- nPar
          } else {
            if (criterion == "AIC") {
              criterionCur <- -2 * mr$loglik + 2 * nPar
              criterionPrev <- -2 * bestModel$loglik + 2 * nBest
            } else {
              criterionCur  <- -2 * mr$loglik + log(length(mr$fitted.values)) * nPar
              criterionPrev <- -2 * bestModel$loglik +
                log(length(bestModel$fitted.values)) * nBest
            }
            if (criterionCur < criterionPrev) {
              bestModel <- mr
              bestChoice <- choices[i]
              nBest <- nPar
            }
          }
        }
      }
      if (criterion == "AIC") {
        bestTab <- bestTab[order(bestTab[, "AIC"]), ]
      } else {
        bestTab <- bestTab[order(bestTab[, "BIC"]), ]
      }
      ## Outputs.
      res$choice <- bestChoice
      res$summaryTab <- bestTab
      tempFile <- tempfile()
      sink(tempFile)
      res$vcovBest <- predict(bestModel, classify = "env",
                              data = TD, vcov = TRUE)$predictions$vcov
      sink(NULL)
      unlink(tempFile)
      colnames(res$vcovBest) <- rownames(res$vcovBest) <- levels(TD$env)
      if (criterion == "AIC") {
        res$AICBest <- min(criterionCur, criterionPrev)
      } else {
        res$BICBest <- min(criterionCur, criterionPrev)
      }
    } else{
      stop("Failed to load 'asreml'.\n")
    }
  }
  return(res)
}
