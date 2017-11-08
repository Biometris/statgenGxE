#' Selects the best variance-covariance model for a set of environments
#'
#' This function selects the best covariance structure for genetic correlations between
#' environments. It fits a range of variance-covariance models to
#' compare (e.g., identity, compound symmetry, diagonal, heterogeneous compound symmetry,
#' first order factor analysis, second order factor analysis, unstructured),
#' and selects the best one using a goodness-of-fit criterion.
#'
#' @param Y A data frame object.
#' @param trait A character string specifying the name for a trait found in \code{Y}.
#' @param genotype A character string specifying the name for the genotypes found in \code{Y}.
#' @param env A character string specifying the name for enviroments found in \code{Y}.
#' @param engine A string specifying the name of a mixed modelling engine to be used.
#' @param criterion A string specifying a goodness-of-fit criterion, i.e., "AIC" or "BIC".
#' @param ... Further arguments to be passed to \code{asreml}.
#'
#' @note If \code{engine="lme4"}, only the compound symmetry model can be fitted.
#'
#' @return A list object consisting of the fitted model objects, a string specifying
#' the best model and its related goodness-of-fit criterion.
#'
#' @examples
#' mydat <- GE.read.csv(system.file("extdata", "F2maize_pheno.csv", package = "RAP"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' names(mydat)=c("env", "genotype","yld")
#' model1 <- GE.VarComp(mydat, trait="yld", genotype="genotype", env="env",
#'                      engine = "lme4", #engine = "asreml",
#'                      criterion = "BIC")
#' model1$BIC
#' model1$choice
#' summary(model1$model[[model1$choice]], nice=TRUE)
#'
#' @export
GE.VarComp <- function(Y,
                       trait,
                       genotype,
                       env,
                       engine = "asreml",
                       criterion = "BIC",
                       ...) {
  # keep NA values
  Y <- droplevels(Y[, c(env, genotype, trait)])
  naEntries <- which(is.na(tapply(X = rep(1, nrow(Y)),
                                  INDEX = Y[, c(genotype, env)],
                                  FUN = identity)),
                     arr.ind = TRUE)
  if (length(naEntries)) {
    tempNames <- names(Y)
    YExt <- data.frame(Y[naEntries[, env], env], Y[naEntries[, genotype], genotype],
                       rep(NA, nrow(naEntries)))
    names(YExt) <- tempNames
    Y <- rbind(Y, YExt)
  }
  #sort data in order of genotype, env
  Y <- Y[order(Y[[genotype]], Y[[env]]), ]
  qvInitial <- function(Y,
                        trait,
                        genotype,
                        env,
                        unitError = NA,
                        vcmodel = c("identity", "cs", "diagonal", "hcs",
                                    "outside", "fa", "fa2", "unstructured"),
                        fixed = NULL,
                        unitFactor = NA, ...) {
    # Replicates '_qvInitial' procedure (S. J. Welham 15/05/09) in GenStat
    # TODO: factanal() ? fa()
    # First, form estimate of unstructured matrix
    # exclude the rows including NA
    X <- na.omit(Y[, c(trait, genotype, env)])
    nEnv <- nlevels(X[, env])
    nGen <- nlevels(X[, genotype])
    # Get fixed df by stealth - in absence of other info, share among environments
    P <- 1
    if (!is.null(fixed)) {
      mr <- asreml::asreml(fixed = fixed, rcov = ~id(units), data = X, ...)
      P <- length(mr$fitted.values) - (1 + mr$nedf)
    }
    # Get number of effects contributing to each sum of squares
    tmpTable <- table(X[, c(env, genotype)])
    nobsEnv <- rowSums(tmpTable)
    Rnobs <- matrix(data = nobsEnv, nrow = nEnv, ncol = nEnv, byrow = TRUE)
    Cnobs <- t(Rnobs)
    Nobs <- Rnobs * (Rnobs - Cnobs < 0) + Cnobs * (Cnobs - Rnobs <= 0)
    if (!is.na(unitError)) {
      # This case is trickier, becuase of partitioning between two random terms,
      # but use of the diag structure is better than nothing!
      weights <- 1 / unitError
      X["weights"] <- weights
      if (!is.null(fixed)) {
        initValues <- asreml::asreml(fixed = fixed,
                                     random = as.formula(paste("~", genotype, ":idh(", env, ")")),
                                     weights = weights, start.values = TRUE, data = X, ...)
        tmp <- initValues$gammas.table
        tmp[, "Constraint"] <- as.character(tmp[,"Constraint"])
        tmp[which(tmp[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
        tmp[, "Constraint"] <- as.factor(tmp[, "Constraint"])
        mr <- asreml::asreml(fixed = fixed,
                             random = as.formula(paste("~", genotype, ":idh(", env, ")")),
                             weights = weights, R.param = tmp, data = X, ...)
        Ores <- matrix(data = mr$coeff$random, nrow = nGen, ncol = nEnv, byrow = TRUE)
      } else {
        initValues <- asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                                     random = as.formula(paste("~", genotype, ":idh(", env, ")")),
                                     weights = weights, start.values = TRUE, data = X, ...)
        tmp <- initValues$gammas.table
        tmp[ ,"Constraint"] <- as.character(tmp[, "Constraint"] )
        tmp[which(tmp[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
        tmp[, "Constraint"] <- as.factor(tmp[, "Constraint"])
        mr <- asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                             random = as.formula(paste("~", genotype, ":idh(", env, ")")),
                             weights = weights, R.param = tmp, data = X, ...)
        Ores <- matrix(data = mr$coeff$random, nrow = nGen, ncol = nEnv, byrow = TRUE)
      }
    } else {
      # This gives correct answers for complete balanced data
      if (!is.null(fixed)) {
        mr <- asreml::asreml(fixed = fixed, data = X, ...)
        residuals <- mr$residuals
      } else {
        mr <- asreml::asreml(fixed = as.formula(paste(trait, "~", env)), data = X, ...)
        residuals <- mr$residuals
      }
      Ores <- tapply(X = residuals, INDEX = list(X[, genotype], X[, env]), FUN = mean)
      if (sum(is.na(Ores)) != 0) {
        Ores[which(is.na(Ores))] <- 0
      }
    }
    Rmat <- Ores
    Evcov <- matrix(nrow = nEnv, ncol = nEnv)
    Evcov <- t(Rmat) %*% Rmat / (Nobs - (P / nEnv))
    # Get off-diagonal elements of Evcov only
    offdiag <- Evcov
    diag(offdiag) <- NA
    # Get off-diagonal elements of cor(Evcov) only
    offCorr <- cor(Evcov)
    diag(offCorr) <- NA
    #Get initial values for each model
    if (vcmodel == "identity") {
      vcInitial <- list(vge = mean(diag(Evcov)))
    } else if (vcmodel == "cs") {
      vg <- mean(offdiag, na.rm = TRUE)
      vge <- mean(diag(Evcov))
      vge <- (vge > vg) * (vge - vg) + (vge <= vg) * 0.1 * vge
      vcInitial <- list(vg = vg, vge = vge)
    } else if (vcmodel == "diagonal") {
      vcInitial <- list(diag = diag(Evcov))
    } else if (vcmodel == "hcs") {
      vg <- mean(offdiag, na.rm = TRUE) / 2
      diag <- diag(Evcov)
      diag <- (diag > vg) * (diag - vg) + (diag <= vg) * 0.1 * diag
      vcInitial <- list(vg = vg, diag = diag)
    } else if (vcmodel == "outside") {
      vg <- mean(offCorr, na.rm = TRUE)
      diag <- diag(Evcov)
      vcInitial <- list(vg = vg, diag = diag)
    } else if (vcmodel == "fa") {
      if (requireNamespace("psych", quietly = TRUE)) {
        factorAnalysis <- try(psych::fa(r = Evcov, nfactors = 1, fm = "mle"), silent = TRUE)
        if (inherits(factorAnalysis, "try-error")) {
          factorAnalysis <- psych::fa(r = Evcov, nfactors = 1, fm = "minres")
        }
        loading <- as.vector(factorAnalysis$loadings)
        comm <- factorAnalysis$comm
        var <- diag(Evcov)
        psi <- var * (1 - comm)
        Smat <- loading %*% t(loading)
        chol <- chol(Smat, pivot = TRUE)
        gamma <- chol[1, ]
        vcInitial <- list(gamma = gamma, psi = psi)
      } else {
        vcInitial <- NULL
        warning("psych package is required but failed to load.\n")
      }
    } else if (vcmodel == "fa2") {
      if (requireNamespace("psych", quietly = TRUE) &&
          requireNamespace("GPArotation", quietly = TRUE)) {
        factorAnalysis <- try(psych::fa(r = Evcov, nfactors = 2, fm = "mle"), silent = TRUE)
        if (inherits(factorAnalysis, "try-error")) {
          factorAnalysis <- psych::fa(r = Evcov, nfactors = 2, fm = "minres")
        }
        loading <- as.vector(factorAnalysis$loadings[, 1])
        comm <- factorAnalysis$comm
        var <- diag(Evcov)
        psi <- var * (1 - comm)
        Smat <- loading %*% t(loading)
        chol <- chol(Smat, pivot = TRUE)
        gamma <- chol[1:2, ]
        vcInitial <- list(gamma = gamma, psi = psi)
      }else{
        vcInitial <- NULL
        warning("psych and GPArotation packages are required but failed to load!\n")
      }
    } else if (vcmodel == "unstructured") {
      vcInitial <- list(Evcov = Evcov)
    }
    return(vcInitial)
  }
  # Main procedure to fit mixed models
  if (engine == "lme4") {
    # Compound symmetry ("cs") only
    if (requireNamespace("lme4", quietly = TRUE)) {
      mr = lme4::lmer(as.formula(paste(trait, "~ ", env, " +(1 |", genotype, ")")),
                      data = Y, ...)
      nPar <- 2
      # Outputs
      res <- new.env()
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
      res <- new.env()
      for (i in 1:length(choices)) {
        if (choices[i] == "identity") {
          mr <- asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                               rcov = as.formula(paste("~", genotype, ":", env)), data = Y, ...)
          mr$call$fixed <- eval(mr$call$fixed)
          mr$call$rcov <- eval(mr$call$rcov)
          res$model[["identity"]] <- mr
          if (!mr$converge) {
            mr$loglik <- -Inf
          }
          nPar <- 1
        } else if (choices[i] == "cs") {
          mr <- asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                               random = as.formula(paste("~", genotype)),
                               rcov = as.formula(paste("~", genotype, ":", env)), data = Y, ...)
          mr$call$fixed <- eval(mr$call$fixed)
          mr$call$random <- eval(mr$call$random)
          mr$call$rcov <- eval(mr$call$rcov)
          res$model[["cs"]] <- mr
          if (!mr$converge) {
            mr$loglik <- -Inf
          }
          nPar <- 2
        } else if (choices[i] == "diagonal") {
          mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                                   rcov = as.formula(paste("~", genotype, ":diag(", env, ")")),
                                   data = Y, ...), silent = TRUE)
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
          nPar <- length(levels(Y[, env]))
        } else if (choices[i] == "hcs") {
          mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                                   random = as.formula(paste("~", genotype)),
                                   rcov = as.formula(paste("~", genotype, ":diag(", env, ")")),
                                   data = Y, ...), silent = TRUE)
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
          nPar <- length(levels(Y[,env])) + 1
        } else if (choices[i] == "outside") {
          initVals <- asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                                     random = as.formula(paste("~", genotype, ":corh(", env, ")")),
                                     start.values = TRUE, data = Y, ...)
          tmpTable <- initVals$gammas.table
          tmp <- tempfile()
          sink(file = tmp)
          tmpValues <- qvInitial(Y = Y, trait = trait, genotype = genotype, env = env,
                                 unitError = NA, vcmodel = "outside",
                                 fixed = as.formula(paste(trait, "~", env)), unitFactor = NA, ...)
          sink()
          unlink(tmp)
          #print(tmpTable)
          tmpTable[, "Value"] <- c(c(tmpValues$vg, tmpValues$diag), 1)
          tmpTable[, "Constraint"] <- as.character(tmpTable[, "Constraint"] )
          tmpTable[which(tmpTable[, "Gamma"]== "R!variance"), "Constraint"] <- "F"
          tmpTable[, "Constraint"] <- as.factor(tmpTable[, "Constraint"])
          #print(tmpTable)
          mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                                   random = as.formula(paste("~", genotype, ":corh(", env,")")),
                                   G.param = tmpTable, R.param = tmpTable, data = Y, ...),
                    silent=TRUE)
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
          nPar <- length(levels(Y[, env])) + 1
        }
        if (nlevels(Y[[env]]) > 4) {
          if (choices[i] == "fa") {
            initVals <- asreml::asreml(fixed = as.formula(paste(trait, "~" ,env)),
                                       random = as.formula(paste("~", genotype, ":fa(", env, ",1)")),
                                       start.values = TRUE, data = Y, ...)
            tmpTable <- initVals$gammas.table
            tmp <- tempfile()
            sink(file = tmp)
            tmpValues <- qvInitial(Y = Y, trait = trait, genotype = genotype,
                                   env = env, unitError = NA,
                                   vcmodel = "fa", fixed = as.formula(paste(trait, "~", env)),
                                   unitFactor = NA, ...)
            sink()
            unlink(tmp)
            if (is.null(tmpValues)) {
              mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                                       random = as.formula(paste("~", genotype, ":fa(", env, ",1)")),
                                       data = Y, ...), silent = TRUE)
            } else {
              tmpTable[, "Value"] <- c(tmpValues$psi, tmpValues$gamma, 1)
              tmpTable[, "Constraint"] <- as.character(tmpTable[, "Constraint"] )
              tmpTable[which(tmpTable[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
              tmpTable[, "Constraint"] <- as.factor(tmpTable[, "Constraint"])
              mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                                       random = as.formula(paste("~", genotype, ":fa(", env, ",1)")),
                                       R.param = tmpTable, G.param = tmpTable, data = Y, ...),
                        silent=TRUE)
            }
            if (inherits(mr, "try-error")){
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
            nPar <- length(levels(Y[, env])) * 2
          }
          if (choices[i] == "fa2"){
            initVals <- asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                                       random = as.formula(paste("~", genotype, ":fa(", env, ",2)")),
                                       start.values = T, data = Y, ...)
            tmpTable <- initVals$gammas.table
            tmp <- tempfile()
            sink(file = tmp)
            tmpValues <- qvInitial(Y = Y, trait = trait, genotype = genotype, env = env,
                                   unitError = NA, vcmodel = "fa2",
                                   fixed = as.formula(paste(trait, "~", env)), unitFactor = NA, ...)
            sink()
            unlink(tmp)
            if (is.null(tmpValues)) {
              mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                                       random = as.formula(paste("~", genotype, ":fa(", env, ",2)")),
                                       data = Y, ...), silent = TRUE)
            } else {
              # Keep loadings of factor 2 away from 0
              tmpValues$gamma[2, which(tmpValues$gamma[2,] < 1e-3)] <- 1e-3
              # Make sure that first entry is 0
              tmpValues$gamma[2, 1] <- 0
              tmpTable[, "Value"] <- c(tmpValues$psi, tmpValues$gamma[1, ], tmpValues$gamma[2, ], 1)
              tmpTable[, "Constraint"] <- as.character(tmpTable[, "Constraint"])
              tmpTable[which(tmpTable[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
              tmpTable[, "Constraint"] <- as.factor(tmpTable[, "Constraint"])
              mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                                       random = as.formula(paste("~", genotype, ":fa(", env, ",2)")),
                                       R.param = tmpTable, G.param = tmpTable, data = Y, ...),
                        silent=TRUE)
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
            nPar <- length(levels(Y[, env])) * 3 - 1
          }
        }
        # Check model
        if (choices[i] == "unstructured"){
          initVals <- asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                                     random = as.formula(paste("~", genotype, ":us(", env, ")")),
                                     start.values = TRUE, data = Y, ...)
          tmpTable <- initVals$gammas.table
          tmp <- tempfile()
          sink(file = tmp)
          tmpValues <- qvInitial(Y = Y, trait = trait, genotype = genotype, env = env,
                                 unitError = NA, vcmodel = "unstructured",
                                 fixed = as.formula(paste(trait, "~", env)), unitFactor = NA, ...)
          sink()
          unlink(tmp)
          tmpValues <- tmpValues$Evcov[upper.tri(tmpValues$Evcov, diag = TRUE)]
          tmpTable[, "Value"] <- c(tmpValues, 1)
          tmpTable[, "Constraint"] <- as.character(tmpTable[, "Constraint"] )
          tmpTable[which(tmpTable[, "Gamma"] == "R!variance"), "Constraint"] <- "F"
          tmpTable[, "Constraint"] <- as.factor(tmpTable[, "Constraint"])
          mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                                   random = as.formula(paste("~", genotype, ":us(", env, ")")),
                                   G.param = tmpTable, R.param = tmpTable, data = Y, ...),
                    silent=TRUE)
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
          nPar <- length(levels(Y[, env])) * (length(levels(Y[, env])) - 1) / 2 +
            length(levels(Y[, env]))
        }
        if (!(nlevels(Y[[env]]) <= 4 && choices[i] %in% c("fa", "fa2"))) {
          bestTab[choices[i], "AIC"] <- -2 * mr$loglik + 2 * nPar
          bestTab[choices[i], "BIC"] <- -2 * mr$loglik + log(length(mr$fitted.values)) * nPar
          bestTab[choices[i], "Deviance"] <- -2 * mr$loglik
          bestTab[choices[i], "NParameters"] <- nPar
          if (i == 1){
            bestModel <- mr
            bestChoice <- choices[i]
            nBest <- nPar
          } else {
            if (criterion == "AIC") {
              criterionCur <- -2 * mr$loglik + 2 * nPar
              criterionPrev <- -2 * bestModel$loglik + 2 * nBest
            } else {
              criterionCur  <- -2 * mr$loglik + log(length(mr$fitted.values)) * nPar
              criterionPrev <- -2 * bestModel$loglik + log(length(bestModel$fitted.values)) * nBest
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
      # Outputs
      res$choice <- bestChoice
      res$summaryTab <- bestTab
      tempFile <- tempfile()
      sink(tempFile)
      res$vcovBest <- predict(bestModel, classify = env, data = Y, vcov = TRUE)$predictions$vcov
      sink(NULL)
      unlink(tempFile)
      colnames(res$vcovBest) <- rownames(res$vcovBest) <- levels(Y[, env])
      if (criterion == "AIC"){
        res$AICBest <- min(criterionCur, criterionPrev)
      }else{
        res$BICBest <- min(criterionCur, criterionPrev)
      }
    } else{
      stop("Failed to load 'asreml'.\n")
    }
  }
  if (engine %in% c("asreml", "lme4")) {
    return(as.list(res))
  } else {
    stop("Please use either 'asreml' or 'lme4' as the mixing modelling engine.\n")
  }
}
