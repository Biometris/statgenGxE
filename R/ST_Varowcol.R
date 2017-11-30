#' Find the best random model for a row-column design (asreml only)
#'
#' This function fits a variety of random and spatial covariance models
#' and selects the best one using a goodness-of-fit criterion.
#'
#' @inheritParams ST.run.model
#' @param tryRep A logical value indicating if 'replicates' are included in the model.
#' Default, \code{TRUE}.
#' @param criterion A string specifies a goodness of fit criterion, i.e., "AIC" or "BIC".
#' @param ... Further arguments to be passed to \code{asreml}.
#'
#' @return an object of class \code{\link{SSA}}.
#'
#' @note This function can only be used if asreml is installed. If \code{trySpatial} is set
#' to "always" or "ifregular", the names for \code{rowCoordinates} and \code{colCoordinates}
#' must be supplied; otherwise, no spatial model will be fitted.
#'
#' @export

ST.Varowcol <- function(TD,
                        trait,
                        covariates = NULL,
                        repId = NULL,
                        rowId = NULL,
                        colId = NULL,
                        tryRep = TRUE,
                        checkId = NULL,
                        rowCoordinates = NULL,
                        colCoordinates = NULL,
                        trySpatial = NULL,
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
  if (!is.null(covariates) && (!is.character(covariates) ||
                              !(all(covariates %in% colnames(TD))))) {
    stop("covariates have to be a columns in TD.\n")
  }
  for (param in c(repId, rowId, colId, rowCoordinates, colCoordinates, checkId)) {
    if (!is.null(param) && (!is.character(param) || length(param) > 1 ||
                            !param %in% colnames(TD))) {
      stop(paste(deparse(param), "has to be NULL or a column in data.\n"))
    }
  }
  if (!is.null(trySpatial) && (!is.character(trySpatial) || length(trySpatial) > 1 ||
                               !trySpatial %in% c("always", "ifregular"))) {
    stop("trySpatial should be NULL, always or ifregular.\n")
  }
  # TODO: Starting values for more complex models? (See some warnings captured from asreml)
  # Run mixed and fixed models using asreml
  # check validity of variable name, trait
  ok <- isValidVariableName(trait)
  trait0 <- trait
  if (!all(ok)) {
    trait0[!ok] <- sapply(X = trait0[!ok], FUN = function(x) {
      paste0("`", x, "`")
    })
  }
  covT <- FALSE
  if (!is.null(covariates)) {
    if (is.character(covariates)) {
      covT <- TRUE
    }
  }
  flag <- 1
  # See if the design is regular
  if (!missing(repId)) {
    reptab <- table(TD[[repId]], TD[[rowId]], TD[[colId]])
  } else {
    reptab <- table(TD[[rowId]], TD[[colId]])
  }
  if (min(reptab) > 1) {
    warning("There must be only one plot at each REPLICATES x ROWS x COLUMNS location.\n
    spatial models will not be tried.\n")
    flag <- 0
  }
  if (min(reptab) == 1 && max(reptab) == 1) {
    regular <- TRUE
  } else {
    regular <- FALSE
  }
  if (is.null(rowCoordinates)) {
    flag <- 0
  }
  if (is.null(colCoordinates)) {
    flag <- 0
  }
  # does not use spatial models
  if (flag == 0) {
    trySpatial <- NULL
  }
  tmp <- tempfile()
  sink(file = tmp)
  # default no spatial models
  if (is.null(trySpatial)) {
    if (tryRep) {
      fixedFormR <- as.formula(paste(trait0, "~", repId,
                                     if (!is.null(checkId)) paste("+", checkId),
                                     if (covT) paste(c("", covariates),
                                                     collapse = "+")))
      mr <- asreml::asreml(fixed = fixedFormR,
                           random = as.formula(paste("~ genotype +", repId, ":",
                                                     rowId, "+", repId, ":", colId)),
                           rcov = ~ units, aom = TRUE, data = TD, ...)
      # constrain variance of the variance components to be fixed as the values in mr
      GParamTmp <- mr$G.param
      tmpPos <- which(names(GParamTmp) == paste(repId, rowId, sep = ":"))
      if (length(tmpPos) > 0) {
        if (!is.null(GParamTmp[[tmpPos]][[repId]]$con)) {
          GParamTmp[[tmpPos]][[repId]]$con <- "F"
        }
      }
      tmpPos <- which(names(GParamTmp) == paste(repId, colId, sep = ":"))
      if (length(tmpPos) > 0) {
        if (!is.null(GParamTmp[[tmpPos]][[repId]]$con)) {
          GParamTmp[[tmpPos]][[repId]]$con <- "F"
        }
      }
      fixedFormF <- as.formula(paste(deparse(fixedFormR), "+ genotype"))
      mf <- asreml::asreml(fixed = fixedFormF,
                           random = as.formula(paste("~", repId, ":", rowId, "+",
                                                     repId, ":", colId)),
                           rcov = ~ units, G.param = GParamTmp, aom = TRUE,
                           data = TD, ...)

    } else {
      fixedFormR <- as.formula(paste(trait0, "~",
                                     if (!is.null(checkId)) paste("+", checkId),
                                     if (covT) paste(c("", covariates),
                                                     collapse = "+")))
      mr <- asreml::asreml(fixed = fixedFormR,
                           random = as.formula(paste("~ genotype +", rowId, "+", colId)),
                           rcov = ~ units, aom = TRUE, data = TD, ...)
      # constrain variance of the variance components to be fixed as the values in mr
      GParamTmp <- mr$G.param
      if (!is.null(GParamTmp[[rowId]][[rowId]]$con)) {
        GParamTmp[[rowId]][[rowId]]$con <- "F"
      }
      if (!is.null(GParamTmp[[colId]][[colId]]$con)) {
        GParamTmp[[colId]][[colId]]$con <- "F"
      }
      fixedFormF <- as.formula(paste(deparse(fixedFormR), "+ genotype"))
      mf <- asreml::asreml(fixed = fixedFormF,
                           random = as.formula(paste("~", rowId, "+", colId)),
                           rcov = ~ units, G.param = GParamTmp, aom = TRUE,
                           data = TD, ...)
    }
    bestModel <- mr
  } else {
    if (regular) {
      if (tryRep) {
        randomChoice <- c(rep(x = c("Identity", "Measurement_error"), each = 3),
                          paste(repId, rowId, sep = ":"), paste(repId, colId, sep = ":"),
                          paste(paste(repId, rowId, sep = ":"), paste(repId, colId, sep = ":"),
                                sep = "+"),
                          paste(paste(repId, rowId, sep = ":"), "Measurement_error", sep = "+"),
                          paste(paste(repId, colId, sep = ":"), "Measurement_error", sep = "+"),
                          paste(paste(repId, rowId, sep = ":"), paste(repId, colId, sep = ":"),
                                "Measurement_error", sep = "+"))
        randomTerm <- c(rep(x = c("NULL", "units"), each = 3),
                        paste(repId, rowId, sep = ":"), paste(repId, colId, sep = ":"),
                        paste(paste(repId, rowId, sep = ":"), paste(repId, colId, sep = ":"),
                              sep = "+"),
                        paste(paste(repId, rowId, sep = ":"), "units", sep = "+"),
                        paste(paste(repId, colId, sep = ":"), "units", sep = "+"),
                        paste(paste(repId, rowId, sep = ":"), paste(repId, colId, sep = ":"),
                              "units", sep = "+"))
      } else {
        randomChoice <- c(rep(x = c("Identity", "Measurement_error"), each = 3), rowId, colId,
                          paste(rowId, colId, sep = "+"), paste(rowId, "Measurement_error",
                                                            sep = "+"),
                          paste(colId, "Measurement_error", sep = "+"),
                          paste(rowId, colId, "Measurement_error", sep = "+"))
        randomTerm <- c(rep(x = c("NULL", "units"), each = 3), rowId, colId,
                        paste(rowId, colId, sep = "+"), paste(rowId, "units", sep = "+"),
                        paste(colId, "units", sep = "+"), paste(rowId, colId, "units", sep = "+"))
      }
      spatialChoice <- rep(x = c("AR1(x)Identity", "Identity(x)AR1", "AR1(x)AR1"), times = 4)
      spatialTerm <- rep(x = c(paste("ar1(", rowCoordinates, "):", colCoordinates),
                               paste(rowCoordinates, ":ar1(", colCoordinates,")"),
                               paste("ar1(", rowCoordinates, "):ar1(", colCoordinates, ")")),
                         times = 4)
      modelChoice <- paste("Random:", randomChoice, "&   Spatial:", spatialChoice)
      for (ii in 1:length(randomTerm)) {
        if (tryRep) {
          fixedFormR <- as.formula(paste(trait0, "~", repId,
                                         if (!is.null(checkId)) paste("+", checkId),
                                         if (covT) paste(c("", covariates),
                                                         collapse = "+")))
          mr <- asreml::asreml(fixed = fixedFormR,
                               random = as.formula(paste("~ genotype +",
                                                         randomTerm[ii])),
                               rcov = as.formula(paste("~", spatialTerm[ii])),
                               aom = TRUE, data = TD, ...)
        } else {
          fixedFormR <- as.formula(paste(trait0, "~",
                                         if (!is.null(checkId)) checkId else "1",
                                         if (covT) paste(c("", covariates),
                                                         collapse = "+")))
          mr <- asreml::asreml(fixed = fixedFormR,
                               random = as.formula(paste("~ genotype +",
                                                         randomTerm[ii])),
                               rcov = as.formula(paste("~", spatialTerm[ii])),
                               aom = TRUE, data = TD, ...)
        }
        if (ii == 1) {
          bestModel <- mr
          bestChoice <- modelChoice[ii]
          bestLoc <- 1
        } else {
          if (criterion == "AIC") {
            criterionCur  <- -2 * mr$loglik + 2 * length(mr$gammas)
            criterionPrev <- -2 * bestModel$loglik + 2 * length(bestModel$gammas)
          } else {
            criterionCur  <- -2 * mr$loglik + log(length(mr$fitted.values)) * length(mr$gammas)
            criterionPrev <- -2 * bestModel$loglik +
              log(length(bestModel$fitted.values)) * length(bestModel$gammas)
          }
          if (criterionCur < criterionPrev) {
            bestModel <- mr
            bestChoice <- modelChoice[ii]
            bestLoc <- ii
          }
        }
      }
    } else {
      if (!is.null(trySpatial) && trySpatial == "always") {
        if (tryRep) {
          randomChoice <- c(rep(x = c("Identity", "Measurement_error"), each = 3),
                            paste(repId, rowId, sep = ":"), paste(repId, colId, sep = ":"),
                            paste(paste(repId, rowId, sep = ":"), paste(repId, colId, sep = ":"),
                                  sep = "+"),
                            paste(paste(repId, rowId, sep = ":"), "Measurement_error", sep = "+"),
                            paste(paste(repId, colId, sep = ":"), "Measurement_error", sep = "+"),
                            paste(paste(repId, rowId, sep = ":"), paste(repId, colId, sep = ":"),
                                  "Measurement_error", sep = "+"))
          randomTerm <- c(rep(x = c("NULL", "units"), each = 3),
                          paste(repId, rowId, sep = ":"), paste(repId, colId, sep = ":"),
                          paste(paste(repId, rowId, sep = ":"), paste(repId, colId, sep = ":"),
                                sep = "+"),
                          paste(paste(repId, rowId, sep = ":"), "units", sep = "+"),
                          paste(paste(repId, colId, sep = ":"), "units", sep = "+"),
                          paste(paste(repId, rowId, sep = ":"), paste(repId, colId, sep = ":"),
                                "units", sep = "+"))
        } else {
          randomChoice <- c(rep(x = c("Identity", "Measurement_error"), each = 3), rowId, colId,
                            paste(rowId, colId, sep = "+"), paste(rowId, "Measurement_error", sep = "+"),
                            paste(colId, "Measurement_error", sep = "+"),
                            paste(rowId, colId, "Measurement_error", sep = "+"))
          randomTerm <- c(rep(x = c("NULL", "units"), each = 3),
                          rowId, colId, paste(rowId, colId, sep = "+"),
                          paste(rowId, "units", sep = "+"), paste(colId, "units", sep = "+"),
                          paste(rowId, colId, "units", sep = "+"))
        }
        spatialChoice <- rep(c("Exponential(x)Identity", "Identity(x)Exponential",
                               "Isotropic exponential"), times = 4)
        spatialTerm <- rep(c(paste0("exp(", rowCoordinates, "):", colCoordinates),
                             paste0(rowCoordinates, ":exp(", colCoordinates, ")"),
                             paste0("iexp(", rowCoordinates, ",", colCoordinates, ")")),
                           times = 4)
        modelChoice <- paste("Random:", randomChoice, "&   Spatial:", spatialChoice)
        for (ii in 1:length(randomTerm)) {
          if (tryRep) {
            fixedFormR <- as.formula(paste(trait0, "~", repId,
                                           if (!is.null(checkId)) paste("+", checkId),
                                           if (covT) paste(c("", covariates),
                                                           collapse = "+")))
            mr <- asreml::asreml(fixed = fixedFormR,
                                 random = as.formula(paste("~ genotype +",
                                                           randomTerm[ii])),
                                 rcov = as.formula(paste("~", spatialTerm[ii])),
                                 aom = TRUE, data = TD, ...)
          } else {
            fixedFormR <- as.formula(paste(trait0, "~",
                                           if (!is.null(checkId)) checkId else "1",
                                           if (covT) paste(c("", covariates),
                                                           collapse = "+")))
            mr <- asreml::asreml(fixed = fixedFormR,
                                 random = as.formula(paste("~ genotype +",
                                                           randomTerm[ii])),
                                 rcov = as.formula(paste("~", spatialTerm[ii])),
                                 aom = TRUE, data = TD, ...)
          }
          if (ii == 1) {
            bestModel <- mr
            bestChoice <- modelChoice[ii]
            bestLoc <- 1
          } else {
            if (criterion == "AIC") {
              criterionCur  <- -2 * mr$loglik + 2 * length(mr$gammas)
              criterionPrev <- -2 * bestModel$loglik + 2 * length(bestModel$gammas)
            } else {
              criterionCur  <- -2 * mr$loglik + log(length(mr$fitted.values)) * length(mr$gammas)
              criterionPrev <- -2 * bestModel$loglik +
                log(length(bestModel$fitted.values)) * length(bestModel$gammas)
            }
            if (criterionCur < criterionPrev) {
              bestModel <- mr
              bestChoice <- modelChoice[ii]
              bestLoc <- ii
            }
          }
        }
      }
    }
    fixedFormF <- as.formula(paste(deparse(fixedFormR), "+ genotype"))
    # constrain variance of the variance components to be fixed as the values in mr
    if (tryRep) {
      GParamTmp = bestModel$G.param
      tmpPos <- which(names(GParamTmp) == paste(repId, rowId, sep = ":"))
      if (length(tmpPos) > 0) {
        if (!is.null(GParamTmp[[tmpPos]][[repId]]$con)) {
          GParamTmp[[tmpPos]][[repId]]$con <- "F"
        }
      }
      tmpPos <- which(names(GParamTmp) == paste(repId, colId, sep = ":"))
      if (length(tmpPos) > 0) {
        if (!is.null(GParamTmp[[tmpPos]][[repId]]$con)) {
          GParamTmp[[tmpPos]][[repId]]$con <- "F"
        }
      }
      mf <- asreml::asreml(fixed = fixedFormF,
                           random = as.formula(paste("~", randomTerm[bestLoc])),
                           rcov = as.formula(paste("~", spatialTerm[ii])),
                           G.param = GParamTmp, aom = TRUE, data = TD, ...)
    } else {
      GParamTmp <- bestModel$G.param
      if (!is.null(GParamTmp[[rowId]][[rowId]]$con)) {
        GParamTmp[[rowId]][[rowId]]$con <- "F"
      }
      if (!is.null(GParamTmp[[colId]][[colId]]$con)) {
        GParamTmp[[colId]][[colId]]$con <- "F"
      }
      mf <- asreml::asreml(fixed = fixedFormF,
                           random = as.formula(paste("~", randomTerm[bestLoc])),
                           rcov = as.formula(paste("~", spatialTerm[ii])),
                           G.param = GParamTmp, aom = TRUE, data = TD, ...)
    }
  }
  # run predict
  if (!is.null(trySpatial)) {
    ii <- bestLoc
  }
  bestModel$call$fixed <- eval(bestModel$call$fixed)
  bestModel$call$random <- eval(bestModel$call$random)
  bestModel$call$rcov <- eval(bestModel$call$rcov)
  mf$call$fixed <- eval(mf$call$fixed)
  mf$call$random <- eval(mf$call$random)
  mf$call$rcov <- eval(mf$call$rcov)
  bestModel = predict(bestModel, classify = "genotype", data = TD)
  if (!is.null(checkId)) {
    mf <- predict(mf, classify = "genotype", vcov = TRUE, data = TD,
                  associate = as.formula(paste0("~", checkId,":genotype")))
  } else {
    mf <- predict(mf, classify = "genotype", vcov = TRUE, data = TD)
  }
  sink()
  unlink(tmp)
  res = createSSA(mMix = bestModel, mFix = mf, data = TD, trait = trait,
                  engine = "asreml")
  return(res)
}
