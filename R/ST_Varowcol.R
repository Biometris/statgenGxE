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
                        covariate = NULL,
                        rep,
                        row,
                        col,
                        tryRep = TRUE,
                        checkId = NA,
                        rowCoordinates = NA,
                        colCoordinates = NA,
                        trySpatial = NA,
                        criterion = "BIC",
                        ...) {
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
  if (!is.null(covariate)) {
    if (is.character(covariate)) {
      covT <- TRUE
    }
  }
  flag <- 1
  # See if the design is regular
  if (!missing(rep)) {
    reptab <- table(TD[[rep]], TD[[row]], TD[[col]])
  } else {
    reptab <- table(TD[[row]], TD[[col]])
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
  if (is.na(rowCoordinates)) {
    flag <- 0
  }
  if (is.na(colCoordinates)) {
    flag <- 0
  }
  # does not use spatial models
  if (flag == 0) {
    trySpatial <- NA
  }
  tmp <- tempfile()
  sink(file = tmp)
  # default no spatial models
  if (is.na(trySpatial)) {
    if (tryRep) {
      if (!is.na(checkId)) {
        mr <- asreml::asreml(fixed = as.formula(paste(trait0, "~", rep, "+", checkId,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"))),
                             random = as.formula(paste("~ genotype +", rep, ":",
                                                       row, "+", rep, ":", col)),
                             rcov = ~units, aom = TRUE, data = TD, ...)
      } else {
        mr <- asreml::asreml(fixed = as.formula(paste(trait0, "~", rep,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"))),
                             random = as.formula(paste("~ genotype +", rep, ":", row,
                                                       "+", rep, ":", col)),
                             rcov = ~units, aom = TRUE, data = TD, ...)
      }
      # constrain variance of the variance components to be fixed as the values in mr
      GParamTmp <- mr$G.param
      tmpPos <- which(names(GParamTmp) == paste(rep, row, sep = ":"))
      if (length(tmpPos) > 0) {
        if (!is.null(GParamTmp[[tmpPos]][[rep]]$con)) {
          GParamTmp[[tmpPos]][[rep]]$con <- "F"
        }
      }
      tmpPos <- which(names(GParamTmp) == paste(rep, col, sep = ":"))
      if (length(tmpPos) > 0) {
        if (!is.null(GParamTmp[[tmpPos]][[rep]]$con)) {
          GParamTmp[[tmpPos]][[rep]]$con <- "F"
        }
      }
      if (!is.na(checkId)) {
        mf <- asreml::asreml(fixed = as.formula(paste(trait0, "~", rep, "+", checkId,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"),
                                                      "+ genotype")),
                             random = as.formula(paste("~", rep, ":", row, "+",
                                                       rep, ":", col)),
                             rcov = ~units, G.param = GParamTmp, aom = TRUE,
                             data = TD, ...)
      } else {
        mf <- asreml::asreml(fixed = as.formula(paste(trait0, "~", rep,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"),
                                                      "+ genotype")),
                             random = as.formula(paste("~", rep, ":", row, "+",
                                                       rep, ":", col)),
                             rcov = ~units, G.param = GParamTmp, aom = TRUE,
                             data = TD, ...)
      }
    } else {
      if (!is.na(checkId)) {
        mr <- asreml::asreml(fixed = as.formula(paste(trait0, "~", checkId,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"))),
                             random = as.formula(paste("~ genotype +", row, "+", col)),
                             rcov = ~units, aom = TRUE, data = TD, ...)
      } else {
        mr <- asreml::asreml(fixed = as.formula(paste(trait0, "~", 1,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"))),
                             random = as.formula(paste("~ genotype +", row, "+", col)),
                             rcov = ~units, aom = TRUE, data = TD, ...)
      }
      # constrain variance of the variance components to be fixed as the values in mr
      GParamTmp <- mr$G.param
      if (!is.null(GParamTmp[[row]][[row]]$con)) {
        GParamTmp[[row]][[row]]$con <- "F"
      }
      if (!is.null(GParamTmp[[col]][[col]]$con)) {
        GParamTmp[[col]][[col]]$con <- "F"
      }
      if (!is.na(checkId)) {
        mf <- asreml::asreml(fixed = as.formula(paste(trait0, "~", checkId,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"),
                                                      "+ genotype")),
                             random = as.formula(paste("~", row, "+", col)),
                             rcov = ~units, G.param = GParamTmp, aom = TRUE,
                             data = TD, ...)
      } else {
        mf <- asreml::asreml(fixed = as.formula(paste(trait0, "~1",
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"),
                                                      "+ genotype")),
                             random = as.formula(paste("~", row, "+", col)),
                             rcov = ~units, G.param = GParamTmp, aom = TRUE,
                             data = TD, ...)
      }
    }
    bestModel <- mr
  } else {
    if (regular) {
      if (tryRep) {
        randomChoice <- c(rep(x = c("Identity", "Measurement_error"), each = 3),
                         paste(rep, row, sep = ":"), paste(rep, col, sep = ":"),
                         paste(paste(rep, row, sep = ":"), paste(rep, col, sep = ":"),
                               sep = "+"),
                         paste(paste(rep, row, sep = ":"), "Measurement_error", sep = "+"),
                         paste(paste(rep, col, sep = ":"), "Measurement_error", sep = "+"),
                         paste(paste(rep, row, sep = ":"), paste(rep, col, sep = ":"),
                               "Measurement_error", sep = "+"))
        randomTerm <- c(rep(x = c("NULL", "units"), each = 3),
                        paste(rep, row, sep = ":"), paste(rep, col, sep = ":"),
                        paste(paste(rep, row, sep = ":"), paste(rep, col, sep = ":"),
                              sep = "+"),
                        paste(paste(rep, row, sep = ":"), "units", sep = "+"),
                        paste(paste(rep, col, sep = ":"), "units", sep = "+"),
                        paste(paste(rep, row, sep = ":"), paste(rep, col, sep = ":"),
                              "units", sep = "+"))
      } else {
        randomChoice <- c(rep(x = c("Identity", "Measurement_error"), each = 3), row, col,
                         paste(row, col, sep = "+"), paste(row, "Measurement_error",
                                                           sep = "+"),
                         paste(col, "Measurement_error", sep = "+"),
                         paste(row, col, "Measurement_error", sep = "+"))
        randomTerm <- c(rep(x = c("NULL", "units"), each = 3), row, col,
                        paste(row, col, sep = "+"), paste(row, "units", sep = "+"),
                        paste(col, "units", sep = "+"), paste(row, col, "units", sep = "+"))
      }
      spatialChoice <- rep(x = c("AR1(x)Identity", "Identity(x)AR1", "AR1(x)AR1"), times = 4)
      spatialTerm <- rep(x = c(paste("ar1(", rowCoordinates, "):", colCoordinates),
                               paste(rowCoordinates, ":ar1(", colCoordinates,")"),
                               paste("ar1(", rowCoordinates, "):ar1(", colCoordinates, ")")),
                         times = 4)
      modelChoice <- paste("Random:", randomChoice, "&   Spatial:", spatialChoice)
      for (ii in 1:length(randomTerm)) {
        if (tryRep) {
          if (!is.na(checkId)) {
            mr <- asreml::asreml(fixed = as.formula(paste(trait0, "~", rep, "+", checkId,
                                                          if (covT) paste(c("", covariate),
                                                                          collapse = "+"))),
                                 random = as.formula(paste("~ genotype +", randomTerm[ii])),
                                 rcov = as.formula(paste("~", spatialTerm[ii])),
                                 aom = TRUE, data = TD, ...)
          } else {
            mr <- asreml::asreml(fixed = as.formula(paste(trait0, "~", rep,
                                                          if (covT) paste(c("", covariate),
                                                                          collapse = "+"))),
                                 random = as.formula(paste("~ genotype +", randomTerm[ii])),
                                 rcov = as.formula(paste("~", spatialTerm[ii])),
                                 aom = TRUE, data = TD, ...)
          }
        } else {
          if (!is.na(checkId)) {
            mr <- asreml::asreml(fixed = as.formula(paste(trait0, "~", checkId,
                                                          if (covT) paste(c("", covariate),
                                                                          collapse = "+"))),
                                 random = as.formula(paste("~ genotype +", randomTerm[ii])),
                                 rcov = as.formula(paste("~", spatialTerm[ii])),
                                 aom = TRUE, data = TD, ...)
          } else {
            mr <- asreml::asreml(fixed = as.formula(paste(trait0, "~", 1,
                                                          if (covT) paste(c("", covariate),
                                                                         collapse = "+"))),
                                 random = as.formula(paste("~ genotype +", randomTerm[ii])),
                                 rcov = as.formula(paste("~", spatialTerm[ii])),
                                 aom = TRUE, data = TD, ...)
          }
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
      if (trySpatial == "always") {
        if (tryRep) {
          randomChoice <- c(rep(x = c("Identity", "Measurement_error"), each = 3),
                           paste(rep, row, sep = ":"), paste(rep, col, sep = ":"),
                           paste(paste(rep, row, sep = ":"), paste(rep, col, sep = ":"),
                                 sep = "+"),
                           paste(paste(rep, row, sep = ":"), "Measurement_error", sep = "+"),
                           paste(paste(rep, col, sep = ":"), "Measurement_error", sep = "+"),
                           paste(paste(rep, row, sep = ":"), paste(rep, col, sep = ":"),
                                 "Measurement_error", sep = "+"))
          randomTerm <- c(rep(x = c("NULL", "units"), each = 3),
                          paste(rep, row, sep = ":"), paste(rep, col, sep = ":"),
                          paste(paste(rep, row, sep = ":"), paste(rep, col, sep = ":"),
                                sep = "+"),
                          paste(paste(rep, row, sep = ":"), "units", sep = "+"),
                          paste(paste(rep, col, sep = ":"), "units", sep = "+"),
                          paste(paste(rep, row, sep = ":"), paste(rep, col, sep = ":"),
                                "units", sep = "+"))
        } else {
          randomChoice <- c(rep(x = c("Identity", "Measurement_error"), each = 3), row, col,
                           paste(row, col, sep = "+"), paste(row, "Measurement_error", sep = "+"),
                           paste(col, "Measurement_error", sep = "+"),
                           paste(row, col, "Measurement_error", sep = "+"))
          randomTerm <- c(rep(x = c("NULL", "units"), each = 3),
                          row, col, paste(row, col, sep = "+"),
                          paste(row, "units", sep = "+"), paste(col, "units", sep = "+"),
                          paste(row, col, "units", sep = "+"))
        }
        spatialChoice <- rep(c("Exponential(x)Identity", "Identity(x)Exponential",
                               "Isotropic exponential"), times = 4)
        spatialTerm <- rep(c(paste("exp(", rowCoordinates, "):", colCoordinates),
                             paste(rowCoordinates, ":exp(", colCoordinates, ")"),
                             paste("iexp(", rowCoordinates, ",", colCoordinates, ")")),
                           times = 4)
        modelChoice <- paste("Random:", randomChoice, "&   Spatial:", spatialChoice)
        for (ii in 1:length(randomTerm)) {
          if (tryRep) {
            if (!is.na(checkId)) {
              mr <- asreml::asreml(fixed = as.formula(paste(trait0, "~", rep, "+", checkId,
                                                            if (covT) paste(c("", covariate),
                                                                            collapse = "+"))),
                                   random = as.formula(paste("~ genotype +", randomTerm[ii])),
                                   rcov = as.formula(paste("~", spatialTerm[ii])),
                                   aom = TRUE, data = TD, ...)
            } else {
              mr <- asreml::asreml(fixed = as.formula(paste(trait0, "~", rep,
                                                            if (covT) paste(c("", covariate),
                                                                            collapse = "+"))),
                                   random = as.formula(paste("~ genotype +", randomTerm[ii])),
                                   rcov = as.formula(paste("~", spatialTerm[ii])),
                                   aom = TRUE, data = TD, ...)
            }
          } else {
            if (!is.na(checkId)) {
              mr <- asreml::asreml(fixed = as.formula(paste(trait0, "~", checkId,
                                                            if (covT) paste(c("", covariate),
                                                                            collapse = "+"))),
                                   random = as.formula(paste("~ genotype +", randomTerm[ii])),
                                   rcov = as.formula(paste("~", spatialTerm[ii])),
                                   aom = TRUE, data = TD, ...)
            } else {
              mr <- asreml::asreml(fixed = as.formula(paste(trait0, "~", 1,
                                                            if (covT) paste(c("", covariate),
                                                                            collapse = "+"))),
                                   random = as.formula(paste("~ genotype +", randomTerm[ii])),
                                   rcov = as.formula(paste("~", spatialTerm[ii])),
                                   aom = TRUE, data = TD, ...)
            }
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
    # constrain variance of the variance components to be fixed as the values in mr
    if (tryRep) {
      GParamTmp = bestModel$G.param
      tmpPos <- which(names(GParamTmp) == paste(rep, row, sep = ":"))
      if (length(tmpPos) > 0) {
        if (!is.null(GParamTmp[[tmpPos]][[rep]]$con)) {
          GParamTmp[[tmpPos]][[rep]]$con <- "F"
        }
      }
      tmpPos <- which(names(GParamTmp) == paste(rep, col, sep = ":"))
      if (length(tmpPos) > 0) {
        if (!is.null(GParamTmp[[tmpPos]][[rep]]$con)) {
          GParamTmp[[tmpPos]][[rep]]$con <- "F"
        }
      }
      if (!is.na(checkId)) {
        mf <- asreml::asreml(fixed = as.formula(paste(trait0, "~", rep, "+", checkId,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"),
                                                      "+ genotype")),
                             random = as.formula(paste("~", randomTerm[bestLoc])),
                             rcov = as.formula(paste("~", spatialTerm[ii])),
                             G.param = GParamTmp, aom = TRUE, data = TD, ...)
      } else {
        mf <- asreml::asreml(fixed = as.formula(paste(trait0, "~", rep,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"),
                                                      "+ genotype")),
                             random = as.formula(paste("~", randomTerm[bestLoc])),
                             rcov = as.formula(paste("~", spatialTerm[ii])),
                             G.param = GParamTmp, aom = TRUE, data = TD, ...)
      }
    } else {
      GParamTmp <- bestModel$G.param
      if (!is.null(GParamTmp[[row]][[row]]$con)) {
        GParamTmp[[row]][[row]]$con <- "F"
      }
      if (!is.null(GParamTmp[[col]][[col]]$con)) {
        GParamTmp[[col]][[col]]$con <- "F"
      }
      if (!is.na(checkId)) {
        mf <- asreml::asreml(fixed = as.formula(paste(trait0, "~", checkId,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"),
                                                      "+ genotype")),
                             random = as.formula(paste("~", randomTerm[bestLoc])),
                             rcov = as.formula(paste("~", spatialTerm[ii])),
                             G.param = GParamTmp, aom = TRUE, data = TD, ...)
      } else {
        mf <- asreml::asreml(fixed = as.formula(paste(trait0, "~1",
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"),
                                                      "+ genotype")),
                             random = as.formula(paste("~", randomTerm[bestLoc])),
                             rcov = as.formula(paste("~", spatialTerm[ii])),
                             G.param = GParamTmp, aom = TRUE, data = TD, ...)
      }
    }
  }
  # run predict
  if (!is.na(trySpatial)) {
    ii <- bestLoc
  }
  bestModel$call$fixed <- eval(bestModel$call$fixed)
  bestModel$call$random <- eval(bestModel$call$random)
  bestModel$call$rcov <- eval(bestModel$call$rcov)
  mf$call$fixed <- eval(mf$call$fixed)
  mf$call$random <- eval(mf$call$random)
  mf$call$rcov <- eval(mf$call$rcov)
  bestModel = predict(bestModel, classify = "genotype", data = TD)
  if (!is.na(checkId)) {
    mf <- predict(mf, classify = "genotype", vcov = TRUE, data = TD,
                  associate = ~as.formula(paste0(checkId,":genotype")))
  } else {
    mf <- predict(mf, classify = "genotype", vcov = TRUE, data = TD)
  }
  sink()
  unlink(tmp)
  res = createSSA(mMix = bestModel, mFix = mf, data = TD, trait = trait,
                  genotype = "genotype", rep = ifelse(tryRep, rep, NULL),
                  engine = "asreml")
  return(res)
}
