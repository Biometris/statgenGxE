#' S3 class megaEnv
#'
#' Function for creating objects of S3 class megaEnv.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param TD A data.frame containing values for the
#' cultivar-superiority measure of Lin and Binns.
#' @param summTab A data.frame containing values for Shukla's stability variance.
#' @param trait A character string indicating the trait that has been analyzed.
#'
#' @seealso \code{\link{predict.megaEnv}}
#'
#' @name megaEnv
NULL

#' @rdname megaEnv
#' @keywords internal
createMegaEnv <- function(TD,
                          summTab,
                          trait) {
  megaEnv <- structure(list(TD = TD,
                            summTab = summTab,
                            trait = trait),
                       class = "megaEnv")
  attr(megaEnv, which = "timestamp") <- Sys.time()
  return(megaEnv)
}

#' @export
print.megaEnv <- function(x,
                          ...) {
  cat("Mega environments based on ", x$trait, "\n\n", sep = "")
  print(x$summTab, row.names = FALSE)
}

#' @export
summary.megaEnv <- function(object,
                            ...) {
  print(object, ...)
}

#' Plot function for class megaEnv
#'
#' Function for creating scatter plots of predicted values in computed mega
#' environments.
#'
#' @param x An object of class megaEnv.
#' @param ... Further arguments to be passed on to underlying plot functions.
#' @param engine A character string specifying the engine used for making the
#' predictions on which the plots are based.
#' @param colorGenoBy A character string indicating a column in \code{TD} by
#' which the genotypes in the scatter plots are colored. If \code{NULL} all
#' genotypes are displayed in black.
#' @param title A character string used a title for the plot.
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a ggtable object is invisibly returned.
#'
#' @examples
#' ## Compute mega environments for TDMaize.
#' geMegaEnv <- gxeMegaEnv(TD = TDMaize, trait = "yld")
#'
#' ## Create a scatter plot of predicted values.
#' plot(geMegaEnv)
#'
#' @family mega environments
#'
#' @export
plot.megaEnv <- function(x,
                         ...,
                         engine = c("lme4", "asreml"),
                         colorGenoBy = NULL,
                         title = paste("Scatterplot of mega environments for",
                                       x$trait),
                         output = TRUE) {
  chkChar(title, len = 1, null = FALSE)
  engine <- match.arg(engine)
  if (!is.null(colorGenoBy)) {
    TDTot <- do.call(rbind, args = x$TD)
    chkCol(column = colorGenoBy, obj = TDTot)
  }
  pred <- predict(x, engine = engine)$predictedValue
  predLong <- reshape(pred, direction = "long",
                      varying = list(megaEnv = colnames(pred)),
                      ids = rownames(pred), idvar = "genotype",
                      timevar = "megaEnv", v.names = "pred")
  predLong[["megaEnv"]] <- factor(predLong[["megaEnv"]],
                                  labels = colnames(pred))
  if (!is.null(colorGenoBy)) {
    predLong <- merge(predLong, TDTot[c("genotype", colorGenoBy)])
  }
  predTD <- createTD(predLong, genotype = "genotype", trial = "megaEnv")
  p <- plot(predTD, plotType = "scatter", traits = "pred",
            colorGenoBy = colorGenoBy, title = title, output = output)
  invisible(p)
}

#' Compute BLUPS based on a set of mega environments
#'
#' This function calculates Best Linear Unbiased Predictors (BLUPS) and
#' associated standard errors based on a set of mega environments.
#'
#' @inheritParams gxeAmmi
#'
#' @param object An object of class megaEnv.
#' @param useYear Should year be used for modeling (as years within
#' trials). If \code{TRUE}, \code{TD} should contain a column "year".
#' @param engine A character string specifying the engine used for modeling.
#' @param ... Further parameters passed to either \code{asreml} or \code{lmer}.
#'
#' @returns A list consisting of two data.frames, \code{predictedValue}
#' containing BLUPs per genotype per mega environment and \code{standardError}
#' containing standard errors for those BLUPs.
#'
#' @examples
#' ## Compute mega environments for TDMaize.
#' geMegaEnv <- gxeMegaEnv(TD = TDMaize, trait = "yld")
#'
#' ## Compute BLUPS and standard errors for those mega environments.
#' megaEnvPred <- predict(geMegaEnv)
#' head(megaEnvPred$predictedValue)
#' head(megaEnvPred$standardError)
#'
#' @family mega environments
#'
#' @export
predict.megaEnv <- function(object,
                            ...,
                            trials = names(object$TD),
                            useYear = FALSE,
                            engine = c("lme4", "asreml")) {
  ## Checks.
  TD <- object$TD
  trait <- object$trait
  trials <- chkTrials(trials, TD)
  TDTot <- Reduce(f = rbind, x = TD[trials])
  if (useYear) {
    chkCol("year", TDTot)
  }
  engine <- match.arg(engine)
  if (length(trials) < 10) {
    warning("One should be cautious with the interpretation of predictions ",
            "for mega environments that are based on less than 10 trials.\n")
  }
  TDTot <- droplevels(TDTot)
  if (engine == "asreml") {
    if (requireNamespace("asreml", quietly = TRUE)) {
      if (useYear) {
        fixForm <- formula(paste0("`", trait, "`~ trial / year"))
        randForm <- formula("~ genotype:us(megaEnv) + genotype:megaEnv:year")
      } else {
        fixForm <- formula(paste0("`", trait, "`~ trial"))
        randForm <- formula("~ genotype:us(megaEnv)")
      }
      mr <- tryCatchExt(asreml::asreml(fixed = fixForm, random = randForm,
                                       data = TDTot, maxiter = 200,
                                       trace = FALSE,
                                       workspace = 160e7, ...))
      if (!is.null(mr$error)  || (!is.null(mr$warning) && !mr$value$converge)) {
        ## In case asreml gave an error return a data.frame with NAs.
        warning("asreml gave an error. Empty data.frame returned.\n")
        predVals <- se <-
          data.frame(matrix(nrow = nlevels(TDTot[["genotype"]]),
                            ncol = nlevels(TDTot[["megaEnv"]]),
                            dimnames = list(levels(TDTot[["genotype"]]),
                                            levels(TDTot[["megaEnv"]]))),
                     check.names = FALSE)
      } else {
        mr <- mr$value
        ## Eval of fixed and random is needed for predictions.
        mr$call$fixed <- eval(mr$call$fixed)
        mr$call$random <- eval(mr$call$random)
        ## Predict and extract BLUPs
        mr <- suppressWarnings(predictAsreml(model = mr,
                                             classify = "genotype:megaEnv",
                                             TD = TDTot))
        if (asreml4()) {
          pVal <- mr$pvals
        } else {
          pVal <- mr$predictions$pvals
        }
        ## If megaEnv consists of numerical values megaEnv will be a numerical
        ## column in pVal instead of the expected factor. This causes a
        ## shift in column order. Therefore set it back to factor.
        pVal[["megaEnv"]] <- factor(pVal[["megaEnv"]],
                                    levels = levels(TDTot[["megaEnv"]]))
        ## Reshape to a data.frame with genotypes in rows and megaEnv in cols.
        predVals <-
          as.data.frame(tapply(X = pVal[["predicted.value"]],
                               INDEX = pVal[, c("genotype", "megaEnv")],
                               FUN = I))
        se <- as.data.frame(tapply(X = pVal[[if (asreml4()) "std.error" else
          "standard.error"]], INDEX = pVal[, c("genotype", "megaEnv")],
          FUN = I))
      }
    } else {
      stop("Failed to load asreml.\n")
    }
  } else if (engine == "lme4") {
    if (useYear) {
      mr <- try(lme4::lmer(formula(paste0("`", trait, "`~ trial / year + ",
                                          "(0 + megaEnv | genotype) + ",
                                          "(0 + megaEnv | genotype:year)")),
                           data = TDTot, ...), silent = TRUE)
    } else {
      mr <- try(lme4::lmer(formula(paste0("`", trait, "`~ trial + ",
                                          "(0 + megaEnv | genotype)")),
                           data = TDTot, ...), silent = TRUE)
    }
    if (inherits(mr, "try-error")) {
      warning("lme4 gave an error. Empty data.frame returned.\n")
      predVals <- se <-
        data.frame(matrix(nrow = nlevels(TDTot[["genotype"]]),
                          ncol = nlevels(TDTot[["megaEnv"]]),
                          dimnames = list(levels(TDTot[["genotype"]]),
                                          levels(TDTot[["megaEnv"]]))),
                   check.names = FALSE)
    } else {
      ## Extract fixed effects needed to compute intercept.
      fixEff = lme4::fixef(mr)
      fixEf = fixEff[grep("trial", names(fixEff))]
      ## Extract random effects for genotypes.
      ranEff = lme4::ranef(mr, drop = TRUE)[["genotype"]]
      ## Compute BLUPs.
      predVals = fixEff[1] + mean(c(fixEf, 0)) + ranEff
      ## Extract seBlups
      predErrs <- attr(lme4::ranef(mr, condVar = TRUE)[["genotype"]],
                       "postVar")
      ## Reshape to a data.frame with genotypes in rows and megaEnv in cols.
      seBlups = t(sqrt(apply(X = predErrs, MARGIN = 3, FUN = diag)))
      se <- data.frame(seBlups, row.names = rownames(predVals),
                       check.names = FALSE)
      colnames(se) <- colnames(predVals) <- levels(TDTot[["megaEnv"]])
    }
  }
  return(list(predictedValue = predVals, standardError = se))
}
