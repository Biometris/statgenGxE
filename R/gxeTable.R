#' Compute BLUPS based on a set of mega environments
#'
#' This function calculates Best Lineair Unbiased Predictors (BLUPS) and
#' associated standard errors based on a set of mega environments.
#'
#' @inheritParams gxeAmmi
#'
#' @param useYear Should year be used for modeling (as years within
#' trials). If \code{TRUE}, \code{TD} should contain a column "year".
#' @param engine A character string specifying the engine used for modeling.
#' Either "lme4" or "asreml".
#' @param ... Further parameters passed to either \code{asreml} or \code{lmer}.
#'
#' @return A list consisting of two data.frames, \code{predictedValue}
#' containing BLUPs per genotype per mega environment and \code{standardError}
#' containing standard errors for those BLUPs.
#'
#' @examples
#' ## Compute mega environments for TDMaize.
#' TDMegaEnv <- gxeMegaEnv(TD = TDMaize, trait = "yld", sumTab = FALSE)
#' ## Compute BLUPS and standard errors for those mega environments.
#' geTab <- gxeTable(TD = TDMegaEnv, trait = "yld")
#'
#' @export
gxeTable <- function(TD,
                     trials = names(TD),
                     trait,
                     useYear = FALSE,
                     engine = c("lme4", "asreml"),
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
  chkCol("megaEnv", TDTot)
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
                                       trace = FALSE, ...))
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
        pVal$megaEnv <- factor(pVal[["megaEnv"]],
                               levels = levels(TDTot[["megaEnv"]]))
        ## Reshape to a data.frame with genotypes in rows and megaEnv in cols.
        predVals <-
          as.data.frame(tapply(X = pVal[["predicted.value"]],
                               INDEX = pVal[, c("genotype", "megaEnv")],
                               FUN = I))
        se <-
          as.data.frame(tapply(X = pVal[[if (asreml4()) "std.error" else
            "standard.error"]], INDEX = pVal[, c("genotype", "megaEnv")],
            FUN = identity))
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
      mr <- try(lme4::lmer(as.formula(paste0("`", trait, "`~ trial + ",
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
