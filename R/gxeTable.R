#' Compute BLUPS and based on a set of mega-environments
#'
#' This function calculates predicted means (BLUPS) and associated standard
#' errors based on a set of mega-environments.
#'
#' @inheritParams gxeAmmi
#'
#' @param useYear Should year be used for modelling (as years within
#' environments). If \code{TRUE} TD should contain a column "year".
#' @param engine A character string specifying the name of the mixed engine to
#' use, either lme4 or asreml.
#' @param ... Other parameters passed to either \code{asreml} or \code{lmer}.
#'
#' @return A list consisting of two data.frames, \code{predictedValue}
#' containing BLUPs per genotype per mega-environment and \code{standardError}
#' containing standard errors for those BLUPs.
#'
#' @examples
#' ## Not run since asreml is used for modeling
#' \dontrun{
#' TDMegaEnv <- gxeMegaEnv(TD = TDMaize, trait = "yld", sumTab = FALSE)
#' geTab <- gxeTable(TD = TDMegaEnv, trait = "yld")
#' }
#' @export
gxeTable <- function(TD,
                     trait,
                     useYear = FALSE,
                     engine = c("asreml", "lme4"),
                     ...) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TD)) {
    stop("trait has to be a column in TD.\n")
  }
  if (useYear && !"year" %in% colnames(TD)) {
    stop("year has to be a column in TD.\n")
  }
  engine <- match.arg(engine)
  if (engine == "asreml") {
    if (requireNamespace("asreml", quietly = TRUE)) {
      ## Create tempfile to suppress asreml output messages.
      tmp <- tempfile()
      sink(file = tmp)
      if (useYear) {
        mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env / year")),
                                 random = as.formula("~ genotype:us(megaEnv) +
                                                     genotype:megaEnv:year"),
                                 data = TD, ...), silent = TRUE)
      } else {
        mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                                 random = as.formula("~ genotype:us(megaEnv)"),
                                 data = TD, ...), silent = TRUE)
      }
      sink()
      unlink(tmp)
      if (inherits(mr, "try-error")) {
        ## In case asreml gave an error return a data.frame with NAs.
        warning("asreml gave an error. Empty data.frame returned.\n")
        predVals <- se <- data.frame(matrix(nrow = nlevels(TD$genotype),
                                            ncol = nlevels(TD$megaEnv),
                                            dimnames = list(levels(TD$genotype),
                                                            levels(TD$megaEnv))),
                                     check.names = FALSE)
      } else {
        ## Eval of fixed is needed since fixed contains a variable term trait.
        mr$call$fixed <- eval(mr$call$fixed)
        ## Predict and extract BLUPs
        mr <- predictAsreml(model = mr, classify = "genotype:megaEnv", TD = TD)
        pVal <- mr$predictions$pvals
        ## If megaEnv consists of numeric values megaEnv will be a numeric
        ## column in pVal instead of the expected factor. This causes a
        ## shift in column order. Therefore set it back to factor.
        pVal$megaEnv <- factor(pVal$megaEnv, levels = levels(TD$megaEnv))
        ## Reshape to a data.frame with genotypes in rows and megaEnv in cols.
        predVals <- as.data.frame(tapply(X = pVal$predicted.value,
                                         INDEX = pVal[, c("genotype", "megaEnv")],
                                         FUN = identity))
        se <- as.data.frame(tapply(X = pVal$standard.error,
                                   INDEX = pVal[, c("genotype", "megaEnv")],
                                   FUN = identity))
      }
    } else {
      stop("Failed to load asreml.\n")
    }
  } else if (engine == "lme4") {
    if (useYear) {
      mr <- try(lme4::lmer(as.formula(paste(trait, "~ env / year +
                                            (0 + megaEnv | genotype) +
                                            (0 + megaEnv | genotype:year)")),
                           data = TD, ...), silent = TRUE)
    } else {
      mr <- try(lme4::lmer(as.formula(paste(trait, "~ env +
                                            (0 + megaEnv | genotype)")),
                           data = TD, ...), silent = TRUE)
    }
    if (inherits(mr, "try-error")) {
      warning("lme4 gave an error. Empty data.frame returned.\n")
      predVals <- se <- data.frame(matrix(nrow = nlevels(TD$genotype),
                                          ncol = nlevels(TD$megaEnv),
                                          dimnames = list(levels(TD$genotype),
                                                          levels(TD$megaEnv))),
                                   check.names = FALSE)
    } else {
      ## Extract fixed effects needed to compute intercept.
      fixEff = lme4::fixef(mr)
      fixEf = fixEff[grep("env", names(fixEff))]
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
      colnames(se) <- colnames(predVals) <- levels(TD$megaEnv)
    }
  }
  return(list(predictedValue = predVals, standardError = se))
}
