#' Compute BLUPS and based on a set of mega-environments
#'
#' This function calculates predicted means (BLUPS) and associated standard
#' errors based on a set of mega-environments.
#'
#' @inheritParams gxeAmmi
#'
#' @param useYear Should year be used for modelling (as years within
#' trials). If \code{TRUE TD} should contain a column "year".
#' @param engine A character string specifying the engine used for modeling.
#' Either "lme4" or "asreml".
#' @param ... Further parameters passed to either \code{asreml} or \code{lmer}.
#'
#' @return A list consisting of two data.frames, \code{predictedValue}
#' containing BLUPs per genotype per mega-environment and \code{standardError}
#' containing standard errors for those BLUPs.
#'
#' @examples
#' ## Compute mega-environments for TDMaize.
#' TDMegaEnv <- gxeMegaEnv(TD = TDMaize, trait = "yld", sumTab = FALSE)
#' ## Compute BLUPS and standard errors for those mega-environments.
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
  if (!is.character(trials) || !all(trials %in% names(TD))) {
    stop("All trials should be in TD.")
  }
  TDTot <- Reduce(f = rbind, x = TD[trials])
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TDTot)) {
    stop("trait has to be a column in TD.\n")
  }
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TDTot)) {
    stop("trait has to be a column in TD.\n")
  }
  if (useYear && !"year" %in% colnames(TDTot)) {
    stop("year has to be a column in TD.\n")
  }
  TDTot$megaEnv <- droplevels(TDTot$megaEnv)
  engine <- match.arg(engine)
  if (engine == "asreml") {
    if (requireNamespace("asreml", quietly = TRUE)) {
      ## Create tempfile to suppress asreml output messages.
      tmp <- tempfile()
      sink(file = tmp)
      if (useYear) {
        mr <- tryCatchExt(asreml::asreml(fixed = as.formula(paste(trait,
                                                                  "~ trial / year")),
                                 random = as.formula("~ genotype:us(megaEnv) +
                                                     genotype:megaEnv:year"),
                                 data = TDTot, ...))
      } else {
        mr <- tryCatchExt(asreml::asreml(fixed = as.formula(paste(trait,
                                                                  "~ trial")),
                                 random = as.formula("~ genotype:us(megaEnv)"),
                                 data = TDTot, ...))
      }
      sink()
      unlink(tmp)
      if (!is.null(mr$warning) || !is.null(mr$error)) {
        ## In case asreml gave an error return a data.frame with NAs.
        warning("asreml gave an error. Empty data.frame returned.\n")
        predVals <- se <- data.frame(matrix(nrow = nlevels(TDTot$genotype),
                                            ncol = nlevels(TDTot$megaEnv),
                                            dimnames = list(levels(TDTot$genotype),
                                                            levels(TDTot$megaEnv))),
                                     check.names = FALSE)
      } else {
        mr <-  mr$value
        ## Eval of fixed is needed since fixed contains a variable term trait.
        mr$call$fixed <- eval(mr$call$fixed)
        ## Predict and extract BLUPs
        mr <- predictAsreml(model = mr, classify = "genotype:megaEnv", TD = TDTot)
        pVal <- mr$predictions$pvals
        ## If megaEnv consists of numeric values megaEnv will be a numeric
        ## column in pVal instead of the expected factor. This causes a
        ## shift in column order. Therefore set it back to factor.
        pVal$megaEnv <- factor(pVal$megaEnv, levels = levels(TDTot$megaEnv))
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
      mr <- try(lme4::lmer(as.formula(paste(trait, "~ trial / year +
                                            (0 + megaEnv | genotype) +
                                            (0 + megaEnv | genotype:year)")),
                           data = TDTot, ...), silent = TRUE)
    } else {
      mr <- try(lme4::lmer(as.formula(paste(trait, "~ trial +
                                            (0 + megaEnv | genotype)")),
                           data = TDTot, ...), silent = TRUE)
    }
    if (inherits(mr, "try-error")) {
      warning("lme4 gave an error. Empty data.frame returned.\n")
      predVals <- se <- data.frame(matrix(nrow = nlevels(TDTot$genotype),
                                          ncol = nlevels(TDTot$megaEnv),
                                          dimnames = list(levels(TDTot$genotype),
                                                          levels(TDTot$megaEnv))),
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
      colnames(se) <- colnames(predVals) <- levels(TDTot$megaEnv)
    }
  }
  return(list(predictedValue = predVals, standardError = se))
}
