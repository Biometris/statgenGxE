#' Forms predicted means (BLUPS) and produces tables based on a set of mega-environments
#'
#' This function calculates predicted means (BLUPS) and associated standard errors based on a
#' set of mega-environments
#'
#' @inheritParams gxeAmmi
#'
#' @param year A character string specifying years within environments.
#' @param engine asreml or lme4
#' @param ... Other parameters passed to either \code{asreml()} or \code{lmer()}.
#'
#' @examples
#' data(TDMaize)
#' myTDMegaEnv <- gxeMegaEnvironment(TD = TDMaize, trait = "yld")
#' geTab <- gxeTable(TD = myTDMegaEnv, trait = "yld")
#'
#' @export
gxeTable <- function(TD,
                     trait,
                     year = NULL,
                     engine = "asreml",
                     ...) {
  if (!trait %in% names(TD)) {
    stop(trait," not found in ", TD)
  }
  if (engine == "asreml") {
    tmp <- tempfile()
    sink(file = tmp)
    if (is.null(year)) {
      mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env")),
                               random = as.formula("~ genotype:us(megaEnv)"),
                               data = TD, ...),
                silent = TRUE)
    } else {
      mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~ env /", year)),
                               random = as.formula(paste("~ genotype:us(megaEnv) +
                                                         genotype:megaEnv:", year)),
                               data = TD, ...),
                silent = TRUE)
    }
    sink()
    if (inherits(mr, "try-error")) {
      genoLevels <- levels(TD$genotype)
      nGenoUnique <- nlevels(TD$genotype)
      megaEnvLevels <- levels(TD$megaEnv)
      nmegaEnvUnique <- nlevels(TD$megaEnv)
      predVals <- se <-  data.frame(matrix(nrow = nGenoUnique, ncol = nmegaEnvUnique),
                                    row.names = genoLevels, check.names = FALSE)
      names(predVals) <- names(se) <- megaEnvLevels
    } else {
      mr$call$fixed <- eval(mr$call$fixed)
      mr$call$random <- eval(mr$call$random)
      mr$call$rcov <- eval(mr$call$rcov)
      mr$call$R.param <- eval(mr$call$R.param)
      mr <- predictAsreml(model = mr, classify = "genotype:megaEnv",
                          TD = TD)
      predictions <- mr$predictions$pvals
      predVals <- tapply(X = predictions$predicted.value,
                         INDEX = predictions[, c("genotype", "megaEnv")],
                         FUN = identity)
      se <- tapply(X = predictions$standard.error,
                   INDEX = predictions[, c("genotype", "megaEnv")],
                   FUN = identity)
    }
    unlink(tmp)
  } else if (engine == "lme4") {
    if (is.null(year)) {
      mr <- try(lme4::lmer(as.formula(paste(trait, "~ env +
                                            (0 + megaEnv | genotype)")),
                           data = TD, ...),
                silent = TRUE)
    } else {
      mr <- try(lme4::lmer(as.formula(paste(trait, "~ env / " , year, "+ (0 +
                                            megaEnv | genotype) + (0 + megaEnv
                                            | genotype:", year, ")")),
                           data = TD, ...), silent = TRUE)
    }
    genoLevels <- levels(TD$genotype)
    nGenoUnique <- nlevels(TD$genotype)
    megaEnvLevels <- levels(TD$megaEnv)
    nmegaEnvUnique <- nlevels(TD$megaEnv)
    if (inherits(mr, "try-error")) {
      predVals <- se <-  data.frame(matrix(nrow = nGenoUnique, ncol = nmegaEnvUnique),
                                    row.names = genoLevels, check.names = FALSE)
      names(predVals) <- names(se) <- megaEnvLevels
    } else {
      # Extract coeffcients mr
      fixEff = lme4::fixef(mr)
      cr = fixEff[grep("env", names(fixEff))]
      ranEff = lme4::ranef(mr, drop = TRUE)[["genotype"]]
      blo = mean(c(cr, 0))
      # Predictions BLUPs
      predVals = fixEff[1] + blo + ranEff
      # Compute seBlups
      seBlups = t(sqrt(apply(X = attr(lme4::ranef(mr, condVar = TRUE)[["genotype"]],
                                      "postVar"),
                             MARGIN = 3, FUN = diag)))
      se <- data.frame(seBlups, row.names = rownames(predVals), check.names = FALSE)
      names(se) <- names(predVals) <- megaEnvLevels
    }
  }
  return(list(predictedValue = predVals, standardError = se))
}
