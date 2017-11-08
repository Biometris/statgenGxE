#' Forms predicted means (BLUPS) and produces tables based on a set of mega-environments
#'
#' This function calculates predicted meanes (BLUPS) and associated standard errors based on a
#' set of mega-environments
#'
#' @param trait A character string specifying a trait column of the data.
#' @param genotype A character string specifying a genotype column of the data.
#' @param env A character string specifying an environment column of the data.
#' @param year A character string specifying years within environments.
#' @param megaEnv A character string specifying the mega-environment factor.
#' @param data A data frame object.
#' @param ... Other parameters passed to either \code{asreml()} or \code{lmer()}.
#'
#' @examples
#' mydat <- GE.read.csv(system.file("extdata", "F2maize_pheno.csv", package = "RAP"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' names(mydat)=c("env", "genotype","yld")
#' Y <- GE.megaEnvironment(Y=mydat, trait="yld", genotype="genotype",
#'                         env="env", megaEnv="megaEnv")
#' GE.GxETable(trait="yld", genotype="genotype", env="env", year=NULL,
#'             megaEnv="megaEnv", data=Y)
#'
#' @import utils
#' @importFrom methods slot
#' @export

GE.GxETable <- function(trait,
                        genotype,
                        env,
                        year = NULL,
                        megaEnv,
                        data,
                        ...) {
  if (!trait %in% names(data)) {
    stop(trait," not found in ", data)
  }
  if (!genotype %in% names(data)) {
    stop(genotype, " not found in ", data)
  }
  if (!env %in% names(data)) {
    stop(env, " not found in ", data)
  }
  if (!megaEnv %in% names(data)) {
    stop(megaEnv, " not found in ", data)
  }
  if (asremlORlme4() == "asreml"){
    if (is.null(year)){
      #      sv <- asreml::asreml(fixed=as.formula(paste(trait,"~",env)),
      #                   random <- as.formula(paste("~",genotype,":us(",megaEnv,")")),
      #                   data=data, start.values=T,...)
      #      sv$R.param$R$variance$con <- "F"
      mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~", env)),
                               random = as.formula(paste("~", genotype, ":us(", megaEnv, ")")),
                               data = data, ...),
                silent = TRUE)
    } else{
      #      sv <- asreml::asreml(fixed=as.formula(paste(trait,"~",env,"/",year)),
      #                   random <- as.formula(paste("~",genotype,":us(",megaEnv,")+",genotype,":",megaEnv,":",year)),
      #                   data=data, start.values=T,...)
      #      sv$R.param$R$variance$con <- "F"
      mr <- try(asreml::asreml(fixed = as.formula(paste(trait, "~", env, "/", year)),
                               random = as.formula(paste("~", genotype, ":us(", megaEnv, ")+",
                                                         genotype, ":", megaEnv, ":", year)),
                               data = data, ...),
                silent = TRUE)
    }
    if (inherits(mr, "try-error")){
      genoLevels <- levels(data[[genotype]])
      nGenoUnique <- nlevels(data[[genotype]])
      megaEnvLevels <- levels(data[[megaEnv]])
      nmegaEnvUnique <- nlevels(data[[megaEnv]])
      predVals <- se <-  data.frame(matrix(nrow = nGenoUnique, ncol = nmegaEnvUnique),
                                    row.names = genoLevels, check.names = FALSE)
      names(predVals) <- names(se) <- megaEnvLevels
    } else {
      mr$call$fixed <- eval(mr$call$fixed)
      mr$call$random <- eval(mr$call$random)
      mr$call$rcov <- eval(mr$call$rcov)
      mr$call$R.param <- eval(mr$call$R.param)
      mr <- predict(mr, classify = paste0(genotype, ":", megaEnv), data = data)
      predictions <- mr$predictions$pvals
      if (!is.factor(predictions[[megaEnv]])) {
        predictions[[megaEnv]] <- as.factor(predictions[[megaEnv]])
      }
      if (!is.factor(predictions[[genotype]])) {
        predictions[[genotype]] <- as.factor(predictions[[genotype]])
      }
      predVals <- tapply(X = predictions$predicted.value,
                         INDEX = predictions[, c(genotype, megaEnv)], FUN = identity)
      se <- tapply(X = predictions$standard.error,
                   INDEX = predictions[, c(genotype, megaEnv)], FUN = identity)
    }
  } else if (asremlORlme4() == "lme4") {
    if (is.null(year)) {
      mr <- try(lme4::lmer(as.formula(paste(trait, "~", env, "+ (0+", megaEnv,
                                            "|", genotype, ")")), data = data, ...),
                silent = TRUE)
    } else {
      mr <- try(lme4::lmer(as.formula(paste(trait, "~", env, "/", year, "+ (0+",
                                            megaEnv, "|", genotype, ") + (0+", megaEnv,
                                            "|", genotype, ":" , year, ")")),
                           data = data, ...), silent = TRUE)
    }
    genoLevels <- levels(data[[genotype]])
    nGenoUnique <- nlevels(data[[genotype]])
    megaEnvLevels <- levels(data[[megaEnv]])
    nmegaEnvUnique <- nlevels(data[[megaEnv]])
    if (inherits(mr, "try-error")){
      predVals <- se <-  data.frame(matrix(nrow = nGenoUnique, ncol = nmegaEnvUnique),
                                    row.names = genoLevels, check.names = FALSE)
      names(predVals) <- names(se) <- megaEnvLevels
    } else {
      # Extract coeffcients mr
      if (class(mr) == 'lmerMod') {
        fe = lme4::fixef(mr)
      }
      if(class(mr) == 'mer') {
        fe = slot(mr, "fixef")
      }
      rr = grep(env, names(fe))
      cr = fe[rr]
      ng = length(unique(slot(mr, "flist")[[genotype]]))
      reff = lme4::ranef(mr, drop = TRUE)[[genotype]]
      blo = mean(c(cr, 0))
      # Predictions BLUPs
      predVals = fe[1] + blo + reff
      # Compute seBlups
      if(class(mr) == 'lmerMod') {
        seBlups = t(sqrt(apply(X = attr(lme4::ranef(mr, condVar = TRUE)[[genotype]], "postVar"),
                               MARGIN = 3, FUN = diag)))
      }
      if(class(mr) == 'mer') {
        seBlups = t(sqrt(apply(X = attr(lme4::ranef(mr, postVar=TRUE)[[genotype]], "postVar"),
                               MARGIN = 3, FUN = diag)))
      }
      se <- data.frame(seBlups, row.names = rownames(predVals), check.names = FALSE)
      names(se) <- names(predVals) <- megaEnvLevels
    }
  } else{
    stop("Either asreml or lme4 is not loaded correctly.")
  }
  result <- vector(mode = "list")
  return(list(predictedValue = predVals, standardError = se))
}
