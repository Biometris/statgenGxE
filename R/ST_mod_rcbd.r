#'  Single trial (ST) modelling for a randomized complete block design (rcbd)
#'
#' Phenotypic traits are analysed by fitting mixed models to obtain estimates of
#' genotypic means and several genetic parameters. Two mixed models are fitted; the first
#' fits genotypes as a random factor to obtain genetic variance components; and the second
#' fits genotypes as fixed to obtain estimates of genotypic means.
#'
#' @inheritParams ST.run.model
#'
#' @param engine A string specifying the name of the mixed modelling engine to use.
#' @param ... Further arguments to be passed to either \code{asreml} or \code{lme4}.
#'
#' @return an object of class \code{\link{SSA}}.
#'
#' @seealso \code{\link{createSSA}}, \code{\link{summary.SSA}}
#'
#' @examples
#' myDat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames = c("Env", "Genotype", "Rep", "Subblock"),
#'                      traitNames = "yield", env = "Env", rowSelect = "HEAT05",
#'                      colSelect = c("Env", "Genotype", "Rep", "Subblock", "yield"))
#' myTD <- createTD(data = myDat, genotype = "Genotype", env = "Env")
#' myModel <- ST.mod.rcbd(TD = myTD, trait = "yield", rep = "Rep",
#'                        engine = "lme4") #engine = "asreml"
#' summary(myModel)
#'
#' @export
ST.mod.rcbd <- function(TD,
                        trait,
                        covariate,
                        rep,
                        checkId,
                        engine,
                        ...) {
  # any check ID
  if (missing(checkId)) {
    checks <- FALSE
    iNames <-c(trait, "genotype", rep)
  } else {
    checks <- checkId %in% colnames(TD)
    iNames <- c(trait, "genotype", rep, checkId)
  }
  # any covariate
  covT <- FALSE
  if (!missing(covariate)) {
    if (is.character(covariate)) {
      covT <- TRUE
      iNames <- c(iNames, covariate)
    }
  }
  #check validility of column names of TD
  TDNames <- names(TD)
  if (all(iNames %in% TDNames)) {
    vNameTest <- isValidVariableName(iNames)
    if (!all(vNameTest)) {
      warning(paste(iNames[!vNameTest], collapse = ",")," not syntactically valid name(s).\n")
    }
  } else {
    stop(paste(iNames[!(iNames %in% TDNames)], collapse = ","), " not found in the names of Y.\n")
  }
  if (engine == "asreml") {
    tmp <- tempfile()
    sink(file = tmp)
    # Run mixed and fixed models using asreml
    if (checks) {
      mr <- asreml::asreml(fixed = as.formula(paste(trait, "~", checkId,
                                                    if (covT) paste(c("", covariate),
                                                                    collapse = "+"))),
                           random = as.formula(paste("~", rep, "+ genotype")),
                           rcov = ~units, aom = TRUE, data = TD, ...)
    } else {
      mr <- asreml::asreml(fixed = as.formula(paste(trait, "~ 1",
                                                    if (covT) paste(c("", covariate),
                                                                    collapse = "+"))),
                           random = as.formula(paste("~", rep, "+ genotype")),
                           rcov = ~units, aom = TRUE, data = TD, ...)
    }
    # constrain variance of the variance components to be fixed as the values in mr
    GParamTmp <- mr$G.param
    GParamTmp[[rep]][[rep]]$con <- "F"
    if (checks) {
      mf <- asreml::asreml(fixed = as.formula(paste(trait, "~", checkId,
                                                    if (covT) paste(c("", covariate),
                                                                    collapse = "+"),
                                                    "+ genotype")),
                           random = as.formula(paste0("~", rep)), rcov = ~ units,
                           G.param = GParamTmp, aom = TRUE, data = TD, ...)
    } else {
      mf <- asreml::asreml(fixed = as.formula(paste(trait, "~1",
                                                    if (covT) paste(c("", covariate),
                                                                    collapse = "+"),
                                                    "+ genotype")),
                           random = as.formula(paste0("~", rep)), rcov = ~units,
                           G.param = GParamTmp, aom = TRUE, data = TD, ...)
    }
    # run predict
    mr$call$fixed <- eval(mr$call$fixed)
    mr$call$random <- eval(mr$call$random)
    mr$call$rcov <- eval(mr$call$rcov)
    mf$call$fixed <- eval(mf$call$fixed)
    mf$call$random <- eval(mf$call$random)
    mf$call$rcov <- eval(mf$call$rcov)
    if (checks) {
      mf <- predict(mf, classify = "genotype", vcov = TRUE, data = TD,
                    associate = as.formula(paste("~", checkId, ":genotype")))
    } else{
      mf <- predict(mf, classify = "genotype", vcov = TRUE, data = TD)
    }
    mr <- predict(mr, classify = "genotype", data = TD)
    mf$call$data <- substitute(TD)
    mr$call$data <- substitute(TD)
    sink()
    unlink(tmp)
  } else if (engine == "lme4") {
    # Run mixed and fixed models using lme4
    if (checks) {
      frm <- as.formula(paste(trait, "~ (1 |", rep, ")+", checkId,
                              if(covT) paste(c("", covariate), collapse = "+"),
                              "+ (1 | genotype)"))
    } else {
      frm <- as.formula(paste(trait, "~ (1 |", rep, ")+ (1 | genotype)"))
    }
    mr <- lme4::lmer(frm, data = TD, ...)
    if (checks) {
      mf <- lm(as.formula(paste(trait, "~", rep, "+", checkId,
                                if(covT) paste(c("", covariate), collapse = "+"),
                                "+ genotype")), data = TD, ...)
    } else {
      mf <- lm(as.formula(paste(trait, "~", rep,
                                if(covT) paste(c("", covariate), collapse = "+"),
                                "+ genotype")), data = TD, ...)
    }
  } else {
    stop("Please use either asreml or lme4 for engine")
  }
  model = createSSA(mMix = mr, mFix = mf, data = TD, trait = trait,
                    genotype = "genotype", rep = rep,
                    design = "rcbd", engine = engine)
  return(model)
}
