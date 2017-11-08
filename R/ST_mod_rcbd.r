#'  Single trial (ST) modelling for a randomized complete block design (rcbd)
#'
#' Phenotypic traits are analysed by fitting mixed models to obtain estimates of
#' genotypic means and several genetic parameters. Two mixed models are fitted; the first
#' fits genotypes as a random factor to obtain genetic variance components; and the second
#' fits genotypes as fixed to obtain estimates of genotypic means.
#' @param Y A data frame object.
#' @param trait A string specifying the trait name.
#' @param covariate A string specifying (a) covariate name(s).
#' @param genotype A string specifying the column name of the genotypes.
#' @param rep A string specifying the column name of the replicates.
#' @param checkId (Optional) a string specifying the column name of the check ID(s).
#' @param engine A string specifying the name of the mixed modelling engine to use.
#' @param ... Further arguments to be passed to either \code{asreml} or \code{lme4}.
#'
#' @return an object of class \code{\link{SSA}}.
#'
#' @examples
#' mydat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames=c("Env","Genotype","Rep"),
#'                      traitNames="yield", env ="Env",rowSelect="HEAT05",
#'                      colSelect=c("Env","Genotype","Rep","yield"))
#' mymodel <- ST.mod.rcbd(Y=mydat, trait="yield",
#'                        genotype="Genotype", rep="Rep",
#'                        engine="lme4") #engine="asreml"
#' summary(mymodel)
#'
#' @export
ST.mod.rcbd <- function(Y,
                        trait,
                        covariate,
                        genotype,
                        rep,
                        checkId,
                        engine,
                        ...) {
  # any check ID
  if (missing(checkId)) {
    checks <- FALSE
    iNames <-c(trait, genotype, rep)
  } else {
    checks <- checkId %in% colnames(Y)
    iNames <- c(trait, genotype, rep, checkId)
  }
  # any covariate
  covT <- FALSE
  if (!missing(covariate)) {
    if (is.character(covariate)) {
      covT <- TRUE
      iNames <- c(iNames, covariate)
    }
  }
  #check validility of column names of Y
  YNames <- names(Y)
  if (all(iNames %in% YNames)) {
    vNameTest <- isValidVariableName(iNames)
    if (!all(vNameTest)) {
      warning(paste(iNames[!vNameTest], collapse = ",")," not syntactically valid name(s).\n")
    }
  } else {
    stop(paste(iNames[!(iNames %in% YNames)], collapse = ","), " not found in the names of Y.\n")
  }
  if (engine == "asreml") {
    # Run mixed and fixed models using asreml
    if (checks) {
      mr <- asreml::asreml(fixed = as.formula(paste(trait, "~", checkId,
                                                    if (covT) paste(c("", covariate), collapse = "+"))),
                           random = as.formula(paste("~", rep, "+", genotype)),
                           rcov = ~units, aom = TRUE, data = Y, ...)
    } else {
      mr <- asreml::asreml(fixed = as.formula(paste(trait, "~ 1",
                                                    if (covT) paste(c("", covariate), collapse = "+"))),
                           random = as.formula(paste("~", rep, "+", genotype)),
                           rcov = ~units, aom = TRUE, data = Y, ...)
    }
    # constrain variance of the variance components to be fixed as the values in mr
    GParamTmp <- mr$G.param
    GParamTmp[[rep]][[rep]]$con <- "F"
    if (checks) {
      mf <- asreml::asreml(fixed = as.formula(paste(trait, "~", checkId,
                                                    if (covT) paste(c("", covariate), collapse = "+"),
                                                    "+", genotype)),
                           random = as.formula(paste0("~", rep)), rcov = ~ units,
                           G.param = GParamTmp, aom = TRUE, data = Y, ...)
    } else {
      mf <- asreml::asreml(fixed = as.formula(paste(trait, "~1",
                                                    if (covT) paste(c("", covariate), collapse = "+"),
                                                    "+", genotype)),
                           random = as.formula(paste0("~", rep)), rcov = ~units,
                           G.param = GParamTmp, aom = TRUE, data = Y, ...)
    }
    # run predict
    mr$call$fixed <- eval(mr$call$fixed)
    mr$call$random <- eval(mr$call$random)
    mr$call$rcov <- eval(mr$call$rcov)
    mf$call$fixed <- eval(mf$call$fixed)
    mf$call$random <- eval(mf$call$random)
    mf$call$rcov <- eval(mf$call$rcov)
    tmp <- tempfile()
    sink(file = tmp)
    if (checks) {
      mf <- predict(mf, classify = genotype, vcov = TRUE, data = Y, associate = ~checkId:genotype)
    } else{
      mf <- predict(mf, classify = genotype, vcov = TRUE, data = Y)
    }
    mr <- predict(mr, classify = genotype, data = Y)
    mf$call$data <- substitute(Y)
    mr$call$data <- substitute(Y)
    sink()
    unlink(tmp)
  } else if (engine == "lme4") {
    # Run mixed and fixed models using lme4
    if (checks) {
      frm <- as.formula(paste(trait, "~ (1 |", rep, ")+", checkId,
                              if(covT) paste(c("", covariate), collapse = "+"),
                              "+ (1 |", genotype, ")"))
    } else {
      frm <- as.formula(paste(trait, "~ (1 |", rep, ")+ (1 |", genotype, ")"))
    }
    mr <- lme4::lmer(frm, data = Y, ...)
    if (checks) {
      mf <- lm(as.formula(paste(trait, "~", rep, "+", checkId,
                                if(covT) paste(c("", covariate), collapse = "+"),
                                "+", genotype)), data = Y, ...)
    } else {
      mf <- lm(as.formula(paste(trait, "~", rep,
                                if(covT) paste(c("", covariate), collapse = "+"),
                                "+", genotype)), data = Y, ...)
    }
  } else {
    stop("Please use either asreml or lme4 for engine")
  }
  model = createSSA(mMix = mr, mFix = mf, data = Y)
  attr(model, "Trait") <- trait
  attr(model, "Design") <- "rcbd"
  attr(model, "Engine") <- engine
  attr(model, "genotype") <- genotype
  attr(model, "rep") <- rep
  return(model)
}
