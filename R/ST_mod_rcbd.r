#' Single trial (ST) modelling for a randomized complete block design (rcbd)
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
#' myModel <- ST.mod.rcbd(TD = myTD, trait = "yield", repId = "Rep",
#'                        engine = "lme4") #engine = "asreml"
#' summary(myModel)
#'
#' @export
ST.mod.rcbd <- function(TD,
                        trait,
                        covariate = NULL,
                        repId = NULL,
                        checkId = NULL,
                        engine,
                        ...) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TD)) {
    stop("trait has to be a column in TD.\n")
  }
  if (!is.null(covariate) && (!is.character(covariate) ||
                              !(all(covariate %in% colnames(TD))))) {
    stop("covariate have to be a columns in TD.\n")
  }
  for (param in c(repId, checkId)) {
    if (!is.null(param) && (!is.character(param) || length(param) > 1 ||
                            !param %in% colnames(TD))) {
      stop(paste(deparse(param), "has to be NULL or a column in data.\n"))
    }
  }
  if (!is.null(engine) && (!is.character(engine) || length(engine) > 1 ||
                           !engine %in% c("asreml", "lme4"))) {
    stop("engine should be asreml or lme4")
  }
  # any check ID
  if (missing(checkId)) {
    checks <- FALSE
    iNames <- c(trait, "genotype", repId)
  } else {
    checks <- checkId %in% colnames(TD)
    iNames <- c(trait, "genotype", repId, checkId)
  }
  # any covariate
  covT <- FALSE
  if (!is.null(covariate)) {
    if (is.character(covariate)) {
      covT <- TRUE
      iNames <- c(iNames, covariate)
    }
  }
  ## Check validility of column names of TD
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
    ## Run mixed and fixed models using asreml
    fixedFormR <- as.formula(paste(trait, "~",
                                   if (checks) checkId else "1",
                                   if (covT) paste(c("", covariate), collapse = "+")))
    mr <- asreml::asreml(fixed = fixedFormR,
                         random = as.formula(paste("~", repId, "+ genotype")),
                         rcov = ~ units, aom = TRUE, data = TD, ...)
    ## Constrain variance of the variance components to be fixed as the values in mr
    GParamTmp <- mr$G.param
    GParamTmp[[repId]][[repId]]$con <- "F"
    fixedFormF <- as.formula(paste(deparse(fixedFormR), "+ genotype"))
    mf <- asreml::asreml(fixed = fixedFormF,
                         random = as.formula(paste0("~", repId)), rcov = ~ units,
                         G.param = GParamTmp, aom = TRUE, data = TD, ...)
    ## run predict
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
    ## Run mixed and fixed models using lme4
    frm <- as.formula(paste(trait, "~ (1 |", repId, ")",
                            if (checks) paste("+", checkId),
                            if (covT) paste(c("", covariate), collapse = "+"),
                            "+ (1 | genotype)"))
    mr <- lme4::lmer(frm, data = TD, ...)
    ffm <- as.formula(paste(trait, "~", repId,
                            if (checks) paste("+", checkId),
                            if (covT) paste(c("", covariate), collapse = "+"),
                            "+ genotype"))
    mf <- lm(ffm, data = TD, ...)
  }
  model = createSSA(mMix = mr, mFix = mf, data = TD, trait = trait,
                    genotype = "genotype", repId = repId,
                    design = "rcbd", engine = engine)
  return(model)
}
