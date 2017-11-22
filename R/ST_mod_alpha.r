#' Single trial (ST) modelling for an incomplete-block design (ibd) or resolvable
#' incomplete-block design (res.ibd)
#'
#' Phenotypic traits are analysed by fitting mixed models to obtain estimates of
#' genotypic means and several genetic parameters. Two mixed models are fitted; the first
#' fits genotypes as a random factor to obtain genetic variance components; and the second
#' fits genotypes as fixed to obtain estimates of genotypic means.
#'
#' @inheritParams ST.run.model
#'
#' @param subDesign A string specifying whether to analyse an incomplete-block design (ibd) or
#' resolvable incomplete-block design (res.ibd).
#' @param engine A string specifying the name of the mixed modelling engine to use, either
#' asreml or lme4.
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
#' myModel <- ST.mod.alpha(TD = myTD, subDesign = "res.ibd", trait = "yield",
#'                         rep = "Rep", subBlock = "Subblock",
#'                         engine = "lme4") #engine = "asreml"
#' summary(myModel)
#'
#' @export
ST.mod.alpha <- function(TD,
                         trait,
                         covariate = NULL,
                         rep = NULL,
                         subBlock = NULL,
                         checkId = NULL,
                         subDesign = NULL,
                         engine,
                         ...) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  subDesigns <- c("ibd", "res.ibd")
  if ((is.null(subDesign) && !attr(TD, "subDesign") %in% subDesigns) ||
      (!is.null(subDesign) && (!is.character(subDesign) || length(subDesign) > 1 ||
                               !subDesign %in% subDesigns))) {
    stop("subDesign should either be an attribute of TD or one of ibd or res.ibd.\n")
  }
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TD)) {
    stop("trait has to be a column in TD.\n")
  }
  if (!is.null(covariate) && (!is.character(covariate) ||
                              !(all(covariate %in% colnames(TD))))) {
    stop("covariate have to be a columns in TD.\n")
  }
  for (param in c(rep, subBlock, checkId)) {
    if (!is.null(param) && (!is.character(param) || length(param) > 1 ||
                            !param %in% colnames(TD))) {
      stop(paste(deparse(param), "has to be NULL or a column in data.\n"))
    }
  }
  if (!is.null(engine) && (!is.character(engine) || length(engine) > 1 ||
                           !engine %in% c("asreml", "lme4", "SpATS"))) {
    stop("engine should be asreml or lme4")
  }
  ## Extract design from TD if needed.
  if (is.null(subDesign)) {
    subDesign <- attr(TD, "design")
  }
  # any check ID
  if (missing(checkId)) {
    checks <- FALSE
    if (is.null(rep)) {
      iNames <- c(trait, "genotype", subBlock)
    } else {
      iNames <- c(trait, "genotype", rep, subBlock)
    }
  } else {
    checks <- checkId %in% colnames(TD)
    if (is.null(rep)) {
      iNames <- c(trait, "genotype", subBlock, checkId)
    } else {
      iNames <- c(trait, "genotype", rep, subBlock, checkId)
    }
  }
  ## any covariate
  covT <- FALSE
  if (!is.null(covariate)) {
    if (is.character(covariate)) {
      covT <- TRUE
      iNames <- c(iNames, covariate)
    }
  }
  ## check validility of column names of TD
  TDNames <- names(TD)
  if (all(iNames %in% TDNames)) {
    vNameTest <- isValidVariableName(iNames)
    if (!all(vNameTest)) {
      warning(paste(iNames[!vNameTest], collapse = ","), " not syntactically valid name(s)")
    }
  } else {
    stop(paste(iNames[!(iNames %in% TDNames)], collapse = ","), " not found in the names of TD")
  }
  if (engine == "SpATS") {
    nSeg <- c(floor(nlevels(TD$col) / 2), floor(nlevels(TD$row) / 2))
    if (subDesign == "res.ibd") {
      mr <- SpATS::SpATS(response = trait, genotype = "genotype",
                         genotype.as.random = TRUE,
                         spatial = ~ SAP(c, r, nseg = nSeg, degree = 3, pord = 2),
                         fixed = as.formula(paste("~", rep,
                                                  if (checks) paste("+", checkId),
                                                  if (covT) paste(c("", covariate),
                                                                  collapse = "+"))),
                         random = ~ rep:subBlock,
                         data = TD,
                         ...)
    } else if (subDesign == "ibd") {

    }
  } else if (engine == "asreml") {
    tmp <- tempfile()
    sink(file = tmp)
    ## Run mixed and fixed models using asreml
    if (subDesign == "res.ibd") {
      fixedFormR <- as.formula(paste(trait, "~", rep,
                                     if (checks) paste("+", checkId),
                                     if (covT) paste(c("", covariate), collapse = "+")))
      mr <- asreml::asreml(fixed = fixedFormR,
                           random = as.formula(paste0("~ genotype +", rep, ":",
                                                      subBlock)),
                           rcov = ~ units, aom = TRUE, data = TD, ...)
      ## Constrain variance of the variance components to be fixed as the values in mr.
      GParamTmp <- mr$G.param
      GParamTmp[[paste0("`", rep, ":", subBlock, "`")]][[rep]]$con <- "F"
      fixedFormF <- as.formula(paste(deparse(fixedFormR), "+ genotype"))
      mf <- asreml::asreml(fixed = fixedFormF,
                           random = as.formula(paste0("~", rep, ":", subBlock)),
                           rcov = ~ units, G.param = GParamTmp, aom = TRUE, data = TD, ...)
    } else if (subDesign == "ibd") {
      fixedFormR <- as.formula(paste(trait, "~",
                                     if (checks) checkId else "1",
                                     if (covT) paste(c("", covariate), collapse = "+")))
      mr <- asreml::asreml(fixed = fixedFormR,
                           random = as.formula(paste0("~ genotype:", subBlock)),
                           rcov = ~ units, aom = TRUE, data = TD, ...)
      ## Constrain variance of the variance components to be fixed as the values in mr.
      GParamTmp <- mr$G.param
      GParamTmp[[subBlock]][[subBlock]]$con <- "F"
      fixedFormF <- as.formula(paste(deparse(fixedFormR), "+ genotype"))
      mf <- asreml::asreml(fixed = fixedFormF,
                           random = as.formula(paste("~", subBlock)), rcov = ~ units,
                           G.param = GParamTmp, aom = TRUE, data = TD, ...)
    }
    ## Run predict.
    mr$call$fixed <- eval(mr$call$fixed)
    mr$call$random <- eval(mr$call$random)
    mr$call$rcov <- eval(mr$call$rcov)
    mf$call$fixed <- eval(mf$call$fixed)
    mf$call$random <- eval(mf$call$random)
    mf$call$rcov <- eval(mf$call$rcov)
    if (checks) {
      mf <- predict(mf, classify = "genotype", vcov = TRUE, data = TD,
                    associate = as.formula(paste("~", checkId, ":genotype")))
    } else {
      mf <- predict(mf, classify = "genotype", vcov = TRUE, data = TD)
    }
    mr <- predict(mr, classify = "genotype", data = TD)
    mf$call$data <- substitute(TD)
    mr$call$data <- substitute(TD)
    sink()
    unlink(tmp)
  } else if (engine == "lme4") {
    ## Run mixed and fixed models using lme4
    if (subDesign == "res.ibd") {
      frm <- as.formula(paste(trait, "~", rep,
                              if (checks) paste("+", checkId),
                              if (covT) paste(c("", covariate), collapse = "+"),
                              "+ (1 | genotype) + (1 | ", rep, ":", subBlock, ")"))
      mr <- lme4::lmer(frm, data = TD, ...)
      ffm <- as.formula(paste(trait, "~", rep,
                              if (checks) paste("+", checkId),
                              if (covT) paste(c("", covariate), collapse = "+"),
                              "+ genotype + (1 | ", rep, ":", subBlock, ")"))
      mf <- lme4::lmer(ffm, data = TD, ...)
    } else if (subDesign == "ibd") {
      frm <- as.formula(paste(trait, "~",
                              if (checks) paste("+", checkId),
                              if (covT) paste(c("", covariate), collapse = "+"),
                              "+ (1 | genotype) + (1 |", subBlock, ")"))
      mr <- lme4::lmer(frm, data = TD, ...)
      ffm <- as.formula(paste(trait, "~",
                              if (checks) paste("+", checkId),
                              if (covT) paste(c("", covariate), collapse = "+"),
                              "+ genotype + (1 |", subBlock, ")"))
      mf <- lme4::lmer(ffm, data = TD, ...)
    }
  }
  model = createSSA(mMix = mr, mFix = mf, data = TD, trait = trait,
                    genotype = "genotype",
                    rep = ifelse(subDesign == "res.ibd", rep, NULL),
                    design = subDesign, engine = engine)
  return(model)
}
