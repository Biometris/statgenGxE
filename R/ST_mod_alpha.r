#' Single trial (ST) modelling for an incomplete-block design (ibd) or resolvable
#' incomplete-block design (res.ibd)
#'
#' Phenotypic traits are analysed by fitting mixed models to obtain estimates of
#' genotypic means and several genetic parameters. Two mixed models are fitted; the first
#' fits genotypes as a random factor to obtain genetic variance components; and the second
#' fits genotypes as fixed to obtain estimates of genotypic means.
#'
#' @param Y A data frame object.
#' @param trait A string specifying the trait name.
#' @param covariate A string specifying (a) covariate name(s).
#' @param genotype A string specifying the column name of the genotypes.
#' @param rep A string specifying the column name of the replicates.
#' @param subBlock A string specifying the column name of the sub-blocks.
#' @param subDesign A string specifying whether to analyse an incomplete-block design (ibd) or
#' resolvable incomplete-block design (res.ibd).
#' @param checkId (optional) a string specifying the column name of the check ID(s).
#' @param engine A string specifying the name of the mixed modelling engine to use.
#' @param ... Further arguments to be passed to either \code{asreml} or \code{lme4}.
#'
#' @return an object of class \code{\link{SSA}}.
#' @examples
#' mydat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames=c("Env","Genotype","Rep","Subblock"),
#'                      traitNames="yield", env ="Env",rowSelect="HEAT05",
#'                      colSelect=c("Env","Genotype","Rep","Subblock","yield"))
#' mymodel <- ST.mod.alpha(Y=mydat, subDesign="res.ibd", trait="yield",
#'                         genotype="Genotype", rep="Rep", subBlock="Subblock",
#'                         engine="lme4") #engine="asreml"
#' summary(mymodel)
#'
#' @export
ST.mod.alpha = function(Y,
                        trait,
                        covariate,
                        genotype,
                        rep,
                        subBlock,
                        subDesign,
                        checkId,
                        engine,
                        ...) {
  # any check ID
  if (missing(checkId)) {
    checks <- FALSE
    if (missing(rep)) {
      iNames <- c(trait, genotype, subBlock)
    } else {
      iNames <- c(trait, genotype, rep, subBlock)
    }
  } else {
    checks <- checkId %in% colnames(Y)
    if (missing(rep)) {
      iNames <- c(trait, genotype, subBlock, checkId)
    } else {
      iNames <- c(trait, genotype, rep, subBlock, checkId)
    }
  }
  # any covariate
  covT <- FALSE
  if (!missing(covariate)) {
    if (is.character(covariate)){
      covT <- TRUE
      iNames <-c(iNames, covariate)
    }
  }
  #check validility of column names of Y
  YNames <- names(Y)
  if (all(iNames %in% YNames)) {
    vNameTest <- isValidVariableName(iNames)
    if(!all(vNameTest)) {
      warning(paste(iNames[!vNameTest], collapse = ","), " not syntactically valid name(s)")
    }
  } else {
    stop(paste(iNames[!(iNames%in% YNames)], collapse = ","), " not found in the names of Y")
  }
  if (engine == "asreml"){
    # Run mixed and fixed models using asreml
    if (subDesign == "res.ibd") {
      if (checks) {
        mr <- asreml::asreml(fixed = as.formula(paste(trait, "~", rep, "+", checkId,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"))),
                             random = as.formula(paste0("~", genotype, "+", rep, ":",
                                                        subBlock)),
                             rcov = ~units, aom = TRUE, data = Y, ...)
      } else {
        mr <- asreml::asreml(fixed = as.formula(paste(trait, "~", rep,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"))),
                             random = as.formula(paste0("~", genotype, "+", rep, ":",
                                                        subBlock)),
                             rcov = ~units, aom = TRUE, data = Y, ...)
      }
      # constrain variance of the variance components to be fixed as the values in mr
      GParamTmp <- mr$G.param
      GParamTmp[[paste0("`", rep, ":", subBlock, "`")]][[rep]]$con <- "F"
      if (checks) {
        mf <- asreml::asreml(fixed = as.formula(paste(trait, "~", rep, "+", checkId,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"),
                                                      "+", genotype)),
                             random = as.formula(paste0("~", rep, ":", subBlock)),
                             rcov = ~units, G.param = GParamTmp, aom = TRUE, data = Y, ...)
      } else {
        mf <- asreml::asreml(fixed = as.formula(paste(trait, "~", rep,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"),
                                                      "+", genotype)),
                             random = as.formula(paste0("~", rep, ":", subBlock)),
                             rcov =~ units, G.param = GParamTmp, aom = TRUE, data = Y, ...)
      }
    } else if (subDesign == "ibd") {
      if (checks) {
        mr <- asreml::asreml(fixed = as.formula(paste(trait, "~", checkId,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"))),
                             random = as.formula(paste0("~", genotype, ":", subBlock)),
                             rcov =~ units, aom = TRUE, data = Y, ...)
      } else {
        mr <- asreml::asreml(fixed = as.formula(paste(trait, "~1",
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"))),
                             random = as.formula(paste0("~", genotype, ":", subBlock)),
                             rcov = ~units, aom = TRUE, data = Y, ...)
      }
      # constrain variance of the variance components to be fixed as the values in mr
      GParamTmp <- mr$G.param
      GParamTmp[[subBlock]][[subBlock]]$con <- "F"
      if (checks) {
        mf <- asreml::asreml(fixed = as.formula(paste(trait, "~", checkId,
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"),
                                                      "+", genotype)),
                             random = as.formula(paste("~", subBlock)), rcov = ~units,
                             G.param = GParamTmp, aom = TRUE, data = Y, ...)
      } else {
        mf <- asreml::asreml(fixed = as.formula(paste(trait, "~1",
                                                      if (covT) paste(c("", covariate),
                                                                      collapse = "+"),
                                                      "+", genotype)),
                             random = as.formula(paste("~", subBlock)), rcov = ~units,
                             G.param = GParamTmp, aom = TRUE, data = Y, ...)
      }
    }
    # run predict
    mr$call$fixed <- eval(mr$call$fixed)
    mr$call$random <- eval(mr$call$random)
    mr$call$rcov <- eval(mr$call$rcov)
    mf$call$fixed <- eval(mf$call$fixed)
    mf$call$random <- eval(mf$call$random)
    mf$call$rcov <- eval(mf$call$rcov)
    tmp <- tempfile()
    sink(file=tmp)
    if (checks) {
      mf <- predict(mf, classify = genotype, vcov = TRUE, data = Y, associate = ~checkId:genotype)
    } else {
      mf <- predict(mf, classify = genotype, vcov = TRUE, data = Y)
    }
    mr <- predict(mr, classify = genotype, data = Y)
    mf$call$data <- substitute(Y)
    mr$call$data <- substitute(Y)
    sink()
    unlink(tmp)
  } else {
    if (engine == "lme4"){
      # Run mixed and fixed models using lme4
      if (subDesign == "res.ibd") {
        if (checks) {
          frm <- as.formula(paste(trait, "~", rep, "+", checkId,
                                  if (covT) paste(c("", covariate), collapse = "+"),
                                  "+ (1 | ", genotype, ")+ (1 | ", rep, ":", subBlock, ")"))
        } else {
          frm <- as.formula(paste(trait, "~", rep, if (covT) paste(c("",covariate), collapse = "+"),
                                  "+ (1 | ", genotype, ")+ (1 | ", rep, ":", subBlock, ")"))
        }
        mr <- lme4::lmer(frm, data = Y, ...)
        if (checks) {
          ffm <- as.formula(paste(trait, "~", rep, "+", checkId,
                                  if (covT) paste(c("", covariate), collapse = "+"),
                                  "+", genotype, "+ (1 | ", rep, ":", subBlock, ")"))
        } else {
          ffm <- as.formula(paste(trait, "~", rep, if(covT) paste(c("", covariate), collapse = "+"),
                                  "+", genotype, "+ (1 | ", rep, ":", subBlock, ")"))
        }
        mf <- lme4::lmer(ffm, data = Y, ...)
      } else if (subDesign == "ibd") {
        if (checks) {
          frm = as.formula(paste(trait, "~ ", checkId,
                                 if (covT) paste(c("", covariate), collapse = "+"),
                                 "+ (1 | ", genotype, ")+ (1 | ", subBlock, ")"))
        } else {
          frm <- as.formula(paste(trait, "~1", if (covT) paste(c("", covariate), collapse = "+"),
                                  "+ (1 |", genotype, ")+ (1 |", subBlock, ")"))
        }
        mr <- lme4::lmer(frm, data = Y, ...)
        if (checks) {
          ffm <- as.formula(paste(trait, "~ ", checkId,
                                  if (covT) paste(c("", covariate), collapse = "+"),
                                  "+", genotype,  "+ (1 | ", subBlock, ")"))
        } else {
          ffm <- as.formula(paste(trait, "~ 1", if (covT) paste(c("", covariate), collapse = "+"),
                                  "+", genotype, "+ (1 | ", subBlock, ")"))
        }
        mf <- lme4::lmer(ffm, data = Y, ...)
      }
    } else {
      stop("Please use either asreml or lme4 for engine")
    }
  }
  model = createSSA(mMix = mr, mFix = mf, data = Y)
  attr(model, "Trait") <- trait
  attr(model, "Design") <- subDesign
  attr(model, "Engine") <- engine
  attr(model, "genotype") <- genotype
  if (subDesign == "res.ibd") {
    attr(model, "rep") <- rep
  }
  return(model)
}
