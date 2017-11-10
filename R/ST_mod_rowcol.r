#' Single trials (ST) models for row-column design (rowcol) or resolvable row-column design (res.rowcol)
#'
#' Phenotypic traits will be analysed by fitting mixed models to obtain estimates of
#' genotypic means and several genetic parameters. Two mixed models are fitted; model \code{a}
#' fits genotypes as a random factor to obtain genetic variance components; model \code{b}
#' fits genotypes as fixed to obtain estimates of genotypic means.
#' @param Y A data frame object.
#' @param engine A string specifying the name of the mixed modelling engine to use.
#' @param trait A string specifying the trait name.
#' @param covariate A string specifying (a) covariate name(s).
#' @param genotype A string specifying the column name of the genotypes.
#' @param rep A string specifying the column name of the replicates.
#' @param row A string specifying the column name of the rows.
#' @param col A string specifying the column name of the columns.
#' @param rowCoordinates A string specifying row coordinates for fitting spatial models.
#' Default, \code{NA}.
#' @param colCoordinates A string specifying col coordinates for fitting spatial models.
#' Default, \code{NA}.
#' @param checkId (Optional) a string specifying the column name of the check ID(s).
#' @param subDesign A string specifying whether to analyse a row-column (rowcol) or
#' resolvable row-column design (res.rowcol).
#' @param trySpatial Whether to try spatial models (always, ifregular); default no spatial models.
#' @param ... Further arguments to be passed to either \code{asreml} or \code{lme4}.
#'
#' @return an object of class \code{\link{SSA}}.
#'
#' @examples
#' mydat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames=c("Env","Genotype","Rep","Subblock","Row","Column"),
#'                      traitNames="yield", env ="Env",rowSelect="HEAT05",
#'                      colSelect=c("Env","Genotype","Rep","Row","Column","yield"))
#' mymodel <- ST.mod.rowcol(Y=mydat, subDesign="res.rowcol", trait="yield",
#'                          genotype="Genotype", rep="Rep", row="Row", col="Column",
#'                          engine="lme4") #engine="asreml"
#' summary(mymodel)
#'
#' @export
ST.mod.rowcol <- function(Y,
                          trait,
                          covariate,
                          genotype,
                          rep,
                          row,
                          col,
                          rowCoordinates = NA,
                          colCoordinates = NA,
                          checkId,
                          subDesign,
                          trySpatial = NA,
                          engine,
                          ...) {
  # any check ID
  if (missing(checkId)) {
    checks <- FALSE
    if (missing(rep)) {
      iNames <-c(trait, genotype, row, col)
    } else {
      iNames <-c(trait, genotype, rep, row, col)
    }
  } else {
    checks <- checkId %in% colnames(Y)
    if (missing(rep)) {
      iNames <-c(trait, genotype, row, col, checkId)
    } else {
      iNames <-c(trait, genotype, rep, row, col, checkId)
    }
  }
  if (!is.na(rowCoordinates)) {
    iNames <- c(iNames, rowCoordinates)
  }
  if (!is.na(colCoordinates)) {
    iNames <- c(iNames, colCoordinates)
  }
  # any covariate
  covT <- FALSE
  if (!missing(covariate)) {
    if (is.character(covariate)) {
      covT <- TRUE
      iNames <-c(iNames, covariate)
    }
  } else {
    covariate <- NULL
  }
  #check validility of column names of Y
  YNames <- names(Y)
  if (all(iNames %in% YNames)) {
    vNameTest <- isValidVariableName(iNames)
    if (!all(vNameTest)) {
      warning(paste(iNames[!vNameTest], collapse = ",")," not syntactically valid name(s).\n")
    }
  } else {
    stop(paste(iNames[!(iNames%in% YNames)], collapse = ","), " not found in the names of Y.\n")
  }
  if (!checks) {
    checkId <- NA
  }
  if (engine == "asreml"){
    if (subDesign == "res.rowcol") {
      model <- ST.Varowcol(Y = Y, trait = trait, covariate = covariate, genotype = genotype,
                         rep = rep, row = row, col = col, tryRep = TRUE,
                         checkId = checkId, rowCoordinates = rowCoordinates,
                         colCoordinates = colCoordinates, trySpatial = trySpatial,
                         criterion = "BIC", ...)
    } else if (subDesign == "rowcol") {
      model <- ST.Varowcol(Y = Y, trait = trait, covariate = covariate, genotype = genotype,
                         rep = rep, row = row, col = col, tryRep = FALSE,
                         checkId = checkId, rowCoordinates = rowCoordinates,
                         colCoordinates = colCoordinates, trySpatial = trySpatial,
                         criterion = "BIC", ...)
    }
    model$mFix$call$data <- substitute(Y)
    model$mMix$call$data <- substitute(Y)
    model$design = subDesign
  } else if (engine == "lme4") {
    if (subDesign  == "res.rowcol") {
      if (checks) {
        frm <- as.formula(paste(trait, "~", rep, "+", checkId,
                                if (covT) paste(c("", covariate), collapse = "+"),
                                "+ (1|", genotype, ") + (1|", rep, ":", row, ") + (1|",
                                rep, ":",col,")"))
      } else {
        frm <- as.formula(paste(trait, "~", rep,
                                if (covT) paste(c("", covariate), collapse = "+"),
                                "+ (1|", genotype, ") + (1|", rep, ":", row,
                                ") + (1|", rep, ":", col, ")"))
      }
      mr <- lme4::lmer(frm, data = Y, ...)
      if (checks) {
        ffm <- as.formula(paste(trait, "~", rep, "+", checkId,
                                if (covT) paste(c("", covariate), collapse = "+"),
                                "+", genotype, "+ (1|", rep, ":", row,
                                ") + (1|", rep, ":", col,")"))
      } else {
        ffm <- as.formula(paste(trait, "~", rep,
                                if (covT) paste(c("", covariate), collapse = "+"),
                                "+", genotype,  "+ (1|", rep, ":", row, ") + (1|",
                                rep, ":", col, ")"))
      }
      mf <- lme4::lmer(ffm, data = Y, ...)
    } else if(subDesign == "rowcol") {
      if (checks) {
        frm <- as.formula(paste(trait, "~", checkId,
                                if(covT) paste(c("", covariate), collapse = "+"),
                                "+ (1|", genotype, ") + (1|", row, ") + (1|", col, ")"))
      } else {
        frm <- as.formula(paste(trait, "~1",
                                if (covT) paste(c("", covariate), collapse = "+"),
                                "+ (1|", genotype, ") + (1|", row, ") + (1|", col, ")"))
      }
      mr <- lme4::lmer(frm, data = Y, ...)
      if (checks) {
        ffm <- as.formula(paste(trait, "~", checkId,
                                if (covT) paste(c("", covariate), collapse = "+"),
                                "+", genotype, "+ (1|", row, ") + (1|", col, ")"))
      } else {
        ffm <- as.formula(paste(trait, "~1",
                                if (covT) paste(c("", covariate), collapse = "+"),
                                "+", genotype, "+ (1|", row, ") + (1|", col, ")"))
      }
      mf <- lme4::lmer(ffm, data = Y, ...)
    }
    model = createSSA(mMix = mr, mFix = mf, data = Y, trait = trait,
                      genotype = genotype,
                      rep = ifelse(subDesign == "res.rowcol", rep, NULL),
                      design = subDesign, engine = engine)
  } else {
    stop("Please use either asreml or lme4 for engine")
  }
  return(model)
}
