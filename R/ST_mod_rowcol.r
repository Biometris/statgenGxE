#' Single trials (ST) models for row-column design (rowcol) or resolvable
#' row-column design (res.rowcol)
#'
#' Phenotypic traits will be analysed by fitting mixed models to obtain estimates of
#' genotypic means and several genetic parameters. Two mixed models are fitted; model \code{a}
#' fits genotypes as a random factor to obtain genetic variance components; model \code{b}
#' fits genotypes as fixed to obtain estimates of genotypic means.
#'
#' @inheritParams ST.run.model
#'
#' @param subDesign A string specifying whether to analyse a row-column (rowcol) or
#' resolvable row-column design (res.rowcol).
#' @param engine A string specifying the name of the mixed modelling engine to use.
#' @param trySpatial Whether to try spatial models (always, ifregular); default no spatial models.
#' @param ... Further arguments to be passed to either \code{asreml} or \code{lme4}.
#'
#' @return an object of class \code{\link{SSA}}.
#'
#' @seealso \code{\link{createSSA}}, \code{\link{summary.SSA}}
#'
#' @examples
#' myDat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames = c("Env", "Genotype", "Rep", "Subblock", "Row", "Column"),
#'                      traitNames = "yield", env = "Env", rowSelect = "HEAT05",
#'                      colSelect = c("Env", "Genotype", "Rep", "Row", "Column", "yield"))
#' myTD <- createTD(data = myDat, genotype = "Genotype", env = "Env")
#' myModel <- ST.mod.rowcol(TD = myTD, subDesign = "res.rowcol", trait = "yield",
#'                          rep = "Rep", row = "Row", col = "Column",
#'                          engine = "lme4") #engine = "asreml"
#' summary(myModel)
#'
#' @export
ST.mod.rowcol <- function(TD,
                          trait,
                          covariate = NULL,
                          rep = NULL,
                          row = NULL,
                          col = NULL,
                          rowCoordinates = NULL,
                          colCoordinates = NULL,
                          checkId = NULL,
                          subDesign = NULL,
                          trySpatial = NULL,
                          engine,
                          ...) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  subDesigns <- c("rowcol", "res.rowcol")
  if ((is.null(subDesign) && (is.null(attr(TD, "subDesign")) ||
                              !attr(TD, "subDesign") %in% subDesigns)) ||
      (!is.null(subDesign) && (!is.character(subDesign) || length(subDesign) > 1 ||
                               !subDesign %in% subDesigns))) {
    stop("subDesign should either be an attribute of TD or one of rowcol or res.rowcol.\n")
  }
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TD)) {
    stop("trait has to be a column in TD.\n")
  }
  if (!is.null(covariate) && (!is.character(covariate) ||
                              !(all(covariate %in% colnames(TD))))) {
    stop("covariate have to be a columns in TD.\n")
  }
  for (param in c(rep, row, col, rowCoordinates, colCoordinates, checkId)) {
    if (!is.null(param) && (!is.character(param) || length(param) > 1 ||
                            !param %in% colnames(TD))) {
      stop(paste(deparse(param), "has to be NULL or a column in data.\n"))
    }
  }
  if (!is.null(trySpatial) && (!is.character(trySpatial) || length(trySpatial) > 1 ||
                               !trySpatial %in% c("always", "ifregular"))) {
    stop("trySpatial should be NULL, always or ifregular.\n")
  }
  if (!is.null(engine) && (!is.character(engine) || length(engine) > 1 ||
                           !engine %in% c("asreml", "lme4"))) {
    stop("engine should be asreml or lme4.\n")
  }
  ## Extract design from TD if needed.
  if (is.null(subDesign)) {
    subDesign <- attr(TD, "design")
  }
  ## any check ID
  if (is.null(checkId)) {
    checks <- FALSE
    if (is.null(rep)) {
      iNames <- c(trait, "genotype", row, col)
    } else {
      iNames <- c(trait, "genotype", rep, row, col)
    }
  } else {
    checks <- checkId %in% colnames(TD)
    if (is.null(rep)) {
      iNames <- c(trait, "genotype", row, col, checkId)
    } else {
      iNames <- c(trait, "genotype", rep, row, col, checkId)
    }
  }
  if (!is.null(rowCoordinates)) {
    iNames <- c(iNames, rowCoordinates)
  }
  if (!is.null(colCoordinates)) {
    iNames <- c(iNames, colCoordinates)
  }
  # any covariate
  covT <- FALSE
  if (!is.null(covariate)) {
    if (is.character(covariate)) {
      covT <- TRUE
      iNames <- c(iNames, covariate)
    }
  } else {
    covariate <- NULL
  }
  #check validility of column names of TD
  TDNames <- names(TD)
  if (all(iNames %in% TDNames)) {
    vNameTest <- isValidVariableName(iNames)
    if (!all(vNameTest)) {
      warning(paste(iNames[!vNameTest], collapse = ",")," not syntactically valid name(s).\n")
    }
  } else {
    stop(paste(iNames[!(iNames %in% TDNames)], collapse = ","), " not found in the names of TD.\n")
  }
  if (!checks) {
    checkId <- NULL
  }
  if (engine == "SpATS") {

  } else if (engine == "asreml") {
    if (subDesign == "res.rowcol") {
      model <- ST.Varowcol(TD = TD, trait = trait, covariate = covariate,
                           rep = rep, row = row, col = col, tryRep = TRUE,
                           checkId = checkId, rowCoordinates = rowCoordinates,
                           colCoordinates = colCoordinates, trySpatial = trySpatial,
                           criterion = "BIC", ...)
    } else if (subDesign == "rowcol") {
      model <- ST.Varowcol(TD = TD, trait = trait, covariate = covariate,
                           rep = rep, row = row, col = col, tryRep = FALSE,
                           checkId = checkId, rowCoordinates = rowCoordinates,
                           colCoordinates = colCoordinates, trySpatial = trySpatial,
                           criterion = "BIC", ...)
    }
    model$mFix$call$data <- substitute(TD)
    model$mMix$call$data <- substitute(TD)
    model$design = subDesign
  } else if (engine == "lme4") {
    if (subDesign  == "res.rowcol") {
      frm <- as.formula(paste(trait, "~", rep,
                              if (checks) paste("+", checkId),
                              if (covT) paste(c("", covariate), collapse = "+"),
                              "+ (1| genotype) + (1|", rep, ":", row, ") + (1|",
                              rep, ":",col,")"))
      mr <- lme4::lmer(frm, data = TD, ...)
      ffm <- as.formula(paste(trait, "~", rep,
                              if (checks) paste("+", checkId),
                              if (covT) paste(c("", covariate), collapse = "+"),
                              "+ genotype + (1|", rep, ":", row,
                              ") + (1|", rep, ":", col,")"))
      mf <- lme4::lmer(ffm, data = TD, ...)
    } else if (subDesign == "rowcol") {
      frm <- as.formula(paste(trait, "~",
                              if (checks) checkId else "1",
                              if (covT) paste(c("", covariate), collapse = "+"),
                              "+ (1| genotype) + (1|", row, ") + (1|", col, ")"))
      mr <- lme4::lmer(frm, data = TD, ...)
      ffm <- as.formula(paste(trait, "~",
                              if (checks) checkId else "1",
                              if (covT) paste(c("", covariate), collapse = "+"),
                              "+ genotype + (1|", row, ") + (1|", col, ")"))
      mf <- lme4::lmer(ffm, data = TD, ...)
    }
    model = createSSA(mMix = mr, mFix = mf, data = TD, trait = trait,
                      genotype = "genotype",
                      rep = ifelse(subDesign == "res.rowcol", rep, NULL),
                      design = subDesign, engine = engine)
  }
  return(model)
}
