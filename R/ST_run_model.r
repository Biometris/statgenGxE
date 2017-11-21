#' Run the single trial mixed model analysis
#'
#' Perform REML analysis given a specific experimental design.
#' This is a wrapper function of \code{\link{ST.mod.alpha}}, \code{\link{ST.mod.rcbd}}
#' and \code{\link{ST.mod.rowcol}}.
#'
#' @param TD An object of class TD.
#' @param design A string specifying the experimental design.
#' @param trait A string specifying the selected trait.
#' @param covariate A string specifying (a) covariate name(s).
#' @param rep A string specifying the column name of the replicates.
#' @param subBlock A string specifying the name of the sub-blocks.
#' @param row A string specifying the column name of the rows.
#' @param col A string specifying the column name of the columns.
#' @param rowCoordinates A string specifying row coordinates for fitting spatial models.
#' Default, \code{NULL}.
#' @param colCoordinates A string specifying col coordinates for fitting spatial models.
#' Default, \code{NULL}.
#' @param checkId (optional) a string specifying the column name of the check ID(s).
#' @param trySpatial Whether to try spatial models ("always", "ifregular"); default no spatial
#' models, i.e., \code{NULL}.
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
#' myModel <- ST.run.model(TD = myTD, design = "res.rowcol", trait = "yield",
#'                         rep = "Rep", row = "Row", col = "Column")
#' names(myModel)
#'
#' @export
ST.run.model = function(TD,
                        design = NULL,
                        trait,
                        covariate = NULL,
                        rep = NULL,
                        subBlock  = NULL,
                        row = NULL,
                        col = NULL,
                        rowCoordinates = NULL,
                        colCoordinates = NULL,
                        checkId = NULL,
                        trySpatial = NULL,
                        ...) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  designs <- c("ibd", "res.ibd", "rcbd", "rowcol", "res.rowcol")
  if ((is.null(design) && !attr(TD, "design") %in% designs) ||
      (!is.null(design) && (!is.character(design) || length(design) > 1 ||
                            !design %in% designs))) {
    stop("design should either be an attribute of TD or one of ibd,
         res.ibd, rcbd, rowcol or res.rowcol.\n")
  }
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TD)) {
    stop("trait has to be a column in TD.\n")
  }
  if (!is.null(covariate) && (!is.character(covariate) ||
      !(all(covariate %in% colnames(TD))))) {
    stop("covariate have to be a columns in TD.\n")
  }
  for (param in c(rep, subBlock, row, col, rowCoordinates,
                  colCoordinates, checkId)) {
    if (!is.null(param) && (!is.character(param) || length(param) > 1 ||
                            !param %in% colnames(TD))) {
      stop(paste(deparse(param), "has to be NULL or a column in data.\n"))
    }
  }
  if (!is.null(trySpatial) && (!is.character(trySpatial) || length(trySpatial) > 1 ||
      !trySpatial %in% c("always", "ifregular"))) {
    stop("trySpatial should be NULL, always or ifregular")
  }
  ## Extract design from TD if needed.
  if (is.null(design)) {
    design <- attr(TD, "design")
  }
  ## Automatically choose a package for mixed modelling.
  engine <- asremlORlme4()
  ## Modelling depending on design.
  if (design == "res.ibd") {
    model <- ST.mod.alpha(TD = TD, trait = trait, covariate = covariate,
                          rep = rep, subBlock = subBlock,
                          subDesign = design, checkId = checkId,
                          engine = engine, ...)
  } else if (design == "ibd") {
    model <- ST.mod.alpha(TD = TD, trait = trait, covariate = covariate,
                          subBlock = subBlock, subDesign = design,
                          checkId = checkId, engine = engine, ...)
  } else if (design == "rcbd") {
    model <- ST.mod.rcbd(TD = TD, trait = trait, covariate = covariate,
                         rep = rep, checkId = checkId,
                         engine = engine, ...)
  } else if (design == "res.rowcol") {
    model <- ST.mod.rowcol(TD = TD, trait = trait, covariate = covariate,
                           rep = rep, row = row, col = col,
                           rowCoordinates = rowCoordinates,
                           colCoordinates = colCoordinates,
                           checkId = checkId, subDesign = design,
                           trySpatial = trySpatial, engine = engine, ...)
  } else if (design == "rowcol") {
    model <- ST.mod.rowcol(TD = TD, trait = trait, covariate = covariate,
                           row = row, col = col, rowCoordinates = rowCoordinates,
                           colCoordinates = colCoordinates, checkId = checkId,
                           subDesign = design, trySpatial = trySpatial,
                           engine = engine, ...)
  }
  return(model)
}
