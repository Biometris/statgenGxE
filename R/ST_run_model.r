#' Run the single trial mixed model analysis
#'
#' Perform REML analysis given a specific experiment design.
#' This is a wrapper funtion of \code{\link{ST.mod.alpha}}, \code{\link{ST.mod.rcbd}}
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
#' Default, \code{NA}.
#' @param colCoordinates A string specifying col coordinates for fitting spatial models.
#' Default, \code{NA}.
#' @param checkId (optional) a string specifying the column name of the check ID(s).
#' @param trySpatial Whether to try spatial models ("always", "ifregular"); default no spatial
#' models, i.e., \code{NA}.
#' @param ... Further arguments to be passed to either \code{asreml} or \code{lme4}.
#'
#' @return an object of class SSA.
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
#'                         rep = "Rep", row = "Row", col = "Column", tryspatial = NA)
#' names(myModel)
#'
#' @export
ST.run.model = function(TD,
                        design,
                        trait,
                        covariate,
                        rep,
                        subBlock,
                        row,
                        col,
                        rowCoordinates = NA,
                        colCoordinates = NA,
                        checkId,
                        trySpatial = NA,
                        ...) {
  #choose a package for mixed modelling
  engine <- asremlORlme4()
  if (design == "res.ibd") {
    model <- ST.mod.alpha(TD = TD, trait = trait, rep = rep, subBlock = subBlock,
                          subDesign = design, checkId = checkId,
                          engine = engine, ...)
  } else if (design == "ibd") {
    model <- ST.mod.alpha(TD = TD, trait = trait, subBlock = subBlock,
                          subDesign = design, checkId = checkId, engine = engine, ...)
  } else if (design == "rcbd") {
    model <- ST.mod.rcbd(TD = TD, trait = trait, rep = rep,
                         checkId = checkId, engine = engine, ...)
  } else if (design == "res.rowcol") {
    model <- ST.mod.rowcol(TD = TD, trait = trait, rep = rep,
                           row = row, col = col, rowCoordinates = rowCoordinates,
                           colCoordinates = colCoordinates, checkId = checkId,
                           subDesign = design, trySpatial = trySpatial, engine = engine, ...)
  } else if (design == "rowcol") {
    model <- ST.mod.rowcol(TD = TD, trait = trait, row = row,
                           col = col, rowCoordinates = rowCoordinates,
                           colCoordinates = colCoordinates, checkId = checkId,
                           subDesign = design, trySpatial = trySpatial, engine = engine, ...)
  }
  return(model)
}
