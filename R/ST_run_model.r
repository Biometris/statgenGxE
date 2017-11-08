#' Run the single trial mixed model analysis
#'
#' Perform REML analysis given a specific experiment design.
#' This is a wrapper funtion of \code{\link[RAP]{ST.mod.alpha}}, \code{\link[RAP]{ST.mod.rcbd}}
#' and \code{\link[RAP]{ST.mod.rowcol}}.
#' @param data A data frame object.
#' @param design A string specifying the experimental design.
#' @param trait A string specifying the selected trait.
#' @param covariate A string specifying (a) covariate name(s).
#' @param genotype A string specifying the column name of the genotypes.
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
#' @return S3 \code{model} object; a list consists of:
#'   \item{mMix}{an asreml or an lme4 class object.}
#'   \item{mFix}{an asreml or an lme4 class object.}
#'   \item{data}{a list of 4 fields, data.frame: data containing traits, environment
#'   and design factors, character: time created, character: meta information, data.frame: checks.}
#' The object also has the following attributes:
#'   \item{attr(*, "Trait")}{character: a name specifying trait;}
#'   \item{attr(*, "Design")}{character: a type of experimental design;}
#'   \item{attr(*, "Engine")}{character: a mixed modelling engine;}
#'   \item{attr(*, "genotype")}{character: a name specifying genotype;}
#'   \item{attr(*, "rep")}{character: a name specifying rep;}
#'
#' @examples
#' mydat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames=c("Env","Genotype","Rep","Subblock","Row","Column"),
#'                      traitNames="yield", env ="Env",rowSelect="HEAT05",
#'                      colSelect=c("Env","Genotype","Rep","Row","Column","yield"))
#' mymodel <- ST.run.model(mydat, design="res.rowcol", trait="yield",
#'                         genotype="Genotype", rep="Rep", row="Row",
#'                         col="Column", trySpatial=NA)
#' names(mymodel)
#'
#' @export
ST.run.model = function(data,
                        design,
                        trait,
                        covariate,
                        genotype,
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
    model <- ST.mod.alpha(Y = data, trait = trait, genotype = genotype, rep = rep,
                          subBlock = subBlock, subDesign = design, checkId = checkId,
                          engine = engine, ...)
  } else if (design == "ibd") {
    model <- ST.mod.alpha(Y = data, trait = trait, genotype = genotype, subBlock = subBlock,
                          subDesign = design, checkId = checkId, engine = engine, ...)
  } else if (design == "rcbd") {
    model <- ST.mod.rcbd(Y = data, trait = trait, genotype = genotype, rep = rep,
                         checkId = checkId, engine = engine, ...)
  } else if (design == "res.rowcol") {
    model <- ST.mod.rowcol(Y = data, trait = trait, genotype = genotype, rep = rep,
                           row = row, col = col, rowCoordinates = rowCoordinates,
                           colCoordinates = colCoordinates, checkId = checkId,
                           subDesign = design, trySpatial = trySpatial, engine = engine, ...)
  } else if (design == "rowcol") {
    model <- ST.mod.rowcol(Y = data, trait = trait, genotype = genotype, row = row,
                           col = col, rowCoordinates = rowCoordinates,
                           colCoordinates = colCoordinates, checkId = checkId,
                           subDesign = design, trySpatial = trySpatial, engine = engine, ...)
  }
  return(model)
}
