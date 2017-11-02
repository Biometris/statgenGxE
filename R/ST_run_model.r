#' Run the single trial mixed model analysis
#'
#' Perform REML analysis given a specific experiment design.
#' This is a wrapper funtion of \code{\link[RAP]{ST.mod.alpha}}, \code{\link[RAP]{ST.mod.rcbd}} and \code{\link[RAP]{ST.mod.rowcol}}.
#' @param data A data frame object.
#' @param design A string specifying the experimental design.
#' @param trait A string specifying the selected trait.
#' @param covariate A string specifying (a) covariate name(s).
#' @param genotype A string specifying the column name of the genotypes.
#' @param rep A string specifying the column name of the replicates.
#' @param subblock A string specifying the name of the sub-blocks.
#' @param row A string specifying the column name of the rows.
#' @param col A string specifying the column name of the columns.
#' @param rowcoordinates A string specifying row coordinates for fitting spatial models. Default, \code{NA}.
#' @param colcoordinates A string specifying col coordinates for fitting spatial models. Default, \code{NA}.
#' @param checkid (optional) a string specifying the column name of the check ID(s).
#' @param tryspatial Whether to try spatial models ("always", "ifregular"); default no spatial models, i.e., \code{NA}.
#' @param ... Further arguments to be passed to either \code{asreml} or \code{lme4}.
#' @return S3 \code{model} object; a list consists of:
#'   \item{mmix}{an asreml or an lme4 class object.}
#'   \item{mfix}{an asreml or an lme4 class object.}
#'   \item{Data}{a list of 4 fields, data.frame: data containing traits, environment and design factors,
#'               character: time created, character: meta information, data.frame: checks.}
#' The object also has the following attributes:
#'   \item{attr(*, "Trait")}{character: a name specifying trait;}
#'   \item{attr(*, "Design")}{character: a type of experimental design;}
#'   \item{attr(*, "Engine")}{character: a mixed modelling engine;}
#'   \item{attr(*, "genotype")}{character: a name specifying genotype;}
#'   \item{attr(*, "rep")}{character: a name specifying rep;}
#' @examples
#' mydat <- ST.read.csv(file.path(path.package("RAP"),"SB_yield.csv"),
#'                      factor.names=c("Env","Genotype","Rep","Subblock","Row","Column"),
#'                      trait.names="yield", env ="Env",rowSelect="HEAT05",
#'                      colSelect=c("Env","Genotype","Rep","Row","Column","yield"))
#' mymodel <- ST.run.model(mydat, design="res.rowcol", trait="yield",
#'                         genotype="Genotype", rep="Rep", row="Row",
#'                         col="Column", tryspatial=NA)
#' names(mymodel)
#'
#' @export
ST.run.model = function(data, design, trait, covariate, genotype, rep, subblock, row, col,
               rowcoordinates=NA, colcoordinates=NA, checkid, tryspatial=NA, ...){

  #choose a package for mixed modelling
  engine <- asremlORlme4()

  if (design == "res.ibd")
    model = ST.mod.alpha(Y=data, trait=trait, genotype=genotype, rep=rep, subblock=subblock, sub.design=design, checkid=checkid, engine=engine, ...)
  if (design == "ibd")
    model = ST.mod.alpha(Y=data, trait=trait, genotype=genotype, subblock=subblock, sub.design=design, checkid=checkid, engine=engine, ...)
  if (design == "rcbd")
    model = ST.mod.rcbd(Y=data, trait=trait, genotype=genotype, rep=rep, checkid=checkid, engine=engine, ...)
  if (design == "res.rowcol")
    model = ST.mod.rowcol(Y=data, trait=trait, genotype=genotype, rep=rep, row=row, col=col, rowcoordinates=rowcoordinates, colcoordinates=colcoordinates, checkid=checkid, sub.design=design, tryspatial=tryspatial, engine=engine, ...)
  if (design == "rowcol")
    model = ST.mod.rowcol(Y=data, trait=trait, genotype=genotype, row=row, col=col, rowcoordinates=rowcoordinates, colcoordinates=colcoordinates, checkid=checkid, sub.design=design, tryspatial=tryspatial, engine=engine, ...)

  model
}
