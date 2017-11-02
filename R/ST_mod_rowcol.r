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
#' @param rowcoordinates A string specifying row coordinates for fitting spatial models. Default, \code{NA}.
#' @param colcoordinates A string specifying col coordinates for fitting spatial models. Default, \code{NA}.
#' @param checkid (Optional) a string specifying the column name of the check ID(s).
#' @param sub.design A string specifying whether to analyse a row-column (rowcol) or resolvable row-column design (res.rowcol).
#' @param tryspatial Whether to try spatial models (always, ifregular); default no spatial models.
#' @param ... Further arguments to be passed to either \code{asreml} or \code{lme4}.
#' @return a list with fields \code{mmix}, \code{mfix} and \code{Data}.
#' @examples
#' mydat <- ST.read.csv(file.path(path.package("RAP"),"SB_yield.csv"),
#'                      factor.names=c("Env","Genotype","Rep","Subblock","Row","Column"),
#'                      trait.names="yield", env ="Env",rowSelect="HEAT05",
#'                      colSelect=c("Env","Genotype","Rep","Row","Column","yield"))
#' mymodel <- ST.mod.rowcol(Y=mydat, sub.design="res.rowcol", trait="yield",
#'                          genotype="Genotype", rep="Rep", row="Row", col="Column",
#'                          engine="lme4") #engine="asreml"
#' summary(mymodel)
#' 
#' @export
ST.mod.rowcol = function(Y, trait, covariate, genotype, rep, row, col, rowcoordinates=NA,colcoordinates=NA, checkid, sub.design, tryspatial=NA, engine,...) {

  # any check ID
  if (missing(checkid)){
    checks <- FALSE
    if (missing(rep))
      inames <-c(trait, genotype, row, col)
    else
      inames <-c(trait, genotype, rep, row, col)
  }else{
    checks <- checkid %in% colnames(Y)
    if (missing(rep))
      inames <-c(trait, genotype, row, col, checkid)
    else
      inames <-c(trait, genotype, rep, row, col, checkid)
  }
  if (!is.na(rowcoordinates))
    inames <- c(inames, rowcoordinates)
  if (!is.na(colcoordinates))
    inames <- c(inames, colcoordinates)
  # any covariate
  covT <- FALSE
  if (!missing(covariate)) {
    if (is.character(covariate)){
      covT <- TRUE
      inames <-c(inames, covariate)
    }
  }else{
      covariate <- NULL
  }
  #check validility of column names of Y
  Y_names <- names(Y)
  if (all(inames%in% Y_names)){
    vname_test <- is_valid_variable_name(inames)
    if(!all(vname_test))
      warning(paste(inames[!vname_test], collapse=",")," not syntactically valid name(s)")
  }else{
    stop(paste(inames[!(inames%in% Y_names)], collapse=","), " not found in the names of Y")
  }

  if (!checks) checkid = NA
  if (engine=="asreml"){
    if (sub.design =="res.rowcol"){
      res <- ST.Varowcol(Y=Y, trait=trait, covariate=covariate, genotype=genotype, rep=rep, row=row, col=col, tryrep = TRUE,
      checkid = checkid, rowcoordinates = rowcoordinates, colcoordinates = colcoordinates, tryspatial = tryspatial, criterion = "BIC",...)
    }
    if (sub.design == "rowcol"){
      res <- ST.Varowcol(Y=Y, trait=trait, covariate=covariate, genotype=genotype, rep=rep, row=row, col=col, tryrep = FALSE,
      checkid = checkid, rowcoordinates = rowcoordinates, colcoordinates = colcoordinates, tryspatial = tryspatial, criterion = "BIC",...)
    }
    res$mfix$call$data <- substitute(Y)
    res$mmix$call$data <- substitute(Y)
  }else{
    if (engine=="lme4"){
        # Run mixed and fixed models using lme4
        suppressPackageStartupMessages(require(lme4))
        if (sub.design =="res.rowcol"){
            frm = as.formula(paste(trait,"~", rep, if(covT) paste(c("",covariate), collapse="+"), "+ (1|", genotype, ") + (1|", rep, ":", row, ") + (1|", rep, ":", col, ")"))
            if (checks) frm = as.formula(paste(trait,"~", rep, "+", checkid, if(covT) paste(c("",covariate), collapse="+"),  "+ (1|", genotype, ") + (1|", rep, ":", row, ") + (1|", rep, ":",col,")"))
            mr = lme4::lmer(frm, data = Y,...)
            ffm = as.formula(paste(trait,"~", rep, if(covT) paste(c("",covariate), collapse="+"), "+", genotype,  "+ (1|", rep, ":", row, ") + (1|", rep, ":", col, ")"))
            if (checks) ffm = as.formula(paste(trait,"~", rep, "+", checkid, if(covT) paste(c("",covariate), collapse="+"), "+", genotype, "+ (1|", rep, ":", row, ") + (1|", rep, ":", col,")"))
            mf = lme4::lmer(ffm, data = Y,...)
        }
        if(sub.design =="rowcol"){
            frm = as.formula(paste(trait,"~1", if(covT) paste(c("",covariate), collapse="+"), "+ (1|", genotype, ") + (1|", row, ") + (1|", col, ")"))
            if (checks) frm = as.formula(paste(trait,"~", checkid, if(covT) paste(c("",covariate), collapse="+"), "+ (1|",genotype, ") + (1|", row, ") + (1|", col,")"))
            mr = lme4::lmer(frm, data = Y,...)
            ffm = as.formula(paste(trait,"~1", if(covT) paste(c("",covariate), collapse="+"), "+", genotype,  "+ (1|", row, ") + (1|", col, ")"))
            if (checks) ffm = as.formula(paste(trait,"~", checkid, if(covT) paste(c("",covariate), collapse="+"), "+", genotype, "+ (1|", row, ") + (1|", col, ")"))
            mf = lme4::lmer(ffm, data = Y,...)
        }
        res = list(mmix = mr, mfix = mf, Data = Y)

    }else{
        stop("Please use either asreml or lme4 for engine")
    }
  }

  attr(res, "Trait") <- trait
  attr(res, "Design") <- sub.design
  attr(res, "Engine") <- engine
  attr(res, "genotype") <- genotype
  if (sub.design=="res.rowcol")
    attr(res, "rep") <- rep

  class(res) <- "SSA"
  res
}
