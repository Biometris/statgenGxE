#' Single trial (ST) modelling for an incomplete-block design (ibd) or resolvable incomplete-block design (res.ibd)
#'
#' Phenotypic traits are analysed by fitting mixed models to obtain estimates of
#' genotypic means and several genetic parameters. Two mixed models are fitted; the first
#' fits genotypes as a random factor to obtain genetic variance components; and the second
#' fits genotypes as fixed to obtain estimates of genotypic means.
#' @param Y A data frame object.
#' @param trait A string specifying the trait name.
#' @param covariate A string specifying (a) covariate name(s).
#' @param genotype A string specifying the column name of the genotypes.
#' @param rep A string specifying the column name of the replicates.
#' @param subblock A string specifying the column name of the sub-blocks.
#' @param sub.design A string specifying whether to analyse an incomplete-block design (ibd) or resolvable incomplete-block design (res.ibd).
#' @param checkid (optional) a string specifying the column name of the check ID(s).
#' @param engine A string specifying the name of the mixed modelling engine to use.
#' @param ... Further arguments to be passed to either \code{asreml} or \code{lme4}.
#' @return a list with fields \code{mmix}, \code{mfix} and \code{Data}.
#' @examples
#' mydat <- ST.read.csv(file.path(path.package("RAP"),"SB_yield.csv"),
#'                      factor.names=c("Env","Genotype","Rep","Subblock"),
#'                      trait.names="yield", env ="Env",rowSelect="HEAT05",
#'                      colSelect=c("Env","Genotype","Rep","Subblock","yield"))
#' mymodel <- ST.mod.alpha(Y=mydat, sub.design="res.ibd", trait="yield",
#'                         genotype="Genotype", rep="Rep", subblock="Subblock",
#'                         engine="lme4") #engine="asreml"
#' summary(mymodel)
#' 
#' @export
ST.mod.alpha = function(Y, trait, covariate, genotype, rep, subblock, sub.design, checkid, engine, ...) {

  # any check ID
  if (missing(checkid)){
    checks <- FALSE
    if (missing(rep))
      inames <-c(trait, genotype, subblock)
    else
      inames <-c(trait, genotype, rep, subblock)
  }else{
    checks <- checkid %in% colnames(Y)
    if (missing(rep))
      inames <-c(trait, genotype, subblock, checkid)
    else
      inames <-c(trait, genotype, rep, subblock, checkid)
  }
  # any covariate
  covT <- FALSE
  if (!missing(covariate)) {
    if (is.character(covariate)){
      covT <- TRUE
      inames <-c(inames, covariate)
    }
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

  if (engine=="asreml"){
    # Run mixed and fixed models using asreml
    require(asreml, quietly=T)
    if (sub.design=="res.ibd"){
        if (checks){
            mr = asreml::asreml(fixed=as.formula(paste(trait,"~", rep, "+", checkid, if(covT) paste(c("",covariate), collapse="+"))),
                 random=as.formula(paste("~",genotype, "+",rep,":",subblock, sep="")),
                 rcov=~units, aom=T, data=Y, ...)
        }else{
            mr = asreml::asreml(fixed=as.formula(paste(trait,"~",rep, if(covT) paste(c("",covariate), collapse="+"))),
            random=as.formula(paste("~",genotype, "+",rep,":",subblock, sep="")),
            rcov=~units, aom=T, data=Y, ...)
        }
        # constrain variance of the variance components to be fixed as the values in mr
        G.param.tmp=mr$G.param
        G.param.tmp[[paste("`",rep,":",subblock,"`",sep="")]][[rep]]$con="F"
        if (checks){
            mf = asreml::asreml(fixed=as.formula(paste(trait,"~",rep,"+",checkid, if(covT) paste(c("",covariate), collapse="+"), "+", genotype)),
                 random=as.formula(paste("~",rep,":",subblock, sep="")), rcov=~units, G.param=G.param.tmp, aom=T, data=Y, ...)
        }else{
            mf = asreml::asreml(fixed=as.formula(paste(trait,"~", rep, if(covT) paste(c("",covariate), collapse="+"), "+", genotype)),
                 random=as.formula(paste("~",rep,":",subblock, sep="")), rcov=~units, G.param=G.param.tmp, aom=T, data=Y, ...)
        }
    }
    if (sub.design=="ibd"){
        if (checks){
            mr = asreml::asreml(fixed=as.formula(paste(trait,"~",checkid, if(covT) paste(c("",covariate), collapse="+"))),
                 random=as.formula(paste("~",genotype,":",subblock, sep="")), rcov=~units, aom=T, data=Y, ...)
        }else{
            mr = asreml::asreml(fixed=as.formula(paste(trait,"~1", if(covT) paste(c("",covariate), collapse="+"))),
                 random=as.formula(paste("~",genotype,":",subblock, sep="")), rcov=~units, aom=T, data=Y, ...)
        }
        # constrain variance of the variance components to be fixed as the values in mr
        G.param.tmp=mr$G.param
        G.param.tmp[[subblock]][[subblock]]$con="F"
        if (checks){
            mf = asreml::asreml(fixed=as.formula(paste(trait,"~", checkid, if(covT) paste(c("",covariate), collapse="+"), "+", genotype)),
                 random=as.formula(paste("~",subblock)), rcov=~units, G.param=G.param.tmp, aom=T, data=Y, ...)
        }else{
            mf = asreml::asreml(fixed=as.formula(paste(trait,"~1", if(covT) paste(c("",covariate), collapse="+"), "+", genotype)),
                 random=as.formula(paste("~",subblock)), rcov=~units, G.param=G.param.tmp, aom=T, data=Y, ...)
        }
    }
    # run predict
    mr$call$fixed <- eval(mr$call$fixed)
    mr$call$random <- eval(mr$call$random)
    mr$call$rcov <- eval(mr$call$rcov)
    mf$call$fixed <- eval(mf$call$fixed)
    mf$call$random <- eval(mf$call$random)
    mf$call$rcov <- eval(mf$call$rcov)

    tmp=tempfile()
    sink(file=tmp)
    if (checks){
      mf = predict(mf, classify=genotype, vcov=T, data=Y, associate=~checkid:genotype)
    }else{
      mf = predict(mf, classify=genotype, vcov=T, data=Y)
    }
    mr = predict(mr, classify=genotype, data=Y)
    mf$call$data <- substitute(Y)
    mr$call$data <- substitute(Y)
    sink()
    unlink(tmp)
  }else{
    if (engine=="lme4"){
        # Run mixed and fixed models using lme4
        suppressPackageStartupMessages(require(lme4))
        if (sub.design=="res.ibd"){
            frm = as.formula(paste(trait,"~", rep, if(covT) paste(c("",covariate), collapse="+"),  "+ (1 | ",genotype, ")+ (1 | ", rep,":",subblock,")"))
            if (checks) frm = as.formula(paste(trait,"~", rep, "+", checkid, if(covT) paste(c("",covariate), collapse="+"), "+ (1 | ",genotype,")+ (1 | ",rep,":",subblock,")"))
            mr = lme4::lmer(frm, data = Y, ...)
            ffm = as.formula(paste(trait,"~", rep, if(covT) paste(c("",covariate), collapse="+"), "+", genotype, "+ (1 | ", rep,":",subblock,")"))
            if (checks) ffm = as.formula(paste(trait,"~", rep, "+", checkid, if(covT) paste(c("",covariate), collapse="+"),  "+",  genotype, "+ (1 | ", rep,":",subblock,")"))
            mf = lme4::lmer(ffm, data = Y, ...)
        }
        if (sub.design=="ibd"){
            frm = as.formula(paste(trait,"~1", if(covT) paste(c("",covariate), collapse="+"),"+ (1 |", genotype, ")+ (1 |", subblock,")"))
            if (checks) frm = as.formula(paste(trait,"~ ",checkid, if(covT) paste(c("",covariate), collapse="+"), "+ (1 | ",genotype, ")+ (1 | ",subblock,")"))
            mr = lme4::lmer(frm, data = Y, ...)
            ffm = as.formula(paste(trait,"~ 1", if(covT) paste(c("",covariate), collapse="+"), "+", genotype, "+ (1 | ", subblock, ")"))
            if (checks) ffm = as.formula(paste(trait,"~ ", checkid, if(covT) paste(c("",covariate), collapse="+"), "+", genotype,  "+ (1 | ",subblock, ")"))
            mf = lme4::lmer(ffm, data = Y, ...)
        }
    }else{
        stop("Please use either asreml or lme4 for engine")
    }
  }

  model = list(mmix = mr, mfix = mf, Data = Y)

  attr(model, "Trait") <- trait
  attr(model, "Design") <- sub.design
  attr(model, "Engine") <- engine
  attr(model, "genotype") <- genotype
  if (sub.design=="res.ibd")
    attr(model, "rep") <- rep

  class(model) <- "SSA"
  model
}
