#' Find the best random model for a row-column design (asreml only)
#'
#' This function fits a variety of random and spatial covariance models
#' and selects the best one using a goodness-of-fit criterion.
#'
#' @param Y A data frame object.
#' @param trait A character string specifying a 'trait' column of the data.
#' @param covariate A string specifying (a) covariate name(s); Default, \code{NULL}.
#' @param genotype A character string specifying a 'genotypes' column of the data.
#' @param rep A string specifying a 'replicates' column of the data.
#' @param row A string specifying a 'rows' column of the data.
#' @param col A string specigying a 'columns' column of the data.
#' @param tryrep A logical value indicating if 'replicates' are included in the model. Default, \code{TRUE}.
#' @param checkid A string specifying a 'checkid' column of the data. Default, \code{NA}.
#' @param rowcoordinates A string specifying row coordinates for fitting spatial models. Default, \code{NA}.
#' @param colcoordinates A string specifying col coordinates for fitting spatial models. Default, \code{NA}.
#' @param tryspatial Whether to try spatial models ("always", "ifregular"); default no spatial models, i.e., \code{NA}.
#' @param criterion A string specifies a goodness of fit criterion, i.e., "AIC" or "BIC".
#' @param ... Further arguments to be passed to \code{asreml}.
#' @return a list with fields \code{mmix}, \code{mfix} and \code{Data}.
#' @note This function can only be used if asreml is installed. If \code{tryspatial} is set to "always" or "ifregular",
#' the names for \code{rowcoordinates} and \code{colcoordinates} must be supplied; otherwise, no spatial model will be fitted.
#'
#' @export

ST.Varowcol <- function(Y, trait, covariate=NULL, genotype, rep, row, col, tryrep = TRUE,
checkid =NA, rowcoordinates = NA, colcoordinates = NA, tryspatial = NA, criterion = "BIC", ...){
  # TODO: Starting values for more complex models? (See some warnings captured from asreml)
  # Run mixed and fixed models using asreml
  require(asreml, quietly=T)

  # check validity of variable name, trait
  ok = is_valid_variable_name(trait)
  TRAIT = trait
  if (!all(ok)) {
    TRAIT[!ok] = sapply(TRAIT[!ok], function(x) paste("`",x,"`",sep=""))
  }
  covT <- FALSE
  if (!is.null(covariate)) {
    if (is.character(covariate)){
      covT <- TRUE
    }
  }
  flag <- 1
  # See if the design is regular
  if (!missing(rep))
    reptab <- with(Y,table(rep,row,col))
  else
    reptab <- with(Y,table(row,col))
  if (min(reptab) > 1) {
    warning("There must be only one plot at each REPLICATES x ROWS x COLUMNS location.\n
    spatial models will not be tried.\n")
    flag <- 0
  }
  if (min(reptab)==1 && max(reptab)==1) regular <- TRUE else regular <- FALSE

  if (is.na(rowcoordinates)) flag <-0
  if (is.na(colcoordinates)) flag <-0

  # does not use spatial models
  if (flag == 0) tryspatial <- NA

  tmp=tempfile()
  sink(file=tmp)
  # default no spatial models
  if (is.na(tryspatial)){
    if (tryrep){
      if (!is.na(checkid)){
        mr <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",rep,"+",checkid, if(covT) paste(c("",covariate), collapse="+"))),
        random=as.formula(paste("~",genotype,"+",rep,":",row,"+",rep,":",col)), rcov=~units, aom=T, data=Y, ...)
      }else{
        mr <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",rep, if(covT) paste(c("",covariate), collapse="+"))),
        random=as.formula(paste("~",genotype,"+",rep,":",row,"+",rep,":",col)), rcov=~units, aom=T, data=Y, ...)
      }

      # constrain variance of the variance components to be fixed as the values in mr
      G.param.tmp=mr$G.param
      tmppos <- which(names(G.param.tmp)==paste(rep,row,sep=":"))
      if (length(tmppos)>0)
        if (!is.null(G.param.tmp[[tmppos]][[rep]]$con))
          G.param.tmp[[tmppos]][[rep]]$con <- "F"
      tmppos <- which(names(G.param.tmp)==paste(rep,col,sep=":"))
      if (length(tmppos)>0)
        if (!is.null(G.param.tmp[[tmppos]][[rep]]$con))
          G.param.tmp[[tmppos]][[rep]]$con <- "F"
      if (!is.na(checkid)){
        mf = asreml::asreml(fixed=as.formula(paste(TRAIT,"~",rep,"+",checkid, if(covT) paste(c("",covariate), collapse="+"),"+",genotype)),
        random=as.formula(paste("~",rep,":",row,"+",rep,":",col)), rcov=~units, G.param=G.param.tmp, aom=T, data=Y,...)
      }else{
        mf = asreml::asreml(fixed=as.formula(paste(TRAIT,"~",rep, if(covT) paste(c("",covariate), collapse="+"), "+",genotype)),
        random=as.formula(paste("~",rep,":",row,"+",rep,":",col)), rcov=~units, G.param=G.param.tmp, aom=T, data=Y,...)
      }
    }else{
      if (!is.na(checkid)){
        mr <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",checkid, if(covT) paste(c("",covariate), collapse="+"))),
        random=as.formula(paste("~",genotype,"+",row,"+",col)), rcov=~units, aom=T, data=Y,...)
      }else{
        mr <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",1, if(covT) paste(c("",covariate), collapse="+"))),
        random=as.formula(paste("~",genotype,"+",row,"+",col)), rcov=~units, aom=T, data=Y,...)
      }

      # constrain variance of the variance components to be fixed as the values in mr
      G.param.tmp=mr$G.param
      if (!is.null(G.param.tmp[[row]][[row]]$con))
        G.param.tmp[[row]][[row]]$con <- "F"
      if (!is.null(G.param.tmp[[col]][[col]]$con))
        G.param.tmp[[col]][[col]]$con <- "F"
      if (!is.na(checkid)){
        mf = asreml::asreml(fixed=as.formula(paste(TRAIT,"~",checkid, if(covT) paste(c("",covariate), collapse="+"),"+",genotype)), random=as.formula(paste("~",row,"+",col)),
        rcov=~units, G.param=G.param.tmp, aom=T, data=Y,...)
      }else{
        mf = asreml::asreml(fixed=as.formula(paste(TRAIT,"~1", if(covT) paste(c("",covariate), collapse="+"), "+", genotype)), random=as.formula(paste("~",row,"+",col)),
        rcov=~units, G.param=G.param.tmp, aom=T, data=Y,...)
      }
    }
    best.model <- mr

  }else{
    if (regular){
      if (tryrep){
        random.choice <-c("Identity", "Identity", "Identity", "Measurement_error","Measurement_error","Measurement_error",
        paste(rep,row,sep=":"), paste(rep,col,sep=":"), paste(paste(rep,row,sep=":"),paste(rep,col,sep=":"),sep="+"),
        paste(paste(rep,row,sep=":"),"Measurement_error",sep="+"), paste(paste(rep,col,sep=":"),"Measurement_error",sep="+"), paste(paste(rep,row,sep=":"),paste(rep,col,sep=":"),"Measurement_error",sep="+"))

        randomterm <- c("NULL","NULL","NULL","units","units","units",
        paste(rep,row,sep=":"), paste(rep,col,sep=":"), paste(paste(rep,row,sep=":"),paste(rep,col,sep=":"),sep="+"),
        paste(paste(rep,row,sep=":"),"units",sep="+"), paste(paste(rep,col,sep=":"),"units",sep="+"), paste(paste(rep,row,sep=":"),paste(rep,col,sep=":"),"units",sep="+"))
      }else{
        random.choice <-c("Identity", "Identity", "Identity", "Measurement_error","Measurement_error","Measurement_error",
        row, col, paste(row,col,sep="+"),
        paste(row,"Measurement_error",sep="+"), paste(col,"Measurement_error",sep="+"), paste(row,col,"Measurement_error",sep="+"))

        randomterm <- c("NULL","NULL","NULL","units","units","units",
        row, col, paste(row,col,sep="+"),
        paste(row,"units",sep="+"), paste(col,"units",sep="+"), paste(row,col,"units",sep="+"))
      }
      spatial.choice <- c("AR1(x)Identity", "Identity(x)AR1", "AR1(x)AR1",
      "AR1(x)Identity", "Identity(x)AR1", "AR1(x)AR1",
      "AR1(x)Identity", "Identity(x)AR1", "AR1(x)AR1",
      "AR1(x)Identity", "Identity(x)AR1", "AR1(x)AR1")

      spatialterm <-c(paste("ar1(",rowcoordinates,"):",colcoordinates), paste(rowcoordinates,":ar1(",colcoordinates,")"), paste("ar1(",rowcoordinates,"):ar1(",colcoordinates,")"),
      paste("ar1(",rowcoordinates,"):",colcoordinates), paste(rowcoordinates,":ar1(",colcoordinates,")"), paste("ar1(",rowcoordinates,"):ar1(",colcoordinates,")"),
      paste("ar1(",rowcoordinates,"):",colcoordinates), paste(rowcoordinates,":ar1(",colcoordinates,")"), paste("ar1(",rowcoordinates,"):ar1(",colcoordinates,")"),
      paste("ar1(",rowcoordinates,"):",colcoordinates), paste(rowcoordinates,":ar1(",colcoordinates,")"), paste("ar1(",rowcoordinates,"):ar1(",colcoordinates,")"))

      model.choice <- paste("Random:", random.choice, "&   Spatial:", spatial.choice)
      for (ii in 1:length(randomterm)){
        if (tryrep){
          if (!is.na(checkid)){
            mr <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",rep,"+",checkid, if(covT) paste(c("",covariate), collapse="+"))),
            random=as.formula(paste("~",genotype,"+",randomterm[ii])), rcov=as.formula(paste("~",spatialterm[ii])), aom=T, data=Y,...)
          }else{
            mr <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",rep, if(covT) paste(c("",covariate), collapse="+"))),
            random=as.formula(paste("~",genotype,"+",randomterm[ii])), rcov=as.formula(paste("~",spatialterm[ii])), aom=T, data=Y,...)
          }
        }else{
          if (!is.na(checkid)){
            mr <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",checkid, if(covT) paste(c("",covariate), collapse="+"))),
            random=as.formula(paste("~",genotype,"+",randomterm[ii])), rcov=as.formula(paste("~",spatialterm[ii])), aom=T, data=Y,...)
          }else{
            mr <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",1, if(covT) paste(c("",covariate), collapse="+"))),
            random=as.formula(paste("~",genotype,"+",randomterm[ii])), rcov=as.formula(paste("~",spatialterm[ii])), aom=T, data=Y,...)
          }
        }
        if (ii==1){
          best.model <- mr
          best.choice <- model.choice[ii]
          best.loc <- 1
        }else{
          if (criterion == "AIC"){
            criterion.cur  <- -2*mr$loglik+2*length(mr$gammas)
            criterion.prev <- -2*best.model$loglik+2*length(best.model$gammas)
          }else{
            criterion.cur  <- -2*mr$loglik+log(length(mr$fitted.values))*length(mr$gammas)
            criterion.prev <- -2*best.model$loglik+log(length(best.model$fitted.values))*length(best.model$gammas)
          }
          if (criterion.cur < criterion.prev){
            rm(best.model)
            best.model <- mr
            best.choice <- model.choice[ii]
            best.loc <- ii
          }
        }
      }

    }else{
      if (tryspatial == "always"){
        if (tryrep){
          random.choice <-c("Identity", "Identity", "Identity", "Measurement_error","Measurement_error","Measurement_error",
          paste(rep,row,sep=":"), paste(rep,col,sep=":"), paste(paste(rep,row,sep=":"),paste(rep,col,sep=":"),sep="+"),
          paste(paste(rep,row,sep=":"),"Measurement_error",sep="+"), paste(paste(rep,col,sep=":"),"Measurement_error",sep="+"), paste(paste(rep,row,sep=":"),paste(rep,col,sep=":"),"Measurement_error",sep="+"))

          randomterm <- c("NULL","NULL","NULL","units","units","units",
          paste(rep,row,sep=":"), paste(rep,col,sep=":"), paste(paste(rep,row,sep=":"),paste(rep,col,sep=":"),sep="+"),
          paste(paste(rep,row,sep=":"),"units",sep="+"), paste(paste(rep,col,sep=":"),"units",sep="+"), paste(paste(rep,row,sep=":"),paste(rep,col,sep=":"),"units",sep="+"))
        }else{
          random.choice <-c("Identity", "Identity", "Identity", "Measurement_error","Measurement_error","Measurement_error",
          row, col, paste(row,col,sep="+"),
          paste(row,"Measurement_error",sep="+"), paste(col,"Measurement_error",sep="+"), paste(row,col,"Measurement_error",sep="+"))

          randomterm <- c("NULL","NULL","NULL","units","units","units",
          row, col, paste(row,col,sep="+"),
          paste(row,"units",sep="+"), paste(col,"units",sep="+"), paste(row,col,"units",sep="+"))
        }
        spatial.choice <- c("Exponential(x)Identity", "Identity(x)Exponential", "Isotropic exponential",
        "Exponential(x)Identity", "Identity(x)Exponential", "Isotropic exponential",
        "Exponential(x)Identity", "Identity(x)Exponential", "Isotropic exponential",
        "Exponential(x)Identity", "Identity(x)Exponential", "Isotropic exponential")

        spatialterm <-c(paste("exp(",rowcoordinates,"):",colcoordinates), paste(rowcoordinates,":exp(",colcoordinates,")"), paste("iexp(",rowcoordinates,",",colcoordinates,")"),
        paste("exp(",rowcoordinates,"):",colcoordinates), paste(rowcoordinates,":exp(",colcoordinates,")"), paste("iexp(",rowcoordinates,",",colcoordinates,")"),
        paste("exp(",rowcoordinates,"):",colcoordinates), paste(rowcoordinates,":exp(",colcoordinates,")"), paste("iexp(",rowcoordinates,",",colcoordinates,")"),
        paste("exp(",rowcoordinates,"):",colcoordinates), paste(rowcoordinates,":exp(",colcoordinates,")"), paste("iexp(",rowcoordinates,",",colcoordinates,")"))

        model.choice <- paste("Random:", random.choice, "&   Spatial:", spatial.choice)

        for (ii in 1:length(randomterm)){
          if (tryrep){
            if (!is.na(checkid)){
              mr <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",rep,"+",checkid, if(covT) paste(c("",covariate), collapse="+"))),
              random=as.formula(paste("~",genotype,"+",randomterm[ii])), rcov=as.formula(paste("~",spatialterm[ii])), aom=T, data=Y,...)
            }else{
              mr <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",rep, if(covT) paste(c("",covariate), collapse="+"))),
              random=as.formula(paste("~",genotype,"+",randomterm[ii])), rcov=as.formula(paste("~",spatialterm[ii])), aom=T, data=Y,...)
            }
          }else{
            if (!is.na(checkid)){
              mr <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",checkid, if(covT) paste(c("",covariate), collapse="+"))),
              random=as.formula(paste("~",genotype,"+",randomterm[ii])), rcov=as.formula(paste("~",spatialterm[ii])), aom=T, data=Y,...)
            }else{
              mr <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",1, if(covT) paste(c("",covariate), collapse="+"))),
              random=as.formula(paste("~",genotype,"+",randomterm[ii])), rcov=as.formula(paste("~",spatialterm[ii])), aom=T, data=Y,...)
            }
          }
          if (ii==1){
            best.model <- mr
            best.choice <- model.choice[ii]
            best.loc <- 1
          }else{
            if (criterion == "AIC"){
              criterion.cur  <- -2*mr$loglik+2*length(mr$gammas)
              criterion.prev <- -2*best.model$loglik+2*length(best.model$gammas)
            }else{
              criterion.cur  <- -2*mr$loglik+log(length(mr$fitted.values))*length(mr$gammas)
              criterion.prev <- -2*best.model$loglik+log(length(best.model$fitted.values))*length(best.model$gammas)
            }
            if (criterion.cur < criterion.prev){
              rm(best.model)
              best.model <- mr
              best.choice <- model.choice[ii]
              best.loc <- ii
            }
          }
        }
      }
    }

    # constrain variance of the variance components to be fixed as the values in mr
    if (tryrep){
      G.param.tmp=best.model$G.param
      tmppos <- which(names(G.param.tmp)==paste(rep,row,sep=":"))
      if (length(tmppos)>0)
        if (!is.null(G.param.tmp[[tmppos]][[rep]]$con))
          G.param.tmp[[tmppos]][[rep]]$con <- "F"
      tmppos <- which(names(G.param.tmp)==paste(rep,col,sep=":"))
      if (length(tmppos)>0)
        if (!is.null(G.param.tmp[[tmppos]][[rep]]$con))
          G.param.tmp[[tmppos]][[rep]]$con <- "F"
      if (!is.na(checkid)){
        mf <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",rep,"+",checkid, if(covT) paste(c("",covariate), collapse="+"),"+",genotype)),
        random=as.formula(paste("~",randomterm[best.loc])), rcov=as.formula(paste("~",spatialterm[ii])), G.param=G.param.tmp, aom=T, data=Y,...)
      }else{
        mf <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",rep, if(covT) paste(c("",covariate), collapse="+"),"+",genotype)),
        random=as.formula(paste("~",randomterm[best.loc])), rcov=as.formula(paste("~",spatialterm[ii])), G.param=G.param.tmp, aom=T, data=Y,...)
      }
    }else{
      G.param.tmp=best.model$G.param
      if (!is.null(G.param.tmp[[row]][[row]]$con))
        G.param.tmp[[row]][[row]]$con <- "F"
      if (!is.null(G.param.tmp[[col]][[col]]$con))
        G.param.tmp[[col]][[col]]$con <- "F"
      if (!is.na(checkid)){
        mf <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~",checkid, if(covT) paste(c("",covariate), collapse="+"),"+",genotype)),
        random=as.formula(paste("~",randomterm[best.loc])), rcov=as.formula(paste("~",spatialterm[ii])), G.param=G.param.tmp, aom=T, data=Y,...)
      }else{
        mf <- asreml::asreml(fixed=as.formula(paste(TRAIT,"~1", if(covT) paste(c("",covariate), collapse="+"), "+",genotype)),
        random=as.formula(paste("~",randomterm[best.loc])), rcov=as.formula(paste("~",spatialterm[ii])), G.param=G.param.tmp, aom=T, data=Y,...)
      }
    }
  }

  # run predict
  if (!is.na(tryspatial)) ii <- best.loc
  best.model$call$fixed <- eval(best.model$call$fixed)
  best.model$call$random <- eval(best.model$call$random)
  best.model$call$rcov <- eval(best.model$call$rcov)

  mf$call$fixed <- eval(mf$call$fixed)
  mf$call$random <- eval(mf$call$random)
  mf$call$rcov <- eval(mf$call$rcov)
  best.model = predict(best.model, classify=genotype, data=Y)
  if (!is.na(checkid)){
    mf = predict(mf, classify=genotype, vcov=T, data=Y, associate=~checkid:genotype)
  }else{
    mf = predict(mf, classify=genotype, vcov=T, data=Y)
  }
  sink()
  unlink(tmp)

  res = list(mmix = best.model, mfix = mf, Data = Y)
  return(res)
}
