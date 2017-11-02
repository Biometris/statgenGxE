#' Forms predicted means (BLUPS) and produces tables based on a set of mega-environments
#'
#' This function calculates predicted meanes (BLUPS) and associated standard errors based on a set of mega-environments
#' @param trait A character string specifying a trait column of the data.
#' @param genotype A character string specifying a genotype column of the data.
#' @param env A character string specifying an environment column of the data.
#' @param year A character string specifying years within environments.
#' @param megaenv A character string specifying the mega-environment factor.
#' @param data A data frame object.
#' @param ... Other parameters passed to either \code{asreml()} or \code{lmer()}.
#' @examples
#' mydat <- GE.read.csv(file.path(path.package("RAP"),"F2maize_pheno.csv"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' names(mydat)=c("env", "genotype","yld") 
#' Y <- GE.megaenvironment(Y=mydat, trait="yld", genotype="genotype", 
#'                         env="env", megaenv="megaenv")
#' GE.GxETable(trait="yld", genotype="genotype", env="env", year=NULL, 
#'             megaenv="megaenv", data=Y)
#'             
#' @export             

GE.GxETable <- function(trait, genotype, env, year=NULL, megaenv, data,...){
  
  if (!trait %in% names(data))
    stop(trait," not found in ", data)
  if (!genotype %in% names(data))
    stop(genotype, " not found in ", data)
  if (!env %in% names(data))
    stop(env, " not found in ", data)
  if (!megaenv %in% names(data))
    stop(megaenv, " not found in ", data)
  
  ok <- require(asreml, quietly = T)
  if (ok){
    if (is.null(year)){
#      sv <- asreml::asreml(fixed=as.formula(paste(trait,"~",env)),
#                   random <- as.formula(paste("~",genotype,":us(",megaenv,")")),
#                   data=data, start.values=T,...)
#      sv$R.param$R$variance$con <- "F"
      mr <- try(asreml::asreml(fixed=as.formula(paste(trait,"~",env)),
                       random=as.formula(paste("~",genotype,":us(",megaenv,")")),
                       data=data,...), silent=TRUE)
    }else{
#      sv <- asreml::asreml(fixed=as.formula(paste(trait,"~",env,"/",year)),
#                   random <- as.formula(paste("~",genotype,":us(",megaenv,")+",genotype,":",megaenv,":",year)),
#                   data=data, start.values=T,...)
#      sv$R.param$R$variance$con <- "F"
      mr <- try(asreml::asreml(fixed=as.formula(paste(trait,"~",env,"/",year)),
                       random=as.formula(paste("~",genotype,":us(",megaenv,")+",genotype,":",megaenv,":",year)),
                       data=data,...), silent=TRUE)
    }
    if (inherits(mr, "try-error")){
      genolevels <- levels(data[[genotype]])
      ngeno.unique <- nlevels(data[[genotype]])
      megaenvlevels <- levels(data[[megaenv]])
      nmegaenv.unique <- nlevels(data[[megaenv]])
      predVals <- se <-  data.frame(matrix(, nrow=ngeno.unique, ncol=nmegaenv.unique), row.names =  genolevels, check.names = FALSE)
      names(predVals) <- names(se) <- megaenvlevels
    }else{
      mr$call$fixed <- eval(mr$call$fixed)
      mr$call$random <- eval(mr$call$random)
      mr$call$rcov <- eval(mr$call$rcov)
      mr$call$R.param <- eval(mr$call$R.param)
      mr <- predict(mr, classify=paste(genotype,":",megaenv,sep=""), data=data)
      predictions <- mr$predictions$pvals
      if (!is.factor(predictions[[megaenv]]))
        predictions[[megaenv]] <- as.factor(predictions[[megaenv]])
      if (!is.factor(predictions[[genotype]]))
        predictions[[genotype]] <- as.factor(predictions[[genotype]])
      predVals <- tapply(X=predictions$predicted.value, INDEX=predictions[,c(genotype, megaenv)], FUN=identity)
      se <- tapply(X=predictions$standard.error, INDEX=predictions[,c(genotype, megaenv)], FUN=identity)
    }
  }else{
    ok <- require(lme4, quietly = T)
    if (ok){
      if (is.null(year)){
        mr <- try(lme4::lmer(as.formula(paste(trait,"~",env, "+ (0+",megaenv,"|",genotype, ")")), data=data, ...), silent=TRUE)
      }else{
        mr <- try(lme4::lmer(as.formula(paste(trait,"~",env,"/",year,
        "+ (0+",megaenv,"|",genotype, ") + (0+",megaenv,"|",genotype,":",year,")")),
        data=data, ...), silent=TRUE)
      }
      genolevels <- levels(data[[genotype]])
      ngeno.unique <- nlevels(data[[genotype]])
      megaenvlevels <- levels(data[[megaenv]])
      nmegaenv.unique <- nlevels(data[[megaenv]])
      if (inherits(mr, "try-error")){
        predVals <- se <-  data.frame(matrix(, nrow=ngeno.unique, ncol=nmegaenv.unique), row.names =  genolevels, check.names = FALSE)
        names(predVals) <- names(se) <- megaenvlevels
      }else{
        # Extract coeffcients mr
        if(class(mr) == 'lmerMod') fe = lme4::fixef(mr)
        if(class(mr) == 'mer') fe = slot(mr,"fixef")
        rr = grep(env, names(fe))
        cr = fe[rr]
        ng = length(unique(slot(mr,"flist")[[genotype]]))
        reff = lme4::ranef(mr,drop=T)[[genotype]]
        blo = mean(c(cr, 0))
        
        # Predictions BLUPs
        predictions.blups = fe[1] + blo + reff
        
        # Compute se.blups
        if(class(mr) == 'lmerMod') se.blups = t(sqrt(apply(attr(lme4::ranef(mr,condVar=T)[[genotype]], "postVar"),3, diag)))
        if(class(mr) == 'mer') se.blups = t(sqrt(apply(attr(lme4::ranef(mr,postVar=T)[[genotype]], "postVar"),3, diag)))
        
        predVals <- predictions.blups
        se <- data.frame(se.blups, row.names = rownames(predVals), check.names = FALSE)
        names(se) <- names(predVals) <- megaenvlevels
      }
    }else{
      stop("Either asreml or lme4 is not loaded correctly.")
    }
  }
  
  result <- new.env()
  result$predicted.value <- predVals
  result$standard.error <- se
  result <- as.list(result)
  result
}
