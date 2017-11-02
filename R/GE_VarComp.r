#' Selects the best variance-covariance model for a set of environments
#'
#' This function selects the best covariance structure for genetic correlations between environments. It fits a
#' range of variance-covariance models to compare (e.g., identity, compound symmetry,
#' diagonal, heterogeneous compound symmetry, first order factor analysis, second order factor analysis, unstructured),
#' and selects the best one using a goodness-of-fit criterion.
#'
#' @param Y A data frame object.
#' @param trait A character string specifying the name for a trait found in \code{Y}.
#' @param genotype A character string specifying the name for the genotypes found in \code{Y}.
#' @param env A character string specifying the name for enviroments found in \code{Y}.
#' @param engine A string specifying the name of a mixed modelling engine to be used.
#' @param criterion A string specifying a goodness-of-fit criterion, i.e., "AIC" or "BIC".
#' @param ... Further arguments to be passed to \code{asreml}.
#' @note If \code{engine="lme4"}, only the compound symmetry model can be fitted.
#' @return A list object consisting of the fitted model objects, a string specifying
#' the best model and its related goodness-of-fit criterion.
#' @examples
#' mydat <- GE.read.csv(file.path(path.package("RAP"),"F2maize_pheno.csv"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' names(mydat)=c("env", "genotype","yld")
#' model1 <- GE.VarComp(mydat, trait="yld", genotype="genotype", env="env",
#'                      engine = "lme4", #engine = "asreml",
#'                      criterion = "BIC")
#' model1$BIC
#' model1$choice
#' summary(model1$model[[model1$choice]], nice=TRUE)
#'
#' @export
GE.VarComp <- function(Y,
                       trait,
                       genotype,
                       env,
                       engine = "asreml",
                       criterion = "BIC",
                       ...) {
  # keep NA values
  Y <- droplevels(Y[,c(env, genotype, trait)])
  na.entries <- which(is.na(tapply(rep(1, nrow(Y)), Y[, c(genotype,env)], identity)),  arr.ind =T)
  if (length(na.entries)){
    tempnames <- names(Y)
    Y.ext <- data.frame(Y[na.entries[,env], env], Y[na.entries[,genotype], genotype], rep(NA, nrow(na.entries)))
    names(Y.ext) <- tempnames
    Y <- rbind(Y, Y.ext)
  }

  #sort data in order of genotype, env
  Y <- Y[order(Y[[genotype]], Y[[env]]),]

  qvinitial <- function(Y, trait, genotype, env, uniterror = NA,
  vcmodel=c("identity", "cs", "diagonal", "hcs", "outside", "fa", "fa2", "unstructured"), fixed=NULL, unitfactor=NA, ...){
    # Replicates '_QVINITIAL' procedure (S. J. Welham 15/05/09) in GenStat
    # TODO: factanal() ? fa()

    # First, form estimate of unstructured matrix

    # exclude the rows including NA
    X <- na.omit(Y[,c(trait,genotype,env)])

    nenv <- nlevels(X[,env])
    ngen <- nlevels(X[,genotype])

    # Get fixed df by stealth - in absence of other info, share among environments
    P <- 1
    if (!is.null(fixed)){
       mr <- asreml::asreml(fixed=fixed, rcov = ~id(units), data=X,...)
       P <- length(mr$fitted.values)-(1+mr$nedf)
    }

    # Get number of effects contributing to each sum of squares
    tmptable <- table(X[,c(env,genotype)])
    nobs_env <- rowSums(tmptable)
    Rnobs <- Cnobs <- matrix(nobs_env, nrow=nenv, ncol=nenv,byrow=T)
    Cnobs <- t(Rnobs)
    Nobs <- Rnobs*(Rnobs-Cnobs<0) + Cnobs*(Cnobs-Rnobs<=0)

    if (!is.na(uniterror)){
      # This case is trickier, becuase of partitioning between two random tems,
      # but use of the diag structure is better than nothing!

      weights <- 1/uniterror
      X["weights"] <- weights
      if (!is.null(fixed)){
        init.values <- asreml::asreml(fixed=fixed,random=as.formula(paste("~",genotype,":idh(",env,")")), weights=weights, start.values=T, data=X, ...)
        tmp <- init.values$gammas.table
        tmp[,"Constraint"] <- as.character(tmp[,"Constraint"] )
        tmp[which(tmp[,"Gamma"]=="R!variance"),"Constraint"] <- "F"
        tmp[,"Constraint"] <- as.factor(tmp[,"Constraint"])
        mr <- asreml::asreml(fixed=fixed, random=as.formula(paste("~",genotype,":idh(",env,")")), weights=weights, R.param=tmp, data=X, ...)
        Ores <- matrix(mr$coeff$random, ngen, nenv, byrow=T)
      }else{
        init.values <- asreml::asreml(fixed=as.formula(paste(trait,"~", env)),random=as.formula(paste("~",genotype,":idh(",env,")")), weights=weights,
        start.values=T, data=X, ...)
        tmp <- init.values$gammas.table
        tmp[,"Constraint"] <- as.character(tmp[,"Constraint"] )
        tmp[which(tmp[,"Gamma"]=="R!variance"),"Constraint"] <- "F"
        tmp[,"Constraint"] <- as.factor(tmp[,"Constraint"])
        mr <- asreml::asreml(fixed=as.formula(paste(trait,"~", env)), random=as.formula(paste("~",genotype,":idh(",env,")")), weights=weights, R.param=tmp,
        data=X, ...)
        Ores <- matrix(mr$coeff$random, ngen, nenv, byrow=T)
      }
    }else{
      # This gives correct answers for complete balanced data
      if (!is.null(fixed)){
         mr <- asreml::asreml(fixed=fixed, data=X, ...)
         residuals <- mr$residuals
      }else{
         mr <- asreml::asreml(fixed=as.formula(paste(trait,"~", env)), data=X, ...)
         residuals <- mr$residuals
      }
      Ores <- tapply(residuals, list(X[,genotype], X[,env]), mean)
      if (sum(is.na(Ores))!=0)  Ores[which(is.na(Ores))]=0
    }
    Rmat <- Ores
    Evcov <- matrix(,nenv,nenv)
    Evcov <- t(Rmat)%*%Rmat/(Nobs-(P/nenv))

    # Get off-diagonal elements of Evcov only
    Offdiag <- Evcov
    diag(Offdiag) <- NA
    # Get off-diagonal elements of cor(Evcov) only
    Offcorr <- cor(Evcov)
    diag(Offcorr) <- NA

    #Get initial values for each model
    if (vcmodel == "identity"){
      vcinitial <- list(vge = mean(diag(Evcov)))
    }
    if (vcmodel == "cs"){
      vg <- mean(Offdiag, na.rm=T)
      vge <- mean(diag(Evcov))
      vge <- (vge>vg)*(vge-vg) + (vge<=vg)*0.1*vge
      vcinitial <- list(vg=vg, vge=vge)
    }
    if (vcmodel == "diagonal"){
      vcinitial <- list(Diag=diag(Evcov))
    }
    if (vcmodel == "hcs"){
      vg <- mean(Offdiag, na.rm=T)/2
      Diag <- diag(Evcov)
      Diag <- (Diag>vg)*(Diag-vg) + (Diag<=vg)*0.1*Diag
      vcinitial <- list(vg=vg, Diag=Diag)
    }
    if (vcmodel == "outside"){
      vg <- mean(Offcorr, na.rm=T)
      Diag <- diag(Evcov)
      vcinitial <- list(vg=vg, Diag=Diag)
    }
    if (vcmodel == "fa"){
      ok <- require(psych, quietly=T)
      if (ok){
        factor.analysis <- try(psych::fa(r=Evcov, nfactors=1, fm="mle"),T)
        if (inherits(factor.analysis,"try-error")){
          factor.analysis <- psych::fa(r=Evcov, nfactors=1, fm="minres")
        }
        Loading <- as.vector(factor.analysis$loadings)
        comm <- factor.analysis$comm
        Var <- diag(Evcov)
        Psi <- Var*(1-comm)
        Smat <- Loading%*%t(Loading)
        Chol <- chol(Smat, pivot = TRUE)
        Gamma <- Chol[1,]
        vcinitial <- list(Gamma=Gamma, Psi=Psi)
      }else{
        vcinitial <- NULL
        warning("psych package is required but failed to load!\n")
      }
    }
    if (vcmodel == "fa2"){
      ok1 <- require(psych, quietly=T)
      ok2 <- require(GPArotation, quietly=T)
      if (ok1&&ok2){
        factor.analysis <- try(psych::fa(r=Evcov, nfactors=2, fm="mle"),T)
        if (inherits(factor.analysis,"try-error")){
          factor.analysis <- psych::fa(r=Evcov, nfactors=2, fm="minres")
        }
        Loading <- as.vector(factor.analysis$loadings[,1])
        comm <- factor.analysis$comm
        Var <- diag(Evcov)
        Psi <- Var*(1-comm)
        Smat <- Loading%*%t(Loading)
        Chol <- chol(Smat, pivot = TRUE)
        Gamma <- Chol[1:2,]
        vcinitial <- list(Gamma=Gamma, Psi=Psi)
      }else{
        vcinitial <- NULL
        warning("psych and GPArotation packages are required but failed to load!\n")
      }
    }
    if (vcmodel == "unstructured"){
      vcinitial <- list(Evcov=Evcov)
    }
    vcinitial
  }

  # Main procedure to fit mixed models
  if (engine == "lme4"){
    # Compound symmetry ("cs") only
    ok <- require(lme4, quietly = T)
    if (ok){
      mr = lme4::lmer(as.formula(paste(trait,"~ ",env," +(1 |", genotype,")")), data = Y, ...)
      n.par <- 2
      # Outputs
      res <- new.env()
      res$model$cs <- mr
      res$choice <- "cs"
      if (criterion == "AIC"){
        res$AIC.best <- -2*as.numeric(logLik(mr))+2*n.par
      }else{
        res$BIC.best <- -2*as.numeric(logLik(mr))+log(length(fitted(mr)))*n.par
      }
    }else{
      stop("Fail to load 'lme4'.\n")
    }
  }
  if (engine == "asreml"){
    ok <- require(asreml, quietly = T)
    if (ok){
      choices <- c("identity", "cs", "diagonal", "hcs", "outside", "fa", "fa2", "unstructured")
      BestTab <- matrix(,nrow=8,ncol=4)
      colnames(BestTab) <- c("AIC", "BIC",	"Deviance", "NParameters")
      rownames(BestTab) <- choices
      res <- new.env()
      for (i in 1:length(choices)){
        if (choices[i] == "identity"){
          mr <- asreml::asreml(fixed=as.formula(paste(trait,"~", env)), rcov = as.formula(paste("~",genotype,":",env)), data=Y, ...)
          mr$call$fixed <- eval(mr$call$fixed)
          mr$call$rcov <- eval(mr$call$rcov)
          res$model[["identity"]] <- mr
          if (!mr$converge) mr$loglik <- -Inf
          n.par <- 1
        }
        if (choices[i] == "cs"){
          mr <- asreml::asreml(fixed=as.formula(paste(trait,"~", env)), random = as.formula(paste("~",genotype)), rcov=as.formula(paste("~",genotype,":",env)), data=Y, ...)
          mr$call$fixed <- eval(mr$call$fixed)
          mr$call$random <- eval(mr$call$random)
          mr$call$rcov <- eval(mr$call$rcov)
          res$model[["cs"]] <- mr
          if (!mr$converge) mr$loglik <- -Inf
          n.par <- 2
        }
        if (choices[i] == "diagonal"){
          mr <- try(asreml::asreml(fixed=as.formula(paste(trait,"~", env)), rcov = as.formula(paste("~",genotype,":diag(",env,")")), data=Y, ...), silent=TRUE)
          if (inherits(mr, "try-error")){
            mr <- vector("list")
            mr$loglik <- -Inf
          }else{
            mr$call$fixed <- eval(mr$call$fixed)
            mr$call$rcov <- eval(mr$call$rcov)
            res$model[["diagonal"]] <- mr
            if (!mr$converge) mr$loglik <- -Inf
          }
          n.par <- length(levels(Y[,env]))
        }
        if (choices[i] == "hcs"){
          mr <- try(asreml::asreml(fixed=as.formula(paste(trait,"~", env)), random = as.formula(paste("~",genotype)), rcov =
          as.formula(paste("~",genotype,":diag(",env,")")), data=Y, ...), silent=TRUE)
          if (inherits(mr, "try-error")){
            mr <- vector("list")
            mr$loglik <- -Inf
          }else{
            mr$call$fixed <- eval(mr$call$fixed)
            mr$call$random <- eval(mr$call$random)
            mr$call$rcov <- eval(mr$call$rcov)
            res$model[["hcs"]] <- mr
            if (!mr$converge) mr$loglik <- -Inf
          }
          n.par <- length(levels(Y[,env]))+1
        }
        if (choices[i] == "outside"){
          init.vals <- asreml::asreml(fixed=as.formula(paste(trait,"~", env)), random = as.formula(paste("~",genotype,":corh(",env,")")), start.values=T, data=Y, ...)
          tmp.table <- init.vals$gammas.table
          tmp=tempfile()
          sink(file=tmp)
          tmp.values <- qvinitial(Y, trait, genotype, env, uniterror = NA, vcmodel="outside", fixed=as.formula(paste(trait,"~", env)), unitfactor=NA, ...)
          sink()
          unlink(tmp)
          #print(tmp.table)
          tmp.table[,"Value"] <-c(c(tmp.values$vg,tmp.values$Diag), 1)
          tmp.table[,"Constraint"] <- as.character(tmp.table[,"Constraint"] )
          tmp.table[which(tmp.table[,"Gamma"]=="R!variance"),"Constraint"] <- "F"
          tmp.table[,"Constraint"] <- as.factor(tmp.table[,"Constraint"])
          #print(tmp.table)
          mr <- try(asreml::asreml(fixed=as.formula(paste(trait,"~",env)), random = as.formula(paste("~",genotype,":corh(",env,")")),
          G.param=tmp.table, R.param=tmp.table, data=Y, ...), silent=TRUE)
          if (inherits(mr, "try-error")){
            mr <- vector("list")
            mr$loglik <- -Inf
          }else{
            mr$call$fixed <- eval(mr$call$fixed)
            mr$call$random <- eval(mr$call$random)
            mr$call$rcov <- eval(mr$call$rcov)
            mr$call$G.param <- eval(mr$call$G.param)
            mr$call$R.param <- eval(mr$call$R.param)
            res$model[["outside"]] <- mr
            if (!mr$converge) mr$loglik <- -Inf
          }
          n.par <- length(levels(Y[,env]))+1
        }

        if (nlevels(Y[[env]])>4){
          if (choices[i] == "fa"){
            init.vals <- asreml::asreml(fixed=as.formula(paste(trait,"~",env)), random = as.formula(paste("~",genotype,":fa(",env,",1)")),
             start.values=T, data=Y, ...)
            tmp.table <- init.vals$gammas.table
            tmp=tempfile()
            sink(file=tmp)
            tmp.values <- qvinitial(Y, trait, genotype, env, uniterror = NA, vcmodel="fa", fixed=as.formula(paste(trait,"~",env)), unitfactor=NA, ...)
            sink()
            unlink(tmp)
            if (is.null(tmp.values)){
               mr <- try(asreml::asreml(fixed=as.formula(paste(trait,"~",env)), random = as.formula(paste("~",genotype,":fa(",env,",1)")),
               data=Y, ...), silent=TRUE)
            }else{
              tmp.table[,"Value"] <-c(tmp.values$Psi,tmp.values$Gamma, 1)
              tmp.table[,"Constraint"] <- as.character(tmp.table[,"Constraint"] )
              tmp.table[which(tmp.table[,"Gamma"]=="R!variance"),"Constraint"] <- "F"
              tmp.table[,"Constraint"] <- as.factor(tmp.table[,"Constraint"])
              mr <- try(asreml::asreml(fixed=as.formula(paste(trait,"~",env)), random = as.formula(paste("~",genotype,":fa(",env,",1)")),
               R.param=tmp.table, G.param=tmp.table,  data=Y, ...), silent=TRUE)
            }
            if (inherits(mr, "try-error")){
              mr <- vector("list")
              mr$loglik <- -Inf
            }else{
              mr$call$fixed <- eval(mr$call$fixed)
              mr$call$random <- eval(mr$call$random)
              mr$call$rcov <- eval(mr$call$rcov)
              mr$call$G.param <- eval(mr$call$G.param)
              mr$call$R.param <- eval(mr$call$R.param)
              res$model[["fa"]] <- mr
              if (!mr$converge) mr$loglik <- -Inf
            }
            n.par <- length(levels(Y[,env]))*2
          }
          if (choices[i] == "fa2"){
            init.vals <- asreml::asreml(fixed=as.formula(paste(trait,"~",env)), random = as.formula(paste("~",genotype,":fa(",env,",2)")),
             start.values=T,  data=Y, ...)
            tmp.table <- init.vals$gammas.table
            tmp=tempfile()
            sink(file=tmp)
            tmp.values <- qvinitial(Y, trait, genotype, env, uniterror = NA, vcmodel="fa2", fixed=as.formula(paste(trait,"~",env)), unitfactor=NA, ...)
            sink()
            unlink(tmp)
            if (is.null(tmp.values)){
               mr <- try(asreml::asreml(fixed=as.formula(paste(trait,"~",env)), random = as.formula(paste("~",genotype,":fa(",env,",2)")),
               data=Y, ...), silent=TRUE)
            }else{
              # Keep loadings of factor 2 away from 0
              tmp.values$Gamma[2,which(tmp.values$Gamma[2,] < 1e-3)] <- 1e-3
              # Make sure that first entry is 0
              tmp.values$Gamma[2,1] <- 0
              tmp.table[,"Value"] <-c(tmp.values$Psi,tmp.values$Gamma[1,],tmp.values$Gamma[2,], 1)
              tmp.table[,"Constraint"] <- as.character(tmp.table[,"Constraint"] )
              tmp.table[which(tmp.table[,"Gamma"]=="R!variance"),"Constraint"] <- "F"
              tmp.table[,"Constraint"] <- as.factor(tmp.table[,"Constraint"])
              mr <- try(asreml::asreml(fixed=as.formula(paste(trait,"~",env)), random = as.formula(paste("~",genotype,":fa(",env,",2)")),
               R.param=tmp.table, G.param=tmp.table,  data=Y, ...), silent=TRUE)
            }
            if (inherits(mr, "try-error")){
              mr <- vector("list")
              mr$loglik <- -Inf
            }else{
              mr$call$fixed <- eval(mr$call$fixed)
              mr$call$random <- eval(mr$call$random)
              mr$call$rcov <- eval(mr$call$rcov)
              mr$call$G.param <- eval(mr$call$G.param)
              mr$call$R.param <- eval(mr$call$R.param)
              res$model[["fa2"]] <- mr
              if (!mr$converge) mr$loglik <- -Inf
            }
            n.par <- length(levels(Y[,env]))*3-1
          }
        }
        # Check model!
        if (choices[i] == "unstructured"){
          init.vals <- asreml::asreml(fixed=as.formula(paste(trait,"~",env)), random = as.formula(paste("~",genotype,":us(",env,")")), start.values=T, data=Y, ...)
          tmp.table <- init.vals$gammas.table
          tmp=tempfile()
          sink(file=tmp)
          tmp.values <- qvinitial(Y, trait, genotype, env, uniterror = NA, vcmodel="unstructured", fixed=as.formula(paste(trait,"~",env)), unitfactor=NA, ...)
          sink()
          unlink(tmp)
          tmp.values <- tmp.values$Evcov[upper.tri(tmp.values$Evcov, diag=T)]
          tmp.table[,"Value"] <-c(tmp.values, 1)
          tmp.table[,"Constraint"] <- as.character(tmp.table[,"Constraint"] )
          tmp.table[which(tmp.table[,"Gamma"]=="R!variance"),"Constraint"] <- "F"
          tmp.table[,"Constraint"] <- as.factor(tmp.table[,"Constraint"])
          mr <- try(asreml::asreml(fixed=as.formula(paste(trait,"~",env)), random = as.formula(paste("~",genotype,":us(",env,")")),
          G.param=tmp.table, R.param=tmp.table, data=Y, ...), silent=TRUE)
          if (inherits(mr, "try-error")){
            mr <- vector("list")
            mr$loglik <- -Inf
          }else{
            mr$call$fixed <- eval(mr$call$fixed)
            mr$call$random <- eval(mr$call$random)
            mr$call$rcov <- eval(mr$call$rcov)
            mr$call$G.param <- eval(mr$call$G.param)
            mr$call$R.param <- eval(mr$call$R.param)
            res$model[["unstructured"]] <- mr
            if (!mr$converge) mr$loglik <- -Inf
          }
          n.par <- length(levels(Y[,env]))*(length(levels(Y[,env]))-1)/2 + length(levels(Y[,env]))
        }
        if (!(nlevels(Y[[env]])<=4 && choices[i] %in% c("fa", "fa2"))){
          BestTab[choices[i], "AIC"] <-  -2*mr$loglik+2*n.par
          BestTab[choices[i], "BIC"] <-  -2*mr$loglik+log(length(mr$fitted.values))*n.par
          BestTab[choices[i], "Deviance"] <- -2*mr$loglik
          BestTab[choices[i], "NParameters"] <- n.par
          if (i==1){
            best.model <- mr
            best.choice <- choices[i]
            n.best <- n.par
          }else{
            if (criterion == "AIC"){
              criterion.cur  <- -2*mr$loglik+2*n.par
              criterion.prev <- -2*best.model$loglik+2*n.best
            }else{
              criterion.cur  <- -2*mr$loglik+log(length(mr$fitted.values))*n.par
              criterion.prev <- -2*best.model$loglik+log(length(best.model$fitted.values))*n.best
            }
            if (criterion.cur < criterion.prev){
              rm(best.model)
              best.model <- mr
              best.choice <- choices[i]
              n.best <- n.par
            }
          }
        }
      }
      if (criterion == "AIC")
        BestTab <- BestTab[order(BestTab[, "AIC"]), ]
      else
        BestTab <- BestTab[order(BestTab[, "BIC"]), ]

      # Outputs
      res$choice <- best.choice
      res$summaryTab <- BestTab
      tempFile <- tempfile()
      sink(tempFile)
      res$vcov.best <- predict(best.model, classify=env, data=Y, vcov=T)$predictions$vcov
      sink(NULL)
      unlink(tempFile)
      colnames(res$vcov.best) <- rownames(res$vcov.best) <- levels(Y[,env])
      if (criterion == "AIC"){
        res$AIC.best <- min(criterion.cur,criterion.prev)
      }else{
        res$BIC.best <- min(criterion.cur,criterion.prev)
      }
    }else{
      stop("Fail to load 'asreml'.\n")
    }
  }
  if (engine=="asreml"||engine=="lme4")
    as.list(res)
  else
    stop("Please use either 'asreml' or 'lme4' as the mixing modelling engine.\n")
}
