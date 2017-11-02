#' Line x Tester Analysis
#'
#' This function does a mixed-model analysis of data from a line-by-tester trial, using either \code{asreml} or \code{lme4}.
#' @param fixed a formula specifying fixed model terms, in addition to the \code{testers} main effect and any control comparisons.
#' @param random a formula specifying random model terms, in addition to the terms involving \code{lines} and \code{testers} that
#' are included automatically.
#' @param lines a character string specifying the line (usually female parent) in \code{data}.
#' @param testers a character string specifying the tester (usually male parent) in \code{data}.
#' @param controls a character string specifying a factor in \code{data}, which distinguishes between control
#' and test (line x tester) genotypes.
#' @param data a data frame object containing data of a line-by-tester trial.
#' @param method a character string specifying the criterion (either \code{"aic"} or \code{"bic"}) to look for the best model
#' if \code{asreml} or \code{lme4} cannot fit the current specifed model; default, \code{"bic"}.
#' @param engine a string specifies the name of a mixed modelling engine (either \code{"asreml"} or \code{"lme4"}); default, \code{"asreml"}.
#' @param na.method.Y a character string (\code{"include"}, \code{"omit"} or \code{"fail"}) specifying how missing data in the response is to be handled.
#' This is applied to the model frame after any subset argument has been used. The default (\code{"include"}) is to estimate missing values;
#' this is necessary in spatial models to preserve the spatial structure. \code{"omit"} deletes observations that contain missing values in the response.
#' Note this is only available when \code{asreml} is used.
#' @param na.method.X a character string (\code{"include"}, \code{"omit"} or \code{"fail"}) specifying how missing data in covariates is to be handled.
#' This is applied to the model frame after any subset argument has been used but before any \code{at()} functions are invoked.
#' The default \code{"fail"} causes an error if there are missing values in any covariate.
#' \code{"omit"} deletes observations that contain missing values in any covariate.
#' Note this is only available when \code{asreml} is used.
#' @param na.action a function that indicates what should happen when the data contain NAs.
#' The default action (\code{na.omit}) strips any observations with any missing values in any variables.
#' This is only available when \code{lme4} is used as a mixed modelling engine.
#' @param recover logical. Whether to try to recover with a simpler random model if REML cannot fit the model; default \code{FALSE}.
#' @param refit logical indicating if objects of class \code{lmerMod} should be refitted with ML before comparing models.
#' The default is \code{TRUE} to prevent the common mistake of inappropriately comparing REML-fitted models with different fixed effects,
#' whose likelihoods are not directly comparable. This option is available for \code{lme4} only.
#' @param verbose logical. Whether to display a summary of the line-by-tester analysis on screen; Default, \code{TRUE}.
#' @param ... other parameters to be passed on either \code{asreml} or \code{lme4}.
#' @return
#' a list consists of the fitted model object (\code{lxtmodel}), a data frame object containing the GCA effects (\code{blups.line}),
#' a data frame object containing SCA effects (\code{blups.line.tester})
#' and a data frame object containing tests for combinability effects (\code{test.combinability}).
#' @examples
#' mydat <- ST.read.csv(file.path(path.package("RAP"),"VLIN-1.csv"),
#'                      factor.names=c("Replicates","Blocks","Controls","Lines","Testers"),
#'                      trait.names="yield", env ="Env")
#' #library(lme4)
#' lxt1 <- ST.lineXtester(fixed=yield~1, random=~Replicates/Blocks, lines="Lines",
#'                        testers="Testers", controls="Controls", data=mydat, engine="lme4")
#' #library(asreml)
#' lxt2 <- ST.lineXtester(fixed=yield~1, random=~Replicates/Blocks, lines="Lines",
#'                        testers="Testers", controls="Controls",
#'                        data=mydat, engine="asreml", maxiter=30)
#' lxt3 <- ST.lineXtester(fixed=yield~1, random=~Replicates/Blocks, lines="Lines",
#'                        testers="Testers", data=mydat)
#'
#' @export
ST.lineXtester <- function(fixed, random, lines, testers, controls, data, method ="bic", engine="asreml",
           na.method.Y = "include", na.method.X = "fail", na.action=na.omit, recover=FALSE, refit = TRUE, verbose=TRUE,...){

  if (missing(fixed))
    stop("a fixed part formula needed.")
  else{
    fterms <- terms(fixed, keep.order = TRUE)
    if (length(attr(fterms, "factors")) == 0L) {
      resp <- as.character(fterms)[2]
      if (attr(fterms, "intercept") == 0L)
        Fterms <- "-1"
      else
        Fterms <- NULL
    }else{
      resp.loc <- attr(fterms, "response")
      resp <- rownames(attr(fterms, "factors"))[resp.loc]
      tFterms <- Fterms <- colnames(attr(fterms, "factors"))
      tFterms <- unlist(strsplit(tFterms,"\\:"))
      tFterms <- gsub("^[[:space:]]","",tFterms)
      tFterms <- gsub("\\s$","",tFterms)
      tFterms <- unique(tFterms)
      if (!all(tFterms %in% names(data)))
        stop(paste(tFterms[!(tFterms %in% names(data))], collapse=","), " not found in column names of data")
      if (attr(fterms, "intercept") == 0L)
        Fterms <- c("-1", Fterms)
    }
    if (!(resp %in% names(data)))
        stop(resp, " not found in column names of data")
  }

  if (!missing(random)){
    rterms <- terms(random, keep.order = TRUE)
    Rterms <- labels(rterms)
  }else{
    Rterms <- NULL
  }

  if (!missing(controls)){
    if (!(controls %in% Rterms) && !(controls %in% Fterms)){
      #add CONTROLS to FIXED model, unless already included in RANDOM
      Fterms <- c(Fterms, controls)
    }
    lterm <- paste(controls, lines, sep=":")
    tterm <- paste(controls, testers, sep=":")
    ltterm<- paste(controls, lines, testers, sep=":")
    ltmodel <- c(lterm, ltterm)
  }else{
    lterm <- lines
    tterm <- testers
    ltterm<- paste(lines, testers, sep=":")
    ltmodel <- c(lterm, ltterm)
  }

  Rterms <- unique(c(Rterms, ltmodel))
  Rterms0 <- Rterms[!(Rterms%in%ltmodel)]

  if (!(tterm %in% Rterms) && !(tterm %in% Fterms)){
    #add TESTERS to FIXED model, unless already included in RANDOM
    Fterms <- c(Fterms, tterm)
  }
  if (engine == "asreml"){
    ok <- require(asreml, quietly = TRUE)
    if (!ok)
      stop("asreml cannot be successfully loaded")
    if (!is.null(Fterms))
      lxtmodel1 <- try(asreml::asreml(fixed=as.formula(paste(resp, paste(Fterms, collapse="+"),sep="~")),
      random=as.formula(paste("~", paste(Rterms, collapse="+"), sep="")), aom=T,
      data=data, na.method.Y = na.method.Y, na.method.X = na.method.X,...), silent = TRUE)
    else
      lxtmodel1 <- try(asreml::asreml(fixed=as.formula(paste(resp, 1, sep="~")),
      random=as.formula(paste("~", paste(Rterms, collapse="+"), sep="")), aom=T,
      data=data, na.method.Y = na.method.Y, na.method.X = na.method.X,...), silent = TRUE)
    flag1 <- inherits(lxtmodel1, "try-error")
    if (!flag1){
      lxtmodel1$call$fixed  <- eval(lxtmodel1$call$fixed)
      lxtmodel1$call$random <- eval(lxtmodel1$call$random)
      lxtmodel1$call$rcov <- eval(lxtmodel1$call$rcov)
      lxtmodel1$call$na.method.X <- na.method.X
      lxtmodel1$call$na.method.Y <- na.method.Y
      if (!lxtmodel1$converge)
        flag1 <- TRUE
    }

    if (flag1){
      message("Model could not be fitted successfully.")
      anysuccess <- FALSE
      if (recover){
        n0 <- length(Rterms0)
        for (ii in 1:n0){
          if (!is.null(Fterms))
            model.cur <- try(asreml::asreml(fixed=as.formula(paste(resp, paste(Fterms, collapse="+"),sep="~")),
            random=as.formula(paste("~", paste(c(Rterms0[1:ii], ltmodel), collapse="+"), sep="")),
            aom=T, data=data, na.method.Y = na.method.Y, na.method.X = na.method.X,...), silent = TRUE)
          else
            model.cur <- try(asreml::asreml(fixed=as.formula(paste(resp, 1, sep="~")),
            random=as.formula(paste("~", paste(c(Rterms0[1:ii], ltmodel), collapse="+"), sep="")),
            aom=T, data=data, na.method.Y = na.method.Y, na.method.X = na.method.X,...), silent = TRUE)

          flag2 <- inherits(model.cur, "try-error")
          if (!flag2){
            model.cur$call$random <- eval(model.cur$call$random)
            model.cur$call$fixed  <- eval(model.cur$call$fixed)
            model.cur$call$na.method.X <- na.method.X
            model.cur$call$na.method.Y <- na.method.Y
            if (!model.cur$converge)
              flag2 <- TRUE
          }

          if (!flag2){
            anysuccess <- TRUE
            if (ii==1){
              lxtmodel2 <- model.cur
            }else{
              if (method == "aic"){
                criterion.cur  <- -2*model.cur$loglik+2*length(model.cur$gammas)
                criterion.prev <- -2*lxtmodel2$loglik+2*length(lxtmodel2$gammas)
              }else{
                criterion.cur  <- -2*model.cur$loglik+log(length(model.cur$fitted.values))*length(model.cur$gammas)
                criterion.prev <- -2*lxtmodel2$loglik+log(length(lxtmodel2$fitted.values))*length(lxtmodel2$gammas)
              }
              if (criterion.cur < criterion.prev)
                lxtmodel2 <- model.cur
            }
          }
        }
        if (!anysuccess)
          message("Model could not be fitted successfully.")
      }
    }

    res <- new.env()

    # BLUPs for lterm & ltterm
    if (!flag1){
      res$lxtmodel <- lxtmodel1
      tmp=tempfile()
      sink(file=tmp)
      blups.line <- predict(lxtmodel1, classify=lterm, only=lterm, data=data)$predictions$pvals
      blups.line.tester <- predict(lxtmodel1, classify=ltterm, only=ltterm, data=data)$predictions$pvals
      lxtmodel00 <- update(lxtmodel1, random=as.formula(paste( "~.",ltmodel[2], sep="-")))
      lxtmodel01 <- update(lxtmodel1, random=as.formula(paste( "~.",paste(ltmodel,collapse="-"), sep="-")))
      sink()
      unlink(tmp)
      deviance2 <- 2*(lxtmodel1$loglik-lxtmodel00$loglik)
      deviance1 <- 2*(lxtmodel00$loglik-lxtmodel01$loglik)
      if (verbose){
        cat("\n")
        print(summary(lxtmodel1))
        cat("\n")
        print(asreml::wald(lxtmodel1))
        cat("\n")
      }
    }else{
      if (anysuccess){
        res$lxtmodel <- lxtmodel2
        tmp=tempfile()
        sink(file=tmp)
        blups.line <- predict(lxtmodel2, classify=lterm, only=lterm, data=data)$predictions$pvals
        blups.line.tester <- predict(lxtmodel2, classify=ltterm, only=ltterm, data=data)$predictions$pvals
        lxtmodel00 <- update(lxtmodel2, random=as.formula(paste( "~.",ltmodel[2], sep="-")))
        lxtmodel01 <- update(lxtmodel2, random=as.formula(paste( "~.",paste(ltmodel,collapse="-"), sep="-")))
        sink()
        unlink(tmp)
        deviance2 <- 2*(lxtmodel2$loglik-lxtmodel00$loglik)
        deviance1 <- 2*(lxtmodel00$loglik-lxtmodel01$loglik)
        if (verbose){
          cat("\n")
          print(summary(lxtmodel2))
          cat("\n")
          print(asreml::wald(lxtmodel2))
          cat("\n")
        }
      }else{
        blups.line <-blups.line.tester <- NULL
        deviance2 <- deviance1 <- NULL
      }
    }
    if (!is.null(deviance2)&&!is.null(deviance1)){
      tlxtmodel <- res$lxtmodel
      G.param.tmp <- tlxtmodel$G.param
      GpNames <- names(G.param.tmp)
      for (ii in 1:length(GpNames)){
        firstName <- unlist(strsplit(GpNames[ii], "\\:"))[1]
        if (!is.null(G.param.tmp[[GpNames[ii]]][[firstName]]$con))
          G.param.tmp[[GpNames[ii]]][[firstName]]$con <- "F"
      }
      fterms <- terms(tlxtmodel$call$fixed)
      fterms <- attr(fterms, "term.labels")
      rterms <- terms(tlxtmodel$call$random)
      rterms <- attr(rterms, "term.labels")

      if (!missing(controls)){
        rterms0 <- rterms[!(rterms %in% c(paste(controls,lines,sep=":"),paste(controls,lines,testers,sep=":")))]
        rform <- as.formula(paste("~",paste(rterms0, collapse="+")))
        fterms <- c(fterms, rterms[rterms %in% c(paste(controls,lines,sep=":"),paste(controls,lines,testers,sep=":"))])
        fform <- as.formula(paste(resp, "~", paste(fterms, collapse="+")))
      }else{
        rterms0 <- rterms[!(rterms %in% c(lines, paste(lines,testers,sep=":")))]
        if (length(rterms0)>0){
          rform <- as.formula(paste("~",paste(rterms0, collapse="+")))
        }else{
          rform <- NULL
        }
        fterms <- c(fterms, rterms[rterms %in% c(lines, paste(lines,testers,sep=":"))])
        fform <- as.formula(paste(resp, "~", paste(fterms, collapse="+")))
      }
      if (!is.null(rform)){
        res$lxtmodel2 <- try(asreml::asreml(fixed=fform, random=rform, G.param = G.param.tmp,
                                            aom=T, data=data, na.method.Y = na.method.Y,
                                            na.method.X = na.method.X,...),
                             silent = TRUE)
      }else{
        res$lxtmodel2 <- try(asreml::asreml(fixed=fform,
                                            aom=T, data=data, na.method.Y = na.method.Y,
                                            na.method.X = na.method.X,...),
                             silent = TRUE)
      }
      if (verbose && !inherits(res$lxtmodel2, "try-error")){
        res$lxtmodel2$call$fixed <- eval(res$lxtmodel2$call$fixed)
        res$lxtmodel2$call$random <- eval(res$lxtmodel2$call$random)
        res$lxtmodel2$call$rcov <- eval(res$lxtmodel2$call$rcov)
        res$lxtmodel2$call$data <- eval(res$lxtmodel2$call$data)
        res$lxtmodel2$call$G.param <- eval(res$lxtmodel2$call$G.param)
        res$lxtmodel2$call$na.method.Y <- eval(res$lxtmodel2$call$na.method.Y)
        res$lxtmodel2$call$na.method.X <- eval(res$lxtmodel2$call$na.method.X)
        cat("\n")
        print(summary(res$lxtmodel2))
        cat("\n")
        print(asreml::wald(res$lxtmodel2))
        cat("\n")
      }
      combTab <- matrix(, nrow=2, ncol=3)
      colnames(combTab) <- c("Chisq", "Chi Df", "Pr(>Chisq)")
      rownames(combTab) <- ltmodel
      combTab[1,"Chisq"] <- deviance1
      combTab[2,"Chisq"] <- deviance2
      combTab[,"Chi Df"] <- 1
      combTab[,"Pr(>Chisq)"] <- 1 - pchisq(c(deviance1, deviance2),1)
      res$test.combinability <- combTab
    }else{
      res$test.combinability <- NULL
    }
    if (verbose && !is.null(blups.line) && !is.null(blups.line.tester)){
      cat("Tests for combinability effects\n\n")
      printCoefmat(combTab, signif.stars = TRUE, P.values = TRUE, has.Pvalue= TRUE, digits = 4)
      if (!missing(controls)){
        # Assume the second level of controls factor specifies test lines;
        # controls factor only contains two levels, (e.g. check & test line).
        testlabl <- levels(data[[controls]])[2]
        cat("\nLines BLUPs (GCA effects)\n\n")
        blups.line2 <- blups.line[blups.line[[controls]]==testlabl, c(lines, "predicted.value", "standard.error")]
        colnames(blups.line2) <- c(lines, "GCA", "s.e.")
        rownames(blups.line2) <- NULL
        print(format(blups.line2, digits = 4), quote = FALSE)

        cat("\nLines:Testers BLUPs (SCA effects)\n\n")
        blups.line.tester2 <- blups.line.tester[blups.line.tester[[controls]]==testlabl, c(lines, testers, "predicted.value", "standard.error")]
        rownames(blups.line.tester2) <- NULL
        blups.line.tester2[[lines]] <- as.factor(blups.line.tester2[[lines]])
        blups.line.tester2[[testers]] <- as.factor(blups.line.tester2[[testers]])
        blups.line.tester.est <- tapply(X=blups.line.tester2[["predicted.value"]], INDEX=blups.line.tester2[,c(lines,testers)], identity)
        colnames(blups.line.tester.est) <- paste("SCA",colnames(blups.line.tester.est), sep=" ")
        blups.line.tester.se  <- tapply(X=blups.line.tester2[["standard.error"]], INDEX=blups.line.tester2[,c(lines,testers)], identity)
        colnames(blups.line.tester.se) <- paste("s.e.",colnames(blups.line.tester.se), sep=" ")
        blups.line.tester3 <- cbind(blups.line.tester.est,blups.line.tester.se)
        nblt <- ncol(blups.line.tester3)
        firsthalf <- 1:(nblt/2)
        secondhalf <- (nblt/2+1):nblt
        oldseq <- newseq <- 1:nblt
        newseq[as.logical(oldseq%%2)] <- firsthalf
        newseq[!as.logical(oldseq%%2)] <- secondhalf
        blups.line.tester3 <- as.data.frame(blups.line.tester3[,newseq])
        blups.line.tester3 <- cbind(rownames(blups.line.tester3),blups.line.tester3)
        rownames(blups.line.tester3) <- NULL
        colnames(blups.line.tester3)[1] <- lines
        print(format(blups.line.tester3, digits = 4), quote = FALSE)
      }else{
        cat("\nLines BLUPs (GCA effects)\n\n")
        blups.line2 <- blups.line[,c(lines,"predicted.value", "standard.error")]
        colnames(blups.line2) <- c(lines, "GCA", "s.e.")
        print(format(blups.line2, digits = 4), quote = FALSE)

        cat("\nLines:Testers BLUPs (SCA effects)\n\n")
        blups.line.tester2 <- blups.line.tester[, c(lines, testers, "predicted.value", "standard.error")]
        rownames(blups.line.tester2) <- NULL
        blups.line.tester2[[lines]] <- as.factor(blups.line.tester2[[lines]])
        blups.line.tester2[[testers]] <- as.factor(blups.line.tester2[[testers]])
        blups.line.tester.est <- tapply(X=blups.line.tester2[["predicted.value"]], INDEX=blups.line.tester2[,c(lines,testers)], identity)
        colnames(blups.line.tester.est) <- paste("SCA",colnames(blups.line.tester.est), sep=" ")
        blups.line.tester.se  <- tapply(X=blups.line.tester2[["standard.error"]], INDEX=blups.line.tester2[,c(lines,testers)], identity)
        colnames(blups.line.tester.se) <- paste("s.e.",colnames(blups.line.tester.se), sep=" ")
        blups.line.tester3 <- cbind(blups.line.tester.est,blups.line.tester.se)
        nblt <- ncol(blups.line.tester3)
        firsthalf <- 1:(nblt/2)
        secondhalf <- (nblt/2+1):nblt
        oldseq <- newseq <- 1:nblt
        newseq[as.logical(oldseq%%2)] <- firsthalf
        newseq[!as.logical(oldseq%%2)] <- secondhalf
        blups.line.tester3 <- as.data.frame(blups.line.tester3[,newseq])
        blups.line.tester3 <- cbind(rownames(blups.line.tester3),blups.line.tester3)
        rownames(blups.line.tester3) <- NULL
        colnames(blups.line.tester3)[1] <- lines
        print(format(blups.line.tester3, digits = 4), quote = FALSE)
      }
    }
  }

  if (engine=="lme4"){
    ok <- require(lme4, quietly = TRUE)
    if (!ok)
      stop("lme4 cannot be successfully loaded")
    Rpart <- paste("(1 | ",Rterms, ")")
    ltpart <- paste("(1 | ",ltmodel, ")")
    if (!is.null(Fterms))
      lxtmodel1 <- try(lme4::lmer(as.formula(paste(resp, paste(c(Fterms, Rpart), collapse="+"),sep="~")),
      data=data, na.action=na.action, ...), silent = TRUE)
    else
      lxtmodel1 <- try(lme4::lmer(as.formula(paste(resp, paste(Rpart, collapse="+"), sep="~")),
      data=data, na.action=na.action,...), silent = TRUE)

    flag1 <- inherits(lxtmodel1, "try-error")
    n0 <- length(Rterms0)
    if (n0>0)
      Rpart0 <- paste("(1 | ",Rterms0, ")")
    anysuccess <- FALSE
    if (flag1 && n0>0){
      message("Model could not be fitted successfully.")
      if (recover){
        for (ii in 1:n0){
          if (!is.null(Fterms))
            model.cur <- try(lme4::lmer(as.formula(paste(resp, paste(c(Fterms,Rpart0[1:ii],ltpart), collapse="+"),sep="~")),
            data=data, na.action=na.action,...), silent = TRUE)
          else
            model.cur <- try(lme4::lmer(as.formula(paste(resp, paste(c(Rpart0[1:ii],ltpart), collapse="+"), sep="~")),
            data=data, na.action=na.action,...), silent = TRUE)

          flag2 <- inherits(model.cur, "try-error")
          if (!flag2){
            anysuccess <- TRUE
            if (ii==1){
              lxtmodel2 <- model.cur
              index <- 1
            }else{
              if (method == "aic"){
                criterion.cur  <- AIC(logLik(model.cur))
                criterion.prev <- AIC(logLik(lxtmodel2))
              }else{
                criterion.cur  <- BIC(logLik(model.cur))
                criterion.prev <- BIC(logLik(lxtmodel2))
              }
              if (criterion.cur < criterion.prev){
                lxtmodel2 <- model.cur
                index <- ii
              }
            }
          }
        }
        if (!anysuccess)
          message("Model could not be fitted successfully.")
      }
    }

    res <- new.env()

    # BLUPs for lterm & ltterm
    if (!flag1){
      res$lxtmodel <- lxtmodel1
      blups.line <- lme4::ranef(lxtmodel1)[[lterm]]
      blups.line.tester <- lme4::ranef(lxtmodel1)[[ltterm]]
      if (n0>0){
        if (!is.null(Fterms))
          fform <- as.formula(paste(resp, paste(c(Fterms,Rpart0,ltmodel), collapse="+"),sep="~"))
        else
          fform <- as.formula(paste(resp, paste(c(Rpart0,ltmodel), collapse="+"),sep="~"))
      }else{
        if (!is.null(Fterms))
          fform <- as.formula(paste(resp, paste(c(Fterms,ltmodel), collapse="+"),sep="~"))
        else
          fform <- as.formula(paste(resp, paste(ltmodel, collapse="+"),sep="~"))
      }
      if (verbose){
        cat("\n")
        print(summary(lxtmodel1))
        cat("\n")
      }
    }else{
      if (anysuccess){
        res$lxtmodel <- lxtmodel2
        blups.line <- lme4::ranef(lxtmodel2)[[lterm]]
        blups.line.tester <- lme4::ranef(lxtmodel2)[[ltterm]]
        if (!is.null(Fterms))
          fform <- as.formula(paste(resp, paste(c(Fterms,Rpart0[1:index],ltmodel), collapse="+"),sep="~"))
        else
          fform <- as.formula(paste(resp, paste(c(Rpart0[1:index],ltmodel), collapse="+"),sep="~"))
        if (verbose){
          cat("\n")
          print(summary(lxtmodel2))
          cat("\n")
        }
      }else{
        blups.line <- blups.line.tester <- NULL
      }
    }
    if (!is.null(blups.line) && !is.null(blups.line.tester)){
      if (n0>0){
        res$lxtmodel2 <- try(lme4::lmer(fform, data=data, na.action=na.action,...), silent = TRUE)
      }else{
        res$lxtmodel2 <- try(lm(fform, data=data, na.action=na.action,...), silent = TRUE)
      }
      if (verbose && !inherits(res$lxtmodel2, "try-error")){
        cat("\n")
        print(summary(res$lxtmodel2))
        cat("\n")
      }
      if (n0>0){
        lxtmodel00 <- update(res$lxtmodel, as.formula(paste( ".~.",ltpart[2], sep="-")), na.action= na.action)
        lxtmodel01 <- update(res$lxtmodel, as.formula(paste( ".~.",paste(ltpart,collapse="-"), sep="-")), na.action= na.action)
        res$test.combinability <- combTab <- anova(lxtmodel01, lxtmodel00, res$lxtmodel, refit=refit)
        rownames(combTab) <-  c("", lterm, ltterm)
        combTab <-  combTab[,c("Chisq", "Chi Df", "Pr(>Chisq)")]
        combTab <- combTab[-1,]
      }
    }else{
      res$test.combinability <- NULL
    }
    if (verbose && !is.null(blups.line) && !is.null(blups.line.tester) && n0>0){
      cat("Tests for combinability effects\n\n")
      printCoefmat(combTab, signif.stars = TRUE, P.values = TRUE, has.Pvalue= TRUE, digits = 4)
      if (!missing(controls)){
        # Assume the second level of controls factor specifies test lines;
        # controls factor only contains two levels, (e.g. check & test line).
        testlabl <- levels(data[[controls]])[2]
        cat("\nLines BLUPs (GCA effects)\n\n")
        rownames.bl <- rownames(blups.line)
        nameM.bl <- sapply(rownames.bl, function(x) unlist(strsplit(x, "\\:")))
        ind.incl <- nameM.bl[1,]==testlabl
        blups.line <- blups.line[ind.incl,, drop=F]
        nameM.bl <- nameM.bl[,ind.incl]
        blups.line2 <- cbind(nameM.bl[2,],blups.line)
        colnames(blups.line2) <- c(lines, "GCA")
        rownames(blups.line2) <- NULL
        print(format(blups.line2, digits = 4), quote = FALSE)

        cat("\nLines:Testers BLUPs (SCA effects)\n\n")
        rownames.blt <- rownames(blups.line.tester)
        nameM.blt <- sapply(rownames.blt, function(x) unlist(strsplit(x, "\\:")))
        ind.incl <- nameM.blt[1,]==testlabl
        blups.line.tester <- blups.line.tester[ind.incl,, drop=F]
        names(blups.line.tester) <- "SCA"
        nameM.blt <- nameM.blt[,ind.incl]
        lineFactor <- as.factor(nameM.blt[2,])
        testerFactor <- as.factor(nameM.blt[3,])
        blups.line.tester[[lines]] <-lineFactor
        blups.line.tester[[testers]] <-testerFactor
        blups.line.tester2 <- tapply(blups.line.tester[["SCA"]], blups.line.tester[,c(lines,testers)], identity)
        print(format(blups.line.tester2, digits = 4), quote = FALSE)
      }else{
        cat("\nLines BLUPs (GCA effects)\n\n")
        blups.line2 <- cbind(rownames(blups.line), blups.line)
        colnames(blups.line2) <- c(lines, "GCA")
        rownames(blups.line2) <- NULL
        print(format(blups.line2, digits = 4), quote = FALSE)

        cat("\nLines:Testers BLUPs (SCA effects)\n\n")
        rownames.blt <- rownames(blups.line.tester)
        nameM.blt <- sapply(rownames.blt, function(x) unlist(strsplit(x, "\\:")))
        names(blups.line.tester) <- "SCA"
        lineFactor <- as.factor(nameM.blt[1,])
        testerFactor <- as.factor(nameM.blt[2,])
        blups.line.tester[[lines]] <-lineFactor
        blups.line.tester[[testers]] <-testerFactor
        blups.line.tester2 <- tapply(blups.line.tester[["SCA"]], blups.line.tester[,c(lines,testers)], identity)
        print(format(blups.line.tester2, digits = 4), quote = FALSE)
      }
    }
  }
  res$blups.line <- blups.line
  res$blups.line.tester <- blups.line.tester
  as.list(res)
}
