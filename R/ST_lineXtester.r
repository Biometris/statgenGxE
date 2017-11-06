#' Line x Tester Analysis
#'
#' This function does a mixed-model analysis of data from a line-by-tester trial,
#' using either \code{asreml} or \code{lme4}.
#' @param fixed a formula specifying fixed model terms, in addition to the \code{testers}
#' main effect and any control comparisons.
#' @param random a formula specifying random model terms, in addition to the terms involving
#' \code{lines} and \code{testers} that are included automatically.
#' @param lines a character string specifying the line (usually female parent) in \code{data}.
#' @param testers a character string specifying the tester (usually male parent) in \code{data}.
#' @param controls a character string specifying a factor in \code{data}, which distinguishes
#' between control and test (line x tester) genotypes.
#' @param data a data frame object containing data of a line-by-tester trial.
#' @param method a character string specifying the criterion (either \code{"aic"} or \code{"bic"})
#' to look for the best model if \code{asreml} or \code{lme4} cannot fit the current specifed
#' model; default, \code{"bic"}.
#' @param engine a string specifies the name of a mixed modelling engine (either \code{"asreml"}
#'  or \code{"lme4"}); default, \code{"asreml"}.
#' @param naMethodY a character string (\code{"include"}, \code{"omit"} or \code{"fail"})
#' specifying how missing data in the response is to be handled.
#' This is applied to the model frame after any subset argument has been used. The
#' default (\code{"include"}) is to estimate missing values;
#' this is necessary in spatial models to preserve the spatial structure. \code{"omit"} deletes
#' observations that contain missing values in the response.
#' Note this is only available when \code{asreml} is used.
#' @param naMethodX a character string (\code{"include"}, \code{"omit"} or \code{"fail"})
#' specifying how missing data in covariates is to be handled.
#' This is applied to the model frame after any subset argument has been used but before
#' any \code{at()} functions are invoked. The default \code{"fail"} causes an error if there are
#' missing values in any covariate.
#' \code{"omit"} deletes observations that contain missing values in any covariate.
#' Note this is only available when \code{asreml} is used.
#' @param na.action a function that indicates what should happen when the data contain NAs.
#' The default action (\code{na.omit}) strips any observations with any missing values in any
#' variables. This is only available when \code{lme4} is used as a mixed modelling engine.
#' @param recover logical. Whether to try to recover with a simpler random model if REML
#' cannot fit the model; default \code{FALSE}.
#' @param refit logical indicating if objects of class \code{lmerMod} should be refitted with
#' ML before comparing models. The default is \code{TRUE} to prevent the common mistake of
#' inappropriately comparing REML-fitted models with different fixed effects, whose likelihoods
#' are not directly comparable. This option is available for \code{lme4} only.
#' @param verbose logical. Whether to display a summary of the line-by-tester analysis on
#' screen; Default, \code{TRUE}.
#' @param ... other parameters to be passed on either \code{asreml} or \code{lme4}.
#' @return
#' a list consists of the fitted model object (\code{lxtmodel}), a data frame object containing
#' the GCA effects (\code{blupsLine}), a data frame object containing SCA effects
#' (\code{blupsLineTester}) and a data frame object containing tests for combinability
#' effects (\code{testCombinability}).
#' @examples
#' mydat <- ST.read.csv(system.file("extdata", "VLIN-1.csv", package = "RAP"),
#'                      factorNames=c("Replicates","Blocks","Controls","Lines","Testers"),
#'                      traitNames="yield", env ="Env")
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
ST.lineXtester <- function(fixed,
                           random,
                           lines,
                           testers,
                           controls,
                           data,
                           method = "bic",
                           engine = "asreml",
                           naMethodY = "include",
                           naMethodX = "fail",
                           na.action = na.omit,
                           recover = FALSE,
                           refit = TRUE,
                           verbose = TRUE,
                           ...) {
  if (missing(fixed)) {
    stop("a fixed part formula needed.")
  } else {
    fTerms <- terms(fixed, keep.order = TRUE)
    if (length(attr(fTerms, "factors")) == 0L) {
      resp <- as.character(fTerms)[2]
      if (attr(fTerms, "intercept") == 0L) {
        fTerms <- "-1"
      } else {
        fTerms <- NULL
      }
    } else {
      respLoc <- attr(fTerms, "response")
      resp <- rownames(attr(fTerms, "factors"))[respLoc]
      tFTerms <- fTerms <- colnames(attr(fTerms, "factors"))
      tFTerms <- unlist(strsplit(tFTerms, "\\:"))
      tFTerms <- gsub(pattern = "^[[:space:]]", replacement = "", x = tFTerms)
      tFTerms <- gsub(pattern = "\\s$", replacement = "", x = tFTerms)
      tFTerms <- unique(tFTerms)
      if (!all(tFTerms %in% names(data))) {
        stop(paste(tFTerms[!(tFTerms %in% names(data))], collapse=","),
             " not found in column names of data")
      }
      if (attr(fTerms, "intercept") == 0L) {
        fTerms <- c("-1", fTerms)
      }
    }
    if (!(resp %in% names(data))) {
      stop(resp, " not found in column names of data")
    }
  }
  if (!missing(random)){
    rTerms <- terms(random, keep.order = TRUE)
    RTerms <- labels(rTerms)
  } else{
    RTerms <- NULL
  }
  if (!missing(controls)) {
    if (!(controls %in% RTerms) && !(controls %in% fTerms)) {
      #add CONTROLS to FIXED model, unless already included in RANDOM
      fTerms <- c(fTerms, controls)
    }
    lTerm <- paste(controls, lines, sep=":")
    tTerm <- paste(controls, testers, sep=":")
    ltTerm<- paste(controls, lines, testers, sep=":")
    ltModel <- c(lTerm, ltTerm)
  } else{
    lTerm <- lines
    tTerm <- testers
    ltTerm <- paste(lines, testers, sep=":")
    ltModel <- c(lTerm, ltTerm)
  }
  RTerms <- unique(c(RTerms, ltModel))
  RTerms0 <- RTerms[!(RTerms %in% ltModel)]
  if (!(tTerm %in% RTerms) && !(tTerm %in% fTerms)) {
    #add TESTERS to FIXED model, unless already included in RANDOM
    fTerms <- c(fTerms, tTerm)
  }
  if (engine == "asreml") {
    if (!requireNamespace("asreml", quietly = TRUE)) {
      stop("asreml cannot be successfully loaded")
    }
    if (!is.null(fTerms)) {
      lxtModel1 <- try(asreml::asreml(fixed = as.formula(paste(resp,
                                                               paste(fTerms, collapse = "+"),
                                                               sep = "~")),
                                      random = as.formula(paste("~", paste(RTerms, collapse = "+"),
                                                                sep="")),
                                      aom = TRUE, data = data, na.method.Y = naMethodY,
                                      na.method.X = naMethodX, ...),
                       silent = TRUE)
    } else {
      lxtModel1 <- try(asreml::asreml(fixed = as.formula(paste(resp, 1, sep = "~")),
                                      random = as.formula(paste("~", paste(RTerms, collapse = "+"),
                                                                sep = "")),
                                      aom = TRUE, data = data, na.method.Y = naMethodY,
                                      na.method.X = naMethodX, ...),
                       silent = TRUE)
    }
    flag1 <- inherits(lxtModel1, "try-error")
    if (!flag1) {
      lxtModel1$call$fixed  <- eval(lxtModel1$call$fixed)
      lxtModel1$call$random <- eval(lxtModel1$call$random)
      lxtModel1$call$rcov <- eval(lxtModel1$call$rcov)
      lxtModel1$call$na.method.X <- naMethodX
      lxtModel1$call$na.method.Y <- naMethodY
      if (!lxtModel1$converge) {
        flag1 <- TRUE
      }
    }
    if (flag1) {
      message("Model could not be fitted successfully.\n")
      anySuccess <- FALSE
      if (recover) {
        n0 <- length(RTerms0)
        for (ii in 1:n0) {
          if (!is.null(fTerms)) {
            modelCur <- try(asreml::asreml(fixed = as.formula(paste(resp,
                                                                    paste(fTerms, collapse = "+"),
                                                                    sep = "~")),
                                           random = as.formula(paste("~",
                                                                     paste(c(RTerms0[1:ii], ltModel),
                                                                           collapse = "+"),
                                                                     sep = "")),
                                           aom = TRUE, data = data, na.method.Y = naMethodY,
                                           na.method.X = naMethodX, ...),
                            silent = TRUE)
          } else {
            modelCur <- try(asreml::asreml(fixed = as.formula(paste(resp, 1, sep = "~")),
                                           random = as.formula(paste("~",
                                                                     paste(c(RTerms0[1:ii], ltModel),
                                                                           collapse = "+"),
                                                                     sep = "")),
                                           aom = TRUE, data = data, na.method.Y = naMethodY,
                                           na.method.X = naMethodX, ...),
                            silent = TRUE)
          }
          flag2 <- inherits(modelCur, "try-error")
          if (!flag2) {
            modelCur$call$random <- eval(modelCur$call$random)
            modelCur$call$fixed <- eval(modelCur$call$fixed)
            modelCur$call$naMethodX <- naMethodX
            modelCur$call$naMethodY <- naMethodY
            if (!modelCur$converge) {
              flag2 <- TRUE
            }
          }
          if (!flag2) {
            anySuccess <- TRUE
            if (ii == 1) {
              lxtModel2 <- modelCur
            } else {
              if (method == "aic") {
                criterionCur <- -2 * modelCur$loglik + 2 * length(modelCur$gammas)
                criterionPrev <- -2 * lxtModel2$loglik + 2 * length(lxtModel2$gammas)
              } else {
                criterionCur <- -2 * modelCur$loglik +
                  log(length(modelCur$fitted.values)) * length(modelCur$gammas)
                criterionPrev <- -2 * lxtModel2$loglik +
                  log(length(lxtModel2$fitted.values)) * length(lxtModel2$gammas)
              }
              if (criterionCur < criterionPrev) {
                lxtModel2 <- modelCur
              }
            }
          }
        }
        if (!anySuccess) {
          message("Model could not be fitted successfully.\n")
        }
      }
    }
    res <- new.env()
    # BLUPs for lTerm & ltTerm
    if (!flag1) {
      res$lxtmodel <- lxtModel1
      tmp <- tempfile()
      sink(file = tmp)
      blupsLine <- predict(lxtModel1, classify = lTerm, only = lTerm,
                           data = data)$predictions$pvals
      blupsLineTester <- predict(lxtModel1, classify = ltTerm, only = ltTerm,
                                 data=data)$predictions$pvals
      lxtModel00 <- update(lxtModel1,
                           random = as.formula(paste("~.", ltModel[2], sep = "-")))
      lxtModel01 <- update(lxtModel1,
                           random = as.formula(paste("~.", paste(ltModel, collapse = "-"), sep = "-")))
      sink()
      unlink(tmp)
      deviance2 <- 2 * (lxtModel1$loglik - lxtModel00$loglik)
      deviance1 <- 2 * (lxtModel00$loglik - lxtModel01$loglik)
      if (verbose){
        cat("\n")
        print(summary(lxtModel1))
        cat("\n")
        print(asreml::wald.asreml(lxtModel1))
        cat("\n")
      }
    } else if (anySuccess) {
      res$lxtmodel <- lxtModel2
      tmp <- tempfile()
      sink(file = tmp)
      blupsLine <- predict(lxtModel2, classify = lTerm, only = lTerm,
                           data = data)$predictions$pvals
      blupsLineTester <- predict(lxtModel2, classify = ltTerm, only = ltTerm,
                                 data = data)$predictions$pvals
      lxtModel00 <- update(lxtModel2,
                           random = as.formula(paste("~.", ltModel[2], sep = "-")))
      lxtModel01 <- update(lxtModel2,
                           random = as.formula(paste("~.", paste(ltModel, collapse = "-"), sep = "-")))
      sink()
      unlink(tmp)
      deviance2 <- 2 * (lxtModel2$loglik - lxtModel00$loglik)
      deviance1 <- 2 * (lxtModel00$loglik - lxtModel01$loglik)
      if (verbose) {
        cat("\n")
        print(summary(lxtModel2))
        cat("\n")
        print(asreml::wald(lxtModel2))
        cat("\n")
      }
    } else {
      blupsLine <- blupsLineTester <- NULL
      deviance2 <- deviance1 <- NULL
    }
    if (!is.null(deviance2) && !is.null(deviance1)) {
      tlxtModel <- res$lxtmodel
      GParamTmp <- tlxtModel$G.param
      GpNames <- names(GParamTmp)
      for (ii in 1:length(GpNames)) {
        firstName <- unlist(strsplit(GpNames[ii], "\\:"))[1]
        if (!is.null(GParamTmp[[GpNames[ii]]][[firstName]]$con)) {
          GParamTmp[[GpNames[ii]]][[firstName]]$con <- "F"
        }
      }
      fTerms <- terms(tlxtModel$call$fixed)
      fTerms <- attr(fTerms, "term.labels")
      rTerms <- terms(tlxtModel$call$random)
      rTerms <- attr(rTerms, "term.labels")
      if (!missing(controls)) {
        rTerms0 <- rTerms[!(rTerms %in% c(paste(controls, lines, sep = ":"),
                                          paste(controls, lines, testers, sep = ":")))]
        rForm <- as.formula(paste("~", paste(rTerms0, collapse = "+")))
        fTerms <- c(fTerms, rTerms[rTerms %in% c(paste(controls, lines, sep = ":"),
                                                 paste(controls, lines, testers, sep = ":"))])
        fForm <- as.formula(paste(resp, "~", paste(fTerms, collapse = "+")))
      } else {
        rTerms0 <- rTerms[!(rTerms %in% c(lines, paste(lines, testers, sep = ":")))]
        if (length(rTerms0) > 0) {
          rForm <- as.formula(paste("~", paste(rTerms0, collapse = "+")))
        } else {
          rForm <- NULL
        }
        fTerms <- c(fTerms, rTerms[rTerms %in% c(lines, paste(lines, testers, sep = ":"))])
        fForm <- as.formula(paste(resp, "~", paste(fTerms, collapse = "+")))
      }
      if (!is.null(rForm)) {
        res$lxtModel2 <- try(asreml::asreml(fixed = fForm, random = rForm, G.param = GParamTmp,
                                            aom = TRUE, data = data, na.method.Y = naMethodY,
                                            na.method.X = naMethodX, ...),
                             silent = TRUE)
      } else {
        res$lxtModel2 <- try(asreml::asreml(fixed = fForm, aom = TRUE, data = data,
                                            na.method.Y = naMethodY, na.method.X = naMethodX, ...),
                             silent = TRUE)
      }
      if (verbose && !inherits(res$lxtModel2, "try-error")) {
        res$lxtModel2$call$fixed <- eval(res$lxtModel2$call$fixed)
        res$lxtModel2$call$random <- eval(res$lxtModel2$call$random)
        res$lxtModel2$call$rcov <- eval(res$lxtModel2$call$rcov)
        res$lxtModel2$call$data <- eval(res$lxtModel2$call$data)
        res$lxtModel2$call$G.param <- eval(res$lxtModel2$call$G.param)
        res$lxtModel2$call$na.method.Y <- eval(res$lxtModel2$call$naMethodY)
        res$lxtModel2$call$na.method.X <- eval(res$lxtModel2$call$naMethodX)
        cat("\n")
        print(summary(res$lxtModel2))
        cat("\n")
        print(asreml::wald.asreml(res$lxtModel2))
        cat("\n")
      }
      combTab <- matrix(nrow = 2, ncol = 3)
      colnames(combTab) <- c("Chisq", "Chi Df", "Pr(>Chisq)")
      rownames(combTab) <- ltModel
      combTab[1, "Chisq"] <- deviance1
      combTab[2, "Chisq"] <- deviance2
      combTab[, "Chi Df"] <- 1
      combTab[, "Pr(>Chisq)"] <- 1 - pchisq(q = c(deviance1, deviance2), df = 1)
      res$testCombinability <- combTab
    } else {
      res$testCombinability <- NULL
    }
    if (verbose && !is.null(blupsLine) && !is.null(blupsLineTester)) {
      cat("Tests for combinability effects\n\n")
      printCoefmat(combTab, signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE, digits = 4)
      if (!missing(controls)) {
        # Assume the second level of controls factor specifies test lines;
        # controls factor only contains two levels, (e.g. check & test line).
        testLabl <- levels(data[[controls]])[2]
        cat("\nLines BLUPs (GCA effects)\n\n")
        blupsLine2 <- blupsLine[blupsLine[[controls]] == testLabl,
                                c(lines, "predicted.value", "standard.error")]
        colnames(blupsLine2) <- c(lines, "GCA", "s.e.")
        rownames(blupsLine2) <- NULL
        print(format(blupsLine2, digits = 4), quote = FALSE)
        cat("\nLines:Testers BLUPs (SCA effects)\n\n")
        blupsLineTester2 <- blupsLineTester[blupsLineTester[[controls]] == testLabl,
                                            c(lines, testers, "predicted.value", "standard.error")]
        rownames(blupsLineTester2) <- NULL
        blupsLineTester2[[lines]] <- as.factor(blupsLineTester2[[lines]])
        blupsLineTester2[[testers]] <- as.factor(blupsLineTester2[[testers]])
        blupsLineTesterEst <- tapply(X = blupsLineTester2[["predicted.value"]],
                                     INDEX = blupsLineTester2[, c(lines, testers)],
                                     FUN = identity)
        colnames(blupsLineTesterEst) <- paste("SCA", colnames(blupsLineTesterEst))
        blupsLineTesterSe <- tapply(X = blupsLineTester2[["standard.error"]],
                                    INDEX = blupsLineTester2[, c(lines, testers)],
                                    FUN = identity)
        colnames(blupsLineTesterSe) <- paste("s.e.", colnames(blupsLineTesterSe))
        blupsLineTester3 <- cbind(blupsLineTesterEst, blupsLineTesterSe)
        nBlt <- ncol(blupsLineTester3)
        firstHalf <- 1:(nBlt / 2)
        secondHalf <- (nBlt / 2 + 1):nBlt
        oldSeq <- newSeq <- 1:nBlt
        newSeq[as.logical(oldSeq %% 2)] <- firstHalf
        newSeq[!as.logical(oldSeq %% 2)] <- secondHalf
        blupsLineTester3 <- as.data.frame(blupsLineTester3[, newSeq])
        blupsLineTester3 <- cbind(rownames(blupsLineTester3), blupsLineTester3)
        rownames(blupsLineTester3) <- NULL
        colnames(blupsLineTester3)[1] <- lines
        print(format(blupsLineTester3, digits = 4), quote = FALSE)
      } else {
        cat("\nLines BLUPs (GCA effects)\n\n")
        blupsLine2 <- blupsLine[, c(lines, "predicted.value", "standard.error")]
        colnames(blupsLine2) <- c(lines, "GCA", "s.e.")
        print(format(blupsLine2, digits = 4), quote = FALSE)
        cat("\nLines:Testers BLUPs (SCA effects)\n\n")
        blupsLineTester2 <- blupsLineTester[, c(lines, testers, "predicted.value", "standard.error")]
        rownames(blupsLineTester2) <- NULL
        blupsLineTester2[[lines]] <- as.factor(blupsLineTester2[[lines]])
        blupsLineTester2[[testers]] <- as.factor(blupsLineTester2[[testers]])
        blupsLineTesterEst <- tapply(X = blupsLineTester2[["predicted.value"]],
                                     INDEX = blupsLineTester2[, c(lines, testers)],
                                     FUN = identity)
        colnames(blupsLineTesterEst) <- paste("SCA",colnames(blupsLineTesterEst))
        blupsLineTesterSe <- tapply(X = blupsLineTester2[["standard.error"]],
                                    INDEX = blupsLineTester2[, c(lines, testers)],
                                    FUN = identity)
        colnames(blupsLineTesterSe) <- paste("s.e.", colnames(blupsLineTesterSe))
        blupsLineTester3 <- cbind(blupsLineTesterEst, blupsLineTesterSe)
        nBlt <- ncol(blupsLineTester3)
        firstHalf <- 1:(nBlt / 2)
        secondHalf <- (nBlt / 2 + 1):nBlt
        oldSeq <- newSeq <- 1:nBlt
        newSeq[as.logical(oldSeq %%2)] <- firstHalf
        newSeq[!as.logical(oldSeq %% 2)] <- secondHalf
        blupsLineTester3 <- as.data.frame(blupsLineTester3[, newSeq])
        blupsLineTester3 <- cbind(rownames(blupsLineTester3), blupsLineTester3)
        rownames(blupsLineTester3) <- NULL
        colnames(blupsLineTester3)[1] <- lines
        print(format(blupsLineTester3, digits = 4), quote = FALSE)
      }
    }
  }
  if (engine == "lme4"){
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("lme4 cannot be successfully loaded")
    }
    rPart <- paste("(1 | ", RTerms, ")")
    ltPart <- paste("(1 | ", ltModel, ")")
    if (!is.null(fTerms)) {
      lxtModel1 <- try(lme4::lmer(as.formula(paste(resp,
                                                   paste(c(fTerms, rPart), collapse = "+"),
                                                   sep = "~")),
                                  data = data, na.action = na.action, ...),
                       silent = TRUE)
    } else {
      lxtModel1 <- try(lme4::lmer(as.formula(paste(resp, paste(rPart, collapse = "+"), sep = "~")),
                                  data = data, na.action = na.action, ...), silent = TRUE)
    }
    flag1 <- inherits(lxtModel1, "try-error")
    n0 <- length(RTerms0)
    if (n0>0) {
      Rpart0 <- paste("(1 | ",RTerms0, ")")
    }
    anySuccess <- FALSE
    if (flag1 && n0 > 0) {
      message("Model could not be fitted successfully.\n")
      if (recover) {
        for (ii in 1:n0) {
          if (!is.null(fTerms)) {
            modelCur <- try(lme4::lmer(as.formula(paste(resp, paste(c(fTerms, Rpart0[1:ii], ltPart),
                                                                    collapse = "+"), sep = "~")),
                                       data = data, na.action = na.action, ...), silent = TRUE)
          } else {
            modelCur <- try(lme4::lmer(as.formula(paste(resp, paste(c(Rpart0[1:ii], ltPart),
                                                                    collapse = "+"), sep = "~")),
                                       data = data, na.action = na.action,...), silent = TRUE)
          }
          flag2 <- inherits(modelCur, "try-error")
          if (!flag2) {
            anySuccess <- TRUE
            if (ii == 1){
              lxtModel2 <- modelCur
              index <- 1
            } else {
              if (method == "aic") {
                criterionCur <- AIC(logLik(modelCur))
                criterionPrev <- AIC(logLik(lxtModel2))
              } else {
                criterionCur  <- BIC(logLik(modelCur))
                criterionPrev <- BIC(logLik(lxtModel2))
              }
              if (criterionCur < criterionPrev) {
                lxtModel2 <- modelCur
                index <- ii
              }
            }
          }
        }
        if (!anySuccess) {
          message("Model could not be fitted successfully.\n")
        }
      }
    }
    res <- new.env()
    # BLUPs for lTerm & ltTerm
    if (!flag1) {
      res$lxtmodel <- lxtModel1
      blupsLine <- lme4::ranef(lxtModel1)[[lTerm]]
      blupsLineTester <- lme4::ranef(lxtModel1)[[ltTerm]]
      if (n0 > 0) {
        if (!is.null(fTerms)) {
          fForm <- as.formula(paste(resp, paste(c(fTerms, Rpart0, ltModel),
                                                collapse = "+"), sep = "~"))
        } else {
          fForm <- as.formula(paste(resp, paste(c(Rpart0, ltModel),
                                                collapse = "+"), sep = "~"))
        }
      } else {
        if (!is.null(fTerms)) {
          fForm <- as.formula(paste(resp, paste(c(fTerms, ltModel),
                                                collapse = "+"), sep = "~"))
        } else {
          fForm <- as.formula(paste(resp, paste(ltModel, collapse = "+"), sep = "~"))
        }
      }
      if (verbose) {
        cat("\n")
        print(summary(lxtModel1))
        cat("\n")
      }
    } else {
      if (anySuccess) {
        res$lxtmodel <- lxtModel2
        blupsLine <- lme4::ranef(lxtModel2)[[lTerm]]
        blupsLineTester <- lme4::ranef(lxtModel2)[[ltTerm]]
        if (!is.null(fTerms)) {
          fForm <- as.formula(paste(resp, paste(c(fTerms, Rpart0[1:index], ltModel),
                                                collapse = "+"), sep = "~"))
        } else {
          fForm <- as.formula(paste(resp, paste(c(Rpart0[1:index], ltModel),
                                                collapse = "+"), sep = "~"))
        }
        if (verbose) {
          cat("\n")
          print(summary(lxtModel2))
          cat("\n")
        }
      } else {
        blupsLine <- blupsLineTester <- NULL
      }
    }
    if (!is.null(blupsLine) && !is.null(blupsLineTester)) {
      if (n0 > 0) {
        res$lxtModel2 <- try(lme4::lmer(fForm, data = data, na.action = na.action,...),
                             silent = TRUE)
      } else {
        res$lxtModel2 <- try(lm(fForm, data = data, na.action = na.action,...),
                             silent = TRUE)
      }
      if (verbose && !inherits(res$lxtModel2, "try-error")) {
        cat("\n")
        print(summary(res$lxtModel2))
        cat("\n")
      }
      if (n0 > 0) {
        lxtModel00 <- update(res$lxtmodel,
                             as.formula(paste(".~.", ltPart[2], sep = "-")), na.action = na.action)
        lxtModel01 <- update(res$lxtmodel,
                             as.formula(paste(".~.", paste(ltPart, collapse = "-"), sep = "-")),
                             na.action = na.action)
        res$testCombinability <- combTab <- anova(lxtModel01, lxtModel00,
                                                  res$lxtmodel, refit = refit)
        rownames(combTab) <- c("", lTerm, ltTerm)
        combTab <- combTab[, c("Chisq", "Chi Df", "Pr(>Chisq)")]
        combTab <- combTab[-1, ]
      }
    } else{
      res$testCombinability <- NULL
    }
    if (verbose && !is.null(blupsLine) && !is.null(blupsLineTester) && n0 > 0) {
      cat("Tests for combinability effects\n\n")
      printCoefmat(combTab, signif.stars = TRUE, P.values = TRUE, has.Pvalue = TRUE, digits = 4)
      if (!missing(controls)) {
        # Assume the second level of controls factor specifies test lines;
        # controls factor only contains two levels, (e.g. check & test line).
        testLabl <- levels(data[[controls]])[2]
        cat("\nLines BLUPs (GCA effects)\n\n")
        rowNamesBl <- rownames(blupsLine)
        nameMBl <- sapply(X = rowNamesBl, FUN = function(x) {
          unlist(strsplit(x, "\\:"))
        })
        indIncl <- nameMBl[1, ] == testLabl
        blupsLine <- blupsLine[indIncl, , drop = FALSE]
        nameMBl <- nameMBl[, indIncl]
        blupsLine2 <- cbind(nameMBl[2, ], blupsLine)
        colnames(blupsLine2) <- c(lines, "GCA")
        rownames(blupsLine2) <- NULL
        print(format(blupsLine2, digits = 4), quote = FALSE)
        cat("\nLines:Testers BLUPs (SCA effects)\n\n")
        rowNamesBlt <- rownames(blupsLineTester)
        nameMBlt <- sapply(X = rowNamesBlt, FUN = function(x) {
          unlist(strsplit(x, "\\:"))
        })
        indIncl <- nameMBlt[1, ] == testLabl
        blupsLineTester <- blupsLineTester[indIncl, , drop = FALSE]
        names(blupsLineTester) <- "SCA"
        nameMBlt <- nameMBlt[, indIncl]
        lineFactor <- as.factor(nameMBlt[2, ])
        testerFactor <- as.factor(nameMBlt[3, ])
        blupsLineTester[[lines]] <- lineFactor
        blupsLineTester[[testers]] <- testerFactor
        blupsLineTester2 <- tapply(X = blupsLineTester[["SCA"]],
                                   INDEX = blupsLineTester[, c(lines, testers)],
                                   FUN = identity)
        print(format(blupsLineTester2, digits = 4), quote = FALSE)
      } else {
        cat("\nLines BLUPs (GCA effects)\n\n")
        blupsLine2 <- cbind(rownames(blupsLine), blupsLine)
        colnames(blupsLine2) <- c(lines, "GCA")
        rownames(blupsLine2) <- NULL
        print(format(blupsLine2, digits = 4), quote = FALSE)
        cat("\nLines:Testers BLUPs (SCA effects)\n\n")
        rowNamesBlt <- rownames(blupsLineTester)
        nameMBlt <- sapply(X = rowNamesBlt, FUN = function(x) {
          unlist(strsplit(x, "\\:"))
        })
        names(blupsLineTester) <- "SCA"
        lineFactor <- as.factor(nameMBlt[1, ])
        testerFactor <- as.factor(nameMBlt[2, ])
        blupsLineTester[[lines]] <- lineFactor
        blupsLineTester[[testers]] <- testerFactor
        blupsLineTester2 <- tapply(X = blupsLineTester[["SCA"]],
                                   INDEX = blupsLineTester[, c(lines, testers)],
                                   FUN = identity)
        print(format(blupsLineTester2, digits = 4), quote = FALSE)
      }
    }
  }
  res$blupsLine <- blupsLine
  res$blupsLineTester <- blupsLineTester
  return(as.list(res))
}
