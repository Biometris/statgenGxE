#' Extracts statistics from the fitted model results
#'
#' This function is to extract and calculate various results such as
#' heritabilities, genotypic means, unit errors etc.
#'
#' @param mixfix A list of model results with fields \code{mMix}, \code{mFix} and \code{data}.
#' @return A list of extracted statistics.
#' @seealso
#' \code{\link{ST.run.model}}, \code{\link{ST.mod.rcbd}}, \code{\link{ST.mod.alpha}} and
#' \code{\link{ST.mod.rowcol}}
#' @examples
#' mydat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames=c("Env","Genotype","Rep","Subblock","Row","Column"),
#'                      traitNames="yield", env ="Env",rowSelect="HEAT05",
#'                      colSelect=c("Env","Genotype","Rep","Row","Column","yield"))
#' mymodel <- ST.run.model(mydat, design="res.rowcol", trait="yield",
#'                         genotype="Genotype", rep="Rep", row="Row",
#'                         col="Column", tryspatial=NA)
#' ##extr <- ST.extract(mymodel)
#' ##str(extr)
#'
#' @importFrom methods slot
#' @export

ST.extract = function(mixfix) {
  # Choose between asreml and mle4
  engine <- attr(mixfix, "Engine")
  trait <- attr(mixfix, "Trait")
  genotype <- attr(mixfix, "genotype")
  rep <- attr(mixfix, "rep")
  # Extract statistics from mixed and fixed models in list mixfix
  mr <- mixfix$mMix
  mf <- mixfix$mFix
  Y <- mixfix$data
  # Use lme4 as an engine for mixed modelling
  if (engine == "lme4") {
    # Extract coeffcients mf
    if (class(mf) == "lmerMod") {
      fe <- lme4::fixef(mf)
    } else if (class(mf) == "mer") {
      fe <- slot(mf, "fixef")
    } else if (class(mf) == "lm") {
      fe <- mf$coeff
    }
    if (!is.null(rep)) {
      rr <- grep(pattern = rep, x = names(fe))
    } else {
      rr <- integer(0)
    }
    rg <- grep(pattern = genotype, x = names(fe))
    rNA <- which(is.na(fe))
    cr <- fe[rr]
    cg <- c(0, fe[rg])
    # Predictions BLUEs
    if (length(cr) == 0) {
      blo <- 0
    } else {
      blo <- mean(c(cr, 0))
    }
    predictions0 <- fe[1] + blo + cg
    tempNames <- as.character(sort(unique(slot(mr, "flist")[[genotype]])))
    tempNames2 <- gsub(pattern = genotype, replacement =  "",
                       x = names(fe)[grep(pattern = genotype, x = names(fe))])
    predictions <- rep(NA, length(tempNames))
    names(predictions) <- tempNames
    predictions[tempNames %in% c(as.character(tempNames[1]), tempNames2)] <- predictions0
    # Compute weights
    if (length(rNA) > 0) {
      rg0 <- rg <- rg[!rg %in% rNA]
    }
    V <- as.matrix(vcov(mf))
    rg <- which(colnames(V) %in% names(fe)[rg])
    Vinv <- try(chol2inv(chol(V)), TRUE)
    if (!inherits(Vinv, "try-error")) {
      w <- diag(Vinv)
    } else {
      w <- diag(solve(V))
    }
    ue <- 1 / w[rg] #according to the literature w[rg], otherwise 1 / w[rg]
    se <- sqrt(diag(V))[rg]
    if (length(rNA) > 0) {
      #insert NA
      uevec <- sevec <- rep(NA, length(fe))
      uevec[rg0] <- ue
      sevec[rg0] <- se
      ue <- uevec[rg0[1]:rg0[length(rg0)]]
      se <- sevec[rg0[1]:rg0[length(rg0)]]
      names(ue) <- names(se) <- names(fe)[rg0[1]:rg0[length(rg0)]]
    } else {
      names(ue) <- names(se) <- names(fe)[rg]
    }
    tempNames <- paste0(genotype,tempNames)
    ue0 <- se0 <- rep(NA, length(tempNames))
    names(ue0) <- names(se0) <- tempNames
    ue0[tempNames %in% names(ue)] <- ue
    se0[tempNames %in% names(se)] <- se
    ue <- ue0
    se <- se0
    #calculate wald test for genotype coeffients
    if (length(rNA) > 0) {
      waldTestGeno <- waldTest(b = fe[-rNA], Sigma = vcov(mf), positions = rg)
    } else {
      waldTestGeno <- waldTest(b = fe, Sigma = vcov(mf), positions = rg)
    }
    #waldpval <- waldTestGeno$result$chi2["P"]
    #calculate Coefficient of Variation
    CV <- 100 * summary(mf)[["sigma"]] / mean(fitted(mf))
    # Extract coeffcients mr
    if (class(mr) == "lmerMod") {
      fe <- lme4::fixef(mr)
    } else if(class(mr) == "mer") {
      fe <- slot(mr, "fixef")
    }
    if (!is.null(rep)) {
      rr <- grep(pattern = rep, x = names(fe))
    } else {
      rr <- integer(0)
    }
    cr <- fe[rr]
    ng <- length(unique(slot(mr, "flist")[[genotype]]))
    rEff <- lme4::ranef(mr, drop = TRUE)[[genotype]]
    # Predictions BLUPs
    if (length(cr) == 0) {
      blo <- 0
    } else {
      blo <- mean(c(cr, 0))
    }
    predictionsBlups <- fe[1] + blo + rEff
    names(predictionsBlups) <- sort(unique(slot(mr, "flist")[[genotype]]))
    # Compute seBlups
    if (class(mr) == "lmerMod") {
      seBlups = sqrt(attr(lme4::ranef(mr, condVar = TRUE)[[genotype]], "postVar"))
    } else if(class(mr) == "mer") {
      seBlups = sqrt(attr(lme4::ranef(mr, postVar = TRUE)[[genotype]], "postVar"))
    }
    seBlups = as.vector(seBlups)
    # Collect results in data frame
    stats = data.frame(names(predictions), predictions, se, predictionsBlups, seBlups, ue)
    names(stats) = c(genotype, "predicted (BLUEs)","s.e. (BLUEs)",
                     "predicted (BLUPs)","s.e. (BLUPs)", "ue")
    # Extract variances
    R = lme4::VarCorr(mr)
    varGen = R[[genotype]]
    varErr = attr(R, "sc") ^ 2
    # for marginal residuals
    if (class(mf) %in% c("lmerMod", "mer")) {
      rdf <- nrow(model.frame(mf)) - length(lme4::fixef(mf))
    } else {
      rdf <- df.residual(mf)
    }
    # Combine results
    result <- list(stats = stats, varGen = varGen, varErr = varErr)
    # Estimating heritability on a line mean basis
    if (attr(mixfix, "Design") %in% c("ibd", "rowcol")) {
      result$heritability <- varGen/(varGen + varErr)
    } else if (attr(mixfix, "Design") %in% c("res.ibd", "res.rowcol", "rcbd")) {
      if (!is.null(rep)) {
        result$heritability <- varGen / (varGen + (varErr / length(unique(Y[[rep]]))))
      } else {
        result$heritability <- varGen / (varGen + varErr)
      }
    }
    result$fitted <- fitted(mf)
    result$resid <- resid(mf)
    result$stdres <- resid(mf, scaled = TRUE)
    if (class(mr) == "mer") {
      result$rmeans <- slot(mr, "mu")
    } else if (class(mr) == "lmerMod") {
      result$rmeans <- lme4::getME(mr, "mu")
    }
    result$ranef <- rEff
    result$model <- attr(mixfix, "Design")
    result$engine <- "lme4"
    result$waldTestGeno <- waldTestGeno
    result$CV <- CV
    result$rdf <- rdf
  } else if (engine == "asreml") {
    # Use asreml as an engine for mixed modelling
    if (class(mf) == "asreml") {
      # Extract coeffcients
      fe <- mf$coe$fixed
      rg <- grep(pattern = genotype, x = names(fe))
      # Predictions
      predictions <- mf$predictions$pvals$predicted.value
      se <- mf$predictions$pvals$standard.error
      names(predictions) <- mf$predictions$pvals[[genotype]]
      rEff <- mr$coe$random[grep(pattern = genotype, x =names(mr$coe$random))]
      # Compute weights
      V <- mf$predictions$vcov
      Vinv <- try(chol2inv(chol(V)), silent = TRUE)
      if (!inherits(Vinv, "try-error")) {
        ue <- 1 / diag(Vinv)
      } else {
        ue <- 1 / diag(solve(V))
      }
      #calculate wald test for genotype coefficients
      df <- mf$nedf
      if (!requireNamespace("asreml", quietly = TRUE)) {
        stop("asreml cannot be successfully loaded.\n")
      }
      wtt <- asreml::wald.asreml(mf, ssType = "conditional", denDF = "numeric")
      pos <- grep(pattern = "genotype", x = row.names(wtt$Wald))
      chi2 <- wtt$Wald$F.con[pos] * wtt$Wald$Df[pos]
      prob <- 1 - pchisq(q = chi2, df = wtt$Wald$Df[pos])
      resWald <- list(chi2 = c(chi2 = chi2, df = wtt$Wald$Df[pos], P = prob),
                      Ftest = c(Fstat = wtt$Wald$F.con[pos],
                                df1 = wtt$Wald$Df[pos],
                                df2 = wtt$Wald$denDF[pos],
                                P = wtt$Wald$Pr[pos]))
      waldTestGeno <- list(result = resWald)
      #waldpval <- waldTestGeno$result$Ftest["P"]
      #calculate Coefficient of Variation
      CV <- 100 * summary(mf)$sigma / mean(fitted(mf))
      # save SED
      sed <- mf$predictions$avsed
      #calculate LSD; significance level (5%)
      lsd <- qt(p = .975, df = df) * sed
    } else if(class(mf) == "lm") {
      # Extract coeffcients
      fe <- mf$coeff
      if (!is.null(rep)) {
        rr <- grep(pattern = rep, x = names(fe))
      } else {
        rr <- integer(0)
      }
      rg <- grep(pattern = genotype, x = names(fe))
      cr <- fe[rr]
      cg <- c(0, fe[rg])
      # Predictions
      if (length(cr) == 0) {
        blo <- 0
      } else {
        blo <- mean(c(cr, 0))
      }
      predictions <- fe[1] + blo + cg
      names(predictions) <- sort(unique(slot(mr, "flist")[[genotype]]))
      rEff <- mr$coe$random[grep(pattern = genotype, x = names(mr$coe$random))]
      # Compute weights
      V <- vcov(mf)
      Vinv <- try(chol2inv(chol(V)), silent = TRUE)
      if (!inherits(Vinv, "try-error")) {
        w <- diag(Vinv)
      } else {
        w <- diag(solve(V))
      }
      ue <- 1 / w[rg] #according to the literature w[rg], otherwise 1 / w[rg]
      ue <- c(ue[1], ue)
      se <- sqrt(diag(V))[rg]
      se <- c(se[1], se)
      #calculate wald test for genotype coeffients
      df <- mf$df.residual
      waldTestGeno <- waldTest(Sigma = V, b = fe, positions = rg, df = df)
      #waldpval <- waldTestGeno$result$Ftest["P"]
      #calculate Coefficient of Variation
      CV <- 100 * summary(mf)$sigma / mean(fitted(mf))
      sed <- NULL
      lsd <- NULL
    }
    # Predictions BLUPs
    predictionsBlups <- mr$predictions$pvals$predicted.value
    seBlups <- mr$predictions$pvals$standard.error
    # Collect results in data frame
    stats <- data.frame(names(predictions), predictions, se, predictionsBlups, seBlups, ue)
    names(stats) <- c(genotype, "predicted (BLUEs)","s.e. (BLUEs)",
                      "predicted (BLUPs)","s.e. (BLUPs)", "ue")
    # Extract variances
    genopos <- grep(pattern = paste0(genotype, "!", genotype, ".var"), x = names(mr$gammas))
    varGen <- mr$gammas[genopos] * mr$sigma2
    resipos <- grep(pattern = "R!variance", x = names(mr$gammas))
    varErr <- mr$gammas[resipos] * mr$sigma2
    # Combine results
    result <- list(stats = stats, varGen = varGen, varErr = varErr)
    # a generalized heritability
    #genopos = grep(genotype, names(mr$coefficients$random))
    #result$heritability =1-mean(mr$vcoeff$random[genopos])/(varGen/mr$sigma2)
    tmpfile <- tempfile()
    sink(file = tmpfile)
    newSed <- predict(mr, classify = genotype, only = genotype, sed = TRUE, data = Y)$predictions$sed
    sink(NULL)
    unlink(tmpfile)
    sedSq <- newSed ^ 2
    result$heritability <- 1 - mean(sedSq[lower.tri(sedSq)]) / (varGen * 2)
    result$fitted <- fitted(mf)
    if (class(mf) == "asreml") {
      mf$call$data <- substitute(Y)
      result$resid <- residuals(mf, type = "response")
      result$stdres <- residuals(mf, type = "stdCond" )
      result$rdf <- mf$nedf
    } else if(class(mf) == "lm") {
      result$resid <- resid(mf)
      result$stdres <- rstandard(mf)
      result$rdf <- df.residual(mf)
    }
    result$rmeans <- fitted(mr)
    result$predictionsSed <- sed
    result$predictionsLsd <- lsd
    result$ranef <- rEff
    result$model <- attr(mixfix, "Design") #"alpha"
    result$engine <- "asreml"
    result$waldTestGeno <- waldTestGeno
    result$CV <- CV
  }
  if (!is.null(rep)) {
    myvec <- tapply(X = Y[, rep], INDEX = Y[, genotype], FUN = nlevels)
    result$minReps <- min(myvec, na.rm = TRUE)
    result$meanReps <- mean(myvec, na.rm = TRUE)
    result$maxReps <- max(myvec, na.rm = TRUE)
  } else {
    result$minReps <- 1
    result$meanReps <- 1
    result$maxReps <- 1
  }
  return(result)
}
