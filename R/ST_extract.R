#' Extracts statistics from the fitted model results
#'
#' This function is to extract and calculate various results such as
#' heritabilities, genotypic means, unit errors etc.
#'
#' @param mixfix A list of model results with fields \code{mmix}, \code{mfix} and \code{Data}.
#' @return A list of extracted statistics.
#' @seealso
#' \code{\link{ST.run.model}}, \code{\link{ST.mod.rcbd}}, \code{\link{ST.mod.alpha}} and \code{\link{ST.mod.rowcol}}
#' @examples
#' mydat <- ST.read.csv(file.path(path.package("RAP"),"SB_yield.csv"),
#'                      factor.names=c("Env","Genotype","Rep","Subblock","Row","Column"),
#'                      trait.names="yield", env ="Env",rowSelect="HEAT05",
#'                      colSelect=c("Env","Genotype","Rep","Row","Column","yield"))
#' mymodel <- ST.run.model(mydat, design="res.rowcol", trait="yield",
#'                         genotype="Genotype", rep="Rep", row="Row",
#'                         col="Column", tryspatial=NA)
#' ## This crashes
#' ## extr <- ST.extract(mymodel)
#' ## str(extr)
#'
#' @importFrom methods slot
#' @export

ST.extract = function(mixfix) {

  # Choose between asreml and mle4
  engine = attr(mixfix, "Engine")
  trait = attr(mixfix,"Trait")
  genotype <- attr(mixfix, "genotype")
  rep <- attr(mixfix, "rep")


  # Extract statistics from mixed and fixed models in list mixfix
  mr = mixfix$mmix
  mf = mixfix$mfix
  Y = mixfix$Data

  # Use lme4 as an engine for mixed modelling
  if (engine=="lme4"){
      suppressPackageStartupMessages(require(lme4))
      # Extract coeffcients mf
      if(class(mf) == 'lmerMod') fe = lme4::fixef(mf)
      if(class(mf) == 'mer') fe = slot(mf,"fixef")
      if(class(mf) == 'lm') fe = mf$coeff
      if (!is.null(rep))
        rr = grep(rep, names(fe))
      else
        rr = integer(0)
      rg = grep(genotype, names(fe))
      rNA= which(is.na(fe))
      cr = fe[rr]
      cg = c(0, fe[rg])

      # Predictions BLUEs
      if (length(cr)==0){
        blo = 0
      }else{
        blo = mean(c(cr, 0))
      }
      predictions0 = fe[1] + blo + cg
      tempNames <- as.character(sort(unique(slot(mr,"flist")[[genotype]])))
      tempNames2 <- gsub(genotype, "",names(fe)[grep(genotype, names(fe))])
      predictions <- rep(NA, length(tempNames))
      names(predictions) <- tempNames
      predictions[tempNames%in%c(as.character(tempNames[1]), tempNames2)] <- predictions0

      # Compute weights
      if (length(rNA)>0){
        rg0 <- rg <- rg[!rg%in%rNA]
      }
      V = as.matrix(vcov(mf))
      rg <- which(colnames(V)%in%names(fe)[rg])
      Vinv=try(chol2inv(chol(V)),T)
      if (!inherits(Vinv,"try-error")){
        w = diag(Vinv)
      }else{
        w = diag(solve(V))
      }
      ue = 1/w[rg] #according to the literature w[rg], otherwise 1 / w[rg]
      se = sqrt(diag(V))[rg]
      if (length(rNA)>0){
        #insert NA
        uevec <- sevec <- rep(NA, length(fe))
        uevec[rg0] <- ue
        sevec[rg0] <- se
        ue <- uevec[rg0[1]:rg0[length(rg0)]]
        se <- sevec[rg0[1]:rg0[length(rg0)]]
        names(ue) <- names(se) <- names(fe)[rg0[1]:rg0[length(rg0)]]
      }else{
        names(ue) <- names(se) <- names(fe)[rg]
      }
      tempNames <- paste0(genotype,tempNames)
      ue0 <- se0 <- rep(NA, length(tempNames))
      names(ue0) <- names(se0) <- tempNames
      ue0[tempNames%in%names(ue)] <- ue
      se0[tempNames%in%names(se)] <- se
      ue = ue0
      se = se0

      #calculate wald test for genotype coeffients
      if (length(rNA)>0){
        waldTest.geno <- wald.test(b = fe[-rNA], Sigma = vcov(mf), positions = rg)
      }else{
        waldTest.geno <- wald.test(b = fe, Sigma = vcov(mf), positions = rg)
      }
      #waldpval <- waldTest.geno$result$chi2["P"]

      #calculate Coefficient of Variation
      CV <- 100* summary(mf)[["sigma"]]/mean(fitted(mf))

      # Extract coeffcients mr
      if(class(mr) == 'lmerMod') fe = lme4::fixef(mr)
      if(class(mr) == 'mer') fe = slot(mr,"fixef")
      if (!is.null(rep))
        rr = grep(rep, names(fe))
      else
        rr = integer(0)
      cr = fe[rr]
      ng = length(unique(slot(mr,"flist")[[genotype]]))
      reff = lme4::ranef(mr,drop=T)[[genotype]]

      # Predictions BLUPs
      if (length(cr)==0){
        blo = 0
      }else{
        blo = mean(c(cr, 0))
      }
      predictions.blups = fe[1] + blo + reff
      names(predictions.blups) = sort(unique(slot(mr,"flist")[[genotype]]))

      # Compute se.blups
      if(class(mr) == 'lmerMod') se.blups = sqrt(attr(lme4::ranef(mr,condVar=T)[[genotype]], "postVar"))
      if(class(mr) == 'mer') se.blups = sqrt(attr(lme4::ranef(mr,postVar=T)[[genotype]], "postVar"))
      se.blups = as.vector(se.blups)

      # Collect results in data frame
      Stats = data.frame(names(predictions), predictions, se, predictions.blups, se.blups, ue)
      names(Stats) = c(genotype, "predicted (BLUEs)","s.e. (BLUEs)", "predicted (BLUPs)","s.e. (BLUPs)", "ue")

      # Extract variances
      R = lme4::VarCorr(mr)
      var.gen = R[[genotype]]
      var.err = attr(R, "sc")^2

      # for marginal residuals
      if(class(mf) == 'lmerMod' || class(mr) == 'mer')
        rdf <- nrow(model.frame(mf))-length(lme4::fixef(mf))
      else
        rdf <- df.residual(mf)


      # Combine results
      Result = list(Stats = Stats, var.gen = var.gen, var.err = var.err)
      # Estimating heritability on a line mean basis
      if (attr(mixfix, "Design")=="ibd"||attr(mixfix, "Design")=="rowcol")
        Result$heritability = var.gen/(var.gen + var.err)
      if (attr(mixfix, "Design")=="res.ibd"||attr(mixfix, "Design")=="res.rowcol"||attr(mixfix, "Design")=="rcbd"){
        if (!is.null(rep))
          Result$heritability = var.gen/(var.gen + (var.err/length(unique(Y[[rep]]))))
        else
          Result$heritability = var.gen/(var.gen + (var.err))
      }
      Result$fitted = fitted(mf)
      Result$resid = resid(mf)
      Result$stdres = resid(mf, scaled=T)
      if(class(mr) == 'mer') Result$rmeans = slot(mr,"mu")
      if(class(mr) == 'lmerMod') Result$rmeans = lme4::getME(mr,"mu")
      Result$ranef = reff
      Result$model = attr(mixfix, "Design")
      Result$engine = "lme4"
      Result$waldTest.geno <- waldTest.geno
      Result$CV <- CV
      Result$rdf <- rdf
  }
  # Use asreml as an engine for mixed modelling
  if (engine=="asreml"){
      if(class(mf) == 'asreml'){
          # Extract coeffcients
          fe = mf$coe$fixed
          rg = grep(genotype, names(fe))

          # Predictions
          predictions = mf$predictions$pvals$predicted.value
          se =  mf$predictions$pvals$standard.error
          names(predictions) = mf$predictions$pvals[[genotype]]
          reff = mr$coe$random[grep(genotype, names(mr$coe$random))]

          # Compute weights
          V = mf$predictions$vcov
          Vinv=try(chol2inv(chol(V)),T)
          if (!inherits(Vinv,"try-error")){
            ue = 1/diag(Vinv)
          }else{
            ue = 1/diag(solve(V))
          }

          #calculate wald test for genotype coefficients
          df <- mf$nedf
		  ok <- require(asreml, quietly = TRUE)
		  if (!ok)
			stop("asreml cannot be successfully loaded")
		  wtt <- asreml::wald.asreml(mf, ssType="conditional", denDF="numeric")
	      pos <- grep("genotype",row.names(wtt$Wald))
		  chi2 <- wtt$Wald$F.con[pos]*wtt$Wald$Df[pos]
		  prob <- 1 - pchisq(chi2, df=wtt$Wald$Df[pos])
          chiobj <- data.frame(chi2=chi2, df=wtt$Wald$Df[pos], P=prob)
		  fobj <- data.frame(Fstat=wtt$Wald$F.con[pos],df1=wtt$Wald$Df[pos],df2=wtt$Wald$denDF[pos],P=wtt$Wald$Pr[pos])
		  reswald <- list(chi2 = c(chi2 = chi2, df = wtt$Wald$Df[pos], P = prob), Ftest = c(Fstat = wtt$Wald$F.con[pos],
            df1 = wtt$Wald$Df[pos], df2 = wtt$Wald$denDF[pos], P = wtt$Wald$Pr[pos]))
		  waldTest.geno <- list(result=reswald)
          #waldpval <- waldTest.geno$result$Ftest["P"]

          #calculate Coefficient of Variation
          CV <- 100* summary(mf)$sigma/mean(fitted(mf))

          # save SED
          sed <- mf$predictions$avsed

          #calculate LSD; significance level (5%)
          lsd <- qt(.975,df)*sed
      }
      if(class(mf) == 'lm') {
          # Extract coeffcients
          fe = mf$coeff
          if (!is.null(rep))
            rr = grep(rep, names(fe))
          else
            rr = integer(0)
          rg = grep(genotype, names(fe))
          cr = fe[rr]
          cg = c(0, fe[rg])

          # Predictions
          if (length(cr)==0){
            blo = 0
          }else{
            blo = mean(c(cr, 0))
          }
          predictions = fe[1] + blo + cg
          names(predictions) = sort(unique(slot(mr,"flist")[[genotype]]))
          reff = mr$coe$random[grep(genotype, names(mr$coe$random))]

          # Compute weights
          V = vcov(mf)
          Vinv=try(chol2inv(chol(V)),T)
          if (!inherits(Vinv,"try-error")){
            w = diag(Vinv)
          }else{
            w = diag(solve(V))
          }
          ue = 1/w[rg] #according to the literature w[rg], otherwise 1 / w[rg]
          ue =  c(ue[1], ue)
          se = sqrt(diag(V))[rg]
          se = c(se[1], se)

          #calculate wald test for genotype coeffients
          df <- mf$df.residual
          waldTest.geno <- wald.test(b = fe, Sigma = V, positions = rg, df = df)
          #waldpval <- waldTest.geno$result$Ftest["P"]

          #calculate Coefficient of Variation
          CV <- 100* summary(mf)$sigma/mean(fitted(mf))

          sed <- NULL
          lsd <- NULL
      }

      # Predictions BLUPs
      predictions.blups = mr$predictions$pvals$predicted.value
      se.blups =  mr$predictions$pvals$standard.error

      # Collect results in data frame
      Stats = data.frame(names(predictions), predictions, se, predictions.blups, se.blups, ue)
      names(Stats) = c(genotype, "predicted (BLUEs)","s.e. (BLUEs)", "predicted (BLUPs)","s.e. (BLUPs)", "ue")

      # Extract variances
      genopos = grep(paste(genotype,"!",genotype,".var", sep=""), names(mr$gammas))
      var.gen = mr$gammas[genopos]*mr$sigma2
      resipos = grep("R!variance", names(mr$gammas))
      var.err = mr$gammas[resipos]*mr$sigma2

      # Combine results
      Result = list(Stats = Stats, var.gen = var.gen, var.err = var.err)
      # a generalized heritability
      #genopos = grep(genotype, names(mr$coefficients$random))
      #Result$heritability =1-mean(mr$vcoeff$random[genopos])/(var.gen/mr$sigma2)
      tmpfile <- tempfile()
      sink(file=tmpfile)
      newsed <- predict(mr, classify=genotype, only=genotype, sed=T, data=Y)$predictions$sed
      sink(NULL)
      unlink(tmpfile)
      sedsq <-newsed^2
      Result$heritability =1-mean(sedsq[lower.tri(sedsq)])/(var.gen*2)
      Result$fitted = fitted(mf)
      if(class(mf) == 'asreml') {
        mf$call$data <- substitute(Y)
        Result$resid = residuals(mf, type="response")
        Result$stdres = residuals(mf, type="stdCond" )
        Result$rdf = mf$nedf
      }
      if(class(mf) == 'lm') {
        Result$resid = resid(mf)
        Result$stdres = rstandard(mf)
        Result$rdf = df.residual(mf)
      }
      Result$rmeans = fitted(mr)
      Result$predictions.sed = sed
      Result$predictions.lsd = lsd
      Result$ranef = reff
      Result$model = attr(mixfix, "Design") #"alpha"
      Result$engine = "asreml"
      Result$waldTest.geno <- waldTest.geno
      Result$CV <- CV
  }
  if (!is.null(rep)){
    myvec <-  tapply(Y[,rep], Y[,genotype], nlevels)
    Result$min.reps <- min(myvec, na.rm = T)
    Result$mean.reps <- mean(myvec, na.rm = T)
    Result$max.reps <- max(myvec, na.rm = T)
  }else{
    Result$min.reps <- 1
    Result$mean.reps <- 1
    Result$max.reps <- 1
  }
  return(Result)
}
