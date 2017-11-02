#' Modified joint regression analysis
#'
#' This function performs a modified joint analysis of data classified by two factors.
#'
#' @param Y A data frame object.
#' @param trait A character string specifying a trait column of the data.
#' @param genotype A character string specifying a genotype column of the data.
#' @param env A character string specifying the envirnoment column of the data.
#' @param maxcycle A integer specifying the maximum number of iterations to be achieved. By default, \code{maxcycle = 15}.
#' @param tol A small positive numerical value specifying convergence tolerance. By default, \code{tol = 0.001}.
#' @param sortBYsens A character string specifying whether the results are to be sorted in an increasing (or decreasing) order of sensitivities.
#' By default, \code{sortBYsens = "ascending"}. Other options are "descending" and NA.
#' @param scatterplot A logical value specifying if a scatterplot of sensitivities is produced.
#' @param lineplot A logical value specifying if a fitted line for each genotype is produced.
#' @param trellisplot A logical value specifying if the trellis plot of the individual genotype slopes is produced.
#' @examples
#' mydat <- GE.read.csv(file.path(path.package("RAP"),"F2maize_pheno.csv"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' names(mydat)=c("env", "genotype","yld")
#' fw.anlysis <- GE.FW(mydat, trait="yld", genotype="genotype", env="env", maxcycle = 15, tol = 0.001,
#'                sortBYsens = "ascending")
#' fw.anlysis
#'
#' @import graphics grDevices
#' @export

GE.FW <- function(Y, trait, genotype, env, maxcycle = 15, tol = 0.001,
                  sortBYsens = c("ascending", "descending", NA),
                  scatterplot = T, lineplot = F, trellisplot = F){
  ## ZZ replicates RJOINT procedure in GenStat
  ## Handling missing values ?na.action = na.exclude?

  if (missing(sortBYsens)){
    sortBYsens <- "ascending"
  }
  nlab <- nlevels(Y[,genotype])
  nenvs <- nlevels(Y[,env])

  # Save Sum Sq & Df
  rdev <- rdf <- rep(NA,5)

  # Initializing...
  # Estimating env effects with the sensitivity beta=1
  model0 <- lm(as.formula(paste(trait, "~-1 +", env, "+", genotype)), data = Y, na.action = na.exclude)
  aov.0 <-  anova(model0)
  tpos <- rownames(aov.0)=="Residuals"
  rdev[2] <- aov.0[["Sum Sq"]][tpos]
  rdf[2]  <- aov.0[["Df"]][tpos]
  coeffs.model0 <- coefficients(model0)
  enveffs0 <-  coeffs.model0[grep(env,names(coeffs.model0))]
  # Adjust env effects to have mean zero
  enveffs0 <- enveffs0 - mean(enveffs0, na.rm = T)
  enveffs <- rep(NA,nrow(Y))
  for (ii in names(enveffs0)){
    tposition <- sapply(Y[[env]],function(x) grepl(x, ii))
    enveffs[tposition] <- enveffs0[ii]
  }
  Y <- cbind(Y,enveffs)
  # Initial values for sensitivity beta
  Y$beta <- rep(1, nlab)
  # Set a relative difference to be large
  maxdiff <- 10000
  # Set iteration to be 1
  iter <- 1
  # Iterate to fit 'y(i,j) = genmean(i)+beta(i)*enveffs(j)'
  while (maxdiff > tol && iter <= maxcycle){
    beta0 <- Y$beta
    # Form variate with current genotype sensitivity relevant to each unit
    model1 <- lm(as.formula(paste(trait, "~-1 +", genotype, "+", genotype,":enveffs")), data = Y, na.action = na.exclude)
    coeffs.model1 <- coefficients(model1)
    # Update beta
    Y$beta <- coeffs.model1[match(paste(genotype,Y[[genotype]],":enveffs",sep=""), names(coeffs.model1))]
    Y$beta <- Y$beta/mean(Y$beta, na.rm = T)
    print(logLik(model1))

    # Form variate with current env means relevant to each unit
    model2 <- lm(as.formula(paste(trait, "~-1 +", genotype, "+", env,":beta")), data = Y, na.action = na.exclude)
    coeffs.model2 <- coefficients(model2)

    # Update enveffs
    Y$enveffs <- coeffs.model2[match(paste(env,Y[[env]],":beta",sep=""), names(coeffs.model2))]
    Y[is.na(Y$enveffs),"enveffs"] <- 0
    Y$enveffs <- Y$enveffs-mean(Y$enveffs)

    # Maximum difference of sensitivities between the succesive iterations
    maxdiff <- max(abs(Y$beta - beta0), na.rm = T)
    if (iter == maxcycle && maxdiff > tol)
      warning(paste('Convergence not achieved in ',iter,' iterations. Tolerance ',
            tol,', criterion at last iteration ',signif(maxdiff,4),'.\n',sep=""))
    iter <- iter + 1
  }


  # Environments
  aov.1 <- anova(model1)
  tpos <- rownames(aov.1)=="Residuals"
  rdev[4] <- aov.1[["Sum Sq"]][tpos]
  rdf[4]  <- aov.1[["Df"]][tpos]

  # Extract total deviance
  modelA <- lm(as.formula(paste(trait, "~", genotype)), data = Y, na.action = na.exclude)
  aov.A <- anova(modelA)
  rdev[5] <- sum(aov.A[["Sum Sq"]])
  rdf[5]  <- sum(aov.A[["Df"]])

  # Fit varieties only for first entry in aov
  modelB <- lm(as.formula(paste(trait, "~-1 +", genotype)), data = Y, na.action = na.exclude)
  aov.B <- anova(modelB)
  tpos <- rownames(aov.B)=="Residuals"
  rdev[1] <- aov.B[["Sum Sq"]][tpos]
  rdf[1]  <- aov.B[["Df"]][tpos]

  # Calculate deviances and d.f.
  rdev[c(3,2,1)] <- rdev[c(2,1,5)]-rdev[c(4,2,1)]
  rdf[c(2,1)] <- rdf[c(1,5)]-rdf[c(2,1)]
  rdf[3] <- rdf[1]
  rdf[4] <- rdf[5]-rdf[1]-rdf[2]-rdf[3]
  # Calculate mean deviances and F statistics
  mdev <- rdev/rdf
  rmdev <- mdev[4]
  devr <- mdev/rmdev
  devr[c(4,5)] <-  NA
  devr[rdf %in% 0] <- NA
  devr[!is.na(devr) & devr <0] <- NA
  fprob <- pf(devr, rdf, rdf[4], lower.tail = FALSE)
  aov.table <- data.frame("Df"=rdf, "Sum Sq"=rdev, "Mean Sq"=mdev, "F value"=devr, "Pr(>F)"=fprob,
   row.names=c(genotype, env, "Sensitivities", "Residual", "Total"), check.names=F)

  # Sensitivity beta(i)
  sens <- Y$beta
  # standard error for env
  sigma.e <- sqrt(diag(vcov(model1))[match(paste(genotype,Y[[genotype]],":enveffs",sep=""), names(coeffs.model1))])
  # Mean for each genotype genmean(i)
  genmean <- coeffs.model1[match(paste(genotype,Y[[genotype]],sep=""), names(coeffs.model1))]
  # residual standard error
  sigma <- sqrt(diag(vcov(model1))[match(paste(genotype,Y[[genotype]],sep=""), names(coeffs.model1))])

  fitted.gen <- fitted(model1)
  resi.gen <- residuals(model1)
  fval <- tapply(fitted.gen, Y[,c(genotype,env)], function(x) mean(x, na.rm = T))
  # mean squared error (MSE) of the trait means for each genotype
  mse <- tapply(resi.gen, Y[,genotype],
                function(x){
                  checkg <- length(x)
                  if (checkg>2){
                    sum(x^2, na.rm = T)/(checkg-2)
                  }else{
                    rep(NA, checkg)
                  }
                }
               )
  # mean estimates for each genotype
  G <- levels(Y[,genotype])
  sens <- tapply(sens, Y[,genotype], function(x) mean(x, na.rm = T))
  sigma.e <- tapply(sigma.e, Y[,genotype], function(x) mean(x, na.rm = T))
  genmean <- tapply(genmean, Y[,genotype], function(x) mean(x, na.rm = T))
  sigma <- tapply(sigma, Y[,genotype], function(x) mean(x, na.rm = T))
  if (sortBYsens == "ascending"){
    order.sens <- order(sens)
    res <- data.frame(G,sens, sigma.e, genmean, sigma, mse, row.names=1:length(sens))[order.sens,]
  }else{
    if (sortBYsens == "descending"){
      order.sens <- order(sens, decreasing = T)
      res <- data.frame(G,sens, sigma.e, genmean, sigma, mse, row.names=1:length(sens))[order.sens,]
    }else{
      res <- data.frame(G,sens, sigma.e, genmean, sigma, mse, row.names=1:length(sens))
    }
  }

  # ANOVA table

  # Environment effects
  match.positions <- match(paste(env,levels(Y[[env]]),":beta",sep=""), names(coeffs.model2))
  Enveffs <- coeffs.model2[match.positions]
  na.position <- is.na(Enveffs)
  Enveffs[na.position] <- 0
  Enveffs <- Enveffs - mean(Enveffs)
  var.Enveffs <- diag(vcov(model2)[match.positions[!na.position], match.positions[!na.position]])
  se.Enveffs <- rep(NA, length(Enveffs))
  names(se.Enveffs) <- names(Enveffs)
  se.Enveffs[names(var.Enveffs)] <- sqrt(var.Enveffs)
  match.positions2 <- match(paste(env,levels(Y[[env]]),":beta",sep=""), names(Enveffs))
  if (!is.null(model1$na.action))
    means.fitted <- tapply(model1$fitted, Y[[env]][-model1$na.action], mean)
  else
    means.fitted <- tapply(model1$fitted, Y[[env]], mean)
  means.fitted <- means.fitted[match.positions2]
  means.rank <- rank(-means.fitted)
  means.names <- names(means.fitted)
  Enveffs.summary <- data.frame(Environment=means.names, Effect=Enveffs,
   s.e.=se.Enveffs, Mean=means.fitted, Rank=means.rank, row.names=NULL)

  # Draw various plots
  if (scatterplot == T){
    if (!all(is.na(mse))){
      xx <- cbind(genmean, mse, sens)
      colnames(xx) <- c("Mean", "m.s.deviation", "Sensitivity")
      pairs(xx, upper.panel = NULL, main = "Finlay & Wilkinson analysis")
    }else{
      plot(x = genmean,y = sens, xlab = "Mean", ylab = "Sensitivity",
      main = "Finlay & Wilkinson analysis")
    }
  }

  if (lineplot == T){
    min.fval <- min(fitted.gen, na.rm = T)
    max.fval <- max(fitted.gen, na.rm = T)
    xeff <- Enveffs
    min.xeff <- min(xeff, na.rm = T)
    max.xeff <- max(xeff, na.rm = T)
    dev.new()
    plot(NA, xlim = c(min.xeff, max.xeff),ylim = c(min.fval,max.fval),
    ylab = trait, xlab = "Environment", xaxt="n")
    axis (1, xeff, levels(Y[,env]), las = 2, cex.axis =.75)
    colour <- 1
    for (i in 1:nlab) {
      if (!is.na(sortBYsens)){
        xfval <- fval[order.sens[i],]
      }else{
        xfval <- fval[i,]
      }
      lines(xeff[order(xeff)], xfval[order(xeff)], col = colour)
      colour <- colour + 1
    }
  }

  if (trellisplot == T){
    trellisdata <- data.frame(genotype=Y[[genotype]], trait=Y[[trait]], fitted=fitted.gen, xeff=rep(Enveffs,nlab))
    if (nlab>64){
      first64 <- levels(Y[[genotype]])[1:64]
      first64 <- Y[[genotype]] %in% first64
      trellisdata <- droplevels(trellisdata[first64,])
    }
    dev.new()
    print(lattice::xyplot(trait + fitted ~ xeff | genotype, data=trellisdata,
                 panel = function(x,y, subscripts) {
                   lattice::panel.xyplot(x,y)
                   lattice::panel.lines(trellisdata$xeff[subscripts], trellisdata$fitted[subscripts])
                 }, as.table=T, subscripts=T, xlab="Environment", ylab=trait,
                 main=paste0("Finlay & Wilkinson analysis for ",trait)))
  }

  info <- data.frame("Response variate"=trait, "Number of genotypes"=nlab, "Number of environments"=nenvs,
  "Convergence criterion"= tol, "Number of iterations"=iter-1, row.names=NULL, check.names=F, stringsAsFactors = FALSE)
  result <- new.env()
  # Returns statistics
  result$estimates <- res
  result$aov.table <- aov.table
  result$EnvEffs <- Enveffs.summary
  result <- as.list(result)
  attr(result, "info") <- info
  result
}
