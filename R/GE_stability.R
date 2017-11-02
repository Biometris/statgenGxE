#' Calculates stability coefficients for genotype-by-environment data
#'
#' This function calculate difference measures of stability, such as the cultivar-superiority measure of Lin & Binns (1988),
#' Shukla's (1972) stability variance and Wricke's (1962) ecovalence.
#'
#' @param Y A data frame object.
#' @param trait A character string specifying a trait column of the data.
#' @param genotype A character string specifying a genotype column of the data.
#' @param env A character string specifying the environment column of the data.
#' @param method A character string specifying (a) measure(s) of stability. By default, \code{method = c("superiority","static","wricke")}.
#' @param superiority.bestmethod A character string specifying the criterion to define the best genotype ("max","min").
#' By default, \code{superiority.bestmethod="max"}.
#' @param sorted A character string specifying whether the results are sorted by increasing (or decreasing) order of stability coefficients.
#' By default, \code{sortBYsens = "descending"}. Other options are \code{"ascending"} and \code{NA}.
#' @param plot A logical value indicating whether to produce a plot of stability against the mean.
#' @return A list of measures of stability.
#' @examples
#' mydat <- GE.read.csv(file.path(path.package("RAP"),"F2maize_pheno.csv"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' names(mydat)=c("env", "genotype","yld") 
#' GE.stability(Y=mydat, trait="yld", genotype="genotype", env="env",
#'              method = "superiority", superiority.bestmethod = "max",
#'              sorted = "descending", plot=TRUE)
#'              
#' @export              
    
GE.stability <- function(Y, trait, genotype, env, method = c("superiority","static","wricke"),
                         superiority.bestmethod = "max", sorted = c("ascending", "descending", NA),
                         plot=T){
  if (missing(sorted)){
    sorted <- "descending"
  }

  if (any(is.na(Y[[trait]]))){
    y0 <- tapply(Y[[trait]], Y[,c(genotype,env)], mean)
    yindex <- tapply(1:nrow(Y), Y[,c(genotype,env)], identity)
    na_yes_no <- is.na(y0)
    # imputaion
    y1 <- RAP.multmissing(y0, maxcycle = 10, na.strings=NA)
    replace.val <- y1[na_yes_no]
    yindex.replace <- yindex[na_yes_no]
    if (is.list(yindex.replace)){
      nyr <- length(yindex.replace)
      for (ii in 1:nyr){
        for (jj in 1:length(Y[yindex.replace[[ii]], trait])){
          if (is.na(Y[yindex.replace[[ii]][jj], trait]))
            Y[yindex.replace[[ii]][jj], trait] <- replace.val[ii]
        }
      }
    }else{
      Y[yindex.replace, trait] <- replace.val
    }
  }
  
  lab <- levels(Y[,genotype])
  nlab <- length(lab)
  envs <- levels(Y[,env])
  nenvs <- length(envs)

  # The centered trait mean for eniroment j
  Ej <- sapply(envs, function(x) mean(Y[which(Y[,env] == x),trait], na.rm = T))

  # Maximum or minimum trait mean among all genotype in jth enviroment
  if (superiority.bestmethod == "max"){
    Mj <- sapply(envs, function(x) max(Y[which(Y[,env] == x),trait], na.rm = T))
  }else{
    Mj <- sapply(envs, function(x) min(Y[which(Y[,env] == x),trait], na.rm = T))
  }

  # The genotype trait mean across environments
  Ei <- sapply(lab, function(x) mean(Y[which(Y[,genotype] == x),trait], na.rm = T))
 
  # The grand mean
  E <- mean(Y[,trait], na.rm = T)
  
  W <- S <- LB <- rep(NA,nlab)
  for ( i in 1:nlab){
    # Observed genotype field response in the enviroment j
    # (averaged across experiment replicates)
    
    Rij <- sapply(envs, function(x) mean(Y[which(Y[,genotype]==lab[i]&Y[,env]==x),trait], na.rm = T))
    pos <- (1:nenvs)[!is.na(Rij)]
    
    # Static measure (Shukla's (1972a) stability variance)
    if ("static" %in% method) {
      if (length(pos)==0){
        S[i] <- NA
      }else{
        S[i] <-  sum(sapply(pos, function(j) (Rij[j]-Ei[i])^2)/(nenvs-1))
      }
    }
    
    # Superiority measure (LIN&BINNS 1988)
    if ("superiority" %in% method){
      if (length(pos)==0){
        LB[i] <- NA
      }else{
        LB[i] <- sum(sapply(pos, function(j) (Rij[j] - Mj[j])^2)/ (2*nenvs))
      }
    }
    
    # Wricke's (1962) ecovalence
    if ("wricke" %in% method){
      if (length(pos)==0){
        W[i] <- NA
      }else{
        W[i] <- sum(sapply(pos, function(j) (Rij[j] - Ei[i] - Ej[j] + E)^2))
      }
    }
  }
  #Calculate trait mean per genotype
  means <- sapply(lab, function(x) mean(Y[Y[,genotype]==x, trait], na.rm=T))
  
  if (plot){
    # Prepare panels
    mflag <- c("static", "superiority", "wricke") %in% method
    if (sum(mflag)!=1){
      if (sum(mflag)==3){
        opar <- par(mfrow = c(2, 2), mar = c(4, 4, 4, 2), oma = c(.5,.5,2,.3))
      }else{
        opar <- par(mfrow=c(1,2), mar = c(4, 4, 4, 2), oma = c(.5,.5,2,.3))
      }
    }else{
      opar <- par(mfrow=c(1,1))
    }
    on.exit(par(opar))
    
    if ("superiority" %in% method){
      plot(x=means, y=LB, xlab="Mean", ylab="Cultivar superiority")
    }
    if ("static" %in% method){
      plot(x=means, y=S, xlab="Mean", ylab="Static stability")
    }
    if ("wricke" %in% method){
      plot(x=means, y=W, xlab="Mean", ylab="Wricke's ecovalence")
    }
    if (sum(mflag)==1){
      title(paste('Stability coefficients for', trait))
    }else{
      mtext(paste('Stability coefficients for', trait), side = 3, outer = TRUE, cex = 1.3, font=2)
    }
  }
    
  if (is.na(sorted)){
    if ("static" %in% method) {
      S.out <- data.frame(lab,S, means, row.names=1:nlab)
      names(S.out) <- c("genotype","static", "mean")
    }
    if ("superiority" %in% method){
      LB.out <- data.frame(lab,LB, means, row.names=1:nlab)
      names(LB.out) <- c("genotype", "superiority", "mean")
    }
    if ("wricke" %in% method){
      W.out <- data.frame(lab,W,means, row.names=1:nlab)
      names(W.out) <- c("genotype", "wricke", "mean")
    }
  }else{
    if (sorted == "ascending"){
      if ("static" %in% method){
        orderS <- order(S)
        S.out <- data.frame(lab,S,means, row.names=1:nlab)[orderS,]
        names(S.out) <- c("genotype","static", "mean")
      }
      if ("superiority" %in% method){
        orderLB <- order(LB)
        LB.out <- data.frame(lab,LB,means, row.names=1:nlab)[orderLB,]
        names(LB.out) <- c("genotype", "superiority", "mean")
      }
      if ("wricke" %in% method){
        orderW <- order(W)
        W.out <- data.frame(lab,W,means, row.names=1:nlab)[orderW,]
        names(W.out) <- c("genotype", "wricke", "mean")
      }
    }else{
      if (sorted == "descending"){
        if ("static" %in% method){
          orderS <- order(S, decreasing = T)
          S.out <- data.frame(lab,S,means, row.names=1:nlab)[orderS,]
          names(S.out) <- c("genotype","static", "mean")
        }
        if ("superiority" %in% method){
          orderLB <- order(LB, decreasing = T)
          LB.out <- data.frame(lab,LB,means, row.names=1:nlab)[orderLB,]
          names(LB.out) <- c("genotype", "superiority", "mean")
        }
        if ("wricke" %in% method){
          orderW <- order(W, decreasing = T)
          W.out <- data.frame(lab,W,means, row.names=1:nlab)[orderW,]
          names(W.out) <- c("genotype", "wricke", "mean")
        }
      }else{
        stop('options for sorted must be one of "ascending", "descending" and NA')
      }
    }
  }
  
  res <- new.env()
  if ("static" %in% method) res$static <- S.out
  if ("superiority" %in% method) res$superiority <- LB.out
  if ("wricke" %in% method) res$wricke <- W.out
  as.list(res)
}
