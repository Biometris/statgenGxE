#' GGE biplot
#'
#' This function draws GGE biplots which are useful for assessing the performance of genotypes in different environments.
#'
#' @param Y A data frame object.
#' @param trait A character string specifying a trait column of the data.
#' @param genotype A character string specifying a genotype column of the data.
#' @param env A character string specifying an environment column of the data.
#' @param nPC An integer specifying the number of principal components used
#' as the multiplicative term of genotype by environment. \code{nPC=2} by default.
#' @param center  A logical value indicating whether the variables
#'  should be shifted to be zero centered.
#' @param scale A logical value indicating whether the variables should
#'  be scaled to have unit variance before the analysis takes
#'  place. The default is \code{FALSE}.
#' @param GGE2plot A logical value determining whether the GGE biplot is drawn.
#' @param scale.biplot The variables are scaled by \code{lambda ^ scale} and the observations are scaled by \code{lambda ^ (1-scale)}
#' where \code{lambda} are the singular values as computed by \code{\link[stats]{princomp}}. Normally \code{0 <= scale <= 1},
#' and a warning will be issued if the specified scale is outside this range.
#' @return A list of two objects, a list with class \code{\link[stats]{prcomp}} and an object of class \code{\link[stats]{anova}}.
#' @examples
#' mydat <- GE.read.csv(file.path(path.package("RAP"),"F2maize_pheno.csv"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' names(mydat)=c("env", "genotype","yld") 
#' GE.GGE(Y=mydat, trait="yld", genotype="genotype", env="env", scale.biplot=1)
#'
#' @export
GE.GGE <- function(Y, trait, genotype, env, nPC=2, center=T, scale=F, GGE2plot=T, scale.biplot=0){

  #drop factor levels
  Y[[genotype]] <- droplevels(Y[[genotype]])
  Y[[env]] <- droplevels(Y[[env]])

  ngeno <- nlevels(Y[[genotype]])
  nenv  <- nlevels(Y[[env]])
  ntrait <- nrow(Y)
      
  #check if the supplied data contains the genotype by environment means
  if (ntrait != ngeno * nenv)
    stop("Only allows the genotype by environment means, \ni.e., one trait value per genotype per enviroment")

  if (any(is.na(Y[[trait]]))){
    y0 <- tapply(Y[[trait]], Y[,c(genotype,env)], identity)
    yindex <- tapply(1:ntrait, Y[,c(genotype,env)], identity)
    na_yes_no <- is.na(y0)
    # imputaion
    y1 <- RAP.multmissing(y0, maxcycle = 10, na.strings=NA)
    replace.val <- y1[na_yes_no]
    Y[yindex[na_yes_no], trait] <- replace.val
  }
  
  # Fit the linear model
  model <- lm(as.formula(paste(trait,"~", env)), data = Y)

  # Use R in-built prcomp
  X <- tapply(resid(model), Y[,c(genotype,env)], identity)
  X <- as.matrix(X)

  pca <- prcomp(x=X,center=center, scale.=scale)
  if (GGE2plot){
    cump2 <- summary(pca)$importance[3,2]
    prop.pc1 <- summary(pca)$importance[2,1]
    prop.pc2 <- summary(pca)$importance[2,2]
    if (scale.biplot==1){
      info <- "environment scaling"
    }else if (scale.biplot==0){
      info <- "genotype scaling"
    }else if (scale.biplot==0.5){
      info <- "symmetric scaling"
    }else{
      info <- paste0(round(cump2*100,1), "%")
    }
    old.par <- par(xpd=NA)
    biplot(pca, scale=scale.biplot, col=c("orange3", "navyblue"),
      main=paste("GGE2 biplot for ",trait," (", info,")",sep=""),
      xlab=paste("PC1 (",round(prop.pc1*100,1),"%)",sep=""),
      ylab=paste("PC2 (",round(prop.pc2*100,1),"%)",sep=""))
    par(old.par)
  }
  
  # ANOVA table for linear model
  a1 <- anova(model)
  tnames <- rownames(a1)
  rownames(a1)[which(tnames=="Residuals")] <- "Interactions"
  
  # Extend the existing ANOVA table
  addTbl <- matrix(NA, nPC+1, 5)
  colnames(addTbl) <- c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
  tnames <- paste("PC", 1:nPC, sep="")
  rownames(addTbl) <- c(tnames, "Residuals")
  
  # Add the df for PC scores and residuals
  ngeno <- nlevels(Y[,genotype])
  nenv  <- nlevels(Y[,env])
  dfPC <- ngeno+nenv-3-(2*(1:nPC-1))
  dfresid <- a1["Interactions", "Df"] - sum(dfPC)
  addTbl[,"Df"] <- c(dfPC, dfresid)
  
  # Add the sum of squaures for PC scores and residuals
  PCAVar <- pca$sdev^2
  totalVar <- sum(pca$sdev^2)
  propVar <- PCAVar/totalVar
  ssPC <- a1["Interactions", "Sum Sq"]*propVar[1:nPC]
  ssresid <- a1["Interactions", "Sum Sq"] - sum(ssPC)
  addTbl[,"Sum Sq"] <- c(ssPC, ssresid)
  
  # Add the mean squaures for PC scores and residuals
  addTbl[,"Mean Sq"] <- addTbl[,"Sum Sq"]/addTbl[,"Df"]
  whichinf <- is.infinite(addTbl[,"Mean Sq"])
  addTbl[,"Mean Sq"][whichinf] <-NA
    
  # Add the F-values for PC scores
  addTbl[-(nPC+1),"F value"] <- addTbl[-(nPC+1),"Mean Sq"]/addTbl["Residuals","Mean Sq"]
  
  # Add the p-value for PC scores
  addTbl[-(nPC+1),"Pr(>F)"] <- 1 - sapply(1:nPC, function(i) pf(q=addTbl[i,"F value"],
  df1=addTbl[i,"Df"], df2=addTbl["Residuals","Df"]))
  
  # ANOVA table for AMMI model
  a0 <- rbind(a1,as.data.frame(addTbl))
  
  res <- new.env()
  res$PCA = pca
  res$ANOVA = a0
  as.list(res)
}
