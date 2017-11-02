#' AMMI analysis
#'
#' This function fits a model which involves the Additive Main effects (i.e. genotype and environment) along
#' with the Multiplicative Interaction effects of principal component analysis (PCA).
#'
#' @param Y A data frame object.
#' @param trait A character string specifying a trait column of the data.
#' @param genotype A character string specifying a genotype column of the data.
#' @param env A character string specifying an environment column of the data.
#' @param nPC An integer specifying the number of principal components used.
#' as the multiplicative term of genotype-by-environment. \code{nPC=2} by default.
#' @param center  A logical value indicating whether the variables
#'  should be shifted to be zero centered.
#' @param scale A logical value indicating whether the variables should
#'  be scaled to have unit variance before the analysis takes
#'  place. The default is \code{FALSE}.
#' @param AMMI2plot A logical value specifying whether an AMMI2 biplot is drawn.
#' @param scale.AMMI2 The variables are scaled by \code{lambda ^ scale} and the observations are scaled by \code{lambda ^ (1-scale)}
#' where \code{lambda} are the singular values as computed by \code{\link[stats]{princomp}}. Normally \code{0 <= scale <= 1},
#' and a warning will be issued if the specified scale is outside this range.
#' @param AMMI1plot A logical value determining whether the AMMI1 biplot (genotypes and environments means vs PC1) is drawn.
#' @param scale.AMMI1 as same as described in \code{scale.AMMI2}.
#' @return A list of three object, a data frame object of environment scores, a data frame object of genotype scores and
#' an object of class \code{\link[stats]{anova}},
#' and a matrix of the fitted values from the AMMI model.
#' @examples
#' mydat <- GE.read.csv(file.path(path.package("RAP"),"F2maize_pheno.csv"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' names(mydat)=c("env", "genotype","yld")
#' GE.AMMI(Y=mydat, trait="yld", genotype="genotype", env="env", nPC = 2,
#'         center = TRUE, scale = FALSE, AMMI2plot = TRUE, scale.AMMI2=1)
#'
#' @import stats graphics grDevices
#' @export

GE.AMMI <- function(Y, trait, genotype, env, nPC=2, center=T, scale=F,
  AMMI2plot=T, scale.AMMI2=1,
  AMMI1plot=F,   scale.AMMI1=1){

  #drop factor levels
  Y[[genotype]] <- droplevels(Y[[genotype]])
  Y[[env]] <- droplevels(Y[[env]])

  ngeno <- nlevels(Y[[genotype]])
  nenv  <- nlevels(Y[[env]])
  ntrait <- nrow(Y)

  # requre number of environments >=3
  if (nenv < 3)
    stop("Requires number of environments greater and equal than 3 for running the AMMI model")

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

  # Descriptive statistics
  env.mean <- tapply(Y[[trait]], Y[[env]], mean)
  geno.mean <- tapply(Y[[trait]], Y[[genotype]], mean)
  overall.mean <- mean(Y[[trait]])

  # Fit the linear model
  model <- lm(as.formula(paste(trait,"~",genotype, "+", env)), data = Y)

  # calculate residuals & fitted values of the linear model
  X <- tapply(resid(model), Y[,c(genotype,env)], identity)
  fittedvals <- tapply(fitted(model), Y[,c(genotype,env)], identity)

  X <- as.matrix(X)
  gnames <- rownames(X)
  if (is.null(gnames)) rownames(X) <- rownames(fittedvals) <- gnames

  # Use R in-built prcomp
  pca <- prcomp(x=X, retx = TRUE,center=center, scale.=scale)
  cump2 <- summary(pca)$importance[3,2]
  prop.pc1 <- summary(pca)$importance[2,1]
  prop.pc2 <- summary(pca)$importance[2,2]
  loadings <- pca$rotation
  scores <- pca$x

  if (AMMI2plot){
    if (scale.AMMI2==1){
      info <- "environment scaling"
    }else if (scale.AMMI2==0){
      info <- "genotype scaling"
    }else if (scale.AMMI2==0.5){
      info <- "symmetric scaling"
    }else{
      info <- paste0(round(cump2*100,1), "%")
    }
    old.par <- par(xpd=NA)
    biplot(pca, scale=scale.AMMI2, col=c("orange3", "navyblue"),
      main=paste("AMMI2 biplot for ",trait," (",info,")",sep=""),
      xlab=paste("PC1 (",round(prop.pc1*100,1),"%)",sep=""),
      ylab=paste("PC2 (",round(prop.pc2*100,1),"%)",sep=""))
    par(old.par)
  }

  if (AMMI1plot){
    # calculate lambda scale
    lam <- pca$sdev[1]
    if(is.null(n <- pca$n.obs)) n <- 1
    lam <- lam * sqrt(n)
    if(scale.AMMI1 < 0 || scale.AMMI1 > 1) warning("'scale' is outside [0, 1]")
    if(scale.AMMI1 != 0) lam <- lam^scale.AMMI1 else lam <- 1

    #biplot(x=cbind(env.mean, loadings[,1]*lam), y=cbind(geno.mean, scores[,1]/lam), var.axes=F,
    #col=c("red", "blue"), xlab="Main effects", ylab="PC1")
    dev.new()
    old.par <- par(xpd=NA)
    plot(1, type = 'n', xlim = range(c(env.mean, geno.mean)), ylim = range(c(loadings[,1]*lam, scores[,1]/lam)),
      xlab = "Main Effects", ylab =paste("PC1 (",round(prop.pc1*100,1),"%)",sep=""),
      main=paste("AMMI1 biplot for ",trait, sep=""))
    points(env.mean, loadings[,1]*lam, "n", col = "navyblue", lwd = 5)
    text(env.mean, loadings[,1]*lam, labels = row.names(env.mean), adj = c(0.5, 0.5), col = "navyblue")
    points(geno.mean, scores[,1]/lam, "n", col = "orange3", lwd = 5)
    text(geno.mean, scores[,1]/lam, labels = row.names(geno.mean), adj = c(0.5, 0.5), col = "orange3")
    par(old.par)
    abline(h = 0, v = overall.mean, lty = 5)
  }

  # calculating the AMMI-estimates per genotype per environment
  mterms <- matrix(0, nrow = ngeno, ncol = nenv)
  for (ii in 1:nPC)
    mterms <- mterms +  outer(scores[,ii], loadings[,ii])
  fitted <- fittedvals + mterms

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
  #res$PCA = pca
  res$Environment_scores <- pca$rotation
  res$Genotype_scores <- pca$x
  res$ANOVA <- a0
  res$fitted <- fitted
  as.list(res)
}







## explicitly use singular value decomposition to compute pca
#  s <- svd(X)
#  U <-s$u
#  L <- s$d
#  V <- s$v
#  LL <- sqrt(diag(L))
#  loadings <- V %*% LL
#  scores <- U %*% LL
#  res <- new.env()
#  res$loadings <- loadings
#  res$scores <- scores
#  as.list(res)
