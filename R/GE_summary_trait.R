#' Summary statistics for a trait
#'
#' This function calculates summary statistics for a trait within and between a set of environments.
#'
#' @param data A data frame object consisting of Genotype-by-Environment data.
#' @param trait A string specifying the column name of a trait column found in \code{data}.
#' @param env A string specifying the column name for the environments found in \code{data}.
#' @param genotype A string specifying the column name of the genotypes found in \code{data}.
#' @param envSelect A character vector of unique names to be selected from the enviroment ID(s). If \code{NA}, all the rows are included.
#' @param NVals A logical value indicating if number of values is included.
#' @param NNVals A logical value indicating if number of non-missing values is included.
#' @param NMiss A logical value indicating if number of missing values is included.
#' @param Mean A logical value indicating if mean is calculated.
#' @param Median A logical value indicating if median is calculated.
#' @param Min A logical value indicating if minimum is calculated.
#' @param Max A logical value indicating if maximum is calculated.
#' @param Range A logical value indicating if range (maximum - minimum) is calculated.
#' @param LowerQ A logical value indicating if (25\%) lower quantile is calculated.
#' @param UpperQ A logical value indicating if (75\%) upper quantile is calculated.
#' @param SD A logical value indicating if standard deviation is calculated.
#' @param SEMean A logical value indicating if standard error of mean is calculated.
#' @param Var A logical value indicating if variance is calculated.
#' @param SEVar A logical value indicating if standard error of variance is calculated.
#' @param CV A logical value indicating if coefficient of variation is calculated.
#' @param Sum A logical value indicating if sum is calculated.
#' @param SumSq A logical value indicating if sum of squares is calculated.
#' @param UncorSumSq A logical value indicating if uncorrected sum of squares is calculated.
#' @param Skew A logical value indicating if skewness is calculated.
#' @param SESkew A logical value indicating if standard error of skewness is calculated.
#' @param Kurt A logical value indicating if kurtosis is calculated.
#' @param SEKurt A logical value indicating if standard error of kurtosis is calculated.
#' @param Boxplot A logical value indicating if boxplot is drawn.
#' @param Histogram A logical value indicating if histogram is drawn.
#' @param Corplot A logical value indicating if a heatmap plot is drawn from the correlation matrix.
#' @param Scatterplot A logical value indicating if a scatter plot is drawn.
#' @param all A logical, if \code{all=TRUE}, all the statistics will be calculated.
#' @examples
#' mydat <- GE.read.csv(file.path(path.package("RAP"),"F2maize_pheno.csv"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' names(mydat)=c("env", "genotype","yld") 
#' GE.summary.trait(data=mydat, trait="yld", env="env", genotype="genotype", 
#'                  envSelect=c("HN96b", "IS92a", "IS94a", "LN96a",
#'                   "LN96b", "NS92a", "SS92a", "SS94a"))
#'
#' @export
GE.summary.trait = function(data, trait, env, genotype, envSelect=NA,
    NVals = F, NNVals = T, NMiss = T,
    Mean = T, Median = T, Min = T, Max = T,
    Range = F, LowerQ = T, UpperQ = T,
    SD = F, SEMean = F, Var = T, SEVar =F,
    CV = F, Sum = F, SumSq = F, UncorSumSq = F,
    Skew = F, SESkew = F, Kurt = F, SEKurt = F,
    Boxplot = F, Histogram = F,
    Corplot = T, Scatterplot= T, all = F){
    
#    load(datname)
#    data = GE.data

    if (!all(is.na(envSelect))){
        data=data[data[[env]]%in%envSelect,]
        data=droplevels(data)
        # check if the number of rows is correct
        if(length(unique(data[[env]]))!= length(envSelect))
          warning("CHECK: Either some selected row(s) are not read in\n or the selected row name(s) do not match the enviroment ID(s)!")
    }else{
        envSelect <- levels(data[[env]])
    }
    
    if (all){
      NVals=NNVals=NMiss=Mean=Median=Min=Max=
      Range=LowerQ=UpperQ=SD=SEMean=Var=SEVar=CV=Sum=SumSq=
      UncorSumSq=Skew=SESkew=Kurt=SEKurt=CorMat=Boxplot=Histogram=Corplot=Scatterplot= TRUE
    }
    
    #a matrix to store the values
    stats=matrix(NA,length(envSelect), 22)
    rownames(stats)=envSelect
    colnames(stats)=c("No. of values", "No. of non-missing values",
    "No. of missing values", "Mean", "Variance",
    "Standard deviation", "Min", "Max",
    "Range (max-min)", "Median", "Lower quartile",
    "Upper quartile", "Sum of values", "Standard error of mean",
    "Standard error of variance", "Coefficient of variation", "Sum of squares",
    "Uncorrected sum of squares", "Skewness", "Standard Error of Skewness",
    "Kurtosis", "Standard Error of Kurtosis")
    for (i in 1:length(envSelect)){
        y <- data[which(data[,env]==envSelect[i]),trait]
        if(NVals) stats[i,1] <- length(y)
        if(NNVals) stats[i,2] <- length(y)- sum(is.na(y))
        if(NMiss) stats[i,3] <- sum(is.na(y))
        if(Mean) stats[i,4] <- mean(y,na.rm=T)
        if(Var) stats[i,5] <- var(y,na.rm=T)
        if(SD) stats[i,6] <- sd(y,na.rm=T)
        if(Min) stats[i,7] <- min(y,na.rm=T)
        if(Max) stats[i,8] <- max(y,na.rm=T)
        if(Range) stats[i,9] <- max(y,na.rm=T)- min(y,na.rm=T)
        if(Median) stats[i,10] <- median(y,na.rm=T)
        if(LowerQ) stats[i,11] <- quantile(y,prob=.25,na.rm=T)
        if(UpperQ) stats[i,12] <- quantile(y,prob=.75,na.rm=T)
        if(Sum) stats[i,13] <- sum(y,na.rm=T)
        if(SEMean) stats[i,14] <- sd(y,na.rm=T)/sqrt(length(y)- sum(is.na(y)))
        if(SEVar)  stats[i,15] <- se.var(y,na.rm=T)
        if(CV) stats[i,16] <- 100*sd(y,na.rm=T)/mean(y,na.rm=T)
        if(SumSq) stats[i,17] <- sum((na.omit(y)-mean(y,na.rm=T))^2)
        if(UncorSumSq) stats[i,18] <- sum(y^2,na.rm=T)
        if(Skew) stats[i,19]  <- skewness(y,na.rm=T)
        if(SESkew) stats[i,20] <- se.skewness(length(na.omit(y)))
        if(Kurt) stats[i,21] <- kurtosis(y,na.rm=T)
        if(SEKurt) stats[i,22] <- se.kurtosis(length(na.omit(y)))

    }

    if(Boxplot){
      f0 <- as.formula(paste(trait,'~',env))
      boxplot(formula = f0 , data = data,
          main = "Boxplot by Enviroment",
          xlab = "Enviroment", ylab = trait)
    }

    if(Histogram){
      f0 <- as.formula(paste('~',trait,"|",env))
      dev.new()
      print(histogram(x = f0 , data = data,
          main = "Histograms by Enviroment",
          xlab = trait))
    }
             
    if (Scatterplot){
      ## put (absolute) correlations on the upper panels,
      ## with size proportional to the correlations.
      panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
      {
          usr <- par("usr"); on.exit(par(usr))
          par(usr = c(0, 1, 0, 1))
          r <- abs(cor(x, y))
          txt <- format(c(r, 0.123456789), digits = digits)[1]
          txt <- paste(prefix, txt, sep = "")
          if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
          text(0.5, 0.5, txt, cex = cex.cor * r)
      }
      ## put histograms on the diagonal
      panel.hist <- function(x, ...)
      {
          usr <- par("usr"); on.exit(par(usr))
          par(usr = c(usr[1:2], 0, 1.5) )
          h <- hist(x, plot = FALSE)
          breaks <- h$breaks; nB <- length(breaks)
          y <- h$counts; y <- y/max(y)
          rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
      }
      X <- tapply(data[,trait], data[, c(genotype,env)], identity)
      if (any(is.na(X))){
        X <- RAP.multmissing(X, maxcycle = 10, na.strings=NA)
      }
      dev.new()
      pairs(X,upper.panel=panel.cor,diag.panel=panel.hist,
      main=paste("Scatterplot matrix and correlations by enviroment:",trait))
    }
    if (Corplot){
      if (is.data.frame(data)&& !missing(env)){
        if (all(c(trait,env)%in%names(data))){
          # create a correlation matrix
          mydat <- tapply(data[,trait], data[, c(genotype,env)], identity)
          if (any(is.na(mydat))){
            mydat <- RAP.multmissing(mydat, maxcycle = 10, na.strings=NA)
          }
          if (is.matrix(mydat)){
            cormat <- cor(mydat)
            cormat[upper.tri(cormat, diag = F)] <- NA
            cormat
          }else{
            stop("Numbers of observations for all environments are not equal.\n")
          } 
        }
      }
      dev.new()
      print(levelplot(cormat[nrow(cormat):1,], xlab="", ylab="",
      main=paste("Correlation of environments for", trait)))
      
    }
    na.col <- apply(stats, 2, function(x) all(is.na(x)))
    stats2 <- stats[,!na.col]
    as.data.frame(stats2)
}
