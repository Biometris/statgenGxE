#' Summary statistics of the trait
#'
#' This function is to calculate summary statistics of all traits.
#'
#' @param data A string path where the data list is saved.
#' @param trait A string (vector) specifying trait name(s).
#' @param NVals A logical value indicating if the number of values is included.
#' @param NObs A logical value indicating if the number of observations is included.
#' @param NMiss A logical value indicating if the number of missing values is included.
#' @param Mean A logical value indicating if the mean is calculated.
#' @param Median A logical value indicating if median is calculated.
#' @param Min A logical value indicating if the minimum is calculated.
#' @param Max A logical value indicating if the maximum is calculated.
#' @param Range A logical value indicating if the range (maximum - minimum) is calculated.
#' @param LowerQ A logical value indicating if the lower (25\%) quantile is calculated.
#' @param UpperQ A logical value indicating if the upper (75\%) quantile is calculated.
#' @param SD A logical value indicating if the standard deviation is calculated.
#' @param SEMean A logical value indicating if the standard error of mean is calculated.
#' @param Var A logical value indicating if the variance is calculated.
#' @param SEVar A logical value indicating if the standard error of variance is calculated.
#' @param CV A logical value indicating if the coefficient of variation is calculated.
#' @param Sum A logical value indicating if the sum is calculated.
#' @param SumSq A logical value indicating if the sum of squares is calculated.
#' @param UncorSumSq A logical value indicating if the uncorrected sum of squares is calculated.
#' @param Skew A logical value indicating if the skewness is calculated.
#' @param SESkew A logical value indicating if the standard error of skewness is calculated.
#' @param Kurt A logical value indicating if the kurtosis is calculated.
#' @param SEKurt A logical value indicating if the standard error of kurtosis is calculated.
#' @param all Logical, if \code{all=TRUE}, all the statistics will be calculated.
#' @param printTable Logical, if \code{printTable=TRUE}, summary statistics will be shown on screen.
#' @return A data frame of summary statistics.
#' @examples
#' mydat <- ST.read.csv(file.path(path.package("RAP"),"SB_yield.csv"),
#'                      factor.names=c("Env","Genotype","Rep","Subblock","Row","Column"),
#'                      trait.names="yield", env ="Env", rowSelect="HEAT06",
#'                      colSelect=c("Env","Genotype","Rep","Row","Column","yield"))
#' ST.summary.trait(data=mydat, trait="yield")
#' 
#' @export
ST.summary.trait = function(data, trait, #env, envSel,
    NVals = FALSE, NObs = TRUE, NMiss = TRUE,
    Mean = TRUE, Median = TRUE, Min = TRUE, Max = TRUE,
    Range = FALSE, LowerQ = TRUE, UpperQ = TRUE,
    SD = FALSE, SEMean = FALSE, Var = TRUE, SEVar = FALSE,
    CV = FALSE, Sum = FALSE, SumSq = FALSE, UncorSumSq = FALSE,
    Skew = FALSE, SESkew = FALSE, Kurt = FALSE, SEKurt = FALSE,
    all = FALSE, printTable = TRUE){

    if (all){
      NVals=NObs=NMiss=Mean=Median=Min=Max=
      Range=LowerQ=UpperQ=SD=SEMean=Var=SEVar=CV=Sum=SumSq=
      UncorSumSq=Skew=SESkew=Kurt=SEKurt= TRUE
    }

    #Data = data[which(data[[env]]==envSel),]
    #a matrix to store the values
    stats=matrix(NA,22,length(trait))
    colnames(stats)=trait
    rownames(stats)=c("Number of values","Number of observations","Number of missing values","Mean","Median",
                    "Min","Max","Range","Lower quartile","Upper quartile","Standard deviation", "Standard error of mean","Variance",
                    "Standard error of variance", "Coefficient of variation", "Sum of values", "Sum of squares",
                    "Uncorrected sum of squares", "Skewness", "Standard Error of Skewness", "Kurtosis", "Standard Error of Kurtosis")
    for (i in 1:length(trait)){
        if(NVals) stats[1,i] = length(data[,trait[i]])
        if(NObs) stats[2,i] = length(na.omit(data[,trait[i]]))
        if(NMiss) stats[3,i] = sum(is.na(data[,trait[i]]))
        if(Mean) stats[4,i] = mean(data[,trait[i]],na.rm=T)
        if(Median) stats[5,i] = median(data[,trait[i]],na.rm=T)
        if(Min) stats[6,i] = min(data[,trait[i]],na.rm=T)
        if(Max) stats[7,i] = max(data[,trait[i]],na.rm=T)
        if(Range) stats[8,i] = stats[7,i] - stats[6,i]
        if(LowerQ) stats[9,i]=quantile(data[,trait[i]],prob=.25,na.rm=T)
        if(UpperQ) stats[10,i]=quantile(data[,trait[i]],prob=.75,na.rm=T)
        if(SD) stats[11,i] = sd(data[,trait[i]],na.rm=T)
        if(SEMean) stats[12,i] = sd(data[,trait[i]],na.rm=T)/sqrt(length(na.omit(data[,trait[i]])))
        if(Var) stats[13,i] = var(data[,trait[i]],na.rm=T)
        if(SEVar) stats[14,i] = se.var(data[,trait[i]],na.rm=T)
        if(CV) stats[15,i] = 100*sd(data[,trait[i]],na.rm=T)/mean(data[,trait[i]],na.rm=T)
        if(Sum) stats[16,i] = sum(data[,trait[i]],na.rm=T)
        if(SumSq) stats[17,i] = sum((na.omit(data[,trait[i]])-mean(data[,trait[i]],na.rm=T))^2)
        if(UncorSumSq) stats[18,i] = sum(data[,trait[i]]^2,na.rm=T)
        if(Skew) stats[19,i] = skewness(data[,trait[i]],na.rm=T)
        if(SESkew) stats[20,i] = se.skewness(length(na.omit(data[,trait[i]])))
        if(Kurt) stats[21,i] = kurtosis(data[,trait[i]],na.rm=T)
        if(SEKurt) stats[22,i] = se.kurtosis(length(na.omit(data[,trait[i]])))

        if (printTable){
          cat("\n")
          if(NVals)
            cat(paste("             Number of values = ",stats[1,i],"\n",sep=""))
          if(NObs)
            cat(paste("       Number of observations = ",stats[2,i],"\n",sep=""))
          if(NMiss)
            cat(paste("     Number of missing values = ",stats[3,i],"\n",sep=""))
          if(Mean)
            cat(paste("                         Mean = ",round(stats[4,i],2),"\n",sep=""))
          if(Median)
            cat(paste("                       Median = ",round(stats[5,i],2),"\n",sep=""))
          if(Min)
            cat(paste("                          Min = ",round(stats[6,i],2),"\n",sep=""))
          if(Max)
            cat(paste("                          Max = ",round(stats[7,i],2),"\n",sep=""))
          if(Range)
            cat(paste("               Range(Max-Min) = ",round(stats[8,i],2),"\n",sep=""))
          if(LowerQ)
            cat(paste("               Lower quartile = ",round(stats[9,i],2),"\n",sep=""))
          if(UpperQ)
            cat(paste("               Upper quartile = ",round(stats[10,i],2),"\n",sep=""))
          if(SD)
            cat(paste("           Standard deviation = ",round(stats[11,i],3),"\n",sep=""))
          if(SEMean)
            cat(paste("       Standard error of mean = ",round(stats[12,i],3),"\n",sep=""))
          if(Var)
            cat(paste("                     Variance = ",round(stats[13,i],3),"\n",sep=""))
          if(SEVar)
            cat(paste("   Standard error of variance = ",round(stats[14,i],3),"\n",sep=""))
          if(CV)
            cat(paste("     Coefficient of variation = ",round(stats[15,i],3),"\n",sep=""))
          if(Sum)
            cat(paste("                Sum of values = ",round(stats[16,i],2),"\n",sep=""))
          if(SumSq)
            cat(paste("               Sum of Squares = ",round(stats[17,i],2),"\n",sep=""))
          if(UncorSumSq)
            cat(paste("   Uncorrected sum of squares = ",round(stats[18,i],2),"\n",sep=""))
          if(Skew)
            cat(paste("                     Skewness = ",round(stats[19,i],3),"\n",sep=""))
          if(SESkew)
            cat(paste("   Standard Error of Skewness = ",round(stats[20,i],3),"\n",sep=""))
          if(Kurt)
            cat(paste("                     Kurtosis = ",round(stats[21,i],3),"\n",sep=""))
          if(SEKurt)
            cat(paste("   Standard Error of Kurtosis = ",round(stats[22,i],3),"\n",sep=""))
          cat("\n")
        }
    }
    invisible(stats)
}
