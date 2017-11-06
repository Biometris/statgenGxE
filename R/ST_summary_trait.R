#' Summary statistics of the trait
#'
#' This function is to calculate summary statistics of all traits.
#'
#' @param data A string path where the data list is saved.
#' @param trait A string (vector) specifying trait name(s).
#' @param nVals A logical value indicating if the number of values is included.
#' @param nObs A logical value indicating if the number of observations is included.
#' @param nMiss A logical value indicating if the number of missing values is included.
#' @param mean A logical value indicating if the mean is calculated.
#' @param median A logical value indicating if median is calculated.
#' @param min A logical value indicating if the minimum is calculated.
#' @param max A logical value indicating if the maximum is calculated.
#' @param range A logical value indicating if the range (maximum - minimum) is calculated.
#' @param lowerQ A logical value indicating if the lower (25\%) quantile is calculated.
#' @param upperQ A logical value indicating if the upper (75\%) quantile is calculated.
#' @param sd A logical value indicating if the standard deviation is calculated.
#' @param seMean A logical value indicating if the standard error of mean is calculated.
#' @param var A logical value indicating if the variance is calculated.
#' @param seVar A logical value indicating if the standard error of variance is calculated.
#' @param CV A logical value indicating if the coefficient of variation is calculated.
#' @param sum A logical value indicating if the sum is calculated.
#' @param sumSq A logical value indicating if the sum of squares is calculated.
#' @param uncorSumSq A logical value indicating if the uncorrected sum of squares is calculated.
#' @param skew A logical value indicating if the skewness is calculated.
#' @param seSkew A logical value indicating if the standard error of skewness is calculated.
#' @param kurt A logical value indicating if the kurtosis is calculated.
#' @param seKurt A logical value indicating if the standard error of kurtosis is calculated.
#' @param all Logical, if \code{all=TRUE}, all the statistics will be calculated.
#' @param printTable Logical, if \code{printTable=TRUE}, summary statistics will be shown on screen.
#' @return A data frame of summary statistics.
#' @examples
#' mydat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames=c("Env","Genotype","Rep","Subblock","Row","Column"),
#'                      traitNames="yield", env ="Env", rowSelect="HEAT06",
#'                      colSelect=c("Env","Genotype","Rep","Row","Column","yield"))
#' ST.summary.trait(data=mydat, trait="yield")
#'
#' @export
ST.summary.trait <- function(data,
                             trait, #env, envSel,
                             nVals = FALSE,
                             nObs = TRUE,
                             nMiss = TRUE,
                             mean = TRUE,
                             median = TRUE,
                             min = TRUE,
                             max = TRUE,
                             range = FALSE,
                             lowerQ = TRUE,
                             upperQ = TRUE,
                             sd = FALSE,
                             seMean = FALSE,
                             var = TRUE,
                             seVar = FALSE,
                             CV = FALSE,
                             sum = FALSE,
                             sumSq = FALSE,
                             uncorSumSq = FALSE,
                             skew = FALSE,
                             seSkew = FALSE,
                             kurt = FALSE,
                             seKurt = FALSE,
                             all = FALSE,
                             printTable = TRUE) {
  if (all) {
    nVals <- nObs <- nMiss <- mean <- median <- min <- max <-
      range <- lowerQ <- upperQ <- sd <- seMean <- var <- seVar <- CV <- sum <- sumSq <-
      uncorSumSq <- skew <- seSkew <- kurt <- seKurt <- TRUE
  }
  #Data = data[which(data[[env]]==envSel),]
  #a matrix to store the values
  stats <- matrix(nrow = 22, ncol = length(trait))
  colnames(stats) <- trait
  rownames(stats) <- c("Number of values","Number of observations","Number of missing values",
                       "mean","median", "min","max","range","Lower quartile","Upper quartile",
                       "Standard deviation", "Standard error of mean","Variance",
                       "Standard error of variance", "Coefficient of variation",
                       "sum of values", "sum of squares", "Uncorrected sum of squares",
                       "Skewness", "Standard Error of Skewness", "Kurtosis",
                       "Standard Error of Kurtosis")
  for (i in 1:length(trait)) {
    if (nVals) stats[1, i] <- length(data[, trait[i]])
    if (nObs) stats[2, i] <- length(na.omit(data[, trait[i]]))
    if (nMiss) stats[3, i] <- sum(is.na(data[, trait[i]]))
    if (mean) stats[4, i] <- mean(data[, trait[i]], na.rm = TRUE)
    if (median) stats[5, i] <- median(data[, trait[i]],na.rm = TRUE)
    if (min) stats[6, i] <- min(data[, trait[i]],na.rm = TRUE)
    if (max) stats[7, i] <- max(data[, trait[i]],na.rm = TRUE)
    if (range) stats[8, i] <- stats[7, i] - stats[6, i]
    if (lowerQ) stats[9, i] <- quantile(data[, trait[i]], prob = .25, na.rm = TRUE)
    if (upperQ) stats[10, i] <- quantile(data[,trait[i]], prob = .75, na.rm = TRUE)
    if (sd) stats[11, i] <- sd(data[, trait[i]], na.rm = TRUE)
    if (seMean) stats[12, i] <- sd(data[, trait[i]], na.rm = TRUE) /
        sqrt(length(na.omit(data[, trait[i]])))
    if (var) stats[13, i] <- var(data[, trait[i]], na.rm = TRUE)
    if (seVar) stats[14, i] <- seVar(data[, trait[i]], na.rm = TRUE)
    if (CV) stats[15, i] <- 100 * sd(data[, trait[i]], na.rm = TRUE) /
        mean(data[, trait[i]], na.rm = TRUE)
    if (sum) stats[16, i] <- sum(data[, trait[i]], na.rm = TRUE)
    if (sumSq) stats[17, i] <- sum((na.omit(data[, trait[i]]) -
                                      mean(data[, trait[i]], na.rm = TRUE))^2)
    if (uncorSumSq) stats[18, i] <- sum(data[, trait[i]] ^ 2, na.rm = TRUE)
    if (skew) stats[19, i] <- skewness(data[, trait[i]], na.rm = TRUE)
    if (seSkew) stats[20, i] <- seSkewness(length(na.omit(data[, trait[i]])))
    if (kurt) stats[21, i] <- kurtosis(data[, trait[i]], na.rm = TRUE)
    if (seKurt) stats[22, i] <- seKurtosis(length(na.omit(data[, trait[i]])))
    if (printTable){
      cat("\n")
      if(nVals)
        cat(paste0("             Number of values = ", stats[1, i], "\n"))
      if(nObs)
        cat(paste0("       Number of observations = ", stats[2, i], "\n"))
      if(nMiss)
        cat(paste0("     Number of missing values = ", stats[3,i ], "\n"))
      if(mean)
        cat(paste0("                         mean = ", round(stats[4, i], 2), "\n"))
      if(median)
        cat(paste0("                       median = ", round(stats[5, i], 2), "\n"))
      if(min)
        cat(paste0("                          min = ", round(stats[6, i], 2), "\n"))
      if(max)
        cat(paste0("                          max = ", round(stats[7, i], 2), "\n"))
      if(range)
        cat(paste0("               range(max-min) = ", round(stats[8, i], 2), "\n"))
      if(lowerQ)
        cat(paste0("               Lower quartile = ", round(stats[9, i], 2), "\n"))
      if(upperQ)
        cat(paste0("               Upper quartile = ", round(stats[10, i], 2), "\n"))
      if(sd)
        cat(paste0("           Standard deviation = ", round(stats[11, i], 3), "\n"))
      if(seMean)
        cat(paste0("       Standard error of mean = ", round(stats[12, i], 3), "\n"))
      if(var)
        cat(paste0("                     Variance = ", round(stats[13, i], 3), "\n"))
      if(seVar)
        cat(paste0("   Standard error of variance = ", round(stats[14, i], 3), "\n"))
      if(CV)
        cat(paste0("     Coefficient of variation = ", round(stats[15, i], 3), "\n"))
      if(sum)
        cat(paste0("                sum of values = ", round(stats[16, i], 2), "\n"))
      if(sumSq)
        cat(paste0("               sum of Squares = ", round(stats[17, i], 2), "\n"))
      if(uncorSumSq)
        cat(paste0("   Uncorrected sum of squares = ", round(stats[18, i], 2), "\n"))
      if(skew)
        cat(paste0("                     Skewness = ", round(stats[19, i], 3), "\n"))
      if(seSkew)
        cat(paste0("   Standard Error of Skewness = ", round(stats[20, i], 3), "\n"))
      if(kurt)
        cat(paste0("                     Kurtosis = ", round(stats[21, i], 3), "\n"))
      if(seKurt)
        cat(paste0("   Standard Error of Kurtosis = ", round(stats[22, i], 3), "\n"))
      cat("\n")
    }
  }
  invisible(stats)
}
