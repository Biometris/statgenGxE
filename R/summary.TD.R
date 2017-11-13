#' Summary statistics of the trait
#'
#' This function is to calculate summary statistics of all traits.
#'
#' @param object An object of class TD.
#' @param ... Further arguments - currently not used.
#' @param trait A string vector specifying trait name(s) to be summarised.
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
#' @param all a logical value indicating if all the statistics should be calculated.
#'
#' @return A data.frame containing the selected summary statistics.
#' @seealso \code{\link{createTD}}
#'
#' @examples
#' myDat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames = c("Env", "Genotype", "Rep", "Subblock"),
#'                      traitNames = "yield", env = "Env", rowSelect = "HEAT05",
#'                      colSelect = c("Env", "Genotype", "Rep", "Subblock", "yield"))
#' myTD <- createTD(data = myDat, genotype = "Genotype", env = "Env")
#' summary(object = myTD, trait = "yield")
#'
#' @export
summary.TD <- function(object,
                       ...,
                       trait,
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
                       all = FALSE) {
  if (all) {
    nVals <- nObs <- nMiss <- mean <- median <- min <- max <-
      range <- lowerQ <- upperQ <- sd <- seMean <- var <- seVar <- CV <- sum <- sumSq <-
      uncorSumSq <- skew <- seSkew <- kurt <- seKurt <- TRUE
  }
  ## Create a matrix to store the values
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
    if (nVals) stats[1, i] <- length(object[, trait[i]])
    if (nObs) stats[2, i] <- length(na.omit(object[, trait[i]]))
    if (nMiss) stats[3, i] <- sum(is.na(object[, trait[i]]))
    if (mean) stats[4, i] <- mean(object[, trait[i]], na.rm = TRUE)
    if (median) stats[5, i] <- median(object[, trait[i]],na.rm = TRUE)
    if (min) stats[6, i] <- min(object[, trait[i]],na.rm = TRUE)
    if (max) stats[7, i] <- max(object[, trait[i]],na.rm = TRUE)
    if (range) stats[8, i] <- stats[7, i] - stats[6, i]
    if (lowerQ) stats[9, i] <- quantile(object[, trait[i]], prob = .25, na.rm = TRUE)
    if (upperQ) stats[10, i] <- quantile(object[,trait[i]], prob = .75, na.rm = TRUE)
    if (sd) stats[11, i] <- sd(object[, trait[i]], na.rm = TRUE)
    if (seMean) stats[12, i] <- sd(object[, trait[i]], na.rm = TRUE) /
        sqrt(length(na.omit(object[, trait[i]])))
    if (var) stats[13, i] <- var(object[, trait[i]], na.rm = TRUE)
    if (seVar) stats[14, i] <- seVar(object[, trait[i]], na.rm = TRUE)
    if (CV) stats[15, i] <- 100 * sd(object[, trait[i]], na.rm = TRUE) /
        mean(object[, trait[i]], na.rm = TRUE)
    if (sum) stats[16, i] <- sum(object[, trait[i]], na.rm = TRUE)
    if (sumSq) stats[17, i] <- sum((na.omit(object[, trait[i]]) -
                                      mean(object[, trait[i]], na.rm = TRUE))^2)
    if (uncorSumSq) stats[18, i] <- sum(object[, trait[i]] ^ 2, na.rm = TRUE)
    if (skew) stats[19, i] <- skewness(object[, trait[i]], na.rm = TRUE)
    if (seSkew) stats[20, i] <- seSkewness(length(na.omit(object[, trait[i]])))
    if (kurt) stats[21, i] <- kurtosis(object[, trait[i]], na.rm = TRUE)
    if (seKurt) stats[22, i] <- seKurtosis(length(na.omit(object[, trait[i]])))
  }
  return(structure(stats,
                   class = c("summary.TD", "table")))
}

#' @export
print.summary.TD <- function(x, ...) {
  for (i in 1:ncol(x)) {
    cat("\nSummary statistics for", colnames(x)[i], "\n\n")
    if(!is.na(x[1, i]))
      cat(paste0("             Number of values = ", x[1, i], "\n"))
    if(!is.na(x[2, i]))
      cat(paste0("       Number of observations = ", x[2, i], "\n"))
    if(!is.na(x[3, i]))
      cat(paste0("     Number of missing values = ", x[3, i], "\n"))
    if(!is.na(x[4, i]))
      cat(paste0("                         mean = ", round(x[4, i], 2), "\n"))
    if(!is.na(x[5, i]))
      cat(paste0("                       median = ", round(x[5, i], 2), "\n"))
    if(!is.na(x[6, i]))
      cat(paste0("                          min = ", round(x[6, i], 2), "\n"))
    if(!is.na(x[7, i]))
      cat(paste0("                          max = ", round(x[7, i], 2), "\n"))
    if(!is.na(x[8, i]))
      cat(paste0("               range(max-min) = ", round(x[8, i], 2), "\n"))
    if(!is.na(x[9, i]))
      cat(paste0("               Lower quartile = ", round(x[9, i], 2), "\n"))
    if(!is.na(x[10, i]))
      cat(paste0("               Upper quartile = ", round(x[10, i], 2), "\n"))
    if(!is.na(x[11, i]))
      cat(paste0("           Standard deviation = ", round(x[11, i], 3), "\n"))
    if(!is.na(x[12, i]))
      cat(paste0("       Standard error of mean = ", round(x[12, i], 3), "\n"))
    if(!is.na(x[13, i]))
      cat(paste0("                     Variance = ", round(x[13, i], 3), "\n"))
    if(!is.na(x[14, i]))
      cat(paste0("   Standard error of variance = ", round(x[14, i], 3), "\n"))
    if(!is.na(x[15, i]))
      cat(paste0("     Coefficient of variation = ", round(x[15, i], 3), "\n"))
    if(!is.na(x[16, i]))
      cat(paste0("                sum of values = ", round(x[16, i], 2), "\n"))
    if(!is.na(x[17, i]))
      cat(paste0("               sum of Squares = ", round(x[17, i], 2), "\n"))
    if(!is.na(x[18, i]))
      cat(paste0("   Uncorrected sum of squares = ", round(x[18, i], 2), "\n"))
    if(!is.na(x[19, i]))
      cat(paste0("                     Skewness = ", round(x[19, i], 3), "\n"))
    if(!is.na(x[20, i]))
      cat(paste0("   Standard Error of Skewness = ", round(x[20, i], 3), "\n"))
    if(!is.na(x[21, i]))
      cat(paste0("                     Kurtosis = ", round(x[21, i], 3), "\n"))
    if(!is.na(x[22, i]))
      cat(paste0("   Standard Error of Kurtosis = ", round(x[22, i], 3), "\n"))
    cat("\n")
  }
}


