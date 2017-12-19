#' S3 class TD
#'
#' Function for creating objects of S3 class TD (Trial Data). The input data is checked
#' and columns are renamed to default column names for ease of further computations.
#' The columns for genotype, env, megaEnv, year, repId, subBlock, rowId, colId and
#' checkId are converted to factor columns, whereas rowCoordinates and colCoordinates
#' are converted to numerical columsn.
#'
#' \code{\link{print}} and \code{\link{summary}} methods are available.
#'
#' @param data a data.frame containing trial data with a least a column for
#' genotype.
#' @param genotype a character string indicating the column in \code{data} that
#' contains genotypes.
#' @param env an optional character string indicating the column in \code{data} that
#' contains environments. If \code{NULL} a default environment will be added.
#' @param megaEnv an optional character string indicating the column in \code{data} that
#' contains megaEnvironments as constructed by \code{\link{gxeMegaEnvironment}}.
#' @param year an optional character string indicating the column in \code{data} that
#' contains years.
#' @param repId an optional character string indicating the column in \code{data} that
#' contains replicates.
#' @param subBlock an optional character string indicating the column in \code{data} that
#' contains sub blocks.
#' @param rowId an optional character string indicating the column in \code{data} that
#' contains field rows.
#' @param colId an optional character string indicating the column in \code{data} that
#' contains field columns.
#' @param rowCoordinates an optional character string indicating the column in
#' \code{data} that contains the rowId coordinates used for fitting spatial models.
#' @param colCoordinates an optional character string indicating the column in
#' \code{data} that contains the column coordinates used for fitting spatial models.
#' @param checkId an optional character string indicating the column in \code{data} that
#' contains the check ID(s).
#' @param design an optional character string indicating the design of the trial. Accepted
#' values are "ibd" (incomplete-block design), "res.ibd" (resolvable incomplete-block design),
#' "rcbd" (randomized complete block design), "rowcol" (rowId-column design) and
#' "res.rowcol" (resolvable rowId-column design).
#' @param x an \code{R} object
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{summary.TD}}
#'
#' @name TD
NULL

#' @rdname TD
#' @export
createTD <- function(data,
                     genotype = "genotype",
                     env = NULL,
                     megaEnv = NULL,
                     year = NULL,
                     repId = NULL,
                     subBlock = NULL,
                     rowId = NULL,
                     colId = NULL,
                     rowCoordinates = NULL,
                     colCoordinates = NULL,
                     checkId = NULL,
                     design = NULL) {
  ## Checks.
  if (missing(data) || !is.data.frame(data)) {
    stop("data has to be a data.frame.\n")
  }
  cols <- colnames(data)
  if (is.null(genotype) || !is.character(genotype) || length(genotype) > 1 ||
      !genotype %in% cols) {
    stop("genotype has to be a column in data.\n")
  }
  for (param in c(env, megaEnv, year, repId, subBlock, rowId, colId,
                  rowCoordinates, colCoordinates, checkId)) {
    if (!is.null(param) && (!is.character(param) || length(param) > 1 || !param %in% cols)) {
      stop(paste(deparse(param), "has to be NULL or a column in data.\n"))
    }
  }
  if (!is.null(design) && (!is.character(design) || length(design) > 1 ||
                           !design %in% c("ibd", "res.ibd", "rcbd", "rowcol", "res.rowcol"))) {
    stop("design has to be NULL or one of ibd, res.ibd, rcbd, rowcol or res.rowcol.\n")
  }
  ## Rename columns.
  renameCols <- c("genotype", "env", "megaEnv", "year", "repId", "subBlock", "rowId",
                  "colId", "rowCoordinates", "colCoordinates", "checkId")
  ## First rename duplicate colums and add duplicated columns to data
  renameFrom <- as.character(sapply(X = renameCols, FUN = function(x) {
    get(x)
  }))
  ## Get duplicate columns
  dupCols <- which(duplicated(renameFrom) & renameFrom != "NULL")
  for (dupCol in dupCols) {
    ## Copy original column as extra column in data for each duplicate.
    tempName <- paste0(".temp", dupCol)
    data[tempName] <- data[, which(colnames(data) == renameFrom[dupCol])]
    ## Add new replacementname to cols and renameFrom.
    cols[length(cols) + 1] <- tempName
    renameFrom[dupCol] <- tempName
  }
  ## Rename columns.
  for (i in 1:length(renameCols)) {
    cols[cols == renameFrom[i]] <- renameCols[i]
  }
  colnames(data) <- cols
  ## Convert columns to factor if neccessary.
  factorCols <-  c("genotype", "env", "megaEnv", "year", "repId", "subBlock", "rowId",
                   "colId", "checkId")
  for (factorCol in factorCols) {
    if (factorCol %in% cols && !is.factor(data[, which(cols == factorCol)])) {
      data[, which(cols == factorCol)] <-
        as.factor(data[, which(cols == factorCol)])
    }
  }
  ## Convert columns to numeric if neccessary.
  numCols <- c("rowCoordinates", "colCoordinates")
  for (numCol in numCols) {
    if (numCol %in% cols && !is.numeric(data[, which(cols == numCol)])) {
      data[, which(cols == numCol)] <-
        as.numeric(levels(data[, which(cols == numCol)]))
    }
  }
  TD <- structure(data,
                  class = c("TD", "data.frame"))
  if (!is.null(design)) {
    attr(x = TD, which = "design") <- design
  }
  return(TD)
}

#' @rdname TD
#' @export
is.TD <- function(x) {
  inherits(x, "TD")
}

#' Summary statistics of the trait
#'
#' This function is to calculate summary statistics of all traits.
#'
#' @param object An object of class TD.
#' @param ... Further arguments - currently not used.
#' @param traits A character vector specifying the name(s) of the traits
#' to be summarised.
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
#' data(TDHeat05)
#' summary(object = TDHeat05, traits = "yield")
#'
#' @export
summary.TD <- function(object,
                       ...,
                       traits,
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
  stats <- matrix(nrow = 22, ncol = length(traits))
  colnames(stats) <- traits
  rownames(stats) <- c("Number of values","Number of observations","Number of missing values",
                       "Mean","Median", "Min","Max","Range","Lower quartile","Upper quartile",
                       "Standard deviation", "Standard error of mean","Variance",
                       "Standard error of variance", "Coefficient of variation",
                       "sum of values", "sum of squares", "Uncorrected sum of squares",
                       "Skewness", "Standard Error of Skewness", "Kurtosis",
                       "Standard Error of Kurtosis")
  for (i in 1:length(traits)) {
    if (nVals) stats[1, i] <- length(object[, traits[i]])
    if (nObs) stats[2, i] <- length(na.omit(object[, traits[i]]))
    if (nMiss) stats[3, i] <- sum(is.na(object[, traits[i]]))
    if (mean) stats[4, i] <- mean(object[, traits[i]], na.rm = TRUE)
    if (median) stats[5, i] <- median(object[, traits[i]], na.rm = TRUE)
    if (min) stats[6, i] <- min(object[, traits[i]], na.rm = TRUE)
    if (max) stats[7, i] <- max(object[, traits[i]], na.rm = TRUE)
    if (range) stats[8, i] <- max(object[, traits[i]], na.rm = TRUE) -
        min(object[, traits[i]], na.rm = TRUE)
    if (lowerQ) stats[9, i] <- quantile(object[, traits[i]], prob = .25, na.rm = TRUE)
    if (upperQ) stats[10, i] <- quantile(object[,traits[i]], prob = .75, na.rm = TRUE)
    if (sd) stats[11, i] <- sd(object[, traits[i]], na.rm = TRUE)
    if (seMean) stats[12, i] <- sd(object[, traits[i]], na.rm = TRUE) /
        sqrt(length(na.omit(object[, traits[i]])))
    if (var) stats[13, i] <- var(object[, traits[i]], na.rm = TRUE)
    if (seVar) stats[14, i] <- seVar(object[, traits[i]], na.rm = TRUE)
    if (CV) stats[15, i] <- 100 * sd(object[, traits[i]], na.rm = TRUE) /
        mean(object[, traits[i]], na.rm = TRUE)
    if (sum) stats[16, i] <- sum(object[, traits[i]], na.rm = TRUE)
    if (sumSq) stats[17, i] <- sum((na.omit(object[, traits[i]]) -
                                      mean(object[, traits[i]], na.rm = TRUE)) ^ 2)
    if (uncorSumSq) stats[18, i] <- sum(object[, traits[i]] ^ 2, na.rm = TRUE)
    if (skew) stats[19, i] <- skewness(object[, traits[i]], na.rm = TRUE)
    if (seSkew) stats[20, i] <- seSkewness(length(na.omit(object[, traits[i]])))
    if (kurt) stats[21, i] <- kurtosis(object[, traits[i]], na.rm = TRUE)
    if (seKurt) stats[22, i] <- seKurtosis(length(na.omit(object[, traits[i]])))
  }
  return(structure(stats,
                   class = c("summary.TD", "table")))
}

#' @export
print.summary.TD <- function(x, ...) {
  for (i in 1:ncol(x)) {
    cat("\nSummary statistics for", colnames(x)[i], "\n\n")
    if (!is.na(x[1, i]))
      cat(paste0("             Number of values = ", x[1, i], "\n"))
    if (!is.na(x[2, i]))
      cat(paste0("       Number of observations = ", x[2, i], "\n"))
    if (!is.na(x[3, i]))
      cat(paste0("     Number of missing values = ", x[3, i], "\n"))
    if (!is.na(x[4, i]))
      cat(paste0("                         Mean = ", round(x[4, i], 2), "\n"))
    if (!is.na(x[5, i]))
      cat(paste0("                       Median = ", round(x[5, i], 2), "\n"))
    if (!is.na(x[6, i]))
      cat(paste0("                          Min = ", round(x[6, i], 2), "\n"))
    if (!is.na(x[7, i]))
      cat(paste0("                          Max = ", round(x[7, i], 2), "\n"))
    if (!is.na(x[8, i]))
      cat(paste0("               Range(max-min) = ", round(x[8, i], 2), "\n"))
    if (!is.na(x[9, i]))
      cat(paste0("               Lower quartile = ", round(x[9, i], 2), "\n"))
    if (!is.na(x[10, i]))
      cat(paste0("               Upper quartile = ", round(x[10, i], 2), "\n"))
    if (!is.na(x[11, i]))
      cat(paste0("           Standard deviation = ", round(x[11, i], 3), "\n"))
    if (!is.na(x[12, i]))
      cat(paste0("       Standard error of mean = ", round(x[12, i], 3), "\n"))
    if (!is.na(x[13, i]))
      cat(paste0("                     Variance = ", round(x[13, i], 3), "\n"))
    if (!is.na(x[14, i]))
      cat(paste0("   Standard error of variance = ", round(x[14, i], 3), "\n"))
    if (!is.na(x[15, i]))
      cat(paste0("     Coefficient of variation = ", round(x[15, i], 3), "\n"))
    if (!is.na(x[16, i]))
      cat(paste0("                sum of values = ", round(x[16, i], 2), "\n"))
    if (!is.na(x[17, i]))
      cat(paste0("               sum of Squares = ", round(x[17, i], 2), "\n"))
    if (!is.na(x[18, i]))
      cat(paste0("   Uncorrected sum of squares = ", round(x[18, i], 2), "\n"))
    if (!is.na(x[19, i]))
      cat(paste0("                     Skewness = ", round(x[19, i], 3), "\n"))
    if (!is.na(x[20, i]))
      cat(paste0("   Standard Error of Skewness = ", round(x[20, i], 3), "\n"))
    if (!is.na(x[21, i]))
      cat(paste0("                     Kurtosis = ", round(x[21, i], 3), "\n"))
    if (!is.na(x[22, i]))
      cat(paste0("   Standard Error of Kurtosis = ", round(x[22, i], 3), "\n"))
    cat("\n")
  }
}






