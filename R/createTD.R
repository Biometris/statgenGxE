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
  ## Create a data.frame with renamed cols to add to TD as an attribute.
  renamed <- data.frame(orig = renameFrom[renameFrom != "NULL"],
                        new = renameCols[renameFrom != "NULL"],
                        stringsAsFactors = FALSE)
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
  attr(TD, "renamedCols") <- renamed
  return(TD)
}

#' @rdname TD
#' @export
is.TD <- function(x) {
  inherits(x, "TD")
}

#' Summarizing objects of class \code{TD}
#'
#' \code{summary} method for class \code{TD}.
#'
#' @param object An object of class TD.
#' @param ... Further arguments - currently not used.
#' @param traits A character vector specifying the name(s) of the traits
#' to be summarised.
#' @param what A character vector indicating which statistics should be computed.\cr
#' If \code{what = "all"} all available statistics are computed.\cr
#' Possible options are\cr
#' \describe{
#' \item{nVals}{the number of values}
#' \item{nObs}{the number of observations}
#' \item{nMiss}{the number of missing values}
#' \item{mean}{the mean}
#' \item{median}{the median}
#' \item{min}{the minimum}
#' \item{max}{the maximum}
#' \item{range}{the range (maximum - minimum)}
#' \item{lowerQ}{the lower (25\%) quantile}
#' \item{upperQ}{the upper (75\%) quantile}
#' \item{sd}{the standard deviation}
#' \item{seMean}{standard error of mean}
#' \item{var}{the variance}
#' \item{seVar}{the standard error of variance}
#' \item{CV}{the coefficient of variation}
#' \item{sum}{the sum}
#' \item{sumSq}{sum of squares}
#' \item{uncorSumSq}{uncorrected sum of squares}
#' \item{skew}{the skewness}
#' \item{seSkew}{the standard error of skewness}
#' \item{kurt}{the kurtosis}
#' \item{seKurt}{the standard error of kurtosis}
#' }
#'
#' @return A data.frame containing the selected summary statistics.
#' @seealso \code{\link{createTD}}
#'
#' @examples
#' summary(object = TDHeat05, traits = "yield")
#'
#' @export
summary.TD <- function(object,
                       ...,
                       traits,
                       what = c("nObs", "nMiss", "mean", "median", "min",
                                "max", "lowerQ", "upperQ", "var")) {
  allWhat <- c("nVals", "nObs", "nMiss", "mean", "median", "min",
               "max", "range", "lowerQ", "upperQ", "sd", "seMean",
               "var", "seVar", "CV", "sum", "sumSq", "uncorSumSq",
               "skew", "seSkew", "kurt", "seKurt")
  if (what[[1]] == "all") {
    what <- allWhat
  }
  whichWhat <- which(allWhat %in% what)
  ## Create a data.frame to store the values
  stats <- matrix(nrow = length(what), ncol = length(traits),
                  dimnames = list(what, traits))
  for (i in 1:length(traits)) {
    if ("nVals" %in% what) {
      stats["nVals", i] <- length(object[, traits[i]])
    }
    if ("nObs" %in% what) {
      stats["nObs", i] <- length(na.omit(object[, traits[i]]))
    }
    if ("nMiss" %in% what) {
      stats["nMiss", i] <- sum(is.na(object[, traits[i]]))
    }
    if ("mean" %in% what) {
      stats["mean", i] <- mean(object[, traits[i]], na.rm = TRUE)
    }
    if ("median" %in% what) {
      stats["median", i] <- median(object[, traits[i]], na.rm = TRUE)
    }
    if ("min" %in% what) {
      stats["min", i] <- min(object[, traits[i]], na.rm = TRUE)
    }
    if ("max" %in% what) {
      stats["max", i] <- max(object[, traits[i]], na.rm = TRUE)
    }
    if ("range" %in% what) {
      stats["range", i] <- max(object[, traits[i]], na.rm = TRUE) -
        min(object[, traits[i]], na.rm = TRUE)
    }
    if ("lowerQ" %in% what) {
      stats["lowerQ", i] <- quantile(object[, traits[i]], prob = .25, na.rm = TRUE)
    }
    if ("upperQ" %in% what) {
      stats["upperQ", i] <- quantile(object[,traits[i]], prob = .75, na.rm = TRUE)
    }
    if ("sd" %in% what) {
      stats["sd", i] <- sd(object[, traits[i]], na.rm = TRUE)
    }
    if ("seMean" %in% what) {
      stats["seMean", i] <- sd(object[, traits[i]], na.rm = TRUE) /
        sqrt(length(na.omit(object[, traits[i]])))
    }
    if ("var" %in% what) {
      stats["var", i] <- var(object[, traits[i]], na.rm = TRUE)
    }
    if ("seVar" %in% what) {
      stats["seVar", i] <- seVar(object[, traits[i]], na.rm = TRUE)
    }
    if ("CV" %in% what) {
      stats["CV", i] <- 100 * sd(object[, traits[i]], na.rm = TRUE) /
        mean(object[, traits[i]], na.rm = TRUE)
    }
    if ("sum" %in% what) {
      stats["sum", i] <- sum(object[, traits[i]], na.rm = TRUE)
    }
    if ("sumSq" %in% what) {
      stats["sumSq", i] <- sum((na.omit(object[, traits[i]]) -
                                  mean(object[, traits[i]], na.rm = TRUE)) ^ 2)
    }
    if ("uncorSumSq" %in% what) {
      stats["uncorSumSq", i] <- sum(object[, traits[i]] ^ 2, na.rm = TRUE)
    }
    if ("skew" %in% what) {
      stats["skew", i] <- skewness(object[, traits[i]], na.rm = TRUE)
    }
    if ("seSkew" %in% what) {
      stats["seSkew", i] <- seSkewness(length(na.omit(object[, traits[i]])))
    }
    if ("kurt" %in% what) {
      stats["kurt", i] <- kurtosis(object[, traits[i]], na.rm = TRUE)
    }
    if ("seKurt" %in% what) {
      stats["seKurt", i] <- seKurtosis(length(na.omit(object[, traits[i]])))
    }
  }
  rownames(stats) <- c("Number of values","Number of observations","Number of missing values",
                       "Mean","Median", "Min","Max","Range","Lower quartile","Upper quartile",
                       "Standard deviation", "Standard error of mean","Variance",
                       "Standard error of variance", "Coefficient of variation",
                       "sum of values", "sum of squares", "Uncorrected sum of squares",
                       "Skewness", "Standard Error of Skewness", "Kurtosis",
                       "Standard Error of Kurtosis")[whichWhat]
  attr(x = stats, which = "whichWhat") <- whichWhat
  return(structure(stats,
                   class = c("summary.TD", "table")))
}

#' @export
print.summary.TD <- function(x, ...) {
  whichWhat <- attr(x, "whichWhat")
  decimals <- c(rep(x = 0, times = 3), rep(x = 2, times = 7), rep(x = 3, times = 5),
                rep(x = 2, times = 3), rep(x = 3, times = 4))[whichWhat]
  maxLength <- max(nchar(rownames(x)))
  for (i in 1:ncol(x)) {
    cat("\nSummary statistics for", colnames(x)[i], "\n\n")
    for (j in 1:nrow(x)) {
      cat(paste0(paste0(rep(x = " ",
                            times = maxLength - nchar(rownames(x)[j]) + 2),
                        collapse = ""),
          rownames(x)[j], "  ", round(x[j, i], decimals[j]), "\n"))
    }
  }
  cat("\n")
}






