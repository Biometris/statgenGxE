#' S3 class TD
#'
#' Function for creating objects of S3 class TD (Trial Data). The input data is
#' checked and columns are renamed to default column names for ease of further
#' computations. The columns for genotype, env, megaEnv, year, repId, subBlock,
#' rowId, colId and checkId are converted to factor columns, whereas
#' rowCoordinates and colCoordinates are converted to numerical columns. One
#' single column can be mapped to multiple defaults, e.g. one column with x
#' coordinates can be mapped to both colId and colCoordinates.\cr
#' Columns other than the default columns, e.g. traits or other covariates
#' will be included in the output unchanged.\cr
#' \code{\link{print}} and \code{\link{summary}} and \code{\link{plot}} methods
#' are available.
#'
#' @param data A data.frame containing trial data with a least a column for
#' genotype.
#' @param genotype An optional character string indicating the column in
#' \code{data} that contains genotypes.
#' @param env An optional character string indicating the column in \code{data}
#' that contains environments.
#' @param megaEnv An optional character string indicating the column in
#' \code{data} that contains mega-environments as constructed by
#' \code{\link{gxeMegaEnv}}.
#' @param year An optional character string indicating the column in \code{data}
#' that contains years.
#' @param repId An optional character string indicating the column in \code{data}
#' that contains replicates.
#' @param subBlock An optional character string indicating the column in
#' \code{data} that contains sub blocks.
#' @param rowId An optional character string indicating the column in
#' \code{data} that contains field rows.
#' @param colId An optional character string indicating the column in
#' \code{data} that contains field columns.
#' @param rowCoordinates An optional character string indicating the column in
#' \code{data} that contains the rowId coordinates used for fitting spatial
#' models.
#' @param colCoordinates An optional character string indicating the column in
#' \code{data} that contains the column coordinates used for fitting spatial
#' models.
#' @param checkId An optional character string indicating the column in
#' \code{data} that contains the check IDs.
#' @param design An optional character string indicating the design of the
#' trial. Either are "ibd" (incomplete-block design), "res.ibd"
#' (resolvable incomplete-block design), "rcbd" (randomized complete block
#' design), "rowcol" (row-column design) or "res.rowcol" (resolvable
#' row-column design).
#'
#' @return An object of class TD, the input data.frame with renamed columns
#' and an attribute \code{renamedCols} containing info on which columns have
#' been renamed.
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{summary.TD}}, \code{\link{plot.TD}}
#'
#' @name TD
NULL

#' @rdname TD
#' @export
createTD <- function(data,
                     genotype = NULL,
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
  for (param in c(genotype, env, megaEnv, year, repId, subBlock, rowId, colId,
                  rowCoordinates, colCoordinates, checkId)) {
    if (!is.null(param) && (!is.character(param) || length(param) > 1 ||
                            !param %in% cols)) {
      stop(paste(deparse(param), "has to be NULL or a column in data.\n"))
    }
  }
  if (!is.null(design)) {
    design <- match.arg(design, choices = c("ibd", "res.ibd", "rcbd",
                                            "rowcol", "res.rowcol"))
  }
  ## Create list of reserved column names for renaming columns.
  renameCols <- c("genotype", "env", "megaEnv", "year", "repId", "subBlock",
                  "rowId", "colId", "rowCoordinates", "colCoordinates",
                  "checkId")
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
  factorCols <-  c("genotype", "env", "megaEnv", "year", "repId", "subBlock",
                   "rowId", "colId", "checkId")
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
        as.numeric(data[, which(cols == numCol)])
    }
  }
  TD <- structure(data,
                  class = c("TD", "data.frame"))
  if (!is.null(design)) {
    attr(x = TD, which = "design") <- design
  }
  if (nrow(renamed) > 0) {
    attr(TD, "renamedCols") <- renamed
  }
  return(TD)
}

#' Summarizing objects of class \code{TD}
#'
#' \code{summary} method for class \code{TD}.
#'
#' @param object An object of class TD.
#' @param ... Further arguments - currently not used.
#' @param traits A character vector specifying the name(s) of the traits
#' to be summarised.
#' @param what A character vector indicating which summary statistics should be
#' computed. If \code{what = "all"} all available statistics are computed.\cr
#' Possible options are\cr
#' \describe{
#' \item{nVals}{The number of values, i.e. non-missing + missing values.}
#' \item{nObs}{The number of non-missing observations.}
#' \item{nMiss}{The number of missing values.}
#' \item{mean}{The mean.}
#' \item{median}{The median.}
#' \item{min}{The minimum.}
#' \item{max}{The maximum.}
#' \item{range}{The range (maximum - minimum).}
#' \item{lowerQ}{The lower (25\%) quantile.}
#' \item{upperQ}{The upper (75\%) quantile.}
#' \item{sd}{The standard deviation.}
#' \item{seMean}{The standard error of mean.}
#' \item{var}{The variance.}
#' \item{seVar}{The standard error of variance.}
#' \item{CV}{The coefficient of variation.}
#' \item{sum}{The sum.}
#' \item{sumSq}{The sum of squares.}
#' \item{uncorSumSq}{The uncorrected sum of squares.}
#' \item{skew}{The skewness.}
#' \item{seSkew}{The standard error of the skewness.}
#' \item{kurt}{The kurtosis.}
#' \item{seKurt}{The standard error of the kurtosis.}
#' \item{all}{All summary statistics.}
#' }
#'
#' @return A table containing the selected summary statistics.
#' @seealso \code{\link{createTD}}
#'
#' @examples
#' ## Summarize TDHeat05.
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
  ## Checks.
  if (!is.character(traits) || !all(traits %in% colnames(object))) {
    stop(paste("All traits should be columns in", substitute(object), ".\n"))
  }
  if (what[[1]] == "all") {
    what <- allWhat
  }
  if (!is.character(what) || sum(what %in% allWhat) == 0) {
    stop("At least one statistic should be chosen.\n")
  }
  whichWhat <- which(allWhat %in% what)
  ## Create a data.frame to store the values".
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
  rownames(stats) <- c("Number of values", "Number of observations",
                       "Number of missing values", "Mean", "Median", "Min",
                       "Max", "Range", "Lower quartile", "Upper quartile",
                       "Standard deviation", "Standard error of mean", "Variance",
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
  decimals <- c(rep(x = 0, times = 3), rep(x = 2, times = 7),
                rep(x = 3, times = 5), rep(x = 2, times = 3),
                rep(x = 3, times = 4))[whichWhat]
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

#' Plot function for class TD
#'
#' Four types of plot can be made. Boxplots and histograms can be made for
#' all objects of class \code{\link{TD}}. In case there is a column \code{"env"}
#' in TD boxplots and histograms will be made per environment.\cr
#' Scatterplots and plots of correlation between environments can only be made
#' if a column \code{"env"} is present in TD and will result in an error if this
#' is not the case.
#'
#' @param x An object of class TD.
#' @param ... Further graphical parameters. For boxplots, histograms and
#' scatterplots all graphical parameters can be changed, for the correlation
#' only \code{"main"} will be used and other parameters are ignored.
#' @param trait A character string specifying the name of the trait to be
#' plotted.
#' @param plotType A character string indicating which plot should be made.
#' Either \code{"box"} for a boxplot, \code{"hist"} for histograms,
#' \code{"scatter"} for scatter plots and correlations or \code{"cor"} for a
#' plot of the correlations between environments.
#'
#' @return One of four plots depending on \code{plotType}.
#'
#' @seealso \code{\link{TD}}
#'
#' @examples
#' ## Create a boxplot for TDMaize.
#' plot(TDMaize, trait = "yld")
#' ## Create a histogram for TDMaize.
#' plot(TDMaize, trait = "yld", plotType = "hist")
#' ## Create a scatter plot for TDMaize.
#' plot(TDMaize, trait = "yld", plotType = "scatter")
#' ## Create a plot of correlations between environments for TDMaize.
#' plot(TDMaize, trait = "yld", plotType = "cor")
#'
#'
#' @importFrom utils modifyList
#' @export
plot.TD <- function(x,
                    ...,
                    trait,
                    plotType = c("box", "hist", "scatter", "cor")) {
  if (!is.character(trait) || !trait %in% colnames(x)) {
    stop(paste("Traits should be a column in", substitute(x), ".\n"))
  }
  plotType <- match.arg(plotType)
  dotArgs <- list(...)
  isEnv <- "env" %in% colnames(x)
  if (plotType == "box") {
    ## Set arguments for boxplot
    if (isEnv) {
      ## Call for simple boxplot is different than boxplot per env.
      ## Therefore diffentiate between those two.
      bpArgs <- list(formula = as.formula(paste(trait, "~ env")),
                     data = x)
    } else {
      bpArgs <- list(x = x[[trait]])
    }
    ## Arguments that are equal for both calls.
    ## show.names needed to display envname even when there is only 1 env.
    bpArgs <- c(bpArgs, list(
      main = paste0("Boxplot", ifelse(isEnv, " by Enviroment", "")),
      xlab = ifelse(isEnv, "Enviroment", ""),
      ylab = trait, show.names = TRUE))
    ## Add and overwrite args with custom args from ...
    fixedArgs <- c("formula", "data", "x")
    bpArgs <- modifyList(bpArgs, dotArgs[!names(dotArgs) %in% fixedArgs])
    do.call(boxplot, args = bpArgs)
  } else if (plotType == "hist") {
    ## Set arguments for histogram
    histArgs <- list(x = as.formula(paste("~", trait, if (isEnv) "| env")),
                     data = x, main = paste0("Histogram",
                                             ifelse(isEnv, " by Enviroment", "")),
                     xlab = trait)
    ## Add and overwrite args with custom args from ...
    fixedArgs <- c("x", "data")
    histArgs <- modifyList(histArgs, dotArgs[!names(dotArgs) %in% fixedArgs])
    do.call(lattice::histogram, args = histArgs)
  } else if (plotType == "scatter") {
    if (isEnv && dplyr::n_distinct(x$env) > 1) {
      ## Function for 'plotting' absolute correlations with text size
      ## proportional to the correlations.
      panelCor <- function(x, y, ...) {
        oldPar <- par(usr = c(0, 1, 0, 1))
        on.exit(par(oldPar))
        ## Compute correlations.
        r <- abs(cor(x, y))
        ## Format and print correlations with proper size.
        txt <- format(c(r, 0.123456789), digits = 2)[1]
        cexCor <- 0.8 / strwidth(txt)
        text(0.5, 0.5, txt, cex = cexCor * r)
      }
      ## Function for plotting histograms on the diagonal.
      panelHist <- function(x, ...) {
        oldPar <- par(usr = c(par("usr")[1:2], 0, 1.5))
        on.exit(par(oldPar))
        ## Create but don't plot histogram.
        h <- hist(x, plot = FALSE)
        ## Rescale.
        breaks <- h$breaks
        nB <- length(breaks)
        y <- h$counts
        y <- y / max(y)
        ## Plot rescaled histogram.
        rect(breaks[-nB], 0, breaks[-1], y, ...)
      }
      ## Compute mean trait value per genotype per environment.
      X <- tapply(X = x[, trait], INDEX = x[, c("genotype", "env")],
                  FUN = mean, na.rm = TRUE)
      ## Set arguments for pairs.
      #diag.panel = panelHist, -- suppressed for now.
      pairsArgs <- list(x = X, upper.panel = panelCor,
                        main = paste("Scatterplot matrix and correlations",
                                     "by enviroment:", trait))
      ## Add and overwrite args with custom args from ...
      fixedArgs <- c("x", "upper.panel")
      pairsArgs <- modifyList(pairsArgs, dotArgs[!names(dotArgs) %in% fixedArgs])
      ## Create scatterplots with absolute correlations on the upper part.
      do.call(pairs, args = pairsArgs)
    } else {
      stop(paste("No column env in data or column env contains only 1 environment.\n",
                 "Scatterplot cannot be made.\n"))
    }
  } else if (plotType == "cor") {
    if (isEnv) {
      ## Compute mean trait value per genotype per environment.
      myDat <- tapply(X = x[, trait], INDEX = x[, c("genotype", "env")],
                      FUN = mean, na.rm = TRUE)
      ## Compute correlation matrix.
      corMat <- cor(myDat, use = "pairwise.complete.obs")
      corMat[upper.tri(corMat, diag = FALSE)] <- NA
      ## Plot correlation matrix.
      ## Set arguments for plotCorMat
      corMatArgs <- list(corMat = corMat,
                         main = paste("Correlations between environments for", trait))
      ## Add and overwrite args with custom args from ...
      fixedArgs <- c("corMat")
      corMatArgs <- modifyList(corMatArgs, dotArgs[!names(dotArgs) %in% fixedArgs])
      do.call(plotCorMat, corMatArgs)
    } else {
      stop("No column env in data. Correlation plot cannot be made.\n")
    }
  }
}
