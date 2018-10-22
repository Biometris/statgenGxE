#' S3 class TD
#'
#' \code{createTD}\cr
#' Function for creating objects of S3 class TD (Trial Data). The function
#' performs the following steps:
#' \itemize{
#' \item{Check input data}
#' \item{Rename columns to default column names - default column names:
#' genotype, trial, loc, year, repId, subBlock, rowCoord, colCoord, rowId,
#' colId, checkId}
#' \item{Convert column types to default column types - rowCoord and colCoord
#' are converted to numeric columns, all other renamed columns to factor
#' columns. Columns other than the default columns, e.g. traits or other
#' covariates will be included in the output unchanged.}
#' \item{Split input data by trial - each trial in the input data will become
#' a list item in the output}
#' \item{Add meta data - the trial meta data are added as attributes to the
#' different output items. The function parameters starting with tr provide
#' the meta data. Their values will be recycled if needed, so by setting a
#' single trDesign all trials will get the same design. The meta data can be
#' changed later on using \code{getMeta} and \code{setMeta}}
#' }
#' \code{addTD}\cr
#' Function for adding extra trial data to an existing object of class TD. The
#' data for the new trials will be added after the data for existing trials. It
#' is possible to add data for an already existing trial, but this will cause
#' multiple items in the output with identical names, which might cause problems
#' later on in the analysis. Therefore a warning will be issued in this
#' case.\cr\cr
#' \code{dropTD}\cr
#' Function for removing data for selected trials from an existing object of
#' class TD.\cr\cr
#' \code{\link{summary.TD}} and \code{\link{plot.TD}} methods are available.
#'
#' @param data A data.frame containing trial data with a least a column for
#' genotype.
#' @param genotype An optional character string indicating the column in
#' \code{data} that contains genotypes.
#' @param trial An optional character string indicating the column in
#' \code{data} that contains trials.
#' @param loc An optional character string indicating the column in
#' \code{data} that contains trial locations.
#' @param year An optional character string indicating the column in \code{data}
#' that contains years.
#' @param repId An optional character string indicating the column in
#' \code{data} that contains replicates.
#' @param subBlock An optional character string indicating the column in
#' \code{data} that contains sub blocks.
#' @param plot An optional character string indicating the column in
#' \code{data} that contains plots. This column will be combined with trial
#' to a single output factor.
#' @param rowCoord An optional character string indicating the column in
#' \code{data} that contains the row coordinates.
#' @param colCoord An optional character string indicating the column in
#' \code{data} that contains the column coordinates.
#' @param rowId An optional character string indicating the column in
#' \code{data} that contains field rows. If not supplied this is assumed to
#' be the same as rowCoord.
#' @param colId An optional character string indicating the column in
#' \code{data} that contains field columns. If not supplied this is assumed to
#' be the same as colCoord.
#' @param checkId An optional character string indicating the column in
#' \code{data} that contains the check IDs.
#' @param trLocation An optional character vector indicating the locations of
#' the trials. This will be used for constructing default names for plots and
#' such. If no locations are provided the trialname will be used as default
#' name.
#' @param trDate An optional date vector indicating the dates of the trials.
#' @param trDesign An optional character vector indicating the designs of the
#' trials. Either "none" (no (known) design), "ibd" (incomplete-block design),
#' "res.ibd" (resolvable incomplete-block design), "rcbd" (randomized complete
#' block design), "rowcol" (row-column design) or "res.rowcol" (resolvable
#' row-column design).
#' @param trLat An optional numerical vector indicating the latitudes of the
#' trials on a scale of -90 to 90.
#' @param trLong An optional numerical vector indicating the longitudes of the
#' trials on a scale of -180 to 180.
#' @param trPlWidth An optional positive numerical vector indicating the
#' widths of the plots.
#' @param trPlLength An optional positive numerical vector indicating the
#' lengths of the plots.
#'
#' @return An object of class TD, a list of data.frames with renamed columns
#' and an attribute \code{renamedCols} containing info on which columns have
#' been renamed. For each unique value of trial the output has a data.frame in
#' the list with the same name as the trial. These data.frames have attributes
#' containing the metadata for the corresponding trial. If there is no column
#' for trial the list will contain one item named after the input data.
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{summary.TD}}, \code{\link{plot.TD}},
#' \code{\link{getMeta}}, \code{\link{setMeta}}
#'
#' @importFrom utils hasName
#'
#' @name TD
NULL

#' @rdname TD
#' @export
createTD <- function(data,
                     genotype = NULL,
                     trial = NULL,
                     loc = NULL,
                     year = NULL,
                     repId = NULL,
                     subBlock = NULL,
                     plot = NULL,
                     rowCoord = NULL,
                     colCoord = NULL,
                     rowId = rowCoord,
                     colId = colCoord,
                     checkId = NULL,
                     trLocation = NULL,
                     trDate = NULL,
                     trDesign = NULL,
                     trLat = NULL,
                     trLong = NULL,
                     trPlWidth = NULL,
                     trPlLength = NULL) {
  ## Save name of original data for naming output.
  dataName <- deparse(substitute(data))
  if (length(dataName) > 1) {
    dataName <- "dat"
  }
  ## Checks.
  if (missing(data) || !is.data.frame(data)) {
    stop("data has to be a data.frame.\n")
  }
  cols <- colnames(data)
  for (param in c(genotype, trial, loc, year, repId, subBlock, plot,
                  rowId, colId, rowCoord, colCoord, checkId)) {
    if (!is.null(param) && (!is.character(param) || length(param) > 1 ||
                            !hasName(data, param))) {
      stop(paste(deparse(param), "has to be NULL or a column in data.\n"))
    }
  }
  checkTDMeta(trDesign = trDesign, trLat = trLat, trLong = trLong,
              trPlWidth = trPlWidth, trPlLength = trPlLength)
  ## Create list of reserved column names for renaming columns.
  renameCols <- c("genotype", "trial", "loc", "year", "repId", "plot",
                  "subBlock", "rowId", "colId", "rowCoord", "colCoord",
                  "checkId")
  ## First rename duplicate colums and add duplicated columns to data
  renameFrom <- as.character(sapply(X = renameCols, FUN = function(x) {
    get(x)
  }))
  ## Create a data.frame with renamed cols to add to TD as an attribute.
  renamed <- data.frame(orig = renameFrom[renameFrom != "NULL"],
                        new = renameCols[renameFrom != "NULL"],
                        stringsAsFactors = FALSE)
  ## Get duplicate columns.
  dupCols <- which(duplicated(renameFrom) & renameFrom != "NULL")
  for (dupCol in dupCols) {
    ## Copy original column as extra column in data for each duplicate.
    tempName <- paste0(".temp", dupCol)
    data[tempName] <- data[, colnames(data) == renameFrom[dupCol]]
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
  factorCols <-  c("genotype", "trial", "loc", "year", "repId", "subBlock",
                   "plot", "rowId", "colId", "checkId")
  for (factorCol in factorCols) {
    if (hasName(data, factorCol)) {
      data[cols == factorCol] <- as.factor(data[, cols == factorCol])
    }
  }
  ## Combine plot and trial into a single factor if both are available.
  ## If trial is not available plot itself was converted to factor in the
  ## previous step.
  if (all(hasName(data, c("trial", "plot")))) {
    data$plot <- interaction(data$trial, data$plot, sep = "_")
  }
  ## Convert columns to numeric if neccessary.
  numCols <- c("rowCoord", "colCoord")
  for (numCol in numCols) {
    if (hasName(data, numCol) && !is.numeric(data[cols == numCol])) {
      data[cols == numCol] <- as.numeric(data[, cols == numCol])
    }
  }
  if (hasName(data, "trial")) {
    listData <- split(x = data, f = droplevels(data$trial))
  } else {
    listData <- setNames(list(data), dataName)
  }
  ## Define meta data to set from input variables.
  meta <- c("trLocation", "trDate", "trDesign", "trLat", "trLong",
            "trPlWidth", "trPlLength")
  ## Expand input values for meta variables to number of trials.
  metaVals <- sapply(X = meta, FUN = function(m) {
    if (!is.null(get(m))) {
      metaVal <- rep(x = get(m), length.out = length(listData))
      if (is.null(names(metaVal)) || !all(hasName(listData, names(metaVal)))) {
        names(metaVal) <- names(listData)
      }
      return(metaVal)
    } else {
      NULL
    }
  }, simplify = FALSE)
  ## Set meta for all trials in data.
  for (tr in names(listData)) {
    for (m in meta) {
      ## Set meta data. Set to NULL if not in input so meta variable is
      ## left out meta data. This to avoid a list of NULL.
      attr(x = listData[[tr]], which = m) <- unname(metaVals[[m]][tr])
    }
    ## Location should always be filled since it is used in plot titles as
    ## well. Use trial name as default value.
    if (is.null(trLocation)) {
      attr(x = listData[[tr]], which = "trLocation") <- tr
    }
    ## Add a list of columns that have been renamed as attribute to TD.
    attr(x = listData[[tr]], which = "renamedCols") <-
      if (nrow(renamed) > 0) renamed else NULL
  }
  TD <- structure(listData,
                  class = c("TD", "list"))
  return(TD)
}

#' @inheritParams createTD
#'
#' @param TD An object of class TD which should be modified.
#'
#' @rdname TD
#' @export
addTD <- function(TD,
                  data,
                  genotype = NULL,
                  trial = NULL,
                  loc = NULL,
                  year = NULL,
                  repId = NULL,
                  subBlock = NULL,
                  plot = NULL,
                  rowCoord = NULL,
                  colCoord = NULL,
                  rowId = rowCoord,
                  colId = colCoord,
                  checkId = NULL,
                  trLocation = NULL,
                  trDate = NULL,
                  trDesign = NULL,
                  trLat = NULL,
                  trLong = NULL,
                  trPlWidth = NULL,
                  trPlLength = NULL) {
  TDNw <- createTD(data = data, genotype = genotype, trial = trial,
                   loc = loc, year = year, repId = repId,
                   subBlock = subBlock, plot = plot, rowCoord = rowCoord,
                   colCoord = colCoord, rowId = rowId, colId = colId,
                   checkId = checkId, trLocation = trLocation, trDate = trDate,
                   trDesign = trDesign, trLat = trLat, trLong = trLong,
                   trPlWidth = trPlWidth, trPlLength = trPlLength)
  dupTrials <- names(TDNw)[names(TDNw) %in% names(TD)]
  if (length(dupTrials) > 0) {
    warning(paste0("The following trials already existed in TD and will be ",
                   "added again: ", paste(dupTrials, collapse = ", "), ".\n"),
            call. = FALSE)
  }
  TDTot <- c(TD, TDNw)
  class(TDTot) <- c("TD", "list")
  return(TDTot)
}

#' @inheritParams addTD
#'
#' @param trials A character vector of trials that should be removed.
#'
#' @rdname TD
#' @export
dropTD <- function(TD,
                   trials) {
  naTrials <- trials[!trials %in% names(TD)]
  if (length(naTrials) > 0) {
    warning(paste0("The following trials are not in TD: ",
                   paste(naTrials, collapse = ", "), ".\n"), call. = FALSE)
  }
  leftTrials <- names(TD)[!names(TD) %in% trials]
  if (length(leftTrials) == 0) {
    warning("All trials have been removed from TD.\n", call. = FALSE)
  }
  return(TD[!names(TD) %in% trials])
}

#' Summarizing objects of class \code{TD}
#'
#' \code{summary} method for class \code{TD}.
#'
#' @param object An object of class TD.
#' @param ... Further arguments - currently not used.
#' @param trial A character vector specifying the trials to be plotted.
#' @param traits A character vector specifying the traits to be summarised.
#' @param what A character vector indicating which summary statistics should be
#' computed. If \code{what = "all"} all available statistics are computed.\cr
#' Possible options are:
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
                       trial = names(object),
                       traits,
                       what = c("nObs", "nMiss", "mean", "median", "min",
                                "max", "lowerQ", "upperQ", "var")) {
  allWhat <- c("nVals", "nObs", "nMiss", "mean", "median", "min",
               "max", "range", "lowerQ", "upperQ", "sd", "seMean",
               "var", "seVar", "CV", "sum", "sumSq", "uncorSumSq",
               "skew", "seSkew", "kurt", "seKurt")
  ## Checks.
  if (!is.character(trial) || length(trial) > 1 || !trial %in% names(object)) {
    stop(paste0("trial should be a single character string in ",
                deparse(substitute(object)), ".\n"))
  }
  trDat <- object[[trial]]
  if (!is.character(traits) || !all(traits %in% colnames(trDat))) {
    stop("All traits should be columns in trial.\n")
  }
  if (what[[1]] == "all") {
    what <- allWhat
  }
  if (!is.character(what) || !all(what %in% allWhat)) {
    stop("At least one statistic should be chosen.\n")
  }
  whichWhat <- which(allWhat %in% what)
  ## Create a matrix to store the values".
  stats <- matrix(nrow = length(what), ncol = length(traits),
                  dimnames = list(what, traits))
  for (i in 1:length(traits)) {
    if ("nVals" %in% what) {
      stats["nVals", i] <- length(trDat[, traits[i]])
    }
    if ("nObs" %in% what) {
      stats["nObs", i] <- length(na.omit(trDat[, traits[i]]))
    }
    if ("nMiss" %in% what) {
      stats["nMiss", i] <- sum(is.na(trDat[, traits[i]]))
    }
    if ("mean" %in% what) {
      stats["mean", i] <- mean(trDat[, traits[i]], na.rm = TRUE)
    }
    if ("median" %in% what) {
      stats["median", i] <- median(trDat[, traits[i]], na.rm = TRUE)
    }
    if ("min" %in% what) {
      stats["min", i] <- min(trDat[, traits[i]], na.rm = TRUE)
    }
    if ("max" %in% what) {
      stats["max", i] <- max(trDat[, traits[i]], na.rm = TRUE)
    }
    if ("range" %in% what) {
      stats["range", i] <- max(trDat[, traits[i]], na.rm = TRUE) -
        min(trDat[, traits[i]], na.rm = TRUE)
    }
    if ("lowerQ" %in% what) {
      stats["lowerQ", i] <- quantile(trDat[, traits[i]], prob = .25, na.rm = TRUE)
    }
    if ("upperQ" %in% what) {
      stats["upperQ", i] <- quantile(trDat[,traits[i]], prob = .75, na.rm = TRUE)
    }
    if ("sd" %in% what) {
      stats["sd", i] <- sd(trDat[, traits[i]], na.rm = TRUE)
    }
    if ("seMean" %in% what) {
      stats["seMean", i] <- sd(trDat[, traits[i]], na.rm = TRUE) /
        sqrt(length(na.omit(trDat[, traits[i]])))
    }
    if ("var" %in% what) {
      stats["var", i] <- var(trDat[, traits[i]], na.rm = TRUE)
    }
    if ("seVar" %in% what) {
      stats["seVar", i] <- seVar(trDat[, traits[i]], na.rm = TRUE)
    }
    if ("CV" %in% what) {
      stats["CV", i] <- 100 * sd(trDat[, traits[i]], na.rm = TRUE) /
        mean(trDat[, traits[i]], na.rm = TRUE)
    }
    if ("sum" %in% what) {
      stats["sum", i] <- sum(trDat[, traits[i]], na.rm = TRUE)
    }
    if ("sumSq" %in% what) {
      stats["sumSq", i] <- sum((na.omit(trDat[, traits[i]]) -
                                  mean(trDat[, traits[i]], na.rm = TRUE)) ^ 2)
    }
    if ("uncorSumSq" %in% what) {
      stats["uncorSumSq", i] <- sum(trDat[, traits[i]] ^ 2, na.rm = TRUE)
    }
    if ("skew" %in% what) {
      stats["skew", i] <- skewness(trDat[, traits[i]], na.rm = TRUE)
    }
    if ("seSkew" %in% what) {
      stats["seSkew", i] <- seSkewness(length(na.omit(trDat[, traits[i]])))
    }
    if ("kurt" %in% what) {
      stats["kurt", i] <- kurtosis(trDat[, traits[i]], na.rm = TRUE)
    }
    if ("seKurt" %in% what) {
      stats["seKurt", i] <- seKurtosis(length(na.omit(trDat[, traits[i]])))
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
                   class = c("summary.TD", "table"),
                   trial = trial))
}

#' @export
print.summary.TD <- function(x, ...) {
  whichWhat <- attr(x, "whichWhat")
  decimals <- c(rep(x = 0, times = 3), rep(x = 2, times = 7),
                rep(x = 3, times = 5), rep(x = 2, times = 3),
                rep(x = 3, times = 4))[whichWhat]
  maxLength <- max(nchar(rownames(x)))
  for (i in 1:ncol(x)) {
    cat(paste("\nSummary statistics for", colnames(x)[i], "in",
              attr(x, "trial"), "\n\n"))
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
#' Plotting function for objects of class TD. Plots either the layout of the
#' different trials within the TD object or locates the trials on a map. Also a
#' boxplot can be made for selected traits per trial and a plot of correlations
#' between trials. A detailed description and optional extra parameters of the
#' different plots is given in the sections below.
#'
#' @section Layout Plot:
#' A layout plot plots the layout of the selected trials (all available trials
#' by default). This plot can only be made for trials that contain both
#' row (\code{rowCoord}) and column(\code{colCoord}) information. If either one
#' of those is missing the trial is skipped with a warning. If blocks
#' (\code{subBlock}) are available for a trial these are indicated in different
#' colors per block, otherwise all plots are colored in grey. If replicates
#' (\code{repId}) are available a black line is plotted between diffent
#' replicates. Missing plots are indicated in white. This can either be single
#' plots in a trial or complete missing columns or rows.\cr
#' Extra parameter options:
#' \itemize{
#' \item{showGeno} {Should individual genotypes be indicated in the plot?
#' Defaults to \code{FALSE}}
#' \item{highlight} {A character vector of genotypes to be highlighted in the
#' plot.}
#' }
#'
#' @section Map Plot:
#' A map is plotted with the locations of the trials in the TD object.
#' Mapping the trials is done based on lattitude and longitude that can be
#' added when creating an object of class TD. Trials without latitude and/or
#' longitude available are skipped with a warning message. The countries in
#' which the trials are located will be plotted on a single map and the
#' location of the trials will be indicated on this map. The actual plot is
#' made using ggplot, but for getting the data for the borders of the countries
#' the maps package is needed.
#'
#' @section Box Plot:
#' Per selected trait a boxplot is created grouped per trial.
#' Extra parameter options:
#' \itemize{
#' \item{groupBy} {A character string indicating a column in \code{TD} by which
#' the boxes in the plot should be grouped. By default the boxes are grouped
#' per trial.}
#' \item{colorBy} {A character string indicating a column in \code{TD} by which
#' the boxes are colored. Coloring will be done within the groups indicated by
#' the \code{groupBy} parameter.}
#' \item{orderBy} {A character vector indicating columns in \code{TD} on which
#' the trials in the plot should be ordered. By default ordering is done
#' alphabetically.}
#' }
#'
#' @section Correlation Plot:
#' Per selected trait a plot is drawn of correlations between the trials in the
#' \code{TD} object. If genotypes are replicated within trials genotypic means are
#' taken before computing correlations.
#' Extra parameter options:
#' \itemize{
#' \item{orderBy} {A character vector indicating columns in \code{TD} on which
#' the trials in the plot should be ordered. By default ordering is done
#' alphabetically.}
#' }
#'
#' @param x An object of class TD.
#' @param ... Extra plot options. Described per plotType in their respective
#' section.
#' @param plotType A single character string indicating which plot should be
#' made. See the sections below for a detailed explanation of the plots.
#' @param trials A character vector indicating the trials to be plotted when
#' plotting field layouts. Only used if \code{plotType} = "layout".
#' @param traits A character vector indicating the traits to be plotted in
#' a boxplot. Only used if \code{plotType} = "box" or "cor".
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returnfed.
#'
#' @export
plot.TD <- function(x,
                    ...,
                    plotType = c("layout", "map", "box", "cor"),
                    trials = names(x),
                    traits = NULL,
                    output = TRUE) {
  ## Maps seems to change graphics parameters without resetting. Do so here.
  if (!is.character(trials) || !all(trials %in% names(x))) {
    stop(paste0("All trials should be in ", deparse(x), ".\n"))
  }
  plotType <- match.arg(plotType)
  dotArgs <- list(...)
  if (plotType == "layout") {
    showGeno <- isTRUE(dotArgs$showGeno)
    highlight <- dotArgs$highlight
    if (!is.null(highlight) && !is.character(highlight)) {
      stop("highlight should be a character vector.\n")
    }
    p <- setNames(vector(mode = "list", length = length(trials)), trials)
    for (trial in trials) {
      trDat <- x[[trial]]
      if (!"rowCoord" %in% colnames(trDat)) {
        warning(paste0("rowCoord should be a column in ", trial, ".\n",
                       "Plot skipped.\n"), call. = FALSE)
        break
      }
      if (!"colCoord" %in% colnames(trDat)) {
        warning(paste0("colCoord should be a column in ", trial, ".\n",
                       "Plot skipped.\n"), call. = FALSE)
        break
      }
      if (length(highlight) > 0) {
        trDat$highlight. <- ifelse(trDat$genotype %in% highlight,
                                   as.character(trDat$genotype), NA)
      }
      trLoc <- attr(trDat, "trLocation")
      ylen <- attr(trDat, "trPlLength")
      xlen <- attr(trDat, "trPlWidth")
      ## Compute aspect for proper depiction of field size. If no information
      ## is available plots are assumed to be square.
      if (is.null(ylen) || is.null(xlen)) {
        aspect <- length(unique(trDat$colCoord)) /
          length(unique(trDat$rowCoord))
      } else {
        aspect <- ylen / xlen
      }
      ## Create data for lines between replicates.
      if ("repId" %in% colnames(trDat)) {
        yMin <- min(trDat$rowCoord)
        yMax <- max(trDat$rowCoord)
        xMin <- min(trDat$colCoord)
        xMax <- max(trDat$colCoord)
        ## Create matrix containing replicates.
        ## First create an empty matrix containing all row/column values
        ## between min and max to assure complete missing rows/columns
        ## are added.
        M <- matrix(nrow = yMax - yMin + 1, ncol = xMax - xMin + 1,
                    dimnames = list(yMin:yMax, xMin:xMax))
        for (i in 1:nrow(trDat)) {
          M[as.character(trDat[i, "rowCoord"]),
            as.character(trDat[i, "colCoord"])] <- trDat[i, "repId"]
        }
        ## Create an imputed version of M for plotting borders around NA values.
        MImp <- M
        MImp[is.na(MImp)] <- nlevels(trDat$repId) + 1
        has.breaks <- function(x) {
          ncol(x) == 2 & nrow(x) > 0
        }
        ## Create a data.frame with positions where the value of rep in the
        ## data changes in vertical direction.
        vertW <- do.call(rbind.data.frame,
                         Filter(f = has.breaks, x = Map(function(i, x) {
                           cbind(y = i, x = which(diff(c(0, x, 0)) != 0))
                         }, 1:nrow(MImp), split(MImp, 1:nrow(MImp)))))
        ## Remove vertical walls that are on the outside bordering an NA value
        ## to prevent drawing of unneeded lines.
        vertW <- vertW[!(vertW$x == 1 & is.na(M[vertW$y, 1])) &
                         !(vertW$x == ncol(M) + 1 &
                             is.na(M[vertW$y, ncol(M)])), ]
        ## Add min row value for plotting in the correct position.
        vertW$y <- vertW$y + yMin - 1
        vertW$x <- vertW$x + xMin - 1
        ## For horizontal walls follow the same procedure as above.
        horW <- do.call(rbind.data.frame,
                        Filter(f = has.breaks, x = Map(function(i, y) {
                          cbind(x = i, y = which(diff(c(0, y, 0)) != 0))
                        }, 1:ncol(MImp), as.data.frame(MImp))))
        horW <- horW[!(horW$y == 1 & is.na(M[1, horW$x])) &
                       !(horW$y == nrow(M) + 1 & is.na(M[nrow(M), horW$x])), ]
        horW$y <- horW$y + yMin - 1
        horW$x <- horW$x + xMin - 1
      }
      ## Create base plot.
      pTr <- ggplot2::ggplot(data = trDat,
                             ggplot2::aes_string(x = "colCoord",
                                                 y = "rowCoord")) +
        ggplot2::coord_fixed(ratio = aspect,
                             xlim = range(trDat$colCoord) + c(-0.5, 0.5),
                             ylim = range(trDat$rowCoord) + c(-0.5, 0.5),
                             clip = "off") +
        ggplot2::theme(panel.background = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ## Move ticks to edge of the plot.
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(),
                                    expand = c(0, 0)) +
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(),
                                    expand = c(0, 0)) +
        ggplot2::ggtitle(trLoc)
      if (hasName(x = trDat, name = "subBlock") &&
          sum(!is.na(trDat$highlight.)) == 0) {
        ## If subblocks are available color tiles by subblock.
        pTr <- pTr + ggplot2::geom_tile(
          ggplot2::aes_string(fill = "subBlock"), color = "grey50")
      } else if (sum(!is.na(trDat$highlight.)) > 0) {
        ## Genotypes to be highlighted get a color.
        ## Everything else the NA color.
        pTr <- pTr + ggplot2::geom_tile(
          ggplot2::aes_string(fill = "highlight."), color = "grey50") +
          ggplot2::labs(fill = "Highlighted") +
          ## Remove NA from scale.
          ggplot2::scale_fill_discrete(na.translate = FALSE)
      } else {
        ## No subblocks and no hightlights so just a single fill color.
        pTr <- pTr + ggplot2::geom_tile(color = "grey50", fill = "pink")
      }
      if (showGeno) {
        ## Add names of genotypes to the center of the tiles.
        pTr <- pTr + ggplot2::geom_text(ggplot2::aes_string(label = "genotype"),
                                        size = 2, check_overlap = TRUE)
      }
      if ("repId" %in% colnames(trDat)) {
        ## Add lines for replicates.
        pTr <- pTr +
          ## Add verical lines as segment.
          ## adding/subtracting 0.5 assures plotting at the borders of
          ## the tiles.
          ggplot2::geom_segment(
            ggplot2::aes_string(x = "x - 0.5", xend = "x - 0.5",
                                y = "y - 0.5", yend = "y + 0.5",
                                linetype = "'replicates'"), data = vertW,
            size = 1) +
          ggplot2::geom_segment(
            ggplot2::aes_string(x = "x - 0.5", xend = "x + 0.5",
                                y = "y - 0.5", yend = "y - 0.5"), data = horW,
            size = 1) +
          ## Just needed for adding a legend entry for replicates.
          ggplot2::scale_linetype_manual("replicates",
                                         values = c("replicates" = "solid"),
                                         name = ggplot2::element_blank())
      }
      p[[trial]] <- pTr
      if (output) {
        plot(pTr)
      }
    }
  } else if (plotType == "map") {
    ## Create a data.frame for plotting trials.
    ## Population has a random value but if left out nothing is plotted.
    locs <- setNames(getMeta(x)[c("trLocation", "trLat", "trLong")],
                     c("name", "lat", "long"))
    locs <- locs[!is.na(locs$lat) & !is.na(locs$long), ]
    if (nrow(locs) == 0) {
      stop(paste("At least one trial should have latitute and longitude",
                 "for plotting on map.\n"))
    }
    longR <- range(locs$long)
    latR <- range(locs$lat)
    ## Create data useable by ggplot geom_polygon.
    mapDat <- ggplot2::map_data("world", xlim = longR, ylim = latR)
    p <- ggplot2::ggplot(mapDat, ggplot2::aes_string(x = "long", y = "lat")) +
      ggplot2::geom_polygon(ggplot2::aes_string(group = "group"),
                            fill = "white", color = "black") +
      ## Add a proper map projection.
      ggplot2::coord_map(clip = "on", xlim = longR + c(-0.1, 0.1) * diff(longR),
                         ylim = latR + c(-0.1, 0.1) * diff(latR)) +
      ## Add trial locations.
      ggplot2::geom_point(data = locs) +
      ggplot2::geom_text(ggplot2::aes_string(label = "name"), data = locs,
                         color = "red", size = 3, nudge_x = 0.01 * diff(longR),
                         nudge_y = 0.04 * diff(latR), check_overlap = TRUE) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.background =
                       ggplot2::element_rect(fill = "steelblue2")) +
      ggplot2::ggtitle("Trial locations")
    if (output) {
      plot(p)
    }
  } else if (plotType == "box") {
    if (is.null(traits) || !is.character(traits)) {
      stop("traits should be a character vector.\n")
    }
    groupBy <- dotArgs$groupBy
    if (!is.null(groupBy) && (!is.character(groupBy) || length(groupBy) > 1)) {
      stop("groupBy should be a single character string.\n")
    }
    if (!is.null(groupBy) && !all(sapply(X = x, FUN = function(trial) {
      hasName(x = trial, name = groupBy)
    }))) {
      stop("groupBy should be a column in TD.\n")
    }
    colorBy <- dotArgs$colorBy
    if (!is.null(colorBy) && (!is.character(colorBy) || length(colorBy) > 1)) {
      stop("colorBy should be a single character string.\n")
    }
    if (!is.null(colorBy) && !all(sapply(X = x, FUN = function(trial) {
      hasName(x = trial, name = colorBy)
    }))) {
      stop("colorBy should be a column in TD.\n")
    }
    orderBy <- dotArgs$orderBy
    if (!is.null(orderBy) && !is.character(orderBy)) {
      stop("orderBy should be a character vector.\n")
    }
    if (!is.null(orderBy) && !all(sapply(X = x, FUN = function(trial) {
      all(hasName(x = trial, name = orderBy))
    }))) {
      stop("All items in orderBy should be columns in TD.\n")
    }
    p <- setNames(vector(mode = "list", length = length(traits)), traits)
    for (trait in traits) {
      ## Create a single data.frame from x with only columns trial and trait.
      ## trail where trait is not measured/available are removed by setting
      ## them to NULL.
      xVar <- if (is.null(groupBy)) "trial" else groupBy
      plotDat <- Reduce(f = rbind, x = lapply(X = x, function(trial) {
        if (!hasName(x = trial, name = trait)) {
          NULL
        } else {
          trial[c(trait, xVar, if (!is.null(colorBy)) colorBy,
                if (!is.null(orderBy)) orderBy)]
        }
      }))
      if (is.null(plotDat)) {
        warning(paste0(trait, " isn't a column in any of the trials.\n",
                       "Plot skipped.\n"), call. = FALSE)
        break
      }
      if (!is.null(orderBy)) {
        ## Reorder levels in trial so plotting is done according to orderBy.
        ## do.call needed since order doesn't accept a vector as input.
        levNw <- unique(plotDat[[xVar]][do.call(order, lapply(orderBy,
                                                            FUN = function(x) {
                                                              plotDat[x]}))])
        plotDat[xVar] <- factor(plotDat[[xVar]], levels = levNw)
      }
      ## Create boxplot.
      pTr <- ggplot2::ggplot(plotDat, ggplot2::aes_string(x = xVar, y = trait,
                                                          fill = colorBy)) +
        ggplot2::geom_boxplot(na.rm = TRUE) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                           vjust = 0.5,
                                                           hjust = 1))
      p[[trait]] <- pTr
      if (output) {
        plot(pTr)
      }
    }
  } else if (plotType == "cor") {
    if (is.null(traits) || !is.character(traits)) {
      stop("traits should be a character vector.\n")
    }
    orderBy <- dotArgs$orderBy
    if (!is.null(orderBy) && !is.character(orderBy)) {
      stop("orderBy should be a character vector.\n")
    }
    if (!is.null(orderBy) && !all(sapply(X = x, FUN = function(trial) {
      all(hasName(x = trial, name = orderBy))
    }))) {
      stop("All items in orderBy should be columns in TD.\n")
    }
    p <- setNames(vector(mode = "list", length = length(traits)), traits)
    for (trait in traits) {
      ## Create a single data.frame from x with only columns trial and trait.
      ## trails where trait is not measured/available are removed by setting
      ## them to NULL.
      plotDat <- Reduce(f = rbind, x = lapply(X = x, FUN = function(trial) {
        if (!hasName(x = trial, name = trait)) {
          NULL
        } else {
          trial[c("genotype", "trial", trait, if (!is.null(orderBy)) orderBy)]
        }
      }))
      if (is.null(plotDat)) {
        warning(paste0(trait, " isn't a column in any of the trials.\n",
                       "Plot skipped.\n"), call. = FALSE)
        break
      }
      if (!is.null(orderBy)) {
        ## Reorder levels in trial so plotting is done according to orderBy.
        ## do.call needed since order doesn't accept a vector as input.
        levNw <- unique(plotDat$trial[do.call(order, lapply(orderBy,
                                                            FUN = function(x) {
                                                              plotDat[x]}))])
        plotDat$trial <- factor(plotDat$trial, levels = levNw)
      }
      ## Create table with values trait per genotype per trial.
      ## If TD already contains BLUEs/BLUPs taking means doesn't do anything
      ## but it is needed for raw data where there can be replicates.
      plotTab <- tapply(plotDat[[trait]],
                        INDEX = list(plotDat$genotype, plotDat$trial),
                        FUN = mean)
      ## Create a correlation matrix.
      corMat <- cor(plotTab, use = "pairwise.complete.obs")
      ## Melt to get the proper format for ggplot.
      meltedCorMat <- reshape2::melt(corMat)
      ## Remove top left of the plot. Only plotting a bottom right triangle.
      ## Diagonal is removed as well.
      meltedCorMat <- meltedCorMat[as.numeric(meltedCorMat$Var1) >
                                     as.numeric(meltedCorMat$Var2), ]
      ## Create plot.
      pTr <- ggplot2::ggplot(data = meltedCorMat,
                             ggplot2::aes_string("Var1", "Var2",
                                                 fill = "value")) +
        ggplot2::geom_raster() +
        ## Create a gradient scale.
        ggplot2::scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                      na.value = "grey", limit = c(-1, 1)) +
        ## Move y-axis to the right for easier reading.
        ggplot2::scale_y_discrete(position = "right") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                           vjust = 1, size = 6,
                                                           hjust = 1)) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6)) +
        ## Center title.
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        ## Remove grid behind empty bit of triangle.
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank()) +
        ## No axis and legend titles.
        ggplot2::labs(x = "", y = "", fill = "") +
        ggplot2::ggtitle(paste("Correlations of environments for", trait)) +
        ## Fix coordinates to get a square sized plot.
        ggplot2::coord_fixed()
      p[[trait]] <- pTr
      if (output) {
        plot(pTr)
      }
    }
  }
  invisible(p)
}

#' Extract metadata from TD objects
#'
#' Function for extracting metadata as a data.frame from objects of class TD.
#' Location, data, design, latitude, longitude, plotWidth and plotLength for
#' all trials in TD will be extracted and return in the form of a data.frame.
#'
#' @param TD An object of class TD.
#'
#' @return A data.frame containing the metadata for all trials in TD.
#'
#' @export
getMeta <- function(TD) {
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be an object of class TD.\n")
  }
  metaVars <- c("trLocation", "trDate", "trDesign", "trLat", "trLong",
                "trPlWidth", "trPlLength")
  meta <- as.data.frame(matrix(nrow = length(TD), ncol = length(metaVars),
                               dimnames = list(names(TD), metaVars)))
  for (mv in metaVars) {
    meta[mv] <- sapply(X = TD, FUN = function(tr) {
      mvTr <- attr(tr, which = mv)
      ## Replace NULL by NA to ensure correct output format for inserting in df.
      if (!is.null(mvTr)) {
        return(mvTr)
      } else {
        return(NA)
      }
    })
  }
  class(meta$trDate) <- "Date"
  return(meta)
}

#' Set metadata of TD objects
#'
#' Function for setting metadata of a TD object for one or more trials
#' simultaneously. Metadata can be set using a data.frame with rownames
#' corresponding to the trials in \code{TD}. The data.frame should contain one
#' or  more of the following columns: trLocation, trDate, trDesign, trLat,
#' trLong, trPlWidth and trPlLength. The values of the metadata of TD
#' will be set to the values in the corresponding column in \code{meta}.
#' Existing values will be overwritten, but \code{NA} will be ignored so
#' setting a value to \code{NA} won't result in accidentally removing it.
#'
#' @param TD An object of class TD.
#' @param meta A data.frame containing metadata.
#'
#' @return The object of class TD with updated metadata.
#'
#' @export
setMeta <- function(TD,
                    meta) {
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be an object of class TD.\n")
  }
  if (missing(meta) || !inherits(meta, "data.frame")) {
    stop("meta should be a data.frame.\n")
  }
  naTr <- rownames(meta)[!rownames(meta) %in% names(TD)]
  if (length(naTr) > 0) {
    warning(paste0("The following trials in meta are not in TD: ",
                   paste(naTr, collapse = ", "), ".\n"), call. = FALSE)
  }
  metaVars <- c("trLocation", "trDate", "trDesign", "trLat", "trLong",
                "trPlWidth", "trPlLength")
  ## Set metadata for trials in meta that are also in TD.
  for (tr in rownames(meta)[rownames(meta) %in% names(TD)]) {
    for (mv in metaVars) {
      mvTr <- meta[tr, mv]
      if (!is.na(mvTr)) {
        chk <- try(do.call(what = checkTDMeta, args = setNames(list(mvTr), mv)),
                   silent = TRUE)
        if (inherits(chk, "try-error")) {
          ## Get message from check function but remove first 8 chars to
          ## prevent having an error text with 3x error in it.
          stop(paste0("\nError for ", tr, ":\n",
                      substring(text = chk, first = 9)))
        }
        attr(TD[[tr]], which = mv) <- mvTr
      }
    }
  }
  return(TD)
}

#' Function for extracting for objects of class TD that keeps class.
#'
#' @keywords internal
`[.TD` <- function(x, i, ...) {
  r <- NextMethod("[")
  attr(r, "class") <- attr(x, "class")
  return(r)
}

#' Helper function for checking metadata structure for TD objects.
#'
#' @keywords internal
checkTDMeta <- function(trLocation = NULL,
                        trDate = NULL,
                        trDesign = NULL,
                        trLat = NULL,
                        trLong = NULL,
                        trPlWidth = NULL,
                        trPlLength = NULL) {
  if (!is.null(trDesign)) {
    trDesign <- match.arg(trDesign, choices = c("none", "ibd", "res.ibd",
                                                "rcbd", "rowcol", "res.rowcol"),
                          several.ok = TRUE)
  }
  if (!is.null(trLat) && (!is.numeric(trLat) || any(abs(trLat) > 90))) {
    stop("trLat should be a numerical vector between -90 and 90.\n",
         call. = FALSE)
  }
  if (!is.null(trLong) && (!is.numeric(trLong) || any(abs(trLong) > 180))) {
    stop("trLat should be a numerical vector between -180 and 180.\n",
         call. = FALSE)
  }
  if (!is.null(trLat) && !is.null(trLong)) {
    locLen <- max(length(trLat), length(trLong))
    ## Check that coordinates point to a proper location so plotting can be done.
    loc <- maps::map.where(x = rep(x = trLong, length.out = locLen),
                           y = rep(x = trLat, length.out = locLen))
    if (length(loc) > 0 && anyNA(loc)) {
      warning(paste("Values for trLat and trLong should all match a known",
                    "land location.\n"),
              call. = FALSE)
    }
  }
  if (!is.null(trPlWidth) && (!is.numeric(trPlWidth) || any(trPlWidth < 0))) {
    stop("trPlWidth should be a positive numerical vector.\n",
         call. = FALSE)
  }
  if (!is.null(trPlLength) && (!is.numeric(trPlLength) || any(trPlLength < 0))) {
    stop("should be a positive numerical vector.\n",
         call. = FALSE)
  }
}
