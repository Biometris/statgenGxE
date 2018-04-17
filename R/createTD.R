#' S3 class TD
#'
#' \code{createTD}\cr
#' Function for creating objects of S3 class TD (Trial Data). The input data is
#' checked and columns are renamed to default column names for ease of further
#' computations. The columns for genotype, trial, megaEnv, year, repId,
#' subBlock, rowId, colId and checkId are converted to factor columns, whereas
#' rowCoord and colCoord are converted to numerical columns. One
#' single column can be mapped to multiple defaults, e.g. one column with x
#' coordinates can be mapped to both colId and colCoord.\cr
#' Columns other than the default columns, e.g. traits or other covariates
#' will be included in the output unchanged.\cr
#' The input data is split by trial. So for an input data frame with three
#' different trials in the column assigned to \code{trial} the output will be
#' generated as a list with three items. In this case metadata describing the
#' trials will be identical.\cr
#' To generate a TD object with different metadata for each trial start by
#' creating a TD object for one trial and then add new trials to it using the
#' addTD function. Alternatively create a TD object with all data and then use
#' \code{\link{setMeta}} to add metadata for all trials at once.\cr\cr
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
#' @param megaEnv An optional character string indicating the column in
#' \code{data} that contains mega-environments as constructed by
#' \code{\link{gxeMegaEnv}}.
#' @param year An optional character string indicating the column in \code{data}
#' that contains years.
#' @param repId An optional character string indicating the column in
#' \code{data} that contains replicates.
#' @param subBlock An optional character string indicating the column in
#' \code{data} that contains sub blocks.
#' @param rowId An optional character string indicating the column in
#' \code{data} that contains field rows.
#' @param colId An optional character string indicating the column in
#' \code{data} that contains field columns.
#' @param rowCoord An optional character string indicating the column in
#' \code{data} that contains the rowId coordinates used for fitting spatial
#' models.
#' @param colCoord An optional character string indicating the column in
#' \code{data} that contains the column coordinates used for fitting spatial
#' models.
#' @param checkId An optional character string indicating the column in
#' \code{data} that contains the check IDs.
#' @param trLocation An optional character string indicating the location of
#' the trial. This will be used for constructing default names for plots and
#' such. If no location is provided the trialname will be used as default name.
#' @param trDate An optional date indicating the date of the trial.
#' @param trDesign An optional character string indicating the design of the
#' trial. Either "ibd" (incomplete-block design), "res.ibd"
#' (resolvable incomplete-block design), "rcbd" (randomized complete block
#' design), "rowcol" (row-column design) or "res.rowcol" (resolvable
#' row-column design).
#' @param trLat An optional numerical value indicating the latitude of the
#' trial on a scale of -90 to 90.
#' @param trLong An optional numerical value indicating the longitude of the
#' trial on a scale of -180 to 180.
#' @param trPlotWidth An optional positive numerical value indicating the
#' width of the plot.
#' @param trPlotLength An optional positive numerical value indicating the
#' length of the plot.
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
#' @seealso \code{\link{summary.TD}}, \code{\link{plot.TD}}
#'
#' @name TD
NULL

#' @rdname TD
#' @export
createTD <- function(data,
                     genotype = NULL,
                     trial = NULL,
                     megaEnv = NULL,
                     year = NULL,
                     repId = NULL,
                     subBlock = NULL,
                     rowId = NULL,
                     colId = NULL,
                     rowCoord = NULL,
                     colCoord = NULL,
                     checkId = NULL,
                     trLocation = NULL,
                     trDate = NULL,
                     trDesign = NULL,
                     trLat = NULL,
                     trLong = NULL,
                     trPlotWidth = NULL,
                     trPlotLength = NULL) {
  ## Save name of original data for naming output.
  dataName <- deparse(substitute(data))
  ## Checks.
  if (missing(data) || !is.data.frame(data)) {
    stop("data has to be a data.frame.\n")
  }
  cols <- colnames(data)
  for (param in c(genotype, trial, megaEnv, year, repId, subBlock, rowId, colId,
                  rowCoord, colCoord, checkId)) {
    if (!is.null(param) && (!is.character(param) || length(param) > 1 ||
                            !param %in% cols)) {
      stop(paste(deparse(param), "has to be NULL or a column in data.\n"))
    }
  }
  checkTDMeta(trDesign = trDesign, trLat = trLat, trLong = trLong,
              trPlotWidth = trPlotWidth, trPlotLength = trPlotLength)
  ## Create list of reserved column names for renaming columns.
  renameCols <- c("genotype", "trial", "megaEnv", "year", "repId", "subBlock",
                  "rowId", "colId", "rowCoord", "colCoord",
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
  factorCols <-  c("genotype", "trial", "megaEnv", "year", "repId", "subBlock",
                   "rowId", "colId", "checkId")
  for (factorCol in factorCols) {
    if (factorCol %in% cols && !is.factor(data[, which(cols == factorCol)])) {
      data[, which(cols == factorCol)] <-
        as.factor(data[, which(cols == factorCol)])
    }
  }
  ## Convert columns to numeric if neccessary.
  numCols <- c("rowCoord", "colCoord")
  for (numCol in numCols) {
    if (numCol %in% cols && !is.numeric(data[, which(cols == numCol)])) {
      data[, which(cols == numCol)] <-
        as.numeric(data[, which(cols == numCol)])
    }
  }
  if ("trial" %in% colnames(data)) {
    listData <- split(x = data, f = droplevels(data$trial))
  } else {
    listData <- setNames(list(data), dataName)
  }
  ## Define meta data to set from input variables.
  meta <- c("trLocation", "trDate", "trDesign", "trLat", "trLong",
            "trPlotWidth", "trPlotLength")
  ## Set meta for all trials in data.
  for (i in seq_along(listData)) {
    for (m in meta) {
      ## Set meta data. Set to NULL if not in input so meta variable is
      ## left out meta data. This to avoid a list of NULL.
      attr(x = listData[[i]], which = m) <-
        if (!is.null(get(m))) get(m) else NULL
    }
    ## Location should always be filled since it is used in plot titles as
    ## well. Use trial name as default value.
    if (is.null(trLocation)) {
      attr(x = listData[[i]], which = "trLocation") <- names(listData)[i]
    }
    ## Add a list of columns that have been renamed as attribute to TD.
    attr(x = listData[[i]], which = "renamedCols") <-
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
                  megaEnv = NULL,
                  year = NULL,
                  repId = NULL,
                  subBlock = NULL,
                  rowId = NULL,
                  colId = NULL,
                  rowCoord = NULL,
                  colCoord = NULL,
                  checkId = NULL,
                  trLocation = NULL,
                  trDate = NULL,
                  trDesign = NULL,
                  trLat = NULL,
                  trLong = NULL,
                  trPlotWidth = NULL,
                  trPlotLength = NULL) {
  TDNw <- createTD(data = data, genotype = genotype, trial = trial,
                   megaEnv = megaEnv, year = year, repId = repId,
                   subBlock = subBlock, rowId = rowId, colId = colId,
                   rowCoord = rowCoord,
                   colCoord = colCoord, checkId = checkId,
                   trLocation = trLocation, trDate = trDate,
                   trDesign = trDesign, trLat = trLat, trLong = trLong,
                   trPlotWidth = trPlotWidth, trPlotLength = trPlotLength)
  dupTrials <- names(TDNw)[names(TDNw) %in% names(TD)]
  if (length(dupTrials) > 0) {
    warning(paste0("The following trials already existed in TD and will be ",
                   "added again: ", paste(dupTrials, collapse = ", "), ".\n"),
            call. = FALSE)
  }
  return(c(TD, TDNw))
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
    stop(paste0(trial, " should be a single character string in ",
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
  ## Create a data.frame to store the values".
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
#' different trials within the TD object or locates the trials on a map.
#' Mapping the trials is done based on lattitude and longitude that can be
#' added when creating an object of class TD. The countries in which the trials
#' are located will be plotted on a single map and the location of the trials
#' will be indicated on this map.
#'
#' @param x An object of class TD.
#' @param ... Further graphical parameters. Fully functional when
#' \code{plotType = "layout"}. For \code{plotType = "layout"} only parameters
#' effecting the plot title will be passed.
#' @param trials A character vector indicating the trials to be plotted.
#' @param plotType A character string indicating which plot should be made.
#' This can be "layout" for a plot of the field layout for the diffent trials
#' in the TD object. This is only possible if the data contains row and column
#' information. Alternatively, for \code{plotType = "map"} a plot will be made
#' depicting the trials on a country map. This is only possible if lattitude
#' and longitude of the trials are available.
#'
#' @export
plot.TD <- function(x,
                    ...,
                    trials = names(x),
                    plotType = c("layout", "map")) {
  ## Maps seems to change graphics parameters without resetting. Do so here.
  #op <- par(no.readonly = TRUE)
  #on.exit(par(op))
  if (!is.character(trials) || !all(trials %in% names(x))) {
    stop(paste0("All trials should be in ", deparse(x), ".\n"))
  }
  plotType <- match.arg(plotType)
  dotArgs <- list(...)
  if (plotType == "layout") {
    for (trial in trials) {
      trDat <- x[[trial]]
      if (!"rowCoord" %in% colnames(trDat)) {
        warning(paste0("rowCoord should be a column in ", trial, ".\n",
                       "Plot skipped."), call. = FALSE)
        break
      }
      if (!"colCoord" %in% colnames(trDat)) {
        warning(paste0("colCoord should be a column in ", trial, ".\n",
                       "Plot skipped."), call. = FALSE)
        break
      }
      trLoc <- attr(trDat, "trLocation")
      ylen <- attr(trDat, "trPlotLength")
      xlen <- attr(trDat, "trPlotWidth")
      ## Compute aspect for proper depiction of field size. If no information
      ## is available plots are assumed to be square.
      if (is.null(ylen) || is.null(xlen)) {
        aspect <- length(unique(trDat$rowCoord)) /
          length(unique(trDat$colCoord))
      } else {
        aspect <- ylen / xlen
      }
      ## Desplot uses lattice for plotting which doesn't plot within a loop.
      ## This is solved by using print.
      plotVar <- ifelse("subBlock" %in% colnames(trDat), "subBlock", "trial")
      plotArgs <- list(form = formula(paste(plotVar, "~ colCoord + rowCoord")),
                       out1 = if ("repId" %in% colnames(trDat)) "repId" else NULL,
                       data = trDat, ticks = TRUE, main = trLoc,
                       aspect = aspect)
      ## Add and overwrite args with custom args from ...
      fixedArgs <- c("form", "data")
      plotArgs <- modifyList(plotArgs, dotArgs[!names(dotArgs) %in% fixedArgs])
      print(do.call(desplot::desplot, args = plotArgs))
    }
  } else if (plotType == "map") {
    ## Create a data.frame for plotting trials.
    ## Population has a random value but if left out nothing is plotted.
    locs <- setNames(getMeta(x)[c("trLocation", "trLat", "trLong")],
                     c("name", "lat", "long"))
    locs <- locs[!is.na(locs$lat) && !is.na(locs$long), ]
    if (nrow(locs) == 0) {
      stop(paste("At leaste one trial should have latitute and longitude",
                 "for plotting on map.\n"))
    } else {
      locs[c("capital", "pop")] <- rep(c(0, 10000), each = nrow(locs))
      ## Use lattitude and longitude to extract trial regions.
      regions <- unique(maps::map.where(x = locs$long, y = locs$lat))
      maps::map(regions = regions)
      maps::map.scale(relwidth = .15, ratio = FALSE, cex = .5)
      maps::map.cities(x = locs, col = seq_along(locs))
      plotArgs <- list(main = "Trial locations")
      plotArgs <- modifyList(plotArgs, dotArgs)
      do.call(title, args = plotArgs)
    }
  }
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
                "trPlotWidth", "trPlotLength")
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
#' trLong, trPlotWidth and trPlotLength. The values of the metadata of TD
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
                "trPlotWidth", "trPlotLength")
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
                        trPlotWidth = NULL,
                        trPlotLength = NULL) {
  if (!is.null(trDesign)) {
    trDesign <- match.arg(trDesign, choices = c("ibd", "res.ibd", "rcbd",
                                                "rowcol", "res.rowcol"))
  }
  if (!is.null(trLat) && (!is.numeric(trLat) || length(trLat) > 1 ||
                          abs(trLat) > 90)) {
    stop("trLat should be a single numerical value between -90 and 90.\n",
         call. = FALSE)
  }
  if (!is.null(trLong) && (!is.numeric(trLong) || length(trLong) > 1 ||
                           abs(trLong) > 180)) {
    stop("trLat should be a single numerical value between -180 and 180.\n",
         call. = FALSE)
  }
  if (!is.null(trLat) && !is.null(trLong)) {
    ## Check that coordinates point to a proper location so plotting can be done.
    loc <- maps::map.where(x = trLong, y = trLat)
    if (length(loc) > 0 && is.na(loc)) {
      warning("Values for trLat and trLong don't match a known land location.\n",
              call. = FALSE)
    }
  }
  if (!is.null(trPlotWidth) && (!is.numeric(trPlotWidth) ||
                                length(trPlotWidth) > 1 || trPlotWidth < 0)) {
    stop("trPlotWidth should be a single positive numerical value.\n",
         call. = FALSE)
  }
  if (!is.null(trPlotLength) && (!is.numeric(trPlotLength) ||
                                 length(trPlotLength) > 1 || trPlotLength < 0)) {
    stop("trPlotLength should be a single positive numerical value.\n",
         call. = FALSE)
  }
}



