#' Identifying outliers in objects of class SSA
#'
#' Function to identify observations with standardized residuals exceeding
#' \code{rLimit}. If not provided \code{rLimit} is computed as
#' \code{qnorm(1 - 0.5 / rDf)} where \code{rDf} are the residual degrees
#' of freedom for the model. This value is then restricted to the interval
#' 2..4. Alternatively a custom limit may be provided.\cr
#' A summary is printed of outliers and observations that have the same
#' value for \code{commonFactors}. The latter ones will be marked as
#' similar to distinguish them from the former ones.
#'
#' @param SSA An object of class \code{\link{SSA}}.
#' @param trial A character string specifying the trial for which outliers
#' should be identified. If \code{NULL} and \code{SSA} contains only one trial
#' that trial is used.
#' @param traits A character vector specifying the names of the traits for
#' which outliers should be identified.
#' @param rLimit A numerical value used for determining when a value is
#' considered an outlier. All observations with standardized residuals
#' exceeding \code{rLimit} will be marked as outliers.
#' @param commonFactors A character vector specifying the names of columns
#' in \code{TD} used for selecting observations that are similar to the
#' outliers. If \code{commonFactors = NULL} only outliers are reported and
#' no similar observations.
#' @param verbose Should the outliers be printed to the console?
#'
#' @return A list with two components:
#' \itemize{
#' \item{indicator - a data.frame containing logical values indicating if the
#' observation is an outlier.}
#' \item{outliers - a data.frame containing the outliers and observations
#' similar to the outliers as defined by \code{commonFactors}}
#' }
#'
#' @examples
#' ## Fit a model using lme4.
#' myModel <- STRunModel(TD = TDHeat05, traits = "yield" ,
#'                       design = "res.rowcol", engine = "lme4")
#' ## Detect outliers in the standardized residuals of the fitted model.
#' outliers <- outlierSSA(SSA = myModel, traits = "yield")
#'
#' @export
outlierSSA <- function(SSA,
                       trial = NULL,
                       traits,
                       rLimit = NULL,
                       commonFactors = NULL,
                       verbose = TRUE) {
  ## Checks.
  if (missing(SSA) || !inherits(SSA, "SSA")) {
    stop("SSA should be a valid object of class SSA.\n")
  }
  if (is.null(trial) && length(SSA) > 1) {
    stop("No trial provided but multiple trials found in SSA object.\n")
  }
  if (!is.null(trial) && (!is.character(trial) || length(trial) > 1 ||
                          !trial %in% names(SSA))) {
    stop("Trial has to be a single character string defining a trial in SSA.\n")
  }
  if (is.null(trial)) {
    trial <- names(SSA)
  }
  if (!is.null(traits) && (!is.character(traits) ||
                           !all(traits %in% colnames(SSA[[trial]]$TD[[trial]])))) {
    stop("Trait has to be a character vector defining columns in TD.\n")
  }
  if (!is.null(commonFactors) && !is.character(commonFactors) &&
      !all(commonFactors %in% colnames(SSA[[trial]]$TD[[trial]]))) {
    stop("commonFactor has to be a character vector defining columns in TD.\n")
  }
  if (!is.null(rLimit) && (!is.numeric(rLimit) || length(rLimit) > 1 ||
                           rLimit < 0)) {
    stop("rLimit should be NULL or a positive numerical value.\n")
  }
  stdRes <- STExtract(SSA, trials = trial, traits = traits,
                      what = "stdRes")[[trial]][["stdRes"]]
  rDf <- STExtract(SSA, trials = trial, traits = traits,
                   what = "rDf")[[trial]][["rDf"]]
  ## Create empty data.frame for storing results.
  indicator <- data.frame(matrix(data = FALSE,
                                 nrow = nrow(SSA[[trial]]$TD[[trial]]),
                                 ncol = length(traits),
                                 dimnames = list(NULL, traits)))
  outTrait <- setNames(vector(mode = "list", length = length(traits)), traits)
  for (trait in traits) {
    ## Compute limit value for residuals.
    if (is.null(rLimit)) {
      rLimit <- min(max(2, qnorm(p = 1 - 0.5 / rDf[trait])), 4)
    }
    datTr <- SSA[[trial]]$TD[[trial]]
    datTr <- datTr[!colnames(datTr) %in% setdiff(traits, trait)]
    ## Compute outliers.
    ## Set missing values to 0 to prevent problems when comparing to rLimit.
    stdRes[is.na(stdRes[[trait]]), trait] <- 0
    outVals <- stdRes[abs(stdRes[[trait]]) > rLimit, trait]
    if (length(outVals > 0)) {
      ## Fill indicator column for current trait.
      indicator[[trait]] <- abs(stdRes[[trait]]) > rLimit
      ## Rename column for easier joining.
      colnames(stdRes)[colnames(stdRes) == trait] <- "res"
      ## Create data.frame with outliers for current trait.
      outTr <- cbind(datTr, stdRes["res"])
      if (!is.null(commonFactors)) {
        ## If commonFactors are given merge to data.
        outTr <- unique(merge(x = outTr,
                              y = outTr[abs(outTr$res) > rLimit,
                                        commonFactors, drop = FALSE],
                              by = commonFactors))
      } else {
        outTr <- outTr[abs(outTr$res) > rLimit, ]
      }
      outTr$similar <- abs(outTr$res) <= rLimit
      outTr$trait <- trait
      ## Rename column trait to value.
      colnames(outTr)[colnames(outTr) == trait] <- "value"
      ## Change order of columns to always display trait, value and res first.
      outTr <- cbind(outTr[, c("trait", "value", "res")],
                     outTr[!colnames(outTr) %in% c("trait", "value", "res")])
      outTrait[[trait]] <- outTr
    }
  }
  ## Create one single outlier data.frame.
  pMat <- Reduce(f = rbind, x = outTrait)
  if (verbose) {
    if (!is.null(pMat)) {
      cat(paste("Large standardized residuals\n\n"))
      print(format(pMat, quote = FALSE))
    } else {
      cat("No large standardized residuals.\n")
    }
  }
  return(list(indicator = indicator, outliers = pMat))
}


