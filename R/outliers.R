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
#'
#' @return A data.frame containing logical values indicating if the
#' observation is an outlier.
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
                       commonFactors = NULL) {
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
                      what = "stdRes")[[trial]]
  rDf <- STExtract(SSA, trials = trial, traits = traits,
                   what = "rDf")[[trial]]
  ## Create empty data.frame for storing results.
  indicator <- data.frame(matrix(data = FALSE,
                                 nrow = nrow(SSA[[trial]]$TD[[trial]]),
                                 ncol = length(traits),
                                 dimnames = list(NULL, traits)))
  outTrait <- vector(mode = "list", length = length(traits)) %>%
    setNames(traits)
  for (trait in traits) {
    ## Set commonFactors to trait when empty so joining is always possible.
    if (is.null(commonFactors)) {
      commonFactors <- trait
    }
    ## Compute limit value for residuals.
    if (is.null(rLimit)) {
      rLimit <- min(max(2, qnorm(p = 1 - 0.5 / rDf[trait])), 4)
    }
    ## Compute outliers.
    outVals <- stdRes[abs(stdRes[[trait]]) > rLimit, trait]
    if (length(outVals > 0)) {
      ## Fill indicator column for current trait.
      indicator[[trait]] <- stdRes[[trait]] %in% outVals
      ## Rename column for easier joining.
      stdRes <- dplyr::rename(stdRes, "res" = !!trait)
      ## Create data.frame with outliers for current trait.
      outTrait[[trait]] <- dplyr::inner_join(SSA[[trial]]$TD[[trial]],
                                             stdRes,
                                             by = setdiff(colnames(stdRes),
                                                          "res")) %>%
        dplyr::semi_join(.[.[["res"]] %in% outVals, commonFactors,
                           drop = FALSE], by = commonFactors) %>%
        ## Add column for similar and column trait with the current trait.
        dplyr::mutate(similar = !.[["res"]] %in% outVals, trait = trait) %>%
        ## Rename column trait to "value"
        dplyr::rename("value" = !!trait) %>%
        ## Change order of columns so always display trait, value and res first.
        dplyr::select(!!rlang::sym("trait"), !!rlang::sym("value"),
                      !!rlang::sym("res"), dplyr::everything())
    }
  }
  ## Create one single outlier matrix.
  pMat <- Reduce(f = rbind, x = outTrait)
  if (!is.null(pMat)) {
    cat(paste("Large standardized residuals\n\n"))
    print(format(pMat, quote = FALSE))
  } else {
    cat(paste("No large standardized residuals.\n"))
  }
  invisible(indicator)
}


