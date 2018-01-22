#' Identifying Outliers
#'
#' Function to identify observations that exceed \code{coef} times the
#' interquartile range. A summary is printed of outliers and observations
#' that have the same value for \code{commonFactors}. The latter ones will
#' be marked as similar to distinguish them from the former ones.
#'
#' @param TD An object of class \code{\link{TD}}.
#' @param traits A character vector specifying the name(s) of the traits for
#' which outliers should be identified.
#' @param coef A numerical value used for determining when a value is
#' considered an outlier. All observations that exceed \code{coef} times
#' the interquartile range will be marked as outliers.
#' @param commonFactors A character vector specifying the names of columns
#' in \code{TD} used for selecting observations that are similar to the
#' outliers. If \code{commonFactors = NULL} only outliers are reported and
#' no similar observations.
#'
#' @return A data.frame containing logical values indicating if the
#' observation is an outlier.
#'
#' @examples
#' outliers <- outlierTD(TD = TDHeat05, traits = "yield")
#'
#' @export
outlierTD <- function(TD,
                      traits,
                      coef = 1.5,
                      commonFactors = NULL) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (is.null(traits) || !is.character(traits) || !all(traits %in% colnames(TD))) {
    stop("traits has to be a vector of columns in TD.\n")
  }
  if (!is.null(commonFactors) && !is.character(commonFactors) &&
      !all(commonFactors %in% colnames(TD))) {
    stop("commonFactor has to be a vector of columns in TD.\n")
  }
  if (is.null(coef) || !is.numeric(coef) || length(coef) > 1 || coef < 0) {
    stop("coef should be a positive numerical value.\n")
  }
  ## Create empty data.frame for storing results.
  indicator <- data.frame(matrix(data = FALSE, nrow = nrow(TD),
                                 ncol = length(traits),
                                 dimnames = list(NULL, traits)))
  outTrait <- vector(mode = "list", length = length(traits)) %>%
    setNames(traits)
  for (trait in traits) {
    ## Set commonFactors to trait when empty so joining is always possible.
    if (is.null(commonFactors)) {
      commonFactors <- trait
    }
    ## Compute outliers.
    outVals <- boxplot.stats(x = TD[[trait]], coef = coef)$out
    if (length(outVals > 0)) {
      ## Fill indicator column for current trait.
      indicator[[trait]] <- TD[[trait]] %in% outVals
      ## Create data.frame with outliers for current trait.
      outTrait[[trait]] <- dplyr::semi_join(TD, TD[TD[[trait]] %in% outVals,
                                                   commonFactors, drop = FALSE],
                                            by = commonFactors) %>%
        ## Add column for similar and column trait with the current trait.
        dplyr::mutate(similar = !.[[trait]] %in% outVals, trait = trait) %>%
        ## Rename column trait to "value"
        dplyr::rename("value" = !!trait) %>%
        ## Change order of columns so always display trait and value first.
        dplyr::select(!!rlang::sym("trait"), !!rlang::sym("value"),
                      dplyr::everything())
    }
  }
  ## Create one single outlier matrix.
  pMat <- Reduce(f = rbind, x = outTrait)
  if (!is.null(pMat)) {
    cat(paste("Observations that exceed", coef, "times the interquartile range\n\n"))
    print(format(pMat, quote = FALSE))
  } else {
    cat(paste("No observations that exceed", coef, "times the interquartile range.\n"))
  }
  invisible(indicator)
}

#' Identifying Outliers
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
#' @inheritParams outlierTD
#'
#' @param SSA An object of class \code{\link{SSA}}.
#' @param rLimit A numerical value used for determining when a value is
#' considered an outlier. All observations with standardized residuals
#' exceeding \code{rLimit} will be marked as outliers.
#'
#' @return A data.frame containing logical values indicating if the
#' observation is an outlier.
#'
#' @examples
#' myModel <- STRunModel(TD = TDHeat05, traits = "yield" ,
#'                       design = "res.rowcol", engine = "lme4")
#' outliers <- outlierSSA(SSA = myModel, traits = "yield")
#'
#' @export
outlierSSA <- function(SSA,
                       traits,
                       rLimit = NULL,
                       commonFactors = NULL) {
  ## Checks.
  if (missing(SSA) || !inherits(SSA, "SSA")) {
    stop("SSA should be a valid object of class SSA.\n")
  }
  if (is.null(traits) || !is.character(traits) || !all(traits %in% colnames(SSA$TD))) {
    stop("traits has to be a vector of columns in TD.\n")
  }
  if (!is.null(commonFactors) && !is.character(commonFactors) &&
      !all(commonFactors %in% colnames(SSA$TD))) {
    stop("commonFactor has to be a vector of columns in TD.\n")
  }
  if (!is.null(rLimit) && (!is.numeric(rLimit) || length(rLimit) > 1 ||
                           rLimit < 0)) {
    stop("rLimit should be NULL or a positive numerical value.\n")
  }
  stdRes <- STExtract(SSA, traits = traits, what = "stdRes")
  rDf <- STExtract(SSA, traits = traits, what = "rDf")
  ## Create empty data.frame for storing results.
  indicator <- data.frame(matrix(data = FALSE, nrow = nrow(SSA$TD),
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
      outTrait[[trait]] <- dplyr::inner_join(SSA$TD, stdRes,
                                             by = setdiff(colnames(stdRes), "res")) %>%
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


