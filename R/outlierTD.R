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
#' outliers.
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
                      commonFactors = "genotype") {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (is.null(traits) || !is.character(traits) || !all(traits %in% colnames(TD))) {
    stop("trait has to be a vector of columns in TD.\n")
  }
  if (is.null(commonFactors) || !is.character(commonFactors) ||
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
    ## Compute outliers.
    outVals <- boxplot.stats(x = TD[[trait]], coef = coef)$out
    if (length(outVals > 0)) {
      ## Fill indicater column for current trait.
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
