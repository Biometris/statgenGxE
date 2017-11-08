#' S3 class SSA
#'
#' Function for creating objects of S3 class Single Site Analysis (SSA).
#'
#' @param mMix a mixed model created using either asreml or lme4
#' @param mFix a fixed model created using either asreml or lme4
#' @param data a data.frame containing the data on which mMix and mFix are based.
#' @param x \code{R} object
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{summary.SSA}}
#'
#' @name SSA
NULL

#' @rdname SSA
#' @export
createSSA <- function(mMix,
                      mFix,
                      data) {
  SSA <- structure(list(mMix = mMix,
                        mFix = mFix,
                        data = data),
                   class = "SSA")
  return(SSA)
}

#' @rdname SSA
#' @export
is.SSA <- function(x) {
  inherits(x, "SSA")
}
