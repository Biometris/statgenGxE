#' S3 class SSA
#'
#' Function for creating objects of S3 class Single Site Analysis (SSA).
#'
#' @param mMix a mixed model created using either asreml or lme4
#' @param mFix a fixed model created using either asreml or lme4
#' @param data an object of class TD containing the data on which mMix and mFix are based.
#' @param trait a character sting indicating the trait for which the analysis is done.
#' @param genotype a character sting indicating genotype column in the data.
#' @param rep a character sting indicating the replicates column in the data.
#' @param design a character string containing the design of the trial.
#' @param engine a character string containing the engine used to do the analysis.
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
                      data,
                      trait = NULL,
                      genotype = NULL,
                      rep = NULL,
                      design = NULL,
                      engine = NULL) {
  SSA <- structure(list(mMix = mMix,
                        mFix = mFix,
                        data = data,
                        trait = trait,
                        genotype = genotype,
                        rep = rep,
                        design = design,
                        engine = engine),
                   class = "SSA")
  return(SSA)
}

#' @rdname SSA
#' @export
is.SSA <- function(x) {
  inherits(x, "SSA")
}
