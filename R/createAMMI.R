#' S3 class AMMI
#'
#' Function for creating objects of S3 class AMMI.
#'
#' \code{\link{print}} and \code{\link{plot}} methods are available.
#'
#' @param envScores a matrix containing environment scores
#' @param genoScores a matrix containing genotypic scores
#' @param importance a data.frame containing the importance of the principal components
#' @param anova a data.frame containing anova scores of the AMMI analysis
#' @param fitted a matrix containing fitted values from the AMMI model
#' @param trait a character value indicating the analysed trait
#' @param envMean a numerical vector containing the means per environment
#' @param genoMean a numerical vector containing the means per genotype
#' @param overallMean a numberical value containing the overall mean
#' @param x an \code{R} object
#'
#' @author Bart-Jan van Rossum
#'
#' @seealso \code{\link{plot.AMMI}}
#'
#' @name AMMI
NULL

#' @rdname AMMI
#' @export
createAMMI <- function(envScores,
                       genoScores,
                       importance,
                       anova,
                       fitted,
                       trait,
                       envMean,
                       genoMean,
                       overallMean) {
  AMMI <- structure(list(envScores = envScores,
                         genoScores = genoScores,
                         importance = importance,
                         anova = anova,
                         fitted = fitted,
                         trait = trait,
                         envMean = envMean,
                         genoMean = genoMean,
                         overallMean = overallMean),
                    class = "AMMI")
  return(AMMI)
}

#' @rdname AMMI
#' @export
is.AMMI <- function(x) {
  inherits(x, "AMMI")
}

#' @export
print.AMMI <- function(x, ...) {
  cat("Principal components",
      "\n====================\n")
  print(x$importance)
  cat("\nAnova",
      "\n=====\n")
  printCoefmat(x$anova)
  cat("\nEnvironment scores",
      "\n==================\n")
  print(x$envScores, ...)
  cat("\nGenotypic scores",
      "\n================\n")
  print(x$genoScores, ..., max.print = 50)
}


