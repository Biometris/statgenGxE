#' Multmissing procedure
#'
#' This function estimates missing values for units in a multivariate data set,
#' using an iterative regression technique.
#'
#' Initial estimates of the missing values in each variate are formed from the
#' variate means using the values for units that have no missing values for
#' any variate. Estimates of the missing values for each variate are then
#' recalculated as the fitted values from the multiple regression of that
#' variate on all the other variates. When all the missing values have been
#' estimated the variate means are recalculated. If any of the means differs
#' from the previous mean by more than a tolerance (the initial standard error
#' divided by 1000) the process is repeated, subject to a maximum number of
#' repetitions defined by \code{maxIter} option. The default maximum number
#' of iterations (10) is usually sufficient when there are few missing
#' values, say two or three. If there are many more, 20 or so, it may be
#' necessary to increase the maximum number of iterations to around 30. The
#' method is similar to that of Orchard & Woodbury (1972), but does not adjust
#' for bias in the variance-covariance matrix as suggested by Beale &
#' Little (1975).
#'
#' @param Y A matrix, data.frame or vector of multivariate data.
#' @param maxIter An integer specifying the maximum number of iterations.
#' @param naStrings A character vector of strings which are to be interpreted
#' as \code{NA} values.
#'
#' @return An object of the same class as the input \code{\link{Y}} with the
#' missing values replaced by their estimates.
#'
#' @references Beale, E.M.L. & Little, R.J.A. (1975). Missing values in
#' multivariate analysis. Journal of the Royal Statistical Society, Series B,
#' 37, 129-145.
#' @references Orchard, T. & Woodbury, M.A. (1972). A missing information principle:
#' theory and applications. In: Proceedings of the 6th Berkeley Symposium in
#' Mathematical Statistics and Probability, Vol I, 697-715.
#'
#' @examples
#' M <- matrix(c("1", "2", "3", NA, "b", "5", "6",
#'               "6", "5", "b", NA, "3", "2", "1"), nrow = 7, ncol = 2)
#' ## Estimate missing values treating "b" as NA.
#' multMissing(M, naStrings = "b")
#'
#' @export
multMissing <- function(Y,
                        maxIter = 10,
                        naStrings = NULL) {
  ## Checks.
  if (!is.matrix(Y) && !is.data.frame(Y) && !is.vector(Y)) {
    stop("Y should be a data.frame, matrix or vector.\n")
  }
  if (!is.numeric(maxIter) || length(maxIter) > 1 ||
      round(maxIter) != maxIter || maxIter < 0) {
    stop("maxIter should be an integer.\n")
  }
  if (!is.null(naStrings) && !is.character(naStrings)) {
    stop("naStrings should be NULL or a character vector.\n")
  }
  if (!is.vector(Y)) {
    inMat <- is.matrix(Y)
    if (is.data.frame(Y)) {
      Y <- as.matrix(Y)
    }
    ## Save row and column names to reassign in the end.
    cNames0 <- colnames(Y)
    rNames0 <- rownames(Y)
    ## Simplify column names for easier processing.
    cNames <- colnames(Y) <- paste0("V", 1:ncol(Y))
  }
  ## set all naStrings to NA.
  if (!is.null(naStrings)) {
    Y[Y %in% naStrings] <- NA
  }
  if (is.vector(Y)) {
    ## Replace missings by mean.
    Y <- as.numeric(Y)
    if (anyNA(Y)) {
      Y[is.na(Y)] <- mean(Y, na.rm = TRUE)
    } else {
      warning("no missing values found.\n", call. = FALSE)
    }
  } else {
    Y <- tryCatchExt(apply(X = Y, MARGIN = 2, FUN = as.numeric))
    if (!is.null(Y$warning)) {
      warning("Warning when converting data to numeric values. Invalid values
              were converted to NA. Make sure to check output.\n",
              call. = FALSE)
    }
    Y <- Y$value
    allNA <- apply(X = Y, MARGIN = 2, FUN = function(x) {
      all(is.na(x))
    })
    if (sum(allNA) != 0) {
      stop("At least one unit must have no missing values.\n")
    }
    ## missing positions (row, col)
    missInd <- which(is.na(Y), arr.ind = TRUE)
    if (nrow(missInd) == 0) {
      warning("no missing values found.\n", call. = FALSE)
    } else {
      ## Extract variates with missing values.
      missVariateInd <- unique(missInd[, 2])
      ## Replace missing values by the means of non-missing values per variable.
      for (i in missVariateInd) {
        Y[is.na(Y[, i]), i] <- mean(Y[, i], na.rm = TRUE)
      }
      ## Create a matrix of weights with values 0/1 for missings/non-missings.
      W <- matrix(data = 1, nrow = nrow(Y), ncol = ncol(Y))
      W[missInd] <- 0
      ## Set initial values for iterative process.
      maxDiff <- Inf
      tol <- sd(Y) / 1000
      iter <- 1
      Y <- as.data.frame(Y)
      while (maxDiff > tol && iter <= maxIter) {
        estPrev <- colMeans(Y)
        ## Estimate missing values.
        for (i in missVariateInd) {
          ## Estimate current column using all other columns.
          formula <- as.formula(paste(cNames[i], "~", paste(cNames[-i],
                                                            collapse = "+")))
          model <- lm(formula = formula, data = Y,
                      weights = W[, i])
          ## Replace original missing values with fitted values.
          fVals <- fitted(model)
          Y[missInd[missInd[, 2] == i, 1], i] <-
            fVals[missInd[missInd[, 2] == i, 1]]
        }
        ## Update values for iterative process.
        maxDiff <- max(abs(colMeans(Y) - estPrev))
        if (iter == maxIter && maxDiff > tol) {
          warning(paste("No convergence achieved after", iter," iterations.
                        Tolerance at last iteration", signif(maxDiff, 4),
                        ".\n"))
        }
        iter <- iter + 1
      }
      ## Reset row and column names to original values.
      colnames(Y) <- cNames0
      rownames(Y) <- rNames0
      if (inMat) {
        Y <- as.matrix(Y)
      }
    }
  }
  return(Y)
}
