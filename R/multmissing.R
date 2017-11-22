#' Multmissing procedure
#'
#' This function is to estimate missing values for units in a multivariate data set,
#' using an iterative regression technique.
#'
#' Initial estimates of the missing values in each variate are formed from the
#' variate means using the values for units that have no missing values for any variate.
#' Estimates of the missing values for each variate are then recalculated as the fitted values
#' from the multiple regression of that variate on all the other variates. When all the missing v
#' alues have been estimated the variate means are recalculated. If any of the means differs from
#' the previous mean by more than a tolerance (the initial standard error divided by 1000) the process
#' is repeated, subject to a maximum number of repetitions defined by the MAXCYCLE option.
#' The default maximum number of iterations (10) is usually sufficient when there are few missing
#' values, say two or three. If there are many more, 20 or so, it may be necessary to increase
#' the maximum number of iterations to around 30. The method is similar to that of
#' Orchard & Woodbury (1972), but does not adjust for bias in the variance-covariance matrix
#' as suggested by Beale & Little (1975).
#'
#' @param Y A data.frame or matrix object of multivariate data.
#' @param maxcycle An integer specifying the maximum allowed number of iterations; default 10.
#' @param na.strings A character vector of strings which are to be interpreted as \code{NA} values.
#'
#' @return A data.frame object of multivarite data with the missing values replaced by their estimates.
#' @references (ed.) R.W. Payne (2011). GenStat Release 14 Reference Manual, Part 3 Procedure
#' Library PL20. VSN International, Hemel Hempstead, UK.\cr
#' Beale, E.M.L. & Little, R.J.A. (1975). Missing values in multivariate analysis.
#' Journal of the Royal Statistical Society, Series B, 37, 129-145.\cr
#' Orchard, T. & Woodbury, M.A. (1972). A missing information principle: theory and applications.
#' In: Proceedings of the 6th Berkeley Symposium in Mathematical Statistics and Probability,
#' Vol I, 697-715.
#'
#' @export
RAP.multmissing <- function(Y,
                            maxcycle = 10,
                            na.strings = NA) {
  ## TODO: ?weights update procedure ?adjust for bias in the variance-covariance matrix
  ##?logistic regression
  if (is.data.frame(Y)) {
    cNames <- names(Y)
    rNames <- rownames(Y)
  } else if (is.matrix(Y)) {
    cNames <- colnames(Y)
    rNames <- rownames(Y)
    if (is.null(cNames)) {
      cNames <- colnames(Y) <- paste0("V", 1:ncol(Y))
    }
  } else if (!is.vector(Y)) {
    stop("'Y' must be a data.frame, matrix or vector object.\n")
  }
  # To avoid that column name can be converted to integers
  cNames0 <- cNames
  cNames <- paste0("V", 1:ncol(Y))
  if (!all(is.na(na.strings))) {
    for (i in 1:length(na.strings)) {
      tmp <- which(Y == na.strings[i])
      Y[tmp] <- NA
    }
  }
  if (is.vector(Y)) {
    Y <- as.numeric(Y)
    missInd <- which(is.na(Y))
    Y[missInd] <- mean(Y, na.rm = TRUE)
  } else {
    Y <- apply(X = Y, MARGIN = 2, FUN = as.numeric)
    allNA <- apply(X = Y, MARGIN = 2, FUN = function(x) {
      all(is.na(x))
    })
    if (sum(allNA) != 0) {
      stop("At least one unit must have no missing values.\n")
    }
    d1 <- nrow(Y)
    d2 <- ncol(Y)
    # missing positions (row, col) [nmiss x 2]
    missInd <- which(is.na(Y), arr.ind = TRUE)
    # weights
    W <- matrix(1, d1, d2)
    W[missInd] <- 0
    if (nrow(missInd) == 0) {
      warning("no missing values found.\n")
    } else {
      missColLoc <- missInd[, 2]
      missRowLoc <- missInd[, 1]
      nmiss <- nrow(missInd)
      missVariateInd <- unique(missColLoc)
      ## initializing...
      # replaces missing values by the means of
      # non-missing values of each variate
      for (j in missVariateInd) {
        tmp <- which(is.na(Y[, j]))
        Y[tmp, j] <- mean(Y[, j], na.rm = TRUE)
      }
      colnames(Y) <- cNames
      estPrev <- rep(Inf, d2)
      estCurr <- colMeans(Y)
      tol <- sd(as.vector(Y)) / 1000
      ## iterative regression
      kk <- 0
      while (max(abs(estCurr - estPrev) - tol) > 0 && kk <= maxcycle) {
        kk <- kk + 1
        for (j in missVariateInd) {
          formula <- as.formula(paste(cNames[j], "~", paste(cNames[-j], collapse = "+")))
          model <- lm(formula = formula, data = as.data.frame(Y), weights = W[, j])
          fVals <- model$fitted.values
          tmp <- which(missColLoc == j)
          Y[missRowLoc[tmp], j] <- fVals[missRowLoc[tmp]]
        }
        # Updates
        estPrev <- estCurr
        estCurr <- colMeans(Y)
      }
    }
  }
  colnames(Y) <- cNames0
  rownames(Y) <- rNames
  return(Y)
}
