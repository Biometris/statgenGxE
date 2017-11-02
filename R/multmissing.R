#' Multmissing procedure
#'
#' This function is to estimate missing values for units in a multivariate data set, using an iterative regression technique.
#'
#' @param Y A data frame or matrix object of multivariate data.
#' @param maxcycle An integer specifying the maximum allowed number of iterations; default 10.
#' @param na.strings A character vector of strings which are to be interpreted as \code{NA} values.
#' @details Initial estimates of the missing values in each variate are formed from the variate means using the values for units that have no missing values for any variate. Estimates of the missing values for each variate are then recalculated as the fitted values from the multiple regression of that variate on all the other variates. When all the missing values have been estimated the variate means are recalculated. If any of the means differs from the previous mean by more than a tolerance (the initial standard error divided by 1000) the process is repeated, subject to a maximum number of repetitions defined by the MAXCYCLE option. The default maximum number of iterations (10) is usually sufficient when there are few missing values, say two or three. If there are many more, 20 or so, it may be necessary to increase the maximum number of iterations to around 30. The method is similar to that of Orchard & Woodbury (1972), but does not adjust for bias in the variance-covariance matrix as suggested by Beale & Little (1975).
#' @return A data.frame object of multivarite data with the missing values replaced by their estimates.
#' @references (ed.) R.W. Payne (2011). GenStat Release 14 Reference Manual, Part 3 Procedure Library PL20. VSN International, Hemel Hempstead, UK.
#'
#' Beale, E.M.L. & Little, R.J.A. (1975). Missing values in multivariate analysis. Journal of the Royal Statistical Society, Series B, 37, 129-145.
#'
#' Orchard, T. & Woodbury, M.A. (1972). A missing information principle: theory and applications. In: Proceedings of the 6th Berkeley Symposium in Mathematical Statistics and Probability, Vol I, 697-715.
#'
#' @export
RAP.multmissing <- function(Y, maxcycle = 10, na.strings=NA){
  ## TODO: ?weights update procedure ?adjust for bias in the variance-covariance matrix ?logistic regression
  
  if(is.data.frame(Y)) {
    Cnames <- names(Y)
    Rnames <- rownames(Y)
  }else{
    if(is.matrix(Y)){
      Cnames <- colnames(Y)
      Rnames <- rownames(Y)
      if(is.null(Cnames)) Cnames <- colnames(Y) <- paste("V",1:ncol(Y),sep="")
    }else{
      if(!is.vector(Y)) stop("'Y' must be a data.frame, matrix or vector object.\n")
    }
  }
  # To avoid that column name can be converted to integers
  Cnames0 <- Cnames
  Cnames <- paste("V",1:ncol(Y),sep="")
  Rnames0 <- Rnames
  
  if(!all(is.na(na.strings))){
    for (i in 1:length(na.strings)){
      tmp <- which(Y==na.strings[i])
      Y[tmp] <- NA
    }
  }
  if (is.vector(Y)){
    Y <- as.numeric(Y)
    miss.ind <- which(is.na(Y))
    Y[miss.ind] <- mean(Y, na.rm=TRUE)
  }else{
    Y <- apply(Y, 2, as.numeric)
    allNA <- apply(Y, 2, function(x) all(is.na(x)))
    if (sum(allNA)!=0) stop("At least one unit must have no missing values.\n")
    d1 <- nrow(Y)
    d2 <- ncol(Y)
    
    # missing positions (row, col) [nmiss x 2]
    miss.ind <- which(is.na(Y), arr.ind=T)
    # weights
    W <- matrix(1,d1,d2)
    W[miss.ind] <- 0
    if (nrow(miss.ind)==0){
      warning("no missing value found \n")
    }else{
      miss.col.loc <- miss.ind[,2]
      miss.row.loc <- miss.ind[,1]
      nmiss <- nrow(miss.ind)
      miss.variate.ind <- unique(miss.col.loc)
          
      ## initializing...
      # replaces missing values by the means of
      # non-missing values of each variate
      for (j in miss.variate.ind){
        tmp <- which(is.na(Y[,j]))
        Y[tmp,j] <- mean(Y[,j], na.rm=TRUE)
      }
      colnames(Y) <- Cnames
      estprev <- rep(Inf, d2)
      estcurr <- colMeans(Y)
      tol <- sd(as.vector(Y))/1000
      
      ## iterative regression
      kk <- 0
      while (max(abs(estcurr-estprev)-tol)>0 &&  kk<=maxcycle){
        kk <- kk+1
        for (j in miss.variate.ind){
          formula <- as.formula(paste(Cnames[j], "~", paste(Cnames[-j],collapse="+")))
          model <- lm(formula=formula, data=as.data.frame(Y), weights=W[,j])
          fvals <- model$fitted.values
          tmp <- which(miss.col.loc==j)
          Y[miss.row.loc[tmp],j] <- fvals[miss.row.loc[tmp]]
        }
        # Updates
        estprev <- estcurr
        estcurr <- colMeans(Y)
      }
    }
  }
  colnames(Y) <- Cnames0
  rownames(Y) <- Rnames0
  Y
}
