#' Calculates stability coefficients for genotype-by-environment data
#'
#' This function calculate difference measures of stability, such as the cultivar-superiority measure of Lin & Binns (1988),
#' Shukla's (1972) stability variance and Wricke's (1962) ecovalence.
#'
#' @param Y A data frame object.
#' @param trait A character string specifying a trait column of the data.
#' @param genotype A character string specifying a genotype column of the data.
#' @param env A character string specifying the environment column of the data.
#' @param method A character string specifying (a) measure(s) of stability. By default, \code{method = c("superiority","static","wricke")}.
#' @param superiorityBestMethod A character string specifying the criterion to define the best genotype ("max","min").
#' By default, \code{superiorityBestMethod="max"}.
#' @param sorted A character string specifying whether the results are sorted by increasing (or decreasing) order of stability coefficients.
#' By default, \code{sortBYsens = "descending"}. Other options are \code{"ascending"} and \code{NA}.
#' @param plot A logical value indicating whether to produce a plot of stability against the mean.
#' @return A list of measures of stability.
#' @examples
#' mydat <- GE.read.csv(file.path(path.package("RAP"),"F2maize_pheno.csv"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' names(mydat)=c("env", "genotype","yld")
#' GE.stability(Y=mydat, trait="yld", genotype="genotype", env="env",
#'              method = "superiority", superiorityBestMethod = "max",
#'              sorted = "descending", plot=TRUE)
#'
#' @export

GE.stability <- function(Y,
                         trait,
                         genotype,
                         env,
                         method = c("superiority", "static", "wricke"),
                         superiorityBestMethod = "max",
                         sorted = c("ascending", "descending", NA),
                         plot = TRUE) {
  if (missing(sorted)) {
    sorted <- "descending"
  }
  if (any(is.na(Y[[trait]]))) {
    y0 <- tapply(Y[[trait]], Y[,c(genotype,env)], mean)
    yIndex <- tapply(X = 1:nrow(Y), INDEX = Y[, c(genotype, env)], FUN = identity)
    na_yes_no <- is.na(y0)
    # imputaion
    y1 <- RAP.multmissing(y0, maxcycle = 10, na.strings = NA)
    replaceVal <- y1[na_yes_no]
    yIndexReplace <- yIndex[na_yes_no]
    if (is.list(yIndexReplace)) {
      nyr <- length(yIndexReplace)
      for (ii in 1:nyr) {
        for (jj in 1:length(Y[yIndexReplace[[ii]], trait])) {
          if (is.na(Y[yIndexReplace[[ii]][jj], trait]))
            Y[yIndexReplace[[ii]][jj], trait] <- replaceVal[ii]
        }
      }
    } else {
      Y[yIndexReplace, trait] <- replaceVal
    }
  }
  lab <- levels(Y[,genotype])
  nLab <- length(lab)
  envs <- levels(Y[, env])
  nEnvs <- length(envs)
  # The centered trait mean for eniroment j
  Ej <- sapply(X = envs, FUN = function(x) {
    mean(Y[which(Y[, env] == x), trait], na.rm = TRUE)
  })
  # Maximum or minimum trait mean among all genotype in jth enviroment
  if (superiorityBestMethod == "max") {
    Mj <- sapply(X = envs, FUN = function(x) {
      max(Y[which(Y[, env] == x), trait], na.rm = TRUE)
    })
  } else {
    Mj <- sapply(X = envs, FUN = function(x) {
      min(Y[which(Y[, env] == x), trait], na.rm = TRUE)
    })
  }
  # The genotype trait mean across environments
  Ei <- sapply(X = lab, FUN = function(x) {
    mean(Y[which(Y[, genotype] == x), trait], na.rm = TRUE)
  })
  # The grand mean
  E <- mean(Y[, trait], na.rm = TRUE)
  W <- S <- LB <- rep(NA, nLab)
  for (i in 1:nLab) {
    # Observed genotype field response in the enviroment j
    # (averaged across experiment replicates)
    Rij <- sapply(X = envs, FUN = function(x) {
      mean(Y[which(Y[, genotype] == lab[i] & Y[, env] == x), trait], na.rm = TRUE)
    })
    pos <- (1:nEnvs)[!is.na(Rij)]
    # Static measure (Shukla's (1972a) stability variance)
    if ("static" %in% method) {
      if (length(pos) == 0) {
        S[i] <- NA
      } else {
        S[i] <- sum(sapply(X = pos, FUN = function(j) {
          (Rij[j] - Ei[i]) ^ 2
        }) / (nEnvs - 1))
      }
    }
    # Superiority measure (LIN&BINNS 1988)
    if ("superiority" %in% method) {
      if (length(pos) == 0){
        LB[i] <- NA
      } else{
        LB[i] <- sum(sapply(X = pos, FUN = function(j) {
          (Rij[j] - Mj[j]) ^ 2
        }) / (2 * nEnvs))
      }
    }
    # Wricke's (1962) ecovalence
    if ("wricke" %in% method){
      if (length(pos)==0) {
        W[i] <- NA
      } else {
        W[i] <- sum(sapply(X = pos, FUN = function(j) {
          (Rij[j] - Ei[i] - Ej[j] + E) ^ 2
        }))
      }
    }
  }
  #Calculate trait mean per genotype
  means <- sapply(X = lab, FUN = function(x) {
    mean(Y[Y[, genotype] == x, trait], na.rm = TRUE)
  })
  if (plot) {
    # Prepare panels
    mFlag <- c("static", "superiority", "wricke") %in% method
    if (sum(mFlag) !=1) {
      if (sum(mFlag) == 3) {
        opar <- par(mfrow = c(2, 2), mar = c(4, 4, 4, 2), oma = c(.5, .5, 2, .3))
      } else{
        opar <- par(mfrow=c(1, 2), mar = c(4, 4, 4, 2), oma = c(.5, .5, 2, .3))
      }
    } else{
      opar <- par(mfrow = c(1, 1))
    }
    on.exit(par(opar))
    if ("superiority" %in% method) {
      plot(x = means, y = LB, xlab = "Mean", ylab = "Cultivar superiority")
    }
    if ("static" %in% method) {
      plot(x = means, y = S, xlab = "Mean", ylab = "Static stability")
    }
    if ("wricke" %in% method) {
      plot(x = means, y = W, xlab = "Mean", ylab = "Wricke's ecovalence")
    }
    if (sum(mFlag) == 1){
      title(paste('Stability coefficients for', trait))
    } else{
      mtext(paste('Stability coefficients for', trait), side = 3, outer = TRUE,
            cex = 1.3, font = 2)
    }
  }
  if (is.na(sorted)) {
    if ("static" %in% method) {
      SOut <- data.frame(lab, S, means, row.names = 1:nLab)
      names(SOut) <- c("genotype","static", "mean")
    }
    if ("superiority" %in% method) {
      LBOut <- data.frame(lab, LB, means, row.names = 1:nLab)
      names(LBOut) <- c("genotype", "superiority", "mean")
    }
    if ("wricke" %in% method) {
      WOut <- data.frame(lab,W,means, row.names=1:nLab)
      names(WOut) <- c("genotype", "wricke", "mean")
    }
  } else {
    if (sorted == "ascending"){
      if ("static" %in% method){
        orderS <- order(S)
        SOut <- data.frame(lab, S, means, row.names = 1:nLab)[orderS, ]
        names(SOut) <- c("genotype", "static", "mean")
      }
      if ("superiority" %in% method){
        orderLB <- order(LB)
        LBOut <- data.frame(lab, LB, means, row.names = 1:nLab)[orderLB, ]
        names(LBOut) <- c("genotype", "superiority", "mean")
      }
      if ("wricke" %in% method){
        orderW <- order(W)
        WOut <- data.frame(lab,W,means, row.names = 1:nLab)[orderW, ]
        names(WOut) <- c("genotype", "wricke", "mean")
      }
    } else if (sorted == "descending") {
      if ("static" %in% method) {
        orderS <- order(S, decreasing = TRUE)
        SOut <- data.frame(lab, S, means, row.names = 1:nLab)[orderS, ]
        names(SOut) <- c("genotype", "static", "mean")
      }
      if ("superiority" %in% method){
        orderLB <- order(LB, decreasing = TRUE)
        LBOut <- data.frame(lab, LB, means, row.names = 1:nLab)[orderLB, ]
        names(LBOut) <- c("genotype", "superiority", "mean")
      }
      if ("wricke" %in% method){
        orderW <- order(W, decreasing = TRUE)
        WOut <- data.frame(lab, W, means, row.names = 1:nLab)[orderW, ]
        names(WOut) <- c("genotype", "wricke", "mean")
      }
    } else {
      stop('options for sorted must be one of "ascending", "descending" and NA')
    }
  }
  res <- vector(mode = "list")
  if ("static" %in% method) res$static <- SOut
  if ("superiority" %in% method) res$superiority <- LBOut
  if ("wricke" %in% method) res$wricke <- WOut
  return(res)
}
