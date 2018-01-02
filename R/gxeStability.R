#' Calculates stability coefficients for genotype-by-environment data
#'
#' This function calculate difference measures of stability, such as the cultivar-superiority
#' measure of Lin & Binns (1988), Shukla's (1972) stability variance and Wricke's (1962)
#' ecovalence.
#'
#' @inheritParams gxeAmmi
#'
#' @param method A character string specifying (a) measure(s) of stability. By default,
#' \code{method = c("superiority","static","wricke")}.
#' @param superiorityBestMethod A character string specifying the criterion to define the best
#' genotype ("max","min"). By default, \code{superiorityBestMethod = "max"}.
#' @param sorted A character string specifying whether the results are sorted by increasing
#' (or decreasing) order of stability coefficients. By default, \code{sortBYsens = "descending"}.
#' Other options are \code{"ascending"} and \code{NA}.
#' @param plot A logical value indicating whether to produce a plot of stability against the mean.
#'
#' @return A list of one to three data.frames containing the stability measures.
#'
#' @references LtN, C. S. aNt BINNS, M. R. 1988. A superiority measure of cultivar performance
#' for cultivar x location data. Can. J. Plant Sci. 68: 193-198\cr
#' Shukla, G.K. 1972. Some statistical aspects of partitioning genotype-environmental
#' components of variability. Heredity 29:237-245\cr
#' Wricke, G. Uber eine method zur erfassung der okologischen streubreit in feldversuchen.
#' Zeitschrift für Pflanzenzucht, v. 47, p. 92-96, 1962\cr
#'
#' @examples
#' gxeStability(TD = TDMaize, trait = "yld", method = "superiority",
#'              superiorityBestMethod = "max", sorted = "descending", plot = TRUE)
#'
#' @export

gxeStability <- function(TD,
                         trait,
                         method = c("superiority", "static", "wricke"),
                         superiorityBestMethod = "max",
                         sorted = c("ascending", "descending", NA),
                         plot = TRUE) {
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (!"env" %in% colnames(TD)) {
    stop("TD should contain a column env to be able to run an AMMI analysis.\n")
  }
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TD)) {
    stop("trait has to be a column in TD.\n")
  }
  if (is.null(method) || !is.character(method) ||
      !method %in% c("superiority", "static", "wricke")) {
    stop("method should be superiority, static and/or wricke.\n")
  }
  if (is.null(superiorityBestMethod) || !is.character(superiorityBestMethod) ||
      length(superiorityBestMethod) > 1 ||
      !superiorityBestMethod %in% c("min", "max")) {
    stop("superiorityBestMethod should be either min or max.\n")
  }
  if (!is.na(sorted) && (!is.character(sorted) || length(sorted) > 1 ||
      !sorted %in% c("ascending", "descending"))) {
    stop("sorted should be one of ascending, descending or NA")
  }
  if (missing(sorted)) {
    sorted <- "descending"
  }
  if (any(is.na(TD[[trait]]))) {
    y0 <- tapply(TD[[trait]], TD[, c("genotype","env")], mean)
    yIndex <- tapply(X = 1:nrow(TD), INDEX = TD[, c("genotype","env")], FUN = identity)
    na_yes_no <- is.na(y0)
    ## imputation
    y1 <- multMissing(y0, maxcycle = 10, na.strings = NA)
    replaceVal <- y1[na_yes_no]
    yIndexReplace <- yIndex[na_yes_no]
    if (is.list(yIndexReplace)) {
      nyr <- length(yIndexReplace)
      for (ii in 1:nyr) {
        for (jj in 1:length(TD[yIndexReplace[[ii]], trait])) {
          if (is.na(TD[yIndexReplace[[ii]][jj], trait]))
            TD[yIndexReplace[[ii]][jj], trait] <- replaceVal[ii]
        }
      }
    } else {
      TD[yIndexReplace, trait] <- replaceVal
    }
  }
  lab <- levels(TD$genotype)
  nGeno <- length(lab)
  envs <- levels(TD$env)
  nEnv <- length(envs)
  ## The centered trait mean for eniroment j
  Ej <- sapply(X = envs, FUN = function(x) {
    mean(TD[which(TD$env == x), trait], na.rm = TRUE)
  })
  ## Maximum or minimum trait mean among all genotype in jth enviroment
  if (superiorityBestMethod == "max") {
    Mj <- sapply(X = envs, FUN = function(x) {
      max(TD[which(TD$env == x), trait], na.rm = TRUE)
    })
  } else {
    Mj <- sapply(X = envs, FUN = function(x) {
      min(TD[which(TD$env == x), trait], na.rm = TRUE)
    })
  }
  ## The genotype trait mean across environments
  Ei <- sapply(X = lab, FUN = function(x) {
    mean(TD[which(TD$genotype == x), trait], na.rm = TRUE)
  })
  ## The grand mean
  E <- mean(TD[, trait], na.rm = TRUE)
  W <- S <- LB <- rep(NA, nGeno)
  for (i in 1:nGeno) {
    ## Observed genotype field response in the enviroment j
    ## (averaged across experiment replicates)
    Rij <- sapply(X = envs, FUN = function(x) {
      mean(TD[which(TD$genotype == lab[i] & TD$env == x), trait], na.rm = TRUE)
    })
    pos <- (1:nEnv)[!is.na(Rij)]
    ## Static measure (Shukla's (1972a) stability variance)
    if ("static" %in% method) {
      if (length(pos) == 0) {
        S[i] <- NA
      } else {
        S[i] <- sum(sapply(X = pos, FUN = function(j) {
          (Rij[j] - Ei[i]) ^ 2
        }) / (nEnv - 1))
      }
    }
    ## Superiority measure (LIN&BINNS 1988)
    if ("superiority" %in% method) {
      if (length(pos) == 0) {
        LB[i] <- NA
      } else{
        LB[i] <- sum(sapply(X = pos, FUN = function(j) {
          (Rij[j] - Mj[j]) ^ 2
        }) / (2 * nEnv))
      }
    }
    ## Wricke's (1962) ecovalence
    if ("wricke" %in% method) {
      if (length(pos) == 0) {
        W[i] <- NA
      } else {
        W[i] <- sum(sapply(X = pos, FUN = function(j) {
          (Rij[j] - Ei[i] - Ej[j] + E) ^ 2
        }))
      }
    }
  }
  ## Calculate trait mean per genotype
  means <- sapply(X = lab, FUN = function(x) {
    mean(TD[TD$genotype == x, trait], na.rm = TRUE)
  })
  if (plot) {
    ## Prepare panels
    mFlag <- c("static", "superiority", "wricke") %in% method
    if (sum(mFlag) != 1) {
      if (sum(mFlag) == 3) {
        oldPar <- par(mfrow = c(2, 2), mar = c(4, 4, 4, 2), oma = c(.5, .5, 2, .3))
      } else {
        oldPar <- par(mfrow = c(1, 2), mar = c(4, 4, 4, 2), oma = c(.5, .5, 2, .3))
      }
    } else{
      oldPar <- par(mfrow = c(1, 1))
    }
    on.exit(par(oldPar))
    if ("superiority" %in% method) {
      plot(x = means, y = LB, xlab = "Mean", ylab = "Cultivar superiority")
    }
    if ("static" %in% method) {
      plot(x = means, y = S, xlab = "Mean", ylab = "Static stability")
    }
    if ("wricke" %in% method) {
      plot(x = means, y = W, xlab = "Mean", ylab = "Wricke's ecovalence")
    }
    if (sum(mFlag) == 1) {
      title(paste('Stability coefficients for', trait))
    } else{
      mtext(paste('Stability coefficients for', trait), side = 3, outer = TRUE,
            cex = 1.3, font = 2)
    }
  }
  if (is.na(sorted)) {
    if ("static" %in% method) {
      SOut <- data.frame(lab, S, means, row.names = 1:nGeno)
      names(SOut) <- c("genotype","static", "mean")
    }
    if ("superiority" %in% method) {
      LBOut <- data.frame(lab, LB, means, row.names = 1:nGeno)
      names(LBOut) <- c("genotype", "superiority", "mean")
    }
    if ("wricke" %in% method) {
      WOut <- data.frame(lab, W, means, row.names = 1:nGeno)
      names(WOut) <- c("genotype", "wricke", "mean")
    }
  } else {
    if ("static" %in% method) {
      orderS <- order(S, decreasing = (sorted == "descending"))
      SOut <- data.frame(lab, S, means, row.names = 1:nGeno)[orderS, ]
      names(SOut) <- c("genotype", "static", "mean")
    }
    if ("superiority" %in% method) {
      orderLB <- order(LB, decreasing = (sorted == "descending"))
      LBOut <- data.frame(lab, LB, means, row.names = 1:nGeno)[orderLB, ]
      names(LBOut) <- c("genotype", "superiority", "mean")
    }
    if ("wricke" %in% method) {
      orderW <- order(W, decreasing = (sorted == "descending"))
      WOut <- data.frame(lab,W,means, row.names = 1:nGeno)[orderW, ]
      names(WOut) <- c("genotype", "wricke", "mean")
    }
  }
  res <- vector(mode = "list")
  if ("static" %in% method) res$static <- SOut
  if ("superiority" %in% method) res$superiority <- LBOut
  if ("wricke" %in% method) res$wricke <- WOut
  return(res)
}
