#' Calculates stability coefficients for genotype-by-environment data
#'
#' This function calculate difference measures of stability, such as the
#' cultivar-superiority measure of Lin & Binns (1988), Shukla's (1972) stability
#' variance and Wricke's (1962) ecovalence.
#'
#' @inheritParams gxeAmmi
#'
#' @param method A character vector specifying (a) measure(s) of stability to be
#' calculated. Options are superiority (cultivar-superiority measure), static
#' (Shukla's stability variance) or wricke (wricke's ecovalence). By default
#' all three measures are calculated.
#' @param bestMethod A character string specifying the criterion to define
#' the best genotype ("max","min").
#' @param sorted A character string specifying the sorting order of the results.
#'
#' @return A list of one to three data.frames containing the stability measures.
#'
#' @references LiN, C. S. and Binns, M. R. 1988. A superiority measure of cultivar
#' performance for cultivar x location data. Can. J. Plant Sci. 68: 193-198
#' @references Shukla, G.K. 1972. Some statistical aspects of partitioning
#' genotype-environmental components of variability. Heredity 29:237-245
#' @references Wricke, G. Uber eine method zur erfassung der okologischen
#' streubreit in feldversuchen. Zeitschrift für Pflanzenzucht,
#' v. 47, p. 92-96, 1962
#'
#' @examples
#' geStab <- gxeStability(TD = TDMaize, trait = "yld", sorted = "descending")
#' report(geStab, outfile = "./testReports/reportStability.pdf")
#'
#' @export
gxeStability <- function(TD,
                         trait,
                         method = c("superiority", "static", "wricke"),
                         bestMethod = c("max", "min"),
                         sorted = c("descending", "ascending", "none")) {
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
  method <- match.arg(method, several.ok = TRUE)
  bestMethod <- match.arg(bestMethod)
  sorted <- match.arg(sorted)
  if (any(is.na(TD[[trait]]))) {
    y0 <- tapply(TD[[trait]], TD[, c("genotype","env")], mean)
    yIndex <- tapply(X = 1:nrow(TD), INDEX = TD[, c("genotype","env")],
                     FUN = identity)
    ## imputate missing values.
    y1 <- multMissing(y0, maxIter = 10)
    replaceVal <- y1[is.na(y0)]
    yIndexReplace <- yIndex[is.na(y0)]
    if (is.list(yIndexReplace)) {
      for (i in 1:length(yIndexReplace)) {
        for (j in 1:length(TD[yIndexReplace[[i]], trait])) {
          if (is.na(TD[yIndexReplace[[i]][j], trait]))
            TD[yIndexReplace[[i]][j], trait] <- replaceVal[i]
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
  ## Compute the centered trait mean per eniroment.
  Ej <- tapply(TD[[trait]], TD[, "env"], mean)
  ## Compute the maximum or minimum trait mean among all genotypes per enviroment.
  if (bestMethod == "max") {
    Mj <- tapply(TD[[trait]], TD[, "env"], max, na.rm = TRUE)
  } else {
    Mj <- tapply(TD[[trait]], TD[, "env"], min, na.rm = TRUE)
  }
  ## Compute the genotype trait mean per environment.
  Ei <- as.vector(tapply(TD[[trait]], TD[, "genotype"], mean, na.rm = TRUE))
  ## Compute the grand mean.
  E <- mean(TD[, trait], na.rm = TRUE)
  ## Create empty vectors for storing output
  W <- S <- LB <- rep(NA, nGeno)
  for (i in 1:nGeno) {
    ## Observed genotype field response in the enviroment j
    ## (averaged across experiment replicates)
    Rij <- sapply(X = envs, FUN = function(x) {
      mean(TD[which(TD$genotype == lab[i] & TD$env == x), trait], na.rm = TRUE)
    })
    pos <- (1:nEnv)[!is.na(Rij)]
    ## Superiority measure (Lin & Binns 1988).
    if ("superiority" %in% method) {
      if (length(pos) > 0) {
        LB[i] <- sum(sapply(X = pos, FUN = function(j) {
          (Rij[j] - Mj[j]) ^ 2
        }) / (2 * nEnv))
      }
    }
    ## Static measure (Shukla's (1972a) stability variance).
    if ("static" %in% method) {
      if (length(pos) > 0) {
        S[i] <- sum(sapply(X = pos, FUN = function(j) {
          (Rij[j] - Ei[i]) ^ 2
        }) / (nEnv - 1))
      }
    }
    ## Wricke's (1962) ecovalence.
    if ("wricke" %in% method) {
      if (length(pos) > 0) {
        W[i] <- sum(sapply(X = pos, FUN = function(j) {
          (Rij[j] - Ei[i] - Ej[j] + E) ^ 2
        }))
      }
    }
  }
  ## Create output data.frames.
  if ("superiority" %in% method) {
    if (sorted == "none") {
      orderLB <- 1:nGeno
    } else {
      orderLB <- order(LB, decreasing = (sorted == "descending"))
    }
    LBOut <- setNames(data.frame(lab, LB, Ei, row.names = 1:nGeno)[orderLB, ],
                      c("genotype", "superiority", "mean"))
  }
  if ("static" %in% method) {
    if (sorted == "none") {
      orderS <- 1:nGeno
    } else {
      orderS <- order(S, decreasing = (sorted == "descending"))
    }
    SOut <- setNames(data.frame(lab, S, Ei, row.names = 1:nGeno)[orderS, ],
                     c("genotype", "static", "mean"))
  }
  if ("wricke" %in% method) {
    if (sorted == "none") {
      orderW <- 1:nGeno
    } else {
      orderW <- order(W, decreasing = (sorted == "descending"))
    }
    WOut <- setNames(data.frame(lab, W, Ei, row.names = 1:nGeno)[orderW, ],
                     c("genotype", "wricke", "mean"))
  }
  res <- createstability(superiority = if ("superiority" %in% method) LBOut,
                         static = if ("static" %in% method) SOut,
                         wricke = if ("wricke" %in% method) WOut,
                         trait = trait)
  return(res)
}
