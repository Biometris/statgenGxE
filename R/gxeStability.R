#' Calculate stability coefficients for genotype-by-environment data
#'
#' This function calculates different measures of stability, the
#' cultivar-superiority measure of Lin & Binns (1988), Shukla's (1972) stability
#' variance and Wricke's (1962) ecovalence.
#'
#' @inheritParams gxeAmmi
#'
#' @param method A character vector specifying the measures of stability to be
#' calculated. Options are superiority (cultivar-superiority measure), static
#' (Shukla's stability variance) or wricke (wricke's ecovalence).
#' @param bestMethod A character string specifying the criterion to define
#' the best genotype. Either \code{"max"} or \code{"min"}.
#' @param sorted A character string specifying the sorting order of the results.
#'
#' @return An object of class \code{\link{stability}}, a list containing:
#' \item{superiority}{A data.frame containing values for the
#' cultivar-superiority measure of Lin and Binns.}
#' \item{static}{A data.frame containing values for Shukla's stability
#' variance.}
#' \item{wricke}{A data.frame containing values for Wricke's ecovalence.}
#' \item{trait}{A character string indicating the trait that has been analyzed.}
#'
#' @references LiN, C. S. and Binns, M. R. 1988. A superiority measure of
#' cultivar performance for cultivar x location data. Can. J. Plant Sci. 68:
#' 193-198
#' @references Shukla, G.K. 1972. Some statistical aspects of partitioning
#' genotype-environmental components of variability. Heredity 29:237-245
#' @references Wricke, G. Uber eine method zur erfassung der okologischen
#' streubreit in feldversuchen. Zeitschrift für Pflanzenzucht,
#' v. 47, p. 92-96, 1962
#'
#' @seealso \code{\link{stability}}, \code{\link{plot.stability}},
#' \code{\link{report.stability}}
#'
#' @examples
#' ## Compute three stability measures for TDMaize.
#' geStab <- gxeStability(TD = TDMaize, trait = "yld")
#' ## Summarize results.
#' summary(geStab)
#' ## Create plot of the computed stability measures against the means.
#' plot(geStab)
#' \dontrun{
#' ## Create a .pdf report summarizing the stability measures.
#' report(geStab, outfile = "./testReports/reportStability.pdf")
#' }
#'
#' ## Compute Wricke's ecovalance for TDMaize with minimal values for yield as
#' ## the best values. Sort results in ascending order.
#' geStab2 <- gxeStability(TD = TDMaize, trait = "yld", method = "wricke",
#'                        bestMethod = "min", sorted = "ascending")
#' summary(geStab2)
#'
#' @export
gxeStability <- function(TD,
                         trials = names(TD),
                         trait,
                         method = c("superiority", "static", "wricke"),
                         bestMethod = c("max", "min"),
                         sorted = c("descending", "ascending", "none")) {
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (!is.character(trials) || !all(trials %in% names(TD))) {
    stop("All trials should be in TD.")
  }
  TDTot <- Reduce(f = rbind, x = TD[trials])
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TDTot)) {
    stop("trait has to be a column in TD.\n")
  }
  if (!"trial" %in% colnames(TDTot)) {
    stop("TD should contain a column trial to be able to calculate stabilities.\n")
  }
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TDTot)) {
    stop("trait has to be a column in TD.\n")
  }
  method <- match.arg(method, several.ok = TRUE)
  bestMethod <- match.arg(bestMethod)
  sorted <- match.arg(sorted)
  ## Remove genotypes that contain only NAs
  allNA <- by(TDTot, TDTot$genotype, FUN = function(x) {
    all(is.na(x[trait]))
  })
  TDTot <- TDTot[!TDTot$genotype %in% names(allNA[allNA]), ]
  if (any(is.na(TDTot[[trait]]))) {
    y0 <- tapply(TDTot[[trait]], TDTot[, c("genotype","trial")], mean)
    yIndex <- tapply(X = 1:nrow(TDTot), INDEX = TDTot[, c("genotype","trial")],
                     FUN = identity)
    ## imputate missing values.
    y1 <- multMissing(y0, maxIter = 10)
    replaceVal <- y1[is.na(y0)]
    yIndexReplace <- yIndex[is.na(y0)]
    if (is.list(yIndexReplace)) {
      for (i in 1:length(yIndexReplace)) {
        for (j in 1:length(TDTot[yIndexReplace[[i]], trait])) {
          if (is.na(TDTot[yIndexReplace[[i]][j], trait]))
            TDTot[yIndexReplace[[i]][j], trait] <- replaceVal[i]
        }
      }
    } else {
      TDTot[yIndexReplace, trait] <- replaceVal
    }
  }
  lab <- levels(TDTot$genotype)
  nGeno <- length(lab)
  trials <- levels(TDTot$trial)
  nTr <- length(trials)
  ## Compute the centered trait mean per eniroment.
  Ej <- tapply(TDTot[[trait]], TDTot[, "trial"], mean)
  ## Compute the max or min trait mean among all genotypes per trial.
  if (bestMethod == "max") {
    Mj <- tapply(TDTot[[trait]], TDTot[, "trial"], max, na.rm = TRUE)
  } else {
    Mj <- tapply(TDTot[[trait]], TDTot[, "trial"], min, na.rm = TRUE)
  }
  ## Compute the genotype trait mean per trial.
  Ei <- as.vector(tapply(TDTot[[trait]], TDTot[, "genotype"], mean,
                         na.rm = TRUE))
  ## Compute the grand mean.
  E <- mean(TDTot[, trait], na.rm = TRUE)
  ## Create empty vectors for storing output
  W <- S <- LB <- rep(NA, nGeno)
  for (i in 1:nGeno) {
    ## Observed genotype field response in the trial j
    ## (averaged across experiment replicates)
    Rij <- sapply(X = trials, FUN = function(x) {
      mean(TDTot[which(TDTot$genotype == lab[i] & TDTot$trial == x), trait],
           na.rm = TRUE)
    })
    pos <- (1:nTr)[!is.na(Rij)]
    ## Superiority measure (Lin & Binns 1988).
    if ("superiority" %in% method) {
      if (length(pos) > 0) {
        LB[i] <- sum(sapply(X = pos, FUN = function(j) {
          (Rij[j] - Mj[j]) ^ 2
        }) / (2 * nTr))
      }
    }
    ## Static measure (Shukla's (1972a) stability variance).
    if ("static" %in% method) {
      if (length(pos) > 0) {
        S[i] <- sum(sapply(X = pos, FUN = function(j) {
          (Rij[j] - Ei[i]) ^ 2
        }) / (nTr - 1))
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
