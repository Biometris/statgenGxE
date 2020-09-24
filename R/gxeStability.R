#' Calculate stability coefficients for genotype-by-environment data
#'
#' This function calculates different measures of stability, the
#' cultivar-superiority measure of Lin & Binns (1988), Shukla's (1972) stability
#' variance and Wricke's (1962) ecovalence.
#'
#' @inheritParams gxeAmmi
#'
#' @param method A character vector specifying the measures of stability to be
#' calculated. Options are "superiority" (cultivar-superiority measure),
#' "static" (Shukla's stability variance) or "wricke" (wricke's ecovalence).
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
#' @references Lin, C. S. and Binns, M. R. 1988. A superiority measure of
#' cultivar performance for cultivar x location data. Can. J. Plant Sci. 68:
#' 193-198
#' @references Shukla, G.K. 1972. Some statistical aspects of partitioning
#' genotype-environmental components of variability. Heredity 29:237-245
#' @references Wricke, G. Uber eine method zur erfassung der okologischen
#' streubreit in feldversuchen. Zeitschrift für Pflanzenzucht,
#' v. 47, p. 92-96, 1962
#'
#' @examples
#' ## Compute three stability measures for TDMaize.
#' geStab <- gxeStability(TD = TDMaize, trait = "yld")
#'
#' ## Summarize results.
#' summary(geStab)
#'
#' ## Create plot of the computed stability measures against the means.
#' plot(geStab)
#'
#' \donttest{
#' ## Create a .pdf report summarizing the stability measures.
#' report(geStab, outfile = tempfile(fileext = ".pdf"))
#' }
#'
#' ## Compute Wricke's ecovalance for TDMaize with minimal values for yield as
#' ## the best values. Sort results in ascending order.
#' geStab2 <- gxeStability(TD = TDMaize, trait = "yld", method = "wricke",
#'                        bestMethod = "min", sorted = "ascending")
#' summary(geStab2)
#'
#' @family stability
#'
#' @importFrom methods getFunction
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
  trials <- chkTrials(trials, TD)
  TDTot <- do.call(rbind, args = TD[trials])
  TDTot <- droplevels(TDTot)
  chkCol(trait, TDTot)
  chkCol("genotype", TDTot)
  method <- match.arg(method, several.ok = TRUE)
  bestMethod <- match.arg(bestMethod)
  sorted <- match.arg(sorted)
  trCol <- "trial"
  chkCol(trCol, TDTot)
  ## Remove genotypes that contain only NAs
  allNA <- by(TDTot, TDTot[["genotype"]], FUN = function(x) {
    all(is.na(x[trait]))
  })
  TDTot <- TDTot[!TDTot[["genotype"]] %in% names(allNA[allNA]), ]
  if (any(is.na(TDTot[[trait]]))) {
    y0 <- tapply(TDTot[[trait]], TDTot[, c("genotype", trCol)], mean)
    ## Actual imputation.
    y1 <- multMissing(y0, maxIter = 50)
    ## Insert imputed values back into original data.
    for (missVal in which(is.na(TDTot[[trait]]))) {
      TDTot[missVal, trait] <- y1[TDTot[missVal, "genotype"],
                                  TDTot[missVal, "trial"]]
    }
  }
  lab <- levels(TDTot[["genotype"]])
  nGeno <- length(lab)
  nTr <- length(trials)
  ## Compute the centered trait mean per environment.
  Ej <- tapply(X = TDTot[[trait]], INDEX = TDTot[[trCol]], FUN = mean,
               na.rm = TRUE)
  ## Compute the max or min trait mean among all genotypes per trial.
  Mj <- tapply(X = TDTot[[trait]], INDEX = TDTot[[trCol]],
               FUN = getFunction(bestMethod), na.rm = TRUE)
  ## Compute the genotype trait mean per trial.
  Ei <- as.numeric(tapply(X = TDTot[[trait]], INDEX = TDTot[["genotype"]],
                          FUN = mean, na.rm = TRUE))
  ## Compute the grand mean.
  E <- mean(TDTot[, trait], na.rm = TRUE)
  ## Create empty vectors for storing output.
  LB <- S <- W <- rep(NA, nGeno)
  for (i in seq_along(lab)) {
    ## Observed genotype field response in the trial j.
    ## (averaged across experiment replicates).
    Rij <- sapply(X = trials, FUN = function(x) {
      mean(TDTot[TDTot[["genotype"]] == lab[i] & TDTot[[trCol]] == x, trait],
           na.rm = TRUE)
    })
    pos <- (1:nTr)[!is.na(Rij)]
    if (length(pos) > 0) {
      ## Superiority measure (Lin & Binns 1988).
      if ("superiority" %in% method) {
        LB[i] <- sum(sapply(X = pos, FUN = function(j) {
          (Rij[j] - Mj[j]) ^ 2
        }) / (2 * nTr))
      }
      ## Static measure (Shukla's (1972a) stability variance).
      if ("static" %in% method) {
        S[i] <- sum(sapply(X = pos, FUN = function(j) {
          (Rij[j] - Ei[i]) ^ 2
        }) / (nTr - 1))
      }
      ## Wricke's (1962) ecovalence.
      if ("wricke" %in% method) {
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
    LBOut <- setNames(data.frame(factor(lab, levels = lab), Ei, LB,
                                 row.names = 1:nGeno)[orderLB, ],
                      c("genotype", "mean", "superiority"))
  }
  if ("static" %in% method) {
    if (sorted == "none") {
      orderS <- 1:nGeno
    } else {
      orderS <- order(S, decreasing = (sorted == "descending"))
    }
    SOut <- setNames(data.frame(factor(lab, levels = lab), Ei, S,
                                row.names = 1:nGeno)[orderS, ],
                     c("genotype", "mean", "static"))
  }
  if ("wricke" %in% method) {
    if (sorted == "none") {
      orderW <- 1:nGeno
    } else {
      orderW <- order(W, decreasing = (sorted == "descending"))
    }
    WOut <- setNames(data.frame(factor(lab, levels = lab), Ei, W,
                                row.names = 1:nGeno)[orderW, ],
                     c("genotype", "mean", "wricke"))
  }
  res <- createStability(superiority = if ("superiority" %in% method) LBOut,
                         static = if ("static" %in% method) SOut,
                         wricke = if ("wricke" %in% method) WOut,
                         TD = TD, trait = trait)
  return(res)
}
