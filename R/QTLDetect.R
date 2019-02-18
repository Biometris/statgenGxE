#' QTL detection
#'
#' This function is essentially a wrapper for \code{\link[qtl]{scanone}} and
#' \code{\link[qtl]{cim}} in the qtl package. Depending on \code{type}, one
#' of these functions is used for (Quantitative Trait Locus) QTL detection.
#' After this is done, from the set of candidate QTLs that are returned, proper
#' peaks are selected by an iterative process using the \code{thr}(eshold) and
#' \code{window} provided. All resulting peaks will have a LOD score above
#' \code{thr} and the distance between pairs of peaks will always be at least
#' the value given as \code{window}.
#'
#' The option \code{thrType} specifies the method for calculating the number of
#' tests to be used as the denominator in a modified Bonferroni correction. The
#' enumeration is given by \code{thrAlpha}.
#' \code{liji} uses the effective number of independent tests, as described by
#' Li & Ji (2005). \code{bonferroni} assumes one independent test at every
#' fixed distance on the genome, defined by the \code{thrDist} (default 4
#' centiMorgans). The threshold is expressed as the P-value on a -log10 scale.
#' By setting {thrType = "fixed"} all the above is ignored and a fixed value is
#' used as threshold instead. This value can be supplied in \code{thrFixed}
#'
#' @note Composite Interval Mapping is not implemented for 4 way crosses.
#'
#' @param cross An object of class cross created by the qtl package.
#' @param trait A character string indicating the trait to be analyzed.
#' @param type A character string indicating the type of QTL detection to be
#' performed. Either "MR" (Marker Response), "SIM" (Simple Interval Mapping)
#' or "CIM" (Composite Interval Mapping)
#' @param step A numerical value indicating the maximum distance (in cM)
#' between positions at which the genotype probabilities are calculated. If
#' \code{step = 0} probabilities are only calculated at the marker locations.
#' @param thrType A character string indicating the algorithm to calculate the
#' lower threshold for the lodscore of the peaks. See details.
#' @param thrAlpha A numerical value between 0 and 1 used as enumerator in
#' calculating the modified Bonferroni correction.
#' @param thrDist A positive numerical value indicating the assumed fixed
#' distance on the genome for independent tests (in cM). Only used
#' when \code{thrType = "bonferroni"}
#' @param thrFixed A numerical value indicating a lower threshold for the
#' lodscore of the peaks. Only used when \code{thrType = "fixed"}
#' @param window A numerical value indicating the window (in cM) used when
#' selecting peaks.
#' @param ... Further parameters to be passed on to underlying functions used
#' for qtl detection, \code{\link[qtl]{scanone}} when \code{type} is "MR"
#' or "SIM" and \code{\link[qtl]{cim}} when \code{type} is "CIM".
#'
#' @return An object of class \code{\link{QTLDet}}, a list containing:
#' \item{scores}{A data.frame containing the lod scores.}
#' \item{peaks}{A data.frame containing the peaks found.}
#' \item{type}{A character string indicating the type of QTL detection
#' performed.}
#' \item{cross}{The object of class cross used for QTL detection.}
#' \item{trait}{A character string indicating the trait for which the detection
#' was done.}
#' \item{info}{A list containing information on the settings used for
#' QTL detection, i.e. step, threshold and window.}
#'
#' @seealso \code{\link[qtl]{scanone}}, \code{\link[qtl]{cim}}
#'
#' @references Broman et al. (2003) R/qtl: QTL mapping in experimental crosses.
#' Bioinformatics 19:889-890
#' @references Cheverud, J.M. (2001). A simple correction for multiple
#' comparisons in interval mapping genome scans. Heredity, 87, 52-58.
#' @references Li, J, & Ji, L. (2005). Adjusting multiple testing in multilocus
#' analyses using the eigenvalues of a correlation matrix. Heredity, 95,
#' 221-227.
#'
#' @examples
#' ## Read the data.
#' F2 <- qtl::read.cross(format="csv",
#'                       file = system.file("extdata", "F2maize_geno.csv",
#'                                         package = "RAP"),
#'                       genotypes = c("AA", "AB", "BB"),
#'                       alleles = c("A", "B"), estimate.map = FALSE)
#' ## Perform a composite interval mapping for detecting QTLs.
#' QTLDet <- QTLDetect(cross = F2, trait = "trait", type = "CIM")
#' ## Summarize results.
#' summary(QTLDet)
#' ## Create a manhattan plot of the results.
#' plot(QTLDet)
#' \dontrun{
#' ## Create a pdf report summarizing the results.
#' report(QTLDet, outfile = "./testReports/reportQTLDectection.pdf")
#' }
#'
#' ## Perform a simple interval mapping for detecting QTLs.
#' ## Choose custom step, threshold and window sizes.
#' QTLDet2 <- QTLDetect(cross = F2, trait = "trait", type = "SIM", step = 15,
#'                     thrType = "fixed", thrFixed = 2.5, window = 50)
#' summary(QTLDet2)
#'
#' @export
QTLDetect <- function(cross,
                      trait,
                      type = c("MR", "SIM", "CIM"),
                      step = 5,
                      thrType = c("liji", "bonferroni", "fixed"),
                      thrAlpha = 0.05,
                      thrDist = 4,
                      thrFixed = 3,
                      window = 30,
                      ...) {
  ## Checks.
  if (missing(cross) || !inherits(cross, "cross")) {
    stop("cross should be an object of class cross\n")
  }
  if (missing(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(cross$pheno)) {
    stop("trait should be a single character string for a trait in cross")
  }
  type <- match.arg(type)
  if (type == "CIM" && "4way" %in% class(cross)) {
    stop("CIM has not been implemented for 4-way crosses.\n")
  }
  if (!is.numeric(step) || length(step) > 1 || step < 0) {
    stop("step should be a single positive numerical value.\n")
  }
  thrType <- match.arg(thrType)
  if (thrType != "fixed") {
    if (!is.numeric(thrAlpha) || length(thrAlpha) > 1 || thrAlpha < 0 ||
        thrAlpha > 1) {
      stop("thrFixed should be a single positive numerical value betwee 0 and 1.\n")
    }
  } else {
    if (!is.numeric(thrFixed) || length(thrFixed) > 1 || thrFixed < 0) {
      stop("thrFixed should be a single positive numerical value.\n")
    }
  }
  if (thrType == "bonferroni") {
    if (!is.numeric(thrDist) || length(thrDist) > 1 || thrDist < 0) {
      stop("thrDist should be a single positive numerical value.\n")
    }
  }
  if (!is.numeric(window) || length(window) > 1 || window < 0) {
    stop("window should be a single positive numerical value.\n")
  }
  ## Calculate genotype probabilities.
  ## Not strictly needed for "MR" but when perferming multiQTLFit after
  ## QTLDetect is has to be done.
  cross <- qtl::calc.genoprob(cross, step = step, error.prob = 0)
  if (type == "MR") {
    ## Perform a marker-based QTL detection.
    scores <- qtl::scanone(cross, pheno.col = trait, method = "mr", ...)
  } else if (type == "SIM") {
    ## Perform a QTL search by Simple Interval Mapping (SIM)
    ## (Haley-Knott regression)
    scores <- qtl::scanone(cross, pheno.col = trait, method = "hk", ...)
  } else if (type == "CIM") {
    ## Perform a QTL search by Composite Interval Mapping (CIM)
    scores <- qtl::cim(cross, pheno.col = trait, n.marcovar = 5,
                       window = window, method = "hk", map.function = "haldane",
                       ...)
  }
  if (thrType == "liji") {
    ## Compute correction values following the algorithm by Li and Ji.
    corVal <- sapply(X = cross$geno, FUN = function(chr) {
      ## Compute correaltions between markers on chromosome.
      mrkCor <- cor(chr$data, use = "complete.obs")
      ## Compute eigenvalues of the correlations and restrict these to 2.
      mrkEig <- eigen(mrkCor, only.values = TRUE)$values
      mrkEig[mrkEig >= 1] <- mrkEig[mrkEig > 1] - floor(mrkEig[mrkEig > 1]) + 1
      sum(mrkEig)
    })
    ## Compute modified Bonferroni threshold.
    thr <- -log10(thrAlpha / sum(corVal))
  } else if (thrType == "bonferroni") {
    ## Compute modified bonferroni threshold using distance.
    thr <- -log10(thrAlpha / ceiling(sum(qtl::chrlen(cross)) / thrDist))
  } else if (thrType == "fixed") {
    thr <- thrFixed
  }
  ## Select set of candidate QTLs to be fitted in a final multiQTL model.
  qtlCand <- scores[scores$lod > thr, ]
  ## Iteratively remove peaks and surrounding qtls from candidates.
  ## Add those to peaks data.
  peaksTot <- qtlCand[FALSE, ]
  while (nrow(qtlCand) > 0) {
    peaks <- summary(qtlCand)
    peaksTot <- rbind(peaksTot, peaks)
    for (i in 1:nrow(peaks)) {
      ## Remove candidates on the same chromosome that are within window.
      qtlCand <- qtlCand[!(qtlCand$chr == peaks[i, "chr"] &
                             abs(qtlCand$pos - peaks[i, "pos"]) < window), ]
    }
  }
  peaksTot <- peaksTot[order(peaksTot$chr, peaksTot$pos), ]
  if (nrow(peaksTot) > 0) {
    ## Add alternative name to peaks.
    ## stringsAsFactors is needed to assure altName is not converted to factor.
    peaksTot <- cbind(altName = paste0("Q", peaksTot$chr, "@", peaksTot$pos),
                      peaksTot, stringsAsFactors = FALSE)
  }
  attr(scores, "marker.covar") <- rownames(peaksTot)
  attr(scores, "marker.covar.pos") <- peaksTot[, c("chr", "pos")]
  info = list(step = step, thrType = thrType, thrAlpha = thrAlpha, thr = thr,
              window = window)
  QTLDet <- createQTLDet(scores = scores, peaks = peaksTot, type = type,
                         cross = cross, trait = trait, info = info)
  return(QTLDet)
}
