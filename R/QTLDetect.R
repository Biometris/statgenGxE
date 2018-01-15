#' Function for QTL detection
#'
#' @param cross An object of class cross created by the qtl package
#' @param trait A character string indicating the trait to be analysed.
#' @param type A character sting indicating the type of QTL detection to be performed.
#' Either "MR" (Marker Response), "SIM" (Simple Interval Mapping) or "CIM"
#' (Composite Interval Mapping)
#' @param thr A numerical value indicating a lower threshold for the lodscore
#' of the peaks.
#' @param window A numerical value indicating the window (in cM) used when
#' selecting peaks.
#' @param ... Other parameters to be passed on to underlying functions.
#'
#' @examples
#' ## Read the data
#' F2 <- qtl::read.cross(format="csv",
#'                       file = system.file("extdata", "F2_maize_practical3_ex2.csv",
#'                       package = "RAP"),
#'                       genotypes = c("AA", "AB", "BB"),
#'                       alleles = c("A", "B"), estimate.map = FALSE)
#' QTLDet <- QTLDetect(cross = F2, trait = "trait", type = "SIM")
#' report(QTLDet, outfile = "./testReports/reportQTLDectection.pdf")
#'
#' @export
QTLDetect <- function(cross,
                      trait,
                      type = c("MR", "SIM", "CIM"),
                      thr = 3,
                      window = 15,
                      ...) {
  type <- match.arg(type)
  if (type == "MR") {
    ## Perform a marker-based QTL detection.
    scores <- qtl::scanone(cross, pheno.col = trait, method = "mr", ...)
  } else if (type == "SIM") {
    ## Calculate genotype probabilities.
    cross <- qtl::calc.genoprob(cross, step = 5, error.prob = 0)
    ## Perform a QTL search by Simple Interval Mapping (SIM)
    ## (Haley-Knott regression)
    scores <- qtl::scanone(cross, pheno.col = trait, method = "hk", ...)
  } else if (type == "CIM") {
    ## Calculate genotype probabilities.
    cross <- qtl::calc.genoprob(cross, step = 5, error.prob = 0)
    ## Perform a QTL search by Simple Interval Mapping (CIM)
    scores <- qtl::cim(cross, pheno.col = trait, n.marcovar = 5, window = 50,
                       method = "hk", map.function = "haldane",
                       ...)
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
      qtlCand <- qtlCand[!(qtlCand$chr == peaks[i, "chr"] &
                             abs(qtlCand$pos - peaks[i, "pos"]) < window), ]
    }
  }
  peaksTot <- peaksTot[order(peaksTot$chr, peaksTot$pos), ]
  attr(scores, "marker.covar") <- rownames(peaksTot)
  attr(scores, "marker.covar.pos") <- peaksTot[, c("chr", "pos")]
  QTLDet <- createQTLDet(scores = scores,
                         peaks = peaksTot,
                         type = type,
                         cross = cross,
                         trait = trait)
  return(QTLDet)
}





