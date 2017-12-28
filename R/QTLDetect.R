#' Function for QTL detection
#'
#' @param cross An object of class cross created by the qtl package
#' @param type A sting indicating the type of QTL detection to be performed.
#' Either "MR" (Marker Response), "SIM" (Simple Interval Mapping) or "CIM"
#' (Composite Interval Mapping)
#'
#'
#' @export
QTLDetect <- function(cross,
                      type = "MR") {

  if (type == "MR") {
    ## Perform a marker-based QTL detection.
    scores <- qtl::scanone(cross, method = "mr")
  } else if (type == "SIM") {
    ## Calculate genotype probabilities.
    cross <- qtl::calc.genoprob(cross, step = 5, error.prob = 0)
    ## Perform a QTL search by Simple Interval Mapping (SIM)
    ## (Haley-Knott regression)
    scores <- qtl::scanone(cross, method = "hk")
  } else if (type == "CIM") {
    ## Calculate genotype probabilities.
    cross <- qtl::calc.genoprob(cross, step = 5, error.prob = 0)
    ## Perform a QTL search by Simple Interval Mapping (CIM)
    scores <- qtl::cim(cross, n.marcovar = 5, window = 50,
                       method = "hk", map.function = "haldane")
  }
  QTLDet <- createQTLDet(scores = scores,
                         type = type,
                         cross = cross)
  return(QTLDet)
}

