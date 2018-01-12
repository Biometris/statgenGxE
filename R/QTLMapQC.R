#' Function for Quality Control of a cross object
#'
#' Function for performing quality control of a cross object.
#'
#' @param cross An object of class cross created by the qtl package
#' @param missMrk A numerical value between 0 and 1 indicating the maximum
#' allowed fraction of missing values per marker. Markers with a fraction of
#' missing values above \code{missMrk} will be removed.
#' @param missInd A numerical value between 0 and 1 indicating the maximum
#' allowed fraction of missing values per individual. Individuals with a fraction of
#' missing values above \code{missMrk} will be removed.
#' @param removeDuplicates Should duplicate markers be removed?
#' @param segDistortion A numberical value between 0 and 1 used a threshold for
#' Mendelian segregation. Markers with a P-value below \code{segDistortion} will
#' be removed.
#' @param recombination A positive numerical value used a threshold for checking
#' recombination between pairs of markers.
#' @param reestimateMap Should the map be reestimated based on the observed markers?
#' @param crossover A numerical value between 0 and 1 indicating the maximum
#' allowed fraction of crossovers per individual. Individuals with a fraction of
#' crossovers above \code{crossover} will be removed.
#'
#' @return A cleaned version of the input cross object after markers and
#' individuals have been removed that are outside the respective thresholds.
#'
#' @seealso \code{\link[qtl]{nmissing}}, \code{\link[qtl]{drop.markers}},
#' \code{\link[qtl]{drop.dupmarkers}}, \code{\link[qtl]{geno.table}},
#' \code{\link[qtl]{nmissing}}, \code{\link[qtl]{checkAlleles}},
#' \code{\link[qtl]{replace.map}}, \code{\link[qtl]{countXO}}
#'
#' @examples
#' # Read the data
#' F2 <- qtl::read.cross(format="csv",
#'                       file = system.file("extdata", "F2_maize_practical3_ex2.csv",
#'                                          package = "RAP"),
#'                       genotypes = c("AA", "AB", "BB"),
#'                       alleles = c("A", "B"), estimate.map = FALSE)
#' # Run quality control
#' F2QC <- QTLMapQC(F2)
#'
#' @export
QTLMapQC <- function(cross,
                     missMrk = 0.05,
                     missInd = 0.05,
                     removeDuplicates = TRUE,
                     segDistortion = 0.001,
                     recombination = 3,
                     reestimateMap = FALSE,
                     crossover = 0.20) {
  ## Extract number of individuals
  nInd <- qtl::nind(cross)
  if (missMrk > 0) {
    ## Calculate pct missing per marker and remover markers with pct above missMrk.
    pctMiss <- qtl::nmissing(cross, what = "mar") / nInd
    dropMrk <- names(pctMiss[pctMiss > missMrk])
    cross <- qtl::drop.markers(cross, markers = dropMrk)
  }
  if (removeDuplicates) {
    ## Remove duplicate markers.
    cross <- qtl::drop.dupmarkers(cross, verbose = FALSE)
  }
  nMrk <- sum(qtl::nmar(cross))
  if (missInd > 0) {
    ## Calculate pct missing per individual and remover individuals with pct above missInd.
    pctMiss <- qtl::nmissing(cross, what = "ind") / nMrk
    dropInd <- which(pctMiss > missInd)
    if (length(dropInd) > 0) {
      cross <- cross[, -dropInd]
    }
  }
  ## Compute the segregation distortion per marker and remove those showing evidence
  ## of distortion.
  if (segDistortion > 0) {
    segDist <- qtl::geno.table(cross)
    dropSegDist <- rownames(segDist[segDist$P.value < segDistortion, ])
    cross <- qtl::drop.markers(cross, markers = dropSegDist)
  }
  ## Estimate recombination frequencies between all pairs of markers.
  crossRec <- qtl::est.rf(cross)
  if (recombination > 0) {
    ## Identify markers which might have been switched and remove those.
    dropRecom <- qtl::checkAlleles(crossRec, threshold = recombination,
                                   verbose = FALSE)
    crossRec <- qtl::drop.markers(crossRec, markers = dropRecom$marker)
  }
  if (reestimateMap) {
    ## Re-estimate map based on the observed markers
    newMap <- qtl::est.map(crossRec, error.prob = 1e-3)
    crossRec <- qtl::replace.map(crossRec, newMap)
  }
  if (crossover > 0) {
    ## Check pct of crossovers per individual and remove individuals with a
    ## pct above crossover
    pctCross <- qtl::countXO(crossRec) / sum(qtl::nmar(crossRec))
    dropCrossInd <- which(pctCross > crossover)
    if (length(dropCrossInd) > 0) {
      crossRec <- crossRec[, -dropCrossInd]
    }
  }
  return(crossRec)
}

#' Report method for class cross
#'
#' A pdf report will be created containing a summary of cross analysis.
#' Simultaneously the same report will be created as a tex file.
#'
#' @param x an object of class cross.
#' @param ... further arguments passed on from other functions - not used yet.
#' @param outfile a character string, the name and location of the output .pdf and .tex
#' file for the report. If \code{NULL} a report will be created in the current working
#' directory.
#'
#' @examples
#' ## Read the data
#' F2 <- qtl::read.cross(format="csv",
#'                       file = system.file("extdata", "F2_maize_practical3_ex2.csv", package = "RAP"),
#'                       genotypes = c("AA", "AB", "BB"),
#'                       alleles = c("A", "B"), estimate.map = FALSE)
#' report(F2, outfile = "./testReports/reportCross.pdf")
#'
#' @export
report.cross <- function(x,
                         ...,
                         outfile = NULL) {
  createReport(x = x, reportName = "crossReport.Rnw",
               outfile = outfile, ...)
}





