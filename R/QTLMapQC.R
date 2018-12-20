#' Quality Control and Cleaning of a cross object
#'
#' Function for performing quality control and cleaning of a cross object.
#' Quality control is done in several subsequent steps:
#' \enumerate{
#' \item{Markers with a fraction of missing values higher than \code{missMrk}
#' are removed.}
#' \item{Duplicate markers are removed.}
#' \item{Individuals with a fraction of missing values higher than \code{missInd}
#' are removed.}
#' \item{Markers that show evidence of segregation distortion (P-value below
#' \code{segDistortion}) are removed.}
#' \item{Markers which might have been switched (with threshold
#' \code{recombination} are removed. See also \code{\link[qtl]{checkAlleles}}).}
#' \item{The map is reestimated based on the observed markers.}
#' \item{Individuals with a fraction of crossovers higher than
#' \code{crossover} are removed.}
#' }
#' Steps 1, 3, 4, 5 and 7 are only performed if their respective threshold
#' values are positive. Setting them to 0 suppresses the corresponding check.\cr
#' Steps 2 and 6 are performed if respectively \code{removeDuplicates} and
#' \code{reestimateMap} are \code{TRUE}.
#'
#' @param cross An object of class cross created by the qtl package.
#' @param missMrk A numerical value between 0 and 1 indicating the maximum
#' allowed fraction of missing values per marker. Markers with a fraction of
#' missing values above \code{missMrk} will be removed.
#' @param missInd A numerical value between 0 and 1 indicating the maximum
#' allowed fraction of missing values per individual. Individuals with a fraction of
#' missing values above \code{missMrk} will be removed.
#' @param removeDuplicates Should duplicate markers be removed?
#' @param segDistortion A numerical value between 0 and 1 used a threshold for
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
#' @references Broman et al. (2003) R/qtl: QTL mapping in experimental crosses.
#' Bioinformatics 19:889-890
#'
#' @examples
#' ## Read the data.
#' F2 <- qtl::read.cross(format="csv",
#'                       file = system.file("extdata", "F2maize_geno.csv",
#'                                          package = "RAP"),
#'                       genotypes = c("AA", "AB", "BB"),
#'                       alleles = c("A", "B"), estimate.map = FALSE)
#' ## Run quality control.
#' F2QC <- QTLMapQC(F2)
#' ## Compare cross object before and after cleaning.
#' summary(F2)
#' summary(F2QC)
#'
#' ## Run quality control: only remove markers with a fraction of missing
#' ## values higher than 0.02
#' F2QC2 <- QTLMapQC(F2, missMrk = 0.02, missInd = 0, removeDuplicates = FALSE,
#'                  segDistortion = 0, recombination = 0, crossover = 0)
#' summary(F2QC2)
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
  ## Checks.
  if (missing(cross) || !inherits(cross, "cross")) {
    stop("cross should be an object of class cross\n")
  }
  if (!is.numeric(missMrk) || length(missMrk) > 1 || missMrk < 0 ||
      missMrk > 1) {
    stop("missMrk should be a single numerical value between 0 and 1.\n")
  }
  if (!is.numeric(missInd) || length(missInd) > 1 || missInd < 0 ||
      missInd > 1) {
    stop("missInd should be a single numerical value between 0 and 1.\n")
  }
  if (!is.numeric(recombination) || length(recombination) > 1 ||
      recombination < 0) {
    stop("recombination should be a single numerical value between 0 and 1.\n")
  }
  if (!is.numeric(crossover) || length(crossover) > 1 || crossover < 0 ||
      crossover > 1) {
    stop("crossover should be a single numerical value between 0 and 1.\n")
  }
  crossType <- setdiff(class(cross), "cross")
  if (!crossType %in% c("f2", "bc", "risib", "riself", "dh", "haploid") &&
      recombination > 0) {
    warning(paste0("Recombination is not possible for cross of type ",
                   crossType, ". Recombination set to 0.\n"), call. = FALSE)
    recombination <- 0
  }
  ## Extract number of individuals
  nInd <- qtl::nind(cross)
  if (missMrk > 0) {
    ## Calculate pct missing per marker and remove markers with pct above missMrk.
    pctMiss <- qtl::nmissing(cross, what = "mar") / nInd
    dropMrk <- names(pctMiss[pctMiss > missMrk])
    cross <- qtl::drop.markers(cross, markers = dropMrk)
  }
  if (removeDuplicates) {
    ## Remove duplicate markers.
    dupMar <- qtl::findDupMarkers(cross)
    cross <- qtl::drop.markers(cross, markers = unlist(dupMar))
  }
  nMrk <- sum(qtl::nmar(cross))
  if (missInd > 0) {
    ## Calculate pct missing per individual and remover individuals with pct
    ## above missInd.
    pctMiss <- qtl::nmissing(cross, what = "ind") / nMrk
    dropInd <- which(pctMiss > missInd)
    if (length(dropInd) > 0) {
      cross <- cross[, -dropInd]
    }
  }
  ## Compute the segregation distortion per marker and remove those showing
  ## evidence of distortion.
  if (segDistortion > 0) {
    segDist <- qtl::geno.table(cross)
    dropSegDist <- rownames(segDist[segDist$P.value < segDistortion, ])
    cross <- qtl::drop.markers(cross, markers = dropSegDist)
  }
  ## Estimate recombination frequencies between all pairs of markers.
  ## Suppress warning generated when there is recombination.
  crossRec <- supprWarn(qtl::est.rf(cross), "Alleles potentially switched")
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
#' A pdf report will be created containing a summary of a cross object.
#' Simultaneously the same report will be created as a tex file.
#'
#' @inheritParams report.AMMI
#'
#' @param x An object of class cross.
#'
#' @examples
#' ## Read the data
#' F2 <- qtl::read.cross(format="csv",
#'                       file = system.file("extdata", "F2maize_geno.csv",
 #'                                         package = "RAP"),
#'                       genotypes = c("AA", "AB", "BB"),
#'                       alleles = c("A", "B"), estimate.map = FALSE)
#' \dontrun{
#' ## Create a report
#' report(F2, outfile = "./testReports/reportCross.pdf")
#' }
#'
#' @export
report.cross <- function(x,
                         ...,
                         outfile = NULL) {
  createReport(x = x, reportName = "crossReport.Rnw", outfile = outfile, ...)
}





