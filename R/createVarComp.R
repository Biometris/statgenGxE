#' S3 class varComp
#'
#' Function for creating objects of S3 class varComp.\cr
#' \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}} and
#' \code{\link{report}} methods are available.
#'
#' @param fitMod A fitted variance components model.
#' @param modDat A data.frame containing the data used in fitting the model.
#'
#' @seealso \code{\link{plot.varComp}}, \code{\link{report.varComp}}
#'
#' @keywords internal
createVarComp <- function(fitMod,
                          modDat,
                          engine) {
  varComp <- structure(list(fitMod = fitMod,
                            modDat = modDat,
                            engine = engine),
                       class = "varComp")
  attr(varComp, which = "timestamp") <- Sys.time()
  return(varComp)
}

#' @export
print.varComp <- function(x, ...) {
  summary(object$fitmod)
}

#' @export
summary.varComp <- function(object, ...) {
  summary(object$fitMod)
}

#' Plot function for class varComp
#'
#' Plot function for class varComp.
#'
#' @param x An object of class varComp
#' @param ... Not used
#' @param output Should the plot be output to the current device? If
#' \code{FALSE} only a list of ggplot objects is invisibly returned.
#'
#' @examples
#' \dontrun{
#'
#' }
#'
#' @export
plot.varComp <- function(x,
                         ...,
                         output = TRUE) {
  NULL
}

#' Report method for class varComp
#'
#' A pdf report will be created containing a summary of an object of class
#' varComp. Simultaneously the same report will be created as a tex
#' file.
#'
#' @inheritParams report.AMMI
#'
#' @param x An object of class varComp.
#'
#' @return A pdf and tex report.
#'
#' @examples
#' \dontrun{

#' }
#' @export
report.varComp <- function(x,
                           ...,
                           outfile = NULL) {
  ## Checks.
  if (nchar(Sys.which("pdflatex")) == 0) {
    stop("An installation of LaTeX is required to create a pdf report.\n")
  }
  createReport(x = x, reportName = "varCompReport.Rnw", outfile = outfile,
               reportPackage = "statgenGxE", ...)
}

#' @export
vc <- function(varComp) {
  if (varComp$engine == "lme4") {
    varcomps <- as.data.frame(lme4::VarCorr(varComp$fitMod))
    rownames(varcomps) <- varcomps[["grp"]]
    varcomps <- varcomps[c((nrow(varcomps)-1):1, nrow(varcomps)),
                         "vcov", drop = FALSE]
    colnames(varcomps) <- "component"
  } else if (varComp$engine == "asreml") {
    varcomps <- summary(varComp$fitMod)$varcomp
    rownames(varcomps)[nrow(varcomps)] <- "Residual"
    varcomps <- varcomps[, "component", drop = FALSE]
  }
  return(varcomps)
}

#' @export
herit <- function(varComp) {
  fitMod <- varComp$fitMod
  modDat <- varComp$modDat
  varcomps <- vc(varComp)
  sigmaG <- varcomps["genotype", "component"]
  sigmaRes <- varcomps["Residual", "component"]
  numerator <- sigmaG
  modTerms <- rownames(varcomps)
  if (varComp$engine == "lme4") {
    modVars <- rownames(attr(x = terms(fitMod, random.only = TRUE),
                             which = "factors"))[-c(1, 2)]

  } else if (varComp$engine == "asreml") {
    modVars <- rownames(attr(x = terms(fitMod$call$random), which = "factors"))[-1]
  }
  for (term in modTerms[-c(1, length(modTerms))]) {
    sigmaTerm <- varcomps[term, "component"]
    termVars <- unlist(strsplit(x = term, split = ":"))[-1]
    numerator <- numerator + sigmaTerm /
      prod(sapply(X = termVars, FUN = function(termVar) {
        nlevels(modDat[[termVar]])}))
  }
  if (length(modVars) > 0) {
    numerator <- numerator + sigmaRes /
      prod(sapply(X = modVars, FUN = function(modVar) {
        nlevels(modDat[[modVar]])}))
  } else {
    numerator <- numerator + sigmaRes
  }
  return(sigmaG / numerator)
}
