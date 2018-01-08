#' Fit Single Trial Mixed Model
#'
#' Perform REML analysis given a specific experimental design.
#' This is a wrapper function of \code{\link{STModSpATS}}, \code{\link{STModLme4}}
#' and \code{\link{STModAsreml}}. See details for the exact models fitted.
#' SpATS is used as a default method when design is rowcol or res.rowcol, lme4
#' for other designs.
#'
#' The actual model fitted depends on the design. For the supported designs the
#' following models are used:
#' \describe{
#' \item{ibd}{trait = genotype + \emph{subBlock} + e}
#' \item{res.ibd}{trait = genotype + \strong{repId} +
#' \emph{repId:subBlock} + e}
#' \item{rcbd}{trait = genotype + \strong{repId} + e}
#' \item{rowcol}{trait = genotype + \emph{rowId} + \emph{colId} + e}
#' \item{res.rowcol}{trait = genotype + \strong{repId} +
#' \emph{repId:rowId} + \emph{repId:colId} + e}
#' }
#' In the above models fixed effects are indicated in bold, random effects in
#' italics. genotype is fitted as fixed or random effect depending on
#' \code{what}.\cr
#' In case \code{useCheckId = TRUE} an extra fixed effect \strong{checkId} is included
#' in the model.\cr
#' Variables in \code{covariates} are fitted as extra fixed effects.\cr\cr
#' When \code{SpATS} is used for modelling an extra spatial term is included in the
#' model. This term is constructed using \code{\link[SpATS]{PSANOVA}} as
#' \code{PSANOVA(colCoordinates, rowCoordinates, nseg = nSeg, nest.div = 2)} where
#' \code{nSeg = (number of columns / 2, number of rows / 2)}. nseg and nest.div can
#' be modified using the \code{control} parameter.
#'
#' @param TD An object of class \code{\link{TD}}.
#' @param design A string specifying the experimental design. One of "ibd"
#' (incomplete block design), "res.ibd" (resolvable incomplete block design),
#' "rcbd" (randomized complete block design), "rowcol" (row column design) or
#' "res.rowcol" (resolvable row column design).
#' @param traits A character vector specifying the selected traits.
#' @param what A character vector specififying whether "genotype" should
#' be fitted as "fixed" or "random". If not specified both models
#' are fitted.
#' @param covariates A string specifying (a) covariate name(s).
#' @param useCheckId Should a checkId be used as a fixed parameter in the model?\cr
#' If \code{TRUE} \code{TD} has to contain a column 'checkId'.
#' @param trySpatial Should spatial models be tried? Spatials models are can only be
#' fitted with SpATS and asreml.
#' @param engine A string specifying the name of the mixed modelling engine to use,
#' either SpATS, lme4 or asreml. For spatial models SpaTS is used as a default, for
#' other models lme4.
#' @param control An optional list with control parameters to be passed to the actual
#' fitting funcions. For now only nSeg and nestDiv are valid and pass a value to
#' nseg and nest.div in \code{\link[SpATS]{PSANOVA}} respectively.
#' @param ... Further arguments to be passed to \code{SpATS}, \code{lme4} or \code{asreml}.
#'
#' @return an object of class \code{\link{SSA}}.
#'
#' @seealso \code{\link{STModSpATS}}, \code{\link{STModLme4}}, \code{\link{STModAsreml}}
#'
#' @references
#' Maria Xose Rodriguez-Alvarez, Martin P. Boer, Fred A. van Eeuwijk, Paul H.C. Eilers (2017).
#' Correcting for spatial heterogeneity in plant breeding experiments with P-splines. Spatial
#' Statistics \url{https://doi.org/10.1016/j.spasta.2017.10.003}\cr
#' Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting Linear
#' Mixed-Effects Models Using lme4. Journal of Statistical Software, 67(1), 1-48.
#' \url{https::/doi:10.18637/jss.v067.i01}.
#'
#' @examples
#' ## Fit model using lme4.
#' myModel1 <- STRunModel(TD = TDHeat05, design = "ibd", traits = "yield", what = "fixed")
#' summary(myModel1)
#' \dontrun{
#' report(myModel1, outfile = "./testReports/reportModelLme4.pdf")
#' }
#'
#' ## Fit model using SpATS.
#' myModel2 <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield",
#'                        what = "fixed")
#' summary(myModel2)
#' \dontrun{
#' report(myModel2, outfile = "./testReports/reportModelSpATS.pdf")
#' }
#'
#' #' ## Fit model using asreml.
#' \dontrun{
#' myModel3 <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield",
#'                        what = "fixed", engine = "asreml")
#' summary(myModel3)
#' report(myModel3, outfile = "./testReports/reportModelAsreml.pdf")
#' }
#'
#' @export
STRunModel = function(TD,
                      design = NULL,
                      traits,
                      what = c("fixed", "random"),
                      covariates = NULL,
                      useCheckId = FALSE,
                      trySpatial = FALSE,
                      engine = NA,
                      control = NULL,
                      ...) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  designs <- c("ibd", "res.ibd", "rcbd", "rowcol", "res.rowcol")
  if ((is.null(design) && (is.null(attr(TD, "design")) ||
                           !attr(TD, "design") %in% designs)) ||
      (!is.null(design) && (!is.character(design) || length(design) > 1 ||
                            !design %in% designs))) {
    stop("design should either be an attribute of TD or one of ibd,
         res.ibd, rcbd, rowcol or res.rowcol.\n")
  }
  ## Extract design from TD if needed.
  if (is.null(design)) {
    design <- attr(TD, "design")
  }
  if (is.null(traits) || !is.character(traits) || !all(traits %in% colnames(TD))) {
    stop("All traits have to be columns in TD.\n")
  }
  what <- match.arg(arg = what, several.ok = TRUE)
  if (!is.null(covariates) && (!is.character(covariates) ||
                               !(all(covariates %in% colnames(TD))))) {
    stop("covariates have to be a columns in TD.\n")
  }
  if (!is.logical(trySpatial) || length(trySpatial) > 1) {
    stop("trySpatial should be a single logical value.\n")
  }
  if (!is.na(engine) && (!is.character(engine) || length(engine) > 1 ||
                         !engine %in% c("asreml", "lme4", "SpATS"))) {
    stop("engine should be SpATS, asreml, lme4.\n")
  }
  if (is.na(engine)) {
    if (design %in% c("rowcol", "res.rowcol")) {
      message("Using SpATS for fitting models.")
      engine = "SpATS"
    } else {
      message("Using lme4 for fitting models.")
      engine = "lme4"
    }
  }
  if (trySpatial && engine == "lme4") {
    warning("Spatial models can only be fitted using SpATS or asreml.\n
            Defaulting to SpATS.")
    engine <- "SpATS"
  }
  ## Columns needed depend on design.
  for (colName in c(if (engine == "SpATS") c("rowCoordinates", "colCoordinates"),
                    if (design %in% c("rowcol", "res.rowcol")) c("rowId", "colId"),
                    if (design %in% c("res.ibd", "res.rowcol", "rcbd")) "repId",
                    if (design %in% c("ibd", "res.ibd")) "subBlock",
                    if (useCheckId) "checkId")) {
    if (!is.null(colName) && (!is.character(colName) || length(colName) > 1 ||
                              !colName %in% colnames(TD))) {
      stop(paste(deparse(colName), "has to be NULL or a column in TD.\n"))
    }
  }
  if (!is.null(control) && !is.list(control)) {
    stop("control has to be null or a list.\n")
  }
  ## Run model depending on engine.
  model <- do.call(what = paste0("STMod", tools::toTitleCase(engine)),
                   args = list(TD = TD, traits = traits,
                               what = what,
                               covariates = covariates,
                               useCheckId = useCheckId,
                               trySpatial = trySpatial,
                               design = design,
                               control = control,
                               ... = ...))
  return(model)
}
