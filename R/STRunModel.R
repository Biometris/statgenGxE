#' Fit Single Trial Mixed Model
#'
#' Perform REML analysis given a specific experimental design.
#' This is a wrapper function of \code{\link{STModSpATS}}, \code{\link{STModLme4}}
#' and \code{\link{STModAsreml}}.
#'
#' @param TD An object of class TD.
#' @param design A string specifying the experimental design. One of "ibd"
#' (incomplete block design), "res.ibd" (resolvable incomplete block design),
#' "rcbd" (randomized complete block design), "rowcol" (row column design) or
#' "res.rowcol" (resolvable row column design).
#' @param traits A character vector specifying the selected traits.
#' @param covariates A string specifying (a) covariate name(s).
#' @param useCheckId Should a checkId be used as a fixed parameter in the model?\cr
#' If \code{TRUE} \code{TD} has to contain a column 'checkId'.
#' @param trySpatial Whether to try spatial models ("always", "ifregular"); default no spatial
#' models, i.e., \code{FALSE}.
#' @param engine A string specifying the name of the mixed modelling engine to use,
#' either SpATS, lme4 or asreml. For spatial models SpaTS is used as a default, for
#' other models lme4.
#' @param ... Further arguments to be passed to \code{SpATS}, \code{lme4} or \code{asreml}.
#'
#' @return an object of class \code{\link{SSA}}.
#'
#' @seealso \code{\link{STModSpATS}}, \code{\link{STModLme4}}, \code{\link{STModAsreml}}
#'
#' @examples
#' ## Load data.
#' data(TDHeat05)
#'
#' ## Fit model using lme4.
#' myModel1 <- STRunModel(TD = TDHeat05, design = "res.rowcol", trait = "yield")
#'
#' ## Fit model using SpATS.
#' myModel2 <- STRunModel(TD = TDHeat05, design = "res.rowcol", trait = "yield",
#'                       engine = "SpATS")
#'
#' @export
STRunModel = function(TD,
                      design = NULL,
                      traits,
                      covariates = NULL,
                      useCheckId = FALSE,
                      trySpatial = FALSE,
                      engine = if (!is.logical(trySpatial)) "SpATS" else "lme4",
                      control = NULL,
                      ...) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  designs <- c("ibd", "res.ibd", "rcbd", "rowcol", "res.rowcol")
  if ((is.null(design) && !attr(TD, "design") %in% designs) ||
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
  if (!is.null(covariates) && (!is.character(covariates) ||
                               !(all(covariates %in% colnames(TD))))) {
    stop("covariates have to be a columns in TD.\n")
  }
  if ((is.logical(trySpatial) && trySpatial) ||
      (is.character(trySpatial) && (length(trySpatial) > 1 ||
                                    !trySpatial %in% c("always", "ifregular")))) {
    stop("trySpatial should be NULL, always or ifregular.\n")
  }
  if (!is.null(engine) && (!is.character(engine) || length(engine) > 1 ||
                           !engine %in% c("asreml", "lme4", "SpATS"))) {
    stop("engine should be SpATS, asreml, lme4.\n")
  }
  if (is.character(trySpatial) && engine == "lme4") {
    warning("Spatial models can only be fitted using SpATS or asreml.\n
            Defaulting to SpATS.")
    engine <- "SpATS"
  }
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
                               covariates = covariates,
                               useCheckId = useCheckId,
                               design = design,
                               control = control,
                               ... = ...))
  return(model)
}
