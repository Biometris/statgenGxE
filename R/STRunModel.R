#' Run the single trial mixed model analysis
#'
#' Perform REML analysis given a specific experimental design.
#' This is a wrapper function of \code{\link{ST.mod.alpha}}, \code{\link{ST.mod.rcbd}}
#' and \code{\link{ST.mod.rowcol}}.
#'
#' @param TD An object of class TD.
#' @param design A string specifying the experimental design.
#' @param trait A string specifying the selected trait.
#' @param covariates A string specifying (a) covariate name(s).
#' @param useCheckId Should a checkId be used as a fixed parameter in the model?\cr
#' If \code{TRUE} \code{TD} has to contain a column 'checkId'.
#' @param trySpatial Whether to try spatial models ("always", "ifregular"); default no spatial
#' models, i.e., \code{NULL}.
#' @param ... Further arguments to be passed to either \code{asreml} or \code{lme4}.
#'
#' @return an object of class \code{\link{SSA}}.
#'
#' @seealso \code{\link{SSA}}, \code{\link{summary.SSA}}
#'
#' @examples
#' data(TDHeat05)
#' myModel <- STRunModel(TD = TDHeat05, design = "res.rowcol", trait = "yield")
#' names(myModel)
#'
#' @export
STRunModel = function(TD,
                      design = NULL,
                      trait,
                      covariates = NULL,
                      useCheckId = FALSE,
                      trySpatial = FALSE,
                      engine = if (trySpatial) "SpATS" else "lme4",
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
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TD)) {
    stop("trait has to be a column in TD.\n")
  }
  if (!is.null(covariates) && (!is.character(covariates) ||
                               !(all(covariates %in% colnames(TD))))) {
    stop("covariates have to be a columns in TD.\n")
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
  # if (trySpatial || (!trySpatial && (!is.character(trySpatial) || length(trySpatial) > 1 ||
  #                              !trySpatial %in% c("always", "ifregular")))) {
  #   stop("trySpatial should be NULL, always or ifregular")
  # }
  ## Extract design from TD if needed.
  if (is.null(design)) {
    design <- attr(TD, "design")
  }
  ## Run model depending on engine.
  model <- do.call(what = paste0("STMod", tools::toTitleCase(engine)),
                   args = list(TD = TD, trait = trait,
                               covariates = covariates,
                               useCheckId = useCheckId,
                               design = design))
  return(model)
}
