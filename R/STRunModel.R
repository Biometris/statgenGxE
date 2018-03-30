#' Fit Single Trial Mixed Model
#'
#' Perform REML analysis given a specific experimental design. This is a wrapper
#' function of \code{\link{STModSpATS}}, \code{\link{STModLme4}} and
#' \code{\link{STModAsreml}}. See details for the exact models fitted. SpATS is
#' used as a default method when design is rowcol or res.rowcol, lme4 for other
#' designs.
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
#' In the above models fixed effects are indicated in \strong{bold}, random
#' effects in \emph{italics}. genotype is fitted as fixed or random effect
#' depending on the value of \code{what}.\cr
#' In case \code{useCheckId = TRUE} an extra fixed effect \strong{checkId} is
#' included in the model.\cr
#' Variables in \code{covariates} are fitted as extra fixed effects.\cr\cr
#' When \code{SpATS} is used for modelling an extra spatial term is included
#' in the model. This term is constructed using the function
#' \code{\link[SpATS]{PSANOVA}} from the SpATS package as
#' \code{PSANOVA(colCoordinates, rowCoordinates, nseg = nSeg, nest.div = 2)}
#' where \code{nSeg = (number of columns / 2, number of rows / 2)}. nseg and
#' nest.div can be modified using the \code{control} parameter.\cr\cr
#' When \code{asreml} is used for modeling and \code{trySpatial} is \code{TRUE}
#' 6 models are fitted with different random term and covariance structure.
#' The best model is determined based on a goodness-of-fit criterion, either
#' AIC or BIC. This can be set using the control parameter \code{criterion},
#' default is AIC.
#' The fitted random terms depend on the structure of the data. If the design
#' has a regular structure, i.e. all replicates appear the same amount of times
#' in the design, the following combinations of random and spatial terms are
#' fitted
#' \itemize{
#' \item{random = NULL, spatial = exp(rowCoordinates):colCoordinates}
#' \item{random = NULL, spatial = rowCoordinates:exp(colCoordinates)}
#' \item{random = NULL, spatial = iexp(rowCoordinates,colCoordinates)}
#' \item{random = repId:rowId, spatial = exp(rowCoordinates):colCoordinates}
#' \item{random = repId:colId, spatial = rowCoordinates:exp(colCoordinates)}
#' \item{random = repId:rowId + repId:colId,
#' spatial = iexp(rowCoordinates,colCoordinates)}
#' }
#' If the design is not regular the following following combinations of random
#' and spatial terms are fitted
#' \itemize{
#' \item{random = NULL, spatial = ar1(rowId):colId}
#' \item{random = NULL, spatial = rowId:ar1(colId)}
#' \item{random = NULL, spatial = ar1(rowId):ar1(colId)}
#' \item{random = repId:rowId, spatial = ar1(rowId):colId}
#' \item{random = repId:colId, spatial = rowId:ar1(colId)}
#' \item{random = repId:rowId + repId:colId, spatial = ar1(rowId):ar1(colId)}
#' }
#' If there are no replicates in the model, in the random parts above, repId is
#' left out.
#'
#' @param TD An object of class \code{\link{TD}}.
#' @param trials A character vector specifying the trials for modeling.
#' @param design A character string specifying the experimental design. Either
#' "ibd" (incomplete block design), "res.ibd" (resolvable incomplete block
#' design), "rcbd" (randomized complete block design), "rowcol" (row column
#' design) or "res.rowcol" (resolvable row column design).
#' @param traits A character vector specifying the traits for modeling.
#' @param what A character vector specifying whether "genotype" should
#' be fitted as "fixed" or "random" effect. If not specified both models
#' are fitted.
#' @param covariates A character vector specifying covariates to be fitted as
#' extra fixed effects in the model.
#' @param useCheckId Should checkId be used as a fixed effect in the model?\cr
#' If \code{TRUE} \code{TD} has to contain a column 'checkId'.
#' @param trySpatial Should spatial models be tried? Spatial models can
#' only be fitted with SpATS and asreml. If SpATS is used for modeling only
#' spatial models can be fitted and trySpatial is always set to \code{TRUE}. If
#' asreml is used fitting spatial models is optional.
#' @param engine A string specifying the name of the mixed modelling engine to
#' use, either SpATS, lme4 or asreml. For spatial models SpaTS is used as a
#' default, for other models lme4.
#' @param control An optional list with control parameters to be passed to the
#' actual fitting funcions. Currently \code{nSeg} and \code{nestDiv} are valid
#' parameters when fitting a model using SpATS. They pass a value to nseg and
#' nest.div in \code{\link[SpATS]{PSANOVA}} respectively. \code{criterion} is a
#' valid parameter when fitting a spatial model using asreml. Use this to pass
#' a goodness-of-fit criterion for comparing different spatial models. See also
#' in details. Other parameters are ignored.
#' @param ... Further arguments to be passed to \code{SpATS}, \code{lme4} or
#' \code{asreml}.
#'
#' @return An object of class \code{\link{SSA}}, a list containing per trial
#' that has been analyzed a list of:
#' \item{mRand}{A list of models with fitted with genotype as random effect.}
#' \item{mFix}{A list of models fitted with genotype as fixed effect.}
#' \item{TD}{An object of class \code{\link{TD}} containing the data on which
#' \code{mRand} and \code{mFix} are based.}
#' \item{traits}{A character vector indicating the traits for which the analysis
#' is done.}
#' \item{design}{A character string containing the design of the trial.
#' (see \code{\link{STRunModel}} for the possible designs).}
#' \item{spatial}{A character string indicating the spatial part of the model.
#' \code{FALSE} if no spatial design has been used.}
#' \item{engine}{A character string containing the engine used for the
#' analysis.}
#' \item{predicted}{A character string indicating the variable that has been
#' predicted.}
#'
#' @seealso \code{\link{STModSpATS}}, \code{\link{STModLme4}},
#' \code{\link{STModAsreml}}
#'
#' @references
#' Maria Xose Rodriguez-Alvarez, Martin P. Boer, Fred A. van Eeuwijk, Paul H.C.
#' Eilers (2017). Correcting for spatial heterogeneity in plant breeding
#' experiments with P-splines. Spatial Statistics
#' \url{https://doi.org/10.1016/j.spasta.2017.10.003}
#' @references
#' Butler, D. G., et al. (2010). Analysis of Mixed Models for S language
#' environments: ASReml-R reference manual. Brisbane, DPI Publications
#' @references
#' Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015). Fitting
#' Linear Mixed-Effects Models Using lme4. Journal of Statistical Software,
#' 67(1), 1-48. \url{https://www.jstatsoft.org/article/view/v067i01/0}.
#'
#' @examples
#' ## Fit model using lme4.
#' myModel1 <- STRunModel(TD = TDHeat05, design = "ibd", traits = "yield",
#'                       what = "fixed")
#' ## Summarize results.
#' summary(myModel1)
#' ## Create base plots of the results.
#' plot(myModel1)
#' \dontrun{
#' ## Create a pdf report summarizing results.
#' report(myModel1, outfile = "./testReports/reportModelLme4.pdf")
#' }
#' ## Fit model using SpATS.
#' myModel2 <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield",
#'                       what = "fixed")
#' summary(myModel2)
#' ## Create spatial plots of the results.
#' plot(myModel2, plotType = "spatial")
#' \dontrun{
#' report(myModel2, outfile = "./testReports/reportModelSpATS.pdf")
#' }
#'
#' ## Fit model using asreml.
#' \dontrun{
#' myModel3 <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield",
#'                       what = "fixed", engine = "asreml")
#' summary(myModel3)
#' report(myModel3, outfile = "./testReports/reportModelAsreml.pdf")
#' }
#' @export
STRunModel = function(TD,
                      trials = names(TD),
                      design = NULL,
                      traits,
                      what = c("fixed", "random"),
                      covariates = NULL,
                      useCheckId = FALSE,
                      trySpatial = FALSE,
                      engine = NA,
                      control = NULL,
                      ...) {
  ## Base check.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  ## Run models depending on engine.
  models <- sapply(X = trials, FUN = function(trial) {
    ## Checks.
    checkOut <- modelChecks(TD = TD, trial = trial, design = design,
                            traits = traits, what = what,
                            covariates = covariates, trySpatial = trySpatial,
                            engine = engine, useCheckId = useCheckId,
                            control = control)
    ## Convert output to variables.
    list2env(x = checkOut, envir = environment())
    model <- do.call(what = paste0("STMod", tools::toTitleCase(engine)),
                     args = list(TD = TD, trial = trial, traits = traits,
                                 what = what, covariates = covariates,
                                 useCheckId = useCheckId,
                                 trySpatial = trySpatial, design = design,
                                 control = control, checks = FALSE,
                                 ... = ...))
    return(model)
  }, simplify = FALSE)
  return(createSSA(models = models))
}

## Helper function for performing checks for single trial modeling.
modelChecks <- function(TD,
                        trial,
                        design,
                        traits,
                        what = c("fixed", "random"),
                        covariates,
                        trySpatial,
                        engine,
                        useCheckId,
                        control) {
  designs <- c("ibd", "res.ibd", "rcbd", "rowcol", "res.rowcol")
  engines <- c("SpATS", "lme4", "asreml")
  if (!is.character(trial) || !trial %in% names(TD)) {
    stop("trial should be in TD.\n")
  }
  if ((is.null(design) && (is.null(attr(TD[[trial]], "trDesign")) ||
                           !attr(TD[[trial]], "trDesign") %in% designs)) ||
      (!is.null(design) && (!is.character(design) || length(design) > 1 ||
                            !design %in% designs))) {
    stop(paste("design should either be an attribute of TD or one of ",
               paste(designs, collapse = ", "), ".\n"))
  }
  ## Extract design from TD if needed.
  if (is.null(design)) {
    design <- attr(TD[[trial]], "trDesign")
  }
  if (is.null(traits) || !is.character(traits)) {
    stop("Traits should be a character vector.\n")
  }
  if (!all(traits %in% colnames(TD[[trial]]))) {
    stop(paste0("All traits should be columns in ", trial, ".\n"))
  }
  what <- match.arg(arg = what, several.ok = TRUE)
  if (!is.null(covariates) && !is.character(covariates)) {
    stop("Covariates should be NULL or a character vector.\n")
  }
  if (!all(covariates %in% colnames(TD[[trial]]))) {
    stop(paste0("All covariates should be columns in ", trial, ".\n"))
  }
  if (!is.logical(trySpatial) || length(trySpatial) > 1) {
    stop("trySpatial should be a single logical value.\n")
  }
  if (!is.na(engine) && (!is.character(engine) || length(engine) > 1 ||
                         !engine %in% engines)) {
    stop(paste0("engine should be one of ", paste(engines, collapse = ", "),
                ".\n"))
  }
  if (is.na(engine)) {
    if (design %in% c("rowcol", "res.rowcol")) {
      message("Using SpATS for fitting models.")
      engine <- "SpATS"
    } else {
      message("Using lme4 for fitting models.")
      engine <- "lme4"
    }
  }
  if (trySpatial && engine == "lme4") {
    warning("Spatial models can only be fitted using SpATS or asreml.\n
              Defaulting to SpATS.", call. = FALSE)
    engine <- "SpATS"
  }
  ## Columns needed depend on design.
  desCols <- c(if (engine == "SpATS") c("rowCoordinates", "colCoordinates"),
               if (design %in% c("rowcol", "res.rowcol")) c("rowId", "colId"),
               if (design %in% c("res.ibd", "res.rowcol", "rcbd")) "repId",
               if (design %in% c("ibd", "res.ibd")) "subBlock",
               if (useCheckId) "checkId")
  for (desCol in desCols) {
    if (!desCol %in% colnames(TD[[trial]])) {
      stop(paste0(desCol, " should be a column in ", trial, ".\n"))
    }
  }
  if (!is.null(control) && !is.list(control)) {
    stop("control has to be NULL or a list.\n")
  }
  return(list(design = design, what = what, engine = engine))
}


