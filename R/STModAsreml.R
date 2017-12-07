#' Fit Single Trial Model using asreml
#'
#' Fit Single Trial Model using asreml
#'
#' @inheritParams STRunModel
#'
#' @examples
#' ## Load data
#' data(TDHeat05)
#'
#' ## Fit model for row column design.
#' STModAsreml_1 <- STModAsreml(TD = TDHeat05, trait = "yield")
#'
#' ## Fit model for row column including replicates.
#' STModAsreml_2 <- STModAsreml(TD = TDHeat05, trait = "yield", design = "res.rowcol")
#'
#' ## Fit model for resolvable incomplete block design.
#' STModAsreml_3 <- STModAsreml(TD = TDHeat05, trait = "yield", design = "res.ibd")
#'
#' @export

STModAsreml <- function(TD,
                        traits,
                        covariates = NULL,
                        useCheckId = FALSE,
                        design = "rowcol",
                        trySpatial = FALSE,
                        ...) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  designs <- c("ibd", "res.ibd", "rcbd", "rowcol", "res.rowcol")
  if ((is.null(design) && !attr(TD, "design") %in% designs) ||
      (!is.null(design) && (!is.character(design) || length(design) > 1 ||
                            !design %in% designs))) {
    stop("design should either be an attribute of TD or one of 'ibd',
         'res.ibd', 'rcbd', 'rowcol' or 'res.rowcol'.\n")
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
  for (colName in c(if (design %in% c("rowcol", "res.rowcol")) c("rowId", "colId"),
                    if (design %in% c("res.ibd", "res.rowcol", "rcbd")) "repId",
                    if (design %in% c("ibd", "res.ibd")) "subBlock",
                    if (useCheckId) "checkId")) {
    if (!is.null(colName) && (!is.character(colName) || length(colName) > 1 ||
                              !colName %in% colnames(TD))) {
      stop(paste(deparse(colName), "has to be NULL or a column in data.\n"))
    }
  }
  if (is.character(trySpatial)) {
    stop("Spatial models not yet implemented for SpATS.\n")
  }
  ## Should repId be used as fixed effect in the model.
  useRepIdFix <- design %in% c("res.ibd", "res.rowcol", "rcbd")
  ## Indicate extra random effects.
  if (design %in% c("ibd", "res.ibd")) {
    randEff <- "subBlock"
  } else if (design %in% c("rowcol", "res.rowcol")) {
    randEff <- c("rowId", "colId")
  } else if (design == "rcbd") {
    randEff <- character()
  }
  ## Construct formula for fixed part.
  fixedForm <- paste("~",
                     if (useRepIdFix) "repId" else "1",
                     if (useCheckId) "+ checkId",
                     if (!is.null(covariates)) paste(c("", covariates),
                                                     collapse = "+"))
  ## Construct formula for random part. Include repId depending on design.
  if (length(randEff) != 0) {
    randomForm <- paste0(if (useRepIdFix) "repId:",
                         paste(randEff, collapse = paste("+", if (useRepIdFix) "repId:")))
  } else {
    randomForm <- character()
  }
  ## Create tempfile to suppress asreml output messages.
  tmp <- tempfile()
  sink(file = tmp)
  ## Create empty base lists.
  mr <- mf <- setNames(vector(mode = "list", length = length(traits)),
                       traits)
  for (trait in traits) {
    ## Fit model with genotype random.
    mrTrait <- asreml::asreml(fixed = as.formula(paste(trait, fixedForm)),
                         random = as.formula(paste("~", randomForm,
                                                   if (length(randomForm) != 0) "+",
                                                   "genotype")),
                         rcov = ~ units, aom = TRUE, data = TD, ...)
    ## Constrain variance of the variance components to be fixed as the values in mr.
    GParamTmp <- mrTrait$G.param
    for (randEf in randEff) {
      ## When there are no replicates the structure is [[randEf]][[randEf]]
      ## otherwise it is [[repId:randEf]][[repId]]
      GParamTmp[[paste0(ifelse(useRepIdFix, "repId:", ""),
                        randEf)]][[ifelse(useRepIdFix, "repId", randEf)]]$con <- "F"
    }
    ## Fit model with genotype fixed.
    if (length(randomForm) != 0) {
      mfTrait <- asreml::asreml(fixed = as.formula(paste(trait, fixedForm,
                                                         "+ genotype")),
                           random = as.formula(paste("~", randomForm)),
                           rcov = ~ units, G.param = GParamTmp, aom = TRUE,
                           data = TD, ...)
    } else {
      mfTrait <- asreml::asreml(fixed = as.formula(paste(trait, fixedForm,
                                                    "+ genotype")),
                           rcov = ~ units, G.param = GParamTmp, aom = TRUE,
                           data = TD, ...)
    }
    ## evaluate call terms in mr and mf so predict can be run.
    mrTrait$call$fixed <- eval(mrTrait$call$fixed)
    mrTrait$call$random <- eval(mrTrait$call$random)
    mrTrait$call$rcov <- eval(mrTrait$call$rcov)
    mfTrait$call$fixed <- eval(mfTrait$call$fixed)
    mfTrait$call$random <- eval(mfTrait$call$random)
    mfTrait$call$rcov <- eval(mfTrait$call$rcov)
    ## Run predict.
    if (useCheckId) {
      mfTrait <- predict(mfTrait, classify = "genotype", vcov = TRUE, data = TD,
                    associate = ~ checkId:genotype)
    } else {
      mfTrait <- predict(mfTrait, classify = "genotype", vcov = TRUE, data = TD)
    }
    mrTrait <- predict(mrTrait, classify = "genotype", data = TD)
    mfTrait$call$data <- substitute(TD)
    mrTrait$call$data <- substitute(TD)
    mr[[trait]] <- mrTrait
    mf[[trait]] <- mfTrait
  }
  sink()
  unlink(tmp)
  ## Construct SSA object.
  model <- createSSA(mMix = mr, mFix = mf, data = TD, traits = traits,
                     design = design, engine = "asreml")
  return(model)
}
