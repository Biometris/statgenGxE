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
#' STModLme4_1 <- STModLme4(TD = TDHeat05, trait = "yield")
#'
#' ## Fit model for row column including replicates.
#' STModLme4_2 <- STModLme4(TD = TDHeat05, trait = "yield", design = "res.rowcol")
#'
#' ## Fit model for resolvable incomplete block design.
#' STModLme4_3 <- STModLme4(TD = TDHeat05, trait = "yield", design = "res.ibd")
#'
#' @export

STModAsreml <- function(TD,
                        trait,
                        covariates = NULL,
                        useCheckId = FALSE,
                        design = "rowcol",
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
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TD)) {
    stop("trait has to be a column in TD.\n")
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
  ## Extract design from TD if needed.
  if (is.null(design)) {
    design <- attr(TD, "design")
  }
  ## Should repId be used as fixed or random effect in the model.
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
  fixedForm <- paste(trait, "~",
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
  ## Fit model with genotype random.
  mr <- asreml::asreml(fixed = as.formula(fixedForm),
                       random = as.formula(paste("~", randomForm,
                                                 if (length(randomForm) != 0) "+",
                                                 "genotype")),
                       rcov = ~ units, aom = TRUE, data = TD, ...)
  ## Constrain variance of the variance components to be fixed as the values in mr.
  GParamTmp <- mr$G.param
  for (randEf in randEff) {
    ## When there are no replicates the structure is [[randEf]][[randEf]]
    ## otherwise it is [[repId:randEf]][[repId]]
    GParamTmp[[paste0(ifelse(useRepIdFix, "repId:", ""),
                      randEf)]][[ifelse(useRepIdFix, "repId", randEf)]]$con <- "F"
  }
  ## Fit model with genotype fixed.
  if (length(randomForm) != 0) {
    mf <- asreml::asreml(fixed = as.formula(paste(fixedForm, "+ genotype")),
                         random = as.formula(paste("~", randomForm)),
                         rcov = ~ units, G.param = GParamTmp, aom = TRUE,
                         data = TD, ...)
  } else {
    mf <- asreml::asreml(fixed = as.formula(paste(fixedForm, "+ genotype")),
                         rcov = ~ units, G.param = GParamTmp, aom = TRUE,
                         data = TD, ...)
  }
  ## evaluate call terms in mr and mf so predict can be run.
  mr$call$fixed <- eval(mr$call$fixed)
  mr$call$random <- eval(mr$call$random)
  mr$call$rcov <- eval(mr$call$rcov)
  mf$call$fixed <- eval(mf$call$fixed)
  mf$call$random <- eval(mf$call$random)
  mf$call$rcov <- eval(mf$call$rcov)
  ## Run predict.
  if (useCheckId) {
    mf <- predict(mf, classify = "genotype", vcov = TRUE, data = TD,
                  associate = ~ checkId:genotype)
  } else {
    mf <- predict(mf, classify = "genotype", vcov = TRUE, data = TD)
  }
  mr <- predict(mr, classify = "genotype", data = TD)
  sink()
  unlink(tmp)
  mf$call$data <- substitute(TD)
  mr$call$data <- substitute(TD)
  ## Construct SSA object.
  model <- createSSA(mMix = mr, mFix = mf, data = TD, trait = trait,
                     design = design, engine = "asreml")
  return(model)
}
