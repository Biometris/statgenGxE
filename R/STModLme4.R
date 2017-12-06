#' Fit Single Trial Model using lme4
#'
#' Fit Single Trial Model using lme4
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

STModLme4 <- function(TD,
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
    randomForm <- paste0("(1 | ", if (useRepIdFix) "repId:",
                         paste(randEff,
                               collapse = paste(") + (1 | ", if (useRepIdFix) "repId:")),
                         ")")
  } else {
    randomForm <- character()
  }
  ## Fit model with genotype random.
  mr <- lme4::lmer(as.formula(paste(fixedForm,
                                    "+ (1 | genotype) ",
                                    if (length(randomForm) != 0) paste("+", randomForm))),
                   data = TD, na.action = na.exclude, ...)
  ## Fit model with genotype fixed.
  ## lme4 cannot handle models without random effect so in that case lm is called.
  if (length(randomForm) != 0) {
    mf <- lme4::lmer(as.formula(paste(fixedForm,
                                      "+ genotype + ", randomForm)),
                     data = TD, na.action = na.exclude, ...)
  } else  {
    mf <- lm(as.formula(paste(fixedForm, "+ genotype")),
             data = TD, na.action = na.exclude, ...)
  }
  ## Construct SSA object.
  model <- createSSA(mMix = mr, mFix = mf, data = TD, trait = trait,
                     design = design, engine = "lme4")
  return(model)
}
