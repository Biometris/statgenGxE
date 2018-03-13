#' Fit Single Trial Model using sommer
#'
#' Fit Single Trial Model using sommer
#'
#' @inheritParams STRunModel
#'
#' @seealso \code{\link{STRunModel}}
#'
#' @keywords internal
STModSommer <- function(TD,
                        traits,
                        what = c("fixed", "random"),
                        covariates = NULL,
                        useCheckId = FALSE,
                        control = NULL,
                        trySpatial = FALSE,
                        design = "rowcol",
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
  if (is.null(traits) || !is.character(traits) ||
      !all(traits %in% colnames(TD))) {
    stop("All traits have to be columns in TD.\n")
  }
  what <- match.arg(arg = what, several.ok = TRUE)
  if (!is.null(covariates) && (!is.character(covariates) ||
                               !(all(covariates %in% colnames(TD))))) {
    stop("covariates have to be columns in TD.\n")
  }
  for (colName in c(if (design %in% c("rowcol", "res.rowcol")) c("rowId",
                                                                 "colId"),
                    if (design %in% c("res.ibd", "res.rowcol", "rcbd")) "repId",
                    if (design %in% c("ibd", "res.ibd")) "subBlock",
                    if (useCheckId) "checkId")) {
    if (!is.null(colName) && (!is.character(colName) || length(colName) > 1 ||
                              !colName %in% colnames(TD))) {
      stop(paste(deparse(colName), "has to be NULL or a column in TD.\n"))
    }
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
                         paste(randEff, collapse =
                                 paste("+", if (useRepIdFix) "repId:")))
  } else {
    randomForm <- character()
  }
  if ("random" %in% what) {
    mr <- sapply(X = traits, FUN = function(trait) {
      ## Fit model with genotype random.
      sommer::mmer2(fixed = as.formula(paste(trait, fixedForm)),
                    random =  as.formula(paste("~", randomForm,
                                               if (length(randomForm) != 0) "+",
                                               "genotype")),
                    data = TD, silent = TRUE)
    }, simplify = FALSE)
  } else {
    mr <- NULL
  }
  if ("fixed" %in% what) {
    ## Fit model with genotype fixed.
    mf <- sapply(X = traits, FUN = function(trait) {
      sommer::mmer2(fixed = as.formula(paste(trait, fixedForm, "+ genotype")),
                    random = as.formula(paste("~", randomForm)),
                    data = TD, silent = TRUE)
    }, simplify = FALSE)
  } else {
    mf <- NULL
  }
  ## Construct SSA object.
  model <- createSSA(mRand = if ("random" %in% what) mr else NULL,
                     mFix = if ("fixed" %in% what) mf else NULL,
                     TD = TD, traits = traits,
                     design = design, engine = "sommer")
  return(model)
}
