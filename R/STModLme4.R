#' Fit Single Trial Model using lme4
#'
#' Fit Single Trial Model using lme4
#'
#' @inheritParams STRunModel
#'
#' @seealso \code{\link{STRunModel}}
#'
#' @keywords internal
STModLme4 <- function(TD,
                      trial = NULL,
                      traits,
                      what = c("fixed", "random"),
                      covariates = NULL,
                      useCheckId = FALSE,
                      control = NULL,
                      trySpatial = FALSE,
                      design = "rowcol",
                      checks = TRUE,
                      ...) {
  ## Base check.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (checks) {
    ## Checks.
    checkOut <- modelChecks(TD = TD, trial = trial, design = design,
                            traits = traits, what = what,
                            covariates = covariates, trySpatial = trySpatial,
                            engine = "lme4", useCheckId = useCheckId,
                            control = control)
    ## Convert output to variables.
    list2env(x = checkOut, envir = environment())
  }
  TDTr <- TD[[trial]]
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
    randomForm <- paste0("(1 | ", if (useRepIdFix) "repId:",
                         paste(randEff,
                               collapse = paste(") + (1 | ", if (useRepIdFix) "repId:")),
                         ")")
  } else {
    randomForm <- character()
  }
  if ("random" %in% what) {
    mr <- sapply(X = traits, FUN = function(trait) {
      ## Fit model with genotype random.
      lme4::lmer(as.formula(paste(trait, fixedForm,
                                  "+ (1 | genotype) ",
                                  if (length(randomForm) != 0) paste("+", randomForm))),
                 data = TDTr, na.action = na.exclude, ...)
    }, simplify = FALSE)
  } else {
    mr <- NULL
  }
  if ("fixed" %in% what) {
    ## Fit model with genotype fixed.
    ## lme4 cannot handle models without random effect so in that case lm is called.
    mf <- sapply(X = traits, FUN = function(trait) {
      if (length(randomForm) != 0) {
        lme4::lmer(as.formula(paste(trait, fixedForm,
                                    "+ genotype + ", randomForm)),
                   data = TDTr, na.action = na.exclude, ...)
      } else  {
        lm(as.formula(paste(trait, fixedForm, "+ genotype")),
           data = TDTr, na.action = na.exclude, ...)
      }}, simplify = FALSE)
  } else {
    mf <- NULL
  }
  spatial <- setNames(rep(FALSE, times = length(traits)), traits)
  ## Construct SSA object.
  return(list(mRand = if ("random" %in% what) mr else NULL,
              mFix = if ("fixed" %in% what) mf else NULL, TD = TD[trial],
              traits = traits, design = design, spatial = spatial,
              engine = "lme4", predicted = "genotype"))
}
