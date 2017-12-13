#' Fit Single Trial Model using asreml
#'
#' Fit Single Trial Model using asreml
#'
#' @inheritParams STRunModel
#'
#' @seealso \code{\link{STRunModel}}
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
                        what = c("fixed", "random"),
                        covariates = NULL,
                        useCheckId = FALSE,
                        design = "rowcol",
                        trySpatial = FALSE,
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
  what <- match.arg(arg = what, choices = c("fixed", "random"), several.ok = TRUE)
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
  ## Create empty base lists.
  mr <- mf <- setNames(vector(mode = "list", length = length(traits)),
                       traits)
  ## Create tempfile to suppress asreml output messages.
  tmp <- tempfile()
  sink(file = tmp)
  for (trait in traits) {
    if ("random" %in% what) {
      ## Fit model with genotype random.
      mrTrait <- asreml::asreml(fixed = as.formula(paste(trait, fixedForm)),
                                random = as.formula(paste("~", randomForm,
                                                          if (length(randomForm) != 0) "+",
                                                          "genotype")),
                                rcov = ~ units, aom = TRUE, data = TD, ...)
      if ("fixed" %in% what) {
        ## Constrain variance of the variance components to be fixed as the values in mr.
        GParamTmp <- mrTrait$G.param
        for (randEf in randEff) {
          ## When there are no replicates the structure is [[randEf]][[randEf]]
          ## otherwise it is [[repId:randEf]][[repId]]
          GParamTmp[[paste0(ifelse(useRepIdFix, "repId:", ""),
                            randEf)]][[ifelse(useRepIdFix, "repId", randEf)]]$con <- "F"
        }
      }
      ## evaluate call terms in mr and mf so predict can be run.
      mrTrait$call$fixed <- eval(mrTrait$call$fixed)
      mrTrait$call$random <- eval(mrTrait$call$random)
      mrTrait$call$rcov <- eval(mrTrait$call$rcov)
      # Run predict.
      ## Asreml has a bug that may throw a warning message:
      ## Abnormal termination
      ## Insufficient workspace - (reset workspace or pworkspace arguments)
      ## This may be avoided by increasing pworkspace, but this doesn't
      ## always work.
      ## If this happens pworkspace is increased in 'small' steps.
      mrTraitP <- tryCatchExt(predict(mrTrait, classify = "genotype",
                                      vcov = TRUE, data = TD))
      pWorkSpace <- 8e6
      while (!is.null(mrTraitP$warning) && pWorkSpace < 10e6) {
        pWorkSpace <- pWorkSpace + 8e6
        mrTraitP <- tryCatchExt(predict(mrTrait, classify = "genotype",
                                        vcov = TRUE, data = TD,
                                        pworkspace = pWorkSpace))
      }
      if (is.null(mrTraitP$warning) && is.null(mrTraitP$error)) {
        mrTrait <- mrTraitP$value
      } else {
        stop(paste("Problem in asreml when running predict. Asreml message:\n",
                   mrTraitP$error$message, "\n",
                   mrTraitP$warning$message, "\n"))
      }
      mrTrait$call$data <- substitute(TD)
      mr[[trait]] <- mrTrait
    }
    if ("fixed" %in% what) {
      ## Fit model with genotype fixed.
      if (!"random" %in% what) {
        GParamTmp <- NULL
      }
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
      mfTrait$call$fixed <- eval(mfTrait$call$fixed)
      mfTrait$call$random <- eval(mfTrait$call$random)
      mfTrait$call$rcov <- eval(mfTrait$call$rcov)
      ## Construct assocForm for use in associate in predict.
      if (useCheckId) {
        assocForm <- as.formula("~ checkId:genotype")
      } else {
        assocForm <- as.formula("~ NULL")
      }
      ## Run predict.
      ## Asreml has a bug that may throw a warning message:
      ## Abnormal termination
      ## Insufficient workspace - (reset workspace or pworkspace arguments)
      ## This may be avoided by increasing pworkspace, but this doesn't
      ## always work.
      ## If this happens pworkspace is increased in 'small' steps.
      mfTraitP <- tryCatchExt(predict(mfTrait, classify = "genotype",
                                      associacte = assocForm, vcov = TRUE,
                                      data = TD))
      pWorkSpace <- 8e6
      while (!is.null(mfTraitP$warning) && pWorkSpace < 10e6) {
        pWorkSpace <- pWorkSpace + 8e6
        mfTraitP <- tryCatchExt(predict(mfTrait, classify = "genotype",
                                        associacte = assocForm, vcov = TRUE,
                                        data = TD, pworkspace = pWorkSpace))
      }
      if (is.null(mfTraitP$warning) && is.null(mfTraitP$error)) {
        mfTrait <- mfTraitP$value
      } else {
       stop(paste("Problem in asreml when running predict. Asreml message:\n",
                  mfTraitP$error$message, "\n",
                  mfTraitP$warning$message, "\n"))
      }
      mfTrait$call$data <- substitute(TD)
      mf[[trait]] <- mfTrait
    }
  }
  sink()
  unlink(tmp)
  ## Construct SSA object.
  model <- createSSA(mRand = if ("random" %in% what) mr else NULL,
                     mFix = if ("fixed" %in% what) mf else NULL,
                     data = TD, traits = traits,
                     design = design, engine = "asreml")
  return(model)
}
